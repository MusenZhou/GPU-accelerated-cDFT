cc************************************************************
cc* THIS SUBROUTINE IS TO CALCULATE THE INTRINSIC POTENTIAL **
CC* USING FFT METHOD.                                       **
CC*                             S.L. ZHAO                   **
CC************************************************************

      SUBROUTINE Compute_Vint

      use, intrinsic::iso_c_binding
      use rho
      use vext
      use pe
      include 'param.h'
      include 'fft.h'
      include 'system.h'
      include 'numbers.h'
      include 'c.h'
      include 'rho.h'
      include 'Vext.h'
      include 'iteration.h'
      include 'pe.h'



*-------------  VARIABLES LOCALES  ------------------------------
      INTEGER i, j, k, l, n, m

      REAL*8 x_com, y_com, z_com, r_nm,rr

      integer nk, k_index, error
      integer m1, m2, m3, ml, mm, mn 
      real*8 x1, y1, z1, x2, y2, z2
      real*8 nr, norm_k, kx, ky, kz, k2

**energy

      complex(C_DOUBLE_COMPLEX),allocatable:: delta_rhoO(:,:,:) ! 2, nfft1, nfft2, nfft3
      complex(C_DOUBLE_COMPLEX),allocatable:: delta_rhoH(:,:,:) ! 2, nfft1, nfft2, nfft3
c      real*8,allocatable:: ele_rho(:,:,:,:) ! 2, nfft1, nfft2, nfft3
      complex(C_DOUBLE_COMPLEX),allocatable:: VkOO(:,:,:)
      complex(C_DOUBLE_COMPLEX),allocatable:: VkOH(:,:,:) ! 2, nfft1, nfft2, nfft3
c      complex(C_DOUBLE_COMPLEX),allocatable:: ele_temp(:,:,:)
      real*8 ckOO_grid, ckOH_grid, ckHH_grid
c      real*8 ckOO_grid_r, ckOH_grid_r, ckHH_grid_r
      

      real*8 DeltaVk
      real*8 rhoO, rhoH
c interpolation
      real*8 GG
      integer n1, n2, n3

*** for output
      CHARACTER*40 file_O, file_H
      CHARACTER NAM*7, CH*3
      CHARACTER*11 chemin

      real*8 ele_cut,ele_n,faitemp

      real*8 NumO, NumH,NumOH,ep1,ep2,tem_rhoO,tem_rhoH 
      real*8 NumH2,NumO2

      real*8 w_pe

CC**************** PART ONE : INPUT density AND direct correlation **************
c   Compute angular grid

      ALLOCATE(delta_rhoO(nfft1,nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION delta_rhoO FAIL!"

      ALLOCATE(delta_rhoH(nfft1,nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION delta_rhoH FAIL!"
     
c      ALLOCATE(ele_temp(nfft1,nfft2,nfft3), stat = error)
c      IF(ERROR.NE.0) STOP "ALLOCATION delta_rhoH FAIL!"
c      ALLOCATE(ele_rho(2,nfft1,nfft2,nfft3), stat = error)
c      IF(ERROR.NE.0) STOP "ALLOCATION delta_rhoH FAIL!"

       ele_cut=lg/2.0-4.0
       ele_tot=0.0D0
       ele_n=0.0D0


       DO  k=1,nfft3
         DO  j=1,nfft2
           DO  i=1,nfft1
            x_com = (i-1)*Lg/FLOAT(nfft1)-lg/2.0
            y_com = (j-1)*Lg/FLOAT(nfft2)-lg/2.0
            z_com = (k-1)*Lg/FLOAT(nfft3)-lg/2.0

            rr = DSQRT(x_com**2 + y_com**2 + z_com**2)
            if(rr.lt.ele_cut)then
              ele_n=ele_n+1.0D0
            endif

c** prepare delta_rhoO and delta_rhoH for Fourier transform purpose

            delta_rhoO(i,j,k)=cmplx(rho_O_old(i,j,k)-1.0d0,0.0)   
            delta_rhoH(i,j,k)=cmplx(rho_H_old(i,j,k)-1.0d0,0.0)  
            ele_rho(i,j,k)=rho_O_old(i,j,k)-rho_H_old(i,j,k)
c            ele_temp(i,j,k)=cmplx(ele_rho(i,j,k),0.0D0)
            ele_tot=ele_tot+ele_rho(i,j,k)

           END DO
          END DO
        END DO


        write(*,*)'rhoO-rhoH =',ele_tot


        Call Compute_ele_pot_3D

        call fftw_3d_for(nfft1,nfft2,nfft3,delta_rhoO,delta_rhoO)
        call fftw_3d_for(nfft1,nfft2,nfft3,delta_rhoH,delta_rhoH)
      
c        call fftw_3d_for(nfft1,nfft2,nfft3,ele_temp,ele_temp)
       
  
CC************************* PART TWO  *****************************

       ALLOCATE(VkOO(nfft1,nfft2,nfft3), stat = error)
       IF(ERROR.NE.0) STOP "ALLOCATION VkOO FAIL!"

       ALLOCATE(VkOH(nfft1,nfft2,nfft3), stat = error)
       IF(ERROR.NE.0) STOP "ALLOCATION VkOH FAIL!"

       IF(IVIN.eq.0) then
         write(*,*) "--- Jump over Vint Calculation --- "
         goto 888
       END IF


c      IF(nb_solute_sites.eq.1) then
c        Call Compute_ele_pot_1D   ! long-range correlation contribution
c      ELSE 
c        Call Compute_ele_pot_3D   ! long-range correlation contribution
c      END IF

c  Loop over k-vectors

      DO n=1, nfft3
         m3 = n
         if ( n .gt. nf3 )m3 = n - nfft3   

        DO m=1, nfft2
          m2 = m
         if ( m .gt. nf2 )m2 = m - nfft2

         DO l=1, nfft1
            m1 = l
           if ( l .gt. nf1 ) m1 = l -  nfft1


         kx=twopi*(m1-1)/Lg
         ky=twopi*(m2-1)/Lg
         kz=twopi*(m3-1)/Lg

         k2 = kx**2 + ky**2 + kz**2
         norm_k = sqrt(k2)

         k_index = int(norm_k/delta_ck)! + 1
         nr = norm_k/delta_ck - k_index

       if(norm_k.lt.0.02066D0)then
c        norm_k=0.02066D0
c        k2=norm_k**2
       endif

       IF(k_index.lt.nb_ck) THEN 

           ckOO_grid = ckOO(k_index)*(one-nr) + ckOO(k_index+1)*nr  !GOT IT
           ckOH_grid = ckOH(k_index)*(one-nr) + ckOH(k_index+1)*nr  !GOT IT
           ckHH_grid = ckHH(k_index)*(one-nr) + ckHH(k_index+1)*nr  !GOT IT

**** long range subtracted
          ckOO_grid =  ckOO_grid
     &      + 4985.964/(k2 + 1.0D0)
c     &         -4836.4/k2 
          ckOH_grid = ckOH_grid
     &      - 2492.982/(k2 + 1.0D0)
c     &         +2418.2/k2       
          ckHH_grid = ckHH_grid
     &      + 1246.491/(k2 + 1.0D0) 
c     &         -1209.1/k2 
       ELSE

**** long range subtracted 

           ckOO_grid =  zero
c     &      + 4836.4*1.0/(k2 + 1.44)
c     &         -4836.4/k2 
          ckOH_grid = zero
c     &      - 2418.2*1.0/(k2 + 1.44)
c     &         +2418.2/k2       
          ckHH_grid = zero
c     &      + 1209.1*1.0/(k2 + 1.44) 
c     &         -1209.1/k2 

       END IF




c   Compute intrinsic potential Vk(l,m,n) in k-space
c        if(norm_k.gt.0.002066D0)then

           VkOO(l,m,n) = ckOO_grid*delta_rhoO(l,m,n)
     &                   + 2.0*ckOH_grid*delta_rhoH(l,m,n)



           VkOH(l,m,n) = ckOH_grid*delta_rhoO(l,m,n) 
     &                   + 2.0*ckHH_grid*delta_rhoH(l,m,n)

c           ele_rho(l,m,n)=-ele_rho(l,m,n)*2418.2/(k2+1.0D-7)

           END DO
          END DO
         END DO   !  end loop over k-vectors


       DEALLOCATE(delta_rhoO)
       DEALLOCATE(delta_rhoH)


        call fftw_3d_bak(nfft1,nfft2,nfft3,VkOO,VkOO)
        call fftw_3d_bak(nfft1,nfft2,nfft3,VkOH,VkOH)

c        call fftw_3d_bak(nfft1,nfft2,nfft3,ele_temp,ele_temp)


888    continue

**** for output only
      
       

       sfe_tot=0.0D0
       detan=0.0D0

        DO k=1, nfft3
         DO j=1,nfft2
          DO i=1,nfft1
           x_com = (i-1)*Lg/FLOAT(nfft1)-lg/2.0
           y_com = (j-1)*Lg/FLOAT(nfft2)-lg/2.0
           z_com = (k-1)*Lg/FLOAT(nfft3)-lg/2.0

           rr = DSQRT(x_com**2 + y_com**2 + z_com**2)
           
    
c           w_pe=1.0D0/((rho_O_old(i,j,k)/1.2)**12+1.0D0)
c           fai_pe(i,j,k)=fai_pe(i,j,k)*w_pe
c     &                  +(1.0D0-w_pe)*ele_rho(1,i,j,k)*rho_0
c        IF(IVin.eq.0) ratio = 0.0d0 
c        IF(IVIn.eq.0) ratio = 0.0d0 


c           if((rr.gt.1.0D-4))then
c           ele_rho(1,i,j,k)=ele_rho(1,i,j,k)
c     &                 +ele_tot*2418.2/rr/4.0/pi*Deltav
c           endif

           VintO(i,j,k) = -REAL(VkOO(i,j,k))*rho_0
     &                  + beta*VextO(i,j,k)
c     &                  +zo*faitemp
     &                  + zo*fai_pe(i,j,k)  ! beta included in fai_pe
c     &                  + zo*REAL(ele_temp(i,j,k))*rho_0
           VintH(i,j,k) = -REAL(VkOH(i,j,k))*rho_0
     &                  + beta*VextH(i,j,k)
c     &                  +zh*faitemp
     &                  + zh*fai_pe(i,j,k)  ! beta included in fai_pe
c     &                  + zh*REAL(ele_temp(i,j,k))*rho_0



           sfe_tot=sfe_tot+(REAL(VkOO(i,j,k))*rho_0-zo*fai_pe(i,j,k))
c     &     - zo*REAL(ele_temp(i,j,k))*rho_0)
     &      *(rho_O_old(i,j,k)+1.0d0)
           sfe_tot=sfe_tot+(REAL(VkOH(i,j,k))*rho_0-zh*fai_pe(i,j,k))
c     &     - zh*REAL(ele_temp(i,j,k))*rho_0)
     &      *(rho_H_old(i,j,k)+1.0d0)*2
 
           detan=detan+(rho_O_old(i,j,k)-1.0d0)
           detan=detan+(rho_H_old(i,j,k)-1.0d0)*2



c           if(VintO(i,j,k).gt.50.0)then
c              VintO(i,j,k)=50.0
c           endif
          END DO
        END DO


       END DO
       sfe_tot=sfe_tot*DeltaV*rho_0*0.5
       detan=detan*DeltaV*rho_0/3.0

       sfe_tot=(sfe_tot-detan)
       write(*,*)"detan=",detan
       write(*,*)"F(HNC)=",sfe_tot*8.31*0.3,"kJ/mol",
     & sfe_tot*8.31*0.3/4.183,"kcal/mol"

c       DEALLOCATE(ele_rho)
       DEALLOCATE(VkOO)
       DEALLOCATE(VkOH)
c       DEALLOCATE(ele_temp)


c       stop "stop in Vint."

       RETURN
       END



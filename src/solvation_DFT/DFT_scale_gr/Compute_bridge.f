c******************************************************
c* THIS PROGRAM IS TO CALCULATE THE BRIDGE TERM.     ** 
c* Bridge = Excess chemical potential difference     **
c*        + Intrinsic HS potential                   **
c*                                                   **
c*                           S.L. Zhao               **              
c* The Fourier transform calculation of HS Excess    **
c* chemical potential cf.                            **
c*   M.G. Knepley et al, JCP, 132, 124101, 2010.     ** 
c* NOTE: Eq.(B14) is not correct in this article.    **
c* Correct expression can be found in Wu's MFMT work.** 
C******************************************************


      SUBROUTINE Compute_bridge
      use, intrinsic::iso_c_binding
      use rho
      use vext
      include 'param.h'
      include 'system.h'
      include 'numbers.h'
      INCLUDE 'rho.h'
      INCLUDE 'fft.h'
      INCLUDE 'bridge.h'
      INCLUDE 'c.h'
      INCLUDE 'Vext.h'
      INCLUDE 'iteration.h'

*-------------  LOCAL VARIABLES   --------------
      INTEGER i, j, k, l, n, m

      REAL*8 MIU_EX, F_rho, FEX_HS, Fint_HS,fb
      real*8 NI(7), XM, YM, ZM, EXT, BX, BY, BZ
    
      real*8 nr, ck

      integer k_index, error
      integer m1, m2, m3, ml, mm, mn
      real*8 norm_k, kx, ky, kz, k2,ff1,ff2

      complex(C_DOUBLE_COMPLEX),allocatable:: delta_rhoO(:,:,:) ! 2, nfft1, nfft2, nfft3
      complex(C_DOUBLE_COMPLEX),allocatable:: VkHS(:,:,:)
      real*8,allocatable:: Br(:,:,:)

      complex(C_DOUBLE_COMPLEX),allocatable:: input(:,:,:)
      complex(C_DOUBLE_COMPLEX),allocatable:: VkN2(:,:,:),VkN3(:,:,:)
  
      complex(C_DOUBLE_COMPLEX),allocatable::VkNV2X(:,:,:),VkNV2Y(:,:,:)
      complex(C_DOUBLE_COMPLEX),allocatable::VkNV2Z(:,:,:)

      real*8 rhoO
c interpolation
      REAL*8 HBOX_X, HBOX_Y, HBOX_Z
      integer n1, n2, n3
      character*40 file_B, file_L
      CHARACTER NAM*7, CH*3
      CHARACTER*11 chemin

       HBOX_X = LG/2.0
       HBOX_Y = LG/2.0
       HBOX_Z = LG/2.0

!       ALLOCATE(VkN2(2,mfft1,mfft2,mfft3), stat=error)
!      IF(error.ne.0) STOP "ALLOCATION VkN2 FAIL !"

!       ALLOCATE(VkN3(2,mfft1,mfft2,mfft3), stat=error)
!      IF(error.ne.0) STOP "ALLOCATION VkN3 FAIL !"
!       ALLOCATE(VkNv2x(2,mfft1,mfft2,mfft3), stat=error)
!      IF(error.ne.0) STOP "ALLOCATION VkNv2x FAIL !"
!       ALLOCATE(VkNv2y(2,mfft1,mfft2,mfft3), stat=error)
!      IF(error.ne.0) STOP "ALLOCATION VkNv2y FAIL !"
!       ALLOCATE(VkNv2z(2,mfft1,mfft2,mfft3), stat=error)
!      IF(error.ne.0) STOP "ALLOCATION VkNv2z FAIL !"

       ALLOCATE(delta_rhoO(nfft1,nfft2,nfft3), stat = error)
       IF(ERROR.NE.0) STOP "ALLOCATION delta_rhoO FAIL!"

       ALLOCATE(input(nfft1,nfft2,nfft3), stat = error)
       IF(ERROR.NE.0) STOP "ALLOCATION delta_rhoO FAIL!"
        
        DO k = 1, nfft3
         DO j = 1, nfft2
            DO i = 1, nfft1

            input(i,j,k) =cmplx(rho_O_old(i,j,k)*rho_0,0.0d0)  !  rho(r)

            delta_rhoO(i,j,k) =cmplx(rho_O_old(i,j,k)-1.0d0,0.0d0)   ! h(r)  real part
                   ! imaginary part 
           END DO
          END DO
        END DO

** test code: density profile is sucessfully transferred. **
c       OPEN(30, file = "../Results/rho_O_Z_test.dat")
c       OPEN(31, file = "../Results/rho_O_Y_test.dat")
c        DO k = 1, nfft3
c         write(30,*) (k-1)*lg/nfft3-lg/2.0, input(1,nf1,nf2,k)
c         write(31,*) (k-1)*lg/nfft3-lg/2.0, input(1,nf1,k,nf3)
c        END DO
c        close(30)
c        close(31)
** test code **


c*** Fourier transform of rho(r) to rho(k)
       
       call fftw_3d_for(nfft1,nfft2,nfft3,input,input)

** test code: Fourier transform is sucessfully implemented. ****
c       CALL fft_backward(input,fftable,fftwork,nfft1,nfft2,nfft3
c     &     ,mfft1,mfft2,mfft3,nfftable,nfftwork) ! input(k)

c       OPEN(30, file = "../Results/Input_densityO.dat")
c       DO i = 1, nfft3
c          write(30,*) (i-1)*lg/nfft3 -lg/2.0, input(1,nf1,nf2,i)
c       END DO
c       CLOSE(30)
** test code ****



       call fftw_3d_for(nfft1,nfft2,nfft3,delta_rhoO,delta_rhoO)

       ALLOCATE(VkN2(nfft1,nfft2,nfft3), stat = error)
       IF(ERROR.NE.0) STOP "ALLOCATION VkN2 FAIL!"
       ALLOCATE(VkN3(nfft1,nfft2,nfft3), stat = error)
       IF(ERROR.NE.0) STOP "ALLOCATION VkN3 FAIL!"
       ALLOCATE(VkNV2X(nfft1,nfft2,nfft3), stat = error)
       IF(ERROR.NE.0) STOP "ALLOCATION VkN2X FAIL!"
       ALLOCATE(VkNV2Y(nfft1,nfft2,nfft3), stat = error)
       IF(ERROR.NE.0) STOP "ALLOCATION VkN2Y FAIL!"
       ALLOCATE(VkNV2Z(nfft1,nfft2,nfft3), stat = error)
       IF(ERROR.NE.0) STOP "ALLOCATION VkN2Z FAIL!"

       ALLOCATE(VkHS(nfft1,nfft2,nfft3), stat = error)
       IF(ERROR.NE.0) STOP "ALLOCATION VkOO FAIL!"
       
c*** In Fourier space, calculated rho(k)* omiga(k)
c*** the expression of omiga(k) can be found in appendix B of reference paper

       DO  n=1,nfft3
          m3 = n
          if ( n .gt. nf3 )m3 = n - nfft3

          DO m=1,nfft2
             m2 = m
             if ( m .gt. nf2 )m2 = m - nfft2

            DO l=1,nfft1
               m1 = l
               if ( l .gt. nf1 ) m1 = l -  nfft1
        
             kx = (float(m1) - 1.)*PI/(HBOX_X)           
             ky = (float(m2) - 1.)*PI/(HBOX_Y)
             kz = (float(m3) - 1.)*PI/(HBOX_Z)

             k2 = sqrt(kx*kx + ky*ky + kz*kz)

c** N2 *
            If(k2.lt.1.d-6)then  ! for k2 = 0
                VkN2(l,m,n) = input(l,m,n)*PI*DP**2    ! k = 0 limit

            else                ! for k2 not equal 0
               VkN2(l,m,n)=input(l,m,n)*4.*PI*DP/2.*sin(DP/2.*k2)/k2

            endif

c** N3 *
            If(k2.lt.1.d-6)then
                VkN3(l,m,n) = input(l,m,n)*PI/6.*DP**3   ! omega = 1/6*pi*DP**3

            else
                VkN3(l,m,n) =
     &              input(l,m,n)*4.*PI/k2**3*(sin(DP/2.*k2)-  ! rho(k)*omega(3)
     &              DP/2.*k2*cos(DP/2.*k2))
            endif

c** NV2X *
            If(k2.lt.1.d-6)then
                VkNV2X(l,m,n) = cmplx(0.0,0.0)

            else
                VkNV2X(l,m,n) =cmplx(0.0,-1.0)*
     &            input(l,m,n)*4.*PI/k2**3*(sin(DP/2.*k2)-
     &            DP/2.*k2*cos(DP/2.*k2))*kx

            END IF

        
c** NV2Y *
            If(k2.lt.1.d-6)then
                VkNV2y(l,m,n) = cmplx(0.0,0.0)

            else
                VkNV2y(l,m,n) =cmplx(0.0,-1.0)*
     &              input(l,m,n)*4.*PI/k2**3*(sin(DP/2.*k2)-
     &              DP/2.*k2*cos(DP/2.*k2))*ky

            endif


c** NV2Z *
            If(k2.lt.1.d-6)then
                VkNV2z(l,m,n) = cmplx(0.0,0.0)
                
            else
                VkNV2z(l,m,n) =cmplx(0.0,-1.0)*
     &              input(l,m,n)*4.*PI/k2**3*(sin(DP/2.*k2)-
     &              DP/2.*k2*cos(DP/2.*k2))*kz

            endif


C*** INTRINSIC POTENTIAL
          norm_k = k2 
          k_index = int(norm_k/delta_ckHS)! + 1
          nr = norm_k/delta_ckHS - k_index

         IF(k_index.lt.nb_ckHS) THEN
           ck = ckHS(k_index)*(one-nr) + ckHS(k_index+1)*nr  !GOT IT
         ELSE
           ck = zero !GOT IT
         END IF

c   Compute intrinsic potential Vk(l,m,n,n_om) in k-space

           VkHS(l,m,n) = delta_rhoO(l,m,n)*ck



           END DO
          END DO
         END DO   !  end loop over k-vectors    
          

         deallocate(delta_rhoO)
         deallocate(input)

c** backward fourier transform, and get weighted density n(r)


        call fftw_3d_bak(nfft1,nfft2,nfft3,VkN2,VkN2)
        

        call fftw_3d_bak(nfft1,nfft2,nfft3,VkN3,VkN3)


        call fftw_3d_bak(nfft1,nfft2,nfft3,VkNV2X,VkNV2X)


        call fftw_3d_bak(nfft1,nfft2,nfft3,VkNV2Y,VkNV2Y)

       
        call fftw_3d_bak(nfft1,nfft2,nfft3,VkNV2Z,VkNV2Z)
**

        call fftw_3d_bak(nfft1,nfft2,nfft3,VkHS,VkHS)

        
        FEX_HS = 0.0d0
c       OPEN(40, file ="../Results/N3_FMT.dat")
        DO k=1, nfft3
         DO j=1, nfft2
          DO i=1, nfft1 
        
cc** NI(1) and NI(2) can be calculated with the knowledge of NI(3)
          NI(3)=real(VkN2(i,j,k)) ! weighted density in real space
          NI(4)=real(VkN3(i,j,k))
          NI(5)=real(VkNV2X(i,j,k))
          NI(6)=real(VkNV2Y(i,j,k))
          NI(7)=real(VkNV2Z(i,j,k))
***** CALCULATE THE EXCESS FREE ENERGY **********
          NI(1)=NI(3)/(PI*DP**2)
          NI(2)=NI(3)/(2.*PI*DP)

        IF (NI(4).LE.1.E-8.OR.(1.-NI(4)).LE.1.D-8) THEN
          F_RHO = 0.0d0
        ELSE
          F_RHO =                      ! PHI(r) in paper
     &  -NI(1)*DLOG(1.-NI(4))+(NI(2)*NI(3)-(NI(5)**2+NI(6)**2+NI(7)**2)
     &  /(2.*PI*DP))/(1.-NI(4))+1./(36.*PI)*(NI(4)*DLOG(1.-NI(4))+
     &  NI(4)**2/(1.-NI(4))**2)*(NI(3)**3-3.*NI(3)*(NI(5)**2+NI(6)**2+
     &  NI(7)**2))/NI(4)**3
        endif
        FEX_HS = FEX_HS + F_RHO*DeltaV
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(NI(4).GT.0.9)then
          NI(4)=0.9
        endif

        IF(NI(4).GT.1.D-6)THEN  ! mark
          VkN2(i,j,k)=cmplx(-LOG(1.-NI(4))/(PI*DP**2)+
     &          NI(3)/(PI*DP*(1.-NI(4)))
     &         +(NI(3)**2-(NI(5)**2+NI(6)**2+NI(7)**2))*(LOG(1.-NI(4))
     &         /NI(4)+1./(1.-NI(4))**2)/(12.*PI*NI(4)),0.0)
        ELSE
          VkN2(i,j,k)=cmplx(0.D0,0.0D0)
        ENDIF


       IF(NI(4).GT.1.D-6)THEN
          VkN3(i,j,k)=cmplx(NI(3)/(PI*(1.-NI(4))*DP**2)
     &         +(NI(3)**2-(NI(5)**2+NI(6)**2+NI(7)**2))/
     &         (2.*PI*DP*(1.-NI(4))**2)
     &        +(NI(3)**3-3.*NI(3)*(NI(5)**2+NI(6)**2+NI(7)**2))
     &        *(-LOG(1.-NI(4))/(18.*PI*NI(4)**3)
     &         -1./(36.*PI*(1.-NI(4))*NI(4)**2)
     &         -(1.-3.*NI(4))/(36.*PI*NI(4)**2*(1.-NI(4))**3)),0.0)

        ELSE
          VkN3(i,j,k)=cmplx(0.D0,0.0D0)
        ENDIF


        IF(NI(4).GT.1.D-6)THEN
          VkNV2X(i,j,k)=cmplx(-NI(5)/(PI*DP*(1.-NI(4)))
     &         -NI(3)*NI(5)*(LOG(1.-NI(4))/NI(4)+
     &         1./(1.-NI(4))**2)/(6.*PI*NI(4)),0.0)
        ELSE
          VkNV2X(i,j,k)=cmplx(0.D0,0.0D0)
        ENDIF

        IF(NI(4).GT.1.D-6)THEN
          VkNV2Y(i,j,k)=cmplx(-NI(6)/(PI*DP*(1.-NI(4)))
     &         -NI(3)*NI(6)*(LOG(1.-NI(4))/NI(4)+
     &         1./(1.-NI(4))**2)/(6.*PI*NI(4)),0.0)
        ELSE
          VkNV2Y(i,j,k)=cmplx(0.D0,0.0D0)
        ENDIF


        IF(NI(4).GT.1.D-6)THEN
          VkNV2Z(i,j,k)=cmplx(-NI(7)/(PI*DP*(1.-NI(4)))
     &         -NI(3)*NI(7)*(LOG(1.-NI(4))/NI(4)+
     &         1./(1.-NI(4))**2)/(6.*PI*NI(4)),0.0)
        ELSE
          VkNV2Z(i,j,k)=cmplx(0.D0,0.0D0)
        ENDIF


cc** Derivatives of excess free energy functional to weight density 
c** imaginary part



c*** test code ***
c        IF(i.eq.nf1.and.j.eq.nf2) then
c         write(40,*) (k-1)*lg/nfft3-lg/2.0, NI(3)
c        END IF
c*** test code ***

          END DO
         END DO
        END DO
c        close(40)


C*** VKN2, VKN3, VKNV2x, VKNV2y, VKNV2Z are derivatives of excess 
c*** chemical potential to weighted density.


        call fftw_3d_for(nfft1,nfft2,nfft3,VkN2,VkN2)

        call fftw_3d_for(nfft1,nfft2,nfft3,VkN3,VkN3)

        call fftw_3d_for(nfft1,nfft2,nfft3,VkNV2X,VkNV2X)

        call fftw_3d_for(nfft1,nfft2,nfft3,VkNV2Y,VkNV2Y)

        call fftw_3d_for(nfft1,nfft2,nfft3,VkNV2Z,VkNV2Z)

c  Loop over k-vectors
        

c*** IN FOURIER SPACE, calculate VKN2 ... * omega ..

       DO  n=1,nfft3
          m3 = n
          if ( n .gt. nf3 )m3 = n - nfft3

          DO m=1,nfft2
             m2 = m
             if ( m .gt. nf2 )m2 = m - nfft2

            DO l=1,nfft1
              m1 = l
              if ( l .gt. nf1 ) m1 = l -  nfft1

        
              kx = (float(m1) - 1.)*PI/(HBOX_X)           
              ky = (float(m2) - 1.)*PI/(HBOX_Y)
              kz = (float(m3) - 1.)*PI/(HBOX_Z)

              k2 = sqrt(kx*kx + ky*ky + kz*kz)


              If(k2.lt.1.d-6)then
                 VkN2(l,m,n) = VkN2(l,m,n)*PI*DP**2
              else
                 VkN2(l,m,n) =
     &              VkN2(l,m,n)*4.*PI*DP/2.*sin(DP/2.*k2)/k2

              endif

              If(k2.lt.1.d-6)then
                 VkN3(l,m,n) = VkN3(l,m,n)*PI/6.*DP**3
              else
                 VkN3(l,m,n) =
     &              VkN3(l,m,n)*4.*PI/k2**3*(sin(DP/2.*k2)-
     &              DP/2.*k2*cos(DP/2.*k2))
              endif

              If(k2.lt.1.d-6)then
                 VkNV2X(l,m,n) = cmplx(0.d0,0.0D0)

              else

                 VkNV2X(l,m,n) =cmplx(0.0,-1.0)*
     &              VkNV2X(l,m,n)*4.*PI/k2**3*(sin(DP/2.*k2)-
     &              DP/2.*k2*cos(DP/2.*k2))*kx

              endif
        
             if(k2.lt.1.d-6)then
                 VkNV2y(l,m,n) =cmplx(0.d0,0.0D0)

                 else

                 VkNV2y(l,m,n) =cmplx(0.0,-1.0)*
     &              VkNV2y(l,m,n)*4.*PI/k2**3*(sin(DP/2.*k2)-
     &              DP/2.*k2*cos(DP/2.*k2))*ky

                 endif


             if(k2.lt.1.d-6)then
                 VkNV2z(l,m,n) =cmplx(0.d0,0.0D0)

                 else
                 VkNV2z(l,m,n) =cmplx(0.0,-1.0)*
     &              VkNV2z(l,m,n)*4.*PI/k2**3*(sin(DP/2.*k2)-
     &              DP/2.*k2*cos(DP/2.*k2))*kz

                 endif

           END DO
          END DO
         END DO   !  end loop over k-vectors       


        call fftw_3d_bak(nfft1,nfft2,nfft3,VkN2,VkN2)

        call fftw_3d_bak(nfft1,nfft2,nfft3,VkN3,VkN3)

        call fftw_3d_bak(nfft1,nfft2,nfft3,VkNV2X,VkNV2X)

        call fftw_3d_bak(nfft1,nfft2,nfft3,VkNV2Y,VkNV2Y)

        call fftw_3d_bak(nfft1,nfft2,nfft3,VkNV2Z,VkNV2Z)

        
        
       ALLOCATE(Br(nfft1,nfft2,nfft3), stat = error)
       IF(ERROR.NE.0) STOP "ALLOCATION VkOO FAIL!"

**** for output only

       FBI = 0.0d0
       Fint_HS = 0.0d0

        Do k=1, nfft3
         Do j=1, nfft2
          Do i=1, nfft1

          Miu_ex = real(VkN2(I,J,K) + VkN3(I,J,K)-
     &        VkNV2X(I,J,K) - VkNV2Y(I,J,K) - VkNV2Z(I,J,K))


           Br(i,j,k) = Miu_ex - CHEM_P + real(VkHS(i,j,k))*rho_0    ! unit: KT

           rhoO = rho_O_old(i,j,k)-one ! hrom

           FBI = FBI -Br(i,j,k)*rho_O_old(i,j,k)

           Fint_HS = Fint_HS + rhoO*REAL(VkHS(i,j,k))
          if(IB.eq.1)then
           VintO(i,j,k) = VintO(i,j,k) + Br(i,j,k)*ratio                
          endif
c        IF(i.eq.nf1.and.j.eq.nf2) then
c        write(42,*) (k-1)*lg/nfft3 - lg/2.0 , Miu_ex
c        write(41,*) (k-1)*lg/nfft3 - lg/2.0 , Miu_ex - CHEM_P
c        END IF

          enddo
         enddo


c        write(44,*) (k-1)*lg/nfft3 - lg/2.0 , VintO(nf1,nf2,k)
c        write(43,*) (k-1)*lg/nfft3 - lg/2.0 , VkHS(1,nf1,nf2,k)*rho_0

        enddo
 
        fbi=fbi*deltaV*rho_0
        fint_hs=fint_hs*0.5*deltaV*rho_0*rho_0

        FB = Fex_HS - FBULK - detan*CHEM_P + Fint_HS+fbi
        sfe_tot=sfe_tot+fb
        write(*,*) "Fex_HS,-Fbulk,-detan*chem_P,Fint_HS,fbi"
        write(*,*) Fex_HS , -FBULK , -detan*CHEM_P , Fint_HS,fbi
        write(*,*)"FB=",fb*8.31*0.3,"kJ/mol"
        write(*,*)"Ftot=",sfe_tot*8.31*0.3,"kJ/mol",
     & sfe_tot*8.31*0.3/4.183,"kcal/mol"
c        CLOSE(41)
c        CLOSE(42)
c        CLOSE(43)
c        CLOSE(44)

        Deallocate(VKHS)
        Deallocate(VkN2)
        Deallocate(VkN3)
        Deallocate(VkNV2X)
        Deallocate(VkNV2Y)
        Deallocate(VkNV2Z)
        Deallocate(Br)





       RETURN
       END

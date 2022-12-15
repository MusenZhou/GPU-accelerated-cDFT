c******************************************************
c* THIS PROGRAM IS TO UPDATE THE SITE DENSITY        ** 
c* PROFILE OF OXYGEN AND HYROGEN IN INHOMOGENEOUS    **
c* SYSTEM.                                           **              
C******************************************************


      SUBROUTINE Compute_new_rho
      use vext
      use rho
	use control
      IMPLICIT NONE

      include 'param.h'
      include 'system.h'
      include 'numbers.h'
      INCLUDE 'Vext.h'
      INCLUDE 'rho.h'
      INCLUDE 'iteration.h'
      include 'bridge.h'

*-------------  LOCAL VARIABLES   --------------
      INTEGER i, j, k, n, m,ic
      REAL*8 IFO, IFH
      REAL*8 VO, VH
      REAL*8 X_COM, Y_COM, Z_COM, RR
C****
      CHARACTER*40 file_O, file_H
      CHARACTER NAM*7, CH*3
      CHARACTER*11 chemin

      integer iup,idown,jup,jdown,kup,kdown
      logical unval
      real*8 wiup,widown,wjup,wjdown,wkup,wkdown

ccccccccccccccccccccccccccccccccccccccccccccccc
c      integer grd_sl(2)
c      real*8 grd_cut(3)
c      grd_sl(1)=2
c      grd_sl(2)=4
c      grd_cut(1)=5.0
c      grd_cut(2)=7.0
c      grd_cut(3)=lg/2-4.0 
ccccccccccccccccccccccccccccccccccccccccccccccc



      IFO = 1.0  ! defaut value
      IFH = 1.0  ! defaut value

      write(*,*)
      write(*,*)"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
      write(*,*)"    ======   NITER =", NITER,"   ======"


c** Calcualte the intrinsic potential
      CALL  Compute_Vint


c** Calculate the bridge term 
      IF((IB.eq.1)) then
        IF(IVIN.eq.1) then
      CALL  Compute_bridge

        END IF
      END IF


	if(one_step)then
	return
	end if

      DO  k=1, nfft3
        DO  j=1, nfft2
           DO  i=1, nfft1

           x_com = (i-1)*Lg/FLOAT(nfft1)-lg/2.0
           y_com = (j-1)*Lg/FLOAT(nfft2)-lg/2.0
           z_com = (k-1)*Lg/FLOAT(nfft3)-lg/2.0

           rr = SQRT(x_com**2 + y_com**2 + z_com**2)

       IF(rr.gt.bulk_cut) then
         rho_O_new(i,j,k) = 1.0d0 !EXP(- VO)*IFO   ! g(r)

         rho_H_new(i,j,k) = 1.0d0 !EXP(- VH)*IFH

         goto 456
       endif
       
       if((VintO(i,j,k).gt.400.0).and.(VintH(i,j,k).gt.400.0))then
         rho_O_new(i,j,k) = 0.0d0 !EXP(- VO)*IFO   ! g(r)

         rho_H_new(i,j,k) = 0.0d0 !EXP(- VH)*IFH

         goto 456
       endif
 
       !do ic=1,nb_fre
        !if((rr.gt.fre_cut(ic)).and.(rr.lt.fre_cut(ic+1)))then
         !if(mod(niter,fre_sl(ic)).ne.0)then
         !rho_O_new(i,j,k) = rho_O_old(i,j,k) !EXP(- VO)*IFO   ! g(r)

         !rho_H_new(i,j,k) = rho_H_old(i,j,k) !EXP(- VH)*IFH
        
         !goto 456
         !endif
        !endif
       !enddo

       do ic=1,nb_grd
        if((rr.gt.grd_cut(ic)).and.(rr.lt.grd_cut(ic+1)))then
          if(mod(i,grd_sl(ic)).ne.0)then
             goto 456
          endif
          if(mod(j,grd_sl(ic)).ne.0)then
             goto 456
          endif
          if(mod(k,grd_sl(ic)).ne.0)then
             goto 456
          endif
        endif
       enddo

c** Calculate the intra-molecule factor 

         CALL  Compute_Intra_Factor(i, j, k, IFH, IFO)

        
!         IF(IFO.gt.20.) IFO = 20.0
!         IF(IFH.gt.20.) IFH = 20.0

c         VO = VintO(i,j,k)
c         IF(VO.lt.-5.0d0) VO = -5.0d0

c         VH = VintH(i,j,k)
c         IF(VH.lt.-5.0d0) VH = -5.0d0

         rho_O_new(i,j,k) = IFO !* EXP(- VO) ! g(r)

         rho_H_new(i,j,k) = IFH !* EXP(- VH)

            
456      continue
c         if(rr.lt.2.0)then
c           write(*,*)i,j,k,VO,VH,IFO,IFH
c         endif
*** test code ***
c         IF(IFO.gt.400.0) then
c         write(*,*) "IFO =", IFO, i,j,k
c         END IF
c         IF(IFH.gt.400.0) then
c         write(*,*) "IFH =", IFH, i,j,k
c         END IF
c**** test code ****


           END DO
         END DO

c** test code **

c       write(32, *) (k-1)*lg/nfft3-lg/2.0, rho_O_new(nf1, nf2, k)
c       write(31, *) (k-1)*lg/nfft3-lg/2.0, rho_H_new(nf1, nf2, k)

c       write(42, *) (k-1)*lg/nfft3-lg/2.0, EXP(-VO) 
c       write(41, *) (k-1)*lg/nfft3-lg/2.0, EXP(-VH)

c** test code **

       END DO

       DO  k=1, nfft3
        DO  j=1, nfft2
           DO  i=1, nfft1

           x_com = (i-1)*Lg/FLOAT(nfft1)-lg/2.0
           y_com = (j-1)*Lg/FLOAT(nfft2)-lg/2.0
           z_com = (k-1)*Lg/FLOAT(nfft3)-lg/2.0

           rr = SQRT(x_com**2 + y_com**2 + z_com**2)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           do ic=1,nb_grd
            if((rr.ge.grd_cut(ic)).and.(rr.le.grd_cut(ic+1)))then
               unval=.false.
               iup=i
               idown=i
               jup=j
               jdown=j
               kup=k
               kdown=k
               if(mod(i,grd_sl(ic)).ne.0)then
                 idown=(i/grd_sl(ic))*grd_sl(ic)
                 iup=idown+grd_sl(ic)
                 wiup=(i-idown)*1.0/(grd_sl(ic)*1.0)
                 widown=1.0-wiup
                 unval=.true.
               endif
               if(mod(j,grd_sl(ic)).ne.0)then
                 jdown=(j/grd_sl(ic))*grd_sl(ic)
                 jup=jdown+grd_sl(ic)
                 wjup=(j-jdown)*1.0/(grd_sl(ic)*1.0)
                 wjdown=1.0-wjup
                 unval=.true.
               endif
               if(mod(k,grd_sl(ic)).ne.0)then
                 kdown=(k/grd_sl(ic))*grd_sl(ic)
                 kup=kdown+grd_sl(ic)
                 wkup=(k-kdown)*1.0/(grd_sl(ic)*1.0)
                 wkdown=1.0-wkup
                 unval=.true.
               endif
               if(unval)then
                 rho_O_new(i,j,k)=
     &  rho_O_new(iup,jup,kup)*wiup*wjup*wkup
     &  +rho_O_new(idown,jup,kup)*widown*wjup*wkup
     &  +rho_O_new(iup,jdown,kup)*wiup*wjdown*wkup
     &  +rho_O_new(iup,jup,kdown)*wiup*wjup*wkdown
     &  +rho_O_new(iup,jdown,kdown)*wiup*wjdown*wkdown
     &  +rho_O_new(idown,jup,kdown)*widown*wjup*wkdown
     &  +rho_O_new(idown,jdown,kup)*widown*wjdown*wkup
     &  +rho_O_new(idown,jdown,kdown)*widown*wjdown*wkdown


                 rho_H_new(i,j,k)=
     &  rho_H_new(iup,jup,kup)*wiup*wjup*wkup
     &  +rho_H_new(idown,jup,kup)*widown*wjup*wkup
     &  +rho_H_new(iup,jdown,kup)*wiup*wjdown*wkup
     &  +rho_H_new(iup,jup,kdown)*wiup*wjup*wkdown
     &  +rho_H_new(iup,jdown,kdown)*wiup*wjdown*wkdown
     &  +rho_H_new(idown,jup,kdown)*widown*wjup*wkdown
     &  +rho_H_new(idown,jdown,kup)*widown*wjdown*wkup
     &  +rho_H_new(idown,jdown,kdown)*widown*wjdown*wkdown

               endif
             endif
            enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


           enddo
         enddo
       enddo



       write(*,*) "  # Update the density profile."
       write(*,*)

c       stop "stop in new-rho"


       RETURN
       END

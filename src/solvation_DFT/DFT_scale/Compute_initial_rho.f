c******************************************************
c* THIS PROGRAM IS TO GENERATE THE INITIAL DENSITY   ** 
c* PROFILE OF OXYGEN AND HYROGEN IN INHOMOGENEOUS    **
c* SYSTEM.                                           **              
C******************************************************


      SUBROUTINE Compute_initial_rho
      use vext
      use rho
	use control
      IMPLICIT NONE

      include 'param.h'
      include 'system.h'
      include 'numbers.h'
      INCLUDE 'Vext.h'
      INCLUDE 'rho.h'
      INCLUDE 'bridge.h'
	INCLUDE 'iteration.h'

*-------------  LOCAL VARIABLES   --------------
      INTEGER i, j, k, n, m
      REAL*8 x_com, y_com, z_com, rr
      integer,PARAMETER:: n_rg=156
      real*8 rg(n_rg),grO(n_rg),grH(n_rg),deta_rg

      integer nr,rnr
	  real*8 d1,d2

c      open(222,file='../../density_Ne/gr_oxygen.dat')
c      open(223,file='../../density_Ne/gr_hydrogen.dat')

c      do i=1,n_rg
c        read(222,*)rg(i),grO(i)
c        read(223,*)rg(i),grH(i)
c      enddo
c      deta_rg=rg(2)-rg(1)

c      close(222)
c      close(223)
	if(in_dens)then


      open(55,file="density_input.dat")
	 do i=1,nfft1
	  do j=1,nfft2
	   do k=1,nfft3

	     x_com = (i-1)*Lg/FLOAT(nfft1)-lg/2.0
           y_com = (j-1)*Lg/FLOAT(nfft2)-lg/2.0
           z_com = (k-1)*Lg/FLOAT(nfft3)-lg/2.0

           rr = SQRT(x_com**2 + y_com**2 + z_com**2) 

	read(55,*)x_com,y_com,z_com,rho_O_old(i,j,k),
     &  rho_H_old(i,j,k)
	     
	    !if(rr.gt.bulk_cut)then
	      !rho_O_old(i,j,k)=1.0
	      !rho_H_old(i,j,k)=1.0
	    !endif
		d1=1.0+in_dens_cut
		d2=1.0-in_dens_cut
	if((rho_O_old(i,j,k).gt.d1).and.(rho_O_old(i,j,k).lt.d2))then
	  rho_O_old(i,j,k)=1.0
	endif
	if((rho_H_old(i,j,k).gt.d1).and.(rho_H_old(i,j,k).lt.d2))then
	  rho_H_old(i,j,k)=1.0
	endif


	   end do
	  end do
	  end do
	  close(55)

	   else

       DO  k=1,nfft3
        DO  j=1,nfft2
          DO  i=1,nfft1

           x_com = (i-1)*Lg/FLOAT(nfft1)-lg/2.0
           y_com = (j-1)*Lg/FLOAT(nfft2)-lg/2.0
           z_com = (k-1)*Lg/FLOAT(nfft3)-lg/2.0

           rr = SQRT(x_com**2 + y_com**2 + z_com**2) 

          IF(VLJO(i,j,k).gt.0.0) then
            rho_O_old(i,j,k) = 0.0 !EXP(-beta* VLJO(i,j,k))   ! g(r)
            rho_H_old(i,j,k) = 0.0 !EXP(-beta* VLJO(i,j,k))   ! same as rho_O
          ELSE
            rho_O_old(i,j,k) = 1.0 !EXP(-beta* VLJO(i,j,k))   ! g(r)
            rho_H_old(i,j,k) = 1.0 !EXP(-beta* VLJO(i,j,k))   ! same as rho_O
          END IF

 
c           IF(rr.lt.(LG/3.0-1.0)) THEN
c           rho_O_old(i,j,k) = EXP(-beta* VLJO(i,j,k))   ! g(r)

c           rho_H_old(i,j,k) = EXP(-beta* VLJO(i,j,k))   ! same as rho_O
c           ELSE
c           rho_O_old(i,j,k) = 1.0 !EXP(-beta* VLJO(i,j,k))   ! g(r)
c           rho_H_old(i,j,k) = 1.0 !EXP(-beta* VLJO(i,j,k))   ! same as rho_O
c           END IF

c            read(55,*)x_com,y_com,z_com,rho_H_old(i,j,k),
c     &  rho_O_old(i,j,k)

c            nr=int((rr-rg(1))/deta_rg)+1
c            rnr=(rr-rg(1))/deta_rg-nr+1
c            if(nr.lt.1)then
c              rho_O_old(i,j,k) = 0.0D0
c              rho_H_old(i,j,k) = 0.0D0
c            elseif(rr.gt.10.0)then
c              rho_O_old(i,j,k) = 1.0D0
c              rho_H_old(i,j,k) = 1.0D0
c            else
c              rho_O_old(i,j,k) = (1.0D0-rnr)*grO(nr)+rnr*grO(nr+1)
c              rho_H_old(i,j,k) = (1.0D0-rnr)*grH(nr)+rnr*grH(nr+1)
c            endif
          

c           read(41,*) rho_O_old(i,j,k)
c           read(42,*) rho_H_old(i,j,k)

c           rho_O_old(i,j,k) = EXP(-beta* VextO(i,j,k))   ! g(r)
c           rho_H_old(i,j,k) = EXP(-beta* VextH(i,j,k))   ! same as rho_O

           END DO
         END DO


       END DO

	end if
c       close(55)

       write(*,*) "  # Generate initial density profile."


       RETURN
       END

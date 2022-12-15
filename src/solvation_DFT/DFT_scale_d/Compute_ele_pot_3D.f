c*********************************************************
c* THIS SUBROUTINE IS TO SOLVE THE 3D POISSON EQUATION. **
C* THE SOLVER IS FROM INTEL INC.                        **
c* /opt/intel/composerxe/mkl/examples/pdepoissonf/source **
C*********************************************************

      Subroutine Compute_ele_pot_3D
      use mkl_poisson
      use mkl_service 
      use rho
      use pe
      use vext
      implicit none

      include 'param.h'
      include 'fft.h'
      include 'system.h'
      include 'numbers.h'
      include 'c.h'
      include 'rho.h'
      include 'Vext.h'
      include 'iteration.h'
      include 'pe.h'
      

      integer ix, iy, iz, i, stat
      integer ipar(128)
      double precision ax, bx, ay, by, az, bz, lx, ly, lz, xi, yi, zi
      double precision cx, cy, cz, c1
      
      double precision q
      type(DFTI_DESCRIPTOR), pointer :: xhandle, yhandle
      character(6) BCtype
      real*8 x_com, y_com, z_com, rr,etot,etemp

      write(*,*)'  ** solve 3D Poisson equation. **'

      ax=-lg/2.
      bx=lg/2.
      ay=-lg/2.
      by=lg/2.
      az=-lg/2.
      bz=lg/2.
      q=0.
cccccccccccccccccccinitialize f ccccccccccccccccccc
      
      do iz=1,nfft3
        do iy=1,nfft2
          do ix=1,nfft1
c*** Note Poisson equation here:  U_xx + U_yy + U_zz = - fai_pe 
            
           x_com = (ix-1)*Lg/FLOAT(nfft1)-lg/2.0
           y_com = (iy-1)*Lg/FLOAT(nfft2)-lg/2.0
           z_com = (iz-1)*Lg/FLOAT(nfft3)-lg/2.0

           rr = SQRT(x_com**2 + y_com**2 + z_com**2)

          IF(rr.lt.bulk_cut) then

c          fai_pe(ix,iy,iz)= -1209.1*rho_0*(2.*zh*rho_H_old(ix,iy,iz)
c     &      +zo*rho_O_old(ix,iy,iz))
c           etemp=ele_rho(1,ix-1,iy,iz)+ele_rho(1,ix+1,iy,iz)
c           etemp=etemp+ele_rho(1,ix,iy-1,iz)+ele_rho(1,ix,iy+1,iz)
c           etemp=etemp+ele_rho(1,ix,iy,iz-1)+ele_rho(1,ix,iy,iz+1)
c           etemp=etemp/6.0
           etemp=ele_rho(ix,iy,iz)
c           fai_pe(ix,iy,iz)=-1209.1*2*rho_0*etemp !spce
c           fai_pe(ix,iy,iz)=-1070.16*2*rho_0*etemp !spc
            fai_pe(ix,iy,iz)=-1246.491*2*rho_0*etemp !tip3p
c           fai_pe(ix,iy,iz)=-1209.1*2*rho_0*ele_rho(1,ix,iy,iz)
          ELSE

          fai_pe(ix,iy,iz)= 0.0d0 
 
          END IF

          enddo
        enddo
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccc
      etot=-1209.1*ele_tot*deltaV*rho_0
      
cccccccccccccccccccbound conditionccccccccccccccccc
      BCtype='DDDDDD'   ! Dirichlet boundary condition at 6 planes

      do iz=1,nfft3
        do iy=1,nfft2
         x_com = ax
         y_com = (iy-1)*Lg/FLOAT(nfft2)-lg/2.0
         z_com = (iz-1)*Lg/FLOAT(nfft3)-lg/2.0

         rr = SQRT(x_com**2 + y_com**2 + z_com**2)

         bd_ax(iy,iz)=0.0!etot/4.0D0/pi/rr
         bd_bx(iy,iz)=0.0!etot/4.0D0/pi/rr

       enddo
      enddo

      do iz=1,nfft3
        do ix=1,nfft1
           x_com = (ix-1)*Lg/FLOAT(nfft1)-lg/2.0
           y_com = ay
           z_com = (iz-1)*Lg/FLOAT(nfft3)-lg/2.0

           rr = SQRT(x_com**2 + y_com**2 + z_com**2)


         bd_ay(ix,iz)=0.0!etot/4.0D0/pi/rr
         bd_by(ix,iz)=0.0!etot/4.0D0/pi/rr
       enddo
      enddo

      do iy=1,nfft2
        do ix=1,nfft1
           x_com = (ix-1)*Lg/FLOAT(nfft1)-lg/2.0
           y_com = (iy-1)*Lg/FLOAT(nfft2)-lg/2.0
           z_com = az

           rr = SQRT(x_com**2 + y_com**2 + z_com**2)


          bd_az(ix,iy)=0.0!etot/4.0D0/pi/rr
          bd_bz(ix,iy)=0.0!etot/4.0D0/pi/rr
        enddo
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,128
        ipar(i)=0
      enddo
ccccccccccccccccccccsolve Possion equationcccccccccc

      call d_init_Helmholtz_3D(ax,bx,ay,by,az,bz,nfft1-1,nfft2-1,
     &                       nfft3-1,BCtype,q,ipar,dpar,stat)

      IF(stat.ne.0) then
       write(*,*) "3D Poisson equation fails to reach solution." 
       stop
      END IF

      
      call d_commit_Helmholtz_3D(fai_pe,bd_ax,bd_bx,bd_ay,bd_by,bd_az,
     &                   bd_bz,xhandle, yhandle, ipar, dpar, stat)

      IF(stat.ne.0) then
       write(*,*) "3D Poisson equation fails to reach solution." 
       stop
      END IF
cc output is also written in: fai_pe


      call d_Helmholtz_3D(fai_pe, bd_ax, bd_bx, bd_ay, bd_by, bd_az,  
     &                  bd_bz,xhandle, yhandle, ipar, dpar, stat)

      IF(stat.ne.0) then
       write(*,*) "3D Poisson equation fails to reach solution." 
       stop
      END IF

cccccccccccccccccccccccccccccccccccccccccccccccccccc

      call free_Helmholtz_3D(xhandle, yhandle, ipar, stat)
      IF(stat.ne.0) then
       write(*,*) "3D Poisson equation fails to reach solution." 
       stop
      END IF
cccccccccccccccccccccccccccccccccccccccccccccccccccc



      call mkl_free_buffers


      return
      end





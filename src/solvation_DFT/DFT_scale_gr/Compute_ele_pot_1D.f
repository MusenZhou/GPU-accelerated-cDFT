c*********************************************************
c* THIS SUBROUTINE IS TO SOLVE THE 3D POISSON EQUATION. **
C* THE SOLVER IS FROM INTEL INC.                        **
c* /opt/intel/composerxe/mkl/examples/pdepoissonf/source **
C*********************************************************

      Subroutine Compute_ele_pot_1D
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
      

      integer ix, iy, iz, nhist, stat
      real*8 ax, bx, ay, by, az, bz, lx, ly, lz, xi, yi, zi
      real*8 cx, cy, cz, c1
      
      real*8 x_com, y_com, z_com, rr, radius
      real*8 gO(mfft3), gH(mfft3)
      real*8 nr
      real*8 got,ght,deta_fai
      write(*,*)'  ** solve 1D Poisson equation. **'

cccccccccccccccccccinitialize f ccccccccccccccccccc
      Open(81, file = "../Results/gO_1D.dat")
      Open(71, file = "../Results/gH_1D.dat")
      
      got=0.0D0
      ght=0.0D0
      do iz=1, nfft3
            
         gO(iz)= rho_O_old(nf1,nf2,iz)*rho_0
         gH(iz)= rho_H_old(nf1,nf2,iz)*rho_0 
 
         if(go(iz).lt.1.0D-10)then
          go(iz)=0.0D0
         endif
  
         if(gh(iz).lt.1.0D-10)then
          gh(iz)=0.0D0
         endif
         
         rr = (iz-1) * LG/Float(nfft3) - LG/2

         got=got+go(iz)*rr*rr
         ght=ght+gh(iz)*rr*rr

      enddo
     
      close(81)
      close(71)
ccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccsolve Possion equationcccccccccc

      Open(82, file = "../Results/ele_pot_1D.dat")
       DO nhist = 1, nf1

         radius = (nhist -1) * LG/Float(nfft1) 

         fai_pe1D(nhist) = 0.0

        IF(radius.gt.1.0d-6) then

         DO iz =nf1+nhist, nfft3

          rr = (iz-1) * LG/Float(nfft3) - LG/2

         fai_pe1D(nhist) = fai_pe1D(nhist)  
     & + gO(iz) * zo *(rr - rr**2/radius)*lg/float(nfft3)
     & + 2.0*gH(iz) * zh *(rr - rr**2/radius)*lg/float(nfft3)


         END DO
    
       END IF
        
        fai_pe1D(nhist) = fai_pe1D(nhist) * 1209.1

      write(82,*) radius,fai_pe1D(nhist)
      END DO

      close(82)

      write(*,*) "  # Obtain Electro-static potential (1D poisson Eqn)."

      Open(83, file = "../Results/ele_pot_Z.dat")

      do iz=1,nfft3
        do iy=1,nfft2
          do ix=1,nfft1

           x_com = (ix-1)*Lg/FLOAT(nfft1)-lg/2.0
           y_com = (iy-1)*Lg/FLOAT(nfft2)-lg/2.0
           z_com = (iz-1)*Lg/FLOAT(nfft3)-lg/2.0

           rr = SQRT(x_com**2 + y_com**2 + z_com**2)

           nhist = INT(rr*nfft3/lg) + 1
           nr = rr*nfft3/lg +1 - nhist

          IF(nhist.lt.nfft3) then

c          fai_pe(ix,iy,iz)=  fai_pe1D(nhist)*(1-nr) 
c     &                     + fai_pe1D(nhist+1)*nr

          ELSE

c          fai_pe(ix,iy,iz)= 0.0d0

          END IF

          enddo
        enddo
c       write(83,*) z_com, fai_pe(nf1,nf2, iz)

      enddo
      close(83)
     
c      stop "stop in poisson equaion 1D"

      return
      end





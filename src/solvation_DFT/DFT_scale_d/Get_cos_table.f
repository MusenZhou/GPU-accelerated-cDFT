c*********************************************************
c* THIS SUBROUTINE IS TO PREPARE THE TABLE OF COS_THETA **
C* AND COS_PHI FOR THE FUTURE USE.                      **
C*                               S.L. ZHAO              **
C*********************************************************

        
       SUBROUTINE Get_cos_table
       IMPLICIT NONE
      
       INCLUDE "numbers.h"
       INCLUDE "Cos_theta.h"
 
cc **** local variable
       INTEGER I, J, K
       real*8 phi,fai(nb_phi)
       real*8 xn, yn, zn, rx, ry, rz
       INTEGER icos, iphi, jphi

       
       DELTA_COS = 2.0/nb_cos
       DELTA_PHI = TWOPI/nb_PHI

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       DO i = 1, nb_cos
c        cos_theta(i) = (i-0.5)*delta_cos - 1.0
c        sin_theta(i) = sqrt(1.0 - cos_theta(i)**2)
c        
c       DO J = 1, nb_phi
c          phi = (J - 0.5) * delta_phi

c          cos_phi(J) = cos(phi)
c          sin_phi(J) = sin(phi)
        
c        DX_thetaPhi(J, I) = Sin_theta(I)*cos_phi(J)
c        DY_thetaPhi(J, I) = Sin_theta(I)*Sin_phi(J)
c        DZ_thetaPhi(J, I) = Cos_theta(I)

c        END DO
c        END DO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccc Guass Legrand integration cccccccccccccccccccc
cccccccccc 5 point cccccccccccccccc
        a12(1)=0.2369269
        a12(2)=0.4786287
        a12(3)=0.5688889
        a12(4)=0.4786287
        a12(5)=0.2369269
ccccccccccccccccccccccccccccccccccc
cccccccccc 6 point cccccccccccccccc
c        a12(1)=0.17132449
c        a12(2)=0.36076157
c        a12(3)=0.46791363
c        a12(6)=0.17132449
c        a12(5)=0.36076157
c        a12(4)=0.46791363
cccccccccccccccccccccccccccccccccccc
cccccccccc 7 point ccccccccccccccccc
c         a12(1)=0.1294849662D0
c         a12(2)=0.2797053915D0
c         a12(3)=0.3818300505D0
c         a12(4)=0.4179591837D0
c         a12(5)=a12(3)
c         a12(6)=a12(2)
c         a12(7)=a12(1)
cccccccccccccccccccccccccccccccccccc
cccccccccc  10 point ccccccccccccccc
c        a12(1)=0.2955242247D0
c        a12(2)=0.2692667193D0
c        a12(3)=0.2190863625D0
c        a12(4)=0.1494513492D0
c        a12(5)=0.0666713443D0
c        a12(6)=a12(5)
c        a12(7)=a12(4)
c        a12(8)=a12(3)
c        a12(9)=a12(2)
c        a12(10)=a12(1)
cccccccccccccccccccccccccccccccccccc
cccccccccc Guass Lobatto 9 ccccccccc
c        a12(1)=0.0277777777D0
c        a12(2)=0.1654953615D0
c        a12(3)=0.2745387125D0
c        a12(4)=0.3464285109D0
c        a12(5)=0.3715192743D0
c        a12(6)=a12(4)
c        a12(7)=a12(3)
c        a12(8)=a12(2)
c        a12(9)=a12(1)
cccccccccccccccccccccccccccccccccccc
        do i=1,nb_phi
          a2pi(i)=a12(i)*pi
        enddo

cccccccccccc  5 ccccccccccccc
        cos_theta(1) = -0.9061798
        cos_theta(2) = -0.5384693
        cos_theta(3) = 0.0
        cos_theta(4) = 0.5384693
        cos_theta(5) = 0.9061798
ccccccccccccccccccccccccccccc
cccccccccccccc  6  cccccccccccc
c         cos_theta(1)=-0.93246951
c         cos_theta(2)=-0.66120939
c         cos_theta(3)=-0.23861919
c         cos_theta(6)=+0.93246951
c         cos_theta(5)=+0.66120939
c         cos_theta(4)=+0.23861919
cccccccccccccccccccccccccccccccc
cccccccccccccc 7 ccccccccccccccc
c         cos_theta(1)=-0.9491079123D0
c         cos_theta(2)=-0.7415311856D0
c         cos_theta(3)=-0.4058451514D0
c         cos_theta(4)=0.0D0
c         cos_theta(5)=-cos_theta(3)
c         cos_theta(6)=-cos_theta(2)
c         cos_theta(7)=-cos_theta(1)
cccccccccccccccccccccccccccccccc
cccccccccccccc 10 cccccccccccccc
c         cos_theta(1)=-0.1488743390D0
c         cos_theta(2)=-0.4333953941D0
c         cos_theta(3)=-0.6794095683D0
c         cos_theta(4)=-0.8650633667D0
c         cos_theta(5)=-0.9739065285D0
c         cos_theta(6)=-cos_theta(5)
c         cos_theta(7)=-cos_theta(4)
c         cos_theta(8)=-cos_theta(3)
c         cos_theta(9)=-cos_theta(2)
c         cos_theta(10)=-cos_theta(1)
cccccccccccccccccccccccccccccccc
cccccccccccc Lobatto 9 cccccccccccc
c         cos_theta(1)=-1.0D0
c         cos_theta(2)=-0.8997579954D0
c         cos_theta(3)=-0.6771862795D0
c         cos_theta(4)=-0.3631174638D0
c         cos_theta(5)=0.0D0
c         cos_theta(6)=-cos_theta(4)
c         cos_theta(7)=-cos_theta(3)
c         cos_theta(8)=-cos_theta(2)
c         cos_theta(9)=-cos_theta(1)
ccccccccccccccccccccccccccccccccccc
c        fai(1)=0.2947449
c        fai(2)=1.449941
c        fai(3)=pi
c        fai(4)=4.8332439
c        fai(5)=5.9884405


        DO i = 1, nb_cos
        sin_theta(i) = sqrt(1.0 - cos_theta(i)**2)
        
       DO J = 1, nb_phi
          phi=(1.0+cos_theta(j))*pi
          cos_phi(J) = cos(phi)
          sin_phi(J) = sin(phi)
        
        DX_thetaPhi(J, I) = Sin_theta(I)*cos_phi(J)*0.9572D0
        DY_thetaPhi(J, I) = Sin_theta(I)*Sin_phi(J)*0.9572D0
        DZ_thetaPhi(J, I) = Cos_theta(I)*0.9572D0

        END DO
        END DO

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c*******  additional table ***************
c      Cos_theta0 = Cos(70.53/180.0*PI)  ! here 70.53 = 180 - 109.47 in SPC/E water model
      Cos_theta0 = Cos(75.48/180.0*PI)  !TIP3P
      Sin_theta0 = DSQRT(1.0 - Cos_theta0**2)

cc test code **
c      write(*,*) "Cos_theta0 =", Cos_theta0
c      write(*,*) "Sin_theta0 =", Sin_theta0
c      stop
cc test code **

c      Cos_theta1 = Cos(109.47/180*PI)  ! here 70.53 = 180 - 109.47 in SPC/E water model
      Cos_theta1 = Cos(104.52/180*PI)
      Sin_theta1 = DSQRT(1.0 - Cos_theta1**2)
ccc*****************************************

        DO icos = 1, nb_cos
         DO iphi = 1, nb_phi
          DO jphi = 1, nb_phi

cc*********
            xn = sin_theta0*cos_phi(jphi)
            yn = sin_theta0*sin_phi(jphi)
            zn = cos_theta0

           call transf_coord(xn, yn, zn, icos, iphi, rx, ry, rz)

            RX_thetaphi0(jphi, iphi, icos) = Rx*0.9572D0
            RY_thetaphi0(jphi, iphi, icos) = Ry*0.9572D0
            RZ_thetaphi0(jphi, iphi, icos) = Rz*0.9572D0


cc*********
            xn = sin_theta1*cos_phi(jphi)
            yn = sin_theta1*sin_phi(jphi)
            zn = cos_theta1

          call transf_coord(xn, yn, zn, icos, iphi, rx, ry, rz)

            RX_thetaphi1(jphi, iphi, icos) = Rx*0.9572D0
            RY_thetaphi1(jphi, iphi, icos) = Ry*0.9572D0
            RZ_thetaphi1(jphi, iphi, icos) = Rz*0.9572D0

         END DO
        END DO
       END DO



       RETURN
       END 



cc ####################################################################
cc this subroutine give the coordinates (x,y,z)  in the old frame
cc givn we know the its coordinates (xn,yn,zn) in the new frame.
cc where new frame is obtained by rotating (theta,phi) of the old frame
cc Viz.                r_old = r_new * M(theta,phi) 
cc HERE M(theta, phi) is the reverse transformation matrix
cc                                 by S.L. Zhao
cc program tested: set (xn, yn, zn) = (0, 0, a), and set phi = 0
cc then  x = a*sin(theta), y = 0, and z =a*cos(theta) 
cc ####################################################################

       subroutine transf_coord(xn, yn, zn, ncos, nphi, x, y, z)
       implicit none
       include 'numbers.h'
       include 'Cos_theta.h'

       real*8 x,y,z
       real*8 xn,yn,zn
       integer ncos, nphi
       real*8 costheta, sintheta, sinphi, cosphi

       costheta = cos_theta(ncos)
       sintheta = sin_theta(ncos)

       sinphi = sin_phi(nphi)
       cosphi = cos_phi(nphi)


       x = costheta*cosphi*xn - sinphi*yn + sintheta*cosphi*zn

       y = costheta*sinphi*xn + cosphi*yn + sintheta*sinphi*zn

       z = -sintheta*xn + costheta*zn


c       x = xn*cos_theta*cos_phi+yn*cos_theta*sin_phi-zn*sin_theta

c       y = costheta*sinphi*xn + cosphi*yn + sintheta*sinphi*zn

c       z = -sintheta*xn + costheta*zn

       return
       end


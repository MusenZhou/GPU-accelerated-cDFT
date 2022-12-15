c******************************************************
c* THIS PROGRAM IS TO CALCULATE THE INTRA-CORRELATED ** 
c* FACTOR for hydrogen.                              **
c*                           S.L. Zhao               **              
c*   Many tables applied in this subroutine are      **
c*   prepared in priori "Get_cos_table.f".           **
C******************************************************


      SUBROUTINE Compute_Intra_Factor(l, m, n, IFH, IFO)
      use vext
      use rho
      IMPLICIT NONE

      include 'param.h'
      include 'system.h'
      include 'numbers.h'
      INCLUDE 'Vext.h'
      INCLUDE 'rho.h'
      INCLUDE 'Cos_theta.h'

*-------------  LOCAL VARIABLES   --------------
      INTEGER l, m, n
      real*8 IFO, IFH,if1,delV
      INTEGER i, j, k
      INTEGER icos, iphi, jphi
      REAL*8 xx, yy, zz, rx, ry, rz, xn, yn, zn
      REAL*8 x_O, y_O, z_O, x_H, y_H, z_H
      REAL*8 BFO, BFH    ! for factor on hydrogen
      REAL*8 BFH1, BFH2  ! for factor on oxygen
      INTEGER INSIDE
      REAL*8 VO, VH,voc,vhc,vtp
      real*8 DeltaX(0:1), DeltaY(0:1), DeltaZ(0:1)
      real*8 GG,ifo1,ifo2,ifh1,ifh2
      integer n1, n2, n3

      voc=VintO(l,m,n)
      vhc=VintH(l,m,n)

       IFH = 0.0
       IFO = 0.0
       if1=0.0
       DO ICOS = 1, nb_cos   ! first loop
          ifo1=0.0
          ifh1=0.0
          DO iphi = 1, nb_phi  ! second loop
           ifo2=0.0
           ifh2=0.0
** intra factor for hydrogen site
c%% first address r' for oxygen, r' = (xx, yy, zz)
c%%     r = (l, m, n) on grid. (l-1)*lg/nfft1, (m-1)*lg/nfft2, (n-1)*lg/nfft3


cc** DX_thetaPhi(iphi, icos) is prepared in priori in "Get_cos_theta.f"

            xx = (l-1)*lg/FLOAT(nfft1) + DX_thetaPhi(iphi,icos) !sin_theta(icos)*cos_phi(iphi)
            yy = (m-1)*lg/FLOAT(nfft2) + DY_thetaPhi(iphi,icos) !Sin_theta(icos)*sin_phi(iphi)
            zz = (n-1)*lg/FLOAT(nfft3) + DZ_thetaPhi(iphi,icos) !cos_theta(icos)

            INSIDE = 1

            IF(xx.gt.lg.or.xx.lt.zero) THEN 
               BFO = 0.0d0     ! Boltzmann factor for oxygen
               INSIDE = 0
               BFH1 = 0.0d0           
            END IF
            IF(yy.gt.lg.or.yy.lt.zero) THEN
               BFO = 0.0d0
               INSIDE = 0
               BFH1 = 0.0d0
            END IF
            IF(zz.gt.lg.or.zz.lt.zero) THEN
               BFO = 0.0d0
               INSIDE = 0
               BFH1 = 0.0d0
            END IF


            IF(INSIDE.EQ.1) THEN
              I = INT(xx*nfft1/lg) + 1
              J = INT(yy*nfft2/lg) + 1
              K = INT(zz*nfft3/lg) + 1

** interpolation 1 *****
             if(I.ge.nfft1) I = nfft1 - 1
             if(J.ge.nfft2) J = nfft2 - 1
             if(K.ge.nfft3) K = nfft3 - 1

            DeltaX(1) = (xx*nfft1/lg + one) - I
            DeltaX(0) = one - DeltaX(1)

            DeltaY(1) = (yy*nfft2/lg + one) - J
            DeltaY(0) = one - DeltaY(1)

            DeltaZ(1) = (zz*nfft3/lg + one) - K
            DeltaZ(0) = one - DeltaZ(1)

             VO = 0.0d0
             VH = 0.0d0

             DO n1 = 0, 1
              DO n2 = 0, 1
               DO n3 = 0, 1

            VO = VO + 
     &      VintO(i+n1,j+n2,k+n3)*deltaX(n1)*deltaY(n2)*deltaZ(n3) 

            VH = VH + 
     &      VintH(i+n1,j+n2,k+n3)*deltaX(n1)*deltaY(n2)*deltaZ(n3) 

               END DO
              END DO
            END DO

c            if(VO.lt.100.0)then
c            BFO = EXP(- VO)
c            else
c            BFO=0.0D0
c            endif
c            if(VH.lt.100.0)then
c            BFH1 = EXP(- VH)
c            else
c            BFH1=0.0D0
c            endif
             BFO=VO
             BFH1=VH
** interpolation 1 *******

c           BFO = EXP(- VintO(i,j,k))
c           BFH1 = EXP(-VintH(i,j,k))

            END IF


cc test code *****
c             IF(BFH1.gt.8.0) then
c               write(*,*) "BFH1 =", BFH1
c               write(*,*) "VintH =", VintH(i,j,k)
c             END IF
cc test code *****



           DO jphi = 1, nb_phi    ! third and last loop
c            delV=a2pi(jphi)*a2pi(iphi)*a12(icos)
cccc for factor on hydrogen site
c%% then address r'' for the other hydrogen ( r is for first hydrogen)
c%%       r'' = r' + (sin_theta0*cos_phi, sin_theta0 * sin_phi, cos_theta0)*M(theta,iphi) 
c%%       r'' = (x_H, y_H, z_H)

c*******************************************
c            xn = sin_theta0*cos_phi(jphi)
c            yn = sin_theta0*sin_phi(jphi)
c            zn = cos_theta0

c          call transf_coord(xn, yn, zn, icos, iphi, rx, ry, rz) ! this part has been stored in a table
c********************************************


           x_H = xx + RX_thetaphi0(jphi, iphi, icos)  ! read from stored table 
           y_H = yy + RY_thetaphi0(jphi, iphi, icos)
           z_H = zz + RZ_thetaphi0(jphi, iphi, icos)


            INSIDE = 1

            IF(x_H.gt.lg.or.x_H.lt.zero) THEN
               BFH = 0.0d0
               INSIDE = 0
            END IF
            IF(y_H.gt.lg.or.y_H.lt.zero) THEN
               BFH = 0.0d0
               INSIDE = 0
            END IF
            IF(z_H.gt.lg.or.z_H.lt.zero) THEN
               BFH = 0.0d0
               INSIDE = 0
            END IF

            IF(INSIDE.EQ.1) THEN
              I = INT(x_H*nfft1/lg) + 1
              J = INT(y_H*nfft2/lg) + 1
              K = INT(z_H*nfft3/lg) + 1

** interpolation 2 *******
             if(I.ge.nfft1) I = nfft1 - 1
             if(J.ge.nfft2) J = nfft2 - 1
             if(K.ge.nfft3) K = nfft3 - 1

            DeltaX(1) = (x_H*nfft1/lg + one) - I
            DeltaX(0) = one - DeltaX(1)

            DeltaY(1) = (y_H*nfft2/lg + one) - J
            DeltaY(0) = one - DeltaY(1)

            DeltaZ(1) = (z_H*nfft3/lg + one) - K
            DeltaZ(0) = one - DeltaZ(1)


             VH = 0.0d0

             DO n1 = 0, 1
              DO n2 = 0, 1
               DO n3 = 0, 1

            VH = VH +
     &      VintH(i+n1,j+n2,k+n3)*deltaX(n1)*deltaY(n2)*deltaZ(n3)

               END DO
              END DO
            END DO

c            if(VH.lt.100.0)then
c            BFH = EXP(- VH)
c            else
c            BFH=0.0D0
c            endif
             BFH=VH
** interpolation 2 **********

c            BFH = EXP(- VintH(i,j,k))

            END IF

           vtp=BFO+BFH+vhc
           if(vtp.lt.10.0)then
            IFH2=IFH2+exp(-vtp)*a2pi(jphi)
           endif

ccc***** for factor on oxygen site ******************
c%% then address r'' for the other hydrogen (r is for oxygen, and r' is for first hydrogen)
c%%       r'' = r + M(theta,iphi)*(sin_theta1*cos_phi, sin_theta1 * sin_phi, cos_theta1)
c%%       r'' = (x_H, y_H, z_H)
ccc**********************************************
c            xn = sin_theta1*cos_phi(jphi)
c            yn = sin_theta1*sin_phi(jphi)
c            zn = cos_theta1

c          call transf_coord(xn, yn, zn, icos, iphi, rx, ry, rz)  ! this part has been stored in a table
ccc**********************************************

            x_H = (l-1)*lg/FLOAT(nfft1) + RX_thetaphi1(jphi,iphi,icos) ! read from stored table
            y_H = (m-1)*lg/FLOAT(nfft2) + RY_thetaphi1(jphi,iphi,icos)
            z_H = (n-1)*lg/FLOAT(nfft3) + RZ_thetaphi1(jphi,iphi,icos) 

            INSIDE = 1

            IF(x_H.gt.lg.or.x_H.lt.zero) THEN
               BFH2 = 0.0d0
               INSIDE = 0
            END IF
            IF(y_H.gt.lg.or.y_H.lt.zero) THEN
               BFH2 = 0.0d0
               INSIDE = 0
            END IF
            IF(z_H.gt.lg.or.z_H.lt.zero) THEN
               BFH2 = 0.0d0
               INSIDE = 0
            END IF

            IF(INSIDE.EQ.1) THEN
              I = INT(x_H*nfft1/lg) + 1
              J = INT(y_H*nfft2/lg) + 1
              K = INT(z_H*nfft3/lg) + 1

** interpolation 3  ***********
             if(I.ge.nfft1) I = nfft1 - 1
             if(J.ge.nfft2) J = nfft2 - 1
             if(K.ge.nfft3) K = nfft3 - 1

            DeltaX(1) = (x_H*nfft1/lg + one) - I
            DeltaX(0) = one - DeltaX(1)

            DeltaY(1) = (y_H*nfft2/lg + one) - J
            DeltaY(0) = one - DeltaY(1)

            DeltaZ(1) = (z_H*nfft3/lg + one) - K
            DeltaZ(0) = one - DeltaZ(1)


             VH = 0.0d0

             DO n1 = 0, 1
              DO n2 = 0, 1
               DO n3 = 0, 1

            VH = VH +
     &      VintH(i+n1,j+n2,k+n3)*deltaX(n1)*deltaY(n2)*deltaZ(n3)

               END DO
              END DO
            END DO

c            if(VH.lt.100)then
c              BFH2 = EXP(- VH)
c            else
c              BFH2=0.0D0
c            endif
            BFH2=VH
** interpolation 3  *********


c           BFH2 = EXP(- VintH(i,j,k))

            END IF
        
         vtp=BFH1+BFH2+voc
         if(vtp.lt.10.0)then
          IFO2=IFO2+exp(-vtp)*a2pi(jphi)
         endif
         
c         write(*,*) "IFO =", IFO
c         write(*,*) "BFH1, BFH2 =", BFH1, BFH2


             END DO      ! end loop on jphi
         ifo1=ifo1+ifo2*a2pi(iphi)
         ifh1=ifh1+ifh2*a2pi(iphi)
           END DO      ! end loop on iphi
         ifo=ifo+ifo1*a12(icos)
         ifh=ifh+ifh1*a12(icos)
        END DO    ! end loop on icos


        IFH = IFH/(8.0*PI**2)
        IFO = IFO/(8.0*PI**2)


*** test code ***
c        IF(IFO.gt.30) then

c        write(*,*) "IFH =", IFH
c        write(*,*) "IFO =", IFO
c        write(*,*) " l, m, n =", l, m, n
c        stop
c        END IF
*** test code ***


       RETURN
       END





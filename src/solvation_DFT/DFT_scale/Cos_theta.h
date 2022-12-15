* Cos_theta

      INTEGER,parameter:: nb_cos=5, nb_phi=5

      REAL*8 delta_cos, delta_phi
      COMMON/DELTA_COS/ delta_cos, delta_phi
     
      REAL*8 cos_theta(nb_cos), cos_phi(nb_phi),
     &       sin_theta(nb_cos), sin_phi(nb_phi)
      common/cos_value/cos_theta, sin_theta, cos_phi, sin_phi

      REAL*8 Cos_theta0, Sin_theta0
      Common/theta0/Cos_theta0, Sin_theta0

      REAL*8 Cos_theta1, Sin_theta1
      Common/theta1/Cos_theta1, Sin_theta1

      REAL*8 DX_thetaPhi(nb_phi, nb_cos)  ! used in calculation of intra-factor
      REAL*8 DY_thetaPhi(nb_phi, nb_cos)
      REAL*8 DZ_thetaPhi(nb_phi, nb_cos)
      Common/D_thetaphi/DX_thetaPhi, DY_thetaPhi, DZ_thetaPhi

      REAL*8 RX_thetaphi0(nb_phi, nb_phi, nb_cos)
      REAL*8 RY_thetaphi0(nb_phi, nb_phi, nb_cos)
      REAL*8 RZ_thetaphi0(nb_phi, nb_phi, nb_cos)
      Common/R_thetaPhi/RX_thetaPhi0, RY_thetaPhi0, RZ_thetaPhi0

      REAL*8 RX_thetaphi1(nb_phi, nb_phi, nb_cos)
      REAL*8 RY_thetaphi1(nb_phi, nb_phi, nb_cos)
      REAL*8 RZ_thetaphi1(nb_phi, nb_phi, nb_cos)
      Common/R_thetaPhi/RX_thetaPhi1, RY_thetaPhi1, RZ_thetaPhi1
      
      real*8 a12(5),a2pi(5)
      common/Guass/a12,a2pi

*  param.h  

c  Grid size
      INTEGER mfft1,mfft2,mfft3
      common/mfft/mfft1,mfft2,mfft3

c  parametres variables
      INTEGER nfft1,nfft2,nfft3, nf1, nf2, nf3
      COMMON/PARAM/nfft1,nfft2,nfft3, nf1, nf2, nf3

c  thermodynamic conditions

      real*8 rho_0,n_0, TEMP
      COMMON/DENS/rho_0,n_0, TEMP



      SUBROUTINE sysdef
      use rho
      use pe
      use vext
	use solute_kind
	use control
      implicit none

      include 'param.h'
      include 'system.h'
      include 'numbers.h'
      include 'rho.h'
      include 'iteration.h'
      include 'bridge.h'
      include 'pe.h'
      include 'Vext.h'



C-------------  VARIABLES LOCALES  ----------------------
      integer n, NREAD,error,m

      integer n_theta, n_phi, n_omega
      real*8 theta, phi


      real*8 EPSO,PIEPSO,QUNIT,MUUNIT,Navo,Boltz

      integer kin_temp(nbmax_solute),typ_mol,ntemp
      real*8 sig_temp(nbmax_solute),eps_temp(nbmax_solute)
      real*8 rtemp
      character(1) chtemp
c--- definition des constantes

      gamma = 1.2254167024d0   !  Gamma(3/4)
      ln2 = 0.69314718056d0    !  ln(2)
      pi=acos(-1.0)
      twopi = 2.0d00*pi
      fourpi = 4.0d00*pi
      EPSO = 8.8542d-12
      PIEPSO = 1.11265d-10
      QUNIT = 1.6d-19
      Boltz = 1.38d-23
      Navo = 6.023d23

c--- unite d energie electrostatique (en kJ/mol)

      QFACT = QUNIT**2*1.0d-3*Navo/(PIEPSO*1d-10)

      TEMP = 300.0! K  system temperature
      n_0 = 0.0332891 ! solvent number density (particile/A**3)


      NREAD = 5    !  fichier de donnees

*----------- Read data from input.dat ----------------------
       
	READ(NREAD,*)
        READ(NREAD,*) BoxSize, weight, IVin 
        READ(NREAD,*)
	READ(NREAD,*) nfft1, nfft2, nfft3
        mfft1=nfft1
        mfft2=nfft2
        mfft3=nfft3
       
      ALLOCATE(rho_O_old(nfft1,nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION rho_O_old FAIL!"
      ALLOCATE(rho_H_old(nfft1,nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION rho_H_old FAIL!"
      ALLOCATE(rho_O_new(nfft1,nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION rho_O_new FAIL!"
      ALLOCATE(rho_H_new(nfft1,nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION rho_H_new FAIL!"

      ALLOCATE(fai_pe(nfft1,nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION fai_pe FAIL!"
      ALLOCATE(bd_ax(nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION bd_ax FAIL!"
      ALLOCATE(bd_bx(nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION bd_bx FAIL!"
      ALLOCATE(bd_ay(nfft1,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION bd_ay FAIL!"
      ALLOCATE(bd_by(nfft1,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION bd_by FAIL!"
      ALLOCATE(bd_az(nfft1,nfft2), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION bd_az FAIL!"
      ALLOCATE(bd_bz(nfft1,nfft2), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION bd_bz FAIL!"
      ALLOCATE(dpar(13*(mfft1+mfft2)/2+9), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION dpar FAIL!"
      ALLOCATE(fai_pe1D(mfft1), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION fai_pe1D FAIL!"
      ALLOCATE(ele_rho(nfft1,nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION ele_rho FAIL!"

      ALLOCATE(VextO(nfft1,nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION VextO FAIL!"
      ALLOCATE(VextH(nfft1,nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION VextH FAIL!"
      ALLOCATE(VLJO(nfft1,nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION VLJO FAIL!"
      ALLOCATE(VintO(nfft1,nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION VintO FAIL!"
      ALLOCATE(VintH(nfft1,nfft2,nfft3), stat = error)
      IF(ERROR.NE.0) STOP "ALLOCATION VintH FAIL!"

       KBT = Boltz*Navo*TEMP*1.0d-3
       beta = one/kBT

       nf1 = nfft1/2
       nf2 = nfft2/2
       nf3 = nfft3/2

       Lg = BoxSize
       Deltav = BoxSize**3/Float(nfft1*nfft2*nfft3)

       write(*,*)"***" 
       write(*,*) 'BoxSize =', BoxSize, "Angstrom."
       write(*,*) 'KBT =', KBT, "KJ/Mol"
       write(*,*) "nfft1= ", nfft1

        READ(NREAD,*)
	READ(NREAD,*) IB, ratio

	Read(NREAD,*)
	Read(NREAD,*)in_dens,one_step,out_dens,scale_md
	Read(NREAD,*)
	Read(NREAD,*)md_peak,rhs_cut,vext_cut,in_dens_cut



         !IVin : option of including Vint: 0 for no, 1 for yes

       rho_0 = n_0 

*------------Read solvent (linear molecule) ---------------
      READ(NREAD,*)
      READ(NREAD,*) nb_solvent_sites
      READ(NREAD,*)
      do n=1,nb_solvent_sites
        READ(NREAD,*) chg_solv(n),sig_solv(n),eps_solv(n),
     &                             x_solv(n),y_solv(n),z_solv(n)
      end do 


      READ(NREAD,*)
*----------- Read correction------------------------------

*----------- Read solute ---------------------------------
      READ(NREAD,*) 
       

      QS = 0.0d0
      READ(NREAD,*) typ_mol
      READ(NREAD,*)
      do n=1,typ_mol
        READ(NREAD,*)kin_temp(n),sig_temp(n),eps_temp(n)
        
      enddo
      READ(NREAD,*)
      READ(NREAD,*) nb_solute_sites
      READ(NREAD,*)
      do n=1,nb_solute_sites
       READ(NREAD,*) ntemp,x_mol(n), y_mol(n), z_mol(n),
     &           rtemp,chg_mol(n)

         x_mol(n) = x_mol(n)-xcor + Lg/2.0
         y_mol(n) = y_mol(n)-ycor + Lg/2.0
         z_mol(n) = z_mol(n)-zcor + Lg/2.0
     
         do m=1,typ_mol
            if(ntemp.eq.kin_temp(m))then
              sig_mol(n)=sig_temp(m)
              eps_mol(n)=eps_temp(m)
            endif
         enddo
         
       write(*,*)n,sig_mol(n),eps_mol(n),chg_mol(n),
     &           x_mol(n),y_mol(n),z_mol(n)
         QS = QS + chg_mol(n)
      end do

      READ(NREAD,*)
      READ(NREAD,*)nb_grd,bulk_cut
      READ(NREAD,*)
      do n=1,nb_grd
        READ(NREAD,*)grd_sl(n),grd_cut(n)
      enddo
      grd_cut(nb_grd+1)=bulk_cut

      READ(NREAD,*)
      READ(NREAD,*)nb_fre
      READ(NREAD,*)
      do n=1,nb_fre
        READ(NREAD,*)fre_sl(n),fre_cut(n)
      enddo
      fre_cut(nb_fre+1)=bulk_cut
	READ(NREAD,*)
	READ(NREAD,*) solute_name
      close(NREAD)
         
        write(*,*) 
      IF(QS.lt.1.0d-6) then
        write(*,*) "  ** NEUTRAL SOLUTE SOLVATED ! **"
      ELSE
        write(*,*) "  ** SOLUTE CARRIES CHARGE VALENCE: ", QS, "**"
      END IF
        write(*,*) 


       write(*,*) "  # HS size for effective brdige (Angstrom) ", DP 

c**** the follow part will be applied in bridge calculation.
        PACKING_F = PI / 6.0 * rho_0 * DP ** 3

! HS BULK CHEMICAL POTENTIAL CALCULATED WITH CARNAHAN-STARLING EOS 
        CHEM_P = 8.0 - 9.0 * PACKING_F + 3.0 * PACKING_F ** 2
        CHEM_P = CHEM_P * PACKING_F / (1.0 - PACKING_F)**3
CC  UNIT: KBT, NAMELY, CHEMICAL POTENTIAL = CHEM_P * KBT
        write(*,*) "   (BULK HS CHEMICAL POTENTIAL/KT = :", CHEM_P,")"
        FBULK = 4.0 * PACKING_F - 3.0 * PACKING_F ** 2
       FBULK = FBULK * RHO_0 * lg**3
       FBULK = FBULK / (1.0d0 - PACKING_F)**2

	 min_de=1E10
	 re_sfe=0.0
      RETURN
      END

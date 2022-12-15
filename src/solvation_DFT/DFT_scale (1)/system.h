*  system.h
*----------

       REAL*8 Lg,Boxsize,DeltaV
	   COMMON/BOX/Lg,Boxsize,DeltaV

      INTEGER nbmax_solute,nbmax_solvent
      PARAMETER(nbmax_solute=100,nbmax_solvent=3)

	   INTEGER nb_solute_sites, nb_solvent_sites
	   REAL*8 x_mol,y_mol,z_mol
           REAL*8 x_solv, y_solv, z_solv
           Real*8 chg_mol,sig_mol,eps_mol
	   REAL*8 chg_solv,sig_solv,eps_solv,d_solv
           REAL*8 Rc,xcor,ycor,zcor

        COMMON/SOLUTE/x_mol(nbmax_solute),y_mol(nbmax_solute),
     &                z_mol(NBMAX_solute),chg_mol(NBMAX_solute),
     &                sig_mol(NBMAX_solute),eps_mol(NBMAX_solute),
     &                Rc,xcor,ycor,zcor,nb_solute_sites
        COMMON/SOLVENT/chg_solv(NBMAX_solvent),sig_solv(NBMAX_solvent),
     &                   eps_solv(NBMAX_solvent),x_solv(NBMAX_solvent),
     &                   y_solv(NBMAX_solvent),z_solv(NBMAX_solvent),
     &                   nb_solvent_sites
 
 
        logical final
        common/out/final
        
        real*8 sfe_tot,detan
        common/sfe/sfe_tot,detan
		
		logical scale_md
		real*8 md_peak
		common/md/scale_md,md_peak
		
		real*8 in_dens_cut
		common/cut_off/in_dens_cut

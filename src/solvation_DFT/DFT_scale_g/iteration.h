* iteration.h

       REAL*8 Torre
       PARAMETER(Torre=0.001)
 
       INTEGER CONVERGE, conv_id
       COMMON/converge/converge, conv_id

       integer IVIN
       common/switchofvint/IVin

       INTEGER NITER
       COMMON/ITERATION NUMBER/NITER

       REAL*8 WEIGHT   ! iteration weight on old density profile
       real*8 weight2
       COMMON/WEIGHT/WEIGHT,weight2

       REAL*8 rho_factor        ! weight on bulk density for iteration trick 
       common/iteration_trick/rho_factor

       REAL*8 RATIO
       COMMON/RATIO_VINT/RATIO

       real*8 srmd_old,srmax_old,fmx,fmx_old
       integer speed
       common/srmd_old/srmd_old,srmax_old,fmx,fmx_old,speed

       real*8 min_de,re_sfe
       common/fd/min_de,re_sfe

       
       integer nb_grd
       integer grd_sl(20)
       real*8 bulk_cut
       real*8 grd_cut(21)
       common/bulkcut/bulk_cut
       common/grid_control/nb_grd
       common/grid_zone/grd_sl,grd_cut
       
       integer nb_fre
       integer fre_sl(20)
       real*8 fre_cut(21)
       common/fre_control/nb_fre
       common/frequency_zone/fre_sl,fre_cut
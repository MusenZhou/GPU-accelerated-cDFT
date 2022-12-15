*c.h

       !real*8 crOO(0:200000), crOH(0:200000), crHH(0:200000)
       !real*8 delta_cr
       !COMMON/DCFs_water/crOO, crOH, crHH, delta_cr

       real*8 crHS(0:2000), delta_crHS
       COMMON/DCF_HS/crHS, delta_crHS

       real*8 ckOO(0:40000), ckOH(0:40000), ckHH(0:40000)
       real*8 delta_ck
       INTEGER nb_ck
       COMMON/DCFs_water/ckOO, ckOH, ckHH, delta_ck, nb_ck

       real*8 ckHS(0:2000)
       real*8 delta_ckHS
       integer nb_ckHS
       COMMON/DCF_HS/ckHS, delta_ckHS, nb_ckHS



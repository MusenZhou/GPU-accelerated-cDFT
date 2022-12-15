* bridge.h

      REAL*8 DP   ! diameter of HS particle
      PARAMETER (DP = 2.86D0)

      INTEGER IB,statB        ! option of including bridge, 1 for including bridge and 0 for not including bridge
      common/option_bridge/IB,statB

      real*8 PACKING_F, CHEM_P,FBULK
      common/HS_bulk_excesschemP/PACKING_F, CHEM_P,FBULK

     

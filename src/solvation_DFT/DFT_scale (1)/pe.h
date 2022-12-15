* pe.h
*---------    
      
      real*8 zh, zo
      PARAMETER(zh=1.00, zo = -2.0)
       
      real*8 QH, QO
      parameter(QH = 0.416d0)
      parameter(QO = -2.0 * QH )

      real*8 QS
      common/Solute charge/QS


      real*8 ele_tot 
      common/elerho/ele_tot
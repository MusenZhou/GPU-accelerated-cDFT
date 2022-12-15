****** donne le temps cpu en secondes ************
      SUBROUTINE temps(ZTIME)
      REAL*8 ZTIME
      dimension T(2)
       time=ETIME(T)
       ZTIME=dble(T(1))
C     ztime=0.0
      RETURN
      END

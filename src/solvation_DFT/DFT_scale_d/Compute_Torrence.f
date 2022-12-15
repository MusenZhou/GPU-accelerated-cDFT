cc*********************************************************c
cc THIS SUBROUTINE IS TO CALCULATE THE DIFFERENCE BETWEEN **
CC NEW DENSITY PROFILE AND OLD DENSITY PROFILE.           **
CC                                S.L. ZHAO               **
CC**********************************************************


       SUBROUTINE COMPUTE_Torrence
       use vext
       use rho
       IMPLICIT NONE

       include 'param.h'
       include 'system.h'
       include 'numbers.h'
       INCLUDE 'Vext.h'
       INCLUDE 'rho.h'
       INCLUDE 'iteration.h'


       REAL*8 srmdO, srmdH    ! squre root of mean deviation
       REAL*8 srmd, srmax
       INTEGER i, j, k
       CHARACTER*40 file_O, file_H 
       CHARACTER NAM*7, CH*3
       CHARACTER*11 chemin
       REAL*8 X_COM, Y_COM, Z_COM, RR
       real*8 NumO, NumH,NumOH,ep1
       integer abs_sp
       real*8 rpo_new,rne_new,rpo_old,rne_old,rO_H,derate
       real*8 temp_o,temp_h,derho
       logical con_trend
       real*8 w_old
       
       real*8 rele_cut
       real*8 w_back,w_cut
       real*8 rr_p,cor
       
       rele_cut=(LG/2.0-1.0)

       con_trend=.true.



cc Iteration judgement

        srmdO = 0.0
        srmdH = 0.0




        !open(40,file = "./rho_O_updated.dat") 
        !open(41,file = "./rho_H_updated.dat") 

        NumO = 0.0d0
        NumH = 0.0d0
        NumOH=0.0D0

        srmax = 0.0d0
        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        DO k = 1, nfft3
         DO j = 1, nfft2
          DO i = 1, nfft1

           x_com = (i-1)*Lg/FLOAT(nfft1)-lg/2.0
           y_com = (j-1)*Lg/FLOAT(nfft2)-lg/2.0
           z_com = (k-1)*Lg/FLOAT(nfft3)-lg/2.0

           rr = DSQRT(x_com**2 + y_com**2 + z_com**2)
          


           srmdO = srmdO + (rho_O_new(i,j,k)-rho_O_old(i,j,k))**2
           srmdH = srmdH + (rho_H_new(i,j,k)-rho_H_old(i,j,k))**2
           if(abs(rho_O_new(i,j,k)-rho_O_old(i,j,k)).gt.srmax)then
             srmax=abs(rho_O_new(i,j,k)-rho_O_old(i,j,k))
             fmx=rho_O_new(i,j,k)-rho_O_old(i,j,k)
             rr_p=rr
           endif
           if(abs(rho_H_new(i,j,k)-rho_H_old(i,j,k)).gt.srmax)then
             srmax=abs(rho_H_new(i,j,k)-rho_H_old(i,j,k))
             fmx=rho_H_new(i,j,k)-rho_H_old(i,j,k)
             rr_p=rr
           endif



          enddo
         enddo
        enddo
        srmdO =sqrt(srmdO/float(nfft1*nfft2*nfft3))
        srmdH =sqrt(srmdH/float(nfft1*nfft2*nfft3))
        srmd = srmdO + srmdH
        write(*,*) '$$ square root mean deviation:',srmd
        write(*,*) '$$ maximum deviation:',fmx,rr_p!srmax

	  if(srmax.lt.min_de)then
	    min_de=srmax
	    re_sfe=sfe_tot
	  endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
c        weight=0.95
         IF(niter.eq.1) then
           srmd_old = srmd + 1.0
           srmax_old=srmax+1.0
           speed = 0
         END IF
         w_cut=0.95
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if(ratio.gt.0.02)then 
c           IVIN=1
           IF((srmd.gt.srmd_old).or.(srmax.gt.srmax_old)) then
             w_old=weight
             weight=1.0D0-(1.0D0-weight)*0.1D0
             con_trend=.false.
           else
             weight=w_cut*0.1+weight*0.9
c             weight=1.0D0-(1.0D0-weight)*weight/w_cut
           endif
           if(weight.ge.0.999)then 
             weight =0.999
           endif
           if(weight.lt.w_cut)then 
             weight =w_cut
           endif
         endif

c           weight=0.95        
c         if(srmax.gt.20.0)then
c           weight=0.999
c         endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        DO k = 1, nfft3
         DO j = 1, nfft2
          DO i = 1, nfft1


c           srmdO = srmdO + (rho_O_new(i,j,k)-rho_O_old(i,j,k))**2
c           srmdH = srmdH + (rho_H_new(i,j,k)-rho_H_old(i,j,k))**2
c           if(abs(rho_O_new(i,j,k)-rho_O_old(i,j,k)).gt.srmax)then
c             srmax=abs(rho_O_new(i,j,k)-rho_O_old(i,j,k))
c           endif
c           if(abs(rho_H_new(i,j,k)-rho_H_old(i,j,k)).gt.srmax)then
c             srmax=abs(rho_H_new(i,j,k)-rho_H_old(i,j,k))
c           endif
           

           x_com = (i-1)*Lg/FLOAT(nfft1)-lg/2.0
           y_com = (j-1)*Lg/FLOAT(nfft2)-lg/2.0
           z_com = (k-1)*Lg/FLOAT(nfft3)-lg/2.0

           rr = DSQRT(x_com**2 + y_com**2 + z_com**2)
                 
c        if(ratio.gt.0.5)then
c           derho=1.0D-7
c         if((rho_O_old(i,j,k).lt.0.1))then
c           rho_O_new(i,j,k) = rho_O_new(i,j,k)*(1.0-weight) 
c     &                      + rho_O_old(i,j,k)*weight
c         else
c           rho_O_new(i,j,k) = (rho_O_new(i,j,k)*rho_O_old(i,j,k)
c     &     +derho*rho_O_old(i,j,k))/(rho_O_new(i,j,k)*weight
c     &     + rho_O_old(i,j,k)*(1.-weight)+derho)
c         endif 
c         if(rho_H_old(i,j,k).lt.0.1)then
c           rho_H_new(i,j,k) = rho_H_new(i,j,k)*(1.0-weight) 
c     &                      + rho_H_old(i,j,k)*weight
c         else
c           rho_H_new(i,j,k) = (rho_H_new(i,j,k)*rho_H_old(i,j,k)
c     &     +derho*rho_H_old(i,j,k))/(rho_H_new(i,j,k)*weight
c     &     + rho_H_old(i,j,k)*(1.-weight)+derho)
c         endif   
c       endif
c         IF(rr.lt.(LG/2.0-4.0)) THEN
c          if(niter.gt.1)then
           rho_O_new(i,j,k) = rho_O_new(i,j,k)*(1.0d0-weight) 
     &                      + rho_O_old(i,j,k)*weight

           rho_H_new(i,j,k) = rho_H_new(i,j,k)*(1.0d0-weight) 
     &                      + rho_H_old(i,j,k)*weight
c          else
c           rho_O_new(i,j,k) = rho_O_new(i,j,k)*0.003 
c     &                      + rho_O_old(i,j,k)*0.997

c           rho_H_new(i,j,k) = rho_H_new(i,j,k)*0.003 
c     &                      + rho_H_old(i,j,k)*0.997
c          endif
c         else
c           rho_O_new(i,j,k) =1.0D0
c           rho_H_new(i,j,k) =1.0D0
c         endif
         



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          if((niter.gt.0))then
          rO_H=1.

          rpo_new=rho_H_new(i,j,k)+rO_H*rho_O_new(i,j,k)
          rne_new=rho_O_new(i,j,k)-rho_H_new(i,j,k)
          rpo_old=rho_H_old(i,j,k)+rO_H*rho_O_old(i,j,k)
          rne_old=rho_O_old(i,j,k)-rho_H_old(i,j,k)


c          derate = 0.1
          derate=1.0/(1.0+rr**2)
	    if(derate.gt.0.1)then
	      derate=0.1
	    endif
          if(rpo_new.gt.(rpo_old+derate))then
             rpo_new=rpo_old+derate
          endif
          if(rpo_new.lt.(rpo_old-derate))then
             rpo_new=rpo_old-derate     
          endif

          if(abs(rne_new).gt.abs(rne_old))then
c           derate = 0.05
           derate=derate*0.5
           if(rne_new.gt.(rne_old+derate))then
             rne_new=rne_old+derate
           endif
           if(rne_new.lt.(rne_old-derate))then
            rne_new=rne_old-derate
           endif
          endif
c          rho_O_old3(i,j,k)=rho_O_old2(i,j,k)
c          rho_H_old3(i,j,k)=rho_H_old2(i,j,k)
c          rho_O_old2(i,j,k)=rho_O_old(i,j,k)
c          rho_H_old2(i,j,k)=rho_H_old(i,j,k)
          rho_O_old(i,j,k)=(rpo_new+rne_new)/(rO_H+1.)
          rho_H_old(i,j,k)=(rpo_new-rO_H*rne_new)/(rO_H+1.)
          else
c          rho_O_old3(i,j,k)=rho_O_old2(i,j,k)
c          rho_H_old3(i,j,k)=rho_H_old2(i,j,k)
c          rho_O_old2(i,j,k)=rho_O_old(i,j,k)
c          rho_H_old2(i,j,k)=rho_H_old(i,j,k)
          rho_O_old(i,j,k)=rho_O_new(i,j,k)
          rho_H_old(i,j,k)=rho_H_new(i,j,k)
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c         ELSE
c           rho_O_old(i,j,k) = 1.0d0 !EXP(-beta* VLJO(i,j,k))   ! g(r)
c           rho_H_old(i,j,k) = 1.0d0 
c         END IF

c        IF(NITER.eq.1) THEN
c          write(30,*) rho_O_old(i,j,k)
c          write(31,*) rho_H_old(i,j,k)
c        END IF

           
c           IF(rr.lt.rele_cut) THEN

c           NumOH = NumOH+rho_O_old(i,j,k)*rho_H_old(i,j,k)
          
c           endif

           END DO
          END DO

** test code ***
c         write(*,*) nf1,nf2,k, rho_O_old(nf1,nf2,k) 
c         write(*,*) nf1,nf2,k, rho_H_old(nf1,nf2,k) 
         
c         write(40,*) (k-1)*Lg/nfft3 - lg/2.0, rho_O_old(nf1,nf2,k) 
c         write(41,*) (k-1)*Lg/nfft3 - lg/2.0, rho_H_old(nf1,nf2,k) 
** test code ***
        END DO
c        close(40)
c        close(41)

        cor=dsqrt(NumH/NumO)
c        IF(NumO.ne.NumH) then

        !DO k = 1, nfft3
        !write(40,*) (k-1)*Lg/nfft3 - lg/2.0, rho_O_old(nf1+1,nf2+1,k) 
        !write(41,*) (k-1)*Lg/nfft3 - lg/2.0, rho_H_old(nf1+1,nf2+1,k)
        !END DO
        !close(40)
        !close(41)

c        END IF


        

c        IF(mod(niter-1, 10).           endif eq.0) THEN
        
c        END IF

        If(srmax.gt.Torre) then
           converge = 0
        ELSE
           ConV_id = 1

c          IF(IVIN.ge.1) THEN  ! including the non-idea part
          IF(ratio.gt.0.99D0) THEN
           converge = 1
          ELSE
           converge = 0
          END IF


          IVIN = IVIN + 1
          niter=1
c          weight = 0.99!99
          ratio = ratio + 0.1d0
          weight=1.0D0-(1.0D0-weight)*0.1D0
        endif
          w_cut=0.2D0+ratio*0.75

         

c       IF(IVIN.ge.1) then
c         IF(srmd.lt.srmd_old) then
c           speed = speed - 1
c         ELSE
c           speed = speed + 1
c         END IF
 
c         abs_sp = ABS(speed)

c         IF(ABS_sp.ge.2) then
c           weight = weight + speed/abs_sp* 9.0/10**(2+ABS_sp)
c         END IF 

c         IF(weight.gt.1.0) then
c           weight = 0.999
c         END IF


          write(*,*) "Picard mixing weight :", weight 

          srmd_old = srmd
          srmax_old = srmax
          fmx_old=fmx


c          de_srmd=srmd-srmd_old
c          if(con_trend)then
          
c          endif

c       END IF


ccccccccccccccccc Update rho cccccccccccccccccccccccc
c       if(con_trend)then
c         DO k = 1, nfft3
c          DO j = 1, nfft2
c           DO i = 1, nfft1
c            rho_O_old(i,j,k)=rho_O_new(i,j,k)
c            rho_H_old(i,j,k)=rho_H_new(i,j,k)
c          enddo
cc         enddo
c         enddo
c       else
c         DO k = 1, nfft3
c          DO j = 1, nfft2
c           DO i = 1, nfft1
c           rho_O_old(i,j,k)=weight*rho_O_old2(i,j,k)
c    &        +0.1D0*(rho_O_old(i,j,k)-w_old*rho_O_old2(i,j,k))
c            rho_H_old(i,j,k)=weight*rho_H_old2(i,j,k)
c     &        +0.1D0*(rho_H_old(i,j,k)-w_old*rho_H_old2(i,j,k))
c           enddo
c          enddo
c         enddo
c       endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
       RETURN
       END


C    *******************************************************************
C    **   LITTLE PROGRAM TO TRANSFER NUMBER TO CHARACTER              **
C    *******************************************************************

          SUBROUTINE CHARTRANS (NUM, CH )
          IMPLICIT NONE

C    *******************************************************************
C    ** TRANSFER NUMBER TO CHARACTERS                                 **
C    *******************************************************************


          INTEGER NUM

          INTEGER   I, J, L
          CHARACTER CH*3


          DO I = 1, 3

             L = MOD ( NUM / ( 10 ** ( I - 1 ) ), 10 )

             J = 4 - I

             CH ( J:J ) = CHAR ( L + 48 )

          ENDDO

          RETURN
          END

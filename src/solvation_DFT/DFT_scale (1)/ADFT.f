      module pe
      real*8,allocatable:: dpar(:)
      real*8,allocatable:: fai_pe(:,:,:)
      real*8,allocatable:: bd_ax(:,:), bd_bx(:,:)
      real*8,allocatable:: bd_ay(:,:), bd_by(:,:)
      real*8,allocatable:: bd_az(:,:), bd_bz(:,:)
      real*8,allocatable:: fai_pe1D(:)
      real*8,allocatable:: ele_rho(:,:,:)
      end module
 
      module rho
      real*8,allocatable:: rho_O_old(:,:,:)
       real*8,allocatable:: rho_H_old(:,:,:)

       real*8,allocatable:: rho_O_new(:,:,:)
       real*8,allocatable:: rho_H_new(:,:,:)
      end module

      module vext
      real*8,allocatable:: VextO(:,:,:),VextH(:,:,:),VLJO(:,:,:)
      real*8,allocatable:: VintO(:,:,:),VintH(:,:,:)
	  real*8 rhs_cut,vext_cut
      end module

      module solute_kind
      Character*100 solute_name 
      end module

     	module control
	Logical in_dens,one_step,out_dens
	end module


      PROGRAM ADFT_PROFILE
      use pe
      use rho
      use vext
	use control
      IMPLICIT NONE

      include 'param.h'
      include 'system.h'
      include 'numbers.h'
      include 'c.h'
      include 'rho.h'
      include 'iteration.h'
      include 'bridge.h'


C*** LOCAL VARIABLE ***
      REAL*8 TIME0, TIME1, TIME2, TIME_TOT 

      final=.false.
      statB=0
c*************************************************
c************** PART 1: READING INPUT  ***********

       CALL sysdef   ! define the calculation system

*===========================


*  Compute external potential

       CALL Compute_Vext
 

*  Read direct correlation functions (DCFs)
       CALL get_DCFs


       CALL Get_cos_table
c*****************************************************
C******* PART 2: INITIAL GUESS OF DENSITY PROFILE ****

       CALL Compute_initial_rho

        write(*,*) "       @@@           @@@           @@@"
        write(*,*)
        IF(IB.eq.1) THEN
        write(*,*) "      && WARNING: BRIDGE INCLUDED! &&"
        ELSE IF(IB.EQ.0) THEN
        write(*,*) "    && WARNING: BRIDGE not INCLUDED! &&"
        END IF
       
        NITER = 1
        ConV_id = 0   ! ideal converge
        Converge = 0

c****************************************************
C******** PART 3: UPDATE THE DENSITY PROFILE  *******
333       CONTINUE    ! iteration loops

       CALL TEMPS(TIME1)
         
       IF(NITER.eq.1) TIME0 = TIME1

c       IF(NITER.gt.1) STOP

       if(mod(niter,10).eq.1)then
         statB=1
       else
         statB=0
       endif

       CALL Compute_new_rho

	if(one_step)then
	goto 444
	end if
c****************************************************
C********** PART 4: ITERATION JUDGEMENT *************

        CALL COMPUTE_Torrence


        CALL TEMPS(TIME2)
	  TIME_TOT=	(TIME2-TIME0)/60.

        IF((Converge.EQ.0).and.(TIME_TOT.lt.360.0)) THEN
          niter=niter+1
          GOTO 333
        ELSE
          write(*,*) 
          write(*,*)  " - - -" 
          write(*,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" 
          write(*,*) "*" 
          write(*,*) "$$ Elapsed time in total:",INT(TIME_TOT),
     & "Mins."
          write(*,*) "*" 
        END IF




c****************************************************
C********* PART 5: OUTPUT FINAL RESULT  ************
444	  CONTINUE
        final=.true.
        call final_output
		
		if(scale_md.and.(.not.onestep))then
		  one_step=.true.
		  goto 333
		endif


        call free_space



      WRITE(*,*)
      WRITE(*,*) "  #  NORMAL END !"
      WRITE(*,*)

      STOP
      END


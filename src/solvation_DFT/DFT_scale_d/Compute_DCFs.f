c************************************************************c
c** THIS SUBROUTINE IS TO INPUT THE DIRECT CORRELATIONS OF **
C** BULK WATER, AND THAT OF BULK HARD SPHERE.              **
C**                              S.L. ZHAO                 **
c** Reference for DCFs of bulk SPC/E water:                **
c**   S.L. Zhao and J. Wu, Molecular Physcs, V109, P2553   **
c**   (2011)                                               **
c** Reference for DCF of bulk hard-sphere:                 **
c**   use Modified FMT.                                    **
c************************************************************


        SUBROUTINE get_DCFS
        IMPLICIT NONE
        
        INCLUDE 'c.h'

* local variables
        INTEGER nhist
        REAL*8 rr

        INTEGER nk
        REAL*8 norm_k

c Input the DCFs of bulk SPC/E water 


c*********** Real space *********

c       open(11,file="DCF/cr_OO_final_highres.dat")
c       open(12,file="DCF/cr_OH_final_highres.dat")
c       open(13,file="DCF/cr_HH_final_highres.dat")

c       open(20,file="../Results/cr_tot_100.dat")
c       open(31,file="../Results/cr_SR.dat")

c       do nhist = 1, 119624
c          read(11,*) rr, crOO(nhist)
c          read(12,*) rr, crOH(nhist)
c          read(13,*) rr, crHH(nhist)
 
c          write(20,*) rr,crOO(nhist)+4.0*crOH(nhist)+4.0*crHH(nhist)

c          If(rr.lt.10.0) then
c          write(31,*) rr, crOO(nhist)
c          ELSE
c          write(31,*) rr, crOO(nhist) + 96.22/rr * 4.0
c          END IF

c       end do
c       CLOSE(11)
c       CLOSE(12)
c       CLOSE(13)

c       CLOSE(20)
c       CLOSE(31)

c       delta_cr = 4.98230286115115550E-002
c       stop "stop in DCFS."
            

c Input the DCFs of bulk hard sphere 

c       open(16, file="DCF/crHS_d265.dat")

c       DO nhist = 0, 2000
c           read(16,*) rr, crHS(nhist)
c       END DO
c       close(16)

c       delta_crHS = 0.01d0

c*********** Fourier space *********
c   total ck
c        delta_ck = 2.10543812536071219E-002
c        nb_ck = 1550

c       open(21,file="DCF/ck_OO_final.dat")
c       open(22,file="DCF/ck_OH_final.dat")
c       open(23,file="DCF/ck_HH_final.dat")
c       open(21,file="DCF/ckOO_convert.dat")
c         open(22,file="DCF/ckOH_convert.dat")
c         open(23,file="DCF/ckHH_convert.dat")

c short range ck

c        delta_ck = 2.43530273437500000E-002
        
c        nb_ck = 1638
c       open(21,file="DCF/ckOO_site_SR1.2.dat")
c       open(22,file="DCF/ckOH_site_SR1.2.dat")
c       open(23,file="DCF/ckHH_site_SR1.2.dat")

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        delta_ck = 0.006185407D0
        
        nb_ck = 6400
       open(21,file="../DCF/ck_OO_TIP3P_SR1.0.dat")
       open(22,file="../DCF/ck_OH_TIP3P_SR1.0.dat")
       open(23,file="../DCF/ck_HH_TIP3P_SR1.0.dat")
c       open(21,file="DCF/cksr_OO.dat")
c       open(22,file="DCF/cksr_OH.dat")
c       open(23,file="DCF/cksr_HH.dat")
       !kapa=5.0
c       delta_ck =4.90873865783214569E-003
c       nb_ck=7000
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       do nk=1, nb_ck

        read(21,*) norm_k, ckOO(nk)
        read(22,*) norm_k, ckOH(nk)
        read(23,*) norm_k, ckHH(nk)
c        ckoo(nk)=ckOO(nk)+4836.4/norm_k**2
c        ckoh(nk)=ckoh(nk)-2418.2/norm_k**2
c        ckhh(nk)=ckhh(nk)+1209.1/norm_k**2
       end do
        ckoo(0)=2*ckoo(1)-ckoo(2)
        ckoh(0)=2*ckoh(1)-ckoh(2)
        ckhh(0)=2*ckhh(1)-ckhh(2)
       close(21)
       close(22)
       close(23)



      delta_ckHS = 0.05d0
      nb_ckHS = 1000

c      nb_ckHS = 600
c      open(21,file="DCF/ck_FMT_d2.85.dat")


c      open(21,file="DCF/ck_FMT_d3.00.dat")
c      open(21,file="DCF/ck_FMT_d2.964.dat")
      open(21,file="../DCF/ck_FMT_d2.86.dat")

       do nk=0, nb_ckHS

        read(21,*) norm_k, ckHS(nk)

       end do
       close(21)

  
       write(*,*) "  # Input the DCFs. " 

        RETURN
        END

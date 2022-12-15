c************************************************************c
c* THIS SUBROUTINE IS TO CALCULATE THE EXTERNAL POTENTIAL  **
C* EXERTED ON EACH SITE DUE TO THE PRESENCE OF SOLUTE.     **
C*                              S.L. ZHAO                  **
C************************************************************


      SUBROUTINE Compute_Vext
      use vext
      implicit none

      include 'param.h'
      include 'fft.h'
      include 'system.h'
      include 'numbers.h'
      include 'Vext.h'

*-------------  LOCAL VARIABLES   --------------
      INTEGER i, j, k, n, m

      REAL*8 x_com, y_com, z_com, x_nm, y_nm, z_nm, x_m, y_m, z_m

      REAL*8 r_nm,r_nm2, r_nm6, sig_nm, sig_nm6, eps_nm, V_LJ, V_coul

      REAL*8 F,uljmin,uljtemp,blr,vadd
     
      logical in_hs


*BEGIN
      open(43, file ="../Results/VextO_Z.dat")
      open(44, file ="../Results/VextH_Z.dat")
c---- Compute the solute-solvent interaction ----

      DO  k=1,nfft3
        DO  j=1,nfft2
         DO  i=1,nfft1

** solvent site position
           x_com = (i-1)*Lg/FLOAT(nfft1)
           y_com = (j-1)*Lg/FLOAT(nfft2)
           z_com = (k-1)*Lg/FLOAT(nfft3)

         do m = 1, 2 !nb_solvent_sites  ! solvent

          V_LJ = zero
          V_coul = zero
          in_hs=.false.
          vadd=0.0D0

            do n = 1, nb_solute_sites  ! solute

c  Compute solute-solvent site distance

                x_nm = x_com - x_mol(n) 
                y_nm = y_com - y_mol(n) 
                z_nm = z_com - z_mol(n) 

                r_nm2 = x_nm**2+y_nm**2+z_nm**2
                r_nm = sqrt(r_nm2)

                if(r_nm.lt.1.0d-1) then
c                  write(*,*) 'ERROR: r_nm=0 for', i, j, k
                  r_nm = 1.0d-1
                  r_nm2 = r_nm**2
c                  vadd=5000.0
                end if
          if((r_nm.lt.1.4d0).and.(chg_mol(n)*chg_solv(m).lt.0.0))then
                  vadd=5000.0
          else
             vadd=0.0
          endif
                r_nm6 = r_nm2*r_nm2*r_nm2

c Lorentz combining rule
                eps_nm = sqrt(eps_mol(n)*eps_solv(m))
                sig_nm = (sig_mol(n)+sig_solv(m))/two

c                if(r_nm.lt.sig_nm)then
c                   V_LJ=800.0D0
c                else
c                   V_LJ=0.0D0
c                endif

c                eps_nm = eps_mol(n)
c                sig_nm = sig_mol(n)

                sig_nm6 = sig_nm*sig_nm*sig_nm
                sig_nm6 = sig_nm6*sig_nm6/r_nm6

                uljtemp=four*eps_nm*(sig_nm6*sig_nm6 - sig_nm6)

                V_LJ = V_LJ + uljtemp
             if(vadd.lt.1.0)then
                V_coul = V_coul + QFACT*chg_mol(n)*chg_solv(m)/r_nm
             endif

            end do   ! end solute
 
            IF(m.eq.1) then
             VextO(i,j,k) = V_LJ + V_coul+vadd
             
             VLJO(i,j,k) = V_LJ 
c             blr=-0.036*QFACT*chg_solv(m)*chg_solv(m)/r_nm
c             VextO(i,j,k)=VextO(i,j,k)+blr
c** test code
c            write(*,*) i,j,k,V_LJ
c** test code

            ELSE if(m.eq.2) then
             VextH(i,j,k) =   V_coul+vadd
c             blr=-0.036*QFACT*chg_solv(m)*chg_solv(m)/r_nm
c             VextH(i,j,k)=VextH(i,j,k)+blr
            END IF

           end do ! end solvent loop


          END DO 
          END DO 

cc  test code **
          write(43,*) (k-1)*lg/nfft3-lg/2.0, Beta*VextO(nf1,nf2,k)
          write(44,*) (k-1)*lg/nfft3-lg/2.0, Beta*VextH(nf1,nf2,k)
cc  test code **

          END DO  ! end the position on grid 
          close(43)
          close(44)

        write(*,*)"  # Calculate external potential."


       RETURN
       END






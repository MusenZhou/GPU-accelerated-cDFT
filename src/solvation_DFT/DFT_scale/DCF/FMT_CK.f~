CC  THIS IS THE FINAL PROGRAM TO CALCULATE THE C(K) FOR HARD SPHERE SYSTEM

        program main
!     AN0,AN1,AN2,AN3:  WEIGHTED DENSITIES FOR BULK
	INTEGER NR   ! kmax = NR * DK
	PARAMETER (NR=1001,PI=3.141592654,NK=1001)
	REAL,PARAMETER::DK=0.05d0  ! real step length 1/A
	REAL SIGMA(2),RHO(2),ETA(2),C(2,2,NR),RADI(2),CK(2,2,NK)
	REAL AN0,AN1,AN2,AN3
	REAL OM0(2),OM1(2),OM2(2),OM3(2),OMV1(2),OMV2(2),S11(NK),S22(NK)
        character*16 file_ckHS

       SIGMA(1)=2.92d0  ! hard sphere size
       file_ckHS = 'ck_FMT_d2.92.dat'  ! output result

!        file_ckHS = 'ck_FMT_d2.65.dat'

       write(*,*) 
       write(*,*) "HS size: ", sigma(1)

	SIGMA(2)=1.0

	RADI(1)=SIGMA(1)/2.0
	RADI(2)=SIGMA(2)/2.0

        RHO(1)= 0.033289050 !0.45/SIGMA(1)**3

 	RHO(2)= 0.0  !0.45/SIGMA(2)**3
 	ETA(1)=RHO(1)*PI*SIGMA(1)**3/6.
 	ETA(2)=RHO(2)*PI*SIGMA(2)**3/6.

C	ETA(1)=0.32805
c	ETA(2)=0.5
C       RHO(1)=ETA(1)*6./(PI*SIGMA(1)**3)
c	RHO(2)=ETA(2)*6./(PI*SIGMA(2)**3)

        WRITE(*,*)"Kmax =", NR*DK
        WRITE(*,*)"RHO =", RHO(1)

	ETAT=ETA(1)+ETA(2)
	RHOT=RHO(1)+RHO(2)

	DR=0.01
	
	AN0=RHO(1)+RHO(2)
	AN1=0.5*(RHO(1)*SIGMA(1)+RHO(2)*SIGMA(2))
	AN2=PI*(RHO(1)*SIGMA(1)**2+RHO(2)*SIGMA(2)**2)
	AN3=PI*(RHO(1)*SIGMA(1)**3+RHO(2)*SIGMA(2)**3)/6.0

        CHI0=1./(1.-AN3)
	CHI1=AN2/(1.-AN3)**2
	CHI2=AN1/(1.-AN3)**2-AN2**2/(12.*PI)
     &       *(2.*ALOG(1.-AN3)/AN3**3
     &       +(2.-5.*AN3+AN3**2)/((1.-AN3)**3*AN3**2))
	CHI3=AN0/(1.-AN3)**2+2.*AN1*AN2/(1.-AN3)**3
     &       +AN2**3/(36.*PI*AN3**4)*(6.*ALOG(1.-AN3)
     &       +(6.-21.*AN3+26.*AN3**2-5.*AN3**3)*AN3/(1.-AN3)**4)
      CHI22=AN2/(6.*PI*AN3**2)*(ALOG(1.-AN3)+AN3/(1.-AN3)**2)

	DO 100 I=1,2
	DO 100 J=1,2
	RI=RADI(I)
	RJ=RADI(J)
	RIJ=RI+RJ
	DO 100 K=1,NR
	R=FLOAT(K-1)*DR
	IF(R.LT.RIJ) THEN
	CALL PYDCF(RI,RJ,R,OV,OS,OR,OX)
	C(I,J,K)=-CHI3*OV-CHI2*OS-CHI1*OR-CHI0
     &         -(CHI22-CHI1/(4.*PI))*OX
	ELSE
	C(I,J,K)=0.0
	ENDIF
100	CONTINUE

        OPEN(10, file='OM.dat')

	DO 200 K=1,NK
	T=FLOAT(K-1)*DK 
	IF (K.EQ.1) THEN
		DO J=1,2
		OM0(J)=1.0
		OM1(J)=SIGMA(J)/2.0
		OM2(J)=PI*SIGMA(J)**2
		OM3(J)=PI*SIGMA(J)**3/6.0
		OMV1(J)=0.0
		OMV2(J)=0.0
		ENDDO
      ELSE
		DO J=1,2
		OM0(J)=2.*SIN(T*SIGMA(J)/2.0)/(T*SIGMA(J))
		OM1(J)=SIN(T*SIGMA(J)/2.0)/T
		OM2(J)=2.*PI*SIGMA(J)*SIN(T*SIGMA(J)/2.0)/T
		OM3(J)=4.*PI*(SIN(T*SIGMA(J)/2.0)
     &           -T*SIGMA(J)/2.0*COS(T*SIGMA(J)/2.0))/T**3
		OMV1(J)=-T*OM3(J)/(2.*PI*SIGMA(J))
		OMV2(J)=-T*OM3(J)
		ENDDO	
	ENDIF

	DO I=1,2
	DO J=1,2
	CK(I,J,K)=-CHI3*OM3(I)*OM3(J)
     &          -CHI2*(OM3(I)*OM2(J)+OM2(I)*OM3(J))
     &          -CHI1*(OM3(I)*OM1(J)+OM1(I)*OM3(J))
     &          -CHI22*(OM2(I)*OM2(J)-OMV2(I)*OMV2(J))
     &          -CHI0*(OM3(I)*OM0(J)+OM0(I)*OM3(J)+
     &                 OM2(I)*OM1(J)+OM1(I)*OM2(J)-
     &                 OMV2(I)*OMV1(J)-OMV1(I)*OMV2(J))

	ENDDO
	ENDDO

!       write(*,*) T,CK(1,1,K) 
!       read(*,*)

       write(10,*) T, OM0(1),OM1(1),OM2(1),OM3(1)
!       write(*,*) T, OM0(1),OM1(1),OM2(1),OM3(1)
!       write(*,*) T, CHI0, CHI1,CHI2,CHI3,CHI22
!       write(*,*) T, CK(1,1,k)
!       read(*,*)


200   CONTINUE 
       close(10)

	OPEN (4,FILE='cr_FMT.dat')	
	OPEN (5,FILE=file_ckHS)
	DO K=1,NR
	R=FLOAT(K-1)*DR/SIGMA(1)
	XK=FLOAT(K-1)*DK!*SIGMA(1)
	S11(K)=1./(1.-RHO(1)*CK(1,1,K))  ! modification here

	S11(K)=CK(1,1,K)

c       write(*,*) XK,S11(K) 

c       read(*,*)
c	S22(K)=1./(1.-RHO(2)*CK(2,2,K))
	S22(K)=CK(1,1,K)/(1.-RHO(1)*CK(1,1,K))
C	WRITE(4,300) R, C(1,1,K),C(2,1,K),C(1,2,K),C(2,2,K),chs1(etaT,r) 
C       WRITE(4,300) XK, S11(K),S22(K),R,C(1,1,K),C(1,2,K),C(2,2,K)
	WRITE(4,*) XK, S11(K),S22(K)
	WRITE(5,*) XK, S11(K)
	ENDDO

	CLOSE(4)
	CLOSE(5)

       write(*,*) "file saved in ", file_ckHS
       write(*,*)

300 	FORMAT(2X,F10.4,2F12.4,2X,F10.6,3F12.3)
	END 

! calculate the overlap volume(OV),the overlap surface area(OS), the mean radius of sphercone (OR)
! between two spheres of radii RI and RJ separated by R
	SUBROUTINE PYDCF(RI,RJ,R,OV,OS,OR,OX)	
	IMPLICIT NONE
	REAL OV,OS,OR,OX,RI,RJ,R,PI,LAMBDA,RIJ,RMIN,X
	PI=3.141592654
	RMIN=MIN(RI,RJ)
	LAMBDA=ABS(RI-RJ)
	RIJ=RI+RJ
	IF (R.LE.LAMBDA) THEN
	OV=4.*PI*RMIN**3/3.
	OS=4.*PI*RMIN**2
	OR=RMIN
	OX=0.0
	ELSEIF ((R.GT.LAMBDA).AND.(R.LE.RIJ)) THEN
	X=-1.5*(RI**2-RJ**2)**2+4.*R*(RI**3+RJ**3)
     &  -3.*R**2*(RI**2+RJ**2)+R**4/2.
	OV=X*PI/(6.*R)

	X=-RIJ*LAMBDA**2+2.*R*(RI**2+RJ**2)-R**2*RIJ
	OS=X*PI/R

	X=-LAMBDA**2/4.+R*RIJ/2.-R**2/4.
	OR=X/R

	OX=R*PI-PI*LAMBDA**2/R
	ELSE
	OV=0.0
	OS=0.0
	OR=0.0
	OX=0.0
	ENDIF
	END

	function chs1(eta,r)
	c1=(1.+eta*(4.+eta*(3.-2.*eta)))/(1.-eta)**4
	c2=(2.-eta+14.*eta**2-6.*eta**3)/(1.-eta)**4+2.*alog(1.-eta)/eta
	c3=(3.+5.*eta*(eta-2.)*(1.-eta))/(1.-eta)**4+3.*alog(1.-eta)/eta
	chs1=-c1+c2*r-c3*r**3
	end



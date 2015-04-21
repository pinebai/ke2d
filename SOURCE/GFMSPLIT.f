	SUBROUTINE GFMSPLIT

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		04-03-2014:	CREATED
C	---------------------------------------------------------------
C	SPLIT THE FIELDS USING A GFM - QUITE AN ORDEAL IN 2D...
C	---------------------------------------------------------------
C	DECLARATIONS!
	INTEGER ib,jb
	INTEGER	im1,ip1,jm1,jp1
	DOUBLE PRECISION spim,spip,spjm,spjp
	DOUBLE PRECISION spimjm,spimjp,spipjm,spipjp
	DOUBLE PRECISION storeA,minangle
	DOUBLE PRECISION rhoA,rhoB,uA,uB,pA,pB,uAt,uBt,ugt
	DOUBLE PRECISION rhoLI,rhoRI,PI,UI,sp1,sp2
	DOUBLE PRECISION ddx1(9,nx1max,nx2max),ddx2(9,nx1max,nx2max)
	DOUBLE PRECISION advcfl

	advcfl = 1.0*min(dx1,dx2)
C	Put copies of states into both domains for now. We will soon
C	overwrite everything in ghost regions...
	DO i=1,nx1
		DO j=1,nx2
			Uprop(1,1,i,j)=rho(i,j)
			Uprop(1,2,i,j)=rho(i,j)*ux1(i,j)
			Uprop(1,3,i,j)=rho(i,j)*ux2(i,j)
			Uprop(1,4,i,j)=E(i,j)
			Uprop(1,5,i,j)=p(i,j)
			Uprop(1,6,i,j)=ux1(i,j)
			Uprop(1,7,i,j)=ux2(i,j)
			Uprop(2,1,i,j)=rho(i,j)
			Uprop(2,2,i,j)=rho(i,j)*ux1(i,j)
			Uprop(2,3,i,j)=rho(i,j)*ux2(i,j)
			Uprop(2,4,i,j)=E(i,j)
			Uprop(2,5,i,j)=p(i,j)
			Uprop(2,6,i,j)=ux1(i,j)
			Uprop(2,7,i,j)=ux2(i,j)
		END DO
	END DO
	DO k=1,bcount
		i=bci(k)
		j=bcj(k)
C	for each type of boundary cell...
		call FINDPAIR(i,j,ib,jb)
		rhoA = rho(i,j)
		pA = p(i,j)
		uA = ux1(i,j)*normx1(i,j) + 
     +			ux2(i,j)*normx2(i,j)
!		uAt = sqrt(ux1(i,j)**2.0 + 
!     +			ux2(i,j)**2.0 - uA**2.0)
		uAt = ux1(i,j)*normx2(i,j)-
     +			ux2(i,j)*normx1(i,j)
		rhoB = rho(ib,jb)
		pB = p(ib,jb)
		uB = ux1(ib,jb)*normx1(ib,jb) + 
     +			ux2(ib,jb)*normx2(ib,jb)
!		uBt = sqrt(ux1(ib,jb)**2.0 + 
!     +			ux2(ib,jb)**2.0 - uB**2.0)
		uBt = ux1(ib,jb)*normx2(ib,jb)-
     +			ux2(ib,jb)*normx1(ib,jb)
c      !!		for viscosity...
		IF(viscflag.EQ.0)THEN
			ugt = ubt
		ELSE
			ugt = uat
		END IF
c	!!!!
		IF (bcell(i,j) .EQ. 2) THEN
			call RIEMANNINTERFACE (
     +			rhoB,pB,uB,gm(1),pc(1),
     +			rhoA,pA,uA,gm(2),pc(2),
     +			rhoLI,rhoRI,UI,PI,sp1,sp2
     +			)
c			write(6,*) rhoLI,rhoRI,uI,pI
c ! comment the next 2 lines for rGFM
c  ! leave in for mGFM(ish)
c ! copying pressures from air, rGFM for all else
			PI = pB
c
			Uprop(1,1,i,j) = rhoLI
			Uprop(1,5,i,j) = PI
			Uprop(1,8,i,j) = UI
			Uprop(1,9,i,j) = ugt
			Uprop(1,6,i,j) = UI*normx1(i,j)
     +				+ugt*normx2(i,j)
			Uprop(1,7,i,j) = UI*normx2(i,j)
     +				-ugt*normx1(i,j)
			Uprop(1,2,i,j)=rhoLI*Uprop(1,6,i,j)
			Uprop(1,3,i,j)=rhoLI*Uprop(1,7,i,j)
			Uprop(1,4,i,j)=(PI+gm(1)*pc(1))/
     +			(gm(1)-1.0) + 0.5*rhoLI*UI**2.0
		ELSE IF (bcell(i,j) .EQ. 1) THEN
			call RIEMANNINTERFACE (
     +			rhoA,pA,uA,gm(1),pc(1),
     +			rhoB,pB,uB,gm(2),pc(2),
     +			rhoLI,rhoRI,UI,PI,sp1,sp2
     +			)
c ! comment the next line for rGFM
c  ! leave in for mGFM(ish)
c ! Copying pressures from air. rGFM for all else...
			PI = pA
c
			Uprop(2,1,i,j) = rhoRI
			Uprop(2,5,i,j) = PI
			Uprop(2,8,i,j) = UI
			Uprop(2,9,i,j) = ugt
			Uprop(2,6,i,j) = UI*normx1(i,j)
     +				+ugt*normx2(i,j)
			Uprop(2,7,i,j) = UI*normx2(i,j)
     +				-ugt*normx1(i,j)
			Uprop(2,2,i,j)=rhoRI*Uprop(2,6,i,j)
			Uprop(2,3,i,j)=rhoRI*Uprop(2,7,i,j)
			Uprop(2,4,i,j)=(PI+gm(2)*pc(2))/
     +			(gm(2)-1.0) + 0.5*rhoRI*UI**2.0
		END IF
	END DO
C	ADVECT PROPERTIES - use upwind difference scheme...
	DO k=1,int(nx1/5)
		DO i=1,nx1
			im1 = max(i-1,1)
			ip1 = min(i+1,nx1)
			DO j=1,nx2
				jm1 = max(j-1,1)
				jp1 = min(j+1,nx2)
				IF (alpha(i,j).GT.0.AND.
     +				bcell(i,j).EQ.0)THEN
C				If in fluid 2 region and beyond bcells
				DO l=1,9
					IF (normx1(i,j).GT.0.0) THEN
C					Left sided x1 deriv
						ddx1(l,i,j) = 
     +						(Uprop(1,l,i,j)-
     +						Uprop(1,l,im1,j))/dx1
					ELSE
C					right sided x1 deriv
						ddx1(l,i,j) = 
     +						(Uprop(1,l,ip1,j)-
     +						Uprop(1,l,i,j))/dx1
					END IF
					IF (normx2(i,j).GT.0.0) THEN
C					Left sided x2 deriv
						ddx2(l,i,j) = 
     +						(Uprop(1,l,i,j)-
     +						Uprop(1,l,i,jm1))/dx2
					ELSE
C					right sided x2 deriv
						ddx2(l,i,j) = 
     +						(Uprop(1,l,i,jp1)-
     +						Uprop(1,l,i,j))/dx2
					END IF
				END DO
				END IF
			END DO
		END DO
		DO i=1,nx1
		DO j=1,nx2
			IF (alpha(i,j).GT.0.AND.
     +			bcell(i,j).EQ.0)THEN
				DO l=1,9
					Uprop(1,l,i,j)=Uprop(1,l,i,j)- 
     +					advcfl*(normx1(i,j)*
     +					ddx1(l,i,j) + normx2(i,j)*
     +					ddx2(l,i,j)/x1(i))
				END DO
c
				Uprop(1,2,i,j) = Uprop(1,1,i,j)*
     +					Uprop(1,6,i,j)
				Uprop(1,3,i,j) = Uprop(1,1,i,j)*
     +					Uprop(1,7,i,j)
				Uprop(1,4,i,j) = (Uprop(1,5,i,j)+
     +					gm(1)*pc(1))/
     +				(gm(1)-1.0) + 0.5*Uprop(1,1,i,j)*
     +				(Uprop(1,6,i,j)**2.0 + 
     +				Uprop(1,7,i,j)**2.0)
			END IF
		END DO
		END DO
		DO i=1,nx1
			im1 = max(i-1,1)
			ip1 = min(i+1,nx1)
			DO j=1,nx2
				jm1 = max(j-1,1)
				jp1 = min(j+1,nx2)
				IF (alpha(i,j).LE.0.AND.
     +				bcell(i,j).EQ.0)THEN
C				If in fluid 1 region and beyond bcells
				DO l=1,9
					IF (normx1(i,j).LE.0.0) THEN
C					Left sided x1 deriv
						ddx1(l,i,j) = 
     +						(Uprop(2,l,i,j)-
     +						Uprop(2,l,im1,j))/dx1
					ELSE
C					right sided x1 deriv
						ddx1(l,i,j) = 
     +						(Uprop(2,l,ip1,j)-
     +						Uprop(2,l,i,j))/dx1
					END IF
					IF (normx2(i,j).LE.0.0) THEN
C					Left sided x2 deriv
						ddx2(l,i,j) = 
     +						(Uprop(2,l,i,j)-
     +						Uprop(2,l,i,jm1))/dx2
					ELSE
C					right sided x2 deriv
						ddx2(l,i,j) = 
     +						(Uprop(2,l,i,jp1)-
     +						Uprop(2,l,i,j))/dx2
					END IF
				END DO
				END IF
			END DO
		END DO
		DO i=1,nx1
		DO j=1,nx2
			IF (alpha(i,j).LE.0.AND.
     +			bcell(i,j).EQ.0)THEN
				DO l=1,9
					Uprop(2,l,i,j)=Uprop(2,l,i,j)+ 
     +					advcfl*(normx1(i,j)*
     +					ddx1(l,i,j) + normx2(i,j)*
     +					ddx2(l,i,j)/x1(i))
				END DO
c
				Uprop(2,2,i,j) = Uprop(2,1,i,j)*
     +				Uprop(2,6,i,j)
				Uprop(2,3,i,j) = Uprop(2,1,i,j)*
     +				Uprop(2,7,i,j)
				Uprop(2,4,i,j) = (Uprop(2,5,i,j)+
     +				gm(2)*pc(2))/
     +				(gm(2)-1.0) + 0.5*Uprop(2,1,i,j)*
     +				(Uprop(2,6,i,j)**2.0 + 
     +				Uprop(2,7,i,j)**2.0)
			END IF
		END DO
		END DO
	END DO
	DO i=1,nx1
		DO j=1,nx2
			IF (alpha(i,j).LE.0)THEN
				ulsx1(i,j) = Uprop(2,6,i,j)
				ulsx2(i,j) = Uprop(2,7,i,j)
			ELSE
				ulsx1(i,j) = Uprop(1,6,i,j)
				ulsx2(i,j) = Uprop(1,7,i,j)
			END IF
c				ulsx1(i,j) = ux1(i,j)
c				ulsx2(i,j) = ux2(i,j)
		END DO
	END DO
	RETURN
	END
C	---------------------------------------------------------------
	SUBROUTINE RIEMANNINTERFACE (rhoL,pL,uL,GAML,PCRITL,rhoR,pR,uR,
     +					GAMR,PCRITR,
     +					rhoLstar,rhoRstar,ustar,pstar,
     +					bL,bR)
	IMPLICIT NONE
C	---------------------------------------------------------------
C	SUBROUTINE TO SOLVE A RIEMANN PROBLEM - HLLC SOLVER
C	---------------------------------------------------------------
C	SIMPLY RETURNS THE STAR STATES
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		15-07-2013:	CREATED
C		25-07-2013:	EMBEDDED IN SPLITFIELDS.f
	DOUBLE PRECISION	rhoL,PL,uL,GAML,PCRITL,rhoR,PR,uR,
     +				GAMR,PCRITR
C
	DOUBLE PRECISION	eL,eR,HL,HR,GL,GR,PHIL,PHIR
	DOUBLE PRECISION	pstar,ustar,rhoLstar,rhoRstar
	DOUBLE PRECISION	ENLstar,ENRstar
	DOUBLE PRECISION	cL,cR
	DOUBLE PRECISION	bL,bR,bM
	DOUBLE PRECISION	ubar,cbar,rhobar,Povrhobar,Pbar,hbar
	DOUBLE PRECISION	PHIbar,Gbar
	DOUBLE PRECISION	Gdebar,PHIdrhobar
c	CALCULATE LEFT AND RIGHT SPEEDS OF SOUND ======================
	cL = sqrt(GAML*(PL+PCRITL)/rhoL)
	cR = sqrt(GAMR*(PR+PCRITR)/rhoR)
C	CALCULATE LEFT AND RIGHT ENERGY, ENTHALPY, etc ================
	eL = (PL+GAML*PCRITL)/((GAML-1.0)*rhoL)
	eR = (PR+GAMR*PCRITR)/((GAMR-1.0)*rhoR)
	HL = eL + PL/rhoL + 0.5*uL**2.0
	HR = eR + PR/rhoR + 0.5*uR**2.0
	GL = GAML-1.0
	GR = GAMR-1.0
	PHIL = (GAML-1.0)*eL
	PHIR = (GAMR-1.0)*eR
C	CALCULATE ROE AVERAGES ((.)bar ) FOR PROPERTIES ===============
	rhobar = sqrt(rhoL*rhoR)
	call ROEAVG (rhoL,rhoR,uL,uR,ubar)
	call ROEAVG (rhoL,rhoR,HL,HR,hbar)
	Povrhobar = (PL/sqrt(rhoL) + PR/sqrt(rhoR))/
     +			(sqrt(rhoL)+sqrt(rhoR)) + 
     +			0.5*((uR-uL)/(sqrt(rhoL)+sqrt(rhoR)))**2
	Pbar = rhobar*Povrhobar
C	WE ARE LOOKING AT 2 FLUIDS! I KNOW THIS BECAUSE WE ALWAYS ARE -
		call ROEAVG (rhoL,rhoR,GL,GR,Gbar)
		call ROEAVG (rhoL,rhoR,PHIL,PHIR,PHIbar)
		cbar = sqrt(PHIbar + Gbar*Povrhobar)
C	CALCULATE STAR STATES AND WAVE SPEEDS =========================
	bL = min(uL-cL,ubar-cbar)
	bR = max(ubar+cbar,uR+cR)
	ustar = (rhoR*uR*(bR-uR) + rhoL*uL*(uL-bL) + PL - PR)/
     +		(rhoR*(bR-uR) + rhoL*(uL-bL))
	bM = ustar
	pstar = PL + rhoL*(uL-bL)*(uL-ustar)
C 	pstar = PR + rhoR*(bR-uR)*(ustar-uR)
	rhoLstar = rhoL*(bL-uL)/(bL-ustar)
	rhoRstar = rhoR*(bR-uR)/(bR-ustar)
	ENLstar = (pstar + GAML*PCRITL)/(GAML-1.0) + 
     +		0.5*rhoLstar*ustar**2
	ENRstar = (pstar + GAMR*PCRITR)/(GAMR-1.0) + 
     +		0.5*rhoRstar*ustar**2
C	APPROXIMATE SOLVER - LINEAR
c	cbar = 0.5*(cL+cR)
c	rhobar = 0.5*(rhoL+rhoR)
c	pstar = 0.5*(PL+PR) - 0.5*(UR-UL)*rhobar*cbar
c	ustar = 0.5*(uL+uR) - 0.5*(PR-PL)/(rhobar*cbar)
c	rhoLstar = rhoL + (uL-ustar)*rhobar/cbar
c	rhoRstar = rhoR + (ustar-uR)*rhobar/cbar
	RETURN
	END


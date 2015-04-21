	SUBROUTINE EVOLVE

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		20-07-2013:	CREATED
C		07-08-2013:	SPHERICALITY BY SOURCE TERMS
C		07-08-2013:	MODIFIED TO WORK ON CONSERVATIVE VARS
C		13-08-2013:	MUSCL SCHEME IS NOW 2-STEP(ie. CORRECT)
C		14-08-2013:	STOPPED SOLVING DUPLICATE RIEMANN PROBS
C		14-08-2013:	ONLY USE NECESSARY GHOST CELLS..
C		13-09-2013:	MUSCL RECONSTRUCTION OF PRIMITIVE VARS
C		25-09-2013:	WENO-Z RECONSTRUCTION OPTION!
C		25-09-2013:	3rd ORDER TVD TIME INTEGRATION OPTION!
C	---------------------------------------------------------------
C	THE 1 PHASE EULER SOLVER. CHOOSES A FIELD TO WORK ON. RECONSTR-
C	UCTS THE FIELD. CALCULATES FLUXES. UPDATES THE VALUES. ADDS THE
C	SOURCE TERMS. UPDATES THE VALUES OF THE CHOSEN FIELD.
C	---------------------------------------------------------------
	INTEGER x1LOW,x1HIGH,x2LOW,x2HIGH
	DOUBLE PRECISION cosij,sinij,sro,sux1,sux2
c	viscous bits!
	INTEGER im1,ip1,jm1,jp1
	DOUBLE PRECISION divu,mu,taux1x1(nx1max,nx2max),
     +		taux1x2(nx1max,nx2max),
     +		taux2x2(nx1max,nx2max),vis(4,nx1max,nx2max)
	DOUBLE PRECISION uij,uip1,uim1,ujp1,ujm1,vij,vip1,vim1,vjp1,
     +			vjm1
	DOUBLE PRECISION dudr,dudth,dvdr,dvdth,costh,sinth
	DOUBLE PRECISION strp,strrho,strT,sT0,smu0,sC
	x1LOW = 1
	x1HIGH = nx1
	x2LOW = 1
	x2HIGH = nx2
C	FOR EACH FLUID!!
	DO k=1,2
C		HOW SHALL WE SWEEP?
C		APPLY BC IN X1 DIRECTION
		DO j=1,nx2
			Uprop(k,1,1,j) = Uprop(k,1,2,j)
			Uprop(k,2,1,j) = -1.0*Uprop(k,2,2,j)
			Uprop(k,3,1,j) = Uprop(k,3,2,j)
			Uprop(k,4,1,j) = Uprop(k,4,2,j)
		END DO
		IF (boundflag .EQ. 1) THEN
			IF (ghostflag .EQ. 0) THEN
				call NLAABC(k)
			ELSE 
				call NLAABCGHOST(k)
			END IF
		ELSE
			DO j=1,nx2
				Uprop(k,1,nx1,j)=Uprop(k,1,nx1-1,j)
				Uprop(k,2,nx1,j)=-1.*Uprop(k,2,nx1-1,j)
				Uprop(k,3,nx1,j)=Uprop(k,3,nx1-1,j)
				Uprop(k,4,nx1,j)=Uprop(k,4,nx1-1,j)
			END DO
		END IF
C		EVOLVE IN X1 DIRECTION
		call SWEEP1D (k,1)
C		APPLY BC IN X2 DIRECTION
		DO i=1,nx1
			Uprop(k,1,i,1) = Uprop(k,1,i,2)
			Uprop(k,2,i,1) = Uprop(k,2,i,2)
			Uprop(k,3,i,1) = -1.0*Uprop(k,3,i,2)
			Uprop(k,4,i,1) = Uprop(k,4,i,2)
C	
			Uprop(k,1,i,nx2) = Uprop(k,1,i,nx2-1)
			Uprop(k,2,i,nx2) = Uprop(k,2,i,nx2-1)
			Uprop(k,3,i,nx2) = -1.0*Uprop(k,3,i,nx2-1)
			Uprop(k,4,i,nx2) = Uprop(k,4,i,nx2-1)
		END DO
C		EVOLVE IN X2 DIRECTION
		call SWEEP1D (k,2)
C		GRAVITY SOURCE TERMS!!
		IF (grav .NE. 0.0) THEN
		DO i=2,nx1-1
			DO j=2,nx2-1
				cosij = cos(x2(j))
				sinij = sin(x2(j))
				sro = Uprop(k,1,i,j)
				sux1 = Uprop(k,2,i,j)/sro
				sux2 = Uprop(k,3,i,j)/sro
				Uprop(k,2,i,j)=Uprop(k,2,i,j)-dt*sro*
     +					grav*cosij
				Uprop(k,3,i,j)=Uprop(k,3,i,j)+dt*sro*
     +					grav*sinij
				Uprop(k,4,i,j)=Uprop(k,4,i,j)-dt*sro*
     +				grav*(cosij*sux1 - sinij*sux2)
			END DO	
		END DO
		END IF
	END DO
	IF (viscflag .EQ. 1) THEN
c	VISCOUS TERMS!!!
	DO k=1,2
		IF (k.EQ.1)THEN
			mu = 2e-5
			sT0 = 273.15
			sC = 110.4
			smu0 = 1.716e-5
		ELSE
			mu = 1.0D-3
		END IF
		DO i=1,nx1
		DO j=1,nx2
			im1 = max(1,i-1)
			ip1 = min(nx1,i+1)
			jm1 = max(1,j-1)
			jp1 = min(nx2,j+1)
c			set costh and sinth
			costh = cos(x2(j))
			sinth = sin(x2(j))
c			set the u,v for here and next door
			uij = Uprop(k,2,i,j)/Uprop(k,1,i,j)
			uip1 = Uprop(k,2,ip1,j)/Uprop(k,1,ip1,j)
			uim1 = Uprop(k,2,im1,j)/Uprop(k,1,im1,j)
			ujp1 = Uprop(k,2,i,jp1)/Uprop(k,1,i,jp1)
			ujm1 = Uprop(k,2,i,jm1)/Uprop(k,1,i,jm1)
			vij = Uprop(k,3,i,j)/Uprop(k,1,i,j)
			vip1 = Uprop(k,3,ip1,j)/Uprop(k,1,ip1,j)
			vim1 = Uprop(k,3,im1,j)/Uprop(k,1,im1,j)
			vjp1 = Uprop(k,3,i,jp1)/Uprop(k,1,i,jp1)
			vjm1 = Uprop(k,3,i,jm1)/Uprop(k,1,i,jm1)
c			calculate the dudr,dudth,dvdr,dvdth
			IF (i.EQ.1)THEN
				dudr = (uip1-uij)/dx1
				dvdr = (vip1-vij)/dx1
			ELSE IF (i.EQ.nx1) THEN
				dudr = (uij-uim1)/dx1
				dvdr = (vij-vim1)/dx1
			ELSE
				dudr = (uip1-uim1)/(2.0*dx1)
				dvdr = (vip1-vim1)/(2.0*dx1)
			END IF
			IF (j.EQ.1)THEN
				dudth = (ujp1-uij)/dx2
				dvdth = (vjp1-vij)/dx2
			ELSE IF (j.EQ.nx2) THEN
				dudth = (uij-ujm1)/dx2
				dvdth = (vih-vjm1)/dx2
			ELSE
				dudth = (ujp1-ujm1)/(2.0*dx2)
				dvdth = (vjp1-vjm1)/(2.0*dx2)
			END IF
c			calculate divergence of velocity field..
			divu = dudr + 2.0*uij/x1(i) + dvdth/x1(i) + 
     +				vij*costh/(x1(i)*sinth)
c			calculate stress tensor components...
			IF (k.EQ.1.AND.alpha(i,j).LE.0)THEN
				strrho = Uprop(k,1,i,j)
				strp = (gm(1)-1.0)*(Uprop(k,4,i,j) - 
     +					0.5*strrho*(uij**2 +vij**2))
				strT = strp/(287.0*strrho)
				mu = smu0*((sT0+sC)/(strT+sC))*
     +				(strT/sT0)**(3.0/2.0)
			END IF
			taux1x1(i,j)=mu*(2.0*dudr - 2.0*divu/3.0)
			taux1x2(i,j)=mu*(dvdr-(vij/x1(i))+dudth/x1(i))
			taux2x2(i,j)=mu*(2.0*(dvdth/x1(i)+uij/x1(i)) -
     +				2.0*divu/3.0)
		END DO
		END DO
		DO i=2,nx1-1
		DO j=2,nx2-1
c			set costh and sinth
			costh = cos(x2(j))
			sinth = sin(x2(j))
c			calculate viscous terms....
			vis(1,i,j) = 0.0
			vis(2,i,j) = (taux1x1(i+1,j)-taux1x1(i-1,j))
     +				/(2.0*dx1) + 2.0*taux1x1(i,j)/x1(i) 
     +				+ (taux1x2(i,j+1)-taux1x2(i,j-1))/
     +				(2.0*dx2*x1(i)) + 
     +				taux1x2(i,j)*costh/(x1(i)*sinth)		
			vis(3,i,j) = (taux1x2(i+1,j)-taux1x2(i-1,j))/
     +				(2.0*dx1) + 2.0*taux1x2(i,j)/x1(i) 
     +				+ (taux2x2(i,j+1)-taux2x2(i,j-1))/
     +				(2.0*dx2*x1(i)) +
     +				taux2x2(i,j)*costh/(x1(i)*sinth)
			vis(4,i,j) = (uip1*taux1x1(i+1,j)+
     +				vip1*taux1x2(i+1,j) - 
     +				uim1*taux1x1(i-1,j)-
     +				vim1*taux1x2(i-1,j))/(2.0*dx1)
     +				+ 2.0*(uij*taux1x1(i,j)+
     +				vij*taux1x2(i,j))/x1(i) +
     +				(ujp1*taux1x2(i,j+1)+
     +				vjp1*taux2x2(i,j+1)-
     +				ujm1*taux1x2(i,j-1)-
     +				vjm1*taux2x2(i,j-1))/(2.0*dx2*x1(i))
     +				+ (uij*taux1x2(i,j)+vij*taux2x2(i,j))
     +				*costh/(x1(i)*sinth)
		END DO
		END DO
		DO l=1,4		
		DO i=2,nx1-1
		DO j=2,nx2-1
			Uprop(k,l,i,j) = Uprop(k,l,i,j) + dt*vis(l,i,j)
			dUmu(l,i,j)=dt*vis(l,i,j)
		END DO
		END DO
		END DO
	END DO
	END IF
	RETURN
	END

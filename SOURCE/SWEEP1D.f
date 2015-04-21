	SUBROUTINE SWEEP1D (mfl,dflag)

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
	INTEGER dflag,mfl
	INTEGER x1LOW,x1HIGH,x2LOW,x2HIGH
	DOUBLE PRECISION tsfrac,fiph(4,nx1max),fjph(4,nx2max)
	DOUBLE PRECISION sr,sux1,sux2,se,sp
	DOUBLE PRECISION S(4),cosoverRsin
	x1LOW = 1
	x1HIGH = nx1
	x2LOW = 1
	x2HIGH = nx2
	tsfrac = 1.0
	IF (dflag.EQ.1) THEN
C	SWEEPS IN THE x1 DIRECTION
	DO j=x2LOW+1,x2HIGH-1
C		CALCULATE FLUXES
		DO i=x1LOW,x1HIGH-1
C			CALCULATE AT RIGHT "CELL WALLS": i,i+1
			call RIEMANNHLLC (Uprop(mfl,1,i,j),
     +			Uprop(mfl,2,i,j),Uprop(mfl,3,i,j),
     +			Uprop(mfl,4,i,j),Uprop(mfl,1,i+1,j),
     +			Uprop(mfl,2,i+1,j),Uprop(mfl,3,i+1,j),
     +			Uprop(mfl,4,i+1,j),
     +			gm(mfl),Pc(mfl),
     +			fiph(1,i),fiph(2,i),fiph(3,i),fiph(4,i))
		END DO
C		USE THE FLUXES TO UPDATE THE PROPERTIES ---------------
		DO i=x1LOW+1,x1HIGH-1
			sr = Uprop(mfl,1,i,j)
			sux1 = Uprop(mfl,2,i,j)/sr
			sux2 = Uprop(mfl,3,i,j)/sr
			se = Uprop(mfl,4,i,j)
			sp=(gm(mfl)-1)*(se-0.5*sr*(sux1*2 + sux2**2))- 
     +				gm(mfl)*pc(mfl)
			S(1) = 2.0*sr*sux1/x1(i)
			S(2) = 2.0*sr*sux1*sux1/x1(i)
			S(3) = 2.0*sr*sux1*sux2/x1(i)
			S(4) = 2.0*sux1*(se+sp)/x1(i)
			DO k=1,4
				Uprop(mfl,k,i,j) = Uprop(mfl,k,i,j) - 
     +				(tsfrac*dt/dx1)*
     +				(fiph(k,i)-fiph(k,i-1)) - 
     +				tsfrac*dt*S(k)
			END DO
		END DO	
	END DO
	ELSE IF (dflag.EQ.2) THEN
C	SWEEPS IN THE x2 DIRECTION
	DO i=x1LOW+1,x1HIGH-1
C 		CALCULATE FLUXES
		DO j=x2LOW,x2HIGH-1
C			CALCULATE AT RIGHT "CELL WALLS": j,j+1
			call RIEMANNHLLC (Uprop(mfl,1,i,j),
     +			Uprop(mfl,3,i,j),Uprop(mfl,2,i,j),
     +			Uprop(mfl,4,i,j),Uprop(mfl,1,i,j+1),
     +			Uprop(mfl,3,i,j+1),Uprop(mfl,2,i,j+1),
     +			Uprop(mfl,4,i,j+1),
     +			gm(mfl),Pc(mfl),
     +			fjph(1,j),fjph(3,j),fjph(2,j),fjph(4,j))
		END DO
C		USE THE FLUXES TO UPDATE THE PROPERTIES ---------------
		DO j=x2LOW+1,x2HIGH-1
			sr = Uprop(mfl,1,i,j)
			sux1 = Uprop(mfl,2,i,j)/sr
			sux2 = Uprop(mfl,3,i,j)/sr
			se = Uprop(mfl,4,i,j)
			sp=(gm(mfl)-1)*(se-0.5*sr*(sux1*2 + sux2**2))- 
     +				gm(mfl)*pc(mfl)
			cosoverRsin = cos(x2(j))/(x1(i)*sin(x2(j)))
			S(1) = sr*sux2*cosoverRsin
			S(2) = sr*sux1*sux2*cosoverRsin
			S(3) = sr*sux2*sux2*cosoverRsin
			S(4) = sux2*(se+sp)*cosoverRsin
			DO k=1,4
				Uprop(mfl,k,i,j) = Uprop(mfl,k,i,j) - 
     +				(tsfrac*dt/dx2)*
     +				(fjph(k,j)-fjph(k,j-1)) -
     +				tsfrac*dt*S(k)
			END DO
		END DO	
	END DO
	END IF
	RETURN
	END

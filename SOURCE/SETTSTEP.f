	SUBROUTINE SETTSTEP
	
	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		20-11-2012:	CREATED
C	---------------------------------------------------------------
C	SUBROUTINE TO CALCULATE THE TIME STEP BASED ON CFL NUMBER, GRID
C	SIZE, MAXIMUM VELOCITY AND MAXIMUM SOUND SPEED
C	---------------------------------------------------------------
	DOUBLE PRECISION strx1,strx2,strmax
C	FIND THE LARGEST SPEEDS ---------------------------------------
	strmax = 0.0
	DO i=1,nx1
		DO j=1,nx2
			IF (alpha(i,j) .LE. 0 ) THEN
c				Here there is a bodge - c = max(c,300).
				lc(i,j) = max(sqrt(gm(1)*(p(i,j)+Pc(1))
     +					/rho(i,j)),00.0)
			ELSE
				lc(i,j) = sqrt(gm(2)*(p(i,j)+Pc(2))
     +					/rho(i,j))
			END IF
			strx1 = (lc(i,j)+abs(ux1(i,j)))/dx1
			strx2 = (lc(i,j)+abs(ux2(i,j)))/(x1(i)*dx2)
			IF((strx1+strx2).GT.strmax) THEN
				strmax = strx1+strx2
c				write(6,*) "dt limit cell = ",i,j
			END IF
		END DO	
	END DO
C	SET THE TIME STEP ---------------------------------------------
	dt = CFL/strmax
c	REMOVE THIS WHEN FINISHED DEBUGGING w.r.t. 1D code
c	dt = 1e-5
	RETURN
	END

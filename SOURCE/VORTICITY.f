	SUBROUTINE VORTICITY
	
	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		04-06-2014:	CREATED
C	---------------------------------------------------------------
C	SUBROUTINE TO CALCULATE THE VORTICITY FIELD!! 
c	A basic calculation, but we should only bother doing it when we 
c	need to write some outputs...
C	---------------------------------------------------------------
	DOUBLE PRECISION dux1dx2,dux2dx1
C	---------------------------------------------------------------
C	BASIC FIRST ORDER DERIVATIVE CALCULATION...
C	---------------------------------------------------------------
	DO i=2,nx1-1
		DO j=2,nx2-1
			dux1dx2 = (ux1(i,j+1)-ux1(i,j-1))/(2.0*dx2)
			dux2dx1 = (ux2(i+1,j)-ux2(i-1,j))/(2.0*dx1)
			IF (i.EQ.2)THEN
			dux2dx1 = (ux2(i+1,j)-ux2(i,j))/dx1
			END IF
			vort(i,j) = (1.0/x1(i))*(ux2(i,j) + 
     +				x1(i)*dux2dx1 - dux1dx2)
		END DO	
	END DO
	RETURN
	END


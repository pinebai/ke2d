	PROGRAM KINGEULER2D
	
	INCLUDE	"commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
	INTEGER qr
	DOUBLE PRECISION store1
	DOUBLE PRECISION dV
C	OPEN FILES ----------------------------------------------------
C	INPUT FILES = 1
C	FIELD OUTPUT FILES = 7,8,9,10,11,12
C	TIME VARYING PROPERTY FILES = 20,21
c	GRID = 30,31,32,33
C	LOG FILE = 5
	OPEN(unit=1,file="init.params",status='old')
	OPEN(unit=7,file="rho.out")
	OPEN(unit=8,file="P.out")
	OPEN(unit=9,file="ux1.out")
	OPEN(unit=10,file="ux2.out")
	OPEN(unit=11,file="E.out")
	OPEN(unit=12,file="A.out")
	OPEN(unit=13,file="ux.out")
	OPEN(unit=14,file="uz.out")
	OPEN(unit=15,file="vort.out")
	OPEN(unit=16,file="IMPmu.out")
	OPEN(unit=20,file="dt.out")
	OPEN(unit=21,file="interface.out")
	OPEN(unit=22,file="pbound.out")	
	OPEN(unit=23,file="pint.out")
	OPEN(unit=24,file="bmass.out")
	OPEN(unit=25,file="bz.out")
	OPEN(unit=26,file="rhobound.out")
	OPEN(unit=27,file="ubound.out")
	OPEN(unit=30,file="x1.out")
	OPEN(unit=31,file="x2.out")
	OPEN(unit=32,file="x.out")
	OPEN(unit=33,file="z.out")
C	SET SOME INITIAL VARIABLES ------------------------------------
	call SETUP
	call SETTSTEP
	call LEVELSET(0)
	tn=1
C	TIME LOOP -----------------------------------------------------
	DO tn = 1,nt
C		PROGRESS COUNTER
		WRITE (6,*) 'Progress: ', 
     +				100.0*float(tn)/float(nt), '%'
C		SET THE TIME STEP
		call SETTSTEP
	!	IF (t.GT. 0.001)THEN !bodge for seeing ghost wave
	!	dt = 1.48D-6
	!	END IF
C		WRITE RESULTS -----------------------------------------
		call DUMPRESULTS
C		SPLIT FIELDS USING THE GFM
		call GFMSPLIT
C		UPDATE THE LEVEL SET!!!
		call LEVELSET(1)
C		FIND FLUXES AND UPDATE IN TIME
		call EVOLVE
C		RECOMBINE FIELDS
		DO i=1,nx1
			DO j=1,nx2
				IF (alpha(i,j) .LT. 0.0) THEN
C		TO DETERMINE MAGNITUDE OF CHANGE...
			dUprop(1,i,j)=Uprop(1,1,i,j)-rho(i,j)
			dUprop(2,i,j)=Uprop(1,2,i,j)-rho(i,j)*ux1(i,j)
			dUprop(3,i,j)=Uprop(1,3,i,j)-rho(i,j)*ux2(i,j)
			dUprop(4,i,j)=Uprop(1,4,i,j)-E(i,j)
				rho(i,j) = Uprop(1,1,i,j)
				ux1(i,j)=Uprop(1,2,i,j)/Uprop(1,1,i,j)
				ux2(i,j)=Uprop(1,3,i,j)/Uprop(1,1,i,j)
				E(i,j) = Uprop(1,4,i,j)
				P(i,j) =  (gm(1)-1)*(E(i,j)-0.5*
     +				rho(i,j)*(ux1(i,j)**2 + ux2(i,j)**2))
     +				-gm(1)*pc(1)
				ELSE
C		TO DETERMINE MAGNITUDE OF CHANGE...
			dUprop(1,i,j)=Uprop(2,1,i,j)-rho(i,j)
			dUprop(2,i,j)=Uprop(2,2,i,j)-rho(i,j)*ux1(i,j)
			dUprop(3,i,j)=Uprop(2,3,i,j)-rho(i,j)*ux2(i,j)
			dUprop(4,i,j)=Uprop(2,4,i,j)-E(i,j)
				rho(i,j) = Uprop(2,1,i,j)
				ux1(i,j)=Uprop(2,2,i,j)/Uprop(2,1,i,j)
				ux2(i,j)=Uprop(2,3,i,j)/Uprop(2,1,i,j)
				E(i,j) = Uprop(2,4,i,j)
				P(i,j) =  (gm(2)-1)*(E(i,j)-0.5*
     +				rho(i,j)*(ux1(i,j)**2 + ux2(i,j)**2))
     +				-gm(2)*pc(2)
				END IF
C		TO DETERMINE MAGNITUDE OF CHANGE
				DO l=1,4				
					IMPmu(l,i,j)=abs(dUmu(l,i,j))/
     +					abs(dUprop(l,i,j))
				END DO
			END DO
		END DO
c		END DO
C		UPDATE TIME -------------------------------------------
		t = t + dt
	END DO
C	DONE, call FINISH
	call FINISH(0)
	STOP
	END

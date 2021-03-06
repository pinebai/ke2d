	SUBROUTINE SETUP

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		20-11-2012:	CREATED
C		05-08-2013:	A MYRIAD OF CHANGES, UNRECORDED
C	---------------------------------------------------------------
C	SUBROUTINE TO READ THE INPUT FILES FOR PROGRAM KINGEULER1D
C	AND TO SET THE INITIAL FIELDS. THIS SUBROUTINE IS A BIT MESSY.
C	---------------------------------------------------------------
	INTEGER	dump_read_flag
	DOUBLE PRECISION storr2,eta
	READ (1,*) nx1
	READ (1,*) nx2
	READ (1,*) nt
	READ (1,*) outfreq
	READ (1,*) x1d
	READ (1,*) x2d
	READ (1,*) CFL
	READ (1,*) coordsno
	READ (1,*) grav
	READ (1,*) viscflag
	READ (1,*) ghostflag
	READ (1,*) dgun
	READ (1,*) dump_read_flag
	READ (1,*) boundflag
	READ (1,*) distflag
	OPEN(unit=50,file="AG.dat")
	READ (50,*) rho1
	READ (50,*) p1
	READ (50,*) ux11
	READ (50,*) ux21
	READ (50,*) gm(1)
	READ (50,*) pc(1)
	READ (50,*) rho2
	READ (50,*) p2
	READ (50,*) ux12
	READ (50,*) ux22
	READ (50,*) gm(2)
	READ (50,*) pc(2)
	READ (50,*) Rint
	IF (dgun .LE. 1D3) THEN
		p2=p2+(dgun-7.7D0)*rho2*9.81D0 !this is ghost stuff
	ELSE 	! to make this work, for G/NG comparisons, set NG 
		p2=p2+(1D-3*dgun-7.7D0)*rho2*9.81D0 !depth to G depth
	END IF ! multiplied by 1000!
C	SET THE INITIAL FIELDS ----------------------------------------
C	CREATE THE MESH
	call MESHGEN
	DO i=1,nx1
		DO j=1,nx2
C		SET THE FLOW PARAMETERS ACCORDING TO WHETHER WE'RE ----
C		INSIDE OR OUTSIDE Rint --------------------------------
C		eta is a disturbance. Set to zero mostly...
			IF(distflag.EQ.0)THEN
				eta = 0.0
			ELSE IF (distflag.EQ.1)THEN
				eta = 1D-3*sin(20*x2(j))
			ELSE IF (distflag.EQ.2)THEN
				eta = 1D-3*cos(20*x2(j))
			END IF
			IF (x1(i).LE. Rint+eta) THEN
				rho(i,j) = rho1
				P(i,j) = p1 - rho(i,j)*grav*z(i,j)
				ux1(i,j) = ux11
				ux2(i,j) = ux21
				E(i,j)=(P(i,j)+Pc(1)*gm(1))/(gm(1)-1)+
     +				0.5D0*rho(i,j)*(ux1(i,j)**2+
     +				ux2(i,j)**2)
			ELSE
				rho(i,j) = rho2
				p(i,j) = p2 - rho(i,j)*grav*z(i,j)
				ux1(i,j) = ux12
				ux2(i,j) = ux22
				E(i,j)=(P(i,j)+Pc(2)*gm(2))/(gm(2)-1)+
     +				0.5D0*rho(i,j)*(ux1(i,j)**2+
     +				ux2(i,j)**2)
			END IF
			alpha(i,j) = x1(i) - Rint - eta
			ulsx1(i,j) = ux1(i,j)
			ulsx2(i,j) = ux2(i,j)
		END DO
	END DO
	DO j=1,nx2
		Pinf(j) = P(nx1,j)
		Pinf0(j) = P(nx1,j)
	END DO
C	Counter for outputting 
	ofc1 = outfreq + 1
	t=0.0
C	Loading from dump?
c	Shall we start from the last output?! Go on then!!
	IF (dump_read_flag.eq.1)THEN
		OPEN(3,file="DUMP.tron")
		READ(3,*) t
		DO i=1,nx1
			DO j=1,nx2
				READ (3,*) rho(i,j),p(i,j),ux1(i,j),
     +					ux2(i,j),E(i,j),alpha(i,j)
				ulsx1(i,j) = ux1(i,j)
				ulsx2(i,j) = ux2(i,j)
			END DO	
		END DO
	END IF
	RETURN
	END
C	---------------------------------------------------------------
	SUBROUTINE MESHGEN

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	CREATE THE MESH BASED ON x1d,x2d,nx1,nx2. DETERMINE CELL FACE
C	AREAS AND CELL VOLUMES.
C	---------------------------------------------------------------
	DOUBLE PRECISION x1min,x2min,x1ph,x1mh,x2ph,x2mh
	x1min = 0.0D0
	x2min = 0.0D0
!!	dx1 = 1.0D0/49.0D0 !use this to compare domain sizes
	dx1 = (x1d-x1min)/float(nx1-1) !use this normally
	dx2 = (x2d-x2min)/float(nx2-2)
	DO i=1,nx1
		x1(i) = x1min + dx1*float(i-1) + 0.5D0*dx1
	END DO
	DO j=1,nx2
		x2(j) = x2min + dx2*float(j-1) - 0.5D0*dx2
	END DO
C	CALCULATE x,z,r,th VALUES FOR EACH NODE, FOR OUTPUTTING...
	IF (coordsno.EQ.2.0) THEN
		DO i=1,nx1
			DO j=1,nx2
				x(i,j) = x1(i)*sin(x2(j))
				z(i,j) = x1(i)*cos(x2(j))
				r(i,j) = x1(i)
				th(i,j) = x2(j)
			END DO
		END DO
	ELSE IF (coordsno.EQ.1.0)THEN
		DO i=1,nx1
			DO j=1,nx2
				x(i,j) = x1(i)
				z(i,j) = x2(j)
				r(i,j) = sqrt(x1(i)**2 + x2(j)**2)
				th(i,j) = acos(x2(j)/r(i,j))
			END DO
		END DO
	ELSE 
		DO i=1,nx1
			DO j=1,nx2
				x(i,j) = x1(i)
				z(i,j) = x2(j)
				r(i,j) = sqrt(x1(i)**2 + x2(j)**2)
				th(i,j) = acos(x2(j)/r(i,j))
			END DO
		END DO
	END IF
	DO i=1,nx1
		DO j=1,nx2
			WRITE(30,*) x1(i)
			WRITE(31,*) x2(j)
			WRITE(32,*) x(i,j)
			WRITE(33,*) z(i,j)
		END DO
	END DO
	CLOSE(30)
	CLOSE(31)
	CLOSE(32)
	CLOSE(33)
	RETURN
	END

	SUBROUTINE LEVELSET(lsv_flag)

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		25-02-2014:	CREATED
C	---------------------------------------------------------------
C	UPDATE A LEVEL SET!!! -----------------------------------------
	INTEGER lsv_flag
	INTEGER	im3,im2,im1,ip1,ip2,ip3,jm3,jm2,jm1,jp1,jp2,jp3
	DOUBLE PRECISION dadx1(nx1max,nx2max),dadx2(nx1max,nx2max)
	DOUBLE PRECISION dadx1m,dadx1p,dadx2m,dadx2p
	DOUBLE PRECISION Sp(nx1max,nx2max),dls(nx1max,nx2max)
	DOUBLE PRECISION GLS(nx1max,nx2max),ga,gb,gc,gd,dtau
	DOUBLE PRECISION s,storex1,storex2,mgrada
	DOUBLE PRECISION store1,spare1,spR
C	PART 1 - ADVECTION OF LEVEL SET 
	DO j=1,nx2
		DO i=4,nx1-3
			IF (ulsx1(i,j) .GT. 0.0) THEN
C				LEFT SIDED DERIVATIVES ----------------
				call WENOZ((alpha(i-2,j)-alpha(i-3,j)),
     +					(alpha(i-1,j)-alpha(i-2,j)),
     +					(alpha(i,j)-alpha(i-1,j)),
     +					(alpha(i+1,j)-alpha(i,j)),
     +					(alpha(i+2,j)-alpha(i+1,j)),
     +					spare1,store1)
				dadx1(i,j) = store1/dx1
!	dadx1(i,j)=(alpha(i,j)-alpha(i-1,j))/dx1 !1st order!
			ELSE
C				RIGHT SIDED DERIVATIVES ---------------
				call WENOZ((alpha(i+3,j)-alpha(i+2,j)),
     +					(alpha(i+2,j)-alpha(i+1,j)),
     +					(alpha(i+1,j)-alpha(i,j)),
     +					(alpha(i,j)-alpha(i-1,j)),
     +					(alpha(i-1,j)-alpha(i-2,j)),
     +					spare1,store1)
				dadx1(i,j) = store1/dx1
!	dadx1(i,j)=(alpha(i+1,j)-alpha(i,j))/dx1 !1st order!
			END IF
		END DO
		DO i=2,3
			IF (ulsx1(i,j) .GT. 0.0) THEN
C				LEFT SIDED DERIVATIVES ----------------
				dadx1(i,j) = (alpha(i,j)-alpha(i-1,j))
     +					/dx1
			ELSE
C				RIGHT SIDED DERIVATIVES ---------------
				dadx1(i,j) = (alpha(i+1,j)-alpha(i,j))
     +					/dx1
			END IF
		END DO
		DO i=nx1-2,nx1-1
			IF (ulsx1(i,j) .GT. 0.0) THEN
C				LEFT SIDED DERIVATIVES ----------------
				dadx1(i,j) = (alpha(i,j)-alpha(i-1,j))
     +					/dx1
			ELSE
C				RIGHT SIDED DERIVATIVES ---------------
				dadx1(i,j) = (alpha(i+1,j)-alpha(i,j))
     +					/dx1
			END IF
		END DO
		dadx1(1,j) = dadx1(2,j)
		dadx1(nx1,j) = dadx1(nx1-1,j)
	END DO
	DO i=1,nx1
		DO j=4,nx2-3
			IF (ulsx2(i,j) .GT. 0.0) THEN
C				LEFT SIDED DERIVATIVES ----------------
				call WENOZ((alpha(i,j-2)-alpha(i,j-3)),
     +					(alpha(i,j-1)-alpha(i,j-2)),
     +					(alpha(i,j)-alpha(i,j-1)),
     +					(alpha(i,j+1)-alpha(i,j)),
     +					(alpha(i,j+2)-alpha(i,j+1)),
     +					spare1,store1)
				dadx2(i,j) = store1/dx2
!	dadx2(i,j) = (alpha(i,j)-alpha(i,j-1))/dx2 !1st order!
			ELSE
C				RIGHT SIDED DERIVATIVES ---------------
				call WENOZ((alpha(i,j+3)-alpha(i,j+2)),
     +					(alpha(i,j+2)-alpha(i,j+1)),
     +					(alpha(i,j+1)-alpha(i,j)),
     +					(alpha(i,j)-alpha(i,j-1)),
     +					(alpha(i,j-1)-alpha(i,j-2)),
     +					spare1,store1)
				dadx2(i,j) = store1/dx2
!	dadx2(i,j) = (alpha(i,j-1)-alpha(i,j))/dx2 !1st order!
			END IF
		END DO
		DO j=2,3
			IF (ulsx2(i,j) .GT. 0.0) THEN
C				LEFT SIDED DERIVATIVES ----------------
				dadx2(i,j) = (alpha(i,j)-alpha(i,j-1))
     +					/dx2
			ELSE
C				RIGHT SIDED DERIVATIVES ---------------
				dadx2(i,j) = (alpha(i,j+1)-alpha(i,j))
     +					/dx2
			END IF
		END DO
		DO j=nx2-2,nx2-1
			IF (ulsx2(i,j) .GT. 0.0) THEN
C				LEFT SIDED DERIVATIVES ----------------
				dadx2(i,j) = (alpha(i,j)-alpha(i,j-1))
     +					/dx2
			ELSE
C				RIGHT SIDED DERIVATIVES ---------------
				dadx2(i,j) = (alpha(i,j+1)-alpha(i,j))
     +					/dx2
			END IF
		END DO
C		CONSTRAINT: theta-derivatives of level set equal zero
c		at theta-boundaries. Physically, this means the bubble
C		doesn't have any corners - which is as it should be!
		dadx2(i,1) = 0.0
		dadx2(i,nx2) = 0.0
	END DO
	IF (lsv_flag.EQ.1)THEN
	DO i=1,nx1
		DO j=1,nx2
			alpha(i,j) = alpha(i,j) - dt*
     +			(ulsx1(i,j)*dadx1(i,j)+ulsx2(i,j)*dadx2(i,j)
     +			/x1(i))
		END DO
	END DO
	j=1
	DO i=1,nx1
		alpha(i,j) = alpha(i,j+1)
	END DO
	i=1
	DO j=1,nx2
		alpha(i,j) = alpha(i+1,j)
	END DO
	j=nx2
	DO i=1,nx1
		alpha(i,j) = alpha(i,j-1)
	END DO
	END IF
C	PART 2: RE-INITIALISATION OF THE LEVEL SET! 
C	Based on the scheme by Russo & Smereka (2000), JCP 163, p51-67
C	Or maybe on Spelt (2005), JCP 207, p389-404.
C	Label all boundary cells!
	call SETBCELL
	DO i=1,nx1
		DO j=1,nx2
C			1. Calculate the sign of the level set at 
c			each node. "Mollified" sign function.
			sp(i,j) = alpha(i,j)/sqrt(alpha(i,j)**2.0 + 
     +				dx1**2.0)
c			2. Set dLS(i,j) to the value of the level set.
			dLS(i,j) = alpha(i,j)
		END DO
	END DO
	dtau = 0.1*x1(1)*dx1*dx2/(x1(1)*dx2+dx1)
	DO k=1,int(2*nx1/5)
C		3. Calculate GLS, upwind disc' of abs(grad(alpha))-1
		DO i=1,nx1
		im1 = max(1,i-1)
		ip1 = min(nx1,i+1)
		DO j=1,nx2
			jm1 = max(1,j-1)
			jp1 = min(nx2,j+1)
			ga = (alpha(i,j)-alpha(im1,j))/dx1
			gb = (alpha(ip1,j)-alpha(i,j))/dx1
			gc = (alpha(i,j)-alpha(i,jm1))/(dx2*x1(i))
			gd = (alpha(i,jp1)-alpha(i,j))/(dx2*x1(i))
			IF (alpha(i,j).GT.0)THEN
				GLS(i,j) = sqrt(max((max(ga,0.0)**2),
     +				(min(gb,0.0)**2))+max((max(gc,0.0)**2),
     +				(min(gd,0.0)**2)))-1.0
			ELSE
				GLS(i,j) = sqrt(max((min(ga,0.0)**2),
     +				(max(gb,0.0)**2))+max((min(gc,0.0)**2),
     +				(max(gd,0.0)**2)))-1.0
			END IF
		END DO
		END DO
C		4. Update the values of dLS
		DO i=1,nx1
		DO j=1,nx2
			IF (bcell(i,j) .EQ. 0)THEN
C			Not a bcell, normal reinitialisation!
				dLS(i,j) = dLS(i,j) - dtau*sp(i,j)*
     +					GLS(i,j)
			ELSE
C			A bcell, special reinitialisation??
c		or just leave it be!!!
			END IF
		END DO
		END DO
	END DO
	IF (lsv_flag.EQ.1)THEN
	DO i=1,nx1
	DO j=1,nx2
		alpha(i,j) = dLS(i,j)
	END DO
	END DO
	j=1
	DO i=1,nx1
		alpha(i,j) = alpha(i,j+1)
	END DO
	i=1
	DO j=1,nx2
		alpha(i,j) = alpha(i+1,j)
	END DO
	j=nx2
	DO i=1,nx1
		alpha(i,j) = alpha(i,j-1)
	END DO
	END IF	
C	PART 3: CALCULATION OF NORMALS 
	call SETBCELL
	DO i=2,nx1-1
		DO j=2,nx2-1
			dadx1(i,j) = (alpha(i+1,j)-alpha(i-1,j))/
     +				(2.0*dx1)
			dadx2(i,j) = (alpha(i,j+1)-alpha(i,j-1))/
     +				(2.0*dx2)
			spare1 = sqrt(dadx1(i,j)**2.0 + 
     +				(dadx2(i,j)/x1(i))**2.0)
			normx1(i,j) = dadx1(i,j)/spare1
			normx2(i,j) = dadx2(i,j)/spare1
		END DO
	END DO
	DO i=2,nx1-1
		normx1(i,1) = normx1(i,2)
		normx1(i,nx2) = normx1(i,nx2-1)
		normx2(i,1) = normx2(i,2)
		normx2(i,nx2) = normx2(i,nx2-1)

	END DO
	DO j=2,nx2-1
		normx1(1,j) = normx1(2,j)
		normx1(nx1,j) = normx1(nx1-1,j)
		normx2(1,j) = normx2(2,j)
		normx2(nx1,j) = normx2(nx1-1,j)
	END DO
	normx1(1,1) = normx1(2,2)
	normx2(1,1) = normx2(2,2)
	normx1(nx1,1) = normx1(nx1-1,2)
	normx2(nx1,1) = normx2(nx1-1,2)
	normx1(1,nx2) = normx1(2,nx2-1)
	normx2(1,nx2) = normx2(2,nx2-1)
	normx1(nx1,nx2) = normx1(nx1-1,nx2-1)
	normx2(nx1,nx2) = normx2(nx1-1,nx2-1)
C	Calculate the bubble volume and hence an effective radius
	spR = 0.0
	DO j=1,nx2
		DO i=1,nx1-1
			IF (alpha(i,j)*alpha(i+1,j).LE.0.0)THEN
				q(j)=i
			END IF
		END DO
		Rint = x1(q(j)) - alpha(q(j),j)/
     +			max(1.0,-1.0*alpha(q(j),j)/dx1)
		IF (alpha(q(j),j).Eq.0.0)THEN
			Rint = x1(q(j))
		ELSE IF (alpha(q(j)+1,j).EQ.0.0)THEN
			Rint = x1(q(j)+1)
		END IF
		spR = spR + (2.0/3.0)*pie*(Rint**3.0)*sin(x2(j))*dx2
	END DO
	Rint = (0.75*spR/pie)**(1.0/3.0)
	RETURN
	END 

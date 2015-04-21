	SUBROUTINE SETBCELL

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		04-03-2014:	CREATED
C	---------------------------------------------------------------
C	SET BCELL VALUES
C	---------------------------------------------------------------
C	DECLARATIONS!
	INTEGER	im1,ip1,jm1,jp1
	DOUBLE PRECISION spim,spip,spjm,spjp
	DOUBLE PRECISION spimjm,spimjp,spipjm,spipjp
	DOUBLE PRECISION storeA
C	Set a flag to say whether a cell is beside the interface...
	bcount = 0
	DO i=2,nx1-1
		DO j=1,nx2
			im1=max(1,i-1)
			ip1=min(nx1,i+1)
			jm1=max(1,j-1)
			jp1=min(nx2,j+1)
			spim = alpha(i,j)*alpha(im1,j)
			spip = alpha(i,j)*alpha(ip1,j)
			spjm = alpha(i,j)*alpha(i,jm1)
			spjp = alpha(i,j)*alpha(i,jp1)
			spimjm = alpha(i,j)*alpha(im1,jm1)
			spimjp = alpha(i,j)*alpha(im1,jp1)
			spipjm = alpha(i,j)*alpha(ip1,jm1)
			spipjp = alpha(i,j)*alpha(ip1,jp1)
			storeA = min(spim,spip,spjm,spjp,spimjm,spimjp,
     +				spipjm,spipjp)
			IF (storeA.GT.0)THEN
				bcell(i,j) = 0
			ELSE
				bcount = bcount + 1
				bci(bcount) = i
				bcj(bcount) = j
				IF (alpha(i,j).LE.0) THEN
					bcell(i,j) = 1
				ELSE
					bcell(i,j) = 2
				END IF
			END IF
		END DO
	END DO
	RETURN
	END

	SUBROUTINE FINDPAIR (iA,jA,iB,jB)

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		04-03-2014:	CREATED
C	---------------------------------------------------------------
C	SPLIT THE FIELDS USING A GFM
C	---------------------------------------------------------------
C	DECLARATIONS!
	INTEGER ia,ja,ib,jb
	INTEGER	im1,ip1,jm1,jp1
	INTEGER tc(2,8)
	DOUBLE PRECISION spim,spip,spjm,spjp
	DOUBLE PRECISION spimjm,spimjp,spipjm,spipjp
	DOUBLE PRECISION nax,naz,nbx,nbz
	DOUBLE PRECISION storeA,minangle
	im1 = max(1,ia-1)
	ip1 = min(nx1,ia+1)
	jm1 = max(1,ja-1)
	jp1 = min(nx2,ja+1)
	tc(1,1) = ip1
	tc(2,1) = ja
	tc(1,2) = ip1
	tc(2,2) = jm1
	tc(1,3) = ia
	tc(2,3) = jm1
	tc(1,4) = im1
	tc(2,4) = jm1
	tc(1,5) = im1
	tc(2,5) = ja
	tc(1,6) = im1
	tc(2,6) = jp1
	tc(1,7) = ia
	tc(2,7) = jp1
	tc(1,8) = ip1
	tc(2,8) = jp1
	minangle = 1e12
C	Loop over all neighbours
	DO k=1,8
C		Check whether the potential neighbour cell is a bound-
c		-ary cell of the right type..
		IF (bcell(tc(1,k),tc(2,k)).NE.0.AND.
     +		bcell(tc(1,k),tc(2,k)).NE.bcell(ia,ja))THEN
C		It's a potential pair. Is it the best choice?
		nax = normx1(ia,ja)*sin(x2(ja)) + 
     +			normx2(ia,ja)*cos(x2(ja))
		naz = normx1(ia,ja)*cos(x2(ja)) - 
     +			normx2(ia,ja)*sin(x2(ja))
		nbx = normx1(tc(1,k),tc(2,k))*sin(x2(tc(2,k))) + 
     +			normx2(tc(1,k),tc(2,k))*cos(x2(tc(2,k)))
		nbz = normx1(tc(1,k),tc(2,k))*cos(x2(tc(2,k))) - 
     +			normx2(tc(1,k),tc(2,k))*sin(x2(tc(2,k)))
		storeA = nax*nbx + naz*nbz
c		storeA = normx1(ia,ja)*normx1(tc(1,k),tc(2,k)) + 
c     +			normx2(ia,ja)*normx2(tc(1,k),tc(2,k))
			IF (abs(1.0-storeA).LE.abs(1.0-minangle))THEN
C			If it's the closest yet then it'll do for now!
				minangle = storeA
				ib = tc(1,k)
				jb = tc(2,k)
			END IF
		END IF
	END DO
C	write(6,*) ia,ja,ib,jb,bcell(ia,ja)
	if (ib.lt.1.or.ib.gt.nx1.or.jb.lt.1.or.jb.gt.nx2)then
		write(6,*) "FIND PAIR BAD MATCH!!!",ia,ja,ib,jb
		call FINISH(2)
	end if
	RETURN
	END

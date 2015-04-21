	SUBROUTINE NLAABCGHOST (mfl)

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		20-07-2013:	CREATED
C		06-08-2013:	MODIFIED TO UPSTREAM WENO-Z SCHEME
C		30-09-2013:	FOR A GENERAL NODE 'ibn' AT RIGHT BOUND
C		??<06-2014:	NOW 2D...		
c		25-06-2014:	INCLUSION OF THE GHOST!!
C	---------------------------------------------------------------
C	THIS SUBROUTINE USES THE STANDARD CHARACTERISTIC BOUNDARY COND-
C	ITION FORMULATION PRESENTED IN Thompson, JCP 89, 439-461(1990)
C	TO APPLY BOUNDARY CONDITIONS AT THE RIGHT EDGE OF THE DOMAIN
C
C	TAKES A CELL INDEX ibn AND USES SOLUTION AT THIS CELL TO APPLY
C	NLAA BC BY PROVIDING NEW VALUES FOR THIS CELL -----------------
C	---------------------------------------------------------------
C	WARNING: The inclusion of the ghost will not work correctly if
C	KE2D is run from a dump file. For simulations including the 
C	ghost, KE2D must be run from scratch.
C	---------------------------------------------------------------
	INTEGER mfl
	DOUBLE PRECISION drhodrB,durdrB,dutdrB,dpdrB
	DOUBLE PRECISION drhodthB,durdthB,dutdthB,dpdthB
	DOUBLE PRECISION cB,L1B,L2B,L3B,L4B
	DOUBLE PRECISION henth
	DOUBLE PRECISION aL1,NA,NB
	DOUBLE PRECISION durdtB,dutdtB,drhodtB,dpdtB
	DOUBLE PRECISION newrho(nx2max),newp(nx2max),
     +			newur(nx2max),newut(nx2max)
c	o=o=o=o=o=o=o=o=o=o=o=o=o
	INTEGER ig
	DOUBLE PRECISION rg,thg,d,taug,ttarget,ref_ghost
	DOUBLE PRECISION ug,pg,dugdrg,dpgdrg,rhog
	DOUBLE PRECISION ugr,ugth,dugrdr,dugrdth,dpgdr,dpgdth,dugthdr
     +		,dugthdth
	DOUBLE PRECISION cos1,cos2,sin1,sin2,gL1,gL2,gL3
	DOUBLE PRECISION dubdt,dhdt
	IF (mfl .EQ. 2)THEN !if looking at water...
c	0. Set f,df,d2f...
	IF (tn .EQ. 1) THEN
		ubold = ux1(nx1,2)
	END IF
	cB = sqrt(gm(2)*(p(nx1,2)+pc(2))/rho(nx1,2))
	henth = (p(nx1,2)-pinf(2))/rho(nx1,2)
	ug = ux1(nx1,2)
	dubdt = (ug-ubold)/dt
	f(tn) = x1(nx1)*x1(nx1)*(ug-henth/cB - 
     +			(ug**2)/(2.0*cB))
	df(tn) = x1(nx1)*(henth + 0.5*ug**2)
	d2f(tn) = -1.0*cb*henth -0.5*cb*ug**2 + x1(nx1)*cB*dubdt
	thist(tn) = t
	ubold = ux1(nx1,2)
	d = dgun
c	o=o=o=o=o=o=o=o=o=o=o=o=o
C	FOR EVERY NODE ON THE nx1 BOUNDARY
	DO j=2,nx2-1
		rob = Uprop(mfl,1,nx1,j)
		robm1 = Uprop(mfl,1,nx1-1,j)
		robm2 = Uprop(mfl,1,nx1-2,j)
		robm3 = Uprop(mfl,1,nx1-3,j)
		robm = Uprop(mfl,1,nx1,j-1)
		robp = Uprop(mfl,1,nx1,j+1)
		ub = Uprop(mfl,2,nx1,j)/Uprop(mfl,1,nx1,j)
		ubm1 = Uprop(mfl,2,nx1-1,j)/Uprop(mfl,1,nx1-1,j)
		ubm2 = Uprop(mfl,2,nx1-2,j)/Uprop(mfl,1,nx1-2,j)
		ubm3 = Uprop(mfl,2,nx1-3,j)/Uprop(mfl,1,nx1-3,j)
		ubm = Uprop(mfl,2,nx1,j-1)/Uprop(mfl,1,nx1,j-1)
		ubp = Uprop(mfl,2,nx1,j+1)/Uprop(mfl,1,nx1,j+1)
		vb = Uprop(mfl,3,nx1,j)/Uprop(mfl,1,nx1,j)
		vbm1 = Uprop(mfl,3,nx1-1,j)/Uprop(mfl,1,nx1-1,j)
		vbm2 = Uprop(mfl,3,nx1-2,j)/Uprop(mfl,1,nx1-2,j)
		vbm3 = Uprop(mfl,3,nx1-3,j)/Uprop(mfl,1,nx1-3,j)
		vbm = Uprop(mfl,3,nx1,j-1)/Uprop(mfl,1,nx1,j-1)
		vbp = Uprop(mfl,3,nx1,j+1)/Uprop(mfl,1,nx1,j+1)
		pb = (gm(mfl)-1.0)*(Uprop(mfl,4,nx1,j)-
     +			0.5*rob*(ub**2+vb**2))-gm(mfl)*pc(mfl)
		pbm1 = (gm(mfl)-1.0)*(Uprop(mfl,4,nx1-1,j)-
     +			0.5*robm1*(ubm1**2+vbm1**2))-gm(mfl)*pc(mfl)
		pbm2 = (gm(mfl)-1.0)*(Uprop(mfl,4,nx1-2,j)-
     +			0.5*robm2*(ubm2**2+vbm2**2))-gm(mfl)*pc(mfl)
		pbm3 = (gm(mfl)-1.0)*(Uprop(mfl,4,nx1-3,j)-
     +			0.5*robm3*(ubm3**2+vbm3**2))-gm(mfl)*pc(mfl)
		pbm = (gm(mfl)-1.0)*(Uprop(mfl,4,nx1,j-1)-
     +			0.5*robm*(ubm**2+vbm**2))-gm(mfl)*pc(mfl)
		pbp = (gm(mfl)-1.0)*(Uprop(mfl,4,nx1,j+1)-
     +			0.5*robp*(ubp**2+vbp**2))-gm(mfl)*pc(mfl)
C		CALCULATE THE LOCAL SPEED OF SOUND --------------------
		cB = sqrt(gm(2)*(pb+pc(2))/rob)
C		FIRST CALCULATE DERIVATIVES OF PRIMITIVE VARIABLES AT -
C		BOUNDARY: 1 SIDED 3rd ORDER ---------------------------
		drhodrB = (-11.0*rob+18.0*robm1-9.0*robm2+2.0*robm3)/
     +			(-6.0*dx1)
		durdrB = (-11.0*ub+18.0*ubm1-9.0*ubm2+2.0*ubm3)/
     +			(-6.0*dx1)
		dutdrB = (-11.0*vb+18.0*vbm1-9.0*vbm2+2.0*vbm3)/
     +			(-6.0*dx1)
		dpdrB = (-11.0*pb+18.0*pbm1-9.0*pbm2+2.0*pbm3)/
     +			(-6.0*dx1)
C		DERIVATIVES IN THE THETA DIRECTION!!
		drhodthB= (robp-robm)/(2.0*dx2)
		durdthB = (ubp-ubm)/(2.0*dx2)
		dutdthB = (vbp-vbm)/(2.0*dx2)
		dpdthB = (pbp-pbm)/(2.0*dx2)
C		ASSUME L4 OUTGOING ------------------------------------
		L4B = (ub+cB)*(dpdrB + rob*cB*durdrB)
c		=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=
c		GHOST!!
c		=<o>=<o>=<o>=<o>=<o>=<o>=<o>=<o>=<o>=<o>=<o>=<o>=<o>=<o
c		1. Find rg,thg
		thg = atan(sin(x2(j))/((2.0*d/x1(nx1))-cos(x2(j))))
		rg = x1(nx1)*sin(x2(j))/sin(thg)
c		2. Caclulate delay...
		taug = (rg-x1(nx1))/cB
c		search through to find the index of f we're looking for
		ttarget = t-taug
		IF (ttarget .LE. 0)THEN
			ig = 1
		ELSE
			DO i=1,tn
				IF (thist(i).LE.ttarget.AND.
     +				thist(i+1).GT.ttarget)THEN
					ig=i
				END IF
			END DO
		END IF
c		3. set ug,pg,dugdrg,dpgdrg,rhog
		rhog=rho2 !CHEAT!
		IF (ghostflag .EQ. 1)THEN
			ref_ghost = 1.0
		ELSE IF (ghostflag .EQ. -1)THEN
			ref_ghost = -1.0
		ELSE 
			ref_ghost = 0.0
		END IF
		ug = ref_ghost*(f(ig)/(rg*rg) + df(ig)/(rg*cb))
		pg = ref_ghost*(df(ig)/rg - 0.5*ug*ug)*rhog + pinf(j)
		dugdrg = ref_ghost*(-2.0*f(ig)/(rg**3) - 
     +			2.0*df(ig)/(rg*rg*cb) - d2f(ig)/(rg*cb*cb))
		dpgdrg = -1.0*ref_ghost*rhog*(ug*dugdrg +df(ig)/(rg*rg)
     +			 + d2f(ig)/(rg*cb))
c		4. Find ugr,ugth,dugrdr,dugrdth,dpgdr,dpgdth,
c		dugthdr,dugthdth
		cos1 = cos(pie-x2(j)-thg)
		sin1 = sin(pie-x2(j)-thg)
		ugr = ug*cos1
		ugth = ug*sin1
		dugrdr = dugdrg*cos1*cos1
		dugrdth = x1(nx1)*dugdrg*sin1*cos1 + ug*sin1
		dugthdr = dugdrg*cos1*sin1
		dugthdth = x1(nx1)*dugdrg*sin1*sin1 - ug*cos1
		dpgdr = dpgdrg*cos1
		dpgdth = x1(nx1)*dpgdrg*sin1
  !!!!!!	Modulate the hydrostatic pressure!!
		pinf(j)=pinf0(j)+ref_ghost*(df(ig)/rg-0.5*ug*ug)*rhog		
c		5. Use the above to set gL1,gL2,gL3
		IF (ttarget .LE. 0.0D0) THEN
			gL1 = 0.0D0
			gL2 = 0.0D0
			gL3 = 0.0D0
		ELSE
			gL1 = (ux1(nx1,j)-cB)*(dpgdr-rhog*cb*dugrdr)
			gL2 = -1.0*ux1(nx1,j)*dpgdr
			gL3 = ux1(nx1,j)*dugthdr
		END IF
C		NLAA BCs!!!
C		+++++++++++++++++++++++++++++++++++++++++++++++++++++++
		henth = (pb-Pinf(j))/rob
C
		aL1 = (rob*(ub-cB)/x1(nx1))*((3.0-ub/cB)*0.5*ub**2.0 
     +		      + 2.0*ub*cB - (1.0+ub/cB)*henth)
     +			+ gL1 !!!!!! =o=o=o=o=
		NA = 1.0 - 0.5*(ub-cB)*ub/(cB**2.0)
		NB = (1.0/(2.0*rob*cB))*(L4B-aL1)
     +			+ vb*durdthB/x1(nx1)+ (vb**2.0)/x1(nx1)
     +			+grav*cos(x2(j)) !
		durdtB = -1.0*NB/NA
		L1B = aL1+rob*ub*(ub-cB)*durdtB/cB
		IF (ub .LT. 0.0) THEN
C			INCOMING
			L2B = ((rob*ub**2.0)/x1(nx1))*(henth/cB - 
     +				2.0*ub + 0.5*(ub**2.0)/cB) + rob*ub*
     +				(1.0-ub/cB)*durdtB
     +				+ gL2 !!! =o=o=o=o=o=
			L3B = 0.0 + gL3 !!=o=o=o=o=o=o
		ELSE
C			OUTGOING --------------------------------------
			L2B = ub*(cB*cB*drhodrB-dpdrB)
			L3B = ub*dutdrB
		END IF
C		DETERMINE THE TIME DERIVATIVES ------------------------
C			NEED TO CHECK SOURCE TERMS!!!
		drhodtB = -1.0*(L2B + 0.5*(L4B+L1B))/cB**2.0
     +			- vb*drhodthB/x1(nx1)
     +			- rob*dutdthB/x1(nx1)
     +			- 2.0*rob*ub/x1(nx1)
     +			- rob*vb/(x1(nx1)*tan(x2(j)))

		dpdtB = -0.5*(L4B+L1B)
     +			-vb*dpdthB/x1(nx1)
     +			-rob*cB*cB*dutdthB/x1(nx1)
     +			- 2.0*rob*ub*cB*cB/x1(nx1)
     +			- rob*vb*cB*cB/(x1(nx1)*tan(x2(j)))

		dutdtB = -1.0*L3B
     + 			- vb*dutdthB/x1(nx1)
     +			- dpdthB/(rob*x1(nx1))
     +			- ub*vb/x1(nx1)
     +			+ grav*sin(x2(j)) 
C		1st ORDER TIME INTEGRATION ----------------------------
		newrho(j) = rob + drhodtB*dt
		newp(j) = pb + dpdtB*dt
		newur(j) = ub + durdtB*dt
		newut(j) = vb + dutdtB*dt
	END DO
	DO j=2,nx2-1
		Uprop(mfl,1,nx1,j) = newrho(j)
		Uprop(mfl,2,nx1,j) = newrho(j)*newur(j)
		Uprop(mfl,3,nx1,j) = newrho(j)*newut(j)
		Uprop(mfl,4,nx1,j) = (newp(j)+gm(2)*pc(2))/(gm(2)-1.0)+ 
     +			0.5*newrho(j)*(newur(j)**2.0 + newut(j)**2.0)
	END DO
	ELSE
		DO j=2,nx2-1
C			STICK IN A REFLECTING BOUNDARY...
			Uprop(mfl,1,nx1,j) = Uprop(mfl,1,nx1-1,j)
			Uprop(mfl,2,nx1,j) = -1.0*Uprop(mfl,2,nx1-1,j)
			Uprop(mfl,3,nx1,j) = Uprop(mfl,3,nx1-1,j)
			Uprop(mfl,4,nx1,j) = Uprop(mfl,4,nx1-1,j)
		END DO
	END IF
	RETURN
	END 

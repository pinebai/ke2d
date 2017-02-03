	SUBROUTINE ROEAVG (roL,roR,ftnL,ftnR,RAVG)

	IMPLICIT NONE
	DOUBLE PRECISION roL,roR,ftnL,ftnR,RAVG
	RAVG = (sqrt(roL)*ftnL + sqrt(roR)*ftnR)/
     +			(sqrt(roL)+sqrt(roR))
	RETURN
	END 
C	---------------------------------------------------------------
 	subroutine riemannhllc(roL,mL,mtL,EL,roR,mR,mtR,ER,
     +	    gam,pcrit,fl1,fl2,fl3,fl4)
	implicit none
   	real roL,mL,mtL,EL
    	real roR,mR,mtR,ER
    	real gam
    	real fl1,fl2,fl3,fl4
    	real uL,PL,uR,PR,utL,utR,ubar,utbar
    	real usqL,usqR,usqbar
    	real hL,hR,hbar
    	real cL,cR,cbar
    	real bL,bR,SL,SR,SM,SX
    	real wL,wR,qL,qR
    	real roX,uX,utX,EX,PX
    	real roS,mS,mtS,ES,PS
    	real scratch

    	logical subsonic
    	! calculate velocities and averages
    	uL = mL/roL
    	uR = mR/roR
    	utL = mtL/roL
    	utR = mtR/roR
    	call roeavg(roL,roR,uL,uR,ubar)
    	call roeavg(roL,roR,utL,utR,utbar)
    	usqL = uL**2.0 + utL**2.0
    	usqR = uR**2.0 + utR**2.0
    	usqbar = ubar**2.0 + utbar**2.0
    	! calculate pressures, enthalpies and sound speeds
    	PL = (gam-1.0)*(EL-0.5*roL*usqL) - gam*pcrit
    	PR = (gam-1.0)*(ER-0.5*roR*usqR) - gam*pcrit
    	hL=(EL+PL)/roL
    	hR=(ER+PR)/roR
    	cL = sqrt(gam*(PL+pcrit)/roL)
    	cR = sqrt(gam*(PR+pcrit)/roR)
    	call roeavg(roL,roR,hL,hR,hbar)
C    	cbar=sqrt((gam-1.0)*(hbar-0.5*usqbar))
C    	bL=sqrt((gam-1.0)/(2.0*gam))
C    	bR=sqrt((gam-1.0)/(2.0*gam))
C    	SL=min(ubar-cbar,uL-bL*cL)
C    	SR=max(ubar+cbar,uR+bR*cR)
	cbar = sqrt((gam-1.0)*(hbar-0.5*(ubar**2+utbar**2)))
	SL = min(uL-cL,ubar-cbar)
	SR = max(ubar+cbar,uR+cR)
    
    
    	! directly calculate fluxes
    	subsonic = .true.
    	if (SL .gt. 0.0) then
    	   fl1 = roL*uL
    	   fl2 = fl1*uL + PL
    	   fl3 = fl1*utL	
    	   fl4 = uL*(PL + EL)
    	   subsonic = .false.
    	end if
    	if (SR .lt. 0.0) then
    	   fl1 = roR*uR
    	   fl2 = fl1*uR + PR
    	   fl3 = fl1*utR
    	   fl4 = uR*(PR + ER)
    	   subsonic = .false.
    	end if
    	if (subsonic) then
    	   qL=roL*uL*(SL-uL)-PL
    	   qR=roR*uR*(SR-uR)-PR
    	   wL=roL*(SL-uL)
    	   wR=roR*(SR-uR)
    	   SM=(qR-qL)/(wR-WL)
    	   if (SM .gt. 0.0) then
    	      roX=roL
    	      uX=uL
    	      utX=utL
    	      EX=EL
    	      PX=PL
    	      SX=SL
    	   else
    	      roX=roR
    	      uX=uR
    	      utX=utR
    	      EX=ER
    	      PX=PR
    	      SX=SR
    	   end if
    	   scratch=SX-SM
    	   roS=roX*(SX-uX)/scratch
    	   pS=roX*(uX-SX)*(uX-SM)+pX
    	   mS=((SX-uX)*roX*uX + (pS-PX))/scratch
    	   mtS=((SX-uX)*roX*utX)/scratch
    	   ES= ((SX-uX)*EX - pX*uX + pS*SM)/scratch
    	   if(SM .gt. 0.0) then
    	      fl1 = roL*uL + SX*(roS-roX)
    	      fl2 = roL*uL*uL + PL + SX*(mS-roX*uX)
    	      fl3 = roL*uL*utL + SX*(mtS-roX*utX)
    	      fl4 = uL*(PL+EL) + SX*(ES-EX)
    	   else
    	      fl1 = roR*uR + SX*(roS-roX)
    	      fl2 = roR*uR*uR + PL + SX*(mS-roX*uX)
    	      fl3 = roR*uR*utR + SX*(mtS-roX*utX)
    	      fl4 = uR*(PR+ER) + SX*(ES-EX)
    	   end if
    	end if
    	return
 	end
C	---------------------------------------------------------------

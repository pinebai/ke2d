	SUBROUTINE WENOZ (fim2,fim1,fi,fip1,fip2,fimh,fiph)

	IMPLICIT NONE
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		??-03-2013:	CREATED
C	---------------------------------------------------------------
C	SUBROUTINE TO CALCULATE FLUXES USING BORGES 5TH ORDER MODIFIED
C	WENO SCHEME (as mentioned in Hu et al 2009)
C	---------------------------------------------------------------
C	INPUT FUNCTION
	DOUBLE PRECISION	fim2,fim1,fi,fip1,fip2
C	OUTPUT WENOZ APPROXIMATION
	DOUBLE PRECISION	fimh,fiph
C	SMOOTHNESS INDICATORS
	DOUBLE PRECISION	beta0,beta1,beta2
	DOUBLE PRECISION	tau5,eps
	DOUBLE PRECISION	beta0z,beta1z,beta2z
C	WEIGHTINGS
	DOUBLE PRECISION	d0,d1,d2
	DOUBLE PRECISION	alpha0z,alpha1z,alpha2z,sumalphaz
	DOUBLE PRECISION	w0z,w1z,w2z
C	SMALL STENCIL INTERPOLATIONS
	DOUBLE PRECISION	fimh0,fimh1,fimh2,fiph0,fiph1,fiph2
C	SET BASIC WEIGHTINGS ------------------------------------------
	d0=0.3
	d1=0.6
	d2=0.1
C	CALCULATE INITIAL SMOOTHNESS INDICATORS -----------------------
	beta0 = (13.0/12.0)*(fim2-2.0*fim1+fi)**2.0 + 
     +		0.25*(fim2-4.0*fim1+3.0*fi)**2.0
	beta1 = (13.0/12.0)*(fim1-2.0*fi+fip1)**2.0 + 
     +		0.25*(fim1-fip1)**2.0
	beta2 = (13.0/12.0)*(fi-2.0*fip1+fip2)**2.0 + 
     +		0.25*(3.0*fi-4.0*fip1+fip2)**2.0
	tau5 = abs(beta0-beta2)
C	CALCULATE BETTER SMOOTHNESS INDICATORS ------------------------
	eps = 1e-40
	beta0z = (beta0+eps)/(beta0+tau5+eps)
	beta1z = (beta1+eps)/(beta1+tau5+eps)
	beta2z = (beta2+eps)/(beta2+tau5+eps)
C	CALCULATE DIFFERENT WEIGHTINGS --------------------------------	
	alpha0z = d0/beta0z
	alpha1z = d1/beta1z
	alpha2z = d2/beta2z
	sumalphaz = alpha0z + alpha1z + alpha2z
C	CALCULATE FINAL WEIGHTINGS ------------------------------------
	w0z = alpha0z/sumalphaz
	w1z = alpha1z/sumalphaz
	w2z = alpha2z/sumalphaz
C	CALCULATE THREE 3 POINT STENCIL ESTIMATES AT i+0.5 ------------
	fiph0 = fim2/3.0 - 7.0*fim1/6.0 + 11.0*fi/6.0
	fiph1 = -1.0*fim1/6.0 + 5.0*fi/6.0 + fip1/3.0
	fiph2 = fi/3.0 + 5.0*fip1/6.0 - fip2/6.0
C	AND A 5th ORDER ESTIMATE USING THE FINAL WEIGHTS --------------
	fiph = w0z*fiph0 + w1z*fiph1 + w2z*fiph2	
C	CALCULATE THREE 3 POINT STENCIL ESTIMATES AT i-0.5 ------------
	fimh0 = -1.0*fim2/6.0 + 5.0*fim1/6.0 + fi/3.0
	fimh1 =	fim1/3.0 + 5.0*fi/6.0 - fip1/6.0
	fimh2 = 11.0*fi/6.0 - 7.0*fip1/6.0 + fip2/3.0
C	AND A 5th ORDER ESTIMATE USING THE FINAL WEIGHTS --------------
	fimh = w0z*fimh0 + w1z*fimh1 + w2z*fimh2
	RETURN
	END

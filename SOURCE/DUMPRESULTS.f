	SUBROUTINE DUMPRESULTS

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		20-11-2012:	CREATED
C		??<08-2013:	MODIFIED FOR FAR FIELD OUTPUT IF REQ
C		??<04-2014:	MODIFIED FOR THE 2D CODE...
C	---------------------------------------------------------------
C	SUBROUTINE TO WRITE THE OUTPUTS IF IT IS THE REQUIRED TIME
C	---------------------------------------------------------------
	INTEGER ofspc
	CHARACTER dmpfile*(14),tron*(5)
	DOUBLE PRECISION sp1,sp2,sp3,sp4
	DOUBLE PRECISION v_cell,m_cell,summz,summ,b_mass,b_z,rplus
C	EITHER WRITE TO OUTPUT OR ADD 1 TO THE OUTPUT COUNTER ---------
	IF (ofc1 .GE. outfreq) THEN
C		IT IS TIME TO WRITE THE OUTPUT
		WRITE (6,*) "WRITING"
c		what is the vorticity?
		call VORTICITY
C		WRITE (6,*) "-------- t = ",t
		write(tron,'(i5.5)') tn
		tron=adjustl(tron)
		dmpfile="DUMP"//tron//".DUMP"
		OPEN(2,file=dmpfile)
		write(2,*) t
		DO i=1,nx1
			DO j=1,nx2
				WRITE(7,*) rho(i,j)
				WRITE(8,*) p(i,j)
				WRITE(9,*) ux1(i,j)
				WRITE(10,*) ux2(i,j)
				WRITE(11,*) E(i,j)
				WRITE(12,*) alpha(i,j)
				WRITE(13,*) ux1(i,j)*sin(x2(j))+
     +					ux2(i,j)*cos(x2(j))
				WRITE(14,*) ux1(i,j)*cos(x2(j))-
     +					ux2(i,j)*sin(x2(j))
				WRITE(15,*) vort(i,j)
				WRITE(16,*) IMPmu(1,i,j),IMPmu(2,i,j),
     +					IMPmu(3,i,j),IMPmu(4,i,j)
				WRITE(2,*) rho(i,j),p(i,j),ux1(i,j),
     +					ux2(i,j),E(i,j),alpha(i,j)
			END DO
		END DO
		close(2)
C		RESET THE OUTPUT COUNTER
		ofc1 = 1
	ELSE
C		ADD 1 TO THE OUTPUT COUNTER
		ofc1 = ofc1 + 1
	END IF
C	Writing other stuff
	WRITE (20,*) t,dt !time step
	WRITE (21,*) t,Rint	! interface location
	sp1 = 0.0
	sp2 = 0.0
	sp3 = 0.0
	sp4 = 0.0
	DO j=2,nx2-1
		sp1 = sp1+p(nx1,j)
		sp2 = sp2+p(q(j),j)*sin(x2(j))*x1(q(j))**2.0
		sp3 = sp3 + sin(x2(j))*x1(q(j))**2.0
		sp4 = sp4+ux1(q(j),j)*sin(x2(j))*x1(q(j))**2.0
	END DO
	sp1 = sp1/float(nx2-2)
!	sp2 = sp2/float(nx2-2)
	sp2 = sp2/sp3	
	WRITE(22,*) t,p(nx1,25) ! boundary pressure
	WRITE(23,*) t,sp2 ! interface pressure
	WRITE(26,*) t,rho(nx1,25)
	WRITE(27,*) t,ux1(nx1,25)
	summz=0.0
	summ=0.0
	DO j=1,nx2
		DO i=1,nx1-1
			IF (alpha(i,j).LT.0)THEN
			IF (alpha(i,j)*alpha(i+1,j).LE.0.0)THEN
				q(j)=i
			END IF
			Rint = x1(q(j)) - alpha(q(j),j)/
     +				max(1.0,-1.0*alpha(q(j),j)/dx1)
			IF (alpha(q(j),j).Eq.0.0)THEN
				Rint = x1(q(j))
			ELSE IF (alpha(q(j)+1,j).EQ.0.0)THEN
				Rint = x1(q(j)+1)
			END IF
			IF (alpha(i,j)*alpha(i+1,j).LE.0.0)THEN
				rplus = Rint
			ELSE
				rplus = x1(i)+0.5*dx1
			END IF
				v_cell = (1.0/3.0)*((rplus)**
     +				3.0 - (x1(i)-0.5*dx1)**3.0)*
     +				(cos(x2(j)-0.5*dx2)-cos(x2(j)+0.5*dx2))
     +				*2.0*pie
				m_cell = rho(i,j)*v_cell
!				m_cell = v_cell
				summz = summz + m_cell*z(i,j)
				summ = summ + m_cell
			END IF
		END DO
	END DO
	b_z = summz/summ
	b_mass = summ
	WRITE(24,*) t,b_mass
	WRITE(25,*) t,b_z
	RETURN
	END

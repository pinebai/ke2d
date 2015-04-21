	SUBROUTINE FINISH(error_number)

	IMPLICIT NONE
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		19-06-2014:	CREATED
C	---------------------------------------------------------------
C	SUBROUTINE TO STOP THE PROGRAM...
c	---------------------------------------------------------------
	INTEGER error_number,i
	WRITE (6,*) "KE2D TERMINATING"
	IF (error_number .EQ. 0) THEN
		WRITE (6,*) "NORMALLY"
	END IF	
	IF (error_number .EQ. 1) THEN
		WRITE (6,*) "ERROR CODE = ",error_number
		WRITE (6,*) "BAD INPUT FILE REQUESTED"
	END IF	
	IF (error_number .EQ. 2) THEN
		WRITE (6,*) "ERROR CODE = ",error_number
		WRITE (6,*) "FAILURE TO CORRECTLY FIND PAIRS"
		WRITE (6,*) "CLOSING TO PREVENT SEGMENTATION FAULT"
	END IF	
	WRITE (6,*) "CLOSING FILES"
	CLOSE(1)
	DO i=7,16
		CLOSE(i)
	END DO
	DO i=20,26
	CLOSE(i)
	END DO
	DO i=30,33
		CLOSE(i)
	END DO
	CLOSE (41)
	WRITE (6,*) "FINISHED. GOOD BYE!"
	STOP
	END

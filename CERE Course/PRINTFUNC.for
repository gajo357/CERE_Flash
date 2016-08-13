      MODULE PRINTFUNC
      
      PUBLIC PRINT_RESULT, PRINT_STABILITY_RESULT
      CONTAINS
      
      
      SUBROUTINE PRINT_RESULT(HEADER, P, T, NCA, Z, X, Y, BETA, 
     &                        COMPS, NITER, IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      CHARACTER(*) HEADER
      INTEGER NCA, NITER
      INTEGER IER
      DOUBLE PRECISION T,P
      DOUBLE PRECISION Z(NCA), X(NCA), Y(NCA)
      CHARACTER*4 COMPS(NCA)
      DOUBLE PRECISION BETA
      
      WRITE(*,*)
      WRITE(*,"(A50)") HEADER
      WRITE(*,3000) "P = ", P, " MPa"
      WRITE(*,3000) "T = ", T, " K"
3000  FORMAT (A, F10.3, A)
      WRITE(*,*)
      WRITE(*,"(A, F10.4)") "CALCULATED BETA ", BETA 
      WRITE(*,*)
      WRITE(*,"(A, I5)") "NUMBER OF ITERATIONS ", NITER
      
      IF (IER .GT. 0) THEN          
          WRITE(*,*)
          WRITE(*,"(A, I5)") "AN ERROR OCCURED ", IER
          GOTO 101
      END IF      
      
      WRITE(*,*)
      WRITE(*,2000) "COMPS", "TOTAL", "VAPOR", "LIQUID"
2000  FORMAT (4(A15))
      
      DO I = 1,NCA
          WRITE(*,1000) COMPS(I), Z(I), Y(I), X(I)
      ENDDO
1000  FORMAT (A15, F15.3, F15.3, F15.3) 
      
101   CONTINUE
      END SUBROUTINE
      
      SUBROUTINE PRINT_STABILITY_RESULT(HEADER, P, T, TRIAL, IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      CHARACTER(*) HEADER
      INTEGER IER
      DOUBLE PRECISION T, P, TRIAL(2)
      
      WRITE(*,*)
      WRITE(*,"(A50)") HEADER
      WRITE(*,"(A, F10.3, A)") "P = ", P, " MPa"
      WRITE(*,"(A, F10.3, A)") "T = ", T, " K"
      WRITE(*,*)
      IF (ABS (TRIAL(1)) .LT. 1E-5) THEN
          WRITE(*,*) "LIQUID TRIAL TRIVIAL"
      ELSE
          WRITE(*,"(A, F10.4, A)") "LIQUID TRIAL", TRIAL(1)
      END IF
      
      IF (ABS (TRIAL(2)) .LT. 1E-5) THEN
          WRITE(*,*) "VAPOR TRIAL TRIVIAL"
      ELSE
          WRITE(*,"(A, F10.4, A)") "VAPOR TRIAL", TRIAL(2)
      END IF
      WRITE(*,*)
      
      END SUBROUTINE
      
      END MODULE PRINTFUNC
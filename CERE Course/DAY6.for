      MODULE DAY6
      
      USE THERMCALC
      USE DAY2
      
      PUBLIC RUN_DAY5
      CONTAINS
      
      SUBROUTINE RUN_DAY6(NCA, T, P, Z, COMPS)
C     NCA      no. of components      (i)
C     Z       feed mole fractions     (i)
C     P, T    PRESSURE, TEMPERATURE   (I)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NCA
      DOUBLE PRECISION T,P
      DOUBLE PRECISION Z(NCA)
      CHARACTER*4 COMPS(NCA)
      
      INTEGER IER
      DOUBLE PRECISION X(NCA), Y(NCA)
      DOUBLE PRECISION TDEW, TBUB, PDEW, PBUB, BETA
            
      BETA = 0
C      CALL GETSATTEMP(NCA, P, Z, BETA, TBUB, X, Y, IER)
C      CALL PRINTRSULT("BUBBLE POINT TEMPERATURE", P, TBUB, IER)
C      CALL GETSATPRES(NCA, T, Z, BETA, PBUB, IER)
C      CALL PRINTRSULT("BUBBLE POINT PRESSURE", PBUB, T, IER)
      
      BETA = 1
C      CALL GETSATTEMP(NCA, P, Z, BETA, TDEW, X, Y, IER)
C      CALL PRINTRSULT("DEW POINT TEMPERATURE", P, TDEW, IER)
C      CALL GETSATPRES(NCA, T, Z, BETA, PDEW, IER)
C      CALL PRINTRSULT("DEW POINT PRESSURE", PDEW, T, IER)
      
      CALL PHAENV(NCA, Z, IER)
      
      END SUBROUTINE

      SUBROUTINE PHAENV(NCA, Z, IER)
      IMPLICIT DOUBLE PRECISION (A-H, O,Z)
C     VAR   VARIABLE TO GO BY
      PARAMETER (NPOINTS=1000)
      INTEGER NCA, VAR, IER
      DOUBLE PRECISION Z(NCA)
      DOUBLE PRECISION PHAT(NPOINTS), PHAP(NPOINTS)
      
      INTEGER IPOINT, INDX(NCA + 2), MAXITER, NITER, TNITER
      DOUBLE PRECISION P, T, BETA, SPEC
      DOUBLE PRECISION ZNORM(NCA), L(NCA), V(NCA)
      DOUBLE PRECISION X(NCA + 2), X_1(NCA + 2)
      DOUBLE PRECISION DXDS(NCA + 2), DXDS_1(NCA + 2), TEMP, DS
      CHARACTER CHOICE
      LOGICAL RETRO
                 
      RETRO = .FALSE.
      IER = 0
      
      SPEC = 0

C     NORMALIZE THE FEED COMPOSITION
      DO I=1,NCA
          ZNORM(I) = Z(I) / SUM(Z)
      ENDDO  
      
C     SMALL INITIAL PRESSURE AND START FROM THE BUBBLE POINT LINE
      P = 0.5
      BETA = 0
C     NUMBER OF ITERATIONS FOR SAT T ALGORITHM
      MAXITER = 3
      CALL GETSATTEMP(NCA, P, ZNORM, BETA, T, L, V, 
     &                IER, MAXITER = MAXITER)
C     SET X FACTOR
      DO I = 1,NCA
          X(I) = LOG(V(I)/L(I))
      ENDDO
      X(NCA + 1) = LOG(T)
      X(NCA + 2) = LOG(P)
      X_1(:) = 0
C     SELECT PRESSURE AS FIRST SPECIFICATION
      VAR = NCA + 2
      DS = 1
      
C     MAX NUMBER OF ITERATIONS FOR NEWTON ALGORITHM
      MAXITER = 50
      TNITER = 4
      IPOINT = 1
      DO WHILE (.TRUE.)
          
          IF (IPOINT .EQ. 1)  GOTO 100
          
C         CALCULATE SENSITIVITIES AND PROCEED CHOOSE THE NEW STEP
300       CONTINUE
          IF (RETRO) GOTO 400
C         GET THE VARIABLE WITH BIGGEST SENSITIVITY
          TEMP = 0
          DO I = 1,NCA + 2              
              IF(ABS(DXDS(I)) .GT. TEMP) THEN
                  VAR = I
                  TEMP = ABS(DXDS(I))
              ENDIF
          ENDDO
          
          IF(NITER .GE. TNITER) THEN
              DS = DS/(1.2*(NITER - TNITER + 1))
          ELSE 
              DS = DS *(1.2*(TNITER - NITER))
          ENDIF
          
C         CHECK IF STEP IS TOO SMALL AND BRAKE
400       IF (DS .LT. 1E-5) GOTO 500
          DO I = 1,NCA + 2
              X(I) = X(I) + DXDS(I)*DS
          ENDDO

C         RUN THE CALCULATION TO FIND TH POINT
100       CALL RUNNEWTON(NCA, ZNORM, BETA, VAR, MAXITER,
     &                    X, DXDS, NITER, IER)
          IF(IER .GT. 0) THEN
              WRITE(*,*) "THE CALCULATION DID NOT CONVERGE ", DS
              X(:) = X_1(:)
              DXDS(:) = DXDS_1(:)
              DS = DS/2
              GOTO 400
          ENDIF
          
C         WE ARE IN THE RETROGRADE AREA, P,T SHOUD DECREASE
          IF(RETRO) THEN
              IF (X(NCA + 1) .GE. X_1(NCA + 1) .OR. 
     &            X(NCA + 2) .GE. X_1(NCA + 2)) THEN 
                  X(:) = X_1(:)
                  DXDS(:) = DXDS_1(:)                     
                  DS = DS/2                  
                  GOTO 400
              ENDIF
          ENDIF          
     
C         IF BOTH P AND T ARE DECREASING, WE ARE IN THE RETROGRADE AREA
          IF (X(NCA + 1) .LT. X_1(NCA + 1) .AND. 
     &        X(NCA + 2) .LT. X_1(NCA + 2)) THEN
              IF (RETRO) THEN
                  GOTO 200
              ELSE                               
                  WRITE(*,*) "ENTERED RETROGRADE CALCULATION", VAR
                  DS = -DS
                  VAR = NCA + 2
                  RETRO = .TRUE.
                  GOTO 200
              ENDIF              
          ENDIF
         
C         POINT IS ACCEPTED
200       CONTINUE
          PHAT(IPOINT) = EXP(X(NCA + 1))
          PHAP(IPOINT) = EXP(X(NCA + 2))
          IPOINT = IPOINT + 1
          X_1(:) = X(:)
          DXDS_1(:) = DXDS(:)
          
1000      FORMAT(A10, F30.20)
          WRITE(*,*) "STEP SOLUTION"
C          DO I=1,NCA
C              WRITE(*,"(A, I2, F30.20)") "lnK", I, X(I)
C          ENDDO
C          WRITE(*,1000) "lnT", X(NCA + 1)
C          WRITE(*,1000) "lnP", X(NCA + 2)
          WRITE(*,*) "VARIABLE ", VAR, " NUMBER ITERATIONS ", NITER
          WRITE(*,1000) "T", PHAT(IPOINT - 1)
          WRITE(*,1000) "P", PHAP(IPOINT - 1)
          WRITE(*,*)
          
C         IF P OR T ARE BELOW MINIMUM, WE ARE DONE
          IF (PHAT(IPOINT - 1) .LT. 0.5 .OR. 
     &        PHAP(IPOINT - 1) .LT. 0.5) THEN
              WRITE(*,*) "TERMINATED"
              GOTO 500
          ENDIF
      ENDDO

500   OPEN(13,FILE='PHAENV.TXT',STATUS='NEW',BLANK='ZERO',ERR=101)
      GOTO 301
      OPEN(13,FILE='PHAENV.TXT',STATUS='OLD',BLANK='ZERO',ERR=101)
301   CONTINUE
      DO I = 1,NPOINTS
          IF (PHAT(I) .LE. 0 .OR. PHAP(I) .LE. 0) GOTO 600
          WRITE(13, "(F30.20, F30.20)") PHAT(I), PHAP(I)
      ENDDO      
600   CLOSE(13)
      
101   CONTINUE
      END SUBROUTINE
      
      SUBROUTINE RUNNEWTON(NCA, Z, BETA, VAR, MAXITER,
     &                    X, DXDS, NITER, IER)
      IMPLICIT DOUBLE PRECISION (A-H, O,Z)
C     VAR   VARIABLE TO GO BY
      INTEGER NCA, VAR, IER, MAXITER
      DOUBLE PRECISION Z(NCA)
      
      INTEGER IPOINT, INDX(NCA + 2)
      DOUBLE PRECISION BETA, SPEC
      DOUBLE PRECISION KFACT(NCA)
      DOUBLE PRECISION L(NCA), V(NCA)
      DOUBLE PRECISION FUG(NCA), FUGL(NCA), FUGV(NCA)
      DOUBLE PRECISION FUGT(NCA), FUGTL(NCA), FUGTV(NCA)
      DOUBLE PRECISION FUGP(NCA), FUGPL(NCA), FUGPV(NCA)
      DOUBLE PRECISION FUGX(NCA, NCA), FUGXL(NCA, NCA), FUGXV(NCA, NCA)
      DOUBLE PRECISION JACOB(NCA + 2, NCA + 2), F(NCA + 2)
      DOUBLE PRECISION X(NCA + 2), DX(NCA + 2)
      DOUBLE PRECISION DXDS(NCA + 2)
      DOUBLE PRECISION MAXERROR
      
      MAXERROR = 1E-6
      
      IER = 0
      
      NITER = 1                        
      DO WHILE (.TRUE.)
C         GET NEW VAPOR LIQUID COMPOSITION
          DO I = 1,NCA
              KFACT(I) = EXP(X(I))
          ENDDO
          CALL GETCOMPOSITION(NCA, Z, BETA, KFACT, L, V)     
          
C         CHECK FUGACITY CEFFICIENTS
          CALL THERMO(EXP(X(NCA + 1)), EXP(X(NCA + 2)), L, FUGL,
     &     FUGT = FUGTL, FUGP = FUGPL, FUGX=FUGXL)
C         CHECK FUGACITY CEFFICIENTS
          CALL THERMO(EXP(X(NCA + 1)), EXP(X(NCA + 2)), V, FUGV,
     &     FUGT = FUGTV, FUGP = FUGPV, FUGX=FUGXV)
          
C         CALCULATE JACOBIAN
          DO I = 1,NCA+2
C             SET dF/dS
              DXDS(I) = 0
                            
C             CALCULATE FUNCTION F
              IF(I .LE. NCA) THEN
                  F(I) = X(I) + FUGV(I) - FUGL(I)
              ELSE IF(I .EQ. NCA + 1) THEN
                  F(I) = 0
                  DO J = 1,NCA
                      F(I) = F(I) + (V(J) - L(J))
                  ENDDO
              ELSE 
                  F(I) = 0
                  SPEC = X(VAR)
              ENDIF
              
              DO J = 1,NCA+2
                  JACOB(I,J) = 0
                  
C                 DF(I)/DlnK(J)   I,J = 1...NCA
                  IF(I .LE. NCA .AND. J .LE. NCA) THEN
                      TEMP = 0
                      IF (I .EQ. J) TEMP = 1.0
                      JACOB(I,J) = TEMP + 
     &                    FUGXV(I, J)*(1 - BETA)*L(J)*V(J)/Z(J) +
     &                    FUGXL(I, J)*BETA*L(J)*V(J)/Z(J)
C                 Df(I)/DlnP      I=1...NCA
                  ELSE IF (I .LE. NCA .AND. J .EQ. NCA + 2) THEN
                      JACOB(I,J) = (FUGPV(I) - FUGPL(I)) * EXP(X(J)) 
C                 Df(I)/DlnT      I=1...NCA
                  ELSE IF (I .LE. NCA .AND. J .EQ. NCA + 1) THEN
                      JACOB(I,J) = (FUGTV(I) - FUGTL(I)) * EXP(X(J)) 
C                 Df(NCA+1)/DlnK(J)
                  ELSE IF (I .EQ. NCA + 1 .AND. J .LE. NCA) THEN
                      JACOB(I,J) = L(J)*V(J)/Z(J)
C                 Df(NCA+1)/DlnP    Df(NCA + 1)/DlnT          
                  ELSE IF (I .EQ. NCA + 1 .AND. J .GT. NCA) THEN
                      JACOB(I,J) = 0
C                 SPEC DERIVATIVE
                  ELSE IF (I .EQ. NCA + 2) THEN
                      IF(VAR .EQ. J) JACOB(I,J) = 1
                  ENDIF                  
              ENDDO
          ENDDO

C         LAST VALUE OF dF/dS EQUALS 1
          DXDS(NCA + 2) = 1
          
C         INVERT JACOBIAN AND GET DELTA X
          CALL LUDEC(NCA + 2, INDX, JACOB)
          CALL BACKSUBST(NCA + 2, INDX, JACOB, F)
          CALL BACKSUBST(NCA + 2, INDX, JACOB, DXDS)
          
C         DELTA X COMES IN THE VALUES OF F
C         CALCULATE NEW X
          DO I = 1,NCA + 2
              DX(I) = F(I)
              IF (DX(I) .NE. DX(I)) THEN
                  IER = 601
                  GOTO 101
              END IF
              
              X(I) = X(I) - DX(I)              
          ENDDO
          IF(EXP(X(NCA + 1)) .LT. 0 .OR.
     &       EXP(X(NCA + 2)) .LE. 0) THEN
              IER = 601
              GOTO 101
          ENDIF          
          
C         CHECK IF CONVERGED
          DO I = 1,NCA + 2
              IF (DX(I) .GE. MAXERROR) GOTO 100
          ENDDO          
          GOTO 200
          
100       CONTINUE
          IF (NITER .GE. MAXITER) THEN
              IER = 701
              GOTO 101
          ENDIF
          NITER = NITER + 1
      ENDDO
      
200   CONTINUE
      
101   CONTINUE
      END SUBROUTINE

      SUBROUTINE GETSATPRES(NCA, T, Z, BETA, PSAT, IER)
C     NCA      no. of components      (i)
C     Z       feed mole fractions     (i)
C     P                    PRESSURE   (I)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NCA
      DOUBLE PRECISION PSAT, T, BETA
      DOUBLE PRECISION Z(NCA), ZNORM(NCA)
      
      INTEGER IER, ITER, CHECKCOMP
      DOUBLE PRECISION KFACT(NCA), X(NCA), Y(NCA)
      DOUBLE PRECISION FUGX(NCA), FUGY(NCA)
      DOUBLE PRECISION FUGPX(NCA), FUGPY(NCA)
      DOUBLE PRECISION F, DF
      DOUBLE PRECISION ERROR, MAXERROR
      DOUBLE PRECISION DP, DY(NCA), DX(NCA)
      DOUBLE PRECISION EP(1000), EK(1000)
      
      IER = 0
      MAXERROR = 1E-5
      
      CHECKCOMP = 1
      
      DO I=1,NCA
          ZNORM(I) = Z(I) / SUM(Z)
      ENDDO      
      
C     GET INITIAL K FACTOR ESTIMATES
      CALL GETINITP(NCA, T, ZNORM, BETA, PSAT)
      CALL WILSON_KFACT(NCA, PSAT, T, KFACT)
      CALL GETCOMPOSITION(NCA, ZNORM, BETA, KFACT, X, Y)
      
      ITER = 1
      DO WHILE(.TRUE.)
C         LOG THE VALUES
          EP(ITER) = PSAT
          EK(ITER) = KFACT(CHECKCOMP)
          
          
C         CHECK FUGACITY CEFFICIENTS
          CALL THERMO(T, PSAT, X, FUGX,
     &     FUGP = FUGPX)
C         CHECK FUGACITY CEFFICIENTS
          CALL THERMO(T, PSAT, Y, FUGY,
     &     FUGP = FUGPY)
          
C         CALCULAT THE K FACTORS
          DO I = 1,NCA
              KFACT(I) = EXP(FUGX(I) - FUGY(I))
          ENDDO
          
C         CALCULTE NEW COMPOSITION
          CALL GETCOMPOSITION(NCA, ZNORM, BETA, KFACT, DX, DY)
          DO I = 1,NCA
              TEMP = DY(I)
              DY(I) = TEMP - Y(I)              
              Y(I) = TEMP
              
              TEMP = DX(I)
              DX(I) = TEMP - X(I)
              X(I) = TEMP
          ENDDO
          
C         CALCULATE F AND DF/DT          
          DF = 0
          F = -1
          DO I = 1,NCA
C              F = F + (Y(I) - X(I))
              
              IF(BETA .GT. 1E-5) THEN
                  F = F + ZNORM(I) / KFACT(I)
                  DF = DF + ZNORM(I) / KFACT(I) * (FUGPY(I) - FUGPX(I))
              ELSE
                  F = F + ZNORM(I) * KFACT(I)
                  DF = DF + ZNORM(I) * KFACT(I) * (FUGPX(I) - FUGPY(I))
              ENDIF
          ENDDO

          DP = -F/DF
C         CHECK FOR ILLEGAL SPOTS
          IF(DP .NE. DP) THEN
              IER = 501
              GOTO 101
          ENDIF          
          PSAT = PSAT + DP
          
C         CHECK CONVERGENCE
          IF(ABS(DP) .LT. MAXERROR) THEN
              DO I=1,NCA
                  IF(ABS(DY(I)) .GE. MAXERROR) GOTO 100
                  IF(ABS(DX(I)) .GE. MAXERROR) GOTO 100
              ENDDO
              
              GOTO 200
          ENDIF
          
100       CONTINUE
          
          ITER = ITER + 1
      ENDDO

200   CONTINUE
      
      
      OPEN(13,FILE='PSATERROR.TXT',STATUS='NEW',BLANK='ZERO',ERR=101)
      GOTO 301
      OPEN(13,FILE='PSATERROR.TXT',STATUS='OLD',BLANK='ZERO',ERR=101)
301   CONTINUE
      DO I = 1,ITER
          EP(I) = LOG(ABS (LOG(EP(I)) - LOG(PSAT)))
          EK(I) = LOG(ABS (LOG(EK(I)) - LOG(KFACT(CHECKCOMP))))
          WRITE(13,'(I10, F30.20, F30.20)') I, EP(I), EK(I)
      ENDDO
      CLOSE(13)
      
101   CONTINUE
      END SUBROUTINE

      SUBROUTINE GETSATTEMP(NCA, P, Z, BETA, TSAT, X, Y, IER, MAXITER)
C     NCA      no. of components      (i)
C     Z       feed mole fractions     (i)
C     P                    PRESSURE   (I)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NCA
      DOUBLE PRECISION TSAT,P, BETA
      DOUBLE PRECISION Z(NCA), ZNORM(NCA)
      INTEGER, OPTIONAL :: MAXITER
      
      INTEGER IER, ITER, CHECKCOMP
      DOUBLE PRECISION KFACT(NCA), X(NCA), Y(NCA)
      DOUBLE PRECISION FUGX(NCA), FUGY(NCA)
      DOUBLE PRECISION FUGTX(NCA), FUGTY(NCA)
      DOUBLE PRECISION F, DF
      DOUBLE PRECISION ERROR, MAXERROR
      DOUBLE PRECISION DT, DY(NCA), DX(NCA)
      DOUBLE PRECISION ET(1000), EK(1000)
      
      IER = 0
      MAXERROR = 1E-5
      
      CHECKCOMP = 1
      
      DO I=1,NCA
          ZNORM(I) = Z(I) / SUM(Z)
      ENDDO      
      
C     GET INITIAL K FACTOR ESTIMATES
      CALL GETINITT(NCA, P, ZNORM, BETA, TSAT)
      CALL WILSON_KFACT(NCA,P,TSAT,KFACT)
      CALL GETCOMPOSITION(NCA, ZNORM, BETA, KFACT, X, Y)
      
      ITER = 1
      DO WHILE(.TRUE.)
C         LOG THE VALUES
          ET(ITER) = TSAT
          EK(ITER) = KFACT(CHECKCOMP)
          
          
C         CHECK FUGACITY CEFFICIENTS
          CALL THERMO(TSAT, P, X, FUGX,
     &     FUGT = FUGTX)
C         CHECK FUGACITY CEFFICIENTS
          CALL THERMO(TSAT, P, Y, FUGY,
     &     FUGT = FUGTY)
          
C         CALCULAT THE K FACTORS
          DO I = 1,NCA
              KFACT(I) = EXP(FUGX(I) - FUGY(I))
          ENDDO
          
C         CALCULTE NEW COMPOSITION
          CALL GETCOMPOSITION(NCA, ZNORM, BETA, KFACT, DX, DY)
          DO I = 1,NCA
              TEMP = DY(I)
              DY(I) = TEMP - Y(I)              
              Y(I) = TEMP
              
              TEMP = DX(I)
              DX(I) = TEMP - X(I)
              X(I) = TEMP
          ENDDO
          
C         CALCULATE F AND DF/DT
          DF = 0
          F = -1
          DO I = 1,NCA
C              F = F +  (Y(I) - X(I))
              
              IF(BETA .GT. 1E-5) THEN
                  F = F + ZNORM(I) / KFACT(I)
                  DF = DF + ZNORM(I) / KFACT(I) * (FUGTY(I) - FUGTX(I))
              ELSE
                  F = F + ZNORM(I) * KFACT(I)
                  DF = DF + ZNORM(I) * KFACT(I) * (FUGTX(I) - FUGTY(I))
              ENDIF
              
          ENDDO
          
C         CALCULATE TEMPERATURE STEP
          DT = -F/DF
          IF( DT .NE. DT) THEN
              IER = 501
              GOTO 101
          ENDIF          
          TSAT = TSAT + DT
          
C         CHECK MAX NUMBER OF ITERATIONS
          IF(PRESENT(MAXITER)) THEN
              IF(ITER .EQ. MAXITER) GOTO 101
          ENDIF
                    
C         CHECK CONVERGENCE
          IF(ABS(DT) .LT. MAXERROR) THEN
              DO I=1,NCA
                  IF(ABS(DY(I)) .GE. MAXERROR) GOTO 100
                  IF(ABS(DX(I)) .GE. MAXERROR) GOTO 100
              ENDDO
              
              GOTO 200
          ENDIF
          
100       CONTINUE
          
          ITER = ITER + 1
      ENDDO

200   CONTINUE
      
      
      OPEN(13,FILE='TSATERROR.TXT',STATUS='NEW',BLANK='ZERO',ERR=401)
      GOTO 301
401   OPEN(13,FILE='TSATERROR.TXT',STATUS='OLD',BLANK='ZERO',ERR=101)      
301   CONTINUE
      DO I = 1,ITER
          ET(I) = LOG(ABS (LOG(ET(I)) - LOG(TSAT)))
          EK(I) = LOG(ABS (LOG(EK(I)) - LOG(KFACT(CHECKCOMP))))
          WRITE(13,'(I10, F30.20, F30.20)') I, ET(I), EK(I)
      ENDDO
      CLOSE(13)
      
101   CONTINUE
      END SUBROUTINE

      SUBROUTINE GETINITT(NCA, P, Z, BETA, T)
C     NCA      no. of components      (i)
C     Z       feed mole fractions     (i)
C     P, T    PRESSURE, TEMPERATURE   (I)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NCA
      DOUBLE PRECISION T,P, BETA
      DOUBLE PRECISION Z(NCA)
      CHARACTER*4 COMPS(NCA)
            
      INTEGER IER
      DOUBLE PRECISION ERROR, MAXERROR
      DOUBLE PRECISION PCALC, TSTEP
      
      MAXERROR = 1E-5
      ERROR = 1
      T = 1E-5
      TSTEP = 100
      
      DO WHILE(ABS(ERROR) .GT. MAXERROR)            
          
          CALL GETINITP(NCA, T, Z, BETA, PCALC)
          ERROR = P - PCALC          
          
          IF(ABS(ERROR) .LT. MAXERROR) THEN
              GOTO 101
          ELSE IF(ERROR .LT. 0) THEN  
              T = T - TSTEP
              TSTEP = TSTEP/2
          ENDIF
          
          T = T + TSTEP
      ENDDO
      
101   CONTINUE
          
      END SUBROUTINE

      SUBROUTINE GETINITP(NCA, T, Z, BETA, P)
C     NCA      no. of components      (i)
C     Z       feed mole fractions     (i)
C     P, T    PRESSURE, TEMPERATURE   (I)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NCA
      DOUBLE PRECISION T,P
      DOUBLE PRECISION TC, PC, OMEGA
      DOUBLE PRECISION Z(NCA)
      CHARACTER*4 COMPS(NCA)
      DOUBLE PRECISION TEMP
      
      INTEGER IER
      
      P = 0
      DO I=1,NCA
          CALL GETCRIT(I,TC,PC,OMEGA)
          
          TEMP = PC*EXP(5.373*(1+OMEGA)*(1-TC/T))
          
          IF(BETA .GT. 1E-5) THEN
              P = P + Z(I)/TEMP
          ELSE
              P = P + Z(I)*TEMP
          ENDIF
          
      ENDDO
      
      IF(BETA .GT. 1E-5) THEN
          P = 1 / P
      ENDIF
      
      END SUBROUTINE

      SUBROUTINE PRINTRSULT(HEADER, P, T, IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) HEADER
      DOUBLE PRECISION P, T
      INTEGER IER
      
      WRITE(*,*)
      WRITE(*,"(A50)") HEADER
      WRITE(*,"(A, F10.3, A)") "P = ", P, " MPa"
      WRITE(*,"(A, F10.3, A)") "T = ", T, " K"
	
      IF (IER .GT. 0) THEN          
          WRITE(*,*)
          WRITE(*,"(A, I5)") "AN ERROR OCCURED ", IER          
      END IF
      
      WRITE(*,*)
	  
      END SUBROUTINE
      
      END MODULE DAY6
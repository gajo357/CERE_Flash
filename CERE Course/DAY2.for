      MODULE DAY2
      
      USE PRINTFUNC
      USE THERMCALC
      
      PUBLIC PR_SOLVE, WILSON_KFACT
      CONTAINS
      
      SUBROUTINE RUN_DAY2(NCA, T, P, Z, COMPS)      
C     NCA      no. of components      (i)
C     Z       feed mole fractions     (i)
C     P, T    PRESSURE, TEMPERATURE   (I)
C     COMPS   LIST OF COMP NAMES      (I)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NCA
      DOUBLE PRECISION T,P
      DOUBLE PRECISION Z(NCA)
      CHARACTER*4 COMPS(NCA)
      
      INTEGER NITER, IER
      DOUBLE PRECISION BETA
      DOUBLE PRECISION X(NCA), Y(NCA), KFACT(NCA)
          
            
      CALL WILSON_KFACT(NCA,P,T,KFACT)
      
      CALL RR_SOLVE(NCA,Z,KFACT,BETA,X,Y)    
      CALL PRINT_RESULT("RR", P, T, 
     &                 NCA, Z, X, Y, BETA, COMPS, NITER, IER)  
      
      CALL SLOPPY_RR_SOLVE(NCA,Z,KFACT,BETA,X,Y)      
      CALL PRINT_RESULT("SLOPPY RR", P, T, 
     &                 NCA, Z, X, Y, BETA, COMPS, NITER, IER)
      
      CALL RR_NEGATIVE_SOLVE(NCA,Z,KFACT,BETA,X,Y)      
      CALL PRINT_RESULT("NEGATIVE FLASH RR", P, T, 
     &                 NCA, Z, X, Y, BETA, COMPS, NITER, IER)
      
      END SUBROUTINE
      
      SUBROUTINE WILSON_KFACT(NCA,P,T,KFACT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     NCA     NUMBER OF COMPONENTS
C     PC      ARRAY OF CRITICAL PRESSURES
C     TC      ARRAY OF CRITICAL TEMPERATURES
C     OMEGA     ACCENTRIC FACTORS
C     P, T    PRESSURE, TEMPERATURE
C     KFACT   LIST OF K-FACTORS (O)
      INTEGER NCA
      DOUBLE PRECISION PC(NCA), TC(NCA), OMEGA(NCA)
      DOUBLE PRECISION P, T
      DOUBLE PRECISION KFACT(NCA)
      
      DO I = 1,NCA
          CALL GETCRIT(I,TC(I),PC(I),OMEGA(I))
      END DO          
      
      DO I = 1,NCA
          KFACT(I) = PC(I) / P *
     &                EXP(5.373 * (1 + OMEGA(I)) * (1 - TC(I)/T)) 
      ENDDO
      
      END SUBROUTINE
      
      SUBROUTINE RR_SOLVE(NCA,Z,KFACT,BETA,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     NCA      no. of components (i)
C     Z       feed mole fractions (i)
C     KFACT   k-factors (i)
C     BETA    vapour fraction (o)
C     X       liquid mole fractions (o)
C     Y       vapour mole fractions (o)
      INTEGER NCA
      DOUBLE PRECISION Z(NCA)
      DOUBLE PRECISION KFACT(NCA)
      DOUBLE PRECISION BETA, DBETA, BETANEW, BETAMIN, BETAMAX
      DOUBLE PRECISION X(NCA)
      DOUBLE PRECISION Y(NCA)
      DOUBLE PRECISION G, DG, DET
      
      BETAMIN = 0
      BETAMAX = 1
      
C     CHECK INITIAL CONDITIONS (G(0) > 0 && G(1) < 0)
      CALL RR_GET_G(NCA, Z, KFACT, BETAMIN, G)
      IF (G .LE. 0) THEN
          BETA = BETAMIN
          GOTO 300
      END IF
      
      CALL RR_GET_G(NCA, Z, KFACT, BETAMAX, G)
      IF (G .GE. 0) THEN
          BETA = BETAMAX
          GOTO 300
      ENDIF
            
C     INITIAL GUESS FOR BETA
      BETA = (BETAMIN + BETAMAX)/2
      
C     JUST TO MAKE THE PROGRAM ENTER THE LOOP
      G = 1
      DBETA = 1
      DO WHILE (ABS(G) .GE. 1E-6 .OR. ABS(DBETA) .GE. 1E-6)
          
          CALL RR_GET_G(NCA, Z, KFACT, BETA, G)
          IF (G .LT. 0) THEN 
              BETAMAX = BETA
          ELSE IF (G .GT. 0) THEN
              BETAMIN = BETA
          ELSE 
              GOTO 100
          END IF
          
          CALL RR_GET_GPRIME(NCA, Z, KFACT, BETA, DG)
          
          BETANEW = BETA - G/DG
          
          IF (BETANEW .GT. BETAMAX .OR. BETANEW .LT. BETAMIN)
     &        BETANEW = (BETAMIN + BETAMAX)/2          
                  
c         WRITE(*,*) (BETA - BETANEW) / DBETA**2
          DBETA = BETA - BETANEW
          BETA = BETANEW                    
      ENDDO
      
C     CALCULATE VAPOR AND LIQUID COMPOSITION
100   DO I = 1,NCA
          CALL RR_DETERMINATOR(BETA, KFACT(I), DET) 
          X(I) = Z(I) / DET
          Y(I) = KFACT(I) * Z(I) / DET
      ENDDO
      
300   CONTINUE
      END SUBROUTINE
      
      SUBROUTINE SLOPPY_RR_SOLVE(NCA,Z,KFACT,BETA,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     NCA      no. of components (i)
C     Z       feed mole fractions (i)
C     KFACT   k-factors (i)
C     BETA    vapour fraction (o)
C     X       liquid mole fractions (o)
C     Y       vapour mole fractions (o)
      INTEGER NCA
      DOUBLE PRECISION Z(NCA), ZNORM(NCA)
      DOUBLE PRECISION KFACT(NCA)
      DOUBLE PRECISION BETA, DBETA, BETANEW, BETAMIN, BETAMAX
      DOUBLE PRECISION X(NCA)
      DOUBLE PRECISION Y(NCA)
      DOUBLE PRECISION G, DG, DET
      DOUBLE PRECISION V(NCA), L(NCA)
      
C     NORMALIZE COMPOSITION
      DO I = 1,NCA
          ZNORM(I) = Z(I) / 100         
      ENDDO
      
      BETAMIN = 0
      BETAMAX = 1
      
      DO I = 1,NCA
          IF (KFACT(I) .GT. 1) THEN
              BETANEW = (KFACT(I) * ZNORM(I) - 1) / (KFACT(I) - 1)
              IF (BETANEW .GT. BETAMIN .AND. BETANEW .LT. BETAMAX) THEN
                  BETAMIN = BETANEW
              END IF
              
          ELSE IF (KFACT(I) .LT. 1) THEN
              BETANEW = (1 - ZNORM(I)) / (1 - KFACT(I))
              IF (BETANEW .LT. BETAMAX .AND. BETANEW .GT. BETAMIN) THEN
                  BETAMAX = BETANEW
              END IF              
          END IF          
      END DO
            
C     INITIAL GUESS FOR BETA
      BETA = (BETAMIN + BETAMAX)/2
      
C     SINGLE NEWTON STEP
      CALL RR_GET_G(NCA, Z, KFACT, BETA, G)
      CALL RR_GET_GPRIME(NCA, Z, KFACT, BETA, DG)
      
      BETANEW = BETA - G/DG
      IF (BETANEW .GT. BETAMAX .OR. BETANEW .LT. BETAMIN)
     &        BETANEW = (BETAMIN + BETAMAX)/2          
      BETA = BETANEW   
      
C     Calculate vapour and liquid flows from
      DO I = 1,NCA
          CALL RR_DETERMINATOR(BETA, KFACT(I), DET) 
          L(I) = ((1 - BETA) * ZNORM(I)) / DET
          V(I) = (BETA * KFACT(I) * ZNORM(I)) / DET          
      ENDDO
C     Calculate the corrected value of the vapour fraction as
      BETA = 0
      DO I = 1,NCA
          BETA = BETA + V(I)
      ENDDO
      
C     CALCULATE VAPOR AND LIQUID COMPOSITION
100   DO I = 1,NCA
          X(I) = L(I) / (1 - BETA) * 100
          Y(I) = V(I) / BETA * 100
      ENDDO
      
300   CONTINUE
      END SUBROUTINE
      
      SUBROUTINE RR_NEGATIVE_SOLVE(NCA,Z,KFACT,BETA,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     NCA      no. of components (i)
C     Z       feed mole fractions (i)
C     KFACT   k-factors (i)
C     BETA    vapour fraction (o)
C     X       liquid mole fractions (o)
C     Y       vapour mole fractions (o)
      INTEGER NCA
      DOUBLE PRECISION Z(NCA)
      DOUBLE PRECISION KFACT(NCA), KMIN, KMAX
      DOUBLE PRECISION BETA, DBETA, BETANEW, BETAMIN, BETAMAX
      DOUBLE PRECISION X(NCA)
      DOUBLE PRECISION Y(NCA)
      DOUBLE PRECISION G, DG, DET
      
      BETAMIN = 0
      BETAMAX = 1
      
      KMIN = 1
      KMAX = 0
      DO I = 1,NCA
          IF (KFACT(I) .GT. KMAX) KMAX = KFACT(I)
          IF (KFACT(I) .LT. KMIN) KMIN = KFACT(I)    
      END DO
            
      IF (KMIN .LT. 1 .AND. KMAX .GT. 1) THEN
          BETAMIN = -1 / (KMAX - 1)
          BETAMAX = 1 / (1 - KMIN)
      END IF
            
C     INITIAL GUESS FOR BETA
      BETA = (BETAMIN + BETAMAX)/2
      
C     JUST TO MAKE THE PROGRAM ENTER THE LOOP
      G = 1
      DBETA = 1
      DO WHILE (ABS(G) .GE. 1E-6 .OR. ABS(DBETA) .GE. 1E-6)
          
          CALL RR_GET_G(NCA, Z, KFACT, BETA, G)
          IF (G .LT. 0) THEN 
              BETAMAX = BETA
          ELSE IF (G .GT. 0) THEN
              BETAMIN = BETA
          ELSE 
              GOTO 100
          END IF
          
          CALL RR_GET_GPRIME(NCA, Z, KFACT, BETA, DG)
          
          BETANEW = BETA - G/DG
          
          IF (BETANEW .GT. BETAMAX .OR. BETANEW .LT. BETAMIN)
     &        BETANEW = (BETAMIN + BETAMAX)/2          
                  
          DBETA = BETA - BETANEW
          BETA = BETANEW                    
      ENDDO
      
100   CALL GETCOMPOSITION(NCA, Z, BETA, KFACT, X, Y)
      
300   CONTINUE
      END SUBROUTINE
      
      SUBROUTINE GETCOMPOSITION(NCA, Z, BETA, KFACT, X, Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NCA
      DOUBLE PRECISION BETA, KFACT(NCA), DET
      DOUBLE PRECISION Z(NCA), X(NCA), Y(NCA)
      
C     CALCULATE VAPOR AND LIQUID COMPOSITION
      DO I = 1,NCA
          CALL RR_DETERMINATOR(BETA, KFACT(I), DET) 
          X(I) = Z(I) / DET
          Y(I) = KFACT(I) * Z(I) / DET
      ENDDO
      
      END SUBROUTINE
      
C     CALCULATES THE G FROM RACHFOR-RICE EQUATION
      SUBROUTINE RR_GET_G(NCA, Z, KFACT, BETA, G)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     NCA      no. of components (i)
C     Z       feed mole fractions (i)
C     KFACT   k-factors (i)
C     BETA    vapour fraction (I)
C     G       GIBBS ENERGY (O)
      INTEGER NCA
      DOUBLE PRECISION Z(NCA)
      DOUBLE PRECISION KFACT(NCA)
      DOUBLE PRECISION BETA, G, DET
      
      G = 0
      DO I = 1,NCA
          CALL RR_DETERMINATOR(BETA, KFACT(I), DET) 
          G = G + Z(I) * (KFACT(I) - 1)/DET
      ENDDO
      
      END SUBROUTINE
      
C     CALCULATES THE G PRIME FROM RACHFOR-RICE EQUATION
      SUBROUTINE RR_GET_GPRIME(NCA, Z, KFACT, BETA, DG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     NCA      no. of components (i)
C     Z       feed mole fractions (i)
C     KFACT   k-factors (i)
C     BETA    vapour fraction (I)
C     G       GIBBS ENERGY (O)
      INTEGER NCA
      DOUBLE PRECISION Z(NCA)
      DOUBLE PRECISION KFACT(NCA)
      DOUBLE PRECISION BETA, DG, DET
      
      DG = 0
      DO I = 1,NCA
          CALL RR_DETERMINATOR(BETA, KFACT(I), DET) 
          DG = DG + 
     &     Z(I) * ((KFACT(I) - 1)/DET)**2
      ENDDO
      
      DG = 0 - DG
      
      END SUBROUTINE
      
      SUBROUTINE RR_DETERMINATOR(BETA, KFACT, DET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION BETA, KFACT, DET      
      
      DET = 1 - BETA + BETA * KFACT
      END SUBROUTINE
      
      END MODULE
      
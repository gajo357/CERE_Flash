      MODULE DAY4
      
      USE PRINTFUNC
      USE THERMCALC
      USE DAY2
      
      PUBLIC RUN_DAY4, CHECKSTABILITY
      CONTAINS
      
      SUBROUTINE RUN_DAY4(NCA, T, P, Z)
C     NCA      no. of components      (i)
C     Z       feed mole fractions     (i)
C     P, T    PRESSURE, TEMPERATURE   (I)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NCA
      DOUBLE PRECISION T,P
      DOUBLE PRECISION Z(NCA)
      
      INTEGER IER
      DOUBLE PRECISION KFACT(NCA)
      DOUBLE PRECISION TRIAL(2)
      
      CALL CHECKSTABILITY(NCA, Z, P, T, KFACT, TRIAL, IER)
      CALL PRINT_STABILITY_RESULT("STABILITY ANALYSIS", 
     &                             P, T, TRIAL, IER)
      
      END SUBROUTINE
      
      SUBROUTINE CHECKSTABILITY(NCA, Z, P, T, KFACT, TRIAL, IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     NCA      no. of components (i)
C     Z       feed mole fractions (i)
C     KFACT   k-factors (i)
C     P, T    PRESSURE, TEMPERATURE
      INTEGER NCA, IER
      DOUBLE PRECISION Z(NCA), ZNORM(NCA), P, T
      DOUBLE PRECISION KFACT(NCA)
      DOUBLE PRECISION WLIGHT(NCA), WHEAVY(NCA)
      DOUBLE PRECISION TRIAL(2)
      DOUBLE PRECISION D(NCA), TMH, TML
      DOUBLE PRECISION FUG(NCA), FUGL(NCA), FUGV(NCA)
      
C     TRIVIAL SOLUTION
      TRIAL(1) = 0
      TRIAL(2) = 0
      
      CALL THERMO(T,P,Z,FUG)
            
C     Get initial K-estimate
      CALL WILSON_KFACT(NCA,P,T,KFACT)
100   DO I = 1,NCA
          ZNORM(I) = Z(I) / 100
C          WLIGHT(I) = KFACT(I) * ZNORM(I)
C          WHEAVY(I) = ZNORM(I) / KFACT(I)
          WLIGHT(I) = ZNORM(I) * EXP(EXP(FUG(I)))
          WHEAVY(I) = ZNORM(I) * EXP(EXP(FUG(I))) / KFACT(I)
      ENDDO
      
C     GET MINIMUM FOR LIGHT AND HEAVY ESTIMATE
      CALL CHECKCONDITION(NCA, ZNORM, WLIGHT, P, T, TML, IER)
      CALL CHECKCONDITION(NCA, ZNORM, WHEAVY, P, T, TMH, IER)
      
C     IF NEGATIVE, THEN RESTART
      CALL THERMO(T,P,WLIGHT,FUGV)
      CALL THERMO(T,P,WHEAVY,FUGL)
      IF (TML .LT. 0 .AND. TMH .LT. 0) THEN
          DO I = 1,NCA
              KFACT(I) = EXP(FUGL(I)) / EXP(FUGV(I))                    
          ENDDO
      ELSE IF (TML .LT. 0) THEN
          DO I = 1,NCA
              KFACT(I) = EXP(FUG(I)) / EXP(FUGV(I))
          ENDDO
      ELSE IF (TMH .LT. 0) THEN
          DO I = 1,NCA
              KFACT(I) = EXP(FUGL(I)) / EXP(FUG(I))
          ENDDO
      ENDIF       
      
      TRIAL(1) = TML
      TRIAL(2) = TMH
      
101   CONTINUE      
      END SUBROUTINE
      
      SUBROUTINE CHECKCONDITION(NCA, Z, W, P, T, TM, IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      INTEGER NCA, IER
      DOUBLE PRECISION Z(NCA), P, T
      DOUBLE PRECISION KFACT(NCA)
      DOUBLE PRECISION BETA
      DOUBLE PRECISION W(NCA)
      DOUBLE PRECISION TRIAL(2)
      DOUBLE PRECISION FUG(NCA), FUGW(NCA)
      DOUBLE PRECISION D(NCA), TM, TM_1, DTM
      DOUBLE PRECISION EMARGIN
      
      EMARGIN = 1E-4
      
      TM = 100000
      
      CALL THERMO(T,P,Z,FUG)
      
      DO WHILE (.TRUE.)
C          Re-evaluate fugacity coefficients from
          CALL THERMO(T,P,W,FUGW)
          
          TM_1 = TM
          TM = 0
          DO I = 1,NCA
              D(I) = FUG(I) + LOG (Z(I))
              TM = TM + W(I) * (LOG (W(I)) + FUGW(I) - D(I) - 1)
          ENDDO
                    
          TM = TM + 1
          
C          IF (TM .LT. 0) THEN
C              GOTO 100
C          ENDIF
          
          DTM = ABS (TM - TM_1)
          IF (DTM .LT. 1E-10) THEN
              GOTO 100
          ENDIF
          
C         NEXT STEP
          DO I = 1,NCA
              W(I) = EXP(D(I) - FUGW(I))
          ENDDO
      ENDDO
      
100   CONTINUE          
      
C     CHECK FOR NON-TRIVIAL SOLUTION
      DO I = 1,NCA
          IF( (Z(I) - W(I))**2 .GT. EMARGIN) THEN
              GOTO 200
          ENDIF
      ENDDO
      
C     IT'S A TRIVIAL SOLUTION
      TM = 0
      GOTO 101
      
C     NON TRIVIAL SOLUTION
200   CONTINUE
     

C     END
101   CONTINUE      
      END SUBROUTINE
      
      
      END MODULE DAY4
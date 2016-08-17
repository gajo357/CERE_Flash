      MODULE DAY7
      
      USE THERMCALC
      USE DAY2
            
      PUBLIC RUN_DAY7
      CONTAINS 
      
      SUBROUTINE RUN_DAY7(NCA, T, P, Z, LIST, COMPS)
C     NCA      no. of components      (i)
C     Z       feed mole fractions     (i)
C     P, T    PRESSURE, TEMPERATURE   (I)
C     LIST    LIST OF COMPONENT INDEXES
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      PARAMETER (NF = 3)
      
      INTEGER NCA, IER, LIST(NCA)
      DOUBLE PRECISION T,P
      DOUBLE PRECISION Z(NCA)
      CHARACTER*4 COMPS(NCA)
      DOUBLE PRECISION BETA(NF), FUGFAC(NCA, NF), PHAZ(NCA, NF)
      DOUBLE PRECISION KFACT(NCA)
      
      CALL WILSON_KFACT(NCA, P, T, KFACT)
      DO K = 1,NF
          DO I = 1,NCA
              IF(K .EQ. 1) THEN
                  FUGFAC(I, K) = 1
              ELSE
                  FUGFAC(I, K) = KFACT(I)
C                 H2S IN METHANE RICH LIQUID PHASE
                  IF(K .EQ. 2 .AND. LIST(I) .EQ. 15) THEN
                      FUGFAC(I, K) = EXP(LOG(KFACT(I)) + 1)
                  END IF
          
C                 METHANE IN H2S RICH LIQUID PHASE
                  IF(K .EQ.3 .AND. LIST(I) .EQ. 1) THEN
                      FUGFAC(I, K) = EXP(LOG(KFACT(I)) + 1)
                  END IF      
              ENDIF
          ENDDO
          
          BETA(K) = 1.0/NF
      ENDDO
      
      
      CALL MULTRACF(NCA, NF, Z, BETA, FUGFAC, PHAZ, IER)
      WRITE(*,"(A20, A20, A20)") "VAPOR", "LIQUID 1", "LIQUID 2"
      WRITE(*,"(F20.15, F20.15, F20.15)") BETA(1), BETA(2), BETA(3)
      
      END SUBROUTINE

      SUBROUTINE MULTRACF (NCA,NF,Z,BETA,FUGFAC, PHAZ, IER)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     NC: NO. OF COMPONENTS
C     NF: NO. OF PHASES
C     BETA: VECTOR OF PHASE FRACTIONS
C     FUGFAC: INVERSE FUGACITY COEFFICIENT ARRAY
      INTEGER NCA, NF, IER
      DOUBLE PRECISION Z(NCA), BETA(NF), FUGFAC(NCA, NF), PHAZ(NCA, NF)
      
      DOUBLE PRECISION Q, G(NF), H(NF, NF)
      DOUBLE PRECISION G_1(NF), E(NCA), QNEW
      DOUBLE PRECISION ALPHA, ALPHANEW, BETANEW(NF), DBETA(NF)
      DOUBLE PRECISION QMARGIN
      INTEGER INDX(NF), KZERO
      LOGICAL PHAACT(NF)
      
      IER = 0
      QMARGIN = 1E-10
      
      ALPHA = 1
      
C     STEP 1: MARK ACTIVE PHASES
      PHAACT(:) = .FALSE.
      DO K = 1,NF
          IF(BETA(K) .GT. 0.0) PHAACT(K) = .TRUE.
      ENDDO
      
C     STEP 2: CALCULATE GRADIENT AND HESSIAN
200   CALL CALCQ(NCA, NF, Z, BETA, FUGFAC, PHAACT, Q, G, H)
C     SAVE THE COPY OF GRADIENT VALUES
      G_1(:) = G(:)
      
C     STEP 3: CALCULATE THE NEWTON STEP, VALUES ARE STORED IN G
      CALL LUDEC(NF, INDX, H)
      CALL BACKSUBST(NF, INDX, H, G)
      DBETA(:) = G(:)
      
C     STEP 4: CALCULATE THE NEW BETA VECTOR
      DO K = 1,NF
          BETANEW(K) = BETA(K) - ALPHA * DBETA(K)
      ENDDO
C     CHECK IF ALL BETANEW VALUES ARE POSITIVE
      DO K = 1,NF
          IF (BETANEW(K) .LT. 0) GOTO 401          
      ENDDO 
C     ALL NEW BETA'S ARE POSITIVE, GO TO STEP 5
      GOTO 500
      
C     FIND ALPHA, RECALCULATE Q
401   ALPHANEW = ALPHA
      CALL FINDALPHA(NF, BETA, BETANEW, DBETA, ALPHANEW, KZERO)
      CALL CALCQ(NCA, NF, Z, BETANEW, FUGFAC, PHAACT, QNEW)
      IF (QNEW .LT. Q + QMARGIN) THEN
C         MARK THE ELEMENT INACTIVE ,SET Q AND BETA, GO TO STEP 2
          PHAACT(KZERO) = .FALSE.
          ALPHA = ALPHANEW
          Q = QNEW
          BETA(:) = BETANEW(:)
          GOTO 200
      ELSE
C         STEP NOT EXCEPTED, HALF THE STEP, RECALCULATE
          ALPHA = ALPHANEW / 2
          DO K = 1,NF
              BETANEW(K) = BETA(K) - ALPHA * DBETA(K)
          ENDDO
      ENDIF      
      
C     FOUND VALID ALPHA
C     STEP 5: CALCULATE NEW Q
500   CALL CALCQ(NCA, NF, Z, BETANEW, FUGFAC, PHAACT, QNEW)
      
C     STEP 6: CHECK NEW Q
      IF (QNEW .LT. Q + QMARGIN) THEN
C         STEP IS ACCEPTED
          Q = QNEW
          BETA(:) = BETANEW(:)
C         CHECK THE CONVERGENCE
          DO K = 1,NF
              IF (ABS(DBETA(K)) .GE. QMARGIN) GOTO 200
          ENDDO          
C         CONVERGED, GO TO STEP 7
          GOTO 700
      ELSE
C         HALF THE ALPHA, GO TO STEP 5
          ALPHA = ALPHA / 2
          DO K = 1,NF
              BETANEW(K) = BETA(K) - ALPHA * DBETA(K)
          ENDDO          
          GOTO 500
      ENDIF

C     STEP 7: CHECK THE GRADIENT VECTOR
700   DO K = 1,NF
          IF(PHAACT(K) .EQ. .FALSE.) THEN
              IF(G_1(K) .LT. 0) THEN
C         IF GRADIENT OF INACTIVE PHASE IS NEGATIVE
C         ACTIVATE THE PHASE AND REDO FROM STEP 2
                  PHAACT(K) = .TRUE.
                  GOTO 200
              ENDIF              
          ENDIF          
      ENDDO
              
C     STEP 8: CALCULATE PHASE COMPOSITIONS
      DO K = 1,NF
          DO I = 1,NCA
              PHAZ(I, K) = Z(I)*BETA(K)
          ENDDO          
      ENDDO
      
      
      END SUBROUTINE 
      
      SUBROUTINE CALCQ(NCA, NF, Z, BETA, FUGFAC, PHAACT, Q, G, H)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     NC: NO. OF COMPONENTS
C     NF: NO. OF PHASES
C     BETA: VECTOR OF PHASE FRACTIONS
C     FUGFAC: INVERSE FUGACITY COEFFICIENT ARRAY
      INTEGER NCA, NF
      DOUBLE PRECISION Z(NCA), BETA(NF), FUGFAC(NCA, NF), E(NCA)
      DOUBLE PRECISION, OPTIONAL :: G(NF), H(NF, NF)
      LOGICAL PHAACT(NF)
                       
      Q = SUM(BETA)
      DO I = 1,NCA
          E(I) = 0
          DO K = 1,NF
              E(I) = E(I) + BETA(K)/FUGFAC(I, K)
          ENDDO
C         Q = SUM(BETA) - SUM(Z * ln(E))
          Q = Q - Z(I) * LOG(E(I))
      ENDDO
      
      IF(PRESENT(G) .AND. PRESENT(H)) THEN
          DO K = 1,NF                            
              IF(PHAACT(K)) THEN
C             CALCULATE GRADIENT
                  G(K) = 1
                  DO I = 1,NCA
                      G(K) = G(K) - Z(I)/(E(I) * FUGFAC(I, K))
                  ENDDO
C             CALCULATE HESSIAN
                  DO L = 1, NF
                      H(K,L) = 0
                      DO I = 1,NCA
                          H(K, L) = H(K, L) + 
     &                    Z(I)/(E(I)**2 * FUGFAC(I, L) * FUGFAC(I, K))
                      ENDDO
                  ENDDO
              ELSE
C             INACTIVE PHASE
                  G(K) = 0
                  DO L = 1, NF
                      H(K, L) = 0
                      H(L, K) = 0
                  ENDDO                  
C                 DECOUPLE PHASE FROM THE SET
                  H(K, K) = 1
              ENDIF
          ENDDO
          
      ENDIF
      
      
      END SUBROUTINE 
      
      SUBROUTINE FINDALPHA(NF, BETA, BETANEW, DBETA, ALPHA, KZERO)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     NF: NO. OF PHASES
C     BETA: VECTOR OF PHASE FRACTIONS
C     FUGFAC: INVERSE FUGACITY COEFFICIENT ARRAY
      INTEGER NF, KZERO
      DOUBLE PRECISION BETA(NF), BETANEW(NF), DBETA(NF), ALPHA
      
C     ONE NEGATIVE BETA HAS TO BECOME 0, OTHERS POSITIVE
      KZERO = 1
402   DO K = KZERO,NF
          IF (BETANEW(K) .LT. 0) THEN
              ALPHA = BETA(K) / DBETA(K)
              KZERO = K
              GOTO 403
          ENDIF
      ENDDO      
      
C     CHECK IF THIS ALPHA DID THE TRICK
403   IF(ALPHA .GE. 1) THEN
          KZERO = KZERO + 1
          GOTO 402
      ENDIF
      
      DO K = 1,NF
          IF((BETA(K) - ALPHA * DBETA(K)) .LT. 0) THEN
              KZERO = KZERO + 1
              GOTO 402   
          ENDIF   
      ENDDO
      
C     ALL BETANEW WILL BE NON-NEGATIVE FOR THIS ALPHA
      DO K = 1,NF
          BETANEW(K) = BETA(K) - ALPHA * DBETA(K)
      ENDDO
      
      END SUBROUTINE
      
      END MODULE DAY7
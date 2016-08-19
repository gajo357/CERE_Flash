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
      
      CALL MULTI(NCA, NF, T, P, Z, LIST, COMP)
      
      END SUBROUTINE

      SUBROUTINE MULTI(NCA, NF, T, P, Z, LIST, COMP)
C     NCA      no. of components      (i)
C     Z       feed mole fractions     (i)
C     P, T    PRESSURE, TEMPERATURE   (I)
C     LIST    LIST OF COMPONENT INDEXES
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      INTEGER NCA, NF, IER, LIST(NCA)
      DOUBLE PRECISION T,P
      DOUBLE PRECISION Z(NCA)
      CHARACTER*4 COMPS(NCA)
      
      DOUBLE PRECISION BETA(NF), BETA0(NF), BETANEW(NF)
      DOUBLE PRECISION FUGFAC(NCA, NF), Y(NCA, NF), FUG(NCA)
      DOUBLE PRECISION FUGFACNEW(NCA, NF), EFUG(NF), EMARGIN
      DOUBLE PRECISION FUGFAC_1(NCA, NF), FUGFAC_2(NCA, NF)
      DOUBLE PRECISION DMARGIN
      INTEGER ITER, CHECKSTEP
      
      EMARGIN = 1E-7
      DMARGIN = 1E-7
      CHECKSTEP = 5
      DO K = 1,NF
          EFUG(K) = 1
      ENDDO      
      
C     INITIAL ESTIMATES
      CALL MULTIWILSON(NCA, NF, T, P, Z, LIST, FUGFAC, BETA)
      
      NITER = 1      
      DO WHILE (.TRUE.)
C         ACCURACY SHOULD BE GREATER AFTER EACH ACCELARATION AND THE FIRST STEP
          IF (MOD(NITER, CHECKSTEP) .EQ. 1) THEN
              DMARGIN = 1E-10
          ELSE
              DMARGIN = 1E-5
          ENDIF
          
C         CALL RR ALGORITHM
          CALL MULTRACF(NCA, NF, Z, BETA, FUGFAC, DMARGIN, Y, IER)
C         SAVE INITIAL BETA VALUES FOR PRINTING
          IF(NITER .EQ. 1) BETA0(:) = BETA(:)
          
C         RECALCULATE FUGACITIES BASED ON CALCULATED COMPOSITIONS
          DO K = 1,NF
              CALL THERMO(T, P, Y(:,K), FUG)
              FUGFACNEW(:,K) = EXP(FUG(:))
          ENDDO
          
C         CHECK FOR CONVERGENCE
          DO I = 1,NCA              
              DO K = 1,NF
                  ETEMP = ABS(
     &                ((FUGFACNEW(I, K)) - (FUGFAC(I, K)))
     &                    / (FUGFAC(I, K)) )
                  IF (ETEMP < EFUG(K)) EFUG(K) = ETEMP
              ENDDO
          ENDDO
          
          DO K = 1,NF
              IF(ABS(EFUG(K)) .GE. EMARGIN) GOTO 100
          ENDDO
          GOTO 200
          
C         PREPARE THE NEXT STEP
100       FUGFAC_2(:,:) = FUGFAC_1(:,:)
          FUGFAC_1(:,:) = FUGFAC(:,:)
          FUGFAC(:,:) = FUGFACNEW(:,:)
          
C         ACCELARATE THE CALCULATION
          IF (MOD(NITER, CHECKSTEP) .EQ. 0) THEN
              CALL FUGACC(NCA, NF, T, P, Z, Y, DMARGIN,
     &                    FUGFAC, FUGFAC_1, FUGFAC_2, FUGFACNEW, BETA)
          ENDIF
                    
          NITER = NITER + 1
      ENDDO

200   CONTINUE
      WRITE(*,*) "NUMBER OF ITERATIONS ", NITER
      WRITE(*,"(A20, A20, A20)") "VAPOR", "LIQUID 1", "LIQUID 2"
      WRITE(*,"(F20.15, F20.15, F20.15)") BETA0(1), BETA0(2), BETA0(3)
      WRITE(*,"(F20.15, F20.15, F20.15)") BETA(1), BETA(2), BETA(3)  
      
101   CONTINUE
      END SUBROUTINE

      SUBROUTINE MULTRACF(NCA,NF,Z,BETA,FUGFAC, DMARGIN, Y, IER)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     NC: NO. OF COMPONENTS
C     NF: NO. OF PHASES
C     BETA: VECTOR OF PHASE FRACTIONS
C     FUGFAC: INVERSE FUGACITY COEFFICIENT ARRAY
      INTEGER NCA, NF, IER
      DOUBLE PRECISION Z(NCA), BETA(NF), FUGFAC(NCA, NF), Y(NCA, NF)
      
      DOUBLE PRECISION Q, G(NF), H(NF, NF)
      DOUBLE PRECISION E(NCA), QNEW
      DOUBLE PRECISION ALPHA, ALPHANEW, BETANEW(NF), DBETA(NF)
      DOUBLE PRECISION QMARGIN, DMARGIN
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
200   CALL CALCQ(NCA, NF, Z, BETA, FUGFAC, PHAACT, E, Q, G, H)
      
C     STEP 3: CALCULATE THE NEWTON STEP, VALUES ARE STORED IN G
      DBETA(:) = G(:)      
      CALL LUDEC(NF, INDX, H)
      CALL BACKSUBST(NF, INDX, H, DBETA)
      
C     STEP 4: CALCULATE THE NEW BETA VECTOR
      DO K = 1,NF
          IF (PHAACT(K) .EQ. .FALSE.) DBETA(K) = 0
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
      CALL CALCQ(NCA, NF, Z, BETANEW, FUGFAC, PHAACT, E, QNEW)
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
500   IF (ALPHA .LT. QMARGIN) GOTO 700
      CALL CALCQ(NCA, NF, Z, BETANEW, FUGFAC, PHAACT, E, QNEW)
      
C     STEP 6: CHECK NEW Q
      IF (QNEW .LT. Q + QMARGIN) THEN
C         STEP IS ACCEPTED
          Q = QNEW
          BETA(:) = BETANEW(:)
C         CHECK THE CONVERGENCE, GO TO STEP 2 IF NOT
          DO K = 1,NF
              IF (ABS(DBETA(K)) .GE. DMARGIN) GOTO 200
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
700   CALL CALCQ(NCA, NF, Z, BETA, FUGFAC, .NOT.PHAACT, E, Q, G, H)
      DO K = 1,NF
          IF(PHAACT(K) .EQ. .FALSE.) THEN
              IF(G(K) .LT. 0) THEN
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
              Y(I, K) = Z(I)/E(I)/FUGFAC(I, K)
          ENDDO          
      ENDDO
      
      
      END SUBROUTINE 
      
      SUBROUTINE CALCQ(NCA, NF, Z, BETA, FUGFAC, PHAACT, E, Q, G, H)
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
          IF((K .NE. KZERO) .AND. 
     &         (BETA(K) - ALPHA * DBETA(K)) .LT. 0) THEN
              KZERO = KZERO + 1
              GOTO 402   
          ENDIF   
      ENDDO
      
C     ALL BETANEW WILL BE NON-NEGATIVE FOR THIS ALPHA
      DO K = 1,NF
          BETANEW(K) = BETA(K) - ALPHA * DBETA(K)
          IF (K .EQ. KZERO) BETANEW(K) = 0
      ENDDO
      
      END SUBROUTINE
      
      SUBROUTINE MULTIWILSON(NCA, NF, T, P, Z, LIST, FUGFAC, BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      INTEGER NCA, NF, LIST(NCA)
      DOUBLE PRECISION T,P
      DOUBLE PRECISION Z(NCA)
      DOUBLE PRECISION BETA(NF), KFACT(NCA)
      DOUBLE PRECISION FUGFAC(NCA, NF)
      
      
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
      END SUBROUTINE
      
      SUBROUTINE FUGACC(NCA, NF, T, P, Z, Y, DMARGIN,
     &     FUGFAC, FUGFAC_1, FUGFAC_2, FUGFACNEW, BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      INTEGER NCA, NF, LIST(NCA), IER
      DOUBLE PRECISION T,P
      DOUBLE PRECISION Z(NCA)
      DOUBLE PRECISION BETA(NF), BETANEW(NF)
      DOUBLE PRECISION FUGFAC(NCA, NF), Y(NCA, NF), FUG(NCA)
      DOUBLE PRECISION FUGFAC_1(NCA, NF), FUGFAC_2(NCA, NF)
      DOUBLE PRECISION FUGFACNEW(NCA, NF)
      DOUBLE PRECISION LAMBDA, L1, L2, D(NCA), D1(NCA)
      DOUBLE PRECISION DMARGIN, G1, G2, DG
      
C     CHECK NEW FUGACITY FACTORS IF THEY DECREASE THE GIBBS ENERGY
C     GIBBS ENERGY BEFORE ACCELARATION
      G1 = 0
      DO I = 1,NCA
	    DO K = 1,NF
	        G1 = G1 + Y(I, K) * 
     &                    (LOG(Y(I, K)) + LOG(FUGFAC(I, K)))
	    ENDDO
      ENDDO
      
C     ACCELARATION      
      DO K = 1,NF
	    L1 = 0
	    L2 = 0
	    DO I = 1,NCA
	        D(I) = LOG(FUGFAC_1(I, K)) - LOG(FUGFAC_2(I, K))
	        D1(I) = LOG(FUGFAC(I, K)) - LOG(FUGFAC_1(I, K))
	  
	        L1 = L1 + D1(I)**2
	        L2 = L2 + D(I) * D1(I)
	    ENDDO
	    LAMBDA = L1 / L2

C         USE THERMCALC CORRECTED K FACTORS              
	    DO I = 1,NCA
	      FUGFACNEW(I, K) = EXP (LOG(FUGFAC(I, K)) + 
     &                           (D1(I) * LAMBDA) / (1 - LAMBDA))
	    ENDDO
      ENDDO

C     GIBBS ENERGY AFTER ACCELARATION
      BETANEW(:) = BETA(:)
      CALL MULTRACF(NCA, NF, Z, BETANEW, FUGFACNEW,
     &                        DMARGIN, Y, IER)
      G2 = 0
      DO K = 1,NF
	    CALL THERMO(T, P, Y(:,K), FUG)
          DO I = 1,NCA
                G2 = G2 + Y(I, K) * 
     &                    (LOG(Y(I, K)) + FUG(I))
          ENDDO

          FUGFACNEW(:,K) = EXP(FUG(:))
      ENDDO

C     CALCULATE THE DELTA GIBBS
      DG = G2 - G1
C     IF DELTA GIBBS IS NEGATIVE, WE ACCEPT ACCELARATION
      IF (DG .LT. 0) THEN
          FUGFAC(:,:) = FUGFACNEW(:,:)
      ENDIF
      
      END SUBROUTINE
      
      END MODULE DAY7
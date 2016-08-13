      MODULE DAY3
      
      USE PRINTFUNC
      USE THERMCALC
      USE DAY2
      USE DAY4
      
      PUBLIC RUN_DAY3, PTFLASH, PTFLASHACC
      CONTAINS
      
      SUBROUTINE RUN_DAY3(NCA, T, P, Z, COMPS)
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
      DOUBLE PRECISION X(NCA), Y(NCA)
      
      
      CALL PTFLASH(NCA, Z, P, T, BETA, X, Y, NITER, IER)     
      CALL PRINT_RESULT("PT FLASH", P, T, 
     &                 NCA, Z, X, Y, BETA, COMPS, NITER, IER)
      
      CALL PTFLASHACC(NCA, Z, P, T, BETA, X, Y, NITER, IER)     
      CALL PRINT_RESULT("PT FLASH ACCELERATED", P, T, 
     &                 NCA, Z, X, Y, BETA, COMPS, NITER, IER)
      
      END SUBROUTINE
      
      SUBROUTINE PTFLASH(NCA, Z, P, T, BETA, X, Y, NITER, IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     NCA      no. of components (i)
C     Z       feed mole fractions (i)
C     KFACT   k-factors (i)
C     P, T    PRESSURE, TEMPERATURE
C     BETA    vapour fraction (o)
C     X       liquid mole fractions (o)
C     Y       vapour mole fractions (o)
      INTEGER NCA, NITER, IER
      DOUBLE PRECISION Z(NCA)
      DOUBLE PRECISION KFACT(NCA), KFACTNEW(NCA)
      DOUBLE PRECISION BETA
      DOUBLE PRECISION PC(NCA), TC(NCA), OMEGA(NCA)
      DOUBLE PRECISION X(NCA)
      DOUBLE PRECISION Y(NCA)
      DOUBLE PRECISION EK, EFUGV, EFUGL, ETEMP
      DOUBLE PRECISION FUGL(NCA), FUGV(NCA), FUGLNEW(NCA), FUGVNEW(NCA)
      DOUBLE PRECISION TRIAL(2)
      
      EK = 1
      EFUGV = 1
      EFUGL = 1
      
C     CHECK THE STABILITY
      CALL CHECKSTABILITY(NCA, Z, P, T, KFACT, TRIAL, IER)
      
C     Get initial K-estimate
C      CALL WILSON_KFACT(NCA,P,T,KFACT)
            
      NITER = 1
      DO WHILE ((EFUGV .GE. 1E-7) .OR. (EFUGL .GE. 1E-7))
C      DO WHILE (EK .GE. 1E-7)
                    
C         SOLVE FOR BETA
          CALL RR_SOLVE(NCA,Z,KFACT,BETA,X,Y)
          
C         Re-evaluate fugacity coefficients from
          CALL THERMO(T,P,X,FUGLNEW)
          CALL THERMO(T,P,Y,FUGVNEW)
C         Calculate new K-factors
          DO I = 1,NCA
              KFACTNEW(I) = EXP(FUGLNEW(I)) / EXP(FUGVNEW(I))
          ENDDO
          
C         Check the new FUGACITIES against the old
          EFUGL = 1
          EFUGV = 1
          DO I = 1,NCA
              ETEMP = ABS(
     &                (EXP(FUGLNEW(I)) - EXP(FUGL(I))) / EXP(FUGL(I)) )
              IF (ETEMP < EFUGL) EFUGL = ETEMP
              
              ETEMP = ABS( 
     &                (EXP(FUGVNEW(I)) - EXP(FUGV(I))) / EXP(FUGV(I)) )
              IF (ETEMP < EFUGV) EFUGV = ETEMP
          ENDDO
          
C         Check the new K against the old
          EK = 1
          DO I = 1,NCA
              ETEMP = ABS( LOG(KFACTNEW(I)) - LOG(KFACT(I)))
              IF (ETEMP < EK) EK = ETEMP
          ENDDO
          
C          WRITE(*, "(I10, E15.8, E15.8, E15.8)") ITER, EFUGL, EFUGV, EK
          
C         ASSIGN THE NEW K FACTORS AND FUGACITIES
          DO I = 1,NCA
              KFACT(I) = KFACTNEW(I)
              
              FUGL(I) = FUGLNEW(I)
              FUGV(I) = FUGVNEW(I)
          ENDDO                   
              
          NITER = NITER + 1
      ENDDO
      
      
      
      END SUBROUTINE

      SUBROUTINE PTFLASHACC(NCA, Z, P, T, BETA, X, Y, NITER, IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     NCA      no. of components (i)
C     Z       feed mole fractions (i)
C     KFACT   k-factors (i)
C     P, T    PRESSURE, TEMPERATURE
C     BETA    vapour fraction (o)
C     X       liquid mole fractions (o)
C     Y       vapour mole fractions (o)
      PARAMETER (CHECKSTEP = 5)
      INTEGER NCA, ITER, NITER, IER
      DOUBLE PRECISION Z(NCA)
      DOUBLE PRECISION KFACT(NCA), KFACTNEW(NCA)
      DOUBLE PRECISION BETA
      DOUBLE PRECISION PC(NCA), TC(NCA), OMEGA(NCA)
      DOUBLE PRECISION X(NCA)
      DOUBLE PRECISION Y(NCA)
      DOUBLE PRECISION EK, ETEMP
      DOUBLE PRECISION KFACTHIST(NCA, CHECKSTEP)
      DOUBLE PRECISION FUG(NCA), FUGL(NCA), FUGV(NCA)
      DOUBLE PRECISION D(NCA), D1(NCA), LAMBDA, L1, L2
      DOUBLE PRECISION OF, OF1, DOF
      DOUBLE PRECISION V(NCA), L(NCA), S1, S2, S3
                  
      EK = 1
      
      CALL THERMO(T,P,Z,FUG)
      
C     Get initial K-estimate
      CALL WILSON_KFACT(NCA,P,T,KFACT)
            
      NITER = 1
      ITER = 1
      DO WHILE (EK .GE. 1E-7)
          
C         SOLVE FOR BETA
          CALL RR_SOLVE(NCA,Z,KFACT,BETA,X,Y)
          
C         Re-evaluate fugacity coefficients from
          CALL THERMO(T,P,X,FUGL)
          CALL THERMO(T,P,Y,FUGV)
C         Calculate new K-factors
          DO I = 1,NCA
              KFACTNEW(I) = EXP(FUGL(I)) / EXP(FUGV(I))
          ENDDO
          
C         Check the new K against the old
          EK = 1
          DO I = 1,NCA
              ETEMP = ABS( LOG(KFACTNEW(I)) - LOG(KFACT(I)))
              IF (ETEMP < EK) EK = ETEMP
          ENDDO
                   
C         ASSIGN THE NEW K FACTORS
          DO I = 1,NCA
              KFACT(I) = KFACTNEW(I)
              
              KFACTHIST(I, ITER) = LOG( KFACT(I) )
          ENDDO
          
C         CORRECT THE K FACTORS
          IF (ITER .EQ. CHECKSTEP) THEN
              L1 = 0
              L2 = 0
              DO I = 1,NCA
                  D(I) = KFACTHIST(I, ITER - 1) - KFACTHIST(I, ITER - 2)
                  D1(I) = KFACTHIST(I, ITER) - KFACTHIST(I, ITER - 1)
                      
                  L1 = L1 + D1(I)**2
                  L2 = L2 + D(I) * D1(I)
              ENDDO
              LAMBDA = L1 / L2
              
C             USE THERMCALC CORRECTED K FACTORS              
              DO I = 1,NCA
                  KFACT(I) = EXP (KFACTHIST(I, ITER) + 
     &                     (D1(I) * LAMBDA) / (1 - LAMBDA))
              ENDDO              
          END IF

C         Calculate vapour and liquid flows from
          S1 = 0
          S2 = 0
          S3 = 0
          DO I = 1,NCA
              CALL RR_DETERMINATOR(BETA, KFACT(I), DET) 
              L(I) = ((1 - BETA) * Z(I) / 100) / DET
              V(I) = (BETA * KFACT(I) * Z(I) / 100) / DET          
              
              S1 = S1 - Z(I)/100 * (LOG(Z(I)/100) + (FUG(I)))
              S2 = S2 + V(I) * (LOG(Y(I)/100) + (FUGV(I)))
              S3 = S3 + L(I) * (LOG(X(I)/100) + (FUGL(I)))
          ENDDO
          DOF = OF
          OF = S1 + S2 + S3
          
          IF ( (OF - DOF) .LT. 0) THEN
              IER = 5
              GOTO 101
          END IF
          
          
          ITER = ITER + 1
          IF (ITER .GT. CHECKSTEP) ITER = 1
          
          NITER = NITER + 1
      ENDDO
            
101   CONTINUE   
      END SUBROUTINE
      
      END MODULE DAY3
      MODULE DAY5
    
      USE THERMCALC
      
      PUBLIC TEST_HELM, TEST_FUGCOEF, TEST_DPRESS, TEST_THERMO
      CONTAINS
      
      SUBROUTINE TEST_HELM()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      PARAMETER (NCA = 7)
      INTEGER MODEL, LIST(NCA)
      
      DOUBLE PRECISION T, P
      DOUBLE PRECISION FEED(NCA), FUG(NCA), ZFACTOR
      DOUBLE PRECISION FUGL(NCA), FUGH(NCA), TFEED(NCA)
      DOUBLE PRECISION ZFACTORL, ZFACTORH
      DOUBLE PRECISION EPSILON, DERROR, F1, F2, F
      
      EPSILON = 1E-5
      DERROR = 1E-8
      
      
      LIST(:NCA) = (/1,2,3,5,7,8,13/)
      
      P = 5
      T = 200
      FEED(:NCA) = (/94.30,2.70,0.74,0.49,0.27,0.10,1.40/)
      
      
      WRITE(*,*)
      WRITE(*,*) "TESTING FIRST ORDER COMPOSITIONAL", 
     &            " DERIVATIVES OF HELMHOLZ FUNCTION"
      WRITE(*,*)
      
      DO MODEL = 1,3
          WRITE(*,*)
          WRITE(*,*) "TESTING MODEL NUMBER ", MODEL
          WRITE(*,*)
          
          CALL INIT(NCA, MODEL, LIST)
          
C         CHECK FUGACITY CEFFICIENTS
          CALL FUGAC(T, P, FEED, ZFACTOR, FUG)
          
          DO J = 1,NCA
              DO I = 1,NCA
                  TFEED(I) = FEED(I)
              ENDDO
              
              TFEED(J) = TFEED(J) - EPSILON
              CALL FUGACCOEF(NCA, T, P, TFEED, F1)
              
              TFEED(J) = TFEED(J) + 2 * EPSILON
              CALL FUGACCOEF(NCA, T, P, TFEED, F2)        
              
              F = (F2 - F1) / (2 * EPSILON)
              
              IF (ABS (F - FUG(J)) .GE. DERROR) THEN
                  WRITE(*,*)
                  WRITE(*,*) "TEST FOR COMP ", J, " FAILED" 
                  WRITE(*,*) "EXP ", FUG(J), ", GOT ", F
                  WRITE(*,*)
              END IF              
          ENDDO          
          
      ENDDO      
      
               
      END SUBROUTINE
      
      SUBROUTINE TEST_FUGCOEF()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      PARAMETER (NCA = 7)
      INTEGER MODEL, LIST(NCA)
      
      DOUBLE PRECISION T, P
      DOUBLE PRECISION FEED(NCA), FUG(NCA), ZFACTOR
      DOUBLE PRECISION FUGL(NCA), FUGH(NCA), TFEED(NCA)
      DOUBLE PRECISION ZFACTORL, ZFACTORH
      DOUBLE PRECISION EPSILON, DERROR, F
      
      EPSILON = 1E-5
      DERROR = 1E-8
      
      
      LIST(:NCA) = (/1,2,3,5,7,8,13/)
      
      P = 5
      T = 200
      FEED(:NCA) = (/94.30,2.70,0.74,0.49,0.27,0.10,1.40/)
      
      
      WRITE(*,*)
      WRITE(*,*) "TESTING FIRST ORDER COMPOSITIONAL", 
     &            " DERIVATIVES OF FUGACITY"
      WRITE(*,*)
      
      DO MODEL = 1,3
          WRITE(*,*)
          WRITE(*,*) "TESTING MODEL NUMBER ", MODEL
          WRITE(*,*)
          
          CALL INIT(NCA, MODEL, LIST)
          
C         CHECK FUGACITY CEFFICIENTS
          CALL FUGAC(T, P, FEED, ZFACTOR, FUG)
          
          DO J = 1,NCA
              DO I = 1,NCA
                  TFEED(I) = FEED(I)
              ENDDO
              
              TFEED(J) = TFEED(J) - EPSILON
              CALL FUGAC(T, P, TFEED, ZFACTORL, FUGL)
              
              TFEED(J) = TFEED(J) + 2 * EPSILON
              CALL FUGAC(T, P, TFEED, ZFACTORH, FUGH)      
              
              F = 0
              DO I = 1,NCA
                  F = F + FEED(I) * (FUGH(I) - FUGL(I)) / (2 * EPSILON)
              ENDDO
              
              IF (ABS (F - 0) .GE. DERROR) THEN
                  WRITE(*,*)
                  WRITE(*,*) "TEST FOR COMP ", J, " FAILED" 
                  WRITE(*,*) "EXP ", 0, ", GOT ", F
                  WRITE(*,*)
              END IF              
          ENDDO          
          
      ENDDO      
      
               
      END SUBROUTINE
      
      SUBROUTINE TEST_DPRESS()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      PARAMETER (NCA = 7)
      INTEGER MODEL, LIST(NCA)
      
      DOUBLE PRECISION T, P
      DOUBLE PRECISION FEED(NCA), FUG(NCA), ZFACTOR
      DOUBLE PRECISION FUGL(NCA), FUGH(NCA), TFEED(NCA)
      DOUBLE PRECISION ZFACTORL, ZFACTORH
      DOUBLE PRECISION EPSILON, DERROR, F, FEXP
      
      EPSILON = 1E-5
      DERROR = 1E-8
      
      
      LIST(:NCA) = (/1,2,3,5,7,8,13/)
      
      P = 5
      T = 200
      FEED(:NCA) = (/94.30,2.70,0.74,0.49,0.27,0.10,1.40/)
      
      
      WRITE(*,*)
      WRITE(*,*) "TESTING FIRST ORDER PRESSURE", 
     &            " DERIVATIVE OF FUGACITY"
      WRITE(*,*)
      
      DO MODEL = 1,3
          WRITE(*,*)
          WRITE(*,*) "TESTING MODEL NUMBER ", MODEL
          WRITE(*,*)
          
          CALL INIT(NCA, MODEL, LIST)
          
          CALL FUGAC(T, P, FEED, ZFACTOR, FUG)
          CALL FUGAC(T, P - EPSILON, FEED, ZFACTORL, FUGL)         
          CALL FUGAC(T, P + EPSILON, FEED, ZFACTORH, FUGH)      
              
          F = 0
          DO I = 1,NCA
              F = F + FEED(I) * (FUGH(I) - FUGL(I)) / (2 * EPSILON)
          ENDDO
              
          FEXP = (ZFACTOR - 1) * SUM(FEED) / P
              
          IF (ABS (F - FEXP) .GE. DERROR) THEN
              WRITE(*,*)
              WRITE(*,*) "EXP ", FEXP, ", GOT ", F
              WRITE(*,*)
          END IF          
      ENDDO      
               
      END SUBROUTINE
      
      SUBROUTINE TEST_THERMO()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      PARAMETER (NCA = 7)
      INTEGER MODEL, LIST(NCA)
      
      DOUBLE PRECISION T, P
      DOUBLE PRECISION FEED(NCA), FUG(NCA), ZFACTOR
      DOUBLE PRECISION FUGT(NCA), FUGP(NCA), FUGX(NCA, NCA), AUX(6)
      
      DOUBLE PRECISION TFEED(NCA)
      DOUBLE PRECISION EPSILON, DERROR, F, F1, F2
      
      EPSILON = 1E-5
      DERROR = 1E-8
      
      
      LIST(:NCA) = (/1,2,3,5,7,8,13/)
      
      P = 5
      T = 200
      FEED(:NCA) = (/94.30,2.70,0.74,0.49,0.27,0.10,1.40/)
      
      
      WRITE(*,*)
      WRITE(*,*) "TESTING THERMO"
      WRITE(*,*)

      DO MODEL = 1,4
          IF (MODEL .EQ. 3) GOTO 100
          
          WRITE(*,*)
          WRITE(*,*) "TESTING MODEL NUMBER ", MODEL
          WRITE(*,*)
          
          CALL INDATA(NCA, MODEL, LIST)
          
C         CHECK FUGACITY CEFFICIENTS
          CALL THERMO(T, P, FEED, FUG,
     &     FUGT = FUGT, FUGP = FUGP, FUGX = FUGX, AUX = AUX)
          
          DO J = 1,NCA
              DO I = 1,NCA
                  TFEED(I) = FEED(I)
              ENDDO
              
              TFEED(J) = TFEED(J) - EPSILON
              CALL THERMOCOEF(NCA, T, P, TFEED, F1)
                            
              TFEED(J) = TFEED(J) + 2 * EPSILON
              CALL THERMOCOEF(NCA, T, P, TFEED, F2)
              
C             HELMHOLZ
              F = (F2 - F1) / (2 * EPSILON)              
              IF (ABS (F - FUG(J)) .GE. DERROR) THEN
                  WRITE(*,*)
                  WRITE(*,*) "HELM TEST FOR COMP ", J, " FAILED" 
                  WRITE(*,*) "EXP ", FUG(J), ", GOT ", F
                  WRITE(*,*)
              END IF              
              
C             COMPOSITIONAL DERIVATIVE
              F = 0
              DO I = 1,NCA
                  F = F + FEED(I) * FUGX(I, J) * SUM(FEED)
              ENDDO
              IF (ABS (F - 0) .GE. DERROR) THEN
                  WRITE(*,*)
                  WRITE(*,*) "COMP DER TEST FOR COMP ", J, " FAILED" 
                  WRITE(*,*) "EXP ", 0, ", GOT ", F
                  WRITE(*,*)
              END IF
          ENDDO
          
C         CHECK PRESSURE DERIVATIVE
          F = 0
          DO I = 1,NCA
              F = F + FEED(I) * FUGP(I)
          ENDDO              
          FEXP = (AUX(1) - 1) * SUM(FEED) / P              
          IF (ABS (F - FEXP) .GE. DERROR) THEN
              WRITE(*,*)
              WRITE(*,*) "PRESSURE DERIVATIVE FAILED"
              WRITE(*,*) "EXP ", FEXP, ", GOT ", F
              WRITE(*,*)
          END IF
          
          
C         CHECK TEMPERATURE DERIVATIVE
          F = 0
          DO I = 1,NCA
              F = F + FEED(I) * FUGT(I) / SUM(FEED)
          ENDDO              
          FEXP = - (AUX(3) / T)              
          IF (ABS (F - FEXP) .GE. DERROR) THEN
              WRITE(*,*)
              WRITE(*,*) "TEMPERATURE DERIVATIVE FAILED"
              WRITE(*,*) "EXP ", FEXP, ", GOT ", F
              WRITE(*,*)
          END IF
          
100   CONTINUE    
      ENDDO      
      
               
      END SUBROUTINE
      
      SUBROUTINE FUGACCOEF(NCA, T, P, FEED, F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      INTEGER J      
      DOUBLE PRECISION T, P, F
      DOUBLE PRECISION FEED(NCA), Z, FUG(NCA)
      
            
      CALL FUGAC(T, P, FEED, Z, FUG)
           
      F = 0
      DO I = 1,NCA
          F = F + FEED(I) * FUG(I)
      ENDDO
      END SUBROUTINE
      
      SUBROUTINE THERMOCOEF(NCA, T, P, FEED, F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      INTEGER J      
      DOUBLE PRECISION T, P, F
      DOUBLE PRECISION FEED(NCA), Z, FUG(NCA)
                  
      CALL THERMO(T, P, FEED, FUG)
           
      F = 0
      DO I = 1,NCA
          F = F + FEED(I) * FUG(I)
      ENDDO
      END SUBROUTINE
      
      END MODULE DAY5
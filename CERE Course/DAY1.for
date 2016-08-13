      MODULE DAY1
      
      USE PRINTFUNC
      USE THERMCALC
      
      PUBLIC RUN_DAY1
      CONTAINS
      
      SUBROUTINE RUN_DAY1(NCA, T, P, Z, COMPS)
C     NCA      no. of components      (i)
C     Z       feed mole fractions     (i)
C     P, T    PRESSURE, TEMPERATURE   (I)
C     COMPS   LIST OF COMP NAMES      (I)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NCA
      DOUBLE PRECISION T,P
      DOUBLE PRECISION Z(NCA)
      CHARACTER*4 COMPS(NCA)
      CHARACTER FTYPE
      
      INTEGER NITER, IER
      DOUBLE PRECISION BETA
      DOUBLE PRECISION X(NCA), Y(NCA), FUG(NCA)
            
      CALL THERMO(T,P,Z,FUG,FLUID='L',FTYPE=FTYPE)
      
      WRITE(*,*)
      WRITE(*,*) "FLUID TYPE ", FTYPE
      WRITE(*,"(A20, A20)") "COMPPONENT", "FUGACITY"
      DO I = 1,NCA
          WRITE(*,"(A20, F20.10)") COMPS(I), FUG(I)
      ENDDO
      WRITE(*,*) 
      
      END SUBROUTINE
      
      END MODULE DAY1
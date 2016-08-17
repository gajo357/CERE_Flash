      PROGRAM CERECourse
      
      USE PRINTFUNC
      USE THERMCALC
      USE DAY1
      USE DAY2
      USE DAY3
      USE DAY4
      USE DAY6
      USE DAY7

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
C      PARAMETER (NCA=7)
      PARAMETER (NCA=5)
      
      INTEGER ICEQ
      DIMENSION LIST(NCA)
      
      DOUBLE PRECISION T,P
      DOUBLE PRECISION Z(NCA)
      CHARACTER CONT
      CHARACTER*4 COMPS(NCA)
      
      IER = 0
      
      ICEQ =  0
      
C      LIST(:NCA) = (/1,2,3,5,7,8,13/)
      LIST(:NCA) = (/1,2,3,14,15/)
                       
      CALL INDATA(NCA,ICEQ,LIST)
      
      DO I = 1,NCA
          CALL GETNAME(I, COMPS(I))
      END DO
      
C      Z(:NCA) = (/94.30,2.70,0.74,0.49,0.27,0.10,1.40/)
      Z(:NCA) = (/0.66, 0.03, 0.01, 0.05, 0.25/)
      
100   CONTINUE
      
      WRITE(*,*)
      WRITE(*,*) "DO YOU WANT TO CONTINUE Y/N"
      READ(*,*) CONT
      IF(CONT.EQ.'N') GOTO 200
      WRITE(*,*) "ENTER T,P"
      READ(*,*) T,P
      WRITE(*,*)
      
C      CALL RUN_DAY1(NCA, T, P, Z, COMPS)
      
C      CALL RUN_DAY2(NCA, T, P, Z, COMPS)
      
C      CALL RUN_DAY3(NCA, T, P, Z, COMPS)
      
C      CALL RUN_DAY4(NCA, T, P, Z) 
      
C      CALL RUN_DAY6(NCA, T, P, Z, COMPS)
      
      CALL RUN_DAY7(NCA, T, P, Z, LIST, COMPS)      
      
      GOTO 100
      
200   CONTINUE
      
      READ(*,*)
      
      end program CERECourse
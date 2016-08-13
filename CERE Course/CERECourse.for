      PROGRAM CERECourse
      
      USE PRINTFUNC
      USE THERMCALC
      USE DAY1
      USE DAY2
      USE DAY3
      USE DAY4

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      PARAMETER (NCA=7)
      
      INTEGER ICEQ
      DIMENSION LIST(NCA)
      
      DOUBLE PRECISION T,P
      DOUBLE PRECISION Z(NCA)
      CHARACTER CONT
      CHARACTER*4 COMPS(NCA)
      
      IER = 0
      
      ICEQ =  0
      
      LIST(:NCA) = (/1,2,3,5,7,8,13/)
                       
      CALL INDATA(NCA,ICEQ,LIST)
      
      DO I = 1,NCA
          CALL GETNAME(I, COMPS(I))
      END DO
      
      Z(:NCA) = (/94.30,2.70,0.74,0.49,0.27,0.10,1.40/)
      
100   CONTINUE
      
      WRITE(*,*)
      WRITE(*,*) "DO YOU WANT TO CONTINUE Y/N"
      READ(*,*) CONT
      IF(CONT.EQ.'N') GOTO 200
      WRITE(*,*) "ENTER T,P"
      READ(*,*) T,P
      WRITE(*,*)
      
      CALL RUN_DAY1(NCA, T, P, Z, COMPS)
      
      CALL RUN_DAY2(NCA, T, P, Z, COMPS)
      
      CALL RUN_DAY3(NCA, T, P, Z, COMPS)
      
      CALL RUN_DAY4(NCA, T, P, Z)      
      
      GOTO 100
      
200   CONTINUE
      
      READ(*,*)
      
      end program CERECourse
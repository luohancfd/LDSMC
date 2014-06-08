      SUBROUTINE TOUCH_BOUND(N,K,JJ,XI,YI,X,Y,X1,Y1,X2,Y2,DTIM,CELL0,STAT)
   !---N---code of molecule
   !---K---cell of the molecule
   !---JJ--num of boundary 
      
      USE CELLINFO
      USE NODEINFO
      USE MOLECS
      USE CALC
      USE GAS
      
      INTEGER :: N,K,JJ,L,I,CELL0,DRIFT,STAT
      REAL*8 :: X,Y,XI,YI,X1,Y1,X2,Y2,DTIM,XINEW,YINEW
      INTEGER :: DIS_DETECT
      
      I=ICELL(3+JJ,K)
      IF( I .EQ. 1 .OR. I .EQ. 2)THEN 
         STAT=1
         DTIM=0.0D00
      ELSE
         CALL INTERSECT_SOL(XI,YI,X,Y,X1,Y1,X2,Y2,XINEW,YINEW)
         DTIM=DTIM-(XINEW-XI)/PV(1,N)
         IF(I .EQ. 3)THEN
            CALL SYM_MOL(N,X,Y,X1,Y1,X2,Y2)
         ELSE IF( I .EQ. 4)THEN
            CALL REFLECT_2D(N,K,JJ)
            X=XINEW+PV(1,N)*DTIM
            Y=YINEW+PV(2,N)*DTIM
         END IF
!加一定的漂移防止分子落到网格边界这种尴尬的情况
         DRIFT=-1
         XI=XINEW+DTIM*(10.0D00**DRIFT)*PV(1,N)
         YI=YINEW+DTIM*(2.0D00**DRIFT)*PV(2,N)
         DO WHILE(DIS_DETECT(XI,YI,K) .EQ. 0)
            DRIFT=DRIFT-1
            XI=XINEW+DTIM*(10.0D00**DRIFT)*PV(1,N)
            YI=YINEW+DTIM*(10.0D00**DRIFT)*PV(2,N)
         END DO
         DTIM=DTIM-DTIM*(10.0D00**DRIFT)
         CELL0=K
      END IF
      
      END SUBROUTINE
      FUNCTION DIS_DETECT(X,Y,I)
    !this function is used to detect whether point(x,y) is in the sampling cell I
      USE CELLINFO
      USE NODEINFO
      IMPLICIT NONE
      
      INTEGER :: I,J,K,DIS_DETECT
      REAL*8 :: X,Y,VEC(2,3),A(3),B(3)
      LOGICAL D(3)
      REAL*8,EXTERNAL :: VEC_PRODUCT
      INTEGER,EXTERNAL :: ADD_ONE
      
      !VEC is the vector from (X,Y) to the node of the cell
      !B is the length of the VEC
      DO J=1,3
         VEC(1,J)=NODES(1,ICELL(J,I))-X
         VEC(2,J)=NODES(2,ICELL(J,I))-Y
         B(J)=DSQRT(VEC(1,J)**2.0D00+VEC(2,J)**2.0D00)
         !nominal VEC for numeric precision
         VEC(1,J)=VEC(1,J)/B(J)
         VEC(2,J)=VEC(2,J)/B(J)
      END DO
         
      DO J=1,3
         A(J)=VEC_PRODUCT(VEC(1,J),VEC(2,J),VEC(1,ADD_ONE(J)),VEC(2,ADD_ONE(J)))
      END DO
      
      IF( (A(1).LT. 0.0D00) .AND. (A(2).LT. 0.0D00) .AND. (A(3).LT. 0.0D00))THEN
         DIS_DETECT=1
      ELSE
         DIS_DETECT=0
      END IF
      
      RETURN
      END FUNCTION
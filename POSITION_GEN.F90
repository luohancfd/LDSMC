      SUBROUTINE POSITION_GEN(I,J,X,Y)
!----I is the code of collision cell in the sampling cell
!----J is the code of sampling cell
!----X is the coordinate X of the molecule
!----Y is the coordinate Y of the molecule
      USE CELLINFO
      USE NODEINFO
      IMPLICIT NONE
      INTEGER :: I,J
      REAL*8 :: X,Y,RANF1,RANF2,x1,x2,x3,y1,y2,y3,TEMP
      REAL*8,EXTERNAL :: VEC_PRODUCT

      x1=NODES(1,ICELL(1,J))
      x2=NODES(1,ICELL(2,J))
      x3=NODES(1,ICELL(3,J))
      y1=NODES(2,ICELL(1,J))
      y2=NODES(2,ICELL(2,J))
      y3=NODES(2,ICELL(3,J))
      
      TEMP=0
!--When temp=0,the molecule lies on the edge of a coolision cell
      DO WHILE (TEMP .EQ. 0)
         CALL RANDOM_NUMBER(RANF1)
         CALL RANDOM_NUMBER(RANF2)
      !--When RANF1&RANF2 =0,the molecule lies on the node of a collision cell
         DO WHILE(RANF1 .EQ. 0 .OR. RANF1 .EQ. 1)
            CALL RANDOM_NUMBER(RANF1)
         END DO
         DO WHILE(RANF2 .EQ. 0 .OR. RANF2 .EQ. 1)
            CALL RANDOM_NUMBER(RANF2)
         END DO
         X=0.5*RANF1*(x2-x1)+0.5*RANF2*(x3-x2)
         Y=0.5*RANF1*(y2-y1)+0.5*RANF2*(y3-y2)
         TEMP=VEC_PRODUCT(x3-x1,y3-y1,X,Y)
      END DO
      
      SELECT CASE(I)
      CASE(1)
         IF (TEMP .LT. 0.0D00)THEN
            X=0.5*(x1+x3)-X
            Y=0.5*(y1+y3)-Y
         ELSE IF(TEMP .GT. 0.0D00)THEN
            X=x1+X
            Y=y1+Y
         ELSE
            WRITE(*,*) "Initiate error!!!"
         END IF
      CASE(2)
         IF (TEMP .LT. 0.0D00)THEN
            X=0.5*(x2+x3)-X
            Y=0.5*(y2+y3)-Y
         ELSE IF(TEMP .GT. 0.0D00)THEN
            X=0.5*(x1+x2)+X
            Y=0.5*(y1+y2)+Y
         ELSE
            WRITE(*,*) "Initiate error!!!"
         END IF
      CASE(3)
         IF (TEMP .LT. 0.0D00)THEN
            X=x3-X
            Y=y3-Y
         ELSE IF(TEMP .GT. 0.0D00)THEN
            X=0.5*(x1+x3)+X
            Y=0.5*(y1+y3)+Y
         ELSE
            WRITE(*,*) "Initiate error!!!"
         END IF
      CASE(4)
         IF (TEMP .LT. 0.0D00)THEN
            X=0.5*(x1+x2)+X
            Y=0.5*(y1+y2)+Y
         ELSE IF(TEMP .GT. 0.0D00)THEN
            X=0.5*(x2+x3)-X
            Y=0.5*(y2+y3)-Y
         ELSE
            WRITE(*,*) "Initiate error!!!"
         END IF
      END SELECT
      END SUBROUTINE
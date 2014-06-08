      LOGICAL FUNCTION INTERSECT_DETECT(X1,Y1,X2,Y2,X3,Y3,X4,Y4)
   !---This function is used to detect whether two segments intersect
   !!!!---this function return a LOGICAL variable--!!!!
      IMPLICIT NONE
      REAL*8 :: X1,Y1,X2,Y2,X3,Y3,X4,Y4
      !REAL*8,EXTERNAL :: VEC_PRODUCT这个也不用了，影响计算效率
      !LOGICAL,EXTERNAL :: ONEDGE_DETECT这个函数不用了，已经改进算法
      
      REAL*8 :: D1,D2,D3,D4
      INTERSECT_DETECT=.false.
      IF( (MIN(X1,X2) .LE. MAX(X3,X4)) .AND. (MIN(X3,X4) .LE. MAX(X1,X2))&  !快速排斥实验
         & .AND. (MIN(Y1,Y2) .LE. MAX(Y3,Y4)) .AND. (MIN(Y3,Y4) .LE. MAX(Y1,Y2)))THEN
         
         !D1=VEC_PRODUCT(X3-X1,Y3-Y1,X2-X1,Y2-Y1)
         D1=(X3-X1)*(Y2-Y1)-(X2-X1)*(Y3-Y1)
         !D2=VEC_PRODUCT(X2-X1,Y2-Y1,X4-X1,Y4-Y1)
         D2=(X2-X1)*(Y4-Y1)-(X4-X1)*(Y2-Y1)
         !D3=VEC_PRODUCT(X1-X3,Y1-Y3,X4-X3,Y4-Y3)
         D3=(X1-X3)*(Y4-Y3)-(X4-X3)*(Y1-Y3)
         !D4=VEC_PRODUCT(X4-X3,Y4-Y3,X2-X3,Y2-Y3)
         D4=(X4-X3)*(Y2-Y3)-(X2-X3)*(Y4-Y3)
      
         IF( (D1*D2 .GE. 0.0D00) .AND. (D3*D4 .GE. 0) )THEN!跨立实验
            INTERSECT_DETECT=.true.
            RETURN
         ELSE
            INTERSECT_DETECT=.false.
            RETURN
         END IF
      ELSE
          INTERSECT_DETECT=.false. 
          RETURN
      END IF
      
      END FUNCTION INTERSECT_DETECT
            
            
            
            
   !        
   !   LOGICAL FUNCTION ONEDGE_DETECT(X1,Y1,X2,Y2,X3,Y3)
   !   IMPLICIT NONE
   !!---this function is used to be compled with function above to detect whether a point
   !!---(X1,Y1) is on the line composed by (X2,Y2),(X3,Y3)
   !!!!!-----this function return a LOGICAL variable--!!!!
   !   REAL*8 :: X1,X2,X3,Y1,Y2,Y3
   !   IF(  ((MIN(X2,X3) .LE. X1) .AND. (X1 .LE. MAX(X2,X3)))  &
   !      &.AND. ((MIN(Y2,Y3) .LE. Y1) .AND. (Y1 .LE. MAX(Y2,Y3))))THEN
   !      ONEDGE_DETECT=.true.
   !   ELSE
   !      ONEDGE_DETECT=.false.
   !   END IF
   !   RETURN
   !   END FUNCTION ONEDGE_DETECT
   !         
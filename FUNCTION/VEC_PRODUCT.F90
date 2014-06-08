      FUNCTION VEC_PRODUCT(X1,Y1,X2,Y2)
   !---This function is used to compute the vectoral product of 2D vectors, which means z is zero for both vectors
   !!!!---this function return a REAL*8 variable--!!!!
   !--(X1,Y1).product.(X2,Y2)
      IMPLICIT NONE
      REAL*8 :: X1,X2,Y1,Y2,VEC_PRODUCT
      
      VEC_PRODUCT=X1*Y2-X2*Y1
      END FUNCTION VEC_PRODUCT
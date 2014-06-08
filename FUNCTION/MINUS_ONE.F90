      FUNCTION MINUS_ONE(I)
      IMPLICIT NONE
      
      INTEGER :: I,MINUS_ONE
      
      SELECT CASE(I)
      CASE(1) 
         MINUS_ONE=3
      CASE DEFAULT  !I=1,2
         MINUS_ONE=I-1
      END SELECT
      END FUNCTION
      
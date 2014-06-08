      FUNCTION ADD_ONE(I)
      IMPLICIT NONE
      
      INTEGER :: I,ADD_ONE
      
      SELECT CASE(I)
      CASE(3) 
         ADD_ONE=1
      CASE DEFAULT  !I=1,2
         ADD_ONE=I+1
      END SELECT
      END FUNCTION
      
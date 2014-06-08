SUBROUTINE LIN2_SOL(A,B,SOL)
   IMPLICIT NONE
   
   REAL*8 :: A(2,2),B(2,1),SOL(2,1)
   
   SOL(1,1)=(B(1,1)*A(2,2)-B(2,1)*A(1,2))/(A(1,1)*A(2,2)-A(2,1)*A(1,2))
   SOL(2,1)=(B(1,1)*A(2,1)-B(2,1)*A(1,1))/(A(1,2)*A(2,1)-A(2,2)*A(1,1))
   
   RETURN
   END SUBROUTINE

   
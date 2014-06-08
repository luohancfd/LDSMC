      SUBROUTINE INTERSECT_SOL(X1,Y1,X2,Y2,X3,Y3,X4,Y4,X,Y)
      !INCLUDE 'link_fnl_shared.h' !change to [use imsl] for CVF6.6
      !use lin_sol_gen_int
      IMPLICIT NONE
      !use to solve the common point on edge (x1,y1---x2,y2) and edge (x3,y3---x4,y4)
      REAL*8 :: X1,Y1,X2,Y2,X3,Y3,X4,Y4,X,Y
      REAL*8 :: A(2,2),B(2,1),SOL(2,1)
      
      A(1,:)=(/ Y3-Y4, X4-X3 /)
      A(2,:)=(/ Y1-Y2 ,X2-X1 /)
      B(1,1)=X4*Y3-X3*Y4
      B(2,1)=X2*Y1-X1*Y2
      
      !CALL D_LIN_SOL_GEN(A,B,SOL)
      CALL LIN2_SOL(A,B,SOL)

      X=SOL(1,1)
      Y=SOL(2,1)
      
      END SUBROUTINE
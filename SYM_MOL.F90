      SUBROUTINE SYM_MOL(N,X,Y,X1,Y1,X2,Y2)
      !INCLUDE 'link_fnl_shared.h' !change to [use imsl] for CVF6.6
      !use lin_sol_gen_int

      USE CELLINFO
      USE MOLECS
      USE NODEINFO
      
      INTEGER :: N
      REAL*8 :: X,Y,X1,X2,Y1,Y2
      REAL*8 :: A(2,2),B(2,1),SOL(2,1)  !solution for symetry boundary
      
      
     !solve the symmetrical velocity
      A(1,:)=(/ Y2-Y1, X1-X2 /)
      A(2,:)=(/ X1-X2, Y1-Y2 /)
      B(1,1)=PV(1,N)*(Y1-Y2)-PV(2,N)*(X1-X2)
      B(2,1)=PV(1,N)*(X1-X2)+PV(2,N)*(Y1-Y2)
      !CALL D_LIN_SOL_GEN(A,B,SOL)
      CALL LIN2_SOL(A,B,SOL)
      PV(1,N)=SOL(1,1)
      PV(2,N)=SOL(2,1)

      
      !solve the symmetrical final position
      A(1,:)=(/ Y2-Y1, X1-X2 /)
      A(2,:)=(/ X1-X2, Y1-Y2 /)
      B(1,1)=2*X1*Y2-2*X2*Y1+X*(Y1-Y2)-Y*(X1-X2)
      B(2,1)=X*(X1-X2)+Y*(Y1-Y2)
      !CALL D_LIN_SOL_GEN(A,B,SOL)
      CALL LIN2_SOL(A,B,SOL)
      X=SOL(1,1)
      Y=SOL(2,1)
      
      
      END SUBROUTINE
   
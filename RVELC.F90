      SUBROUTINE RVELC(U,V,VMP)
      !
      USE CONST
      !
      IMPLICIT NONE
      !
      !--generates two random velocity components U and V in an equilibrium
      !--gas with most probable speed VMP
      REAL(KIND=8) :: U,V,VMP,A,B,RANF
      !
      CALL RANDOM_NUMBER(RANF)
      A=DSQRT(-DLOG(RANF))
      CALL RANDOM_NUMBER(RANF)
      B=DPI*RANF
      U=A*DSIN(B)*VMP
      V=A*DCOS(B)*VMP
      RETURN
      !
      END SUBROUTINE RVELC
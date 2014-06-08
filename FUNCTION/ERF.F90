FUNCTION ERF(S)
!
!--evaluates the error function of S
!
   IMPLICIT NONE
   REAL*8 :: ERF,S,B,C,D,T
   B=ABS(S)
   IF (B < 4.) THEN
      C=EXP(-B*B)
      T=1./(1.+0.3275911*B)
      D=1.-(0.254829592*T-0.284496736*T*T+1.421413741*T*T*T- &
         1.453152027*T*T*T*T+1.061405429*T*T*T*T*T)*C
   ELSE
      D=1.
   END IF
   IF (S < 0.) D=-D
   ERF=D
   RETURN
END FUNCTION ERF
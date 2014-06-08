FUNCTION GAM(X)
!
!--calculates the Gamma function of X.
!
   IMPLICIT NONE
   REAL*8 :: A,Y,X,GAM
   A=1.
   Y=X
   IF (Y < 1.) THEN
      A=A/Y
   ELSE
      Y=Y-1.
      DO WHILE (Y >= 1.)
         A=A*Y
         Y=Y-1.
      END DO
   END IF
   GAM=A*(1.-0.5748646*Y+0.9512363*Y**2-0.6998588*Y**3+  &
      0.4245549*Y**4-0.1010678*Y**5)
   !
   RETURN
   !
END FUNCTION GAM
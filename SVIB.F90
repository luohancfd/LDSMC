      SUBROUTINE SVIB(L,TEMP,IVIB,K)
      !
      !--sets a typical vibrational state at temp. TEMP of mode K of species L
      !
      USE GAS
      !
      IMPLICIT NONE
      !
      INTEGER :: K,L,N
      REAL(KIND=8) :: TEMP,RANF
      INTEGER :: IVIB
      !
      CALL RANDOM_NUMBER(RANF)
      N=-DLOG(RANF)*TEMP/SPVM(1,K,L)                 !eqn(11.24)
      !--the state is truncated to an integer
      IVIB=N
      END
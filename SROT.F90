      SUBROUTINE SROT(L,TEMP,ROTE)
      !
      !--sets a typical rotational energy ROTE of species L
      !
      USE GAS
      USE CONST
      !
      IMPLICIT NONE
      !
      INTEGER :: I,L
      REAL(KIND=8) :: A,B,ROTE,ERM,TEMP,RANF
      IF (ISPR(1,L).EQ.2) THEN
   !--rotational degree=2, distribution function:-1/(kT)*exp[-rote/(kT)]
        CALL RANDOM_NUMBER(RANF)
        ROTE=-DLOG(RANF)*BOLTZ*TEMP
      ELSE
        A=0.5D00*ISPR(1,L)-1.D00
        I=0
        DO WHILE (I == 0)
          CALL RANDOM_NUMBER(RANF)
          ERM=RANF*10.D00
      !--rotational degree range from 0~10
      !--there is an energy cut-off at 10 kT
          B=((ERM/A)**A)*DEXP(A-ERM)
          CALL RANDOM_NUMBER(RANF)
          IF (B > RANF) I=1
        END DO
        ROTE=ERM*BOLTZ*TEMP
      END IF
      !
      RETURN
      !
      END SUBROUTINE SROT
      SUBROUTINE COLLISION
      !
      !--calculate the collisions
      !

      USE MOLECS
      USE CELLINFO
      USE GAS
      USE TRANSC
      USE SAMPLES
      USE CALC
      IMPLICIT NONE

      INTEGER :: I
      INTEGER :: J,K,NSEL,M1,M2,MPOS,MT
      REAL*8 :: DTC,ASEL,SEPSMIN,SEPS,RANF
      INTEGER :: ISM
      !--I,J,K,L,M,N---loop number
      !--II,JJ,KK------loop control
      !---ISM----------Select method---0 for transient sub cell 1 for virtual sub cell
      !---M1-----------selected molecule 1
      !---M2-----------selected molecule 2
      !---MT-----------temporal selected molecule
      !---MPOS---------selected molecule position in ICREF
      !---SEPS---------square of current collision pair separation
      !---SEPSMIN------square of current minimum molecule separation
      
      !IPCP=0 this is done in molecules move

      ALLOCATE(IPCP(NM))
      IPCP=0

      !$OMP PARALLEL DEFAULT(SHARED) 
      !$OMP DO SCHEDULE(STATIC) PRIVATE(I,J,K,DTC,ISM,ASEL,NSEL,MPOS,M1,M2,SEPS,SEPSMIN,MT,RANF) REDUCTION(+:TOTCOL,TCOL,COLLS,CLSEP)
      DO I=1,NCCELLS
         IF (FTIME-CCELL(5,I) > CCELL(3,I)) THEN
            DTC=2.0D00*CCELL(3,I)

            IF(ICCELL(2,I) .GT. 1) CCELL(5,I)=CCELL(5,I)+DTC

            IF( ICCELL(2,I) .GT. 30)THEN
               ISM=0!transient sub cell
            ELSE
               ISM=1!virtual sub cell
            END IF
            !NTC method
            ASEL=0.5D00*ICCELL(2,I)*(ICCELL(2,I)-1)*&
               &FNUM*CCELL(4,I)*DTC/CCELL(1,I)+CCELL(2,I)
            IF(ASEL .LT. 0) ASEL=0.0D00
            NSEL=ASEL
            CCELL(2,I)=ASEL-NSEL

            IF( (NSEL .GT. 0) .AND. (ICCELL(2,I) .GT. 1))THEN
               IF(ISM .EQ. 0)THEN
                  CALL TRANSIENT_SUBCELL_GEN(I)
                  DO J=1,NSEL
                     !select first molecule
                     CALL RANDOM_NUMBER(RANF)
                     MPOS=INT(RANF*DFLOAT(ICCELL(2,I))-0.0001)+ICCELL(1,I)+1
                     M1=ICREF(MPOS)
                     !Select second molecule
                     CALL TRANSIENT_SUBCELL_SELECT(M1,M2)
                     !Collision
                     CALL MOLECULES_COL(M1,M2)
                  END DO
                  CALL TRANSIENT_SUBCELL_DEL
               ELSE
                  DO J=1,NSEL
                     !select first molecule
                     CALL RANDOM_NUMBER(RANF)
                     MPOS=INT(RANF*DFLOAT(ICCELL(2,I))-0.0001)+ICCELL(1,I)+1
                     M1=ICREF(MPOS)
                     !select second molecule
                     SEPSMIN=(MAXCELL(2)*0.5D00)**2.0D00
                     DO K=1,ICCELL(2,I)
                        MPOS=ICCELL(1,I)+K
                        MT=ICREF(MPOS)
                        IF((MT .NE. M1) .AND. (MT .NE. IPCP(M1)))THEN
                           SEPS=(PX(MT)-PX(M1))**2.0D00+(PY(MT)-PY(M1))**2.0D00
                           IF( SEPS .LT. SEPSMIN)THEN
                              M2=MT
                              SEPSMIN=SEPS
                           END IF
                        END IF
                     END DO
                     !collision
                     CALL MOLECULES_COL(M1,M2)
                  END DO
               END IF
            END IF   !end nsel possible
         END IF
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
      DEALLOCATE(IPCP)
      END SUBROUTINE
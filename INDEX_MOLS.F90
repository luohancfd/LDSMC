      SUBROUTINE INDEX_MOLS
      !
      !--index the molecules to the collision cells
      !
      USE MOLECS
      USE CALC
      USE CELLINFO
      !
      IMPLICIT NONE
      !
      INTEGER :: N,M,K
      !
      !--N,M,K working integer
      !
      ICCELL(2,:)=0
      !
      IF (NM .NE. 0) THEN
         DO N=1,NM
            M=IPCCELL(N)
            ICCELL(2,M)=ICCELL(2,M)+1
         END DO
         !
         M=0
         DO N=1,NCCELLS
            ICCELL(1,N)=M
            M=M+ICCELL(2,N)
            ICCELL(2,N)=0
         END DO
         !
         DO N=1,NM
            M=IPCCELL(N)
            ICCELL(2,M)=ICCELL(2,M)+1
            K=ICCELL(1,M)+ICCELL(2,M)
            ICREF(K)=N
         END DO
         !
      END IF
      !
          
      RETURN
      !
      END SUBROUTINE INDEX_MOLS
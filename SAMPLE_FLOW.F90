SUBROUTINE SAMPLE_FLOW
!
!--sample the flow properties
!
USE GAS
USE MOLECS
USE SAMPLES
USE CONST
USE CALC
USE CELLINFO
!
IMPLICIT NONE
!
INTEGER :: NC,NCC,LS,N,M,K,L,I,KV
REAL(KIND=8) :: A,TE,TT,WF
!
!--NC the sampling cell number
!--NCC the collision cell number
!--LS the species code
!--N,M,K working integers
!--TE total translational energy
!
NSAMP=NSAMP+1
WRITE (*,*) 'Sample',NSAMP
WRITE (9,*) NM,'Mols. at sample',NSAMP
!

DO N=1,NM
  NCC=IPCCELL(N) 
  NC=ICCELL(3,NCC)
  EMOLS(NC)=EMOLS(NC)+1.0D00
  LS=IPSP(N)
  CS(1,NC,LS)=CS(1,NC,LS)+1.0D00     
  DO M=1,3
     CS(M+1,NC,LS)=CS(M+1,NC,LS)+PV(M,N)
     CS(M+4,NC,LS)=CS(M+4,NC,LS)+PV(M,N)**2.0D00
  END DO
  IF (MMRM > 0) CS(8,NC,LS)=CS(8,NC,LS)+PROT(N)
  IF (MMVM > 0) THEN
      IF (ISPV(LS).GT.0) THEN
        DO K=1,ISPV(LS)
          CS(K+8,NC,LS)=CS(K+8,NC,LS)+DFLOAT(IPVIB(K,N))*BOLTZ*SPVM(1,K,LS)
        END DO
      END IF
  END IF
END DO

!
!
RETURN
!
END SUBROUTINE SAMPLE_FLOW
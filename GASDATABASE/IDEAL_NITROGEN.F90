SUBROUTINE IDEAL_NITROGEN
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
MSP=1
MMRM=1
MMVM=0
MNRE=0
MTBP=0
MNSR=0
MEX=0
MMEX=0
!
CALL ALLOCATE_GAS
!
SP(1,1)=4.17D-10
SP(2,1)=273.
SP(3,1)=0.74
SP(4,1)=1.
SP(5,1)=4.65D-26
ISPR(1,1)=2
ISPR(2,1)=0
SPR(1,1)=5.
FSP(1)=1.0D00

!
RETURN
END SUBROUTINE IDEAL_NITROGEN
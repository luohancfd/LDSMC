      SUBROUTINE READ_DATA
      USE GAS
      USE CALC
      USE BOUNDINFO
      IMPLICIT NONE
      
      INTEGER :: N,I
      REAL*8 :: TEMP
      OPEN (3,FILE='DS1VD.TXT')
      OPEN(10,FILE="INIT.TXT",ACTION='READ')
   !The gas type
      !GASCODE=3
      READ(10,*)
      READ(10,*)
      READ(10,*) GASCODE
      
      SELECT CASE (GASCODE)
      CASE(1)   
         WRITE (3,*) 'Hard sphere gas'
         CALL HARD_SPHERE
      CASE(2) 
         WRITE (3,*) 'Argon'
         CALL ARGON
      CASE(3)   
          WRITE (3,*) 'Nitrogen'
          CALL IDEAL_NITROGEN
      CASE(5) 
          WRITE (3,*) 'Ideal air'
          CALL IDEAL_AIR
      CASE(4)
         WRITE (3,*) 'Reak Nitrogen'
         CALL REAL_NITROGEN
      CASE DEFAULT
          WRITE(*,*) "UNSUPPORTTED GAS"
      END SELECT
   !--collision control
      NNC=1
     !IF (NNC == 0) WRITE (*,*) ' Collision partners are selected randomly from the collision cell  '
     !IF (NNC == 1) WRITE (*,*) ' Nearest neighbor collisions are employed '
   !--control of time step   
      IMTS=1
      IF (IMTS == 0) WRITE (3,*) ' The move time step is uniform over the cells  '
      IF (IMTS == 1) WRITE (3,*) ' The move time step can vary over the cells '
   !--property of fluid region   
      !FND=3.0E19
      READ(10,*)
      READ(10,*)FND
      WRITE (3,*) '    The stream number density is',FND
      
      !FTMP=100.0D00
      READ(10,*)
      READ(10,*)FTMP
      WRITE (3,*) '    The stream temperature is',FTMP
      
      !FVTMP=100.0D00
      READ(10,*)
      READ(10,*)FVTMP
      WRITE (3,*) '    The stream vibrational temperature is',FVTMP
   !-------------------------following is important, they are about boundary property--------
   !--Inflow-------------------------
      !VFX=1000.0D00
      READ(10,*)
      READ(10,*) VFX
      WRITE (3,*) '    The stream velocity in the x direction is',VFX    

      !VFY=0.0D00
      READ(10,*)
      READ(10,*) VFY
      WRITE (3,*) '    The stream velocity in the y direction is',VFY  
   !---solid bound, which collided with molecule--------------------
      !---following data get from Padilla's graduate thesis
      !---ASSESSMENT OF GAS-SURFACE INTERACTION MODELS FOR COMPUTATION OF RAREFIED HYPERSONIC FLOWS
      !---Page 142, simulation of Apollo6 reentry
      !---这里设置为不依赖于分子种类，实际可以改变
      !---另外本程序认为所有固壁面都一样，以后也可以更改
      ALLOCATE(SSSEG_INT(9,MSP))
      SSSEG_INT(1,:)=0.0D00
      SSSEG_INT(2,:)=0.0D00
      SSSEG_INT(3,:)=0.0D00
      READ(10,*)
      READ(10,*) TEMP
      SSSEG_INT(4,:)=TEMP
      READ(10,*)
      READ(10,*) TEMP
      SSSEG_INT(5,:)=TEMP
      READ(10,*)
      READ(10,*) TEMP
      SSSEG_INT(6,:)=TEMP
      READ(10,*)
      READ(10,*) TEMP
      SSSEG_INT(7,:)=TEMP
      READ(10,*)
      READ(10,*) TEMP
      SSSEG_INT(8,:)=TEMP
      READ(10,*)
      READ(10,*) TEMP
      SSSEG_INT(9,:)=TEMP
      READ(10,*)
      READ(10,*) NWALL
      ALLOCATE(SSEG_INT(5,NWALL))
      READ(10,*)
      DO I=1,NWALL
         READ(10,*) SSEG_INT(1,I)
      END DO
      SSEG_INT(2,:)=0.0D00
      SSEG_INT(3,:)=0.0D00
      SSEG_INT(4,:)=0.0D00
      SSEG_INT(5,:)=0.0D00
      !time control
      CPDTM=0.2D00
      
      TPDTM=0.5D00 !这个大小的设置还不确定
      
      
      !MOLSC=8  !number of simulated molecules in the smallest collision cell
      READ(10,*)
      READ(10,*) MOLSC
      MOLSC=4*MOLSC!number of simulated molecules in the smallest samping cell
      
      !sample control
      SAMPRAT=10
      WRITE (3,*) ' The number if time steps in a sampling interval is ',SAMPRAT
      OUTRAT=10
      WRITE (3,*) ' The number of smpling intervals in an output interval is ',OUTRAT
      
      ISF=0
      IF (ISF == 0) WRITE (3,*) ' The sampling is for an eventual steady flow '
      IF (ISF == 1) WRITE (3,*) ' The sampling is for a continuing unsteady flow '
      CLOSE(3)
      CLOSE(10)
      ICENS=1
   END SUBROUTINE
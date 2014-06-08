      SUBROUTINE OUTPUT_RESULT
      USE CELLINFO
      USE NODEINFO
      USE SAMPLES
      USE CALC
      USE MOLECS
      USE CONST
      USE OUTPUT
      USE BOUNDINFO
      USE GAS
      IMPLICIT NONE
      
      INTEGER :: OCELL,OBOUND,OWALL,OINTERVAL,OEDGE,I,J,K,L
      REAL*8 :: COOR_CELL(2,3),A,UU,DOF,SVDF,VDF,VDOFM,TVIBM,EVIBM,VEL,SUM(12),SUMS(7,2)
      REAL*8 :: OAREA,OPAREA,VG,TS,S,X,Y,TTG,TRG
      REAL*8,ALLOCATABLE :: TV(:,:),VDOF(:),TVIB(:)
      CHARACTER*20 :: FILENAME
      
      INTEGER,EXTERNAL :: ADD_ONE
!----A------------real variable
!----COOR_CELL----coordinate of the sampling cell
!----DOF----------degree of freedom
!----EVIBM--------vibrational energy of mixture       
!----OCELL--------sampling cell code
!----OBOUND-------code of BC
!----OWALL--------code of wall( means solid surface)
!----OINTERVAL----code of interval
!----OEDGE--------code of edge
!----OAREA--------area of the edge
!----SUM(1)-------sum of molecule number
!----SVDF---------the sum of the degrees of freedom of species L
!----TV(K,L)------the temperature of vibrational mode K of species L
!----UU-----------sum of the square of velocity
!----VDOFM--------vibrational degrees of freedom of mixture                    
!----VG-----------slip velocity
!----VEL----------flow speed
!----TTG----------the translational temperature just above the surface      
!----TRG----------the rotational energy just above the surface              
!----TS-----------temperature of surface
!----S------------dis to start point




      
      
      !calculate flow-field properties
      ALLOCATE(VAR(23,NCELLS),VARSP(11,NCELLS,MSP),&
         &TV(MMVM,MSP),DF(NCELLS,MMVM,MSP),VDOF(MSP),TVIB(MSP))
      VAR=0.0D00;DF=0.0D00; VARSP=0.0D00
      
      DO OCELL=1,NCELLS
         !get coordinate of sampling cell
         DO I=1,3
            DO J=1,2
               COOR_CELL(J,I)=NODES(J,ICELL(I,OCELL))
            END DO
         END DO
         !initiate summary variable
         SUM=0.0D00
         A=FNUM/(NSAMP*CELL(1,OCELL))
         !set the coordinate of the center of sampling cell
         DO I=1,2   !I=1 for x, I=2 for y
            DO J=1,3
               VAR(I,OCELL)=VAR(I,OCELL)+COOR_CELL(I,J)
            END DO
            VAR(I,OCELL)=VAR(I,OCELL)/3.0D00
         END DO
         DO I=1,MSP
            SUM(1)=SUM(1)+CS(1,OCELL,I)  !sum of molecule number
            SUM(2)=SUM(2)+CS(1,OCELL,I)*SP(5,I) !sum of mass
            DO J=1,3
               SUM(2+J)=SUM(2+J)+CS(1+J,OCELL,I)*SP(5,I)! sum of kinetic
               IF (CS(1,OCELL,I) > 0.1) THEN
                  VARSP(J+1,OCELL,I)=CS(J+4,OCELL,I)/CS(1,OCELL,I)
                  !--VARSP(2,3,4 are temporarily the mean of the squares of the velocities
                  VARSP(J+8,OCELL,I)=CS(J+1,OCELL,I)/CS(1,OCELL,I)
                  !--VARSP(9,10,11 are temporarily the mean of the velocities
               END IF
            END DO
            SUM(10)=SUM(10)+SP(5,I)*CS(5,OCELL,I)  !sum of mass*vx^2
            SUM(11)=SUM(11)+SP(5,I)*CS(6,OCELL,I)  !sum of mass*vy^2
            SUM(12)=SUM(12)+SP(5,I)*CS(7,OCELL,I)  !sum of mass*vz^2
            IF (CS(1,OCELL,I) > 0.5) THEN
               SUM(7)=SUM(7)+CS(5,OCELL,I)+CS(6,OCELL,I)+CS(7,OCELL,I)
            END IF

            IF (ISPR(1,I) > 0) THEN
               SUM(8)=SUM(8)+CS(8,OCELL,I)         !sum of rotational energy
               SUM(9)=SUM(9)+CS(1,OCELL,I)*ISPR(1,I)! sum of rotational degree
            END IF
         END DO
         SUM(6)=SUM(10)+SUM(11)+SUM(12)!sum of 2*translational energy
         
         DO I=1,MSP
            VARSP(1,OCELL,I)=0.0D00
            VARSP(6,OCELL,I)=0.0D00
            VARSP(7,OCELL,I)=0.0D00
            VARSP(8,OCELL,I)=0.0D00
            VARSP(5,OCELL,I)=0.0D00
            IF (SUM(1) > 1.0D-6) THEN
               VARSP(1,OCELL,I)=CS(1,OCELL,I)/SUM(1)
               IF ((ISPR(1,I) > 0).AND.(CS(1,OCELL,I) > 0.5)) THEN
                  VARSP(6,OCELL,I)=(2.0D00/BOLTZ)*CS(8,OCELL,I)/(ISPR(1,I)*CS(1,OCELL,I))
               END IF
            END IF
            DO J=1,3
               !--VARSP(2,3,4 are temporarily the mean of the squares of the velocities
               !--VARSP(9,10,11 are temporarily the mean of overall velocities
               VARSP(J+1,OCELL,I)=(SP(5,I)/BOLTZ)*(VARSP(J+1,OCELL,I)-VARSP(J+8,OCELL,I)**2)
               VARSP(5,OCELL,I)=VARSP(5,OCELL,I)+VARSP(J+1,OCELL,I)
            END DO
            VARSP(5,OCELL,I)=VARSP(5,OCELL,I)/3.0D00
         END DO

         IF( SUM(1) .GT. 0.0D00)THEN
            VAR(3,OCELL)=SUM(1)*A  !number density
            VAR(4,OCELL)=SUM(2)*A  !density
            VAR(5,OCELL)=SUM(3)/SUM(2) !x velocity
            VAR(6,OCELL)=SUM(4)/SUM(2) !y velocity
            VAR(7,OCELL)=SUM(5)/SUM(2) !z velocity
            UU=VAR(5,OCELL)**2.0D00+VAR(6,OCELL)**2.0D00+VAR(7,OCELL)**2.0D00
            VAR(8,OCELL)=DABS(SUM(6)-SUM(2)*UU)/SUM(1)/(BOLTZ*3.0D00) !translational temperature
            VAR(20,OCELL)=(ABS(SUM(10)-SUM(2)*VAR(5,OCELL)**2.0D00))/(BOLTZ*SUM(1))! x translational temperature
            VAR(21,OCELL)=(ABS(SUM(11)-SUM(2)*VAR(6,OCELL)**2.0D00))/(BOLTZ*SUM(1))! y translational temperature
            VAR(22,OCELL)=(ABS(SUM(12)-SUM(2)*VAR(7,OCELL)**2.0D00))/(BOLTZ*SUM(1))! z translational temperature
            VAR(19,OCELL)=VAR(3,OCELL)*BOLTZ*VAR(8,OCELL)            !--scalar pressure (now (from V3) based on the translational temperature)
            
            !rotational variable
            IF (SUM(9) > 0.0D00) THEN
               VAR(9,OCELL)=(2.0D00/BOLTZ)*SUM(8)/SUM(9)    !--rotational temperature
            ELSE
               VAR(9,OCELL)=0.0D00
            END IF
            VAR(10,OCELL)=0.0D00  !vibration default
            DOF=3.0D00+SUM(9)/SUM(1) !degree of freedom=rot+trans
            
            !temperature based on translational energy and rotational energy
            VAR(11,OCELL)=(3.0D00*VAR(8,OCELL)+(SUM(9)/SUM(1))*VAR(9,OCELL))/DOF
            
            !vibrational variable
            IF (MMVM > 0) THEN
               DO I=1,MSP
                  SVDF=0.0D00   !--SVDF is the sum of the degrees of freedom of species I
                  VDOF(I)=0.0D00 !sum of i specie vib DOF
                  IF (ISPV(I).GT.0) THEN
                     DO J=1,ISPV(I)
                        IF (CS(J+8,OCELL,I) < BOLTZ) THEN !vibrational energy ~= 0
                           TV(J,I)=0.0D00
                           DF(OCELL,J,I)=0.0d00
                        ELSE
                           TV(J,I)=SPVM(1,J,I)/DLOG(1.+BOLTZ*SPVM(1,J,I)*CS(1,OCELL,I)/CS(J+8,OCELL,I)) !11.32
                           DF(OCELL,J,I)=2.0D00*(SPVM(1,J,I)/TV(J,I))/(DEXP(SPVM(1,J,I)/TV(J,I))-1.0D00)
                        END IF
                        VDOF(I)=VDOF(I)+DF(OCELL,J,I)
                        SVDF=SVDF+DF(OCELL,J,I)
                     END DO
                     TVIB(I)=0.0D00
                     DO J=1,ISPV(I)
                        IF (SVDF > 1.0E-6) THEN
                           TVIB(I)=TVIB(I)+TV(J,I)*DF(OCELL,J,I)/SVDF
                        ELSE
                           TVIB(I)=VAR(8,OCELL)   !SVDF~=0, DF,TV~=0 in case of 0/0=NaN
                        END IF
                     END DO
                  ELSE
                     TVIB(I)=0.0D00  
                     VDOF(I)=0.0D00
                  END IF
                  VARSP(7,OCELL,I)=TVIB(I)
                  VARSP(8,OCELL,I)=VDOF(I)  !temporary assignment
               END DO
               VDOFM=0.0D00
               EVIBM=0.0D00
               DO I=1,MSP
                  VDOFM=VDOFM+VDOF(I)*CS(1,OCELL,I)/SUM(1)
                  EVIBM=EVIBM+0.5*BOLTZ*VDOF(I)*TVIB(I)*CS(1,OCELL,I)/SUM(1)
               END DO
               IF( VDOFM .LT. 1.0E-6)THEN !VDOF~=0
                  TVIBM=VAR(8,OCELL)
               ELSE
                  TVIBM=(2.0D00/BOLTZ)*EVIBM/VDOFM
               END IF
               VAR(10,OCELL)=TVIBM
               DOF=3.0D00+SUM(9)/SUM(1)+VDOFM!--DOF is the number of degrees of freedom               
               VAR(11,OCELL)=(3.0D00*VAR(8,OCELL)+(SUM(9)/SUM(1))*VAR(9,OCELL)+VDOFM*TVIBM)/DOF
               !--the overall temperature now includes vibration, see eqn (11.12)
            END IF
            
            DO I=1,MSP
               !--set the overall temperature for the individual species
               VARSP(8,OCELL,I)=(3.0D00*VARSP(5,OCELL,I)+ISPR(1,I)*VARSP(6,OCELL,I)+     &
                  VARSP(8,OCELL,I)*VARSP(7,OCELL,I))/(3.0D00+ISPR(1,I)+VARSP(8,OCELL,I))
               !--convert the species velocity components to diffusion velocities
               DO J=1,3
                  VARSP(J+8,OCELL,I)=VARSP(J+8,OCELL,I)-VAR(J+4,OCELL)
               END DO
            END DO
            
            VEL=DSQRT(VAR(5,OCELL)**2.0D00+VAR(6,OCELL)**2.0D00+VAR(7,OCELL)**2.0D00)
            VAR(12,OCELL)=VEL/DSQRT((DOF+2.0D00)/DOF*BOLTZ/(SUM(2)/SUM(1))*VAR(11,OCELL)) !Mach number
            VAR(13,OCELL)=SUM(1)/NSAMP
            IF (COLLS(OCELL).GT.2.) THEN
               VAR(14,OCELL)=(FTIME-TISAMP)/(COLLS(OCELL)*2.0D00/(EMOLS(OCELL)/NSAMP))
               !--mean collision time
               !这个与DS2一样但是公式略有变形，做如下解释（我自己的解释，没有找到对应公式）
               !FTIME-TISAMP summary of time 总时间
               !EMOLC(OCELL)/NSAMP 每个sampling cell中平均的分子数
               !COLLS(OCELL)*2/（EMOLC(OCELL)/NSAMP）平均每个分子所存在的碰撞对的数目
               !如现有A\B\C分子，A分子，和B、C分别发生过一次碰撞，则A所存在的碰撞对数目为2,B\C分别为1
               !继而VAR(14,OCELL)表示平均而言每个分子每两次相邻碰撞的间隙即碰撞时间
               VAR(15,OCELL)=0.92132*SQRT(ABS(SUM(7)/SUM(1)-UU))*VAR(14,OCELL)
               !--mean free path (based on r.m.s speed with correction factor for equilib.)
               !--0.92132修正参数不知从何而来
            ELSE
               VAR(14,OCELL)=1.E10
               VAR(15,OCELL)=1.41E18/VAR(3,OCELL)
               !--m.f.p set by nominal values
            END IF
            IF (ICENS == 1) VAR(23,OCELL)=SQRT(0.33333*(((VAR(20,OCELL)/VAR(8,OCELL))-1.)**2.0D00+       &
               ((VAR(21,OCELL)/VAR(8,OCELL))-1.)**2.0D00+((VAR(22,OCELL)/VAR(8,OCELL))-1.)**2.0D00))
            IF (ICENS == 2) VAR(23,OCELL)=VAR(20,OCELL)/VAR(21,OCELL)-1.
            IF (ICENS == 3) VAR(23,OCELL)=VAR(9,OCELL)/VAR(8,OCELL)
            IF (ICENS == 4) VAR(23,OCELL)=VAR(10,OCELL)/VAR(8,OCELL)

            IF (COLLS(OCELL) > 0.0D00) THEN
               VAR(16,OCELL)=(CLSEP(OCELL)/COLLS(OCELL))/VAR(15,OCELL)
            ELSE
               VAR(16,OCELL)=0.0D00
            END IF
         ELSE
            DO I=3,19
               VAR(I,OCELL)=0.0D00
            END DO
         END IF
         VAR(17,OCELL)=VEL
         IF (SUM(1) > 0.0D00) THEN
            VAR(18,OCELL)=ATAN2(VAR(6,OCELL),VAR(5,OCELL))*180.0D00/PI
         END IF
      END DO
      
      CALL POST
      !calculate surface property
      ALLOCATE(VARS(30,MAXVAL(IWALL(0,:)),NWALL))
      DO OWALL=1,NWALL
         WRITE(FILENAME,*) OWALL
         OPEN(10,FILE=".\result\LDSMC_WALL"//TRIM(ADJUSTL(FILENAME))//".dat")
101      FORMAT ('VARIABLES = "S","X","Y","Number Flux","Vertical_P","Paralle_P",&
                  "Z_P","Total Incident E","Total Reflect E","Net Heat","Temperature","VG","TTG","TRG"')
         WRITE(10,101)
         S=0.0D00
         OPAREA=0.0D00
         DO OINTERVAL=1,IWALL(0,OWALL)
            SUMS=0.0D00
            OEDGE=IWALL(OINTERVAL,OWALL)
            OCELL=IBEDGE4(1,OEDGE)
            K=IBEDGE4(2,OEDGE)
            L=ADD_ONE(K)
            COOR_CELL(1,K)=NODES(1,ICELL(K,OCELL))
            COOR_CELL(2,K)=NODES(2,ICELL(K,OCELL))
            COOR_CELL(1,L)=NODES(1,ICELL(L,OCELL))
            COOR_CELL(2,L)=NODES(2,ICELL(L,OCELL))
            OAREA=(COOR_CELL(1,K)-COOR_CELL(1,L))**2.0D00+&
               &(COOR_CELL(2,K)-COOR_CELL(2,L))**2.0D00
            OAREA=DSQRT(OAREA)            
            S=S+0.5D00*OPAREA+0.5D00*OAREA
            OPAREA=OAREA            
            
            X=0.5D00*(COOR_CELL(1,K)+COOR_CELL(1,L))
            Y=0.5D00*(COOR_CELL(2,K)+COOR_CELL(2,L))
            
            A=FNUM/(OAREA*(FTIME-TISAMP))
            IF(MMVM .GT. 0)THEN
               L=7
            ELSE IF(MMRM .GT. 0)THEN
               L=6
            ELSE
               L=5
            END IF
            
            DO I=1,MSP
               IF(ISPV(I) .GT. 0 )THEN
                  L=7
               ELSE IF(ISPR(1,I) .GT. 0)THEN
                  L=6
               ELSE
                  L=5
               END IF
               DO J=1,L
                  DO K=1,2
                     SUMS(J,K)=SUMS(J,K)+CSS(J,OEDGE,I,K)
                  END DO
               END DO
            END DO
            
            TS=SSEG(1,OEDGE)
            IF (CSSS(1,OEDGE) > 0.5) THEN
               VG=CSSS(3,OEDGE)/CSSS(2,OEDGE)
               TTG=(CSSS(4,OEDGE)-CSSS(2,OEDGE)*VG*VG)/(CSSS(1,OEDGE)*3.0D00*BOLTZ)-TS
               IF (CSSS(6,OEDGE) > 0.0D00) THEN
                  TRG=(2.0D00/BOLTZ)*(CSSS(5,OEDGE)/CSSS(6,OEDGE))-TS
               ELSE
                  TRG=0
               END IF
            ELSE
               VG=0.0D00
               TTG=0.0D00
               TRG=0.0D00
            END IF
            VARS(1,OINTERVAL,OWALL)=SUMS(1,1)
            VARS(2,OINTERVAL,OWALL)=SUMS(1,2)
            VARS(3,OINTERVAL,OWALL)=SUMS(1,1)*A
            VARS(4,OINTERVAL,OWALL)=SUMS(1,2)*A
            VARS(5,OINTERVAL,OWALL)=SUMS(2,1)*A
            VARS(6,OINTERVAL,OWALL)=SUMS(2,2)*A
  !---shear stress>0, the direction of shear stress is the same as the direction of edge
  !------------------(the direction of edge is the clockwise route of sampling cell)
            VARS(7,OINTERVAL,OWALL)=SUMS(3,1)*A
            VARS(8,OINTERVAL,OWALL)=SUMS(3,2)*A
            VARS(9,OINTERVAL,OWALL)=SUMS(4,1)*A
            VARS(10,OINTERVAL,OWALL)=SUMS(4,2)*A
            VARS(11,OINTERVAL,OWALL)=SUMS(5,1)*A
            VARS(26,OINTERVAL,OWALL)=SUMS(5,1)*A
            VARS(12,OINTERVAL,OWALL)=SUMS(5,2)*A
            VARS(27,OINTERVAL,OWALL)=SUMS(5,2)*A
            IF( MMRM .GT. 0)THEN
               VARS(13,OINTERVAL,OWALL)=SUMS(6,1)*A
               VARS(14,OINTERVAL,OWALL)=SUMS(6,2)*A
               VARS(24,OINTERVAL,OWALL)=(SUMS(6,1)-SUMS(6,2))*A
               VARS(26,OINTERVAL,OWALL)=VARS(26,OINTERVAL,OWALL)+SUMS(6,1)*A
               VARS(27,OINTERVAL,OWALL)=VARS(27,OINTERVAL,OWALL)+SUMS(6,2)*A
            END IF
            IF( MMVM .GT. 0)THEN
               VARS(15,OINTERVAL,OWALL)=SUMS(7,1)*A
               VARS(16,OINTERVAL,OWALL)=SUMS(7,2)*A
               VARS(25,OINTERVAL,OWALL)=(SUMS(7,1)-SUMS(7,2))*A
               VARS(26,OINTERVAL,OWALL)=VARS(26,OINTERVAL,OWALL)+SUMS(7,1)*A
               VARS(27,OINTERVAL,OWALL)=VARS(27,OINTERVAL,OWALL)+SUMS(7,2)*A
            END IF
            !下面三个我也不明确，主要是对CSSS不理解
            VARS(17,OINTERVAL,OWALL)=VG
            VARS(18,OINTERVAL,OWALL)=TTG
            VARS(19,OINTERVAL,OWALL)=TRG
            VARS(20,OINTERVAL,OWALL)=(SUMS(2,1)+SUMS(2,2))*A                       
            VARS(21,OINTERVAL,OWALL)=(SUMS(3,1)+SUMS(3,2))*A
            VARS(22,OINTERVAL,OWALL)=(SUMS(4,1)+SUMS(4,2))*A
            VARS(23,OINTERVAL,OWALL)=(SUMS(5,1)-SUMS(5,2))*A
            VARS(28,OINTERVAL,OWALL)=VARS(26,OINTERVAL,OWALL)-VARS(27,OINTERVAL,OWALL)
            VARS(29,OINTERVAL,OWALL)=TS
            WRITE (10,"(13E14.5)") S,X,Y,VARS(3,OINTERVAL,OWALL),VARS(20,OINTERVAL,OWALL),VARS(21,OINTERVAL,OWALL),&
               &VARS(22,OINTERVAL,OWALL),VARS(26,OINTERVAL,OWALL),&
               &VARS(27,OINTERVAL,OWALL),VARS(28,OINTERVAL,OWALL),&
               &TS,VG,TTG,TRG
         END DO
         CLOSE(10)
      END DO
      
     
      DEALLOCATE(VAR,VARSP,TV,DF,VDOF,TVIB,VARS)      
            
      

      
      
      
      
      
      END SUBROUTINE
         
         


            
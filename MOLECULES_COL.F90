      SUBROUTINE MOLECULES_COL(M1,M2)
      USE MOLECS
      USE CONST
      USE CALC
      USE CELLINFO
      USE GAS
      USE SAMPLES
      IMPLICIT NONE
      
      INTEGER :: M1,M2,I,J,K,II,KV,N,OCCELL,OCELL,NSP,MAXLEV,IV
      REAL*8 :: A,B,C,D,OC,SD,VRR,VR,CVR,SEP,ECT,EVIB,ERM,ECC,COLT,&
         &PROB,VCM(3),VRCP(3),VRC(3),RM1,RM2,ZV,RANF
      INTEGER :: LS,MS,KS,JS,NLOOP
!---I,J,K--------loop number
!---II-----------loop control
!---OCCELL-------collision cell of molecule M1 and M2
!---OCELL--------sampling cell of molecule M1 and M2
!---VRC----------components of the relative velocity
!---VRCP---------components of the post collision relative velocity
!---VRR----------square of relative velocity
!---VR-----------relative velocity
!---VCM----------velocity of center
!---CVR----------cross-section*relative velocity
!---SEP----------collision pair separation
!---ECT----------relative translational energy
!---EVIB---------vibrational energy
!---ERM----------rotational energy
!---ECC----------energy to be divided
!---MAXLEV-------maximum vibrational level
!---COLT---------collision temperature=temp+vib temp
!---IV-----------vibrational level
!---LS,MS--------specie of M1 and M2
!---RM1,RMM2-----molecule mass parameters for M1 and M2
!---NLOOP--------NLOOP


      OCCELL=IPCCELL(M1)
      OCELL=ICCELL(3,OCCELL)
      DO I=1,3
         VRC(I)=PV(I,M1)-PV(I,M2)
      END DO
      VRR=VRC(1)**2.0D00+VRC(2)**2.0D00+VRC(3)**2.0D00
      
      IF( VRR .GT. 1.0E-6)THEN
         !different molecule
         VR=DSQRT(VRR)
         IF ( MSP .EQ. 1)THEN
            CVR=VR*CXSS*((2.0D00*BOLTZ*SP(2,1)/(RMAS*VRR))**(SP(3,1)-0.5D00))*RGFS
            IF(CVR .GT. CCELL(4,OCCELL))  CCELL(4,IPCCELL(M1))=CVR
            
            CALL RANDOM_NUMBER(RANF)
            IF( RANF .LT. CVR/CCELL(4,OCCELL))THEN
               !collision occur
               TOTCOL=TOTCOL+1.0D00
               TCOL(1,1)=TCOL(1,1)+2.D00
               COLLS(OCELL)=COLLS(OCELL)+1.D000
               SEP=DSQRT((PX(M1)-PX(M2))**2.0D00+(PY(M1)-PY(M2))**2.0D00)
               CLSEP(OCELL)=CLSEP(OCELL)+SEP   
               
               IF ((ISPR(1,1) > 0)) THEN
                  ECT=0.5D00*RMAS*VRR
                  DO NSP=1,2
                     !--consider the molecules in turn
                     IF (NSP .EQ. 1) THEN
                        K=M1
                     ELSE
                        K=M2
                     END IF
                     
                     !--deal with vibrational energy
                     IF (MMVM > 0) THEN
                        IF (ISPV(1) > 0) THEN
                           DO KV=1,ISPV(1)
                              EVIB=DFLOAT(IPVIB(KV,K))*BOLTZ*SPVM(1,KV,1)
                              ECC=ECT+EVIB
                              IF (SPVM(3,KV,1) > 0.) THEN
                                 !-- quantizing of COLT 
                                 MAXLEV=ECC/(BOLTZ*SPVM(1,KV,1))
                                 COLT=DFLOAT(MAXLEV)*SPVM(1,KV,1)/(3.5-SPM(3,1,1))  !(5.42)
                                 B=SPVM(4,KV,1)/SPVM(3,KV,1)    !Tdiss/Tref          !
                                 A=SPVM(4,KV,1)/COLT       !Tdiss/T
                                 ZV=(A**SPM(3,1,1))*(SPVM(3,KV,1)*(B**(-SPM(3,1,1))))**&
                                    &(((A**0.3333333D00)-1.D00)/((B**0.33333D00)-1.D00))
                              ELSE
                                 ZV=SPVM(2,KV,1)
                              END IF
                              !
                              CALL RANDOM_NUMBER(RANF)
                              IF (1.0D00/ZV > RANF) THEN
                                 !vibrational energy excited
                                 II=0
                                 DO WHILE (II == 0)
                                    CALL RANDOM_NUMBER(RANF)
                                    IV=RANF*(MAXLEV+0.99999D00)
                                    IPVIB(KV,K)=IV
                                    EVIB=IV*BOLTZ*SPVM(1,KV,1)
                                    IF (EVIB < ECC) THEN
                                       PROB=(1.0D00-EVIB/ECC)**(1.5D00-SPM(3,KV,1))
                                       !--PROB is the probability ratio of eqn (5.61)
                                       CALL RANDOM_NUMBER(RANF)
                                       IF (PROB > RANF) II=1
                                    END IF
                                 END DO
                                 ECT=ECC-EVIB
                              END IF
                           END DO
                        END IF
                     END IF

                     !--deal with rotational energy
                     IF (ISPR(1,1).GT.0) THEN
                        IF (ISPR(2,1) == 0) THEN
                           B=1.0D00/SPR(1,1)
                        ELSE
                           COLT=ECC/((2.5-SP(3,1))*BOLTZ)
                           B=1.0D00/(SPR(1,1)+SPR(2,1)*COLT+SPR(3,1)*COLT*COLT)
                        END IF
                        CALL RANDOM_NUMBER(RANF)
                        IF (B > RANF) THEN
                           ECC=ECT+PROT(K)
                           IF (ISPR(1,1) .EQ. 2) THEN
                              CALL RANDOM_NUMBER(RANF)
                              ERM=1.0D00-RANF**(1.0D00/(2.5D00-SP(3,1)))  !eqn(5.46)
                           ELSE
                              CALL LBS(0.5D00*DFLOAT(ISPR(1,1))-1.D00,1.5D00-SPM(3,1,1),ERM)
                           END IF
                           PROT(K)=ERM*ECC
                           ECT=ECC-PROT(K)
                        END IF
                     END IF
                  END DO
                  VR=DSQRT(2.D00*ECT/SPM(1,1,1))
               END IF
         !----end of L-B model
               
               DO I=1,3
                  VCM(I)=0.5D00*(PV(I,M1)+PV(I,M2))
               END DO
               
               IF(DABS(SP(4,1)-1.0D00) .LT. 1.0E-3)THEN
         !----SP(4,1)=1---use VHS model
                  CALL RANDOM_NUMBER(RANF)
                  B=2.0D00*RANF-1.0D00
         !----B is the cosine of a random elevation angle
                  A=DSQRT(1.0D00-B*B)
                  VRCP(1)=B*VR
                  CALL RANDOM_NUMBER(RANF)
                  C=2.D00*PI*RANF
         !----C is a random azimuth angle
                  VRCP(2)=A*DCOS(C)*VR
                  VRCP(3)=A*DSIN(C)*VR
               ELSE
                  CALL RANDOM_NUMBER(RANF)
                  B=2.D00*(RANF**SP(4,1))-1.D00
         !--B is the cosine of the deflection angle for the VSS model (eqn (11.8)
                  A=SQRT(1.D00-B*B)
                  CALL RANDOM_NUMBER(RANF)
                  C=2.D00*PI*RANF
                  OC=DCOS(C)
                  SD=DSIN(C)
                  D=SQRT(VRC(2)**2+VRC(3)**2)
                  IF (D .GT. 1.0E-6) THEN                                                  
                     VRCP(1)=B*VRC(1)+A*SD*D                                            
                     VRCP(2)=B*VRC(2)+A*(VR*VRC(3)*OC-VRC(1)*VRC(2)*SD)/D               
                     VRCP(3)=B*VRC(3)-A*(VR*VRC(2)*OC+VRC(1)*VRC(3)*SD)/D               
                  ELSE                                                                  
                     VRCP(1)=B*VRC(1)                                                   
                     VRCP(2)=A*OC*VRC(1)                                                
                     VRCP(3)=A*SD*VRC(1)                                                
                  END IF
               END IF
               
               DO I=1,3
                  PV(I,M1)=VCM(I)+0.5D00*VRCP(I)
                  PV(I,M2)=VCM(I)-0.5D00*VRCP(I)
               END DO
               
               IPCP(M1)=M2
               IPCP(M2)=M1
            END IF
         ELSE
            LS=IPSP(M1)
            MS=IPSP(M2)
            CVR=VR*SPM(2,LS,MS)*((2.D00*BOLTZ*SPM(5,LS,MS)/(SPM(1,LS,MS)*VRR))**(SPM(3,LS,MS)-0.5D00))*SPM(6,LS,MS)
            IF (CVR > CCELL(4,OCCELL)) CCELL(4,OCCELL)=CVR
            CALL RANDOM_NUMBER(RANF)
            IF (RANF < CVR/CCELL(4,OCCELL)) THEN
               TOTCOL=TOTCOL+1.D00
               TCOL(LS,MS)=TCOL(LS,MS)+1.D00
               TCOL(MS,LS)=TCOL(MS,LS)+1.D00
               COLLS(OCELL)=COLLS(OCELL)+1.D00
               SEP=DSQRT((PX(M1)-PX(M2))**2.0D00+(PY(M1)-PY(M2))**2.0D00)
               CLSEP(OCELL)=CLSEP(OCELL)+SEP
               RM1=SPM(1,LS,MS)/SP(5,MS)
               RM2=SPM(1,LS,MS)/SP(5,LS)
               DO I=1,3
                  VCM(I)=RM1*PV(I,M1)+RM2*PV(I,M2)
               END DO
         !---未来化学反应的添加点
               IF ((ISPR(1,LS) > 0).OR.(ISPR(1,MS) > 0)) THEN
                  ECT=0.5D00*SPM(1,LS,MS)*VRR
                  DO NSP=1,2
                     IF (NSP == 1) THEN
                        K=M1 ; KS=LS ; JS=MS
                     ELSE
                        K=M2 ; KS=MS ; JS=LS
                     END IF
                     !vibrational energy
                     IF (MMVM > 0) THEN
                        IF (ISPV(KS) > 0) THEN
                           DO KV=1,ISPV(KS)
                              EVIB=DFLOAT(IPVIB(KV,K))*BOLTZ*SPVM(1,KV,KS)
                              ECC=ECT+EVIB
                              MAXLEV=ECC/(BOLTZ*SPVM(1,KV,KS))
                              IF (SPVM(3,KV,KS) > 0.) THEN
                                 !--note quantizing of COLT in the following statements
                                 COLT=(DFLOAT(MAXLEV)*SPVM(1,KV,KS))/(3.5D00-SPM(3,KS,JS))        !--quantized collision temperature
                                 B=SPVM(4,KV,KS)/SPVM(3,KV,KS)    !Tdiss/Tref
                                 A=SPVM(4,KV,KS)/COLT       !Tdiss/T
                                 ZV=(A**SPM(3,KS,JS))*(SPVM(2,KV,KS)*(B**(-SPM(3,KS,JS))))**(((A**0.3333333D00)-1.D00)/((B**0.33333D00)-1.D00))
                              ELSE
                                 ZV=SPVM(2,KV,KS)
                              END IF
                              !解离判断点IF(MAXLEV*SPVM(1,KV,KS) .LT. SPVM(4,KV,SV) .OR. (JCI .EQ. 0))THEN
                              CALL RANDOM_NUMBER(RANF)
                              IF (1.D00/ZV > RANF)THEN
                                 II=0
                                 NLOOP=0
                                 DO WHILE (II .EQ. 0)
                                    NLOOP=NLOOP+1
                                    if (NLOOP > 100) then
                                       write (*,*) "Dealing L-B of molecule",K,'NLOOP=',NLOOP
                                    end if
                                    CALL RANDOM_NUMBER(RANF)
                                    IV=RANF*(MAXLEV+0.99999999D00)
                                    IPVIB(KV,K)=IV
                                    EVIB=DFLOAT(IV)*BOLTZ*SPVM(1,KV,KS)
                                    IF (EVIB < ECC) THEN
                                       PROB=(1.D00-EVIB/ECC)**(1.5D00-SPM(3,KS,JS))
                                       !--PROB is the probability ratio of eqn (5.61)
                                       CALL RANDOM_NUMBER(RANF)
                                       IF (PROB > RANF) II=1
                                    END IF
                                 END DO
                                 ECT=ECC-EVIB
                                 !ELSE
                                 !    !--reflects the infinity of levels beyond the dissociation limit
                                 !    IPV(KV,K)=MAXLEV
                                 !    EVIB=MAXLEV*BOLTZ*SPV(KV,KS)
                                 !    ECT=ECC-EVIB
                                 !END IF
                              END IF
                           END DO
                        END IF
                     END IF
                     
                      !--now rotation of this molecule     
                     IF (ISPR(1,KS).GT.0) THEN
                        IF ((ISPR(2,KS) .EQ. 0).AND.(ISPR(2,JS) .EQ. 0))THEN
                           B=1.0D00/SPM(7,KS,JS)
                        ELSE
                           COLT=ECC/((2.5-SPM(3,KS,JS))*BOLTZ)
                           B=1.0D00/(SPM(7,KS,JS)+SPM(8,KS,JS)*COLT+SPM(9,KS,JS)*COLT*COLT)
                        END IF
                        CALL RANDOM_NUMBER(RANF)
                        IF (B > RANF) THEN
                           ECC=ECT+PROT(K)
                           IF (ISPR(1,KS) .EQ. 2) THEN
                              CALL RANDOM_NUMBER(RANF)
                              ERM=1.0D00-RANF**(1./(2.5-SPM(3,KS,JS)))  !eqn(5.46)
                           ELSE
                              CALL LBS(0.5D00*ISPR(1,KS)-1.0D00,1.5D00-SPM(3,KS,JS),ERM)
                           END IF
                           PROT(K)=ERM*ECC
                           ECT=ECC-PROT(K)
                        END IF
                     END IF   
                  END DO
                  IF (ECT < 0.) THEN
                     WRITE (9,*) 'NEG ECT 1',ECT,' SPECIES',LS,MS
                     ECT=1.E-26
                  END IF
                  VR=DSQRT(2.0D00*ECT/SPM(1,LS,MS))
               END IF
               IF (ABS(SPM(10,LS,MS)-1.0D00) < 0.001) THEN
                  !--use the VHS logic
                  CALL RANDOM_NUMBER(RANF)
                  B=2.0D00*RANF-1.D00
                  !--B is the cosine of a random elevation angle
                  A=DSQRT(1.D00-B*B)
                  VRCP(1)=B*VR
                  CALL RANDOM_NUMBER(RANF)
                  C=2.0D00*PI*RANF
                  !--C is a random azimuth angle
                  VRCP(2)=A*DCOS(C)*VR
                  VRCP(3)=A*DSIN(C)*VR
               ELSE
                  !--use the VSS logic
                  CALL RANDOM_NUMBER(RANF)
                  B=2.D00*(RANF**SP(4,1))-1.D00
                  !--B is the cosine of the deflection angle for the VSS model (eqn (11.8)
                  A=SQRT(1.D00-B*B)
                  CALL RANDOM_NUMBER(RANF)
                  C=2.D00*PI*RANF
                  OC=DCOS(C)
                  SD=DSIN(C)
                  D=SQRT(VRC(2)**2+VRC(3)**2)
                  VRCP(1)=B*VRC(1)+A*SD*D
                  VRCP(2)=B*VRC(2)+A*(VR*VRC(3)*OC-VRC(1)*VRC(2)*SD)/D
                  VRCP(3)=B*VRC(3)-A*(VR*VRC(2)*OC+VRC(1)*VRC(3)*SD)/D
               END IF
               DO I=1,3
                  PV(I,M1)=VCM(I)+0.5D00*VRCP(I)
                  PV(I,M2)=VCM(I)-0.5D00*VRCP(I)
               END DO

               IPCP(M1)=M2
               IPCP(M2)=M1
               !
            END IF  !End the collision procedure
         END IF     !end gas mixture or single gas
      ELSE
         !same molecule
         CCELL(2,OCCELL)=CCELL(2,OCCELL)+1.0D00
      END IF
                  
      END SUBROUTINE
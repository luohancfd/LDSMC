      SUBROUTINE SET_INITIAL_STATE
      USE CELLINFO
      USE NODEINFO
      USE GAS
      USE CONST
      USE CALC
      USE MOLECS
      USE SAMPLES
      IMPLICIT NONE
      
      REAL*8 :: A,B,AA,BB,C,SN,NMI,NRMI,DMOM(3),DENG,ctheta,x1,y1,x2,y2,FVEL,RANF
      REAL*8,ALLOCATABLE :: VB(:,:),ROTE(:)
      INTEGER :: I,J,K,L,II,JJ,KK,NSET,KN,ERROR
      INTEGER*8 :: M,N,O,P
      REAL*8,EXTERNAL :: GAM,ERF
      INTEGER,EXTERNAL :: ADD_ONE
!--ctheta---COSINE the angle between the direction of flow and the vector normal to inflow boundary
!--NRMI--the initial number of real molecules
!--NMI--the initial number of molecules
!--DMOM(N) N=1,2,3 for x,y and z momentum sums of initial molecules
!--VB alternative sets of velocity components
!--NSET the alternative set numbers in the setting of exact initial state
!NSET用于初始化分子的时候生成多组可用数据，择优
!--DENG use to describe the degree of the initial gas diverging from equilibrate state
!--FVEL--flow velocity
      DENG=0.0D00
      NSET=2
      ALLOCATE(VB(3,NSET),ROTE(NSET))
      NRMI=0.0D00
      DMOM(:)=0
      



!---Set the overall parameters
      DO I=1,NCELLS
         NRMI=NRMI+CELL(1,I)
      END DO
      NRMI=NRMI*FND

      NMI=0.0D00
      DO I=1,NCELLS
         NMI=NMI+CEILING(MOLSC/MINCELL(1)*CELL(1,I))
      END DO
      FNUM=NRMI/NMI
      NRMI=CEILING(NRMI,KIND=8)
      NMI=CEILING(NMI,KIND=8)
!---Set the information on the molecular species 
!---Set the information for cells
!---this was finished in preprocessing   


!---Set the information on the molecular species
      A=0.D00
      B=0.D00
      DO L=1,MSP
          A=A+SP(5,L)*FSP(L)
          B=B+(3.+ISPR(1,L))*FSP(L)
          VMP(L)=SQRT(2.D00*BOLTZ*FTMP/SP(5,L))
          IF (L == 1) THEN
             VMPM=VMP(L)
          ELSE
             IF (VMP(L) > VMPM) VMPM=VMP(L)
          END IF
      END DO
      WRITE (*,*) 'Maximum thermal velocity =',VMPM 
      FDEN=A*FND
      FPR=FND*BOLTZ*FTMP
      FMA=VFX/DSQRT((B/(B+2.D00))*BOLTZ*FTMP/A) 
      DO L=1,MSP
         DO M=1,MSP
            SPM(4,L,M)=0.5D00*(SP(1,L)+SP(1,M))
            SPM(3,L,M)=0.5D00*(SP(3,L)+SP(3,M))
            SPM(5,L,M)=0.5D00*(SP(2,L)+SP(2,M))
            SPM(1,L,M)=SP(5,L)*(SP(5,M)/(SP(5,L)+SP(5,M)))
            SPM(2,L,M)=0.25D00*PI*(SP(1,L)+SP(1,M))**2.0D00
            AA=2.5D00-SPM(3,L,M)
            A=GAM(AA)
            SPM(6,L,M)=1.D00/A
            SPM(10,L,M)=0.5D00*(SP(4,L)+SP(4,M))
            SPM(7,L,M)=(SPR(1,L)+SPR(1,M))*0.5D00
            IF ((ISPR(2,L) .EQ. 1).OR.(ISPR(2,M) .EQ. 1)) THEN
               SPM(8,L,M)=(SPR(2,L)+SPR(2,M))*0.5D00
               SPM(9,L,M)=(SPR(3,L)+SPR(3,M))*0.5D00
            END IF
         END DO
      END DO
      IF (MSP == 1) THEN
         RMAS=SPM(1,1,1)
         CXSS=SPM(2,1,1)
         RGFS=SPM(6,1,1)
      END IF
      
      DO L=1,MSP
         CR(L)=0.!collision rate
         DO M=1,MSP
            CR(L)=CR(L)+2.D00*SPI*SPM(4,L,M)**2*FND*FSP(M)*(FTMP/SPM(5,L,M))** &
               (1.-SPM(3,L,M))*DSQRT(2.*BOLTZ*SPM(5,L,M)/SPM(1,L,M))
         END DO
      END DO
      
      A=0.D00
      FP=0.D00
      DO L=1,MSP
         A=A+FSP(L)*CR(L)
         FP=FP+FSP(L)*(2./SPI)*VMP(L)/CR(L)  !平均运动速度乘以碰撞时间
      END DO
      CTM=1.D00/A
      WRITE (*,*) 'Approximate collision time in the stream is',CTM
      WRITE (*,*) 'Approximate mean free path in the stream is',FP 
!----Set initial time step
      DTM=CTM*CPDTM !CPDTM=0.2
      FVEL=DSQRT(VFX**2.0D00+VFY**2.0D00)
      IF (FVEL .GT. VMPM) THEN
         A=FVEL
         IF (FVEL > 3.0D00*VMPM) A=3.0D00*VMPM  
         A=0.5D00*(MINCELL(2)*0.5D00)/A
      ELSE
         A=0.5D00*(MINCELL(2)*0.5D00)/VMPM
      END IF
      IF( DTM >A) DTM=A
      WRITE (*,*) 'The initial value of the overall time step is',DTM
      DTSAMP=SAMPRAT*DTM
      DTOUT=OUTRAT*DTSAMP
      TSAMP=DTSAMP
      TOUT=DTOUT
      
      FTIME=0.D00
      TLIM=1.D20   !流域时间
      TOTMOV=0.D00
      TOTCOL=0.D00
      TCOL=0.D00
      ENTMASS=0.0D00
!--initialise cell quantities associated with collisions

      DO N=1,NCCELLS
         CCELL(3,N)=DTM/2.D00
         CCELL(4,N)=2.D00*VMPM*SPM(2,1,1) !a random value
         CALL RANDOM_NUMBER(RANF)
         CCELL(2,N)=RANF
         CCELL(5,N)=0.D00
      END DO
!--set the entry quantities
!!----ATTENTION---!!Whether a boundary edge is inflow boundary or out flow boundary depends on the theta( defined below),but NOT boundary definition
!!----ATTENTION---!!when theta  < 90 dgree（ctheta>0), inflow boundary, theta >90 (ctheta<0)outflow boundary            
      ALLOCATE( ENTR(6,MSP,NBEDGES(2)+NBEDGES(3)))
      ENTR=0.0D00
      DO I=1,NBEDGES(2)
         II=IBEDGE1(2,I)  !ii is the edge code of the cell
         JJ=IBEDGE1(1,I)  !jj is the code of cell
         C=0.0D00 !c is the wide of the edge
         DO KK=1,2
            C=C+(NODES(KK,ICELL(II,JJ))-NODES(KK,ICELL(ADD_ONE(II),JJ)))**2.0D00
         END DO
         C=DSQRT(C)
         DO L=1,MSP  
            ctheta=(ENTRYIN(1,I)*VFX+ENTRYIN(2,I)*VFY)/DSQRT(VFX**2.0D00+VFY*2.0D00)
            IF( DABS(ctheta) .LT. 1.0E-8) ctheta=0.0D00
            SN=DSQRT(VFX**2.0D00+VFY**2.0D00)*ctheta/VMP(L)  !Formula 12.4 in Bird's book
            AA=SN
            A=1.D00+DERF(AA)
            BB=DEXP(-SN**2.0D00)
            ENTR(3,L,I)=SN
            ENTR(4,L,I)=SN+DSQRT(SN**2.0D00+2.0D00)
            ENTR(5,L,I)=0.5D00*(1.D00+SN*(2.D00*SN-ENTR(4,L,I)))
            ENTR(6,L,I)=3.0D00*VMP(L)
            B=BB+SPI*SN*A
            !following A is the wide of edge            
            ENTR(1,L,I)=(FND*FSP(L)*VMP(L))*B*C/(FNUM*2.D00*SPI) !the formula (4.22)
            ENTR(2,L,I)=0.D00
         END DO
      END DO
      
      DO I=1,NBEDGES(3)
         II=IBEDGE2(2,I)  !ii is the edge code of the cell
         JJ=IBEDGE2(1,I)  !jj is the code of cell
         C=0.0D00 !c is the wide of the edge
         DO KK=1,2
            C=C+(NODES(KK,ICELL(II,JJ))-NODES(KK,ICELL(ADD_ONE(II),JJ)))**2.0D00
         END DO
         C=DSQRT(C)
         DO L=1,MSP  
            ctheta=(ENTRYOUT(1,I)*VFX+ENTRYOUT(2,I)*VFY)/DSQRT(VFX**2.0D00+VFY**2.0D00)
            IF( DABS(ctheta) .LT. 1.0E-8) ctheta=0.0D00
            SN=DSQRT(VFX**2.0D00+VFY**2.0D00)*ctheta/VMP(L)  !Formula 12.4 in Bird's book
            AA=SN
            A=1.D00+DERF(AA)
            BB=DEXP(-SN**2.0D00)
            ENTR(3,L,I+NBEDGES(2))=SN
            ENTR(4,L,I+NBEDGES(2))=SN+DSQRT(SN**2.0D00+2.0D00)
            ENTR(5,L,I+NBEDGES(2))=0.5D00*(1.D00+SN*(2.D00*SN-ENTR(4,L,I+NBEDGES(2))))
            ENTR(6,L,I+NBEDGES(2))=3.0D00*VMP(L)
            B=BB+SPI*SN*A
            ENTR(1,L,I+NBEDGES(2))=(FND*FSP(L)*VMP(L))*B*C/(FNUM*2.D00*SPI) !the formula (4.22)
            ENTR(2,L,I+NBEDGES(2))=0.D00
         END DO
      END DO
!--Set the initial molecules and sample
      MNM=IDINT(1.2D00*NMI)
      ALLOCATE (PX(MNM),PY(MNM),PTIM(MNM),IPCCELL(MNM),IPSP(MNM),ICREF(MNM),  &
         PV(3,MNM),STAT=ERROR)
      IF(MMVM .GT. 0)THEN
         ALLOCATE( IPVIB(MMVM,MNM))
         IPVIB=0
      END IF
      IF(MMRM .GT. 0)THEN
         ALLOCATE(PROT(MNM))
         PROT=0.0D00
      END IF
      PX=0.0D00;PY=0.0D00;PTIM=0.0D00
      IPCCELL=0;IPSP=0;ICREF=0
      
      
      NM=0
      DO I=1,MSP
  !--cycle species
            DO J=1,NCELLS
         !---cycle between cells
               L=CEILING(MOLSC*CELL(1,J)/MINCELL(1)/4.0D00*FSP(I))  
               DO K=1,4
            !----cycle in collision cell in a specific sampling cell
                  DO II=1,L
                  !---cycle between the molecules in a specifc collision cell
                     NM=NM+1
                     N=N+1
                     CALL POSITION_GEN(K,J,PX(NM),PY(NM))
                     PTIM(NM)=0.0D00
                     IPSP(NM)=I
                     IPCCELL(NM)=K+ICELL(8,J)
                     DO NSET=1,2
                        DO KK=1,3
                  !----cycle X,Y,Z direction
                           CALL RVELC(A,B,VMP(I))
                !---following code aims to make the sum of momentum approximates zero
                           IF( DMOM(KK) < 0.D00) THEN
                              IF( A<B) THEN
                                 BB=B
                              ELSE
                                 BB=A
                              END IF
                           ELSE
                              IF(A<B)THEN
                                 BB=A
                              ELSE
                                 BB=B
                              END IF
                           END IF
                           VB(KK,NSET)=BB         
                        END DO
                        IF (ISPR(1,I) > 0) CALL SROT(I,FTMP,ROTE(NSET))
                     END DO
                     A=(0.5D00*SP(5,I)*(VB(1,1)**2+VB(2,1)**2+VB(3,1)**2)+ROTE(1))/&
                        &(0.5D00*BOLTZ*FTMP)-3.D00-DFLOAT(ISPR(1,I))
                     B=(0.5D00*SP(5,I)*(VB(1,2)**2+VB(2,2)**2+VB(3,2)**2)+ROTE(2))/&
                        &(0.5D00*BOLTZ*FTMP)-3.D00-DFLOAT(ISPR(1,I))
               !-----follwing code aims to choose a more equilibrate state
                     IF(DENG < 0.D00)THEN
                        IF(A<B)THEN
                           KN=2
                        ELSE
                           KN=1
                        END IF
                     ELSE
                        IF(A<B)THEN
                           KN=1
                        ELSE
                           KN=2
                        END IF
                     END IF
                     DO KK=1,3
                        PV(KK,NM)=VB(KK,KN)
                        DMOM(KK)=DMOM(KK)+VB(KK,KN)*SP(5,I)
                     END DO
                     PV(1,NM)=PV(1,NM)+VFX
                     PV(2,NM)=PV(2,NM)+VFY
                     IF (KN == 1) DENG=DENG+A  
                     IF (KN == 2) DENG=DENG+B
                     IF (MMVM > 0) THEN
                       IF (ISPV(I) > 0) THEN
                         DO KK=1,ISPV(I)
                           CALL SVIB(I,FVTMP,IPVIB(KK,NM),KK)
                         END DO
                       END IF
                     END IF
                  END DO
               END DO
            END DO
      END DO
      WRITE (*,*) 'The initial number of molecules is',NM
      
      CALL INITIALISE_SAMPLES(0)

      deallocate(VB,ROTE)
     END SUBROUTINE
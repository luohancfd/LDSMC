      SUBROUTINE MOLECULES_ENTR
      USE MOLECS
      USE CALC
      USE CELLINFO
      USE GAS
      USE NODEINFO
      USE SAMPLES
      
      IMPLICIT NONE
      
      INTEGER :: I,J,K,L,M,N,NENT,CELL0,CELL1,STAT,TEMP,DRIFT,N0
      REAL*8 :: A,B,C,AA,BB,CC,U,VN,VP,XI,YI,X,Y,DX,DY,DTIM,RANF
      REAL*8,EXTERNAL ::VEC_PRODUCT
      INTEGER, EXTERNAL :: ADD_ONE,SEARCH_CCELL,DIS_DETECT
      !INFLOW BOUNDARY
      TEMP=0
      N0=NM
      DO I=1,NBEDGES(2)
         DO J=1,MSP
            A=ENTR(1,J,I)*DTM+ENTR(2,J,I)
            NENT=A
            TEMP=TEMP+NENT
            ENTR(2,J,I)=A-NENT
            IF (NENT .GT. 0) THEN
               DO K=1,NENT
                  STAT=1
                  DO WHILE( STAT .EQ. 1)
                     DO WHILE (NM .GE. MNM)
                        CALL EXTEND_MNM(1.001D00)
                     END DO
                     NM=NM+1
                     AA=DMAX1(0.D00,ENTR(3,J,I)-3.D00)   !这里得到的是分子速度/最几热运动速度的采样区间
                     BB=DMAX1(3.D00,ENTR(3,J,I)+3.D00)
                     A=-1.0D00
                     CALL RANDOM_NUMBER(RANF)
                     DO WHILE(A .LE. RANF)
                        CALL RANDOM_NUMBER(RANF)
                        B=AA+(BB-AA)*RANF
                        !B is the velocity/most probable thermal velocit of molecule normal to the boundary
                        U=B-ENTR(3,J,I)    !thermal velocity/most probable velocity
                        A=(2.D00*B/ENTR(4,J,I))*DEXP(ENTR(5,J,I)-U*U)
                     END DO
                     VN=B*VMP(J)
                     CALL RVELC(VP,PV(3,NM),VMP(J))
                     VP=VP+VEC_PRODUCT(ENTRYIN(1,I),ENTRYIN(2,I),VFX,VFY)
                     PV(1,NM)=VN*ENTRYIN(1,I)-VP*ENTRYIN(2,I)
                     PV(2,NM)=VN*ENTRYIN(2,I)+VP*ENTRYIN(1,I)
                     IF (ISPR(1,J) .GT. 0) CALL SROT(J,FTMP,PROT(NM))
                     IF (MMVM > 0) THEN
                        DO L=1,ISPV(J)
                           CALL SVIB(J,FVTMP,IPVIB(L,NM),L)
                        END DO
                     END IF
                     IPSP(NM)=J
                     CELL0=IBEDGE1(1,I)
                     !--advance the molecule into the flow
                     CALL RANDOM_NUMBER(RANF)
                     XI=RANF*(NODES(1,ICELL(ADD_ONE(IBEDGE1(2,I)),IBEDGE1(1,I)))-&
                        &NODES(1,ICELL(IBEDGE1(2,I),IBEDGE1(1,I))))
                     XI=XI+NODES(1,ICELL(IBEDGE1(2,I),IBEDGE1(1,I)))
                     YI=RANF*(NODES(2,ICELL(ADD_ONE(IBEDGE1(2,I)),IBEDGE1(1,I)))-&
                        &NODES(2,ICELL(IBEDGE1(2,I),IBEDGE1(1,I))))
                     YI=YI+NODES(2,ICELL(IBEDGE1(2,I),IBEDGE1(1,I)))
                     !The initial position (XI,YI) is on the boundary edge
                     !in order to obviate possible error in SEAR_CELL
                     !generate a little drift
                     CALL RANDOM_NUMBER(RANF)
                     DTIM=RANF*DTM                   
                     DRIFT=-3
                     X=XI+DTIM*PV(1,NM)*(10.0D00)**DRIFT
                     Y=YI+DTIM*PV(2,NM)*(10.0D00)**DRIFT
                     DO WHILE(DIS_DETECT(X,Y,CELL0) .EQ. 0)
                        DRIFT=DRIFT-1
                        X=XI+DTIM*PV(1,NM)*(10.0D00)**DRIFT
                        Y=YI+DTIM*PV(2,NM)*(10.0D00)**DRIFT
                     END DO
                     XI=X; YI=Y
                     X=XI+DTIM*PV(1,NM)*(1.0D00-10.0D00**DRIFT)
                     Y=YI+DTIM*PV(2,NM)*(1.0D00-10.0D00**DRIFT)                     
                     ENTMASS=ENTMASS+SP(5,J)
                     N=NM
                     CALL SEARCH_CELL(N,XI,YI,X,Y,CELL0,CELL1,DTIM,STAT)
                     IF( STAT .EQ. 0)THEN
                        PX(NM)=X
                        PY(NM)=Y
                        IPCCELL(NM)= SEARCH_CCELL(X,Y,CELL1)+ICELL(8,CELL1)
                     ELSE
                        L=IPSP(NM)
                        ENTMASS=ENTMASS-SP(5,L)  
                        NM=NM-1
                     END IF
                     !当分子移出时，NM->NM-1,STAT=1,则以上操作重新进行，重建NM分子
                     !当分子保留时，CELL1返回所在CELL，STAT=0跳出这个do while
                  END DO
               END DO
            END IF
         END DO
      END DO
      !OUT FLOW BOUNDARY
      DO I=1,NBEDGES(3)
         DO J=1,MSP
            A=ENTR(1,J,I+NBEDGES(2))*DTM+ENTR(2,J,I+NBEDGES(2))
            NENT=A
            ENTR(2,J,I+NBEDGES(2))=A-NENT
            IF (NENT .GT. 0) THEN
               DO K=1,NENT
                  STAT=1
                  DO WHILE( STAT .EQ. 1)
                     DO WHILE (NM .GE. MNM)
                        CALL EXTEND_MNM(1.001D00)
                     END DO
                     NM=NM+1
                     AA=DMAX1(0.D00,ENTR(3,J,I+NBEDGES(2))-3.D00)   !这里得到的是分子速度/最几热运动速度的采样区间
                     BB=DMAX1(3.D00,ENTR(3,J,I+NBEDGES(2))+3.D00)
                     A=-1.0D00
                     CALL RANDOM_NUMBER(RANF)
                     DO WHILE(A .LE. RANF)
                        CALL RANDOM_NUMBER(RANF)
                        B=AA+(BB-AA)*RANF
                        !B is the velocity/most probable thermal velocit of molecule normal to the boundary
                        U=B-ENTR(3,J,I+NBEDGES(2))    !thermal velocity/most probable velocity
                        A=(2.D00*B/ENTR(4,J,I+NBEDGES(2)))*DEXP(ENTR(5,J,I+NBEDGES(2))-U*U)
                     END DO
                     VN=B*VMP(J)
                     CALL RVELC(VP,PV(3,NM),VMP(J))
                     VP=VP+VEC_PRODUCT(ENTRYOUT(1,I),ENTRYOUT(2,I),VFX,VFY)
                     PV(1,NM)=VN*ENTRYOUT(1,I)-VP*ENTRYOUT(2,I)
                     PV(2,NM)=VN*ENTRYOUT(2,I)+VP*ENTRYOUT(1,I)
                     IF (ISPR(1,J) .GT. 0) CALL SROT(J,FTMP,PROT(NM))
                     IF (MMVM > 0) THEN
                        DO L=1,ISPV(J)
                           CALL SVIB(J,FVTMP,IPVIB(L,NM),L)
                        END DO
                     END IF
                     IPSP(NM)=J
                     CELL0=IBEDGE2(1,I)
                     !--advance the molecule into the flow
                     CALL RANDOM_NUMBER(RANF)
                     XI=RANF*(NODES(1,ICELL(ADD_ONE(IBEDGE2(2,I)),IBEDGE2(1,I)))-&
                        &NODES(1,ICELL(IBEDGE2(2,I),IBEDGE2(1,I))))
                     XI=XI+NODES(1,ICELL(IBEDGE2(2,I),IBEDGE2(1,I)))
                     YI=RANF*(NODES(2,ICELL(ADD_ONE(IBEDGE2(2,I)),IBEDGE2(1,I)))-&
                        &NODES(2,ICELL(IBEDGE2(2,I),IBEDGE2(1,I))))
                     YI=YI+NODES(2,ICELL(IBEDGE2(2,I),IBEDGE2(1,I)))
                     !The initial position (XI,YI) is on the boundary edge
                     !in order to obviate possible error in SEAR_CELL
                     !generate a little drift
                     CALL RANDOM_NUMBER(RANF)
                     DTIM=RANF*DTM
                     DRIFT=-3
                     X=XI+DTIM*PV(1,NM)*(10.0D00)**DRIFT
                     Y=YI+DTIM*PV(2,NM)*(10.0D00)**DRIFT
                     DO WHILE(DIS_DETECT(X,Y,CELL0) .EQ. 0)
                        DRIFT=DRIFT-1
                        X=XI+DTIM*PV(1,NM)*(10.0D00)**DRIFT
                        Y=YI+DTIM*PV(2,NM)*(10.0D00)**DRIFT
                     END DO
                     XI=X; YI=Y
                     X=XI+DTIM*PV(1,NM)*(1.0D00-10.0D00**DRIFT)
                     Y=YI+DTIM*PV(2,NM)*(1.0D00-10.0D00**DRIFT)
                     ENTMASS=ENTMASS+SP(5,J)
                     N=NM
                     CALL SEARCH_CELL(N,XI,YI,X,Y,CELL0,CELL1,DTIM,STAT)
                     !lout useless
                     IF( STAT .EQ. 0)THEN
                        PX(NM)=X
                        PY(NM)=Y
                        IPCCELL(NM)= SEARCH_CCELL(X,Y,CELL1)+ICELL(8,CELL1)
                     ELSE
                        L=IPSP(NM)
                        ENTMASS=ENTMASS-SP(5,L)  
                        NM=NM-1
                     END IF
                     !当分子移出时，NM->NM-1,STAT=1,则以上操作重新进行，重建NM分子
                     !当分子保留时，CELL1返回所在CELL，STAT=0跳出这个do while
                  END DO
               END DO
            END IF
         END DO
      END DO
      WRITE(*,*)  "FINISH Molecule Enter","Enter=",NM-N0
      WRITE (9,*) "FINISH Molecule Enter","Enter=",NM-N0
   END SUBROUTINE
   
      
                  
                  

                  
            
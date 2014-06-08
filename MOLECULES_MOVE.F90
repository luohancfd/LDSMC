      SUBROUTINE MOLECULES_MOVE
!----The subroutine is mained to move molecule and calculate the collision between the molecule and the surface
!----It also moves new molecules in at the inflow boundary and remove molecules out at the outflow boundary
      USE CELLINFO
      USE MOLECS
      USE CALC
      USE CONST
      USE NODEINFO
      USE OUTPUT
      USE GAS
      USE SAMPLES
      !$ USE OMP_LIB
      IMPLICIT NONE
      INTEGER :: CELL1,CELL0,CCELL0,CCELL1,N,lout,I,J,L,thread,chunk
      INTEGER,ALLOCATABLE :: STAT(:)
      REAL*8 :: DX,DY,XI,YI,X,Y,DTIM,TI
      INTEGER,EXTERNAL :: SEARCH_CCELL

!--TI initial cell time
!--DTIM time interval for the move
!--LOUT--number of molecule touch bound
!--CCELL0--the initial collision cell which molecule was
!--CELL0--then initial sampling cell which molecule was
!--NARC--the sum of arc region


      allocate(STAT(NM))
      !thread=OMP_get_max_threads()
      !write(*,*) "Thread=",thread

      lout=0
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(CELL0,CCELL0,CELL1,CCELL1,DTIM,TI,XI,DX,YI,DY,X,Y,L) 

      !$OMP DO REDUCTION(+:CSS,CSSS,ENTMASS,TOTMOV,lout)
      DO N=1,NM
         IF(N .EQ. 47999)THEN
            I=0
         END IF
         
         CCELL0=IPCCELL(N)
         CELL0=ICCELL(3,CCELL0)
         IF (IMTS == 0) DTIM=DTM
         IF (IMTS == 1) DTIM=2.D00*CCELL(3,CCELL0)
                  
         IF (FTIME-PTIM(N) > 0.5*DTIM) THEN
            TI=PTIM(N)
            PTIM(N)=TI+DTIM
            TOTMOV=TOTMOV+1
            
            XI=PX(N)
            DX=PV(1,N)*DTIM
            X=XI+DX
            
            YI=PY(N)
            DY=PV(2,N)*DTIM
            Y=YI+DY
            
            CALL SEARCH_CELL(N,XI,YI,X,Y,CELL0,CELL1,DTIM,STAT(N))
            
            IF(STAT(N) .EQ. 0)THEN
               PX(N)=X
               PY(N)=Y
               CCELL1=SEARCH_CCELL(X,Y,CELL1)
               IPCCELL(N)=CCELL1+ICELL(8,CELL1)
            ElSE
               lout=lout+1
               L=IPSP(N)
               ENTMASS=ENTMASS-SP(5,L)  
            END IF
         END IF
      END DO
      !$OMP END DO
      !$OMP END PARALLEL

      

      
      
      WRITE(*,"(A12,I20,A1,I20,A18)") "Finsh moving",lout,'/',NM,"molecules move out"
      WRITE (9,"(A12,I20,A1,I20,A18)") "Finsh moving",lout,'/',NM,"molecules move out"
      

      I=1;J=NM
      
      IF (MMRM > 0 .AND. MMVM >0) THEN
         DO WHILE( lout .GT. 0)
            IF(STAT(I) .EQ. 1)THEN
               DO WHILE((STAT(J) .EQ. 1) .AND. (I .LT. J))
                  J=J-1
                  lout=lout-1
                  NM=NM-1
               END DO
               IF( I .ne. J)THEN
                  PX(I)=PX(J)
                  PY(I)=PY(J)
                  DO N=1,3
                     PV(N,I)=PV(N,J)
                  END DO
                  IPCCELL(I)=IPCCELL(J)
                  IPSP(I)=IPSP(J)
                  PROT(I)=PROT(J)
                  DO N=1,MMVM
                     IPVIB(N,I)=IPVIB(N,J)
                  END DO
                  PTIM(I)=PTIM(J)
               END IF
               lout=lout-1
               nm=nm-1
               J=J-1
            END IF
            I=I+1
         END DO
      ELSE IF(MMRM > 0)THEN
         DO WHILE( lout .GT. 0)
            IF(STAT(I) .EQ. 1)THEN
               DO WHILE((STAT(J) .EQ. 1) .AND. (I .LT. J))
                  J=J-1
                  lout=lout-1
                  NM=NM-1
               END DO
               IF( I .ne. J)THEN
                  PX(I)=PX(J)
                  PY(I)=PY(J)
                  DO N=1,3
                     PV(N,I)=PV(N,J)
                  END DO
                  IPCCELL(I)=IPCCELL(J)
                  IPSP(I)=IPSP(J)
                  PROT(I)=PROT(J)
                  PTIM(I)=PTIM(J)
               END IF
               lout=lout-1
               nm=nm-1
               J=J-1
            END IF
            I=I+1
         END DO
      ELSE
         DO WHILE( lout .GT. 0)
            IF(STAT(I) .EQ. 1)THEN
               DO WHILE((STAT(J) .EQ. 1) .AND. (I .LT. J))
                  J=J-1
                  lout=lout-1
                  NM=NM-1
               END DO
               IF( I .ne. J)THEN
                  PX(I)=PX(J)
                  PY(I)=PY(J)
                  DO N=1,3
                     PV(N,I)=PV(N,J)
                  END DO
                  IPCCELL(I)=IPCCELL(J)
                  IPSP(I)=IPSP(J)
                  PTIM(I)=PTIM(J)
               END IF
               lout=lout-1
               nm=nm-1
               J=J-1
            END IF
            I=I+1
         END DO
      END IF
      DEALLOCATE(STAT)
   
   END SUBROUTINE
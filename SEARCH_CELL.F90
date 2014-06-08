      SUBROUTINE SEARCH_CELL(N,XI,YI,X,Y,CELL0,CELL1,DTIM,STAT)
      USE NODEINFO
      USE CELLINFO

      IMPLICIT NONE
      !search the cell of the molecule according to XI,YI,X,Y,CELL0
      !if the molecule move out ,remove it,and copy NM molecule to N
      
      INTEGER :: N,CELL0,CELL1,STAT
      !N---code of molecule
      !CELL0----initial sampling cell
      !CELL1----final sampling cell
      !DTIM-----time reaming
      !STAT-----STAT of the molecule, 0 still in the sampling cell 1 removed
      INTEGER :: K0,K,KK,JJEND,JJSTART,JJSTEP,JJ
      REAL*8 :: XI,YI,X,Y,DTIM,X1,Y1,X2,Y2
      LOGICAL,EXTERNAL :: INTERSECT_DETECT
      INTEGER,EXTERNAL :: ADD_ONE,MINUS_ONE
      INTEGER :: N1,N2,STAT2
      
      STAT2=0 !STAT2=0 if molecule's position is founded
      STAT=0  !STAT =1 if molecule is removed
      DO WHILE( STAT2 .EQ. 0)
         K0=CELL0
         K=CELL0

         JJSTART=1    !JJSTART is the beginning of the loop of edge
         JJEND=3     ! JJEND is the end of the loop
         JJSTEP=1    !JSTEP is the loop step
         KK=1
         !KK=0 keep in current cell
         !  =1 go into neighboring cell
         !  =2 touch bound

         
         DO WHILE( KK .EQ. 1)
            KK=0
            DO JJ=JJSTART,JJEND,JJSTEP
               N2=JJ+1
               IF(N2 .EQ. 4) N2=1
               N2=ICELL(N2,K)
               N1=ICELL(JJ,K)
               X1=NODES(1,N1)
               Y1=NODES(2,N1)
               X2=NODES(1,N2)
               Y2=NODES(2,N2)
               IF(INTERSECT_DETECT(X1,Y1,X2,Y2,XI,YI,X,Y))THEN
                  IF(ICELL(3+JJ,K).EQ.0)THEN
                     KK=1
                     EXIT
                  ELSE
                     KK=2
                     EXIT
                  END IF
               END IF
            END DO
            IF( KK .EQ. 0)THEN
               CELL1=K
               DTIM=0.0D00
               STAT2=1
            ELSE IF( KK .EQ. 1)THEN
               K0=K
               K=ICELL(8+JJ,K)
               JJSTART=1
               DO WHILE( ICELL(JJ,K0) .NE. ICELL(JJSTART,K))
                  JJSTART=JJSTART+1
               END DO
               JJEND=JJSTART+1
               IF(JJEND .EQ. 4) JJEND=1
               JJSTEP=JJEND-JJSTART
            ELSE IF( KK .EQ. 2)THEN
               CALL TOUCH_BOUND(N,K,JJ,XI,YI,X,Y,X1,Y1,X2,Y2,DTIM,CELL0,STAT)
               IF(STAT .EQ. 1) RETURN
            END IF
         END DO
      END DO
      
      END SUBROUTINE
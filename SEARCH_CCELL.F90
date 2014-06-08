      FUNCTION SEARCH_CCELL(X,Y,II)
   !!!!---this function return a INTEGER variable--!!!!
   !this function is used to detect the collision cell where the molecule is
   !----X,Y------the coordinate of the molecule
   !----II-------the cell

      USE CELLINFO
      USE NODEINFO
      USE GAS
      IMPLICIT NONE
      
      
      REAL*8 :: X,Y,x1,x2,x3,y1,y2,y3,x10,y10,x20,y20,x30,y30
      INTEGER :: II,I,K
      INTEGER :: SEARCH_CCELL
      !REAL*8,EXTERNAL :: VEC_PRODUCT !Å×Æúµô£¬Ì«ºÄÄÚ´æ
      LOGICAL d1,d2,d3
   !---X---coordiante x of molecule
   !---Y---coordinate y of molecule
   !---II--code of cell
      
   !---x10,y10,x20,y20,x30,y30-----the coordinate of the cell,clockwise
       x10=NODES(1,ICELL(1,II))
       x20=NODES(1,ICELL(2,II))
       x30=NODES(1,ICELL(3,II))
       y10=NODES(2,ICELL(1,II))
       y20=NODES(2,ICELL(2,II))
       y30=NODES(2,ICELL(3,II))
   !---x1,y1,x2,y2,x3,y3-----the coordinate of the middle collision cell
       x1=(x20+x30)*0.5D00
       x2=(x10+x30)*0.5D00
       x3=(x10+x20)*0.5D00
       y1=(y20+y30)*0.5D00
       y2=(y10+y30)*0.5D00
       y3=(y10+y20)*0.5D00

       
       
       !d1=(VEC_PRODUCT(x1-X,y1-Y,x2-X,y2-Y) .GT. 0)
       d1=(((x1-X)*(y2-Y)-(x2-X)*(y1-Y)) .GT. 0.0D00)
       !d2=(VEC_PRODUCT(x2-X,y2-Y,x3-X,y3-Y) .GT. 0)
       d2=(((x2-X)*(y3-Y)-(x3-X)*(y2-Y)) .GT. 0.0D00)
       !d3=(VEC_PRODUCT(x3-X,y3-Y,x1-X,y1-Y) .GT. 0)
       d3=(((x3-X)*(y1-Y)-(x1-X)*(y3-Y)) .GT. 0.0D00)
       
       IF( II .EQ. 17422)THEN
          II=17422
       END IF
       
       SELECT CASE(d1)
       CASE(.true.)
          SELECT CASE(d2)
          CASE(.true.)
             SELECT CASE(d3)
             CASE(.true.)
                WRITE(*,*) "ERROR!"
                OPEN(525,FILE="DIAG_SEARCH.TXT")
                WRITE(525,*)'SEARCH CELL ERROR',II
                CLOSE(525)
                STOP
             CASE(.false.)
                !SEARCH_CCELL=2
                WRITE(*,*) "ERROR!"
               OPEN(525,FILE="DIAG_SEARCH.TXT")
                WRITE(525,*)'SEARCH CELL ERROR',II
                CLOSE(525)
                STOP
             END SELECT
          CASE(.false.)
             SELECT CASE(d3)
             CASE(.true.)
                !SEARCH_CCELL=6
                WRITE(*,*) "ERROR!"
                OPEN(525,FILE="DIAG_SEARCH.TXT")
                WRITE(525,*)'SEARCH CELL ERROR',II
                CLOSE(525)
                STOP
             CASE(.false.)
                !SEARCH_CCELL=1
                SEARCH_CCELL=3
             END SELECT
          END SELECT
       CASE(.false.)
           SELECT CASE(d2)
          CASE(.true.)
             SELECT CASE(d3)
             CASE(.true.)
                !SEARCH_CCELL=4
                WRITE(*,*) "ERROR!"
                 OPEN(525,FILE="DIAG_SEARCH.TXT")
                WRITE(525,*)'SEARCH CELL ERROR',II
                CLOSE(525)
                STOP
             CASE(.false.)
                !SEARCH_CCELL=3
                SEARCH_CCELL=1
             END SELECT
          CASE(.false.)
             SELECT CASE(d3)
             CASE(.true.)
                !SEARCH_CCELL=5
                SEARCH_CCELL=2
             CASE(.false.)
                !SEARCH_CCELL=7
                SEARCH_CCELL=4
             END SELECT
          END SELECT
       END SELECT
       RETURN
       END FUNCTION SEARCH_CCELL
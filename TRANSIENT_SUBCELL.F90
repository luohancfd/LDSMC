        MODULE TRANSC
        IMPLICIT NONE
        !only used in subroutine collision and TRANSIENT_SUBCELL for transient cell
        INTEGER,SAVE :: NITSC,NDIVX,NDIVY
        INTEGER,DIMENSION(:,:,:),ALLOCATABLE,SAVE:: ITSC
        INTEGER,DIMENSION(:,:),ALLOCATABLE,SAVE :: IPTSC,NTSC
        !$OMP THREADPRIVATE(NITSC,ITSC,IPTSC,NTSC,NDIVX,NDIVY)


        !---NTSC(Y,X)----------number of molecule in (X,Y) transient sub cell
        !---ITSC(I,Y,X)---- transient sub cell index
        !--------X----the index in x direction
        !--------Y----the index in y direction
        !--------I----the code of molecule in transient sub cell
        !-------------I=1,NITSC
        !---IPTSC(I,M)-----molecule to transient cell
        !-------I=1----x index of molecule M
        !-------I=2----y index of molecule M
        !-------I=3----index of molecule in transient cell
        !---NDIVX--------number of transient sub cell in x direction
        !---NDIVY--------number of transient sub cell in y direction
        CONTAINS
            SUBROUTINE TRANSIENT_SUBCELL_GEN(I)
            USE NODEINFO
            USE CELLINFO
            USE MOLECS


            IMPLICIT NONE

            INTEGER :: I
            INTEGER :: II,K,L,J,XC,YC
            INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: ITSC_BACKUP
            REAL*8 :: x1,x2,x3,y1,y2,y3
            REAL*8 :: x10,x20,x30,y10,y20,y30
            REAL*8 :: XI,YI,X,Y,DX,DY
            !---I---------------code of collision cell
            !---II--------------code of sampling cell
            !---NDIVX-----------number of transient sub cell in x direction
            !---NDIVY-----------number of transient sub cell in y direction
            !---ITSC_BACKUP-----BACKUP OF ITSC,in case the I of ITSC is not enough
            !---XI,YI-----------the minimum boundary of transient sub cell
            !---X ,Y -----------the maximum boundary of transient sub cell
            II=ICCELL(3,I)
            !---x1,y1,x2,y2,x3,y3-----the coordinate of the cell,clockwise
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
            SELECT CASE(I-ICELL(8,II))
            CASE(1)
                x1=x10; y1=y10
            CASE(2)
                x2=x20; y2=y20
            CASE(3)
                x3=x30; y3=y30
            END SELECT

            XI=DMIN1(x1,x2,x3); YI=DMIN1(y1,y2,y3)
            X=DMAX1(x1,x2,x3) ; Y=DMAX1(y1,y2,y3)

            !  for triangular cell, num_molecule<0.5*ndivx*ndivy
            !  this is used to ensure there is nearly a molecule in every transient sub cell

            !  the following aims to make DX=DY
            DX=DSQRT((X-XI)*(Y-YI)/(2.0D00*DFLOAT(ICCELL(2,I))))
            DY=DX
            NDIVX=(X-XI)/DX
            DX=(X-XI)/DFLOAT(NDIVX)
            NDIVY=(Y-YI)/DY
            DY=(Y-YI)/DFLOAT(NDIVY)

            NITSC=2
            ALLOCATE(ITSC(NITSC,NDIVY,NDIVX),NTSC(NDIVY,NDIVX),IPTSC(3,NM))
            ITSC=0; NTSC=0
            L=ICCELL(1,I)
            DO J=1+L,ICCELL(2,I)+L
                K=ICREF(J)  !K: code of molecule
                XC=INT((PX(K)-XI-DX*1D-3)/DX)+1  !-DX*1D-3 aims to obviate PX(K)=X
                YC=INT((PY(K)-YI-DY*1D-3)/DY)+1
                IPTSC(1,K)=XC
                IPTSC(2,K)=YC
                NTSC(YC,XC)=NTSC(YC,XC)+1
                IF(NTSC(YC,XC) .GT. NITSC)THEN
                    !increase the index code of molecule in the transient sub cell
                    ALLOCATE(ITSC_BACKUP(NITSC,NDIVY,NDIVX))
                    ITSC_BACKUP=ITSC
                    DEALLOCATE(ITSC)
                    NITSC=NITSC+1
                    ALLOCATE(ITSC(NITSC,NDIVY,NDIVX))
                    ITSC(NITSC:NITSC,:,:)=0
                    ITSC(1:NITSC-1,:,:)=ITSC_BACKUP
                    DEALLOCATE(ITSC_BACKUP)
                END IF
                ITSC(NTSC(YC,XC),YC,XC)=K
                IPTSC(3,K)=NTSC(YC,XC)
            END DO
            RETURN
            END SUBROUTINE
        
            SUBROUTINE TRANSIENT_SUBCELL_SELECT(M1,M2)
            USE MOLECS
            USE CELLINFO

            IMPLICIT NONE

            INTEGER :: M1,M2
            !----M1-----already selected molecule for collsion
            !----M2-----the molecule to be selected
            !----DEL1---the first deleted molecule
            !----DEL2---the second deleted molecule
            INTEGER :: I,J,K,L,JJ,JEND,X,Y,TSCLINK_LEN
            INTEGER :: M1TX,M1TY,MPTX,MPTY
            INTEGER,ALLOCATABLE :: TSCLINK(:)

            REAL*8 :: RANF
            INTEGER :: MP,M1SC,MPSC
            !----I,J,K------------loop number
            !----JJ,JEND----loop control
            !----M1TX,M1TY--------molecule 1 transient sub cell X,Y coordinate
            !----M2T--------------temporal molecule 2
            !----X,Y--------------current
            !----TSCLINK----------the link of transient sub cells which has molecule
            !----TSCLINK_LEN------link length
            !-------------------------------------------------------------------------
            !----NTSC_BACKUP,ITSC_BACKUP is used to back up NTSC,ITSC
            !----NUM_DEL---------if M1 and IPCP(M1) in the same transient sub cell =1 else =2



            !----Delete M1 and IPCP(M1) from transient sub cell grid
            M1TX=IPTSC(1,M1)
            M1TY=IPTSC(2,M1)
            M1SC=IPTSC(3,M1)
            
            MP=IPCP(M1)
            MPTX=0;MPTY=0;MPSC=0
            IF(MP .NE. 0)THEN
                MPTX=IPTSC(1,MP)
                MPTY=IPTSC(2,MP)
                MPSC=IPTSC(3,MP)
            END IF
            
          

            !-----start to select molecule 2
            M2=0
            JJ=0  !layer of transient sub cell molecules from which molecule 2 is chosen
            DO WHILE ( M2 .EQ. 0)
                JEND=8*JJ
                IF(JJ .EQ. 0) JEND=1
                ALLOCATE(TSCLINK(JEND))
                TSCLINK=0
                TSCLINK_LEN=0
                !----loop between the transient sub cells in the same layer in order to create TSCLINK
                IF( JJ .EQ. 0)THEN
                    IF(ITSC(1,M1TY,M1TX) .NE. 0)THEN
                        TSCLINK_LEN=1
                        TSCLINK(1)=1
                    END IF
                ELSE
                    DO J=1,JEND
                        !------find X,Y index for transient sub cell
                        IF( J .LE. JJ*2+1)THEN
                            X=M1TX-JJ+J-1
                            Y=M1TY-JJ
                        ELSE IF( J .LE. JJ*4+1)THEN
                            X=M1TX+JJ
                            Y=M1TY-3*JJ+J-1
                        ELSE IF( J .LE. JJ*6+1) THEN
                            X=M1TX+5*JJ-J+1
                            Y=M1TY+JJ
                        ELSE
                            X=M1TX-JJ
                            Y=M1TY+7*JJ-J+1
                        END IF
                        !-------------------------
                        IF( (Y>= 1) .AND. (X>=1) .AND. (Y<=NDIVY).AND. (X<=NDIVX))THEN
                            IF(ITSC(1,Y,X) .NE. 0)THEN
                                TSCLINK_LEN=TSCLINK_LEN+1
                                TSCLINK(TSCLINK_LEN)=J
                            END IF
                        END IF
                    END DO
                END IF

                !--------select second molecule
                IF( TSCLINK_LEN .EQ. 0)THEN
                    JJ=JJ+1   !no molecule in this layer, go to next layer
                    DEALLOCATE(TSCLINK)
                ELSE
                    CALL RANDOM_NUMBER(RANF)
                    L=INT(TSCLINK_LEN*RANF-0.0001)+1 !code in link
800                 J=TSCLINK(L)                     !code of transient sub cell
                    IF( J .LE. JJ*2+1)THEN
                        X=M1TX-JJ+J-1
                        Y=M1TY-JJ
                    ELSE IF( J .LE. JJ*4+1)THEN
                        X=M1TX+JJ
                        Y=M1TY-3*JJ+J-1
                    ELSE IF( J .LE. JJ*6+1) THEN
                        X=M1TX+5*JJ-J+1
                        Y=M1TY+JJ
                    ELSE
                        X=M1TX-JJ
                        Y=M1TY+7*JJ-J+1
                    END IF
                    
                    IF( JJ .EQ. 0)THEN
                       IF( NTSC(Y,X) .EQ. 2)THEN
                          IF(MPTX .NE. M1TX)THEN
                             M2=ITSC(3-M1SC,M1TY,M1TX)
                          ELSE
                             JJ=JJ+1   !no molecule in this layer, go to next layer
                             DEALLOCATE(TSCLINK)
                          END IF
                       ELSE IF( NTSC(Y,X) .GT. 2)THEN
                          DO WHILE(M2 .EQ. 0)
                             CALL RANDOM_NUMBER(RANF)
                             K=INT(NTSC(Y,X)*RANF-0.0001)+1  !code in transient sub cell
                             M2=ITSC(K,Y,X)                !code of molecule
                             IF((M2 .EQ. M1) .OR. (M2 .EQ. MP))THEN
                                M2=0
                             END IF
                          END DO
                       ELSE
                          JJ=JJ+1   !no molecule in this layer, go to next layer
                          DEALLOCATE(TSCLINK)
                       END IF
                    ELSE
                       K=INT(NTSC(Y,X)*RANF-0.0001)+1  !code in transient sub cell
                       M2=ITSC(K,Y,X)                !code of molecule
                       IF( M2 .EQ. MP)THEN
                          M2=0
                          IF(TSCLINK_LEN .NE. 1)THEN
                             L=L+1
                             IF(L .GT. TSCLINK_LEN)THEN
                                L=L-2
                             END IF
                             GOTO 800
                          END IF
                          JJ=JJ+1   !no molecule in this layer, go to next layer
                          DEALLOCATE(TSCLINK)
                       END IF
                    END IF     
                END IF
            END DO

       
            RETURN
            END SUBROUTINE

            SUBROUTINE TRANSIENT_SUBCELL_DEL
            IMPLICIT NONE
            DEALLOCATE(ITSC,IPTSC,NTSC)

            NITSC=0
            RETURN
            END SUBROUTINE
        END MODULE




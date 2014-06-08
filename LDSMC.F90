        PROGRAM LDSMC
        USE NODEINFO
        USE CELLINFO
        USE BOUNDINFO
        USE CALC
        USE CONST
        USE GAS
        USE MOLECS
        IMPLICIT NONE
        
        INTEGER :: ERROR,IRUN,SRUN,MOVS,SAMPS
        REAL*8 :: A
   !--MOVS the number of moves since last sample
   !--SAMPS the number of samples since last output
        OPEN (9,FILE='DIAG.TXT',FORM='FORMATTED',STATUS='REPLACE')
        WRITE (9,*,IOSTAT=ERROR)
        IF (ERROR .NE. 0) THEN
           WRITE (*,*) 'Stop the LDSMC.EXE that is already running and try again'
           STOP
        ELSE
           WRITE (9,*) 'File DIAG.TXT has been opened'
        END IF
        
        WRITE (*,*) 'Send 1 to continue the current run, or'
        WRITE (*,*) 'send 2 to start a new run :-'
        READ (*,*) IRUN
		
        IF (IRUN .EQ. 1) THEN
           WRITE(*,*) 'Send 1 to continue current sample'
           WRITE(*,*) 'Send 2 to start new sample'
           READ(*,*) SRUN
        END IF
        
        IF (IRUN .EQ. 2) WRITE (9,*) 'Starting a new run'
        
        IF( IRUN .EQ. 1)THEN
           CALL READ_RESTART
           IF(SRUN .EQ. 1)THEN
              WRITE(9,*) 'Continue current run and continue current sample'
           ELSE
              CALL INITIALISE_SAMPLES(1)
              WRITE(9,*) 'Start new sample and continue current sample'
           END IF
           TSAMP=FTIME+DTSAMP
           TOUT=FTIME+DTOUT
        END IF
        
        IF( IRUN .EQ. 2)THEN
           !Read initial parameter
           !Deal with mesh
           CALL READ_DATA
           
           CALL PREPROCESSING
           
           !
           CALL SET_INITIAL_STATE
           WRITE (*,*) 'INITIAL STATE SET'
           WRITE (9,*) 'INITIAL STATE SET'
           !CALL WRITE_RESTART
        END IF
        MOVS=0
        SAMPS=0
        
        TOTCOLI=TOTCOL
        TOTMOVI=TOTMOV
        
        !CALL SYSTEM_CLOCK(COUNT=MCOMPTIMI)

        DO WHILE (FTIME < TLIM)
           FTIME=FTIME+DTM
           WRITE (9,*) 'TIME',FTIME,' Number of molecule',NM,' COLLS',TOTCOL
           WRITE (*,*) 'TIME',FTIME,' Number of molecule',NM,' COLLS',TOTCOL
           
           CALL MOLECULES_MOVE
           
           CALL MOLECULES_ENTR
           
           CALL INDEX_MOLS
           
           CALL COLLISION

           IF ((DABS(FTIME - TSAMP) .LT. DTM*0.5D00) .OR. ((ISF .EQ. 1) .AND. (MOVS .GT. 10))) THEN
              CALL SAMPLE_FLOW
              TSAMP=TSAMP+DTSAMP

              MOVS=0
              SAMPS=SAMPS+1
              
              IF((DABS(FTIME - TOUT) .LT. DTM*0.5D00).OR.((ISF == 1).AND.(SAMPS > 50))) THEN
                 WRITE (9,*) 'TIME',FTIME,"OUT------------------------"
                 WRITE (*,*) 'TIME',FTIME,"OUT------------------------"

                 SAMPS=0                 
                 CALL OUTPUT_RESULT
                 TOUT=FTIME+DTOUT
                 

                 CALL WRITE_RESTART

              END IF               

           END IF
        END DO
        
   end program
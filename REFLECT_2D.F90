      SUBROUTINE REFLECT_2D(N,OCELL,OEDGE)
      USE CELLINFO
      USE MOLECS
      USE SAMPLES
      USE GAS
      
      
      USE CONST
      
      INTEGER :: LS,NBE,I,J,L,INVIB,N,K,OCELL,OEDGE
      REAL*8 :: UN,UP,ETI,ERI,EVI,VMPS,RANF
      REAL*8,EXTERNAL :: VEC_PRODUCT
      INTEGER,EXTERNAL :: ADD_ONE
      !-----OCELL---cell code
      !-----OEDGE---edge code
      !-----LS---molecule specie code
      !-----NBE--code of boundary edge
      !-----UN---normal velocity to the edge
      !-----UP---parallel velocity to the edge
      !-----ETI.ERI,EVI,EEI incident trans, rot, vib and total energy
      !-----VMPS most probable molecular velocity at surface temperature
      LS=IPSP(N)
      NBE=IS(OEDGE,OCELL)


      
      !---boundary is always at the right side of the boundary edge
      UN=PV(1,N)*SURFACEVEC(1,NBE)+PV(2,N)*SURFACEVEC(2,NBE)
      UP=PV(2,N)*SURFACEVEC(1,NBE)-PV(1,N)*SURFACEVEC(2,NBE)
      !---UP>0, the parallel velocity is in the same direction as the direction of boundary edge
      !---UP<0, the parallel velocity is in the opposite direction as the direction of boundary edge

      CSS(1,NBE,LS,1)=CSS(1,NBE,LS,1)+1.D00
      CSS(2,NBE,LS,1)=CSS(2,NBE,LS,1)-SP(5,LS)*UN 
      CSS(3,NBE,LS,1)=CSS(3,NBE,LS,1)-SP(5,LS)*UP 
      CSS(4,NBE,LS,1)=CSS(4,NBE,LS,1)-SP(5,LS)*PV(3,N)

      ETI=0.5*SP(5,LS)*(PV(1,N)**2.0D00+PV(2,N)**2.0D00+PV(3,N)**2.0D00)
      CSS(5,NBE,LS,1)=CSS(5,NBE,LS,1)+ETI

      IF (ISPR(1,LS) > 0) THEN
         CSS(6,NBE,LS,1)=CSS(6,NBE,LS,1)+PROT(N)
      END IF

      EVI=0
      IF (MMVM > 0) THEN
         INVIB=0
         IF (ISPV(LS) > 0) THEN
            CALL RANDOM_NUMBER(RANF)
            IF (RANF < SSSEG(8,LS,NBE)) INVIB=1
            IF (INVIB ==1) THEN
               DO K=1,ISPV(LS)
                  EVI=EVI+IPVIB(K,N)*BOLTZ*SPVM(1,K,LS)
                  CSS(7,NBE,LS,1)=CSS(7,NBE,LS,1)+IPVIB(K,N)*BOLTZ*SPVM(1,K,LS)
               END DO
            END IF
         END IF
      END IF
      
      CSSS(1,NBE)=CSSS(1,NBE)+1.0D00/DABS(UN)
      CSSS(2,NBE)=CSSS(2,NBE)+SP(5,LS)/DABS(UN)
      CSSS(3,NBE)=CSSS(3,NBE)+SP(5,LS)*UP/DABS(UN)
      CSSS(4,NBE)=CSSS(4,NBE)+SP(5,LS)*&
         &(UN**2.0D00+UP**2.0D00+PV(3,N)**2.0D00)/DABS(UN)
      IF (ISPR(1,LS) > 0) THEN
        CSSS(5,NBE)=CSSS(5,NBE)+PROT(N)/DABS(UN)
        CSSS(6,NBE)=CSSS(6,NBE)+ISPR(1,LS)/DABS(UN)
      END IF
      
      CALL RANDOM_NUMBER(RANF)
      IF( RANF .GT. SSSEG(9,LS,NBE))THEN   !diffuse but not specular reflection
         VMPS=SQRT(2.0D00*BOLTZ*SSEG(1,NBE)/SP(5,LS))
         IF (SSSEG(4,LS,NBE) < 0.) THEN
            !--diffuse reflection
            CALL RANDOM_NUMBER(RANF)
            UN=DSQRT(-LOG(RANF))*VMPS
            CALL RVELC(UP,PV(3,N),VMPS)
            IF (ISPR(1,LS) > 0) THEN
               CALL RANDOM_NUMBER(RANF)
               IF (SSSEG(7,LS,NBE) > RANF) THEN
                  CALL SROT(LS,SSEG(1,NBE),PROT(N))
               END IF
            END IF
         ELSE
            !CLL reflection
            CALL CLL(N,LS,UN,UP,NBE,SSEG(1,NBE),VMPS)
         END IF
         
         PV(1,N)=UN*SURFACEVEC(1,NBE)-UP*SURFACEVEC(2,NBE)
         PV(2,N)=UN*SURFACEVEC(2,NBE)+UP*SURFACEVEC(1,NBE)
         IF (MMVM > 0) THEN
          IF (ISPV(LS) > 0) THEN
             IF (INVIB .EQ. 1) THEN
                DO I=1,ISPV(LS)
                   CALL SVIB(LS,SSEG(1,NBE),IPVIB(I,N),I)
                END DO
             END IF
          END IF
         END IF
      ELSE  !specular reflection
         UN=-UN
         PV(1,N)=UN*SURFACEVEC(1,NBE)-UP*SURFACEVEC(2,NBE)
         PV(2,N)=UN*SURFACEVEC(2,NBE)+UP*SURFACEVEC(1,NBE)
      END IF
      
      CSS(1,NBE,LS,2)=CSS(1,NBE,LS,2)+1.00D00
      CSS(2,NBE,LS,2)=CSS(2,NBE,LS,2)+SP(5,LS)*UN
      CSS(3,NBE,LS,2)=CSS(3,NBE,LS,2)+SP(5,LS)*UP
      CSS(4,NBE,LS,2)=CSS(4,NBE,LS,2)+SP(5,LS)*PV(3,N)
      CSS(5,NBE,LS,2)=CSS(5,NBE,LS,2)+(UP**2.0D00+UN**2.0D00+PV(3,N)**2.0D00)*0.5D00*SP(5,LS)
      IF(ISPR(1,LS) .GT. 0)THEN
         CSS(6,NBE,LS,2)=CSS(6,NBE,LS,2)+PROT(N)
      END IF
      
      
      IF (MMVM > 0) THEN
         IF (ISPV(LS).GT.0) THEN
            IF (INVIB == 1) THEN
               DO K=1,ISPV(LS)
                  CSS(7,NBE,LS,2)=CSS(7,NBE,LS,2)+IPVIB(K,N)*BOLTZ*SPVM(1,K,LS)
               END DO
            END IF
         END IF
      END IF
      
      IF(UN .LT. 1.00D-6) UN=1.00D-6
      CSSS(1,NBE)=CSSS(1,NBE)+1.0D00/UN
      CSSS(2,NBE)=CSSS(2,NBE)+SP(5,LS)/UN
      CSSS(3,NBE)=CSSS(3,NBE)+SP(5,LS)*UP/UN
      CSSS(4,NBE)=CSSS(4,NBE)+SP(5,LS)*(UP**2.0D00+UN**2.0D00+PV(3,N)**2.0D00)/UN
      IF (ISPR(1,LS) > 0) THEN
        CSSS(5,NBE)=CSSS(5,NBE)+PROT(N)/UN
        CSSS(6,NBE)=CSSS(6,NBE)+ISPR(1,LS)/UN
      END IF   
      END SUBROUTINE
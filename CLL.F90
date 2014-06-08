      SUBROUTINE CLL(N,LS,UN,UP,NBE,STEMP,VM)
      !
      !--CLL reflection of molecule N of species LS, normal and parallel velocity components UN,UP,
      !      from boundary edge NBE at temperature STEMP and most probable speed VM
      !
      USE MOLECS
      USE CELLINFO
      USE CONST
      USE GAS
      !
      IMPLICIT NONE
      !
      INTEGER :: N,LS,NBE
      REAL*8 :: UN,UP,STEMP,VM,VNI,UPI,WPI,ANG,VPI,ALPHAN,ALPHAT,ALPHAI,R,TH,UM,VN,VP,WP,OM,CTH,A
      REAL*8 :: RANF
      !
      VNI=UN/VM
      UPI=UP/VM
      WPI=PV(3,N)/VM
      ANG=DATAN2(WPI,UPI)
      VPI=DSQRT(UPI*UPI+WPI*WPI)
      !--VNI is the normalized incident normal vel. component (always +ve)
      !--VPI is the normalized incident tangential vel. comp. in int. plane
      !--ANG is the angle between the interaction plane and the x or y axis
      !
      !--first the normal component
      ALPHAN=SSSEG(5,LS,NBE)
      CALL RANDOM_NUMBER(RANF)
      R=DSQRT(-ALPHAN*DLOG(RANF))
      CALL RANDOM_NUMBER(RANF)
      TH=DPI*RANF
      UM=DSQRT(1.0D00-ALPHAN)*VNI  !average velocity
      VN=DSQRT(R*R+UM*UM+2.0D00*R*UM*DCOS(TH))
      !--VN is the normalized magnitude of the reflected normal vel. comp.
      !----from eqns (14.3)
      !
      !--then the tangential component
      ALPHAT=SSSEG(6,LS,NBE)*(2.0D00-SSSEG(6,LS,NBE))
      CALL RANDOM_NUMBER(RANF)
      R=DSQRT(-ALPHAT*DLOG(RANF))
      CALL RANDOM_NUMBER(RANF)
      TH=DPI*RANF
      UM=DSQRT(1.0D00-ALPHAT)*VPI
      VP=UM+R*DCOS(TH)
      WP=R*DSIN(TH)
      !--VP,WP are the normalized reflected tangential vel. components in and
      !----normal to the interaction plane, from eqns (14.4) and (14.5)
      UN=VN*VM
      UP=(VP*DCOS(ANG)-WP*DSIN(ANG))*VM
      PV(3,N)=(VP*DSIN(ANG)+WP*DCOS(ANG))*VM
      
      
      !Internal energy
      IF (ISPR(1,LS) > 0) THEN
         !--set CLL rotational energy by analogy with normal vel. component
         ALPHAI=SSSEG(7,LS,NBE)
         OM=DSQRT(PROT(N)*(1.0D00-ALPHAI)/(BOLTZ*STEMP))
         IF (ISPR(1,LS).EQ.2) THEN
            CALL RANDOM_NUMBER(RANF)
            R=DSQRT(-ALPHAI*DLOG(RANF))
            CALL RANDOM_NUMBER(RANF)
            CTH=DCOS(DPI*RANF)
         ELSE
            !--for polyatomic case, apply acceptance-rejection based on eqn (14.6)
            A=0.0D00
            DO WHILE (A < RANF)
               CALL RANDOM_NUMBER(RANF)
               R=4.0D00*RANF !For f=r^2exp(-r^2), f(r>4)~=0, 这里最小也要取3，可以更大
               A=R**2.0D00*DEXP(-R**2.0D00)/DEXP(-1.0D00)
               CALL RANDOM_NUMBER(RANF)
            END DO
            R=R*DSQRT(ALPHAI)
            CALL RANDOM_NUMBER(RANF)
            CTH=2.0D00*RANF-1.0D00  !COSINE uniform distributed
         END IF
         PROT(N)=BOLTZ*STEMP*(R*R+OM*OM+2.0D00*R*OM*CTH)
      END IF
      !
      RETURN
      !
      END SUBROUTINE CLL
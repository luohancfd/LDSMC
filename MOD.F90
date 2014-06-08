        MODULE NODEINFO
        !  Coordinate X,Y for node
        INTEGER,SAVE :: NNODES
        REAL*8,DIMENSION(:,:),ALLOCATABLE,SAVE :: NODES

        !--NODES(M,NNODES)-information of node
        !-----M=1--x coordinate of nodes
        !-----M=2--y coordinate of nodes
        !--NNODES-----Number of nodes
        END MODULE NODEINFO

        MODULE CELLINFO
        IMPLICIT NONE
        !  Number of cell
        INTEGER,SAVE :: NCELLS,NCCELLS,NBEDGES(5)
        INTEGER,DIMENSION(:,:),ALLOCATABLE,SAVE :: ICELL,ICCELL,IBEDGE1,IBEDGE2,IBEDGE3,IBEDGE4,IS
        REAL*8,SAVE :: MAXCELL(2),MINCELL(2)
        REAL*8,DIMENSION(:,:),ALLOCATABLE,SAVE :: CELL,CCELL,SSEG,ENTRYIN,ENTRYOUT,SURFACEVEC
        REAL*8,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: SSSEG

        !--NCELLS the number of sampling cells,
        !--ICELL(M,CELL) integer information about cell
        !-----M=1,2,3--the three node of the cell
        !-----M=4--the type of first edge
        !---------=0 Stream, No boundary
        !---------=1 Inflow
        !---------=2 Outflow
        !---------=3 SymmetryPlane
        !---------=4 Solid Wall
        !-----M=5--the type of second edge
        !---------=0 Stream, No boundary
        !---------=1 Inflow
        !---------=2 Outflow
        !---------=3 SymmetryPlane
        !---------=4 Solid Wall
        !-----M=6--the type of third
        !---------=0 Stream, No boundary
        !---------=1 Inflow
        !---------=2 Outflow
        !---------=3 SymmetryPlane
        !---------=4 Solid Wall
        !-----M=7--sum of boundary
        !-----M=8--number of collision cells preceding those in sampling cell N
        !-----M=9--the cell which has common edge(composed by node 1 and node 2) to cell N,0 for NULL
        !-----M=10--the cell which has common edge(composed by node 3 and node 2) to cell N,0 for NULL
        !-----M=11--the cell which has common edge(composed by node 1 and node 3) to cell N,0 for NULL
        !--CELL(M,cell) information of cell
        !-----M=1--the volume of cell
        !--MAXCELL(M) information of maxium value
        !-----M=1--max volume of the cell
        !-----M=2--max edge of the cell
        !--MINCELL(M) information of minium value
        !-----M=1--min volume of the cell
        !-----M=2--min edge of the cell
        !--NCCELLS the number of collision cells
        !--CCELL(M,N) information on collision cell N
        !    M=1 volume
        !    M=2 remainder in collision counting
        !    M=3 half the allowed time step
        !    M=4 maximum value of product of cross-section and relative velocity
        !    M=5 collision cell time
        !--ICCELL(M,N) integer information on collision cell N
        !    M=1 the (start address -1) in ICREF of molecules in collision cell N
        !    M=2 the number of molecules in collision cell N
        !    M=3 the sampling cell in which the collision cell lies
        !----------------following is mainly about boundary edge
        !--NBEDGES(M)---the number of different type of boundary edge
        !------M=1----the number of the sum
        !------M=2----the sum of inflow edges  ENTRYIN
        !------M=3----the sum of outflow edges  ENTRYOUT
        !------M=4----the sum of symmetrical surface edges
        !------M=5----the sum of solid surface edges  SSEG
        !--IS(JJ,N) the boundary edge code of JJ edge in cell N, 0 for null
        !---------for different type of boundary return different type of boundary code
        !--IBEDGE1(M,N)----N=1,NBEDGES(2)  for INFLOW
        !-------M=1---the code of cell which contains the solid surface edge
        !-------M=2---the edge code of the cell
        !--IBEDGE2(M,N)----N=1,NBEDGES(3)  for OUTFLOW
        !-------M=1---the code of cell which contains the solid surface edge
        !-------M=2---the edge code of the cell
        !--IBEDGE3(M,N)----N=1,NBEDGES(4)  for SYMETRICAL BOUNDARY
        !-------M=1---the code of cell which contains the solid surface edge
        !-------M=2---the edge code of the cell
        !--IBEDGE4(M,N)----N=1,NBEDGES(5)  for SOLID SURFACE
        !-------M=1---the code of cell which contains the solid surface edge
        !-------M=2---the edge code of the cell
        !--ENTRYIN(M,N)----N=1,NBEDGES(2) for INFLOW boundary
        !-------(ENTRYIN(1,N),ENTRYIN(2,N))---the normal vector of the boundary
        !--ENTRYOUT(M,N)---N=1,NBEDGES(3) for OUTFLOW boundary
        !-------(ENTRYOUT(1,N),ENTRYOUT(2,N))---the normal vector of the boundary,the direction is on the right of the edge
        !--SURFACEVEC(2,N)-N=1,NBEDGES(5)-----the normal vector of the surface,the direction is on the right of the edge
        !--SSEG(M,N) data on solid surface edge N  Solid Surface EdGe=SSEG
        !---------M=1 temperature
        !---------M=2 emissivity (needed for adiabatic element only)
        !---------M=3 in-plane velocity of surface
        !---------M=4 normal-to-plane velocity of surface
        !---------M=5 the calculated temperature (adiabatic cases only)
        !--SSSEG(M,L,N) species L dependent data on solid surface edge N
        !---------M=1 fraction of species  adsorbed
        !---------M=2 initially outgas mass flux of species L
        !---------    finally the number entering per unit time
        !---------M=3 remainder
        !---------M=4 -1.,1. for diffuse, CLL reflection
        !---------M=5 normal energy accommodation coefficient (CLL only) alpha_n
        !---------ATTENTION : alpha is energy accommodation coefficient, thegma is velocity(momentum) accommodation coefficient
        !---------alpha=sigma*(2-sigma)
        !---------1-sigma=sqrt(1-alpha)
        !---------M=6 tangential momentum accommodation coefficient (CLL only) thegma_t
        !---------M=7 rotational accommodation coefficient
        !---------M=8 vibrational accommodation coefficient
        !---------M=9 fraction of specular reflection
        END MODULE

        MODULE BOUNDINFO
        !!ATTENTION!!this module is only used in the preprocessing,output
        IMPLICIT NONE
        INTEGER :: NBOUND
        INTEGER,SAVE :: NWALL
        INTEGER,ALLOCATABLE :: IWBOUND(:)
        INTEGER,ALLOCATABLE,SAVE :: IWALL(:,:)
        !------NBOUND------number of BC
        !------NWALL-------number of solid surface
        !------IWBOUND-----the wall number of N BC
        !------IWALL(0,N)-------the total number of edge in N wall
        !------IWALL(1~*,N)-----the edge code of N wall
        TYPE :: BC
            integer :: type_bound             !type number of bound
            character*20 :: name       !name of boundary condition
            integer :: type_points       !type of grids composing BC,1 for point range,2 for point range
            integer :: num_points        !number of points in BC
            !----------------variables below are useless here,just for function storage-----------
            integer :: index_normal    !index number for patch normal,a vector,USELESS=0
            integer :: size_normal     !size of NormalList,USELESS=0
            integer :: datatype_normal !storage type of NormalList,RealSingle or RealDouble or NULL=0,USELESS
            integer :: ndataset
            !----------------variables above useless here------------------------------------------
            integer,allocatable :: points(:) !PointsList compose BC
        end type BC
        type(BC),allocatable :: bound(:)
        end module

        MODULE CALC
        IMPLICIT NONE
        ! Initial Condition for flow
        INTEGER,SAVE :: MOLSC,IMTS,NNC,ISF,ICENS
        REAL*8,SAVE :: FNUM,DTM,CPDTM,TPDTM,SAMPRAT,OUTRAT,TSAMP,TOUT,DTSAMP,DTOUT,FTIME,TLIM,TOTMOV, &
            TOTCOL,TOTCOLI,TOTMOVI,ENTMASS
        REAL*8,DIMENSION(:,:),ALLOCATABLE,SAVE :: SSSEG_INT,TCOL,SSEG_INT
        !INTEGER(KIND=8) MCOMPTIM,MCOMPTIMI,MCOMPTIMF
        !--ISF 0,1 for steady, unsteady flow sampling
        !--ICENS 1 for pressure tensor anisotropy as the non-equilibrium parameter
        !        2 for Tx/Ty-1 as the non-equilibrium parameter
        !        3 for Trot/Ttrans-1 as the non-equilibrium parameter
        !        4 for Tvib/Ttrans-1 as the non-equilibrium parameter
        !--MOLSC the numbet of molecules per sampling cell
        !--IMTS 0 for uniform move time steps, 1 for time steps that vary over the cells
        !--NNC 0 to select collision partner randomly from collision cell, 1 for nearest neighbor collisions
        !--SAMPRAT the number of time steps between samplings
        !--SAMPRAT the number of time steps between samplings
        !--FNUM the number of real molecules represented by each simulated molecule
        !--CPDTM the maximum number of collisions per time step (standard 0.2)
        !--TPDTM the maximum number of sampling cell transit times of the flow per time step 这个大小的设置还不确定
        !--RANF random value
        !--FTIME the flow time
        !--TLIM the time at which the calculation stops
        !--TOTMOV total molecule moves
        !--TOTCOL total collisions
        !--TCOL species dependent collision counter
        !--SSSEG_INT(M,L)-----this is used to define solid boundary properties,
        !                     the value will be copied to SSSEG since all solid boundary edge all the same
        !                     这个未来还可以得到扩展
        !                     L is the specie of the molecule
        !---------M=1 fraction of species  adsorbed
        !---------M=2 initially outgas mass flux of species L
        !---------    finally the number entering per unit time
        !---------M=3 remainder
        !---------M=4 -1.,1. for diffuse, CLL reflection
        !---------M=5 normal energy accommodation coefficient (CLL only)
        !---------M=6 tangential momentum accommodation coefficient (CLL only)
        !---------M=7 rotational accommodation coefficient
        !---------M=8 vibrational accommodation coefficient
        !---------M=9 fraction of specular reflection
        !--SSEG_INT(M,NWALL) -----------this data is used to define solid boundary properties,
        !                     values will be copied to SSEG since they have the same properties
        !                     Could be extended in the future
        !                     ATTENTION! The diffuse adiabatic option is NOT applicable now, moving boundary is also NOT applicable
        !                                thus M=1 for negative, M=2, M=3, M=5 is just reserved variable but USELESS
        !---------M=1 temperature, diffuse if positive
        !---------                 diffuse adiabatic insulated if negative (abs is minimum temp)
        !---------                 diffuse adiabatic perfectly conducting
        !---------M=2 emissivity (needed for adiabatic element only)
        !---------M=3 in-plane velocity of surface
        !---------M=4 normal-to-plane velocity of surface
        !---------M=5 the calculated temperature (adiabatic cases only)


        END MODULE

        MODULE GAS
        !--declares the variables associated with the molecular species and the stream definition
        IMPLICIT NONE
        REAL(KIND=8),SAVE :: RMAS,CXSS,RGFS,VMPM,FDEN,FPR,FMA,FP,CTM
        REAL(KIND=8),SAVE :: FND,FTMP,FVTMP,VFX,VFY
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:),SAVE :: CR,FSP,VMP
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:),SAVE :: SP,SPR
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:),SAVE :: SPM,SPVM,ENTR
        INTEGER,SAVE :: MSP,MMVM,MMRM,MNRE,MNSR,MTBP,GASCODE,MMEX,MEX
        INTEGER, ALLOCATABLE, DIMENSION(:),SAVE :: ISPV
        INTEGER, ALLOCATABLE, DIMENSION(:,:),SAVE :: ISPR
        INTEGER, ALLOCATABLE, DIMENSION(:,:,:),SAVE :: ISPVM
        !
        !--MSP the number of molecular species
        !--MMVM the maximum number of vibrational modes of any species
        !--MEX number of exchange or chain reactions
        !--MMEX the maximum number of exchange reactions involving the same precollision pair of molecules
        !--MMRM 0 if gas is completely monatomic, 1 if some species have rotation
        !--MNRE the number of gas phase chemical reactions (old scheme)
        !--MNSR the number oF surface reactions
        !--MTBP the number of third body tables
        !--SP(1,L) the reference diameter of species L
        !--SP(2,L) the reference temperature of species L
        !--SP(3,L) the viscosity-temperature power law of species L
        !--SP(4,L) the reciprocal of the VSS scattering parameter
        !--SP(5,L) molecular mass of species L
        !--ISPR(1,L) number of rotational degrees of freedom of species L
        !--ISPR(2,L) 0,1 for constant, polynomial rotational relaxation collision number
        !--SPR(1,L) constant rotational relaxation collision number of species L
        !           or the constant in a second order polynomial in temperature
        !--SPR(2,L) the coefficient of temperature in the polynomial
        !--SPR(3,L) the coefficient of temperature squared in the polynomial
        !--SPM(1,L,M) the reduced mass for species L,M
        !--SPM(2,L,M) the reference collision cross-section for species L,M
        !--SPM(3,L,M) the mean value of the viscosity-temperature power law
        !--SPM(4,L,M) the reference diameter for L,M collisions
        !--SPM(5,L,M) the reference temperature for species L,M
        !--SPM(6,L,M) reciprocal of the gamma function of (5/2-w) for species L,M
        !--SPM(7,L,M) rotational relaxation collision number for species L,M, or const in polynomial
        !--SPM(8,L,M) coefficient of temperature in polynomial for rotational collision number
        !--SPM(9,L,M) coefficient of temperature squared in this polynomial
        !--SPM(10,L,M) reciprocal of VSS scattering parameter
        !--ISPV(L) the number of vibrational modes
        !--SPVM(1,K,L) the characteristic vibrational temperature
        !--SPVM(2,K,L) constant Zv, or reference Zv for mode K
        !--SPVM(3,K,L) -1. for constant Zv, or reference temperature
        !--SPVM(4,K,L) the characteristic dissociation temperature
        !--ISPVM(1,K,L) the species code of the first dissociation product
        !--ISPVM(2,K,L) the species code of the second dissociation product
        !--THBP(N,L) third body efficiency in table N of species L
        !--RMAS reduced mass for single species case
        !--CXSS reference cross-section for single species case
        !--RGFS reciprocal of gamma function for single species case
        !--for the following, J=1 for the reference gas and/or the minimum x boundary, J=2 for the secondary sream at maximum x boundary
        !--FND stream or reference gas number density
        !--FTMP stream temperature
        !--FVTMP the vibrational temperature in the freestream
        !--VFX  the x velocity components of the stream
        !--VFY the y velocity component in the stream
        !--FSP(N)) fraction of species N in the stream
        !--FMA stream Mach number
        !--VMP(N) most probable molecular velocity of species N at FTMP
        !--VMPM the maximum value of VMP in stream 1
        !--ENTR(M,L,K) entry/removal information for species L at K edge
        !---------L=1,MSP;K=1,NBEDGES(2)+NBEDGES(3)
        !    M=1 number of molecule enter per time
        !    M=2 remainder
        !    M=3 speed ratio
        !    M=4 first constant
        !    M=5 second constant
        !    M=6 the maxinum normal velocity component in the removal zone (> XREM)
        !--CTM approximate mean collision time in stream (exact for simple gas)
        !--FP approximate mean free path
        !--FPR stream 1 pressure
        !--FMA stream 1 Mach number
        !--RMAS reduced mass for single species case
        !--CXSS reference cross-section for single species case
        !--RGFS reciprocal of gamma function for single species case
        !--CR(L) collision rate of species L

        END MODULE GAS

        MODULE MOLECS
        !
        !--declares the variables associated with the molecules
        !
        IMPLICIT NONE
        !
        INTEGER, ALLOCATABLE, DIMENSION(:),SAVE :: IPCCELL,IPSP,ICREF,IPCP
        INTEGER, ALLOCATABLE, DIMENSION(:,:),SAVE :: IPVIB
        INTEGER ,SAVE:: NM,MNM
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:),SAVE :: PV
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:),SAVE :: PX,PY,PTIM,PROT
        !
        !--PX(N) position (x coordinate) of molecule N
        !--PY(N) position (y coordinate) of molecule N
        !--PTIM(N) molecule time
        !--PROT(N) rotational energy
        !--PV(1-3,N) u,v,w velocity components
        !--IPSP(N) the molecular species
        !--IPCCELL(N) the collision cell number
        !--ICREF the cross-reference array (molecule numbers in order of collision cells)
        !------给一个基于collision cell的位置编号，返回一个分子编号
        !--IPCP(N) the code number of the last collision partner of molecule
        !--IPVIB(K,N) level of vibrational mode K of molecule N
        !---IPTSC(I,M)-----molecule to transient cell USE FOR transient cell method in COLLISION SUBROUTINE
        !-------I=1----x index of molecule M
        !-------I=2----y index of molecule M
        !-------I=3----the index of molecule in the transient cell
        !--NM number of molecules
        !--MNM the maximum number of molecules
        !
        END MODULE MOLECS

        MODULE CONST
        IMPLICIT NONE
        REAL*8,PARAMETER :: PI=3.1415926535897932D00
        REAL*8,PARAMETER :: DPI=6.283185307179586D00
        REAL*8,PARAMETER :: SPI=1.772453850905516D00
        REAL*8,PARAMETER :: BOLTZ=1.380658D-23
        END MODULE

        MODULE OUTPUT
        IMPLICIT NONE
        REAL*8,DIMENSION(:,:),ALLOCATABLE,SAVE :: VAR
        REAL*8,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: VARSP,VARS,DF

        !--VAR(M,N) the flowfield properties in samping cell N
        !--M=1 the x coordinate
        !--M=2 the y coordinate
        !--M=3 neutral and ion number density
        !--M=4 density
        !--M=5 u velocity component
        !--M=6 v velocity component
        !--M=7 w velocity component
        !--M=8 translational temperature
        !--M=9 rotational temperature
        !--M=10 vibrational temperature
        !--M=11 temperature
        !--M=12 Mach number
        !--M=13 molecules per cell
        !--M=14 mean collision time
        !--M=15 mean free path
        !--M=16 ratio (mean collisional separation) / (mean free path) or VAR(16,N)
        !--M=17 flow speed
        !--M=18 flow angle
        !--M=19 scalar pressure nkT
        !--M=20 x component of translational temperature TTX
        !--M=21 y component of translational temperature TTY
        !--M=22 z component of translational temperature TTZ
        !--M=23 nonequilibrium parameter
        !--     if ICENS=1, rms pressure tensor asymmetry
        !--     if ICENS=2, Tx/Ty-1
        !--     if ICENS=3, Trot/Ttrans
        !--     if ICENS=4, Tvib/Ttrans-1
        !
        !--VARSP(M,N,L) the flowfield properties for species L in cell N
        !--M=1 the fraction
        !--M=2 the temperature component in the x direction
        !--M=3 the temperature component in the y direction
        !--M=4 the temperature component in the z direction
        !--M=5 the translational temperature
        !--M=6 the rotational temperature
        !--M=7 the vibrational temperature
        !--M=8 the temperature
        !--M=9 the x component of the diffusion velocity
        !--M=10 the y component of the diffusion velocity
        !--M=11 the z component of the diffusion velocity

        !--VARS(N,M,L) surface property N on interval M of boundary L


        !--N=1 the incident sample
        !--N=2 the reflected sample
        !
        !--N=3 the incident number flux
        !--N=4 the reflected number flux
        !
        !--N=5 the incident pressure
        !--N=6 the reflected pressure
        !
        !--N=7 the incident parallel shear tress
        !--N=8 the reflected parallel shear stress
        !
        !--N=9 the incident normal-to-plane shear stress
        !--N=10 the reflected normal shear stress
        !
        !--N=11 the incident translational heat flux
        !--N=12 the reflected translational heat flux
        !
        !--N=13 the incident rotational heat flux
        !--N=14 the reflected rotational heat flux
        !
        !--N=15 the incident vibrational heat flux
        !--N=16 the reflected vibrational heat flux
        !
        !--N=17 slip velocity
        !--N=18 temperature slip
        !--N=19 rotational temperature slip
        !--N=20 the net pressure
        !--N=21 the net parallel in-plane shear
        !--N=22 the net parallel normal-to-plane shear
        !--N=23 the net translational energy flux
        !--N=24 the net rotational heat flux
        !--N=25 the net vibrational heat flux
        !--N=26 total incident heat transfer
        !--N=27 total reflected heat transfer
        !--N=28 net heat transfer
        !--N=29 surface temperature
        !--DF(N,M,L) the average d.o.f. of vibrational mode M of species L in sampling cell N

        END MODULE


        MODULE SAMPLES
        IMPLICIT NONE
        INTEGER,SAVE :: NSAMP
        REAL,SAVE :: TISAMP
        REAL*8,DIMENSION(:),ALLOCATABLE,SAVE :: EMOLS
        REAL*8,DIMENSION(:,:),ALLOCATABLE,SAVE :: CSSS
        REAL*8,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: CS
        REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE,SAVE :: CSS
        REAL*8,DIMENSION(:),ALLOCATABLE,SAVE :: COLLS,CLSEP

        !--EMOLS(N) sum of molecules in collision cells that are set to be in sampling cell N
        !--NSAMP-----sampling time
        !--TISAMP the time at which the sampling was last reset
        !--CS(1,N,L) sampled number of species L in cell N
        !--CS(2,N,L), CS(3,N,L), CS(4,N,L) sampled sum of u, v, w
        !--CS(5,N,L), CS(6,N,L), CS(7,N,L) sampled sum of u*u, v*v, w*w
        !--CS(8,N,L) sampled sum of rotational energy of species L in cell N
        !--CS(8+K,N,L) sampled sum of vibrational energy of species L in cell N
        !              K is the mode
        !---CSS(N,J,K,L)---sample information for boundary edge
        !----------J-------then code of boundary edge
        !----------K-------the specie of interacted molecule
        !----------L-------=1 for incident molecule, =2 for reflected molecule
        !----------N=1-----number sum of molecules of species K interacted with edge J
        !----------N=2-----normal momentum sum to edge
        !----------N=3-----parallel (in plane) momentum sum to edge
        !----------N=4-----parallel (normal to plane) momentum sum to edge
        !----------N=5-----translational energy sum to edge
        !----------N=6-----rotational energy sum to edge
        !----------N=7-----vibrational energy sum to edge
        !---CSSS(N,J)-------目前没搞明白是啥
        !----------N=1---sum (over incident AND reflected molecules) of 1/normal vel. component
        !----------N=2---similar sum of molecular mass / normal vel. component
        !----------N=3---similar sum of molecular mass * parallel vel. component / normal vel. component
        !----------N=4---similar sum of molecular mass * speed squared / normal vel. component
        !----------N=5---similar sum of rotational energy / normal vel. component
        !----------N=6---similar sum of rotational degrees of freedom /normal velocity component
        !---COLLS(N) total number of collisions in sampling cell N
        !---CLSEP(N) sum of collision pair separation in cell N
        END MODULE


! * ------------------------------------------------------------------------- *
! * CGNS - CFD General Notation System (http://www.cgns.org)                  *
! * CGNS/MLL - Mid-Level Library header file                                  *
! * Please see cgnsconfig.h file for this local installation configuration    *
! * ------------------------------------------------------------------------- *
!
! * ------------------------------------------------------------------------- *
!
!  This software is provided 'as-is', without any express or implied warranty.
!  In no event will the authors be held liable for any damages arising from
!  the use of this software.
!
!  Permission is granted to anyone to use this software for any purpose,
!  including commercial applications, and to alter it and redistribute it
!  freely, subject to the following restrictions:
!
!  1. The origin of this software must not be misrepresented; you must not
!     claim that you wrote the original software. If you use this software
!     in a product, an acknowledgment in the product documentation would be
!     appreciated but is not required.
!
!  2. Altered source versions must be plainly marked as such, and must not
!     be misrepresented as being the original software.
!
!  3. This notice may not be removed or altered from any source distribution.
!
! * ------------------------------------------------------------------------- *
!
! Fortran version of cgnslib.h
    integer :: CG_BUILD_64BIT
    parameter (CG_BUILD_64BIT = 1)
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      modes for cgns file                                            *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    integer*4 :: CG_MODE_READ, CG_MODE_WRITE, CG_MODE_MODIFY
    parameter (CG_MODE_READ   = 0)
    parameter (CG_MODE_WRITE  = 1)
    parameter (CG_MODE_MODIFY = 2)
!* legacy code support
    integer*4 :: MODE_READ, MODE_WRITE, MODE_MODIFY
    parameter (MODE_READ   = 0)
    parameter (MODE_WRITE  = 1)
    parameter (MODE_MODIFY = 2)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      file types                                                     *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    integer*4 :: CG_FILE_NONE, CG_FILE_ADF, CG_FILE_HDF5
    integer*4 :: CG_FILE_ADF2
    parameter (CG_FILE_NONE = 0)
    parameter (CG_FILE_ADF  = 1)
    parameter (CG_FILE_HDF5 = 2)
    parameter (CG_FILE_ADF2 = 3)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      some error code                                                *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    integer*4 :: CG_OK, CG_ERROR, CG_NODE_NOT_FOUND
    integer*4 :: CG_INCORRECT_PATH, CG_NO_INDEX_DIM
    parameter (CG_OK             = 0)
    parameter (CG_ERROR          = 1)
    parameter (CG_NODE_NOT_FOUND = 2)
    parameter (CG_INCORRECT_PATH = 3)
    parameter (CG_NO_INDEX_DIM   = 4)
!* legacy code support
    integer*4 :: ALL_OK, ERROR, NODE_NOT_FOUND, INCORRECT_PATH
    parameter (ALL_OK         = 0)
    parameter (ERROR          = 1)
    parameter (NODE_NOT_FOUND = 2)
    parameter (INCORRECT_PATH = 3)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Dimensional Units                                              *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    integer*4 :: CG_UserDefined, CG_Null
    parameter (CG_Null = 0)
    parameter (CG_UserDefined = 1)
!* legacy code support
    integer*4 :: Null, UserDefined
    parameter (Null = 0)
    parameter (UserDefined = 1)
    
    integer*4 :: Kilogram, Gram, Slug, PoundMass
    character(32) :: MassUnitsName(0:5)
    parameter (Kilogram  = 2)
    parameter (Gram      = 3)
    parameter (Slug      = 4)
    parameter (PoundMass = 5)
    
    integer*4 :: Meter, Centimeter, Millimeter
    integer*4 :: Foot, Inch
    character(32) :: LengthUnitsName(0:6)
    parameter (Meter      = 2)
    parameter (Centimeter = 3)
    parameter (Millimeter = 4)
    parameter (Foot       = 5)
    parameter (Inch       = 6)
    
    integer*4 :: Second
    character(32) :: TimeUnitsName(0:2)
    parameter (Second = 2)
    
    integer*4 :: Kelvin, Celsius, Rankine, Fahrenheit
    character(32) :: TemperatureUnitsName(0:5)
    parameter (Kelvin     = 2)
    parameter (Celsius    = 3)
    parameter (Rankine    = 4)
    parameter (Fahrenheit = 5)
    
    integer*4 :: Degree, Radian
    character(32) :: AngleUnitsName(0:3)
    parameter (Degree = 2)
    parameter (Radian = 3)
    
    integer*4 :: Ampere, Abampere, Statampere, Edison, auCurrent
    character(32) :: ElectricCurrentUnitsName(0:6)
    parameter (Ampere     = 2)
    parameter (Abampere   = 3)
    parameter (Statampere = 4)
    parameter (Edison     = 5)
    parameter (auCurrent  = 6)
    
    integer*4 :: Mole, Entities, StandardCubicFoot, StandardCubicMeter
    character(32) :: SubstanceAmountUnitsName(0:5)
    parameter (Mole               = 2)
    parameter (Entities           = 3)
    parameter (StandardCubicFoot  = 4)
    parameter (StandardCubicMeter = 5)
    
    integer*4 :: Candela, Candle, Carcel, Hefner, Violle
    character(32) :: LuminousIntensityUnitsName(0:6)
    parameter (Candela = 2)
    parameter (Candle  = 3)
    parameter (Carcel  = 4)
    parameter (Hefner  = 5)
    parameter (Violle  = 6)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Data Class                                                     *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    integer*4 :: Dimensional, NormalizedByDimensional
    integer*4 :: NormalizedByUnknownDimensional
    integer*4 :: NondimensionalParameter, DimensionlessConstant
    character(32) :: DataClassName(0:6)
    parameter (Dimensional                    = 2)
    parameter (NormalizedByDimensional        = 3)
    parameter (NormalizedByUnknownDimensional = 4)
    parameter (NondimensionalParameter        = 5)
    parameter (DimensionlessConstant          = 6)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Grid Location                                                  *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: Vertex, CellCenter, FaceCenter
    integer*4 :: IFaceCenter, JFaceCenter, KFaceCenter, EdgeCenter
    character(32) :: GridLocationName(0:8)
    parameter (Vertex      = 2)
    parameter (CellCenter  = 3)
    parameter (FaceCenter  = 4)
    parameter (IFaceCenter = 5)
    parameter (JFaceCenter = 6)
    parameter (KFaceCenter = 7)
    parameter (EdgeCenter  = 8)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Grid Connectivity Types                                        *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: Overset, Abutting, Abutting1to1
    character(32) :: GridConnectivityTypeName(0:4)
    parameter (Overset      = 2)
    parameter (Abutting     = 3)
    parameter (Abutting1to1 = 4)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Point Set Types                                                *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: PointList, PointListDonor
    integer*4 :: PointRange, PointRangeDonor
    integer*4 :: ElementRange, ElementList, CellListDonor
    character(32) :: PointSetTypeName(0:8)
    parameter (PointList       = 2)
    parameter (PointListDonor  = 3)
    parameter (PointRange      = 4)
    parameter (PointRangeDonor = 5)
    parameter (ElementRange    = 6)
    parameter (ElementList     = 7)
    parameter (CellListDonor   = 8)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Governing Equations and Physical Models Types                  *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: FullPotential, Euler
    integer*4 :: NSLaminar, NSTurbulent
    integer*4 :: NSLaminarIncompressible
    integer*4 :: NSTurbulentIncompressible
    character(32) :: GoverningEquationsTypeName(0:7)
    parameter (FullPotential             = 2)
    parameter (Euler                     = 3)
    parameter (NSLaminar                 = 4)
    parameter (NSTurbulent               = 5)
    parameter (NSLaminarIncompressible   = 6)
    parameter (NSTurbulentIncompressible = 7)
    
!** Any model type will accept both ModelTypeNull and ModelTypeUserDefined.
!** The following models will accept these values as vaild...
!**
!** GasModel_t: Ideal, VanderWaals, CaloricallyPerfect, ThermallyPerfect,
!**    ConstantDensity, RedlichKwong
!**
!** ViscosityModel_t: Constant, PowerLaw, SutherlandLaw
!**
!** ThermalConductivityModel_t: PowerLaw, SutherlandLaw, ConstantPrandtl
!**
!** TurbulenceModel_t: Algebraic_BaldwinLomax, Algebraic_CebeciSmith,
!**    HalfEquation_JohnsonKing, OneEquation_BaldwinBarth,
!**    OneEquation_SpalartAllmaras, TwoEquation_JonesLaunder,
!**    TwoEquation_MenterSST,TwoEquation_Wilcox
!**
!** TurbulenceClosure_t: EddyViscosity, ReynoldsStress,
!**    ReynoldsStressAlgebraic
!**
!** ThermalRelaxationModel_t: Frozen, ThermalEquilib, ThermalNonequilib
!**
!** ChemicalKineticsModel_t: Frozen, ChemicalEquilibCurveFit,
!**    ChemicalEquilibMinimization, ChemicalNonequilib
!**
!** EMElectricFieldModel_t: Voltage, Interpolated, Constant, Frozen
!**
!** EMMagneticFieldModel_t: Interpolated, Constant, Frozen
!**
!** EMConductivityModel_t: Constant, Frozen, Equilibrium_LinRessler,
!**                        Chemistry_LinRessler
    
    integer*4 :: Ideal, VanderWaals
    integer*4 :: Constant
    integer*4 :: PowerLaw, SutherlandLaw
    integer*4 :: ConstantPrandtl
    integer*4 :: EddyViscosity, ReynoldsStress
    integer*4 :: ReynoldsStressAlgebraic
    integer*4 :: Algebraic_BaldwinLomax, Algebraic_CebeciSmith
    integer*4 :: HalfEquation_JohnsonKing, OneEquation_BaldwinBarth
    integer*4 :: OneEquation_SpalartAllmaras
    integer*4 :: TwoEquation_JonesLaunder
    integer*4 :: TwoEquation_MenterSST, TwoEquation_Wilcox
    integer*4 :: CaloricallyPerfect, ThermallyPerfect
    integer*4 :: ConstantDensity, RedlichKwong
    integer*4 :: Frozen, ThermalEquilib, ThermalNonequilib
    integer*4 :: ChemicalEquilibCurveFit
    integer*4 :: ChemicalEquilibMinimization
    integer*4 :: ChemicalNonequilib
    integer*4 :: EMElectricField, EMMagneticField, Voltage
    integer*4 :: Interpolated
    integer*4 :: EMConductivity, Equilibrium_LinRessler
    integer*4 :: Chemistry_LinRessler
    character(32) :: ModelTypeName(0:35)
    
    parameter (Ideal                       = 2)
    parameter (VanderWaals                 = 3)
    parameter (Constant                    = 4)
    parameter (PowerLaw                    = 5)
    parameter (SutherlandLaw               = 6)
    parameter (ConstantPrandtl             = 7)
    parameter (EddyViscosity               = 8)
    parameter (ReynoldsStress              = 9)
    parameter (ReynoldsStressAlgebraic     = 10)
    parameter (Algebraic_BaldwinLomax      = 11)
    parameter (Algebraic_CebeciSmith       = 12)
    parameter (HalfEquation_JohnsonKing    = 13)
    parameter (OneEquation_BaldwinBarth    = 14)
    parameter (OneEquation_SpalartAllmaras = 15)
    parameter (TwoEquation_JonesLaunder    = 16)
    parameter (TwoEquation_MenterSST       = 17)
    parameter (TwoEquation_Wilcox          = 18)
    parameter (CaloricallyPerfect          = 19)
    parameter (ThermallyPerfect            = 20)
    parameter (ConstantDensity             = 21)
    parameter (RedlichKwong                = 22)
    parameter (Frozen                      = 23)
    parameter (ThermalEquilib              = 24)
    parameter (ThermalNonequilib           = 25)
    parameter (ChemicalEquilibCurveFit     = 26)
    parameter (ChemicalEquilibMinimization = 27)
    parameter (ChemicalNonequilib          = 28)
    parameter (EMElectricField             = 29)
    parameter (EMMagneticField             = 30)
    parameter (EMConductivity              = 31)
    parameter (Voltage                     = 32)
    parameter (Interpolated                = 33)
    parameter (Equilibrium_LinRessler      = 34)
    parameter (Chemistry_LinRessler        = 35)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Boundary Condition Types                                       *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: BCAxisymmetricWedge, BCDegenerateLine
    integer*4 :: BCDegeneratePoint
    integer*4 :: BCDirichlet, BCExtrapolate, BCFarfield, BCGeneral
    integer*4 :: BCInflow, BCInflowSubsonic,  BCInflowSupersonic
    integer*4 :: BCNeumann
    integer*4 :: BCOutflow, BCOutflowSubsonic, BCOutflowSupersonic
    integer*4 :: BCSymmetryPlane, BCSymmetryPolar
    integer*4 :: BCTunnelInflow, BCTunnelOutflow
    integer*4 :: BCWall, BCWallInviscid, BCWallViscous
    integer*4 :: BCWallViscousHeatFlux, BCWallViscousIsothermal
    integer*4 :: FamilySpecified
    character(32) :: BCTypeName(0:25)
    parameter (BCAxisymmetricWedge     = 2)
    parameter (BCDegenerateLine        = 3)
    parameter (BCDegeneratePoint       = 4)
    parameter (BCDirichlet             = 5)
    parameter (BCExtrapolate           = 6)
    parameter (BCFarfield              = 7)
    parameter (BCGeneral               = 8)
    parameter (BCInflow                = 9)
    parameter (BCInflowSubsonic        = 10)
    parameter (BCInflowSupersonic      = 11)
    parameter (BCNeumann               = 12)
    parameter (BCOutflow               = 13)
    parameter (BCOutflowSubsonic       = 14)
    parameter (BCOutflowSupersonic     = 15)
    parameter (BCSymmetryPlane         = 16)
    parameter (BCSymmetryPolar         = 17)
    parameter (BCTunnelInflow          = 18)
    parameter (BCTunnelOutflow         = 19)
    parameter (BCWall                  = 20)
    parameter (BCWallInviscid          = 21)
    parameter (BCWallViscous           = 22)
    parameter (BCWallViscousHeatFlux   = 23)
    parameter (BCWallViscousIsothermal = 24)
    parameter (FamilySpecified         = 25)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Data types                                                     *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: Integer, RealSingle, RealDouble, Character
    integer*4 :: LongInteger
    character(32) :: DataTypeName(0:6)
    parameter (Integer     = 2)
    parameter (RealSingle  = 3)
    parameter (RealDouble  = 4)
    parameter (Character   = 5)
    parameter (LongInteger = 6)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      BCData_t types                                                 *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: Dirichlet, Neumann
    character(32) :: BCDataTypeName(0:3)
    parameter (Dirichlet = 2)
    parameter (Neumann   = 3)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Element types                                                  *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: NODE, BAR_2, BAR_3, TRI_3, TRI_6
    integer*4 :: QUAD_4, QUAD_8, QUAD_9
    integer*4 :: TETRA_4, TETRA_10, PYRA_5, PYRA_14
    integer*4 :: PENTA_6, PENTA_15, PENTA_18
    integer*4 :: HEXA_8, HEXA_20, HEXA_27
    integer*4 :: MIXED, PYRA_13, NGON_n, NFACE_n
    integer*4 :: BAR_4, TRI_9, TRI_10
    integer*4 :: QUAD_12, QUAD_16
    integer*4 :: TETRA_16, TETRA_20
    integer*4 :: PYRA_21, PYRA_29, PYRA_30
    integer*4 :: PENTA_24, PENTA_38, PENTA_40
    integer*4 :: HEXA_32, HEXA_56, HEXA_64
    character(32) :: ElementTypeName(0:39)
    parameter (NODE     =  2)
    parameter (BAR_2    =  3)
    parameter (BAR_3    =  4)
    parameter (TRI_3    =  5)
    parameter (TRI_6    =  6)
    parameter (QUAD_4   =  7)
    parameter (QUAD_8   =  8)
    parameter (QUAD_9   =  9)
    parameter (TETRA_4  = 10)
    parameter (TETRA_10 = 11)
    parameter (PYRA_5   = 12)
    parameter (PYRA_14  = 13)
    parameter (PENTA_6  = 14)
    parameter (PENTA_15 = 15)
    parameter (PENTA_18 = 16)
    parameter (HEXA_8   = 17)
    parameter (HEXA_20  = 18)
    parameter (HEXA_27  = 19)
    parameter (MIXED    = 20)
    parameter (PYRA_13  = 21)
    parameter (NGON_n   = 22)
    parameter (NFACE_n  = 23)
    parameter (BAR_4    = 24)
    parameter (TRI_9    = 25)
    parameter (TRI_10   = 26)
    parameter (QUAD_12  = 27)
    parameter (QUAD_16  = 28)
    parameter (TETRA_16 = 29)
    parameter (TETRA_20 = 30)
    parameter (PYRA_21  = 31)
    parameter (PYRA_29  = 32)
    parameter (PYRA_30  = 33)
    parameter (PENTA_24 = 34)
    parameter (PENTA_38 = 35)
    parameter (PENTA_40 = 36)
    parameter (HEXA_32  = 37)
    parameter (HEXA_56  = 38)
    parameter (HEXA_64  = 39)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Zone types                                                     *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: Structured, Unstructured
    character(32) :: ZoneTypeName(0:3)
    parameter (Structured   =  2)
    parameter (Unstructured =  3)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Rigid Grid Motion types                                        *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: ConstantRate, VariableRate
    character(32) :: RigidGridMotionTypeName(0:3)
    parameter (ConstantRate = 2)
    parameter (VariableRate = 3)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Arbitrary Grid Motion types                                    *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: NonDeformingGrid, DeformingGrid
    character(32) :: ArbitraryGridMotionTypeName(0:3)
    parameter (NonDeformingGrid = 2)
    parameter (DeformingGrid = 3)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Simulation type                                                *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: TimeAccurate, NonTimeAccurate
    character(32) :: SimulationTypeName(0:3)
    parameter (TimeAccurate = 2)
    parameter (NonTimeAccurate = 3)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      BC Property types                                              *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: Generic
    character(32) :: WallFunctionTypeName(0:2)
    parameter (Generic = 2)
    
    integer*4 :: BleedArea, CaptureArea
    character(32) :: AreaTypeName(0:3)
    parameter (BleedArea = 2)
    parameter (CaptureArea = 3)
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Grid Connectivity Property types                               *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    integer*4 :: AverageAll, AverageCircumferential, AverageRadial
    integer*4 :: AverageI, AverageJ, AverageK
    character(32) :: AverageInterfaceTypeName(0:7)
    parameter (AverageAll = 2)
    parameter (AverageCircumferential = 3)
    parameter (AverageRadial = 4)
    parameter (AverageI = 5)
    parameter (AverageJ = 6)
    parameter (AverageK = 7)
    
! For portability to Linux Absoft, all data statements were moved after the
! variables and parametres declarations
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Dimensional Units                                              *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    data MassUnitsName /'Null','UserDefined','Kilogram','Gram',     & &
    'Slug','PoundMass'/
    data LengthUnitsName / 'Null', 'UserDefined',                   & &
    'Meter','Centimeter','Millimeter','Foot','Inch'/
    
    data TimeUnitsName /'Null','UserDefined','Second'/
    
    data TemperatureUnitsName /'Null','UserDefined',                & &
    'Kelvin','Celsius','Rankine','Fahrenheit'/
    
    data AngleUnitsName /'Null','UserDefined','Degree','Radian'/
    
    data ElectricCurrentUnitsName /'Null', 'UserDefined', 'Ampere', & &
    'Abampere', 'Statampere', 'Edison', 'a.u.'/
    
    data SubstanceAmountUnitsName /'Null', 'UserDefined', 'Mole',   & &
    'Entities', 'StandardCubicFoot', 'StandardCubicMeter'/
    
    data LuminousIntensityUnitsName /'Null', 'UserDefined',         & &
    'Candela', 'Candle', 'Carcel', 'Hefner', 'Violle'/
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Data Class                                                     *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    data DataClassName / 'Null','UserDefined',                      & &
    'Dimensional','NormalizedByDimensional',                & &
    'NormalizedByUnknownDimensional',                       & &
    'NondimensionalParameter','DimensionlessConstant'/
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Grid Location                                                  *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data GridLocationName / 'Null','UserDefined',                   & &
    'Vertex','CellCenter','FaceCenter','IFaceCenter',       & &
    'JFaceCenter','KFaceCenter','EdgeCenter' /
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Grid Connectivity Types                                        *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data GridConnectivityTypeName / 'Null','UserDefined',           & &
    'Overset','Abutting','Abutting1to1'/
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Point Set Types                                                *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data PointSetTypeName / 'Null','UserDefined',                   & &
    'PointList','PointListDonor',                           & &
    'PointRange','PointRangeDonor',                         & &
    'ElementRange','ElementList','CellListDonor'/
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Governing Equations and Physical Models Types                  *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data GoverningEquationsTypeName / 'Null','UserDefined',         & &
    'FullPotential','Euler', 'NSLaminar', 'NSTurbulent',    & &
    'NSLaminarIncompressible', 'NSTurbulentIncompressible'/
    
    data ModelTypeName / 'Null','UserDefined',                      & &
    'Ideal','VanderWaals', 'Constant','PowerLaw',             & &
    'SutherlandLaw','ConstantPrandtl','EddyViscosity',        & &
    'ReynoldsStress','ReynoldsStressAlgebraic',               & &
    'Algebraic_BaldwinLomax','Algebraic_CebeciSmith',         & &
    'HalfEquation_JohnsonKing','OneEquation_BaldwinBarth',    & &
    'OneEquation_SpalartAllmaras','TwoEquation_JonesLaunder', & &
    'TwoEquation_MenterSST','TwoEquation_Wilcox',             & &
    'CaloricallyPerfect', 'ThermallyPerfect',                 & &
    'ConstantDensity', 'RedlichKwong', 'Frozen',              & &
    'ThermalEquilib', 'ThermalNonequilib',                    & &
    'ChemicalEquilibCurveFit', 'ChemicalEquilibMinimization', & &
    'ChemicalNonequilib', 'EMElectricField',                  & &
    'EMMagneticField', 'EMConductivity', 'Voltage',           & &
    'Interpolated', 'Equilibrium_LinRessler',                 & &
    'Chemistry_LinRessler'/
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Boundary Condition Types                                       *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data BCTypeName / 'Null','UserDefined',                         & &
    'BCAxisymmetricWedge','BCDegenerateLine',               & &
    'BCDegeneratePoint','BCDirichlet','BCExtrapolate',      & &
    'BCFarfield','BCGeneral','BCInflow','BCInflowSubsonic', & &
    'BCInflowSupersonic','BCNeumann','BCOutflow',           & &
    'BCOutflowSubsonic','BCOutflowSupersonic',              & &
    'BCSymmetryPlane','BCSymmetryPolar','BCTunnelInflow',   & &
    'BCTunnelOutflow','BCWall','BCWallInviscid',            & &
    'BCWallViscous','BCWallViscousHeatFlux',                & &
    'BCWallViscousIsothermal','FamilySpecified' /
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Data types                                                     *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data DataTypeName / 'Null','UserDefined',                       & &
    'Integer','RealSingle','RealDouble','Character',        & &
    'LongInteger' /
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      BCData_t types                                                 *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data BCDataTypeName / 'Null','UserDefined',                     & &
    'Dirichlet', 'Neumann' /
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Element types                                                  *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data ElementTypeName / 'Null','UserDefined',                    & &
    'NODE', 'BAR_2', 'BAR_3', 'TRI_3', 'TRI_6',                 & &
    'QUAD_4', 'QUAD_8', 'QUAD_9', 'TETRA_4', 'TETRA_10',        & &
    'PYRA_5', 'PYRA_14', 'PENTA_6', 'PENTA_15',                 & &
    'PENTA_18', 'HEXA_8', 'HEXA_20', 'HEXA_27', 'MIXED',        & &
    'PYRA_13', 'NGON_n', 'NFACE_n',                             & &
    'BAR_4', 'TRI_9', 'TRI_10',                                 & &
    'QUAD_12', 'QUAD_16',                                       & &
    'TETRA_16', 'TETRA_20',                                     & &
    'PYRA_21', 'PYRA_29', 'PYRA_30',                            & &
    'PENTA_24', 'PENTA_38', 'PENTA_40',                         & &
    'HEXA_32', 'HEXA_56', 'HEXA_64' /
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Zone types                                                     *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data ZoneTypeName / 'Null','UserDefined',                       & &
    'Structured', 'Unstructured' /
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Rigid Grid Motion types                                        *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data RigidGridMotionTypeName / 'Null','UserDefined',            & &
    'ConstantRate', 'VariableRate' /
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Arbitrary Grid Motion types                                    *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data ArbitraryGridMotionTypeName / 'Null','UserDefined',        & &
    'NonDeformingGrid', 'DeformingGrid' /
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Simulation type                                                *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data SimulationTypeName / 'Null','UserDefined',                 & &
    'TimeAccurate', 'NonTimeAccurate' /
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      BC Property types                                              *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data WallFunctionTypeName / 'Null','UserDefined',               & &
    'Generic' /
    
    data AreaTypeName / 'Null','UserDefined',                       & &
    'BleedArea', 'CaptureArea' /
    
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!*      Grid Connectivity Property types                               *
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
    
    data AverageInterfaceTypeName / 'Null','UserDefined',           & &
    'AverageAll', 'AverageCircumferential', 'AverageRadial',   & &
    'AverageI', 'AverageJ', 'AverageK' /

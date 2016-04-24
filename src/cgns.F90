MODULE cgns_types

IMPLICIT NONE

  PUBLIC

  INTEGER, PARAMETER :: CG_BUILD_64BIT=1

  INTEGER, PARAMETER :: CGSIZE_T=8
  INTEGER, PARAMETER :: CGLONG_T=8
  INTEGER, PARAMETER :: CGID_T=8

END MODULE



MODULE cgio

  USE :: cgns_types, ONLY: CGSIZE_T

IMPLICIT NONE

  PUBLIC

  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_MODE_READ=0
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_MODE_WRITE=1
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_MODE_MODIFY=2

  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_FILE_NONE=0
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_FILE_ADF=1
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_FILE_HDF5=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_FILE_ADF2=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_FILE_PHDF5=4

  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_MAX_DATATYPE_LENGTH=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_MAX_DIMENSIONS=12
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_MAX_NAME_LENGTH=32
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_MAX_LABEL_LENGTH=32
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_MAX_VERSION_LENGTH=32
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_MAX_DATE_LENGTH=32
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_MAX_ERROR_LENGTH=80
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_MAX_LINK_DEPTH=100
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_MAX_FILE_LENGTH=1024
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_MAX_LINK_LENGTH=4096

  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_NONE=0
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_BAD_CGIO=-1
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_MALLOC=-2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_FILE_MODE=-3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_FILE_TYPE=-4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_NULL_FILE=-5
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_TOO_SMALL=-6
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_NOT_FOUND=-7
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_NULL_PATH=-8
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_NO_MATCH=-9
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_FILE_OPEN=-10
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_READ_ONLY=-11
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_NULL_STRING=-12
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_BAD_OPTION=-13
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_FILE_RENAME=-14
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_TOO_MANY=-15
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_DIMENSIONS=-16
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_BAD_TYPE=-17
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CGIO_ERR_NOT_HDF5=-18


  PRIVATE :: CGSIZE_T

END MODULE cgio



!* ------------------------------------------------------------------------- *
!* CGNS - CFD General Notation System (http://www.cgns.org)                  *
!* CGNS/MLL - Mid-Level Library header file                                  *
!* Please see cgnsconfig.h file for this local installation configuration    *
!* ------------------------------------------------------------------------- *
!
!* ------------------------------------------------------------------------- *
!
!This software is provided 'as-is', without any express or implied warranty.
!In no event will the authors be held liable for any damages arising from
!the use of this software.
!
!Permission is granted to anyone to use this software for any purpose,
!including commercial applications, and to alter it and redistribute it
!freely, subject to the following restrictions:
!
!1. The origin of this software must not be misrepresented; you must not
!claim that you wrote the original software. If you use this software
!in a product, an acknowledgment in the product documentation would be
!appreciated but is not required.
!
!2. Altered source versions must be plainly marked as such, and must not
!be misrepresented as being the original software.
!
!3. This notice may not be removed or altered from any source distribution.
!
!* ------------------------------------------------------------------------- *
!
MODULE cgnslib

  USE :: cgns_types, ONLY: CGSIZE_T

IMPLICIT NONE

  PUBLIC


  ! Fortran version of cgnslib.h
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_BUILD_64BIT=1
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      modes for cgns file                                            *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_MODE_READ=0
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_MODE_WRITE=1
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_MODE_MODIFY=2
  !* legacy code support
  INTEGER(KIND=CGSIZE_T), PARAMETER :: MODE_READ=0
  INTEGER(KIND=CGSIZE_T), PARAMETER :: MODE_WRITE=1
  INTEGER(KIND=CGSIZE_T), PARAMETER :: MODE_MODIFY=2

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      file types                                                     *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_FILE_NONE=0
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_FILE_ADF=1
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_FILE_HDF5=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_FILE_ADF2=3

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      some error code                                                *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_OK=0
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_ERROR=1
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_NODE_NOT_FOUND=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_INCORRECT_PATH=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_NO_INDEX_DIM=4
  !* legacy code support
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ALL_OK=0
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ERROR=1
  INTEGER(KIND=CGSIZE_T), PARAMETER :: NODE_NOT_FOUND=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: INCORRECT_PATH=3

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Dimensional Units                                              *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_Null=0
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CG_UserDefined=1
  !* legacy code support
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Null=0
  INTEGER(KIND=CGSIZE_T), PARAMETER :: UserDefined=1

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Kilogram=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Gram=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Slug=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PoundMass=5

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Meter=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Centimeter=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Millimeter=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Foot=5
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Inch=6

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Second=2

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Kelvin=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Celsius=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Rankine=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Fahrenheit=5

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Degree=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Radian=3

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Ampere=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Abampere=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Statampere=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Edison=5
  INTEGER(KIND=CGSIZE_T), PARAMETER :: auCurrent=6

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Mole=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Entities=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: StandardCubicFoot=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: StandardCubicMeter=5

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Candela=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Candle=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Carcel=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Hefner=5
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Violle=6

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Data Class                                                     *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Dimensional=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: NormalizedByDimensional=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: NormalizedByUnknownDimensional=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: NondimensionalParameter=5
  INTEGER(KIND=CGSIZE_T), PARAMETER :: DimensionlessConstant=6

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Grid Location                                                  *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Vertex=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CellCenter=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: FaceCenter=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: IFaceCenter=5
  INTEGER(KIND=CGSIZE_T), PARAMETER :: JFaceCenter=6
  INTEGER(KIND=CGSIZE_T), PARAMETER :: KFaceCenter=7
  INTEGER(KIND=CGSIZE_T), PARAMETER :: EdgeCenter=8

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Grid Connectivity Types                                        *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Overset=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Abutting=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Abutting1to1=4

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Point Set Types                                                *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: PointList=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PointListDonor=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PointRange=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PointRangeDonor=5
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ElementRange=6
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ElementList=7
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CellListDonor=8

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Governing Equations and Physical Models Types                  *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: FullPotential=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Euler=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: NSLaminar=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: NSTurbulent=5
  INTEGER(KIND=CGSIZE_T), PARAMETER :: NSLaminarIncompressible=6
  INTEGER(KIND=CGSIZE_T), PARAMETER :: NSTurbulentIncompressible=7

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
  !**                             Chemistry_LinRessler


  INTEGER(KIND=CGSIZE_T), PARAMETER :: Ideal=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: VanderWaals=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Constant=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PowerLaw=5
  INTEGER(KIND=CGSIZE_T), PARAMETER :: SutherlandLaw=6
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ConstantPrandtl=7
  INTEGER(KIND=CGSIZE_T), PARAMETER :: EddyViscosity=8
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ReynoldsStress=9
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ReynoldsStressAlgebraic=10
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Algebraic_BaldwinLomax=11
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Algebraic_CebeciSmith=12
  INTEGER(KIND=CGSIZE_T), PARAMETER :: HalfEquation_JohnsonKing=13
  INTEGER(KIND=CGSIZE_T), PARAMETER :: OneEquation_BaldwinBarth=14
  INTEGER(KIND=CGSIZE_T), PARAMETER :: OneEquation_SpalartAllmaras=15
  INTEGER(KIND=CGSIZE_T), PARAMETER :: TwoEquation_JonesLaunder=16
  INTEGER(KIND=CGSIZE_T), PARAMETER :: TwoEquation_MenterSST=17
  INTEGER(KIND=CGSIZE_T), PARAMETER :: TwoEquation_Wilcox=18
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CaloricallyPerfect=19
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ThermallyPerfect=20
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ConstantDensity=21
  INTEGER(KIND=CGSIZE_T), PARAMETER :: RedlichKwong=22
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Frozen=23
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ThermalEquilib=24
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ThermalNonequilib=25
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ChemicalEquilibCurveFit=26
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ChemicalEquilibMinimization=27
  INTEGER(KIND=CGSIZE_T), PARAMETER :: ChemicalNonequilib=28
  INTEGER(KIND=CGSIZE_T), PARAMETER :: EMElectricField=29
  INTEGER(KIND=CGSIZE_T), PARAMETER :: EMMagneticField=30
  INTEGER(KIND=CGSIZE_T), PARAMETER :: EMConductivity=31
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Voltage=32
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Interpolated=33
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Equilibrium_LinRessler=34
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Chemistry_LinRessler=35

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Boundary Condition Types                                       *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCAxisymmetricWedge=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCDegenerateLine=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCDegeneratePoint=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCDirichlet=5
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCExtrapolate=6
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCFarfield=7
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCGeneral=8
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCInflow=9
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCInflowSubsonic=10
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCInflowSupersonic=11
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCNeumann=12
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCOutflow=13
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCOutflowSubsonic=14
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCOutflowSupersonic=15
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCSymmetryPlane=16
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCSymmetryPolar=17
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCTunnelInflow=18
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCTunnelOutflow=19
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCWall=20
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCWallInviscid=21
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCWallViscous=22
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCWallViscousHeatFlux=23
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BCWallViscousIsothermal=24
  INTEGER(KIND=CGSIZE_T), PARAMETER :: FamilySpecified=25

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Data types                                                     *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Integer=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: RealSingle=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: RealDouble=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Character=5
  INTEGER(KIND=CGSIZE_T), PARAMETER :: LongInteger=6

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      BCData_t types                                                 *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Dirichlet=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Neumann=3

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Element types                                                  *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: NODE=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BAR_2=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BAR_3=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: TRI_3=5
  INTEGER(KIND=CGSIZE_T), PARAMETER :: TRI_6=6
  INTEGER(KIND=CGSIZE_T), PARAMETER :: QUAD_4=7
  INTEGER(KIND=CGSIZE_T), PARAMETER :: QUAD_8=8
  INTEGER(KIND=CGSIZE_T), PARAMETER :: QUAD_9=9
  INTEGER(KIND=CGSIZE_T), PARAMETER :: TETRA_4=10
  INTEGER(KIND=CGSIZE_T), PARAMETER :: TETRA_10=11
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PYRA_5=12
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PYRA_14=13
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PENTA_6=14
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PENTA_15=15
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PENTA_18=16
  INTEGER(KIND=CGSIZE_T), PARAMETER :: HEXA_8=17
  INTEGER(KIND=CGSIZE_T), PARAMETER :: HEXA_20=18
  INTEGER(KIND=CGSIZE_T), PARAMETER :: HEXA_27=19
  INTEGER(KIND=CGSIZE_T), PARAMETER :: MIXED=20
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PYRA_13=21
  INTEGER(KIND=CGSIZE_T), PARAMETER :: NGON_n=22
  INTEGER(KIND=CGSIZE_T), PARAMETER :: NFACE_n=23
  INTEGER(KIND=CGSIZE_T), PARAMETER :: BAR_4=24
  INTEGER(KIND=CGSIZE_T), PARAMETER :: TRI_9=25
  INTEGER(KIND=CGSIZE_T), PARAMETER :: TRI_10=26
  INTEGER(KIND=CGSIZE_T), PARAMETER :: QUAD_12=27
  INTEGER(KIND=CGSIZE_T), PARAMETER :: QUAD_16=28
  INTEGER(KIND=CGSIZE_T), PARAMETER :: TETRA_16=29
  INTEGER(KIND=CGSIZE_T), PARAMETER :: TETRA_20=30
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PYRA_21=31
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PYRA_29=32
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PYRA_30=33
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PENTA_24=34
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PENTA_38=35
  INTEGER(KIND=CGSIZE_T), PARAMETER :: PENTA_40=36
  INTEGER(KIND=CGSIZE_T), PARAMETER :: HEXA_32=37
  INTEGER(KIND=CGSIZE_T), PARAMETER :: HEXA_56=38
  INTEGER(KIND=CGSIZE_T), PARAMETER :: HEXA_64=39

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Zone types                                                     *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Structured=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: Unstructured=3

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Rigid Grid Motion types                                        *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: ConstantRate=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: VariableRate=3

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Arbitrary Grid Motion types                                    *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: NonDeformingGrid=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: DeformingGrid=3

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Simulation type                                                *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: TimeAccurate=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: NonTimeAccurate=3

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      BC Property types                                              *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: Generic=2

  INTEGER(KIND=CGSIZE_T), PARAMETER :: BleedArea=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: CaptureArea=3

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Grid Connectivity Property types                               *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  INTEGER(KIND=CGSIZE_T), PARAMETER :: AverageAll=2
  INTEGER(KIND=CGSIZE_T), PARAMETER :: AverageCircumferential=3
  INTEGER(KIND=CGSIZE_T), PARAMETER :: AverageRadial=4
  INTEGER(KIND=CGSIZE_T), PARAMETER :: AverageI=5
  INTEGER(KIND=CGSIZE_T), PARAMETER :: AverageJ=6
  INTEGER(KIND=CGSIZE_T), PARAMETER :: AverageK=7

  ! For portability to Linux Absoft, all data statements were moved after the
  ! variables and parametres declarations

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Dimensional Units                                              *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  CHARACTER(LEN=32), DIMENSION(0:5), PARAMETER :: MassUnitsName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Kilogram                        ', &
                                  'Gram                            ', &
                                  'Slug                            ', &
                                  'PoundMass                       ' /)
  CHARACTER(LEN=32), DIMENSION(0:6), PARAMETER :: LengthUnitsName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Meter                           ', &
                                  'Centimeter                      ', &
                                  'Millimeter                      ', &
                                  'Foot                            ', &
                                  'Inch                            ' /)

  CHARACTER(LEN=32), DIMENSION(0:2), PARAMETER :: TimeUnitsName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Second                          ' /)

  CHARACTER(LEN=32), DIMENSION(0:5), PARAMETER :: TemperatureUnitsName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Kelvin                          ', &
                                  'Celsius                         ', &
                                  'Rankine                         ', &
                                  'Fahrenheit                      ' /)

  CHARACTER(LEN=32), DIMENSION(0:3), PARAMETER :: AngleUnitsName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Degree                          ', &
                                  'Radian                          ' /)

  CHARACTER(LEN=32), DIMENSION(0:6), PARAMETER :: ElectricCurrentUnitsName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Ampere                          ', &
                                  'Abampere                        ', &
                                  'Statampere                      ', &
                                  'Edison                          ', &
                                  'a.u.                            ' /)

  CHARACTER(LEN=32), DIMENSION(0:5), PARAMETER :: SubstanceAmountUnitsName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Mole                            ', &
                                  'Entities                        ', &
                                  'StandardCubicFoot               ', &
                                  'StandardCubicMeter              ' /)

  CHARACTER(LEN=32), DIMENSION(0:6), PARAMETER :: LuminousIntensityUnitsName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Candela                         ', &
                                  'Candle                          ', &
                                  'Carcel                          ', &
                                  'Hefner                          ', &
                                  'Violle                          ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Data Class                                                     *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  CHARACTER(LEN=32), DIMENSION(0:6), PARAMETER :: DataClassName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Dimensional                     ', &
                                  'NormalizedByDimensional         ', &
                                  'NormalizedByUnknownDimensional  ', &
                                  'NondimensionalParameter         ', &
                                  'DimensionlessConstant           ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Grid Location                                                  *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:8), PARAMETER :: GridLocationName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Vertex                          ', &
                                  'CellCenter                      ', &
                                  'FaceCenter                      ', &
                                  'IFaceCenter                     ', &
                                  'JFaceCenter                     ', &
                                  'KFaceCenter                     ', &
                                  'EdgeCenter                      ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Grid Connectivity Types                                        *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:4), PARAMETER :: GridConnectivityTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Overset                         ', &
                                  'Abutting                        ', &
                                  'Abutting1to1                    ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Point Set Types                                                *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:8), PARAMETER :: PointSetTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'PointList                       ', &
                                  'PointListDonor                  ', &
                                  'PointRange                      ', &
                                  'PointRangeDonor                 ', &
                                  'ElementRange                    ', &
                                  'ElementList                     ', &
                                  'CellListDonor                   ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Governing Equations and Physical Models Types                  *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:7), PARAMETER :: GoverningEquationsTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'FullPotential                   ', &
                                  'Euler                           ', &
                                  'NSLaminar                       ', &
                                  'NSTurbulent                     ', &
                                  'NSLaminarIncompressible         ', &
                                  'NSTurbulentIncompressible       ' /)

  CHARACTER(LEN=32), DIMENSION(0:35), PARAMETER :: ModelTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Ideal                           ', &
                                  'VanderWaals                     ', &
                                  'Constant                        ', &
                                  'PowerLaw                        ', &
                                  'SutherlandLaw                   ', &
                                  'ConstantPrandtl                 ', &
                                  'EddyViscosity                   ', &
                                  'ReynoldsStress                  ', &
                                  'ReynoldsStressAlgebraic         ', &
                                  'Algebraic_BaldwinLomax          ', &
                                  'Algebraic_CebeciSmith           ', &
                                  'HalfEquation_JohnsonKing        ', &
                                  'OneEquation_BaldwinBarth        ', &
                                  'OneEquation_SpalartAllmaras     ', &
                                  'TwoEquation_JonesLaunder        ', &
                                  'TwoEquation_MenterSST           ', &
                                  'TwoEquation_Wilcox              ', &
                                  'CaloricallyPerfect              ', &
                                  'ThermallyPerfect                ', &
                                  'ConstantDensity                 ', &
                                  'RedlichKwong                    ', &
                                  'Frozen                          ', &
                                  'ThermalEquilib                  ', &
                                  'ThermalNonequilib               ', &
                                  'ChemicalEquilibCurveFit         ', &
                                  'ChemicalEquilibMinimization     ', &
                                  'ChemicalNonequilib              ', &
                                  'EMElectricField                 ', &
                                  'EMMagneticField                 ', &
                                  'EMConductivity                  ', &
                                  'Voltage                         ', &
                                  'Interpolated                    ', &
                                  'Equilibrium_LinRessler          ', &
                                  'Chemistry_LinRessler            ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Boundary Condition Types                                       *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:25), PARAMETER :: BCTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'BCAxisymmetricWedge             ', &
                                  'BCDegenerateLine                ', &
                                  'BCDegeneratePoint               ', &
                                  'BCDirichlet                     ', &
                                  'BCExtrapolate                   ', &
                                  'BCFarfield                      ', &
                                  'BCGeneral                       ', &
                                  'BCInflow                        ', &
                                  'BCInflowSubsonic                ', &
                                  'BCInflowSupersonic              ', &
                                  'BCNeumann                       ', &
                                  'BCOutflow                       ', &
                                  'BCOutflowSubsonic               ', &
                                  'BCOutflowSupersonic             ', &
                                  'BCSymmetryPlane                 ', &
                                  'BCSymmetryPolar                 ', &
                                  'BCTunnelInflow                  ', &
                                  'BCTunnelOutflow                 ', &
                                  'BCWall                          ', &
                                  'BCWallInviscid                  ', &
                                  'BCWallViscous                   ', &
                                  'BCWallViscousHeatFlux           ', &
                                  'BCWallViscousIsothermal         ', &
                                  'FamilySpecified                 ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Data types                                                     *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:6), PARAMETER :: DataTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Integer                         ', &
                                  'RealSingle                      ', &
                                  'RealDouble                      ', &
                                  'Character                       ', &
                                  'LongInteger                     ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      BCData_t types                                                 *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:3), PARAMETER :: BCDataTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Dirichlet                       ', &
                                  'Neumann                         ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Element types                                                  *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:39), PARAMETER :: ElementTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'NODE                            ', &
                                  'BAR_2                           ', &
                                  'BAR_3                           ', &
                                  'TRI_3                           ', &
                                  'TRI_6                           ', &
                                  'QUAD_4                          ', &
                                  'QUAD_8                          ', &
                                  'QUAD_9                          ', &
                                  'TETRA_4                         ', &
                                  'TETRA_10                        ', &
                                  'PYRA_5                          ', &
                                  'PYRA_14                         ', &
                                  'PENTA_6                         ', &
                                  'PENTA_15                        ', &
                                  'PENTA_18                        ', &
                                  'HEXA_8                          ', &
                                  'HEXA_20                         ', &
                                  'HEXA_27                         ', &
                                  'MIXED                           ', &
                                  'PYRA_13                         ', &
                                  'NGON_n                          ', &
                                  'NFACE_n                         ', &
                                  'BAR_4                           ', &
                                  'TRI_9                           ', &
                                  'TRI_10                          ', &
                                  'QUAD_12                         ', &
                                  'QUAD_16                         ', &
                                  'TETRA_16                        ', &
                                  'TETRA_20                        ', &
                                  'PYRA_21                         ', &
                                  'PYRA_29                         ', &
                                  'PYRA_30                         ', &
                                  'PENTA_24                        ', &
                                  'PENTA_38                        ', &
                                  'PENTA_40                        ', &
                                  'HEXA_32                         ', &
                                  'HEXA_56                         ', &
                                  'HEXA_64                         ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Zone types                                                     *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:3), PARAMETER :: ZoneTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Structured                      ', &
                                  'Unstructured                    ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Rigid Grid Motion types                                        *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:3), PARAMETER :: RigidGridMotionTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'ConstantRate                    ', &
                                  'VariableRate                    ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Arbitrary Grid Motion types                                    *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:3), PARAMETER :: ArbitraryGridMotionTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'NonDeformingGrid                ', &
                                  'DeformingGrid                   ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Simulation type                                                *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:3), PARAMETER :: SimulationTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'TimeAccurate                    ', &
                                  'NonTimeAccurate                 ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      BC Property types                                              *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:2), PARAMETER :: WallFunctionTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'Generic                         ' /)

  CHARACTER(LEN=32), DIMENSION(0:3), PARAMETER :: AreaTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'BleedArea                       ', &
                                  'CaptureArea                     ' /)

  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  !*      Grid Connectivity Property types                               *
  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *

  CHARACTER(LEN=32), DIMENSION(0:7), PARAMETER :: AverageInterfaceTypeName= (/ &
                                  'Null                            ', &
                                  'UserDefined                     ', &
                                  'AverageAll                      ', &
                                  'AverageCircumferential          ', &
                                  'AverageRadial                   ', &
                                  'AverageI                        ', &
                                  'AverageJ                        ', &
                                  'AverageK                        ' /)

  PRIVATE :: CGSIZE_T

END MODULE cgnslib

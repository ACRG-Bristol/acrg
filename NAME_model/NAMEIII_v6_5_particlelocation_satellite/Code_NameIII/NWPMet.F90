! Module:  NWP Met Module

Module NWPMetModule

! This is a met module designed to read gridded NWP met data.

! It uses input from files of NWP met data and a file of topographic height data.
!
! The coords stored in NWPMet%C are
! Horizontal: 1: Coord system for horizontal grids used for met data.
! Vertical:   1: Coord system of vertical grids used for met data.
!             2: m agl.
!             3: Pa.
!
! No grids are stored in NWPMet%C.
!
! It only works in a calendar time frame. $$ add test for this.
!
! The status of the optional met module features is as follows (where ???? here =
! NWP):
!     SetUpCoordsEtc_????Met routine - used
!     SetUp????Met_CoordsEtc routine - used
!     Reset????Met routine           - used.
!
! Time labelling of met data:
! Each set of met data must be labelled with the time.
!
! Varying met data between cases:
! This is possible with multiple sets of met files.
! The cases step through the sets of met files.
!
! Varying met with FixedMet = true by incrementing FixedMetTime: $$
! This is allowed.

! Notes on the format of Name II met files and associated topography files.
!
! These are unformatted files containing a series of 2-d (horizontal) fields. 3-d
! fields are read in as a series of 2-d (horizontal) fields. Each field contains a
! header record followed, if the field isn't missing, by nY records of nX 32-bit reals
! giving the values on the nX by nY grid on which the field is stored. The header
! record contains 64 32-bit reals of which the following are significant here:
!     Header(1)  :: If not > 0, this means the field is missing.
!     Header(14) :: Forecast period in hours. $$ Add value to met read message?.
!     Header(34) ::
!     Header(35) ::
!     Header(36) ::
!     Header(37) ::
! The coord systems and grids used and the order in which the fields appear in the
! file need to be specified in instances of the types NWPMetDefn_ and NWPMetDefn2_
! respectively. This information is read in from the main set of headed input files
! in 'NWP Met Definitions' and 'NWP Met File Structure' blocks.
!
! The associated topography file is an unformatted file containing a single 2-d
! (horizontal) field of topographic height in the Name II format described above,
! but without the header record.

! Notes on the format of PP met files and associated topography files.
!
! These are unformatted files containing a series of 2-d (horizontal) fields. 3-d
! fields are read in as a series of 2-d (horizontal) fields. Each field contains a
! header record followed, if the field isn't missing, by a single record of nX*nY
! 32-bit reals giving the values on the nX by nY grid on which the field is stored.
! The header record contains 45 32-bit integers and 19 32-bit reals of which the
! following are significant here:
!     $$ currently none used. This should probably change (e.g. test for missing
!        fields - or change comments above)
! The coord systems and grids used and the order in which the fields appear in the
! file need to be specified in instances of the types NWPMetDefn_ and NWPMetDefn2_
! respectively. This information is read in from the main set of headed input files
! in 'NWP Met Definitions' and 'NWP Met File Structure' blocks.
!
! The associated topography file is an unformatted file containing a single 2-d
! (horizontal) field of topographic height in the PP format described above.

! Notes on coord systems and grids for met and topography and on allowed and required
! met fields.
!
! There are a number of restrictions on the possible coord systems and grids and the
! fields that can be and are required to be input. These are tested for in the code.
! If a field can be omitted in a met file structure definition (held in an
! NWPMetDefn2_ instance and input in an 'NWP Met File Structure' block), then the
! field can also be missing in the met file even if it is given in a met file
! structure definition.
!
! The topography data must be on the main horizontal grid used for the met data (i.e.
! the grid used for everything with the possible exception of the u and v components
! of wind and surface stress).

! $$ list restrictions on coords/grids here and possible/required fields with coords,
! grids, units (inc topog units) and special values (e.g. cloud base when no cloud).
! horizontal wind and stress cpts orientated with grid.
! Document restrictions in input documentation.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use CommonMetModule
Use TimerModule
Use NWPUtilitiesModule

Use OpenMPModule, Only : IOMaxLookahead,         &
                         OpenMPOpts_,            &
                         OpenMPMaxThreads

!$ Use omp_lib

# ifdef GRIBsupport
  ! If support for GRIB is requested, the code references external routines in the GRIB_API module, which are
  ! linked dynamically at run time from their shared object libraries. The GRIB_API is an ECMWF package for
  ! encoding/decoding gridded meteorological fields in GRIB-1 and GRIB-2 formats. For further details, see
  ! the ReadMe file for optional third-party software used by NAME.
  Use grib_api
# endif

# ifdef NetCDFsupport
  ! If support for NetCDF is requested, the code references external routines in a NetCDF API, which is
  ! linked dynamically at run time from shared object libraries. For further details, see the ReadMe file
  ! for optional third-party software used by NAME.
  Use netcdf
# endif

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: MA_NWPFlow             !} Codes for the module attributes. These indicate that the module can supply
Public :: MA_NWPCloud            !} flow, cloud, rain, soil moisture, canopy height and canopy information.
Public :: MA_NWPRain             !} Exactly what is needed/supplied for each attribute is noted
                                 !} in the Needed3d
Public :: MA_NWPSoilMoisture     !} and Needed2d arrays.
Public :: MA_NWPCanopyHeight     !}
Public :: MA_NWPCanopy           !}
Public :: NWPMetDefn_            ! An NWP met definition (part 1 - basic definitions).
Public :: NWPMetDefn2_           ! An NWP met definition (part 2 - met file
                                 ! structure definition).
Public :: MetDefns_              ! A collection of met definitions.
Public :: NWPMet_                ! The state of an NWP met module instance.
Public :: InitMetDefns           ! Initialises a collection of met definitions.
Public :: InitNWPMetDefn         ! Initialises an NWP met definition (part 1).
Public :: AddNWPMetDefn          ! Adds an NWP met definition (part 1) to a collection
                                 ! of met definitions.
Public :: FindNWPMetDefnIndex    ! Finds the index of an NWP met definition (part 1).
Public :: InitNWPMetDefn2        ! Initialises an NWP met definition (part 2).
Public :: AddNWPMetDefn2         ! Adds an NWP met definition (part 2) to a collection
                                 ! of met definitions.
Public :: FindNWPMetDefn2Index   ! Finds the index of an NWP met definition (part 2).
Public :: InitNWPMet             ! Initialises the state of an NWP met module
                                 ! instance.
Public :: SetUpCoordsEtc_NWPMet  ! Sets up Coords and Grids by adding any extra coords
                                 ! and grids which NWPMet wants to define.
Public :: SetUpNWPMet_CoordsEtc  ! Sets up NWPMet using information from EtaDefns,
                                 ! Coords and Grids.
Public :: PrepareForUpdateNWPMet ! Prepares for updating an instance of the NWP met module.
Public :: UpdateNWPMet           ! Updates an instance of the NWP met module.
Public :: ResetNWPMet            ! Resets the state of an NWP met module instance for
                                 ! a new realisation.
Public :: ReadNWPMet             ! } These subroutines have to be public so that the IO thread can call them
Public :: ProcessNWPMet          ! } in IOThreadModule.F90
Public :: NWPMetModuleTimerSummary
Public :: NWPMetModuleTimerInitialise
Public :: EulerianModelGlobal

!-------------------------------------------------------------------------------------------------------------

Logical, Save :: EulerianModelGlobal  ! $$ This global variable should be removed in due course

!-------------------------------------------------------------------------------------------------------------

! Codes for the module attributes.
Integer, Parameter :: MA_NWPFlow         = 1 !} Codes for the module attributes. These indicate that the
Integer, Parameter :: MA_NWPCloud        = 2 !} module can supply flow, cloud, rain, soil moisture, canopy
Integer, Parameter :: MA_NWPRain         = 3 !} height and canopy information. Exactly what is needed/supplied
Integer, Parameter :: MA_NWPSoilMoisture = 4 !} for each attribute is noted in the Needed3d and Needed2d
Integer, Parameter :: MA_NWPCanopyHeight = 5 !} arrays.
Integer, Parameter :: MA_NWPCanopy       = 6 !}

! NWP met field names and information on which fields can be missing:
Integer,                  Parameter :: nFieldNames3d = 11
Integer,                  Parameter :: nFieldNames2d = 20
Character(MaxCharLength), Parameter :: FieldNames3d(nFieldNames3d) =   &
                                       (/                              &
                                         'wind (u-cpt)              ', & ! 1
                                         'wind (v-cpt)              ', & ! 2
                                         'wind (w-cpt)              ', & ! 3
                                         'temperature (K)           ', & ! 4
                                         'specific humidity         ', & ! 5
                                         'pressure (Pa)             ', & ! 6
                                         'cloud liquid water (kg/kg)', & ! 7
                                         'cloud ice (kg/kg)         ', & ! 8
                                         'canopy height (m)         ', & ! 9
                                         'canopy water (kg/m^2)     ', & !10
                                         'stomatal conductance (m/s)'  & !11
                                       /)
Character(MaxCharLength), Parameter :: FieldNames2d(nFieldNames2d) =          &
                                       (/                                     &
                                         'roughness length                 ', & !  1
                                         'surface stress (u-cpt) (N/m^2)   ', & !  2
                                         'surface stress (v-cpt) (N/m^2)   ', & !  3
                                         'surface sensible heat flux       ', & !  4
                                         'boundary layer depth             ', & !  5
                                         'sea level pressure (Pa)          ', & !  6
                                         'convective cloud amount (0-1)    ', & !  7
                                         'convective cloud base            ', & !  8
                                         'convective cloud top             ', & !  9
                                         'high cloud amount (0-1)          ', & ! 10
                                         'medium cloud amount (0-1)        ', & ! 11
                                         'low cloud amount (0-1)           ', & ! 12
                                         'convective rain rate (kg/(m^2 s))', & ! 13
                                         'convective snow rate (kg/(m^2 s))', & ! 14
                                         'dynamic rain rate (kg/(m^2 s))   ', & ! 15
                                         'dynamic snow rate (kg/(m^2 s))   ', & ! 16
                                         'dummy                            ', & ! 17
                                         'dummy (u-grid)                   ', & ! 18
                                         'dummy (v-grid)                   ', & ! 19
                                         'soil moisture in layer (kg/m^2)  '  & ! 20
                                       /)
Integer,                  Parameter :: ExtrapolationLevelsBelow(nFieldNames3d) = &
                                       (/                                        &
                                         3,                                      & ! 1
                                         3,                                      & ! 2
                                         3,                                      & ! 3
                                         3,                                      & ! 4
                                         2,                                      & ! 5
                                         3,                                      & ! 6
                                         2,                                      & ! 7
                                         2,                                      & ! 8
                                         0,                                      & ! 9
                                         0,                                      & !10
                                         0                                       & !11
                                       /)
Integer,                  Parameter :: ExtrapolationLevelsAbove(nFieldNames3d) = &
                                       (/                                        &
                                         1,                                      & ! 1
                                         1,                                      & ! 2
                                         1,                                      & ! 3
                                         1,                                      & ! 4
                                         0,                                      & ! 5
                                         Huge(ExtrapolationLevelsAbove)/3,       & ! 6
                                         0,                                      & ! 7
                                         0,                                      & ! 8
                                         0,                                      & ! 9
                                         4,                                      & !10
                                         4                                       & !11
                                       /)
Integer,                  Parameter :: InterpolationLevels(nFieldNames3d) = &
                                       (/                                   &
                                         2,                                 & ! 1
                                         2,                                 & ! 2
                                         2,                                 & ! 3
                                         2,                                 & ! 4
                                         0,                                 & ! 5
                                         Huge(InterpolationLevels)/3,       & ! 6
                                         0,                                 & ! 7
                                         0,                                 & ! 8
                                         0,                                 & ! 9
                                         0,                                 & !10
                                         0                                  & !11
                                       /)
Integer,                  Parameter :: PersistTimes3d(nFieldNames3d) = &
                                       (/                              &
                                         0,                            & ! 1
                                         0,                            & ! 2
                                         0,                            & ! 3
                                         0,                            & ! 4
                                         Huge(PersistTimes3d)/3,       & ! 5
                                         Huge(PersistTimes3d)/3,       & ! 6
                                         Huge(PersistTimes3d)/3,       & ! 7
                                         Huge(PersistTimes3d)/3,       & ! 8
                                         0,                            & ! 9
                                         0,                            & !10
                                         0                             & !11
                                       /)
Logical,                  Parameter :: AllowFixUpToDefault3d(nFieldNames3d) = &
                                       (/                                     &
                                         .false.,                             & ! 1
                                         .false.,                             & ! 2
                                         .false.,                             & ! 3
                                         .false.,                             & ! 4
                                         .true.,                              & ! 5
                                         .true.,                              & ! 6
                                         .true.,                              & ! 7
                                         .true.,                              & ! 8
                                         .false.,                             & ! 9
                                         .false.,                              & !10
                                         .false.                               & !11
                                       /)
Integer,                  Parameter :: PersistTimes2d(nFieldNames2d) = &
                                       (/                              &
                                         Huge(PersistTimes2d)/3,       & !  1
                                         0,                            & !  2
                                         0,                            & !  3
                                         0,                            & !  4
                                         0,                            & !  5
                                         Huge(PersistTimes2d)/3,       & !  6
                                         0,                            & !  7
                                         0,                            & !  8
                                         0,                            & !  9
                                         0,                            & ! 10
                                         0,                            & ! 11
                                         0,                            & ! 12
                                         0,                            & ! 13
                                         0,                            & ! 14
                                         0,                            & ! 15
                                         0,                            & ! 16
                                         0,                            & ! 17
                                         0,                            & ! 18
                                         0,                            & ! 19
                                         0                             & ! 20
                                       /)
Logical,                  Parameter :: AllowFixUpToDefault2d(nFieldNames2d) = &
                                       (/                                     &
                                         .true.,                              & !  1
                                         .false.,                             & !  2
                                         .false.,                             & !  3
                                         .false.,                             & !  4
                                         .true.,                              & !  5
                                         .true.,                              & !  6
                                         .true.,                              & !  7
                                         .true.,                              & !  8
                                         .true.,                              & !  9
                                         .false.,                             & ! 10
                                         .false.,                             & ! 11
                                         .false.,                             & ! 12
                                         .true.,                              & ! 13
                                         .true.,                              & ! 14
                                         .false.,                             & ! 15
                                         .false.,                             & ! 16
                                         .true.,                              & ! 17
                                         .true.,                              & ! 18
                                         .true.,                              & ! 19
                                         .false.                              & ! 20
                                       /)
Integer,                  Parameter :: nAttribs = 6
Logical,                  Parameter :: Needed3d(nAttribs, nFieldNames3d) =  &
                                       Reshape(                             &
                                         (/                                 &
!                                          Flow    Cloud   Rain    SoilMoisture   CanopyHeight    Canopy
                                           .true., .true., .true., .false.,    .false.,     .false., & ! 1
                                           .true., .true., .true., .false.,    .false.,     .false., & ! 2
                                           .true., .true., .true., .false.,    .false.,     .false., & ! 3
                                           .true., .true., .true., .false.,    .false.,     .false., & ! 4
                                           .true., .true., .true., .false.,    .false.,     .false., & ! 5
                                           .true., .true., .true., .false.,    .false.,     .false., & ! 6
                                           .true., .true., .true., .false.,    .false.,     .false., & ! 7
                                           .true., .true., .true., .false.,    .false.,     .false., & ! 8
                                           .false.,.false.,.false.,.false.,    .true.,      .false., & ! 9
                                           .false.,.false.,.false.,.false.,    .false.,     .true.,  & !10
                                           .false.,.false.,.false.,.false.,    .false.,     .true.   & !11
                                         /),                                &
                                         (/ 6, 11 /)                         &
                                       )
Logical,                  Parameter :: Needed2d(nAttribs, nFieldNames2d) =     &
                                       Reshape(                                &
                                         (/                                    &
!                                          Flow     Cloud    Rain     SoilMoisture   CanopyHeight    Canopy
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & !  1
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & !  2
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & !  3
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & !  4
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & !  5
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & !  6
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & !  7
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & !  8
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & !  9
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & ! 10
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & ! 11
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & ! 12
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & ! 13
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & ! 14
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & ! 15
                                           .true.,  .true.,  .true.,  .false.,    .false.,    .false., & ! 16
                                           .false., .false., .false., .false.,    .false.,    .false., & ! 17
                                           .false., .false., .false., .false.,    .false.,    .false., & ! 18
                                           .false., .false., .false., .false.,    .false.,    .false., & ! 19
                                           .false., .false., .false., .true.,     .false.,    .false.  & ! 20
                                         /),                                   &
                                         (/ 6, 20 /)                           &
                                       )
! nFieldNames3d            :: Number of 3-d fields.
! nFieldNames2d            :: Number of 2-d (horizontal) fields.
! FieldNames3d             :: Names of 3-d fields.
! FieldNames2d             :: Names of 2-d (horizontal) fields.
! ExtrapolationLevelsBelow :} Variables controlling the vertical extrapolation and interpolation of 3-d met
! ExtrapolationLevelsAbove :} fields and the persistence in time or use of defaults for 3-d and 2-d met
! InterpolationLevels      :} fields. The variables indicate the number of levels which can be extrapolated
! PersistTimes3d           :} below the bottom level present or above the top level present, the number of
! AllowFixUpToDefault3d    :} consecutive levels which can be interpolated, the number of met times over which
! PersistTimes2d           :} persistence is allowed, and whether a fix up to a default value is allowed.
! AllowFixUpToDefault2d    :}
! nAttribs                 :: Number of attributes.
! Needed3d                 :} Indicates which 3-d and 2-d met fields are needed (possibly fixed up) for the
! Needed2d                 :} module to be valid for each of its attributes.

! Fix ups are applied in the order (i) vertical extrapolation/interpolation, (ii) persistence, (iii) default
! values. Extrapolated/interpolated values can be persisted. Default values cannot be persisted, although this
! does not matter except where the default is not a constant. Default values for some variables (pressure and
! boundary layer depth) do depend on other variables and so are not constant (here 'default value used' really
! means 'default estimation method applied'). The other variables are fixed up before calculating pressure
! and/or boundary layer depth (and pressure is, indirectly, one of the 'other variables' for boundary layer
! depth).

! The code currently assumes no 2-d fields need fixing up in any way (persistence or fix up to default) other
! than roughness length, boundary layer depth and sea level pressure. Hence PersistTimes2d and
! AllowFixUpToDefault2d must be zero and false except possibly for roughness length, boundary layer depth and
! sea level pressure (and except for the dummy fields where AllowFixUpToDefault2d must be true).
! $$ add tests for this? alter?

! The requirements for pressure are minimal. The most that is needed is one level of pressure which is used to
! estimate the surface pressure, with the whole pressure field then being reconstructed from surface pressure
! and temperature. ExtrapolationLevelsAbove and InterpolationLevels must be set large to prevent the model
! thinking there is inadequate data where the data are not required. ExtrapolationLevelsBelow controls how
! high a level can be used to estimate the surface pressure, PersistTimes3d controls how long surface pressure
! can be fixed up using persistence, and AllowFixUpToDefault3d controls whether sea level pressure can be used
! to estimate the surface pressure when no appropriate pressure levels are available. AllowFixUpToDefault2d
! controls whether sea level pressure can be fixed up to a default (and then possibly used to estimate the
! surface pressure). Provided surface pressure can be estimated, all pressure levels are computed from surface
! pressure and temperature (and the settings of ExtrapolationLevelsBelow etc do not affect this).
! $$ add tests that these restrictions are satisfied?

! The code will work with all the file types in the ADG met file archive when these files have no missing
! data, provided ExtrapolationLevelsBelow is at least (1, 1, 1, 1, 2, 0, 2, 2, 0, 0, 0) and provided
! AllowFixUpToDefault is true for roughness length and boundary layer depth. In practice these files have some
! missing data and what is needed to cope with missing data is not clear.

! Codes for NWP met fields:
Integer, Parameter :: F_U                   =  1 !} Codes giving the location in the arrays
Integer, Parameter :: F_V                   =  2 !} FieldNames3d and FieldNames2d of the
Integer, Parameter :: F_W                   =  3 !} fields.
Integer, Parameter :: F_T                   =  4 !}
Integer, Parameter :: F_Q                   =  5 !}
Integer, Parameter :: F_PAsRead             =  6 !}
Integer, Parameter :: F_TotalOrDynCloudWater=  7 !}
Integer, Parameter :: F_TotalOrDynCloudIce  =  8 !}
Integer, Parameter :: F_CanopyHeight        =  9 !}
Integer, Parameter :: F_CanopyWater         = 10 !}
Integer, Parameter :: F_StomataConduct      = 11 !}
Integer, Parameter :: F_Z0                  =  1 !}
Integer, Parameter :: F_UStress             =  2 !}
Integer, Parameter :: F_VStress             =  3 !}
Integer, Parameter :: F_HeatFlux            =  4 !}
Integer, Parameter :: F_H                   =  5 !}
Integer, Parameter :: F_PSeaLevel           =  6 !}
Integer, Parameter :: F_ConCloud            =  7 !}
Integer, Parameter :: F_ConCloudBase        =  8 !}
Integer, Parameter :: F_ConCloudTop         =  9 !}
Integer, Parameter :: F_TotalOrDynHCloud    = 10 !}
Integer, Parameter :: F_TotalOrDynMCloud    = 11 !}
Integer, Parameter :: F_TotalOrDynLCloud    = 12 !}
Integer, Parameter :: F_ConRain             = 13 !}
Integer, Parameter :: F_ConSnow             = 14 !}
Integer, Parameter :: F_DynRain       = 15 !}
Integer, Parameter :: F_DynSnow       = 16 !}
Integer, Parameter :: F_Dummy         = 17 !}
Integer, Parameter :: F_DummyU        = 18 !}
Integer, Parameter :: F_DummyV        = 19 !}
Integer, Parameter :: F_SoilMoisture  = 20 !}

!-------------------------------------------------------------------------------------------------------------

Type :: NWPMetDefn_ ! An NWP met definition (part 1 - basic definitions). Together
                    ! with part 2 this provides information on the nature of the NWP
                    ! met module met data.
  Character(MaxCharLength)     :: Name
  Character(MaxCharLength)     :: BinaryFormat
  Character(MaxCharLength)     :: FileType
  Character(MaxCharLength)     :: Prefix
  Character(MaxCharLength)     :: Suffixs(MaxNWPMetDefn2s)
  Character(MaxFileNameLength) :: TopogFile
  Logical                      :: DayPerFile
  Type(Time_)                  :: T0
  Type(Time_)                  :: Dt
  Logical                      :: NextHeatFlux
  Logical                      :: NextPrecip
  Logical                      :: NextCloud
  Character(MaxCharLength)     :: MetDefn2Names(MaxNWPMetDefn2s)
  Character(MaxCharLength)     :: ZCoordNameW
  Character(MaxCharLength)     :: ZCoordNameCl
  Character(MaxCharLength)     :: HGridName
  Character(MaxCharLength)     :: HGridNameU
  Character(MaxCharLength)     :: HGridNameV
  Character(MaxCharLength)     :: ZGridName
  Character(MaxCharLength)     :: ZGridNameUV
  Character(MaxCharLength)     :: ZGridNameW
  Character(MaxCharLength)     :: ZGridNameP
  Real(Std)                    :: SigUUM
  Real(Std)                    :: TauUUM
  Integer                      :: nMetDefn2Names
  ! Name           :: Name of met definition (part 1 - basic definitions).
  ! BinaryFormat   :: Binary format of met and topography files.
  ! FileType       :: Type of met and topography files (Name II, PP or GRIB).
  ! Prefix         :: Prefix for names of met files.
  ! Suffixs        :: Suffixes for names of met files.
  ! TopogFile      :: Topography file.
  ! DayPerFile     :: Indicates each met file contains a whole day rather than a single
  !                   time.
  ! T0             :: Reference time for the first met field.
  ! Dt             :: Time interval between fields.
  ! NextHeatFlux   :: Indicates its best to use the next time step for heat flux and ustar rather than
  !                   interpolating in time.
  ! NextPrecip     :: Indicates its best to use the next time step for precipitation rather than interpolating
  !                   in time.
  ! NextCloud      :: Indicates its best to use the next time step for cloud rather than interpolating
  !                   in time.
  ! MetDefn2Names  :: Names of the NWP met definition (part 2 - met file structure definition).
  ! ZCoordNameW    :: Name of vertical coord system of which the input vertical velocity
  !                   is the rate of change (internally vertical velocity is expressed
  !                   as rate of change of height (m agl)).
  ! ZCoordNameCl   :: Name of vertical coord system used for input cloud base and top
  !                   (internally cloud base and top are expressed as height above
  !                   ground (m)).
  ! HGridName      :: Name of main horizontal grid.
  ! HGridNameU     :: Name of horizontal grid for U.
  ! HGridNameV     :: Name of horizontal grid for V.
  ! ZGridName      :: Name of main vertical grid.
  ! ZGridNameUV    :: Name of vertical grid for U and V.
  ! ZGridNameW     :: Name of vertical grid for W.
  ! ZGridNameP     :: Name of vertical grid for input P (internally P is stored on the
  !                   main vertical grid).
  ! SigUUM         :: Default velocity variance for unresolved mesoscale motions
  ! TauUUM         :: Default Lagrangian timescale for unresolved mesoscale motions
  ! nMetDefn2Names :: Number of MetDefn2Names (i.e. part 2s) and suffixes for part 1 with Name
End Type NWPMetDefn_

!-------------------------------------------------------------------------------------------------------------

Type :: NWPMetDefn2_ ! An NWP met definition (part 2 - met file structure definition).
                     ! Together with part 1 this provides information on the nature of
                     ! the NWP met module met data.
  Character(MaxCharLength) :: Name
  Integer                  :: nFields
  Character(MaxCharLength) :: FieldNames(MaxNWPMetFields)
  Integer                  :: LowestLevels(MaxNWPMetFields)
  Integer                  :: HighestLevels(MaxNWPMetFields)
  Logical                  :: Top(MaxNWPMetFields)
  Integer                  :: FieldCodes(MaxNWPMetFields)
  Logical                  :: ThreeD(MaxNWPMetFields)
  Logical                  :: Total(MaxNWPMetFields)
  Character(MaxCharLength) :: NCFieldNames(MaxNWPMetFields)
  ! Name          :: Name of met definition (part 2 - met file structure definition).
  ! nFields       :: Number of fields in the met file.
  ! FieldNames    :: For each of the fields in the met file, the name of the corresponding NAME III field.
  ! LowestLevels  :: For each of the fields in the met file, the lowest NAME III model level.
  ! HighestLevels :: For each of the fields in the met file, the highest NAME III model level. Defined only if
  !                  Top is false.
  ! Top           :: For each of the fields in the met file, indicates that the highest NAME III model level
  !                  is the top level of the appropriate grid.
  ! FieldCodes    :: For each of the fields in the met file, the NWP field code (i.e. the code according to
  !                  the NWP model supplying the data). For Name II and PP met files this is the 'stash' code.
  !                  For GRIB met files this is the 'unique parameter identifier (paramId)'.
  ! ThreeD        :: For each of the fields in the met file, indicates that the field is part of a 3-d NWP
  !                  field. Note this means it has an NWP level index associated with it (which needn't match
  !                  the NAME III level index).
  ! Total         :: For each of the fields in the met file, indicates that the field is a total (dyn + conv) 
  !                  cloud field (rather than dynamic).
  ! NCFieldNames  :: For each of the fields in the met file, the NetCDF name of the variable.

  ! Note that here we distinguish between: (i) a 'NAME III field' - a complete 2-d or 3-d field as used within
  ! NAME III, (ii) an 'NWP field' - a complete 2-d or 3-d field as used within the NWP model supplying the
  ! data, and (iii) a 'field in the met file' - a 2-d (horizontal) field or a consecutive set of (horizontal)
  ! sections of a 3-d field.

  ! Here 'consecutive' means in the sense of their location in the met file if the format requires this and in
  ! the sense of their (NAME III) level numbers (with level numbers increasing).

  ! $$ could support consecutive decreasing level numbers?

  ! Name II and PP met files require that the order of fields in the met file matches that in NWPMetDefn2_.
  ! This constraint is not required for GRIB met files.

  ! $$ Relax this by reading files differently?

  ! Although here we refer to a field in the met file as a single entity, in fact the met file formats have a
  ! basic unit of data which is a single 2-d field or a single 2-d section of a 3-d field. All formats support
  ! the possibility of some of these basic units being missing in the met file.

  ! The same NAME III field name or the same NWP field code can appear in more than one field in the met file.
  ! This is useful for example if (i) the NAME III field's levels are not all consecutive in the met file and
  ! the format requires that the met file order matches the entries in NWPMetDefn2_, or (ii) the NAME III
  ! fields levels correspond to more than one NWP field code or vice versa (e.g. 10m wind, screen or surface
  ! temperature, surface pressure, or surface stress may be different fields to the 3-d wind, temperature,
  ! pressure and stress fields), or (iii) there is more than one 2-d NWP field or 2-d section of a 3-d NWP
  ! field which can be used for a given 2-d NAME III field or 2-d section of a 3-d NAME III field, or (iv)
  ! there is more than one 2-d NAME III field or 2-d section of a 3-d NAME III field which can use a given 2-d
  ! NWP field or 2-d section of a 3-d NWP field.

  ! Note however that where a section of a 3-d NAME III field uses a section of a 3-d NWP field, the
  ! corresponding level numbers are determined by the grid.

  ! $$ extend to allow a space separated list of NWP level indices to be specified for each field in the met
  ! file. This would make grid indices redundant or could have 'use grid indices' as an option
  ! (could then relax test that grids must have indices defined).

  ! ThreeD must be false for NAME III fields which are 2-d, but doesn't have to be true for NAME III fields
  ! which are 3-d. This is useful for example if a 2-d NWP field such as 10m wind is used as part of a 3-d
  ! NAME III field.

  ! $$ extend to allow 3-D NWP field to fill 2-D NAME III field (e.g. 3-d stress field used to generate
  ! surface stress) - but will need to give NWP level number in metdefn2. If above extension done too, can
  ! remove ThreeD (no level given implies 2-d). Also useful to allow more than one NWP level when reading into
  ! dummy array. Currently we sometimes lie in the met defn files, pretending an NWP field is 2-d when its
  ! actually a level of a 3-d field (e.g. 10m Q, Cl Water and Cl Ice in V3, H, H2001 mesoscale files). This
  ! works for ordered reads but wont work for GRIB if more than one field can meet spec.

  ! If the NWP field is 2-d then LowestLevels and HighestLevels must be equal. If the NAME III field is 2-D,
  ! then LowestLevels and HighestLevels must equal 1.

  ! $$ add tests that this is so. could set to zero instead? blank in met defn file? But note wish to read
  ! several levels into dummy field (see above).

  ! Where more than one possible NWP field can be used for a given NAME III field (or a 2-d section of a 3-d
  ! NAME III field), this can be requested by specifying more than one field in the met file that can meet the
  ! need. In this case the last NWP field in the order defined in NWPMetDefn2_ which is present in the met
  ! file and which can meet the need will be used.

  ! The same 2-d NWP field (or the same 2-d section of a 3-d NWP field) should not appear twice in the met
  ! file where this can be avoided. This is enforced by the following tests when the met files are read:
  !   Direct access indexed read (e.g. GRIB): check the same NWP field code / level combination doesn't occur
  !   twice in the file.
  !   Sequential read in order not necessarily aligned with NWPMetDefn2_: check the same NWP field code /
  !   level combination doesn't occur twice in the file.
  !   Sequential read in same order as NWPMetDefn2_: here there are some situations which may require the same
  !   field to appear twice and so no checks are made.

  ! $$ add tests required for this. May need to check ensemble index/time too.
  ! $$ for 'Sequential read in same order as NWPMetDefn2_' could possibly check any duplicate field is
  !    the same.

End Type NWPMetDefn2_

!-------------------------------------------------------------------------------------------------------------

Type :: MetDefns_ ! A collection of met definitions.
  Integer            :: nNWPMetDefns                  ! Number of NWP met definitions (part 1 - basic
                                                      ! definitions).
  Type(NWPMetDefn_)  :: NWPMetDefns(MaxNWPMetDefns)   ! NWP met definitions (part 1 - basic definitions).
  Integer            :: nNWPMetDefn2s                 ! Number of NWP met definitions (part 2 - met file
                                                      ! structure definition).
  Type(NWPMetDefn2_) :: NWPMetDefn2s(MaxNWPMetDefn2s) ! NWP met definitions (part 2 - met file structure
                                                      ! definition).
End Type MetDefns_

!-------------------------------------------------------------------------------------------------------------

Type :: NWPMet_ ! The state of an NWP met module instance.

  ! CommonMet_:
  Type(CommonMet_) :: C ! Part of met state common to all flow modules.

  ! Validity flags for each time level.
  Logical :: ValidAttribs(nAttribs, 3) ! Indicates whether, for each time level, the module is valid for each
                                       ! of its attributes.

  ! Input variables:
  Type(NWPMetDefn_)                     :: MetDefn                    ! NWP met definition (part 1 - basic
                                                                      ! definitions).
  Type(NWPMetDefn2_)                    :: MetDefn2s(MaxNWPMetDefn2s) ! NWP met definition (part 2 - met file
                                                                      ! structure definition).
  Type(ShortTime_)                      :: Dt                         ! Time interval between fields.
  Real(Std)                             :: HMin                       !} Minimum and maximum boundary layer
  Real(Std)                             :: HMax                       !} depth to be imposed.
  Logical                               :: UseNWPH                    ! Indicates that the NWP boundary layer
                                                                      ! depth is to be used if its available.
  Real(Std)                             :: SigUUM                     ! Velocity variance of unresolved mesoscale motions.
  Real(Std)                             :: TauUUM                     ! Lagrangian timescale of unresolved mesoscale motions.
  Real(Std)                             :: SigU2HPlus                 !} Free tropospheric velocity variances
  Real(Std)                             :: SigW2HPlus                 !} (horizontal and vertical).
  Real(Std)                             :: TauUHPlus                  !] Free tropospheric Lagrangian timescales
  Real(Std)                             :: TauWHPlus                  !] (horzontal and vertical).
  Integer                               :: MetFolderDefnType          ! Indicates the way the met folder(s) are
                                                                      ! defined - i.e. which of MetFolder,
                                                                      ! MetFolderStem or MetFolders is used.
                                                                      !   0 = the single met folder is valid,
                                                                      !   1 = the met folder stem is valid,
                                                                      !   2 = the array of met folders is valid.
                                                                      !   3 = the single ensemble met folder is valid.
  Character(MaxFileNameLength)          :: MetFolder                  ! Folder containing met data.
  Character(MaxFileNameLength)          :: MetFolderStem              ! Stem of folders containing met data.
  Character(MaxFileNameLength), Pointer :: MetFolders(:)              ! Array of folders containing met data.
  Character(MaxFileNameLength)          :: EnsembleMetFolder          ! Folder containing ensemble of met data.
  Character(MaxFileNameLength)          :: TopogFolder                ! Folder containing topography data.
  Character(MaxFileNameLength)          :: RestoreMetScript           ! File containing script for restoring met data.
                                                                      ! Blank indicates met files will not be restored.
  Logical                               :: DeleteMet                  ! Indicates met files will be deleted after use.

  ! Quantities derived from MetDefn2s and Grids:
  Integer :: iField(MaxNWPMetDefn2s,MaxNWPMetFields) ! Indices in FieldNames3d and FieldNames2d of
                                                     ! fields in the met file.
  Logical :: ThreeD(MaxNWPMetDefn2s,MaxNWPMetFields) ! Indicates if the NAME III field is 3-D for each of the fields in the
                                                     ! met file.
  Integer :: k1(MaxNWPMetDefn2s,MaxNWPMetFields)     ! Lowest NAME III model level of each of the fields in the
                                                     ! met file (1 for 2-D fields).
  Integer :: k2(MaxNWPMetDefn2s,MaxNWPMetFields)     ! Highest NAME III model level of each of the fields in
                                                     ! the met file (1 for 2-D fields).
  Integer :: kMax(nFieldNames3d)                     ! Highest NAME III model level of each of the fields
                                                     ! listed in FieldNames3d.
  Logical :: TotalCloudFlag                          ! Indicates that the non-convective cloud (liquid water / ice, 
                                                     ! base / height) are total cloud values

  ! Grid and coord indices:
  Integer :: iHCoord   !} Indices of the main horizontal and vertical coord systems
  Integer :: iZCoord   !} (as used for the main grids).
  Integer :: iZCoordZ  !] Indices of m agl and Pa coord systems.
  Integer :: iZCoordP  !]
  Integer :: iZCoordW  !} Other coord system indices.
  Integer :: iZCoordCl !}
  Integer :: iHGrid    !] Grid indices.
  Integer :: iHGridU   !]
  Integer :: iHGridV   !]
  Integer :: iZGrid    !]
  Integer :: iZGridUV  !]
  Integer :: iZGridW   !]
  Integer :: iZGridP   !]

  ! Interpolation coefficients:
  Type(GHCoeffs_)         :: GHCoeffs
  Type(GHCoeffs_)         :: GHCoeffsU
  Type(GHCoeffs_)         :: GHCoeffsV
  Type(GHCoeffs_)         :: GHCoeffsdX
  Type(GHCoeffs_)         :: GHCoeffsdY
  Type(ZCoeffs_), Pointer :: ZCoeffsToW(:)
  Type(ZCoeffs_), Pointer :: ZCoeffsUVToW(:)
  Type(ZCoeffs_), Pointer :: ZCoeffsToP(:)
  ! GHCoeffs     :: Coefficients for interpolating from main horizontal grid to main
  !                 horizontal grid.
  ! GHCoeffsU    :: Coefficients for interpolating from horizontal grid for U to main
  !                 horizontal grid.
  ! GHCoeffsV    :: Coefficients for interpolating from horizontal grid for V to main
  !                 horizontal grid.
  ! GHCoeffsdX   :: Coefficients for interpolating x derivatives from main horizontal
  !                 grid to main horizontal grid.
  ! GHCoeffsdY   :: Coefficients for interpolating y derivatives from main horizontal
  !                 grid to main horizontal grid.
  ! ZCoeffsToW   :: Coefficients for interpolatingfrom main vertical grid to levels
  !                 of vertical grid for W.
  ! ZCoeffsUVToW :: Coefficients for interpolating from vertical grid for U and V to
  !                 levels of vertical grid for W.
  ! ZCoeffsToP   :: Coefficients for interpolating from main vertical grid to levels
  !                 of vertical grid for P.

  ! Allocatable arrays:
  Real(Std), Allocatable :: U            (:,:,:,:)
  Real(Std), Pointer :: V            (:,:,:,:)
  Real(Std), Pointer :: W            (:,:,:,:)
  Real(Std), Pointer :: WNoTFix      (:,:,:,:)
  Real(Std), Pointer :: EtaDot       (:,:,:,:)
  Real(Std), Pointer :: EtaDotNoTFix (:,:,:,:)
  Real(Std), Pointer :: T            (:,:,:,:)
  Real(Std), Pointer :: Theta        (:,:,:,:)
  Real(Std), Pointer :: Q            (:,:,:,:)
  Real(Std), Pointer :: PAsRead      (:,:,:,:)
  Real(Std), Pointer :: P            (:,:,:,:)
  Real(Std), Pointer :: FLPa         (:,:,:)
  Real(Std), Pointer :: Rho          (:,:,:,:)
  Real(Std), Pointer :: Z            (:,:,:,:)
  Real(Std), Pointer :: Topog        (:,:,:)
  Real(Std), Pointer :: Z0           (:,:,:)
  Real(Std), Pointer :: UStress      (:,:,:)
  Real(Std), Pointer :: VStress      (:,:,:)
  Real(Std), Pointer :: UStar        (:,:,:)
  Real(Std), Pointer :: HeatFlux     (:,:,:)
  Real(Std), Pointer :: H            (:,:,:)
  Real(Std), Pointer :: PSeaLevel    (:,:,:)
  Real(Std), Pointer :: ConCloud     (:,:,:)
  Real(Std), Pointer :: ConCloudBase (:,:,:)
  Real(Std), Pointer :: ConCloudTop  (:,:,:)
  Real(Std), Pointer :: ConCloudBasePa (:,:,:)
  Real(Std), Pointer :: ConCloudTopPa  (:,:,:)
  Real(Std), Pointer :: Cloud        (:,:,:)
  Real(Std), Pointer :: Cloud3d      (:,:,:,:)
  Real(Std), Pointer :: TotalOrDynCloudWater(:,:,:,:)
  Real(Std), Pointer :: TotalOrDynCloudIce  (:,:,:,:)
  Real(Std), Pointer :: TotalOrDynCloudBase (:,:,:)
  Real(Std), Pointer :: TotalOrDynCloudTop  (:,:,:)
  Real(Std), Pointer :: TotalOrDynHCloud    (:,:,:)
  Real(Std), Pointer :: TotalOrDynMCloud    (:,:,:)
  Real(Std), Pointer :: TotalOrDynLCloud    (:,:,:)
  Real(Std), Pointer :: ConPpt       (:,:,:)
  Real(Std), Pointer :: ConRain      (:,:,:)
  Real(Std), Pointer :: ConSnow      (:,:,:)
  Real(Std), Pointer :: DynPpt       (:,:,:)
  Real(Std), Pointer :: DynRain      (:,:,:)
  Real(Std), Pointer :: DynSnow      (:,:,:)
  Real(Std), Pointer :: SoilMoisture (:,:,:)
  Real(Std), Pointer :: CanopyHeight  (:,:,:,:)
  Real(Std), Pointer :: CanopyWater   (:,:,:,:)
  Real(Std), Pointer :: StomataConduct(:,:,:,:)
  ! U             :} Wind (u, v and w components). U and V are aligned with the X and
  ! V             :} Y coords in the coord system named ZCoordName (but are in m/s),
  ! W             :} while W is read in as rate of change of the coord system named
  !                  ZCoordNameW, but is stored internally as rate of change of height
  !                  (m agl).
  ! WNoTFix       :: W array without the d/dt correction term (used only when input vertical velocities are
  !                  pressure-based).
  ! EtaDot        :: Eta dot array (1/s). That is, the rate of change of the eta coord.
  ! EtaDotNoTFix  :: Eta dot array (1/s) without the d/dt correction term (used only when input vertical
  !                  velocities are pressure-based).
  ! T             :: Temperature (K).
  ! Theta         :: Potential temperature.
  ! Q             :: Specific humidity.
  ! PAsRead       :: Pressure as read from the met file (Pa).
  ! P             :: Pressure (Pa).
  ! FLPa          :: Pressure of the level nearer to the freezing point (Pa).
  ! Rho           :: Density.
  ! Z             :: Height above ground.
  ! Topog         :: Topography (has time dimension for convenience).
  ! Z0            :: Roughness length.
  ! UStress       :} Surface stress (u and v components) (N/m^2).
  ! VStress       :}
  ! UStar         :: Friction velocity.
  ! HeatFlux      :: Surface sensible heat flux.
  ! H             :: Boundary layer depth.
  ! PSeaLevel     :: Sea level pressure (Pa).
  ! ConCloud      :: Convective cloud amount (fraction).
  ! ConCloudBase  :} Convective cloud base and top. These are read in as values in the
  ! ConCloudTop   :} coord system named ZCoordNameCl, but are stored internally as
  !                  height (m agl). Values < 0 indicate no cloud.
  ! ConCloudBasePa :: Convective cloud base in pressure coordinate
  ! ConCloudTopPa  :: Convective cloud top in pressure coordinate
  ! Cloud          :: Total cloud amount (fraction).
  ! Cloud3d        :: 3-d total cloud amount (fraction).
  ! TotalOrDynCloudWater :] Total (TotalCloudFlag true) or dynamic (TotalCloudFlag false) cloud 
  ! TotalOrDynCloudIce   :] liquid water (kg/kg) and cloud ice (kg/kg). 
  ! TotalOrDynCloudBase  :} Total or dynamic cloud base and top (m agl). Values < 0 indicate no cloud.
  ! TotalOrDynCloudTop   :}
  ! TotalOrDynHCloud     :] Total or dynamic high/medium/low cloud amount (fraction).
  ! TotalOrDynMCloud     :]
  ! TotalOrDynLCloud     :]
  ! ConPpt        :: Convective precipitation rate.
  ! ConRain       :: Convective rain rate.
  ! ConSnow       :: Convective snow rate.
  ! DynPpt        :: Dynamic precipitation rate.
  ! DynRain       :: Dynamic rain rate.
  ! DynSnow       :: Dynamic snow rate.
  ! SoilMoisture  :: Soil moisture.
  ! CanopyHeight  :: Canopy Height (m) (5 plant types)
  ! CanopyWater   :: Canopy Water (kg/m^2) (5 plant types)
  ! StomataConduct:: Stomatal conductance (m/s) (5 plant types)
  ! Note that for TotalOrDynCloudWater and TotalOrDynCloudIce we don't distinguish
  ! between mass per unit mass of dry air and mass per unit total mass. Also ConPpt,
  ! ConRain, ConSnow, DynPpt, DynRain and DynSnow are (if read in) read in in units
  ! of kg/(m^2 s) but stored internally in units of mm of water per hour. $$
  ! $$ could save memory by avoiding need for PAsRead?

  ! Information on allocatable arrays:
  Logical                  :: SpaceAllocated
  Integer                  :: New
  Integer                  :: Old
  Integer                  :: Prefetch
  Type(ShortTime_)         :: OldTime
  Type(ShortTime_)         :: NewTime
  Type(ShortTime_)         :: PrefetchTime
  Logical,         Pointer :: FieldPresent3d(:, :)
  Integer                  :: LowestPresent(nFieldNames3d)
  Integer                  :: HighestPresent(nFieldNames3d)
  Integer,         Pointer :: kLower(:, :)
  Integer,         Pointer :: kUpper(:, :)
  Integer,         Pointer :: FieldPersist3d(:, :)
  Logical                  :: FieldPresent2d(nFieldNames2d)
  Integer                  :: FieldPersist2d(nFieldNames2d)
  ! SpaceAllocated :: Indicates arrays allocated and topography read in.
  ! New            :} Indices of latest and latest but one sets of met data. 0
  ! Old            :} indicates data invalid because no data read in, errors occurred
  !                   during reading the data, or data is for the wrong time.
  ! OldTime        :] Time of latest and latest but one sets of met data.
  ! NewTime        :]
  ! FieldPresent3d :: Indicates which levels of the 3-d fields have been successfully read in.
  ! LowestPresent  :: Lowest level at which the field is present. Huge/3 values indicate no such level.
  ! HighestPresent :: Highest level at which the field is present. -Huge/3 values indicate no such level.
  ! kLower         :} Lower and higher levels used in interpolating the field to a given level.
  ! kUpper         :}
  ! FieldPersist3d :: The number of met times for which persistence has been used for the levels of the 3-d
  !                   fields. -1 indicates a fix up to a default value has been applied.
  ! FieldPresent2d :: Indicates which 2-d fields have been successfully read in.
  ! FieldPersist2d :: The number of met times for which persistence has been used for the 2-d fields. -1
  !                   indicates a fix up to a default value has been applied.

  ! $$ note persistence counters should go to restart file.

  ! Miscellaneous:
  Character(MaxFileNameLength) :: OldMetFiles(MaxNWPMetDefn2s) ! Name of the last opened met files. Blank
                                                               ! indicates undefined.
  Type(OpenMPOpts_) :: OpenMPOpts   ! OpenMP options

End Type NWPMet_

!-------------------------------------------------------------------------------------------------------------

! Skipping of MetData files
Integer, Parameter :: MaxSkip = 1 ! MaxSkip = maximal number of files we allow to be skipped

!-------------------------------------------------------------------------------------------------------------

Type(Timer_), Save :: WorkerNWPMetReadTimer(0:OpenMPMaxThreads-1), &
                      WorkerNWPMetProcessTimer(0:OpenMPMaxThreads-1)
Logical            :: TimersInitialised = .false.

!-------------------------------------------------------------------------------------------------------------

Interface Operator(==) ! Equality of NWP met definitions.
  Module Procedure NWPMetDefnEq
  Module Procedure NWPMetDefn2Eq
End Interface

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

! Initialise used in this module

Subroutine NWPMetModuleTimerInitialise(TimerOpts, OpenMPOpts)

  Implicit None

  Type(TimerOpts_),  Intent(In) :: TimerOpts
  Type(OpenMPOpts_), Intent(In) :: OpenMPOpts
  Integer                       :: i

  If (.not.TimersInitialised) Then
    Do i = 0, OpenMPOpts%nParticleThreads-1
      Call TimerCreate(WorkerNWPMetReadTimer(i)   , "WorkerThread"    // &
                                                    Trim(Int2Char(i)) // &
                                                    "_NWPMetRead",       &
                                                    TimerOpts)
      Call TimerCreate(WorkerNWPMetProcessTimer(i), "WorkerThread"    // &
                                                    Trim(Int2Char(i)) // &
                                                    "_NWPMetProcess",    &
                                                    TimerOpts)
    End Do
    TimersInitialised = .true.
  End If

End Subroutine NWPMetModuleTimerInitialise

!-------------------------------------------------------------------------------------------------------------

! Output timer summary information

Subroutine NWPMetModuleTimerSummary(OpenMPOpts)

  Implicit None

  Type(OpenMPOpts_), Intent(In) :: OpenMPOpts

  Integer :: i

  If (TimersInitialised) Then
    Do i = 0, OpenMPOpts%nParticleThreads-1
      Call TimerWriteSummary(WorkerNWPMetReadTimer(i))
      Call TimerWriteSummary(WorkerNWPMetProcessTimer(i))
    End Do
  End If

End Subroutine NWPMetModuleTimerSummary

!-------------------------------------------------------------------------------------------------------------

Function InitMetDefns() Result(MetDefns)
! Initialises a collection of met definitions.

  Implicit None
  ! Function result:
  Type(MetDefns_) :: MetDefns ! Initialised collection of met definitions.

  MetDefns%nNWPMetDefns  = 0
  MetDefns%nNWPMetDefn2s = 0

End Function InitMetDefns

!-------------------------------------------------------------------------------------------------------------

Function InitNWPMetDefn(                                                         &
           Name, BinaryFormat, FileType, Prefix, Suffixs, TopogFile, DayPerFile, &
           T0, Dt, NextHeatFlux, NextPrecip, NextCloud,                          &
           MetDefn2Names,                                                        &
           ZCoordNameW, ZCoordNameCl,                                            &
           HGridName, HGridNameU,  HGridNameV,                                   &
           ZGridName, ZGridNameUV, ZGridNameW, ZGridNameP,                       &
           SigUUM, TauUUM                                                        &
         )                                                                       &
Result(NWPMetDefn)
! Initialises an NWP met definition (part 1 - basic definitions).

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name             ! Name of met definition (part 1 - basic definitions).
  Character(*), Intent(In) :: BinaryFormat     ! Binary format of met and topography files.
  Character(*), Intent(In) :: FileType         ! Type of met and topography files (Name II, PP or GRIB).
  Character(*), Intent(In) :: Prefix           ! Prefix for names of met files.
  Character(MaxCharLength), Intent(In) :: Suffixs(:)       ! Suffix for names of met files.
                                                           ! $$ changed from * to support Intel compiler on Cray
  Character(*), Intent(In) :: TopogFile        ! Topography file.
  Logical,      Intent(In) :: DayPerFile       ! Indicates each met file contains a whole
                                               ! day rather than a single time.
  Type(Time_),  Intent(In) :: T0               ! Reference time for met fields.
  Type(Time_),  Intent(In) :: Dt               ! Time interval between fields.
  Logical,      Intent(In) :: NextHeatFlux     ! Indicates its best to use the next time step for heat flux and
                                               ! ustar rather than interpolating in time.
  Logical,      Intent(In) :: NextPrecip       ! Indicates its best to use the next time step for precipitation
                                               ! rather than interpolating in time.
  Logical,      Intent(In) :: NextCloud        ! Indicates its best to use the next time step for cloud
                                               ! rather than interpolating in time.
  Character(MaxCharLength), Intent(In) :: MetDefn2Names(:) ! Name of NWP met definition (part 2 - met
                                               ! file structure definition).
                                               ! $$ changed from * to support Intel compiler on Cray
  Character(*), Intent(In) :: ZCoordNameW      ! Name of vertical coord system of which
                                               ! the input vertical velocity is the rate
                                               ! of change (internally vertical velocity
                                               ! is expressed as rate of change of height
                                               ! above ground (m)).
  Character(*), Intent(In) :: ZCoordNameCl     ! Name of vertical coord system used for
                                               ! input cloud base and top (internally
                                               ! cloud base and top are expressed as
                                               ! height (m agl)).
  Character(*), Intent(In) :: HGridName        ! Name of main horizontal grid.
  Character(*), Intent(In) :: HGridNameU       ! Name of horizontal grid for U.
  Character(*), Intent(In) :: HGridNameV       ! Name of horizontal grid for V.
  Character(*), Intent(In) :: ZGridName        ! Name of main vertical grid.
  Character(*), Intent(In) :: ZGridNameUV      ! Name of vertical grid for U and V.
  Character(*), Intent(In) :: ZGridNameW       ! Name of vertical grid for W.
  Character(*), Intent(In) :: ZGridNameP       ! Name of vertical grid for input P
                                               ! (internally P is stored on the main
                                               ! vertical grid).
  Real(Std),    Intent(In) :: SigUUM           ! Default velocity variance for unresolved mesoscale motions.
  Real(Std),    Intent(In) :: TauUUM           ! Default Lagrangian timescale for unresolved mesoscale motions.
  ! Function result:
  Type(NWPMetDefn_) :: NWPMetDefn ! The NWP met definition (part 1 - basic
                                  ! definitions).

  ! Locals:
  Integer :: iMetDefn2Names    ! Loop index.

  If (Len_Trim(Name) == 0) Then
    Call Message('ERROR in InitNWPMetDefn: Name is blank', 3)
  End If
  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message(                                           &
           'ERROR in InitNWPMetDefn: Name is given as "' // &
           Trim(Name)                                    // &
           '" and is too long',                             &
           3                                                &
         )
  End If

  If (Len_Trim(FileType) == 0) Then
    Call Message('ERROR in InitNWPMetDefn: File Type is blank', 3)
  End If
  If (Len_Trim(FileType) > MaxCharLength) Then
    Call Message(                                                 &
           'ERROR in InitNWPMetDefn:  File Type is given as "' // &
           Trim(FileType)                                      // &
           '" and is too long',                                   &
           3                                                      &
         )
  End If

  If (Len_Trim(Prefix) > MaxCharLength) Then
    Call Message(                                             &
           'ERROR in InitNWPMetDefn: Prefix is given as "' // &
           Trim(Prefix)                                    // &
           '" and is too long',                               &
           3                                                  &
         )
  End If

  If (Size(Suffixs) /= Size(MetDefn2Names)) Then
    Call Message(                                                &
           'ERROR in InitNWPMetDefn: Size of Suffixs array '  // &
           'is different to the size of MetDefn2Names array',    &
           3                                                     &
         )
  End If

  If (Len_Trim(TopogFile) == 0) Then
    Call Message('ERROR in InitNWPMetDefn: TopogFile is blank', 3)
  End If
  If (Len_Trim(TopogFile) > MaxCharLength) Then
    Call Message(                                                &
           'ERROR in InitNWPMetDefn: TopogFile is given as "' // &
           Trim(TopogFile)                                    // &
           '" and is too long',                                  &
           3                                                     &
         )
  End If


  Do iMetDefn2Names = 1,Size(MetDefn2Names)
    If (Len_Trim(Suffixs(iMetDefn2Names)) > MaxCharLength) Then
      Call Message(                                 &
             'ERROR in InitNWPMetDefn: Element ' // &
             Trim(Int2Char(iMetDefn2Names))      // &
             ' of Suffixs is given as "'         // &
             Trim(Suffixs(iMetDefn2Names))       // &
             '" and is too long',                   &
             3                                      &
           )
    End If

    If (Len_Trim(MetDefn2Names(iMetDefn2Names)) == 0) Then
      Call Message(                                 &
             'ERROR in InitNWPMetDefn: Element ' // &
             Trim(Int2Char(iMetDefn2Names))      // &
             ' of MetDefn2Names is blank',          &
             3                                      &
           )
    End If
    If (Len_Trim(MetDefn2Names(iMetDefn2Names)) > MaxCharLength) Then
      Call Message(                                 &
             'ERROR in InitNWPMetDefn: Element ' // &
             Trim(Int2Char(iMetDefn2Names))      // &
             ' of MetDefn2Names is given as "'   // &
             Trim(MetDefn2Names(iMetDefn2Names)) // &
             '" and is too long',                   &
             3                                      &
           )
    End If
  End Do

  If (Len_Trim(ZCoordNameW) == 0) Then
    Call Message('ERROR in InitNWPMetDefn: ZCoordNameW is blank', 3)
  End If
  If (Len_Trim(ZCoordNameW) > MaxCharLength) Then
    Call Message(                                                  &
           'ERROR in InitNWPMetDefn: ZCoordNameW is given as "' // &
           Trim(ZCoordNameW)                                    // &
           '" and is too long',                                    &
           3                                                       &
         )
  End If

  If (Len_Trim(ZCoordNameCl) == 0) Then
    Call Message('ERROR in InitNWPMetDefn: ZCoordNameCl is blank', 3)
  End If
  If (Len_Trim(ZCoordNameCl) > MaxCharLength) Then
    Call Message(                                                   &
           'ERROR in InitNWPMetDefn: ZCoordNameCl is given as "' // &
           Trim(ZCoordNameCl)                                    // &
           '" and is too long',                                     &
           3                                                        &
         )
  End If

  If (Len_Trim(HGridName) == 0) Then
    Call Message('ERROR in InitNWPMetDefn: HGridName is blank', 3)
  End If
  If (Len_Trim(HGridName) > MaxCharLength) Then
    Call Message(                                                &
           'ERROR in InitNWPMetDefn: HGridName is given as "' // &
           Trim(HGridName)                                    // &
           '" and is too long',                                  &
           3                                                     &
         )
  End If

  If (Len_Trim(HGridNameU) == 0) Then
    Call Message('ERROR in InitNWPMetDefn: HGridNameU is blank', 3)
  End If
  If (Len_Trim(HGridNameU) > MaxCharLength) Then
    Call Message(                                                 &
           'ERROR in InitNWPMetDefn: HGridNameU is given as "' // &
           Trim(HGridNameU)                                    // &
           '" and is too long',                                   &
           3                                                      &
         )
  End If

  If (Len_Trim(HGridNameV) == 0) Then
    Call Message('ERROR in InitNWPMetDefn: HGridNameV is blank', 3)
  End If
  If (Len_Trim(HGridNameV) > MaxCharLength) Then
    Call Message(                                                 &
           'ERROR in InitNWPMetDefn: HGridNameV is given as "' // &
           Trim(HGridNameV)                                    // &
           '" and is too long',                                   &
           3                                                      &
         )
  End If

  If (Len_Trim(ZGridName) == 0) Then
    Call Message('ERROR in InitNWPMetDefn: ZGridName is blank', 3)
  End If
  If (Len_Trim(ZGridName) > MaxCharLength) Then
    Call Message(                                                &
           'ERROR in InitNWPMetDefn: ZGridName is given as "' // &
           Trim(ZGridName)                                    // &
           '" and is too long',                                  &
           3                                                     &
         )
  End If

  If (Len_Trim(ZGridNameUV) == 0) Then
    Call Message('ERROR in InitNWPMetDefn: ZGridNameUV is blank', 3)
  End If
  If (Len_Trim(ZGridNameUV) > MaxCharLength) Then
    Call Message(                                                  &
           'ERROR in InitNWPMetDefn: ZGridNameUV is given as "' // &
           Trim(ZGridNameUV)                                    // &
           '" and is too long',                                    &
           3                                                       &
         )
  End If

  If (Len_Trim(ZGridNameW) == 0) Then
    Call Message('ERROR in InitNWPMetDefn: ZGridNameW is blank', 3)
  End If
  If (Len_Trim(ZGridNameW) > MaxCharLength) Then
    Call Message(                                                 &
           'ERROR in InitNWPMetDefn: ZGridNameW is given as "' // &
           Trim(ZGridNameW)                                    // &
           '" and is too long',                                   &
           3                                                      &
         )
  End If

  If (Len_Trim(ZGridNameP) == 0) Then
    Call Message('ERROR in InitNWPMetDefn: ZGridNameP is blank', 3)
  End If
  If (Len_Trim(ZGridNameP) > MaxCharLength) Then
    Call Message(                                                 &
           'ERROR in InitNWPMetDefn: ZGridNameP is given as "' // &
           Trim(ZGridNameP)                                    // &
           '" and is too long',                                   &
           3                                                      &
         )
  End If

  NWPMetDefn%Name           = Name
  NWPMetDefn%DayPerFile     = DayPerFile
  NWPMetDefn%Prefix         = Prefix
  NWPMetDefn%TopogFile      = TopogFile
  NWPMetDefn%NextHeatFlux   = NextHeatFlux
  NWPMetDefn%NextPrecip     = NextPrecip
  NWPMetDefn%NextCloud      = NextCloud
  NWPMetDefn%ZCoordNameW    = ZCoordNameW
  NWPMetDefn%ZCoordNameCl   = ZCoordNameCl
  NWPMetDefn%HGridName      = HGridName
  NWPMetDefn%HGridNameU     = HGridNameU
  NWPMetDefn%HGridNameV     = HGridNameV
  NWPMetDefn%ZGridName      = ZGridName
  NWPMetDefn%ZGridNameUV    = ZGridNameUV
  NWPMetDefn%ZGridNameW     = ZGridNameW
  NWPMetDefn%ZGridNameP     = ZGridNameP
  NWPMetDefn%nMetDefn2Names = Size(MetDefn2Names)
  NWPMetDefn%SigUUM         = SigUUM
  NWPMetDefn%TauUUM         = TauUUM

  NWPMetDefn%Suffixs      (1:NWPMetDefn%nMetDefn2Names)       = Suffixs(:)
  NWPMetDefn%MetDefn2Names(1:NWPMetDefn%nMetDefn2Names)       = MetDefn2Names(:)

  ! Binary Format.
  If (Len_Trim(BinaryFormat) > MaxCharLength) Then
    Call Message(                                                   &
           'ERROR in InitNWPMetDefn: BinaryFormat is given as "' // &
           Trim(BinaryFormat)                                    // &
           '" and is too long',                                     &
           3                                                        &
         )
  End If
  NWPMetDefn%BinaryFormat = BinaryFormat
  ! Native format.
  If (NWPMetDefn%BinaryFormat == ' ') NWPMetDefn%BinaryFormat = 'NATIVE'
# ifdef CompaqPCCompiler
    ! Check supported formats. ! $$ Compaq fortran can support more formats
    If (                                            &
      NWPMetDefn%BinaryFormat /= 'BIG_ENDIAN' .and. &
      NWPMetDefn%BinaryFormat /= 'NATIVE'           &
    ) Then
      Call Message(                                                                 &
             'ERROR in InitNWPMetDefn: '                                         // &
             'only BIG_ENDIAN and NATIVE supported under Compaq NT PC compiler',    &
             3                                                                      &
           )
    End If
# endif
# ifdef IntelLinCompiler
    ! Check supported formats.
    ! $$ currently different formats obtained by environmental variable.
    !    Value of BinaryFormat is irrelevant and so isn't checked here.
# endif
# ifdef sun
    ! Check supported formats.
    ! $$ currently different formats obtained by compile line xfilebyteorder option.
    !    Value of BinaryFormat is irrelevant and so isn't checked here.
    !    Can run on SPARC or Intel hardware.
    !    SUN SPARC / Unix native format is BIG_ENDIAN, Intel is little endian.
# endif

  ! File Type.
  If (                                     &
    .not.(FileType .CIEq. 'Name II') .and. &
    .not.(FileType .CIEq. 'PP     ') .and. &
    .not.(FileType .CIEq. 'NetCDF ') .and. &
    .not.(FileType .CIEq. 'GRIB   ')       &
  ) Then
    Call Message(                                            &
           'FATAL ERROR: NWP met file type is given as "' // &
           Trim(FileType)                                 // &
           '" and has not been recognised',                  &
           3                                                 &
         )
  End If
# ifdef GRIBsupport
# else
    If (FileType .CIEq. 'GRIB') Then
      Call Message(                                                                                  &
             'FATAL ERROR: This version of the code does not support GRIB format for NWP met files', &
             3                                                                                       &
           )
    End If
# endif
# ifdef NetCDFsupport
# else
    If (FileType .CIEq. 'NetCDF') Then
      Call Message(                                                                                    &
             'FATAL ERROR: This version of the code does not support NetCDF format for NWP met files', &
             3                                                                                         &
           )
    End If
# endif
  NWPMetDefn%FileType = FileType

  ! T0.
  NWPMetDefn%T0 = T0

  ! Dt.
  If (Dt <= ZeroTime()) Then
    Call Message('ERROR in InitNWPMetDefn: Dt is <= 0', 3)
  End If
  NWPMetDefn%Dt = Dt

End Function InitNWPMetDefn

!-------------------------------------------------------------------------------------------------------------

Subroutine AddNWPMetDefn(NWPMetDefn, MetDefns)
! Adds an NWP met definition (part 1 - basic definitions) to a collection of met
! definitions.

  Implicit None
  ! Argument list:
  Type(NWPMetDefn_), Intent(In)    :: NWPMetDefn ! The NWP met definition (part 1 -
                                                 ! basic definitions).
  Type(MetDefns_),   Intent(InOut) :: MetDefns   ! The collection of met definitions.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, MetDefns%nNWPMetDefns
    If (NWPMetDefn%Name .CIEq. MetDefns%NWPMetDefns(i)%Name) Then
      If (NWPMetDefn == MetDefns%NWPMetDefns(i)) Then
        Return
      Else
        Call Message(                                                            &
               'ERROR in adding NWP Met Definition "'                         // &
               Trim(NWPMetDefn%Name)                                          // &
               '": a different definition with the same name already exists',    &
               3                                                                 &
             )
      End If
    End If
  End Do

  If (MetDefns%nNWPMetDefns >= MaxNWPMetDefns) Then
    Call Message(                                        &
           'ERROR in adding an NWP Met Definition: '  // &
           'there are too many NWP Met Definitions ',    &
           3                                             &
         )
  End If

  MetDefns%nNWPMetDefns                       = MetDefns%nNWPMetDefns + 1
  MetDefns%NWPMetDefns(MetDefns%nNWPMetDefns) = NWPMetDefn

End Subroutine AddNWPMetDefn

!-------------------------------------------------------------------------------------------------------------

Function NWPMetDefnEq(NWPMetDefn1, NWPMetDefn2)
! Tests for equality of NWP met definitions (part 1 - basic definitions).

  Implicit None
  ! Argument list:
  Type(NWPMetDefn_), Intent(In) :: NWPMetDefn1 !} The two NWP met definitions (part 1
  Type(NWPMetDefn_), Intent(In) :: NWPMetDefn2 !} - basic definitions).
  ! Function result:
  Logical :: NWPMetDefnEq ! Indicates if NWP met definitions are equal.
  ! Locals:
  Integer :: i ! Loop index.

  NWPMetDefnEq = (NWPMetDefn1%Name          .CIEq. NWPMetDefn2%Name          ) .and. &
                 (NWPMetDefn1%BinaryFormat  .CIEq. NWPMetDefn2%BinaryFormat  ) .and. &
                 (NWPMetDefn1%FileType      .CIEq. NWPMetDefn2%FileType      ) .and. &
                 (NWPMetDefn1%Prefix        .CIEq. NWPMetDefn2%Prefix        ) .and. &
                 (NWPMetDefn1%TopogFile     .CIEq. NWPMetDefn2%TopogFile     ) .and. &
                 (NWPMetDefn1%DayPerFile    .eqv.  NWPMetDefn2%DayPerFile    ) .and. &
                 (NWPMetDefn1%T0              ==   NWPMetDefn2%T0            ) .and. &
                 (NWPMetDefn1%Dt              ==   NWPMetDefn2%Dt            ) .and. &
                 (NWPMetDefn1%NextHeatFlux  .eqv.  NWPMetDefn2%NextHeatFlux  ) .and. &
                 (NWPMetDefn1%NextPrecip    .eqv.  NWPMetDefn2%NextPrecip    ) .and. &
                 (NWPMetDefn1%NextCloud     .eqv.  NWPMetDefn2%NextCloud     ) .and. &
                 (NWPMetDefn1%ZCoordNameW   .CIEq. NWPMetDefn2%ZCoordNameW   ) .and. &
                 (NWPMetDefn1%ZCoordNameCl  .CIEq. NWPMetDefn2%ZCoordNameCl  ) .and. &
                 (NWPMetDefn1%HGridName     .CIEq. NWPMetDefn2%HGridName     ) .and. &
                 (NWPMetDefn1%HGridNameU    .CIEq. NWPMetDefn2%HGridNameU    ) .and. &
                 (NWPMetDefn1%HGridNameV    .CIEq. NWPMetDefn2%HGridNameV    ) .and. &
                 (NWPMetDefn1%ZGridName     .CIEq. NWPMetDefn2%ZGridName     ) .and. &
                 (NWPMetDefn1%ZGridNameUV   .CIEq. NWPMetDefn2%ZGridNameUV   ) .and. &
                 (NWPMetDefn1%ZGridNameW    .CIEq. NWPMetDefn2%ZGridNameW    ) .and. &
                 (NWPMetDefn1%ZGridNameP    .CIEq. NWPMetDefn2%ZGridNameP    ) .and. &
                 (NWPMetDefn1%nMetDefn2Names  ==   NWPMetDefn2%nMetDefn2Names)

  Do i = 1, Min(NWPMetDefn1%nMetDefn2Names, NWPMetDefn2%nMetDefn2Names)
    NWPMetDefnEq = NWPMetDefnEq                                                       .and. &
                   (NWPMetDefn1%Suffixs(i)       .CIEq. NWPMetDefn2%Suffixs(i)      ) .and. &
                   (NWPMetDefn1%MetDefn2Names(i) .CIEq. NWPMetDefn2%MetDefn2Names(i))
  End Do

End Function NWPMetDefnEq

!-------------------------------------------------------------------------------------------------------------

Function FindNWPMetDefnIndex(Name, MetDefns)
! Finds the index of an NWP met definition (part 1 - basic definitions).

  Implicit None
  ! Argument list:
  Character(*),    Intent(In) :: Name     ! Name of NWP met definition (part 1 - basic
                                          ! definitions).
  Type(MetDefns_), Intent(In) :: MetDefns ! Collection of met definitions.
  ! Function result:
  Integer :: FindNWPMetDefnIndex ! Index of NWP met definition (part 1 - basic
                                 ! definitions).
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, MetDefns%nNWPMetDefns
    If (Name .CIEq. MetDefns%NWPMetDefns(i)%Name) Then
      FindNWPMetDefnIndex = i
      Return
    End If
  End Do

  Call Message(                                 &
         'FATAL ERROR: NWP Met Definition "' // &
         Trim(Name)                          // &
         '" not found',                         &
         3                                      &
       )

End Function FindNWPMetDefnIndex

!-------------------------------------------------------------------------------------------------------------

Function InitNWPMetDefn2(                                                                              &
           Name, FieldNames, LowestLevels, HighestLevels, Top, FieldCodes, ThreeD, Total, NCFieldNames &
         )                                                                                             &
Result(NWPMetDefn2)
! Initialises an NWP met definition (part 2 - met file structure definition).

  Implicit None
  ! Argument list:
  Character(*),              Intent(In) :: Name
  Character(MaxTokenLength), Intent(In) :: FieldNames(:)    ! $$ changed from * to support Intel compiler on Cray
  Integer,                   Intent(In) :: LowestLevels(:)
  Integer,                   Intent(In) :: HighestLevels(:)
  Logical,                   Intent(In) :: Top(:)
  Integer,                   Intent(In) :: FieldCodes(:)
  Logical,                   Intent(In) :: ThreeD(:)
  Logical,                   Intent(In) :: Total(:)
  Character(MaxTokenLength), Intent(In) :: NCFieldNames(:)  ! $$ changed from * to support Intel compiler on Cray
  ! Name          :: Name of met definition (part 2 - met file structure definition).
  ! FieldNames    :: For each of the fields in the met file, the name of the corresponding NAME III field.
  ! LowestLevels  :: For each of the fields in the met file, the lowest NAME III model level.
  ! HighestLevels :: For each of the fields in the met file, the highest NAME III model level. Needs to be
  !                  defined only if Top is false.
  ! Top           :: For each of the fields in the met file, indicates that the highest NAME III model level
  !                  is the top level of the appropriate grid.
  ! FieldCodes    :: For each of the fields in the met file, the NWP field code (i.e. the code according to
  !                  the NWP model supplying the data). For Name II and PP met files these are the 'stash'
  !                  codes. For GRIB met files these are the 'unique parameter identifiers (paramId)'.
  ! ThreeD        :: For each of the fields in the met file, indicates that the field is part of a 3-d NWP
  !                  field. Note this means it has an NWP level index associated with it (which needn't match
  !                  the NAME III level index).
  ! Total         :: For each of the fields in the met file, indicates that the field is a total (dyn + conv) 
  !                  cloud field.
  ! NCFieldNames  :: For each of the fields in the met file, the NetCDF name of the variable.
  ! Function result:
  Type(NWPMetDefn2_) :: NWPMetDefn2 ! The NWP met definition (part 2 - met file structure definition).
  ! Locals:
  Integer :: i ! Loop index.

  If (Len_Trim(Name) == 0) Then
    Call Message(                                                            &
           'ERROR in initialising an NWP Met File Structure Definition: ' // &
           'the NWP Met File Structure Definition name is blank',            &
           3                                                                 &
         )
  End If
  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message(                                                            &
           'ERROR in initialising an NWP Met File Structure Definition: ' // &
           'the NWP Met File Structure Definition name is given as '      // &
           Trim(Name)                                                     // &
           ' and is too long',                                               &
           3                                                                 &
         )
  End If

  If (                                           &
    Size(FieldNames) /= Size(LowestLevels)  .or. &
    Size(FieldNames) /= Size(HighestLevels) .or. &
    Size(FieldNames) /= Size(Top)           .or. &
    Size(FieldNames) /= Size(FieldCodes)    .or. &
    Size(FieldNames) /= Size(ThreeD)        .or. &
    Size(FieldNames) /= Size(Total)         .or. &
    Size(FieldNames) /= Size(NCFieldNames)       &
  ) Then
    Call Message(                                                           &
           'UNEXPECTED ERROR in initialising an NWP Met File Structure ' // &
           'Definition: the array sizes are inconsistent',                  &
           3                                                                &
         )
  End If
  If (Size(FieldNames) > MaxNWPMetFields) Then
    Call Message(                                                &
           'ERROR in initialising an NWP Met File Structure ' // &
           'Definition: the number of fields is too large',      &
           3                                                     &
         )
  End If

  Do i = 1, Size(FieldNames)
    If (Len_Trim(FieldNames(i)) == 0) Then
      Call Message(                                                            &
             'ERROR in initialising an NWP Met File Structure Definition: ' // &
             'the met file structure definition "'                          // &
             Trim(Name)                                                     // &
             '" contains an empty field name (field '                       // &
             Trim(Int2Char(i))                                              // &
             ')',                                                              &
             3                                                                 &
           )
    End If
    If (Len_Trim(FieldNames(i)) > MaxCharLength) Then
      Call Message(                                                            &
             'ERROR in initialising an NWP Met File Structure Definition: ' // &
             'the met file structure definition "'                          // &
             Trim(Name)                                                     // &
             '" contains an invalid field name (field '                     // &
             Trim(Int2Char(i))                                              // &
             ' has the name '                                               // &
             Trim(FieldNames(i))                                            // &
             ' which is too long)',                                            &
             3                                                                 &
           )
    End If
    If (Len_Trim(NCFieldNames(i)) > MaxCharLength) Then
      Call Message(                                                            &
             'ERROR in initialising an NWP Met File Structure Definition: ' // &
             'the met file structure definition "'                          // &
             Trim(Name)                                                     // &
             '" contains an invalid NC field name (field '                  // &
             Trim(Int2Char(i))                                              // &
             ' has the name '                                               // &
             Trim(NCFieldNames(i))                                          // &
             ' which is too long)',                                            &
             3                                                                 &
           )
    End If
  End Do

  NWPMetDefn2%Name                              = Name
  NWPMetDefn2%nFields                           = Size(FieldNames)
  NWPMetDefn2%FieldNames   (1:Size(FieldNames)) = FieldNames   (1:Size(FieldNames))
  NWPMetDefn2%LowestLevels (1:Size(FieldNames)) = LowestLevels (1:Size(FieldNames))
  NWPMetDefn2%HighestLevels(1:Size(FieldNames)) = HighestLevels(1:Size(FieldNames))
  NWPMetDefn2%Top          (1:Size(FieldNames)) = Top          (1:Size(FieldNames))
  NWPMetDefn2%FieldCodes   (1:Size(FieldNames)) = FieldCodes   (1:Size(FieldNames))
  NWPMetDefn2%ThreeD       (1:Size(FieldNames)) = ThreeD       (1:Size(FieldNames))
  NWPMetDefn2%Total        (1:Size(FieldNames)) = Total        (1:Size(FieldNames))
  NWPMetDefn2%NCFieldNames (1:Size(FieldNames)) = NCFieldNames (1:Size(FieldNames))

End Function InitNWPMetDefn2

!-------------------------------------------------------------------------------------------------------------

Subroutine AddNWPMetDefn2(NWPMetDefn2, MetDefns)
! Adds an NWP met definition (part 2 - met file structure definition) to a
! collection of met definitions.

  Implicit None
  ! Argument list:
  Type(NWPMetDefn2_), Intent(In)    :: NWPMetDefn2 ! The NWP met definition (part 2 - met file structure
                                                   ! definition).
  Type(MetDefns_),    Intent(InOut) :: MetDefns    ! The collection of met definitions.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, MetDefns%nNWPMetDefn2s
    If (NWPMetDefn2%Name .CIEq. MetDefns%NWPMetDefn2s(i)%Name) Then
      If (NWPMetDefn2 == MetDefns%NWPMetDefn2s(i)) Then
        Return
      Else
        Call Message(                                                          &
               'ERROR in adding NWP Met File Structure Definition "'        // &
               Trim(NWPMetDefn2%Name)                                       // &
               '": a different NWP Met File Structure Definition with the ' // &
               'same name already exists',                                     &
               3                                                               &
             )
      End If
    End If
  End Do

  If (MetDefns%nNWPMetDefn2s >= MaxNWPMetDefn2s) Then
    Call Message(                                                      &
           'ERROR in adding an NWP Met File Structure Definition: ' // &
           'there are too many NWP Met File Structure Definitions',    &
           3                                                           &
         )
  End If

  MetDefns%nNWPMetDefn2s                        = MetDefns%nNWPMetDefn2s + 1
  MetDefns%NWPMetDefn2s(MetDefns%nNWPMetDefn2s) = NWPMetDefn2

End Subroutine AddNWPMetDefn2

!-------------------------------------------------------------------------------------------------------------

Function NWPMetDefn2Eq(NWPMetDefn1, NWPMetDefn2)
! Tests for equality of NWP met definitions (part 2 - met file structure definitions).

  Implicit None
  ! Argument list:
  Type(NWPMetDefn2_), Intent(In) :: NWPMetDefn1 !} The two NWP met definitions (part 2
  Type(NWPMetDefn2_), Intent(In) :: NWPMetDefn2 !} - met file structure definitions).
  ! Function result:
  Logical :: NWPMetDefn2Eq ! Indicates if NWP met definitions are equal.
  ! Locals:
  Integer :: i ! Loop index.

  NWPMetDefn2Eq = (NWPMetDefn1%Name    .CIEq. NWPMetDefn2%Name)   .and. &
                   NWPMetDefn1%nFields   ==   NWPMetDefn2%nFields
  Do i = 1, Min(NWPMetDefn1%nFields, NWPMetDefn2%nFields)
    NWPMetDefn2Eq =                                                            &
      NWPMetDefn2Eq                                                      .and. &
      (NWPMetDefn1%FieldNames(i)    .CIEq. NWPMetDefn2%FieldNames(i)   ) .and. &
      (NWPMetDefn1%LowestLevels(i)    ==   NWPMetDefn2%LowestLevels(i) ) .and. &
      (NWPMetDefn1%HighestLevels(i)   ==   NWPMetDefn2%HighestLevels(i)) .and. &
      (NWPMetDefn1%Top(i)           .eqv.  NWPMetDefn2%Top(i)          ) .and. &
      (NWPMetDefn1%FieldCodes(i)      ==   NWPMetDefn2%FieldCodes(i)   ) .and. &
      (NWPMetDefn1%ThreeD(i)        .eqv.  NWPMetDefn2%ThreeD(i)       ) .and. &
      (NWPMetDefn1%Total(i)         .eqv.  NWPMetDefn2%Total(i)        ) .and. &
      (NWPMetDefn1%NCFieldNames(i)  .CIEq. NWPMetDefn2%NCFieldNames(i) )
  End Do

End Function NWPMetDefn2Eq

!-------------------------------------------------------------------------------------------------------------

Function FindNWPMetDefn2Index(Name, MetDefns)
! Finds the index of an NWP met definition (part 2 - met file structure definition).

  Implicit None
  ! Argument list:
  Character(*),    Intent(In) :: Name     ! Name of NWP met definition (part 2 - met
                                          ! file structure definition).
  Type(MetDefns_), Intent(In) :: MetDefns ! Collection of met definitions.
  ! Function result:
  Integer :: FindNWPMetDefn2Index ! Index of NWP met definition (part 2 - met file
                                  ! structure definition).
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, MetDefns%nNWPMetDefn2s
    If (Name .CIEq. MetDefns%NWPMetDefn2s(i)%Name) Then
      FindNWPMetDefn2Index = i
      Return
    End If
  End Do

  Call Message(                                                &
         'FATAL ERROR: NWP Met File Structure Definition "' // &
         Trim(Name)                                         // &
         '" not found',                                        &
         3                                                     &
       )

End Function FindNWPMetDefn2Index

!-------------------------------------------------------------------------------------------------------------

Function InitNWPMet(                    &
           MetName,                     &
           FixedMet,                    &
           UpdateOnDemand,              &
           NWPMetDefn, NWPMetDefn2s,    &
           HMin, HMax, UseNWPH,         &
           SigUUM, TauUUM,              &
           SigU2HPlus, SigW2HPlus,      &
           TauUHPlus, TauWHPlus,        &
           MetFolder,                   &
           MetFolderStem,               &
           MetFolders,                  &
           EnsembleMetFolder,           &
           TopogFolder,                 &
           RestoreMetScript, DeleteMet, &
           OpenMPOpts                   &
         )                              &
Result(NWPMet)
! Initialises the state of an NWP met module instance.

  Implicit None
  ! Argument list:
  Character(*),       Intent(In) :: MetName           ! Name of met module instance.
  Logical,            Intent(In) :: FixedMet          ! Indicates that the met is fixed.
  Logical,            Intent(In) :: UpdateOnDemand    ! Indicates the met module instance is to be updated
                                                      ! using update-on-demand.
  Type(NWPMetDefn_),  Intent(In) :: NWPMetDefn        ! NWP met definition (part 1 - basic definitions).
  Type(NWPMetDefn2_), Intent(In) :: NWPMetDefn2s(:)   ! NWP met definition (part 2 - met file structure
                                                      ! definition).
  Real(Std),          Intent(In) :: HMin              !} Minimum and maximum boundary layer depth to be
  Real(Std),          Intent(In) :: HMax              !} imposed.
  Logical,            Intent(In) :: UseNWPH           ! Indicates that the NWP boundary layer depth is to be
                                                      ! used if its available.
  Real(Std),          Intent(In) :: SigUUM            ! Velocity variance of unresolved mesoscale motions.
  Real(Std),          Intent(In) :: TauUUM            ! Lagrangian timescale of unresolved mesoscale motions.
  Real(Std),          Intent(In) :: SigU2HPlus        !] Free tropospheric velocity variances (horizontal
  Real(Std),          Intent(In) :: SigW2HPlus        !] and vertical).
  Real(Std),          Intent(In) :: TauUHPlus         !] Free tropospheric Lagrangian timescales
  Real(Std),          Intent(In) :: TauWHPlus         !] (horizontal and vertical).
  Character(*),       Intent(In) :: MetFolder         ! Folder containing met data.
  Character(*),       Intent(In) :: MetFolderStem     ! Stem of folders containing met data.
  Character(MaxTokenLength), Intent(In) :: MetFolders(:)  ! Array of folders containing met data.
                                                          ! $$ changed from * to support Intel compiler on Cray
  Character(*),       Intent(In) :: EnsembleMetFolder ! Folder containing ensemble of met data.
  Character(*),       Intent(In) :: TopogFolder       ! Folder containing topography data.
  Character(*),       Intent(In) :: RestoreMetScript  ! File containing script for restoring met data. Blank
                                                      ! indicates met files will not be restored.
  Logical,            Intent(In) :: DeleteMet         ! Indicates met files will be deleted after use.
  Type(OpenMPOpts_),  Intent(In) :: OpenMPOpts        ! OpenMP options
  ! Function result:
  Type(NWPMet_) :: NWPMet ! Initialised state of an NWP met module instance.
  ! Locals:
  Integer :: L ! Trimed length of MetFolder(s) or TopogFolder.
  Integer :: i ! Loop index.

  NWPMet%C = InitCommonMet('NWP Met', MetName, FixedMet, UpdateOnDemand)

  NWPMet%MetDefn = NWPMetDefn

  NWPMet%Dt = Time2ShortTime(NWPMetDefn%Dt)

  NWPMet%MetDefn2s(1:NWPMet%MetDefn%nMetDefn2Names) = NWPMetDefn2s(:)

  If (HMin > HMax) Then
    Call Message('FATAL ERROR: HMin > HMax', 3)
  End If
  NWPMet%HMin       = HMin
  NWPMet%HMax       = HMax
  NWPMet%UseNWPH    = UseNWPH
  NWPMet%SigUUM     = SigUUM
  NWPMet%TauUUM     = TauUUM
  NWPMet%SigU2HPlus = SigU2HPlus
  NWPMet%SigW2HPlus = SigW2HPlus
  NWPMet%TauUHPlus  = TauUHPlus
  NWPMet%TauWHPlus  = TauWHPlus
  NWPMet%OpenMPOpts = OpenMPOpts

  If (                             &
    MetFolder         == ' ' .and. &
    MetFolderStem     == ' ' .and. &
    Size(MetFolders)  == 0   .and. &
    EnsembleMetFolder == ' '       &
  ) Then
    NWPMet%MetFolderDefnType = 0
  Else If (                        &
    MetFolder         /= ' ' .and. &
    MetFolderStem     == ' ' .and. &
    Size(MetFolders)  == 0   .and. &
    EnsembleMetFolder == ' '       &
  ) Then
    NWPMet%MetFolderDefnType = 0
  Else If (                        &
    MetFolder         == ' ' .and. &
    MetFolderStem     /= ' ' .and. &
    Size(MetFolders)  == 0   .and. &
    EnsembleMetFolder == ' '       &
  ) Then
    NWPMet%MetFolderDefnType = 1
  Else If (                        &
    MetFolder         == ' ' .and. &
    MetFolderStem     == ' ' .and. &
    Size(MetFolders)  /= 0   .and. &
    EnsembleMetFolder == ' '       &
  ) Then
    NWPMet%MetFolderDefnType = 2
  Else If (                        &
    MetFolder         == ' ' .and. &
    MetFolderStem     == ' ' .and. &
    Size(MetFolders)  == 0   .and. &
    EnsembleMetFolder /= ' '       &
  ) Then
    NWPMet%MetFolderDefnType = 3
  Else
    Call Message(                                                              &
           'FATAL ERROR in reading the NWP met module instance "'           // &
           Trim(MetName)                                                    // &
           '" from the input file(s): at most one of '                      // &
           'Met Folder, Met Folder Stem and Met Folders must be specified',    &
           3                                                                   &
         )
  End If

  Select Case (NWPMet%MetFolderDefnType)

    ! Single met folder.
    Case (0)

      L = Len_Trim(MetFolder)
      If (L == 0) Then
        NWPMet%MetFolder = '.\'
      Else If (                                                  &
        (MetFolder(L:L) == '\' .or. MetFolder(L:L) == '/') .and. &
        L <= MaxFileNameLength                                   &
      ) Then
        NWPMet%MetFolder = MetFolder
      Else If (                                                 &
        MetFolder(L:L) /= '\' .and. MetFolder(L:L) /= '/' .and. &
        L <= MaxFileNameLength - 1                              &
      ) Then
        NWPMet%MetFolder = Trim(MetFolder) // '\'
      Else
        Call Message('FATAL ERROR: Met folder name is too long', 3)
      End If

    ! Met folder stem.
    Case (1)

      L = Len_Trim(MetFolderStem)  ! Note MetFolderStem is not blank, so L > 0.
      If (L <= MaxFileNameLength) Then
        NWPMet%MetFolderStem = MetFolderStem
      Else
        Call Message('FATAL ERROR: Met folder stem is too long', 3)
      End If

    ! Array of met folders.
    Case (2)

      Allocate (NWPMet%MetFolders(Size(MetFolders)))
      Do i = 1, Size(MetFolders)
        L = Len_Trim(MetFolders(i))
        If (L == 0) Then
          NWPMet%MetFolders(i) = '.\' ! disallow blank? $$
        Else If (                                                          &
          (MetFolders(i)(L:L) == '\' .or. MetFolders(i)(L:L) == '/') .and. &
          L <= MaxFileNameLength                                           &
        ) Then
          NWPMet%MetFolders(i) = MetFolders(i)
        Else If (                                                         &
          MetFolders(i)(L:L) /= '\' .and. MetFolders(i)(L:L) /= '/' .and. &
          L <= MaxFileNameLength - 1                                      &
        ) Then
          NWPMet%MetFolders(i) = Trim(MetFolders(i)) // '\'
        Else
          Call Message(                                                    &
                 'FATAL ERROR in reading the NWP met module instance "' // &
                 Trim(MetName)                                          // &
                 '" from the input file(s): element '                   // &
                 Int2Char(i)                                            // &
                 ' in the array of met folders is given as "'           // &
                 Trim(MetFolders(i))                                    // &
                 '" and is too long',                                      &
                 3                                                         &
               )
        End If
      End Do

    ! Single ensemble met folder.
    Case (3)

      L = Len_Trim(EnsembleMetFolder) ! Note EnsembleMetFolder is not blank, so L > 0.
      If (                                                                       &
        (EnsembleMetFolder(L:L) == '\' .or. EnsembleMetFolder(L:L) == '/') .and. &
        L <= MaxFileNameLength                                                   &
      ) Then
        NWPMet%EnsembleMetFolder = EnsembleMetFolder
      Else If (                                                                 &
        EnsembleMetFolder(L:L) /= '\' .and. EnsembleMetFolder(L:L) /= '/' .and. &
        L <= MaxFileNameLength - 1                                              &
      ) Then
        NWPMet%EnsembleMetFolder = Trim(EnsembleMetFolder) // '\'
      Else
        Call Message('FATAL ERROR: Ensemble met folder is too long', 3)
      End If

  End Select

  L = Len_Trim(TopogFolder)
  If (L == 0) Then
    NWPMet%TopogFolder = '.\'
  Else If (                                                      &
    (TopogFolder(L:L) == '\' .or. TopogFolder(L:L) == '/') .and. &
    L <= MaxFileNameLength                                       &
  ) Then
    NWPMet%TopogFolder = TopogFolder
  Else If (                                                     &
    TopogFolder(L:L) /= '\' .and. TopogFolder(L:L) /= '/' .and. &
    L <= MaxFileNameLength - 1                                  &
  ) Then
    NWPMet%TopogFolder = Trim(TopogFolder) // '\'
  Else
    Call Message('ERROR: Topography folder name too long', 3)
  End If
  NWPMet%TopogFolder = NWPMet%TopogFolder

  If (Len_Trim(RestoreMetScript) > MaxFileNameLength) Then
    Call Message(                                                  &
           'ERROR in initialising an NWP met module instance: ' // &
           'Restore Met Script is given as '                    // &
           Trim(RestoreMetScript)                               // &
           ' and is too long',                                     &
           3                                                       &
         )
  End If
  NWPMet%RestoreMetScript = RestoreMetScript
  NWPMet%DeleteMet        = DeleteMet

  NWPMet%SpaceAllocated = .false.

  NWPMet%New = 0
  NWPMet%Old = 0

  NWPMet%OldMetFiles(:) = ' '

End Function InitNWPMet

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpCoordsEtc_NWPMet(NWPMet, Coords, Grids)
! Sets up Coords and Grids by adding any extra coords and grids which NWPMet wants to
! define.

  Implicit None
  ! Argument list:
  Type(NWPMet_), Intent(In)    :: NWPMet ! State of an NWP met module instance.
  Type(Coords_), Intent(InOut) :: Coords ! Collection of coord systems.
  Type(Grids_),  Intent(InOut) :: Grids  ! Collection of grids.
  ! Locals:
  Type(ZCoord_) :: ZCoord ! Coord system to be added.

  ! Construct coord system to be added to Coords.
  ZCoord = ZCoord_m_agl()

  ! Add coord system to Coords.
  Call AddZCoord(ZCoord, Coords)

  ! Construct coord system to be added to Coords.
  ZCoord = ZCoord_Pa()

  ! Add coord system to Coords.
  Call AddZCoord(ZCoord, Coords)

End Subroutine SetUpCoordsEtc_NWPMet

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpNWPMet_CoordsEtc(EtaDefns, Coords, Grids, MetEnsembleSize, NWPMet)
! Sets up NWPMet using information from EtaDefns, Coords and Grids. ! $$ not quite true now MetEnsembleSize
! added.

  Implicit None
  ! Argument list:
  Type(EtaDefns_), Intent(In)           :: EtaDefns        ! Collection of eta definitions.
  Type(Coords_),   Intent(In)           :: Coords          ! Collection of coord systems.
  Type(Grids_),    Intent(In),   Target :: Grids           ! Collection of grids.
  Integer,         Intent(In)           :: MetEnsembleSize ! Size of the met ensemble (i.e. number of met
                                                           ! realisations).
  Type(NWPMet_),   Intent(InOut)        :: NWPMet          ! State of an NWP met module instance.
  ! Locals:
  Character(MaxCharLength)         :: ZCoordName1             !} Names of the vertical coord
  Character(MaxCharLength)         :: ZCoordName2             !} systems to add to NWPMet%C.
  Character(MaxCharLength)         :: ZCoordName3             !}
  Type(HGrid_),            Pointer :: HGrid                   !] Abbreviations for grids.
  Type(HGrid_),            Pointer :: HGridU                  !]
  Type(HGrid_),            Pointer :: HGridV                  !]
  Type(ZGrid_),            Pointer :: ZGrid                   !]
  Type(ZGrid_),            Pointer :: ZGridUV                 !]
  Type(ZGrid_),            Pointer :: ZGridW                  !]
  Type(ZGrid_),            Pointer :: ZGridP                  !]
  Logical                          :: FirstCloud              ! Flag to denote first cloud field
  Logical                          :: Error                   ! Error flags - indicates fatal error
                                                              ! which however will not be made
                                                              ! fatal just yet in order to detect
                                                              ! other possible errors.
  Real(Std)                        :: Scale                   ! Tolerance for computing grid
                                                              ! offset / grid spacing.
  Integer                          :: Offset                  ! 1, 0 and -1 indicate a grid offset
                                                              ! of 1/2, 0 and -1/2 times the grid
                                                              ! spacing.
  Logical                          :: OffsetError             ! Indicates the grid offset is not
                                                              ! acceptable.
  Integer                          :: i                       !} Loop indices.
  Integer                          :: j                       !}
  Integer                          :: k                       !}

  ! Initialise Error.
  Error = .false.

  ! Find indices of coords and grids.
  NWPMet%iZCoordZ  = FindZCoordIndex('m agl',                     Coords)
  NWPMet%iZCoordP  = FindZCoordIndex('Pa',                        Coords)
  NWPMet%iZCoordW  = FindZCoordIndex(NWPMet%MetDefn%ZCoordNameW,  Coords)
  NWPMet%iZCoordCl = FindZCoordIndex(NWPMet%MetDefn%ZCoordNameCl, Coords)
  NWPMet%iHGrid    = FindHGridIndex(NWPMet%MetDefn%HGridName,   Grids)
  NWPMet%iHGridU   = FindHGridIndex(NWPMet%MetDefn%HGridNameU,  Grids)
  NWPMet%iHGridV   = FindHGridIndex(NWPMet%MetDefn%HGridNameV,  Grids)
  NWPMet%iZGrid    = FindZGridIndex(NWPMet%MetDefn%ZGridName,   Grids)
  NWPMet%iZGridUV  = FindZGridIndex(NWPMet%MetDefn%ZGridNameUV, Grids)
  NWPMet%iZGridW   = FindZGridIndex(NWPMet%MetDefn%ZGridNameW,  Grids)
  NWPMet%iZGridP   = FindZGridIndex(NWPMet%MetDefn%ZGridNameP,  Grids)
  NWPMet%iHCoord   = FindHCoordIndex(Grids%HGrids(NWPMet%iHGrid)%HCoordName, Coords)
  NWPMet%iZCoord   = FindZCoordIndex(Grids%ZGrids(NWPMet%iZGrid)%ZCoordName, Coords)

  ! Set up the names of vertical coord systems to add to NWPMet%C.
  ZCoordName1 = Grids%ZGrids(NWPMet%iZGrid)%ZCoordName
  ZCoordName2 = 'm agl'
  ZCoordName3 = 'Pa'

  ! Add coord systems to NWPMet%C.
  Call AddCoordsToCommonMet(                              &
         1, (/ Grids%HGrids(NWPMet%iHGrid)%HCoordName /), &
         3, (/ ZCoordName1, ZCoordName2, ZCoordName3 /),  &
         NWPMet%C                                         &
       )

  ! Check number of coord systems in NWPMet.
  If (NWPMet%C%nHCoords /= 1) Then
    Call Message(                                               &
           'UNEXPECTED ERROR in SetUpNWPMet_CoordsEtc: '     // &
           'error in numbering of horizontal coord systems ' // &
           'used by the NWP met module instance "'           // &
           Trim(NWPMet%C%MetName)                            // &
           '"',                                                 &
           4                                                    &
         )
  End If
  If (NWPMet%C%nZCoords /= 3) Then
    Call Message(                                             &
           'UNEXPECTED ERROR in SetUpNWPMet_CoordsEtc: '   // &
           'error in numbering of vertical coord systems ' // &
           'used by the NWP met module instance "'         // &
           Trim(NWPMet%C%MetName)                          // &
           '"',                                               &
           4                                                  &
         )
  End If

  ! Set up abbreviations for grids.
  HGrid   => Grids%HGrids(NWPMet%iHGrid  )
  HGridU  => Grids%HGrids(NWPMet%iHGridU )
  HGridV  => Grids%HGrids(NWPMet%iHGridV )
  ZGrid   => Grids%ZGrids(NWPMet%iZGrid  )
  ZGridUV => Grids%ZGrids(NWPMet%iZGridUV)
  ZGridW  => Grids%ZGrids(NWPMet%iZGridW )
  ZGridP  => Grids%ZGrids(NWPMet%iZGridP )

  ! Check horizontal are grids regular with non-zero spacing and more than one point
  ! in each direction.
  If (                                           &
    HGrid%Unstructured .or.                      &
    HGrid%Variable     .or.                      &
    HGrid%dX == 0.0    .or. HGrid%dY == 0.0 .or. &
    HGrid%nX <= 1      .or. HGrid%nY <= 1        &
  ) Then
    Call Message(                                                            &
           'ERROR: The horizontal grid named "'                           // &
           Trim(HGrid%Name)                                               // &
           '" which is specified as "H-Grid" in the '                     // &
           '"NWP Met Definition" named "'                                 // &
           Trim(NWPMet%MetDefn%Name)                                      // &
           '" has variable or zero spacing or less than 2 points in one ' // &
           '(at least) direction.',                                          &
           2                                                                 &
         )
    Error = .true.
  End If
  If (                                             &
    HGridU%Unstructured .or.                       &
    HGridU%Variable     .or.                       &
    HGridU%dX == 0.0    .or. HGridU%dY == 0.0 .or. &
    HGridU%nX <= 1      .or. HGridU%nY <= 1        &
  ) Then
    Call Message(                                                            &
           'ERROR: The horizontal grid named "'                           // &
           Trim(HGridU%Name)                                              // &
           '" which is specified as "H-Grid - U" in the '                 // &
           '"NWP Met Definition" named "'                                 // &
           Trim(NWPMet%MetDefn%Name)                                      // &
           '" has variable or zero spacing or less than 2 points in one ' // &
           '(at least) direction.',                                          &
           2                                                                 &
         )
    Error = .true.
  End If
  If (                                             &
    HGridV%Unstructured .or.                       &
    HGridV%Variable     .or.                       &
    HGridV%dX == 0.0    .or. HGridV%dY == 0.0 .or. &
    HGridV%nX <= 1      .or. HGridV%nY <= 1        &
  ) Then
    Call Message(                                                            &
           'ERROR: The horizontal grid named "'                           // &
           Trim(HGridV%Name)                                              // &
           '" which is specified as "H-Grid - V" in the '                 // &
           '"NWP Met Definition" named "'                                 // &
           Trim(NWPMet%MetDefn%Name)                                      // &
           '" has variable or zero spacing or less than 2 points in one ' // &
           '(at least) direction.',                                          &
           2                                                                 &
         )
    Error = .true.
  End If

  ! Check horizontal grids for consistency of coord system used.
  If (.not.(HGridU%HCoordName .CIEq. HGrid%HCoordName)) Then
    Call Message(                                                          &
           'ERROR: The horizontal grids named "'                        // &
           Trim(HGrid%Name)                                             // &
           '" and "'                                                    // &
           Trim(HGridU%Name)                                            // &
           '" which are specified as "H-Grid" and "H-Grid - U" in the ' // &
           '"NWP Met Definition" named "'                               // &
           Trim(NWPMet%MetDefn%Name)                                    // &
           '" use different horizontal coordinate systems.',               &
           2                                                               &
         )
    Error = .true.
  End If
  If (.not.(HGridV%HCoordName .CIEq. HGrid%HCoordName)) Then
    Call Message(                                                          &
           'ERROR: The horizontal grids named "'                        // &
           Trim(HGrid%Name)                                             // &
           '" and "'                                                    // &
           Trim(HGridV%Name)                                            // &
           '" which are specified as "H-Grid" and "H-Grid - V" in the ' // &
           '"NWP Met Definition" named "'                               // &
           Trim(NWPMet%MetDefn%Name)                                    // &
           '" use different horizontal coordinate systems.',               &
           2                                                               &
         )
    Error = .true.
  End If

  ! Check horizontal grids for consistency of spacing, wrapping and number of points.
  If (HGridU%dX /= HGrid%dX .or. HGridU%dY /= HGrid%dY) Then
    Call Message(                                                          &
           'ERROR: The horizontal grids named "'                        // &
           Trim(HGrid%Name)                                             // &
           '" and "'                                                    // &
           Trim(HGridU%Name)                                            // &
           '" which are specified as "H-Grid" and "H-Grid - U" in the ' // &
           '"NWP Met Definition" named "'                               // &
           Trim(NWPMet%MetDefn%Name)                                    // &
           '" have different spacing.',                                    &
           2                                                               &
         )
    Error = .true.
  End If
  If (HGridU%Wrap .neqv. HGrid%Wrap) Then
    Call Message(                                                          &
           'ERROR: The horizontal grids named "'                        // &
           Trim(HGrid%Name)                                             // &
           '" and "'                                                    // &
           Trim(HGridU%Name)                                            // &
           '" which are specified as "H-Grid" and "H-Grid - U" in the ' // &
           '"NWP Met Definition" named "'                               // &
           Trim(NWPMet%MetDefn%Name)                                    // &
           '" have different values of "Wrap".',                           &
           2                                                               &
         )
    Error = .true.
  End If
  If (HGridU%Wrap .and. HGrid%Wrap .and. HGridU%nX /= HGrid%nX) Then
    Call Message(                                                             &
           'ERROR: The horizontal grids named "'                           // &
           Trim(HGrid%Name)                                                // &
           '" and "'                                                       // &
           Trim(HGridU%Name)                                               // &
           '" which are specified as "H-Grid" and "H-Grid - U" in the '    // &
           '"NWP Met Definition" named "'                                  // &
           Trim(NWPMet%MetDefn%Name)                                       // &
           '" are both "wrapped" in the X-direction but have a different ' // &
           'number of points.',                                               &
           2                                                                  &
         )
    Error = .true.
  End If
  If (HGridV%dX /= HGrid%dX .or. HGridV%dY /= HGrid%dY) Then
    Call Message(                                                          &
           'ERROR: The horizontal grids named "'                        // &
           Trim(HGrid%Name)                                             // &
           '" and "'                                                    // &
           Trim(HGridV%Name)                                            // &
           '" which are specified as "H-Grid" and "H-Grid - V" in the ' // &
           '"NWP Met Definition" named "'                               // &
           Trim(NWPMet%MetDefn%Name)                                    // &
           '" have different spacing.',                                    &
           2                                                               &
         )
    Error = .true.
  End If
  If (HGridV%Wrap .neqv. HGrid%Wrap) Then
    Call Message(                                                          &
           'ERROR: The horizontal grids named "'                        // &
           Trim(HGrid%Name)                                             // &
           '" and "'                                                    // &
           Trim(HGridV%Name)                                            // &
           '" which are specified as "H-Grid" and "H-Grid - V" in the ' // &
           '"NWP Met Definition" named "'                               // &
           Trim(NWPMet%MetDefn%Name)                                    // &
           '" have different values of "Wrap".',                           &
           2                                                               &
         )
    Error = .true.
  End If
  If (HGridV%Wrap .and. HGrid%Wrap .and. HGridV%nX /= HGrid%nX) Then
    Call Message(                                                             &
           'ERROR: The horizontal grids named "'                           // &
           Trim(HGrid%Name)                                                // &
           '" and "'                                                       // &
           Trim(HGridV%Name)                                               // &
           '" which are specified as "H-Grid" and "H-Grid - V" in the '    // &
           '"NWP Met Definition" named "'                                  // &
           Trim(NWPMet%MetDefn%Name)                                       // &
           '" are both "wrapped" in the X-direction but have a different ' // &
           'number of points.',                                               &
           2                                                                  &
         )
    Error = .true.
  End If

  ! Check offset of U grid in X direction (note that the U and V horizontal grids must
  ! be offset from the main horizontal grid by grid spacing times 1/2, 0 or -1/2).
  OffsetError = .false.
  Scale = Max(Abs(HGridU%X0/HGrid%dX), Abs(HGrid%X0/HGrid%dX), 0.5) * &
          10.0 * Epsilon(1.0_Std)
  If (Abs((HGridU%X0 - HGrid%X0)/HGrid%dX) < Scale) Then
    Offset = 0
    If (HGridU%nX /= HGrid%nX) OffsetError = .true.
  Else If (Abs((HGridU%X0 - HGrid%X0)/HGrid%dX - 0.5) < Scale) Then
    Offset = 1
    If (HGridU%nX /= HGrid%nX .and. HGridU%nX /= HGrid%nX - 1) OffsetError = .true.
  Else If (Abs((HGridU%X0 - HGrid%X0)/HGrid%dX + 0.5) < Scale) Then
    Offset = - 1
    If (HGridU%nX /= HGrid%nX .and. HGridU%nX /= HGrid%nX + 1) OffsetError = .true.
  Else
    OffsetError = .true.
  End If
  If (OffsetError) Then
    Call Message(                                                               &
           'ERROR: The horizontal grid named "'                              // &
           Trim(HGridU%Name)                                                 // &
           '" which is specified as "H-Grid - U" in the '                    // &
           '"NWP Met Definition" named "'                                    // &
           Trim(NWPMet%MetDefn%Name)                                         // &
           '" is neither aligned with the main grid nor offset by +/- grid ' // &
           'spacing / 2 at one (at least) of the edges in the X-direction.',    &
           2                                                                    &
         )
    Error = .true.
  End If

  ! Check offset of U grid in Y direction (note that the U and V horizontal grids must
  ! be offset from the main horizontal grid by grid spacing times 1/2, 0 or -1/2).
  OffsetError = .false.
  Scale = Max(Abs(HGridU%Y0/HGrid%dY), Abs(HGrid%Y0/HGrid%dY), 0.5) * &
          10.0 * Epsilon(1.0_Std)
  If (Abs((HGridU%Y0 - HGrid%Y0)/HGrid%dY) < Scale) Then
    Offset = 0
    If (HGridU%nY /= HGrid%nY) OffsetError = .true.
  Else If (Abs((HGridU%Y0 - HGrid%Y0)/HGrid%dY - 0.5) < Scale) Then
    Offset = 1
    If (HGridU%nY /= HGrid%nY .and. HGridU%nY /= HGrid%nY - 1) OffsetError = .true.
  Else If (Abs((HGridU%Y0 - HGrid%Y0)/HGrid%dY + 0.5) < Scale) Then
    Offset = - 1
    If (HGridU%nY /= HGrid%nY .and. HGridU%nY /= HGrid%nY + 1) OffsetError = .true.
  Else
    OffsetError = .true.
  End If
  If (OffsetError) Then
    Call Message(                                                               &
           'ERROR: The horizontal grid named "'                              // &
           Trim(HGridU%Name)                                                 // &
           '" which is specified as "H-Grid - U" in the '                    // &
           '"NWP Met Definition" named "'                                    // &
           Trim(NWPMet%MetDefn%Name)                                         // &
           '" is neither aligned with the main grid nor offset by +/- grid ' // &
           'spacing / 2 at one (at least) of the edges in the Y-direction.',    &
           2                                                                    &
         )
    Error = .true.
  End If

  ! Check offset of V grid in X direction (note that the U and V horizontal grids must
  ! be offset from the main horizontal grid by grid spacing times 1/2, 0 or -1/2).
  OffsetError = .false.
  Scale = Max(Abs(HGridV%X0/HGrid%dX), Abs(HGrid%X0/HGrid%dX), 0.5) * &
          10.0 * Epsilon(1.0_Std)
  If (Abs((HGridV%X0 - HGrid%X0)/HGrid%dX) < Scale) Then
    Offset = 0
    If (HGridV%nX /= HGrid%nX) OffsetError = .true.
  Else If (Abs((HGridV%X0 - HGrid%X0)/HGrid%dX - 0.5) < Scale) Then
    Offset = 1
    If (HGridV%nX /= HGrid%nX .and. HGridV%nX /= HGrid%nX - 1) OffsetError = .true.
  Else If (Abs((HGridV%X0 - HGrid%X0)/HGrid%dX + 0.5) < Scale) Then
    Offset = - 1
    If (HGridV%nX /= HGrid%nX .and. HGridV%nX /= HGrid%nX + 1) OffsetError = .true.
  Else
    OffsetError = .true.
  End If
  If (OffsetError) Then
    Call Message(                                                               &
           'ERROR: The horizontal grid named "'                              // &
           Trim(HGridV%Name)                                                 // &
           '" which is specified as "H-Grid - V" in the '                    // &
           '"NWP Met Definition" named "'                                    // &
           Trim(NWPMet%MetDefn%Name)                                         // &
           '" is neither aligned with the main grid nor offset by +/- grid ' // &
           'spacing / 2 at one (at least) of the edges in the X-direction.',    &
           2                                                                    &
         )
    Error = .true.
  End If

  ! Check offset of V grid in Y direction (note that the U and V horizontal grids must
  ! be offset from the main horizontal grid by grid spacing times 1/2, 0 or -1/2).
  OffsetError = .false.
  Scale = Max(Abs(HGridV%Y0/HGrid%dY), Abs(HGrid%Y0/HGrid%dY), 0.5) * &
          10.0 * Epsilon(1.0_Std)
  If (Abs((HGridV%Y0 - HGrid%Y0)/HGrid%dY) < Scale) Then
    Offset = 0
    If (HGridV%nY /= HGrid%nY) OffsetError = .true.
  Else If (Abs((HGridV%Y0 - HGrid%Y0)/HGrid%dY - 0.5) < Scale) Then
    Offset = 1
    If (HGridV%nY /= HGrid%nY .and. HGridV%nY /= HGrid%nY - 1) OffsetError = .true.
  Else If (Abs((HGridV%Y0 - HGrid%Y0)/HGrid%dY + 0.5) < Scale) Then
    Offset = - 1
    If (HGridV%nY /= HGrid%nY .and. HGridV%nY /= HGrid%nY + 1) OffsetError = .true.
  Else
    OffsetError = .true.
  End If
  If (OffsetError) Then
    Call Message(                                                               &
           'ERROR: The horizontal grid named "'                              // &
           Trim(HGridV%Name)                                                 // &
           '" which is specified as "H-Grid - V" in the '                    // &
           '"NWP Met Definition" named "'                                    // &
           Trim(NWPMet%MetDefn%Name)                                         // &
           '" is neither aligned with the main grid nor offset by +/- grid ' // &
           'spacing / 2 at one (at least) of the edges in the Y-direction.',    &
           2                                                                    &
         )
    Error = .true.
  End If

  ! Check vertical grids have height increasing with index. Note this also ensures
  ! there are at least two grid points.
  ! $$ support decreasing, or different behaviour for different grids?
  If (ZGrid%IncreasingHeight /= 1) Then
    Call Message(                                                   &
           'ERROR: The vertical grid named "'                    // &
           Trim(ZGrid%Name)                                      // &
           '" which is specified as "Z-Grid" in the '            // &
           '"NWP Met Definition" named "'                        // &
           Trim(NWPMet%MetDefn%Name)                             // &
           '" does not have height increasing with grid level.',    &
           2                                                        &
         )
    Error = .true.
  End If
  If (ZGridUV%IncreasingHeight /= 1) Then
    Call Message(                                                   &
           'ERROR: The vertical grid named "'                    // &
           Trim(ZGridUV%Name)                                    // &
           '" which is specified as "Z-Grid - UV" in the '       // &
           '"NWP Met Definition" named "'                        // &
           Trim(NWPMet%MetDefn%Name)                             // &
           '" does not have height increasing with grid level.',    &
           2                                                        &
         )
    Error = .true.
  End If
  If (ZGridW%IncreasingHeight /= 1) Then
    Call Message(                                                   &
           'ERROR: The vertical grid named "'                    // &
           Trim(ZGridW%Name)                                     // &
           '" which is specified as "Z-Grid - W" in the '        // &
           '"NWP Met Definition" named "'                        // &
           Trim(NWPMet%MetDefn%Name)                             // &
           '" does not have height increasing with grid level.',    &
           2                                                        &
         )
    Error = .true.
  End If
  If (ZGridP%IncreasingHeight /= 1) Then
    Call Message(                                                   &
           'ERROR: The vertical grid named "'                    // &
           Trim(ZGridP%Name)                                     // &
           '" which is specified as "Z-Grid - P" in the '        // &
           '"NWP Met Definition" named "'                        // &
           Trim(NWPMet%MetDefn%Name)                             // &
           '" does not have height increasing with grid level.',    &
           2                                                        &
         )
    Error = .true.
  End If

  ! Check vertical grids for consistency of coord system used.
  If (.not.(ZGridUV%ZCoordName .CIEq. ZGrid%ZCoordName)) Then
    Call Message(                                                           &
           'ERROR: The vertical grids named "'                           // &
           Trim(ZGrid%Name)                                              // &
           '" and "'                                                     // &
           Trim(ZGridUV%Name)                                            // &
           '" which are specified as "Z-Grid" and "Z-Grid - UV" in the ' // &
           '"NWP Met Definition" named "'                                // &
           Trim(NWPMet%MetDefn%Name)                                     // &
           '" use different vertical coordinate systems.',                  &
           2                                                                &
         )
    Error = .true.
  End If
  If (.not.(ZGridW%ZCoordName .CIEq. ZGrid%ZCoordName)) Then
    Call Message(                                                          &
           'ERROR: The vertical grids named "'                          // &
           Trim(ZGrid%Name)                                             // &
           '" and "'                                                    // &
           Trim(ZGridW%Name)                                            // &
           '" which are specified as "Z-Grid" and "Z-Grid - W" in the ' // &
           '"NWP Met Definition" named "'                               // &
           Trim(NWPMet%MetDefn%Name)                                    // &
           '" use different vertical coordinate systems.',                 &
           2                                                               &
         )
    Error = .true.
  End If
  If (.not.(ZGridP%ZCoordName .CIEq. ZGrid%ZCoordName)) Then
    Call Message(                                                          &
           'ERROR: The vertical grids named "'                          // &
           Trim(ZGrid%Name)                                             // &
           '" and "'                                                    // &
           Trim(ZGridP%Name)                                            // &
           '" which are specified as "Z-Grid" and "Z-Grid - P" in the ' // &
           '"NWP Met Definition" named "'                               // &
           Trim(NWPMet%MetDefn%Name)                                    // &
           '" use different vertical coordinate systems.',                 &
           2                                                               &
         )
    Error = .true.
  End If

  ! Check vertical grids have index array.
  If (.not. ZGrid%UseiZ) Then
    Call Message(                                          &
           'ERROR: The vertical grid named "'           // &
           Trim(ZGrid%Name)                             // &
           '" which is specified as "Z-Grid" in the '   // &
           '"NWP Met Definition" named "'               // &
           Trim(NWPMet%MetDefn%Name)                    // &
           '" does not have an index array specified.',    &
           2                                               &
         )
    Error = .true.
  End If
  If (.not. ZGridUV%UseiZ) Then
    Call Message(                                             &
           'ERROR: The vertical grid named "'              // &
           Trim(ZGrid%Name)                                // &
           '" which is specified as "Z-Grid - UV" in the ' // &
           '"NWP Met Definition" named "'                  // &
           Trim(NWPMet%MetDefn%Name)                       // &
           '" does not have an index array specified.',       &
           2                                                  &
         )
    Error = .true.
  End If
  If (.not. ZGridW%UseiZ) Then
    Call Message(                                            &
           'ERROR: The vertical grid named "'             // &
           Trim(ZGrid%Name)                               // &
           '" which is specified as "Z-Grid - W" in the ' // &
           '"NWP Met Definition" named "'                 // &
           Trim(NWPMet%MetDefn%Name)                      // &
           '" does not have an index array specified.',      &
           2                                                 &
         )
    Error = .true.
  End If
  If (.not. ZGridP%UseiZ) Then
    Call Message(                                            &
           'ERROR: The vertical grid named "'             // &
           Trim(ZGrid%Name)                               // &
           '" which is specified as "Z-Grid - P" in the ' // &
           '"NWP Met Definition" named "'                 // &
           Trim(NWPMet%MetDefn%Name)                      // &
           '" does not have an index array specified.',      &
           2                                                 &
         )
    Error = .true.
  End If

  ! Fatal error.
  If (Error) Then
    Call Message(                                        &
           'FATAL ERROR: the grids specified in the ' // &
           '"NWP Met Definition" named "'             // &
           Trim(NWPMet%MetDefn%Name)                  // &
           '" are not useable.',                         &
           3                                             &
         )
  End If

  ! Calculate interpolation coefficients.
  Call InitGHCoeffs(NWPMet%GHCoeffs  )
  Call InitGHCoeffs(NWPMet%GHCoeffsU )
  Call InitGHCoeffs(NWPMet%GHCoeffsV )
  Call InitGHCoeffs(NWPMet%GHCoeffsdX)
  Call InitGHCoeffs(NWPMet%GHCoeffsdY)
  Call GetGHCoeffs(HGrid , HGrid, ' ', NWPMet%GHCoeffs  )
  Call GetGHCoeffs(HGridU, HGrid, ' ', NWPMet%GHCoeffsU )
  Call GetGHCoeffs(HGridV, HGrid, ' ', NWPMet%GHCoeffsV )
  Call GetGHCoeffs(HGrid , HGrid, 'X', NWPMet%GHCoeffsdX)
  Call GetGHCoeffs(HGrid , HGrid, 'Y', NWPMet%GHCoeffsdY)
  Allocate(NWPMet%ZCoeffsToW  (ZGridW%nZ))
  Allocate(NWPMet%ZCoeffsUVToW(ZGridW%nZ))
  Allocate(NWPMet%ZCoeffsToP  (ZGridP%nZ))
  Do k = 1, ZGridW%nZ
    Call GetZCoeffs(ZGridW%Z(k), ZGrid,   NWPMet%ZCoeffsToW  (k))
    Call GetZCoeffs(ZGridW%Z(k), ZGridUV, NWPMet%ZCoeffsUVToW(k))
  End Do
  Do k = 1, ZGridP%nZ
    Call GetZCoeffs(ZGridP%Z(k), ZGrid, NWPMet%ZCoeffsToP(k))
  End Do

  ! Set up kMax.
  Do i = 1, nFieldNames3d
    Select Case (i)
      Case (F_T, F_Q, F_TotalOrDynCloudWater, F_TotalOrDynCloudIce)
        NWPMet%kMax(i) = ZGrid%nZ
      Case (F_U, F_V)
        NWPMet%kMax(i) = ZGridUV%nZ
      Case (F_W)
        NWPMet%kMax(i) = ZGridW%nZ
      Case (F_PAsRead)
        NWPMet%kMax(i) = ZGridP%nZ
      Case (F_CanopyHeight, F_CanopyWater, F_StomataConduct)
        NWPMet%kMax(i) = 5   ! $$ better not to hardwire this long term.
      Case Default
        Call Message(                                                 &
               'UNEXPECTED ERROR in SetUpNWPMet_CoordsEtc: '       // &
               'there is a problem with the codes for the fields',    &
               4                                                      &
             )
    End Select
  End Do

  ! Allocate FieldPresent3d, kLower, kUpper and FieldPersist3d.
  Allocate(                                            &
    NWPMet%FieldPresent3d(                             &
      Max(ZGrid%nZ, ZGridUV%nZ, ZGridW%nZ, ZGridP%nZ), &
      nFieldNames3d                                    &
    )                                                  &
  )
  Allocate(                                            &
    NWPMet%kLower(                                     &
      Max(ZGrid%nZ, ZGridUV%nZ, ZGridW%nZ, ZGridP%nZ), &
      nFieldNames3d                                    &
    )                                                  &
  )
  Allocate(                                            &
    NWPMet%kUpper(                                     &
      Max(ZGrid%nZ, ZGridUV%nZ, ZGridW%nZ, ZGridP%nZ), &
      nFieldNames3d                                    &
    )                                                  &
  )
  Allocate(                                            &
    NWPMet%FieldPersist3d(                             &
      Max(ZGrid%nZ, ZGridUV%nZ, ZGridW%nZ, ZGridP%nZ), &
      nFieldNames3d                                    &
    )                                                  &
  )

  ! Initialise FirstCloud, Error, TotalCloudFlag, FieldPersist3d and FieldPersist2d.
  FirstCloud                  = .true.
  Error                       = .false.
  NWPMet%TotalCloudFlag       = .false.
  NWPMet%FieldPersist3d(:, :) = 0
  NWPMet%FieldPersist2d(:)    = 0

  Do k = 1, NWPMet%MetDefn%nMetDefn2Names
  
    ! Determine the field index and dimension of each field in each NWPMetDefn2s.
    Do i = 1, NWPMet%MetDefn2s(k)%nFields
      NWPMet%iField(k,i) = 0
      Do j = 1, nFieldNames3d
        If (NWPMet%MetDefn2s(k)%FieldNames(i) .CIEq. FieldNames3d(j)) Then
          NWPMet%iField(k,i) = j
          NWPMet%ThreeD(k,i) = .true.
          
          ! Determine TotalCloudFlag
          Select Case (j)
            Case (                    &
              F_TotalOrDynCloudWater, &
              F_TotalOrDynCloudIce    &
            )
              If (FirstCloud) Then
                NWPMet%TotalCloudFlag = NWPMet%MetDefn2s(k)%Total(i)
                FirstCloud            = .false.
              Else
                If (NWPMet%TotalCloudFlag .neqv. NWPMet%MetDefn2s(k)%Total(i)) Then
                  Call Message(                                                                &
                         'FATAL ERROR: Inconsistent field qualifiers for dynamic or total ' // &
                         'cloud fields in the NWP Met File Structure Definitions '          // &
                         'for the NWP Met Definition: '                                     // &
                         Trim(NWPMet%MetDefn%Name),                                            &
                         3                                                                     &
                       )
                End If
              End If
          End Select
        End If
      End Do
      Do j = 1, nFieldNames2d
        If (NWPMet%MetDefn2s(k)%FieldNames(i) .CIEq. FieldNames2d(j)) Then
          NWPMet%iField(k,i) = j
          NWPMet%ThreeD(k,i) = .false.
          
          ! Determine TotalCloudFlag   $$ No dyn / total cloud fields at all?
          Select Case (j)
            Case (                &
              F_TotalOrDynHCloud, &
              F_TotalOrDynMCloud, &
              F_TotalOrDynLCloud  &
            )
              If (FirstCloud) Then
                NWPMet%TotalCloudFlag = NWPMet%MetDefn2s(k)%Total(i)
                FirstCloud            = .false.
              Else
                If (NWPMet%TotalCloudFlag .neqv. NWPMet%MetDefn2s(k)%Total(i)) Then
                  Call Message(                                                                &
                         'FATAL ERROR: Inconsistent field qualifiers for dynamic or total ' // &
                         'cloud fields in the NWP Met File Structure Definitions '          // &
                         'for the NWP Met Definition: '                                     // &
                         Trim(NWPMet%MetDefn%Name),                                            &
                         3                                                                     &
                       )
                End If
              End If
          End Select
        End If
      End Do
      If (NWPMet%iField(k,i) == 0) Then
        Call Message(                                          &
               'ERROR: the met file structure definition "' // &
               Trim(NWPMet%MetDefn2s(k)%Name)                // &
               '" contains an invalid field name "'         // &
               Trim(NWPMet%MetDefn2s(k)%FieldNames(i))       // &
               '".',                                           &
               2                                               &
             )
        Error = .true.
      End If
    End Do

    ! Set up k1 and k2.
    Do i = 1, NWPMet%MetDefn2s(k)%nFields
      If (NWPMet%ThreeD(k,i) .and. NWPMet%iField(k,i) /= 0) Then
        NWPMet%k1(k,i) = NWPMet%MetDefn2s(k)%LowestLevels(i)
        If (NWPMet%MetDefn2s(k)%Top(i)) Then
          NWPMet%k2(k,i) = NWPMet%kMax(NWPMet%iField(k,i))
        Else
          NWPMet%k2(k,i) = NWPMet%MetDefn2s(k)%HighestLevels(i)
        End If
        If (NWPMet%k2(k,i) > NWPMet%kMax(NWPMet%iField(k,i))) Then
          Call Message(                                                    &
                 'ERROR: The highest level of the '                     // &
                 Trim(Int2Char(i)) // Int2Ordinal(i)                    // &
                 ' field ('                                             // &
                 Trim(NWPMet%MetDefn2s(k)%FieldNames(i))                 // &
                 ') in the met file structure definition "'             // &
                 Trim(NWPMet%MetDefn2s(k)%Name)                          // &
                 '" is greater than the number of levels in the grid.',    &
                 2                                                         &
               )
          Error = .true.
        End If
        If (NWPMet%k2(k,i) < NWPMet%k1(k,i)) Then
          Call Message(                                        &
                 'ERROR: The highest level of the '         // &
                 Trim(Int2Char(i)) // Int2Ordinal(i)        // &
                 ' field ('                                 // &
                 Trim(NWPMet%MetDefn2s(k)%FieldNames(i))     // &
                 ') in the met file structure definition "' // &
                 Trim(NWPMet%MetDefn2s(k)%Name)              // &
                 '" is less than the lowest level',            &
                 2                                             &
               )
          Error = .true.
        End If
        If (NWPMet%k1(k,i) < 1) Then
          Call Message(                                        &
                 'ERROR: The lowest level of the '          // &
                 Trim(Int2Char(i)) // Int2Ordinal(i)        // &
                 ' field ('                                 // &
                 Trim(NWPMet%MetDefn2s(k)%FieldNames(i))     // &
                 ') in the met file structure definition "' // &
                 Trim(NWPMet%MetDefn2s(k)%Name)              // &
                 '" is less than 1',                           &
                 2                                             &
               )
          Error = .true.
        End If
      Else
        NWPMet%k1(k,i) = 1
        NWPMet%k2(k,i) = 1
      End If
    End Do

  ! Fatal error.
    If (Error) Then
      Call Message(                                                &
             'FATAL ERROR: the met file structure definition "' // &
             Trim(NWPMet%MetDefn2s(k)%Name)                     // &
             '" is not useable.',                                  &
             3                                                     &
           )
    End If
  End Do

  ! Check consistent with MetEnsembleSize.
  If (NWPMet%MetFolderDefnType == 0) Then
    If (MetEnsembleSize > 1) Then
      Call Message(                                                                      &
             'FATAL ERROR in reading the NWP met module instance "'                   // &
             Trim(NWPMet%C%MetName)                                                   // &
             '" from the input file(s): using a single met folder for multiple met '  // &
             'realisations is currently an unsupported option',                          &
             3                                                                           &
           )
    End If
  Else If (NWPMet%MetFolderDefnType == 2) Then
    If (MetEnsembleSize > Size(NWPMet%MetFolders)) Then
      Call Message(                                                                      &
             'FATAL ERROR in reading the NWP met module instance "'                   // &
             Trim(NWPMet%C%MetName)                                                   // &
             '" from the input file(s): using a single met folder for multiple met '  // &
             'realisations is currently an unsupported option',                          &
             3                                                                           &
           )
    End If
  End If

End Subroutine SetUpNWPMet_CoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForUpdateNWPMet( &
             Coords, Grids,        &
             iCase, iMetCase,      &
             Time,                 &
             TValid, UpdateNow,    &
             NWPMet,               &
             Units                 &
           )
! Prepares for updating an instance of the NWP met module.

! This routine must set TValid and UpdateNow but must not alter the validity of the NWP met module instance.
! $$ this comment should be with generic definition.

  Implicit None
  ! Argument List:
  Type(Coords_), Intent(In)           :: Coords
  Type(Grids_),  Intent(In),   Target :: Grids
  Integer,       Intent(In)           :: iCase
  Integer,       Intent(In)           :: iMetCase
  Type(Time_),   Intent(In)           :: Time
  Type(Time_),   Intent(Out)          :: TValid
  Logical,       Intent(Out)          :: UpdateNow
  Type(NWPMet_), Intent(InOut)        :: NWPMet
  Type(Units_),  Intent(InOut)        :: Units
  ! Coords    :: Collection of coord systems.
  ! Grids     :: Collection of grids.
  ! iCase     :: Number of case.
  ! iMetCase  :: Number of the met realisation in the met ensemble.
  ! Time      :: Time for which the NWP met module instance might be updated.
  ! TValid    :: Earliest time that the validity (overall or for any single attribute) of the met module
  !              instance might change, assuming the met module instance is updated now. The value is that
  !              determined at the end of this routine (the actual time may be later).
  ! UpdateNow :: Indicates the met module instance must be updated now (even if update-on-demand is
  !              specified). If set, TValid need not be set to any particular time.
  ! NWPMet    :: State of an NWP met module instance.
  ! Units     :: Collection of information on input/output unit numbers.
  ! Locals:
  Type(Time_) :: MetTime1 !} Times of the met data required.
  Type(Time_) :: MetTime2 !}

  ! Calculate times of met data.
  MetTime1 = Round(Time, NWPMet%MetDefn%T0, NWPMet%MetDefn%Dt, Up = .false.)
             ! Could name files using a time zone $$
             ! Currently files named in UTC.
  MetTime2  = MetTime1 + NWPMet%MetDefn%Dt

  ! Set validity flags.
  If (NWPMet%C%FixedMet) Then
    TValid = InfFutureTime()
  Else
    TValid = MetTime2
  End If
  UpdateNow = .false.

End Subroutine PrepareForUpdateNWPMet

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateNWPMet(      &
             Coords, Grids,   &
             iCase, iMetCase, &
             Time,            &
             NWPMet,          &
             Units            &
           )
! Updates an instance of the NWP met module.


  Implicit None
  ! Argument List:
  Type(Coords_), Intent(In)           :: Coords   ! Collection of coord systems.
  Type(Grids_),  Intent(In),   Target :: Grids    ! Collection of grids.
  Integer,       Intent(In)           :: iCase    ! Number of case.
  Integer,       Intent(In)           :: iMetCase ! Number of the met realisation in the met ensemble.
  Type(Time_),   Intent(In)           :: Time     ! Time for which the NWP met module instance is to be
                                                  ! updated.
  Type(NWPMet_), Intent(InOut)        :: NWPMet   ! State of an NWP met module instance.
  Type(Units_),  Intent(InOut)        :: Units    ! Collection of information on input/output unit numbers.
  ! Locals:
  Type(Time_)              :: MetTime1        !} Times and short times of the met data
  Type(Time_)              :: MetTime2        !} required.
  Type(Time_)              :: MetTime3        !}
  Type(ShortTime_)         :: SMetTime1       !}
  Type(ShortTime_)         :: SMetTime2       !}
  Type(ShortTime_)         :: SMetTime3       !}
  Logical                  :: ReadOneTimeOnly ! Indicates that only one set of met data is to be read
  Logical                  :: Error           ! Indicates that an error occurred in ReadNWPMet.
  Logical                  :: Skip            ! Indicates that a file should be skipped
  Logical                  :: AllowSkip       ! Indicated whether we allow skipping of a particular file
  Type(HGrid_),    Pointer :: HGrid           !} Abbreviations for grids.
  Type(HGrid_),    Pointer :: HGridU          !}
  Type(HGrid_),    Pointer :: HGridV          !}
  Type(ZGrid_),    Pointer :: ZGrid           !}
  Type(ZGrid_),    Pointer :: ZGridUV         !}
  Type(ZGrid_),    Pointer :: ZGridW          !}
  Type(ZGrid_),    Pointer :: ZGridP          !}
  Logical                  :: DoIfBlock       ! Indicates an if-block is to be executed
                                              ! where the logical test can't be evaluated
                                              ! in-line because of risk of error in testing
                                              ! equality of times.
  Integer                  :: i               ! Loop index.
  Integer                  :: nSkipped        ! Number of skipped MetData files
  Integer                  :: SaveOldIndex    ! Temporary used for buffer reordering
  Integer                  :: NBuffers        ! Number of buffers that we are using

  Logical                  :: first

  Type(OpenMPOpts_)        :: OpenMPOpts

  Integer                  :: nThread

  ! Read OpenMP options for this specific NWPMet module
  OpenMPOpts = NWPMet%OpenMPOpts

  ! Set up abbreviations for grids.
  HGrid   => Grids%HGrids(NWPMet%iHGrid  )
  HGridU  => Grids%HGrids(NWPMet%iHGridU )
  HGridV  => Grids%HGrids(NWPMet%iHGridV )
  ZGrid   => Grids%ZGrids(NWPMet%iZGrid  )
  ZGridUV => Grids%ZGrids(NWPMet%iZGridUV)
  ZGridW  => Grids%ZGrids(NWPMet%iZGridW )
  ZGridP  => Grids%ZGrids(NWPMet%iZGridP )

  If (OpenMPOpts%ParallelMetRead) Then ! add an extra buffer for prefetching input data
    NBuffers=3
  Else ! use the original scheme
    NBuffers=2
  End If

  ! Allocate arrays.
  If (.not.NWPMet%SpaceAllocated) Then
    Allocate(NWPMet%U                   (HGridU%nX, HGridU%nY, ZGridUV%nZ, NBuffers))
    Allocate(NWPMet%V                   (HGridV%nX, HGridV%nY, ZGridUV%nZ, NBuffers))
    Allocate(NWPMet%W                   (HGrid%nX,  HGrid%nY,  ZGridW%nZ,  NBuffers))
    Allocate(NWPMet%T                   (HGrid%nX,  HGrid%nY,  ZGrid%nZ,   NBuffers))
    Allocate(NWPMet%Theta               (HGrid%nX,  HGrid%nY,  ZGrid%nZ,   NBuffers))
    Allocate(NWPMet%Q                   (HGrid%nX,  HGrid%nY,  ZGrid%nZ,   NBuffers))
    Allocate(NWPMet%PAsRead             (HGrid%nX,  HGrid%nY,  ZGridP%nZ,  NBuffers))
    Allocate(NWPMet%P                   (HGrid%nX,  HGrid%nY,  ZGrid%nZ,   NBuffers))
    Allocate(NWPMet%FLPa                (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%Rho                 (HGrid%nX,  HGrid%nY,  ZGrid%nZ,   NBuffers))
    Allocate(NWPMet%Z                   (HGrid%nX,  HGrid%nY,  ZGrid%nZ,   NBuffers))
    Allocate(NWPMet%Topog               (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%Z0                  (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%UStress             (HGridU%nX, HGridU%nY,             NBuffers))
    Allocate(NWPMet%VStress             (HGridV%nX, HGridV%nY,             NBuffers))
    Allocate(NWPMet%UStar               (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%HeatFlux            (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%H                   (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%PSeaLevel           (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%ConCloud            (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%ConCloudBase        (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%ConCloudTop         (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%ConCloudBasePa      (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%ConCloudTopPa       (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%Cloud               (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%Cloud3d             (HGrid%nX,  HGrid%nY,  ZGrid%nZ,   NBuffers))
    Allocate(NWPMet%TotalOrDynCloudWater(HGrid%nX,  HGrid%nY,  ZGrid%nZ,   NBuffers))
    Allocate(NWPMet%TotalOrDynCloudIce  (HGrid%nX,  HGrid%nY,  ZGrid%nZ,   NBuffers))
    Allocate(NWPMet%TotalOrDynCloudBase (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%TotalOrDynCloudTop  (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%TotalOrDynHCloud    (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%TotalOrDynMCloud    (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%TotalOrDynLCloud    (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%ConPpt              (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%ConRain             (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%ConSnow             (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%DynPpt              (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%DynRain             (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%DynSnow             (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%SoilMoisture        (HGrid%nX,  HGrid%nY,              NBuffers))
    Allocate(NWPMet%CanopyHeight        (HGrid%nX,  HGrid%nY,  5,          NBuffers))
    Allocate(NWPMet%CanopyWater         (HGrid%nX,  HGrid%nY,  5,          NBuffers))
    Allocate(NWPMet%StomataConduct      (HGrid%nX,  HGrid%nY,  5,          NBuffers))
    ! Allocate additional array for storing EtaDot (applied only when the Eulerian advection scheme is used).
    ! $$ Replace EulerianModelGlobal with DispOpts%EulerianModel
    If (EulerianModelGlobal) Then
      Allocate(NWPMet%EtaDot        (HGrid%nX,  HGrid%nY,  ZGridW%nZ,  NBuffers))
    End If
    ! Allocate additional array for storing the W array without the d/dt correction term (used only when input
    ! vertical velocities are pressure-based) and the corresponding EtaDot array without the d/dt correction
    ! term when additionally using the Eulerian advection scheme.
    Select Case (Coords%ZCoords(NWPMet%iZCoordW)%CoordType) ! $$ set up abbreviation?
      Case (Z_P, Z_PAsZ, Z_PAsEta)
        Allocate(NWPMet%WNoTFix(HGrid%nX, HGrid%nY, ZGridW%nZ, NBuffers))
        ! $$ Replace EulerianModelGlobal with DispOpts%EulerianModel
        If (EulerianModelGlobal) Then
          Allocate(NWPMet%EtaDotNoTFix(HGrid%nX,  HGrid%nY,  ZGridW%nZ,  NBuffers))
        End If
    End Select
    NWPMet%SpaceAllocated = .true.
    Call ReadTopog(Grids, NWPMet, Units)
  End If

  ! Calculate times of met data.
  MetTime1 = Round(Time, NWPMet%MetDefn%T0, NWPMet%MetDefn%Dt, Up = .false.)
             ! Could name files using a time zone $$
             ! Currently files named in UTC.
  MetTime2  = MetTime1 + NWPMet%MetDefn%Dt
  SMetTime1 = Time2ShortTime(MetTime1)
  SMetTime2 = Time2ShortTime(MetTime2)
  ! Reset NWPMet%Dt to the Dt in the MetDefinitions as it might have changed due to missing data
  NWPMet%Dt = Time2ShortTime(NWPMet%MetDefn%Dt)

  If (OpenMPOpts%ParallelMetRead) Then ! Determine the time for data prefetch
    MetTime3 = MetTime2 + NWPMet%MetDefn%Dt
    SMetTime3 = Time2ShortTime(MetTime3)
  End If

  ! Only read one set of met data for FixedMet runs where the time of fixed met coincides with a data time.
  ReadOneTimeOnly = (Time == MetTime1) .and. NWPMet%C%FixedMet

  !RF Do we need to set prefetch to false if ReadOneTimeOnly is true?

  ! Note could save some time multi-case simulations by avoiding reading met in cases
  ! where (i) the met is already stored (other than previous met stored in (.., New)
  ! which we do check for) and (ii) we already know the reading will give an error.
  ! However this would make the code rather complex.
  ! Could also avoid reading two times if fixed met and if fixed met time hits a data
  ! time. $$

  If (NWPMet%New == 0) Then
    DoIfBlock = .true.
  Else
    DoIfBlock = SMetTime1 /= NWPMet%NewTime
  End If

  ! set variable first to .false. by default
  first=.false.

  If (DoIfBlock) Then

    first = .true. ! there is no valid prefetched data

    ! Mark met data in Old time level as invalid (otherwise ReadNWPMet will assume Old
    ! time level data is for the previous time and may use it to fix-up missing data).
    NWPMet%Old = 0

    ! Read Time1 met data into New time level.
    NWPMet%New     = 1
    NWPMet%NewTime = SMetTime1
    ! Do not allow skipping of first MetData
    AllowSkip = .false.
    ! Read Data
    nThread = 0
    !$ nThread = omp_get_thread_num()
    Call TimerOn(WorkerNWPMetReadTimer(nThread))
    Call ReadNWPMet(MetTime1, iCase, iMetCase, Coords, Grids, AllowSkip, Error, Skip, NWPMet, Units)
    Call TimerOff(WorkerNWPMetReadTimer(nThread))
    Call TimerWriteLast(WorkerNWPMetReadTimer(nThread))
    ! Process Data
    If (Error) Go To 9
    Call TimerOn(WorkerNWPMetProcessTimer(nThread))
    Call ProcessNWPMet(Coords, Grids, NWPMet)
    Call TimerOff(WorkerNWPMetProcessTimer(nThread))
    Call TimerWriteLast(WorkerNWPMetProcessTimer(nThread))

    ! Now we can set the initial value of NWPMet%Old appropriately
    NWPMet%Old = 3

    If (OpenMPOpts%ParallelMetRead) Then
      ! also set up initial prefetch buffer
      NWPMet%Prefetch = 2
    End If

  End If

  ! now swap buffers
  If (OpenMPOpts%ParallelMetRead) Then ! extra prefetch buffer

    !RF For fine grain sync here is where we would wait
    !RF until the IOThread data is available

    ! 3 Buffers Save=Old, Old->New, New->Prefetch, Prefetch->Save

    ! Keep a reference to the Old time buffer index
    SaveOldIndex    = NWPMet%Old

    ! Move Time1 met data into Old time level.
    NWPMet%Old     = NWPMet%New
    NWPMet%OldTime = NWPMet%NewTime

    ! Move Time2 met data from prefetch into New time level
    NWPMet%New=NWPMet%Prefetch
    NWPMet%NewTime = SMetTime2

    If (.not.(first) .and. (NWPMet%PrefetchTime.ne.NWPMet%NewTime)) Then
      Stop "Error, Prefetch buffer does not have expected timestamp"
    End If

    ! Point prefetch at unused Old buffer and set new prefetch time
    NWPMet%Prefetch=SaveOldIndex

    ! If the prefetch time is beyond the end of the run then do not set up the prefetch info
    ! Round IOMaxLookahead up to time of next MetData
    If (SMetTime3<=Time2ShortTime(Round(ShortTime2Time(IOMaxLookahead), &
        NWPMet%MetDefn%T0,NWPMet%MetDefn%Dt,up = .true.))) Then

      ! mark our prefetch read request to the IO thread
      NWPMet%C%prefetch=.true.
      NWPMet%C%iCase=iCase
      NWPMet%C%iMetCase=iMetCase
      NWPMet%C%MetTime=MetTime3
      NWPMet%C%PrefetchBufferIndex=NWPMet%Prefetch
      NWPMet%C%NewBufferIndex=NWPMet%New

    End If

    !RF For fine grain sync here is where we would tell
    ! the IO thread(s) to prefetch the next data

  Else ! 2 buffers Save=Old, Old->New, New->Save

    ! Move Time1 met data into Old time level.
    NWPMet%Old     = NWPMet%New
    NWPMet%OldTime = NWPMet%NewTime

    ! Swap Time2 met data pointer
    NWPMet%New     = 3 - NWPMet%Old
    NWPMet%NewTime = SMetTime2

  End if

  ! Skip second read of met data for FixedMet runs where the time of fixed met coincides with a data time
  ! Here the met at the New time level is set equal to the met at the Old time level.
  ! $$ Note that this may not give precisely identical results to the case where both sets of met data
  !    are read, e.g. different values for heat flux and u-star would be used when the NextHeatFlux flag
  !    is set. Also the d/dt correction term for pressure-based vertical velocity is not applied (see below).
  If (ReadOneTimeOnly) Then

    NWPMet%U                   (:, :, :, NWPMet%New) = NWPMet%U                   (:, :, :, NWPMet%Old)
    NWPMet%V                   (:, :, :, NWPMet%New) = NWPMet%V                   (:, :, :, NWPMet%Old)
    NWPMet%W                   (:, :, :, NWPMet%New) = NWPMet%W                   (:, :, :, NWPMet%Old)
    NWPMet%T                   (:, :, :, NWPMet%New) = NWPMet%T                   (:, :, :, NWPMet%Old)
    NWPMet%Theta               (:, :, :, NWPMet%New) = NWPMet%Theta               (:, :, :, NWPMet%Old)
    NWPMet%Q                   (:, :, :, NWPMet%New) = NWPMet%Q                   (:, :, :, NWPMet%Old)
    NWPMet%PAsRead             (:, :, :, NWPMet%New) = NWPMet%PAsRead             (:, :, :, NWPMet%Old)
    NWPMet%P                   (:, :, :, NWPMet%New) = NWPMet%P                   (:, :, :, NWPMet%Old)
    NWPMet%Rho                 (:, :, :, NWPMet%New) = NWPMet%Rho                 (:, :, :, NWPMet%Old)
    NWPMet%Z                   (:, :, :, NWPMet%New) = NWPMet%Z                   (:, :, :, NWPMet%Old)
    NWPMet%Topog               (:, :,    NWPMet%New) = NWPMet%Topog               (:, :,    NWPMet%Old)
    NWPMet%Z0                  (:, :,    NWPMet%New) = NWPMet%Z0                  (:, :,    NWPMet%Old)
    NWPMet%UStress             (:, :,    NWPMet%New) = NWPMet%UStress             (:, :,    NWPMet%Old)
    NWPMet%VStress             (:, :,    NWPMet%New) = NWPMet%VStress             (:, :,    NWPMet%Old)
    NWPMet%UStar               (:, :,    NWPMet%New) = NWPMet%UStar               (:, :,    NWPMet%Old)
    NWPMet%HeatFlux            (:, :,    NWPMet%New) = NWPMet%HeatFlux            (:, :,    NWPMet%Old)
    NWPMet%H                   (:, :,    NWPMet%New) = NWPMet%H                   (:, :,    NWPMet%Old)
    NWPMet%PSeaLevel           (:, :,    NWPMet%New) = NWPMet%PSeaLevel           (:, :,    NWPMet%Old)
    NWPMet%ConCloud            (:, :,    NWPMet%New) = NWPMet%ConCloud            (:, :,    NWPMet%Old)
    NWPMet%ConCloudBase        (:, :,    NWPMet%New) = NWPMet%ConCloudBase        (:, :,    NWPMet%Old)
    NWPMet%ConCloudTop         (:, :,    NWPMet%New) = NWPMet%ConCloudTop         (:, :,    NWPMet%Old)
    NWPMet%Cloud               (:, :,    NWPMet%New) = NWPMet%Cloud               (:, :,    NWPMet%Old)
    NWPMet%Cloud3d             (:, :, :, NWPMet%New) = NWPMet%Cloud3d             (:, :, :, NWPMet%Old)
    NWPMet%TotalOrDynCloudWater(:, :, :, NWPMet%New) = NWPMet%TotalOrDynCloudWater(:, :, :, NWPMet%Old)
    NWPMet%TotalOrDynCloudIce  (:, :, :, NWPMet%New) = NWPMet%TotalOrDynCloudIce  (:, :, :, NWPMet%Old)
    NWPMet%TotalOrDynCloudBase (:, :,    NWPMet%New) = NWPMet%TotalOrDynCloudBase (:, :,    NWPMet%Old)
    NWPMet%TotalOrDynCloudTop  (:, :,    NWPMet%New) = NWPMet%TotalOrDynCloudTop  (:, :,    NWPMet%Old)
    NWPMet%TotalOrDynHCloud    (:, :,    NWPMet%New) = NWPMet%TotalOrDynHCloud    (:, :,    NWPMet%Old)
    NWPMet%TotalOrDynMCloud    (:, :,    NWPMet%New) = NWPMet%TotalOrDynMCloud    (:, :,    NWPMet%Old)
    NWPMet%TotalOrDynLCloud    (:, :,    NWPMet%New) = NWPMet%TotalOrDynLCloud    (:, :,    NWPMet%Old)
    NWPMet%ConPpt              (:, :,    NWPMet%New) = NWPMet%ConPpt              (:, :,    NWPMet%Old)
    NWPMet%ConRain             (:, :,    NWPMet%New) = NWPMet%ConRain             (:, :,    NWPMet%Old)
    NWPMet%ConSnow             (:, :,    NWPMet%New) = NWPMet%ConSnow             (:, :,    NWPMet%Old)
    NWPMet%DynPpt              (:, :,    NWPMet%New) = NWPMet%DynPpt              (:, :,    NWPMet%Old)
    NWPMet%DynRain             (:, :,    NWPMet%New) = NWPMet%DynRain             (:, :,    NWPMet%Old)
    NWPMet%DynSnow             (:, :,    NWPMet%New) = NWPMet%DynSnow             (:, :,    NWPMet%Old)
    NWPMet%SoilMoisture        (:, :,    NWPMet%New) = NWPMet%SoilMoisture        (:, :,    NWPMet%Old)
    NWPMet%CanopyHeight        (:, :, :, NWPMet%New) = NWPMet%CanopyHeight        (:, :, :, NWPMet%Old)
    NWPMet%CanopyWater         (:, :, :, NWPMet%New) = NWPMet%CanopyWater         (:, :, :, NWPMet%Old)
    NWPMet%StomataConduct      (:, :, :, NWPMet%New) = NWPMet%StomataConduct      (:, :, :, NWPMet%Old)
    ! Set EtaDot when the Eulerian advection scheme is being used.
    ! $$ Replace EulerianModelGlobal with DispOpts%EulerianModel
    If (EulerianModelGlobal) Then
      NWPMet%EtaDot(:, :, :, NWPMet%New) = NWPMet%EtaDot(:, :, :, NWPMet%Old)
    End If
    Select Case (Coords%ZCoords(NWPMet%iZCoordW)%CoordType) ! $$ set up abbreviation?
      Case (Z_P, Z_PAsZ, Z_PAsEta)
        NWPMet%WNoTFix(:, :, :, NWPMet%New) = NWPMet%WNoTFix(:, :, :, NWPMet%Old)
        ! Set EtaDotNoTFix when the Eulerian advection scheme is being used.
        ! $$ Replace EulerianModelGlobal with DispOpts%EulerianModel
        If (EulerianModelGlobal) Then
          NWPMet%EtaDotNoTFix(:, :, :, NWPMet%New) = NWPMet%EtaDotNoTFix(:, :, :, NWPMet%Old)
        End If
    End Select

  Else

! If we update the first MetData, always read and process it with the Worker thread
! If we don't use parallelIO ALWAYS read and process with the Worker thread
    If (First .or. .not.(OpenMPOpts%ParallelMetRead)) Then ! Read the Met Data (and process it)
      nThread = 0
      !$ nThread = omp_get_thread_num()
      Call TimerOn(WorkerNWPMetReadTimer(nThread))
      nSkipped = 0
      Do
! Do not allow skipping MetData if we use parallel IO
        AllowSkip = ((nSkipped < MaxSkip) .and. (.not. OpenMPOpts%ParallelMetRead))
        Call ReadNWPMet(MetTime2, iCase, iMetCase, Coords, Grids, AllowSkip, Error, Skip, NWPMet, Units)
        If (Skip) Then
          nSkipped = nSkipped + 1
! If this is an interpolation file, increase the time by Dt and try to read next file
          MetTime2  = MetTime2 + NWPMet%MetDefn%Dt
          SMetTime2 = Time2ShortTime(MetTime2)
          NWPMet%Dt = NWPMet%Dt + Time2ShortTime(NWPMet%MetDefn%Dt)
          SMetTime2 = Time2ShortTime(MetTime2)
          NWPMet%NewTime = SMetTime2
        Else
          Exit
        End If
      End Do
      Call TimerOff(WorkerNWPMetReadTimer(nThread))
      Call TimerWriteLast(WorkerNWPMetReadTimer(nThread))

      If (Error) Go To 9

! If we use parallel IO, we will only enter this region when updating the first MetData
      Call TimerOn(WorkerNWPMetProcessTimer(nThread))
      Call ProcessNWPMet(Coords, Grids, NWPMet)
      Call TimerOff(WorkerNWPMetProcessTimer(nThread))
      Call TimerWriteLast(WorkerNWPMetProcessTimer(nThread))

    Else
      Error = NWPMet%C%Error ! determine if prefetch read returned an error
    End If

! post structure in NWPMet instance needs iCase,iMetCase,MetTime3,Error
! coords (IN), grids (IN), mets (INOUT)(for NWPMet) and Units (INOUT) are global data
! structures so it gets these by arg passing. mets changes over time but coords and grids are fixed
! Units will need an lock/atomic around the unit increment
! May need an error per buffer
! Need to change Error test under here to test this flag
! modify ReadNWPMet to pass buffer that we are writing to explicitely as currently it is part of NWPMet.

    If (Error) Go To 9

! Only process following data with worker thread if we use parallelIO and don't want to process
! data with IOThread
    If (                                          &
      (.not. first)                         .and. &
      (OpenMPOpts%ParallelMetRead)          .and. &
      (.not. OpenMPOpts%ParallelMetProcess)       &
    ) Then
      nThread = 0
      !$ nThread = omp_get_thread_num()
      Call TimerOn(WorkerNWPMetProcessTimer(nThread))
      Call ProcessNWPMet(Coords, Grids, NWPMet)
      Call TimerOff(WorkerNWPMetProcessTimer(nThread))
      Call TimerWriteLast(WorkerNWPMetProcessTimer(nThread))
    End If

    ! Apply d/dt correction term needed to complete the conversion of vertical velocity from a pressure-based
    ! vertical velocity to the rate of change of height above ground in metres. This needs to be done
    ! separately from ProcessNWPMet because it requires met information at two times. $$ not true that it
    ! needs to be done separately - can test NWPMet%Old /= 0 to avoid doing this on first pass through
    ! ProcessNWPMet.
    Select Case (Coords%ZCoords(NWPMet%iZCoordW)%CoordType)
      Case (Z_P, Z_PAsZ, Z_PAsEta)
        Call CalcZDotCorrection(Coords, Grids, NWPMet)
    End Select

  End If


  ! Set individual attribute validity flags.
  If (ReadOneTimeOnly) Then
    Do i = 1, nAttribs
      NWPMet%C%ValidAttribs(i) = NWPMet%ValidAttribs(i, 1)
    End Do
  Else
    Do i = 1, nAttribs
      NWPMet%C%ValidAttribs(i) = NWPMet%ValidAttribs(i, NWPMet%New) .and. NWPMet%ValidAttribs(i, NWPMet%Old)
    End Do
  End If
  If (.not. Any(NWPMet%C%ValidAttribs(1:nAttribs))) Error = .true.

  If (Error) Go To 9

  ! Set validity flags.
  NWPMet%C%Valid = .true.
  If (NWPMet%C%FixedMet) Then
    NWPMet%C%TValid = InfFutureTime()
  Else
    NWPMet%C%TValid = MetTime2
  End If

  Return
  ! Errors.
9 Continue

  Do i = 1, nAttribs
    NWPMet%C%ValidAttribs(i) = .false.
  End Do

  NWPMet%New = 0

  NWPMet%C%Valid = .false.
  If (NWPMet%C%FixedMet) Then
    NWPMet%C%TValid = InfFutureTime()
  Else
    NWPMet%C%TValid = MetTime2
  End If

End Subroutine UpdateNWPMet

!-------------------------------------------------------------------------------------------------------------

Subroutine ResetNWPMet(NWPMet)
! Resets the state of an NWP met module instance for a new realisation.

  Implicit None
  ! Argument list:
  Type(NWPMet_), Intent(InOut) :: NWPMet ! State of an NWP met module instance.

  ! Reset validity flags for met data (e.g. for a new realisation in an ensemble of met data).
  NWPMet%New = 0
  NWPMet%Old = 0

End Subroutine ResetNWPMet

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadTopog(Grids, M, Units)
! Reads the topography data.

  Implicit None
  ! Argument list:
  Type(Grids_),  Intent(In),   Target :: Grids ! Collection of grids.
  Type(NWPMet_), Intent(InOut)        :: M     ! State of an NWP met module instance.
  Type(Units_),  Intent(InOut)        :: Units ! Collection of information on
                                               ! input/output unit numbers.
  ! Locals:
  Type(HGrid_), Pointer    :: HGrid           ! Abbreviation for grid.
  Integer                  :: Unit            ! Input/output unit number (for Name II and PP files).
  Integer                  :: iFile           ! File identifier          (for a GRIB file).
  Integer                  :: IOStat          ! Error code for open statement.
  Integer                  :: Error           ! Error code for allocate statement.
  Integer                  :: Status          ! Status code from GRIB_API subroutines.
  Character(MaxCharLength) :: ErrMessage      ! Error message from GRIB_API subroutines.
  Integer                  :: nMessages       ! Number of GRIB messages in topography file.
  Integer                  :: gribID          ! Identifier for GRIB message.
  Integer                  :: nX              !} Size of horizontal grid encoded in GRIB message.
  Integer                  :: nY              !}
  Integer                  :: numberOfValues  ! Size of data values array encoded in GRIB message.
  Real(Std), Allocatable   :: LocalField(:)   ! Temporary array for reading data values in GRIB message.
  Integer                  :: IHeader(45)     !} Header record.
  Real(Std)                :: RHeader(19)     !}
  Integer                  :: i               !] Loop indices.
  Integer                  :: j               !]
  Integer                  :: k               ! Index counter.
  Integer                  :: varid           ! NetCDF variable ID returned from nf90_inq_varid
  Integer                  :: start(3)        ! The start and count arrays for reading data
  Integer                  :: count(3)        ! using netCDF library routines.

  ! Set up abbreviation for grid.
  HGrid => Grids%HGrids(M%iHGrid)

  ! Open topography file.

  ! GRIB topography file.
  If (M%MetDefn%FileType .CIEq. 'GRIB') Then

#   ifdef GRIBsupport

      ! $$ Why not read like grib met file - avoids need to open and can check is actually
      !    topog - need to input
      !    topog field code. Input topog field code would be useful for PP files too.

      ! Open topography file using special routine in the GRIB_API.
      !
      ! Use the GRIB_API subroutine 'grib_open_file(iFile, filename, mode, status)'.
      !   Action: opens a GRIB file for binary access and sets up a file identifier for use in
      !           subsequent GRIB_API calls on the file.
      !   Arguments:
      !     iFile    - unique identifier of the opened file to be used in subsequent GRIB_API calls
      !                (note this is not the same as, or compatible with, a Fortran unit number)
      !     filename - name of GRIB file to be opened
      !     mode     - 'r' for read only
      !     status   - return code: GRIB_SUCCESS on a successful execution, or other integer value
      !                               on error (use 'grib_get_error_string' for error message).
      Call grib_open_file(                                                            &
             iFile,                                                                   &
             Trim(ConvertFileName(Trim(M%TopogFolder) // Trim(M%MetDefn%TopogFile))), &
             'r',                                                                     &
             Status                                                                   &
           )
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                                     &
               'ERROR: Unable to open topography file '                                // &
               Trim(ConvertFileName(Trim(M%TopogFolder) // Trim(M%MetDefn%TopogFile))) // &
               ' (GRIB_API returns with the error message "'                           // &
               Trim(ErrMessage)                                                        // &
               '")',                                                                      &
               3                                                                          &
             )
      End If

#   endif

  Else If (M%MetDefn%FileType .CIEq. 'NetCDF') Then

#   ifdef NetCDFsupport

      ! Open existing NetCDF topography file for access
      ! nf90_nowrite specifies: open the dataset with read-only access,
      !                         buffering and caching accesses for efficiency.
      Status = nf90_open(                                                                 &
                 Trim(ConvertFileName(Trim(M%TopogFolder) // Trim(M%MetDefn%TopogFile))), &
                 nf90_nowrite,                                                            &
                 Unit                                                                     &
               )

      If (Status /= nf90_noerr) Then
          Call Message(                                                                     &
                 'ERROR: Unable to open the NetCDF topography file '                     // &
                 Trim(ConvertFileName(Trim(M%TopogFolder) // Trim(M%MetDefn%TopogFile))) // &
                 ' (NetCDF returns with the error message "'                             // &
                 Trim(nf90_strerror(Status))                                             // &
                 '")',                                                                      &
                 3                                                                          &
               )
      End If

#   endif

  ! Name II or PP topography file.
  Else

    Call GetNewUnit(Unit, Units) !$$ use openfile once have sorted out how to handle 'convert ='
#   ifdef UseConvert
      Open (                                                                               &
        Unit    = Unit,                                                                    &
        File    = Trim(ConvertFileName(Trim(M%TopogFolder) // Trim(M%MetDefn%TopogFile))), &
        Status  = 'Old',                                                                   &
        Form    = 'Unformatted',                                                           &
        Action  = 'Read',                                                                  &
        Convert = Trim(M%MetDefn%BinaryFormat),                                            &
        IOStat  = IOStat                                                                   &
      )
#   else
      Open(                                                                               &
        Unit   = Unit,                                                                    &
        File   = Trim(ConvertFileName(Trim(M%TopogFolder) // Trim(M%MetDefn%TopogFile))), &
        Status = 'Old',                                                                   &
        Form   = 'Unformatted',                                                           &
        Action = 'Read',                                                                  &
        IOStat = IOStat                                                                   &
      )
#   endif
    If (IOStat /= 0) Then
      Call Message(                                                                      &
             'ERROR: Unable to open topography file '                                 // &
             Trim(ConvertFileName(Trim(M%TopogFolder) // Trim(M%MetDefn%TopogFile))),    &
             3                                                                           &
           )
    End If

  End If

  ! Read topography.

  ! GRIB topography file.
  If (M%MetDefn%FileType .CIEq. 'GRIB') Then

#   ifdef GRIBsupport

      ! Check that topography file contains (at least) one GRIB message.
      !
      ! Use the GRIB_API subroutine 'grib_count_in_file(iFile, n, status)'.
      !   Action: counts the number of GRIB messages in a file.
      !   Arguments:
      !     iFile    - unique identifier of a file opened with 'grib_open_file'
      !     n        - number of GRIB messages in file
      !     status   - return code: GRIB_SUCCESS on a successful execution, or other integer value
      !                               on error (use 'grib_get_error_string' for error message).
      Call grib_count_in_file(iFile, nMessages, Status)
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                                     &
               'ERROR: problem determining number of fields in the topography file '   // &
               Trim(ConvertFileName(Trim(M%TopogFolder) // Trim(M%MetDefn%TopogFile))) // &
               ' (GRIB_API returns with the error message "'                           // &
               Trim(ErrMessage)                                                        // &
               '")',                                                                      &
               3                                                                          &
             )
      End If
      If (nMessages <= 0) Then
        Call Message(                                                                     &
               'ERROR: the topography file '                                           // &
               Trim(ConvertFileName(Trim(M%TopogFolder) // Trim(M%MetDefn%TopogFile))) // &
               ' appears to not contain any GRIB messages',                               &
               3                                                                          &
             )
      Else If (nMessages > 1) Then
        Call Message(                                                                     &
               'WARNING: the topography file '                                         // &
               Trim(ConvertFileName(Trim(M%TopogFolder) // Trim(M%MetDefn%TopogFile))) // &
               ' appears to contain more than one GRIB message'                        // &
               ' - only the first field present will be read',                            &
               2                                                                          &
             )
      End If
      
      ! Read GRIB message from topography file.
      !
      ! Use the GRIB_API subroutine 'grib_new_from_file(iFile, gribID, status)'.
      !   Action: reads a GRIB message, one per call, from a file into memory.
      !   Arguments:
      !     iFile    - unique identifier of a file opened with 'grib_open_file'
      !     gribID   - unique identifier of GRIB message for use by subsequent GRIB_API calls
      !     status   - return code: GRIB_SUCCESS on a successful execution, GRIB_END_OF_FILE
      !                               when called at end of file, or other integer value
      !                               on error (use 'grib_get_error_string' for error message).
      Call grib_new_from_file(iFile, gribID, Status)
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                                     &
               'ERROR: problem reading field in the topography file '                  // &
               Trim(ConvertFileName(Trim(M%TopogFolder) // Trim(M%MetDefn%TopogFile))) // &
               ' (GRIB_API returns with the error message "'                           // &
               Trim(ErrMessage)                                                        // &
               '")',                                                                      &
               3                                                                          &
             )
      End If

      ! $$ Check topography code?

      ! Check that topography field has grid size consistent with the met defn.
      !
      ! Use the GRIB_API subroutine 'grib_get(gribID, key, value, status)'.
      !   Action: gets the value of a key from a GRIB message loaded in memory.
      !   Arguments:
      !     gribID   - unique identifier of GRIB message
      !     key      - name of key
      !     value    - value of the key
      !     status   - return code: GRIB_SUCCESS on a successful execution, or other integer value
      !                               on error (use 'grib_get_error_string' for error message).
      Call grib_get(gribID, 'Ni', nX, Status)
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                         &
               'ERROR: unable to interpret grid size in topography field ' // &
               '(GRIB_API returns with the error message "'                // &
               Trim(ErrMessage)                                            // &
               '")',                                                          &
               3                                                              &
             )
      End If
      Call grib_get(gribID, 'Nj', nY, Status)
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                         &
               'ERROR: unable to interpret grid size in topography field ' // &
               '(GRIB_API returns with the error message "'                // &
               Trim(ErrMessage)                                            // &
               '")',                                                          &
               3                                                              &
             )
      End If
      
      If (nX /= HGrid%nX .or. nY /= HGrid%nY) Then
        Call Message(                                                                                  &
               'ERROR: grid size of the topography field is inconsistent with the NWP met definition', &
               3                                                                                       &
             )
      End If

      ! Actual size of the data values array in the GRIB message
      !
      ! Use the GRIB_API subroutine 'grib_get_size(gribID, key, size, status)'.
      !   Action: gets the size of an array key from a GRIB message loaded in memory.
      !   Arguments:
      !     gribID   - unique identifier of GRIB message
      !     key      - name of array key
      !     size     - size of the array key
      !     status   - return code: GRIB_SUCCESS on a successful execution, or other integer value
      !                               on error (use 'grib_get_error_string' for error message).
      Call grib_get_size(gribID, 'values', numberOfValues, Status)
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                                  &
               'ERROR: unable to interpret size of the data array in GRIB message ' // &
               '(GRIB_API returns with the error message "'                         // &
               Trim(ErrMessage)                                                     // &
               '")',                                                                   &
               3                                                                       &
             )
      End If
      
      If (numberOfValues /= HGrid%nX * HGrid%nY) Then
        Call Message(                                                              &
               'ERROR: incorrect number of data values in topography field array', &
               3                                                                   &
             )
      End If

      ! Allocate temporary array to store topography field.
      Allocate(LocalField(numberOfValues), Stat = Error)
      If (Error /= 0) Then
        Call Message('FATAL ERROR in ReadTopog: unable to allocate array to store topography field', 3)
      End If

      ! Get data values of topography field.
      Call grib_get(gribID, 'values', LocalField, Status)
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                     &
               'ERROR: unable to read data array in topography field ' // &
               '(GRIB_API returns with the error message "'            // &
               Trim(ErrMessage)                                        // &
               '")',                                                      &
               3                                                          &
             )
      End If

      ! Release the GRIB message from memory.
      !
      ! Use the GRIB_API subroutine 'grib_release(gribID, status)'.
      !   Action: free the memory of the GRIB message specified by gribID.
      !   Arguments:
      !     gribID   - unique identifier of GRIB message to be unloaded
      !     status   - return code: GRIB_SUCCESS on a successful execution, or other integer value
      !                               on error (use 'grib_get_error_string' for error message).
      Call grib_release(gribID, Status)
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                     &
               'ERROR: unable to release memory for the GRIB message ' // &
               '(GRIB_API returns with the error message "'            // &
               Trim(ErrMessage)                                        // &
               '")',                                                      &
               3                                                          &
             )
      End If

      ! Put data into topography array.
      k = 0
      Do j = 1, HGrid%nY
      Do i = 1, HGrid%nX
        k = k + 1
        M%Topog(i, j, 1) = LocalField(k)
      End Do
      End Do

      ! Deallocate temporary array.
      Deallocate(LocalField)
      
#   endif

  ! NetCDF topography file
  Else If (M%MetDefn%FileType .CIEq. 'NetCDF') Then

#   ifdef NetCDFsupport

      ! $$ Hard coded variable name to be read from the topography file (use a command line
      ! editing tool such as "nco ncrename" to change the name of the variable in the
      ! topography file to suit this NAME code).
      ! $$ Longer term the topography could be read in the same way as other met fields.
      
      Status = nf90_inq_varid(Unit, 'topography', varid)
      If (Status /= nf90_noerr) Then
        Call Message(                                                                 &
                 'ERROR: NetCDF_API returns with the error message "'              // &
                 Trim(nf90_strerror(Status))                                       // &
                 '"',                                                                 &
                 3                                                                    &
               )
      End If

      ! Indices specifying the relative location and dimensions of the topography field
      ! in the NetCDF file.
      count = (/ HGrid%nX, HGrid%nY, 1 /)
      start = (/ 1, 1, 1 /)

      ! Read topography.
      Status = nf90_get_var(Unit, varid, M%Topog, start, count)
      If (Status /= nf90_noerr) Then
        Call Message(                                                                 &
                 'ERROR: NetCDF_API returns with the error message "'              // &
                 Trim(nf90_strerror(Status))                                       // &
                 '"',                                                                 &
                 3                                                                    &
               )
      End If

#   endif

  ! PP topography file.
  Else If (M%MetDefn%FileType .CIEq. 'PP') Then

    Read (Unit = Unit, Err = 9, End = 10) IHeader, RHeader
    Read (Unit = Unit, Err = 9, End = 10) M%Topog(:, :, 1)

  ! Name II topography file.
  Else

    Do j = 1, HGrid%nY
      Read (Unit = Unit, Err = 9, End = 10) (M%Topog(i, j, 1), i = 1, HGrid%nX)
    End Do

  End If

  ! Topography for the first time step is the same as for the second.
  M%Topog(:, :, 2) = M%Topog(:, :, 1)

  If (M%OpenMPOpts%ParallelMetRead) Then
    ! also the case for prefetching
    M%Topog(:, :, 3) = M%Topog(:, :, 1)
  End if

  ! Close topography file.

  ! GRIB topography file.
  If (M%MetDefn%FileType .CIEq. 'GRIB') Then

#   ifdef GRIBsupport

      ! Close topography file using special routine in the GRIB_API.
      !
      ! Use the GRIB_API subroutine 'grib_close_file(iFile, status)'.
      !   Action: closes a GRIB file previously opened with the subroutine 'grib_open_file'.
      !   Arguments:
      !     iFile    - unique identifier of the GRIB file to be closed
      !     status   - return code: GRIB_SUCCESS on a successful execution, or other integer value
      !                               on error (use 'grib_get_error_string' for error message).
      Call grib_close_file(iFile, Status)
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                          &
               'ERROR: Unable to close topography file '    // &
               '(GRIB_API returns with the error message "' // &
               Trim(ErrMessage)                             // &
               '")',                                           &
               3                                               &
             )
      End If

#   endif

  Else If (M%MetDefn%FileType .CIEq. 'NetCDF') Then

#   ifdef NetCDFsupport

      ! Close the NetCDF topography file
      Status = nf90_close(Unit)
      If (Status /= nf90_noerr) Then
        Call Message(                                             &
               'ERROR: Unable to close the NetCDF topog file ' // &
               '(NetCDF returns with the error message "'      // &
               Trim(nf90_strerror(Status))                     // &
               '")',                                              &
               3                                                  &
             )
      End If

#   endif

  ! Name II or PP topography file.
  Else

    Call CloseUnit(Unit, Units)

  End If

  Return

  ! Read error.
9 Continue

  Call Message('ERROR reading topography', 3)

10 Continue

   Call Message('ERROR reading topography. Premature end of file', 3)

End Subroutine ReadTopog

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadNWPMet(Time, iCase, iMetCase, Coords, Grids, Allowskip, Error, Skip, M, Units, IdxIn)

! Reads NWP met data for a single time.

! Note this routine reads data into the M%New time level of the data arrays.
!
! Non-fatal errors are issued for file opening errors, read errors and inconsistent file headers. Warnings are
! issued for missing data which is expected in the met definition, except where more than one 2-d field or 2-d
! (horizontal) section of a 3-d field in the met definition can fulfill a given role, and where one but not
! all of these are present.

! $$ Currently if module is invalid, those values which were obtained cannot be used later for persistence fix
!    ups. Is this sensible? (it may be, but may need to be reviewed with attribute dep validity).

  Implicit None
  ! Argument List:
  Type(Time_),   Intent(In)            :: Time     ! Time for which NWP met data is to be read.
  Integer,       Intent(In)            :: iCase    ! Number of case.
  Integer,       Intent(In)            :: iMetCase ! Index of the met realisation in the met ensemble.
  Type(Coords_), Intent(In),    Target :: Coords   ! Collection of coord systems.
  Type(Grids_),  Intent(In),    Target :: Grids    ! Collection of grids.
  Logical,       Intent(In)            :: AllowSkip ! Allow skipping of this file?
  Logical,       Intent(Out)           :: Error    ! Indicates an error occured in this routine.
  Logical,       Intent(InOut)         :: Skip     ! Indicates that the MetData file could not be found
                                                   ! and should be skipped
  Type(NWPMet_), Intent(InOut), Target :: M        ! State of an NWP met module instance.
  Type(Units_),  Intent(InOut)         :: Units    ! Collection of information on input/output unit numbers.
  Integer, Intent(In), Optional :: IdxIn

  ! Local types:
  Type :: F2_ ! 2-d (horizontal) field (with time dimension too).
    Real(Std), Pointer :: P(:,:,:) ! 2-d (horizontal) field.
  End Type F2_
  Type :: F3_ ! 3-d field (with time dimension too).
    Real(Std), Pointer :: P(:,:,:,:) ! 3-d field.
  End Type F3_
  ! Locals:
  Real(Std),allocatable,Target         :: Dummy(:,:,:)
  Real(Std),allocatable,Target         :: DummyU(:,:,:)
  Real(Std),allocatable,Target         :: DummyV(:,:,:)

  Type(F3_)                            :: F3(nFieldNames3d)
  Type(F2_)                            :: F2(nFieldNames2d)
  Integer                              :: nT
  Integer                              :: iMetCaseL
  Character(MaxFileNameLength)         :: MetFolderPart1
  Character(MaxFileNameLength)         :: MetFolderPart2
  Character(MaxFileNameLength)         :: MetFiles(MaxNWPMetDefn2s)
  Character(MaxFileNameLength)         :: MetFileFail
  Integer                              :: iMetFile
  Logical                              :: Exist
  Integer                              :: IOStat
  Integer                              :: ForecastStep
  Logical                              :: InconsistentHeader
  Logical                              :: InconsistentValidityTime
  Logical                              :: Missing
  Logical                              :: NewFile
  Integer                              :: MetFileUnits(MaxNWPMetDefn2s)
 ! Type(HCoord_),               Pointer :: HCoord
  Type(HGrid_),                Pointer :: HGrid
  Type(HGrid_),                Pointer :: HGridU
  Type(HGrid_),                Pointer :: HGridV
  Type(ZGrid_),                Pointer :: ZGrid
  Type(ZGrid_),                Pointer :: ZGridUV
  Type(ZGrid_),                Pointer :: ZGridW
  Type(ZGrid_),                Pointer :: ZGridP
  Integer                              :: i
  Integer                              :: j
  Integer                              :: k
  Integer                              :: iT
  Integer                              :: iField
  Integer                              :: iNWPFieldLevel
  Integer                              :: iNameIIIFieldLevel
  Integer                              :: iAttrib
  Character(MaxCharLength)             :: CharTime
  Integer                              :: idx
  Integer                              :: nDummyBuffers
  Integer                              :: indexML
  Integer                              :: indexSL
  Character(MaxCharLength)             :: ParameterKeysML
  Character(MaxCharLength)             :: ParameterKeysSL
  Integer                              :: NCLevelIndex     ! Physical or pseudo/soil level number obtained
                                                           ! from a counter in the loop over Metdefn2 levels
                                                           ! for indexing NetCDF variables to be read in.
                                                           ! $$ this is not really the ideal way of handling
                                                           !    pseudo levels (but would require some code
                                                           !    restructuring to do better here).
  Integer                              :: Status
  Character(MaxCharLength)             :: ErrMessage
  ! Dummy              :} Dummy 2-d (horizontal) fields (with time dimension too) with horizontal dimensions
  ! DummyU             :} corresponding to the main, U and V grids.
  ! DummyV             :}
  ! F3                 :] 3-d and 2-d (horizontal) fields (with time dimension too) used as alias for fields
  ! F2                 :] in M and for fields Dummy, DummyU and DummyV.
  ! nT                 :: Index of the required time in the met file.
  ! iMetCaseL          :: Local copy of iMetCase.
  ! MetFolderPart1     :} Met folder parts 1 and 2.
  ! MetFolderPart2     :}
  ! MetFiles           :: Met files.
  ! MetFileFail        :: Met file which cannot be opened
  ! iMetFile           :: Loop variable
  ! Exist              :: Indicates that a file exists.
  ! IOStat             :: Error code for open statement.
  ! ForecastStep       :: Forecast step of a met field in hours (used to check that all fields at a given time
  !                       have the same f/c step).
  ! InconsistentHeader :: Indicates the header for the field is inconsistent with what is expected.
  ! InconsistentValidityTime :: Indicates the validity time in the header is inconsistent with what is
  !                             expected (note this does not necessarily imply the met file is not valid).
  ! Missing            :: Indicates a field is missing from the met file.
  ! NewFile            :: Indicates a newly opened file is being considered.
  ! MetFileUnits       :: Input/output unit numbers.
  ! HCoord             :} Abbreviations for coords, grids.
  ! HGrid              :}
  ! HGridU             :}
  ! HGridV             :}
  ! ZGrid              :}
  ! ZGridUV            :}
  ! ZGridW             :}
  ! ZGridP             :}
  ! i                  :] Loop indices.
  ! j                  :]
  ! k                  :]
  ! iT                 :]
  ! iField             :]
  ! iNWPFieldLevel     :: Index of NWP field level.
  ! iNameIIIFieldLevel :: Index of Name III field level.
  ! iAttrib            :: Index of attribute.
  ! CharTime           :: Character version of Time.
  ! idx                :: 
  ! nDummyBuffers      :: Number of Dummy, DummyU and DummyV buffers.
  ! indexML            :: Identifier for index providing model-level access to a GRIB met file.
  ! indexSL            :: Identifier for index providing single-level access to a GRIB met file.
  ! ParameterKeysML    :: Comma-separated list of keys for creating GRIB file index incl. model level.
  ! ParameterKeysSL    :: Comma-separated list of keys for creating GRIB file index excl. model level.
  ! Status             :: Status code from GRIB_API subroutines.
  ! ErrMessage         :: Error message from GRIB_API subroutines.

  If (M%OpenMPOpts%ParallelMetRead) Then
    nDummyBuffers = 3
  Else
    nDummyBuffers = 2
  End If

  Allocate(Dummy(Grids%HGrids(M%iHGrid)%nX,Grids%HGrids(M%iHGrid)%nY,nDummyBuffers))
  Allocate(DummyU(Grids%HGrids(M%iHGridU)%nX,Grids%HGrids(M%iHGridU)%nY,nDummyBuffers))
  Allocate(DummyV(Grids%HGrids(M%iHGridV)%nX,Grids%HGrids(M%iHGridV)%nY,nDummyBuffers))

  If (Present(IdxIn)) Then
    idx=IdxIn
  Else
    idx=M%New
  End If

  ! 1) Set up error flag and abbreviations for grids.

  ! Set error flag to false.
  Error = .false.
  ! Set skip flag to false
  Skip = .false.

  ! Set up abbreviations for coords and grids.
  Associate(                             &
    HCoord  => Coords%HCoords(M%iHCoord) &
  )
  HGrid   => Grids%HGrids(M%iHGrid  )
  HGridU  => Grids%HGrids(M%iHGridU )
  HGridV  => Grids%HGrids(M%iHGridV )
  ZGrid   => Grids%ZGrids(M%iZGrid  )
  ZGridUV => Grids%ZGrids(M%iZGridUV)
  ZGridW  => Grids%ZGrids(M%iZGridW )
  ZGridP  => Grids%ZGrids(M%iZGridP )

  ! 2) Sort out met files.

  ! Derive met folder name.
  If (M%MetFolderDefnType == 0) Then
    MetFolderPart1 = M%MetFolder
    MetFolderPart2 = ' '
  Else If (M%MetFolderDefnType == 1) Then
    MetFolderPart1 = M%MetFolderStem
    MetFolderPart2 = Trim(Int2Char(iMetCase - 1)) // '\' !$$ Currently set up so that the first met case is
                                                         !   the control forecast in subfolder 0, followed
                                                         !   by the perturbed forecasts 1, 2, 3, ...
                                                         !   Better to control this via input file.
  Else If (M%MetFolderDefnType == 2) Then
    MetFolderPart1 = M%MetFolders(iMetCase)
    MetFolderPart2 = ' '
  Else If (M%MetFolderDefnType == 3) Then
    MetFolderPart1 = M%EnsembleMetFolder
    MetFolderPart2 = ' '
  Else
    Call Message('UNEXPECTED FATAL ERROR in ReadNWPMet', 4)
  End If
  MetFolderPart1 = ConvertFileName(MetFolderPart1)
  MetFolderPart2 = ConvertFileName(MetFolderPart2)

  ! Set iMetCaseL. This should be -1 unless the met case number is to be used in finding the met realisation
  ! in a met file with several met realisations.
  If (M%MetFolderDefnType == 0) Then
    iMetCaseL = -1
  Else If (M%MetFolderDefnType == 1) Then
    iMetCaseL = -1
  Else If (M%MetFolderDefnType == 2) Then
    iMetCaseL = -1
  Else If (M%MetFolderDefnType == 3) Then
    iMetCaseL = iMetCase - 1              !$$ Currently set up so that the first met case is the control
                                          !   forecast, followed by the perturbed forecasts 1, 2, 3, ...
                                          !   Better to control this via input file.
  End If

  ! Derive met file name and nT.
  Do iMetFile = 1, M%MetDefn%nMetDefn2Names
    If (M%MetDefn%DayPerFile) Then
      CharTime = FileNameTime(Time)
      nT = (Char2Int(CharTime(9:10)) * 3600) / Int(ShortTime2RealTime(M%Dt)) + 1
      ! $$ robust? interaction with T0 & dt. assumes midnight UTC is a met time (and hence dt divides 1 day)
      MetFiles(iMetFile) =                   &
        Trim(M%MetDefn%Prefix)            // &
        CharTime(7:8)                     // &
        CharTime(5:6)                     // &
        CharTime(3:4)                     // &
        '.'                               // &
        Trim(M%MetDefn%Suffixs(iMetFile))
    Else
      nT = 1
      MetFiles(IMetFile) =                   &
        Trim(M%MetDefn%Prefix)            // &
        Trim(FileNameTime(Time))          // &
        '.'                               // &
        Trim(M%MetDefn%Suffixs(iMetFile))
    End If

    ! Delete previous file.
    If (M%DeleteMet .and. M%OldMetFiles(iMetFile) /= MetFiles(iMetFile) .and. M%OldMetFiles(iMetFile) /= ' ') Then
      Inquire(File = Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(M%OldMetFiles(iMetFile)), Exist = Exist)
      If (Exist) Then
        Call GetNewUnit(MetFileUnits(iMetFile), Units) !$$ to use openfile,
                                                       ! need to add return-on-error option to openfile, with options
                                                       ! to give ERROR
                                                       ! or WARNING instead of FATAL ERROR.
                                                       ! Until then should trap error here and give error or warning message.
        ! Note - should be ok for GRIB files and no need to use 'convert' here as not reading data from file? $$
        Open (                                                                                                           &
          Unit   = MetFileUnits(iMetFile),                                                                               &
          File   = Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(M%OldMetFiles(iMetFile)))), &
          Status = 'Old',                                                                                                &
          Form   = 'Unformatted',                                                                                        &
          IOStat = IOStat                                                                                                &
        ) 
        Call CloseUnit(Unit = MetFileUnits(iMetFile), Units = Units, DeleteFile = .true.)
      End If
    End If

    ! Restore met file.
    If (M%RestoreMetScript /= ' ') Then
      Inquire(File = Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)), Exist = Exist)
      If (.not.Exist) Then
        Call SubmitSystemCommand(                               &
               Trim(ConvertFileName(M%RestoreMetScript))     // &
               ' '                                           // &
               Trim(MetFolderPart1) // Trim(MetFolderPart2)  // &
               ' '                                           // &
               Trim(MetFiles(iMetFile))                         &
           )
      End If
    End If

    ! Store met file name.
    M%OldMetFiles(iMetFile) = MetFiles(iMetFile)

    ! Open met file (GRIB files are 'opened' by creating an index using the relevant GRIB_API routine).
  
    ! GRIB met file.
    If (M%MetDefn%FileType .CIEq. 'GRIB') Then

      ! Check that file exists
      Inquire (                                                                                                  &
        File  = Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)))), &
        Exist = Exist                                                                                            &
      )
      If (.not. Exist) Then
        If (.not. AllowSkip) Then
          Call Message(                                                                                              &
                 "ERROR: Can't open met file "                                                                    // &
                 Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)))),    &
                 2                                                                                                   &
               )
        End If
        Skip = Allowskip
        Error = .true.
        MetFileFail = MetFiles(iMetFile)
        Go To 9
      End If

#     ifdef GRIBsupport

        ! Note that indexML and indexSL need to be assigned some arbitrary (but distinct) positive values
        ! to avoid a crash in the GRIB_API routines under certain circumstances.
        indexML = 1
        indexSL = 2
      
        ! Create two indexes for accessing i) model-level and ii) single-level fields from the GRIB met file.
        ! Two separate file indexes are created here so that the level index does not need to be specified for
        ! single-level (surface/standard-level) fields. The level index is often set to zero for such fields,
        ! but this cannot be guaranteed for all met files that might be read.
        !
        ! Use the GRIB_API subroutine 'grib_index_create(indexID, filename, keys, status)'.
        !   Action: creates an index for referencing the messages in a GRIB file using the supplied keys.
        !   Arguments:
        !     indexID  - unique identifier of the index created
        !     filename - name of the file of GRIB messages to be indexed
        !     keys     - comma-separated list of keys for creating the index
        !     status   - return code: GRIB_SUCCESS on a successful execution, or other integer value
        !                               on error (use 'grib_get_error_string' for error message).
      
        ! Set indexing keys - should include ensemble member when the met case number is to be used in
        ! finding the met realisation in a met file with several met realisations.
        If (iMetCaseL /= -1) Then
          ParameterKeysML = 'validityDate,validityTime,paramId,perturbationNumber,level'
          ParameterKeysSL = 'validityDate,validityTime,paramId,perturbationNumber'
        Else
          ParameterKeysML = 'validityDate,validityTime,paramId,level'
          ParameterKeysSL = 'validityDate,validityTime,paramId'
        End If
      
        ! Create index providing model-level access.
        Call grib_index_create(                                                                                 &
               indexML,                                                                                         &
               Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)))), &
               ParameterKeysML, Status                                                                          &
             )
      
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                                                             &
                 'ERROR: Unable to create index for model-level access to the GRIB met file '                    // &
                 Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)))) // &
                 ' (GRIB_API returns with the error message "'                                                   // &
                 Trim(ErrMessage)                                                                                // &
                 '")',                                                                                              &
                 2                                                                                                  &
               )
          Error = .true.
          Go To 9
        End If

        ! Create index providing single-level access.
        Call grib_index_create(                                                                                 &
               indexSL,                                                                                         &
               Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)))), &
               ParameterKeysSL, Status                                                                          &
             )
      
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                                                             &
                 'ERROR: Unable to create index for single-level access to the GRIB met file '                   // &
                 Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)))) // &
                 ' (GRIB_API returns with the error message "'                                                   // &
                 Trim(ErrMessage)                                                                                // &
                 '")',                                                                                              &
                 2                                                                                                  &
               )
          Error = .true.
          Go To 9
        End If

#     endif

    Else If (M%MetDefn%FileType .CIEq. 'NetCDF') Then

#     ifdef NetCDFsupport

      ! Open the NetCDF file using MetFileUnits to hold ncid
      ! nf90_nowrite specifies: open the dataset with read-only access,
      !                         buffering and caching accesses for efficiency.
      Status = nf90_open(                                                                                         &
                 Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)))), &
                 nf90_nowrite,                                                                                    &
                 MetFileUnits(iMetFile)                                                                           &
               )

      If (Status /= nf90_noerr) Then
        If (.not.AllowSkip) Then
          Call Message(                                                                                             &
                 'ERROR: Unable to open the NetCDF met file '                                                    // &
                 Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)))) // &
                 ' (NetCDF returns with the error message "'                                                     // &
                 Trim(nf90_strerror(Status))                                                                     // &
                 '")',                                                                                              &
                 2                                                                                                  &
               )
        End If
        Skip = Allowskip
        Error = .true.
        MetFileFail = MetFiles(iMetFile)
        Go To 9
      End If

#     endif

    ! Name II or PP met file.
    Else

      Call GetNewUnit(MetFileUnits(iMetFile), Units) !$$ use openfile once have sorted out how to handle 'convert ='
#     ifdef UseConvert
        Open (                                                                                                       &
          Unit    = MetFileUnits(iMetFile),                                                                          &
          File    = Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)))), &
          Status  = 'Old',                                                                                           &
          Form    = 'Unformatted',                                                                                   &
          Action  = 'Read',                                                                                          &
          Convert = Trim(M%MetDefn%BinaryFormat),                                                                    &
          IOStat  = IOStat                                                                                           &
        )
#     else
        Open (                                                                                                      &
          Unit   = MetFileUnits(iMetFile),                                                                          &
          File   = Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)))), &
          Status = 'Old',                                                                                           &
          Form   = 'Unformatted',                                                                                   &
          Action = 'Read',                                                                                          &
          IOStat = IOStat                                                                                           &
        )
#     endif
      If (IOStat /= 0) Then
        Call CloseUnit(MetFileUnits(iMetFile), Units)
        If (.not. AllowSkip) Then
          Call Message(                                                                                              &
                 "ERROR: Can't open met file "                                                                    // &
                 Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)))),    &
                 2                                                                                                   &
               )
        End If
        Skip = Allowskip
        Error = .true.
        MetFileFail = MetFiles(iMetFile)
        Go To 9
      End If

    End If
  End Do

  ! 3) Set up equivalences to enable looping over fields and initialise FieldPresent3d, FieldPresent2d and
  ! NewFile.

  F3(F_U                   )%P => M%U
  F3(F_V                   )%P => M%V
  F3(F_W                   )%P => M%W
  F3(F_T                   )%P => M%T
  F3(F_Q                   )%P => M%Q
  F3(F_PAsRead             )%P => M%PAsRead
  F3(F_TotalOrDynCloudWater)%P => M%TotalOrDynCloudWater
  F3(F_TotalOrDynCloudIce  )%P => M%TotalOrDynCloudIce
  F3(F_CanopyHeight        )%P => M%CanopyHeight
  F3(F_CanopyWater         )%P => M%CanopyWater
  F3(F_StomataConduct      )%P => M%StomataConduct

  F2(F_Z0                 )%P => M%Z0
  F2(F_UStress            )%P => M%UStress
  F2(F_VStress            )%P => M%VStress
  F2(F_HeatFlux           )%P => M%HeatFlux
  F2(F_H                  )%P => M%H
  F2(F_PSeaLevel          )%P => M%PSeaLevel
  F2(F_ConCloud           )%P => M%ConCloud
  F2(F_ConCloudBase       )%P => M%ConCloudBase
  F2(F_ConCloudTop        )%P => M%ConCloudTop
  F2(F_TotalOrDynHCloud   )%P => M%TotalOrDynHCloud
  F2(F_TotalOrDynMCloud   )%P => M%TotalOrDynMCloud
  F2(F_TotalOrDynLCloud   )%P => M%TotalOrDynLCloud
  F2(F_ConRain            )%P => M%ConRain
  F2(F_ConSnow            )%P => M%ConSnow
  F2(F_DynRain            )%P => M%DynRain
  F2(F_DynSnow            )%P => M%DynSnow
  F2(F_Dummy              )%P => Dummy
  F2(F_DummyU             )%P => DummyU
  F2(F_DummyV             )%P => DummyV
  F2(F_SoilMoisture       )%P => M%SoilMoisture

  M%FieldPresent3d(:, :) = .false.
  M%FieldPresent2d(:)    = .false.

  LoopOverMetFiles: Do iMetFile = 1, M%MetDefn%nMetDefn2Names

    NewFile = .true.

    ! 4) Read data (note missing data is ignored at times which don't correspond to the time of interest).

    ! $$ Note currently all files are read in prescribed order. Could change for PP met files and
    !    possibly others.

    ! Reading fields in the order they are in M%MetDefn2s(iMetFile).
    If (.true.) Then

      LoopOverTimesInFile: Do iT = 1, nT

        ! Initialise forecast step variable to indicate that it is not yet set.
        ForecastStep = -1

        LoopOverFieldsAtGivenTime: Do iField = 1, M%MetDefn2s(iMetFile)%nFields

          ! 3-d fields.
          NCLevelIndex = 1  ! Initialise NetCDF physical/pseudo/soil level number for both 2-d and 3-d fields

          If (M%ThreeD(iMetFile,iField)) Then

            LoopOverLevels: Do k = M%k1(iMetFile,iField), M%k2(iMetFile,iField)

              Select Case (M%iField(iMetFile,iField))
                ! U, V grids: U, V .
                Case (F_U, F_V)
                  iNWPFieldLevel = ZGridUV%iZ(k)
                ! W grid: W.
                Case (F_W)
                  iNWPFieldLevel = ZGridW%iZ(k)
                ! P grid: P.
                Case (F_PAsRead)
                  iNWPFieldLevel = ZGridP%iZ(k)
                ! Main grid: all other fields.
                Case Default
                  iNWPFieldLevel = ZGrid%iZ(k)
              End Select

              Select Case (M%iField(iMetFile,iField))
                ! U grid: U.
                Case (F_U)
                  Call Read2dField(                                                                      &
                         NewFile,                                                                        &
                         M%MetDefn%FileType,                                                             &
                         Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)),       &
                         indexML, indexSL, MetFileUnits(iMetFile),                                       &
                         M%MetDefn2s(iMetFile)%FieldCodes(iField), M%MetDefn2s(iMetFile)%ThreeD(iField), &
                         M%MetDefn2s(iMetFile)%NCFieldNames(iField), NCLevelIndex,                       &
                         iNWPFieldLevel, iMetCaseL, Time,                                                &
                         HCoord, HGridU,                                                                 &
                         F3(M%iField(iMetFile,iField))%P(:, :, k, idx),                                  &
                         Error, InconsistentHeader, InconsistentValidityTime, Missing,                   &
                         ForecastStep                                                                    &
                       )
                ! V grid: V.
                Case (F_V)
                  Call Read2dField(                                                                      &
                         NewFile,                                                                        &
                         M%MetDefn%FileType,                                                             &
                         Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)),       &
                         indexML, indexSL, MetFileUnits(iMetFile),                                       &
                         M%MetDefn2s(iMetFile)%FieldCodes(iField), M%MetDefn2s(iMetFile)%ThreeD(iField), &
                         M%MetDefn2s(iMetFile)%NCFieldNames(iField), NCLevelIndex,                       &
                         iNWPFieldLevel, iMetCaseL, Time,                                                &
                         HCoord, HGridV,                                                                 &
                         F3(M%iField(iMetFile,iField))%P(:, :, k, idx),                                  &
                         Error, InconsistentHeader, InconsistentValidityTime, Missing,                   &
                         ForecastStep                                                                    &
                       )
                ! Main grid: all other fields.
                Case Default
                  Call Read2dField(                                                                      &
                         NewFile,                                                                        &
                         M%MetDefn%FileType,                                                             &
                         Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)),       &
                         indexML, indexSL, MetFileUnits(iMetFile),                                       &
                         M%MetDefn2s(iMetFile)%FieldCodes(iField), M%MetDefn2s(iMetFile)%ThreeD(iField), &
                         M%MetDefn2s(iMetFile)%NCFieldNames(iField), NCLevelIndex,                       &
                         iNWPFieldLevel, iMetCaseL, Time,                                                &
                         HCoord, HGrid,                                                                  &
                         F3(M%iField(iMetFile,iField))%P(:, :, k, idx),                                  &
                         Error, InconsistentHeader, InconsistentValidityTime, Missing,                   &
                         ForecastStep                                                                    &
                     )
              End Select

              ! Read errors.
              If (Error) Then
                Call Message(                                            &
                       'ERROR: read error occured in reading field '  // &
                       Trim(M%MetDefn2s(iMetFile)%FieldNames(iField)) // &
                       ' from the met file '                          // &
                       Trim(MetFiles(iMetFile))                       // &
                       ' for the '                                    // &
                       Trim(Int2Char(k)) // Int2Ordinal(k)            // &
                       ' (Name III) vertical level and the '          // &
                       Trim(Int2Char(iT)) // Int2Ordinal(iT)          // &
                       ' time level in the file',                        &
                       2                                                 &
                     )
                Go To 9
              End If

              ! Inconsistent header.
              If (InconsistentHeader .and. iT == nT) Then
                Call ControlledMessage(                                                                      &
                       'The field '                                                                       // &
                       Trim(M%MetDefn2s(iMetFile)%FieldNames(iField))                                     // &
                       ' from met file '                                                                  // &
                       Trim(MetFiles(iMetFile))                                                           // &
                       ' for the '                                                                        // &
                       Trim(Int2Char(k)) // Int2Ordinal(k)                                                // &
                       ' (Name III) vertical level and the '                                              // &
                       Trim(Int2Char(iT)) // Int2Ordinal(iT)                                              // &
                       ' time level in the file has a header which is inconsistent with the information ' // &
                       'given in the NWP Met Definition '                                                 // &
                       Trim(M%MetDefn%Name)                                                               // &
                       '.',                                                                                  &
                       MessageControls    = GlobalMessageControls,                                           &
                       MessageControlName = 'Inconsistent NWP header',                                       &
                       ErrorCode          = 2                                                                &
                     ) ! $$ Set Error = .true. and goto 9? or downgrade to warning
              End If

              ! Inconsistent validity time.
              If (InconsistentValidityTime) Then
                Call ControlledMessage(                                                                      &
                       'The field '                                                                       // &
                       Trim(M%MetDefn2s(iMetFile)%FieldNames(iField))                                     // &
                       ' from met file '                                                                  // &
                       Trim(MetFiles(iMetFile))                                                           // &
                       ' for the '                                                                        // &
                       Trim(Int2Char(k)) // Int2Ordinal(k)                                                // &
                       ' (Name III) vertical level and the '                                              // &
                       Trim(Int2Char(iT)) // Int2Ordinal(iT)                                              // &
                       ' time level in the file has a header validity time which is inconsistent with '   // &
                       'the time in the met file name.',                                                     &
                       MessageControls    = GlobalMessageControls,                                           &
                       MessageControlName = 'Inconsistent validity time',                                    &
                       ErrorCode          = 2                                                                &
                     ) ! $$ Set Error = .true. and goto 9? or downgrade to warning
              End If

              ! FieldPresent3d.
              If (iT == nT .and. .not.Missing) M%FieldPresent3d(k, M%iField(iMetFile,iField)) = .true.

              ! NewFile.
              NewFile = .false.

              NCLevelIndex = NCLevelIndex + 1 ! Increment NetCDF level number

            End Do LoopOverLevels

          ! 2-d fields.
          Else

            Select Case (M%iField(iMetFile,iField))
              ! U grid: UStress.
              Case (F_UStress, F_DummyU)
                Call Read2dField(                                                                      &
                       NewFile,                                                                        &
                       M%MetDefn%FileType,                                                             &
                       Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)),       &
                       indexML, indexSL, MetFileUnits(iMetFile),                                       &
                       M%MetDefn2s(iMetFile)%FieldCodes(iField), M%MetDefn2s(iMetFile)%ThreeD(iField), &
                       M%MetDefn2s(iMetFile)%NCFieldNames(iField), 1,                                  &
                       0, iMetCaseL, Time, HCoord, HGridU,                                             &
                       F2(M%iField(iMetFile,iField))%P(:, :, idx),                                     &
                       Error, InconsistentHeader, InconsistentValidityTime, Missing,                   &
                       ForecastStep                                                                    &
                     )
              ! V grid: VStress.
              Case (F_VStress, F_DummyV)
                Call Read2dField(                                                                      &
                       NewFile,                                                                        &
                       M%MetDefn%FileType,                                                             &
                       Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)),       &
                       indexML, indexSL, MetFileUnits(iMetFile),                                       &
                       M%MetDefn2s(iMetFile)%FieldCodes(iField), M%MetDefn2s(iMetFile)%ThreeD(iField), &
                       M%MetDefn2s(iMetFile)%NCFieldNames(iField), 1,                                  &
                       0, iMetCaseL, Time, HCoord, HGridV,                                             &
                       F2(M%iField(iMetFile,iField))%P(:, :, idx),                                     &
                       Error, InconsistentHeader, InconsistentValidityTime, Missing,                   &
                       ForecastStep                                                                    &
                     )
              ! Main grid: all other fields.
              Case Default
                Call Read2dField(                                                                      &
                       NewFile,                                                                        &
                       M%MetDefn%FileType,                                                             &
                       Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile)),       &
                       indexML, indexSL, MetFileUnits(iMetFile),                                       &
                       M%MetDefn2s(iMetFile)%FieldCodes(iField), M%MetDefn2s(iMetFile)%ThreeD(iField), &
                       M%MetDefn2s(iMetFile)%NCFieldNames(iField), 1,                                  &
                       0, iMetCaseL, Time, HCoord, HGrid,                                              &
                       F2(M%iField(iMetFile,iField))%P(:, :, idx),                                     &
                       Error, InconsistentHeader, InconsistentValidityTime, Missing,                   &
                       ForecastStep                                                                    &
                     )
            End Select

            ! Read errors.
            If (Error) Then
              Call Message(                                            &
                     'ERROR: read error occured in reading field '  // &
                     Trim(M%MetDefn2s(iMetFile)%FieldNames(iField)) // &
                     ' from the met file '                          // &
                     Trim(MetFiles(iMetFile))                       // &
                     ' for the '                                    // &
                     Trim(Int2Char(iT)) // Int2Ordinal(iT)          // &
                     ' time level in the file',                        &
                     2                                                 &
                   )
              Go To 9
            End If

            ! Inconsistent header.
            If (InconsistentHeader .and. iT == nT) Then
              Call ControlledMessage(                                                                      &
                     'The field '                                                                       // &
                     Trim(M%MetDefn2s(iMetFile)%FieldNames(iField))                                     // &
                     ' from met file '                                                                  // &
                     Trim(MetFiles(iMetFile))                                                           // &
                     ' for the '                                                                        // &
                     Trim(Int2Char(iT)) // Int2Ordinal(iT)                                              // &
                     ' time level in the file has a header which is inconsistent with the information ' // &
                     'given in the NWP Met Definition '                                                 // &
                     Trim(M%MetDefn%Name)                                                               // &
                     '.',                                                                                  &
                     MessageControls    = GlobalMessageControls,                                           &
                     MessageControlName = 'Inconsistent NWP header',                                       &
                     ErrorCode          = 2                                                                &
                   ) ! $$ Set Error = .true. and goto 9? or downgrade to warning
            End If

            ! Inconsistent validity time.
            If (InconsistentValidityTime) Then
              Call ControlledMessage(                                                                      &
                     'The field '                                                                       // &
                     Trim(M%MetDefn2s(iMetFile)%FieldNames(iField))                                     // &
                     ' from met file '                                                                  // &
                     Trim(MetFiles(iMetFile))                                                           // &
                     ' for the '                                                                        // &
                     Trim(Int2Char(iT)) // Int2Ordinal(iT)                                              // &
                     ' time level in the file has a header validity time which is inconsistent with '   // &
                     'the time in the met file name.',                                                     &
                     MessageControls    = GlobalMessageControls,                                           &
                     MessageControlName = 'Inconsistent validity time',                                    &
                     ErrorCode          = 2                                                                &
                   ) ! $$ Set Error = .true. and goto 9? or downgrade to warning
            End If

            ! FieldPresent2d.
            If (iT == nT .and. .not.Missing) M%FieldPresent2d(M%iField(iMetFile,iField)) = .true.

            ! NewFile.
            NewFile = .false.

          End If

        End Do LoopOverFieldsAtGivenTime

      End Do LoopOverTimesInFile

    ! Reading fields in the order they are in the file.
    Else

    End If

    ! Close met file (or delete its index if a GRIB file).
  
    ! GRIB met file.
    If (M%MetDefn%FileType .CIEq. 'GRIB') Then

#     ifdef GRIBsupport

        ! Delete the two indexes created for the GRIB met file.
        !
        ! Use the GRIB_API subroutine 'grib_index_release(indexID, status)'.
        !   Action: deletes the index created from a GRIB file using 'grib_index_create'.
        !   Arguments:
        !     indexID  - unique identifier of an index created from a GRIB file
        !     status   - return code: GRIB_SUCCESS on a successful execution, or other integer value
        !                               on error (use 'grib_get_error_string' for error message).
        Call grib_index_release(indexML, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                                                  &
                 'ERROR: Unable to delete index created for model-level access to the GRIB met file ' // &
                 '(GRIB_API returns with the error message "'                                         // &
                 Trim(ErrMessage)                                                                     // &
                 '")',                                                                                   &
                 3                                                                                       &
               )
          Go To 9
        End If

        Call grib_index_release(indexSL, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                                                   &
                 'ERROR: Unable to delete index created for single-level access to the GRIB met file ' // &
                 '(GRIB_API returns with the error message "'                                          // &
                 Trim(ErrMessage)                                                                      // &
                 '")',                                                                                    &
                 3                                                                                        &
               )
          Go To 9
        End If

#     endif

    Else If (M%MetDefn%FileType .CIEq. 'NetCDF') Then

#     ifdef NetCDFsupport

        ! Close the NetCDF Met file with ID stored in MetFileUnits(iMetFile)
        Status = nf90_close(MetFileUnits(iMetFile))
        If (Status /= nf90_noerr) Then
          Call Message(                                           &
                 'ERROR: Unable to close the NetCDF met file ' // &
                 '(NetCDF returns with the error message "'    // &
                 Trim(nf90_strerror(status))                   // &
                 '")',                                            &
                 3                                                &
               )
          Go To 9
        End If

#     endif

    ! Name II or PP met file.
    Else

      Call CloseUnit(MetFileUnits(iMetFile), Units)

    End If

    ! 5) Warnings for missing fields which are expected to be present according to MetDefn2s. Note this is done
    ! here rather than as each field is read, to avoid producing warnings when more than one possible NWP field
    ! can be used for a given NAME III field and some but not all of these NWP fields are missing.

    Do iField = 1, M%MetDefn2s(iMetFile)%nFields

      ! 3-d fields.
      If (M%ThreeD(iMetFile,iField)) Then

        Do k = M%k1(iMetFile,iField), M%k2(iMetFile,iField)

          If (.not. M%FieldPresent3d(k, M%iField(iMetFile,iField))) Then
            Call ControlledMessage(                                            &
                   'Field '                                              //    &
                   Trim(M%MetDefn2s(iMetFile)%FieldNames(iField))        //    &
                   ' missing from met file '                             //    &
                   Trim(MetFiles(iMetFile))                              //    &
                   ' for the '                                           //    &
                   Trim(Int2Char(k)) // Int2Ordinal(k)                   //    &
                   ' (Name III) vertical level and the '                 //    &
                   Trim(Int2Char(nT)) // Int2Ordinal(nT)                 //    &
                   ' time level in the file',                                  &
                   MessageControls    = GlobalMessageControls,                 &
                   MessageControlName = 'Missing NWP field',                   &
                   ErrorCode          = 1                                      &
                 )
          End If

        End Do

      ! 2-d fields.
      Else

        If (.not. M%FieldPresent2d(M%iField(iMetFile,iField))) Then
          Call ControlledMessage(                                         &
                 'Field '                                           //    &
                 Trim(M%MetDefn2s(iMetFile)%FieldNames(iField))     //    &
                 ' missing from met file '                          //    &
                 Trim(MetFiles(iMetFile))                           //    &
                 ' for the '                                        //    &
                 Trim(Int2Char(nT)) // Int2Ordinal(nT)              //    &
                 ' time level in the file',                               &
                 MessageControls    = GlobalMessageControls,              &
                 MessageControlName = 'Missing NWP field',                &
                 ErrorCode          = 1                                   &
               )
        End If

      End If

    End Do

  End Do LoopOverMetFIles


  ! 6) Check adequate data present. Note special treatment of pressure and sea level pressure.
  ! [The requirements for pressure are minimal. The most that is needed is one level of pressure which is used
  ! to estimate the surface pressure, with the whole pressure field then being reconstructed from surface
  ! pressure and temperature. ExtrapolationLevelsAbove and InterpolationLevels must be set large to prevent
  ! the model thinking there is inadequate data where the data are not required. ExtrapolationLevelsBelow
  ! controls how high a level can be used to estimate the surface pressure, PersistTimes3d controls how long
  ! surface pressure can be fixed up using persistence, and AllowFixUpToDefault3d controls whether sea level
  ! pressure can be used to estimate the surface pressure when no appropriate pressure levels are available.
  ! AllowFixUpToDefault2d controls whether sea level pressure can be fixed up to a default (and then possibly
  ! used to estimate the surface pressure). Provided surface pressure can be estimated, all pressure levels
  ! are computed from surface pressure and temperature (and the settings of ExtrapolationLevelsBelow etc do
  ! not affect this).]

  ! Initialise validity flags.
  M%ValidAttribs(1:nAttribs, idx) = .true.

  ! LowestPresent and HighestPresent.
  M%LowestPresent(:)  = Huge(M%LowestPresent)/3
  M%HighestPresent(:) = -Huge(M%HighestPresent)/3
  Do iField = 1, nFieldNames3d
    Do k = 1, M%kMax(iField)
      If (M%FieldPresent3d(k, iField)) Then
        M%LowestPresent(iField) = k
        Exit
      End If
    End Do
    Do k = M%kMax(iField), 1, -1
      If (M%FieldPresent3d(k, iField)) Then
        M%HighestPresent(iField) = k
        Exit
      End If
    End Do
  End Do

  ! kLower and kUpper.
  Do iField = 1, nFieldNames3d
    Do k = M%LowestPresent(iField), M%HighestPresent(iField) - 1
      If (M%FieldPresent3d(k, iField)) Then
        iNameIIIFieldLevel = k
      Else
        M%kLower(k, iField) = iNameIIIFieldLevel
      End If
    End Do
    Do k = M%HighestPresent(iField), M%LowestPresent(iField) + 1, -1
      If (M%FieldPresent3d(k, iField)) Then
        iNameIIIFieldLevel = k
      Else
        M%kUpper(k, iField) = iNameIIIFieldLevel
      End If
    End Do
  End Do

  ! 3-d fields.
  Do iField = 1, nFieldNames3d
    Do k = 1, M%kMax(iField)
      ! Avoid checking for pressure above bottom level.
      If (iField == F_PAsRead .and. k == 2) Exit
      ! Cycle if quantity present or can be fixed up. Note the test for FieldPersist3d >= 0 which prevents
      ! default values being persisted.
      If (M%FieldPresent3d(k, iField)) Cycle
      If (                                                                    &
        k < M%LowestPresent(iField)                                     .and. &
        k >= M%LowestPresent(iField) - ExtrapolationLevelsBelow(iField)       &
      ) Cycle
      If (                                                                     &
        k > M%HighestPresent(iField)                                     .and. &
        k <= M%HighestPresent(iField) + ExtrapolationLevelsAbove(iField)       &
      ) Cycle
      If (                                                                                 &
        k > M%LowestPresent(iField)                                                  .and. &
        k < M%HighestPresent(iField)                                                 .and. &
        M%kUpper(k, iField) - M%kLower(k, iField) <= InterpolationLevels(iField) + 1       &
      ) Cycle
      If (                                                         &
        M%FieldPersist3d(k, iField) < PersistTimes3d(iField) .and. &
        M%FieldPersist3d(k, iField) >= 0                     .and. &
        M%Old /= 0                                                 &
      ) Cycle
      If (AllowFixUpToDefault3d(iField)) Cycle
      Do iAttrib = 1, nAttribs
        If (Needed3d(iAttrib, iField)) M%ValidAttribs(iAttrib, idx) = .false.
      End Do
!      ! Set Error and write error message.
!      Error = .true.
!      Call Message(                                              &
!             'ERROR: Field '                                  // &
!             Trim(FieldNames3d(iField))                       // &
!             ' missing from met file '                        // &
!             Trim(MetFile)                                    // &
!             ' at (Name III) level '                          // &
!             Trim(Int2Char(k))                                // &
!             ' (and is required at this level and cannot be ' // &
!             'fixed up from other available information).',      &
!             2                                                   &
!           )
    End Do
  End Do

  ! 2-d fields.
  Do iField = 1, nFieldNames2d
    ! Cycle if quantity present or can be fixed up. Note the test for FieldPersist2d >= 0 which prevents
    ! default values being persisted.
    If (M%FieldPresent2d(iField)) Cycle
    If (                                                      &
      M%FieldPersist2d(iField) < PersistTimes2d(iField) .and. &
      M%FieldPersist2d(iField) >= 0                     .and. &
      M%Old /= 0                                              &
    ) Cycle
    If (AllowFixUpToDefault2d(iField)) Cycle
    Do iAttrib = 1, nAttribs
      If (Needed2d(iAttrib, iField)) M%ValidAttribs(iAttrib, idx) = .false.
    End Do
!    ! Set Error and write error message.
!    Error = .true.
!    Call Message(                                            &
!           'ERROR: Field '                                // &
!           Trim(FieldNames2d(iField))                     // &
!           ' missing from met file '                      // &
!           Trim(MetFile)                                  // &
!           ' (and is required and cannot be '             // &
!           'fixed up from other available information).',    &
!           2                                                 &
!         )
  End Do

  ! 6) Message.

  Do iMetFile = 1, M%MetDefn%nMetDefn2Names
    Call Message(                                                                      &
           'NWP met data for '                                                      // &
           Trim(Time2Char(Time, .false., 0, .true.))                                // &
           ' read from file '                                                       // &
           Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFiles(iMetFile))    &
         )
  End Do


  Return

  ! 7) Errors.

9 Continue

  If (Skip) Then
    Call Message(                                             &
           'WARNING: Cannot find met data for met module ' // &
           Trim(M%C%MetModName)                            // &
           '.'                                             // &
           Trim(M%C%MetName)                               // &
           ' skipping file '                               // &
           Trim(MetFolderPart1) // Trim(MetFolderPart2)    // &
           Trim(MetFileFail),                                 &
           1                                                  &
          )
  Else
    Call Message(                                          &
           'ERROR: in reading met data for met module ' // &
           Trim(M%C%MetModName)                         // &
           '.'                                          // &
           Trim(M%C%MetName),                              &
           2                                               &
         )
  End If

  M%ValidAttribs(1:nAttribs, idx) = .false.

  End Associate
  
End Subroutine ReadNWPMet

!-------------------------------------------------------------------------------------------------------------

Subroutine Read2dField(                                                    &
             NewFile,                                                      &
             FileType,                                                     &
             FileName,                                                     &
             indexML, indexSL, Unit,                                       &
             FieldCode, ThreeD, NCFieldName, NCLevelIndex,                 &
             iFieldLevel, iMetCase, Time,                                  &
             HCoord, HGrid,                                                &
             Field,                                                        &
             Error, InconsistentHeader, InconsistentValidityTime, Missing, &
             ForecastStep                                                  &
           )
! Reads in a 2-d (horizontal) field of met data from Name II, PP or GRIB met files.

  Implicit None
  ! Argument List:
  Logical,       Intent(In)         :: NewFile       ! Indicates a newly opened file is being considered.
  Character(*),  Intent(In)         :: FileType      ! Type of met file (Name II, PP, GRIB or NetCDF).
  Character(*),  Intent(In)         :: FileName      ! Name of met file (used in checking Name II met files).
  Integer,       Intent(In), Target :: indexML       ! Identifier for index providing model-level access to
                                                     ! a GRIB met file.
  Integer,       Intent(In), Target :: indexSL       ! Identifier for index providing single-level access to
                                                     ! a GRIB met file.
  Integer,       Intent(In)         :: Unit          ! Input/output unit number of met file (used for Name II
                                                     ! and PP met files).
  Integer,       Intent(In)         :: FieldCode     ! NWP code of the field. For Name II and PP met files
                                                     ! this is the 'stash' code. For GRIB met files this is
                                                     ! the 'unique parameter identifier (paramId)'.
  Logical,       Intent(In)         :: ThreeD        ! Indicates that the 2-d field is a section of a 3-d NWP
                                                     ! field.
  Character(*),  Intent(In)         :: NCFieldName   ! The NetCDF name of the variable (NetCDF files only).
  Integer,       Intent(In)         :: NCLevelIndex  ! Level "count" for the netcdf variable indexing. Use
                                                     ! instead of iFieldLevel.
                                                     ! $$ this is not really the ideal way of handling
                                                     !    pseudo levels (but would require some code
                                                     !    restructuring to do better here).
  Integer,       Intent(In)         :: iFieldLevel   ! Index of field level for sections of 3-d NWP fields.
                                                     ! For GRIB met files, the value of iFieldLevel determines
                                                     ! whether the model-level file index (iFieldLevel /= 0)
                                                     ! or single-level file index (iFieldLevel = 0) is used.
  Integer,       Intent(In)         :: iMetCase      ! Index of the met realisation (used with GRIB met files
                                                     ! only).
  Type(Time_),   Intent(In)         :: Time          ! Time of the met field.
  Type(HCoord_), Intent(In)         :: HCoord        ! Coord system in which HGrid is defined.
  Type(HGrid_),  Intent(In)         :: HGrid         ! Grid on which the field is stored.
  Real(Std),     Intent(Out)        :: Field(:, :)   ! 2-d (horizontal) field.
  Logical,       Intent(Out)        :: Error         ! Indicates a read error.
  Logical,       Intent(Out)        :: InconsistentHeader ! Indicates the header for the field is inconsistent
                                                          ! with what is expected.
  Logical,       Intent(Out)        :: InconsistentValidityTime ! Indicates the validity time from the field
                                                                ! header is inconsistent with what is expected
                                                                ! (note this does not necessarily imply the
                                                                ! met file is not valid, e.g. an existing met
                                                                ! file may have been duplicated using
                                                                ! a new time in the filename).
  Logical,       Intent(Out)        :: Missing       ! Indicates a missing field.
  Integer,       Intent(InOut)      :: ForecastStep  ! Forecast step of the met field in hours (used to check
                                                     ! that all fields at a given time have the
                                                     ! same f/c step).
  ! Locals:
  Integer,                 Save :: IHeader(45)                      !} Header record for PP met files.
  !$OMP THREADPRIVATE(IHeader)
  Real(Std),               Save :: RHeader(19)                      !}
  !$OMP THREADPRIVATE(RHeader)
  Logical,                 Save :: MissingL                         ! Indicates the last field was missing
  !$OMP THREADPRIVATE(MissingL)
                                                                    ! for PP met files.
  Real(Std)                     :: Header(64)                       !} Header record for Name II met files.
  Real(Std)                     :: Scale                            ! Tolerance for checking headers.
  Real(Std)                     :: LocalField(HGrid%nX, HGrid%nY)   ! $$ needed on Intel for Read
                                                                    !    to work correctly.
  Real(Std)                     :: LocalField2(HGrid%nX * HGrid%nY) ! Local array for reading met field.
  Integer                       :: PackingAccuracy                  ! WGDOS packing accuracy
  Integer                       :: lengthOfPackedData               ! Size of packed WGDOS data in 32 bit words
  Integer(I64), Allocatable     :: PackedField(:)                   ! Temporary 64 bit integer field for
                                                                    ! storing WGDOS packed data.
  Integer                       :: ErrorCode                        ! Error code from Allocate PackedField.
  Character(MaxCharLength)      :: CharTime                         ! Character version of Time.
  Integer                       :: i                                !} Loop indices.
  Integer                       :: j                                !}
  Integer                       :: k                                ! Counter.
  Integer,                 Save :: IOStat                           ! Error code for read statement.
  !$OMP THREADPRIVATE(IOStat)
  Integer                       :: ForecastStepL                    ! Local copy of the forecast step.
  Integer, Pointer              :: indexID                          ! Identifier for a GRIB file index.
  Integer                       :: gribID                           ! Identifier for a GRIB message.
  Integer                       :: Status                           ! Status code from GRIB_API subroutines.
  Character(MaxCharLength)      :: ErrMessage                       ! Error message from GRIB_API subroutines.
  Integer                       :: nX                               !} Size of field encoded in GRIB message.
  Integer                       :: nY                               !}
  Integer                       :: numberOfValues                   ! Size of data array in GRIB message.
  Character(MaxCharLength)      :: StepUnits                        ! Time unit in GRIB header.
  Character(MaxCharLength)      :: StepType                         !} F/C step information in GRIB header.
  Integer                       :: StartStep                        !}
  Integer                       :: EndStep                          !}
  Integer                       :: varid                            ! NetCDF variable ID returned from nf90_inq_varid
  Integer                       :: start(3)                         ! The start and count arrays for reading data
  Integer                       :: count(3)                         ! using netCDF library routines.

  Real(Std), Parameter :: CoordTolerance = 0.001   ! tolerance (in degrees) when checking coord
                                                   ! values in pp headers

  ! $$ Saved variables should really be returned.

  Error                    = .false.
  InconsistentHeader       = .false.
  InconsistentValidityTime = .false.
  Missing                  = .false.

  ! GRIB met file.
  If (FileType .CIEq. 'GRIB') Then

#   ifdef GRIBsupport

      ! Validity time for the met field as a character string (format YYYYMMDDHHMM).
      CharTime = FileNameTime(Time) ! $$ more efficient to do this before calling this routine?

      ! Select GRIB file index to use for accessing this field.
      If (iFieldLevel == 0) Then
        ! access as single-level field (i.e. surface or standard level)
        indexID => indexSL
      Else
        ! access as model-level field
        indexID => indexML
      End If
      
      ! Select GRIB message in the indexed met file that contains the required field.
      !
      ! Use the GRIB_API subroutine 'grib_index_select(indexID, key, value, status)'.
      !   Action: selects the message subset with key==value in the index.
      !   Arguments:
      !     indexID  - unique identifier of an index created from a GRIB met file
      !     key      - key to be selected
      !     value    - value of the key to select
      !     status   - return code: GRIB_SUCCESS on a successful execution, or other integer value
      !                               on error (use 'grib_get_error_string' for error message).

      ! - select on date
      Call grib_index_select(indexID, 'validityDate', Char2Int(CharTime(1:8)), Status)
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                                &
               'ERROR: Unable to select messages in the GRIB met file on DATE = ' // &
               CharTime(1:8)                                                      // &
               ' (GRIB_API returns with the error message "'                      // &
               Trim(ErrMessage)                                                   // &
               '")',                                                                 &
               2                                                                     &
             )
        Go To 9
      End If

      ! - select on time
      Call grib_index_select(indexID, 'validityTime', Char2Int(CharTime(9:12)), Status)
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                                &
               'ERROR: Unable to select messages in the GRIB met file on TIME = ' // &
               CharTime(9:12)                                                     // &
               ' (GRIB_API returns with the error message "'                      // &
               Trim(ErrMessage)                                                   // &
               '")',                                                                 &
               2                                                                     &
             )
        Go To 9
      End If

      ! - select on parameter
      Call grib_index_select(indexID, 'paramId', FieldCode, Status)
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                                     &
               'ERROR: Unable to select messages in the GRIB met file on PARAMETER = ' // &
               Trim(Int2Char(FieldCode))                                               // &
               ' (GRIB_API returns with the error message "'                           // &
               Trim(ErrMessage)                                                        // &
               '")',                                                                      &
               2                                                                          &
             )
        Go To 9
      End If

      ! - select on level (when using model-level index)
      If (Associated(indexID, Target = indexML)) Then
        Call grib_index_select(indexID, 'level', iFieldLevel, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                                 &
                 'ERROR: Unable to select messages in the GRIB met file on LEVEL = ' // &
                 Trim(Int2Char(iFieldLevel))                                         // &
                 ' (GRIB_API returns with the error message "'                       // &
                 Trim(ErrMessage)                                                    // &
                 '")',                                                                  &
                 2                                                                      &
               )
          Go To 9
        End If
      End If

      ! - select on ensemble number (when the met case number is being used in the index)
      !   note that control forecast = 0
      If (iMetCase /= -1) Then
        Call grib_index_select(indexID, 'perturbationNumber', iMetCase, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                                  &
                 'ERROR: Unable to select messages in the GRIB met file on NUMBER = ' // &
                 Trim(Int2Char(iMetCase))                                             // &
                 ' (GRIB_API returns with the error message "'                        // &
                 Trim(ErrMessage)                                                     // &
                 '")',                                                                   &
                 2                                                                       &
               )
          Go To 9
        End If
      End If

      ! Read the selected GRIB message into memory.
      !
      ! Use the GRIB_API subroutine 'grib_new_from_index(indexID, gribID, status)'.
      !   Action: loads a GRIB message selected from an index into memory.
      !           note that all the index key values must be selected.
      !   Arguments:
      !     indexID  - unique identifier of an index created from a GRIB met file
      !     gribID   - unique identifier of GRIB message for use by subsequent GRIB_API calls
      !     status   - return code: GRIB_SUCCESS on a successful execution, GRIB_END_OF_INDEX
      !                               when no further matching messages are present in
      !                               the index, or other integer value on error
      !                               (use 'grib_get_error_string' for error message).
      Call grib_new_from_index(indexID, gribID, Status)
      If (Status /= GRIB_SUCCESS .and. Status /= GRIB_END_OF_INDEX) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                          &
               'ERROR: problem reading NWP met field with parameter ID = '  // &
               Trim(Int2Char(FieldCode))                                    // &
               ' on met model level '                                       // &
               Trim(Int2Char(iFieldLevel))                                  // &
               ' (GRIB_API returns with the error message "'                // &
               Trim(ErrMessage)                                             // &
               '")',                                                           &
               2                                                               &
             )
        Go To 9
      End If

      ! GRIB product missing from met file.
      If (Status == GRIB_END_OF_INDEX) Then

        ! Currently fixing up a missing field if its field code is '-999999' (sets all field values to zero)
        ! $$ More specific fix up for missing met fields?
        If (FieldCode == -999999) Then
          Field(:, :) = 0.0
        Else
          Missing = .true.
        End If

      ! GRIB product present in met file.
      Else

        ! Check that met field has the expected grid size.
        !
        ! Use the GRIB_API subroutine 'grib_get(gribID, key, value, status)'.
        !   Action: gets the value of a key from a GRIB message loaded in memory.
        !   Arguments:
        !     gribID   - unique identifier of GRIB message
        !     key      - name of key
        !     value    - value of the key
        !     status   - return code: GRIB_SUCCESS on a successful execution, or other integer value
        !                               on error (use 'grib_get_error_string' for error message).
        Call grib_get(gribID, 'Ni', nX, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                  &
                 'ERROR: unable to interpret grid size in met field ' // &
                 '(GRIB_API returns with the error message "'         // &
                 Trim(ErrMessage)                                     // &
                 '")',                                                   &
                 3                                                       &
               )
          Go To 9
        End If
        Call grib_get(gribID, 'Nj', nY, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                  &
                 'ERROR: unable to interpret grid size in met field ' // &
                 '(GRIB_API returns with the error message "'         // &
                 Trim(ErrMessage)                                     // &
                 '")',                                                   &
                 3                                                       &
               )
          Go To 9
        End If
        
        If (nX /= HGrid%nX .or. nY /= HGrid%nY) Then
          Call Message(                                                                           &
                 'ERROR: grid size of the met field is inconsistent with the NWP met definition', &
                 3                                                                                &
               )
          Go To 9
        End If

        ! Check size of the data values array in the GRIB message
        !
        ! Use the GRIB_API subroutine 'grib_get_size(gribID, key, size, status)'.
        !   Action: gets the size of an array key from a GRIB message loaded in memory.
        !   Arguments:
        !     gribID   - unique identifier of GRIB message
        !     key      - name of array key
        !     size     - size of the array key
        !     status   - return code: GRIB_SUCCESS on a successful execution, or other integer value
        !                               on error (use 'grib_get_error_string' for error message).
        Call grib_get_size(gribID, 'values', numberOfValues, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                                  &
                 'ERROR: unable to interpret size of the data array in GRIB message ' // &
                 '(GRIB_API returns with the error message "'                         // &
                 Trim(ErrMessage)                                                     // &
                 '")',                                                                   &
                 3                                                                       &
               )
          Go To 9
        End If
        
        If (numberOfValues /= HGrid%nX * HGrid%nY) Then
          Call Message(                                                            &
                 'ERROR: incorrect size of data values array in the GRIB message', &
                 3                                                                 &
               )
          Go To 9
        End If

        ! Get data values of met field.
        Call grib_get(gribID, 'values', LocalField2, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                            &
                 'ERROR: unable to read data values array in the GRIB message ' // &
                 '(GRIB_API returns with the error message "'                   // &
                 Trim(ErrMessage)                                               // &
                 '")',                                                             &
                 3                                                                 &
               )
          Go To 9
        End If

        ! Put data into 2-d met field array.
        k = 0
        Do j = 1, HGrid%nY
        Do i = 1, HGrid%nX
          k = k + 1
          Field(i, j) = LocalField2(k)
        End Do
        End Do

        ! Check time unit is 'hours' in the GRIB header.
        Call grib_get(gribID, 'stepUnits', StepUnits, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                          &
                 'ERROR: unable to interpret step time unit in GRIB message ' // &
                 '(GRIB_API returns with the error message "'                 // &
                 Trim(ErrMessage)                                             // &
                 '")',                                                           &
                 3                                                               &
               )
          Go To 9
        End If
        If (StepUnits /= 'h') Then
          Call Message(                                                   &
                 'ERROR: time unit of the forecast step is not in hours', &
                 3                                                        &
               )
          InconsistentHeader = .true.
        End If

        ! Check forecast step in the GRIB header is consistent with previous fields read for this met time.
        Call grib_get(gribID, 'stepType', StepType, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                     &
                 'ERROR: unable to interpret step type in GRIB message ' // &
                 '(GRIB_API returns with the error message "'            // &
                 Trim(ErrMessage)                                        // &
                 '")',                                                      &
                 3                                                          &
               )
          Go To 9
        End If
        Call grib_get(gribID, 'startStep', StartStep, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                      &
                 'ERROR: unable to interpret start step in GRIB message ' // &
                 '(GRIB_API returns with the error message "'             // &
                 Trim(ErrMessage)                                         // &
                 '")',                                                       &
                 3                                                           &
               )
          Go To 9
        End If
        Call grib_get(gribID, 'endStep', EndStep, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                    &
                 'ERROR: unable to interpret end step in GRIB message ' // &
                 '(GRIB_API returns with the error message "'           // &
                 Trim(ErrMessage)                                       // &
                 '")',                                                     &
                 3                                                         &
               )
          Go To 9
        End If
        If (StepType == 'instant') Then
          ! Instantaneous field valid at StartStep
          ForecastStepL = StartStep
        Else If (StepType == 'avg') Then
          ! Time-averaged field calculated between StartStep and EndStep
          ForecastStepL = EndStep
        Else
          ! Unsupported field type
          ForecastStepL = -1
          InconsistentHeader = .true.
        End If

        If (ForecastStep == -1) Then
          ! Set forecast step using the first valid field read at this met time
          ForecastStep = ForecastStepL
        Else If (ForecastStep /= ForecastStepL) Then
          ! Latest field has an inconsistent forecast step
          InconsistentHeader = .true.
        End If

        ! $$ More checking of headers, etc. for consistency with what is expected.

        ! Release the GRIB message from memory.
        !
        ! Use the GRIB_API subroutine 'grib_release(gribID, status)'.
        !   Action: free the memory of the GRIB message specified by gribID.
        !   Arguments:
        !     gribID   - unique identifier of GRIB message to be unloaded
        !     status   - return code: GRIB_SUCCESS on a successful execution, or other integer value
        !                               on error (use 'grib_get_error_string' for error message).
        Call grib_release(gribID, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                     &
                 'ERROR: unable to release memory for the GRIB message ' // &
                 '(GRIB_API returns with the error message "'            // &
                 Trim(ErrMessage)                                        // &
                 '")',                                                      &
                 3                                                          &
               )
          Go To 9
        End If

        ! Check that there are no further GRIB messages with the specified key selection.
        Call grib_new_from_index(indexID, gribID, Status)
        If (Status /= GRIB_END_OF_INDEX) Then
          Call Message(                                                          &
                 'WARNING: more than one NWP met field with parameter ID = '  // &
                 Trim(Int2Char(FieldCode))                                    // &
                 ' on met model level '                                       // &
                 Trim(Int2Char(iFieldLevel))                                  // &
                 ' has been found in the met file. '                          // &
                 'Only the first occurrence will be used.',                      &
                 2                                                               &
               )
        End If
        Call grib_release(gribID, Status)
        If (Status /= GRIB_SUCCESS) Then
          Call grib_get_error_string(Status, ErrMessage)
          Call Message(                                                     &
                 'ERROR: unable to release memory for the GRIB message ' // &
                 '(GRIB_API returns with the error message "'            // &
                 Trim(ErrMessage)                                        // &
                 '")',                                                      &
                 3                                                          &
               )
          Go To 9
        End If
      End If

#   endif

  ! NetCDF met file.
  Else If (FileType .CIEq. 'NetCDF') Then

#   ifdef NetCDFsupport

    If (Trim(NCFieldName) .CIEq. '') Then
      ! Ignore any variables which do not have an "NC Field Name" defined.
      ! $$ Could possibly trap this as an error instead here?
    Else
      ! Get variable ID (varid) using the "NC Field Name" attribute.
      Status = nf90_inq_varid(Unit, Trim(NCFieldName), varid)
      If (Status /= nf90_noerr) Then
        Call Message(                                                                 &
                 'ERROR: NetCDF_API returns with the error message "'              // &
                 Trim(nf90_strerror(Status))                                       // &
                 '"',                                                                 &
                 2                                                                    &
               )
        Go To 9
      End If

      ! Set indices specifying the relative location and dimensions of the slice in
      ! the NetCDF variable. 3-d variables are read in as a series of 2-d fields -
      ! lat-lon for a specific level - just like pp fields.
      count = (/ HGrid%nX, HGrid%nY, 1 /)
      start = (/ 1, 1, NCLevelIndex /)

      ! Read 2-d variable from NetCDF file into temporary variable LocalField.
      Status = nf90_get_var(Unit, varid, LocalField, start, count)
      If (Status /= nf90_noerr) Then
        Call Message(                                                                 &
                 'ERROR: NetCDF_API returns with the error message "'              // &
                 Trim(nf90_strerror(Status))                                       // &
                 '"',                                                                 &
                 2                                                                    &
               )
        Go To 9
      End If

      Field = LocalField

    End If

#   endif

  ! PP met file.
  Else If (FileType .CIEq. 'PP') Then

    If (NewFile) MissingL = .false.

    ! Read header.
    If (.not.MissingL) Then
      Read (Unit = Unit, Err = 9, IOStat = IOStat) IHeader, RHeader
    End If

    ! Field missing.
    ! Notes:
    ! IHeader(42) = Code for the field (stash code).
    ! IHeader(33) = Level number.
    If (IHeader(42) /= FieldCode .or. (ThreeD .and. IHeader(33) /= iFieldLevel) .or. IOStat < 0) Then

      MissingL = .true.
      Missing  = .true.

    ! Field not missing.
    Else

      MissingL = .false.

      ! Check header. Note code for the field (stash code) and level number already checked when testing for
      ! missing fields above.
      ! Notes:
      ! IHeader( 1) = } Year, month, day, hour, and minute of time 1. Generally time 1 is the time the data is
      ! IHeader( 2) = } valid for. Here we assume this is the case, but see PP file documentation for other
      ! IHeader( 3) = } possibilities.
      ! IHeader( 4) = }
      ! IHeader( 5) = }
      ! IHeader( 7) = ] Year, month, day, hour, and minute of time 2. Generally time 2 is the time of the data
      ! IHeader( 8) = ] used to calculte the field (e.g. the analysis time from which a forecast is produced).
      ! IHeader( 9) = ] Here we assume this is the case, but see PP file documentation for other
      ! IHeader(10) = ] possibilities. $$ actually not yet checked at all.
      ! IHeader(11) = ]
      ! IHeader(14) = Forecast time in hours.
      ! IHeader(18) = } Number of y and x points in the grid for the field.
      ! IHeader(19) = }
      ! RHeader(11) = } Latitude and longitude of coord system pole.
      ! RHeader(12) = }
      ! RHeader(14) = Zero-th y point of main grid (i.e. first y point - dy).
      ! RHeader(15) = |dy| (is | | right? $$)
      ! RHeader(16) = Zero-th x point of main grid (i.e. first x point - dx).
      ! RHeader(17) = |dx| (is | | right? $$)

      ! Validity time.
      If (                              & ! $$ would be better to maintain private encapsulation of
                                          ! time structure
        IHeader(1) /= Time%Year    .or. & ! - but not immediately clear how best to do this. Time zones?
        IHeader(2) /= Time%Month   .or. &
        IHeader(3) /= Time%Day     .or. &
        IHeader(4) /= Time%Hour    .or. &
        IHeader(5) /= Time%Minute       &
      ) Then
        InconsistentValidityTime = .true.
      End If

      ! Check forecast step is consistent with previous fields read for this met time.
      ForecastStepL = IHeader(14)
      If (ForecastStep == -1) Then
        ! Set forecast step using the first valid field read at this met time.
        ForecastStep = ForecastStepL
      Else If (ForecastStep /= ForecastStepL) Then
        ! Latest field has an inconsistent forecast step.
        InconsistentHeader = .true.
      End If

      ! Check grid size nX, nY and grid spacing dX, dY.
      ! $$ grid origin currently not checked to avoid potential ambiguities e.g. with multiples of 2 pi
      InconsistentHeader = InconsistentHeader      .or. &
                           HGrid%nX /= IHeader(19) .or. &
                           HGrid%nY /= IHeader(18)
      InconsistentHeader = InconsistentHeader                                                           .or. &
                           Abs((HGrid%dX * HCoord%Unit(1) * 180.0 / Pi) - RHeader(17)) > CoordTolerance .or. &
                           Abs((HGrid%dY * HCoord%Unit(2) * 180.0 / Pi) - RHeader(15)) > CoordTolerance

      ! $$ More checks of headers in PP met files?

      ! Read field.

      If (IHeader(21) == 1) Then   ! Process WGDOS packed PP data

        ! Allocate temporary array of 64-bit integers, PackedField, with space to store
        ! (LBLREC - LBEXT) 32-bit words. This explains the division by 2 below. Note that
        ! LBLREC is the total length of the data record, which includes both the data values
        ! of the field itself and possibly 'extra data' (though PP files for NAME do not
        ! generally contain any 'extra data'), LBEXT is the length of this extra data, and
        ! so the difference gives the length of the field array.

        lengthOfPackedData = (IHeader(15) - IHeader(20)) / 2

        Allocate(PackedField(lengthOfPackedData), Stat = ErrorCode)

        If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate array for PackedField', 3)

        Read (Unit = Unit, Err = 9, End = 9) PackedField

        Status=0

        ! Get the packing accuracy of the WGDOS data from the real header
        PackingAccuracy = Int(RHeader(6))

        ! Unpack PackedField data array into LocalField dimensioned by no. of points in
        ! each row IHeader(19) by the no. of rows in a field IHeader(18).
        Call Xpnd_all(LocalField, PackedField, lengthOfPackedData, IHeader(19), IHeader(18), &
                      PackingAccuracy, RHeader(18), Status, ErrMessage)

        Deallocate(PackedField)

         If (Status /= 0) Then
           Call Message(                                                                      &
                 'ERROR:  Xpnd_all (WGDOS unpacking) returns with the error message "'     // &
                 Trim(ErrMessage)                                                          // &
                 '" for stash code '                                                       // &
                 Trim(Int2Char(IHeader(42))),                                                 &
                 3                                                                            &
               )
           Go To 9
         End If

      Else

        Read (Unit = Unit, Err = 9, End = 9) LocalField

      End If

      Field = LocalField


    End If

  ! Name II met file.
  Else If (FileType .CIEq. 'Name II') Then

    ! Read header.
    Read (Unit = Unit, Err = 9, End = 9) Header

    ! Field missing.
    If (Int(Header(1)) <= 0) Then

      Missing = .true.

    ! Field not missing.
    Else

      ! Check header.
      ! Notes:
      ! In the header the horizontal data grid can be defined indirectly as a rectangular subset of another
      ! grid (called here the main grid).
      ! Header( 1) = } Year, month, day, hour, and minute of time 1. Header(1) <= 0 indicates a missing field.
      ! Header( 2) = } Generally time 1 is the time the data is valid for. Here we assume this is the case,
      ! Header( 3) = } but see PP file documentation for other possibilities.
      ! Header( 4) = }
      ! Header( 5) = }
      ! Header( 7) = ] Year, month, day, hour, and minute of time 2. Generally time 2 is the time of the data
      ! Header( 8) = ] used to calculte the field (e.g. the analysis time from which a forecast is produced).
      ! Header( 9) = ] Here we assume this is the case, but see PP file documentation for other possibilities.
      ! Header(10) = ] $$ actually not yet checked at all.
      ! Header(11) = ]
      ! Header(14) = Forecast time in hours.
      ! Header(42) = Code for the field (stash code).
      ! Header(33) = Level number.
      ! Header(18) = } Number of y and x points in the main grid.
      ! Header(19) = }
      ! Header(34) = ] Indices in the main grid of the first and last x points and the first and last y points
      ! Header(35) = ] of the 'rectangular subset' grid. Values <= 0 indicate the main grid is used. If first
      ! Header(36) = ] index > last index the grid wraps. First index = last index isn't allowed.
      ! Header(37) = ]
      ! Header(56) = } Latitude and longitude of coord system pole.
      ! Header(57) = }
      ! Header(59) = First y point of main grid - dy.
      ! Header(60) = |dy| (is ! ! right? $$)
      ! Header(61) = First x point of main grid - dx.
      ! Header(62) = |dx| (is ! ! right? $$)
      ! $$ make more complete checks (in particular sign dx, dy)?
      ! $$ check header defns assumed are correct.
      ! $$ these checks (at present) assume that grid is lat-long (degrees) and Header values are always in
      !    lat-long degrees.
      ! $$ Try to encapsulate some of the coord changes in the coord module?

      ! Note we don't check stash and level number for missing fields since these are not always correct (at
      ! least level isn't).

      ! Validity time.
      If (                             &
        Header(1) /= Time%Year    .or. &
        Header(2) /= Time%Month   .or. &
        Header(3) /= Time%Day     .or. &
        Header(4) /= Time%Hour    .or. &
        Header(5) /= Time%Minute       &
      ) Then
        InconsistentValidityTime = .true.
      End If

      ! Check forecast step is consistent with previous fields read for this met time.
      ForecastStepL = Header(14)
      If (ForecastStep == -1) Then
        ! Set forecast step using the first valid field read at this met time
        ForecastStep = ForecastStepL
      Else If (ForecastStep /= ForecastStepL) Then
        ! Latest field has an inconsistent forecast step
        InconsistentHeader = .true.
      End If

      ! Stash code and level number.
      If (NInt(Header(42)) /= FieldCode .or. (ThreeD .and. NInt(Header(33)) /= iFieldLevel)) Then
        ! Fix up to suppress errors for a known problem in our met files. This involves H regional met files
        ! (including cut down versions of these files) between 1/1/99 and 29/3/99 inclusive and after
        ! 1/1/2001. (Note this includes H regional files after 1/1/2001 but excludes H2001 regional files
        ! which cover the same period.) This avoids a lot of error messages for a known problem with these
        ! files. These files have, on the main vertical grid, NWP level 18 not 19. (Note level 18 also appears
        ! in H2001 global and regional files but, being consistent, this has been treated by altering the
        ! levels in the met definitions to 18.) $$
        i = Scan(FileName, '.', .true.)
        i = Min(i + 1, Len(FileName))
        j = Min(i + 3, Len(FileName))
        If (                                                                                            &
          NInt(Header(33)) == 18                                                                  .and. &
          iFieldLevel      == 19                                                                  .and. &
          FileName(i : j) == 'REGH'                                                               .and. &
          FileName(Min(j + 1, Len(FileName)) : Min(j + 1, Len(FileName))) /= '2'                  .and. &
          (                                                                                             &
            (NInt(Header(1)) == 1999 .and. NInt(Header(2)) <= 2)                             .or.       &
            (NInt(Header(1)) == 1999 .and. NInt(Header(2)) == 3 .and. NInt(Header(3)) <= 29) .or.       &
            (NInt(Header(1)) >= 2001)                                                                   &
          )                                                                                             &
        ) Then
        Else
          InconsistentHeader = .true.
        End If
      End If

      ! Headers 34-37.
      If (Header(34) <= 0 .and. Header(35) <= 0 .and. Header(36) <= 0 .and. Header(37) <= 0) Then

        Header(34) = 1.0
        Header(35) = Header(19)
        Header(36) = 1.0
        Header(37) = Header(18)

      Else If (                                                                             &
        Header(34) > 0 .and. Header(35) > 0 .and. Header(36) > 0 .and. Header(37) > 0 .and. &
        Header(34) <= Header(19) .and. Header(35) <= Header(19)                       .and. &
        Header(36) <= Header(18) .and. Header(37) <= Header(18)                             &
      ) Then

        If (Header(34) > Header(35)) Then
          ! $$ check 'main grid' joins up
          Header(35) = Header(35) + Header(19)
        Else If (Header(34) == Header(35)) Then
          InconsistentHeader = .true.
        End If

        If (Header(36) > Header(37)) Then
          ! $$ check 'main grid' joins up
          Header(37) = Header(37) + Header(18)
        Else If (Header(36) == Header(37)) Then
          InconsistentHeader = .true.
        End If

        ! nX, nY.
        InconsistentHeader = InconsistentHeader                                .or. &
                             HGrid%nX /= Int(Header(35)) - Int(Header(34)) + 1 .or. &
                             HGrid%nY /= Int(Header(37)) - Int(Header(36)) + 1

      Else

        InconsistentHeader = .true.

      End If

      ! nX, nY.
   !   InconsistentHeader = InconsistentHeader                                .or. &
   !                        HGrid%nX /= Int(Header(35)) - Int(Header(34)) + 1 .or. &
   !                        HGrid%nY /= Int(Header(37)) - Int(Header(36)) + 1

      ! dX.
   !   Scale              = 4.0 * Abs(HGrid%dX * HCoord%Unit(1) * 180.0 / Pi) + &
   !                        1.0 * Header(62)
   !   Scale              = Scale * 10.0 * Epsilon(1.0_Std)
   !   InconsistentHeader = InconsistentHeader .or.                                               &
   !                        Abs(Abs(HGrid%dX * HCoord%Unit(1) * 180.0 / Pi) - Header(62)) > Scale

      ! dY.
   !   Scale              = 4.0 * Abs(HGrid%dY * HCoord%Unit(2) * 180.0 / Pi) + &
   !                        1.0 * Header(60)
   !   Scale              = Scale * 10.0 * Epsilon(1.0_Std)
   !   InconsistentHeader = InconsistentHeader .or.                                               &
   !                        Abs(Abs(HGrid%dY * HCoord%Unit(2) * 180.0 / Pi) - Header(60)) > Scale

      ! X0.
   !   Scale              = 4.0 * Abs(             HGrid%X0 * HCoord%Unit(1) * 180.0 / Pi) + &
   !                        4.0 * Abs(Header(34) * HGrid%dX * HCoord%Unit(1) * 180.0 / Pi) + &
   !                        3.0 * Abs(HCoord%Origin(1)                       * 180.0 / Pi) + &
   !                        1.0 * Abs(Header(61))
   !   Scale              = Scale * 10.0 * Epsilon(1.0_Std)
   !   InconsistentHeader = InconsistentHeader .or.                                                    &
   !                        Abs(                                                                       &
   !                          ((HGrid%X0 - Header(34) * HGrid%dX) * HCoord%Unit(1) + HCoord%Origin(1)) &
   !                          * 180.0 / Pi - Header(61)                                                &
   !                        ) > Scale

      ! Y0.
   !   Scale              = 4.0 * Abs(             HGrid%Y0 * HCoord%Unit(2) * 180.0 / Pi) + &
   !                        4.0 * Abs(Header(36) * HGrid%dY * HCoord%Unit(2) * 180.0 / Pi) + &
   !                        3.0 * Abs(HCoord%Origin(2)                       * 180.0 / Pi) + &
   !                        1.0 * Abs(Header(59))
   !   Scale              = Scale * 10.0 * Epsilon(1.0_Std)
   !   InconsistentHeader = InconsistentHeader .or.                                                    &
   !                        Abs(                                                                       &
   !                          ((HGrid%Y0 - Header(36) * HGrid%dY) * HCoord%Unit(2) + HCoord%Origin(2)) &
   !                          * 180.0 / Pi - Header(59)                                                &
   !                        ) > Scale

      ! Read field.
      Do j = 1, HGrid%nY
        Read (Unit = Unit, Err = 9, End = 9) (Field(i, j), i = 1, HGrid%nX)
      End Do

    End If

  Else

    Call Message('UNEXPECTED FATAL ERROR in Read2dField: invalid met file type', 4)

  End If

  Return

  ! Read error.
9 Continue

  Error = .true.

End Subroutine Read2dField

!-------------------------------------------------------------------------------------------------------------

Subroutine ProcessNWPMet(Coords, Grids, M, OldIdx, NewIdx)
! Processes the met data.

  Implicit None
  ! Argument List:
  Type(Coords_),  Intent(In),   Target :: Coords ! Collection of coord systems.
  Type(Grids_),   Intent(In),   Target :: Grids  ! Collection of grids.
  Type(NWPMet_),  Intent(InOut)        :: M      ! State of an NWP met module instance.
  ! Locals:
  Type(ZCoord_), Pointer :: ZCoord                !} Abbreviations for coords and
  Type(ZCoord_), Pointer :: ZCoordZ               !} grids.
  Type(ZCoord_), Pointer :: ZCoordP               !}
  Type(ZCoord_), Pointer :: ZCoordCl              !}
  Type(HGrid_),  Pointer :: HGrid                 !}
  Type(HGrid_),  Pointer :: HGridU                !}
  Type(HGrid_),  Pointer :: HGridV                !}
  Type(ZGrid_),  Pointer :: ZGrid                 !}
  Type(ZGrid_),  Pointer :: ZGridUV               !}
  Type(ZGrid_),  Pointer :: ZGridW                !}
  Type(ZGrid_),  Pointer :: ZGridP                !}
  Real(Std)              :: UStress               !] U and V components of the surface
  Real(Std)              :: VStress               !] stress interpolated to the main
                                                  !] grid.
  Integer                :: i                     !} Loop indices.
  Integer                :: j                     !}
  Integer                :: k                     !}
  Real(Std)              :: PLowestLevelByPGround ! Estimate of P(lowest level present) / P(ground).
  Real(Std)              :: PGroundByPSea         ! Estimate of P(ground) / P(sea).
  Integer                :: LowestLevel           ! Index of lowest level for which P is present.
  Integer, Intent(In),Optional      :: OldIdx, NewIdx        ! Index of the data to be
                                                  ! processed. This can be the
                                                  ! index of the Prefetch buffer if parallel IO is
                                                  ! used.
  Real(Std)              :: MTH                   !} Extra variables for use in vertical
  Real(Std)              :: Htemp                 !} windspeed mass consistency correction.
  Real(Std)              :: H1                    !} Currently not in use
  Real(Std)              :: H2                    !}
  Real(Std)              :: HMax                  !}
  Real(Pos)              :: Point(2)              !} 
  Integer                :: MOld, MNew            
  Real(Std)              :: MinTempDiff           ! Temporary variables for calculation of pressure level nearer to freezing point 
  Real(Std)              :: TempDiffToFL          ! Temporary variables for calculation of pressure level nearer to freezing point 
  Integer                :: fl                    ! Temporary variables for calculation of pressure level nearer to freezing point 

  If (Present(NewIdx)) Then
    MNew = NewIdx
  Else
    MNew = M%New
  End If

  If (Present(OldIdx)) Then
    MOld = OldIdx
  Else
    MOld = M%Old
  End If

  ! 1) Set up abbreviations for coords and grids.

  ZCoord   => Coords%ZCoords(M%iZCoord  )
  ZCoordZ  => Coords%ZCoords(M%iZCoordZ )
  ZCoordP  => Coords%ZCoords(M%iZCoordP )
  ZCoordCl => Coords%ZCoords(M%iZCoordCl)
  HGrid    => Grids%HGrids(M%iHGrid  )
  HGridU   => Grids%HGrids(M%iHGridU )
  HGridV   => Grids%HGrids(M%iHGridV )
  ZGrid    => Grids%ZGrids(M%iZGrid  )
  ZGridUV  => Grids%ZGrids(M%iZGridUV)
  ZGridW   => Grids%ZGrids(M%iZGridW )
  ZGridP   => Grids%ZGrids(M%iZGridP )

  ! 2) Set missing values and update FieldPersist3d/2d (except for P and H which are treated below). It is
  ! assumed that, of the 2-d fields, at most roughness length, boundary layer depth and sea level pressure
  ! need fixing. Note surface velocities are set to zero even if they are not missing.

  ! $$ Allow more 2-d fields to be missing?
  ! $$ could use better extrapolation/interpolation - e.g. log, use z values not indices.
  ! $$ w persistence not good as the persisted value will go through the coord conversion again! (but w
  !    persistence not currently switched on). Either fix or decide w persistence not allowed.
  ! $$ default fix-ups not ideal for T, Q.
  ! $$ general consistency issues of fixed up data - e.g. cloud water/ice fixup could interact badly with
  !    H/M/L cloud.
  ! $$ Q, CloudWater and CloudIce fix ups are (more or less) correct for settings of ExtrapolationLevelsBelow
  !    etc specified above. However the code needs generalising.
  ! $$ Treatment of P is not quite consistent with comments given at start of module (with defns of
  !    ExtrapolationLevelsBelow etc. Should change this. Its better to use a particular p level or sea level
  !    p consistently across whole domain to avoid discontinuities.
  ! $$ Treatment of bl depth doesn't include possibility of persistence.
  ! $$ Worthwhile altering the fix-up blocks below so we can loop over the various fields (as in ReadNWPMet)?

  ! Roughness length: default fix up depends on topography. It is 0.1m for land (topog > 0) and 0.0001m for
  ! sea (topog = 0).


  If (M%FieldPresent2d(F_Z0)) Then

    M%FieldPersist2d(F_Z0) = 0

  Else If (                                             &
    M%FieldPersist2d(F_Z0) < PersistTimes2d(F_Z0) .and. &
    M%FieldPersist2d(F_Z0) >= 0                   .and. &
    MOld /= 0                                          &
  ) Then

    M%Z0(:, :, MNew) = M%Z0(:, :, MOld)

    M%FieldPersist2d(F_Z0) = M%FieldPersist2d(F_Z0) + 1

  Else If (AllowFixUpToDefault2d(F_Z0)) Then

    Do j = 1, HGrid%nY
    Do i = 1, HGrid%nX
      If (M%Topog(i, j, 1) > 0) Then
         M%Z0(i, j, MNew) = 0.1
      Else
         M%Z0(i, j, MNew) = 0.0001
      End If
    End Do
    End Do

    M%FieldPersist2d(F_Z0) = -1

  End If

  ! Sea level pressure: default fix up is PAt0km.
  If (M%FieldPresent2d(F_PSeaLevel)) Then

    M%FieldPersist2d(F_PSeaLevel) = 0

  Else If (                                                           &
    M%FieldPersist2d(F_PSeaLevel) < PersistTimes2d(F_PSeaLevel) .and. &
    M%FieldPersist2d(F_PSeaLevel) >= 0                          .and. &
    MOld /= 0                                                        &
  ) Then

    M%PSeaLevel(:, :, MNew) = M%PSeaLevel(:, :, MOld)

    M%FieldPersist2d(F_PSeaLevel) = M%FieldPersist2d(F_PSeaLevel) + 1

  Else If (AllowFixUpToDefault2d(F_PSeaLevel)) Then

    M%PSeaLevel(:, :, MNew) = PAt0km

    M%FieldPersist2d(F_PSeaLevel) = -1

  End If

  ! Convective cloud amount: default fix up is 0.
  If (M%FieldPresent2d(F_ConCloud)) Then

    M%FieldPersist2d(F_ConCloud) = 0

  Else If (                                                         &
    M%FieldPersist2d(F_ConCloud) < PersistTimes2d(F_ConCloud) .and. &
    M%FieldPersist2d(F_ConCloud) >= 0                         .and. &
    MOld /= 0                                                      &
  ) Then

    M%ConCloud(:, :, MNew) = M%ConCloud(:, :, MOld)

    M%FieldPersist2d(F_ConCloud) = M%FieldPersist2d(F_ConCloud) + 1

  Else If (AllowFixUpToDefault2d(F_ConCloud)) Then

    M%ConCloud(:, :, MNew) = 0.0

    M%FieldPersist2d(F_ConCloud) = -1

  End If

  ! Convective cloud base: default fix up is -Huge.
  If (M%FieldPresent2d(F_ConCloudBase)) Then

    M%FieldPersist2d(F_ConCloudBase) = 0

  Else If (                                                                 &
    M%FieldPersist2d(F_ConCloudBase) < PersistTimes2d(F_ConCloudBase) .and. &
    M%FieldPersist2d(F_ConCloudBase) >= 0                             .and. &
    MOld /= 0                                                              &
  ) Then

    M%ConCloudBase(:, :, MNew) = M%ConCloudBase(:, :, MOld)

    M%FieldPersist2d(F_ConCloudBase) = M%FieldPersist2d(F_ConCloudBase) + 1

  Else If (AllowFixUpToDefault2d(F_ConCloudBase)) Then

    M%ConCloudBase(:, :, MNew) = -Huge(1.0) / 3.0

    M%FieldPersist2d(F_ConCloudBase) = -1

  End If

  ! Convective cloud top: default fix up is -Huge.
  If (M%FieldPresent2d(F_ConCloudTop)) Then

    M%FieldPersist2d(F_ConCloudTop) = 0

  Else If (                                                               &
    M%FieldPersist2d(F_ConCloudTop) < PersistTimes2d(F_ConCloudTop) .and. &
    M%FieldPersist2d(F_ConCloudTop) >= 0                            .and. &
    MOld /= 0                                                            &
  ) Then

    M%ConCloudTop(:, :, MNew) = M%ConCloudTop(:, :, MOld)

    M%FieldPersist2d(F_ConCloudTop) = M%FieldPersist2d(F_ConCloudTop) + 1

  Else If (AllowFixUpToDefault2d(F_ConCloudTop)) Then

    M%ConCloudTop(:, :, MNew) = -Huge(1.0) / 3.0

    M%FieldPersist2d(F_ConCloudTop) = -1

  End If

  ! Convective rain rate: default fix up is 0.
  If (M%FieldPresent2d(F_ConRain)) Then

    M%FieldPersist2d(F_ConRain) = 0

  Else If (                                                       &
    M%FieldPersist2d(F_ConRain) < PersistTimes2d(F_ConRain) .and. &
    M%FieldPersist2d(F_ConRain) >= 0                        .and. &
    MOld /= 0                                                    &
  ) Then

    M%ConRain(:, :, MNew) = M%ConRain(:, :, MOld)

    M%FieldPersist2d(F_ConRain) = M%FieldPersist2d(F_ConRain) + 1

  Else If (AllowFixUpToDefault2d(F_ConRain)) Then

    M%ConRain(:, :, MNew) = 0.0

    M%FieldPersist2d(F_ConRain) = -1

  End If

  ! Convective snow rate: default fix up is 0.
  If (M%FieldPresent2d(F_ConSnow)) Then

    M%FieldPersist2d(F_ConSnow) = 0

  Else If (                                                       &
    M%FieldPersist2d(F_ConSnow) < PersistTimes2d(F_ConSnow) .and. &
    M%FieldPersist2d(F_ConSnow) >= 0                        .and. &
    MOld /= 0                                                    &
  ) Then

    M%ConSnow(:, :, MNew) = M%ConSnow(:, :, MOld)

    M%FieldPersist2d(F_ConSnow) = M%FieldPersist2d(F_ConSnow) + 1

  Else If (AllowFixUpToDefault2d(F_ConSnow)) Then

    M%ConSnow(:, :, MNew) = 0.0

    M%FieldPersist2d(F_ConSnow) = -1

  End If

  ! Wind (u-cpt): default fix up is zero.
  Do k = 2, ZGridUV%nZ

    If (M%FieldPresent3d(k, F_U)) Then

      M%FieldPersist3d(k, F_U) = 0

    Else If (k < M%LowestPresent(F_U) .and. k >= M%LowestPresent(F_U) - ExtrapolationLevelsBelow(F_U)) Then

      M%U(:, :, k, MNew) = (Real(k - 1) / Real(M%LowestPresent(F_U) - 1)) * &
                            M%U(:, :, M%LowestPresent(F_U), MNew)

      M%FieldPersist3d(k, F_U) = 0

    Else If (k > M%HighestPresent(F_U) .and. k <= M%HighestPresent(F_U) + ExtrapolationLevelsAbove(F_U)) Then

      M%U(:, :, k, MNew) = M%U(:, :, M%HighestPresent(F_U), MNew)

      M%FieldPersist3d(k, F_U) = 0

    Else If (                                                                   &
      k > M%LowestPresent(F_U)                                            .and. &
      k < M%HighestPresent(F_U)                                           .and. &
      M%kUpper(k, F_U) - M%kLower(k, F_U) <= InterpolationLevels(F_U) + 1       &
    ) Then

      M%U(:, :, k, MNew) = (Real(k - M%kLower(k, F_U)) / Real(M%kUpper(k, F_U) - M%kLower(k, F_U))) * &
                            M%U(:, :, M%kUpper(k, F_U), MNew) +                                       &
                            (Real(M%kUpper(k, F_U) - k) / Real(M%kUpper(k, F_U) - M%kLower(k, F_U))) * &
                            M%U(:, :, M%kLower(k, F_U), MNew)

      M%FieldPersist3d(k, F_U) = 0

    Else If (                                              &
      M%FieldPersist3d(k, F_U) < PersistTimes3d(F_U) .and. &
      M%FieldPersist3d(k, F_U) >= 0                  .and. &
      MOld /= 0                                           &
    ) Then

      M%U(:, :, k, MNew) = M%U(:, :, k, MOld)

      M%FieldPersist3d(k, F_U) = M%FieldPersist3d(k, F_U) + 1

    Else If (AllowFixUpToDefault3d(F_U)) Then

      M%U(:, :, k, MNew) = 0.0

      M%FieldPersist3d(k, F_U) = -1

    End If

  End Do

  M%U(:, :, 1, MNew) = 0.0

  ! Wind (v-cpt): default fix up is zero.
  Do k = 2, ZGridUV%nZ

    If (M%FieldPresent3d(k, F_V)) Then

      M%FieldPersist3d(k, F_V) = 0

    Else If (k < M%LowestPresent(F_V) .and. k >= M%LowestPresent(F_V) - ExtrapolationLevelsBelow(F_V)) Then

      M%V(:, :, k, MNew) = (Real(k - 1) / Real(M%LowestPresent(F_V) - 1)) * &
                            M%V(:, :, M%LowestPresent(F_V), MNew)

      M%FieldPersist3d(k, F_V) = 0

    Else If (k > M%HighestPresent(F_V) .and. k <= M%HighestPresent(F_V) + ExtrapolationLevelsAbove(F_V)) Then

      M%V(:, :, k, MNew) = M%V(:, :, M%HighestPresent(F_V), MNew)

      M%FieldPersist3d(k, F_V) = 0

    Else If (                                                                   &
      k > M%LowestPresent(F_V)                                            .and. &
      k < M%HighestPresent(F_V)                                           .and. &
      M%kUpper(k, F_V) - M%kLower(k, F_V) <= InterpolationLevels(F_V) + 1       &
    ) Then

      M%V(:, :, k, MNew) = (Real(k - M%kLower(k, F_V)) / Real(M%kUpper(k, F_V) - M%kLower(k, F_V))) * &
                            M%V(:, :, M%kUpper(k, F_V), MNew) +                                       &
                            (Real(M%kUpper(k, F_V) - k) / Real(M%kUpper(k, F_V) - M%kLower(k, F_V))) * &
                            M%V(:, :, M%kLower(k, F_V), MNew)

      M%FieldPersist3d(k, F_V) = 0

    Else If (                                              &
      M%FieldPersist3d(k, F_V) < PersistTimes3d(F_V) .and. &
      M%FieldPersist3d(k, F_V) >= 0                  .and. &
      MOld /= 0                                           &
    ) Then

      M%V(:, :, k, MNew) = M%V(:, :, k, MOld)

      M%FieldPersist3d(k, F_V) = M%FieldPersist3d(k, F_V) + 1

    Else If (AllowFixUpToDefault3d(F_V)) Then

      M%V(:, :, k, MNew) = 0.0

      M%FieldPersist3d(k, F_V) = -1

    End If

  End Do

  M%V(:, :, 1, MNew) = 0.0

  ! Wind (w-cpt): default fix up is zero.
  Do k = 2, ZGridW%nZ

    If (M%FieldPresent3d(k, F_W)) Then

      M%FieldPersist3d(k, F_W) = 0

    Else If (k < M%LowestPresent(F_W) .and. k >= M%LowestPresent(F_W) - ExtrapolationLevelsBelow(F_W)) Then

      M%W(:, :, k, MNew) = (Real(k - 1) / Real(M%LowestPresent(F_W) - 1)) * &
                            M%W(:, :, M%LowestPresent(F_W), MNew)

      M%FieldPersist3d(k, F_W) = 0

    Else If (k > M%HighestPresent(F_W) .and. k <= M%HighestPresent(F_W) + ExtrapolationLevelsAbove(F_W)) Then

      M%W(:, :, k, MNew) = M%W(:, :, M%HighestPresent(F_W), MNew)

      M%FieldPersist3d(k, F_W) = 0

    Else If (                                                                   &
      k > M%LowestPresent(F_W)                                            .and. &
      k < M%HighestPresent(F_W)                                           .and. &
      M%kUpper(k, F_W) - M%kLower(k, F_W) <= InterpolationLevels(F_W) + 1       &
    ) Then

      M%W(:, :, k, MNew) = (Real(k - M%kLower(k, F_W)) / Real(M%kUpper(k, F_W) - M%kLower(k, F_W))) * &
                            M%W(:, :, M%kUpper(k, F_W), MNew) +                                       &
                            (Real(M%kUpper(k, F_W) - k) / Real(M%kUpper(k, F_W) - M%kLower(k, F_W))) * &
                            M%W(:, :, M%kLower(k, F_W), MNew)

      M%FieldPersist3d(k, F_W) = 0

    Else If (                                              &
      M%FieldPersist3d(k, F_W) < PersistTimes3d(F_W) .and. &
      M%FieldPersist3d(k, F_W) >= 0                  .and. &
      MOld /= 0                                           &
    ) Then

      M%W(:, :, k, MNew) = M%W(:, :, k, MOld)

      M%FieldPersist3d(k, F_W) = M%FieldPersist3d(k, F_W) + 1

    Else If (AllowFixUpToDefault3d(F_W)) Then

      M%W(:, :, k, MNew) = 0.0

      M%FieldPersist3d(k, F_W) = -1

    End If

  End Do

  M%W(:, :, 1, MNew) = 0.0

  ! Temperature: default fix up is TKAtTCEq0.
  Do k = 1, ZGrid%nZ

    If (M%FieldPresent3d(k, F_T)) Then

      M%FieldPersist3d(k, F_T) = 0

    Else If (k < M%LowestPresent(F_T) .and. k >= M%LowestPresent(F_T) - ExtrapolationLevelsBelow(F_T)) Then

      M%T(:, :, k, MNew) = M%T(:, :, M%LowestPresent(F_T), MNew)

      M%FieldPersist3d(k, F_T) = 0

    Else If (k > M%HighestPresent(F_T) .and. k <= M%HighestPresent(F_T) + ExtrapolationLevelsAbove(F_T)) Then

      M%T(:, :, k, MNew) = M%T(:, :, M%HighestPresent(F_T), MNew)

      M%FieldPersist3d(k, F_T) = 0

    Else If (                                                                   &
      k > M%LowestPresent(F_T)                                            .and. &
      k < M%HighestPresent(F_T)                                           .and. &
      M%kUpper(k, F_T) - M%kLower(k, F_T) <= InterpolationLevels(F_T) + 1       &
    ) Then

      M%T(:, :, k, MNew) = (Real(k - M%kLower(k, F_T)) / Real(M%kUpper(k, F_T) - M%kLower(k, F_T))) * &
                            M%T(:, :, M%kUpper(k, F_T), MNew) +                                       &
                            (Real(M%kUpper(k, F_T) - k) / Real(M%kUpper(k, F_T) - M%kLower(k, F_T))) * &
                            M%T(:, :, M%kLower(k, F_T), MNew)

      M%FieldPersist3d(k, F_T) = 0

    Else If (                                              &
      M%FieldPersist3d(k, F_T) < PersistTimes3d(F_T) .and. &
      M%FieldPersist3d(k, F_T) >= 0                  .and. &
      MOld /= 0                                           &
    ) Then

      M%T(:, :, k, MNew) = M%T(:, :, k, MOld)

      M%FieldPersist3d(k, F_T) = M%FieldPersist3d(k, F_T) + 1

    Else If (AllowFixUpToDefault3d(F_T)) Then

      M%T(:, :, k, MNew) = TKAtTCEq0

      M%FieldPersist3d(k, F_T) = -1

    End If

  End Do

  ! Canopy water, fix up pseudo z levels 2 to 5 for global met data (1 grid box average value only)
  Do k = 2, M%kMax(F_CanopyWater)

    If (.not.(M%FieldPresent3d(k, F_CanopyWater)) .and. k > M%HighestPresent(F_CanopyWater) .and.            &
          k <= M%HighestPresent(F_CanopyWater) + ExtrapolationLevelsAbove(F_CanopyWater)) Then

      M%CanopyWater(:, :, k, MNew) = M%CanopyWater(:, :, M%HighestPresent(F_CanopyWater), MNew)

    End If

  End Do

  ! stomatal conductance, fix up pseudo z levels 2 to 5 for global met data (1 grid box average value only)
  Do k = 2, M%kMax(F_StomataConduct)

    If (.not.(M%FieldPresent3d(k, F_StomataConduct)) .and. k > M%HighestPresent(F_StomataConduct) .and.      &
          k <= M%HighestPresent(F_StomataConduct) + ExtrapolationLevelsAbove(F_StomataConduct)) Then

      M%StomataConduct(:, :, k, MNew) = M%StomataConduct(:, :, M%HighestPresent(F_StomataConduct), MNew)

    End If

  End Do


  ! Specific humidity: at level 3 and above, fix-up to Old time level or to 0; at level 2 and below, fix-up to
  ! level above.
  Do k = ZGrid%nZ, 3, -1
    If (.not.M%FieldPresent3d(k, F_Q)) Then
      If (MOld /= 0) Then
        M%Q(:, :, k, MNew) = M%Q(:, :, k, MOld)
      Else
        M%Q(:, :, k, MNew) = 0.0
      End If
    End If
  End Do
  Do k = 2, 1, -1
    If (.not.M%FieldPresent3d(k, F_Q)) Then
      M%Q(:, :, k, MNew) = M%Q(:, :, k + 1, MNew)
    End If
  End Do

  ! Total or dynamic cloud liquid water: at level 3 and above, fix-up to Old time level or to 0; 
  ! at level 2 and below, fix-up to level above.
  Do k = ZGrid%nZ, 3, -1
    If (.not.M%FieldPresent3d(k, F_TotalOrDynCloudWater)) Then
      If (MOld /= 0) Then
        M%TotalOrDynCloudWater(:, :, k, MNew) = M%TotalOrDynCloudWater(:, :, k, MOld)
      Else
        M%TotalOrDynCloudWater(:, :, k, MNew) = 0.0
      End If
    End If
  End Do
  Do k = 2, 1, -1
    If (.not.M%FieldPresent3d(k, F_TotalOrDynCloudWater)) Then
      M%TotalOrDynCloudWater(:, :, k, MNew) = M%TotalOrDynCloudWater(:, :, k + 1, MNew)
    End If
  End Do

  ! Total or dynamic cloud ice: at level 3 and above, fix-up to Old time level or to 0; 
  ! at level 2 and below, fix-up to level above.
  Do k = ZGrid%nZ, 3, -1
    If (.not.M%FieldPresent3d(k, F_TotalOrDynCloudIce)) Then
      If (MOld /= 0) Then
        M%TotalOrDynCloudIce(:, :, k, MNew) = M%TotalOrDynCloudIce(:, :, k, MOld)
      Else
        M%TotalOrDynCloudIce(:, :, k, MNew) = 0.0
      End If
    End If
  End Do
  Do k = 2, 1, -1
    If (.not.M%FieldPresent3d(k, F_TotalOrDynCloudIce)) Then
      M%TotalOrDynCloudIce(:, :, k, MNew) = M%TotalOrDynCloudIce(:, :, k + 1, MNew)
    End If
  End Do

  ! 3) Z, P, Rho, Theta.

  ! LowestLevel.
  LowestLevel = M%LowestPresent(F_PAsRead) ! $$ could use righthand side directly at expense of readability

  Select Case (ZCoord%CoordType)

    ! (i) P based (terrain following) coord system.
    Case (Z_PAsEta)

      ! (a) Surface pressure.

      ! Surface pressure available.
      If (M%FieldPresent3d(1, F_PAsRead)) Then

        M%P(:, :, 1, MNew) = M%PAsRead(:, :, 1, MNew)

      ! PSeaLevel available but no PAsRead values available. Here we estimate
      ! P(ground) / P(sea level) and use P(sea level) to estimate surface pressure.
      Else If (M%FieldPresent2d(F_PSeaLevel) .and. .not. LowestLevel <= ZGridP%nZ) Then

        ! P(ground) = P(sea level) * Exp( - Topog * Gravity / (GasConstant * T)).
        M%P(:, :, 1, MNew) = M%PSeaLevel(:, :, MNew) *            &
                              Exp(                                  &
                                - M%Topog(:, :, MNew) * Gravity /  &
                                (GasConstant * M%T(:, :, 1, MNew)) &
                              )

      ! PAsRead values available but PSeaLevel not available. Here we estimate
      ! P(lowest level present) / P(ground) and use P(lowest level present) to
      ! estimate surface pressure.
      Else If (LowestLevel <= ZGridP%nZ .and. .not. M%FieldPresent2d(F_PSeaLevel)) Then

        ! We estimate P(lowest level present) / P(ground) by assuming surface pressure
        ! is the ICAO standard atmosphere value of P(sea level) with a topography
        ! correction. This is fine near the ground where eta = pressure / surface
        ! pressure and it is dangerous to attempt a more accurate calculation higher
        ! up because the eta-p relation may become insensitive to surface pressure.
        Do j = 1, HGrid%nY
        Do i = 1, HGrid%nX
          PLowestLevelByPGround =                                        &
            PBasedToPBased(                                              &
              ZCoordIn  = ZCoord,                                        &
              ZCoordOut = ZCoordP,                                       &
              PS        = PAt0km * Exp(                                  &
                                     - M%Topog(i, j, MNew) * Gravity /  &
                                     (GasConstant * M%T(i, j, 1, MNew)) &
                                   ),                                    &
              ZIn       = ZGridP%Z(LowestLevel)                          &
            )                                                            &
            / PAt0km  ! $$ this looks wrong - check.
          M%P(i, j, 1, MNew) = M%PAsRead(i, j, LowestLevel, MNew) / &
                                PLowestLevelByPGround
        End Do
        End Do

      ! PSeaLevel and PAsRead values available. Here we compare estimates of
      ! P(lowest level present) / P(ground) and P(ground) / P(sea level) to decide
      ! whether to use P(sea level) or P(lowest level present) in estimating
      ! P(ground).
      Else If (M%FieldPresent2d(F_PSeaLevel) .and. LowestLevel <= ZGridP%nZ) Then

        Do j = 1, HGrid%nY
        Do i = 1, HGrid%nX
          PLowestLevelByPGround =                                        &
            PBasedToPBased(                                              &
              ZCoordIn  = ZCoord,                                        &
              ZCoordOut = ZCoordP,                                       &
              PS        = PAt0km * Exp(                                  &
                                     - M%Topog(i, j, MNew) * Gravity /  &
                                     (GasConstant * M%T(i, j, 1, MNew)) &
                                   ),                                    &
              ZIn       = ZGridP%Z(LowestLevel)                          &
            )                                                            &
            / PAt0km
          PGroundByPSea  = Exp(                                  &
                             - M%Topog(i, j, MNew) * Gravity /  &
                             (GasConstant * M%T(i, j, 1, MNew)) &
                           )
          If (PLowestLevelByPGround < PGroundByPSea) Then
            M%P(i, j, 1, MNew) = M%PSeaLevel(i, j, MNew) * PGroundByPSea
          Else
            M%P(i, j, 1, MNew) = M%PAsRead(i, j, LowestLevel, MNew) &
                                  / PLowestLevelByPGround
          End If
        End Do
        End Do

      ! PSeaLevel and PAsRead values not available. Here we estimate P(ground) /
      ! P(sea level) and use the ICAO standard atmosphere value of P(sea level) to
      ! estimate surface pressure.
      Else

        ! P(ground) = PICAO(sea level) * Exp( - Topog * Gravity / (GasConstant * T)).
        M%P(:, :, 1, MNew) = PAt0km *                              &
                              Exp(                                  &
                                - M%Topog(:, :, MNew) * Gravity /  &
                                (GasConstant * M%T(:, :, 1, MNew)) &
                              )

      End If

      ! (b) Pressure above the surface.

      Do k = 2, ZGrid%nZ
      Do j = 1, HGrid%nY
      Do i = 1, HGrid%nX
        M%P(i, j, k, MNew) = PBasedToPBased(                    &
                                ZCoordIn  = ZCoord,              &
                                ZCoordOut = ZCoordP,             &
                                PS        = M%P(i, j, 1, MNew), &
                                ZIn       = ZGrid%Z(k)           &
                              )
      End Do
      End Do
      End Do

      ! (c) Rho = P / (GasConstant * T).
      M%Rho(:, :, :, MNew) = M%P(:, :, :, MNew) /               &
                              (GasConstant * M%T(:, :, :, MNew))

      ! (d) Theta = T * (PRef / P) ** (GasConstant/Cp).
      M%Theta(:, :, :, MNew) = M%T(:, :, :, MNew) *                            &
                                (PRef / M%P(:, :, :, MNew)) ** (GasConstant/Cp)


      ! (e) Z: Z2 - Z1 = GasConstant * T * Log(P1/P2) / Gravity.
      M%Z(:, :, 1, MNew) = 0.0
      Do k = 2, ZGrid%nZ
        M%Z(:, :, k, MNew) =                                     &
          M%Z(:, :, k - 1, MNew) +                               &
          GasConstant *                                           &
          (M%T(:, :, k - 1, MNew) + M%T(:, :, k, MNew)) * 0.5 * &
          Log(M%P(:, :, k - 1, MNew) / M%P(:, :, k, MNew)) /    &
          Gravity
      End Do

    ! (ii) Z based (terrain following) coord system.
    Case (Z_AboveGround, Z_ZAsEta)

      ! (a) Z.
      Do k = 1, ZGrid%nZ
      Do j = 1, HGrid%nY
      Do i = 1, HGrid%nX
        M%Z(i, j, k, MNew) = ZBasedToZBased(                     &
                                ZCoordIn  = ZCoord,               &
                                ZCoordOut = ZCoordZ,              &
                                Topog     = M%Topog(i, j, MNew), &
                                ZIn       = ZGrid%Z(k)            &
                              )
      End Do
      End Do
      End Do

      ! (b) Pressure / surface pressure:
      ! P2 = P1 * Exp( - (Z2 - Z1) * Gravity / (GasConstant * T)).
      M%P(:, :, 1, MNew) = 1.0
      Do k = 2, ZGrid%nZ
        M%P(:, :, k, MNew) = M%P(:, :, k - 1, MNew) *                          &
          Exp(                                                                   &
            - (M%Z(:, :, k, MNew) - M%Z(:, :, k - 1, MNew)) * Gravity /        &
            (GasConstant * (M%T(:, :, k, MNew) + M%T(:, :, k - 1, MNew)) * 0.5 &
            )                                                                    &
          )
      End Do

      ! (c) Surface pressure.

      ! Surface pressure available.
      If (M%FieldPresent3d(1, F_PAsRead)) Then

        M%P(:, :, 1, MNew) = M%PAsRead(:, :, 1, MNew)

      ! PSeaLevel available but no PAsRead values available. Here we estimate
      ! P(ground) / P(sea level) and use P(sea level) to estimate surface pressure.
      Else If (M%FieldPresent2d(F_PSeaLevel) .and. .not. LowestLevel <= ZGridP%nZ) Then

        ! P(ground) = P(sea level) * Exp( - Topog * Gravity / (GasConstant * T)).
        M%P(:, :, 1, MNew) = M%PSeaLevel(:, :, MNew) *            &
                              Exp(                                  &
                                - M%Topog(:, :, MNew) * Gravity /  &
                                (GasConstant * M%T(:, :, 1, MNew)) &
                              )

      ! PAsRead values available but PSeaLevel not available. Here we estimate
      ! P(lowest level present) / P(ground) and use P(lowest level present) to
      ! estimate surface pressure.
      Else If (LowestLevel <= ZGridP%nZ .and. .not. M%FieldPresent2d(F_PSeaLevel)) Then

        Do j = 1, HGrid%nY
        Do i = 1, HGrid%nX
          PLowestLevelByPGround = GInterpXYZ(                  &
                                    M%P(:, :, :, MNew), i, j, &
                                    M%GHCoeffs,                &
                                    M%ZCoeffsToP(LowestLevel)  &
                                  )
          M%P(i, j, 1, MNew) = M%PAsRead(i, j, LowestLevel, MNew) / &
                                PLowestLevelByPGround
        End Do
        End Do

      ! PSeaLevel and PAsRead values available. Here we compare estimates of
      ! P(lowest level present) / P(ground) and P(ground) / P(sea level) to decide
      ! whether to use P(sea level) or P(lowest level present) in estimating
      ! P(ground).
      Else If (M%FieldPresent2d(F_PSeaLevel) .and. LowestLevel <= ZGridP%nZ) Then

        Do j = 1, HGrid%nY
        Do i = 1, HGrid%nX
          PLowestLevelByPGround = GInterpXYZ(                  &
                                    M%P(:, :, :, MNew), i, j, &
                                    M%GHCoeffs,                &
                                    M%ZCoeffsToP(LowestLevel)  &
                                  )
          PGroundByPSea  = Exp(                                  &
                             - M%Topog(i, j, MNew) * Gravity /  &
                             (GasConstant * M%T(i, j, 1, MNew)) &
                           )
          If (PLowestLevelByPGround < PGroundByPSea) Then
            M%P(i, j, 1, MNew) = M%PSeaLevel(i, j, MNew) * PGroundByPSea
          Else
            M%P(i, j, 1, MNew) = M%PAsRead(i, j, LowestLevel, MNew) &
                                  / PLowestLevelByPGround
          End If
        End Do
        End Do

      ! PSeaLevel and PAsRead values not available. Here we estimate P(ground) /
      ! P(sea level) and use the ICAO standard atmosphere value of P(sea level) to
      ! estimate surface pressure.
      Else

        ! P(ground) = PICAO(sea level) * Exp( - Topog * Gravity / (GasConstant * T)).
        M%P(:, :, 1, MNew) = PAt0km *                              &
                              Exp(                                  &
                                - M%Topog(:, :, MNew) * Gravity /  &
                                (GasConstant * M%T(:, :, 1, MNew)) &
                              )

      End If

      ! (d) Pressure above the surface:
      ! P(level k) = P(level 1) * (P(level k) / P(level 1)).
      Do k = 2, ZGrid%nZ
        M%P(:, :, k, MNew) = M%P(:, :, 1, MNew) * M%P(:, :, k, MNew)
      End Do

      ! (e) Rho = P / (GasConstant * T).
      M%Rho(:, :, :, MNew) = M%P(:, :, :, MNew) /               &
                              (GasConstant * M%T(:, :, :, MNew))

      ! (f) Theta = T * (PRef / P) ** (GasConstant/Cp).
      M%Theta(:, :, :, MNew) = M%T(:, :, :, MNew) *                            &
                                (PRef / M%P(:, :, :, MNew)) ** (GasConstant/Cp)

    Case Default

      Call Message(                                   &
             'UNEXPECTED ERROR in ProcessNWPMet: ' // &
             'unexpected coordinate system type',     &
             4                                        &
           )

  End Select

  ! 4) UStar, H and W.

  ! Calculate UStar.
  Do j = 1, HGrid%nY
  Do i = 1, HGrid%nX
    UStress = GInterpXY(M%UStress(:, :, MNew), i, j, M%GHCoeffsU)
    VStress = GInterpXY(M%VStress(:, :, MNew), i, j, M%GHCoeffsV)
    M%UStar(i, j, MNew) = Sqrt(Sqrt(UStress**2 + VStress**2)/M%Rho(i, j, 1, MNew))
  End Do
  End Do

  ! Calculate H.
  If (.not.M%UseNWPH .or. .not.M%FieldPresent2d(F_H)) Then
    Call CalcH(Coords, Grids, M, MNew)
  End If

  ! Apply minimum/maximum boundary layer depth of HMin/HMax.
  M%H(:, :, MNew) = Max(M%H(:, :, MNew), M%HMin)
  M%H(:, :, MNew) = Min(M%H(:, :, MNew), M%HMax)

  ! Calculate eta dot when the Eulerian advection scheme is being used.
  ! $$ Replace EulerianModelGlobal with DispOpts%EulerianModel
  If (EulerianModelGlobal) Call CalcZAboveGroundDot(Coords, Grids, M, MNew, .true.)

  ! Convert W to rate of change of height above ground (m/s).
  Call CalcZAboveGroundDot(Coords, Grids, M, MNew, .false.)

  ! 5) Cloud and Rain.

  ! Convert from kg/(m^2 s) to mm/hr.
  M%DynRain(:, :, MNew) = M%DynRain(:, :, MNew) * 3600.0
  M%DynSnow(:, :, MNew) = M%DynSnow(:, :, MNew) * 3600.0
  M%ConRain(:, :, MNew) = M%ConRain(:, :, MNew) * 3600.0
  M%ConSnow(:, :, MNew) = M%ConSnow(:, :, MNew) * 3600.0

  ! Sum rain and snow rates to give precipitation rates.
  M%DynPpt(:, :, MNew) = M%DynRain(:, :, MNew) + M%DynSnow(:, :, MNew)
  M%ConPpt(:, :, MNew) = M%ConRain(:, :, MNew) + M%ConSnow(:, :, MNew)

  ! fixups for possible bad um met (originally in NAME routine Met) $$
  M%DynPpt(:, :, MNew) = Max(M%DynPpt(:, :, MNew), 0.0)
  M%ConPpt(:, :, MNew) = Max(M%ConPpt(:, :, MNew), 0.0)

  ! Initialise values of ConCloudBasePa and ConCloudTopPa to ConCloudBase and ConCloudTop
  ! in order to capture all negative values corresponding to no ConCloud.
  M%ConCloudBasePa(:, :, MNew) = M%ConCloudBase(:, :, MNew)
  M%ConCloudTopPa(:, :, MNew) = M%ConCloudTop(:, :, MNew)
  ! Convert convective cloud base and top to Pa and to m. 
  ! $$ conversion uses instantaneous pressure and temperature
  ! $$ should perhaps be mean pressure and temperature in line with mean
  ! $$ cloud
  Do j = 1, HGrid%nY
  Do i = 1, HGrid%nX
    If (M%ConCloudBase(i, j, MNew) >= 0.0) Then
      M%ConCloudBasePa(i, j, MNew) =  ConvertZColumn(                &
                                      Coords, Grids, M, i, j,     &
                                      ZCoordCl, ZCoordP,          &
                                      M%ConCloudBase(i, j, MNew), &
                                      MNew                        &
                                    )
      M%ConCloudBase(i, j, MNew) = ConvertZColumn(                &
                                      Coords, Grids, M, i, j,     &
                                      ZCoordCl, ZCoordZ,          &
                                      M%ConCloudBase(i, j, MNew), &
                                      MNew                        &
                                    )
      ! Fix up if cloud base calculated to be below ground due to
      ! use of potentially incorrect pressure / temperature (cloud 
      ! base may not be interpolated in time whereas pressure / 
      ! temperature interpolated)
      M%ConCloudBase(i, j, MNew) = Max(                           &
                                      M%ConCloudBase(i, j, MNew), &
                                      0.0                         &
                                    )
      ! If M%ConCloudBase(i, j, MNew) > 0 then also calculate the nearest level to the Freezing Level (in pressure)
      minTempDiff = 100.0                                
      Do k = 1, ZGrid%nZ
         TempDiffToFL = abs( M%T(i, j, k, MNew) - 273.15) 
         If (TempDiffToFL < minTempDiff) then
            MinTempDiff = TempDiffToFL
            fl=k
         End If
      End Do
      M%FLPa(i, j, MNew) = M%P(i, j, fl, MNew)
    End If
    If (M%ConCloudTop(i, j, MNew) >= 0.0) Then
      M%ConCloudTopPa(i, j, MNew) = ConvertZColumn(               &
                                     Coords, Grids, M, i, j,    &
                                     ZCoordCl, ZCoordP,         &
                                     M%ConCloudTop(i, j, MNew), &
                                     MNew                       &
                                   )
      M%ConCloudTop(i, j, MNew) = ConvertZColumn(               &
                                     Coords, Grids, M, i, j,    &
                                     ZCoordCl, ZCoordZ,         &
                                     M%ConCloudTop(i, j, MNew), &
                                     MNew                       &
                                   )
            
      ! Unlikely fix up if cloud top calculated to be below ground 
      ! due to use of potentially incorrect pressure / temperature 
      ! (cloud base may not be interpolated in time whereas pressure / 
      ! temperature interpolated). Fix up to 10m to ensure not equal 
      ! to cloud base. 
      If (M%ConCloudTop(i, j, MNew) <= 0.0) Then
        M%ConCloudTop(i, j, MNew) = Max(                          &
                                       M%ConCloudTop(i, j, MNew), &
                                       10.0                       &
                                     )
      End If
    End If
  End Do
  End Do

  ! Calculate cloud information (Cloud3d, Cloud, TotalOrDynCloudBase
  ! and TotalOrDynCloudTop).

  Call CalcCloudInfo(Coords, Grids, M, MNew)
  
!  Code for correcting vertical windspeeds for mass consistency

!  MTH=ZCoord%ModelTopHeight
!  Do i=2, HGrid%nX
!    Do j=2, HGrid%nY - 1 
!      Point = [HGrid%X(i),HGrid%Y(j)]
!      Call MetricCoeffs(HCoord, Point, HMax, H1, H2)
!      Do k=2, ZGrid%nZ
!        M%W(i, j, k, MNew) = -((M%U(i, j, k, MNew) - M%U(i-1, j, k, MNew)) / H1  &
!                           + (M%V(i, j, k, MNew) - M%V(i, j-1, k, MNew)) / H2 )  &
!                           * (M%Z(i, j, k, MNew) - M%Z(i, j, k-1, MNew))         &
!                           + M%W(i, j, k-1, MNew)
!      End Do
!    End Do
!  End Do


End Subroutine ProcessNWPMet

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcH(Coords, Grids, M, NewIdx)
! Calculates the boundary layer depth at model grid points.

! The boundary layer depth is calculated using both a Richardson number method and a
! dry adiabat method. The maximum of the two values is taken.

  Implicit None
  ! Argument list:
  Type(Coords_), Intent(In),   Target :: Coords ! Collection of coord systems.
  Type(Grids_),  Intent(In),   Target :: Grids  ! Collection of grids.
  Type(NWPMet_), Intent(InOut)        :: M      ! State of an NWP met module instance.
  Integer,       Intent(In)           :: NewIdx
  Integer                             :: MNew
  ! Local parameters:
  Real(Std), Parameter :: RiCrit   = 1.3 ! Critical Richardson number.
  Real(Std), Parameter :: AddDay   = 0.5 ! Daytime temperature offset.
  Real(Std), Parameter :: AddNight = 0.5 ! Nighttime temperature offset.
  ! Locals:
  Real(Std)                  :: ThetaS   ! Surface potential temperature plus offset.
  Real(Std)                  :: Ri       ! Richardson number.
  Real(Std)                  :: HRi      ! Boundary layer depth using the Richardson
                                         ! number method.
  Real(Std)                  :: HAdiabat ! Boundary layer depth using the dry adiabat
                                         ! method.
  Real(Std),     Allocatable :: T(:)     ! Temperature profile.
  Real(Std),     Allocatable :: Theta(:) ! Potential temperature profile.
  Real(Std),     Allocatable :: U(:)     ! Wind (u-cpt) profile.
  Real(Std),     Allocatable :: V(:)     ! Wind (v-cpt) profile.
  Type(HGrid_),  Pointer     :: HGrid    !} Abbreviations for grids.
  Type(ZGrid_),  Pointer     :: ZGrid    !}
  Type(ZGrid_),  Pointer     :: ZGridUV  !}
  Integer                    :: i        !] Loop indices.
  Integer                    :: j        !]
  Integer                    :: k        !]

  MNew = NewIdx

  ! Set up abbreviations for grids.
  HGrid   => Grids%HGrids(M%iHGrid  )
  ZGrid   => Grids%ZGrids(M%iZGrid  )
  ZGridUV => Grids%ZGrids(M%iZGridUV)

  Allocate(T    (ZGrid  %nZ))
  Allocate(Theta(ZGrid  %nZ))
  Allocate(U    (ZGridUV%nZ))
  Allocate(V    (ZGridUV%nZ))

  Do j = 1, HGrid%nY
  Do i = 1, HGrid%nX

    ! Compute wind and potential temperature profiles.
    Do k = 2, ZGrid%nZ
      T(k)     = M%T(i, j, k, MNew)
      Theta(k) = M%Theta(i, j, k, MNew)
    End Do
    Do k = 2, ZGridUV%nZ
      U(k)     = GInterpXY(M%U(:, :, k, MNew), i, j, M%GHCoeffsU)
      V(k)     = GInterpXY(M%V(:, :, k, MNew), i, j, M%GHCoeffsV)
    End Do

    ! Initialise HRi and HAdiabat to top of ZGrid.
    HRi      = M%Z(i, j, ZGrid%nZ, MNew)
    HAdiabat = HRi

    ! Calculate the boundary layer depth by the Richardson number method.
    ! Need to revise to do new dynamics properly $$ U, V staggering is different
    Do k = 3, ZGrid%nZ - 1
      Ri = (U(k + 1) - U(k)) ** 2 + (V(k + 1) - V(k)) ** 2
      Ri = Max(  (  (T(k) + T(k + 1))/2.0   ) * Ri, 1.E-5)
      Ri = Gravity * (Theta(k + 1) - Theta(k)) *             &
           (M%Z(i, j, k + 1, MNew) - M%Z(i, j, k, MNew)) / &
           Ri
      If (Ri > RiCrit) Then
        HRi = M%Z(i, j, k, MNew)
        Exit
      End If
    End Do

    ! Calculate the boundary layer depth by the dry adiabat method, adding the
    ! appropriate offset to the surface temperature. If the surface is colder than the
    ! 2nd level, then it is assumed to be night time (or stable).
    If (Theta(2) < Theta(3)) Then
      ThetaS = Theta(3) + AddNight
    Else
      ThetaS = Theta(3) + AddDay
    End If
    Do k = 3, ZGrid%nZ - 1
      If (Theta(k + 1) > ThetaS) Then
        HAdiabat = M%Z(i, j, k, MNew) +                             &
                   (M%Z(i, j, k + 1, MNew) - M%Z(i, j, k, MNew)) * &
                   (ThetaS       - Theta(k)) /                       &
                   (Theta(k + 1) - Theta(k))
        Exit
      End If
    End Do

    ! Select highest boundary layer depth.
    M%H(i, j, MNew) = Max(HAdiabat, HRi)

  End Do
  End Do

  DeAllocate(T    )
  DeAllocate(Theta)
  DeAllocate(U    )
  DeAllocate(V    )

End Subroutine CalcH

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcZAboveGroundDot(Coords, Grids, M, NewIdx, CalcEtaDot)
! Calculates rate of change of height above ground in metres following the mean flow
! from the input vertical velocity (which may well be defined in a different coord
! system).

  Implicit None
  ! Argument List:
  Type(Coords_), Intent(In),    Target :: Coords      ! Collection of coord systems.
  Type(Grids_),  Intent(In),    Target :: Grids       ! Collection of grids.
  Type(NWPMet_), Intent(InOut), Target :: M           ! State of an NWP met module instance.
  Integer,       Intent(In)            :: NewIdx
  Logical,       Intent(In)            :: CalcEtaDot  ! If true calculates eta dot; if false calculates z_agl

  ! Locals:
  Real(Std)                :: HMax            ! Max metric coefficient.
  Real(Std)                :: dX              ! Distance in longitudinal direction in
                                              ! metres for unit change in X coord.
  Real(Std)                :: dY              ! Distance in latitudinal direction in
                                              ! metres for unit change in Y coord.
  Real(Std)                :: dZdX(3)         ! d Z(AboveGround) / d (X(m),Y(m),Z(w)) with dX(m) and dY(m)
                                              ! denoting dX and dY in metres and Z(w) being the height in the
                                              ! coord system for which the input vertical velocity is Z(w)Dot.
  Real(Std)                :: dEtadX(3)       ! d Eta / d (X(m),Y(m),Z(w)).
  Real(Std)                :: PStar           ! Surface pressure in Pa.
  Real(Std)                :: dPStardX(2)     ! Horizontal gradient of surface pressure.
  Type(HCoord_),   Pointer :: HCoord          !] Abbreviations for coords, grids.
  Type(ZCoord_),   Pointer :: ZCoord          !]
  Type(ZCoord_),   Pointer :: ZCoordW         !]
  Type(ZCoord_),   Pointer :: ZCoordZ         !]
  Type(HGrid_),    Pointer :: HGrid           !]
  Type(ZGrid_),    Pointer :: ZGrid           !]
  Type(ZGrid_),    Pointer :: ZGridW          !]
  Type(GHCoeffs_), Pointer :: GHCoeffs        !} Interpolation coefficients for
  Type(GHCoeffs_), Pointer :: GHCoeffsU       !} various grids (including for
  Type(GHCoeffs_), Pointer :: GHCoeffsV       !} obtaining derivatives of fields).
  Type(GHCoeffs_), Pointer :: GHCoeffsdX      !}
  Type(GHCoeffs_), Pointer :: GHCoeffsdY      !}
  Type(ZCoeffs_),  Pointer :: ZCoeffsToW(:)   !}
  Type(ZCoeffs_),  Pointer :: ZCoeffsUVToW(:) !}
  Type(ZCoeffs_)           :: dZCoeffsToWdZ   !}
  Integer                  :: i               !] Loop indices.
  Integer                  :: j               !]
  Integer                  :: k               !]
  Integer                  :: MNew
  Real(Std)                :: Temp(2)         ! Temporary array needed to avoid array-copying warnings on
                                              ! Linux Intel compiler with -C option.
                                              ! $$ review whether still needed.
  MNew = NewIdx

  ! Set up abbreviations for coords, grids and interpolation coefficients.
  HCoord       => Coords%HCoords(M%iHCoord )
  ZCoord       => Coords%ZCoords(M%iZCoord )
  ZCoordW      => Coords%ZCoords(M%iZCoordW)
  If (CalcEtaDot) Then
    ZCoordZ      => Coords%ZCoords(M%iZCoord)
  Else
    ZCoordZ      => Coords%ZCoords(M%iZCoordZ)
  EndIf
  HGrid        => Grids%HGrids(M%iHGrid )
  ZGrid        => Grids%ZGrids(M%iZGrid )
  ZGridW       => Grids%ZGrids(M%iZGridW)
  GHCoeffs     => M%GHCoeffs
  GHCoeffsU    => M%GHCoeffsU
  GHCoeffsV    => M%GHCoeffsV
  GHCoeffsdX   => M%GHCoeffsdX
  GHCoeffsdY   => M%GHCoeffsdY
  ZCoeffsToW   => M%ZCoeffsToW
  ZCoeffsUVToW => M%ZCoeffsUVToW

  ! Notes on derivatives: All differences are evaluated as 2 - 1, with the denominator
  ! dX, dY, dZ or dEta chosen consistently.

  ! Note that vertical velocity has already been set to zero at the surface and is skipped here.

  ! $$ does this routine implicitly assume that HGrid%dX,dY = 1 ?

  Select Case (ZCoordW%CoordType)

    ! Input vertical velocity is rate of change of a height-based coordinate.
    Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)

      Do i = 1, HGrid%nX
      Do j = 1, HGrid%nY

        ! Calculate dX, dY.
        Temp = (/ HGrid%X(i), HGrid%Y(j) /)
        Call MetricCoeffs(HCoord, Temp, HMax, dX, dY)

        Do k = 2, ZGridW%nZ

          ! Calculate d Z(AboveGround) / d (X(m), Y(m), Z(w)) with dX(m) and dY(m)
          ! denoting dX and dY in metres and Z(w) being the height in the coord system
          ! for which the input vertical velocity is Z(w) Dot.
          ! Note the last if-block will work for all cases, but the other if-blocks
          ! offer a slightly faster execution by not evaluating irrelevant quantities.
          If (ZCoordW%CoordType == Z_AboveSea .and. .not.CalcEtaDot) Then
            Temp = (/                                                       &
                     GInterpXY(M%Topog(:, :, MNew), i, j, GHCoeffsdX) / dX, &
                     GInterpXY(M%Topog(:, :, MNew), i, j, GHCoeffsdY) / dY  &
                   /)
            Call CalcdZdXZBased(     &
                   ZCoordZ, ZCoordW, &
                   1.0,              &
                   1.0,              &
                   Temp,             &
                   dZdX              &
                 )
          Else If (ZCoordW%CoordType == Z_AboveGround .and. .not.CalcEtaDot) Then
            Temp = (/ 1.0, 1.0 /)
            Call CalcdZdXZBased(     &
                   ZCoordZ, ZCoordW, &
                   1.0,              &
                   1.0,              &
                   Temp,             &
                   dZdX              &
                 )
          Else
            Temp = (/                                                       &
                     GInterpXY(M%Topog(:, :, MNew), i, j, GHCoeffsdX) / dX, &
                     GInterpXY(M%Topog(:, :, MNew), i, j, GHCoeffsdY) / dY  &
                   /)
            Call CalcdZdXZBased(                                                  &
                   ZCoordZ, ZCoordW,                                              &
                   GInterpXYZ(M%Z(:, :, :, MNew), i, j, GHCoeffs, ZCoeffsToW(k)), &
                   M%Topog(i, j, MNew),                                           &
                   Temp,                                                          &
                   dZdX                                                           &
                 )
          End If

          ! Calculate Z(AboveGround) Dot as d Z(AboveGround) / d (X(m), Y(m), Z(w)) * (U, V, Z(w) Dot).
          If (CalcEtaDot) Then
            M%EtaDot(i, j, k, MNew) =                                             &
              dZdX(1) *                                                           &
              GInterpXYZ(M%U(:, :, :, MNew), i, j, GHCoeffsU, ZCoeffsUVToW(k)) +  &
              dZdX(2) *                                                           &
              GInterpXYZ(M%V(:, :, :, MNew), i, j, GHCoeffsV, ZCoeffsUVToW(k)) +  &
              dZdX(3) *                                                           &
              M%W(i, j, k, MNew)
          Else
            M%W(i, j, k, MNew) =                                                  &
              dZdX(1) *                                                           &
              GInterpXYZ(M%U(:, :, :, MNew), i, j, GHCoeffsU, ZCoeffsUVToW(k)) +  &
              dZdX(2) *                                                           &
              GInterpXYZ(M%V(:, :, :, MNew), i, j, GHCoeffsV, ZCoeffsUVToW(k)) +  &
              dZdX(3) *                                                           &
              M%W(i, j, k, MNew)
          End If

        End Do

      End Do
      End Do

    ! Input vertical velocity is rate of change of a pressure-based coordinate.
    ! Note that the d/dt tendencies cannot be evaluated here (because they depend on met information
    ! at two times and will also be recalculated at the next met update). These additional terms are
    ! therefore treated later as a correction WNoTFix --> W of the value WNoTFix (the W array without the
    ! d/dt correction term) calculated here.
    Case (Z_P, Z_PAsZ, Z_PAsEta)

      Do i = 1, HGrid%nX
      Do j = 1, HGrid%nY

        ! Calculate dX, dY.
        Temp = (/ HGrid%X(i), HGrid%Y(j) /)
        Call MetricCoeffs(HCoord, Temp, HMax, dX, dY)

        Do k = 2, ZGridW%nZ

          ! Calculate d Z(AboveGround) / d (X(m), Y(m), Z(w)) with dX(m) and dY(m)
          ! denoting dX and dY in metres and Z(w) being the height in the coord system
          ! for which the input vertical velocity is Z(w) Dot. This vector is computed
          ! as the product
          !
          !  d Z(AboveGround) / d (X(m), Y(m), Eta) * d (X(m), Y(m), Eta) /  d (X(m), Y(m), Z(w))
          !
          ! where Eta is the vertical coord system of the main grids (i.e. model levels).
          ! When Z(w) = Eta, the second term here is trivial.

          ! Calculate vertical interpolation coefficients for calculating d Z(AboveGround) / d Eta.
          ! Note that we use values on ZGrid either side of ZGridW levels.
          ! $$ Asymmetric differences if W grid = main grid
          dZCoeffsToWdZ    =  ZCoeffsToW(k)
          dZCoeffsToWdZ%Z1 =  1.0
          dZCoeffsToWdZ%Z2 = -1.0

          ! Calculate d Z(AboveGround) / d (X(m), Y(m), Eta).
          dZdX(1) = GInterpXYZ(M%Z(:, :, :, MNew), i, j, GHCoeffsdX, ZCoeffsToW(k)) &
                    / dX
          dZdX(2) = GInterpXYZ(M%Z(:, :, :, MNew), i, j, GHCoeffsdY, ZCoeffsToW(k)) &
                    / dY
          dZdX(3) = GInterpXYZ(M%Z(:, :, :, MNew), i, j, GHCoeffs,   dZCoeffsToWdZ) &
                    / (ZGrid%Z(ZCoeffsToW(k)%iZ2) - ZGrid%Z(ZCoeffsToW(k)%iZ1))

          ! Calculate d Z(AboveGround) / d (X(m), Y(m), Z(w)).
          ! Note that dZdX is precisely this term already when Z(w) = Eta.
          ! $$ Note this assumes ZCoord (i.e. eta) is pressure based when ZCoordW is.
          If (M%iZCoordW /= M%iZCoord) Then
            ! Calculate surface pressure and the horizontal gradient of surface pressure.
            ! $$ Independent of height - move outside k loop?
            PStar       = GInterpXY(M%P(:, :, 1, MNew), i, j, GHCoeffs  )
            dPStardX(1) = GInterpXY(M%P(:, :, 1, MNew), i, j, GHCoeffsdX) / dX
            dPStardX(2) = GInterpXY(M%P(:, :, 1, MNew), i, j, GHCoeffsdY) / dY

            ! Calculate transformation coefficients d Eta / d(X, Y, Z(w)).
            Call CalcdZdXPBased(ZCoord, ZCoordW, ZGridW%Z(k), PStar, dPStardX, dEtadX)

            ! Calculate revised dZdX terms.
            dZdX(1) = dZdX(1) + dZdX(3) * dEtadX(1)
            dZdX(2) = dZdX(2) + dZdX(3) * dEtadX(2)
            dZdX(3) =           dZdX(3) * dEtadX(3)
          End If

          ! Calculate Z(AboveGround) Dot as d Z(AboveGround) / d (X(m), Y(m), Z(w)) *
          ! (U, V, Z(w) Dot).
          M%W(i, j, k, MNew) =                                                 &
            dZdX(1) *                                                          &
            GInterpXYZ(M%U(:, :, :, MNew), i, j, GHCoeffsU, ZCoeffsUVToW(k)) + &
            dZdX(2) *                                                          &
            GInterpXYZ(M%V(:, :, :, MNew), i, j, GHCoeffsV, ZCoeffsUVToW(k)) + &
            dZdX(3) *                                                          &
            M%W(i, j, k, MNew)

          ! Save calculated values of W into the WNoTFix array.
          M%WNoTFix(i, j, k, MNew) = M%W(i, j, k, MNew)

        End Do

      End Do
      End Do

    Case Default

      Call Message(                                         &
             'UNEXPECTED ERROR in CalcZAboveGroundDot: ' // &
             'unknown coordinate system type',              &
             4                                              &
           )

  End Select

  ! Might be a good idea to damp ZAboveGroundDot to zero at the ground. $$
  ! Also adjust surface U, V values (although these are not normally input)

End Subroutine CalcZAboveGroundDot

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcZDotCorrection(Coords, Grids, M)
! Calculates d/dt correction term needed to complete the conversion of vertical velocity from a pressure-based
! vertical velocity to the rate of change of height above ground in metres.

  Implicit None
  ! Argument List:
  Type(Coords_), Intent(In),    Target :: Coords ! Collection of coord systems.
  Type(Grids_),  Intent(In),    Target :: Grids  ! Collection of grids.
  Type(NWPMet_), Intent(InOut), Target :: M      ! State of an NWP met module instance.
  ! Locals:
  Real(Std)                :: dZdEtaNew       !} d Z(AboveGround) / d Eta at the current and previous met
  Real(Std)                :: dZdEtaOld       !} update times.
  Real(Std)                :: dZdt            ! d Z(AboveGround) / d t (at constant Eta).
  Real(Std)                :: dEtadX(3)       ! First component stores d Eta / dt at constant Z(w); other
                                              ! components not explicitly used here.
  Real(Std)                :: PStar           ! Surface pressure in Pa.
  Real(Std)                :: dPStardX(2)     ! First component is surface pressure tendency, second component
                                              ! is not used.
  Real(Std)                :: WCorrNew        !} Correction to vertical velocity at the current and previous
  Real(Std)                :: WCorrOld        !} met update times.
  Type(HCoord_),   Pointer :: HCoord          !] Abbreviations for coords, grids.
  Type(ZCoord_),   Pointer :: ZCoord          !]
  Type(ZCoord_),   Pointer :: ZCoordW         !]
  Type(ZCoord_),   Pointer :: ZCoordZ         !]
  Type(HGrid_),    Pointer :: HGrid           !]
  Type(ZGrid_),    Pointer :: ZGrid           !]
  Type(ZGrid_),    Pointer :: ZGridW          !]
  Type(GHCoeffs_), Pointer :: GHCoeffs        !} Interpolation coefficients for various grids (including for
  Type(GHCoeffs_), Pointer :: GHCoeffsU       !} obtaining derivatives of fields).
  Type(GHCoeffs_), Pointer :: GHCoeffsV       !}
  Type(GHCoeffs_), Pointer :: GHCoeffsdX      !}
  Type(GHCoeffs_), Pointer :: GHCoeffsdY      !}
  Type(ZCoeffs_),  Pointer :: ZCoeffsToW(:)   !}
  Type(ZCoeffs_),  Pointer :: ZCoeffsUVToW(:) !}
  Type(ZCoeffs_)           :: dZCoeffsToWdZ   !}
  Integer                  :: i               !] Loop indices.
  Integer                  :: j               !]
  Integer                  :: k               !]

  ! Set up abbreviations for coords, grids and interpolation coefficients.
  HCoord       => Coords%HCoords(M%iHCoord )
  ZCoord       => Coords%ZCoords(M%iZCoord )
  ZCoordW      => Coords%ZCoords(M%iZCoordW)
  ZCoordZ      => Coords%ZCoords(M%iZCoordZ)
  HGrid        => Grids%HGrids(M%iHGrid )
  ZGrid        => Grids%ZGrids(M%iZGrid )
  ZGridW       => Grids%ZGrids(M%iZGridW)
  GHCoeffs     => M%GHCoeffs
  GHCoeffsU    => M%GHCoeffsU
  GHCoeffsV    => M%GHCoeffsV
  GHCoeffsdX   => M%GHCoeffsdX
  GHCoeffsdY   => M%GHCoeffsdY
  ZCoeffsToW   => M%ZCoeffsToW
  ZCoeffsUVToW => M%ZCoeffsUVToW

  ! Notes on derivatives: All differences are evaluated as 2 - 1, with the denominator
  ! dX, dY, dZ or dEta chosen consistently.

  ! Note that the d/dt tendencies cannot be evaluated in the usual met processing step because they
  ! depend on the change in met between two times. Each tendency is also recalculated at the next met
  ! update. These corrections are therefore applied as a modification of WNoTFix, the W array without the
  ! d/dt correction term.

  ! Note that vertical velocity has already been set to zero at the surface and is skipped here.

  Do i = 1, HGrid%nX
  Do j = 1, HGrid%nY

    Do k = 2, ZGridW%nZ

      ! Calculate d Z(AboveGround) / d t on constant (X(m), Y(m), Z(w)) surface, with dX(m) and dY(m)
      ! denoting dX and dY in metres and Z(w) being the height in the coord system for which the
      ! input vertical velocity is Z(w) Dot. When Z(w) = Eta, the value can be computed directly;
      ! otherwise it is constructed as the product
      !
      !  d Z(AboveGround) / d (X(m), Y(m), Eta, t) * d (X(m), Y(m), Eta, t) / (X(m), Y(m), Z(w), t)
      !  d (X(m), Y(m), Eta, t) /  d t
      !
      ! where Eta is the vertical coord system of the main grids (i.e. model levels).
      ! Note that, in the second term, the derivatives are calculated at constant (X(m), Y(m), Z(w))
      ! and only the Eta and t derivatives are non-zero.

      ! Calculate vertical interpolation coefficients for calculating d Z(AboveGround) / d Eta.
      ! Note that we use values on ZGrid either side of ZGridW levels.
      ! $$ Asymmetric differences if W grid = main grid
      dZCoeffsToWdZ    =  ZCoeffsToW(k)
      dZCoeffsToWdZ%Z1 =  1.0
      dZCoeffsToWdZ%Z2 = -1.0

      ! Calculate d Z(AboveGround) / d Eta at current and previous met update times.
      dZdEtaNew = GInterpXYZ(M%Z(:, :, :, M%New), i, j, GHCoeffs,   dZCoeffsToWdZ) &
                  / (ZGrid%Z(ZCoeffsToW(k)%iZ2) - ZGrid%Z(ZCoeffsToW(k)%iZ1))
      dZdEtaOld = GInterpXYZ(M%Z(:, :, :, M%Old), i, j, GHCoeffs,   dZCoeffsToWdZ) &
                  / (ZGrid%Z(ZCoeffsToW(k)%iZ2) - ZGrid%Z(ZCoeffsToW(k)%iZ1))
      ! Calculate d Z(AboveGround) / d t (at constant Eta).
      dZdt = (GInterpXYZ(M%Z(:, :, :, M%New), i, j, GHCoeffs, ZCoeffsToW(k)) -  &
              GInterpXYZ(M%Z(:, :, :, M%Old), i, j, GHCoeffs, ZCoeffsToW(k)))   &
             / ShortTime2RealTime(M%Dt)

      ! Calculate d Z(AboveGround) / d t on constant (X(m), Y(m), Z(w)) surface.
      ! Note that dZdt is precisely this term already when Z(w) = Eta.
      If (M%iZCoordW /= M%iZCoord) Then
        ! Calculate surface pressure and surface pressure tendency.
        ! $$ Independent of height - move outside k loop?
        PStar       =  GInterpXY(M%P(:, :, 1, M%New), i, j, GHCoeffs)
        dPStardX(1) = (GInterpXY(M%P(:, :, 1, M%New), i, j, GHCoeffs) - &
                       GInterpXY(M%P(:, :, 1, M%Old), i, j, GHCoeffs))  &
                      / ShortTime2RealTime(M%Dt)
        dPStardX(2) = 0.0 ! not needed

        ! Calculate transformation coefficient d Eta / dt (at constant Z(w)).
        ! Note: using same routine as for the horizontal gradients but with the
        ! substitution of 'x' by 't'. Only use dEtadX(1) component.
        Call CalcdZdXPBased(ZCoord, ZCoordW, ZGridW%Z(k), PStar, dPStardX, dEtadX)

        ! Calculate correction terms.
        WCorrNew = dZdEtaNew * dEtadX(1) + dZdt
        WCorrOld = dZdEtaOld * dEtadX(1) + dZdt
      Else
        ! Calculate correction terms.
        WCorrNew = dZdt
        WCorrOld = dZdt
      End If

      ! Update W arrays at current and previous met update times.
      M%W(i, j, k, M%New) = M%WNoTFix(i, j, k, M%New) + WCorrNew
      M%W(i, j, k, M%Old) = M%WNoTFix(i, j, k, M%Old) + WCorrOld

    End Do

  End Do
  End Do

End Subroutine CalcZDotCorrection

!-------------------------------------------------------------------------------------------------------------

Function ConvertZColumn(Coords, Grids, M, i, j, ZCoordIn, ZCoordOut, ZIn, NewIdx) Result(ZOut)
! Converts vertical coords between coord systems in a column on the main horizontal
! grid.

  Implicit None
  ! Argument list:
  Type(Coords_), Intent(In),  Target :: Coords    ! Collection of coord systems.
  Type(Grids_),  Intent(In),  Target :: Grids     ! Collection of grids.
  Type(NWPMet_), Intent(In)          :: M         ! State of an NWP met module
                                                  ! instance.
  Integer,       Intent(In)          :: i         !} Grid indices on the main
  Integer,       Intent(In)          :: j         !} horizontal grid.
  Type(ZCoord_), Intent(In)          :: ZCoordIn  ! The coord system to convert from.
  Type(ZCoord_), Intent(In)          :: ZCoordOut ! The coord system to convert to.
  Real(Std),     Intent(In)          :: ZIn       ! Height in coord system ZCoordIn.
  Integer,       Intent(In)          :: NewIdx
  Integer                            :: MNew
  ! Function result:
  Real(Std) :: ZOut ! Height in coord system ZCoordOut.
  ! Locals:
  Type(ZCoord_), Pointer :: ZCoordZ !} Abbreviations for coords.
  Type(ZCoord_), Pointer :: ZCoordP !}
  Real(Std)              :: ZZ      !} Height in coord systems ZCoordZ and ZCoordP.
  Real(Std)              :: ZP      !}

  MNew = NewIdx

  ! Set up abbreviations for coords.
  ZCoordZ => Coords%ZCoords(M%iZCoordZ)
  ZCoordP => Coords%ZCoords(M%iZCoordP)

  Select Case (ZCoordIn%CoordType)

    Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)

      Select Case (ZCoordOut%CoordType)

        Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)

          ZOut = ZBasedToZBased(ZCoordIn, ZCoordOut, M%Topog(i, j, MNew), ZIn)

        Case (Z_P, Z_PAsZ, Z_PAsEta)

          ZZ   = ZBasedToZBased(ZCoordIn, ZCoordZ, M%Topog(i, j, MNew), ZIn)
          ZP   = ZToPColumn(Grids, M, i, j, ZZ, MNew)
          ZOut = PBasedToPBased(ZCoordP, ZCoordOut, M%P(i, j, 1, MNew), ZP)

        Case Default

          Call Message('UNEXPECTED ERROR in ConvertZColumn: unknown coordinate system type', 4)

      End Select

    Case (Z_P, Z_PAsZ, Z_PAsEta)

      Select Case (ZCoordOut%CoordType)

        Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)

          ZP   = PBasedToPBased(ZCoordIn, ZCoordP, M%P(i, j, 1, MNew), ZIn)
          ZZ   = PToZColumn(Grids, M, i, j, ZP, MNew)
          ZOut = ZBasedToZBased(ZCoordZ, ZCoordOut, M%Topog(i, j, MNew), ZZ)

        Case (Z_P, Z_PAsZ, Z_PAsEta)

          ZOut = PBasedToPBased(ZCoordIn, ZCoordOut, M%P(i, j, 1, MNew), ZIn)

        Case Default

          Call Message('UNEXPECTED ERROR in ConvertZColumn: unknown coordinate system type', 4)

      End Select

    Case Default

      Call Message('UNEXPECTED ERROR in ConvertZColumn: unknown coordinate system type', 4)

  End Select

End Function ConvertZColumn

!-------------------------------------------------------------------------------------------------------------

Function PToZColumn(Grids, M, i, j, PIn, NewIdx) Result(ZOut)
! Converts pressure to height above ground within a column of data on the main
! horizontal grid.

  Implicit None
  ! Argument List:
  Type(Grids_),  Intent(In), Target :: Grids ! Collection of grids.
  Type(NWPMet_), Intent(In)         :: M     ! State of an NWP met module instance.
  Integer,       Intent(In)         :: i     !} Grid indices on the main horizontal
  Integer,       Intent(In)         :: j     !} grid.
  Real(Std),     Intent(In)         :: PIn   ! Pressure.
  Integer,       Intent(In)         :: NewIdx
  Integer                           :: MNew
  ! Function result:
  Real(Std) :: ZOut ! Height above ground.
  ! Locals:
  Integer               :: kMin  !} Indices of levels bracketing location.
  Integer               :: kMax  !}
  Integer               :: k     ! Index of a level between kMin and kMax.
  Type(ZGrid_), Pointer :: ZGrid ! Abbreviation for grid.

  MNew = NewIdx

  ! Set up abbreviation for grid.
  ZGrid => Grids%ZGrids(M%iZGrid)

  ! Find levels surrounding point.
  kMax = ZGrid%nZ + 1
  kMin = 0
  Do
    k = kMax - kMin
    If (k <= 1) Exit
    k = kMin + k/2
    If ((PIn - M%P(i, j, k, MNew)) * ZGrid%IncreasingHeight < 0.0) Then
      kMin = k
    Else
      kMax = k
    End If
  End Do

  ! Interpolate. Above and below ZGrid we assume constant temperature.
  If (kMin == ZGrid%nZ) Then
    ZOut = M%Z(i, j, kMin, MNew) +                       &
           Log(M%P(i, j, kMin, MNew) / PIn) *            &
           GasConstant * M%T(i, j, kMin, MNew) / Gravity
  Else If (kMax == 1) Then
    ZOut = M%Z(i, j, kMax, MNew) +                       &
           Log(M%P(i, j, kMax, MNew) / PIn) *            &
           GasConstant * M%T(i, j, kMax, MNew) / Gravity
  Else
    ZOut =  M%Z(i, j, kMin, MNew) +                           &
           (M%Z(i, j, kMax, MNew) - M%Z(i, j, kMin, MNew)) *  &
           (PIn                    - M%P(i, j, kMin, MNew)) / &
           (M%P(i, j, kMax, MNew) - M%P(i, j, kMin, MNew))
  End If

End Function PToZColumn

!-------------------------------------------------------------------------------------------------------------

Function ZToPColumn(Grids, M, i, j, ZIn, NewIdx) Result(POut)
! Converts height above ground to pressure within a column of data on the main
! horizontal grid.

  Implicit None
  ! Argument List:
  Type(Grids_),  Intent(In), Target :: Grids ! Collection of grids.
  Type(NWPMet_), Intent(In)         :: M     ! State of an NWP met module instance.
  Integer,       Intent(In)         :: i     !} Grid indices on the main horizontal
  Integer,       Intent(In)         :: j     !} grid.
  Real(Std),     Intent(In)         :: ZIn   ! Height above ground.
  Integer,       Intent(In)         :: NewIdx
  Integer                           :: MNew
  ! Function result:
  Real(Std) :: POut ! Pressure.
  ! Locals:
  Integer               :: kMin  !} Indices of levels bracketing location.
  Integer               :: kMax  !}
  Integer               :: k     ! Index of a level between kMin and kMax.
  Type(ZGrid_), Pointer :: ZGrid ! Abbreviation for grid.

  MNew = NewIdx

  ! Set up abbreviation for grid.
  ZGrid => Grids%ZGrids(M%iZGrid)

  ! Find levels surrounding point.
  kMax = ZGrid%nZ + 1
  kMin = 0
  Do
    k = kMax - kMin
    If (k <= 1) Exit
    k = kMin + k/2
    If ((ZIn - M%Z(i, j, k, MNew)) * ZGrid%IncreasingHeight > 0.0) Then
      kMin = k
    Else
      kMax = k
    End If
  End Do

  ! Interpolate. Above and below ZGrid we assume constant temperature.
  If (kMin == ZGrid%nZ) Then
    POut = M%P(i, j, kMin, MNew) *                       &
           Exp(                                          &
             - Gravity * (ZIn - M%Z(i, j, kMin, MNew)) / &
             (GasConstant * M%T(i, j, kMin, MNew))       &
           )
  Else If (kMax == 1) Then
    POut = M%P(i, j, kMax, MNew) *                       &
           Exp(                                          &
             - Gravity * (ZIn - M%Z(i, j, kMax, MNew)) / &
             (GasConstant * M%T(i, j, kMax, MNew))       &
           )
  Else
    POut =  M%P(i, j, kMin, MNew) +                           &
           (M%P(i, j, kMax, MNew) - M%P(i, j, kMin, MNew)) *  &
           (ZIn                    - M%Z(i, j, kMin, MNew)) / &
           (M%Z(i, j, kMax, MNew) - M%Z(i, j, kMin, MNew))
  End If

End Function ZToPColumn

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcCloudInfo(Coords, Grids, M, NewIdx)
! Calculates (i) TotalOrDynCloudBase and TotalOrDynCloudTop from TotalOrDynCloudWater,
! TotalOrDynCloudIce, TotalOrDynHCloud, TotalOrDynMCloud and TotalOrDynLCloud and 
! (ii) Cloud3d, Cloud from TotalOrDynCloudWater, TotalOrDynCloudIce, TotalOrDynHCloud,
!  TotalOrDynMCloud and TotalOrDynLCloud and, if TotalCloudFlag false, ConCloud, 
! ConCloudBase and ConCloudTop.

  Implicit None
  ! Argument list:
  Type(Coords_), Intent(In),   Target :: Coords ! Collection of coord systems.
  Type(Grids_),  Intent(In),   Target :: Grids  ! Collection of grids.
  Type(NWPMet_), Intent(InOut)        :: M      ! State of an NWP met module instance.
  Integer,       Intent(In)           :: NewIdx
  Integer                             :: MNew

  ! Local parameters:
  Real(Std), Parameter :: LCloudBotEta  =     0.990 ! Low    cloud bottom (eta, um4).
  Real(Std), Parameter :: LCloudTopEta  =     0.800 ! Low    cloud top    (eta, um4).
  Real(Std), Parameter :: MCloudTopEta  =     0.500 ! Medium cloud top    (eta, um4).
  Real(Std), Parameter :: HCloudTopEta  =     0.150 ! High   cloud top    (eta, um4).
  Real(Std), Parameter :: LCloudBotZSea =   111.0   ! Low    cloud bottom (zsea, um5).
  Real(Std), Parameter :: LCloudTopZSea =  1949.0   ! Low    cloud top    (zsea, um5).
  Real(Std), Parameter :: MCloudTopZSea =  5574.0   ! Medium cloud top    (zsea, um5).
  Real(Std), Parameter :: HCloudTopZSea = 13608.0   ! High   cloud top    (zsea, um5).
  ! $$ these should all be inputs ?? coord systems always ZCoord ??
  ! Locals:
  Real(Std)              :: LCloudBot ! Low    cloud bottom in coord system ZCoord.
  Real(Std)              :: LCloudTop ! Low    cloud top    in coord system ZCoord.
  Real(Std)              :: MCloudTop ! Medium cloud top    in coord system ZCoord.
  Real(Std)              :: HCloudTop ! High   cloud top    in coord system ZCoord.
  Integer                :: i         !} Loop indices.
  Integer                :: j         !}
  Integer                :: k         !}
  Type(ZCoord_), Pointer :: ZCoord    !] Abbreviations for coords and grids.
  Type(HGrid_),  Pointer :: HGrid     !]
  Type(ZGrid_),  Pointer :: ZGrid     !]

  MNew = NewIdx

  ! Set up abbreviations for coords and grids.
  ZCoord => Coords%ZCoords(M%iZCoord)
  HGrid  => Grids%HGrids(M%iHGrid)
  ZGrid  => Grids%ZGrids(M%iZGrid)

  ! LCloudBot, LCloudTop, MCloudTop and HCloudTop. $$
  If (ZCoord%CoordType == Z_PAsEta) Then ! UM 4
    LCloudBot = LCloudBotEta
    LCloudTop = LCloudTopEta
    MCloudTop = MCloudTopEta
    HCloudTop = HCloudTopEta
  Else If (ZCoord%CoordType == Z_ZAsEta) Then ! UM 5
    LCloudBot = LCloudBotZSea / ZCoord%ModelTopHeight
    LCloudTop = LCloudTopZSea / ZCoord%ModelTopHeight
    MCloudTop = MCloudTopZSea / ZCoord%ModelTopHeight
    HCloudTop = HCloudTopZSea / ZCoord%ModelTopHeight
  Else
    Call Message('Cloud fraction routine requires z- or p-eta coord', 3)
  End If

  ! Cloud (total if M%TotalCloudFlag set to true or temporarily dynamic if M%TotalCloudFlag is false)
  ! $$ Note that if M%TotalCloudFlag is true then M%Cloud temporarily only includes cloud up to HCloudTop. 
  M%Cloud(:, :, MNew) = Max(                               & ! OK? $$
                           M%TotalOrDynLCloud(:, :, MNew), &
                           M%TotalOrDynMCloud(:, :, MNew), &
                           M%TotalOrDynHCloud(:, :, MNew)  &
                         )

  ! TotalOrDynCloudBase and TotalOrDynCloudTop. Set to < 0 if no cloud. Otherwise
  ! set to extreme values to prepare for calculating them.
  Do j = 1, HGrid%nY
  Do i = 1, HGrid%nX
    If (M%Cloud(i, j, MNew) > 0.0) Then
      M%TotalOrDynCloudBase(i, j, MNew) = Huge(1.0_Std)
      M%TotalOrDynCloudTop (i, j, MNew) = 0.0
    Else
      M%TotalOrDynCloudBase(i, j, MNew) = -1.0
      M%TotalOrDynCloudTop (i, j, MNew) = -1.0
    End If
  End Do
  End Do

  ! Cloud3d (total if M%TotalCloudFlag set to true or temporarily dynamic if M%TotalCloudFlag is false), 
  ! TotalOrDynCloudBase and TotalOrDynCloudTop.
  If (ZCoord%IncreasingUpwards) Then

    Do k = 1, ZGrid%nZ
    Do j = 1, HGrid%nY
    Do i = 1, HGrid%nX

      If (                                                      &
        M%TotalCloudFlag                                  .and. & ! For total cloud above HCloudTop
        ZGrid%Z(k) >=  HCloudTop                          .and. & ! < <= > >= etc $$
        (M%TotalOrDynCloudWater(i, j, k, MNew) > 0.0      .or.  &
         M%TotalOrDynCloudIce  (i, j, k, MNew) > 0.0)     .and. &
         M%ConCloudBase        (i, j, MNew) >= 0.0        .and. &
         M%ConCloudTop         (i, j, MNew) >= 0.0        .and. &
         M%Z(i, j, k, MNew) >= M%ConCloudBase(i, j, MNew) .and. & 
         M%Z(i, j, k, MNew) <= M%ConCloudTop (i, j, MNew)       &
      ) Then

        M%Cloud3d(i, j, k, MNew)          = M%ConCloud(i, j, MNew)
        M%TotalOrDynCloudBase(i, j, MNew) = Min(                                  &
                                               M%TotalOrDynCloudBase(i, j, MNew), &
                                               M%ConCloudBase(i, j, MNew)         &
                                             )
        M%TotalOrDynCloudTop (i, j, MNew) = Max(                                 &
                                               M%TotalOrDynCloudTop(i, j, MNew), &
                                               M%ConCloudTop(i, j, MNew)         &
                                             )
        M%Cloud(i, j, MNew) = Max(M%Cloud(i, j, MNew), M%ConCloud(i, j, MNew))

      Else If (                                            &
        ZGrid%Z(k) <  HCloudTop                      .and. & 
        ZGrid%Z(k) >= MCloudTop                      .and. &
        M%TotalOrDynHCloud(i, j, MNew) > 0.0         .and. &
        (M%TotalOrDynCloudWater(i, j, k, MNew) > 0.0 .or.  &
         M%TotalOrDynCloudIce  (i, j, k, MNew) > 0.0)      &
      ) Then

        M%Cloud3d(i, j, k, MNew)          = M%TotalOrDynHCloud(i, j, MNew)
        M%TotalOrDynCloudBase(i, j, MNew) = Min(                                  &
                                               M%TotalOrDynCloudBase(i, j, MNew), &
                                               M%Z(i, j, k, MNew)                 &
                                             )
        M%TotalOrDynCloudTop (i, j, MNew) = Max(                                 &
                                               M%TotalOrDynCloudTop(i, j, MNew), &
                                               M%Z(i, j, k, MNew)                &
                                             )

      Else If (                                            &
        ZGrid%Z(k) <  MCloudTop                      .and. &
        ZGrid%Z(k) >= LCloudTop                      .and. &
        M%TotalOrDynMCloud(i, j, MNew) > 0.0         .and. &
        (M%TotalOrDynCloudWater(i, j, k, MNew) > 0.0 .or.  &
         M%TotalOrDynCloudIce  (i, j, k, MNew) > 0.0)      &
      ) Then

        M%Cloud3d(i, j, k, MNew)          = M%TotalOrDynMCloud(i, j, MNew)
        M%TotalOrDynCloudBase(i, j, MNew) = Min(                                  &
                                               M%TotalOrDynCloudBase(i, j, MNew), &
                                               M%Z(i, j, k, MNew)                 &
                                             )
        M%TotalOrDynCloudTop (i, j, MNew) = Max(                                 &
                                               M%TotalOrDynCloudTop(i, j, MNew), &
                                               M%Z(i, j, k, MNew)                &
                                             )

      Else If (                                            &
        ZGrid%Z(k) <  LCloudTop                      .and. &
        ZGrid%Z(k) >= LCloudBot                      .and. &
        M%TotalOrDynLCloud(i, j, MNew) > 0.0         .and. &
        (M%TotalOrDynCloudWater(i, j, k, MNew) > 0.0 .or.  &
         M%TotalOrDynCloudIce  (i, j, k, MNew) > 0.0)      &
      ) Then

        M%Cloud3d(i, j, k, MNew)          = M%TotalOrDynLCloud(i, j, MNew)
        M%TotalOrDynCloudBase(i, j, MNew) = Min(                                  &
                                               M%TotalOrDynCloudBase(i, j, MNew), &
                                               M%Z(i, j, k, MNew)                 &
                                             )
        M%TotalOrDynCloudTop (i, j, MNew) = Max(                                 &
                                               M%TotalOrDynCloudTop(i, j, MNew), &
                                               M%Z(i, j, k, MNew)         &
                                             )

      Else

        M%Cloud3d(i, j, k, MNew)          = 0.0

      End If

    End Do
    End Do
    End Do

  Else

    Do k = 1, ZGrid%nZ
    Do j = 1, HGrid%nY
    Do i = 1, HGrid%nX

      If (                                                      &
        M%TotalCloudFlag                                  .and. & ! For total cloud above HCloudTop
        ZGrid%Z(k) <  HCloudTop                           .and. & ! < <= > >= etc $$
        (M%TotalOrDynCloudWater(i, j, k, MNew) > 0.0      .or.  &
         M%TotalOrDynCloudIce  (i, j, k, MNew) > 0.0)     .and. &
         M%ConCloudBase        (i, j, MNew) >= 0.0        .and. &
         M%ConCloudTop         (i, j, MNew) >= 0.0        .and. &
         M%Z(i, j, k, MNew) >= M%ConCloudBase(i, j, MNew) .and. & 
         M%Z(i, j, k, MNew) <= M%ConCloudTop (i, j, MNew)       &
      ) Then

        M%Cloud3d(i, j, k, MNew)          = M%ConCloud(i, j, MNew)
        M%TotalOrDynCloudBase(i, j, MNew) = Min(                                  &
                                               M%TotalOrDynCloudBase(i, j, MNew), &
                                               M%ConCloudBase(i, j, MNew)         &
                                             )
        M%TotalOrDynCloudTop (i, j, MNew) = Max(                                 &
                                               M%TotalOrDynCloudTop(i, j, MNew), &
                                               M%ConCloudTop(i, j, MNew)         &
                                             )
        M%Cloud(i, j, MNew) = Max(M%Cloud(i, j, MNew), M%ConCloud(i, j, MNew))
        
      Else If (                                                 &
        ZGrid%Z(k) >= HCloudTop                      .and. &
        ZGrid%Z(k) <  MCloudTop                      .and. &
        M%TotalOrDynHCloud(i, j, MNew) > 0.0         .and. &
        (M%TotalOrDynCloudWater(i, j, k, MNew) > 0.0 .or.  &
         M%TotalOrDynCloudIce  (i, j, k, MNew) > 0.0)      &
      ) Then

        M%Cloud3d(i, j, k, MNew)          = M%TotalOrDynHCloud(i, j, MNew)
        M%TotalOrDynCloudBase(i, j, MNew) = Min(                                  &
                                               M%TotalOrDynCloudBase(i, j, MNew), &
                                               M%Z(i, j, k, MNew)                 &
                                             )
        M%TotalOrDynCloudTop (i, j, MNew) = Max(                                 &
                                               M%TotalOrDynCloudTop(i, j, MNew), &
                                               M%Z(i, j, k, MNew)                &
                                             )

      Else If (                                            &
        ZGrid%Z(k) >= MCloudTop                      .and. &
        ZGrid%Z(k) <  LCloudTop                      .and. &
        M%TotalOrDynMCloud(i, j, MNew) > 0.0         .and. &
        (M%TotalOrDynCloudWater(i, j, k, MNew) > 0.0 .or.  &
         M%TotalOrDynCloudIce  (i, j, k, MNew) > 0.0)      &
      ) Then

        M%Cloud3d(i, j, k, MNew)          = M%TotalOrDynMCloud(i, j, MNew)
        M%TotalOrDynCloudBase(i, j, MNew) = Min(                                  &
                                               M%TotalOrDynCloudBase(i, j, MNew), &
                                               M%Z(i, j, k, MNew)                 &
                                             )
        M%TotalOrDynCloudTop (i, j, MNew) = Max(                                 &
                                               M%TotalOrDynCloudTop(i, j, MNew), &
                                               M%Z(i, j, k, MNew)                &
                                             )

      Else If (                                            &
        ZGrid%Z(k) >= LCloudTop                      .and. &
        ZGrid%Z(k) <  LCloudBot                      .and. &
        M%TotalOrDynLCloud(i, j, MNew) > 0.0         .and. &
        (M%TotalOrDynCloudWater(i, j, k, MNew) > 0.0 .or.  &
         M%TotalOrDynCloudIce  (i, j, k, MNew) > 0.0)      &
      ) Then

        M%Cloud3d(i, j, k, MNew)          = M%TotalOrDynLCloud(i, j, MNew)
        M%TotalOrDynCloudBase(i, j, MNew) = Min(                                  &
                                               M%TotalOrDynCloudBase(i, j, MNew), &
                                               M%Z(i, j, k, MNew)                 &
                                             )
        M%TotalOrDynCloudTop (i, j, MNew) = Max(                                 &
                                               M%TotalOrDynCloudTop(i, j, MNew), &
                                               M%Z(i, j, k, MNew)                &
                                             )

      Else

        M%Cloud3d(i, j, k, MNew)          = 0.0

      End If

    End Do
    End Do
    End Do

  End If

  ! Add in convective cloud component to Cloud3d.and Cloud if M%TotalCloudFlag is false 
  Do j = 1, HGrid%nY
  Do i = 1, HGrid%nX

    If (                                         &
       (M%ConCloudBase(i, j, MNew) <  0.0  .and. &
        M%ConCloudTop (i, j, MNew) >= 0.0) .or.  &
       (M%ConCloudTop (i, j, MNew) <  0.0  .and. &
        M%ConCloudBase(i, j, MNew) >= 0.0)       &
       ) Then
      Call Message('WARNING: Convective cloud base and top inconsistent', 1)
    End If

  End Do
  End Do

  If (.not. M%TotalCloudFlag) Then

    M%Cloud(:, :, MNew) = 1.0 - (1.0 - M%ConCloud(:, :, MNew)) * (1.0 - M%Cloud(:, :, MNew))


    Do k = 1, ZGrid%nZ
    Do j = 1, HGrid%nY
    Do i = 1, HGrid%nX

      If (                                                      &
         M%ConCloudBase(i, j, MNew) >= 0.0                .and. &
         M%ConCloudTop (i, j, MNew) >= 0.0                .and. &
         M%Z(i, j, k, MNew) >= M%ConCloudBase(i, j, MNew) .and. & ! < <= > >= etc $$
         M%Z(i, j, k, MNew) <= M%ConCloudTop (i, j, MNew)       &
         ) Then

        M%Cloud3d(i, j, k, MNew) = 1.0 - (1.0 - M%ConCloud(i, j, MNew)) * (1.0 - M%Cloud3d(i, j, k, MNew))

      End If

    End Do
    End Do
    End Do
  
  End If
  
  
End Subroutine CalcCloudInfo

!-------------------------------------------------------------------------------------------------------------

End Module NWPMetModule

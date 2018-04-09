! Module:  Ancillary Met Module

Module AncillaryMetModule

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
! It only works in a calendar time frame. $$ add test for this
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

# ifdef GRIBsupport
  ! If support for GRIB is requested, the code references external routines in the GRIB_API module, which are
  ! linked dynamically at run time from their shared object libraries. The GRIB_API is an ECMWF package for
  ! encoding/decoding gridded meteorological fields in GRIB-1 and GRIB-2 formats.
  Use grib_api
# endif

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: MA_AncillarySoil             !} Codes for the module attributes. These indicate that the module can
Public :: MA_AncillaryLandUse          !} supply soil, land use, soil moisture, leaf area index (LAI) and
Public :: MA_AncillarySoilMoisture     !} canopy height information. Exactly what is needed/supplied for
Public :: MA_AncillaryLAI              !} each attribute is noted in the Needed array.
Public :: MA_AncillaryCanopyHeight     !}
Public :: AncillaryMetDefn_            ! An NWP met definition (part 1 - basic definitions).
Public :: AncillaryMetDefn2_           ! An NWP met definition (part 2 - met file
                                       ! structure definition).
Public :: AncillaryMetDefns_           ! A collection of met definitions.
Public :: AncillaryMet_                ! The state of an NWP met module instance.
Public :: InitAncillaryMetDefns        ! Initialises a collection of met definitions.
Public :: InitAncillaryMetDefn         ! Initialises an NWP met definition (part 1).
Public :: AddAncillaryMetDefn          ! Adds an NWP met definition (part 1) to a collection
                                       ! of met definitions.
Public :: FindAncillaryMetDefnIndex    ! Finds the index of an NWP met definition (part 1).
Public :: InitAncillaryMetDefn2        ! Initialises an NWP met definition (part 2).
Public :: AddAncillaryMetDefn2         ! Adds an NWP met definition (part 2) to a collection
                                       ! of met definitions.
Public :: FindAncillaryMetDefn2Index   ! Finds the index of an NWP met definition (part 2).
Public :: InitAncillaryMet             ! Initialises the state of an NWP met module
                                       ! instance.
Public :: SetUpCoordsEtc_AncillaryMet  ! Sets up Coords and Grids by adding any extra coords
                                       ! and grids which NWPMet wants to define.
Public :: SetUpAncillaryMet_CoordsEtc  ! Sets up NWPMet using information from EtaDefns,
                                       ! Coords and Grids.
Public :: PrepareForUpdateAncillaryMet ! Prepares for updating an instance of the ancillary met module.
Public :: UpdateAncillaryMet           ! Updates an instance of the NWP met module.
Public :: ResetAncillaryMet            ! Resets the state of an NWP met module instance for
                                       ! a new realisation.

!-------------------------------------------------------------------------------------------------------------

! Codes for the module attributes.
Integer, Parameter :: MA_AncillarySoil         = 1 !} Codes for the module attributes. These indicate that the
Integer, Parameter :: MA_AncillaryLandUse      = 2 !} module can supply soil, land use, soil moisture, leaf
Integer, Parameter :: MA_AncillarySoilMoisture = 3 !} area index (LAI) and canopy height information. Exactly
Integer, Parameter :: MA_AncillaryLAI          = 4 !} what is needed/supplied for each attribute is noted in
Integer, Parameter :: MA_AncillaryCanopyHeight = 5 !} the Needed array.


! NWP met field names and information on which fields can be missing:
Integer,                  Parameter :: nFieldNames = 10
Character(MaxCharLength), Parameter :: FieldNames(nFieldNames) =                       &
                                       (/                                              &
                                         'land use fractions (9 land use types)     ', & ! 1
                                         'clay mass fraction                        ', & ! 2
                                         'silt mass fraction                        ', & ! 3
                                         'sand mass fraction                        ', & ! 4
                                         'soil particle mass fractions (6 size bins)', & ! 5
                                         'soil moisture in layer (kg/m^2)           ', & ! 6
                                         'leaf area index (5 plant land use types)  ', & ! 7
                                         'canopy height (m) (5 plant types)         ', & ! 8
                                         'land fraction                             ', & ! 9
                                         'dummy                                     '  & ! 10
                                       /)
Integer,                  Parameter :: FieldDimensions(nFieldNames) = & ! $$ eventually input these
                                       (/                             & ! $$ currently 0 = 2-d
                                         9,                           & ! 1
                                         0,                           & ! 2
                                         0,                           & ! 3
                                         0,                           & ! 4
                                         6,                           & ! 5
                                         0,                           & ! 6
                                         5,                           & ! 7
                                         5,                           & ! 8
                                         0,                           & ! 9
                                         0                            & ! 10
                                       /)
Integer,                  Parameter :: PersistTimes(nFieldNames) = &
                                       (/                          &
                                         0,                        & ! 1
                                         0,                        & ! 2
                                         0,                        & ! 3
                                         0,                        & ! 4
                                         0,                        & ! 5
                                         0,                        & ! 6
                                         0,                        & ! 7
                                         0,                        & ! 8
                                         0,                        & ! 9
                                         0                         & ! 10
                                       /)
Logical,                  Parameter :: AllowFixUpToDefault(nFieldNames) = &
                                       (/                                 &
                                         .false.,                         & ! 1
                                         .false.,                         & ! 2
                                         .false.,                         & ! 3
                                         .false.,                         & ! 4
                                         .false.,                         & ! 5
                                         .false.,                         & ! 6
                                         .false.,                         & ! 7
                                         .false.,                         & ! 8
                                         .false.,                         & ! 9
                                         .true.                           & ! 10
                                       /)
Integer,                  Parameter :: nAttribs = 5
Logical,                  Parameter :: Needed(nAttribs, nFieldNames) =    &
                                       Reshape(                           &
                                         (/                               &
!                                          Soil     LandUse  SoilMoisture  LAI       CanopyHeight
                                           .false., .true.,  .false.,      .false.,  .false.,     & ! 1
                                           .true.,  .false., .false.,      .false.,  .false.,     & ! 2
                                           .true.,  .false., .false.,      .false.,  .false.,     & ! 3
                                           .true.,  .false., .false.,      .false.,  .false.,     & ! 4
                                           .true.,  .false., .false.,      .false.,  .false.,     & ! 5
                                           .false., .false., .true.,       .false.,  .false.,     & ! 6
                                           .false., .false., .false.,      .true.,   .false.,     & ! 7
                                           .false., .false., .false.,      .false.,  .true.,      & ! 8
                                           .false., .true.,  .false.,      .false.,  .false.,     & ! 9
                                           .false., .false., .false.,      .false.,  .false.      & ! 10
                                         /),                                                      &
                                         (/ 5, 10 /)                                              &
                                       )
! nFieldNames2d            :: Number of 2-d (horizontal) fields.
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
! data, provided ExtrapolationLevelsBelow is at least (1, 1, 1, 1, 2, 0, 2, 2) and provided
! AllowFixUpToDefault is true for roughness length and boundary layer depth. In practice these files have some
! missing data and what is needed to cope with missing data is not clear.

! Codes for NWP met fields:
Integer, Parameter :: F_LandUseFracs  = 1 !} Codes giving the location in the array
Integer, Parameter :: F_ClayMassFrac  = 2 !} FieldNames of the
Integer, Parameter :: F_SiltMassFrac  = 3 !} fields.
Integer, Parameter :: F_SandMassFrac  = 4 !}
Integer, Parameter :: F_SoilMassFracs = 5 !}
Integer, Parameter :: F_SoilMoisture  = 6 !}
Integer, Parameter :: F_LAI           = 7 !}
Integer, Parameter :: F_CanopyHeight  = 8 !}
Integer, Parameter :: F_LandFrac      = 9 !}

!-------------------------------------------------------------------------------------------------------------

Type :: AncillaryMetDefn_ ! An NWP met definition (part 1 - basic definitions). Together
                    ! with part 2 this provides information on the nature of the NWP
                    ! met module met data.
  Character(MaxCharLength)     :: Name
  Character(MaxCharLength)     :: BinaryFormat
  Character(MaxCharLength)     :: FileType
  Character(MaxCharLength)     :: Prefix
  Character(MaxCharLength)     :: Suffix
  Type(Time_)                  :: T0
  Type(Time_)                  :: Dt
  Logical                      :: Monthly
  Logical                      :: NextHeatFlux
  Logical                      :: NextPrecip
  Character(MaxCharLength)     :: MetDefn2Name
  Character(MaxCharLength)     :: HGridName
  ! Name         :: Name of met definition (part 1 - basic definitions).
  ! BinaryFormat :: Binary format of met and topography files.
  ! FileType     :: Type of met and topography files (Name II, PP or GRIB).
  ! Prefix       :: Prefix for names of met files.
  ! Suffix       :: Suffix for names of met files.
  ! T0           :: Reference time for the first met field.
  ! Dt           :: Time interval between fields.
  ! NextHeatFlux :: Indicates its best to use the next time step for heat flux and ustar rather than
  !                 interpolating in time.
  ! NextPrecip   :: Indicates its best to use the next time step for precipitation rather than interpolating
  !                 in time.
  ! MetDefn2Name :: Name of the NWP met definition (part 2 - met file structure definition).
  ! HGridName    :: Name of main horizontal grid.
End Type AncillaryMetDefn_

!-------------------------------------------------------------------------------------------------------------

Type :: AncillaryMetDefn2_ ! An NWP met definition (part 2 - met file structure definition).
                     ! Together with part 1 this provides information on the nature of
                     ! the NWP met module met data.
  Character(MaxCharLength) :: Name
  Integer                  :: nFields
  Character(MaxCharLength) :: FieldNames(MaxAncillaryMetFields)
  Integer                  :: LowestLevels(MaxAncillaryMetFields)
  Integer                  :: HighestLevels(MaxAncillaryMetFields)
  Logical                  :: Top(MaxAncillaryMetFields)
  Integer                  :: FieldCodes(MaxAncillaryMetFields)
  Logical                  :: ThreeD(MaxAncillaryMetFields)
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
  ! $$ for 'Sequential read in same order as NWPMetDefn2_' could possibly check any
  ! duplicate field is the same.

End Type AncillaryMetDefn2_

!-------------------------------------------------------------------------------------------------------------

Type :: AncillaryMetDefns_ ! A collection of met definitions.
  Integer                  :: nNWPMetDefns                  ! Number of NWP met definitions (part 1 - basic
                                                            ! definitions).
  Type(AncillaryMetDefn_)  :: NWPMetDefns(MaxNWPMetDefns)   ! NWP met definitions
                                                            ! (part 1 - basic definitions).
  Integer                  :: nNWPMetDefn2s                 ! Number of NWP met definitions (part 2 - met file
                                                            ! structure definition).
  Type(AncillaryMetDefn2_) :: NWPMetDefn2s(MaxNWPMetDefn2s) ! NWP met definitions (part 2 - met file structure
                                                            ! definition).
End Type AncillaryMetDefns_

!-------------------------------------------------------------------------------------------------------------

Type :: AncillaryMet_ ! The state of an NWP met module instance.

  ! CommonMet_:
  Type(CommonMet_) :: C ! Part of met state common to all flow modules.

  ! Validity flags for each time level.
  Logical :: ValidAttribs(nAttribs, 2) ! Indicates whether, for each time level, the module is valid for each
                                       ! of its attributes.

  ! Input variables:
  Type(AncillaryMetDefn_)              :: MetDefn           ! NWP met definition (part 1 - basic
                                                             ! definitions).
  Type(AncillaryMetDefn2_)             :: MetDefn2          ! NWP met definition (part 2 - met file structure
                                                             ! definition).
  Type(ShortTime_)                      :: Dt                ! Time interval between fields.
  Integer                               :: MetFolderDefnType ! Indicates the way the met folder(s) are defined
                                                             ! - i.e. which of MetFolder, MetFolderStem or
                                                             ! MetFolders is used.
                                                             !   0 = the single met folder is valid,
                                                             !   1 = the met folder stem is valid,
                                                             !   2 = the array of met folders is valid.
                                                             !   3 = the single ensemble met folder is valid.
  Character(MaxFileNameLength)          :: MetFolder         ! Folder containing met data.
  Character(MaxFileNameLength)          :: MetFolderStem     ! Stem of folders containing met data.
  Character(MaxFileNameLength), Pointer :: MetFolders(:)     ! Array of folders containing met data.
  Character(MaxFileNameLength)          :: EnsembleMetFolder ! Folder containing ensemble of met data.

  ! Quantities derived from MetDefn2 and Grids:
  Integer :: iField(MaxAncillaryMetFields) ! Indices in FieldNames of
                                           ! fields in the met file.
  Logical :: ThreeD(MaxAncillaryMetFields) ! Indicates if the NAME III field is 3-D
                                           ! for each of the fields in the
                                           ! met file.
  Integer :: k1(MaxAncillaryMetFields)     ! Lowest NAME III model level of each of the fields in the
                                           ! met file (1 for 2-D fields).
  Integer :: k2(MaxAncillaryMetFields)     ! Highest NAME III model level of each of the fields in
                                           ! the met file (1 for 2-D fields).
  Integer :: kMax(nFieldNames)             ! Highest NAME III model level of each of the fields listed in
                                           ! FieldNames.

  ! Grid and coord indices:
  Integer :: iHCoord ! Indices of the main horizontal coord system.
  Integer :: iHGrid  ! Grid indices.

  ! Allocatable arrays:
  Real(Std), Pointer :: LandUseFracs (:,:,:,:)
  Real(Std), Pointer :: ClayMassFrac (:,:,:,:)
  Real(Std), Pointer :: SiltMassFrac (:,:,:,:)
  Real(Std), Pointer :: SandMassFrac (:,:,:,:)
  Real(Std), Pointer :: SoilMassFracs(:,:,:,:)
  Real(Std), Pointer :: SoilMoisture (:,:,:,:)
  Real(Std), Pointer :: LAI          (:,:,:,:)
  Real(Std), Pointer :: CanopyHeight (:,:,:,:)
  Real(Std), Pointer :: LandFrac     (:,:,:,:)
  ! LandUseFracs  ::
  ! ClayMassFrac  ::
  ! SiltMassFrac  ::
  ! SandMassFrac  ::
  ! SoilMassFracs ::
  ! SoilMoisture  ::
  ! LAI           ::
  ! CanopyHeight  ::
  ! Note that LandUseFracs, ClayMassFrac, SiltMassFrac, SandMassFrac, SoilMassFracs, SoilMoisture, LAI
  ! and CanopyHeight are defined to be negative over the sea.

  ! LandFrac   ::
  ! Note that LandFrac is 0 over sea and 1 over land
  ! (and in future could be between 0 and 1 for coastal grid points)

  ! Information on allocatable arrays:
  Logical                  :: SpaceAllocated
  Integer                  :: New
  Integer                  :: Old
  Type(ShortTime_)         :: OldTime
  Type(ShortTime_)         :: NewTime
  Logical,         Pointer :: FieldPresent(:, :)
  Integer,         Pointer :: FieldPersisted(:, :)
  ! SpaceAllocated :: Indicates arrays allocated and topography read in.
  ! New            :} Indices of latest and latest but one sets of met data. 0
  ! Old            :} indicates data invalid because no data read in, errors occurred
  !                   during reading the data, or data is for the wrong time.
  ! OldTime        :] Time of latest and latest but one sets of met data.
  ! NewTime        :]
  ! FieldPresent   :: Indicates which levels of the 3-d fields have been successfully read in.
  ! FieldPersisted :: The number of met times for which persistence has been used for the
  ! levels of the 3-d fields. -1
  !                   indicates a fix up to a default value has been applied.
  ! Present2d      :: Indicates which 2-d fields have been successfully read in.
  ! Persist2d      :: The number of met times for which persistence has been used for the 2-d
  !                   fields. -1 indicates a fix up to a default value has been applied.

  ! $$ note persistence counters should go to restart file.

End Type AncillaryMet_

!-------------------------------------------------------------------------------------------------------------

Interface Operator(==) ! Equality of NWP met definitions.
  Module Procedure NWPMetDefnEq
  Module Procedure NWPMetDefn2Eq
End Interface

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitAncillaryMetDefns() Result(MetDefns)
! Initialises a collection of met definitions.

  Implicit None
  ! Function result:
  Type(AncillaryMetDefns_) :: MetDefns ! Initialised collection of met definitions.

  MetDefns%nNWPMetDefns  = 0
  MetDefns%nNWPMetDefn2s = 0

End Function InitAncillaryMetDefns

!-------------------------------------------------------------------------------------------------------------

Function InitAncillaryMetDefn(                           &
           Name, BinaryFormat, FileType, Prefix, Suffix, &
           T0, Dt, Monthly, NextHeatFlux, NextPrecip,    &
           MetDefn2Name,                                 &
           HGridName                                     &
         )                                               &
Result(NWPMetDefn)
! Initialises an NWP met definition (part 1 - basic definitions).

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name         ! Name of met definition (part 1 - basic definitions).
  Character(*), Intent(In) :: BinaryFormat ! Binary format of met and topography files.
  Character(*), Intent(In) :: FileType     ! Type of met and topography files (Name II, PP or GRIB).
  Character(*), Intent(In) :: Prefix       ! Prefix for names of met files.
  Character(*), Intent(In) :: Suffix       ! Suffix for names of met files.
  Type(Time_),  Intent(In) :: T0           ! Reference time for met fields.
  Type(Time_),  Intent(In) :: Dt           ! Time interval between fields.
  Logical,      Intent(In) :: Monthly      !
  Logical,      Intent(In) :: NextHeatFlux ! Indicates its best to use the next time step for heat flux and
                                           ! ustar rather than interpolating in time.
  Logical,      Intent(In) :: NextPrecip   ! Indicates its best to use the next time step for precipitation
                                           ! rather than interpolating in time.
  Character(*), Intent(In) :: MetDefn2Name ! Name of NWP met definition (part 2 - met
                                           ! file structure definition).
  Character(*), Intent(In) :: HGridName    ! Name of main horizontal grid.
  ! Function result:
  Type(AncillaryMetDefn_) :: NWPMetDefn ! The NWP met definition (part 1 - basic
                                        ! definitions).

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

  If (Len_Trim(Suffix) > MaxCharLength) Then
    Call Message(                                             &
           'ERROR in InitNWPMetDefn: Suffix is given as "' // &
           Trim(Suffix)                                    // &
           '" and is too long',                               &
           3                                                  &
         )
  End If

  If (Len_Trim(MetDefn2Name) == 0) Then
    Call Message('ERROR in InitNWPMetDefn: MetDefn2Name is blank', 3)
  End If
  If (Len_Trim(MetDefn2Name) > MaxCharLength) Then
    Call Message(                                                   &
           'ERROR in InitNWPMetDefn: MetDefn2Name is given as "' // &
           Trim(MetDefn2Name)                                    // &
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

  NWPMetDefn%Name         = Name
  NWPMetDefn%Prefix       = Prefix
  NWPMetDefn%Suffix       = Suffix
  NWPMetDefn%NextHeatFlux = NextHeatFlux
  NWPMetDefn%NextPrecip   = NextPrecip
  NWPMetDefn%MetDefn2Name = MetDefn2Name
  NWPMetDefn%HGridName    = HGridName

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
  NWPMetDefn%FileType = FileType

  ! T0.
  NWPMetDefn%T0 = T0

  ! Dt.
  NWPMetDefn%Monthly = Monthly
  If (Monthly) Then
    NWPMetDefn%Dt = ZeroTime()
  Else
    If (Dt <= ZeroTime()) Then
      Call Message('ERROR in InitNWPMetDefn: Dt is <= 0', 3)
    End If
    NWPMetDefn%Dt = Dt
  End If

End Function InitAncillaryMetDefn

!-------------------------------------------------------------------------------------------------------------

Subroutine AddAncillaryMetDefn(NWPMetDefn, MetDefns)
! Adds an NWP met definition (part 1 - basic definitions) to a collection of met
! definitions.

  Implicit None
  ! Argument list:
  Type(AncillaryMetDefn_), Intent(In)    :: NWPMetDefn ! The NWP met definition (part 1 -
                                                        ! basic definitions).
  Type(AncillaryMetDefns_),   Intent(InOut) :: MetDefns   ! The collection of met definitions.
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

End Subroutine AddAncillaryMetDefn

!-------------------------------------------------------------------------------------------------------------

Function NWPMetDefnEq(NWPMetDefn1, NWPMetDefn2)
! Tests for equality of NWP met definitions (part 1 - basic definitions).

  Implicit None
  ! Argument list:
  Type(AncillaryMetDefn_), Intent(In) :: NWPMetDefn1 !} The two NWP met definitions (part 1
  Type(AncillaryMetDefn_), Intent(In) :: NWPMetDefn2 !} - basic definitions).
  ! Function result:
  Logical :: NWPMetDefnEq ! Indicates if NWP met definitions are equal.

  NWPMetDefnEq = (NWPMetDefn1%Name         .CIEq. NWPMetDefn2%Name        ) .and. &
                 (NWPMetDefn1%BinaryFormat .CIEq. NWPMetDefn2%BinaryFormat) .and. &
                 (NWPMetDefn1%FileType     .CIEq. NWPMetDefn2%FileType    ) .and. &
                 (NWPMetDefn1%Prefix       .CIEq. NWPMetDefn2%Prefix      ) .and. &
                 (NWPMetDefn1%Suffix       .CIEq. NWPMetDefn2%Suffix      ) .and. &
                 (NWPMetDefn1%T0             ==   NWPMetDefn2%T0          ) .and. &
                 (NWPMetDefn1%Dt             ==   NWPMetDefn2%Dt          ) .and. &
                 (NWPMetDefn1%Monthly      .eqv.  NWPMetDefn2%Monthly     ) .and. &
                 (NWPMetDefn1%NextHeatFlux .eqv.  NWPMetDefn2%NextHeatFlux) .and. &
                 (NWPMetDefn1%NextPrecip   .eqv.  NWPMetDefn2%NextPrecip  ) .and. &
                 (NWPMetDefn1%MetDefn2Name .CIEq. NWPMetDefn2%MetDefn2Name) .and. &
                 (NWPMetDefn1%HGridName    .CIEq. NWPMetDefn2%HGridName   )

End Function NWPMetDefnEq

!-------------------------------------------------------------------------------------------------------------

Function FindAncillaryMetDefnIndex(Name, MetDefns)
! Finds the index of an NWP met definition (part 1 - basic definitions).

  Implicit None
  ! Argument list:
  Character(*),    Intent(In) :: Name     ! Name of NWP met definition (part 1 - basic
                                          ! definitions).
  Type(AncillaryMetDefns_), Intent(In) :: MetDefns ! Collection of met definitions.
  ! Function result:
  Integer :: FindAncillaryMetDefnIndex ! Index of NWP met definition (part 1 - basic
                                 ! definitions).
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, MetDefns%nNWPMetDefns
    If (Name .CIEq. MetDefns%NWPMetDefns(i)%Name) Then
      FindAncillaryMetDefnIndex = i
      Return
    End If
  End Do

  Call Message(                                 &
         'FATAL ERROR: NWP Met Definition "' // &
         Trim(Name)                          // &
         '" not found',                         &
         3                                      &
       )

End Function FindAncillaryMetDefnIndex

!-------------------------------------------------------------------------------------------------------------

Function InitAncillaryMetDefn2(                                                         &
           Name, FieldNames, LowestLevels, HighestLevels, Top, FieldCodes, ThreeD &
         )                                                                        &
Result(NWPMetDefn2)
! Initialises an NWP met definition (part 2 - met file structure definition).

  Implicit None
  ! Argument list:
  Character(*),              Intent(In) :: Name
  Character(MaxTokenLength), Intent(In) :: FieldNames(:)  ! $$ changed from * to support Intel compiler on Cray 
  Integer,                   Intent(In) :: LowestLevels(:)
  Integer,                   Intent(In) :: HighestLevels(:)
  Logical,                   Intent(In) :: Top(:)
  Integer,                   Intent(In) :: FieldCodes(:)
  Logical,                   Intent(In) :: ThreeD(:)
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
  ! Function result:
  Type(AncillaryMetDefn2_) :: NWPMetDefn2 ! The NWP met definition (part 2 - met file structure definition).
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
    Size(FieldNames) /= Size(ThreeD)             &
  ) Then
    Call Message(                                                           &
           'UNEXPECTED ERROR in initialising an NWP Met File Structure ' // &
           'Definition: the array sizes are inconsistent',                  &
           3                                                                &
         )
  End If
  If (Size(FieldNames) > MaxAncillaryMetFields) Then
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
  End Do

  NWPMetDefn2%Name                              = Name
  NWPMetDefn2%nFields                           = Size(FieldNames)
  NWPMetDefn2%FieldNames   (1:Size(FieldNames)) = FieldNames   (1:Size(FieldNames))
  NWPMetDefn2%LowestLevels (1:Size(FieldNames)) = LowestLevels (1:Size(FieldNames))
  NWPMetDefn2%HighestLevels(1:Size(FieldNames)) = HighestLevels(1:Size(FieldNames))
  NWPMetDefn2%Top          (1:Size(FieldNames)) = Top          (1:Size(FieldNames))
  NWPMetDefn2%FieldCodes   (1:Size(FieldNames)) = FieldCodes   (1:Size(FieldNames))
  NWPMetDefn2%ThreeD       (1:Size(FieldNames)) = ThreeD       (1:Size(FieldNames))

End Function InitAncillaryMetDefn2

!-------------------------------------------------------------------------------------------------------------

Subroutine AddAncillaryMetDefn2(NWPMetDefn2, MetDefns)
! Adds an NWP met definition (part 2 - met file structure definition) to a
! collection of met definitions.

  Implicit None
  ! Argument list:
  Type(AncillaryMetDefn2_), Intent(In)    :: NWPMetDefn2 ! The NWP met definition (part 2 - met file structure
                                                   ! definition).
  Type(AncillaryMetDefns_),    Intent(InOut) :: MetDefns    ! The collection of met definitions.
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

End Subroutine AddAncillaryMetDefn2

!-------------------------------------------------------------------------------------------------------------

Function NWPMetDefn2Eq(NWPMetDefn1, NWPMetDefn2)
! Tests for equality of NWP met definitions (part 2 - met file structure definitions).

  Implicit None
  ! Argument list:
  Type(AncillaryMetDefn2_), Intent(In) :: NWPMetDefn1 !} The two NWP met definitions (part 2
  Type(AncillaryMetDefn2_), Intent(In) :: NWPMetDefn2 !} - met file structure definitions).
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
      (NWPMetDefn1%ThreeD(i)        .eqv.  NWPMetDefn2%ThreeD(i)       )
  End Do

End Function NWPMetDefn2Eq

!-------------------------------------------------------------------------------------------------------------

Function FindAncillaryMetDefn2Index(Name, MetDefns)
! Finds the index of an NWP met definition (part 2 - met file structure definition).

  Implicit None
  ! Argument list:
  Character(*),    Intent(In) :: Name     ! Name of NWP met definition (part 2 - met
                                          ! file structure definition).
  Type(AncillaryMetDefns_), Intent(In) :: MetDefns ! Collection of met definitions.
  ! Function result:
  Integer :: FindAncillaryMetDefn2Index ! Index of NWP met definition (part 2 - met file
                                  ! structure definition).
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, MetDefns%nNWPMetDefn2s
    If (Name .CIEq. MetDefns%NWPMetDefn2s(i)%Name) Then
      FindAncillaryMetDefn2Index = i
      Return
    End If
  End Do

  Call Message(                                                &
         'FATAL ERROR: NWP Met File Structure Definition "' // &
         Trim(Name)                                         // &
         '" not found',                                        &
         3                                                     &
       )

End Function FindAncillaryMetDefn2Index

!-------------------------------------------------------------------------------------------------------------

Function InitAncillaryMet(             &
           MetName,                    &
           FixedMet,                   &
           UpdateOnDemand,             &
           NWPMetDefn, NWPMetDefn2,    &
       !    HMin, HMax, UseNWPH,        &
       !    SigUUM, TauUUM,             &
       !    SigU2HPlus, SigW2HPlus,     &
       !    TauUHPlus, TauWHPlus,       &
           MetFolder,                  &
           MetFolderStem,              &
           MetFolders,                 &
           EnsembleMetFolder           &
       !    TopogFolder,                &
       !    RestoreMetScript, DeleteMet &
         )                             &
Result(NWPMet)
! Initialises the state of an NWP met module instance.

  Implicit None
  ! Argument list:
  Character(*),              Intent(In) :: MetName           ! Name of met module instance.
  Logical,                   Intent(In) :: FixedMet          ! Indicates that the met is fixed.
  Logical,                   Intent(In) :: UpdateOnDemand    ! Indicates the met module instance is to be
                                                             ! updated using update-on-demand.
  Type(AncillaryMetDefn_),   Intent(In) :: NWPMetDefn        ! NWP met definition (part 1 - basic definitions).
  Type(AncillaryMetDefn2_),  Intent(In) :: NWPMetDefn2       ! NWP met definition (part 2 - met file structure
                                                             ! definition).
  Character(*),              Intent(In) :: MetFolder         ! Folder containing met data.
  Character(*),              Intent(In) :: MetFolderStem     ! Stem of folders containing met data.
  Character(MaxTokenLength), Intent(In) :: MetFolders(:)     ! Array of folders containing met data.
                                                             ! $$ changed from * to support Intel compiler on Cray
  Character(*),              Intent(In) :: EnsembleMetFolder ! Folder containing ensemble of met data.
  ! Function result:
  Type(AncillaryMet_) :: NWPMet ! Initialised state of an NWP met module instance.
  ! Locals:
  Integer :: L ! Trimed length of MetFolder(s) or TopogFolder.
  Integer :: i ! Loop index.

  NWPMet%C = InitCommonMet('Ancillary Met', MetName, FixedMet, UpdateOnDemand)

  NWPMet%MetDefn = NWPMetDefn

  NWPMet%MetDefn2 = NWPMetDefn2

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

  NWPMet%SpaceAllocated = .false.

  NWPMet%New = 0
  NWPMet%Old = 0

End Function InitAncillaryMet

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpCoordsEtc_AncillaryMet(NWPMet, Coords, Grids)
! Sets up Coords and Grids by adding any extra coords and grids which NWPMet wants to
! define.

  Implicit None
  ! Argument list:
  Type(AncillaryMet_), Intent(In)    :: NWPMet ! State of an NWP met module instance.
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

End Subroutine SetUpCoordsEtc_AncillaryMet

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpAncillaryMet_CoordsEtc(EtaDefns, Coords, Grids, MetEnsembleSize, NWPMet)
! Sets up NWPMet using information from EtaDefns, Coords and Grids. ! $$ not quite true now MetEnsembleSize
! added.

  Implicit None
  ! Argument list:
  Type(EtaDefns_), Intent(In)           :: EtaDefns        ! Collection of eta definitions.
  Type(Coords_),   Intent(In)           :: Coords          ! Collection of coord systems.
  Type(Grids_),    Intent(In),   Target :: Grids           ! Collection of grids.
  Integer,         Intent(In)           :: MetEnsembleSize ! Size of the met ensemble (i.e. number of met
                                                           ! realisations).
  Type(AncillaryMet_),   Intent(InOut)        :: NWPMet          ! State of an NWP met module instance.
  ! Locals:
  Character(MaxCharLength)         :: ZCoordName1 !} Names of the vertical coord
  Character(MaxCharLength)         :: ZCoordName2 !} systems to add to NWPMet%C.
  Character(MaxCharLength)         :: ZCoordName3 !}
  Type(HGrid_),            Pointer :: HGrid       !] Abbreviations for grids.
  Logical                          :: Error       ! Error flag - indicates fatal error
                                                  ! which however will not be made
                                                  ! fatal just yet in order to detect
                                                  ! other possible errors.
  Real(Std)                        :: Scale       ! Tolerance for computing grid
                                                  ! offset / grid spacing.
  Integer                          :: Offset      ! 1, 0 and -1 indicate a grid offset
                                                  ! of 1/2, 0 and -1/2 times the grid
                                                  ! spacing.
  Logical                          :: OffsetError ! Indicates the grid offset is not
                                                  ! acceptable.
  Integer                          :: i           !} Loop indices.
  Integer                          :: j           !}
  Integer                          :: k           !}

  ! Initialise Error.
  Error = .false.

  ! Find indices of coords and grids.
  NWPMet%iHGrid  = FindHGridIndex(NWPMet%MetDefn%HGridName,                Grids )
  NWPMet%iHCoord = FindHCoordIndex(Grids%HGrids(NWPMet%iHGrid)%HCoordName, Coords)

  ! Add coord systems to NWPMet%C.
  Call AddCoordsToCommonMet(                              &
         1, (/ Grids%HGrids(NWPMet%iHGrid)%HCoordName /), &
         0, (/ ' ' /),                                    &
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
  If (NWPMet%C%nZCoords /= 0) Then
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
  HGrid   => Grids%HGrids(NWPMet%iHGrid)

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

  ! Check vertical grids have index array. $$

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

  ! Set up kMax.
  Do i = 1, nFieldNames
    NWPMet%kMax(i) = Max(FieldDimensions(i), 1) ! Note FieldDimensions = 0 is code for 2-D field.
  End Do

  ! Allocate FieldPresent and FieldPersisted.
  Allocate(                                                  &
    NWPMet%FieldPresent(MaxVal(NWPMet%kMax(:)), nFieldNames) &
  )
  Allocate(                                                    &
    NWPMet%FieldPersisted(MaxVal(NWPMet%kMax(:)), nFieldNames) &
  )

  ! Initialise Error and FieldPersisted.
  Error                       = .false.
  NWPMet%FieldPersisted(:, :) = 0

  ! Determine the field index and dimension of each field in NWPMetDefn2.
  Do i = 1, NWPMet%MetDefn2%nFields
    NWPMet%iField(i) = 0
    Do j = 1, nFieldNames
      If (NWPMet%MetDefn2%FieldNames(i) .CIEq. FieldNames(j)) Then
        NWPMet%iField(i) = j
        NWPMet%ThreeD(i) = FieldDimensions(j) /= 0
      End If
    End Do
    If (NWPMet%iField(i) == 0) Then
      Call Message(                                          &
             'ERROR: the met file structure definition "' // &
             Trim(NWPMet%MetDefn2%Name)                   // &
             '" contains an invalid field name "'         // &
             Trim(NWPMet%MetDefn2%FieldNames(i))          // &
             '".',                                           &
             2                                               &
           )
      Error = .true.
    End If
  End Do

  ! Set up k1 and k2.
  Do i = 1, NWPMet%MetDefn2%nFields
    If (NWPMet%ThreeD(i) .and. NWPMet%iField(i) /= 0) Then
      NWPMet%k1(i) = NWPMet%MetDefn2%LowestLevels(i)
      If (NWPMet%MetDefn2%Top(i)) Then
        NWPMet%k2(i) = NWPMet%kMax(NWPMet%iField(i))
      Else
        NWPMet%k2(i) = NWPMet%MetDefn2%HighestLevels(i)
      End If
      If (NWPMet%k2(i) > NWPMet%kMax(NWPMet%iField(i))) Then
        Call Message(                                                    &
               'ERROR: The highest level of the '                     // &
               Trim(Int2Char(i)) // Int2Ordinal(i)                    // &
               ' field ('                                             // &
               Trim(NWPMet%MetDefn2%FieldNames(i))                    // &
               ') in the met file structure definition "'             // &
               Trim(NWPMet%MetDefn2%Name)                             // &
               '" is greater than the number of levels in the grid.',    &
               2                                                         &
             )
        Error = .true.
      End If
      If (NWPMet%k2(i) < NWPMet%k1(i)) Then
        Call Message(                                        &
               'ERROR: The highest level of the '         // &
               Trim(Int2Char(i)) // Int2Ordinal(i)        // &
               ' field ('                                 // &
               Trim(NWPMet%MetDefn2%FieldNames(i))        // &
               ') in the met file structure definition "' // &
               Trim(NWPMet%MetDefn2%Name)                 // &
               '" is less than the lowest level',            &
               2                                             &
             )
        Error = .true.
      End If
      If (NWPMet%k1(i) < 1) Then
        Call Message(                                        &
               'ERROR: The lowest level of the '          // &
               Trim(Int2Char(i)) // Int2Ordinal(i)        // &
               ' field ('                                 // &
               Trim(NWPMet%MetDefn2%FieldNames(i))        // &
               ') in the met file structure definition "' // &
               Trim(NWPMet%MetDefn2%Name)                 // &
               '" is less than 1',                           &
               2                                             &
             )
        Error = .true.
      End If
    Else
      NWPMet%k1(i) = 1
      NWPMet%k2(i) = 1
    End If
  End Do

  ! Fatal error.
  If (Error) Then
    Call Message(                                                &
           'FATAL ERROR: the met file structure definition "' // &
           Trim(NWPMet%MetDefn2%Name)                         // &
           '" is not useable.',                                  &
           3                                                     &
         )
  End If

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

End Subroutine SetUpAncillaryMet_CoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForUpdateAncillaryMet( &
             Coords, Grids,              &
             iCase, iMetCase,            &
             Time,                       &
             TValid, UpdateNow,          &
             NWPMet,                     &
             Units                       &
           )
! Prepares for updating an instance of the ancillary met module.

! This routine must set TValid and UpdateNow but must not alter the validity of the ancillary met module
! instance.
! $$ this comment should be with generic definition.

  Implicit None
  ! Argument List:
  Type(Coords_),       Intent(In)           :: Coords
  Type(Grids_),        Intent(In),   Target :: Grids
  Integer,             Intent(In)           :: iCase
  Integer,             Intent(In)           :: iMetCase
  Type(Time_),         Intent(In)           :: Time
  Type(Time_),         Intent(Out)          :: TValid
  Logical,             Intent(Out)          :: UpdateNow
  Type(AncillaryMet_), Intent(InOut)        :: NWPMet
  Type(Units_),        Intent(InOut)        :: Units
  ! Coords       :: Collection of coord systems.
  ! Grids        :: Collection of grids.
  ! iCase        :: Number of case.
  ! iMetCase     :: Number of the met realisation in the met ensemble.
  ! Time         :: Time for which the ancillary met module instance might be updated.
  ! TValid       :: Earliest time that the validity (overall or for any single attribute) of the met module
  !                 instance might change, assuming the met module instance is updated now. The value is that
  !                 determined at the end of this routine (the actual time may be later).
  ! UpdateNow    :: Indicates the met module instance must be updated now (even if update-on-demand is
  !                 specified). If set, TValid need not be set to any particular time.
  ! AncillaryMet :: State of an ancillary met module instance.
  ! Units        :: Collection of information on input/output unit numbers.
  ! Locals:
  Type(Time_) :: MetTime1 !} Times of the met data required.
  Type(Time_) :: MetTime2 !}

  ! Calculate times of met data.
  If (NWPMet%MetDefn%Monthly) Then
    MetTime1 = RoundToMonth(Time,     Up = .false., Strictly = .false.)
    MetTime2 = RoundToMonth(MetTime1, Up = .true.,  Strictly = .true. )
  Else If (IsInfFuture(NWPMet%MetDefn%Dt)) Then
    MetTime1 = InfPastTime()
    MetTime2 = InfFutureTime()
  Else
    MetTime1 = Round(Time, NWPMet%MetDefn%T0, NWPMet%MetDefn%Dt, Up = .false.)
               ! Could name files using a time zone $$
               ! Currently files named in UTC.
    MetTime2 = MetTime1 + NWPMet%MetDefn%Dt
  End If

  ! Set validity flags.
  If (NWPMet%C%FixedMet) Then
    TValid = InfFutureTime()
  Else
    TValid = MetTime2
  End If
  UpdateNow = .false.

End Subroutine PrepareForUpdateAncillaryMet

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateAncillaryMet( &
             Coords, Grids,    &
             iCase, iMetCase,  &
             Time,             &
             NWPMet,           &
             Units             &
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
  Type(AncillaryMet_), Intent(InOut)        :: NWPMet   ! State of an NWP met module instance.
  Type(Units_),  Intent(InOut)        :: Units    ! Collection of information on input/output unit numbers.
  ! Locals:
  Type(Time_)              :: MetTime1  !} Times and short times of the met data
  Type(Time_)              :: MetTime2  !} required.
  Type(ShortTime_)         :: SMetTime1 !}
  Type(ShortTime_)         :: SMetTime2 !}
  Logical                  :: Error     ! Indicates that an error occurred in ReadNWPMet.
  Type(HGrid_),    Pointer :: HGrid     ! Abbreviations for grids.
  Logical                  :: DoIfBlock ! Indicates an if-block is to be executed
                                        ! where the logical test can't be evaluated
                                        ! in-line because of risk of error in testing
                                        ! equality of times.
  Integer                  :: i         ! Loop index.
  Integer                  :: k         ! Loop index

  ! Set up abbreviations for grids.
  HGrid => Grids%HGrids(NWPMet%iHGrid)

  ! Allocate arrays.
  If (.not.NWPMet%SpaceAllocated) Then
    Allocate(NWPMet%LandUseFracs (HGrid%nX, HGrid%nY, NWPMet%kMax(F_LandUseFracs),  2))
    Allocate(NWPMet%ClayMassFrac (HGrid%nX, HGrid%nY, NWPMet%kMax(F_ClayMassFrac),  2))
    Allocate(NWPMet%SiltMassFrac (HGrid%nX, HGrid%nY, NWPMet%kMax(F_SiltMassFrac),  2))
    Allocate(NWPMet%SandMassFrac (HGrid%nX, HGrid%nY, NWPMet%kMax(F_SandMassFrac),  2))
    Allocate(NWPMet%SoilMassFracs(HGrid%nX, HGrid%nY, NWPMet%kMax(F_SoilMassFracs), 2))
    Allocate(NWPMet%SoilMoisture (HGrid%nX, HGrid%nY, NWPMet%kMax(F_SoilMoisture),  2))
    Allocate(NWPMet%LAI          (HGrid%nX, HGrid%nY, NWPMet%kMax(F_LAI),           2))
    Allocate(NWPMet%CanopyHeight (HGrid%nX, HGrid%nY, NWPMet%kMax(F_CanopyHeight),  2))
    Allocate(NWPMet%LandFrac     (HGrid%nX, HGrid%nY, NWPMet%kMax(F_LandFrac),   2))
    NWPMet%SpaceAllocated = .true.
  End If

  ! Calculate times of met data.
  If (NWPMet%MetDefn%Monthly) Then
    MetTime1 = RoundToMonth(Time,     Up = .false., Strictly = .false.)
    MetTime2 = RoundToMonth(MetTime1, Up = .true.,  Strictly = .true. )
  Else If (IsInfFuture(NWPMet%MetDefn%Dt)) Then
    MetTime1 = InfPastTime()
    MetTime2 = InfFutureTime()
  Else
    MetTime1 = Round(Time, NWPMet%MetDefn%T0, NWPMet%MetDefn%Dt, Up = .false.)
               ! Could name files using a time zone $$
               ! Currently files named in UTC.
    MetTime2 = MetTime1 + NWPMet%MetDefn%Dt
  End If
  SMetTime1 = Time2ShortTime(MetTime1)
  SMetTime2 = Time2ShortTime(MetTime2)
  NWPMet%Dt = SMetTime2 - SMetTime1

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

  ! Required so that monthly ancillaries are reread at start of each month (00Z on 1st of month).
  DoIfBlock = DoIfBlock .or. (NWPMet%MetDefn%Monthly .and. Time == MetTime1)

  If (DoIfBlock) Then

    ! Mark met data in Old time level as invalid (otherwise ReadNWPMet will assume Old
    ! time level data is for the previous time and may use it to fix-up missing data).
    NWPMet%Old = 0

    ! Read Time1 met data into New time level.
    NWPMet%New     = 1
    NWPMet%NewTime = SMetTime1
    Call ReadNWPMet(MetTime1, iCase, iMetCase, Coords, Grids, Error, NWPMet, Units)
    If (Error) Go To 9

  End If

  ! Move Time1 met data into Old time level.
  NWPMet%Old     = NWPMet%New
  NWPMet%OldTime = NWPMet%NewTime

  ! Read Time2 met data into New time level.
  ! Note that for monthly data and infinite dt data we simply duplicate the data. Monthly data is regarded as
  ! representative of the month and infinite dt data is valid for all times.
  ! $$ note we could set up monthly to interpolate in NWPFlow, possibly as an option. Then would
  ! have to read two sets of data and alter file name calculation for backwards runs too. Could signal this
  ! with 'Last' (like 'NextHeatFlux') - the absence of Last and Next would indicate interpolation. Could also
  ! allow Yearly and "Yearly and Monthly" - useful for e.g. long climate change sims.
  ! See also "$$" in NWPFlow.
  NWPMet%New     = 3 - NWPMet%Old
  NWPMet%NewTime = SMetTime2
  If (NWPMet%MetDefn%Monthly) Then
    NWPMet%LandUseFracs (:, :, :, NWPMet%New) = NWPMet%LandUseFracs (:, :, :, NWPMet%Old)
    NWPMet%ClayMassFrac (:, :, :, NWPMet%New) = NWPMet%ClayMassFrac (:, :, :, NWPMet%Old)
    NWPMet%SiltMassFrac (:, :, :, NWPMet%New) = NWPMet%SiltMassFrac (:, :, :, NWPMet%Old)
    NWPMet%SandMassFrac (:, :, :, NWPMet%New) = NWPMet%SandMassFrac (:, :, :, NWPMet%Old)
    NWPMet%SoilMassFracs(:, :, :, NWPMet%New) = NWPMet%SoilMassFracs(:, :, :, NWPMet%Old)
    NWPMet%SoilMoisture (:, :, :, NWPMet%New) = NWPMet%SoilMoisture (:, :, :, NWPMet%Old)
    NWPMet%LAI          (:, :, :, NWPMet%New) = NWPMet%LAI          (:, :, :, NWPMet%Old)
    NWPMet%CanopyHeight (:, :, :, NWPMet%New) = NWPMet%CanopyHeight (:, :, :, NWPMet%Old)
    NWPMet%LandFrac     (:, :, :, NWPMet%New) = NWPMet%LandFrac     (:, :, :, NWPMet%Old)
    NWPMet%ValidAttribs(:, NWPMet%New)        = NWPMet%ValidAttribs(:, NWPMet%Old)
  Else If (IsInfFuture(NWPMet%MetDefn%Dt)) Then
    NWPMet%LandUseFracs (:, :, :, NWPMet%New) = NWPMet%LandUseFracs (:, :, :, NWPMet%Old)
    NWPMet%ClayMassFrac (:, :, :, NWPMet%New) = NWPMet%ClayMassFrac (:, :, :, NWPMet%Old)
    NWPMet%SiltMassFrac (:, :, :, NWPMet%New) = NWPMet%SiltMassFrac (:, :, :, NWPMet%Old)
    NWPMet%SandMassFrac (:, :, :, NWPMet%New) = NWPMet%SandMassFrac (:, :, :, NWPMet%Old)
    NWPMet%SoilMassFracs(:, :, :, NWPMet%New) = NWPMet%SoilMassFracs(:, :, :, NWPMet%Old)
    NWPMet%SoilMoisture (:, :, :, NWPMet%New) = NWPMet%SoilMoisture (:, :, :, NWPMet%Old)
    NWPMet%LAI          (:, :, :, NWPMet%New) = NWPMet%LAI          (:, :, :, NWPMet%Old)
    NWPMet%CanopyHeight (:, :, :, NWPMet%New) = NWPMet%CanopyHeight (:, :, :, NWPMet%Old)
    NWPMet%LandFrac     (:, :, :, NWPMet%New) = NWPMet%LandFrac     (:, :, :, NWPMet%Old)
    NWPMet%ValidAttribs(:, NWPMet%New)        = NWPMet%ValidAttribs(:, NWPMet%Old)
  Else If (.not. IsInfFuture(NWPMet%MetDefn%Dt)) Then
    Call ReadNWPMet(MetTime2, iCase, iMetCase, Coords, Grids, Error, NWPMet, Units)
  End If

  ! Adjust land use fractions for sea points (values < 0) and coastal grid points (0 < land fraction < 1)
  Do k = 1, NWPMet%kMax(F_LandUseFracs)
    NWPMet%LandUseFracs(:, :, k, :) = Max(0.0, NWPMet%LandUseFracs(:, :, k, :)) * NWPMet%LandFrac(:, :, 1, :)
  End Do
  ! sea treated as inland water
  NWPMet%LandUseFracs(:, :, 7, :) = NWPMet%LandUseFracs(:, :, 7, :) + 1.0 - NWPMet%LandFrac(:, :, 1, :)

  ! Set individual attribute validity flags.
  Do i = 1, nAttribs
    NWPMet%C%ValidAttribs(i) = NWPMet%ValidAttribs(i, 1) .and. NWPMet%ValidAttribs(i, 2)
  End Do
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

End Subroutine UpdateAncillaryMet

!-------------------------------------------------------------------------------------------------------------

Subroutine ResetAncillaryMet(NWPMet)
! Resets the state of an NWP met module instance for a new realisation.

  Implicit None
  ! Argument list:
  Type(AncillaryMet_), Intent(InOut) :: NWPMet ! State of an NWP met module instance.

  ! Reset validity flags for met data (e.g. for a new realisation in an ensemble of met data).
  NWPMet%New = 0
  NWPMet%Old = 0

End Subroutine ResetAncillaryMet

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadNWPMet(Time, iCase, iMetCase, Coords, Grids, Error, M, Units)
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
  Type(Time_),   Intent(In)           :: Time     ! Time for which NWP met data is to be read.
  Integer,       Intent(In)           :: iCase    ! Number of case.
  Integer,       Intent(In)           :: iMetCase ! Index of the met realisation in the met ensemble.
  Type(Coords_), Intent(In),   Target :: Coords   ! Collection of coord systems.
  Type(Grids_),  Intent(In),   Target :: Grids    ! Collection of grids.
  Logical,       Intent(Out)          :: Error    ! Indicates an error occured in this routine.
  Type(AncillaryMet_), Intent(InOut)        :: M        ! State of an NWP met module instance.
  Type(Units_),  Intent(InOut)        :: Units    ! Collection of information on input/output unit numbers.
  ! Local types:
  Type :: F2_ ! 2-d (horizontal) field (with time dimension too).
    Real(Std), Pointer :: P(:,:,:) ! 2-d (horizontal) field.
  End Type F2_
  Type :: F3_ ! 3-d field (with time dimension too).
    Real(Std), Pointer :: P(:,:,:,:) ! 3-d field.
  End Type F3_
  ! Locals:
  Real(Std),                   Target  :: Dummy(Grids%HGrids(M%iHGrid)%nX, Grids%HGrids(M%iHGrid)%nY, 2)
  Type(F3_)                            :: F3(nFieldNames)
  Integer                              :: iMetCaseL
  Character(MaxFileNameLength)         :: MetFolderPart1
  Character(MaxFileNameLength)         :: MetFolderPart2
  Character(MaxFileNameLength)         :: MetFile
  Logical                              :: Exist
  Integer                              :: IOStat
  Logical                              :: InconsistentHeader
  Logical                              :: Missing
  Logical                              :: NewFile
  Integer                              :: Unit
  Type(HCoord_),               Pointer :: HCoord
  Type(HGrid_),                Pointer :: HGrid
  Integer                              :: i
  Integer                              :: j
  Integer                              :: k
  Integer                              :: iField
  Integer                              :: iNWPFieldLevel
  Integer                              :: iNameIIIFieldLevel
  Integer                              :: iAttrib
  Character(MaxCharLength)             :: CharTime
  Integer                              :: indexML
  Integer                              :: indexSL
  Character(MaxCharLength)             :: ParameterKeysML
  Character(MaxCharLength)             :: ParameterKeysSL
  Integer                              :: Status
  Character(MaxCharLength)             :: ErrMessage
  ! Dummy              :} Dummy 2-d (horizontal) fields (with time dimension too) with horizontal dimensions
  ! DummyU             :} corresponding to the main, U and V grids.
  ! DummyV             :}
  ! F3                 :] 3-d and 2-d (horizontal) fields (with time dimension too) used as alias for fields
  ! F2                 :] in M and for fields Dummy, DummyU and DummyV.
  ! iMetCaseL          :: Local copy of iMetCase.
  ! MetFolderPart1     :} Met folder parts 1 and 2.
  ! MetFolderPart2     :}
  ! MetFile            :: Met file.
  ! Exist              :: Indicates that a file exists.
  ! IOStat             :: Error code for open statement.
  ! InconsistentHeader :: Indicates the header for the field is inconsistent with what is expected.
  ! Missing            :: Indicates a field is missing from the met file.
  ! NewFile            :: Indicates a newly opened file is being considered.
  ! Unit               :: Input/output unit number.
  ! HCoord             :} Abbreviations for coords, grids.
  ! HGrid              :}
  ! i                  :] Loop indices.
  ! j                  :]
  ! k                  :]
  ! iField             :]
  ! iNWPFieldLevel     :: Index of NWP field level.
  ! iNameIIIFieldLevel :: Index of Name III field level.
  ! iAttrib            :: Index of attribute.
  ! CharTime           :: Character version of Time.
  ! indexML            :: Identifier for index providing model-level access to a GRIB met file.
  ! indexSL            :: Identifier for index providing single-level access to a GRIB met file.
  ! ParameterKeysML    :: Comma-separated list of keys for creating GRIB file index incl. model level.
  ! ParameterKeysSL    :: Comma-separated list of keys for creating GRIB file index excl. model level.
  ! Status             :: Status code from GRIB_API subroutines.
  ! ErrMessage         :: Error message from GRIB_API subroutines.

  ! 1) Set up error flag and abbreviations for grids.

  ! Set error flag to false.
  Error = .false.

  ! Set up abbreviations for coords and grids.
  HCoord  => Coords%HCoords(M%iHCoord)
  HGrid   => Grids%HGrids(M%iHGrid  )

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

  ! Set iMetCaseL. This should be -1 unless the met case number is to be used in finding the
  ! met realisation in
  ! a met file with several met realisations.
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

  ! Derive met file name (note that, for monthly data, the files are labelled in effect with the start of the
  ! month, with "start" here meaning in chronological time not run direction time - hence the use of
  ! IsBackwards).
  If (M%MetDefn%Monthly) Then
    If (IsBackwards()) Then
      CharTime = FileNameTime(RoundToMonth(Time, Up = .true., Strictly = .true.))
    Else
      CharTime = FileNameTime(Time)
    End If
    MetFile =                     &
      Trim(M%MetDefn%Prefix)   // &
      CharTime(5:6)            // &
      '.'                      // &
      Trim(M%MetDefn%Suffix)
  Else If (IsInfFuture(M%MetDefn%Dt)) Then
    MetFile =                     &
      Trim(M%MetDefn%Prefix)   // &
      '.'                      // &
      Trim(M%MetDefn%Suffix)
  Else
    MetFile =                     &
      Trim(M%MetDefn%Prefix)   // &
      Trim(FileNameTime(Time)) // &
      '.'                      // &
      Trim(M%MetDefn%Suffix)
  End If

  ! Open met file (GRIB files are 'opened' by creating an index using the relevant GRIB_API routine).
  
  ! GRIB met file.
  If (M%MetDefn%FileType .CIEq. 'GRIB') Then

    ! Check that file exists
    Inquire (                                                                                       &
      File  = Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile))), &
      Exist = Exist                                                                                 &
    )
    If (.not. Exist) Then
      Call Message(                                                                                   &
             "ERROR: Can't open met file "                                                         // &
             Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile))),    &
             2                                                                                        &
           )
      Error = .true.
      Go To 9
    End If

#   ifdef GRIBsupport

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
      Call grib_index_create(                                                                      &
             indexML,                                                                              &
             Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile))), &
             ParameterKeysML, Status                                                               &
           )
      
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                                                  &
               'ERROR: Unable to create index for model-level access to the GRIB met file '         // &
               Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile))) // &
               ' (GRIB_API returns with the error message "'                                        // &
               Trim(ErrMessage)                                                                     // &
               '")',                                                                                   &
               2                                                                                       &
             )
        Error = .true.
        Go To 9
      End If

      ! Create index providing single-level access.
      Call grib_index_create(                                                                      &
             indexSL,                                                                              &
             Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile))), &
             ParameterKeysSL, Status                                                               &
           )
      
      If (Status /= GRIB_SUCCESS) Then
        Call grib_get_error_string(Status, ErrMessage)
        Call Message(                                                                                  &
               'ERROR: Unable to create index for single-level access to the GRIB met file '        // &
               Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile))) // &
               ' (GRIB_API returns with the error message "'                                        // &
               Trim(ErrMessage)                                                                     // &
               '")',                                                                                   &
               2                                                                                       &
             )
        Error = .true.
        Go To 9
      End If

#   endif

  ! Name II or PP met file.
  Else

    Call GetNewUnit(Unit, Units) !$$ use openfile once have sorted out how to handle 'convert ='
#   ifdef UseConvert
      Open (                                                                                            &
        Unit    = Unit,                                                                                 &
        File    = Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile))), &
        Status  = 'Old',                                                                                &
        Form    = 'Unformatted',                                                                        &
        Action  = 'Read',                                                                               &
        Convert = Trim(M%MetDefn%BinaryFormat),                                                         &
        IOStat  = IOStat                                                                                &
      )
#   else
      Open (                                                                                           &
        Unit   = Unit,                                                                                 &
        File   = Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile))), &
        Status = 'Old',                                                                                &
        Form   = 'Unformatted',                                                                        &
        Action = 'Read',                                                                               &
        IOStat = IOStat                                                                                &
      )
#   endif
    If (IOStat /= 0) Then
      Call CloseUnit(Unit, Units)
      Call Message(                                                                                   &
             "ERROR: Can't open met file "                                                         // &
             Trim(ConvertFileName(Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile))),    &
             2                                                                                        &
           )
      Error = .true.
      Go To 9
    End If

  End If

  ! 3) Set up equivalences to enable looping over fields and initialise FieldPresent and NewFile.

  F3(F_LandUseFracs )%P => M%LandUseFracs
  F3(F_ClayMassFrac )%P => M%ClayMassFrac
  F3(F_SiltMassFrac )%P => M%SiltMassFrac
  F3(F_SandMassFrac )%P => M%SandMassFrac
  F3(F_SoilMassFracs)%P => M%SoilMassFracs
  F3(F_SoilMoisture )%P => M%SoilMoisture
  F3(F_LAI          )%P => M%LAI
  F3(F_CanopyHeight )%P => M%CanopyHeight
  F3(F_LandFrac     )%P => M%LandFrac

  M%FieldPresent(:, :) = .false.

  NewFile = .true.

  ! 4) Read data (note missing data is ignored at times which don't correspond to the time of interest).

  ! $$ Note currently all files are read in prescribed order. Could change for PP met
  ! $$ files and possibly others.

  ! Reading fields in the order they are in M%MetDefn2.
  If (.true.) Then

    LoopOverFieldsAtGivenTime: Do iField = 1, M%MetDefn2%nFields

      LoopOverLevels: Do k = M%k1(iField), M%k2(iField)

     !   If (.not. M%ThreeD(iField)) Then ! $$ note this is what NWPMet does (incorrectly)
        If (.not. M%MetDefn2%ThreeD(iField)) Then
          iNWPFieldLevel = 0
        Else
          iNWPFieldLevel = k ! $$ here need to use index array (part of grids in NWPMet).
        End If

        Select Case (M%iField(iField))
          ! U grid: U.
      !    Case (F_U)
      !      Call Read2dField(                                                                &
      !             NewFile,                                                                  &
      !             M%MetDefn%FileType,                                                       &
      !             Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile),            &
      !             indexML, indexSL, Unit,                                                   &
      !             M%MetDefn2%FieldCodes(iField), M%MetDefn2%ThreeD(iField), iNWPFieldLevel, &
      !             iMetCaseL, Time,                                                          &
      !             HCoord, HGridU,                                                           &
      !             F3(M%iField(iField))%P(:, :, k, M%New),                                   &
      !             Error, InconsistentHeader, Missing                                        &
      !           )
          ! Main grid: all other fields.
          Case Default
            Call Read2dField(                                                                &
                   NewFile,                                                                  &
                   M%MetDefn%FileType,                                                       &
                   Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile),            &
                   indexML, indexSL, Unit,                                                   &
                   M%MetDefn2%FieldCodes(iField), M%MetDefn2%ThreeD(iField), iNWPFieldLevel, &
                   iMetCaseL, Time,                                                          &
                   HCoord, HGrid,                                                            &
                   F3(M%iField(iField))%P(:, :, k, M%New),                                   &
                   Error, InconsistentHeader, Missing                                        &
                 )
        End Select

        ! Read errors.
        If (Error) Then
          If (M%ThreeD(iField)) Then
            Call Message(                                           &
                   'ERROR: read error occured in reading field ' // &
                   Trim(M%MetDefn2%FieldNames(iField))           // &
                   ' from the met file '                         // &
                   Trim(MetFile)                                 // &
                   ' for the '                                   // &
                   Trim(Int2Char(k)) // Int2Ordinal(k)           // &
                   ' (Name III) vertical level',                    &
                   2                                                &
                 )
          Else
            Call Message(                                           &
                   'ERROR: read error occured in reading field ' // &
                   Trim(M%MetDefn2%FieldNames(iField))           // &
                   ' from the met file '                         // &
                   Trim(MetFile),                                   &
                   2                                                &
                 )
          End If
          Go To 9
        End If

        ! Inconsistent header.
        If (InconsistentHeader) Then
          If (M%ThreeD(iField)) Then
            Call ControlledMessage(                                               &
                   'The field '                                                // &
                   Trim(M%MetDefn2%FieldNames(iField))                         // &
                   ' from met file '                                           // &
                   Trim(MetFile)                                               // &
                   ' for the '                                                 // &
                   Trim(Int2Char(k)) // Int2Ordinal(k)                         // &
                   ' (Name III) vertical level '                               // &
                   ' has a header which is inconsistent with the information ' // &
                   'given in the NWP Met Definition '                          // &
                   Trim(M%MetDefn%Name)                                        // &
                   '.',                                                           &
                   MessageControls    = GlobalMessageControls,                    &
                   MessageControlName = 'Inconsistent NWP header',                &
                   ErrorCode          = 2                                         &
                 ) ! $$ Set Error = .true. and goto 9? or downgrade to warning
          Else
            Call ControlledMessage(                                               &
                   'The field '                                                // &
                   Trim(M%MetDefn2%FieldNames(iField))                         // &
                   ' from met file '                                           // &
                   Trim(MetFile)                                               // &
                   ' has a header which is inconsistent with the information ' // &
                   'given in the NWP Met Definition '                          // &
                   Trim(M%MetDefn%Name)                                        // &
                   '.',                                                           &
                   MessageControls    = GlobalMessageControls,                    &
                   MessageControlName = 'Inconsistent NWP header',                &
                   ErrorCode          = 2                                         &
                 ) ! $$ Set Error = .true. and goto 9? or downgrade to warning
          End If
        End If

        ! FieldPresent.
        M%FieldPresent(k, M%iField(iField)) = .not.Missing

        ! NewFile.
        NewFile = .false.

      End Do LoopOverLevels

    End Do LoopOverFieldsAtGivenTime

  ! Reading fields in the order they are in the file.
  Else

  End If

  ! Close met file (or delete its index if a GRIB file).
  
  ! GRIB met file.
  If (M%MetDefn%FileType .CIEq. 'GRIB') Then

#   ifdef GRIBsupport

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

#   endif

  ! Name II or PP met file.
  Else

    Call CloseUnit(Unit, Units)

  End If

  ! 5) Warnings for missing fields which are expected to be present according to MetDefn2. Note this is done
  ! here rather than as each field is read, to avoid producing warnings when more than one possible NWP field
  ! can be used for a given NAME III field and some but not all of these NWP fields are missing.

  Do iField = 1, M%MetDefn2%nFields

    Do k = M%k1(iField), M%k2(iField)

      If (.not. M%FieldPresent(k, M%iField(iField))) Then

        If (M%ThreeD(iField)) Then

          Call ControlledMessage(                            &
                 'Field '                              //    &
                 Trim(M%MetDefn2%FieldNames(iField))   //    &
                 ' missing from met file '             //    &
                 Trim(MetFile)                         //    &
                 ' for the '                           //    &
                 Trim(Int2Char(k)) // Int2Ordinal(k)   //    &
                 ' (Name III) vertical level',               &
                 MessageControls    = GlobalMessageControls, &
                 MessageControlName = 'Missing NWP field',   &
                 ErrorCode          = 1                      &
               )

        Else

          Call ControlledMessage(                            &
                 'Field '                              //    &
                 Trim(M%MetDefn2%FieldNames(iField))   //    &
                 ' missing from met file '             //    &
                 Trim(MetFile),                              &
                 MessageControls    = GlobalMessageControls, &
                 MessageControlName = 'Missing NWP field',   &
                 ErrorCode          = 1                      &
               )

        End If

      End If

    End Do

  End Do

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
  M%ValidAttribs(1:nAttribs, M%New) = .true.

  ! 3-d fields.
  Do iField = 1, nFieldNames
    Do k = 1, M%kMax(iField)
      ! Cycle if quantity present or can be fixed up. Note the test for FieldPersisted >= 0 which
      ! prevents default
      ! values being persisted.
      If (M%FieldPresent(k, iField)) Cycle
      If (                                                       &
        M%FieldPersisted(k, iField) < PersistTimes(iField) .and. &
        M%FieldPersisted(k, iField) >= 0                   .and. &
        M%Old /= 0                                               &
      ) Cycle
      If (AllowFixUpToDefault(iField)) Cycle
      Do iAttrib = 1, nAttribs
        If (Needed(iAttrib, iField)) M%ValidAttribs(iAttrib, M%New) = .false.
      End Do
!      ! Set Error and write error message.
!      Error = .true.
!      If (FieldDimensions(iField) /= 0) Then ! $$ note no convenient logical set up for this yet
!        Call Message(                                              &
!               'ERROR: Field '                                  // &
!               Trim(FieldNames3d(iField))                       // &
!               ' missing from met file '                        // &
!               Trim(MetFile)                                    // &
!               ' at (Name III) level '                          // &
!               Trim(Int2Char(k))                                // &
!               ' (and is required at this level and cannot be ' // &
!               'fixed up from other available information).',      &
!               2                                                   &
!             )
!      Else
!        Call Message(                                            &
!               'ERROR: Field '                                // &
!               Trim(FieldNames2d(iField))                     // &
!               ' missing from met file '                      // &
!               Trim(MetFile)                                  // &
!               ' (and is required and cannot be '             // &
!               'fixed up from other available information).',    &
!               2                                                 &
!             )
!      End If
    End Do
  End Do

  ! 6) Message.

  If (M%MetDefn%Monthly) Then
    Call Message(                                                           &
           'Ancillary met data '                                         // &
           'read from file '                                             // &
           Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile)    &
         )
  Else If (IsInfFuture(M%MetDefn%Dt)) Then
    Call Message(                                                           &
           'Ancillary met data '                                         // &
           'read from file '                                             // &
           Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile)    &
         )
  Else If (.not. IsInfFuture(M%MetDefn%Dt)) Then
    Call Message(                                                           &
           'Ancillary met data for '                                     // &
           Trim(Time2Char(Time, .false., 0, .true.))                     // &
           ' read from file '                                            // &
           Trim(MetFolderPart1) // Trim(MetFolderPart2) // Trim(MetFile)    &
         )
  End If

  Return

  ! 7) Errors.

9 Continue

  Call Message(                                         &
         'ERROR in reading met data for met module ' // &
         Trim(M%C%MetModName)                        // &
         '.'                                         // &
         Trim(M%C%MetName),                             &
         2                                              &
       )
  M%ValidAttribs(1:nAttribs, M%New) = .false.

End Subroutine ReadNWPMet

!-------------------------------------------------------------------------------------------------------------

Subroutine Read2dField(                                      &
             NewFile,                                        &
             FileType,                                       &
             FileName,                                       &
             indexML, indexSL, Unit,                         &
             FieldCode, ThreeD, iFieldLevel, iMetCase, Time, &
             HCoord, HGrid,                                  &
             Field,                                          &
             Error, InconsistentHeader, Missing              &
           )
! Reads in a 2-d (horizontal) field of met data from Name II, PP or GRIB met files.

  Implicit None
  ! Argument List:
  Logical,       Intent(In)         :: NewFile       ! Indicates a newly opened file is being considered.
  Character(*),  Intent(In)         :: FileType      ! Type of met file (Name II, PP or GRIB).
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
  Logical,       Intent(Out)        :: Missing       ! Indicates a missing field.
  ! Locals:
  Integer,                 Save :: IHeader(45)                        !} Header record for PP met files.
  Real(Std),               Save :: RHeader(19)                        !}
  Logical,                 Save :: MissingL                           ! Indicates the last field was missing
                                                                      ! for PP met files.
  Real(Std)                     :: Header(64)                         !} Header record for Name II met files.
  Real(Std)                     :: Scale                              ! Tolerance for checking headers.
  Real(Std)                     :: LocalField(HGrid%nX, HGrid%nY)     ! $$ needed on Intel for Read
                                                                      ! to work correctly.
  Real(Std)                     :: LocalField2(HGrid%nX * HGrid%nY)   ! Local array for reading met field.
  Character(MaxCharLength)      :: CharTime                           ! Character version of Time.
  Integer                       :: i                                  !} Loop indices.
  Integer                       :: j                                  !}
  Integer                       :: k                                  ! Counter.
  Integer, Pointer              :: indexID                            ! Identifier for a GRIB file index.
  Integer                       :: gribID                             ! Identifier for a GRIB message.
  Integer                       :: Status                             ! Status code from GRIB_API subroutines.
  Character(MaxCharLength)      :: ErrMessage                         ! Error message from GRIB_API subroutines.
  Integer                       :: nX                                 !} Size of field encoded in GRIB message.
  Integer                       :: nY                                 !}
  Integer                       :: numberOfValues                     ! Size of data array in GRIB message.
  Character(MaxCharLength)      :: StepUnits                          ! Time unit in GRIB header.
  Character(MaxCharLength)      :: StepType                           !} F/C step information in GRIB header.
  Integer                       :: StartStep                          !}
  Integer                       :: EndStep                            !}

  ! $$ Saved variables should really be returned.

  Error              = .false.
  InconsistentHeader = .false.
  Missing            = .false.

  ! GRIB Edition 1 met file.
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

  ! PP met file.
  Else If (FileType .CIEq. 'PP') Then

    If (NewFile) MissingL = .false.

    ! Read header.
    If (.not.MissingL) Then
      Read (Unit = Unit, Err = 9, End = 9) IHeader, RHeader
    End If

    ! Field missing.
    If (IHeader(42) /= FieldCode .or. (ThreeD .and. IHeader(33) /= iFieldLevel)) Then

      MissingL = .true.
      Missing  = .true.

    ! Field not missing.
    Else

      MissingL = .false.
      ! $$ More checks of headers in Name II & PP met files (inc Time)? Possibly make
      ! optional (controlled from
      ! input) or just (non-fatal) error or warning (currently inconsistent header generates an error)

      ! Read field.
      Read (Unit = Unit, Err = 9, End = 9) LocalField
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

      ! Checks for MetDefn being consistent with header information.
      ! Notes:
      ! In the header the horizontal data grid can be defined indirectly as a rectangular subset of another
      ! grid (called here the main grid).
      ! Header( 1) = } Year, month, day, hour, and minute of time 1. Header(1) <= 0 indicates a missing field.
      ! Header( 2) = } Generally time 1 is the time the data is valid for - but see PP file documentation for
      ! Header( 3) = } other possibilities.
      ! Header( 4) = }
      ! Header( 5) = }
      ! Header( 7) = ] Year, month, day, hour, and minute of time 2. Generally time 2 is the time of the data
      ! Header( 8) = ] used to calculte the field (e.g. the analysis time from which a forecast is produced) -
      ! Header( 9) = ] but see PP file documentation for other possibilities.
      ! Header(10) = ]
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
      ! $$ these checks assume that grid is lat-long (degrees) and Header values are always in
      ! $$ lat-long degrees.
      ! $$ Try to encapsulate some of the coord changes in the coord module?

      ! Note we don't check stash and level number for missing fields since these are not
      ! always correct (at least level isn't).

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

End Module AncillaryMetModule

! Module:  Radar Rainfall Met Module

Module RadarMetModule
!.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use CommonMetModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: RadarMetDefn_
Public  :: RadarMetDefn2_
Public  :: RadarMetDefns_
Public  :: RadarMet_
Public  :: InitRadarMetDefn
Public  :: InitRadarMetDefn2
Public  :: InitRadarMetDefns
Public  :: InitRadarMet
Public  :: AddRadarMetDefn
Public  :: AddRadarMetDefn2
Public  :: FindRadarMetDefnIndex
Public  :: FindRadarMetDefn2Index
Public  :: SetUpRadarMet_CoordsEtc
Public  :: PrepareForUpdateRadarMet ! Prepares for updating an instance of the radar met module.
Public  :: UpdateRadarMet

!-------------------------------------------------------------------------------------------------------------

! Radar met field names.
Integer,                  Parameter :: nFieldNames2d = 2
Character(MaxCharLength), Parameter :: FieldNames2d(nFieldNames2d) =          &
                                       (/                                     &
                                         'dummy                     ',        & ! 1
                                         'precipitation rate (mm/hr)'         & ! 2
                                       /)
! nFieldNames2d            :: Number of 2-d (horizontal) fields.
! FieldNames2d             :: Names of 2-d (horizontal) fields.

! Codes for radar met fields:
Integer, Parameter :: F_Dummy           =  1 !} Codes giving the location of each field
Integer, Parameter :: F_Ppt             =  2 !} in the array FieldNames2d.

! Skipping of radar data files.
Integer, Parameter :: MaxSkip = 1  ! maximal number of files we allow to be skipped

!-------------------------------------------------------------------------------------------------------------

Type :: RadarMetDefn_ ! A radar met definition (part 1 - basic definition).
  Character(MaxCharLength)     :: Name
  Character(MaxCharLength)     :: BinaryFormat
  Character(MaxCharLength)     :: FileType
  Character(MaxCharLength)     :: Prefix
  Character(MaxCharLength)     :: Suffixs(MaxRadarMetDefn2s)
  Logical                      :: DayPerFile
  Type(Time_)                  :: T0
  Type(Time_)                  :: Dt
  Logical                      :: NextPrecip
  Integer                      :: nMetDefn2Names
  Character(MaxCharLength)     :: MetDefn2Names(MaxRadarMetDefn2s)
  Character(MaxCharLength)     :: HGridName
  ! Name           :: Name of met definition (part 1 - basic definition).
  ! BinaryFormat   :: Binary format of radar met files.
  ! FileType       :: Type of radar met files (currently only 'nimrod' format is supported).
  ! Prefix         :: Prefix for names of met files.
  ! Suffixs        :: Suffixes for names of met files.
  ! DayPerFile     :: Indicates each met file contains a whole day rather than a single time.
  ! T0             :: Reference time for the first radar met field.
  ! Dt             :: Time interval between fields.
  ! NextPrecip     :: Indicates its best to use the next time step for precipitation rather than interpolating
  !                   in time.
  ! nMetDefn2Names :: Number of MetDefn2Names (i.e. part 2s) and suffixes for file names.
  ! MetDefn2Names  :: Names of the radar met definitions (part 2 - met file structure definitions).
  ! HGridName      :: Name of main horizontal grid.
End Type RadarMetDefn_

!-------------------------------------------------------------------------------------------------------------

Type :: RadarMetDefn2_ ! A radar met definition (part 2 - met file structure definition).
                       ! Together with part 1 this provides information on the nature of
                       ! the radar met module met data.
  Character(MaxCharLength) :: Name
  Integer                  :: nFields
  Character(MaxCharLength) :: FieldNames(MaxRadarMetFields)
  Integer                  :: LowestLevels(MaxRadarMetFields)
  Integer                  :: HighestLevels(MaxRadarMetFields)
  Logical                  :: Top(MaxRadarMetFields)
  Integer                  :: FieldCodes(MaxRadarMetFields)
  Logical                  :: ThreeD(MaxRadarMetFields)
  Logical                  :: Total(MaxRadarMetFields)
  ! Name          :: Name of met definition (part 2 - met file structure definition).
  ! nFields       :: Number of fields in the met file.
  ! FieldNames    :: For each of the fields in the met file, the name of the corresponding NAME III field.
  ! LowestLevels  :: For each of the fields in the met file, the lowest NAME III model level.
  ! HighestLevels :: For each of the fields in the met file, the highest NAME III model level. Defined only if
  !                  Top is false.
  ! Top           :: For each of the fields in the met file, indicates that the highest NAME III model level
  !                  is the top level of the appropriate grid.
  ! FieldCodes    :: For each of the fields in the met file, the radar field code 
  ! ThreeD        :: For each of the fields in the met file, indicates that the field is part of a 3-d radar
  !                  field. Note this means it has a level index associated with it (which needn't match
  !                  the NAME III level index).
  ! Total         :: For each of the fields in the met file, indicates that the field is a total (dyn + conv) 
  !                  cloud field (rather than dynamic) $$ for consistency with NWP met but not used for radar met

  ! Currently radar met is only intended to read a single 2-d precipitation field. However the code has
  ! been designed to give more flexible capabilities similar to those for the NWP met module.

  ! Note that here we distinguish between: (i) a 'NAME III field' - a complete 2-d or 3-d field as used within
  ! NAME III, (ii) a 'Radar field' - a complete 2-d  field supplying the data, and (iii) a 'field in the 
  ! met file' - a 2-d (horizontal) field or a consecutive set of (horizontal) sections of a 3-d field.

  ! Here 'consecutive' means in the sense of their location in the met file if the format requires this and in
  ! the sense of their (NAME III) level numbers (with level numbers increasing).

  ! Although here we refer to a field in the met file as a single entity, in fact the met file formats have a
  ! basic unit of data which is a single 2-d field or a single 2-d section of a 3-d field. All formats support
  ! the possibility of some of these basic units being missing in the met file.

  ! If the radar field is 2-d then LowestLevels and HighestLevels must be equal. If the NAME III field is 2-D,
  ! then LowestLevels and HighestLevels must equal 1.

  ! $$ add tests that this is so. could set to zero instead? blank in met defn file? But note wish to read
  ! several levels into dummy field (see above).

  ! The same 2-d radar field (or the same 2-d section of a 3-d radar field) should not appear twice in the met
  ! file where this can be avoided. This is enforced by the following tests when the met files are read:
  !   Sequential read in order not necessarily aligned with RadarMetDefn2_: check the same radar field code /
  !   level combination doesn't occur twice in the file.
  !   Sequential read in same order as RadarMetDefn2_: here there are some situations which may require the same
  !   field to appear twice and so no checks are made.

End Type RadarMetDefn2_

!-------------------------------------------------------------------------------------------------------------

Type :: RadarMetDefns_ ! A collection of radar met definitions.
  Integer              :: nRadarMetDefns                    ! Number of radar met definitions  (part 1
                                                            ! - basic definitions).
  Type(RadarMetDefn_)  :: RadarMetDefns(MaxRadarMetDefns)   ! Radar met definitions (part 1
                                                            ! - basic definitions).
  Integer              :: nRadarMetDefn2s                   ! Number of radar met definitions (part 2
                                                            ! - met file structure definitions).
  Type(RadarMetDefn2_) :: RadarMetDefn2s(MaxRadarMetDefn2s) ! Radar met definitions (part 2
                                                            ! - met file structure definitions).
End Type RadarMetDefns_

!-------------------------------------------------------------------------------------------------------------

Type :: RadarMet_  ! Information describing the state of an instance of the radar met module.

  ! CommonMet_:
  Type(CommonMet_) :: C   ! The part of the met state common to all met modules

  ! Validity flags for each time level.
  Logical :: Valid(2)     ! Indicates whether met module is valid for each time level

  ! Input variables:
  Type(RadarMetDefn_)          :: MetDefn                      ! Radar met definition (part 1 - basic
                                                               ! definition).
  Type(RadarMetDefn2_)         :: MetDefn2s(MaxRadarMetDefn2s) ! Radar met definition (part 2 - met file
                                                               ! structure definitions).
  Type(ShortTime_)             :: Dt                           ! Time interval between fields.
  Character(MaxFileNameLength) :: MetFolder                    ! Location of radar rainfall files
  Character(MaxFileNameLength) :: RestoreMetScript             ! File containing script for restoring met data.
                                                               ! Blank indicates met files will not be restored.
  Logical                      :: DeleteMet                    ! Indicates met files will be deleted after use.

  ! Quantities derived from MetDefn2s and Grids:
  Integer :: iField(MaxRadarMetDefn2s,MaxRadarMetFields) ! Indices in FieldNames3d and FieldNames2d of
                                                         ! fields in the met file.
  Integer :: k1(MaxRadarMetDefn2s,MaxRadarMetFields)     ! Lowest NAME III model level of each of the fields in
                                                         ! the met file (1 for 2-D fields).
  Integer :: k2(MaxRadarMetDefn2s,MaxRadarMetFields)     ! Highest NAME III model level of each of the fields in
                                                         ! the met file (1 for 2-D fields).
                                                         ! $$ Note that radar met is only intended to read
                                                         !    2-d fields at present, although there is no
                                                         !    reason why 3-d fields could not be read too.

  ! Grid and coord indices:
  Integer :: iHGrid         ! Index of horizontal grid used by radar fields.
  Integer :: iHCoord        ! Index of coord system of horizontal grid used by radar fields.
  
  ! Interpolation coefficients:
  Type(GHCoeffs_) :: GHCoeffs    ! Coefficients for interpolating from horizontal grid to horizontal grid
  Type(GHCoeffs_) :: GHCoeffsdX  ! Coefficients for interpolating x derivatives from grid to grid
  Type(GHCoeffs_) :: GHCoeffsdY  ! Coefficients for interpolating y derivatives from grid to grid

  ! Allocatable arrays:
  Real(Std), Pointer :: Ppt(:,:,:)   ! Precipitation array

  ! Information on allocatable arrays:
  Logical          :: SpaceAllocated                !  Indicates whether the arrays have been allocated.
  Integer          :: NewData                       !} Indices of latest and latest but one sets of met data.
  Integer          :: OldData                       !} 0 indicates data invalid.
  Type(ShortTime_) :: NewTime                       !] Time of latest and latest but one sets of met data.
  Type(ShortTime_) :: OldTime                       !]
  Logical          :: FieldPresent2d(nFieldNames2d) ! Indicates which 2-d fields have been successfully read in.

  ! Miscellaneous:
  Character(MaxFileNameLength) :: OldMetFiles(MaxRadarMetDefn2s) ! Name of the last opened met files. Blank
                                                                 ! indicates undefined.

End Type RadarMet_

!-------------------------------------------------------------------------------------------------------------

Interface Operator(==) ! Equality of radar met definitions.
  Module Procedure RadarMetDefnEq
  Module Procedure RadarMetDefn2Eq
End Interface

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitRadarMetDefn(                                            &
           Name, BinaryFormat, FileType, Prefix, Suffixs, DayPerFile, &
           T0, Dt, NextPrecip, MetDefn2Names, HGridName               &
         )                                                            &
Result(RadarMetDefn)
! Initialises a radar met definition (part 1 - basic definition).

  Implicit None
  ! Argument list:
  Character(*),             Intent(In) :: Name             ! Name of met definition (part 1 - basic definition).
  Character(*),             Intent(In) :: BinaryFormat     ! Binary format of radar met files.
  Character(*),             Intent(In) :: FileType         ! Type of radar met files (currently only 'nimrod'
                                                           ! format is supported).
  Character(*),             Intent(In) :: Prefix           ! Prefix for names of met files.
  Character(MaxCharLength), Intent(In) :: Suffixs(:)       ! Suffixes for names of met files.
                                                           ! $$ changed from * to support Intel compiler on Cray
  Logical,                  Intent(In) :: DayPerFile       ! Indicates each met file contains a whole day
                                                           ! rather than a single time.
  Type(Time_),              Intent(In) :: T0               ! Reference time for the first radar met field.
  Type(Time_),              Intent(In) :: Dt               ! Time interval between fields.
  Logical,                  Intent(In) :: NextPrecip       ! Indicates its best to use the next time step for
                                                           ! precipitation rather than interpolating in time.
  Character(MaxCharLength), Intent(In) :: MetDefn2Names(:) ! Names of radar met definitions (part 2 - met file
                                                           ! structure definitions).
                                                           ! $$ changed from * to support Intel compiler on Cray
  Character(*),             Intent(In) :: HGridName        ! Name of main horizontal grid.

  ! Function result:
  Type(RadarMetDefn_) :: RadarMetDefn ! The radar met definition (part 1 - basic definition).

  ! Locals:
  Integer :: iMetDefn2Names    ! Loop index.

  ! Name.
  If (Len_Trim(Name) == 0) Then
    Call Message('ERROR in InitRadarMetDefn: Name is blank', 3)
  End If
  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message(                                             &
           'ERROR in InitRadarMetDefn: Name is given as "' // &
           Trim(Name)                                      // &
           '" and is too long',                               &
           3                                                  &
         )
  End If
  RadarMetDefn%Name = Name

  ! Binary Format.
  If (Len_Trim(BinaryFormat) > MaxCharLength) Then
    Call Message(                                                      &
           'ERROR in InitRadarMetDefn: Binary Format is given as "' // &
           Trim(BinaryFormat)                                       // &
           '" and is too long',                                        &
           3                                                           &
         )
  End If
  RadarMetDefn%BinaryFormat = BinaryFormat

  ! File Type.
  If (Len_Trim(FileType) == 0) Then
    Call Message('ERROR in InitRadarMetDefn: File Type is blank', 3)
  End If
  If (Len_Trim(FileType) > MaxCharLength) Then
    Call Message(                                                   &
           'ERROR in InitRadarMetDefn:  File Type is given as "' // &
           Trim(FileType)                                        // &
           '" and is too long',                                     &
           3                                                        &
         )
  End If

  If (.not.(FileType .CIEq. 'Nimrod')) Then
    Call Message(                                                            &
           'ERROR in InitRadarMetDefn: radar met file type is given as "' // &
           Trim(FileType)                                                 // &
           '" and has not been recognised',                                  &
           3                                                                 &
         )
  End If
  RadarMetDefn%FileType = FileType

  ! File Prefix.
  If (Len_Trim(Prefix) > MaxCharLength) Then
    Call Message(                                               &
           'ERROR in InitRadarMetDefn: Prefix is given as "' // &
           Trim(Prefix)                                      // &
           '" and is too long',                                 &
           3                                                    &
         )
  End If
  RadarMetDefn%Prefix = Prefix

  ! File Suffixes and Met File Structure Definitions.
  If (Size(Suffixs) /= Size(MetDefn2Names)) Then
    Call Message(                                                                     &
           'ERROR in InitRadarMetDefn: the number of supplied Suffixs should be '  // &
           'the same as the number of Radar Met File Structure Definitions',          &
           3                                                                          &
         )
  End If

  Do iMetDefn2Names = 1,Size(MetDefn2Names)
    If (Len_Trim(Suffixs(iMetDefn2Names)) > MaxCharLength) Then
      Call Message(                                   &
             'ERROR in InitRadarMetDefn: Element ' // &
             Trim(Int2Char(iMetDefn2Names))        // &
             ' of Suffixs is given as "'           // &
             Trim(Suffixs(iMetDefn2Names))         // &
             '" and is too long',                     &
             3                                        &
           )
    End If

    If (Len_Trim(MetDefn2Names(iMetDefn2Names)) == 0) Then
      Call Message(                                                 &
             'ERROR in InitRadarMetDefn: the name of the '       // &
             Trim(Int2Char(iMetDefn2Names))                      // &
             '-th Radar Met File Structure Definition is blank',    &
             3                                                      &
           )
    End If
    If (Len_Trim(MetDefn2Names(iMetDefn2Names)) > MaxCharLength) Then
      Call Message(                                                      &
             'ERROR in InitRadarMetDefn: the name of the '            // &
             Trim(Int2Char(iMetDefn2Names))                           // &
             '-th Radar Met File Structure Definition is given as "'  // &
             Trim(MetDefn2Names(iMetDefn2Names))                      // &
             '" and is too long',                                        &
             3                                                           &
           )
    End If
  End Do
  RadarMetDefn%nMetDefn2Names = Size(MetDefn2Names)

  RadarMetDefn%Suffixs      (1:RadarMetDefn%nMetDefn2Names) = Suffixs(:)
  RadarMetDefn%MetDefn2Names(1:RadarMetDefn%nMetDefn2Names) = MetDefn2Names(:)

  ! H-Grid Name.
  If (Len_Trim(HGridName) == 0) Then
    Call Message('ERROR in InitRadarMetDefn: HGridName is blank', 3)
  End If
  If (Len_Trim(HGridName) > MaxCharLength) Then
    Call Message(                                                  &
           'ERROR in InitRadarMetDefn: HGridName is given as "' // &
           Trim(HGridName)                                      // &
           '" and is too long',                                    &
           3                                                       &
         )
  End If
  RadarMetDefn%HGridName = HGridName

  ! DayPerFile.
  RadarMetDefn%DayPerFile = DayPerFile

  ! NextPrecip.
  RadarMetDefn%NextPrecip = NextPrecip

  ! T0.
  RadarMetDefn%T0 = T0

  ! Dt.
  If (Dt <= ZeroTime()) Then
    Call Message('ERROR in InitRadarMetDefn: radar update time step is not positive', 3)
  End If
  RadarMetDefn%Dt = Dt

End Function InitRadarMetDefn

!-------------------------------------------------------------------------------------------------------------

Function InitRadarMetDefn2(                                                              &
           Name, FieldNames, LowestLevels, HighestLevels, Top, FieldCodes, ThreeD, Total &
         )                                                                               &
Result(RadarMetDefn2)
! Initialises a radar met definition (part 2 - met file structure definition).

  Implicit None
  ! Argument list:
  Character(*),              Intent(In) :: Name
  Character(MaxTokenLength), Intent(In) :: FieldNames(:)   ! $$ changed from * to support Intel compiler on Cray
  Integer,                   Intent(In) :: LowestLevels(:)
  Integer,                   Intent(In) :: HighestLevels(:)
  Logical,                   Intent(In) :: Top(:)
  Integer,                   Intent(In) :: FieldCodes(:)
  Logical,                   Intent(In) :: ThreeD(:)
  Logical,                   Intent(In) :: Total(:)
  ! Name          :: Name of met definition (part 2 - met file structure definition).
  ! FieldNames    :: For each of the fields in the met file, the name of the corresponding NAME III field.
  ! LowestLevels  :: For each of the fields in the met file, the lowest NAME III model level.
  ! HighestLevels :: For each of the fields in the met file, the highest NAME III model level. Needs to be
  !                  defined only if Top is false.
  ! Top           :: For each of the fields in the met file, indicates that the highest NAME III model level
  !                  is the top level of the appropriate grid.
  ! FieldCodes    :: For each of the fields in the met file, the radar field code.
  ! ThreeD        :: For each of the fields in the met file, indicates that the field is part of a 3-d radar
  !                  field. Note this means it has an level index associated with it (which needn't match
  !                  the NAME III level index).
  ! Total         :: For each of the fields in the met file, indicates that the field is a total (dyn + conv) 
  !                  cloud field $$ for consistency with NWP met but not used for radar met
  ! Function result:
  Type(RadarMetDefn2_) :: RadarMetDefn2 ! The radar met definition (part 2 - met file structure definition).
  ! Locals:
  Integer :: i       ! Loop index.
  Integer :: nFields ! Number of fields in the met file structure definition.

  ! Name.
  If (Len_Trim(Name) == 0) Then
    Call Message(                                                             &
           'ERROR in initialising a Radar Met File Structure Definition: ' // &
           'the Radar Met File Structure Definition name is blank',           &
           3                                                                  &
         )
  End If
  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message(                                                             &
           'ERROR in initialising a Radar Met File Structure Definition: ' // &
           'the Radar Met File Structure Definition name is given as "'    // &
           Trim(Name)                                                      // &
           '" and is too long',                                               &
           3                                                                  &
         )
  End If
  RadarMetDefn2%Name = Name

  ! Check on number of fields.
  nFields = Size(FieldNames)
  If (nFields > MaxRadarMetFields) Then
    Call Message(                                                 &
           'ERROR in initialising a Radar Met File Structure ' // &
           'Definition: the number of fields is too large',       &
           3                                                      &
         )
  End If
  RadarMetDefn2%nFields = nFields

  If (                                    &
    (Size(LowestLevels)  /= nFields) .or. &
    (Size(HighestLevels) /= nFields) .or. &
    (Size(Top)           /= nFields) .or. &
    (Size(FieldCodes)    /= nFields) .or. &
    (Size(ThreeD)        /= nFields) .or. &
    (Size(Total)         /= nFields)      &
  ) Then
    Call Message(                                                            &
           'UNEXPECTED ERROR in initialising a Radar Met File Structure ' // &
           'Definition: the array sizes are inconsistent',                   &
           3                                                                 &
         )
  End If

  ! Field names.
  Do i = 1, nFields
    If (Len_Trim(FieldNames(i)) == 0) Then
      Call Message(                                                             &
             'ERROR in initialising a Radar Met File Structure Definition: ' // &
             'the met file structure definition "'                           // &
             Trim(Name)                                                      // &
             '" contains an empty field name (field '                        // &
             Trim(Int2Char(i))                                               // &
             ')',                                                               &
             3                                                                  &
           )
    End If
    If (Len_Trim(FieldNames(i)) > MaxCharLength) Then
      Call Message(                                                             &
             'ERROR in initialising a Radar Met File Structure Definition: ' // &
             'the met file structure definition "'                           // &
             Trim(Name)                                                      // &
             '" contains an invalid field name (field '                      // &
             Trim(Int2Char(i))                                               // &
             ' has the name "'                                               // &
             Trim(FieldNames(i))                                             // &
             '" which is too long)',                                            &
             3                                                                  &
           )
    End If
  End Do

  RadarMetDefn2%FieldNames   (1:nFields) = FieldNames   (1:nFields)
  RadarMetDefn2%LowestLevels (1:nFields) = LowestLevels (1:nFields)
  RadarMetDefn2%HighestLevels(1:nFields) = HighestLevels(1:nFields)
  RadarMetDefn2%Top          (1:nFields) = Top          (1:nFields)
  RadarMetDefn2%FieldCodes   (1:nFields) = FieldCodes   (1:nFields)
  RadarMetDefn2%ThreeD       (1:nFields) = ThreeD       (1:nFields)
  RadarMetDefn2%Total        (1:nFields) = Total        (1:nFields)

End Function InitRadarMetDefn2

!-------------------------------------------------------------------------------------------------------------

Function InitRadarMetDefns() Result(MetDefns)
! Initialises a collection of radar met definitions.

  Implicit None
  ! Function result:
  Type(RadarMetDefns_) :: MetDefns ! Initialised collection of met definitions.

  MetDefns%nRadarMetDefns  = 0
  MetDefns%nRadarMetDefn2s = 0

End Function InitRadarMetDefns

!-------------------------------------------------------------------------------------------------------------

Function InitRadarMet(                               &
           MetName,                                  &
           FixedMet,                                 &
           UpdateOnDemand,                           &
           RadarMetDefn, RadarMetDefn2s,             &
           MetFolder,                                &
           RestoreMetScript,                         &
           DeleteMet                                 &
         )
! Initialises an instance of the radar rainfall met module.

  Implicit None
  ! Argument list:
  Character(*),         Intent(In) :: MetName           ! Name of met module instance.
  Logical,              Intent(In) :: FixedMet
  Logical,              Intent(In) :: UpdateOnDemand    ! Indicates the met module instance is to be updated
                                                        ! using update-on-demand.
  Type(RadarMetDefn_),  Intent(In) :: RadarMetDefn      ! Radar met definition (part 1 - basic definition)
  Type(RadarMetDefn2_), Intent(In) :: RadarMetDefn2s(:) ! Radar met definition (part 2 - met file structure
                                                        ! definitions).
  Character(*),         Intent(In) :: MetFolder         ! Location of radar rainfall files.
  Character(*),         Intent(In) :: RestoreMetScript  ! File containing script for restoring met data.
                                                        ! Blank indicates met files will not be restored.
  Logical,              Intent(In) :: DeleteMet         ! Indicates met files will be deleted after use.
  ! Function result:
  Type(RadarMet_) InitRadarMet ! Initialised instance of the radar rainfall met module.
  ! Locals:
  Type(RadarMet_) RadarMet     ! Local copy of instance of the radar rainfall met module.

  ! Initialise the common part of the met module instance
  RadarMet%C = InitCommonMet('Radar Met', MetName, FixedMet, UpdateOnDemand)

  ! Met definition
  RadarMet%MetDefn = RadarMetDefn
  RadarMet%MetDefn2s(1:RadarMet%MetDefn%nMetDefn2Names) = RadarMetDefn2s(:)

  ! Met folder name
  If (Len_Trim(MetFolder) > MaxFileNameLength) Then
    Call Message(                                                                      &
           'ERROR in InitRadarMet: the name of the radar met folder is given as "'  // &
           Trim(MetFolder)                                                          // &
           '" and is too long',                                                        &
           3                                                                           &
         )
  End If
  RadarMet%MetFolder = MetFolder

  ! Restore met script
  If (Len_Trim(RestoreMetScript) > MaxFileNameLength) Then
    Call Message(                                                                        &
           'ERROR in InitRadarMet: the name of the restore met script is given as "'  // &
           Trim(RestoreMetScript)                                                     // &
           '" and is too long',                                                          &
           3                                                                             &
         )
  End If
  RadarMet%RestoreMetScript = RestoreMetScript
  RadarMet%DeleteMet        = DeleteMet

  ! Initialise met update step
  RadarMet%Dt = Time2ShortTime(RadarMetDefn%Dt)
  
  ! Initialise met module validity
  RadarMet%Valid(:) = .false.

  ! Precipitation array is initially unallocated
  RadarMet%SpaceAllocated = .false.

  RadarMet%OldData = 0
  RadarMet%NewData = 0

  InitRadarMet = RadarMet

End Function InitRadarMet

!-------------------------------------------------------------------------------------------------------------

Function RadarMetDefnEq(RadarMetDefn1, RadarMetDefn2)
! Tests for equality of radar met definitions (part 1 - basic definition).

  Implicit None
  ! Argument list:
  Type(RadarMetDefn_), Intent(In) :: RadarMetDefn1 !} The two radar met definitions (part 1
  Type(RadarMetDefn_), Intent(In) :: RadarMetDefn2 !} - basic definition).
  ! Function result:
  Logical :: RadarMetDefnEq ! Indicates if radar met definitions are equal.
  ! Locals:
  Integer :: i ! Loop index.

  RadarMetDefnEq = (RadarMetDefn1%Name          .CIEq. RadarMetDefn2%Name          ) .and. &
                   (RadarMetDefn1%BinaryFormat  .CIEq. RadarMetDefn2%BinaryFormat  ) .and. &
                   (RadarMetDefn1%FileType      .CIEq. RadarMetDefn2%FileType      ) .and. &
                   (RadarMetDefn1%Prefix        .CIEq. RadarMetDefn2%Prefix        ) .and. &
                   (RadarMetDefn1%DayPerFile    .eqv.  RadarMetDefn2%DayPerFile    ) .and. &
                   (RadarMetDefn1%T0              ==   RadarMetDefn2%T0            ) .and. &
                   (RadarMetDefn1%Dt              ==   RadarMetDefn2%Dt            ) .and. &
                   (RadarMetDefn1%NextPrecip    .eqv.  RadarMetDefn2%NextPrecip    ) .and. &
                   (RadarMetDefn1%HGridName     .CIEq. RadarMetDefn2%HGridName     ) .and. &
                   (RadarMetDefn1%nMetDefn2Names  ==   RadarMetDefn2%nMetDefn2Names)

  Do i = 1, Min(RadarMetDefn1%nMetDefn2Names, RadarMetDefn2%nMetDefn2Names)
    RadarMetDefnEq = RadarMetDefnEq                                                       .and. &
                   (RadarMetDefn1%Suffixs(i)       .CIEq. RadarMetDefn2%Suffixs(i)      ) .and. &
                   (RadarMetDefn1%MetDefn2Names(i) .CIEq. RadarMetDefn2%MetDefn2Names(i))
  End Do

End Function RadarMetDefnEq

!-------------------------------------------------------------------------------------------------------------

Function RadarMetDefn2Eq(RadarMetDefn1, RadarMetDefn2)
! Tests for equality of radar met definitions (part 2 - met file structure definition).

  Implicit None
  ! Argument list:
  Type(RadarMetDefn2_), Intent(In) :: RadarMetDefn1 !} The two radar met definitions (part 2
  Type(RadarMetDefn2_), Intent(In) :: RadarMetDefn2 !} - met file structure definition).
  ! Function result:
  Logical :: RadarMetDefn2Eq ! Indicates if radar met definitions are equal.
  ! Locals:
  Integer :: i ! Loop index.

  RadarMetDefn2Eq = (RadarMetDefn1%Name    .CIEq. RadarMetDefn2%Name   ) .and. &
                    (RadarMetDefn1%nFields   ==   RadarMetDefn2%nFields)

  Do i = 1, Min(RadarMetDefn1%nFields, RadarMetDefn2%nFields)
    RadarMetDefn2Eq =                                                              &
      RadarMetDefn2Eq                                                        .and. &
      (RadarMetDefn1%FieldNames(i)    .CIEq. RadarMetDefn2%FieldNames(i)   ) .and. &
      (RadarMetDefn1%LowestLevels(i)    ==   RadarMetDefn2%LowestLevels(i) ) .and. &
      (RadarMetDefn1%HighestLevels(i)   ==   RadarMetDefn2%HighestLevels(i)) .and. &
      (RadarMetDefn1%Top(i)           .eqv.  RadarMetDefn2%Top(i)          ) .and. &
      (RadarMetDefn1%FieldCodes(i)      ==   RadarMetDefn2%FieldCodes(i)   ) .and. &
      (RadarMetDefn1%ThreeD(i)        .eqv.  RadarMetDefn2%ThreeD(i)       ) .and. &
      (RadarMetDefn1%Total(i)         .eqv.  RadarMetDefn2%Total(i)        )
  End Do

End Function RadarMetDefn2Eq

!-------------------------------------------------------------------------------------------------------------

Subroutine AddRadarMetDefn(RadarMetDefn, MetDefns)
! Adds a radar met definition (part 1 - basic definition) to a collection of radar met definitions.

  Implicit None
  ! Argument list:
  Type(RadarMetDefn_),    Intent(In)    :: RadarMetDefn ! The radar met definition (part 1 - basic definition)
  Type(RadarMetDefns_),   Intent(InOut) :: MetDefns     ! The collection of radar met definitions
  ! Locals:
  Integer :: i ! Loop index

  Do i = 1, MetDefns%nRadarMetDefns
    If (RadarMetDefn%Name .CIEq. MetDefns%RadarMetDefns(i)%Name) Then
      If (RadarMetDefn == MetDefns%RadarMetDefns(i)) Then
        Return
      Else
        Call Message(                                                            &
               'ERROR in adding Radar Met Definition "'                       // &
               Trim(RadarMetDefn%Name)                                        // &
               '": a different definition with the same name already exists',    &
               3                                                                 &
             )
      End If
    End If
  End Do

  If (MetDefns%nRadarMetDefns >= MaxRadarMetDefns) Then
    Call Message(                                         &
           'ERROR in adding a Radar Met Definition: '  // &
           'there are too many Radar Met Definitions ',   &
           3                                              &
         )
  End If

  MetDefns%nRadarMetDefns                         = MetDefns%nRadarMetDefns + 1
  MetDefns%RadarMetDefns(MetDefns%nRadarMetDefns) = RadarMetDefn

End Subroutine AddRadarMetDefn

!-------------------------------------------------------------------------------------------------------------

Subroutine AddRadarMetDefn2(RadarMetDefn2, MetDefns)
! Adds a radar met definition (part 2 - met file structure definition) to a collection of radar met definitions.

  Implicit None
  ! Argument list:
  Type(RadarMetDefn2_), Intent(In)    :: RadarMetDefn2 ! The radar met definition (part 2 - met file structure
                                                       ! definition)
  Type(RadarMetDefns_), Intent(InOut) :: MetDefns      ! The collection of radar met definitions
  ! Locals:
  Integer :: i ! Loop index

  Do i = 1, MetDefns%nRadarMetDefn2s
    If (RadarMetDefn2%Name .CIEq. MetDefns%RadarMetDefn2s(i)%Name) Then
      If (RadarMetDefn2 == MetDefns%RadarMetDefn2s(i)) Then
        Return
      Else
        Call Message(                                                            &
               'ERROR in adding Radar Met File Structure Definition "'        // &
               Trim(RadarMetDefn2%Name)                                       // &
               '": a different Radar Met File Structure Definition with the ' // &
               'same name already exists',                                       &
               3                                                                 &
             )
      End If
    End If
  End Do

  If (MetDefns%nRadarMetDefn2s >= MaxRadarMetDefn2s) Then
    Call Message(                                                       &
           'ERROR in adding a Radar Met File Structure Definition: ' // &
           'there are too many Radar Met File Structure Definitions',   &
           3                                                            &
         )
  End If

  MetDefns%nRadarMetDefn2s                          = MetDefns%nRadarMetDefn2s + 1
  MetDefns%RadarMetDefn2s(MetDefns%nRadarMetDefn2s) = RadarMetDefn2

End Subroutine AddRadarMetDefn2

!-------------------------------------------------------------------------------------------------------------

Function FindRadarMetDefnIndex(Name, MetDefns)
! Finds the index of a radar met definition (part 1 - basic definition).

  Implicit None
  ! Argument list:
  Character(*),         Intent(In) :: Name     ! Name of radar met definition (part 1 - basic definition).
  Type(RadarMetDefns_), Intent(In) :: MetDefns ! Collection of radar met definitions.
  ! Function result:
  Integer :: FindRadarMetDefnIndex ! Index of radar met definition (part 1 - basic definition).
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, MetDefns%nRadarMetDefns
    If (Name .CIEq. MetDefns%RadarMetDefns(i)%Name) Then
      FindRadarMetDefnIndex = i
      Return
    End If
  End Do

  Call Message(                                                                &
         'FATAL ERROR in FindRadarMetDefnIndex: the Radar Met Definition "' // &
         Trim(Name)                                                         // &
         '" has not been found',                                               &
         3                                                                     &
       )

End Function FindRadarMetDefnIndex

!-------------------------------------------------------------------------------------------------------------

Function FindRadarMetDefn2Index(Name, MetDefns)
! Finds the index of a radar met definition (part 2 - met file structure definition).

  Implicit None
  ! Argument list:
  Character(*),         Intent(In) :: Name     ! Name of radar met definition (part 2 - met file
                                               ! structure definition).
  Type(RadarMetDefns_), Intent(In) :: MetDefns ! Collection of radar met definitions.
  ! Function result:
  Integer :: FindRadarMetDefn2Index ! Index of radar met definition (part 2 - met file structure definition).
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, MetDefns%nRadarMetDefn2s
    If (Name .CIEq. MetDefns%RadarMetDefn2s(i)%Name) Then
      FindRadarMetDefn2Index = i
      Return
    End If
  End Do

  Call Message(                                                                                &
         'FATAL ERROR in FindRadarMetDefn2Index: the Radar Met File Structure Definition "' // &
         Trim(Name)                                                                         // &
         '" has not been found',                                                               &
         3                                                                                     &
       )

End Function FindRadarMetDefn2Index

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpRadarMet_CoordsEtc(EtaDefns, Coords, Grids, MetEnsembleSize, RadarMet)
! Sets up RadarMet using information from EtaDefns, Coords and Grids. ! $$ not quite true now MetEnsembleSize
! added.

  Implicit None
  ! Argument list:
  Type(EtaDefns_), Intent(In)           :: EtaDefns        ! Collection of eta definitions.
  Type(Coords_),   Intent(In)           :: Coords          ! Collection of coord systems.
  Type(Grids_),    Intent(In),   Target :: Grids           ! Collection of grids.
  Integer,         Intent(In)           :: MetEnsembleSize ! Size of the met ensemble (i.e. number of met
                                                           ! realisations).
  Type(RadarMet_), Intent(InOut)        :: RadarMet        ! State of a radar met module instance.
  ! Locals:
  Character(MaxCharLength) :: ZCoordName  ! Name of the vertical coord system
  Type(HGrid_),    Pointer :: HGrid       ! Abbreviation for horizontal grid
  Logical                  :: Error       ! Error flag indicates fatal error which however will not be made
                                          ! fatal just yet in order to detect other possible errors.
  Integer                  :: i           !} Loop indices.
  Integer                  :: j           !}
  Integer                  :: k           !}

  ! Initialise Error.
  Error = .false.

  ! Find indices of coords and grids.
  RadarMet%iHGrid  = FindHGridIndex(RadarMet%MetDefn%HGridName, Grids)
  RadarMet%iHCoord = FindHCoordIndex(Grids%HGrids(RadarMet%iHGrid)%HCoordName, Coords)

  ZCoordName = ''
  
  ! Add coord systems to RadarMet%C.
  Call AddCoordsToCommonMet(                                &
         1, (/ Grids%HGrids(RadarMet%iHGrid)%HCoordName /), &
         0, (/ ZCoordName /),                               &
         RadarMet%C                                         &
       )

  ! Check number of coord systems in RadarMet.
  If (RadarMet%C%nHCoords /= 1) Then
    Call Message(                                               &
           'UNEXPECTED ERROR in SetUpRadarMet_CoordsEtc: '   // &
           'error in numbering of horizontal coord systems ' // &
           'used by the radar met module instance "'         // &
           Trim(RadarMet%C%MetName)                          // &
           '"',                                                 &
           4                                                    &
         )
  End If
  If (RadarMet%C%nZCoords /= 0) Then
    Call Message(                                             &
           'UNEXPECTED ERROR in SetUpRadarMet_CoordsEtc: ' // &
           'error in numbering of vertical coord systems ' // &
           'used by the radar met module instance "'       // &
           Trim(RadarMet%C%MetName)                        // &
           '"',                                               &
           4                                                  &
         )
  End If

  ! Set up abbreviation for horizontal grid.
  HGrid => Grids%HGrids(RadarMet%iHGrid)

  ! Check horizontal grids are regular with non-zero spacing and more than one point in each direction.
  If (                                           &
    HGrid%Unstructured .or.                      &
    HGrid%Variable     .or.                      &
    HGrid%dX == 0.0    .or. HGrid%dY == 0.0 .or. &
    HGrid%nX <= 1      .or. HGrid%nY <= 1        &
    ) Then
    Call Message(                                                             &
           'ERROR in SetUpRadarMet_CoordsEtc: The horizontal grid named "' // &
           Trim(HGrid%Name)                                                // &
           '" which is specified as "H-Grid" in the '                      // &
           '"Radar Met Definition" named "'                                // &
           Trim(RadarMet%MetDefn%Name)                                     // &
           '" has variable or zero spacing or less than 2 points in one '  // &
           '(at least) direction',                                            &
           2                                                                  &
         )
    Error = .true.
  End If

  ! Fatal error.
  If (Error) Then
    Call Message(                                        &
           'FATAL ERROR: the grids specified in the ' // &
           '"Radar Met Definition" named "'           // &
           Trim(RadarMet%MetDefn%Name)                // &
           '" are not useable',                          &
           3                                             &
         )
  End If

  ! Calculate interpolation coefficients.
  Call InitGHCoeffs(RadarMet%GHCoeffs  )
  Call InitGHCoeffs(RadarMet%GHCoeffsdX)
  Call InitGHCoeffs(RadarMet%GHCoeffsdY)
  Call GetGHCoeffs(HGrid , HGrid, ' ', RadarMet%GHCoeffs  )
  Call GetGHCoeffs(HGrid , HGrid, 'X', RadarMet%GHCoeffsdX)
  Call GetGHCoeffs(HGrid , HGrid, 'Y', RadarMet%GHCoeffsdY)

  Do k = 1, RadarMet%MetDefn%nMetDefn2Names
    ! Determine the field index and dimension of each field in each RadarMetDefn2.
    Do i = 1, RadarMet%MetDefn2s(k)%nFields
      RadarMet%iField(k,i) = 0
      Do j = 1, nFieldNames2d
        If (RadarMet%MetDefn2s(k)%FieldNames(i) .CIEq. FieldNames2d(j)) Then
          RadarMet%iField(k,i) = j
        End If
      End Do

      If (RadarMet%iField(k,i) == 0) Then
        Call Message(                                          &
               'ERROR: the met file structure definition "' // &
               Trim(RadarMet%MetDefn2s(k)%Name)             // &
               '" contains an invalid field name "'         // &
               Trim(RadarMet%MetDefn2s(k)%FieldNames(i))    // &
               '"',                                            &
               2                                               &
             )
        Error = .true.
      End If
    End Do

    ! Set up k1 and k2.
    Do i = 1, RadarMet%MetDefn2s(k)%nFields
      RadarMet%k1(k,i) = 1
      RadarMet%k2(k,i) = 1
    End Do

  ! Fatal error.
    If (Error) Then
      Call Message(                                                &
             'FATAL ERROR: the met file structure definition "' // &
             Trim(RadarMet%MetDefn2s(k)%Name)                   // &
             '" is not useable',                                   &
             3                                                     &
           )
    End If
  End Do

End Subroutine SetUpRadarMet_CoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForUpdateRadarMet( &
             Coords, Grids,          &
             iCase, iMetCase,        &
             Time,                   &
             TValid, UpdateNow,      &
             RadarMet,               &
             Units                   &
           )
! Prepares for updating an instance of the radar met module.

! This routine must set TValid and UpdateNow but must not alter the validity of the radar met module instance.
! $$ this comment should be with generic definition.

! $$ iCase and iMetCase are passed in to this routine for consistency with NWP routines but are not currently used.

  Implicit None
  ! Argument List:
  Type(Coords_),   Intent(In)           :: Coords
  Type(Grids_),    Intent(In),   Target :: Grids
  Integer,         Intent(In)           :: iCase
  Integer,         Intent(In)           :: iMetCase
  Type(Time_),     Intent(In)           :: Time
  Type(Time_),     Intent(Out)          :: TValid
  Logical,         Intent(Out)          :: UpdateNow
  Type(RadarMet_), Intent(InOut)        :: RadarMet
  Type(Units_),    Intent(InOut)        :: Units
  ! Coords    :: Collection of coord systems.
  ! Grids     :: Collection of grids.
  ! iCase     :: Number of case.
  ! iMetCase  :: Number of the met realisation in the met ensemble.
  ! Time      :: Time for which the radar met module instance might be updated.
  ! TValid    :: Earliest time that the validity (overall or for any single attribute) of the met module
  !              instance might change, assuming the met module instance is updated now. The value is that
  !              determined at the end of this routine (the actual time may be later).
  ! UpdateNow :: Indicates the met module instance must be updated now (even if update-on-demand is
  !              specified). If set, TValid need not be set to any particular time.
  ! RadarMet  :: State of a radar met module instance.
  ! Units     :: Collection of information on input/output unit numbers.
  ! Locals:
  Type(Time_) :: MetTime !

  MetTime = Round(Time, RadarMet%MetDefn%T0, RadarMet%MetDefn%Dt, Up = .false.)

  MetTime = MetTime + RadarMet%MetDefn%Dt

  If (RadarMet%C%FixedMet) Then
    TValid = InfFutureTime()
  Else
    TValid = MetTime
  End If
  UpdateNow = .false.

End Subroutine PrepareForUpdateRadarMet

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateRadarMet(    &
             Coords, Grids,   &
             iCase, iMetCase, &
             Time,            &
             RadarMet,        &
             Units            &
           )
! Updates an instance of the radar rainfall met module.

! $$ iCase and iMetCase are passed in to this routine for consistency with NWP routines but are not currently used.

  Implicit None
  ! Argument list:
  Type(Coords_),   Intent(In)    :: Coords   ! Collection of coord systems.
  Type(Grids_),    Intent(In)    :: Grids    ! Collection of grids.
  Integer,         Intent(In)    :: iCase    ! Number of case.
  Integer,         Intent(In)    :: iMetCase ! Number of the met realisation in the met ensemble.
  Type(Time_),     Intent(In)    :: Time     ! Time for which the radar rainfall met module instance
                                             ! is to be updated.
  Type(RadarMet_), Intent(InOut) :: RadarMet ! State of a radar rainfall met module instance.
  Type(Units_),    Intent(InOut) :: Units    ! Collection of information on input/output unit numbers.
  ! Locals:
  Type(Time_)      :: MetTime1        !} Times and short times of the radar met data required.
  Type(Time_)      :: MetTime2        !}
  Type(ShortTime_) :: SMetTime1       !}
  Type(ShortTime_) :: SMetTime2       !}
  Logical          :: ReadOneTimeOnly ! Indicates that only one set of radar met data is to be read.
  Logical          :: DoIfBlock       ! Indicates if-block is to be executed to read two sets of radar data.
  Logical          :: Error           ! Error flag.
  Integer          :: nSkipped        ! Number of skipped radar met data files
  Logical          :: Skip            ! Indicates that a file should be skipped
  Logical          :: AllowSkip       ! Indicated whether we allow skipping of a particular file

  Error = .false.

  ! Allocate precipitation array.
  If (.not.RadarMet%SpaceAllocated) Then
    Allocate ( RadarMet%Ppt(Grids%HGrids(RadarMet%iHGrid)%nX, Grids%HGrids(RadarMet%iHGrid)%nY, 2) )
    RadarMet%SpaceAllocated = .true.
  End If

  ! Calculate times of radar rainfall data.
  MetTime1 = Round(Time, RadarMet%MetDefn%T0, RadarMet%MetDefn%Dt, Up = .false.)
             ! $$ Could name files using a time zone - currently files named in UTC.
  MetTime2  = MetTime1 + RadarMet%MetDefn%Dt
  SMetTime1 = Time2ShortTime(MetTime1)
  SMetTime2 = Time2ShortTime(MetTime2)

  ! Only read one set of radar data for FixedMet runs where the time of fixed met coincides with a data time.
  ReadOneTimeOnly = (Time == MetTime1) .and. RadarMet%C%FixedMet

  If (RadarMet%NewData == 0) Then
    DoIfBlock = .true.
  Else
    DoIfBlock = SMetTime1 /= RadarMet%NewTime
  End If

  If (DoIfBlock) Then
    
    ! Mark met data in Old time level as invalid.
    RadarMet%OldData = 0

    ! Read Time1 met data into New time level.
    RadarMet%NewData = 1
    RadarMet%NewTime = SMetTime1

    ! Do not allow skipping of first radar data
    AllowSkip = .false.

    Call ReadRadarMet(MetTime1, iCase, iMetCase, Coords, Grids, AllowSkip, Error, Skip, RadarMet, Units)

    If (Error) Go To 9

    ! Set met data as valid at this time
    RadarMet%Valid(RadarMet%NewData) = .true.

    ! Now set the initial value of RadarMet%OldData appropriately
    RadarMet%OldData = 3

  End If

  ! Move Time1 met data into Old time level.
  RadarMet%OldData = RadarMet%NewData
  RadarMet%OldTime = RadarMet%NewTime

  ! Swap Time2 met data pointer
  RadarMet%NewData = 3 - RadarMet%OldData
  RadarMet%NewTime = SMetTime2

  ! Skip second read of radar data for FixedMet runs where the time of fixed met coincides with a data time
  ! Here the rainfall at the New time level is set equal to the rainfall at the Old time level.
  If (ReadOneTimeOnly) Then

    RadarMet%Ppt(:, :, RadarMet%NewData) = RadarMet%Ppt(:, :, RadarMet%OldData)

  Else

    nSkipped = 0

    Do
      ! $$ Do not allow skipping MetData if we use parallel IO
      AllowSkip = (nSkipped < MaxSkip)

      Call ReadRadarMet(MetTime2, iCase, iMetCase, Coords, Grids, AllowSkip, Error, Skip, RadarMet, Units)

      If (Skip) Then
        nSkipped = nSkipped + 1
        ! Increase the time by Dt and try to read next file
        MetTime2  = MetTime2 + RadarMet%MetDefn%Dt
        SMetTime2 = Time2ShortTime(MetTime2)
        RadarMet%Dt      = RadarMet%Dt + Time2ShortTime(RadarMet%MetDefn%Dt)
        RadarMet%NewTime = SMetTime2
      Else
        Exit
      End If

    End Do

    If (Error) Go To 9

    ! Set met data as valid at this time
    RadarMet%Valid(RadarMet%NewData) = .true.

  End If

  ! Check validity of radar met data
  If (ReadOneTimeOnly) Then
    Error = .not.RadarMet%Valid(1)
  Else
    Error = .not.(RadarMet%Valid(RadarMet%OldData) .and. RadarMet%Valid(RadarMet%NewData))
  End If

  If (Error) Go To 9

  ! Set validity flags.
  RadarMet%C%Valid = .true.
  If (RadarMet%C%FixedMet) Then
    RadarMet%C%TValid = InfFutureTime()
  Else
    RadarMet%C%TValid = MetTime2
  End If

  Return

  ! Errors.
9 Continue

  RadarMet%Valid(RadarMet%NewData) = .false.
  RadarMet%NewData = 0
  RadarMet%C%Valid = .false.
  If (RadarMet%C%FixedMet) Then
    RadarMet%C%TValid = InfFutureTime()
  Else
    RadarMet%C%TValid = MetTime2
  End If

End Subroutine UpdateRadarMet

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadRadarMet(Time, iCase, iMetCase, Coords, Grids, AllowSkip, Error, Skip, M, Units, IdxIn)

! Reads radar met data for a single time.

! Note this routine reads data into the M%New time level of the data arrays.
!
! Non-fatal errors are issued for file opening errors, read errors and inconsistent file headers. Warnings are
! issued for missing data which is expected in the met definition, except where more than one 2-d field or 2-d
! (horizontal) section of a 3-d field in the met definition can fulfill a given role, and where one but not
! all of these are present.

! $$ Currently if module is invalid, those values which were obtained cannot be used later for persistence fix
!    ups. Is this sensible? (it may be, but may need to be reviewed with attribute dep validity).

! $$ iCase and iMetCase are passed in to this routine for consistency with NWP routines but are not currently used.

  Implicit None
  ! Argument List:
  Type(Time_),     Intent(In)           :: Time      ! Time for which radar met data is to be read.
  Integer,         Intent(In)           :: iCase     ! Number of case.
  Integer,         Intent(In)           :: iMetCase  ! Index of the met realisation in the met ensemble.
  Type(Coords_),   Intent(In),   Target :: Coords    ! Collection of coord systems.
  Type(Grids_),    Intent(In),   Target :: Grids     ! Collection of grids.
  Logical,         Intent(In)           :: AllowSkip ! Allow skipping of this update time
  Logical,         Intent(Out)          :: Error     ! Indicates an error occured in this routine.
  Logical,         Intent(InOut)        :: Skip      ! Indicates that the radar met data could not be found
                                                     ! and should be skipped
  Type(RadarMet_), Intent(InOut)        :: M         ! State of a radar met module instance.
  Type(Units_),    Intent(InOut)        :: Units     ! Collection of information on input/output unit numbers.
  Integer,         Intent(In), Optional :: IdxIn

  ! Local types:
  Type :: F2_ ! 2-d (horizontal) field (with time dimension too).
    Real(Std), Pointer :: P(:,:,:) ! 2-d (horizontal) field.
  End Type F2_

  ! Locals:
  Real(Std), Allocatable, Target       :: Dummy(:,:,:)

  Type(F2_)                            :: F2(nFieldNames2d)
  Integer                              :: nT
  Character(MaxFileNameLength)         :: MetFolder
  Character(MaxFileNameLength)         :: MetFiles(MaxRadarMetDefn2s)
  Character(MaxFileNameLength)         :: MetFileFail
  Integer                              :: iMetFile
  Logical                              :: Exist
  Integer                              :: IOStat
  Logical                              :: InconsistentHeader
  Logical                              :: InconsistentValidityTime
  Logical                              :: Missing
  Logical                              :: NewFile
  Integer                              :: MetFileUnits(MaxRadarMetDefn2s)
  Type(HCoord_),               Pointer :: HCoord
  Type(HGrid_),                Pointer :: HGrid
  Integer                              :: i
  Integer                              :: j
  Integer                              :: k
  Integer                              :: iT
  Integer                              :: iField
  Integer                              :: iAttrib
  Character(MaxCharLength)             :: CharTime
  Integer                              :: idx
  Integer                              :: nDummyBuffers
  Integer                              :: Status
  Character(MaxCharLength)             :: ErrMessage
  ! F2                 :: 2-d (horizontal) fields used as alias for fields in M and for field Dummy.
  ! nT                 :: Index of the required time in the met file.
  ! MetFolder          :: Directory containing radar files
  ! MetFiles           :: Met files.
  ! MetFileFail        :: Met file which cannot be opened
  ! iMetFile           :: Loop variable
  ! Exist              :: Indicates that a file exists.
  ! IOStat             :: Error code for open statement.
  ! InconsistentHeader :: Indicates the header for the field is inconsistent with what is expected.
  ! InconsistentValidityTime :: Indicates the validity time in the header is inconsistent with what is
  !                             expected (note this does not necessarily imply the met file is not valid).
  ! Missing            :: Indicates a field is missing from the met file.
  ! NewFile            :: Indicates a newly opened file is being considered.
  ! MetFileUnits       :: Input/output unit numbers.
  ! HCoord             :} Abbreviations for coords, grids.
  ! HGrid              :}
  ! i                  :] Loop indices.
  ! j                  :]
  ! k                  :]
  ! iT                 :]
  ! iField             :]
  ! iAttrib            :: Index of attribute.
  ! CharTime           :: Character version of Time.
  ! idx                :: 
  ! Status             :: Status code from GRIB_API subroutines.
  ! ErrMessage         :: Error message from GRIB_API subroutines.

  nDummyBuffers = 2
  Allocate(Dummy(Grids%HGrids(M%iHGrid)%nX,Grids%HGrids(M%iHGrid)%nY,nDummyBuffers))

  If (Present(IdxIn)) Then
    idx=IdxIn
  Else
    idx=M%NewData
  End If

  ! 1) Set up error flag and abbreviations for grids.

  ! Set error flag to false.
  Error = .false.
  ! Set skip flag to false
  Skip = .false.

  ! Set up abbreviations for coords and grids.
  HCoord  => Coords%HCoords(M%iHCoord)
  HGrid   => Grids%HGrids(M%iHGrid)

  ! 2) Sort out radar met files.

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
!        '.'                               // &
        Trim(M%MetDefn%Suffixs(iMetFile))
    Else
      nT = 1
      MetFiles(iMetFile) =                   &
        Trim(M%MetDefn%Prefix)            // &
        Trim(FileNameTime(Time))          // &
!        '.'                               // &
        Trim(M%MetDefn%Suffixs(iMetFile))
    End If

    ! Delete previous file.
    If (M%DeleteMet .and. M%OldMetFiles(iMetFile) /= MetFiles(iMetFile) .and. M%OldMetFiles(iMetFile) /= ' ') Then
      Inquire(File = Trim(M%MetFolder) // Trim(M%OldMetFiles(iMetFile)), Exist = Exist)

      If (Exist) Then
        Call GetNewUnit(MetFileUnits(iMetFile), Units) !$$ to use openfile,
                                                       ! need to add return-on-error option to openfile, with options
                                                       ! to give ERROR or WARNING instead of FATAL ERROR.
                                                       ! Until then should trap error here and give error or warning message.
        Open (                                                                                &
          Unit   = MetFileUnits(iMetFile),                                                    &
          File   = Trim(ConvertFileName(Trim(M%MetFolder) // Trim(M%OldMetFiles(iMetFile)))), &
          Status = 'Old',                                                                     &
          Form   = 'Unformatted',                                                             &
          IOStat = IOStat                                                                     &
        ) 

        Call CloseUnit(Unit = MetFileUnits(iMetFile), Units = Units, DeleteFile = .true.)
      End If
    End If

    ! Restore met file.
    If (M%RestoreMetScript /= ' ') Then
      Inquire(File = Trim(M%MetFolder) // Trim(MetFiles(iMetFile)), Exist = Exist)
      If (.not.Exist) Then
        Call SubmitSystemCommand(                               &
               Trim(ConvertFileName(M%RestoreMetScript))     // &
               ' '                                           // &
               Trim(M%MetFolder)                             // &
               ' '                                           // &
               Trim(MetFiles(iMetFile))                         &
           )
      End If
    End If

    ! Store met file name.
    M%OldMetFiles(iMetFile) = MetFiles(iMetFile)

    ! Open met file 
    Call GetNewUnit(MetFileUnits(iMetFile), Units) !$$ use openfile once have sorted out how to handle 'convert ='
#   ifdef UseConvert
      Open (                                                                            &
        Unit    = MetFileUnits(iMetFile),                                               &
        File    = Trim(ConvertFileName(Trim(M%MetFolder) // Trim(MetFiles(iMetFile)))), &
        Status  = 'Old',                                                                &
        Form    = 'Unformatted',                                                        &
        Action  = 'Read',                                                               &
        Convert = Trim(M%MetDefn%BinaryFormat),                                         &
        IOStat  = IOStat                                                                &
      )
#   else
      Open (                                                                            &
        Unit   = MetFileUnits(iMetFile),                                                &
        File   = Trim(ConvertFileName(Trim(M%MetFolder)  // Trim(MetFiles(iMetFile)))), &
        Status = 'Old',                                                                 &
        Form   = 'Unformatted',                                                         &
        Action = 'Read',                                                                &
        IOStat = IOStat                                                                 &
      )
#   endif
    If (IOStat /= 0) Then
      Call CloseUnit(MetFileUnits(iMetFile), Units)
      If (.not. AllowSkip) Then
        Call Message(                                                                   &
               'ERROR: Cannot open met file '                                        // &
               Trim(ConvertFileName(Trim(M%MetFolder) // Trim(MetFiles(iMetFile)))),    &
               2                                                                        &
             )
      End If
      Skip = AllowSkip
      Error = .true.
      MetFileFail = MetFiles(iMetFile)
      Go To 9
    End If

  End Do

  ! 3) Set up equivalences to enable looping over fields and initialise FieldPresent2d and NewFile.

  F2(F_Dummy          )%P => Dummy
  F2(F_Ppt            )%P => M%Ppt

  M%FieldPresent2d(:) = .false.

  LoopOverMetFiles: Do iMetFile = 1, M%MetDefn%nMetDefn2Names

    NewFile = .true.

    ! 4) Read data (note missing data is ignored at times which don't correspond to the time of interest).

    ! $$ Note currently all files are read in prescribed order.

    ! Reading fields in the order they are in M%MetDefn2s(iMetFile).
    If (.true.) Then

      LoopOverTimesInFile: Do iT = 1, nT

        LoopOverFieldsAtGivenTime: Do iField = 1, M%MetDefn2s(iMetFile)%nFields

          ! 2-d fields.
          Call ReadRadarField(                                                                   &
                 NewFile,                                                                        &
                 M%MetDefn%FileType,                                                             &
                 Trim(M%MetFolder) // Trim(MetFiles(iMetFile)),                                  &
                 MetFileUnits(iMetFile),                                                         &
                 M%MetDefn2s(iMetFile)%FieldCodes(iField),                                       &
                 Time, HCoord, HGrid,                                                            &
                 F2(M%iField(iMetFile,iField))%P(:, :, idx),                                     &
                 Error, InconsistentHeader, InconsistentValidityTime, Missing                    &
               )

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
                     'given in the Radar Met Definition '                                               // &
                     Trim(M%MetDefn%Name)                                                               // &
                     '.',                                                                                  &
                     MessageControls    = GlobalMessageControls,                                           &
                     MessageControlName = 'Inconsistent radar met header',                                 &
                     ErrorCode          = 2                                                                &
                   ) ! $$ Set Error = .true. and goto 9? or downgrade to warning
            End If

            ! Inconsistent validity time.
            If (InconsistentValidityTime .and. iT == nT) Then
              Call ControlledMessage(                                                                      &
                     'The field '                                                                       // &
                     Trim(M%MetDefn2s(iMetFile)%FieldNames(iField))                                     // &
                     ' from met file '                                                                  // &
                     Trim(MetFiles(iMetFile))                                                           // &
                     ' for the '                                                                        // &
                     Trim(Int2Char(iT)) // Int2Ordinal(iT)                                              // &
                     ' time level in the file has a header validity time which is inconsistent with '   // &
                     'the time in the met file name',                                                      &
                     MessageControls    = GlobalMessageControls,                                           &
                     MessageControlName = 'Inconsistent radar validity time',                              &
                     ErrorCode          = 2                                                                &
                   ) ! $$ Set Error = .true. and goto 9? or downgrade to warning
            End If

            ! FieldPresent2d.
            If (iT == nT .and. .not.Missing) M%FieldPresent2d(M%iField(iMetFile,iField)) = .true.

            ! NewFile.
            NewFile = .false.

        End Do LoopOverFieldsAtGivenTime

      End Do LoopOverTimesInFile

    ! Reading fields in the order they are in the file.
    Else

    End If

    ! Close met file 
  
    Call CloseUnit(MetFileUnits(iMetFile), Units)

    ! 5) Warnings for missing fields which are expected to be present according to MetDefn2s. Note this is done
    ! here rather than as each field is read, to avoid producing warnings when more than one possible radar field
    ! can be used for a given NAME III field and some but not all of these radar fields are missing.

    Do iField = 1, M%MetDefn2s(iMetFile)%nFields

      ! 2-d fields.

      If (.not. M%FieldPresent2d(M%iField(iMetFile,iField))) Then
        Call ControlledMessage(                                      &
               'Field '                                           // &
               Trim(M%MetDefn2s(iMetFile)%FieldNames(iField))     // &
               ' missing from met file '                          // &
               Trim(MetFiles(iMetFile))                           // &
               ' for the '                                        // &
               Trim(Int2Char(nT)) // Int2Ordinal(nT)              // &
               ' time level in the file',                            &
               MessageControls    = GlobalMessageControls,           &
               MessageControlName = 'Missing radar met field',       &
               ErrorCode          = 1                                &
             )

        ! If all required fields are not present in the file then return as error.
        Error = .true.
        Skip = AllowSkip
        Go To 9

      End If

    End Do

  End Do LoopOverMetFiles

  ! 6) Message.
  Do iMetFile = 1, M%MetDefn%nMetDefn2Names
    Call Message(                                           &
           'Radar met data for '                         // &
           Trim(Time2Char(Time, .false., 0, .true.))     // &
           ' read from file '                            // &
           Trim(M%MetFolder) // Trim(MetFiles(iMetFile))    &
         )
  End Do

  Return

9 Continue

  If (Skip) Then
    Call Message(                                             &
           'WARNING: Cannot find met data for met module ' // &
           Trim(M%C%MetModName)                            // &
           '.'                                             // &
           Trim(M%C%MetName)                               // &
           ' skipping to next update step',                   &
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

End Subroutine ReadRadarMet

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadRadarField(                                                 &
             NewFile,                                                      &
             FileType, FileName, Unit,                                     &
             FieldCode, Time,                                              &
             HCoord, HGrid,                                                &
             Field,                                                        &
             Error, InconsistentHeader, InconsistentValidityTime, Missing  &
           )
! Reads in a 2-d (horizontal) field of radar data
! $$ Currently supports only 'nimrod' type files.
!    Support for GRIB could be added if/when the postprocessing system moves to GRIB.

  Implicit None
  ! Argument List:
  Logical,       Intent(In)  :: NewFile                  ! Indicates a newly opened file is being considered.
  Character(*),  Intent(In)  :: FileType                 ! Type of radar met file ('nimrod' format).
  Character(*),  Intent(In)  :: FileName                 ! Name of radar met file (not currently used).
  Integer,       Intent(In)  :: Unit                     ! Input/output unit number of radar met file.
  Integer,       Intent(In)  :: FieldCode                ! Radar code of the field.
  Type(Time_),   Intent(In)  :: Time                     ! Time of the field.
  Type(HCoord_), Intent(In)  :: HCoord                   ! Coord system in which HGrid is defined.
  Type(HGrid_),  Intent(In)  :: HGrid                    ! Grid on which the field is defined.
  Real(Std),     Intent(Out) :: Field(:, :)              ! 2-d (horizontal) field.
  Logical,       Intent(Out) :: Error                    ! Indicates a read error.
  Logical,       Intent(Out) :: InconsistentHeader       ! Indicates the header for the field is inconsistent
                                                         ! with what is expected.
  Logical,       Intent(Out) :: InconsistentValidityTime ! Indicates the validity time from the field
                                                         ! header is inconsistent with what is expected
                                                         ! (note this does not necessarily imply the
                                                         ! met file is not valid, e.g. an existing met
                                                         ! file may have been duplicated using
                                                         ! a new time in the filename).
  Logical,       Intent(Out) :: Missing                  ! Indicates a missing field.

  ! Locals:
  Integer(I16)                  :: LocalField(HGrid%nX, HGrid%nY)   ! Local copy of met field.
  Integer                       :: i, j                             ! Loop indices.
  Integer,                 Save :: IOStat                           ! Error code for read statement.
  !$OMP THREADPRIVATE(IOStat)

  ! Headers for Nimrod formatted files
  Integer(I16),            Save :: IHeader(31)         ! General integer header record
  !$OMP THREADPRIVATE(IHeader)
  Real(Std),               Save :: RHeader(28)         ! General real header record
  !$OMP THREADPRIVATE(RHeader)
  Integer(I16),            Save :: IHeaderSpecific(51) ! Data specific integer header record
  !$OMP THREADPRIVATE(IHeaderSpecific)
  Real(Std),               Save :: RHeaderSpecific(45) ! Data specific real header record
  !$OMP THREADPRIVATE(RHeaderSpecific)
  Character*56,            Save :: CHeader             ! Character header entries
  !$OMP THREADPRIVATE(CHeader)

  Logical,                 Save :: MissingL            ! Indicates the last expected field was missing
                                                       ! (so do not read another header section).
  !$OMP THREADPRIVATE(MissingL)

  Real(Std), Parameter :: CoordTolerance = 0.001       ! tolerance (in degrees) when checking coord values

  ! $$ Saved variables should really be returned.

  Error                    = .false.
  InconsistentHeader       = .false.
  InconsistentValidityTime = .false.
  Missing                  = .false.

  ! Nimrod met file.
  If (FileType .CIEq. 'Nimrod') Then

    If (NewFile) MissingL = .false.

    ! Read header.
    If (.not.MissingL) Then
      Read (Unit = Unit, Err = 9, IOStat = IOStat)                       &
            IHeader, RHeader, RHeaderSpecific, CHeader, IHeaderSpecific
    End If

    ! Field missing.
    If (IHeader(19) /= FieldCode .or. IOStat < 0) Then

      MissingL = .true.
      Missing  = .true.

    ! Field present.
    Else

      MissingL = .false.

      ! Check header. Note code for the field (stash code) and level number already checked when testing for
      ! missing fields above.
      ! Notes:
      ! IHeader( 1) = } Year, month, day, hour, minute and second of the Validity Time of the field.
      ! IHeader( 2) = } For data with a time-period of validity (e.g. precip accumulation over one hour),
      ! IHeader( 3) = } this is the end of the time-period.
      ! IHeader( 4) = }
      ! IHeader( 5) = }
      ! IHeader( 6) = }
      ! IHeader(12) = Data type ( = 0 for real, 1 for integer, 2 for byte).
      ! IHeader(13) = Number of bytes for each data element (1, 2, or 4).
      ! IHeader(15) = Horizontal grid type (0=NG, 1=lat/lon, 2=space view, 3=polar stereographic,
      !                                     4=UTM32 (EuroPP), 5=Rotated Lat Lon, 6=other).
      ! IHeader(16) = Number of rows in field.
      ! IHeader(17) = Number of columns in field.
      ! IHeader(19) = Field code number.
      ! IHeader(26) = Period of interest for accumulation, average or probability (minutes).
      ! IHeader(31) = Time averaging (LBPROC) ( = 128 for accumulation or average).
      ! RHeader( 3) = Northing or latitude or start line of first row of data (metres for NG, degrees for PS grids).
      ! RHeader( 4) = Interval between rows (metres or degrees).
      ! RHeader( 5) = Easting or longitude or start pixel of first point of first row of data (metres or degrees).
      ! RHeader( 6) = Interval between columns (metres or degrees).

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

      ! Check grid size nX, nY
      InconsistentHeader = InconsistentHeader      .or. &
                           HGrid%nX /= IHeader(17) .or. &
                           HGrid%nY /= IHeader(16)

      ! $$ More checks of headers in nimrod files?

      ! Read field.
      If (IHeader(12) == 1 .and. IHeader(13) == 2) Then ! 2Byte Integer data type
        Read (Unit = Unit, Err = 9, End = 9, IOStat=IOStat) LocalField

        Do i = 1, HGrid%nX
          Do j = 1, HGrid%nY
            If (LocalField(i,j) < 0.0) Then
              LocalField(i,j) = 0.0
              Call Message('Warning: negative precipitation - set to zero', 1)
            End If
          End Do
        End Do

        ! Nimrod fields are scaled by a factor of 32.0
        Field(:,:) = LocalField(:,:) * (1.0/32.0)

      Else

        Call Message('FATAL ERROR in ReadRadarField: Unexpected data type', 4)

      End If

    End If

  Else

    Call Message('UNEXPECTED FATAL ERROR in ReadRadarField: unsupported met file type', 4)

  End If

  Return

  ! Read error.
9 Continue

  Error = .true.

End Subroutine ReadRadarField

!-------------------------------------------------------------------------------------------------------------

End Module RadarMetModule

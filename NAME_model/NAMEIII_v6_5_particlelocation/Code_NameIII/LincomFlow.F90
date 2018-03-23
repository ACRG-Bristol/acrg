! Module:  Lincom Flow Module

Module LINCOMFlowModule

! This is a flow module designed interface with LINCOM, a linear flow model

! It uses input from a single instance of the LINCOM met module.
!
! It does not use any input from other flow modules (and hence uses no update subset).
!
! It uses the topography and surface pressure information supplied by the LINCOM met
! module.
!
! It has attributes Update, Flow.
!
! The status of the optional flow module features is as follows (where ???? here =
! LINCOM):
!     SetUpCoordsEtc_????Flow routine     - used
!     SetUp????Flow_MetsCoordsEtc routine - used
!     ????FlowCoordIndices routine        - used
!     Get????Flow routine                 - used (with flow memory argument)
!     ????ReflectCoordIndices             - used
!     ????Reflect routine                 - used.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowAndFlowProfileModule
Use MetsModule
Use CommonFlowModule

Use CommonMetModule ! $$
Use NWPMetModule    ! $$
! This module also uses LINCOM, but this is referenced via a system call, not via module routines.

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: LINCOMFlow_                   ! The state of an LINCOM flow module instance.
Public  :: InitLINCOMFlow                ! Initialises the state of an LINCOM flow module
                                         ! instance.
Public  :: SetUpCoordsEtc_LINCOMFlow     ! Sets up Coords and Grids by adding any extra
                                         ! coords and grids which LINCOMFlow wants to
                                         ! define.
Public  :: SetUpLINCOMFlow_MetsCoordsEtc ! Sets up LINCOMFlow using information from
                                         ! EtaDefns, Coords, Grids and Mets.
Public  :: PrepareForUpdateLINCOMFlow    ! Prepares for updating an instance of the LINCOM flow module.
Public  :: LINCOMFlowReqs                ! Specifies what information the flow module
                                         ! instance wants from the other flow module
                                         ! instances.
Public  :: UpdateLINCOMFlow              ! Updates an instance of the LINCOM flow module.
Public  :: LINCOMFlowCoordIndices        ! Returns indices of coord systems in which
                                         ! positions need to be specified when using
                                         ! GetLINCOMFlow.
Public  :: GetLINCOMFlow                 ! Gets flow information from an LINCOM flow module
                                         ! instance.
Public  :: LINCOMReflectCoordIndices     ! Returns indices of coord systems in which
                                         ! positions need to be specified when using
                                         ! LINCOMReflect.
Public  :: LINCOMReflect                 ! Reflects position of particle or puff centroid
                                         ! using an LINCOM flow module instance.

!-------------------------------------------------------------------------------------------------------------

Type :: LINCOMFlow_ ! The state of an LINCOM flow module instance.
  Type(CommonFlow_)            :: C         ! Part of flow state common to all flow modules.

  Character(MaxFileNameLength) ::  TerrainFile   ! Name of terrain data file.
  Character(MaxFileNameLength) ::  RoughnessFile ! Name of roughness data file.
  Character(MaxCharLength)     ::  HGrid         ! Horizontal data LINCOM grid.
  Character(MaxCharLength)     ::  ZGrid         ! Vertical data LINCOM grid.
  Logical             ::  TerrainPert            ! Flag to calc pert caused by terrain.
  Logical             ::  RoughnessPert          ! Flag to calc pert caused by roughness.
  Logical             ::  BackgroundWindAdded    ! Flag to add background wind speed to pert.
  Integer             ::  nXCompGrid             ! X computational grid must be 2^n
                                                 ! and larger than HGrid%X.
  Integer             ::  nYCompGrid             ! Y computational grid must be 2^n
                                                 ! and larger than HGrid%Y.
  Integer             ::  iHCoord
  Integer             ::  iZCoord
  Integer             ::  iHGrid
  Integer             ::  iZGrid
  Integer             ::  iHGridFlowForLINCOM   ! grid used to transfer flow into LINCOM
  Real(Std), Pointer  ::  Velocity(:,:,:,:,:)
  Real(Std), Pointer  ::  UStar(:,:,:)
  Real(Std), Pointer  ::  Z0(:,:,:)
  Real(Std), Pointer  ::  Topog(:,:,:)
  Logical             ::  SpaceAllocated
  Type(Flow_)         ::  Flow(1,1,250,2)        ! $$ Better to allocate these arays
  Type(ProfileData_)  ::  ProfileData(1,1,250,2) ! in LINCOMFlowReqs. In particular 250
                                                 ! is a guess - the old MaxZPoints before
                                                 ! this was removed.
  Character(MaxFileNameLength) :: LincomExe              ! Name of Lincom executable.
  Integer          :: New
  Integer          :: Old
  Type(ShortTime_) :: OldTime
  Type(ShortTime_) :: NewTime
  ! New            :} Indices of latest and latest but one sets of met data. 0
  ! Old            :} indicates data invalid because no data read in, errors occurred
  !                   during reading the data, or data is for the wrong time.
  ! OldTime        :] Time of latest and latest but one sets of met data.
  ! NewTime        :]
End Type LINCOMFlow_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitLINCOMFlow(     &
           FlowName,         &
           TerrainFile,      &
           RoughnessFile,    &
           HGrid,            &
           ZGrid,            &
           TerrainPert,      &
           RoughnessPert,    &
           nXCompGrid,       &
           nYCompGrid,       &
           LincomExe,        &
           DomainName,       &
           UpdateSubsetName, &
           FixedMet,         &
           UpdateOnDemand    &
         )                   &
         Result (LINCOMFlow)

! Initialises the state of an LINCOM flow module instance.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: FlowName         ! Name of flow module instance.
  Character(*), Intent(In) :: TerrainFile      ! Name of terrain file.
  Character(*), Intent(In) :: RoughnessFile    ! Name of roughness file.
  Character(*), Intent(In) :: HGrid
  Character(*), Intent(In) :: ZGrid
  Logical,      Intent(In) :: TerrainPert      ! Calculate perturbation caused by terrain.
  Logical,      Intent(In) :: RoughnessPert    ! Calculate perturbation caused by roughness.
  Integer,      Intent(In) :: nXCompGrid       ! X Computational grid 2^n.
  Integer,      Intent(In) :: nYCompGrid       ! Y Computational grid 2^n.
  Character(*), Intent(In) :: LincomExe        ! Name of Lincom executable.
  Character(*), Intent(In) :: DomainName       ! Name of the domain of the flow module instance.
  Character(*), Intent(In) :: UpdateSubsetName !
  Logical,      Intent(In) :: FixedMet         ! Indicates that the met is fixed.
  Logical,      Intent(In) :: UpdateOnDemand   ! Indicates the flow module instance is to be updated using
                                               ! update-on-demand.
  ! Function result:
  Type(LINCOMFlow_)        :: LINCOMFlow       ! Initialised state of an LINCOM flow
                                               ! module instance.

!  If (.not. IsCalendar()) Then
!    Call Message(                                           &
!           'Error in InitLINCOMFlow: The LINCOM flow '   // &
!           'module only works in a calendar time frame',    &
!           3                                                &
!         )
!  End If

  LINCOMFlow%TerrainFile          = TerrainFile
  LINCOMFlow%RoughnessFile        = RoughnessFile
  LINCOMFlow%HGrid                = HGrid
  LINCOMFlow%ZGrid                = ZGrid
  LINCOMFlow%TerrainPert          = TerrainPert
  LINCOMFlow%RoughnessPert        = RoughnessPert
  LINCOMFlow%BackgroundWindAdded  = .true. ! $$ BackgroundWindAdded is hard wired as
                                           ! true as this is the only sensible option.
                                           ! Eventually could remove BackgroundWindAdded
                                           ! completely.
  LINCOMFlow%nXCompGrid           = nXCompGrid
  LINCOMFlow%nYCompGrid           = nYCompGrid
  LINCOMFlow%LincomExe            = LincomExe

  LINCOMFlow%SpaceAllocated       = .FALSE.
  LINCOMFlow%New                  = 0
  LINCOMFlow%Old                  = 0

  LINCOMFlow%C = InitCommonFlow(                   &
                   FlowModName    = 'LINCOM Flow', &
                   FlowName       = FlowName,      &
                   nMets          = 0,             &
                   MetModNames    = (/ ' ' /),     &
                   MetNames       = (/ ' ' /),     &
                   DomainName     = DomainName,    &
                   nAttribs       = 2,             &
                   AttribNames    = (/             &
                                      'Update',    &
                                      'Flow  '     &
                                    /),            &
                   FixedMet       = FixedMet,      &
                   UpdateOnDemand = UpdateOnDemand &
                 )

  Call AddUpdateSubsetsToCommonFlow(1, (/ UpdateSubsetName /), LINCOMFlow%C)

End Function InitLINCOMFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpCoordsEtc_LINCOMFlow(LINCOMFlow, Coords, Grids)
! Sets up Coords and Grids by adding any extra coords and grids which LINCOMFlow wants
! to define.

  Implicit None
  ! Argument list:
  Type(LINCOMFlow_), Intent(In)    :: LINCOMFlow  ! State of an LINCOM flow module instance.
  Type(Coords_),     Intent(InOut) :: Coords      ! Collection of coord systems.
  Type(Grids_),      Intent(InOut) :: Grids       ! Collection of grids.

  ! Locals:
  Type(HGrid_) :: HGrid
  Integer      :: iHGrid
  Real  ::  X0
  Real  ::  Y0

  iHGrid = FindHGridIndex(LINCOMFlow%HGrid,Grids)

  X0 = Grids%HGrids(iHGrid)%X0 +                                        &
      (Grids%HGrids(iHGrid)%nX - 1) * Grids%HGrids(iHGrid)%dX / 2.0

  Y0 = Grids%HGrids(iHGrid)%Y0 +                                        &
      (Grids%HGrids(iHGrid)%nY - 1) * Grids%HGrids(iHGrid)%dY / 2.0


  HGrid = InitHGrid(                                      &
            Name       = 'Flow for LINCOM',               &
            HCoordName = Grids%HGrids(iHGrid)%HCoordName, &
            Wrap       = .false.,                         &
            nX         = 1,                               &
            nY         = 1,                               &
            dX         = 1.0,                             &
            dY         = 1.0,                             &
            X0         = X0,                              &
            Y0         = Y0                               &
          )

  Call AddHGrid(HGrid, Grids)



End Subroutine SetUpCoordsEtc_LINCOMFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpLINCOMFlow_MetsCoordsEtc(EtaDefns, Coords, Grids, Mets, LINCOMFlow)
! Sets up LINCOMFlow using information from EtaDefns, Coords, Grids and Mets.

  Implicit None
  ! Argument list:
  Type(EtaDefns_), Intent(In)           :: EtaDefns   ! Collection of eta definitions.
  Type(Coords_),   Intent(In)           :: Coords     ! Collection of coord systems.
  Type(Grids_),    Intent(In)           :: Grids      ! Collection of grids.
  Type(Mets_),     Intent(In),   Target :: Mets       ! Set of met module instance
                                                      ! states.
  Type(LINCOMFlow_),  Intent(InOut)     :: LINCOMFlow ! State of an LINCOM flow module
                                                      ! instance.


  LINCOMFlow%iHGrid = FindHGridIndex(LINCOMFlow%HGrid,Grids)
  LINCOMFlow%iZGrid = FindZGridIndex(LINCOMFlow%ZGrid,Grids)
  LINCOMFlow%iHGridFlowForLINCOM = FindHGridIndex('Flow for LINCOM',Grids)

  LINCOMFlow%iHCoord = FindHCoordIndex(Grids%HGrids(LINCOMFlow%iHGrid)%HCoordName,Coords)
  LINCOMFlow%iZCoord = FindZCoordIndex(Grids%ZGrids(LINCOMFlow%iZGrid)%ZCoordName,Coords)


End Subroutine SetUpLINCOMFlow_MetsCoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForUpdateLINCOMFlow( &
             Coords, Grids, Domains,   &
             Mets,                     &
             iCase,                    &
             Time,                     &
             TValid, UpdateNow,        &
             LINCOMFlow,               &
             Units                     &
           )
! Prepares for updating an instance of the LINCOM flow module.

! This routine must set TValid and UpdateNow but must not alter the validity of the LINCOM flow module
! instance.
! $$ this comment should be with generic definition.

  Implicit None
  ! Argument list:
  Type(Coords_),     Intent(In)           :: Coords
  Type(Grids_),      Intent(In)           :: Grids
  Type(Domains_),    Intent(In),   Target :: Domains
  Type(Mets_),       Intent(In)           :: Mets
  Integer,           Intent(In)           :: iCase
  Type(Time_),       Intent(In)           :: Time
  Type(Time_),       Intent(Out)          :: TValid
  Logical,           Intent(Out)          :: UpdateNow
  Type(LINCOMFlow_), Intent(InOut)        :: LINCOMFlow
  Type(Units_),      Intent(InOut)        :: Units
  ! Coords     :: Collection of coord systems.
  ! Grids      :: Collection of grids.
  ! Domains    :: Collection of domains.
  ! Mets       :: Set of met module instance states.
  ! iCase      :: Number of case.
  ! Time       :: Time for which the LINCOM flow module instance might be updated.
  ! TValid     :: Earliest time that the validity (overall or for any single attribute) of the flow module
  !               instance might change, except perhaps for a change caused by a change in the validity of the
  !               met and flow module instances acting as data sources, assuming the flow module instance is
  !               updated now. The value is that determined at the end of this routine (the actual time may be
  !               later).
  ! UpdateNow  :: Indicates the flow module instance must be updated now (even if update-on-demand is
  !               specified). If set, TValid need not be set to any particular time.
  ! LINCOMFlow :: State of a LINCOM flow module instance.
  ! Units      :: Collection of information on input/output unit numbers.
  ! Locals:
  Type(Domain_), Pointer :: Domain ! Abbreviation for domain.

  ! Abbreviations.
  Domain => Domains%Domains(LINCOMFlow%C%iDomain)

  ! Validity due to time domain.
  If (Time < StartTimeOfDomain(Domain)) Then
    TValid = StartTimeOfDomain(Domain)
  Else If (Time >= EndTimeOfDomain(Domain)) Then
    TValid = InfFutureTime()
  Else
    TValid = TMin(EndTimeOfDomain(Domain), Time + Char2Time('01:00', Interval = .true.))
  End If
  UpdateNow = .false.

End Subroutine PrepareForUpdateLINCOMFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine LINCOMFlowReqs(                                              &
             LINCOMFlow, Time,                                          &
             WantFlowField,     WantCloudField,     WantRainField,      &
             iFlowUpdateSubset, iCloudUpdateSubset, iRainUpdateSubset,  &
             FlowField,         CloudField,         RainField           &
           )
! Specifies what information the flow module instance wants from the other flow module
! instances.

  Implicit None
  ! Argument list:
  Type(LINCOMFlow_),  Intent(In), Target  :: LINCOMFlow
  Type(Time_),        Intent(In)          :: Time
  Logical,            Intent(Out)         :: WantFlowField
  Logical,            Intent(Out)         :: WantCloudField
  Logical,            Intent(Out)         :: WantRainField
  Integer,            Intent(Out)         :: iFlowUpdateSubset
  Integer,            Intent(Out)         :: iCloudUpdateSubset
  Integer,            Intent(Out)         :: iRainUpdateSubset
  Type(FlowField_),   Intent(Out)         :: FlowField
  Type(CloudField_),  Intent(Out)         :: CloudField
  Type(RainField_),   Intent(Out)         :: RainField

  ! LINCOMFlow         :: State of an LINCOM flow module instance.
  ! Time               :: Current time in non-fixed met cases and the time of the met
  !                       in fixed met cases.
  ! WantFlowField      :} Indicates the flow module instance wants information for
  ! WantCloudField     :} various attributes from other flow module instances.
  ! WantRainField      :}
  ! iFlowUpdateSubset  :] Index of update subsets for various attributes in the array
  ! iCloudUpdateSubset :] of update subset information in the common part of the flow
  ! iRainUpdateSubset  :] module instance.
  ! FlowField          :} Field of information for various attributes used here to
  ! CloudField         :} specify what information the flow module instance wants from
  ! RainField          :} the other flow module instances.


  ! Flow requirements.
  WantFlowField           = .true.
  iFlowUpdateSubset       = 1
  FlowField%C%Valid       = .false. ! not needed, but convenient for consistency

  FlowField%C%iHGrid      = LINCOMFlow%iHGridFlowForLINCOM
  FlowField%C%iZGrid      = LINCOMFlow%iZGrid

  FlowField%C%Dt          = Char2Time('01:00', Interval = .true.)
  FlowField%C%UseTwoTimes = .true.
  FlowField%Flow          => LINCOMFlow%Flow
  FlowField%ProfileData   => LINCOMFlow%ProfileData


  ! Cloud requirements.
  WantCloudField = .false.
  ! Avoid compiler warning.
  iCloudUpdateSubset = 0
  CloudField%C%Valid = .false.

  ! Rain requirements.
  WantRainField = .false.
  ! Avoid compiler warning.
  iRainUpdateSubset = 0
  RainField%C%Valid = .false.

End Subroutine LINCOMFlowReqs

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateLINCOMFlow(                    &
             Coords, Grids, Domains,            &
             Mets,                              &
             iCase,                             &
             Time,                              &
             FlowField, CloudField, RainField,  &
             LINCOMFlow,                        &
             Units                              &
           )
! Updates an instance of the LINCOM flow module.

  Implicit None
  ! Argument list:
  Type(Coords_),     Intent(In)           :: Coords     ! Collection of coord systems.
  Type(Grids_),      Intent(In),   Target :: Grids      ! Collection of grids.
  Type(Domains_),    Intent(In),   Target :: Domains    ! Collection of domains.
  Type(Mets_),       Intent(In)           :: Mets       ! Set of met module instance states.
  Integer,           Intent(In)           :: iCase      ! Number of case.
  Type(Time_),       Intent(In)           :: Time       ! Time for which the LINCOM flow module
                                                        ! instance is to be updated.
  Type(FlowField_),  Intent(In)           :: FlowField  !} Field of information for various
  Type(CloudField_), Intent(In)           :: CloudField !} attributes needed by the flow
  Type(RainField_),  Intent(In)           :: RainField  !} module instance.
  Type(LINCOMFlow_), Intent(InOut)        :: LINCOMFlow ! State of an LINCOM flow module
                                                        ! instance.
  Type(Units_),      Intent(InOut)        :: Units      ! Collection of information on
                                                        ! input/output unit numbers.
  ! Locals:
  Type(HGrid_),  Pointer       :: HGrid
  Type(ZGrid_),  Pointer       :: ZGrid
  Type(Domain_), Pointer       :: Domain    ! Abbreviation for domain.
  Logical                      :: Valid     !} Validity due to time domain.
  Type(Time_)                  :: TValid    !}
  Logical                      :: MValid    !] Validity due to met.
  Type(Time_)                  :: MTValid   !]
  Integer                      :: IOStat    ! Error code for open statement.
  Integer                      :: Ix, Iy
  Character(MaxFileNameLength) :: FileName
  Integer                      :: Unit
  Logical                      :: DoIfBlock ! Indicates an if-block is to be executed
                                            ! where the logical test can't be evaluated
                                            ! in-line because of risk of error in testing
                                            ! equality of times.

  If (.not.FlowField%C%Valid) Then
    Call Message(                                                      &
           'ERROR updating flow module '                            // &
           Trim(LINCOMFlow%C%FlowName)                              // &
           ': no valid flow information available to drive Lincom',    &
           2                                                           &
         )
    LINCOMFlow%C%Valid                = .false.
    LINCOMFlow%C%ValidAttribParams(:) = .false.
    LINCOMFlow%C%ValidityUnimprovable = LINCOMFlow%C%Valid .or. LINCOMFlow%C%FixedMet
    LINCOMFlow%C%TValid               = InfFutureTime()
    Return
  End If

  HGrid  => Grids%HGrids(LINCOMFlow%iHGrid)
  ZGrid  => Grids%ZGrids(LINCOMFlow%iZGrid)

  Domain => Domains%Domains(LINCOMFlow%C%iDomain)

  If (.NOT. LINCOMFlow%SpaceAllocated) Then
    Allocate(LINCOMFlow%Velocity(Grids%HGrids(LINCOMFlow%iHGrid)%nX,      &
                      Grids%HGrids(LINCOMFlow%iHGrid)%nY,                 &
                      Grids%ZGrids(LINCOMFlow%iZGrid)%nZ,                 &
                      2,3))
    Allocate(LINCOMFlow%UStar(Grids%HGrids(LINCOMFlow%iHGrid)%nX,         &
                      Grids%HGrids(LINCOMFlow%iHGrid)%nY, 2))
    Allocate(LINCOMFlow%Z0(Grids%HGrids(LINCOMFlow%iHGrid)%nX,            &
                      Grids%HGrids(LINCOMFlow%iHGrid)%nY, 2))
    Allocate(LINCOMFlow%Topog(Grids%HGrids(LINCOMFlow%iHGrid)%nX,         &
                      Grids%HGrids(LINCOMFlow%iHGrid)%nY, 2))

    LINCOMFlow%SpaceAllocated = .TRUE.

! read in topography file
    FileName = Trim(LINCOMFlow%TerrainFile)

    Unit = OpenFile(                               &
             File            = FileName,           &
             Units           = Units,              &
             Status          = 'Old',              &
             Action          = 'Read',             &
             FileDescription = 'LINCOM topog file' &
           )

    Read(Unit,*)
    Read(Unit,*)
    Read(Unit,*)
    Read(Unit,*)
    Read(Unit,*)
    Read(Unit,*)

    Do Iy = 1, Grids%HGrids(LINCOMFlow%iHGrid)%nY
      Read(Unit,*)(LINCOMFlow%Topog(Ix,Iy,1),                      &
                    Ix = 1, Grids%HGrids(LINCOMFlow%iHGrid)%nX)
    EndDo

    Call CloseUnit(Unit, Units)

    LINCOMFlow%Topog(:,:,2) = LINCOMFlow%Topog(:,:,1)

  EndIf

  If (LINCOMFlow%New == 0) Then
    DoIfBlock = .false.
  Else
    DoIfBlock = Time2ShortTime(Time) == LINCOMFlow%NewTime
  End If

  If (DoIfBlock) Then

    ! Move Time1 met data into Old time level.
    LINCOMFlow%Old     = LINCOMFlow%New
    LINCOMFlow%OldTime = LINCOMFlow%NewTime

    ! Read Time2 met data into New time level.
    LINCOMFlow%New     = 3 - LINCOMFlow%Old
    LINCOMFlow%NewTime = Time2ShortTime(FlowField%C%Time2)

    Call RunLINCOM(                               &
               Mets,                              &
               Coords, Grids,                     &
               iCase,                             &
               Time,                              &
               FlowField, CloudField, RainField,  &
               LINCOMFlow,                        &
               Units, 2                           &
             )
  Else

    ! Mark met data in Old time level as invalid (otherwise ReadFieldsFile will assume
    ! Old time level data is for the previous time and may use it to fix-up missing
    ! data).
    LINCOMFlow%Old = 0

    ! Read Time1 met data into New time level.
    LINCOMFlow%New     = 1
    LINCOMFlow%NewTime = Time2ShortTime(Time)
    Call RunLINCOM(                               &
               Mets,                              &
               Coords, Grids,                     &
               iCase,                             &
               Time,                              &
               FlowField, CloudField, RainField,  &
               LINCOMFlow,                        &
               Units, 1                           &
             )

! write(6,*)LINCOMFlow%Old,LINCOMFlow%New
! write(6,*)LINCOMFlow%Velocity(60,60,2,LINCOMFlow%New,1)

    ! Move Time1 met data into Old time level.
    LINCOMFlow%Old     = LINCOMFlow%New
    LINCOMFlow%OldTime = LINCOMFlow%NewTime

    ! Read Time2 met data into New time level.
    LINCOMFlow%New     = 3 - LINCOMFlow%Old
    LINCOMFlow%NewTime = Time2ShortTime(FlowField%C%Time2)
    Call RunLINCOM(                               &
               Mets,                              &
               Coords, Grids,                     &
               iCase,                             &
               Time,                              &
               FlowField, CloudField, RainField,  &
               LINCOMFlow,                        &
               Units, 2                           &
             )
! write(6,*)LINCOMFlow%Old,LINCOMFlow%New
! write(6,*)LINCOMFlow%Velocity(60,60,2,LINCOMFlow%Old,1),LINCOMFlow%Velocity(60,60,2,LINCOMFlow%New,1)

  End If

  ! Validity due to time domain.
  Valid = StartTimeOfDomain(Domain) <= Time .and. Time < EndTimeOfDomain(Domain)
  If (Valid) Then
    TValid = TMin(FlowField%C%Time2, EndTimeOfDomain(Domain))
  Else
    TValid = StartTimeOfDomain(Domain)
    If (Time >= TValid) TValid = InfFutureTime()
  End If

  ! Overall validity.
  LINCOMFlow%C%Valid                     = Valid
  LINCOMFlow%C%ValidAttribParams(:)      = .false.
  LINCOMFlow%C%ValidAttribParams(A_Flow) = Valid
  LINCOMFlow%C%ValidityUnimprovable      = Valid .or. LINCOMFlow%C%FixedMet
  LINCOMFlow%C%TValid                    = TValid

End Subroutine UpdateLINCOMFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine RunLINCOM(                           &
             Mets,                              &
             Coords, Grids,                     &
             iCase,                             &
             Time,                              &
             FlowField, CloudField, RainField,  &
             LINCOMFlow,                        &
             Units, ITForFlow                   &
           )
! Updates an instance of the LINCOM flow module.

  Implicit None
  ! Argument list:
  Type(Mets_),       Intent(In)    :: Mets        ! Set of met module instance states.
  Type(Coords_),     Intent(In)    :: Coords      ! Collection of coord systems.

  Type(Grids_),      Intent(In), Target    :: Grids       ! Collection of grids.

  Integer,           Intent(In)    :: iCase       ! Number of case.
  Type(Time_),       Intent(In)    :: Time        ! Time for which the LINCOM flow module
                                                  ! instance is to be updated.
  Type(FlowField_),  Intent(In)    :: FlowField   !} Field of information for various
  Type(CloudField_), Intent(In)    :: CloudField  !} attributes needed by the flow
  Type(RainField_),  Intent(In)    :: RainField   !} module instance.
  Type(LINCOMFlow_), Intent(InOut) :: LINCOMFlow  ! State of an LINCOM flow module
                                                  ! instance.
  Type(Units_),      Intent(InOut) :: Units       ! Collection of information on
                                                  ! input/output unit numbers.
  Integer,           Intent(In)    :: ITForFlow
  ! Locals:
  Integer                :: IOStat   ! Error code for open statement.
  Integer                :: I, J, K
  Integer                :: Ix, Iy, Iz, In, ILevel
  Character(20)          :: FileName
  Character(1)           :: VelComp(3)
  Integer                :: Unit
  Real                   :: Speed
  Real                   :: Direction

  Real                   :: dZdX(3)
  Real                   :: Max_dZdX(2)
  Integer                :: Max_dZdX_X(2)
  Integer                :: Max_dZdX_Y(2)

  Type(HGrid_), Pointer  :: HGrid
  Type(ZGrid_), Pointer  :: ZGrid
  Type(HGrid_), Pointer  :: HGridForLincom

  HGrid          => Grids%HGrids(LINCOMFlow%iHGrid)
  ZGrid          => Grids%ZGrids(LINCOMFlow%iZGrid)
  HGridForLincom => Grids%HGrids(LINCOMFlow%iHGridFlowForLINCOM)


! generate LINCOM input file L1.INP

  Unit = OpenFile(                               &
           File            = 'l1f.inp',          &
           Units           = Units,              &
           Status          = 'Replace',          &
           FileDescription = 'LINCOM input file' &
         )

  Write(Unit,*)' "LINCOM Input File" '
  Write(Unit,*)'@'
  Write(Unit,*)'.FALSE.           FDebug        Flag for debbuging'

  If (LINCOMFlow%TerrainPert .EQV. .TRUE.) Then
    Write(Unit,*)'.TRUE.            FTerra        Flag for terrain pertubations'
  Else
    Write(Unit,*)'.FALSE.           FTerra        Flag for terrain pertubations'
  EndIf

  If (LINCOMFlow%RoughnessPert .EQV. .TRUE.) Then
    Write(Unit,*)'.TRUE.            FRough        Flag for roughness pertubations'
  Else
    Write(Unit,*)'.FALSE.           FRough        Flag for roughness pertubations'
  EndIf

  Write(Unit,*)'.FALSE.           FWtRgh       Flag for Water Roughness calculation'

  If (LINCOMFlow%BackgroundWindAdded .EQV. .TRUE.) Then
    Write(Unit,*)'.TRUE.            FBackg        Flag for adding background wind'
  Else
    Write(Unit,*)'.FALSE.           FBackg        Flag for adding background wind'
  EndIf

  Write(Unit,*)'.FALSE.           IDeriv       Flag for velocity derivative calculation'

  Write(Unit,*)'@'
  Write(Unit,*)'0               LinFld   0: Fields, 1: line, 2: both, -1: no output.'
  Write(Unit,*)'@'

  Write(Unit,*) Grids%HGrids(LINCOMFlow%iHGrid)%nX, Grids%HGrids(LINCOMFlow%iHGrid)%nY   &
                , '        MaxIT, MaxJT    Size of data domain'

  Write(Unit,*) Int(Grids%HGrids(LINCOMFlow%iHGrid)%X0),                                 &
                  Int(Grids%HGrids(LINCOMFlow%iHGrid)%Y0)                                &
                    ,'        XminT, YminT    x, y minimum for terrain data'

  Write(Unit,*) Grids%HGrids(LINCOMFlow%iHGrid)%dX, Grids%HGrids(LINCOMFlow%iHGrid)%dY   &
                , '    XInc, YInc      x, y grid point spacing'

  Write(Unit,*) '49.85         aLat     Latitude of area  [degree]'


  Write(Unit,*)'@'
  Write(Unit,*) "'",Trim(LINCOMFlow%TerrainFile),"'"                                     &
                   ,'     TerFil    Terrain height input file'
  Write(Unit,*)'1.0              TerFac    Terrain height multiplier'
  Write(Unit,*)'@'
  Write(Unit,*) "'",Trim(LINCOMFlow%RoughnessFile),"'"                                   &
                   ,'     RghFil    Roughness data input file'

  Write(Unit,*)'1.0              RghFac    Roughness multiplier'

  If (LINCOMFlow%RoughnessPert .EQV. .TRUE.) Then
    Write(Unit,*)'0.0              RghAdd    Roughness added (m)'
  Else
    Write(Unit,*)LINCOMFlow%ProfileData(1,1,1,1)%z0,'         RghAdd    Roughness added (m)'
  EndIf

  Write(Unit,*)'@                          (z0 if FRough=.False.)'
  Write(Unit,*)'@'
  Write(Unit,*) "'",'wval_ftc.grd',"'",'    WFtFil   Water fetch file'
  Write(Unit,*)'@'

!$$ this ensures that LINCOM uses the 10m wind for it's calcualtions
!$$ however this is not very robust or well thought out.

  ILevel = Grids%ZGrids(LINCOMFlow%iZGrid)%nZ
  Do I = 1, Grids%ZGrids(LINCOMFlow%iZGrid)%nZ
    If (Grids%ZGrids(LINCOMFlow%iZGrid)%Z(I) > 9.1) Then
      ILevel = I
      Exit
    End If
  End Do

  Speed = SQRT( LINCOMFlow%Flow(1,1,2,ITForFlow)%U(1)**2 + LINCOMFlow%Flow(1,1,2,ITForFlow)%U(2)**2 )
  Direction = Atan2ZeroTest(LINCOMFlow%Flow(1,1,2,ITForFlow)%U(2), LINCOMFlow%Flow(1,1,2,ITForFlow)%U(1))  &
              * 180.0 / Pi
  Direction = 270.0 - Direction

  If(Direction <   0.0) Direction = Direction + 360.0
  If(Direction > 360.0) Direction = Direction - 360.0

  Write(Unit,*)' 3               NType    Type of specified wind'
  Write(Unit,*)'@    Speed    Dir   Height      x          y'
  Write(Unit,*) Speed, Direction, Grids%ZGrids(LINCOMFlow%iZGrid)%Z(ILevel), HGridForLincom%X0, &
                HGridForLincom%Y0

  Write(Unit,*)'@'
  Write(Unit,*)'1                NPoint  Number of points of interest'
  Write(Unit,*)'@  Xsp         Ysp  X,Y-position of the point  [m]'
  Write(Unit,*) HGridForLincom%X0, HGridForLincom%Y0


  Write(Unit,*)'@'
  Write(Unit,*)Grids%ZGrids(LINCOMFlow%iZGrid)%nZ-1,'         Nlev      Number of levels'
! level heights
! level 1 is ground level and is not used by LINCOM.
  Do I = 2, Grids%ZGrids(LINCOMFlow%iZGrid)%nz
    Write(Unit,*)Grids%ZGrids(LINCOMFlow%iZGrid)%Z(i)
  EndDo

  Call CloseUnit(Unit, Units)

! system call to run LINCOM.

  Call SubmitSystemCommand(Trim(ConvertFileName(LINCOMFlow%LincomExe)) // ' > LincomOut.txt')

! read LINCOM output
! read in velocity fields

  VelComp = (/ 'U', 'V', 'W' /)

! loop through velocity component files

  Do J = 1,3                                       !  velocity component
    Do I = 2,Grids%ZGrids(LINCOMFlow%iZGrid)%nZ    !  number of LINCOM vertical levels

      FileName = VelComp(J)//Trim(Int2Char(I-1, Justify = 'L', FormatString = 'I2.2'))//'.GRD'

      Unit = OpenFile(                                &
               File            = FileName,            &
               Units           = Units,               &
               Status          = 'Old',               &
               Action          = 'Read',              &
               FileDescription = 'LINCOM output file' &
             )
      Read(Unit,*)
      Read(Unit,*)
      Read(Unit,*)
      Read(Unit,*)
      Read(Unit,*)
      Read(Unit,*)

      Iz = I
      In = J
      Do Iy = 1, Grids%HGrids(LINCOMFlow%iHGrid)%nY
        Read(Unit,*)(LINCOMFlow%Velocity(Ix,Iy,Iz,LINCOMFlow%New,In),            &
                      Ix = 1, Grids%HGrids(LINCOMFlow%iHGrid)%nX)
      EndDo
      Call CloseUnit(Unit, Units)

    EndDo
    LINCOMFlow%Velocity(                    &
      1:Grids%HGrids(LINCOMFlow%iHGrid)%nX, &
      1:Grids%HGrids(LINCOMFlow%iHGrid)%nY, &
      1, LINCOMFlow%New, In                 &
    ) = 0.0
  EndDo

! change W velocity to rate of change above ground

      Do i = 1, HGrid%nX
      Do j = 1, HGrid%nY
!
! Calculate gradients
        If (i == 1) Then
          dZdX(1) = - ( LINCOMFlow%Topog(i+1,j,1) - LINCOMFlow%Topog(i,j,1) ) / HGrid%dX
        Else If (i == HGrid%nX) Then
          dZdX(1) = - ( LINCOMFlow%Topog(i,j,1) - LINCOMFlow%Topog(i-1,j,1) ) / HGrid%dX
        Else
          dZdX(1) = - ( LINCOMFlow%Topog(i+1,j,1) - LINCOMFlow%Topog(i-1,j,1) ) / (2.0 * HGrid%dX)
        EndIf
!
        If (j == 1) Then
          dZdX(2) = - ( LINCOMFlow%Topog(i,j+1,1) - LINCOMFlow%Topog(i,j,1) ) / HGrid%dY
        Else If (j == HGrid%nY) Then
          dZdX(2) = - ( LINCOMFlow%Topog(i,j,1) - LINCOMFlow%Topog(i,j-1,1) ) / HGrid%dY
        Else
          dZdX(2) = - ( LINCOMFlow%Topog(i,j+1,1) - LINCOMFlow%Topog(i,j-1,1) ) / (2.0 * HGrid%dY)
        EndIf
!
        dZdX(3) = 1.0

        Do k = 1, ZGrid%nZ
! Calculate Z(AboveGround) Dot as d Z(AboveGround) / d (X(m), Y(m), Z(w)) *
! (U, V, Z(w) Dot).
!
          LINCOMFlow%Velocity(i,j,k,LINCOMFlow%New,3) =                     &
                 dZdX(1) * LINCOMFlow%Velocity(i,j,k,LINCOMFlow%New,1) +    &
                 dZdX(2) * LINCOMFlow%Velocity(i,j,k,LINCOMFlow%New,2) +    &
                 dZdX(3) * LINCOMFlow%Velocity(i,j,k,LINCOMFlow%New,3)
        End Do
!
! write out topog and W velocity information

        Max_dZdX(:) = 0.0
        
        If(dZdX(1) > Max_dZdX(1)) Then
          Max_dZdX(1) = dZdX(1)
          Max_dZdX_X(1) = i
          Max_dZdX_Y(1) = j
        EndIf
        If(dZdX(2) > Max_dZdX(2)) Then
          Max_dZdX(2) = dZdX(2)
          Max_dZdX_X(2) = i
          Max_dZdX_Y(2) = j
        EndIf

      End Do
      End Do

!! loop through velocity component files
!  Do J = 1,3                                       !  velocity component
!    Do I = 2,Grids%ZGrids(LINCOMFlow%iZGrid)%nZ    !  number of LINCOM vertical levels
!
!      FileName = VelComp(J)//Trim(Int2Char(I-1, Justify = 'L', &
!                    FormatString = 'I2.2'))//'.GRD1'
!      Unit = OpenFile(                              &
!               File   = FileName,                   &
!               Units  = Units,                      &
!               Status = 'Replace'                   &
!             )
!      Write(Unit,*)'100.0 '
!      Write(Unit,*)'100.0 '
!      Write(Unit,*)'100.0 '
!      Write(Unit,*)'100.0 '
!      Write(Unit,*)'100.0 '
!      Write(Unit,*)' '
!      Iz = I
!      In = J
!      Do Iy = 1, Grids%HGrids(LINCOMFlow%iHGrid)%nY
!        Write(Unit,*)(LINCOMFlow%Velocity(Ix,Iy,Iz,LINCOMFlow%New,In),            &
!                      Ix = 1, Grids%HGrids(LINCOMFlow%iHGrid)%nX)
!      Write(Unit,*)' '
!      EndDo
!      Call CloseUnit(Unit, Units)
!    EndDo
!  EndDo
!  stop


 !     Write(6,*)'MaxW = ',  MAXVAL(LINCOMFlow%Velocity(:,:,:,LINCOMFlow%New,3)),    &
 !               'MaxLoc = ',MAXLOC(LINCOMFlow%Velocity(:,:,:,LINCOMFlow%New,3))
 !     Write(6,*)'MinW = ',  MINVAL(LINCOMFlow%Velocity(:,:,:,LINCOMFlow%New,3)),    &
 !               'MinLoc = ',MINLOC(LINCOMFlow%Velocity(:,:,:,LINCOMFlow%New,3))
 !     Write(6,*)'Ave = ',SUM(ABS(LINCOMFlow%Velocity(:,:,:,LINCOMFlow%New,3)))/(HGrid%nX*HGrid%nY)
 !     Write(6,*)'Max_dZdX(1) = ',Max_dZdX(1),'Max_dZdX(2) = ',Max_dZdX(2)
 !     Write(6,*) Max_dZdX_X(1),Max_dZdX_Y(1),Max_dZdX_X(2),Max_dZdX_Y(2)


! Read in UStar and Z0 from LINCOM.
! Only done if roughness perturbation calculated.
!
! Ustar LINCOM output file
  If (LINCOMFlow%RoughnessPert .EQV. .TRUE.) Then
    FileName = 'USTARP.GRD'
    Unit = OpenFile(                                &
             File            = FileName,            &
             Units           = Units,               &
             Status          = 'Old',               &
             Action          = 'Read',              &
             FileDescription = 'LINCOM output file' &
           )
    Read(Unit,*)
    Read(Unit,*)
    Read(Unit,*)
    Read(Unit,*)
    Read(Unit,*)
    Read(Unit,*)

    Do Iy = 1, Grids%HGrids(LINCOMFlow%iHGrid)%nY
      Read(Unit,*)(LINCOMFlow%UStar(Ix,Iy,LINCOMFlow%New),                      &
                    Ix = 1, Grids%HGrids(LINCOMFlow%iHGrid)%nX)
    EndDo

    Call CloseUnit(Unit, Units)

! roughness LINCOM input file
    FileName = Trim(LINCOMFlow%RoughnessFile)
    Unit = OpenFile (                               &
             File            = FileName,            &
             Units           = Units,               &
             Status          = 'Old',               &
             Action          = 'Read',              &
             FileDescription = 'LINCOM output file' &
           )
    Read(Unit,*)
    Read(Unit,*)
    Read(Unit,*)
    Read(Unit,*)
    Read(Unit,*)
    Read(Unit,*)

    Do Iy = 1, Grids%HGrids(LINCOMFlow%iHGrid)%nY
      Read(Unit,*)(LINCOMFlow%Z0(Ix,Iy,LINCOMFlow%New),                         &
                   Ix = 1, Grids%HGrids(LINCOMFlow%iHGrid)%nX)
    EndDo

    Call CloseUnit(Unit, Units)

  EndIf


End Subroutine RunLINCOM

!-------------------------------------------------------------------------------------------------------------

Subroutine LINCOMFlowCoordIndices(LINCOMFlow, nHCoords, iHCoords, nZCoords, iZCoords)
! Returns indices of coord systems in which positions need to be specified when using
! GetLINCOMFlow.

  Implicit None
  ! Argument list:
  Type(LINCOMFlow_), Intent(In)  :: LINCOMFlow        ! State of an LINCOM flow module
                                                      ! instance.
  Integer,        Intent(Out) :: nHCoords             !} Number and values of indices
  Integer,        Intent(Out) :: iHCoords(MaxHCoords) !} of horizontal and vertical
  Integer,        Intent(Out) :: nZCoords             !} coord systems.
  Integer,        Intent(Out) :: iZCoords(MaxHCoords) !}


  ! Horizontal coords needed:
  ! 1) Coord system for horizontal grids used for met data.
  nHCoords    = 1
  iHCoords(1) = LINCOMFlow%iHCoord

  ! Vertical coords needed:
  ! 1) m agl.
  nZCoords    = 1
  iZCoords(1) = LINCOMFlow%iZCoord

End Subroutine LINCOMFlowCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine GetLINCOMFlow(              &
             Coords, Grids,            &
             LINCOMFlow,               &
             Moisture, Inhomog, Homog, &
             Time, Position,           &
             Flow, ProfileData         &
           )
! Gets flow information from an LINCOM flow module instance.

  Implicit None
  ! Argument list:
  Type(Coords_),        Intent(In)             :: Coords
  Type(Grids_),         Intent(In)             :: Grids
  Type(LINCOMFlow_),    Intent(In)             :: LINCOMFlow
  Logical,              Intent(In)             :: Moisture
  Logical,              Intent(In)             :: Inhomog
  Logical,              Intent(In)             :: Homog
  Type(ShortTime_),     Intent(In)             :: Time
  Type(Position_),      Intent(In)             :: Position
  Type(Flow_),          Intent(Out)            :: Flow
  Type(ProfileData_),   Intent(Out),  Optional :: ProfileData
  ! Coords      :: Collection of coord systems.
  ! Grids       :: Collection of grids.
  ! LINCOMFlow  :: State of an LINCOM flow module instance.
  ! Moisture    :: Indicates Q is required.
  ! Inhomog     :: Indicates inhomogeneous quantities are required.
  ! Homog       :: Indicates homogeneous quantities are required.
  ! Time        :: Current time in non-fixed met cases and the time of the met in
  !                fixed met cases.
  ! Position    :: Coords of the point in various coord systems in Coords, with
  !                flags to indicate whether the values are valid.
  ! Flow        :: Flow information.
  ! ProfileData :: Data needed to construct idealised analytic mean flow and
  !                turbulence profiles.
  ! Locals:
  Real(Std)          :: X            !} Horizontal location in coord system of
  Real(Std)          :: Y            !} horizontal grids used for met data.
  Real(Std)          :: Eta          ! Height in coord system of vertical grids used
                                     ! for met data.
  Real(Std)          :: Z            ! Height (m agl).
  Type(ProfileData_) :: ProfileDataL ! Local copy of ProfileData.
  Logical            :: FastRun      ! Indicates fast, but less precise, calculation
                                     ! required.
  Real(Std)          :: TauMin       ! Minimum timescale.
  Real(Std)          :: ZMin         ! Minimum height for turbulence quantities.
  Type(Flow_)        :: Flow2
  Type(TCoeffs_)     :: TCoeffs     !]

  FastRun = .false. ! $$ This should be input.
  TauMin  = 20.0    ! $$ This should be input.
  ZMin    =  0.0    ! $$ This should be input (and should be non-zero).

  Flow%iHCoord = LINCOMFlow%iHCoord
  Flow%iZCoord = LINCOMFlow%iZCoord

  ! $$ use ps systems near poles.

! Coords of location.
  X   = Position%XY(1, LINCOMFlow%iHCoord)
  Y   = Position%XY(2, LINCOMFlow%iHCoord)
  Z   = Position%Z (   LINCOMFlow%iZCoord)

! Set up other quantities required in ProfileDataL.

  ProfileDataL = LINCOMFlow%ProfileData(1,1,1,1)

! boundary layer height
!  Flow%H = ProfileDataL%H

  ! Calculate the following parts of Flow and ProfileDataL:
  ! Flow: U, dUdX, T, Theta, dThetadX, Q, P, Rho, dRhodX, Topog and dTopogdX.
  ! ProfileDataL: Z0, Phi0, UStar, RecipLMO, H and WStar.

  Call LincomFlowInfoFromMet(                          &
         Coords, Grids,                                &
         LINCOMFlow,                                   &
         FastRun,                                      &
         Time, X, Y, Z,                                &
         Flow, ProfileDataL, TCoeffs                   &
       )

  ! Set up other quantities required in Flow.
  Flow%dUdT(:)     = 0.0 ! ! $$ Also remember must be in chronological time direction - use IsBackwards
  Flow%ZS          = 0.0
  Flow%SmoothTopog = .true.

  ! 'Inhomogeneous turbulence quantities', 'inhomogeneous eddy diffusivities',
  ! 'homogeneous turbulence quantities', 'unresolved mesoscale motion quantities', 'boundary layer
  ! characteristics' and 'coord systems' parts of Flow.
  Call TurbProfiles(                 &
         Z           = Z,            &
         ProfileData = LINCOMFlow%ProfileData(1,1,1,1), &
         Inhomog     = Inhomog,      &
         Homog       = Homog,        &
         UEqV        = .true.,       & ! Use sig u \= sig v $$
         TauMin      = TauMin,       &
         ZMin        = ZMin,         &
         Flow        = Flow          &
       )
  Call TurbProfiles(                 &
         Z           = Z,            &
         ProfileData = LINCOMFlow%ProfileData(1,1,1,2), &
         Inhomog     = Inhomog,      &
         Homog       = Homog,        &
         UEqV        = .true.,       & ! Use sig u \= sig v $$
         TauMin      = TauMin,       &
         ZMin        = ZMin,         &
         Flow        = Flow2         &
       )

  ! Interpolate.
  ! $$ Only some quantities interpolated. Other quantities?
!  Flow%dUdT(:) = (Flow2%U(:) - Flow%U(:)) / 3600.0 ! $$
  Flow%SigUU(:)      = Flow%SigUU(:)    * TCoeffs%T2 + Flow2%SigUU(:)    * TCoeffs%T1
  Flow%dSigUUdX(:)   = Flow%dSigUUdX(:) * TCoeffs%T2 + Flow2%dSigUUdX(:) * TCoeffs%T1
  Flow%TauUU(:)      = Flow%TauUU(:)    * TCoeffs%T2 + Flow2%TauUU(:)    * TCoeffs%T1
  Flow%Eps           = Flow%Eps         * TCoeffs%T2 + Flow2%Eps         * TCoeffs%T1
  Flow%H             = Flow%H           * TCoeffs%T2 + Flow2%H           * TCoeffs%T1
  Flow%ZInterface(1) = Flow%H
  Flow%dTauUUdZ(:)   = Flow%dTauUUdZ(:) * TCoeffs%T2 + Flow2%dTauUUdZ(:) * TCoeffs%T1


  ! Time step limits.
  Flow%MaxDZ = Huge(Flow%MaxDZ)
  If (Inhomog) Then
    Flow%MaxDZ = Min(Flow%MaxDZ, Flow%H)
    If (Flow%dSigUUdX(3) /= 0.0) Then
      Flow%MaxDZ = Min(                                    &
                     Flow%MaxDZ,                           &
                     Abs(Flow%SigUU(3) / Flow%dSigUUdX(3)) &
                   )
    End If
    If (Flow%dTauUUdZ(3) /= 0.0) Then
      Flow%MaxDZ = Min(                                    &
                     Flow%MaxDZ,                           &
                     Abs(Flow%TauUU(3) / Flow%dTauUUdZ(3)) &
                   )
    End If
    If (Flow%dKdX(3) /= 0.0) Then
      Flow%MaxDZ = Min(                            &
                     Flow%MaxDZ,                   &
                     Flow%K(3) / Abs(Flow%dKdX(3)) &
                   )
    End If
  End If

  ! Puff size limits.
  Flow%DeltaI = Flow%H/4.0 ! $$

  ! ProfileData.
  If (Present(ProfileData)) ProfileData = ProfileDataL

End Subroutine GetLINCOMFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine LINCOMReflectCoordIndices(                           &
             LINCOMFlow, nHCoords, iHCoords, nZCoords, iZCoords &
           )
! Returns indices of coord systems in which positions need to be specified when using
! LINCOMReflect.

  Implicit None
  ! Argument list:
  Type(LINCOMFlow_), Intent(In)  :: LINCOMFlow           ! State of a LINCOM flow
                                                         ! module instance.
  Integer,           Intent(Out) :: nHCoords             !} Number and values of
  Integer,           Intent(Out) :: iHCoords(MaxHCoords) !} indices of horizontal and
  Integer,           Intent(Out) :: nZCoords             !} vertical coord systems.
  Integer,           Intent(Out) :: iZCoords(MaxHCoords) !}

  ! Horizontal coords needed: none
  nHCoords    = 0
  iHCoords(1) = 0 ! avoids compiler warning

  ! Vertical coords needed:
  ! 1) m agl.
  nZCoords    = 1
  iZCoords(1) = LINCOMFlow%iZCoord

End Subroutine LINCOMReflectCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine LINCOMReflect(                   &
             Coords, LINCOMFlow,            &
             Vel, OldPosition, Position, U, &
             Reflected                      &
           )
! Reflects position of particle or puff centroid using an LINCOM flow module instance.

  Implicit None
  ! Argument list:
  Type(Coords_),     Intent(In)    :: Coords
  Type(LINCOMFlow_), Intent(In)    :: LINCOMFlow
  Logical,           Intent(In)    :: Vel
  Type(Position_),   Intent(In)    :: OldPosition
  Type(Position_),   Intent(InOut) :: Position
  Real(Std),         Intent(InOut) :: U(3)
  Logical,           Intent(Out)   :: Reflected
  ! Coords      :: Collection of coord systems.
  ! LINCOMFlow  :: State of an LINCOM flow module instance.
  ! Vel         :: Indicates a dispersion model with velocity memory is being used.
  ! OldPosition :: Old particle position.
  ! Position    :: Current particle position.
  ! U           :: Current particle velocity. Defined only if Vel is true.
  ! Reflected   :: Indicates a reflection has occurred.
  ! Locals:
  Integer :: iZCoordZ ! Index of m agl coord system.

  iZCoordZ = LINCOMFlow%iZCoord

  Reflected = .false.

  ! Velocity memory.
  If (Vel) Then

    If (Position%Z(iZCoordZ) < 0.0) Then
      U(3)                               = - U(3)
      Position%Z(iZCoordZ)               = - Position%Z(iZCoordZ)
      Position%ZValid(1:Coords%nZCoords) = .false.
      Position%ZValid(iZCoordZ)          = .true.
      Position%RhoValid                  = .false.
      Reflected                          = .true.
    End If

  ! No velocity memory.
  Else

    If (Position%Z(iZCoordZ) < 0.0) Then
      Position%Z(iZCoordZ)               = - Position%Z(iZCoordZ)
      Position%ZValid(1:Coords%nZCoords) = .false.
      Position%ZValid(iZCoordZ)          = .true.
      Position%RhoValid                  = .false.
      Reflected                          = .true.
    End If

  End If

End Subroutine LINCOMReflect

!-------------------------------------------------------------------------------------------------------------

Subroutine LincomFlowInfoFromMet(                          &
             Coords, Grids,                                &
             M,                                            &
             FastRun,                                      &
             Time, X, Y, Z,                                &
             Flow, ProfileData, TCoeffs                    &
           )
! Extracts flow information at a particular location from the met fields.

  Implicit None
  ! Argument List:
  Type(Coords_),        Intent(In),    Target :: Coords
  Type(Grids_),         Intent(In),    Target :: Grids
  Type(LINCOMFlow_),    Intent(In)            :: M
  Logical,              Intent(In)            :: FastRun
  Type(ShortTime_),     Intent(In)            :: Time
  Real(Std),            Intent(In)            :: X
  Real(Std),            Intent(In)            :: Y
  Real(Std),            Intent(In)            :: Z
  Type(Flow_),          Intent(Out)           :: Flow
  Type(ProfileData_),   Intent(Out)           :: ProfileData
  Type(TCoeffs_),       Intent(Out)           :: TCoeffs     !]
 ! Coords         :: Collection of coord systems.
 ! Grids          :: Collection of grids.
 ! M              :: State of a LINCOM met module instance.
 ! Moisture       :: Indicates Q is required.
 ! FastRun        :: Indicates fast, but less precise, calculation required.
 ! Time           :: Current time in non-fixed met cases and the time of the met in
 !                   fixed met cases.
 ! X              :} Horizontal location.
 ! Y              :}
 ! Eta            :: Height in coord system of vertical grid used for met data.
 ! Flow           :: Flow information.
 ! ProfileData    :: Data needed to construct idealised analytic mean flow and
 !                   turbulence profiles.
 ! Locals:
  Real(Std)               :: dX           !} Grid separation in longitudinal and
  Real(Std)               :: dY           !} latitudinal direction in metres.
  Real(Std)               :: dZ           ! Height change between levels.
  Real(Std)               :: dZdX(2)      ! Horizontal height gradient (at constant
                                          ! Eta).
  Real(Std)               :: U0           ! X component of near surface wind.
  Real(Std)               :: V0           ! Y component of near surface wind.
  Type(HCoord_),  Pointer :: HCoord       !} Abbreviations for coords and grids.
  Type(ZCoord_),  Pointer :: ZCoord       !}
  Type(HGrid_),   Pointer :: HGrid        !}
  Type(ZGrid_),   Pointer :: ZGrid        !}
  Type(HCoeffs_)          ::  HCoeffs     !] Interpolation coefficients for various
  Type(HCoeffs_)          :: dHCoeffsdX   !] grids (including for obtaining
  Type(HCoeffs_)          :: dHCoeffsdY   !] derivatives of fields). The pointers are
  Type(ZCoeffs_)          ::  ZCoeffs     !] met datasets have time-averaged values of
  Type(ZCoeffs_)          :: dZCoeffsdZ   !] HeatFlux and UStar).


  ! Notes on derivatives: All differences are evaluated as 2 - 1, with the denominator
  ! dX, dY, dZ or dEta chosen consistently.

  ! 1) Set up abbreviations for coords, grids and interpolation coefficients.

  HCoord  => Coords%HCoords(M%iHCoord)
  ZCoord  => Coords%ZCoords(M%iZCoord)
  HGrid   => Grids%HGrids(M%iHGrid  )
  ZGrid   => Grids%ZGrids(M%iZGrid  )
  
  ProfileData%Canopy = .false.

  ! 2) Calculate interpolation coefficients.

! Calculate temporal interpolation coefficients.
    Call GetTCoeffs(                                 &
           Time,                                     &
           M%OldTime, Int(ShortTime2RealTime(M%NewTime-M%OldTime)), &
           M%Old, M%New,                             &
           TCoeffs                                   &
         )

  ! Calculate vertical interpolation coefficients.
  Call GetZCoeffs(Z, ZGrid, ZCoeffs)
  dZCoeffsdZ    =  ZCoeffs
  dZCoeffsdZ%Z1 =  1.0
  dZCoeffsdZ%Z2 = -1.0

  ! Calculate horizontal interpolation coefficients.
    Call GetHCoeffs(X, Y, HGrid, HCoeffs)
  dHCoeffsdX    =  HCoeffs
  dHCoeffsdY    =  HCoeffs
  dHCoeffsdX%X1 =  1.0
  dHCoeffsdX%X2 = -1.0
  dHCoeffsdY%Y1 =  1.0
  dHCoeffsdY%Y2 = -1.0


  ! 3) Extract data from arrays.

  ! 3-D data: U, V, W, T, Q, P.
  Call InterpXYZT(HCoeffs, ZCoeffs, TCoeffs, M%Velocity(:,:,:,:,1), Flow%U(1))
  Call InterpXYZT(HCoeffs, ZCoeffs, TCoeffs, M%Velocity(:,:,:,:,2), Flow%U(2))
  Call InterpXYZT(HCoeffs, ZCoeffs, TCoeffs, M%Velocity(:,:,:,:,3), Flow%U(3))

! support for gradients might need adding

! 2-D data: U0, V0, TS, PS, Z0, UStar, HeatFlux, H and Topog.
!  If (FastRun) Then
! Example lines for fast 'nearest' grid point approach.
!If (M%RoughnessPert .EQV. .TRUE.) Then
!    ProfileData%Z0    = M%Z0      (HCoeffs%iXNear, HCoeffs%iYNear, TCoeffs%iTNear )
!    ProfileData%UStar = M%UStar   (HCoeffs%iXNear, HCoeffs%iYNear, TCoeffs%iTNear)
!EndIf
!!

!  Else
!    Call InterpXYT( HCoeffs,    TCoeffs,  M%T(:, :, 2, :), TS               )
!    Call InterpXYT( HCoeffs,    TCoeffs,  M%P(:, :, 1, :), PS               )
!    Call InterpXYT( HCoeffs,    TCoeffs0, M%HeatFlux,      HeatFlux         )
!    Call InterpXYT( HCoeffs,    TCoeffs,  M%H,             ProfileData%H    )

    Call InterpXYT( HCoeffs,   TCoeffs,  M%Velocity(:, :, 2, :, 1), U0       )
    Call InterpXYT( HCoeffs,   TCoeffs,  M%Velocity(:, :, 2, :, 2), V0       )

    If (M%RoughnessPert .EQV. .TRUE.) Then
      Call InterpXYT( HCoeffs,    TCoeffs,    M%Z0,           ProfileData%Z0 )
      Call InterpXYT( HCoeffs,    TCoeffs, M%UStar,        ProfileData%UStar )
    EndIf

    Call InterpXYT( HCoeffs,     TCoeffs,  M%Topog,         Flow%Topog       )
    Call InterpXYT( dHCoeffsdX,  TCoeffs,  M%Topog,         Flow%dTopogdX(1) )
    Call InterpXYT( dHCoeffsdY,  TCoeffs,  M%Topog,         Flow%dTopogdX(2) )

!  End If


!  ! Surface wind direction, Monin-Obukhov length and WStar.

  ProfileData%Phi0 = ATan2ZeroTest(V0, U0)

!  RhoS = PS / (GasConstant * TS)
!  If (ProfileData%UStar == 0.0) Then
!    ProfileData%RecipLMO = - Sign(200000.0, HeatFlux)
!  Else
!    ProfileData%RecipLMO = - (VK * Gravity * HeatFlux) /             &
!                           (RhoS * Cp * TS * ProfileData%UStar ** 3)
!  End If
!  If (ProfileData%RecipLMO >= 0.0) Then
!    ProfileData%WStar = 0.0
!  Else
!    ProfileData%WStar = (ProfileData%H * Gravity * HeatFlux / (RhoS * Cp * TS)) &
!                        ** (1.0/3.0)
!  End If


!  Call MetricCoeffs(HCoord, (/ X, Y /), dX, dY)
!  dX = dX * HGrid%dX
!  dY = dY * HGrid%dY
!
!  Flow%dTopogdX(1) = Flow%dTopogdX(1) / dX
!  Flow%dTopogdX(2) = Flow%dTopogdX(2) / dY


End Subroutine LincomFlowInfoFromMet

!-------------------------------------------------------------------------------------------------------------

End Module LINCOMFlowModule



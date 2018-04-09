! Module: Single Site Met Module

Module SingleSiteMetModule

! This is a met module designed to read single site met data.

! It uses input from a file or files of single site met data.
!
! The coords stored in SingleSite%C are
! Horizontal: 1: Cartesian coord system used to define wind directions.
! Vertical:   1: m agl.
!
! No grids are stored in SingleSite%C.
!
! It works in a calendar or relative time frame. $$ SingleSite.F90 needs amending for
! 360-day year - solar elevation etc
!
! The status of the optional met module features is as follows (where ???? here =
! SingleSite):
!     SetUpCoordsEtc_????Met routine - used
!     SetUp????Met_CoordsEtc routine - used
!     Reset????Met routine           - used.     $$ could avoid by storing data ?
!
! Time labelling of met data:
! 1) Calendar time frames:
!    Each met data record must be labelled with the time.
! 2) Relative time frame:
!    The first record is at time zero and subsequent records are at fixed intervals.
!    Any time in the file is used in pre-processing (e.g. to calculate solar elevation
!    etc) but does not affect the time associated with the record.
! IgnoreFixedMetTime causes any value of FixedMetTime to be ignored - instead the
! first met record in the file is used (for the first case).
! For non-sequential data, FixedMet and IgnoreFixedMetTime must be true.
!
! Varying met data between cases:    $$
! Here there are two possibilities:
! 1) Multiple met files for multiple cases:
!    Here the cases step through the list of met files.
! 2) Multiple met records for multiple cases:
!    Here the cases step through the met records in the met file.
!    For this option, FixedMet and IgnoreFixedMetTime must be true.
!
! Varying met with FixedMet = true by incrementing FixedMetTime: $$
! IgnoreFixedMetTime must be false (hence the met data must be sequential and the
! multiple-met-records-for-multiple-cases option can't be used).
!
! $$ consider restartability - don't rely on knowing where one is in the met file. Use
! iCase and don't assume this starts at 1 and increments?

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowAndFlowProfileModule
Use SingleSiteModule
Use CommonMetModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: SingleSiteMet_                ! The state of a single site met module
                                         ! instance.
Public  :: InitSingleSiteMet             ! Initialises the state of a single site met
                                         ! module instance.
Public  :: SetUpCoordsEtc_SingleSiteMet  ! Sets up Coords and Grids by adding any extra
                                         ! coords and grids which SingleSiteMet wants
                                         ! to define.
Public  :: SetUpSingleSiteMet_CoordsEtc  ! Sets up SingleSiteMet using information from
                                         ! EtaDefns, Coords and Grids.
Public  :: PrepareForUpdateSingleSiteMet ! Prepares for updating an instance of the single site met module.
Public  :: UpdateSingleSiteMet           ! Updates an instance of the single site met
                                         ! module.
Public  :: ResetSingleSiteMet            ! Resets the state of a SingleSite met module
                                         ! instance for a new realisation.

!-------------------------------------------------------------------------------------------------------------

Type :: SingleSiteMet_ ! The state of a single site met module instance.

  ! CommonMet_:
  Type(CommonMet_) :: C ! The part of the met state common to all met modules.

  ! Input variables:
  Logical                      :: Sequential         ! Indicates that the met data
                                                     ! records are sequential.
  Type(Time_)                  :: Dt                 ! Time interval of met data (for
                                                     ! the case of sequential data).
  Logical                      :: IgnoreFixedMetTime ! Indicates that FixedMetTime is
                                                     ! to be ignored.
  Type(Time_)                  :: TSample            ! Sampling time for computing
                                                     ! unresolved mesoscale motions.
  Real(Std)                    :: Long               !} Longitude and latitude of the
  Real(Std)                    :: Lat                !} met site/dispersion area.
  Real(Std)                    :: Z                  !] Height of wind measurement,
  Real(Std)                    :: Z0                 !] roughness length, modified
  Real(Std)                    :: Z0D                !] Priestley-Taylor (moisture
  Real(Std)                    :: Alpha              !] availability) parameter,
  Real(Std)                    :: AlphaD             !] surface albedo, and maximum
  Real(Std)                    :: R                  !] value of 1/LMO (with values >
  Real(Std)                    :: RD                 !] 1E6 indicating no maximum).
  Real(Std)                    :: RecipLMOMax        !] Absence of suffix D indicates
  Real(Std)                    :: RecipLMOMaxD       !] value at met site; presence of
                                                     !] suffix D indicates value in
                                                     !] dispersion area. For all these
                                                     !] variables -999.0 indicates
                                                     !] value not specified.
  Real(Std)                    :: PrecipFactor       ! Precipitation correction factor
                                                     ! between met site and dispersion
                                                     ! area. -999.0 indicates value
                                                     ! not specified.
  Logical                      :: Representative     ! Indicates that no distinction
                                                     ! is to be made between met site
                                                     ! and dispersion area properties.
  Real(Std)                    :: SigUUM             ! Velocity variance of unresolved mesoscale motions.
  Real(Std)                    :: TauUUM             ! Lagrangian timescale of unresolved mesoscale motions.
  Real(Std)                    :: SigU2HPlus         !} Free tropospheric velocity variances (horizontal and
  Real(Std)                    :: SigW2HPlus         !} vertical).
  Real(Std)                    :: TauUHPlus          !] Free tropospheric Lagrangian timescales (horizontal
  Real(Std)                    :: TauWHPlus          !] and vertical).
  Character(MaxFileNameLength) :: MetFile            ! Met file.
  Character(MaxFileNameLength) :: MetLimitsFile      ! File containing limits within
                                                     ! which the met variables must
                                                     ! lie to be regarded as sensible.

  ! $$ Currently the priority order for Z, Z0, Alpha, R, RecipLMOMax is as follows:
  ! Use value given here if available. If not then use value from met file if
  ! available. If not then use default value (default available for Alpha and R only).
  ! Its not clear that this is the most sensible order.
  ! Note Z, Z0 must be given somewhere, but RecipLMOMax can be missing.
  ! The above agrees with ADMS.
  ! $$ add default for RecipLMOMax - won't change results but logic simpler.

  ! Internal quantities needed by the single site module:
  Type(Site_)         :: Site         !} Site characteristics at the met site and in
  Type(Site_)         :: SiteD        !} the dispersion area (including any
                                      !} information on the recent met which the
                                      !} single site module wishes to store).
  Type(SSMet_)        :: MinMet       !] Limits within which the met variables must
  Type(SSMet_)        :: MaxMet       !] lie to be regarded as sensible.
  Type(MetFileState_) :: MetFileState ! The state of the met file.

  ! Counters:
  Integer   :: iRecord    ! Index of the met record most recently read from the met
                          ! file.
  Integer   :: NInad      !} Cumulative number of met records and frequency of
  Integer   :: NCalm      !} occasions with inadequate data, with calm conditions, and
  Integer   :: NTotal     !} with neither inadequate data nor calm conditions.
  Real(Std) :: FInad      !}
  Real(Std) :: FCalm      !}
  Real(Std) :: FTotal     !} $$ how to define/treat for multi-case runs
  Logical   :: NoMoreData ! Indicates no more data available in the met file.

  ! Processed data:
  Type(MetFr_)       :: MetFr        ! Met frequency data.
                                     ! $$ need to consider how to treat this
  Type(ProfileData_) :: ProfileData1 !} Data needed to construct idealised analytic
  Type(ProfileData_) :: ProfileData2 !} mean flow and turbulence profiles using the
  Real(Std)          :: Cloud1       !} flow and flow profile module, together with
  Real(Std)          :: Cloud2       !} cloud amount (fraction), precipitation rate
  Real(Std)          :: Ppt1         !} (mm/hr), the data time, and validity
  Real(Std)          :: Ppt2         !} information. 1 and 2 refer to the two time
  Type(Time_)        :: Time1        !} levels.
  Type(Time_)        :: Time2        !}
  Logical            :: Valid1       !}
  Logical            :: Valid2       !}
  Real(Std)          :: ActualDt     ! Time2 - Time1 in seconds. This can be > Dt if data missing.
                                     ! Other reasons it might not = Dt? $$
                                     ! Consider whether should be a time structure. $$

End Type SingleSiteMet_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitSingleSiteMet(                  &
           MetName,                          &
           HCoordName,                       &
           FixedMet,                         &
           UpdateOnDemand,                   &
           Sequential,                       &
           IgnoreFixedMetTime,               &
           TSample,                          &
           Long, Lat,                        &
           Z, Z0, Z0D, Alpha, AlphaD,        &
           R, RD, RecipLMOMax, RecipLMOMaxD, &
           PrecipFactor,                     &
           Representative,                   &
           SigUUM, TauUUM,                   &
           SigU2HPlus, SigW2HPlus,           &
           TauUHPlus, TauWHPlus,             &
           MetFile,  MetLimitsFile           &
         )                                   &
Result(SingleSiteMet)
! Initialises an instance of the single site met module.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: MetName            ! Name of met module instance.
  Character(*), Intent(In) :: HCoordName         ! $$
  Logical,      Intent(In) :: FixedMet           ! Indicates that the met is fixed.
  Logical,      Intent(In) :: UpdateOnDemand     ! Indicates the met module instance is to be updated using
                                                 ! update-on-demand.
  Logical,      Intent(In) :: Sequential         ! Indicates that the met data records
                                                 ! are sequential.
  Logical,      Intent(In) :: IgnoreFixedMetTime ! Indicates that FixedMetTime is to
                                                 ! be ignored.
  Type(Time_),  Intent(In) :: TSample            ! Sampling time for computing
                                                 ! unresolved mesoscale motions.
  Real(Std),    Intent(In) :: Long               !} Longitude and latitude of the met
  Real(Std),    Intent(In) :: Lat                !} site/dispersion area.
  Real(Std),    Intent(In) :: Z                  !] Height of wind measurement,
  Real(Std),    Intent(In) :: Z0                 !] roughness length, modified
  Real(Std),    Intent(In) :: Z0D                !] Priestley-Taylor (moisture
  Real(Std),    Intent(In) :: Alpha              !] availability) parameter, surface
  Real(Std),    Intent(In) :: AlphaD             !] albedo, and maximum value of 1/LMO
  Real(Std),    Intent(In) :: R                  !] (with values > 1E6 indicating no
  Real(Std),    Intent(In) :: RD                 !] maximum). Absence of suffix D
  Real(Std),    Intent(In) :: RecipLMOMax        !] indicates value at met site;
  Real(Std),    Intent(In) :: RecipLMOMaxD       !] presence of suffix D indicates
                                                 !] value in dispersion area. For all
                                                 !] these variables -999.0 indicates
                                                 !] value not specified.
  Real(Std),    Intent(In) :: PrecipFactor       ! Precipitation correction factor
                                                 ! between met site and dispersion
                                                 ! area. -999.0 indicates value not
                                                 ! specified.
  Logical,      Intent(In) :: Representative     ! Indicates that no distinction is to
                                                 ! be made between met site and
                                                 ! dispersion area properties.
  Real(Std),    Intent(In) :: SigUUM             ! Velocity variance of unresolved mesoscale motions.
  Real(Std),    Intent(In) :: TauUUM             ! Lagrangian timescale of unresolved mesoscale motions.
  Real(Std),    Intent(In) :: SigU2HPlus         !} Free tropospheric velocity variances (horizontal and
  Real(Std),    Intent(In) :: SigW2HPlus         !} vertical).
  Real(Std),    Intent(In) :: TauUHPlus          !] Free tropospheric Lagrangian timescales (horizontal and
  Real(Std),    Intent(In) :: TauWHPlus          !] vertical).
  Character(*), Intent(In) :: MetFile            ! Met file.
  Character(*), Intent(In) :: MetLimitsFile      ! File containing limits within which
                                                 ! the met variables must lie to be
                                                 ! regarded as sensible.
  ! Function result:
  Type(SingleSiteMet_) :: SingleSiteMet ! Initialised state of a single site met
                                        ! module instance.
  ! Locals:
!  Logical :: Error ! Indicates an error in calling InitMetInput. $$ del

  ! Note checks on values of TSample, Long, Lat, Z, Z0, Z0D, Alpha, AlphaD, R, RD,
  ! RecipLMOMax, RecipLMOMaxD, PrecipFactor (including whether they are missing or
  ! not) are handled by the single site module.
  ! $$ not necessarily comprehensively (e.g. TSample which in particular might be
  ! better handled here)

  ! Sequential, FixedMet and IgnoreFixedMetTime.
  If (.not. Sequential) Then
    If (.not. FixedMet .or. .not. IgnoreFixedMetTime) Then
      Call Message(                                              &
             'Error in InitSingleSiteMet: For non-seq data, ' // &
             'FixedMet and IgnoreFixedMetTime must be true',     &
             3                                                   &
           )
    End If
  End If

  ! MetFile and MetLimitsFile.
  If (Len_Trim(MetFile) > MaxFileNameLength) Then
    Call Message(                                           &
           'Error in InitSingleSiteMet: met file name (' // &
           Trim(MetFile)                                 // &
           ') is too long',                                 &
           3                                                &
         )
  End If
  If (Len_Trim(MetFile) <= 0) Then
    Call Message('Error in InitSingleSiteMet: met file name is blank', 3)
  End If
  If (Len_Trim(MetLimitsFile) > MaxFileNameLength) Then
    Call Message(                                                  &
           'Error in InitSingleSiteMet: met limits file name (' // &
           Trim(MetFile)                                        // &
           ') is too long',                                        &
           3                                                       &
         )
  End If
  If (Len_Trim(MetFile) <= 0) Then
    Call Message('Error in InitSingleSiteMet: met limits file name is blank', 3)
  End If
  ! Single Site Module expects string of length 32.
  If (Len_Trim(MetLimitsFile) > 32) Then
    Call Message(                                                  &
           'Error in InitSingleSiteMet: met limits file name (' // &
           Trim(MetFile)                                        // &
           ') is too long',                                        &
           3                                                       &
         )
  End If

  SingleSiteMet%C = InitCommonMet('Single Site Met', MetName, FixedMet, UpdateOnDemand)

  Call AddCoordsToCommonMet(               &
         1, (/ HCoordName /), 0, (/ ' '/), &
         SingleSiteMet%C                   &
       )

  SingleSiteMet%Sequential         = Sequential
  SingleSiteMet%Dt                 = Char2Time('01:00', Interval = .true.) ! $$ 1 hour hard wired.
  SingleSiteMet%IgnoreFixedMetTime = IgnoreFixedMetTime
  SingleSiteMet%TSample            = TSample
  SingleSiteMet%Long               = Long
  SingleSiteMet%Lat                = Lat
  SingleSiteMet%Z                  = Z
  SingleSiteMet%Z0                 = Z0
  SingleSiteMet%Z0D                = Z0D
  SingleSiteMet%Alpha              = Alpha
  SingleSiteMet%AlphaD             = AlphaD
  SingleSiteMet%R                  = R
  SingleSiteMet%RD                 = RD
  SingleSiteMet%RecipLMOMax        = RecipLMOMax
  SingleSiteMet%RecipLMOMaxD       = RecipLMOMaxD
  SingleSiteMet%PrecipFactor       = PrecipFactor
  SingleSiteMet%Representative     = Representative
  SingleSiteMet%SigUUM             = SigUUM
  SingleSiteMet%TauUUM             = TauUUM
  SingleSiteMet%SigU2HPlus         = SigU2HPlus
  SingleSiteMet%SigW2HPlus         = SigW2HPlus
  SingleSiteMet%TauUHPlus          = TauUHPlus
  SingleSiteMet%TauWHPlus          = TauWHPlus
  SingleSiteMet%MetFile            = ConvertFileName(MetFile)       ! $$ Conversion redundant once
                                                                    ! $$ singlesite.f90
  SingleSiteMet%MetLimitsFile      = ConvertFileName(MetLimitsFile) ! $$ uses openfile routine ??

 ! $$ del
 ! Call InitMetInput(                                                     &
 !        Long, Lat,                                                      &
 !        ShortTime2RealTime(Time2ShortTime(TSample))/3600.0,             &
 !        Z, Z0, Z0D, Alpha, AlphaD, R, RD, RecipLMOMax, RecipLMOMaxD,    &
 !        PrecipFactor,                                                   &
 !        SingleSiteMet%MetFile, SingleSiteMet%MetLimitsFile(1:32),       & ! $$
 !        Sequential, Representative,                                     &
 !        SingleSiteMet%MetFr,                                            &
 !        SingleSiteMet%iRecord,                                          &
 !        SingleSiteMet%NInad, SingleSiteMet%NCalm, SingleSiteMet%NTotal, &
 !        SingleSiteMet%FInad, SingleSiteMet%FCalm, SingleSiteMet%FTotal, &
 !        Error,                                                          &
 !        SingleSiteMet%MinMet, SingleSiteMet%MaxMet,                     &
 !        SingleSiteMet%Site,   SingleSiteMet%SiteD,                      &
 !        SingleSiteMet%MetFileState                                      &
 !      )

 ! If (Error) Then
 !   Call Message('Error in InitSingleSiteMet: Unable to initialise single site module', 3)
 ! End If

  ! NoMoreData is false until we know there is no more data.
  SingleSiteMet%NoMoreData = .false.

  SingleSiteMet%Valid1 = .false.
  SingleSiteMet%Valid2 = .false.

End Function InitSingleSiteMet

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpCoordsEtc_SingleSiteMet(SingleSiteMet, Coords, Grids)
! Sets up Coords and Grids by adding any extra coords and grids which SingleSiteMet
! wants to define.

  Implicit None
  ! Argument list:
  Type(SingleSiteMet_), Intent(In)    :: SingleSiteMet ! State of a single site met
                                                       ! module instance.
  Type(Coords_),        Intent(InOut) :: Coords        ! Collection of coord systems.
  Type(Grids_),         Intent(InOut) :: Grids         ! Collection of grids.
  ! Locals:
  Type(ZCoord_) :: ZCoord ! Coord system to be added.

  ! Construct coord system to be added to Coords.
  ZCoord = ZCoord_m_agl()

  ! Add coord system to Coords.
  Call AddZCoord(ZCoord, Coords)

End Subroutine SetUpCoordsEtc_SingleSiteMet

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpSingleSiteMet_CoordsEtc(EtaDefns, Coords, Grids, MetEnsembleSize, SingleSiteMet)
! Sets up SingleSiteMet using information from EtaDefns, Coords and Grids.
! $$ not quite true now MetEnsembleSize added.

! $$ Add tests as appropriate on MetEnsembleSize.

  Implicit None
  ! Argument list:
  Type(EtaDefns_),      Intent(In)    :: EtaDefns        ! Collection of eta definitions. !$$ not used
  Type(Coords_),        Intent(In)    :: Coords          ! Collection of coord systems.
  Type(Grids_),         Intent(In)    :: Grids           ! Collection of grids.
  Integer,              Intent(In)    :: MetEnsembleSize ! Size of the met ensemble (i.e. number of met
                                                         ! realisations).
  Type(SingleSiteMet_), Intent(InOut) :: SingleSiteMet   ! SingleSite met module instance states.

  ! Add coord systems to SingleSiteMet%C.
  Call AddCoordsToCommonMet( &
         0, (/ ' ' /),       &
         1, (/ 'm agl' /),   &
         SingleSiteMet%C     &
       )

  ! Check number of coord systems in SingleSiteMet.
  If (SingleSiteMet%C%nHCoords /= 1) Then
    Call Message(                                                  &
           'Unexpected error in SetUpSingleSiteMet_CoordsEtc: ' // &
           'error in numbering of horizontal coord systems '    // &
           'used by the single site met module instance "'      // &
           Trim(SingleSiteMet%C%MetName)                        // &
           '"',                                                    &
           4                                                       &
         )
  End If
  If (SingleSiteMet%C%nZCoords /= 1) Then
    Call Message(                                                  &
           'Unexpected error in SetUpSingleSiteMet_CoordsEtc: ' // &
           'error in numbering of vertical coord systems '      // &
           'used by the single site met module instance "'      // &
           Trim(SingleSiteMet%C%MetName)                        // &
           '"',                                                    &
           4                                                       &
         )
  End If

End Subroutine SetUpSingleSiteMet_CoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForUpdateSingleSiteMet( &
             Coords, Grids,               &
             iCase, iMetCase,             &
             Time,                        &
             TValid, UpdateNow,           &
             SingleSiteMet,               &
             Units                        &
           )
! Prepares for updating an instance of the single site met module.

! This routine must set TValid and UpdateNow but must not alter the validity of the single site met module
! instance.
! $$ this comment should be with generic definition.

  Implicit None
  ! Argument List:
  Type(Coords_),        Intent(In)           :: Coords
  Type(Grids_),         Intent(In),   Target :: Grids
  Integer,              Intent(In)           :: iCase
  Integer,              Intent(In)           :: iMetCase
  Type(Time_),          Intent(In)           :: Time
  Type(Time_),          Intent(Out)          :: TValid
  Logical,              Intent(Out)          :: UpdateNow
  Type(SingleSiteMet_), Intent(InOut)        :: SingleSiteMet
  Type(Units_),         Intent(InOut)        :: Units
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! iCase         :: Number of case.
  ! iMetCase      :: Number of the met realisation in the met ensemble.
  ! Time          :: Time for which the single site met module instance might be updated.
  ! TValid        :: Earliest time that the validity (overall or for any single attribute) of the met module
  !                  instance might change, assuming the met module instance is updated now. The value is that
  !                  determined at the end of this routine (the actual time may be later).
  ! UpdateNow     :: Indicates the met module instance must be updated now (even if update-on-demand is
  !                  specified). If set, TValid need not be set to any particular time.
  ! SingleSiteMet :: State of a single site met module instance.
  ! Units         :: Collection of information on input/output unit numbers.
  ! Locals:
  Type(Time_) :: MetTime1 !} Times of the met data required.
  Type(Time_) :: MetTime2 !}

  TValid    = ReferenceTime()
  UpdateNow = .true.

End Subroutine PrepareForUpdateSingleSiteMet

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateSingleSiteMet( &
             Coords, Grids,     &
             iCase, iMetCase,   &
             Time,              &
             SingleSiteMet,     &
             Units              &
           )
! Updates an instance of the single site met module.

  Implicit None
  ! Argument list:
  Type(Coords_),        Intent(In)    :: Coords        ! Collection of coord systems.
  Type(Grids_),         Intent(In)    :: Grids         ! Collection of grids.
  Integer,              Intent(In)    :: iCase         ! Number of case.
  Integer,              Intent(In)    :: iMetCase      ! Number of the met realisation in the met ensemble.
  Type(Time_),          Intent(In)    :: Time          ! Time for which the single site met module instance is
                                                       ! to be updated.
  Type(SingleSiteMet_), Intent(InOut) :: SingleSiteMet ! State of a single site met module instance.
  Type(Units_),         Intent(InOut) :: Units         ! Collection of information on input/output unit
                                                       ! numbers.

  ! If FixedMet and IgnoreFixedMetTime are true, use next record in met file and use
  ! it for times 1 and 2.
  If (SingleSiteMet%C%FixedMet .and. SingleSiteMet%IgnoreFixedMetTime) Then

    Call ReadSS(SingleSiteMet)
    SingleSiteMet%ProfileData1 = SingleSiteMet%ProfileData2
    SingleSiteMet%Cloud1       = SingleSiteMet%Cloud2
    SingleSiteMet%Ppt1         = SingleSiteMet%Ppt2
    SingleSiteMet%Time1        = SingleSiteMet%Time2
    SingleSiteMet%Valid1       = SingleSiteMet%Valid2

  ! Otherwise, read met file until find valid met records either side of the time Time
  ! or after the time Time (the latter will only occur if the met file doesn't have
  ! any valid met records before the time Time) or until no more data is available.
  ! Note that:
  ! (i)   before first updating, Valid1, Valid2 and NoMoreData are false,
  ! (ii)  before subsequent updatings and before the end of the met file has been
  !       reached, Valid1 and Valid2 are true and NoMoreData is false,
  ! (iii) after the end of the met file has been reached Valid2 is false and
  !       NoMoreData is true.
  ! $$ Note need two records even for fixed met.
  Else

    ! If first update, get two valid records (unless no more data available).
    If (                                  &
      .not.SingleSiteMet%Valid1     .and. &
      .not.SingleSiteMet%Valid2     .and. &
      .not.SingleSiteMet%NoMoreData       &
    ) Then

      Do While (.not.SingleSiteMet%Valid2)
        Call ReadSS(SingleSiteMet)
        If (SingleSiteMet%NoMoreData) Exit
      End Do
      If (.not.SingleSiteMet%NoMoreData) Then
        SingleSiteMet%ProfileData1 = SingleSiteMet%ProfileData2
        SingleSiteMet%Cloud1       = SingleSiteMet%Cloud2
        SingleSiteMet%Ppt1         = SingleSiteMet%Ppt2
        SingleSiteMet%Time1        = SingleSiteMet%Time2
        SingleSiteMet%Valid1       = SingleSiteMet%Valid2
        SingleSiteMet%Valid2       = .false.
        Do While (.not.SingleSiteMet%Valid2)
          Call ReadSS(SingleSiteMet)
          If (SingleSiteMet%NoMoreData) Exit
        End Do
      End If

    End If

    ! Get two valid records either side of the time Time or after the time Time
    ! (unless no more data available).
    If (                                  &
      SingleSiteMet%Valid1          .and. &
      SingleSiteMet%Valid2          .and. &
      .not.SingleSiteMet%NoMoreData       &
    ) Then

      ! Note times in met file assumed to be in order - add check?
      OuterLoop: Do While (SingleSiteMet%Time2 <= Time)
        SingleSiteMet%ProfileData1 = SingleSiteMet%ProfileData2
        SingleSiteMet%Cloud1       = SingleSiteMet%Cloud2
        SingleSiteMet%Ppt1         = SingleSiteMet%Ppt2
        SingleSiteMet%Time1        = SingleSiteMet%Time2
        SingleSiteMet%Valid1       = SingleSiteMet%Valid2
        SingleSiteMet%Valid2       = .false.
        Do While (.not.SingleSiteMet%Valid2)
          Call ReadSS(SingleSiteMet)
          If (SingleSiteMet%NoMoreData) Exit OuterLoop
        End Do
      End Do OuterLoop

    End If

    If (.not.SingleSiteMet%NoMoreData) Then
      If (SingleSiteMet%Time1 > Time) Then
        Call Message(                                                     &
               'ERROR: First valid record in Single Site met dataset ' // &
               Trim(SingleSiteMet%MetFile)                             // &
               ' is too late',                                            &
               2                                                          &
             )
      End If
    End If

  End If

  SingleSiteMet%ActualDt = ShortTime2RealTime(                                                         &
                             Time2ShortTime(SingleSiteMet%Time2) - Time2ShortTime(SingleSiteMet%Time1) &
                           )

  ! $$ Max interval between time1 and time 2 hardwired at 03:00

  ! $$ This routine will be v complex for backwards dispersion. Better
  ! to set up array of ProfileData's etc? Also converge code structure with NWPMet.

  ! Update validity information.
  If (SingleSiteMet%Valid1 .and. SingleSiteMet%Valid2) Then
    If (                                                                  &
      (                                                                   &
        Char2Time('03:00', Interval = .true.)                             &
        >=                                                                &
        SingleSiteMet%Time2 - SingleSiteMet%Time1                         &
      )                                                                   &
      .and.                                                               &
      (                                                                   &
        SingleSiteMet%Time1 <= Time                                       &
        .or.                                                              &
        (SingleSiteMet%C%FixedMet .and. SingleSiteMet%IgnoreFixedMetTime) &
      )                                                                   &
    ) Then
      SingleSiteMet%C%Valid = .true.
    Else
      SingleSiteMet%C%Valid = .false.
    End If
  Else
    SingleSiteMet%C%Valid = .false.
  End If

  If (SingleSiteMet%C%FixedMet .or. SingleSiteMet%NoMoreData) Then
    SingleSiteMet%C%TValid = InfFutureTime()
  Else If (SingleSiteMet%Time1 <= Time) Then
    SingleSiteMet%C%TValid = SingleSiteMet%Time2
  Else
    SingleSiteMet%C%TValid = SingleSiteMet%Time1
  End If

End Subroutine UpdateSingleSiteMet

!-------------------------------------------------------------------------------------------------------------

Subroutine ResetSingleSiteMet(SingleSiteMet)
! Resets the state of a SingleSite met module instance for a new realisation.

  Implicit None
  ! Argument list:
  Type(SingleSiteMet_), Intent(InOut) :: SingleSiteMet ! State of a single site met
                                                       ! module instance.
  ! Locals:
  Logical :: Error ! Indicates an error in calling InitMetInput.

  ! may need to step through to right place or avoid resetting if multi-cases
  ! are to step though the met records $$

  ! note this is called before case 1, so will call InitMetInput before case 1.
  ! Logically, should rename 'reset' as 'prepare' $$

  Call InitMetInput(                                                       &
         SingleSiteMet%Long, SingleSiteMet%Lat,                            &
         ShortTime2RealTime(Time2ShortTime(SingleSiteMet%TSample))/3600.0, &
         SingleSiteMet%Z,                                                  &
         SingleSiteMet%Z0,          SingleSiteMet%Z0D,                     &
         SingleSiteMet%Alpha,       SingleSiteMet%AlphaD,                  &
         SingleSiteMet%R,           SingleSiteMet%RD,                      &
         SingleSiteMet%RecipLMOMax, SingleSiteMet%RecipLMOMaxD,            &
         SingleSiteMet%PrecipFactor,                                       &
         SingleSiteMet%MetFile, SingleSiteMet%MetLimitsFile(1:32),         & ! $$
         SingleSiteMet%Sequential, SingleSiteMet%Representative,           &
         SingleSiteMet%MetFr,                                              &
         SingleSiteMet%iRecord,                                            &
         SingleSiteMet%NInad, SingleSiteMet%NCalm, SingleSiteMet%NTotal,   &
         SingleSiteMet%FInad, SingleSiteMet%FCalm, SingleSiteMet%FTotal,   &
         Error,                                                            &
         SingleSiteMet%MinMet, SingleSiteMet%MaxMet,                       &
         SingleSiteMet%Site,   SingleSiteMet%SiteD,                        &
         SingleSiteMet%MetFileState                                        &
       )

  If (Error) Then
    Call Message('FATAL ERROR: Unable to prepare single site met module', 3)
  End If

  SingleSiteMet%Valid1  = .false.
  SingleSiteMet%Valid2  = .false.

End Subroutine ResetSingleSiteMet

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadSS(SingleSiteMet)
! Reads a record of met data and uses it to update the Profile2, Cloud2, Ppt2, Valid2
! and Time2 parts of SingleSiteMet.

  Implicit None
  ! Argument list:
  Type(SingleSiteMet_), Intent(InOut) :: SingleSiteMet ! State of a single site met
                                                       ! module instance.
  ! Locals:
  Type(SSMet_) :: MetAsRead      ! A record of met data as read from the met file.
  Type(SSMet_) :: ProcMet        ! A record of met data as processed by the single
                                 ! site module.
  Logical      :: NoMoreData     ! Indicates no more data in the met file.
  Logical      :: Error          ! Indicates an error in calling UpdateMetInput.
  Logical      :: InadequateData ! Indicates inadequate data in the record of met
                                 ! data.
  Logical      :: Calm           ! Indicates calm conditions in the record of met
                                 ! data.
  Real(Std)    :: Zone           !} Time zone information for the met data.
  Integer      :: ZoneHour       !}
  Integer      :: ZoneMinute     !}

  ! Call UpdateMetInput.
  Call UpdateMetInput(                                                     &
         ShortTime2RealTime(Time2ShortTime(SingleSiteMet%TSample))/3600.0, &
         SingleSiteMet%PrecipFactor,                                       &
         SingleSiteMet%Sequential, SingleSiteMet%Representative,           &
         MetAsRead, ProcMet,                                               &
         SingleSiteMet%MetFr,                                              &
         SingleSiteMet%iRecord,                                            &
         SingleSiteMet%NInad, SingleSiteMet%NCalm, SingleSiteMet%NTotal,   &
         SingleSiteMet%FInad, SingleSiteMet%FCalm, SingleSiteMet%FTotal,   &
         Error, NoMoreData, InadequateData, Calm,                          &
         SingleSiteMet%MinMet, SingleSiteMet%MaxMet,                       &
         SingleSiteMet%Site,   SingleSiteMet%SiteD,                        &
         SingleSiteMet%MetFileState                                        &
       )

  ! SingleSiteMet%NoMoreData. Note this indicates no more data available in the met
  ! file - this might be due to an error rather than the end of the met file being
  ! reached.
  SingleSiteMet%NoMoreData = NoMoreData .or. Error

  ! SingleSiteMet%Valid2.
  SingleSiteMet%Valid2 = .not.(NoMoreData .or. Error .or. InadequateData)

  If (NoMoreData) Then
    Call Message(                                               &
           'ERROR: No more data in Single Site met dataset ' // &
           Trim(SingleSiteMet%MetFile),                         &
           2                                                    &
         )
    Return
  End If

  ! Relative time frame.
  ! $$ Probably more logical to put calander time frame bit here too. Needs to be here for relative time frame
  ! otherwise data gets out of sync when data isn't valid.
  If (.not. IsCalendar()) Then

    If (SingleSiteMet%iRecord == 1) Then
      SingleSiteMet%Time2 = ReferenceTime()
    Else
      SingleSiteMet%Time2 = SingleSiteMet%Time2 + SingleSiteMet%dT
    End If

  EndIf

  If (Error) Then
    Call Message(                                         &
           'ERROR in reading Single Site met dataset ' // &
           Trim(SingleSiteMet%MetFile),                   &
           2                                              &
         )
    Return
  End If
  If (InadequateData) Return

  ! Update ProfileData2, Cloud2 and Ppt2. $$ warnings needed?
  Call SSMetToProfileData2(ProcMet, SingleSiteMet)
  If (ProcMet%Cl == -999.0) Then
    SingleSiteMet%Cloud2 = 0.0
  Else
    SingleSiteMet%Cloud2 = ProcMet%Cl / 8.0
  End If
  If (ProcMet%P == -999.0) Then
    SingleSiteMet%Ppt2 = 0.0
  Else
    SingleSiteMet%Ppt2 = ProcMet%P
  End If

  ! Calendar time frames.
  If (IsCalendar()) Then

    If (                          &
      ProcMet%Year == -999.0 .or. &
      ProcMet%Day  == -999.0 .or. &
      ProcMet%Hour == -999.0      &
    ) Then
      Call Message(                                                &
             'Error in ReadSS: met data must be time-labelled ' // &
             'for calculations using a calendar time frame',       &
             3                                                     &
           )
    Else
      If (ProcMet%TimeZone == -999.0) Then
        Call ControlledMessage(                                                &
               'Time zone not given with met data in single site met file ' // &
               '- local time based on longitude assumed.',                     &
               MessageControls    = GlobalMessageControls,                     &
               MessageControlName = 'No time zone in single site met file',    &
               ErrorCode          = 1                                          &
             )
        Zone = SingleSiteMet%Site%Long / 15.0
      Else
        Zone = ProcMet%TimeZone
      End If
      ZoneHour   = Int(Zone)
      ZoneMinute = Int((Zone - ZoneHour) * 60.0)
      ! Note ZoneHour and ZoneMinute, if both non-zero, have the same sign.
      If (ZoneHour == 0 .and. ZoneMinute == 0) Then
        SingleSiteMet%Time2 = Char2Time(                             &
                                '1/1/'                            // &
                                Trim(Int2Char(Int(ProcMet%Year))) // &
                                ' '                               // &
                                Trim(Int2Char(Int(ProcMet%Hour))) // &
                                ':00'                                &
                              )
      Else If (Zone > 0) Then
        SingleSiteMet%Time2 = Char2Time(                             &
                                '1/1/'                            // &
                                Trim(Int2Char(Int(ProcMet%Year))) // &
                                ' '                               // &
                                Trim(Int2Char(Int(ProcMet%Hour))) // &
                                ':00 UTC+'                        // &
                                Trim(Int2Char(Abs(ZoneHour)))     // &
                                ':'                               // &
                                Trim(Int2Char(Abs(ZoneMinute)))      &
                              )
      Else If (Zone < 0) Then
        SingleSiteMet%Time2 = Char2Time(                             &
                                '1/1/'                            // &
                                Trim(Int2Char(Int(ProcMet%Year))) // &
                                ' '                               // &
                                Trim(Int2Char(Int(ProcMet%Hour))) // &
                                ':00 UTC-'                        // &
                                Trim(Int2Char(Abs(ZoneHour)))     // &
                                ':'                               // &
                                Trim(Int2Char(Abs(ZoneMinute)))      &
                              )
      End If
      SingleSiteMet%Time2 = SingleSiteMet%Time2 +                                &
                            Char2Time(                                           &
                              Trim(Int2Char(Int(ProcMet%Day) - 1)) // 'd 00:00', &
                              Interval = .true.                                  &
                            )
    End If

  EndIf

End Subroutine ReadSS

!-------------------------------------------------------------------------------------------------------------

Subroutine SSMetToProfileData2(SSMet, SingleSiteMet)
! Calculates SingleSiteMet%ProfileData2 from SSMet and from coordinate, unresolved mesoscale motions and 
! free tropospheric turbulence information within SingleSiteMet.

  Implicit None
  ! Argument list:
  Type(SSMet_),         Intent(In)            :: SSMet
  Type(SingleSiteMet_), Intent(InOut), Target :: SingleSiteMet
  ! SSMet         :: A record of met data as processed by the single site module.
  ! SingleSiteMet :: State of a single site met module instance.
  ! Locals:
  Type(Flow_)                 :: Flow        ! Flow properties at boundary layer top plus some bulk boundary
                                             ! layer properties.
  Real(Std)                   :: Speed       ! Wind speed at boundary layer top.
  Real(Std)                   :: dThetadZ    ! d(potential temperature)/dz just above the boundary layer.
  Type(ProfileData_), Pointer :: ProfileData ! Abbreviation for SingleSiteMet%ProfileData2.

  ProfileData => SingleSiteMet%ProfileData2

  ! Note directions in ProfileData are anticlockwise from x-direction in radians (in [-pi,pi]),
  ! directions in SSMet are 'from bearing' in degrees (in [0,360)).
  ! $$ clarify conventions elsewhere and trap +/-Pi for DeltaPhi which wont work in MeanFlowProfiles

  ! Z0. Note the single site module returns the met site Z0 in SSMet%Z0 if this is the right Z0 to use and the
  ! dispersion area Z0 in SSMet%Z0D if this is the right Z0 to use.
  If (SSMet%Z0 /= -999.0 .and. SSMet%Z0D /= -999.0) Then
    Call Message('UNEXPECTED FATAL ERROR in SSMetToProfileData2', 4)
  Else If (SSMet%Z0 /= -999.0) Then
    ProfileData%Z0 = SSMet%Z0
  Else If (SSMet%Z0D /= -999.0) Then
    ProfileData%Z0 = SSMet%Z0D
  Else
    Call Message('UNEXPECTED FATAL ERROR in SSMetToProfileData2', 4)
  End If

  ProfileData%Phi0       = (- 90.0 - SSMet%Phi0) * Pi / 180.0
  ProfileData%PhiG       = (- 90.0 - SSMet%PhiG) * Pi / 180.0
  ProfileData%DeltaPhi   = - SSMet%DeltaPhi * Pi / 180.0
  Do While (ProfileData%Phi0 > Pi)
    ProfileData%Phi0 = ProfileData%Phi0 - 2.0*Pi
  End Do
  Do While (ProfileData%Phi0 < -Pi)
    ProfileData%Phi0 = ProfileData%Phi0 + 2.0*Pi
  End Do
  Do While (ProfileData%PhiG > Pi)
    ProfileData%PhiG = ProfileData%PhiG - 2.0*Pi
  End Do
  Do While (ProfileData%PhiG < -Pi)
    ProfileData%PhiG = ProfileData%PhiG + 2.0*Pi
  End Do
  Do While (ProfileData%DeltaPhi > Pi)
    ProfileData%DeltaPhi = ProfileData%DeltaPhi - 2.0*Pi
  End Do
  Do While (ProfileData%DeltaPhi < -Pi)
    ProfileData%DeltaPhi = ProfileData%DeltaPhi + 2.0*Pi
  End Do

  ProfileData%UStar      = SSMet%UStar
  ProfileData%UG         = SSMet%UG
  ProfileData%RecipLMO   = SSMet%RecipLMO
  ProfileData%H          = SSMet%H
  ProfileData%WStar      = SSMet%WStar

  ! $$ Use TSample in estimating unresolved mesoscale motions?.
  ProfileData%SigUUM     = SingleSiteMet%SigUUM
  ProfileData%TauUUM     = SingleSiteMet%TauUUM

  ProfileData%SigU2HPlus = SingleSiteMet%SigU2HPlus
  ProfileData%SigW2HPlus = SingleSiteMet%SigW2HPlus
  ProfileData%TauUHPlus  = SingleSiteMet%TauUHPlus
  ProfileData%TauWHPlus  = SingleSiteMet%TauWHPlus

  ProfileData%ZT         = 1.22  ! $$ hardwired.
  ProfileData%T0         = SSMet%T0K
  ProfileData%WT         = SSMet%FTheta0 / (1.225 * Cp) ! $$ rho consistent with ss module
                                                        ! (but not with rest of nameiii); cp not
  If (ProfileData%UStar /= 0.0) Then
    ProfileData%TStar    = - ProfileData%WT / ProfileData%UStar
  Else If (ProfileData%WT > 0.0) Then
    ProfileData%TStar    = - 200000.0
  Else
    ProfileData%TStar    = 200000.0
  End If
  ProfileData%Q0         = SSMet%Q0
  ProfileData%WQ         = SSMet%LambdaE / (1.225 * CalcLatentHeat(ProfileData%T0))
                                                        ! $$ consistency with ss module and rest of code
  If (ProfileData%UStar /= 0.0) Then
    ProfileData%QStar    = - ProfileData%WQ / ProfileData%UStar
  Else If (ProfileData%WQ > 0.0) Then
    ProfileData%QStar    = - 200000.0
  Else
    ProfileData%QStar    = 200000.0
  End If
  ProfileData%RHHPlus    = SSMet%RHU
  ProfileData%dRHdZHPlus = SSMet%dRHdZU
  ProfileData%PS         = PAt0km

  ProfileData%iHCoord = SingleSiteMet%C%iHCoords(1) ! (clarify nature of coords $$)
  ProfileData%iZCoord = SingleSiteMet%C%iZCoords(1)

  ProfileData%MeanFlow = .true.
  ProfileData%Turb     = .true.

  ! Get values just below boundary layer top.
  Call MeanFlowProfiles(                 &
         Z              = ProfileData%H, &
         ProfileData    = ProfileData,   &
         Moisture       = .false.,       &
         HMinus         = .true.,        &
         Speed          = Speed,         &
         Flow           = Flow           &
       )

  ProfileData%UG        = Speed ! $$ removing bl top jump in speed (but not in direction).
  ProfileData%THPlus    = Flow%T + SSMet%DeltaTheta
  ProfileData%PHPlus    = Flow%P
  dThetadZ              = (SSMet%NU ** 2) * Flow%T / Gravity
  ProfileData%dTdZHPlus = dThetadZ * (Flow%P / PRef) ** (GasConstant/Cp) - Gravity/Cp

End Subroutine SSMetToProfileData2

!-------------------------------------------------------------------------------------------------------------

End Module SingleSiteMetModule

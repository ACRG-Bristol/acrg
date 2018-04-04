! Module:  Radar Flow Module

Module RadarFlowModule
!.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowAndFlowProfileModule
Use RadarMetModule
Use MetsModule
Use CommonFlowModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: RadarFlow_                !
Public  :: InitRadarFlow             !
Public  :: SetUpRadarFlow_MetsCoordsEtc
Public  :: PrepareForUpdateRadarFlow ! Prepares for updating an instance of the radar flow module.
Public  :: RadarFlowReqs             !
Public  :: UpdateRadarFlow           !
Public  :: RadarRainCoordIndices     !
Public  :: GetRadarRain              !

!-------------------------------------------------------------------------------------------------------------

Type :: RadarFlow_ !
  Type (CommonFlow_)         :: C        ! The part of the flow state common to all flow modules.
  Type (RadarMet_),  Pointer :: RadarMet ! Stored radar rainfall data.
End Type RadarFlow_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitRadarFlow(          &
           FlowName,             &
           MetModName, MetName,  &
           DomainName, FixedMet, &
           UpdateOnDemand        &
         )
! Initialises an instance of the flow state.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: FlowName
  Character(*), Intent(In) :: MetModName
  Character(*), Intent(In) :: MetName
  Character(*), Intent(In) :: DomainName     ! Name of the domain of the flow module instance.
  Logical,      Intent(In) :: FixedMet
  Logical,      Intent(In) :: UpdateOnDemand ! Indicates the flow module instance is to be updated using
                                             ! update-on-demand.
  ! Function result:
  Type(RadarFlow_) InitRadarFlow !
  ! Locals:
  Type(RadarFlow_)         :: RadarFlow               !
  Integer                  :: nAttribs                !
  Character(MaxCharLength) :: Attribs(MaxFlowAttribs) !
  Character(MaxCharLength) :: MetModNames(MaxMets)    !
  Character(MaxCharLength) :: MetNames(MaxMets)       !

  ! Met module instances used.
  If (Len_Trim(MetModName) > MaxCharLength .or. Len_Trim(MetModName) == 0) Then
    Call Message('Error', 3)
  End If
  If (Len_Trim(MetName) > MaxCharLength .or. Len_Trim(MetName) == 0) Then
    Call Message('Error', 3)
  End If
  MetModNames(1) = MetModName
  MetNames(1)    = MetName
  If (MetModNames(1) /= 'Radar Met') Then
    Call Message('Error in InitRadarFlow', 3)
  End If

  ! Attributes.
  nAttribs   = 2
  Attribs(1) = 'Update'
  Attribs(2) = 'Rain'

  RadarFlow%C = InitCommonFlow(                   &
                  'Radar Flow', FlowName,         &
                  1, MetModNames, MetNames,       &
                  DomainName,                     &
                  nAttribs, Attribs,              &
                  FixedMet       = FixedMet,      &
                  UpdateOnDemand = UpdateOnDemand &
                )

  InitRadarFlow = RadarFlow

End Function InitRadarFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpRadarFlow_MetsCoordsEtc(EtaDefns, Coords, Grids, Mets, RadarFlow)
! Sets up RadarFlow using information from EtaDefns, Coords, Grids and Mets.

  Implicit None
  ! Argument list:
  Type(EtaDefns_),  Intent(In)    :: EtaDefns   ! Collection of eta definitions.
  Type(Coords_),    Intent(In)    :: Coords     ! Collection of coord systems.
  Type(Grids_),     Intent(In)    :: Grids      ! Collection of grids.
  Type(Mets_),      Intent(In)    :: Mets       ! Set of met module instance states.
  Type(RadarFlow_), Intent(InOut) :: RadarFlow  ! State of a radar flow module instance.
  ! Locals:
  Integer                  :: iMetMod     ! Met module index.
  Integer                  :: iMet        ! Met module instance index.
  Character(MaxCharLength) :: HCoordName1 !} Names of coord systems to add to flow module instance.
  Character(MaxCharLength) :: ZCoordName1 !}

  ! Find index of met module instance.
  Call FindMetIndex(                                          &
         RadarFlow%C%MetModNames(1), RadarFlow%C%MetNames(1), &
         Mets,                                                &
         iMetMod, iMet                                        &
       )

  ! Set up the names of the coord systems to add to RadarFlow%C.
  ! $$ should do for all met modules used or check coords (and grids?) the same
  HCoordName1 = Mets%C(iMetMod, iMet)%P%HCoordNames(1)
  ZCoordName1 = ''

  ! Add coord systems from met module instance to RadarFlow%C.
  Call AddCoordsToCommonFlow(                        &
         1, (/ HCoordName1 /), 0, (/ ZCoordName1 /), &
         RadarFlow%C                                 &
       )

  ! Check number of coord systems in RadarFlow.
  If (RadarFlow%C%nHCoords /= 1) Then
    Call Message(                                                  &
           'Unexpected error in SetUpRadarFlow_MetsCoordsEtc: ' // &
           'error in numbering of horizontal coord systems '    // &
           'used by the radar flow module instance "'           // &
           Trim(RadarFlow%C%FlowName)                           // &
           '"',                                                    &
           3                                                       &
         )
  End If
  If (RadarFlow%C%nZCoords /= 0) Then
    Call Message(                                                  &
           'Unexpected error in SetUpRadarFlow_MetsCoordsEtc: ' // &
           'error in numbering of vertical coord systems '      // &
           'used by the radar flow module instance "'           // &
           Trim(RadarFlow%C%FlowName)                           // &
           '"',                                                    &
           3                                                       &
         )
  End If

End Subroutine SetUpRadarFlow_MetsCoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForUpdateRadarFlow( &
             Coords, Grids, Domains,  &
             Mets,                    &
             iCase,                   &
             Time,                    &
             TValid, UpdateNow,       &
             RadarFlow,               &
             Units                    &
           )
! Prepares for updating an instance of the radar flow module.

! This routine must set TValid and UpdateNow but must not alter the validity of the radar flow module
! instance.
! $$ this comment should be with generic definition.

  Implicit None
  ! Argument list:
  Type(Coords_),    Intent(In)           :: Coords
  Type(Grids_),     Intent(In)           :: Grids
  Type(Domains_),   Intent(In),   Target :: Domains
  Type(Mets_),      Intent(In)           :: Mets
  Integer,          Intent(In)           :: iCase
  Type(Time_),      Intent(In)           :: Time
  Type(Time_),      Intent(Out)          :: TValid
  Logical,          Intent(Out)          :: UpdateNow
  Type(RadarFlow_), Intent(InOut)        :: RadarFlow
  Type(Units_),     Intent(InOut)        :: Units
  ! Coords    :: Collection of coord systems.
  ! Grids     :: Collection of grids.
  ! Domains   :: Collection of domains.
  ! Mets      :: Set of met module instance states.
  ! iCase     :: Number of case.
  ! Time      :: Time for which the radar flow module instance might be updated.
  ! TValid    :: Earliest time that the validity (overall or for any single attribute) of the flow module
  !              instance might change, except perhaps for a change caused by a change in the validity of the
  !              met and flow module instances acting as data sources, assuming the flow module instance is
  !              updated now. The value is that determined at the end of this routine (the actual time may be
  !              later).
  ! UpdateNow :: Indicates the flow module instance must be updated now (even if update-on-demand is
  !              specified). If set, TValid need not be set to any particular time.
  ! RadarFlow :: State of a radar flow module instance.
  ! Units     :: Collection of information on input/output unit numbers.
  ! Locals:
  Type(Domain_), Pointer :: Domain ! Abbreviation for domain.

  ! Abbreviations.
  Domain => Domains%Domains(RadarFlow%C%iDomain)

  ! Validity due to time domain.
  If (Time < StartTimeOfDomain(Domain)) Then
    TValid = StartTimeOfDomain(Domain)
  Else If (Time >= EndTimeOfDomain(Domain)) Then
    TValid = InfFutureTime()
  Else
    TValid = EndTimeOfDomain(Domain)
  End If
  UpdateNow = .false.

End Subroutine PrepareForUpdateRadarFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine RadarFlowReqs(                                              &
             RadarFlow, Time,                                          &
             WantFlowField,     WantCloudField,     WantRainField,     &
             iFlowUpdateSubset, iCloudUpdateSubset, iRainUpdateSubset, &
             FlowField,         CloudField,         RainField          &
           )
! .

  Implicit None
  ! Argument list:
  Type(RadarFlow_),  Intent(In)  :: RadarFlow         !
  Type(Time_),       Intent(In)  :: Time              !
  Logical,           Intent(Out) :: WantFlowField
  Logical,           Intent(Out) :: WantCloudField
  Logical,           Intent(Out) :: WantRainField
  Integer,           Intent(Out) :: iFlowUpdateSubset
  Integer,           Intent(Out) :: iCloudUpdateSubset
  Integer,           Intent(Out) :: iRainUpdateSubset
  Type(FlowField_),  Intent(Out) :: FlowField
  Type(CloudField_), Intent(Out) :: CloudField
  Type(RainField_),  Intent(Out) :: RainField

  WantFlowField = .false.
  ! Avoid compiler warning.
  iFlowUpdateSubset = 0
  FlowField%C%Valid = .false.

  WantCloudField = .false.
  ! Avoid compiler warning.
  iCloudUpdateSubset = 0
  CloudField%C%Valid = .false.

  WantRainField = .false.
  ! Avoid compiler warning.
  iRainUpdateSubset = 0
  RainField%C%Valid = .false.

  ! $$ Need con/dyn ppt ratio from nwp to partition rain and assume 50/50 if nwp
  ! values both zero (this is what NAME II does)

End Subroutine RadarFlowReqs

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateRadarFlow(                    &
             Coords, Grids, Domains,           &
             Mets,                             &
             iCase,                            &
             Time,                             &
             FlowField, CloudField, RainField, &
             RadarFlow,                        &
             Units                             &
           )
! Updates an instance of the flow state.

  Implicit None
  ! Argument list:
  Type(Coords_),            Intent(In)           :: Coords     !
  Type(Grids_),             Intent(In)           :: Grids      ! Collection of grids.
  Type(Domains_),           Intent(In),   Target :: Domains    ! Collection of domains.
  Type(Mets_),      Target, Intent(In)           :: Mets       !
  Integer,                  Intent(In)           :: iCase      !
  Type(Time_),              Intent(In)           :: Time       ! Time.
  Type(FlowField_),         Intent(In)           :: FlowField  !
  Type(CloudField_),        Intent(In)           :: CloudField !
  Type(RainField_),         Intent(In)           :: RainField  !
  Type(RadarFlow_),         Intent(InOut)        :: RadarFlow  ! Flow state to be updated.
  Type(Units_),             Intent(InOut)        :: Units      !
  ! Locals:
  Type(Domain_), Pointer :: Domain  ! Abbreviation for domain.
  Logical                :: Valid   !} Validity due to time domain.
  Type(Time_)            :: TValid  !}
  Logical                :: MValid  !] Validity due to met.
  Type(Time_)            :: MTValid !]

  ! Abbreviations.
  Domain => Domains%Domains(RadarFlow%C%iDomain)

  ! Validity due to time domain.
  Valid = StartTimeOfDomain(Domain) <= Time .and. Time < EndTimeOfDomain(Domain)
  If (Valid) Then
    TValid = EndTimeOfDomain(Domain)
  Else
    TValid = StartTimeOfDomain(Domain)
    If (Time >= TValid) TValid = InfFutureTime()
  End If

  ! Validity due to met.
  Call MetValid(Mets,                    &
                RadarFlow%C%iMetMods(1), &
                RadarFlow%C%iMets(1),    &
                Time, MValid, MTValid)

  ! Update met information.
  RadarFlow%RadarMet => Mets%RadarMets(RadarFlow%C%iMets(1))

  ! Overall validity.
  RadarFlow%C%Valid                     = Valid .and. MValid
  RadarFlow%C%ValidAttribParams(:)      = .false.
  RadarFlow%C%ValidAttribParams(A_Rain) = RadarFlow%C%Valid
  RadarFlow%C%ValidityUnimprovable      = RadarFlow%C%Valid .or. RadarFlow%C%FixedMet
  RadarFlow%C%TValid                    = TValid
  If (RadarFlow%C%Valid) RadarFlow%C%TValid = TMin(RadarFlow%C%TValid, MTValid)

End Subroutine UpdateRadarFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine RadarRainCoordIndices(RadarFlow, nHIndices, HIndex, nZIndices, ZIndex)
! .

  Implicit None
  ! Argument list:
  Type(RadarFlow_), Intent(In)  :: RadarFlow          !
  Integer,          Intent(Out) :: nHIndices          !
  Integer,          Intent(Out) :: HIndex(MaxHCoords) !
  Integer,          Intent(Out) :: nZIndices          !
  Integer,          Intent(Out) :: ZIndex(MaxHCoords) !

  nHIndices = 1
  HIndex(1) = RadarFlow%C%iHCoords(1)
  nZIndices = 0
  ZIndex(1) = RadarFlow%C%iZCoords(1)

End Subroutine RadarRainCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine GetRadarRain(Coords, Grids, RadarFlow, Time, Position, Rain)
! Extracts rain information at a particular location from the radar met fields.

  Implicit None
  ! Argument list:
  Type(Coords_),    Intent(In)           :: Coords
  Type(Grids_),     Intent(In),  Target  :: Grids
  Type(RadarFlow_), Intent(In),  Target  :: RadarFlow
  Type(ShortTime_), Intent(In)           :: Time
  Type(Position_),  Intent(In)           :: Position
  Type(Rain_),      Intent(Out)          :: Rain
  ! Locals:
  Real(Std)                :: X          !} Horizontal location in coord system of horizontal grids
  Real(Std)                :: Y          !} used for met data.
  Real(Std)                :: Ppt        ! Precipitation at the specified location and time.
  Type(RadarMet_), Pointer :: M          ! Abbreviation for RadarFlow%RadarMet.
  Type(HGrid_),    Pointer :: HGrid      ! Abbreviation for grid.
  Type(HCoeffs_)           :: HCoeffs    !} Interpolation coefficients.
  Type(TCoeffs_)           :: TCoeffs    !}
  Type(TCoeffs_)           :: TCoeffsEnd !}

  ! 1) Set up abbreviations for met and horizontal grid

  M     => RadarFlow%RadarMet

  HGrid => Grids%HGrids(M%iHGrid)

  ! 2) Calculate interpolation coefficients.

  ! Calculate temporal interpolation coefficients.
  Call GetTCoeffs(                                 &
         Time,                                     &
         M%OldTime, Int(ShortTime2RealTime(M%Dt)), &
         M%OldData, M%NewData,                     &
         TCoeffs                                   &
       )

  ! Calculate time interpolation coefficients for interpolating to a time at the end of the radar met time
  ! interval (note end here means in chronological time not run direction time - hence the use of
  ! IsBackwards).  The "If NextPrecip" is for efficiency as TCoeffsEnd is not used otherwise.
  If (M%MetDefn%NextPrecip) Then

      TCoeffsEnd = TCoeffs
      If (IsBackwards()) Then
        TCoeffsEnd%T1     = 0.0
        TCoeffsEnd%T2     = 1.0
      Else
        TCoeffsEnd%T1     = 1.0
        TCoeffsEnd%T2     = 0.0
      End If

  End If

  ! Calculate horizontal interpolation coefficients.

  X = Position%XY(1, RadarFlow%C%iHCoords(1))
  Y = Position%XY(2, RadarFlow%C%iHCoords(1))

  Call GetHCoeffs(X, Y, HGrid, HCoeffs)

  ! 3) Extract data from radar met array.

  ! If field refers to "NextPrecip" then interpolate in time to the end of the radar met time interval;
  ! otherwise use linear interpolation to the specified time. Spatial interpolation is bi-linear to the
  ! specified location.

  If (M%MetDefn%NextPrecip) Then
    Call InterpXYT(HCoeffs, TCoeffsEnd, M%Ppt, Ppt)
  Else
    Call InterpXYT(HCoeffs, TCoeffs,    M%Ppt, Ppt)
  End If

! $$ Split precipitation evenly into convective and dynamic components.

  Rain%DynPpt = 0.5 * Ppt
  Rain%ConPpt = 0.5 * Ppt

End Subroutine GetRadarRain

!-------------------------------------------------------------------------------------------------------------

End Module RadarFlowModule

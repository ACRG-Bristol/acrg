! Module: Prototype Flow Module

Module PrototypeFlowModule

! This is a simple prototype flow module.
! It requires input from a single instance of the prototype met module.
! It has attributes Convert and Flow.
! The index of the coords needed for GetFlow and mean advect ($$) is 1.
! It assumes no topography and the pressure structure of the ICAO standard atmosphere.
! It returns default values for all bulk boundary layer quantities other than H. The
! defaults are Z0 = 0.1, UStar = 0.4, HeatFlux = 0.0 and RecipLMO = 0.0. The values of
! these quantities are not used for dispersion per se.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowAndFlowProfileModule
Use PrototypeMetModule
Use MetsModule
Use CommonFlowModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: PrototypeFlow_                !
Public  :: InitPrototypeFlow             !
Public  :: PrepareForUpdatePrototypeFlow ! Prepares for updating an instance of the prototype flow module.
Public  :: PrototypeFlowReqs             !
Public  :: UpdatePrototypeFlow           !
Public  :: PrototypeConvertCoordIndices  !
Public  :: PrototypeConvertToZ           !
Public  :: PrototypeFlowCoordIndices     !
Public  :: GetPrototypeFlow              !
Public  :: PrototypeReflectCoordIndices  ! Returns indices of coord systems in which
                                         ! positions need to be specified when using
                                         ! PrototypeReflect.
Public  :: PrototypeReflect              !

!-------------------------------------------------------------------------------------------------------------

Type :: PrototypeFlow_ !
  Type (CommonFlow_) C   ! The part of the flow state common to all flow modules.
  Type (Met_)        Met ! Stored met data.
End Type PrototypeFlow_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitPrototypeFlow(        &
           FlowName,               &
           MetModName, MetName,    &
           HCoordName, ZCoordName, &
           DomainName, FixedMet,   &
           UpdateOnDemand          &
         )
! Initialises an instance of the flow state.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: FlowName
  Character(*), Intent(In) :: MetModName
  Character(*), Intent(In) :: MetName
  Character(*), Intent(In) :: HCoordName
  Character(*), Intent(In) :: ZCoordName
  Character(*), Intent(In) :: DomainName     ! Name of the domain of the flow module instance.
  Logical,      Intent(In) :: FixedMet
  Logical,      Intent(In) :: UpdateOnDemand ! Indicates the flow module instance is to be updated using
                                             ! update-on-demand.
  ! Function result:
  Type(PrototypeFlow_) InitPrototypeFlow !
  ! Locals:
  Type(PrototypeFlow_) PrototypeFlow  !
  Integer              nAttribs            !
  Character(MaxCharLength) Attribs(MaxFlowAttribs) !
  Character(MaxCharLength) HCoordNames(MaxHCoords)
  Character(MaxCharLength) ZCoordNames(MaxZCoords)
  Character(MaxCharLength) MetModNames(MaxMets)
  Character(MaxCharLength) MetNames(MaxMets)

  ! Met module instances used.
  If (Len_Trim(MetModName) > MaxCharLength .or. Len_Trim(MetModName) == 0) Then
    Call Message('Error', 3)
  End If
  If (Len_Trim(MetName) > MaxCharLength .or. Len_Trim(MetName) == 0) Then
    Call Message('Error', 3)
  End If
  MetModNames(1) = MetModName
  MetNames(1)    = MetName
  If (MetModNames(1) /= 'Prototype Met') Then
    Call Message('Error in InitPrototypeFlow', 3)
  End If

  ! Attributes.
  nAttribs   = 3
  Attribs(1) = 'Update'
  Attribs(2) = 'Convert'
  Attribs(3) = 'Flow'

  If (Len_Trim(HCoordName) > MaxCharLength .or. Len_Trim(HCoordName) == 0) Then
    Call Message('Error', 3)
  End If
  If (Len_Trim(ZCoordName) > MaxCharLength .or. Len_Trim(ZCoordName) == 0) Then
    Call Message('Error', 3)
  End If
  HCoordNames(1) = HCoordName
  ZCoordNames(1) = ZCoordName

  PrototypeFlow%C = InitCommonFlow(                   &
                      'Prototype Flow', FlowName,     &
                      1, MetModNames, MetNames,       &
                      DomainName,                     &
                      nAttribs, Attribs,              &
                      FixedMet       = FixedMet,      &
                      UpdateOnDemand = UpdateOnDemand &
                    )
  Call AddCoordsToCommonFlow(            &
         1, HCoordNames, 1, ZCoordNames, &
         PrototypeFlow%C                 &
       )

  InitPrototypeFlow = PrototypeFlow

End Function InitPrototypeFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForUpdatePrototypeFlow( &
             Coords, Grids, Domains,      &
             Mets,                        &
             iCase,                       &
             Time,                        &
             TValid, UpdateNow,           &
             PrototypeFlow,               &
             Units                        &
           )
! Prepares for updating an instance of the prototype flow module.

! This routine must set TValid and UpdateNow but must not alter the validity of the prototype flow module
! instance.
! $$ this comment should be with generic definition.

  Implicit None
  ! Argument list:
  Type(Coords_),        Intent(In)           :: Coords
  Type(Grids_),         Intent(In)           :: Grids
  Type(Domains_),       Intent(In),   Target :: Domains
  Type(Mets_),          Intent(In)           :: Mets
  Integer,              Intent(In)           :: iCase
  Type(Time_),          Intent(In)           :: Time
  Type(Time_),          Intent(Out)          :: TValid
  Logical,              Intent(Out)          :: UpdateNow
  Type(PrototypeFlow_), Intent(InOut)        :: PrototypeFlow
  Type(Units_),         Intent(InOut)        :: Units
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! Domains       :: Collection of domains.
  ! Mets          :: Set of met module instance states.
  ! iCase         :: Number of case.
  ! Time          :: Time for which the prototype flow module instance might be updated.
  ! TValid        :: Earliest time that the validity (overall or for any single attribute) of the flow module
  !                  instance might change, except perhaps for a change caused by a change in the validity of
  !                  the met and flow module instances acting as data sources, assuming the flow module
  !                  instance is updated now. The value is that determined at the end of this routine (the
  !                  actual time may be later).
  ! UpdateNow     :: Indicates the flow module instance must be updated now (even if update-on-demand is
  !                  specified). If set, TValid need not be set to any particular time.
  ! PrototypeFlow :: State of a prototype flow module instance.
  ! Units         :: Collection of information on input/output unit numbers.
  ! Locals:
  Type(Domain_), Pointer :: Domain ! Abbreviation for domain.

  ! Abbreviations.
  Domain => Domains%Domains(PrototypeFlow%C%iDomain)

  ! Validity due to time domain.
  If (Time < StartTimeOfDomain(Domain)) Then
    TValid = StartTimeOfDomain(Domain)
  Else If (Time >= EndTimeOfDomain(Domain)) Then
    TValid = InfFutureTime()
  Else
    TValid = EndTimeOfDomain(Domain)
  End If
  UpdateNow = .false.

End Subroutine PrepareForUpdatePrototypeFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine PrototypeFlowReqs(                                          &
             PrototypeFlow, Time,                                      &
             WantFlowField,     WantCloudField,     WantRainField,     &
             iFlowUpdateSubset, iCloudUpdateSubset, iRainUpdateSubset, &
             FlowField,         CloudField,         RainField          &
           )
! .

  Implicit None
  ! Argument list:
  Type(PrototypeFlow_), Intent(In)  :: PrototypeFlow !
  Type(Time_),          Intent(In)  :: Time
  Logical,              Intent(Out) :: WantFlowField
  Logical,              Intent(Out) :: WantCloudField
  Logical,              Intent(Out) :: WantRainField
  Integer,              Intent(Out) :: iFlowUpdateSubset
  Integer,              Intent(Out) :: iCloudUpdateSubset
  Integer,              Intent(Out) :: iRainUpdateSubset
  Type(FlowField_),     Intent(Out) :: FlowField
  Type(CloudField_),    Intent(Out) :: CloudField
  Type(RainField_),     Intent(Out) :: RainField

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

End Subroutine PrototypeFlowReqs

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdatePrototypeFlow(                &
             Coords, Grids, Domains,           &
             Mets,                             &
             iCase,                            &
             Time,                             &
             FlowField, CloudField, RainField, &
             PrototypeFlow,                    &
             Units                             &
           )
! Updates an instance of the flow state.

  Implicit None
  ! Argument list:
  Type(Coords_),        Intent(In)           :: Coords        !
  Type(Grids_),         Intent(In)           :: Grids         ! Collection of grids.
  Type(Domains_),       Intent(In),   Target :: Domains       ! Collection of domains.
  Type(Mets_),          Intent(In)           :: Mets          !
  Integer,              Intent(In)           :: iCase         !
  Type(Time_),          Intent(In)           :: Time          ! Time.
  Type(FlowField_),     Intent(In)           :: FlowField
  Type(CloudField_),    Intent(In)           :: CloudField    !
  Type(RainField_),     Intent(In)           :: RainField     !
  Type(PrototypeFlow_), Intent(InOut)        :: PrototypeFlow ! Flow state to be updated.
  Type(Units_),         Intent(InOut)        :: Units         !
  ! Locals:
  Type(Domain_), Pointer :: Domain  ! Abbreviation for domain.
  Logical                :: Valid   !} Validity due to time domain.
  Type(Time_)            :: TValid  !}
  Logical                :: MValid  !] Validity due to met.
  Type(Time_)            :: MTValid !]

  ! Abbreviations.
  Domain => Domains%Domains(PrototypeFlow%C%iDomain)

  ! Validity due to time domain.
  Valid = StartTimeOfDomain(Domain) <= Time .and. Time < EndTimeOfDomain(Domain)
  If (Valid) Then
    TValid = EndTimeOfDomain(Domain)
  Else
    TValid = StartTimeOfDomain(Domain)
    If (Time >= TValid) TValid = InfFutureTime()
  End If

  ! Validity due to met.
  Call MetValid(Mets,                        &
                PrototypeFlow%C%iMetMods(1), &
                PrototypeFlow%C%iMets(1),    &
                Time, MValid, MTValid)

  ! Update met information.
  PrototypeFlow%Met = Mets%PrototypeMets(PrototypeFlow%C%iMets(1))%Met

  ! Overall validity.
  PrototypeFlow%C%Valid                        = Valid .and. MValid
  PrototypeFlow%C%ValidAttribParams(:)         = .false.
  PrototypeFlow%C%ValidAttribParams(A_Convert) = PrototypeFlow%C%Valid
  PrototypeFlow%C%ValidAttribParams(A_Flow   ) = PrototypeFlow%C%Valid
  PrototypeFlow%C%ValidityUnimprovable         = PrototypeFlow%C%Valid .or. PrototypeFlow%C%FixedMet
  PrototypeFlow%C%TValid                       = TValid
  If (PrototypeFlow%C%Valid) PrototypeFlow%C%TValid = TMin(PrototypeFlow%C%TValid, MTValid)

End Subroutine UpdatePrototypeFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine PrototypeConvertCoordIndices(PrototypeFlow, nHIndices, HIndex)
! .

  Implicit None
  ! Argument list:
  Type(PrototypeFlow_), Intent(In)  :: PrototypeFlow      !
  Integer,              Intent(Out) :: nHIndices          !
  Integer,              Intent(Out) :: HIndex(MaxHCoords) !

  nHIndices = 0
  HIndex(1) = PrototypeFlow%C%iHCoords(1)

End Subroutine PrototypeConvertCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine PrototypeConvertToZ(            &
             PrototypeFlow, Coords, Grids, &
             ZIndex, Time, Position        &
           )
! Converts vertical coords between coord systems.

  Implicit None
  ! Argument list:
  Type(PrototypeFlow_),         Intent(In)    :: PrototypeFlow !
  Type(Coords_),        Target, Intent(In)    :: Coords        !
  Type(Grids_),                 Intent(In)    :: Grids
  Integer,                      Intent(In)    :: ZIndex        !
  Type(ShortTime_),             Intent(In)    :: Time
  Type(Position_),              Intent(InOut) :: Position
  ! Locals:
  Type(ZCoord_), Pointer :: ZCoordIn  ! Abreviation.
  Type(ZCoord_), Pointer :: ZCoordOut ! Abreviation.
  Type(ZCoord_)          :: ZCoord1   !
  Type(ZCoord_)          :: ZCoord2   !
  Real(Std)              :: ZIn       !
  Real(Std)              :: ZOut      !
  Real(Std)              :: Z1        !
  Real(Std)              :: Z2        !
  Integer                :: ZUseIndex !
  Integer                :: i         !

  Do i = 1, Coords%nZCoords ! replace by efficient selection
    If (Position%ZValid(i)) Then
      ZUseIndex = i
      Exit
    End If
  End Do

  ZCoordIn  => Coords%ZCoords(ZUseIndex)
  ZCoordOut => Coords%ZCoords(ZIndex)
  ZIn       =  Position%Z(ZUseIndex)

  Select Case (ZCoordIn%CoordType)

    Case (Z_AboveGround, Z_AboveSea)

      Select Case (ZCoordOut%CoordType)
        Case (Z_AboveGround, Z_AboveSea)
          ZOut    = ZBasedToZBased(ZCoordIn, ZCoordOut, 0.0, ZIn)
        Case (Z_P, Z_PAsZ, Z_PAsEta)
          ZCoord1 = ZCoord_m_asl()
          ZCoord2 = InitZCoord('PAsZ', Z_PAsZ, 1.0)
          Z1      = ZBasedToZBased    (ZCoordIn, ZCoord1,   0.0,    ZIn)
          Z2      = AboveSeaToPAsZICAO(ZCoord1,  ZCoord2,           Z1 )
          ZOut    = PBasedToPBased    (ZCoord2,  ZCoordOut, PAt0km, Z2 )
        Case Default
          Call Message('Error in PrototypeConvertToZ', 4)
      End Select

    Case (Z_P, Z_PAsZ, Z_PAsEta)

      Select Case (ZCoordOut%CoordType)
        Case (Z_AboveGround, Z_AboveSea)
          ZCoord1 = InitZCoord('Z_PAsZ', Z_PAsZ, 1.0)
          ZCoord2 = ZCoord_m_asl()
          Z1      = PBasedToPBased    (ZCoordIn, ZCoord1,   PAt0km, ZIn)
          Z2      = PAsZToAboveSeaICAO(ZCoord1,  ZCoord2,           Z1 )
          ZOut    = ZBasedToZBased    (ZCoord2,  ZCoordOut, 0.0,    Z2 )
        Case (Z_P, Z_PAsZ, Z_PAsEta)
          ZOut    = PBasedToPBased(ZCoordIn, ZCoordOut, 0.0, ZIn)
        Case Default
          Call Message('Error in PrototypeConvertToZ', 4)
      End Select

    Case Default

      Call Message('Error in PrototypeConvertToZ', 4)

  End Select

  Position%Z(ZIndex) = ZOut

End Subroutine PrototypeConvertToZ

!-------------------------------------------------------------------------------------------------------------

Subroutine PrototypeFlowCoordIndices(PrototypeFlow,                        &
                                     nHIndices, HIndex, nZIndices, ZIndex)
! .

  Implicit None
  ! Argument list:
  Type(PrototypeFlow_), Intent(In)  :: PrototypeFlow      !
  Integer,              Intent(Out) :: nHIndices          !
  Integer,              Intent(Out) :: HIndex(MaxHCoords) !
  Integer,              Intent(Out) :: nZIndices          !
  Integer,              Intent(Out) :: ZIndex(MaxHCoords) !

  nHIndices = 1
  HIndex(1) = PrototypeFlow%C%iHCoords(1)
  nZIndices = 1
  ZIndex(1) = PrototypeFlow%C%iZCoords(1)

End Subroutine PrototypeFlowCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine GetPrototypeFlow(Coords, Grids, PrototypeFlow, &
                            Moisture,       &
                            Inhomog, Homog, &
                            Time, Position, &
                            Flow, ProfileData)
! .

  Implicit None
  ! Argument list:
  Type(Coords_),        Intent(In)  :: Coords
  Type(Grids_),         Intent(In)  :: Grids
  Logical,              Intent(In)  :: Moisture      ! Indicates Q is required.
  Logical,              Intent(In)  :: Inhomog       ! Indicates inhomogeneous
                                                     ! quantities are required.
  Logical,              Intent(In)  :: Homog         ! Indicates homogeneous
                                                     ! quantities are required.
  Type(ShortTime_),     Intent(In)  :: Time
  Type(Position_),      Intent(In)  :: Position      !
  Type(PrototypeFlow_), Intent(In)  :: PrototypeFlow !
  Type(Flow_),          Intent(Out) :: Flow          !
  Type(ProfileData_), Optional, Intent(Out) :: ProfileData
  ! Locals:
  Real(Std) ZeroPoint !
  Real(Std) XYL(2)
  Real(Std) ZL
  Real(Std) dUdZ
  Real(Std) dVdZ
  Real(Std) dDeltaIdZ !
  Real(Std) dTauWDZ

  Logical Vel

  XYL = Position%XY(:,PrototypeFlow%C%iHCoords(1))
  ZL  = Position%Z(PrototypeFlow%C%iZCoords(1))

  Flow%iHCoord = PrototypeFlow%C%iHCoords(1)
  Flow%iZCoord = PrototypeFlow%C%iZCoords(1)

  Call Interp(PrototypeFlow%Met%U,                          &
              ZL, PrototypeFlow%Met%Z00,                    &
              PrototypeFlow%Met%H, PrototypeFlow%Met%HPlus, &
              Flow%U(1), dUdZ)
  Call Interp(PrototypeFlow%Met%V,                          &
              ZL, PrototypeFlow%Met%Z00,                    &
              PrototypeFlow%Met%H, PrototypeFlow%Met%HPlus, &
              Flow%U(2), dVdZ)
  Flow%U(3) = 0.0

  Call Interp(PrototypeFlow%Met%DeltaI,                     &
              ZL, PrototypeFlow%Met%Z00,                    &
              PrototypeFlow%Met%H, PrototypeFlow%Met%HPlus, &
              Flow%DeltaI, dDeltaIdZ)

!  Vel = .false. ! $$ Temp fixup (particle%Vel not available)
                ! Now (July 03) changed so both blocks are calculated - at least
                ! we get siguu values returned, but needs rationalisation

!  If (Vel) Then

    Call Interp(PrototypeFlow%Met%SigW2,                      &
                ZL, PrototypeFlow%Met%Z00,                    &
                PrototypeFlow%Met%H, PrototypeFlow%Met%HPlus, &
                Flow%SigUU(3), Flow%dSigUUdX(3))
    Call Interp(PrototypeFlow%Met%TauW,                       &
                ZL, PrototypeFlow%Met%Z00,                    &
                PrototypeFlow%Met%H, PrototypeFlow%Met%HPlus, &
                Flow%TauUU(3), dTauWdZ)
    Flow%SigUU(1) = 1.0 ! $$
    Flow%SigUU(2) = 1.0 ! $$
    Flow%TauUU(1) = 1000.0 ! $$
    Flow%TauUU(2) = 1000.0 ! $$
    Flow%Sk       = 0.0 ! $$

    Flow%MaxDZ = Huge(Flow%MaxDZ)
    If (PrototypeFlow%Met%TauW%H > PrototypeFlow%Met%TauW%HPlus) Then
      ZeroPoint = (PrototypeFlow%Met%HPlus - PrototypeFlow%Met%H)*            &
                  PrototypeFlow%Met%TauW%HPlus/                               &
                  (PrototypeFlow%Met%TauW%H - PrototypeFlow%Met%TauW%HPlus) + &
                  PrototypeFlow%Met%HPlus
      If (ZL <= PrototypeFlow%Met%HPlus) Then
        Flow%MaxDZ = Min(ZeroPoint - ZL, Flow%MaxDZ)
      Else
        Flow%MaxDZ = Min(ZL - 2.0*PrototypeFlow%Met%HPlus + ZeroPoint, Flow%MaxDZ)
      End If
    End If
    If (PrototypeFlow%Met%SigW2%H > PrototypeFlow%Met%SigW2%HPlus) Then
      ZeroPoint = (PrototypeFlow%Met%HPlus - PrototypeFlow%Met%H)*              &
                  PrototypeFlow%Met%SigW2%HPlus/                                &
                  (PrototypeFlow%Met%SigW2%H - PrototypeFlow%Met%SigW2%HPlus) + &
                  PrototypeFlow%Met%HPlus
      If (ZL <= PrototypeFlow%Met%HPlus) Then
        Flow%MaxDZ = Min(ZeroPoint - ZL, Flow%MaxDZ)
      Else
        Flow%MaxDZ = Min(ZL - 2.0*PrototypeFlow%Met%HPlus + ZeroPoint, Flow%MaxDZ)
      End If
    End If

    Flow%MaxDZ = Min(                                    &
                   Flow%MaxDZ,                           &
                   Flow%TauUU(3) / Abs(dTauWdZ),         &
                   Flow%SigUU(3) / Abs(Flow%dSigUUdX(3)) &
                 )

!  Else

    Call Interp(PrototypeFlow%Met%KZ,                         &
                ZL, PrototypeFlow%Met%Z00,                    &
                PrototypeFlow%Met%H, PrototypeFlow%Met%HPlus, &
                Flow%K(3), Flow%dKdX(3))
    Flow%K(1) = 1.0
    Flow%K(2) = 1.0

    Flow%MaxDZ = Huge(Flow%MaxDZ)
    If (PrototypeFlow%Met%KZ%H > PrototypeFlow%Met%KZ%HPlus) Then
      ZeroPoint = (PrototypeFlow%Met%HPlus - PrototypeFlow%Met%H)*        &
                  PrototypeFlow%Met%KZ%HPlus/                             &
                  (PrototypeFlow%Met%KZ%H - PrototypeFlow%Met%KZ%HPlus) + &
                  PrototypeFlow%Met%HPlus
      If (ZL <= PrototypeFlow%Met%HPlus) Then
        Flow%MaxDZ = Min(ZeroPoint - ZL, Flow%MaxDZ)
      Else
        Flow%MaxDZ = Min(ZL - 2.0*PrototypeFlow%Met%HPlus + ZeroPoint, Flow%MaxDZ)
      End If
    End If

    If (Flow%dKdX(3) /= 0.0) Flow%MaxDZ = Min(Flow%MaxDZ, Flow%K(3) / Abs(Flow%dKdX(3)))

!  End If

  Flow%Z0       = 0.1
  Flow%UStar    = 0.4
  Flow%RecipLMO = 0.0
  Flow%H        = PrototypeFlow%Met%H

  Flow%dUdT(:)     = 0.0 ! $$

  Flow%T = TAt0km ! $$ Avoids failure in calculating relative humidity for foot and
  Flow%P = PAt0km ! mouth decay.

  If (Present(ProfileData)) Then
    ProfileData%WT       = 0.0
    ProfileData%MeanFlow = .false.
    ProfileData%Turb     = .false.
    ProfileData%Canopy   = .false.
  End If

End Subroutine GetPrototypeFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine PrototypeReflectCoordIndices(                           &
             PrototypeFlow, nHCoords, iHCoords, nZCoords, iZCoords &
           )
! Returns indices of coord systems in which positions need to be specified when using
! PrototypeReflect.

  Implicit None
  ! Argument list:
  Type(PrototypeFlow_), Intent(In)  :: PrototypeFlow        ! State of a prototype
                                                            ! flow module instance.
  Integer,              Intent(Out) :: nHCoords             !} Number and values of
  Integer,              Intent(Out) :: iHCoords(MaxHCoords) !} indices of horizontal
  Integer,              Intent(Out) :: nZCoords             !} and vertical coord
  Integer,              Intent(Out) :: iZCoords(MaxHCoords) !} systems.

  ! Horizontal coords needed: none
  nHCoords    = 0
  iHCoords(1) = 0 ! avoids compiler warning

  ! Vertical coords needed:
  ! 1) m agl.
  nZCoords    = 1
  iZCoords(1) = PrototypeFlow%C%iZCoords(1)

End Subroutine PrototypeReflectCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine PrototypeReflect(                &
             Coords, PrototypeFlow,         &
             Vel, OldPosition, Position, U, &
             Reflected                      &
           )
! .

  Implicit None
  ! Argument list:
  Type(Coords_),        Intent(In)    :: Coords     ! Collection of coord systems.
  Type(PrototypeFlow_), Intent(In)    :: PrototypeFlow !
  Logical,              Intent(In)    :: Vel           !
  Type(Position_),      Intent(In)    :: OldPosition
  Type(Position_),      Intent(InOut) :: Position
  Real(Std),            Intent(InOut) :: U(3)          !
  Logical,              Intent(Out)   :: Reflected
  ! Locals:
  Integer :: iZCoordZ ! Index of m agl coord system.

  iZCoordZ = PrototypeFlow%C%iZCoords(1)

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

End Subroutine PrototypeReflect

!-------------------------------------------------------------------------------------------------------------

Subroutine Interp(Levels, Z, Z00, H, HPlus, Value, dValuedZ)
! .

  Implicit None
  ! Argument list:
  Type(Levels_), Intent(In)  :: Levels   !
  Real(Std),     Intent(In)  :: Z        !
  Real(Std),     Intent(In)  :: Z00      !
  Real(Std),     Intent(In)  :: H        !
  Real(Std),     Intent(In)  :: HPlus    !
  Real(Std),     Intent(Out) :: Value    !
  Real(Std),     Intent(Out) :: dValuedZ !

  If (Z > HPlus) Then
    Value    = Levels%HPlus
    dValuedZ = 0.0
  Else If (Z > H) Then
    dValuedZ = (Levels%HPlus - Levels%H)/(HPlus - H)
    Value    = Levels%H + (Z - H)*dValuedZ
  Else
    dValuedZ = (Levels%H - Levels%Z)/H
    Value    = Levels%Z + Max(Z, Z00)*dValuedZ
    If (Z <= Z00) dValuedZ = 0.0
  End If

End Subroutine Interp

!-------------------------------------------------------------------------------------------------------------

End Module PrototypeFlowModule

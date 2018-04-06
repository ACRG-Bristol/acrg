! Module:  Building Flow Module

Module BuildingFlowModule

! This module calculates the flow around buildings.
! It doesn't require any input from met modules.
! It has attribute Flow.
! The index of the coords needed for GetFlow and mean advect ($$) is 1.
! It assumes no topography and the pressure structure of the ICAO standard atmosphere.

! H/Z Grid(1) is for getting flow info from other flow modules.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowAndFlowProfileModule
Use BuildingModule
Use MetsModule
Use CommonFlowModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: BuildingFlow_                   !
Public  :: InitBuildingFlow                !
Public  :: SetUpCoordsEtc_BuildingFlow     !
Public  :: SetUpBuildingFlow_MetsCoordsEtc
Public  :: PrepareForUpdateBuildingFlow    ! Prepares for updating an instance of the building flow module.
Public  :: BuildingFlowReqs                !
Public  :: UpdateBuildingFlow              !
Public  :: BuildingFlowCoordIndices        !
Public  :: GetBuildingFlow                 !
Public  :: BuildingReflectCoordIndices     ! Returns indices of coord systems in which
                                           ! positions need to be specified when using
                                           ! BuildingReflect.
Public  :: BuildingReflect                 !

!-------------------------------------------------------------------------------------------------------------

Type :: BuildingFlow_ !
  Type (CommonFlow_) :: C ! The part of the flow state common to all flow modules.
  Type (Flow_State)  :: flow_state_building
  Real(Std)          :: Xb
  Real(Std)          :: Yb
  Real(Std)          :: BuildingOrientation ! Acute angle of building length
                                            ! to X axis
  Real(Std)          :: WindDirection ! Angle wind is blowing to in degrees
                                      ! anticlockwise from x-axis.
  Type(Flow_)        :: Flow(1,1,1,1)
  Type(ProfileData_) :: ProfileData(1,1,1,1)
End Type BuildingFlow_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitBuildingFlow(              &
           FlowName,                    &
           HCoordName, ZCoordName,      &
           DomainName,                  &
           FixedMet,                    &
           UpdateOnDemand,              &
           Lb, Wb, Hb,                  &
           Xb, Yb, BuildingOrientation, &
           UpdateSubsetName             &
         )
! Initialises an instance of the flow state.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: FlowName
  Character(*), Intent(In) :: HCoordName
  Character(*), Intent(In) :: ZCoordName
  Character(*), Intent(In) :: DomainName           ! Name of the domain of the flow module instance.
  Logical,      Intent(In) :: FixedMet
  Logical,      Intent(In) :: UpdateOnDemand       ! Indicates the flow module instance is to be updated using
                                                   ! update-on-demand.
  Real(Std),    Intent(In) :: Lb
  Real(Std),    Intent(In) :: Wb
  Real(Std),    Intent(In) :: Hb
  Real(Std),    Intent(In) :: Xb
  Real(Std),    Intent(In) :: Yb
  Real(Std),    Intent(In) :: BuildingOrientation
  Character(*), Intent(In) :: UpdateSubsetName
  ! Function result:
  Type(BuildingFlow_) :: InitBuildingFlow !
  ! Locals:
  Type(BuildingFlow_)      :: BuildingFlow            !
  Integer                  :: nAttribs                !
  Character(MaxCharLength) :: Attribs(MaxFlowAttribs) !
  Character(MaxCharLength) :: HCoordNames(MaxHCoords)
  Character(MaxCharLength) :: ZCoordNames(MaxZCoords)
  Character(MaxCharLength) :: MetModNames(MaxMets)
  Character(MaxCharLength) :: MetNames(MaxMets)
  Logical                  :: UseUpdateSubset

  ! Attributes.
  nAttribs   = 2
  Attribs(1) = 'Update'
  Attribs(2) = 'Flow'

  If (Len_Trim(HCoordName) > MaxCharLength .or. Len_Trim(HCoordName) == 0) Then
    Call Message('Error', 3)
  End If
  If (Len_Trim(ZCoordName) > MaxCharLength .or. Len_Trim(ZCoordName) == 0) Then
    Call Message('Error', 3)
  End If
  HCoordNames(1) = HCoordName
  ZCoordNames(1) = ZCoordName

  If (Len_Trim(UpdateSubsetName) > MaxCharLength) Then
    Call Message('Error', 3)
  End If
  If (Len_Trim(UpdateSubsetName) == 0) Then
    UseUpdateSubset = .false.    ! must be true $$.
  Else
    UseUpdateSubset = .true.
  End If

  BuildingFlow%C = InitCommonFlow(                   &
                     'Building Flow', FlowName,      &
                     0, MetModNames, MetNames,       &
                     DomainName,                     &
                     nAttribs, Attribs,              &
                     FixedMet       = FixedMet,      &
                     UpdateOnDemand = UpdateOnDemand &
                   )

  Call AddCoordsToCommonFlow(            &
         1, HCoordNames, 1, ZCoordNames, &
         BuildingFlow%C                  &
       )

  Call AddUpdateSubsetsToCommonFlow( &
         1, (/ UpdateSubsetName /),  &
         BuildingFlow%C              &
       )

! Initialise building dimensions
  BuildingFlow%flow_state_building%cuboid%hb = Hb
  BuildingFlow%flow_state_building%cuboid%wb = Wb
  BuildingFlow%flow_state_building%cuboid%lb = Lb

! Initialise (x,y) coordinates of centre of building
  BuildingFlow%Xb = Xb
  BuildingFlow%Yb = Yb

! Initialise building orientation
  BuildingFlow%BuildingOrientation = BuildingOrientation

  InitBuildingFlow = BuildingFlow

End Function InitBuildingFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpCoordsEtc_BuildingFlow(BuildingFlow, Coords, Grids)
!

  Implicit None
  ! Argument list:
  Type(BuildingFlow_), Intent(In)    :: BuildingFlow
  Type(Coords_),       Intent(InOut) :: Coords
  Type(Grids_),        Intent(InOut) :: Grids
  ! Locals:
  Type(HGrid_) :: HGrid
  Type(ZGrid_) :: ZGrid

  HGrid = InitHGrid(                                    &
            Name       = 'Flow for buildings',          &
            HCoordName = BuildingFlow%C%HCoordNames(1), &
            Wrap       = .false.,                       &
            nX         = 1,                             &
            nY         = 1,                             &
            dX         = 1.0,                           &
            dY         = 1.0,                           &
            X0         = BuildingFlow%XB,               &
            Y0         = BuildingFlow%YB                &
          )

  Call AddHGrid(HGrid, Grids)

  ZGrid = InitZGrid(                                    &
            Name       = 'Flow for buildings',          &
            ZCoordName = BuildingFlow%C%ZCoordNames(1), &
            nZ         = 1,                             &
            dZ         = 1.0,                           &
            Z0         = 0.0                            &
          ) ! Need more points for when profile data invalid $$

  Call AddZGrid(ZGrid, Grids)

End Subroutine SetUpCoordsEtc_BuildingFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpBuildingFlow_MetsCoordsEtc(EtaDefns, Coords, Grids, Mets, BuildingFlow)
! Sets up BuildingFlow using information from EtaDefns, Coords, Grids and Mets.

  Implicit None
  ! Argument list:
  Type(EtaDefns_), Intent(In)           :: EtaDefns ! Collection of eta definitions.
  Type(Coords_),   Intent(In)           :: Coords   ! Collection of coord systems.
  Type(Grids_),    Intent(In)           :: Grids    ! Collection of grids.
  Type(Mets_),     Intent(In),   Target :: Mets     ! Set of met module instance
                                                    ! states.
  Type(BuildingFlow_), Intent(InOut)    :: BuildingFlow  ! State of a building flow module
                                                    ! instance.

  ! Add grids from met module instance to BuildingFlow%C.
  Call AddGridsToCommonFlow(            &
         1, (/ 'Flow for buildings' /), &
         1, (/ 'Flow for buildings' /), &
         BuildingFlow%C                 &
       )

  ! Check first grid.
  If (BuildingFlow%C%nHGrids /= 1) Then
    Call Message('Error in SetUpCoordsEtc_BuildingFlow', 3)
  End If

  ! Check first grid.
  If (BuildingFlow%C%nZGrids /= 1) Then
    Call Message('Error in SetUpCoordsEtc_BuildingFlow', 3)
  End If

End Subroutine SetUpBuildingFlow_MetsCoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForUpdateBuildingFlow( &
             Coords, Grids, Domains,     &
             Mets,                       &
             iCase,                      &
             Time,                       &
             TValid, UpdateNow,          &
             BuildingFlow,               &
             Units                       &
           )
! Prepares for updating an instance of the building flow module.

! This routine must set TValid and UpdateNow but must not alter the validity of the building flow module
! instance.
! $$ this comment should be with generic definition.

  Implicit None
  ! Argument list:
  Type(Coords_),       Intent(In)           :: Coords
  Type(Grids_),        Intent(In)           :: Grids
  Type(Domains_),      Intent(In),   Target :: Domains
  Type(Mets_),         Intent(In)           :: Mets
  Integer,             Intent(In)           :: iCase
  Type(Time_),         Intent(In)           :: Time
  Type(Time_),         Intent(Out)          :: TValid
  Logical,             Intent(Out)          :: UpdateNow
  Type(BuildingFlow_), Intent(InOut)        :: BuildingFlow
  Type(Units_),        Intent(InOut)        :: Units
  ! Coords       :: Collection of coord systems.
  ! Grids        :: Collection of grids.
  ! Domains      :: Collection of domains.
  ! Mets         :: Set of met module instance states.
  ! iCase        :: Number of case.
  ! Time         :: Time for which the building flow module instance might be updated.
  ! TValid       :: Earliest time that the validity (overall or for any single attribute) of the flow module
  !                 instance might change, except perhaps for a change caused by a change in the validity of
  !                 the met and flow module instances acting as data sources, assuming the flow module
  !                 instance is updated now. The value is that determined at the end of this routine (the
  !                 actual time may be later).
  ! UpdateNow    :: Indicates the flow module instance must be updated now (even if update-on-demand is
  !                 specified). If set, TValid need not be set to any particular time.
  ! BuildingFlow :: State of a building flow module instance.
  ! Units        :: Collection of information on input/output unit numbers.
  ! Locals:
  Type(Domain_), Pointer :: Domain ! Abbreviation for domain.

  ! Abbreviations.
  Domain => Domains%Domains(BuildingFlow%C%iDomain)

  ! Validity due to time domain.
  If (Time < StartTimeOfDomain(Domain)) Then
    TValid = StartTimeOfDomain(Domain)
  Else If (Time >= EndTimeOfDomain(Domain)) Then
    TValid = InfFutureTime()
  Else
    TValid = EndTimeOfDomain(Domain)
  End If
  UpdateNow = .false.

End Subroutine PrepareForUpdateBuildingFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine BuildingFlowReqs(                                           &
             BuildingFlow, Time,                                       &
             WantFlowField,     WantCloudField,     WantRainField,     &
             iFlowUpdateSubset, iCloudUpdateSubset, iRainUpdateSubset, &
             FlowField,         CloudField,         RainField          &
           )
! .

  Implicit None
  ! Argument list:
  Type(BuildingFlow_), Intent(In), Target :: BuildingFlow !
  Type(Time_),         Intent(In)  :: Time
  Logical,             Intent(Out) :: WantFlowField
  Logical,             Intent(Out) :: WantCloudField
  Logical,             Intent(Out) :: WantRainField
  Integer,             Intent(Out) :: iFlowUpdateSubset ! Note this is index in
  ! BuildingFlow%C%iUpdateSubsets (to ensure subset acceptable).
  Integer,             Intent(Out) :: iCloudUpdateSubset
  Integer,             Intent(Out) :: iRainUpdateSubset
  Type(FlowField_),    Intent(Out) :: FlowField    ! FlowField type containing
                                                   ! information on what flow
                                                   ! information is required, but
                                                   ! with no flow information.
  Type(CloudField_),   Intent(Out) :: CloudField
  Type(RainField_),    Intent(Out) :: RainField

  WantFlowField           = .true.
  iFlowUpdateSubset       = 1
  FlowField%C%Valid       = .false. ! not needed, but convenient for consistency

  FlowField%C%iHGrid      = BuildingFlow%C%iHGrids(1)
  FlowField%C%iZGrid      = BuildingFlow%C%iZGrids(1)

  FlowField%C%Dt          = Char2Time('01:00', Interval = .true.)
  FlowField%C%UseTwoTimes = .false.
  FlowField%Flow        => BuildingFlow%Flow
  FlowField%ProfileData => BuildingFlow%ProfileData

  WantCloudField = .false.
  ! Avoid compiler warning.
  iCloudUpdateSubset = 0
  CloudField%C%Valid = .false.

  WantRainField = .false.
  ! Avoid compiler warning.
  iRainUpdateSubset = 0
  RainField%C%Valid = .false.

End Subroutine BuildingFlowReqs

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateBuildingFlow(                 &
             Coords, Grids, Domains,           &
             Mets,                             &
             iCase,                            &
             Time,                             &
             FlowField, CloudField, RainField, &
             BuildingFlow,                     &
             Units                             &
           )
! Updates an instance of the flow state.

  Implicit None
  ! Argument list:
  Type(Coords_),       Intent(In),   Target :: Coords       !
  Type(Grids_),        Intent(In),   Target :: Grids        ! Collection of grids.
  Type(Domains_),      Intent(In),   Target :: Domains      ! Collection of domains.
  Type(Mets_),         Intent(In)           :: Mets         !
  Integer,             Intent(In)           :: iCase        !
  Type(Time_),         Intent(In)           :: Time         ! Time.
  Type(FlowField_),    Intent(In)           :: FlowField
  Type(CloudField_),   Intent(In)           :: CloudField   !
  Type(RainField_),    Intent(In)           :: RainField    !
  Type(BuildingFlow_), Intent(InOut)        :: BuildingFlow ! Flow state to be updated.
  Type(Units_),        Intent(InOut)        :: Units        !
  ! Locals:
  Type(Domain_), Pointer :: Domain  ! Abbreviation for domain.
  Logical                :: Valid   !} Validity due to time domain.
  Type(Time_)            :: TValid  !}
  Logical                :: MValid  !] Validity due to met.
  Type(Time_)            :: MTValid !]

  ! Abbreviations.
  Domain => Domains%Domains(BuildingFlow%C%iDomain)

  ! Validity due to time domain.
  Valid = StartTimeOfDomain(Domain) <= Time .and. Time < EndTimeOfDomain(Domain)
  If (Valid) Then
    TValid = EndTimeOfDomain(Domain)
  Else
    TValid = StartTimeOfDomain(Domain)
    If (Time >= TValid) TValid = InfFutureTime()
  End If

  BuildingFlow%flow_state_building%Coords => Coords
  BuildingFlow%flow_state_building%Grids  => Grids

! Initialise met parameters
  BuildingFlow%flow_state_building%met%k = 0.4   ! von Karman constant

  If (                                                 &
    .not. FlowField%ProfileData(1,1,1,1)%MeanFlow .or. &
    .not. FlowField%ProfileData(1,1,1,1)%Turb          &
  ) Then
    Call Message('Error: Buildings currently requires analytic profile', 3)
  End If

  BuildingFlow%flow_state_building%met%u_star = FlowField%ProfileData(1,1,1,1)%UStar
  BuildingFlow%flow_state_building%met%h      = FlowField%ProfileData(1,1,1,1)%H
  BuildingFlow%flow_state_building%met%z_0    = FlowField%ProfileData(1,1,1,1)%Z0
  BuildingFlow%WindDirection                  = FlowField%ProfileData(1,1,1,1)%Phi0
  BuildingFlow%flow_state_building%FlowField  = FlowField

! Modify orientation of building

  BuildingFlow%flow_state_building%vortices%building_phi =             &
    BuildingFlow%BuildingOrientation - BuildingFlow%WindDirection - Pi

  Do While ( Abs(BuildingFlow%flow_state_building%vortices%building_phi) >= Pi/4.0)
    If ( BuildingFlow%flow_state_building%vortices%building_phi < -Pi/4.0 ) Then
      BuildingFlow%flow_state_building%vortices%building_phi =          &
        BuildingFlow%flow_state_building%vortices%building_phi + Pi/2.0
    Else
      BuildingFlow%flow_state_building%vortices%building_phi =          &
        BuildingFlow%flow_state_building%vortices%building_phi - Pi/2.0
    End If
  End Do

! Initialise dimensions of building effects region
  Call building_parameters(BuildingFlow%flow_state_building)

! Initialise constants for 'dividing streamline' and correlation time scale
  BuildingFlow%flow_state_building%constants%pi = pi
  Call calculate_constants(BuildingFlow%flow_state_building, ' ')

  Call vortex_axis_position(BuildingFlow%flow_state_building)

  ! determine axis position of vortices
  Call vortex_parameters(BuildingFlow%flow_state_building, ' ')

  ! Overall validity.
  BuildingFlow%C%Valid                     = Valid
  BuildingFlow%C%ValidAttribParams(:)      = .false.
  BuildingFlow%C%ValidAttribParams(A_Flow) = Valid
  BuildingFlow%C%ValidityUnimprovable      = Valid .or. BuildingFlow%C%FixedMet
  BuildingFlow%C%TValid                    = TValid

End Subroutine UpdateBuildingFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine BuildingFlowCoordIndices(                            &
             BuildingFlow, nHIndices, HIndex, nZIndices, ZIndex &
           )
! .

  Implicit None
  ! Argument list:
  Type(BuildingFlow_), Intent(In)  :: BuildingFlow       !
  Integer,             Intent(Out) :: nHIndices          !
  Integer,             Intent(Out) :: HIndex(MaxHCoords) !
  Integer,             Intent(Out) :: nZIndices          !
  Integer,             Intent(Out) :: ZIndex(MaxHCoords) !

  nHIndices = 1
  HIndex(1) = BuildingFlow%C%iHCoords(1)
  nZIndices = 1
  ZIndex(1) = BuildingFlow%C%iZCoords(1)

End Subroutine BuildingFlowCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine GetBuildingFlow(            &
             Coords, Grids,            &
             BuildingFlow,             &
             Moisture,                 &
             Inhomog, Homog,           &
             Time, Position,           &
             Flow, ProfileData         &
           )
! .

  Implicit None
  ! Argument list:
  Type(Coords_),       Intent(In)            :: Coords
  Type(Grids_),        Intent(In)            :: Grids
  Logical,             Intent(In)            :: Moisture     ! Indicates Q is required.
  Logical,             Intent(In)            :: Inhomog      ! Indicates inhomogeneous
                                                             ! quantities are required.
  Logical,             Intent(In)            :: Homog        ! Indicates homogeneous
                                                             ! quantities are required.
  Type(ShortTime_),    Intent(In)            :: Time         !
  Type(Position_),     Intent(In)            :: Position     !
  Type(BuildingFlow_), Intent(In)            :: BuildingFlow !
  Type(Flow_),         Intent(Out)           :: Flow         !
  Type(ProfileData_),  Intent(Out), Optional :: ProfileData
  ! Locals:
  Real(Std)        :: ZL, XL, YL
  Real(Std)        :: XRot(3)  ! Position rotated into frame aligned with wind
  Type(turbulence) :: building_turbulence

  Real(Std) :: Speed
  Real(Std) :: Theta
  Real(Std) :: enhancement(3)
  Real(Std) :: UVortices(3)
  Real(Std) :: URot(3)

  Flow%iHCoord = BuildingFlow%C%iHCoords(1)
  Flow%iZCoord = BuildingFlow%C%iZCoords(1)

  XL  = Position%XY(1, BuildingFlow%C%iHCoords(1))
  YL  = Position%XY(2, BuildingFlow%C%iHCoords(1))
  ZL  = Position%Z (   BuildingFlow%C%iZCoords(1))

  ! Translate and rotate particle position
  XRot(1) = (YL - BuildingFlow%Yb) * Sin(BuildingFlow%WindDirection) + &
            (XL - BuildingFlow%Xb) * Cos(BuildingFlow%WindDirection)
  XRot(2) = (YL - BuildingFlow%Yb) * Cos(BuildingFlow%WindDirection) - &
            (XL - BuildingFlow%Xb) * Sin(BuildingFlow%WindDirection)
  XRot(3) = ZL

  ! Call mean flow advection
  ! Fix up for particle outside effects region. $$ (could use Vortices for continuity)
  If (XRot(1) <= - BuildingFlow%flow_state_building%cuboid%lu -      &
                   BuildingFlow%flow_state_building%cuboid%lb/2 .or. &
      XRot(1) >=   BuildingFlow%flow_state_building%cuboid%lw +      &
                   BuildingFlow%flow_state_building%cuboid%lb/2 .or. &
      XRot(2) >=   BuildingFlow%flow_state_building%cuboid%we/2 .or. &
      XRot(2) <= - BuildingFlow%flow_state_building%cuboid%we/2 .or. &
      XRot(3) >=   BuildingFlow%flow_state_building%cuboid%he) Then
    URot = u_mean(XRot(3), BuildingFlow%flow_state_building, ' ')
  Else
    call building_vortices(XRot, BuildingFlow%flow_state_building, UVortices)
    ! Dt fixed at 1.0
    Call mean(XRot, URot, 1.0, BuildingFlow%flow_state_building, ' ')
    URot = URot + UVortices
  End If

  ! Rotate velocity back
  Flow%U(1) = URot(1) * Cos(BuildingFlow%WindDirection) - &
              URot(2) * Sin(BuildingFlow%WindDirection)
  Flow%U(2) = URot(1) * Sin(BuildingFlow%WindDirection) + &
              URot(2) * Cos(BuildingFlow%WindDirection)
  Flow%U(3) = URot(3)

  Call turbulence_building(XRot, BuildingFlow%flow_state_building, &
                           building_turbulence, ' ', enhancement)

  Flow%SigUU(1) = building_turbulence%sig2(1) ! rotate $$
  Flow%SigUU(2) = building_turbulence%sig2(2)
  Flow%SigUU(3) = building_turbulence%sig2(3)

!  Flow%dSigUUdX(1) = building_turbulence%dsig2_dx(1)
!  Flow%dSigUUdX(2) = building_turbulence%dsig2_dx(2)
  Flow%dSigUUdX(3) = building_turbulence%dsig2_dx(3)

  Flow%TauUU(1) = building_turbulence%tau(1)
  Flow%TauUU(2) = building_turbulence%tau(2)
  Flow%TauUU(3) = building_turbulence%tau(3)

  Flow%K(1) = Flow%SigUU(1) * Flow%TauUU(1)
  Flow%K(2) = Flow%SigUU(2) * Flow%TauUU(2)
  Flow%K(3) = Flow%SigUU(3) * Flow%TauUU(3)
  Flow%dKdX(3) = Flow%dSigUUdX(3) * Flow%TauUU(3) + &
               Flow%SigUU(3) * building_turbulence%dTauUdX(3)

  Flow%MaxDZ = Huge(Flow%MaxDZ)
 ! If (Flow%dSigUUdX(3) /= 0.0) Flow%MaxDZ = Min(Flow%MaxDZ, Abs(Flow%SigUU(3)/Flow%dSigUUdX(3)))
 ! If (Flow%dTauWDZ  /= 0.0) Flow%MaxDZ = Min(Flow%MaxDZ, Abs(Flow%TauUU(3)/Flow%dTauWDZ))
 ! Flow%MaxDZ = Min(Flow%MaxDZ, SingleSiteFlow%SingleSiteMet%ProcMet%H)
  Flow%MaxDZ = Min(Flow%MaxDZ, Flow%K(3) / Abs(Flow%dKdX(3)))

  Flow%H = BuildingFlow%flow_state_building%met%h

  !$$ need to calculate DeltaI for puffs $$
  Flow%DeltaI = Flow%H/4.0

  Flow%dUdT(:)     = 0.0 ! ! $$ Also remember must be in chronological time direction - use IsBackwards

  If (Present(ProfileData)) Then
    ProfileData%MeanFlow = .false.
    ProfileData%Turb     = .false.
    ProfileData%Canopy   = .false.
  End If

End Subroutine GetBuildingFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine BuildingReflectCoordIndices(                           &
             BuildingFlow, nHCoords, iHCoords, nZCoords, iZCoords &
           )
! Returns indices of coord systems in which positions need to be specified when using
! BuildingReflect.

  Implicit None
  ! Argument list:
  Type(BuildingFlow_), Intent(In)  :: BuildingFlow         ! State of a building flow
                                                           ! module instance.
  Integer,             Intent(Out) :: nHCoords             !} Number and values of
  Integer,             Intent(Out) :: iHCoords(MaxHCoords) !} indices of horizontal
  Integer,             Intent(Out) :: nZCoords             !} and vertical coord
  Integer,             Intent(Out) :: iZCoords(MaxHCoords) !} systems.

  ! Horizontal coords needed: none
  nHCoords    = 1
  iHCoords(1) = BuildingFlow%C%iHCoords(1)

  ! Vertical coords needed:
  ! 1) m agl.
  nZCoords    = 1
  iZCoords(1) = BuildingFlow%C%iZCoords(1)

End Subroutine BuildingReflectCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine BuildingReflect(                 &
             Coords, BuildingFlow,          &
             Vel, OldPosition, Position, U, &
             Reflected                      &
           )
! .

  Implicit None
  ! Argument list:
  Type(Coords_),       Intent(In)    :: Coords       ! Collection of coord systems.
  Type(BuildingFlow_), Intent(In)    :: BuildingFlow !
  Logical,             Intent(In)    :: Vel      !
  Type(Position_),     Intent(In)    :: OldPosition
  Type(Position_),     Intent(InOut) :: Position
  Real(Std),           Intent(InOut) :: U(3)          !
  Logical,             Intent(Out)   :: Reflected
  ! Locals:
  Real(Std) :: X(3)
  Real(Std) :: XOld(3)
  Real(Std) :: XRot(3)
  Real(Std) :: URot(3)
  Real(Std) :: XOldRot(3)
  Real(Std) :: XTemp(3)  ! Temp copy of XRot

  X    = Position2X(Coords, Position,    BuildingFlow%C%iHCoords(1), BuildingFlow%C%iZCoords(1))
     ! $$ could make X a pointer
  XOld = Position2X(Coords, OldPosition, BuildingFlow%C%iHCoords(1), BuildingFlow%C%iZCoords(1))
     ! $$ could make X a pointer

  ! Translate and rotate particle position
  XRot(1) = (X(2) - BuildingFlow%Yb) * Sin(BuildingFlow%WindDirection) + &
            (X(1) - BuildingFlow%Xb) * Cos(BuildingFlow%WindDirection)
  XRot(2) = (X(2) - BuildingFlow%Yb) * Cos(BuildingFlow%WindDirection) - &
            (X(1) - BuildingFlow%Xb) * Sin(BuildingFlow%WindDirection)
  XRot(3) = X(3)

  XOldRot(1) = (XOld(2) - BuildingFlow%Yb) * Sin(BuildingFlow%WindDirection) + &
               (XOld(1) - BuildingFlow%Xb) * Cos(BuildingFlow%WindDirection)
  XOldRot(2) = (XOld(2) - BuildingFlow%Yb) * Cos(BuildingFlow%WindDirection) - &
               (XOld(1) - BuildingFlow%Xb) * Sin(BuildingFlow%WindDirection)
  XOldRot(3) = XOld(3)

  ! Rotate velocity
  URot(1) = U(2) * Sin(BuildingFlow%WindDirection) + &
            U(1) * Cos(BuildingFlow%WindDirection)
  URot(2) = U(2) * Cos(BuildingFlow%WindDirection) - &
            U(1) * Sin(BuildingFlow%WindDirection)
  URot(3) = U(3)

  XTemp(:) = XRot(:)

  ! Call reflection subroutine
  Call reflection(BuildingFlow%flow_state_building%cuboid, XOldRot, &
                  XRot, URot, .False.)

  Reflected = Any(XTemp /= XRot)

  ! Rotate velocity back
  U(1) = URot(1) * Cos(BuildingFlow%WindDirection) - &
         URot(2) * Sin(BuildingFlow%WindDirection)
  U(2) = URot(1) * Sin(BuildingFlow%WindDirection) + &
         URot(2) * Cos(BuildingFlow%WindDirection)
  U(3) = URot(3)

! Rotate and translate position back

  X(1) = XRot(1) * Cos(BuildingFlow%WindDirection) - &
         XRot(2) * Sin(BuildingFlow%WindDirection)
  X(1) = X(1) + BuildingFlow%Xb
  X(2) = XRot(1) * Sin(BuildingFlow%WindDirection) + &
         XRot(2) * Cos(BuildingFlow%WindDirection)
  X(2) = X(2) + BuildingFlow%Yb
  X(3) = XRot(3)

  Position = X2Position(Coords, X, BuildingFlow%C%iHCoords(1), BuildingFlow%C%iZCoords(1))

End Subroutine BuildingReflect

!-------------------------------------------------------------------------------------------------------------

End Module BuildingFlowModule

! Module:  Eulerian interface

Module EulerianInterfaceModule

! This module provides information for the Eulerian solver

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowsModule, Only: Flows_, FlowMemory_,                            &
                       A_Flow, A_Cloud, A_Rain,                        &
                       ConvertToZ, GetAttrib, ResetFlowMemory,         &
                       Mets_,                                          &
                       Flow_, Cloud_, Rain_, Surface_, Soil_, Plant_,  &
                       ConvertFlow, ReverseFlow, CalcdZdZ
Use NWPMetModule
Use SemiLagrangianModule, Only: sl_init
Use SpeciesModule
Use SizeDistModule
Use ParticleModule, Only: AgentDecay

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: EulerianField_
Public :: SetUpEulerianGrid
Public :: SetUpEulerianArrays
Public :: EulerianFlowField
Public :: EulerianSedimentationVelocity
Public :: EulerianVerticalDiffusion
Public :: EulerianDeposition
Public :: EulerianAgentDecay

!-------------------------------------------------------------------------------------------------------------

Type :: EulerianField_ ! The state of an Eulerian flow field and associated grid

  Integer                :: iHCoord
  Integer                :: iZCoord
  Integer                :: iHGrid
  Integer                :: iZGrid
  Integer                :: iZGridBoundary
  Integer                :: iOld
  Integer                :: iNew
  Type(HCoord_), Pointer :: HCoord
  Type(HGrid_),  Pointer :: HGrid
  Type(ZGrid_),  Pointer :: ZGrid
  Type(ZGrid_),  Pointer :: ZGridBoundary
  Real(Std),     Pointer :: ZAboveGround(:,:,:)        => Null() ! $$ could move these initialisations
  Real(Std),     Pointer :: U(:,:,:,:)                 => Null() !    to SetUpEulerianGrid or use 
  Real(Std),     Pointer :: V(:,:,:,:)                 => Null() !    allocatable arrays instead of pointers.
  Real(Std),     Pointer :: W(:,:,:,:)                 => Null()
  Real(Std),     Pointer :: ModifiedW(:,:,:,:)         => Null()
  Real(Std),     Pointer :: Density(:,:,:,:)           => Null()
  Real(Std),     Pointer :: NewSource(:,:,:,:)         => Null()
  Real(Con),     Pointer :: Concentration(:,:,:,:,:)   => Null()
  Real(Std),     Pointer :: Diffusivity(:,:,:,:)       => Null()
  Real(Std),     Pointer :: Temperature(:,:,:,:)       => Null()
  Real(Std),     Pointer :: Humidity(:,:,:,:)          => Null()
  Real(Std),     Pointer :: dZdZ(:,:,:,:)              => Null()
  Real(Std),     Pointer :: Vd(:,:,:,:)                => Null()
  Real(Std),     Pointer :: Lambda(:,:,:,:)            => Null()
  Real(Std),     Pointer :: WSedAtGround(:,:,:,:)      => Null()
  Real(Std),     Pointer :: DryDepositionRate(:,:,:)   => Null()
  Real(Std),     Pointer :: WetDepositionRate(:,:,:)   => Null()
  Real(Std),     Pointer :: TotalDepositionRate(:,:,:) => Null()
  Logical                :: InitialiseFields

End Type EulerianField_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpEulerianGrid(EulerianField, Coords, Grids, Specieses)

  Implicit None

  ! Argument list:
  Type(EulerianField_), Intent(InOut)        :: EulerianField
  Type(Coords_),        Intent(In),   Target :: Coords
  Type(Grids_),         Intent(In),   Target :: Grids
  Type(Specieses_),     Intent(In)           :: Specieses

  ! Locals:
  Integer   :: i
  Integer   :: nX
  Integer   :: nY
  Integer   :: nZ
  Real(Std) :: dLambda
  Real(Std) :: dPhi
  Real(Std) :: Lx
  Real(Std) :: Ly
  Real(Std) :: Y0
  Real(Std) :: Tolerance

  ! Grids
  i = FindHGridIndex('Eulerian HGrid', Grids)
  EulerianField%iHGrid = i
  EulerianField%HGrid => Grids%HGrids(i)
  i = FindZGridIndex('Eulerian ZGrid', Grids)
  EulerianField%iZGrid = i
  EulerianField%ZGrid => Grids%ZGrids(i)
  i = FindZGridIndex('Eulerian ZGridBoundary', Grids)
  EulerianField%iZGridBoundary = i
  EulerianField%ZGridBoundary => Grids%ZGrids(i)

  If (EulerianField%ZGridBoundary%nZ /= EulerianField%ZGrid%nZ) Then
    Call Message('ERROR in SetUpEulerianField: number of Eulerian grid boundary levels ' // &
                 'does not agree with the main Eulerian grid levels',3)
  End If

  nX = EulerianField%HGrid%nX
  nY = EulerianField%HGrid%nY
  nZ = EulerianField%ZGrid%nZ

  i = FindZCoordIndex(EulerianField%ZGrid%ZCoordName, Coords)
  EulerianField%iZCoord = i

  ! Horizontal coordinates
  i = FindHCoordIndex(EulerianField%HGrid%HCoordName, Coords)
  EulerianField%iHCoord = i
  EulerianField%HCoord => Coords%HCoords(i)

  ! Horizontal grid spacing
  dLambda = EulerianField%HGrid%dX * EulerianField%HCoord%Unit(1)
  dPhi    = EulerianField%HGrid%dY * EulerianField%HCoord%Unit(2)
  
  !$$ Note Hcoord must be lat-long based and Zcoord must be height-based eta or possibly height above ground.

  ! Limited area or global domain?
  ! Global: EulerianField%HGrid%Wrap = .True.
  ! Local:  EulerianField%HGrid%Wrap = .False.

  ! Horizontal domain size
  If ( EulerianField%HGrid%Wrap ) Then
    Lx = dLambda * Real(nX)
  Else
    Lx = dLambda * Real(nX - 1)
  End If
  Y0 = EulerianField%HGrid%Y0 * EulerianField%HCoord%Unit(2) + EulerianField%HCoord%Origin(2)
  Ly = dPhi * Real(nY - 1)

  Tolerance = 1.e-3
  If ( EulerianField%HGrid%Wrap .And. (Abs(Lx - 2.0 * Pi) > Tolerance .Or. Abs(Ly - Pi) > Tolerance) ) Then
    Call Message('FATAL ERROR: if a "wrapped" grid is used for Eulerian advection, it must be global.', 3)
  Else If (                                                                &
            .Not.EulerianField%HGrid%Wrap .And.                            &
            (Abs(Lx - 2.0 * Pi) < Tolerance .Or. Abs(Ly - Pi) < Tolerance) &
           ) Then
    Call Message(                                                                     &
           'FATAL ERROR: if a non-"wrapped" grid is used for Eulerian advection, ' // &
           'it must be limited in longitude and latitude.',                           &
           3                                                                          &
         )
  End If
  
  ! Initialise semi-Lagrangian advection scheme
  If (Specieses%AdvectedFields) Call sl_init(nX,                         &
                                             nY,                         &
                                             nZ,                         &
                                             EulerianField%ZGrid%Z,      &
                                             Lx,                         &
                                             Ly,                         &
                                             Y0,                         &
                                             dLambda,                    &
                                             dPhi,                       &
                                             EulerianField%HGrid%Wrap)

  EulerianField%InitialiseFields = .True.

End Subroutine SetUpEulerianGrid

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpEulerianArrays(EulerianField, Specieses, Diffusion, WetDep, AgentDecay)

  Implicit None

  ! Argument list:
  Type(EulerianField_), Intent(InOut) :: EulerianField
  Type(Specieses_),     Intent(In)    :: Specieses
  Logical,              Intent(In)    :: Diffusion
  Logical,              Intent(In)    :: WetDep
  Logical,              Intent(In)    :: AgentDecay

  ! Locals:
  Integer   :: nTracers
  Integer   :: nX
  Integer   :: nY
  Integer   :: nZ

  nX = EulerianField%HGrid%nX
  nY = EulerianField%HGrid%nY
  nZ = EulerianField%ZGrid%nZ

  ! Number of tracers
  nTracers = Specieses%nFields

  ! Allocate concentration field
  If (nTracers > 0) Then
    Allocate(EulerianField%Concentration(-1:nX + 2, -1:nY + 2, nZ, nTracers, 2))
    EulerianField%Concentration(:,:,:,:,:) = 0.0
  End If
  EulerianField%iNew = 1
  EulerianField%iOld = 2

  If (Specieses%AdvectedFields) Then ! $$ currently non advected fields have no dep, sedimentation etc
                                     !    (but could have separate flags).

    ! Allocate velocity field
    Allocate(EulerianField%Density(-1:nX + 2, -1:nY + 2, EulerianField%ZGrid%nZ, 2))
    EulerianField%Density(:,:,:,:) = 0.0

    Allocate(EulerianField%U(-1:nX + 2, -1:nY + 2, nZ, 2))
    EulerianField%U(:,:,:,:) = 0.0

    Allocate(EulerianField%V(-1:nX + 2, -1:nY + 2, nZ, 2))
    EulerianField%V(:,:,:,:) = 0.0

    Allocate(EulerianField%W(-1:nX + 2, -1:nY + 2, nZ, 2))
    EulerianField%W(:,:,:,:) = 0.0

    ! Allocate Eulerian source
    Allocate(EulerianField%NewSource(-1:nX + 2, -1:nY + 2, nZ, nTracers))
    EulerianField%NewSource(:,:,:,:) = 0.0

    ! Allocate total deposition rate
    Allocate(EulerianField%TotalDepositionRate(nX, nY, nTracers))
    EulerianField%TotalDepositionRate(:,:,:) = 0.0

    ! Allocate diffusivity and 'height-above-ground' variable
    If ( Diffusion ) Then
      Allocate(EulerianField%Diffusivity(nX, nY, nZ, 2))
      EulerianField%Diffusivity(:,:,:,:) = 0.0
      Allocate(EulerianField%ZAboveGround(nX, nY, nZ))
      EulerianField%ZAboveGround(:,:,:) = 0.0
      Allocate(EulerianField%Vd(nX, nY, nTracers, 2))
      EulerianField%Vd(:,:,:,:) = 0.0
      Allocate(EulerianField%DryDepositionRate(nX, nY, nTracers))
      EulerianField%DryDepositionRate(:,:,:) = 0.0
    End If

    If ( WetDep ) Then
      Allocate(EulerianField%Lambda(nX, nY, nZ, nTracers))
      EulerianField%Lambda(:,:,:,:) = 0.0
      Allocate(EulerianField%WetDepositionRate(nX, nY, nTracers))
      EulerianField%WetDepositionRate(:,:,:) = 0.0
    End If

    ! Allocate humidity for agent decay
    If ( AgentDecay ) Allocate(EulerianField%Humidity(nX, nY, nZ, 2))

    ! Allocate temperature if sedimenting particles present or agent decay
    If ( Specieses%AdvectedSedimentingFields .Or. AgentDecay ) Then
      Allocate(EulerianField%Temperature(nX, nY, nZ, 2))
    End If

    ! Allocate appropriate fields if sedimenting particles present
    If ( Specieses%AdvectedSedimentingFields ) Then
      Allocate(EulerianField%dZdZ(nX, nY, nZ, 2))
      Allocate(EulerianField%ModifiedW(-1:nX + 2, -1:nY + 2, nZ, 2))
      Allocate(EulerianField%WSedAtGround(nX, nY, nTracers, 2))
      EulerianField%WSedAtGround(:,:,:,:) = 0.0
    End If
  
  End If

End Subroutine SetUpEulerianArrays

!-------------------------------------------------------------------------------------------------------------

Subroutine EulerianFlowField(EulerianField, Specieses, SizeDists, Time, Coords, Grids, Domains, &
                             Mets, Flows, FlowMemory, Units, Diffusion, DryDep, WetDep, AgentDecay)

  Implicit None

  ! ArgumentList:
  Type(EulerianField_), Intent(InOut)      :: EulerianField
  Type(Specieses_),     Intent(In)         :: Specieses
  Type(SizeDists_),     Intent(In), Target :: SizeDists
  Type(Time_),          Intent(In)         :: Time
  Type(Coords_),        Intent(In)         :: Coords
  Type(Grids_),         Intent(In)         :: Grids
  Type(Domains_),       Intent(In)         :: Domains
  Type(Mets_),          Intent(InOut)      :: Mets
  Type(Flows_),         Intent(InOut)      :: Flows
  Type(FlowMemory_),    Intent(InOut)      :: FlowMemory
  Type(Units_),         Intent(InOut)      :: Units
  Logical,              Intent(In)         :: Diffusion
  Logical,              Intent(In)         :: DryDep
  Logical,              Intent(In)         :: WetDep
  Logical,              Intent(In)         :: AgentDecay

  ! Locals:
  Type(SizeDist_), Pointer :: SizeDist
  Type(Position_)          :: Position
  Type(Flow_)              :: Flow
  Type(Cloud_)             :: Cloud
  Type(Rain_)              :: Rain
  Type(Plant_)             :: Plant
  Type(Soil_)              :: Soil
  Type(Surface_)           :: Surface
  Integer                  :: i
  Integer                  :: j
  Integer                  :: k
  Integer                  :: nX
  Integer                  :: nY
  Integer                  :: nZ
  Integer                  :: iNew
  Integer                  :: iTracer
  Integer                  :: iSpecies
  Integer                  :: iSpeciesUse
  Integer                  :: ErrorCode
  Real(Std)                :: X(3)
  Real(Std)                :: TempVar(3)
  Real(Std)                :: dZdZ
  Real(Std)                :: Diameter

  iNew = EulerianField%iNew
  
  nX = EulerianField%HGrid%nX
  nY = EulerianField%HGrid%nY
  nZ = EulerianField%ZGrid%nZ
  
  ! Interpolate meteorological fields to target time
  Do i = 1, nX
    Do j = 1, nY
      Do k = 1, nZ

        Call ResetFlowMemory(Flows, FlowMemory)

        X = (/ EulerianField%HGrid%X(i), EulerianField%HGrid%Y(j), EulerianField%ZGrid%Z(k) /)
        Position = X2Position(Coords, X, EulerianField%HGrid%iHCoord, EulerianField%ZGrid%iZCoord)

        ! Get flow info.
        Call GetAttrib(                              &
               iAttribParam  = A_Flow,               &
               Coords        = Coords,               &
               Grids         = Grids,                &
               Domains       = Domains,              &
               Moisture      = .false.,              &
               Inhomog       = .true.,               &
               Homog         = .false.,              &
               Time          = Time2ShortTime(Time), &
               AnyTravelTime = .true.,               &
               TravelTime    = ZeroShortTime(),      &
               Position      = Position,             &
               Units         = Units,                &
               Mets          = Mets,                 &
               Flows         = Flows,                &
               FlowMemory    = FlowMemory,           &
               Flow          = Flow,                 &
               Cloud         = Cloud,                &
               Rain          = Rain,                 &
               Surface       = Surface,              &
               Soil          = Soil,                 &
               Plant         = Plant,                &
               ErrorCode     = ErrorCode             &
             )
        If (ErrorCode > 0) Then
          Write (*,*) 'Error in getting flow velocity', ErrorCode
          Stop
       End If
       !If (IsBackwards()) Call ReverseFlow(Flow)  !$ Allow backward dispersion in Eulerian scheme

       If ( EulerianField%InitialiseFields ) Then
         Call ConvertToZ(                                                &
                Coords, Grids, Domains,                                  &
                Flow%iZCoord,                                            &  ! Assumes Flow%iZCoord is
                Time2ShortTime(Time), .true., ZeroShortTime(), Position, &  ! height-above-ground index
                Units, Mets, Flows,                                      &
                FlowMemory, ErrorCode                                    &
              )
       End If

       EulerianField%Density(i, j, k, iNew) = Flow%Rho
       EulerianField%U(i, j, k, iNew)       = Flow%U(1)/EarthRadius   ! Scaled by EarthRadius here:
       EulerianField%V(i, j, k, iNew)       = Flow%U(2)/EarthRadius   ! required by SL advection scheme
       EulerianField%W(i, j, k, iNew)       = Flow%EtaDot

       If ( Diffusion ) Then
         EulerianField%Diffusivity(i, j, k, iNew) = Flow%K(3)
         If ( EulerianField%InitialiseFields ) Then
           TempVar(:) = Position2X(Coords, Position, Flow%iHCoord, Flow%iZCoord)
           EulerianField%ZAboveGround(i, j, k) = TempVar(3)
         Else If ( DryDep .And. k == 1 ) Then
           Do iTracer = 1, Specieses%nFields
             iSpecies = Specieses%iField2Species(iTracer)
             ! $$ Remove once land use scheme working for Eulerian model
             If (Specieses%Specieses(iSpecies)%LandUseDryDep) &
               Call Message('FATAL ERROR: land use scheme not yet working for Eulerian advection.', 3)
             If (Specieses%iField2SizeDist(iTracer) == 0) Then
               Diameter = 0.0
             Else
               iSpeciesUse = Specieses%iSpecies2SpeciesUses(iSpecies)
               SizeDist => SizeDists%SizeDists(Specieses%SpeciesUseses(iSpeciesUse)%iSizeDist)
               Diameter = SizeDist%RepresentativeDiameter(Specieses%iField2Size(iTracer))
             End If
               EulerianField%Vd(i, j, iTracer, iNew) = CalcVd(Specieses%Specieses(iSpecies),            &
                                                              Flow,                                     &
                                                              Surface,                                  &
                                                              Plant,                                    &
                                                              Time,                                     &
                                                              Cloud3D = 0.0,                            &
                                                              WSed = 0.0,                               &
                                                              Diameter = Diameter,                      &
                                                              Z = EulerianField%ZAboveGround(i, j, 1),  &
                                                              Zs = EulerianField%ZAboveGround(i, j, 1), &
                                                              YLatLong = 0.0                            &
                                                       )
                                                       ! WSed : sedimentation not treated via Vd on fields
                                                       ! $$ YLatLong - Will need to calculate this properly
                                                       !    if using LandUseDryDep
           End Do
         End If
       End If

       If ( .Not.EulerianField%InitialiseFields .And. WetDep ) Then
         ! Get cloud info.
         Call GetAttrib(                              &
                iAttribParam  = A_Cloud,              &
                Coords        = Coords,               &
                Grids         = Grids,                &
                Domains       = Domains,              &
                Moisture      = .false.,              &
                Inhomog       = .true.,               &
                Homog         = .false.,              &
                Time          = Time2ShortTime(Time), &
                AnyTravelTime = .true.,               &
                TravelTime    = ZeroShortTime(),      &
                Position      = Position,             &
                Units         = Units,                &
                Mets          = Mets,                 &
                Flows         = Flows,                &
                FlowMemory    = FlowMemory,           &
                Flow          = Flow,                 &
                Cloud         = Cloud,                &
                Rain          = Rain,                 &
                Surface       = Surface,              &
                Soil          = Soil,                 &
                Plant         = Plant,                &
                ErrorCode     = ErrorCode             &
              )
         If (ErrorCode > 0) Then
           Write (*,*) 'Error in getting cloud information', ErrorCode
           Stop
         End If
         ! Get rain info.
         Call GetAttrib(                              &
                iAttribParam  = A_Rain,               &
                Coords        = Coords,               &
                Grids         = Grids,                &
                Domains       = Domains,              &
                Moisture      = .false.,              &
                Inhomog       = .true.,               &
                Homog         = .false.,              &
                Time          = Time2ShortTime(Time), &
                AnyTravelTime = .true.,               &
                TravelTime    = ZeroShortTime(),      &
                Position      = Position,             &
                Units         = Units,                &
                Mets          = Mets,                 &
                Flows         = Flows,                &
                FlowMemory    = FlowMemory,           &
                Flow          = Flow,                 &
                Cloud         = Cloud,                &
                Rain          = Rain,                 &
                Surface       = Surface,              &
                Soil          = Soil,                 &
                Plant         = Plant,                &
                ErrorCode     = ErrorCode             &
              )
         If (ErrorCode > 0) Then
           Write (*,*) 'Error in getting rain information', ErrorCode
           Stop
         End If
         Do iTracer = 1, Specieses%nFields
           iSpecies = Specieses%iField2Species(iTracer)
           iSpeciesUse = Specieses%iSpecies2SpeciesUses(iSpecies)
           ! If Lambda depends on diameter need something here like for dry dep $$
           EulerianField%Lambda(i, j, k, iTracer) = CalcWetScavCoeff(                       &
                                                      Specieses%Specieses(iSpecies),        &
                                                      Cloud,                                &
                                                      Rain,                                 &
                                                      EulerianField%ZAboveGround(i, j, k),  &
                                                      Flow%T                                &
                                                    )
         End Do
       End If

       If ( AgentDecay ) EulerianField%Humidity(i, j, k, iNew) = Flow%Q

       If ( Specieses%AdvectedSedimentingFields .Or. AgentDecay) Then 
         EulerianField%Temperature(i, j, k, iNew) = Flow%T
       End If

       If ( Specieses%AdvectedSedimentingFields ) Then 
         Call CalcdZdZ(                                &
                Coords        = Coords,                &
                Grids         = Grids,                 &
                Domains       = Domains,               &
                iZCoord1      = EulerianField%iZcoord, &
                iZCoord2      = Flow%iZCoord,          &
                Time          = Time2ShortTime(Time),  &
                AnyTravelTime = .True. ,               &
                TravelTime    = ZeroShortTime(),       &
                Position      = Position,              &
                Units         = Units,                 &
                Mets          = Mets,                  &
                Flows         = Flows,                 &
                FlowMemory    = FlowMemory,            &
                dZdZ          = dZdZ,                  &
                ErrorCode     = ErrorCode              &
              )
         EulerianField%dZdZ(i, j, k, iNew) = dZdZ
       End If

      End Do
    End Do
  End Do

  ! Fill halo
  If ( EulerianField%HGrid%Wrap ) Then
    Call FillGlobalHalo(EulerianField%U(:,:,:,iNew))
    Call FillGlobalHalo(EulerianField%V(:,:,:,iNew))
    Call FillGlobalHalo(EulerianField%W(:,:,:,iNew))
    Call FillGlobalHalo(EulerianField%Density(:,:,:,iNew))
    !Call FillGlobalHalo(EulerianField%NewSource(:,:,:,iNew))
    ! $$ currently not used. If changes, need to write to restart file too.
    ! $$ Note global halo for concentration handled in SemiLagrangian.F90
    !    - would be good to treat local and global haloes the same way.
  Else
    Call FillLocalHalo(EulerianField%U(:,:,:,iNew))
    Call FillLocalHalo(EulerianField%V(:,:,:,iNew))
    Call FillLocalHalo(EulerianField%W(:,:,:,iNew))
    Call FillLocalHalo(EulerianField%Density(:,:,:,iNew))
    Call LocalAreaBoundaryConditions(EulerianField%Concentration(:,:,:,:,iNew))
  End If

  EulerianField%InitialiseFields = .False.

End Subroutine EulerianFlowField

!-------------------------------------------------------------------------------------------------------------

Subroutine FillGlobalHalo(Field)

  Implicit None

  ! Argument list:
  Real(Std), Intent(InOut) :: Field(-1:,-1:,:)

  ! Locals
  Integer :: nX
  Integer :: nY

  nX = Size(Field, 1)
  nX = nX - 4
  nY = Size(Field, 2)
  nY = nY - 4

  Field(  1:nX,      0,:) = Field(  1:nX,      1,:)
  Field(  1:nX,     -1,:) = Field(  1:nX,      2,:)
  Field(  1:nX, nY + 1,:) = Field(  1:nX,     nY,:)
  Field(  1:nX, nY + 2,:) = Field(  1:nX, nY - 1,:)
  Field(     0,      :,:) = Field(    nX,      :,:)
  Field(    -1,      :,:) = Field(nX - 1,      :,:)
  Field(nX + 1,      :,:) = Field(     1,      :,:)
  Field(nX + 2,      :,:) = Field(     2,      :,:)

End Subroutine FillGlobalHalo

!-------------------------------------------------------------------------------------------------------------

Subroutine FillGlobalHaloConc(Field)

  Implicit None

  ! Argument list:
  Real(Con), Intent(InOut) :: Field(-1:,-1:,:)

  ! Locals
  Integer :: nX
  Integer :: nY

  nX = Size(Field, 1)
  nX = nX - 4
  nY = Size(Field, 2)
  nY = nY - 4

  Field(  1:nX,      0,:) = Field(  1:nX,      1,:)
  Field(  1:nX,     -1,:) = Field(  1:nX,      2,:)
  Field(  1:nX, nY + 1,:) = Field(  1:nX,     nY,:)
  Field(  1:nX, nY + 2,:) = Field(  1:nX, nY - 1,:)
  Field(     0,      :,:) = Field(    nX,      :,:)
  Field(    -1,      :,:) = Field(nX - 1,      :,:)
  Field(nX + 1,      :,:) = Field(     1,      :,:)
  Field(nX + 2,      :,:) = Field(     2,      :,:)

End Subroutine FillGlobalHaloConc

!-------------------------------------------------------------------------------------------------------------

Subroutine FillLocalHalo(Field)

  Implicit None

  ! Argument list:
  Real(Std), Intent(InOut) :: Field(-1:,-1:,:)

  ! Locals
  Integer :: nX
  Integer :: nY

  nX = Size(Field, 1)
  nX = nX - 4
  nY = Size(Field, 2)
  nY = nY - 4

  ! Reflect field at boundary 
  ! (strictly line in between boundary and first halo point)
  ! $$ may want to revisit this $$
  Field(  1:nX,      0,:) = Field(   1:nX,      1, :)
  Field(  1:nX,     -1,:) = Field(   1:nX,      2, :)
  Field(  1:nX, nY + 1,:) = Field(   1:nX,     nY, :)
  Field(  1:nX, nY + 2,:) = Field(   1:nX, nY - 1, :)
  Field(     0,      :,:) = Field(      1,      :, :)
  Field(    -1,      :,:) = Field(      2,      :, :)
  Field(nX + 1,      :,:) = Field(     nX,      :, :)
  Field(nX + 2,      :,:) = Field( nX - 1,      :, :)

End Subroutine FillLocalHalo

!-------------------------------------------------------------------------------------------------------------

Subroutine LocalAreaBoundaryConditions(Concentration)

  Implicit None

  ! Argument list:
  Real(Con), Intent(InOut) :: Concentration(-1:, -1:, :, :)

  ! Locals:
  Integer :: nX
  Integer :: nY

  nX = Size(Concentration, 1)
  nX = nX - 4
  nY = Size(Concentration, 2)
  nY = nY - 4

  ! Absorbing boundary condition  
  Concentration(     -1:1, -1:nY + 2, :, :) = 0.0
  Concentration(nX:nX + 2, -1:nY + 2, :, :) = 0.0
  Concentration(        :,      -1:1, :, :) = 0.0
  Concentration(        :, nY:nY + 2, :, :) = 0.0

  !$$ zero beyond (rather than on) boundary? (so outflow gives realistic boundary value; but inflow may then
  !   be a problem?)
  !$$ need to consider how it interacts with sl_advect and depart_pnt and what this does when back traj
  !   is outside domain.
  !$$ Constant flux boundary condition ???

End Subroutine LocalAreaBoundaryConditions

!-------------------------------------------------------------------------------------------------------------

Subroutine EulerianSedimentationVelocity(iTracer, Density, Diameter, EulerianField)

  Implicit None
 
  ! Argument list:
  Integer,              Intent(In)    :: iTracer
  Real(Std),            Intent(In)    :: Density       ! Particle density.
  Real(Std),            Intent(In)    :: Diameter      ! Particle diameter.
  Type(EulerianField_), Intent(InOut) :: EulerianField
  ! $$ Particle density: Could precalculate in SetUpEulerianField and store in EulerianField?
 
  ! Locals: 
  Integer   :: ix
  Integer   :: iy
  Integer   :: iz
  Integer   :: iOld
  Integer   :: iNew
  Real(Std) :: Pressure
  Real(Std) :: WSed
  Real(Std) :: ParticleShape
  Integer   :: ShapeSchemeCode

  ParticleShape = 1.0                 ! $$ in due course the Eulerian model should allow non-spherical shapes
  ShapeSchemeCode = ShapeScheme_White ! 

  iOld = EulerianField%iOld
  iNew = EulerianField%iNew

  Do iz = 1, EulerianField%ZGrid%nZ
    Do iy = 1, EulerianField%HGrid%nY
      Do ix = 1, EulerianField%HGrid%nX
        
        Pressure = EulerianField%Density(ix, iy, iz, iOld) * GasConstant  &
                   * EulerianField%Temperature(ix, iy, iz, iOld)
        
        WSed = TerminalVelocity(Diameter,                                     &
                                Density,                                      &
                                ParticleShape,                                &
                                ShapeSchemeCode,                              &
                                EulerianField%Temperature(ix, iy, iz, iOld),  &
                                Pressure,                                     &
                                EulerianField%Density(ix, iy, iz, iOld))
        
        If ( iz == 1 ) EulerianField%WSedAtGround(ix, iy, iTracer, iOld) = WSed
        
        EulerianField%ModifiedW(ix, iy, iz, iOld) =                                       &
          EulerianField%W(ix, iy, iz, iOld) - WSed * EulerianField%dZdZ(ix, iy, iz, iOld)
        
        Pressure = EulerianField%Density(ix, iy, iz, iNew) * GasConstant *    &
                     EulerianField%Temperature(ix, iy, iz, iNew)
        
        WSed = TerminalVelocity(Diameter,                                     &
                                Density,                                      &
                                ParticleShape,                                &
                                ShapeSchemeCode,                              &
                                EulerianField%Temperature(ix, iy, iz, iNew),  &
                                Pressure,                                     &
                                EulerianField%Density(ix, iy, iz, iNew))
        
        If ( iz == 1 ) EulerianField%WSedAtGround(ix, iy, iTracer, iNew) = WSed
        
        EulerianField%ModifiedW(ix, iy, iz, iNew) =                                       &
          EulerianField%W(ix, iy, iz, iNew) - WSed * EulerianField%dZdZ(ix, iy, iz, iNew)
        
      End Do
    End Do
  End Do
  If ( EulerianField%HGrid%Wrap ) Then
    Call FillGlobalHalo(EulerianField%ModifiedW(:,:,:,iOld))
    Call FillGlobalHalo(EulerianField%ModifiedW(:,:,:,iNew))
  Else
    Call FillLocalHalo(EulerianField%ModifiedW(:,:,:,iOld))
    Call FillLocalHalo(EulerianField%ModifiedW(:,:,:,iNew))
  End If

End Subroutine EulerianSedimentationVelocity

!-------------------------------------------------------------------------------------------------------------

Subroutine EulerianVerticalDiffusion(OldConcentration, NewConcentration, OldK, NewK, OldDensity, NewDensity, &
                                     TimeStep, nTracers, nX, nY, nZ, Z, OldVd, NewVd)
 
! $$ need to account for ambient density?

  Implicit None

  ! Argument list:
  Real(Con),      Intent(In)  :: OldConcentration(-1:, -1:, :, :)
  Real(Con),      Intent(Out) :: NewConcentration(-1:, -1:, :, :)
  Real(Std),      Intent(In)  :: OldK(:,:,:)
  Real(Std),      Intent(In)  :: NewK(:,:,:)
  Real(Std),      Intent(In)  :: OldDensity(-1:, -1:, :)
  Real(Std),      Intent(In)  :: NewDensity(-1:, -1:, :)
  Real(Std),      Intent(In)  :: TimeStep
  Integer,        Intent(In)  :: nTracers
  Integer,        Intent(In)  :: nX
  Integer,        Intent(In)  :: nY
  Integer,        Intent(In)  :: nZ
  Real(Std),      Intent(In)  :: Z(:,:,:)
  Real(Std),      Intent(In)  :: OldVd(:,:,:)
  Real(Std),      Intent(In)  :: NewVd(:,:,:)

  ! Locals:
  Integer   :: ix
  Integer   :: iy
  Integer   :: iz
  Integer   :: iTracer
  Real(Std) :: DeltaZLower
  Real(Std) :: DeltaZUpper
  Real(Std) :: DeltaZMidPoint
  Real(Std) :: OldKLower
  Real(Std) :: OldKUpper
  Real(Std) :: NewKLower
  Real(Std) :: NewKUpper
  Real(Std) :: LDiagonal
  Real(Std) :: LOffDiagonal
  Real(Std) :: U(nZ - 1)
  Real(Std) :: B
  Real(Std) :: OldConcentrationLower
  Real(Std) :: OldConcentrationUpper
  Real(Std) :: Z0
  Real(Std) :: Y(nZ)
  Real(Std) :: ZnExtra

  ! Crank-Nicolson using Crout factorisation to solve tridiagonal system A x = B
  Do iTracer = 1, nTracers
    Do iy = 1, nY
      Do ix = 1, nX
        
        ! Lower boundary condition
        Z0 = 2.0 * Z(ix, iy, 1) - Z(ix, iy, 2)     ! Level below ground

        ! Solve L y = B
        DeltaZUpper = Z(ix, iy, 2) - Z(ix, iy, 1)
        DeltaZMidPoint = 0.5 * (Z(ix, iy, 2) - Z0)

        OldKUpper = 0.5 * (OldK(ix, iy, 1) * OldDensity(ix, iy, 1) + OldK(ix, iy, 2) * OldDensity(ix, iy, 2))
        NewKUpper = 0.5 * (NewK(ix, iy, 1) * NewDensity(ix, iy, 1) + NewK(ix, iy, 2) * NewDensity(ix, iy, 2))

        LDiagonal = 1.0 + 0.5 * TimeStep * (NewKUpper/(DeltaZUpper * NewDensity(ix, iy, 1))  &
                      + NewVd(ix, iy, iTracer))/DeltaZMidPoint
        U(1) = -0.5 * TimeStep * NewKUpper/(DeltaZMidPoint * DeltaZUpper * NewDensity(ix, iy, 2))
        U(1) = U(1)/LDiagonal

        OldConcentrationUpper = OldConcentration(ix, iy, 2, iTracer)
        B = (1.0 - 0.5 * TimeStep * (OldKUpper/(DeltaZUpper * OldDensity(ix, iy, 1)) +       &
              OldVd(ix, iy, iTracer))/DeltaZMidPoint) * OldConcentration(ix, iy, 1, iTracer) &
            + 0.5 * TimeStep * OldKUpper *                                                   &
              OldConcentrationUpper/(DeltaZUpper * DeltaZMidPoint * OldDensity(ix, iy, 2))
        Y(1) = B/LDiagonal

        ! Intermediate values
        Do iz = 2, nZ - 1

          DeltaZLower = Z(ix, iy, iz) - Z(ix, iy, iz - 1)
          DeltaZUpper = Z(ix, iy, iz + 1) - Z(ix, iy, iz)
          DeltaZMidPoint = 0.5 * (Z(ix, iy, iz + 1) - Z(ix, iy, iz - 1))
          OldKLower  = 0.5 * (OldK(ix, iy, iz) * OldDensity(ix, iy, iz)    &
                       + OldK(ix, iy, iz - 1) * OldDensity(ix, iy, iz - 1))
          OldKUpper  = 0.5 * (OldK(ix, iy, iz) * OldDensity(ix, iy, iz)    &
                       + OldK(ix, iy, iz + 1) * OldDensity(ix, iy, iz + 1))
          NewKLower  = 0.5 * (NewK(ix, iy, iz) * NewDensity(ix, iy, iz)    &
                       + NewK(ix, iy, iz - 1) * NewDensity(ix, iy, iz - 1))
          NewKUpper  = 0.5 * (NewK(ix, iy, iz) * NewDensity(ix, iy, iz)    &
                       + NewK(ix, iy, iz + 1) * NewDensity(ix, iy, iz + 1))

          LOffDiagonal = -0.5 * TimeStep * NewKLower/(DeltaZMidPoint * DeltaZLower * NewDensity(ix, iy, iz - 1))
          LDiagonal    =  1.0 + 0.5 * TimeStep * (NewKUpper/DeltaZUpper + NewKLower/DeltaZLower)  &
                          /(DeltaZMidPoint * NewDensity(ix, iy, iz))                              &
                          - LOffDiagonal * U(iz - 1)

          U(iz) = -0.5 * TimeStep * NewKUpper/(DeltaZMidPoint * DeltaZUpper * NewDensity(ix, iy, iz + 1))
          U(iz) = U(iz)/LDiagonal

          OldConcentrationLower = OldConcentration(ix, iy, iz - 1, iTracer)
          OldConcentrationUpper = OldConcentration(ix, iy, iz + 1, iTracer)
          B = (1.0 - 0.5 * TimeStep * (OldKUpper/DeltaZUpper + OldKLower/DeltaZLower)                          &
              /(DeltaZMidPoint * OldDensity(ix, iy, iz))) * OldConcentration(ix, iy, iz, iTracer)              &
              + 0.5 * TimeStep * (OldKUpper * OldConcentrationUpper/(DeltaZUpper * OldDensity(ix, iy, iz + 1)) &
              + OldKLower * OldConcentrationLower/(DeltaZLower * OldDensity(ix, iy, iz - 1)) )/DeltaZMidPoint
          Y(iz) = (B - LOffDiagonal * Y(iz - 1))/LDiagonal

        End Do

        ! Upper boundary condition: zero flux b.c.
        ZnExtra = 2.0 * Z(ix, iy, nZ) - Z(ix, iy, nZ - 1)       ! Level above top of domain
        DeltaZLower = Z(ix, iy, nZ) - Z(ix, iy, nZ - 1)
        DeltaZUpper = ZnExtra - Z(ix, iy, nZ)
        DeltaZMidPoint = 0.5 * (ZnExtra - Z(ix, iy, nZ - 1))

        OldKLower  = 0.5 * (OldK(ix, iy, iz) * OldDensity(ix, iy, nZ)    &
                     + OldK(ix, iy, nZ - 1) * OldDensity(ix, iy, nZ - 1))
        NewKLower  = 0.5 * (NewK(ix, iy, iz) * NewDensity(ix, iy, nZ)    &
                     + NewK(ix, iy, nZ - 1) * NewDensity(ix, iy, nZ - 1))

        LOffDiagonal = -0.5 * TimeStep * NewKLower/(DeltaZMidPoint * DeltaZLower * NewDensity(ix, iy, nZ - 1))
        LDiagonal = 1.0 + 0.5 * TimeStep * NewKLower/(DeltaZLower * DeltaZMidPoint * NewDensity(ix, iy, nZ))  &
                    - LOffDiagonal * U(nZ - 1)

        OldConcentrationLower = OldConcentration(ix, iy, nZ - 1, iTracer)
        B = (1.0 - 0.5 * TimeStep * OldKLower/(DeltaZLower * DeltaZMidPoint * OldDensity(ix, iy, nZ))) &
            * OldConcentration(ix, iy, nZ, iTracer)                                                    &
            + 0.5 * TimeStep * OldKLower                                                               &
            * OldConcentrationLower/(DeltaZLower * DeltaZMidPoint * OldDensity(ix, iy, nZ - 1))
        Y(nZ) = (B - LOffDiagonal * Y(nZ - 1))/LDiagonal

        ! Solve Ux = z  (x = NewConcentration)
        NewConcentration(ix, iy, nZ, iTracer) = Y(nZ)
        Do iz = nZ - 1, 1, -1
          NewConcentration(ix, iy, iz, iTracer) = Y(iz) - U(iz) * NewConcentration(ix, iy, iz + 1, iTracer)
        End Do
     
      End Do
    End Do
  End Do 

End Subroutine EulerianVerticalDiffusion

!-------------------------------------------------------------------------------------------------------------

Subroutine EulerianDeposition(TimeStep, iTracer, DryDep, WetDep, Diffusion, EulerianField)

  Implicit None

  ! Argument list:
  Real(Std),            Intent(In)    :: TimeStep
  Integer,              Intent(In)    :: iTracer
  Logical,              Intent(In)    :: DryDep
  Logical,              Intent(In)    :: WetDep
  Logical,              Intent(In)    :: Diffusion
  Type(EulerianField_), Intent(InOut) :: EulerianField

  ! Locals: 
  Integer :: ix
  Integer :: iy
  Integer :: iz
  Integer :: nX
  Integer :: nY
  Integer :: nZ
  Integer :: iOld
  Integer :: iNew

  iOld = EulerianField%iOld
  iNew = EulerianField%iNew
  nX = EulerianField%HGrid%nX
  nY = EulerianField%HGrid%nY
  nZ = EulerianField%ZGrid%nZ

  ! Dry deposition rate 
  If ( DryDep .And. Diffusion ) Then
    EulerianField%DryDepositionRate(:, :, iTracer) = EulerianField%TotalDepositionRate(:, :, iTracer) +                 &
                                                     0.5 * (  EulerianField%Concentration(1:nX, 1:nY, 1, iTracer, iOld) &
                                                              * EulerianField%Vd(:, :, iTracer, iOld)                   &
                                                            + EulerianField%Concentration(1:nX, 1:nY, 1, iTracer, iNew) &
                                                              * EulerianField%Vd(:, :, iTracer, iNew)                   &
                                                           )
    EulerianField%TotalDepositionRate(:,:, iTracer) = EulerianField%DryDepositionRate(:, :, iTracer)
  End If

  ! Wet deposition 
  If ( WetDep ) Then
    Do ix = 1, nX
      Do iy = 1, nY
        EulerianField%WetDepositionRate(ix, iy, iTracer) = 0.0
        Do iz = 1, nZ  
          ! Calculate wet deposition rate for output 
          EulerianField%WetDepositionRate(ix, iy, iTracer)                          &
            = EulerianField%WetDepositionRate(ix, iy, iTracer) +                    &
              EulerianField%Concentration(ix, iy, iz, iTracer, iOld)                &
              * (1.0 - Exp(-TimeStep * EulerianField%Lambda(ix, iy, iz, iTracer)))  &
              * EulerianField%ZGridBoundary%AvZ(iz)/TimeStep        
          ! Change in concentration due to wet deposition 
          EulerianField%Concentration(ix, iy, iz, iTracer, iNew)         &
            = EulerianField%Concentration(ix, iy, iz, iTracer, iNew) *   &
              Exp(-TimeStep * EulerianField%Lambda(ix, iy, iz, iTracer))
        End Do
        EulerianField%TotalDepositionRate(ix, iy, iTracer)                                                        &
          = EulerianField%TotalDepositionRate(ix, iy, iTracer) + EulerianField%WetDepositionRate(ix, iy, iTracer)
      End Do
    End Do
  
  End If

End Subroutine EulerianDeposition

!-------------------------------------------------------------------------------------------------------------

Subroutine EulerianAgentDecay(Time, TimeStep, EulerianField, Specieses, Coords, Grids, Domains, &
                              Mets, Flows, FlowMemory, Units)

  Implicit None

  ! Argument list:
  Type(Time_),          Intent(In)    :: Time
  Real(Std),            Intent(In)    :: TimeStep
  Type(EulerianField_), Intent(InOut) :: EulerianField
  Type(Specieses_),     Intent(In)    :: Specieses
  Type(Coords_),        Intent(In)    :: Coords
  Type(Grids_),         Intent(In)    :: Grids
  Type(Domains_),       Intent(In)    :: Domains
  Type(Mets_),          Intent(InOut) :: Mets
  Type(Flows_),         Intent(InOut) :: Flows
  Type(FlowMemory_),    Intent(InOut) :: FlowMemory
  Type(Units_),         Intent(InOut) :: Units

  ! Locals: 
  Integer         :: ix
  Integer         :: iy
  Integer         :: iz
  Integer         :: nX
  Integer         :: nY
  Integer         :: nZ
  Integer         :: iNew
  Integer         :: ErrorCode
  Type(Position_) :: Position
  Type(Flow_)     :: Flow
  Type(Cloud_)    :: Cloud
  Type(Rain_)     :: Rain
  Type(Plant_)    :: Plant
  Type(Soil_)     :: Soil
  Type(Surface_)  :: Surface
  Real(Std)       :: X(3)
  Real(Std)       :: Pressure
  Real(Std)       :: RelHumidity
  Real(Std)       :: ZenithAngle
  Real(Std)       :: DummyConcentration(Specieses%nFields)
  Real(Std)       :: Temperature

  nX = EulerianField%HGrid%nX
  nY = EulerianField%HGrid%nY
  nZ = EulerianField%ZGrid%nZ
  iNew = EulerianField%iNew

  Do ix = 1, nX
    Do iy = 1, nY
      Do iz = 1, nZ

        X = (/ EulerianField%HGrid%X(ix), EulerianField%HGrid%Y(iy), EulerianField%ZGrid%Z(iz) /)
        Position = X2Position(Coords, X, EulerianField%HGrid%iHCoord, EulerianField%ZGrid%iZCoord)
        ZenithAngle = CalcZenithAngle(Coords, Time, Position)

        Call GetAttrib(                              &
               iAttribParam  = A_Rain,               &
               Coords        = Coords,               &
               Grids         = Grids,                &
               Domains       = Domains,              &
               Moisture      = .false.,              &
               Inhomog       = .true.,               &
               Homog         = .false.,              &
               Time          = Time2ShortTime(Time), &
               AnyTravelTime = .true.,               &
               TravelTime    = ZeroShortTime(),      &
               Position      = Position,             &
               Units         = Units,                &
               Mets          = Mets,                 &
               Flows         = Flows,                &
               FlowMemory    = FlowMemory,           &
               Flow          = Flow,                 &
               Cloud         = Cloud,                &
               Rain          = Rain,                 &
               Surface       = Surface,              &
               Soil          = Soil,                 &
               Plant         = Plant,                &
               ErrorCode     = ErrorCode             &
             )

        Temperature = EulerianField%Temperature(ix, iy, iz, iNew)
        Pressure = EulerianField%Density(ix, iy, iz, iNew) * GasConstant * Temperature
        RelHumidity = CalcRH(EulerianField%Humidity(ix, iy, iz, iNew), Temperature, Pressure)
        DummyConcentration(:) = Real(EulerianField%Concentration(ix, iy, iz, :, iNew), Std)
        Call AgentDecay(                    &
                        Time,               &
                        DummyConcentration, &
                        Specieses,          &
                        TimeStep,           &
                        Rain,               &
                        RelHumidity,        &
                        ZenithAngle,        &
                        Temperature,        &
                        Field = .True.      &
                       )
        EulerianField%Concentration(ix, iy, iz, :, iNew) = Real(DummyConcentration(:), Con)

      End Do
    End Do
  End Do

End Subroutine EulerianAgentDecay

!-------------------------------------------------------------------------------------------------------------

End Module EulerianInterfaceModule

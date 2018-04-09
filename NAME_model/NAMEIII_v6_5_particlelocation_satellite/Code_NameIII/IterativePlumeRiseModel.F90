! Module: Iterative Plume Rise Model

Module IterativePlumeRiseModel

! This module contains the iterative plume-rise model

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use SourceModule, Only: Source_
Use FlowsModule,  Only: Flows_, FlowMemory_, A_Flow, GetAttrib, ResetFlowMemory, &
                        Mets_, Flow_, Cloud_, Rain_, Surface_, Soil_, Plant_

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: IterativePlumeRiseVars_
Public :: FlowProfilesAtSource
Public :: IteratePlumeRiseModel

!-------------------------------------------------------------------------------------------------------------

! Module global parameters
Real(P64), Parameter :: dtFactor      = 1.e-3
Real(P64), Parameter :: TRef          = 273.0
Real(P64), Parameter :: RAir          = Real(GasConstant, P64)
Real(P64), Parameter :: gAcceleration = Real(Gravity, P64)
Real(P64), Parameter :: c_pa          = Real(Cp, P64)

Real(P64), Parameter :: Lp = 2.5e6         !! $$ get from elsewhere?

! From here on down eventually move to user-specified values
Real(P64), Parameter :: Alpha         = 0.1
Real(P64), Parameter :: Beta          = 0.5
Real(P64), Parameter :: m             = 1.5
                              
!-------------------------------------------------------------------------------------------------------------

Type :: IterativePlumeRiseVars_

  Integer                :: nZ
  Real(P64)              :: HObserved
  Real(P64), Allocatable :: Pressure(:)
  Real(P64), Allocatable :: Density(:)
  Real(P64), Allocatable :: Qv(:) 
  Real(P64), Allocatable :: Temperature(:)
  Real(P64), Allocatable :: Theta(:)
  Real(P64), Allocatable :: U(:) 
  Real(P64), Allocatable :: V(:) 
  Real(P64), Allocatable :: Z(:)

End Type IterativePlumeRiseVars_

!-------------------------------------------------------------------------------------------------------------

Type :: AmbientProfiles_
  Real(P64) :: U
  Real(P64) :: UOld
  Real(P64) :: V
  Real(P64) :: VOld
  Real(P64) :: dUdt
  Real(P64) :: dVdt
  Real(P64) :: Qv
  Real(P64) :: Pressure
  Real(P64) :: Density
  Real(P64) :: Temperature
End Type AmbientProfiles_

!-------------------------------------------------------------------------------------------------------------

Type :: PlumeProperties_
  Real(P64) :: Up
  Real(P64) :: Vp
  Real(P64) :: Wp
  Real(P64) :: UHat
  Real(P64) :: T
  Real(P64) :: RhoColumn
  Real(P64) :: Qm0
  Real(P64) :: F0
  Real(P64) :: GasFraction
  Real(P64) :: EntrainmentVelocity
  Real(P64) :: Ql
  Real(P64) :: QlOld
  Real(P64) :: dQldt
End Type PlumeProperties_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine FlowProfilesAtSource(Time, Source, Coords, Grids, Domains, &
                                Mets, Flows, FlowMemory, Units,       &
                                IterativePlumeRiseVars) 

  Implicit None

  ! Argument list:
  Type(ShortTime_),              Intent(In)           :: Time
  Type(Source_),                 Intent(InOut)        :: Source
  Type(Coords_),                 Intent(In)           :: Coords     
  Type(Grids_),                  Intent(In),   Target :: Grids
  Type(Domains_),                Intent(In)           :: Domains
  Type(Mets_),                   Intent(InOut)        :: Mets
  Type(Flows_),                  Intent(InOut)        :: Flows
  Type(FlowMemory_),             Intent(InOut)        :: FlowMemory
  Type(Units_),                  Intent(InOut)        :: Units
  Type(IterativePlumeRiseVars_), Intent(Out)          :: IterativePlumeRiseVars

  ! Locals:
  Integer                 :: k
  Integer                 :: nZ
  Real(Std)               :: X(3)
  Type(ZGrid_),   Pointer :: ZGrid
  Type(Position_)         :: Position
  Type(Flow_)             :: Flow
  Type(Cloud_)            :: Cloud
  Type(Rain_)             :: Rain
  Type(Plant_)            :: Plant
  Type(Soil_)             :: Soil
  Type(Surface_)          :: Surface
  Integer                 :: ErrorCode

  ZGrid => Grids%ZGrids(Source%iZGridIterativePlumeRise)
  nZ = ZGrid%nZ
  IterativePlumeRiseVars%nZ = nZ
  If ( .Not.Allocated(IterativePlumeRiseVars%Pressure) ) Then
    Allocate(IterativePlumeRiseVars%Pressure(nZ))
    Allocate(IterativePlumeRiseVars%Density(nZ))
    Allocate(IterativePlumeRiseVars%Qv(nZ))
    Allocate(IterativePlumeRiseVars%Temperature(nZ))
    Allocate(IterativePlumeRiseVars%Theta(nZ))
    Allocate(IterativePlumeRiseVars%U(nZ))
    Allocate(IterativePlumeRiseVars%V(nZ))
    Allocate(IterativePlumeRiseVars%Z(nZ))
  End If

  ! Interpolate meteorological fields to source location at target time
  Do k = 1, nZ

    Call ResetFlowMemory(Flows, FlowMemory)

    IterativePlumeRiseVars%Z(k) = ZGrid%Z(k)

    X = (/ Source%X(1), Source%X(2), ZGrid%Z(k) /)
    Position = X2Position(Coords, X, Source%iHCoord, Source%iZCoord)

    ! Get flow info.
    Call GetAttrib(                             & 
               iAttribParam  = A_Flow,          &
               Coords        = Coords,          &
               Grids         = Grids,           &
               Domains       = Domains,         &
               Moisture      = .true.,          &
               Inhomog       = .true.,          &
               Homog         = .false.,         &
               Time          = Time,            &
               AnyTravelTime = .true.,          &
               TravelTime    = ZeroShortTime(), &
               Position      = Position,        &
               Units         = Units,           &
               Mets          = Mets,            &
               Flows         = Flows,           &
               FlowMemory    = FlowMemory,      &
               Flow          = Flow,            &
               Cloud         = Cloud,           &
               Rain          = Rain,            &
               Surface       = Surface,         &
               Soil          = Soil,            &
               Plant         = Plant,           &
               ErrorCode     = ErrorCode        &
             )
        If (ErrorCode > 0) Then
          Write (*,*) 'Error in getting flow velocity', ErrorCode
          Stop
       End If

       IterativePlumeRiseVars%Pressure(k)    = Real(Flow%P,     P64)
       IterativePlumeRiseVars%Density(k)     = Real(Flow%Rho,   P64)
       IterativePlumeRiseVars%Qv(k)          = Real(Flow%Q,     P64)
       IterativePlumeRiseVars%Temperature(k) = Real(Flow%T,     P64)
       IterativePlumeRiseVars%Theta(k)       = Real(Flow%Theta, P64)
       IterativePlumeRiseVars%U(k)           = Real(Flow%U(1),  P64)
       IterativePlumeRiseVars%V(k)           = Real(Flow%U(2),  P64)

  End Do

End Subroutine FlowProfilesAtSource

!-------------------------------------------------------------------------------------------------------------

Subroutine IteratePlumeRiseModel(Source, IterativePlumeRiseVars, QmOriginal, QmFinal)

  Implicit None

  ! Argument list:
  Type(Source_),                 Intent(In)    :: Source
  Type(IterativePlumeRiseVars_), Intent(InOut) :: IterativePlumeRiseVars
  Real(Std),                     Intent(Out)   :: QmOriginal
  Real(Std),                     Intent(Out)   :: QmFinal                 

  ! Locals: 
  Integer              :: it
  Integer,   Parameter :: NIterations = 100
  Real(P64), Parameter :: Tol = 1.e-2
  Real(P64)            :: Qm
  Real(P64)            :: Qm0
  Real(P64)            :: Difference
  Real(P64)            :: Difference0
  Real(P64)            :: ZMax
  Real(P64)            :: ZMax0
  Real(P64)            :: ZMaxZeroWind
  Real(P64)            :: ZFinal
  Real(P64)            :: ZEq
  Real(P64)            :: Radius
  Real(P64)            :: SummitElevation
  Real(P64)            :: HObserved
  Real(P64)            :: BuoyancyFrequency
  Real(P64)            :: F0
  Real(Std)            :: GasFraction0
  
  ! Initial gas fraction
  GasFraction0 = Source%GasFraction0 

  ! Ground level of source 
  SummitElevation = Real(Source%SummitElevation, P64)

  ! Target height  
  HObserved = Real(Source%PlumeRiseHeight, P64)
  IterativePlumeRiseVars%HObserved = HObserved

  ! Initial mass flux estimate 
  If ( Source%SpeciesNames(1) .CIEq. 'VOLCANIC_ASH' ) Then
    Qm0 = 141.0 * ((HObserved - SummitElevation)/1000.0)**4.15    ! NB: use Mastin formula for volcanic ash
  Else
    Call CalcBuoyancyFrequency(IterativePlumeRiseVars, SummitElevation, BuoyancyFrequency)
    F0  = HObserved**4 * Alpha**2 * BuoyancyFrequency**3/3.42
    Qm0 = F0 * TRef/(gAcceleration * (Real(Source%Temperature, P64) - IterativePlumeRiseVars%Temperature(1)))
    ! NB: assumes Z(1) is ground level
  End If
  QmOriginal = Real(Qm0, Std)

  ! Initial rise height 
  Call IntegratePlumeEquations(Qm0, ZMaxZeroWind, ZEq, Radius, IterativePlumeRiseVars, Source, ZeroWind = .True.)
  Call IntegratePlumeEquations(Qm0, ZMax,         ZEq, Radius, IterativePlumeRiseVars, Source, ZeroWind = .False.)
  ZMax = Min(ZMaxZeroWind, ZMax + Radius)

  Difference0 = ZMax - HObserved
  ZMax0 = ZMax
  If ( Abs(Difference0) < Tol ) Then
    If ( Source%SpeciesNames(1) .CIEq. 'VOLCANIC_ASH' ) Then
      QmFinal = (1.0 - GasFraction0) * Real(Qm0, Std)  
    Else
      QmFinal = Real(Qm0, Std)  
    End If
    Return
  End If

  ! Calculate value of Qm for which Difference changes sign
  Do it = 1, NIterations
    Qm = Qm0 * (HObserved/ZMax)**4 
    Call IntegratePlumeEquations(Qm, ZMaxZeroWind, ZEq, Radius, IterativePlumeRiseVars, Source, ZeroWind = .True.)
    Call IntegratePlumeEquations(Qm, ZMax,         ZEq, Radius, IterativePlumeRiseVars, Source, ZeroWind = .False.)
    ZMax = Min(ZMaxZeroWind, ZMax + Radius)
    Difference = ZMax - HObserved
    If ( Difference0 * Difference < 0.0 ) Then 
      Exit
    Else
      Difference0 = Difference
      Qm0 = Qm
      ZMax0 = ZMax
    End If 
  End Do

  If ( Abs(ZMax - HObserved) < Tol ) Then
    QmFinal = Real(Qm, Std)
    ZFinal = ZMax
    Return
  Else
    If ( ZMax0 < ZMax ) Then
      Call Bisection(Qm0, Qm, ZMax0, ZEq, NIterations, IterativePlumeRiseVars, Source)
      If ( Source%SpeciesNames(1) .CIEq. 'VOLCANIC_ASH' ) Then
        QmFinal = (1.0 - GasFraction0) * Real(Qm, Std)   
      Else
        QmFinal = Real(Qm0, Std)  
      End If
      ZFinal = ZMax0
    Else
      Call Bisection(Qm, Qm0, ZMax, ZEq, NIterations, IterativePlumeRiseVars, Source)
      If ( Source%SpeciesNames(1) .CIEq. 'VOLCANIC_ASH' ) Then
        QmFinal = (1.0 - GasFraction0) * Real(Qm, Std)   
      Else
        QmFinal = Real(Qm0, Std)  
      End If
      ZFinal = ZMax
    End If
  End If

End Subroutine IteratePlumeRiseModel

!-------------------------------------------------------------------------------------------------------------

 Subroutine Bisection(Qm0, Qm, ZMax, ZEq, NIterations, IterativePlumeRiseVars, Source)
  
  Implicit None

  ! Argument list:
  Real(P64),                     Intent(InOut) :: Qm0
  Real(P64),                     Intent(InOut) :: Qm
  Real(P64),                     Intent(InOut) :: ZMax
  Real(P64),                     Intent(Out)   :: ZEq
  Integer,                       Intent(In)    :: NIterations
  Type(IterativePlumeRiseVars_), Intent(In)    :: IterativePlumeRiseVars
  Type(Source_),                 Intent(In)    :: Source

  ! Locals:
  Integer              :: it
  Real(P64)            :: P
  Real(P64)            :: Difference
  Real(P64)            :: Difference0
  Real(P64)            :: ZMaxZeroWind
  Real(P64)            :: Radius
  Real(P64), Parameter :: Tol1 = 1.e-2
  Real(P64), Parameter :: Tol2 = 1.e-20

  Difference0 = ZMax - IterativePlumeRiseVars%HObserved
 
  it = 1
  Do While (it <= NIterations) 
    P = Qm0 + (Qm - Qm0)/2.0
    Call IntegratePlumeEquations(P, ZMaxZeroWind, ZEq, Radius, IterativePlumeRiseVars, Source, ZeroWind = .True.)
    Call IntegratePlumeEquations(P, ZMax,         ZEq, Radius, IterativePlumeRiseVars, Source, ZeroWind = .False.)
    ZMax = Min(ZMaxZeroWind, ZMax + Radius)
    Difference = ZMax - IterativePlumeRiseVars%HObserved
    If ( Abs(Difference) < Tol1 .Or. (Qm - Qm0)/2.0 < Tol2 ) Then
      Qm = P
      Exit  
    Else If ( Difference0 * Difference > 0.0 ) Then
      Qm0 = P
      Difference0 = Difference
    Else
      Qm = P
    End If
    it = it + 1
  End Do       
  If ( it > NIterations ) Then
    Write (*,*) 'Error: bisection method failed to converge'
    ! Need to add error code
  End If

 End Subroutine Bisection

!-------------------------------------------------------------------------------------------------------------

 Subroutine IntegratePlumeEquations(Qm, ZMax, ZEq, Radius, IterativePlumeRiseVars, Source, ZeroWind)

  Implicit None

  ! Argument list:
  Real(P64),                     Intent(In)  :: Qm
  Real(P64),                     Intent(Out) :: ZMax
  Real(P64),                     Intent(Out) :: ZEq
  Real(P64),                     Intent(Out) :: Radius
  Type(IterativePlumeRiseVars_), Intent(In)  :: IterativePlumeRiseVars
  Type(Source_),                 Intent(In)  :: Source
  Logical,                       Intent(In)  :: ZeroWind

  ! Locals:
  Integer,   Parameter   :: NVariables = 6
  Integer                :: it
  Integer                :: j
  Real(P64), Parameter   :: Tol = 1.e-15
  Real(P64)              :: P(NVariables)
  Real(P64)              :: Q(NVariables)
  Real(P64)              :: K(4, NVariables)
  Real(P64)              :: Z
  Real(P64)              :: dt
  Real(P64)              :: OldW
  Real(P64)              :: OldDensity
  Real(P64)              :: OldZ
  Real(P64)              :: BuoyancyFrequency
  Logical                :: CalculateLevelNeutralBuoyancy
  Type(AmbientProfiles_) :: AmbientProfiles
  Type(PlumeProperties_) :: PlumeProperties

  ! Variables
  ! P(1) = Qm (Mass flux)
  ! P(2) = Mw (Vertical momentum flux)
  ! P(3) = E  (Enthalpy flux)
  ! P(4) = Qt (Qv)
  ! P(5) = Mu (Horizontal momentum flux: x direction)
  ! P(6) = Mv (Horizontal momentum flux: y direction)

  ! Initialise variables
  Z = Real(Source%SummitElevation, P64)
  OldZ = Z

  Call AmbientProfileAtZ(Z, AmbientProfiles, IterativePlumeRiseVars)
  If ( ZeroWind ) Then 
    AmbientProfiles%U = 0.0
    AmbientProfiles%V = 0.0
  End If
  AmbientProfiles%UOld = AmbientProfiles%U
  AmbientProfiles%VOld = AmbientProfiles%V
  AmbientProfiles%dUdt = 0.0                 ! $$ Estimate from w du/dz using point below source? 
  AmbientProfiles%dVdt = 0.0

  Call InitialisePlumeVariables(Qm, AmbientProfiles, P, PlumeProperties, Source)

  OldW = PlumeProperties%Wp
  OldDensity = PlumeProperties%RhoColumn
  Call EntrainmentVelocity(AmbientProfiles, PlumeProperties)

  CalculateLevelNeutralBuoyancy = .False.

  ! Time step (1/w)
  Call CalcBuoyancyFrequency(IterativePlumeRiseVars, Z, BuoyancyFrequency)
  dt = dtFactor * PlumeProperties%F0**0.25 * BuoyancyFrequency**(-0.75)/Real(Source%FlowVelocity, P64)

  ! Integrate plume equations
  it = 1
  Do While (.True.) 

    Do j = 1, NVariables
      K(1,j) = dt * F(j, P, AmbientProfiles, PlumeProperties)
    End Do

    Q(:) = P(:) + 0.5 * K(1,:)
    Do j = 1, NVariables
      K(2,j) = dt * F(j, Q, AmbientProfiles, PlumeProperties)
    End Do

    Q(:) = P(:) + 0.5 * K(2,:)
    Do j = 1, NVariables
      K(3,j) = dt * F(j, Q, AmbientProfiles, PlumeProperties)
    End Do

    Q(:) = P(:) + K(3,:)
    Do j = 1, NVariables
      K(4,j) = dt * F(j, Q, AmbientProfiles, PlumeProperties)
    End Do

    P(:) = P(:) + (K(1,:) + 2.0 * K(2,:) + 2.0 * K(3,:) + K(4,:))/6.0

    PlumeProperties%Up = AmbientProfiles%U + P(5)/P(1)
    PlumeProperties%Vp = AmbientProfiles%V + P(6)/P(1)
    PlumeProperties%Wp = P(2)/P(1)
    Z = Z + PlumeProperties%Wp * dt

    Call AmbientProfileAtZ(Z, AmbientProfiles, IterativePlumeRiseVars)
    If ( ZeroWind ) Then 
      AmbientProfiles%U = 0.0
      AmbientProfiles%V = 0.0
    End If
    Call ColumnThermodynamics(P, dt, AmbientProfiles, PlumeProperties, Real(Source%Density, P64), Real(Source%GasFraction0, P64))

    If ( Abs(OldW - PlumeProperties%Wp) < Tol .Or. PlumeProperties%Wp < 0.0 ) Then 
      Exit
    Else
      Call AmbientGradients(dt, AmbientProfiles)
      Call EntrainmentVelocity(AmbientProfiles, PlumeProperties)
      OldW = PlumeProperties%Wp
      PlumeProperties%QlOld = PlumeProperties%Ql
      AmbientProfiles%UOld  = AmbientProfiles%U
      AmbientProfiles%VOld  = AmbientProfiles%V
      it = it + 1
    End If
    Radius = Sqrt( P(1)/(Pi * PlumeProperties%RhoColumn * PlumeProperties%UHat) )

    ! Calculate level of neutral buoyancy
    If ( .Not.CalculateLevelNeutralBuoyancy .And. PlumeProperties%RhoColumn < AmbientProfiles%Density ) Then
      CalculateLevelNeutralBuoyancy = .True.
    End If
    If ( CalculateLevelNeutralBuoyancy .And. PlumeProperties%RhoColumn > AmbientProfiles%Density ) Then 
      ZEq = (Z - OldZ) * (AmbientProfiles%Density - OldDensity)/( PlumeProperties%RhoColumn - OldDensity) + OldZ
      CalculateLevelNeutralBuoyancy = .False.
    End If
    OldDensity = PlumeProperties%RhoColumn
    OldZ = Z

  End Do
  ZMax = Z
 
 End Subroutine IntegratePlumeEquations

!--------------------------------------------------------------------------------------------------------

 Function F(j, P, AmbientProfiles, PlumeProperties)

   Implicit None

   ! Argument list:
   Integer,                Intent(In)         :: j
   Real(P64),              Intent(In)         :: P(:)
   Type(AmbientProfiles_), Intent(In), Target :: AmbientProfiles
   Type(PlumeProperties_), Intent(In), Target :: PlumeProperties 

   ! Function result:
   Real(P64) :: F

   ! Locals: 
   Real(P64)          :: E
   Real(P64)          :: RHSEnergyEquation
   Real(P64)          :: WindSpeed
   Real(P64)          :: PlumeRadius
   Real(P64), Pointer :: RhoColumn
   Real(P64), Pointer :: RhoAmbient
   Real(P64), Pointer :: Pressure
   Real(P64), Pointer :: Up
   Real(P64), Pointer :: Vp
   Real(P64), Pointer :: Wp
   Real(P64), Pointer :: UHat
   Real(P64), Pointer :: q_v

   RhoColumn   => PlumeProperties%RhoColumn
   Up          => PlumeProperties%Up
   Vp          => PlumeProperties%Vp
   Wp          => PlumeProperties%Wp
   UHat        => PlumeProperties%UHat
   Pressure    => AmbientProfiles%Pressure
   RhoAmbient  => AmbientProfiles%Density
   q_v         => AmbientProfiles%Qv

   WindSpeed   = Sqrt(AmbientProfiles%U**2 + AmbientProfiles%V**2)
   PlumeRadius = Sqrt( P(1)/(Pi * RhoColumn * UHat) )

   E = 2.0 * Pi * PlumeRadius * PlumeProperties%EntrainmentVelocity * UHat * RhoAmbient * &
       Sqrt( Max(1.0, RhoColumn/RhoAmbient)) 

   RHSEnergyEquation = ((1.0 - q_v) * c_pa + q_v * c_pv) * AmbientProfiles%Temperature * E  &
                       - gAcceleration * RhoAmbient * Wp * P(1)/RhoColumn                         &
                       + (Lp - 273.0 * (c_pv - c_pw) ) * PlumeProperties%dQldt
 
   If ( j == 1 ) Then                                                           ! Mass flux 
     F = E                                                     
   Else If ( j == 2 ) Then                                                      ! Vertical momentum flux
     F = (RhoAmbient - RhoColumn) * gAcceleration * Pi * PlumeRadius**2 * UHat
   Else If ( j == 3 ) Then                                                      ! Energy
     F = RHSEnergyEquation
   Else If ( j == 4 ) Then                                                      ! Total moisture flux
     F = E * AmbientProfiles%Qv     
   Else If ( j == 5 ) Then
     F = -P(1) * AmbientProfiles%dUdt                                           ! Horizontal momentum flux: x direction
   Else If ( j == 6 ) Then
     F = -P(1) * AmbientProfiles%dVdt                                           ! Horizontal momentum flux: y direction
   End If

 End Function

!----------------------------------------------------------------------------------------------

 Subroutine EntrainmentVelocity(AmbientProfiles, PlumeProperties)

   Implicit None

   ! Argument list:
   Type(AmbientProfiles_), Intent(In)            :: AmbientProfiles
   Type(PlumeProperties_), Intent(InOut), Target :: PlumeProperties 

   ! Locals: 
   Real(P64)          :: CoefficientParallel
   Real(P64)          :: DeltaUParallel
   Real(P64)          :: DeltaUPerpendicular
   Real(P64)          :: DeltaUPerpendicularX
   Real(P64)          :: DeltaUPerpendicularY
   Real(P64)          :: DeltaUPerpendicularZ
   Real(P64), Pointer :: UHat

   UHat => PlumeProperties%UHat

   UHat = Sqrt(PlumeProperties%Up**2 + PlumeProperties%Vp**2 + PlumeProperties%Wp**2)

   CoefficientParallel = (PlumeProperties%Up - AmbientProfiles%U) * PlumeProperties%Up + &
                         (PlumeProperties%Vp - AmbientProfiles%V) * PlumeProperties%Vp + &
                         PlumeProperties%Wp * PlumeProperties%Wp
   CoefficientParallel = CoefficientParallel/UHat**2
   DeltaUParallel = Abs(CoefficientParallel) * UHat
   DeltaUPerpendicularX = PlumeProperties%Up - AmbientProfiles%U - CoefficientParallel * PlumeProperties%Up
   DeltaUPerpendicularY = PlumeProperties%Vp - AmbientProfiles%V - CoefficientParallel * PlumeProperties%Vp
   DeltaUPerpendicularZ = PlumeProperties%Wp - CoefficientParallel * PlumeProperties%Wp
   DeltaUPerpendicular  = Sqrt(DeltaUPerpendicularX**2 + DeltaUPerpendicularY**2 + DeltaUPerpendicularZ**2)
   PlumeProperties%EntrainmentVelocity = ( (Alpha * DeltaUParallel)**m + (Beta * DeltaUPerpendicular)**m )**(1.0/m)

 End Subroutine 

!----------------------------------------------------------------------------------------------

 Subroutine InitialisePlumeVariables(Qm, AmbientProfiles, P, PlumeProperties, Source)

   Implicit None

   ! Argument list:
   Real(P64),              Intent(In)          :: Qm
   Type(AmbientProfiles_), Intent(In)          :: AmbientProfiles
   Real(P64),              Intent(Out)         :: P(:)
   Type(PlumeProperties_), Intent(Out), Target :: PlumeProperties 
   Type(Source_),          Intent(In)          :: Source

   ! Locals: 
   Real(P64), Pointer :: RhoColumn
   Real(P64)          :: n_v
   Real(P64)          :: q_v
   Real(P64)          :: RhoGas
   Real(P64)          :: BulkGasConstant
   Real(P64)          :: Qv0
   Real(P64)          :: W0
   Real(P64)          :: T0
   Real(P64)          :: RhoSolid
   Real(P64)          :: GasFraction0
   Real(P64)          :: WaterVapourFraction0

   RhoColumn => PlumeProperties%RhoColumn
   
   W0                   = Real(Source%FlowVelocity,         P64)
   T0                   = Real(Source%Temperature,          P64)
   RhoSolid             = Real(Source%Density,              P64)
   WaterVapourFraction0 = Real(Source%WaterVapourFraction0, P64)
   GasFraction0         = Real(Source%GasFraction0,         P64) 

   ! Set initial mass flux
   P(1) = Qm
   PlumeProperties%Qm0 = Qm

   ! Set initial momentum flux
   P(2) = Qm * W0

   ! Initial water vapour flux
   Qv0 = WaterVapourFraction0 * Qm
   P(4) = Qv0

   ! Initial enthalpy flux
   n_v = Qv0/Qm
   P(3) = (1.0 - GasFraction0) * c_ps * T0 * Qm + (GasFraction0 - n_v) * c_pa * T0 * Qm + n_v * c_pv * T0 * Qm

   ! Initial horizontal momentum fluxes
   P(5) = -AmbientProfiles%U * Qm
   P(6) = -AmbientProfiles%V * Qm

   ! Initial plume velocities
   PlumeProperties%Up = AmbientProfiles%U + P(5)/P(1)
   PlumeProperties%Vp = AmbientProfiles%V + P(6)/P(1)
   PlumeProperties%Wp = P(2)/P(1)
   PlumeProperties%UHat = Sqrt(PlumeProperties%Up**2 + PlumeProperties%Vp**2 + PlumeProperties%Wp**2)
 
   ! Initial plume density
   q_v = n_v/GasFraction0
   BulkGasConstant = q_v * RVapour + (1.0 - q_v) * RAir     
   RhoGas = AmbientProfiles%Pressure/(BulkGasConstant * T0)
   RhoColumn = GasFraction0/RhoGas + (1.0 - GasFraction0)/RhoSolid   
   RhoColumn = 1.0/RhoColumn
   
   ! Initial temperature and gradient
   PlumeProperties%T = T0

   ! Initial buoyancy flux: assumes all the heat is in the gas phase
   PlumeProperties%F0 = gAcceleration * Qm * (T0 - AmbientProfiles%Temperature)/TRef

   ! Initial gas fraction
   PlumeProperties%GasFraction = GasFraction0

   ! Initial liquid water flux
   PlumeProperties%Ql = 0.0
   PlumeProperties%QlOld = 0.0
   PlumeProperties%dQldt = 0.0

 End Subroutine InitialisePlumeVariables

!----------------------------------------------------------------------------------------------
 
 Subroutine AmbientProfileAtZ(Z, AmbientProfiles, IterativePlumeRiseVars)

   Implicit None

   ! Argument list:
   Real(P64),                     Intent(In)  :: Z
   Type(AmbientProfiles_),        Intent(Out) :: AmbientProfiles
   Type(IterativePlumeRiseVars_), Intent(In)  :: IterativePlumeRiseVars

   ! Locals: 
   Integer   :: iz
   Real(P64) :: ZLevel_iz
   Real(P64) :: ZLevel_iz_1

   ! Determine Z level for interpolation
   Do iz = 1, IterativePlumeRiseVars%nZ - 1
     If ( Z >= IterativePlumeRiseVars%Z(iz) .And. Z <= IterativePlumeRiseVars%Z(iz + 1) ) Exit
   End Do
   If ( iz >= IterativePlumeRiseVars%nZ ) Then
     Write (*,*) 'Error: Z outside range'
     ! Error code 
   End If
   ZLevel_iz   = IterativePlumeRiseVars%Z(iz) 
   ZLevel_iz_1 = IterativePlumeRiseVars%Z(iz + 1)

   ! Interpolate to Z 
   Call Interpolate(ZLevel_iz, ZLevel_iz_1, Z, IterativePlumeRiseVars%U(iz),                 &
                    IterativePlumeRiseVars%U(iz + 1), AmbientProfiles%U)
   Call Interpolate(ZLevel_iz, ZLevel_iz_1, Z, IterativePlumeRiseVars%V(iz),                 &
                    IterativePlumeRiseVars%V(iz + 1), AmbientProfiles%V)
   Call Interpolate(ZLevel_iz, ZLevel_iz_1, Z, IterativePlumeRiseVars%Qv(iz),                &
                    IterativePlumeRiseVars%Qv(iz + 1), AmbientProfiles%Qv)
   Call Interpolate(ZLevel_iz, ZLevel_iz_1, Z, IterativePlumeRiseVars%Pressure(iz),          &
                    IterativePlumeRiseVars%Pressure(iz + 1), AmbientProfiles%Pressure)
   Call Interpolate(ZLevel_iz, ZLevel_iz_1, Z, IterativePlumeRiseVars%Density(iz),           &
                    IterativePlumeRiseVars%Density(iz + 1), AmbientProfiles%Density)
   Call Interpolate(ZLevel_iz, ZLevel_iz_1, Z, IterativePlumeRiseVars%Temperature(iz),       &
                    IterativePlumeRiseVars%Temperature(iz + 1), AmbientProfiles%Temperature)

 End Subroutine AmbientProfileAtZ

!----------------------------------------------------------------------------------------------

 Subroutine AmbientGradients(dt, AmbientProfiles)

   Implicit None

   ! Argument list:
   Real(P64),              Intent(In)    :: dt
   Type(AmbientProfiles_), Intent(InOut) :: AmbientProfiles

   ! Calculate gradients
   AmbientProfiles%dUdt = (AmbientProfiles%U - AmbientProfiles%UOld)/dt
   AmbientProfiles%dVdt = (AmbientProfiles%V - AmbientProfiles%VOld)/dt

 End Subroutine AmbientGradients

!----------------------------------------------------------------------------------------------

 Subroutine Interpolate(X1, X2, X, Y1, Y2, Y)

   Implicit None

   ! Argument list: 
   Real(P64), Intent(In)  :: X1 
   Real(P64), Intent(In)  :: X2
   Real(P64), Intent(In)  :: X
   Real(P64), Intent(In)  :: Y1
   Real(P64), Intent(In)  :: Y2
   Real(P64), Intent(Out) :: Y

   ! Interpolate Y to required value of X
   Y = (Y2 - Y1) * (X - X1)/(X2 - X1) + Y1

 End Subroutine Interpolate

!----------------------------------------------------------------------------------------

 Function Rs(T, Pressure)

   Implicit None

   ! Argument list:
   Real(P64), Intent(In) :: T
   Real(P64), Intent(In) :: Pressure

   ! Function result:
   Real(P64) :: Rs

   Rs = 3.8 * Exp(17.65 * (T - 273.15)/(T - 29.65))/(Pressure/100.0) !NB: pressure in hPa/mbar
 
 End Function Rs

!----------------------------------------------------------------------------------------------

 Subroutine ColumnThermodynamics(P, dt, AmbientProfiles, PlumeProperties, RhoSolid, GasFraction0)

   Implicit None

   ! Argument list:
   Real(P64),              Intent(In)            :: P(:)
   Real(P64),              Intent(In)            :: dt
   Type(AmbientProfiles_), Intent(In),    Target :: AmbientProfiles
   Type(PlumeProperties_), Intent(InOut), Target :: PlumeProperties
   Real(P64),              Intent(In)            :: RhoSolid
   Real(P64),              Intent(In)            :: GasFraction0

   ! Locals: 
   Real(P64)          :: n_t
   Real(P64)          :: n_l
   Real(P64)          :: n_a
   Real(P64)          :: n_v
   Real(P64)          :: q_v
   Real(P64)          :: q_l
   Real(P64)          :: c_pc
   Real(P64)          :: RhoGas
   Real(P64)          :: BulkGasConstant
   Real(P64), Pointer :: Pressure
   Real(P64), Pointer :: RhoColumn
   Real(P64), Pointer :: T   
   Real(P64), Pointer :: GasFraction
   Real(P64), Pointer :: Ql

   Pressure    => AmbientProfiles%Pressure
   RhoColumn   => PlumeProperties%RhoColumn
   T           => PlumeProperties%T
   GasFraction => PlumeProperties%GasFraction
   Ql          => PlumeProperties%Ql

   ! Mass fractions
   n_t = P(4)/P(1)
   n_a = 1.0 - n_t - (1.0 - GasFraction0) * PlumeProperties%Qm0/P(1)
   n_l = Max(0.0, n_t - n_a * Rs(T, Pressure))
   n_v = n_t - n_l
   GasFraction = n_a + n_v
   q_v = n_v/GasFraction
   q_l = n_l/GasFraction

   ! Liquid water flux
   Ql = P(1) * n_l

   ! Specific heat capacity of plume
   c_pc = n_a * c_pa + n_v * c_pv + n_l * c_pw + (1.0 - n_a - n_v - n_l) * c_ps

   ! Plume temperature
   T = P(3)/(P(1) * c_pc)
 
   ! Bulk gas constant for plume
   BulkGasConstant = q_v * RVapour + (1.0 - q_v) * RAir
     
   ! Density of gas constituent of plume
   RhoGas = Pressure/(BulkGasConstant * T)

   ! Density of eruption column
   RhoColumn = GasFraction/RhoGas + (1.0 - n_a - n_v - n_l)/RhoSolid + n_l/RhoWater
   RhoColumn = 1.0/RhoColumn

   ! Liquid water flux gradient
   PlumeProperties%dQldt = (Ql - PlumeProperties%QlOld)/dt

 End Subroutine ColumnThermodynamics

!----------------------------------------------------------------------------------------------

 Subroutine CalcBuoyancyFrequency(IterativePlumeRiseVars, SummitElevation, BuoyancyFrequency)

   Implicit None

   ! Argument list:
   Type(IterativePlumeRiseVars_), Intent(In)  :: IterativePlumeRiseVars
   Real(P64),                     Intent(In)  :: SummitElevation
   Real(P64),                     Intent(Out) :: BuoyancyFrequency

   ! Locals: 
   Integer   :: mH
   Integer   :: nH
   Integer   :: iz
   Real(P64) :: ZLevel(IterativePlumeRiseVars%nZ)
   Real(P64) :: ThetaProfile(IterativePlumeRiseVars%nZ)
   Real(P64) :: ThetaGradient

   ThetaProfile(:) = IterativePlumeRiseVars%Theta(:)
   ZLevel(:)       = IterativePlumeRiseVars%Z(:)

   mH = 1
   nH = IterativePlumeRiseVars%nZ
   Do iz = 1, IterativePlumeRiseVars%nZ
     If ( ZLevel(iz) > SummitElevation ) Then
       mH = Min(iz - 1, 1)
       Exit
     End If
   End Do
   Do iz = mH, IterativePlumeRiseVars%nZ
     If ( ZLevel(iz) > IterativePlumeRiseVars%HObserved ) Then
       nH = iz
       Exit
     End If
   End Do

   ! Buoyancy frequency     
   !  - using least squares fit (B&F p.476)
   ThetaGradient = ( Real(nH - mH + 1) * Sum(ThetaProfile(mH:nH) * ZLevel(mH:nH)) -      &
                     Sum(ZLevel(mH:nH)) * Sum(ThetaProfile(mH:nH)) ) /                   &
                   ( Real(nH - mH + 1) * Sum(ZLevel(mH:nH)**2) - Sum(ZLevel(mH:nH))**2 )
   BuoyancyFrequency = Sqrt(gAcceleration * ThetaGradient/TRef)

 End Subroutine CalcBuoyancyFrequency

!----------------------------------------------------------------------------------------------

End Module IterativePlumeRiseModel

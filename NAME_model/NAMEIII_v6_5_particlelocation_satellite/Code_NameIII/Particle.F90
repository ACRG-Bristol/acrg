! Module:  Particle Module

Module ParticleModule

! This module provides code to treat particles.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowsModule, Only: Flows_, Flow_, Cloud_, Rain_, Surface_, Plant_
Use SizeDistModule
Use SpeciesModule
Use SourceModule
Use OutputModule, Only: OutputOpts_, Reqs_, Results_, CalciTA, OutputPPInfos
Use PlumeRiseModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: Convection_None
Public :: Convection_Old
Public :: Convection_New
Public :: Particle_
Public :: Extra_
Public :: ParticleActive
Public :: ParticleReallyActive
Public :: UseParticleForOutput
Public :: MarkParticle
Public :: ParticleTime
Public :: OldParticleTime
Public :: ParticleMaxTime
Public :: ParticleMinTime
Public :: OldParticleMinTime
Public :: ParticleInstTime
Public :: TravelTime
Public :: OldTravelTime
Public :: InactivateParticle
Public :: InitParticle
Public :: InitParticle2
Public :: CalcRandomWalkDt
Public :: EvolveParticle
Public :: XHomogSdeTerms
Public :: XInhomogSdeTerms
Public :: XUSdeTerms
Public :: InitW
Public :: ReflectW
Public :: UpdateTime
Public :: DiffusionProcessReflect
Public :: RadioactiveDecay
Public :: AgentDecay
Public :: DryDeposition
Public :: WetDeposition
Public :: DeepConvectionNew
Public :: DeepConvectionOld
Public :: ParticleEulerianFields
Public :: ParticleConc
Public :: ParticleCloudGamma
Public :: ParticleNumbers
Public :: ParticleNumbersBySpecies
Public :: ParticleMeanTravelTime
Public :: ParticleMass
Public :: ParticleMeanZ
Public :: ParticleSigZ2
Public :: ParticleXStats
Public :: ParticleInfo

!-------------------------------------------------------------------------------------------------------------

! Codes for convection scheme.
Integer, Parameter :: Convection_None = 0 ! No deep convection scheme used.
Integer, Parameter :: Convection_Old  = 1 ! Old deep convection scheme used.
Integer, Parameter :: Convection_New  = 2 ! New deep convection scheme used.

!-------------------------------------------------------------------------------------------------------------

Type :: Particle_ ! Information defining a particle.
  ! If changing this type, remember restart read and write routines.
  Integer          :: iUP          ! Unique index of particle (i.e. an index unique to the particle
                                   ! which is not
                                   ! reused despite particle recycling). Should be Integer(8) $$
  Real(Std)        :: X(3)         ! Particle position.
  Real(Std)        :: XOld(3)      ! Particle position at previous time step.
  Type(ShortTime_) :: T            ! Travel time.
  Type(ShortTime_) :: TOld         ! Travel time at previous time step.
  Type(ShortTime_) :: T0           ! Release time.
  Logical          :: NeedToSetVel ! Indicates that the particle velocity needs initialising. (e.g. after deep
                                   ! convection).
  Logical          :: Active       ! Indicates particle is active.
  Integer          :: Marked       ! Non-zero values indicate particle is marked to be inactivated at the end
                                   ! of the sync time.
                                   ! 1 = contribute to output
                                   ! 2 = don't contrib to output
                                   ! Note Active is only set to false after output at the end of the sync time
                                   ! (e.g. to avoid re-use)
                                   ! 1 & 2 imply no more time stepping of particle.
  ! Note particles and puffs are inactivated in 3 ways:
  ! 1) Particle/puff has Active set to false and is returned to free particle/puff stack. This
  !    generally only occurs at the end of the sync time because:
  !    a) the order in which puffs are time stepped must not be disrupted by altering the puff-tree
  !    b) the particle/puff may need to contribute to the output
  !   $$ other uses (particles which don't contrib to output) may be possible however and beneficial
  !      to memory needs
  ! 2) Particle/puff has Marked set to 1 or 2. It then does nothing (no more time-stepping) except
  !    possibly contribute to output, followed by being inactivated completely as in (1)
  ! 3) (puffs only). Puff has Active set to false due to being replaced by two split puffs. It remains
  !    in the puff tree and outside the free puff stack.
  ! Note particle/puffs in the free particle/puff stack can be reused for release, particle splitting,
  ! and puff splitting. Hence once a particle/puff is inactivated completely as in (1) one cannot
  ! assume it will survive till the end of the sync step.

  Integer          :: iHCoord           !
  Integer          :: iZCoord           !
  Integer          :: iHCoordOld        ! Coords for XOld
  Integer          :: iZCoordOld        !
  Integer          :: iSource           ! Index of the source the particle was released from.
  Real(Std)        :: H                 ! Boundary layer depth (in metres).
  Real(Std)        :: WSed              ! Gravitational settling velocity (for particulates only).
  Real(Std)        :: Diameter          ! Diameter (for particulates only).
  Real(Std)        :: Density           ! Density (for particulates only).
  Real(Std)        :: ParticleShape     ! Particle Shape (for particulates only).
  Integer          :: ShapeSchemeCode   ! Shape Scheme (for particulates only). 
End Type Particle_

!-------------------------------------------------------------------------------------------------------------

Type :: Extra_ ! A set of extra information for a particle/puff.

! Note that VelMem false  => Skew false, U = 0, SigUU unpredictable
!           MVelMem false => UM = 0
!           FMass = 0     => FM, FH, FMass0, B0Old = 0, and U/T/Theta/PAOld unpredictable
! (its easiest to sometimes set SigUU and U/T/Theta/PAOld, but if Extra(0) is used these may be incorrect)
! VelMem false       } allows Extra(0) to be used.
! MVelMem false      }
! Inhomog false      }
! FMass = 0          }
! TPlus = TMinus = 0 }
! TInst = -infinity  }

  Real(Std)        :: U(3)      ! Particle turbulent velocity.
  Real(Std)        :: UM(2)     ! Particle velocity for unresolved mesoscale motions.

  Real(Std)        :: SigUU(3)  ! Velocity variances at particle position at previous time step.

  Type(ShortTime_) :: TPlus     !} Time of leading and trailing 'time-edges' of particle relative to the
  Type(ShortTime_) :: TMinus    !} nominal time. (TMinus is <=0).
  Type(ShortTime_) :: TInst     ! Output travel time beyond which no time spread is assumed.

  Logical          :: Skew      ! Indicates particle uses skew model
  Logical          :: VelMem    ! Indicates particle uses (x,u)-Markov model
  Logical          :: Inhomog   !
  Logical          :: MVelMem   ! Indicates particle uses (x,u)-Markov model for unresolved mesoscale motions.

  Real(Std)        :: FMass     !} Fluxes of mass, momentum, heat, and mass with entrainment due to plume rise
  Real(Std)        :: FM(3)     !} only. However, after a particle is initialised by InitParticle and before
  Real(Std)        :: FH        !} InitVelocity is called, these hold the volume flow rate, flow velocity,
  Real(Std)        :: FMass0    !} temperature and volume flow rate of the source
  Real(Std)        :: B0Old     ! Plume radius from entrainment due to plume rise only at previous time.
  Real(Std)        :: UAOLD(3)  ! Ambient velocity at particle position at previous time step.
  Real(Std)        :: TAOLD     ! Ambient temperature at particle position at previous time step.
  Real(Std)        :: ThetaAOLD ! Ambient potential temperature at particle position at previous timestep.
  Real(Std)        :: PAOLD     ! Ambient pressure at particle position at previous time step.
End Type Extra_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function ParticleActive(Particle)
! Indicates if particle is active.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In) :: Particle ! Particle.
  ! Function result:
  Logical :: ParticleActive ! Indicates if particle is active.

  ParticleActive = Particle%Active

End Function ParticleActive

!-------------------------------------------------------------------------------------------------------------

Function ParticleReallyActive(Particle)
! Indicates if particle is active.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In) :: Particle ! Particle.
  ! Function result:
  Logical :: ParticleReallyActive ! Indicates if particle is active and not 'marked'.

  ParticleReallyActive = Particle%Active .and. Particle%Marked == 0

End Function ParticleReallyActive

!-------------------------------------------------------------------------------------------------------------

Function UseParticleForOutput(Particle)
! .

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In) :: Particle ! Particle.
  ! Function result:
  Logical :: UseParticleForOutput ! .

  UseParticleForOutput = Particle%Active .and. Particle%Marked <= 1

End Function UseParticleForOutput

!-------------------------------------------------------------------------------------------------------------

Subroutine MarkParticle(UseForOutput, Particle)
! Marks a particle or particle part of a puff to be inactivated at the end of the sync time.

  Implicit None
  ! Argument list:
  Logical,         Intent(In)    :: UseForOutput ! Indicates particle/puff is to be used
                                                 ! in output calculations.
  Type(Particle_), Intent(InOut) :: Particle     ! Particle to be killed.


  ! Mark particle.
  If (UseForOutput) Then
    Particle%Marked = 1
  Else
    Particle%Marked = 2
  End If

End Subroutine MarkParticle

!-------------------------------------------------------------------------------------------------------------

Function ParticleTime(Particle)
! Gives time of particle.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In) :: Particle ! Particle.
  ! Function result:
  Type(ShortTime_) :: ParticleTime ! Time of particle.

  ParticleTime = Particle%T0 + Particle%T

End Function ParticleTime

!-------------------------------------------------------------------------------------------------------------

Function OldParticleTime(Particle)
! Gives time of particle at previous time step.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In) :: Particle ! Particle.
  ! Function result:
  Type(ShortTime_) :: OldParticleTime ! Time of particle.

  OldParticleTime = Particle%T0 + Particle%TOld

End Function OldParticleTime

!-------------------------------------------------------------------------------------------------------------

Function ParticleMaxTime(Particle, Extra)
! Gives max time of particle, accounting for time spread.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In) :: Particle ! Particle.
  Type(Extra_),    Intent(In) :: Extra    ! .
  ! Function result:
  Type(ShortTime_) :: ParticleMaxTime ! Time of particle.

  ParticleMaxTime = Particle%T0 + Particle%T + Extra%TPlus

End Function ParticleMaxTime

!-------------------------------------------------------------------------------------------------------------

Function ParticleMinTime(Particle, Extra)
! Gives min time of particle, accounting for time spread.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In) :: Particle ! Particle.
  Type(Extra_),    Intent(In) :: Extra    ! .
  ! Function result:
  Type(ShortTime_) :: ParticleMinTime ! Time of particle.

  ParticleMinTime = Particle%T0 + Particle%T + Extra%TMinus

End Function ParticleMinTime

!-------------------------------------------------------------------------------------------------------------

Function OldParticleMinTime(Particle, Extra)
! Gives min time of particle at previous time step, accounting for time spread.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In) :: Particle ! Particle.
  Type(Extra_),    Intent(In) :: Extra    ! .
  ! Function result:
  Type(ShortTime_) :: OldParticleMinTime ! Time of particle.

  OldParticleMinTime = Particle%T0 + Particle%TOld + Extra%TMinus

End Function OldParticleMinTime

!-------------------------------------------------------------------------------------------------------------

Function ParticleInstTime(Particle, Extra)
! Gives output time of particle at and beyond which time spread is ignored.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In) :: Particle ! Particle.
  Type(Extra_),    Intent(In) :: Extra    ! .
  ! Function result:
  Type(ShortTime_) :: ParticleInstTime ! Time of particle.

  ParticleInstTime = Particle%T0 + Extra%TInst

End Function ParticleInstTime

!-------------------------------------------------------------------------------------------------------------

Function TravelTime(Particle)
! Gives travel time of particle.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In) :: Particle ! Particle.
  ! Function result:
  Type(ShortTime_) :: TravelTime ! Travel time of particle.

  TravelTime = Particle%T

End Function TravelTime

!-------------------------------------------------------------------------------------------------------------

Function OldTravelTime(Particle)
! Gives travel time of particle at previous time step.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In) :: Particle ! Particle.
  ! Function result:
  Type(ShortTime_) :: OldTravelTime ! Old travel time of particle.

  OldTravelTime = Particle%TOld

End Function OldTravelTime

!-------------------------------------------------------------------------------------------------------------

Subroutine InactivateParticle(                        &
             iP,                                      &
!! AJM
             Coords,DispState_iHCoordLatLong,         &
!! AJM
             Particles,                               &
             nParticles, FreeParticleStack,           &
             nParticleExtras, FreeParticleExtraStack, &
             iParticleExtras                          &
           )
! Kills a single particle and returns it to the free particle stack.

  Implicit None
  ! Argument list:
  Integer,         Intent(In)    :: iP                        ! Index of particle to be killed.
  Type(Particle_), Intent(InOut) :: Particles(:)              !

!!!! AJM
  Type(Coords_),   Intent(In)     :: Coords
  Type(Position_)                 :: Position    !
  Integer,         Intent(In)     :: DispState_iHCoordLatLong
!!!! AJM

  Integer,         Intent(InOut) :: nParticles                !
  Integer,         Intent(InOut) :: FreeParticleStack(:)      !
  Integer,         Intent(InOut) :: nParticleExtras           !
  Integer,         Intent(InOut) :: FreeParticleExtraStack(:) !
  Integer,         Intent(InOut) :: iParticleExtras(:)        !

!!!! AJM
!  Real(Std)        :: x         ! Particle position.
!  Real(Std)        :: y         ! Particle position.
!  Real(Std)        :: z         ! Particle position.
!  Type(ShortTime_) :: t         ! Travel time.
!  Character*12 xstr,ystr,zstr,hCoordstr,zCoordstr
!  Integer          :: iHCoord      !
!  Integer          :: iZCoord      !

  Real(Std)                 :: XLatLong(3) ! Position in Lat Long coord system
  Real(Std) :: IPartAge,IPartTime

!!!! AJM

  
  ! Kill particle.
  Particles(iP)%Active = .false.


!!!! AJM !!!

!  Call Message('De-activate a particle', 0)
!  x=Particles(iP)%X(1)
!  y=Particles(iP)%X(2)
!  z=Particles(iP)%X(3)
!  t=Particles(iP)%T
!  iHCoord=Particles(iP)%iHCoord
!  iZCoord=Particles(iP)%iZCoord
!  xstr=Std2Char(x,12)
!  ystr=Std2Char(y,12)
!  zstr=Std2Char(z,12)
!  hCoordstr=Int2Char(iHCoord,12)
!  zCoordstr=Int2Char(iZCoord,12)
!  Call Message(xstr//' '//ystr//' '//zstr//' '//hCoordstr//' '//zCoordstr, 0)

! Calculate particle position in lat long
  Position = X2Position(Coords, Particles(iP)%X, Particles(iP)%iHCoord, Particles(iP)%iZCoord)

! Changes Position to Position in LatLong
  Call ConvertToH(Coords, DispState_iHCoordLatLong, Position)

! Changes from Position type to X(1:2)
  XLatLong = Position2X(Coords, Position, DispState_iHCoordLatLong, Particles(iP)%iZCoord)

! Calculate Travel Time of Particle in Seconds
  IPartAge = ShortTime2RealTime(TravelTime(Particles(iP)))

!!!!
! This does not work - maybe negative time because backwards?
! IPartTime = ShortTime2RealTime(ParticleTime(Particles(iP)))
!!!!

  WRITE(ParticleEndLocationFileUnit,'(2F7.2,I6,F8.2,I3)')XLatLong(1),XLatLong(2),NINT(XLatLong(3)),&
       IPartAge/3600.,Particles(iP)%iSource

!!!! AJM !!!

  ! Return particle to the stack of free particles and reset number of active
  ! particles.
  FreeParticleStack(nParticles) = iP
  nParticles = nParticles - 1

  ! Return extra to the stack of free extras and reset number of active extras
  ! (except for extra zero which is shared by all cheap particles).
  If (iParticleExtras(iP) > 0) Then
    FreeParticleExtraStack(nParticleExtras) = iParticleExtras(iP)
    nParticleExtras = nParticleExtras - 1
  End If

End Subroutine InactivateParticle

!-------------------------------------------------------------------------------------------------------------

Subroutine InitParticle(                                     &
             iUP,                                            &
             Time,                                           &
             SkewTime, VelMemTime, InhomogTime, MVelMemTime, &
             Coords,                                         &
             Specieses, SizeDists,                           &
             iSource, Source, H3,                            &
             Diameter, Density, ParticleShape, X, dX,        &
             Particle, Extra, iP                             &
           )
! Initialises a particle.

  Implicit None
  ! Argument list:
  Integer,          Intent(In)  :: iUP           !
  Type(ShortTime_), Intent(In)  :: Time          !
  Type(ShortTime_), Intent(In)  :: SkewTime      !
  Type(ShortTime_), Intent(In)  :: VelMemTime    !
  Type(ShortTime_), Intent(In)  :: InhomogTime   !
  Type(ShortTime_), Intent(In)  :: MVelMemTime   !
  Type(Coords_),    Intent(In)  :: Coords        !
  Type(Specieses_), Intent(In)  :: Specieses     !
  Type(SizeDists_), Intent(In)  :: SizeDists     ! $$ needed?
  Integer,          Intent(In)  :: iSource       !
  Type(Source_),    Intent(In)  :: Source        !
  Real(Std),        Intent(In)  :: H3            !
  Real(Std),        Intent(In)  :: Diameter      !
  Real(Std),        Intent(In)  :: ParticleShape !
  Real(Std),        Intent(In)  :: Density       !
  Real(Std),        Intent(In)  :: X(3)          !
  Real(Std),        Intent(In)  :: dX(3)         !
  Type(Particle_),  Intent(Out) :: Particle      !
  Type(Extra_),     Intent(Out) :: Extra         !
  Integer,          Intent(In)  :: iP            ! Particle index
  ! Locals:
  Real(Std) :: R               !
  Real(Std) :: RX              !
  Real(Std) :: RY              !
  Real(Std) :: RZ              !
  Real(Std) :: HMax            ! Max metric coefficient.
  Real(Std) :: H1              !
  Real(Std) :: H2              !
  Real(Std) :: LocalH1         !
  Real(Std) :: LocalH2         !
  Real(Std) :: Factor          !
  Integer   :: SourceDim       ! Number of elements of dX(:) which are non-zero.
  Real(Std) :: GeometricFactor !
  Integer   :: i               !
  Real(Std) :: PX1             !
  Real(Std) :: PX2             !

  Particle%iUP = iUP

  Particle%iSource = iSource

  Extra%Skew    = SkewTime    > ZeroShortTime()
  Extra%VelMem  = VelMemTime  > ZeroShortTime()
  Extra%Inhomog = InhomogTime > ZeroShortTime()
  Extra%MVelMem = MVelMemTime > ZeroShortTime()

 ! Particle%iHCoord = Source%iHCoord
 ! Particle%iZCoord = Source%iZCoord
    Particle%iHCoord  = FindHCoordIndex(Source%HCoordName,  Coords)
    Particle%iZCoord  = FindZCoordIndex(Source%ZCoordName,  Coords)

  Call MetricCoeffs(                       &
         Coords%HCoords(Particle%iHCoord), &
         X(1:2),                           &
         HMax, H1, H2                      &
       )

  ! Move source shape bits to Sources.F90? $$
  ! Use dHCoordName? $$

  Do

    ! Ellipsoid sources.
    If ( Source%Shape .CIEq. 'Ellipsoid' ) Then

      ! Generate uniform distribution on dX aligned coord system.
      If (Source%TopHat) Then

        Do
          If (dX(1) == 0.0) Then
            RX = 0.0
          Else
            Call GetRandomNumber(RX, iP)
            RX = RX - 0.5
          End If
          If (dX(2) == 0.0) Then
            RY = 0.0
          Else
            Call GetRandomNumber(RY, iP)
            RY = RY - 0.5
          End If
          If (dX(3) == 0.0) Then
            RZ = 0.0
          Else
            Call GetRandomNumber(RZ, iP)
            RZ = RZ - 0.5
          End If
          If (RX**2 + RY**2 + RZ**2 <= 0.25) Exit
        End Do
        Particle%X(1) = dX(1) * RX
        Particle%X(2) = dX(2) * RY
        Particle%X(3) = dX(3) * RZ

      ! Generate Gaussian distribution on dX aligned coord system.
      Else

        SourceDim = 0
        Do i = 1, 3
          If (dX(i) > 0.0) SourceDim = SourceDim + 1
        End Do
        GeometricFactor = Sqrt( 1.0 / ( 4.0 * (Real(SourceDim) + 2.0) ) )
        R = Gauss(iP)
        Particle%X(1) = GeometricFactor * dX(1) * R
        R = Gauss(iP)
        Particle%X(2) = GeometricFactor * dX(2) * R
        R = Gauss(iP)
        Particle%X(3) = GeometricFactor * dX(3) * R

      End If

    ! Cylindroid sources.
    Else If ( Source%Shape .CIEq. 'Cylindroid' ) Then

      ! NB: the height of the cylindroid is aligned with the vertical

      ! Generate uniform distribution on dX aligned coord system.
      If (Source%TopHat) Then

        Do
          If (dX(1) == 0.0) Then
            RX = 0.0
          Else
            Call GetRandomNumber(RX, iP)
            RX = RX - 0.5
          End If
          If (dX(2) == 0.0) Then
            RY = 0.0
          Else
            Call GetRandomNumber(RY, iP)
            RY = RY - 0.5
          End If
          If (RX**2 + RY**2 <= 0.25) Exit
        End Do
        If (dX(3) == 0.0) Then
          RZ = 0.0
        Else
          Call GetRandomNumber(RZ, iP)
          RZ = RZ - 0.5
        End If
        Particle%X(1) = dX(1) * RX  
        Particle%X(2) = dX(2) * RY
        Particle%X(3) = dX(3) * RZ          

      ! Generate Gaussian distribution on dX aligned coord system.
      Else

        SourceDim = 0
        Do i = 1, 2
          If (dX(i) > 0.0) SourceDim = SourceDim + 1
        End Do
        GeometricFactor = Sqrt( 1.0 / ( 4.0 * (Real(SourceDim) + 2.0) ) )
        R = Gauss(iP)
        Particle%X(1) = GeometricFactor * dX(1) * R
        R = Gauss(iP)
        Particle%X(2) = GeometricFactor * dX(2) * R
        GeometricFactor = Sqrt(1.0/12.0)
        R = Gauss(iP)
        Particle%X(3) = GeometricFactor * dX(3) * R   

      End If

    ! Cuboid sources.
    Else If ( Source%Shape .CIEq. 'Cuboid' ) Then

      ! Generate uniform distribution on dX aligned coord system.
      If (Source%TopHat) Then

        Call GetRandomNumber(R, iP)
        Particle%X(1) = dX(1) * (R - 0.5)
        Call GetRandomNumber(R, iP)
        Particle%X(2) = dX(2) * (R - 0.5)
        Call GetRandomNumber(R, iP)
        Particle%X(3) = dX(3) * (R - 0.5)

      ! Generate Gaussian distribution on dX aligned coord system.
      Else

        R = Gauss(iP)
        Particle%X(1) = Sqrt(1.0/12.0) * dX(1) * R
        R = Gauss(iP)
        Particle%X(2) = Sqrt(1.0/12.0) * dX(2) * R
        R = Gauss(iP)
        Particle%X(3) = Sqrt(1.0/12.0) * dX(3) * R

      End If

    Else
      Call Message('UNEXPECTED FATAL ERROR in InitParticle regarding source shape.', 4)
    End If

    ! Rotate particle to true coordinate system and apply metric coefficients.
    PX1 = Particle%X(1)
    PX2 = Particle%X(2)
    Particle%X(1) = PX1 * Cos(Source%Angle) - PX2 * Sin(Source%Angle)
    Particle%X(2) = PX2 * Cos(Source%Angle) + PX1 * Sin(Source%Angle)
    If (Source%dHMetres) Then
      Particle%X(1) = Particle%X(1) / H1
      Particle%X(2) = Particle%X(2) / H2
    End If
    If (Source%dZMetres) Then
      Particle%X(3) = Particle%X(3) / H3
    End If
    Particle%X(1) = X(1) + Particle%X(1)
    Particle%X(2) = X(2) + Particle%X(2)
    Particle%X(3) = X(3) + Particle%X(3)

    Particle%XOld       = Particle%X
    Particle%iHCoordOld = Particle%iHCoord
    Particle%iZCoordOld = Particle%iZCoord

    If (Source%UniformArea) Then
      Call MetricCoeffs(                       &
             Coords%HCoords(Particle%iHCoord), &
             Particle%X(1:2),                  &
             HMax, LocalH1, LocalH2            &
           )
      Call GetRandomNumber(R, iP)
      If (R <= LocalH1 * LocalH2 / HMax**2) Exit
    Else
      Exit
    End If

  End Do

  Particle%T          = ZeroShortTime()
  Particle%TOld       = Particle%T
  Particle%T0         = Time
  Particle%Active     = .true.
  Particle%Marked = 0
  Particle%NeedToSetVel = .false.

  ! Particles currently always instantaneous.
  Extra%TPlus  = ZeroShortTime()
  Extra%TMinus = ZeroShortTime()
  Extra%TInst  = InfPastShortTime(Interval = .true.)

  If (Source%PlumeRise) Then
    Extra%FMass  = Source%VolumeFlowRate
    Extra%FM(1)  = 0.0
    Extra%FM(2)  = 0.0
    Extra%FM(3)  = Source%FlowVelocity
    Extra%FH     = Source%Temperature
    Extra%FMass0 = Source%VolumeFlowRate
    Extra%B0Old  = dX(1)/2.0  ! If dX not metres need correction here $$ (puff too)
  Else
    Extra%FMass  = 0.0
  End If

  ! Initialise gravitational settling velocity for solid particle (from tables in TDN
  ! 262)
  If (Source%Particulate) Then

    If (.false.) Then ! $$ del in due course
    ! $$ Above diameter = 100? Better to have cts formula. Use density.
    ! $$ Apply to puffs too.
    If (Diameter <= 0.1) Then
      Particle%WSed = 0.0
    Else If ( Diameter > 0.1 .And. Diameter <= 0.3 ) Then
      Particle%WSed = 2.076E-6
    Else If ( Diameter > 0.3 .And. Diameter <= 1.0 ) Then
      Particle%WSed = 2.076E-5
    Else If ( Diameter > 1.0 .And. Diameter <= 3.0 ) Then
      Particle%WSed = 2.076E-4
    Else If ( Diameter > 3.0 .And. Diameter <= 10.0 ) Then
      Particle%WSed = 2.076E-3
    Else If ( Diameter > 10.0 .And. Diameter <= 30.0 ) Then
      Particle%WSed = 2.076E-2
    Else If ( Diameter > 30.0 .And. Diameter <= 100.0 ) Then
      Particle%WSed = 1.908E-1
    End If
    End If

    Particle%Diameter = Diameter
    Particle%Density  = Density
    Particle%ParticleShape = ParticleShape
    Particle%ShapeSchemeCode = Source%ShapeSchemeCode


  Else

    If (.false.) Then ! $$ del in due course
      Particle%WSed = 0.0
    End If
    Particle%Diameter = 0.0
    Particle%Density  = 0.0

  End If

End Subroutine InitParticle

!-------------------------------------------------------------------------------------------------------------

Subroutine InitParticle2(Flow, Particle, Extra, iP)
! .

  Implicit None
  ! Argument list:
  Type(Flow_),     Intent(In)    :: Flow     !
  Type(Particle_), Intent(InOut) :: Particle !
  Type(Extra_),    Intent(InOut) :: Extra    !
  Integer,         Intent(In)    :: iP       ! Particle index
  ! Locals:
  Real(Std) :: R       !
  Real(Std) :: Rho_p   !
  Real(Std) :: Theta_p !

  If (Extra%VelMem) Then

    R = Gauss(iP)
    Extra%U(1) = Sqrt(Flow%SigUU(1)) * R
    R = Gauss(iP)
    Extra%U(2) = Sqrt(Flow%SigUU(2)) * R
    Extra%U(3) = InitW(Flow, 1.0, Extra%Skew, iP)

    Extra%SigUU(1) = Flow%SigUU(1)
    Extra%SigUU(2) = Flow%SigUU(2)
    Extra%SigUU(3) = Flow%SigUU(3)
    Extra%UAOLD(1) = Flow%U(1)
    Extra%UAOLD(2) = Flow%U(2)
    Extra%UAOLD(3) = Flow%U(3)
    Extra%TAOLD    = Flow%T
    Extra%ThetaAOLD= Flow%Theta
    Extra%PAOLD    = Flow%P

  End If

  If (Extra%MVelMem) Then
    Extra%UM(1) = Sqrt(Flow%SigUUM) * Gauss(iP)
    Extra%UM(2) = Sqrt(Flow%SigUUM) * Gauss(iP)
  Else
    Extra%UM(:) = 0.0
  End If

  If (Extra%FMass /= 0.0) Then ! test on real dangerous - hence else block
    ! Plume density and potential temperature.
    Rho_p   = Flow%P/(GasConstant * Extra%FH)
    Theta_p = Extra%FH * (PRef/Flow%P)**(GasConstant/Cp)
    ! Plume fluxes.
    Extra%FMass  = Rho_p * Extra%FMass
    Extra%FM(1)  = (Extra%FM(1) - Flow%U(1)) * Extra%FMass
    Extra%FM(2)  = (Extra%FM(2) - Flow%U(2)) * Extra%FMass
    Extra%FM(3)  = (Extra%FM(3) - Flow%U(3)) * Extra%FMass
    Extra%FH     = Cp * (Theta_p - Flow%Theta) * Extra%FMass
    Extra%FMass0 = Rho_p * Extra%FMass0
  Else
    Extra%FMass  = 0.0
    Extra%FM(1)  = 0.0
    Extra%FM(2)  = 0.0
    Extra%FM(3)  = 0.0
    Extra%FH     = 0.0
    Extra%FMass0 = 0.0
  End If

  Particle%H = Flow%H

End Subroutine InitParticle2

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcRandomWalkDt(Flow, Particle, Extra, Dt)
! Code to compute random walk time step.

  Implicit None
  ! Argument list:
  Type(Flow_),     Intent(In)  :: Flow     !
  Type(Particle_), Intent(In)  :: Particle !
  Type(Extra_),    Intent(In)  :: Extra    !
  Real(Std),       Intent(Out) :: Dt       !
  ! Locals:
  Integer :: i

  ! Start with a huge value and then impose various upper limits. The huge value must
  ! not be so huge that it can't be converted to a short time.

  Dt = Float(Huge(i))*50.0

  If (Extra%VelMem) Then

!    If (Flow%MaxDZ /= Huge(Flow%MaxDZ)) Dt = Min(Dt, Flow%MaxDZ / Sqrt(Flow%SigUU(3)))
    If (Flow%MaxDZ /= Huge(Flow%MaxDZ)) Then
      Dt = Min(Dt, Flow%MaxDZ / Max(Abs(Extra%U(3)),Sqrt(Flow%SigUU(3))))
    End If
    Dt = Min(Dt, Flow%TauUU(3))
    Dt = 0.1_Std*Dt

  Else

    If (Extra%Inhomog) Then
      If (Flow%MaxDZ /= Huge(Flow%MaxDZ)) Dt = Min(Dt, Flow%MaxDZ**2 / (2.0 * Flow%K(3)))
      Dt = 0.03_Std*Dt
    Else
      If (Flow%MaxDZ /= Huge(Flow%MaxDZ)) Dt = Min(Dt, Flow%MaxDZ**2 / (2.0 * Flow%HK(3)))
      Dt = 0.03_Std*Dt
    End If

  End If

  If (Extra%FMass > 0.0) Then
    If (Flow%MaxDZ /= Huge(Flow%MaxDZ)) Then
      Dt = Min(Dt, 0.03_Std * Flow%MaxDZ / Abs(Flow%U(3) + Extra%FM(3)/Extra%FMass))
      ! $$ trap possible zero denominator - +0 probably OK - presumably Abs(x) can't be -0.
    End If
!   Dt = 0.03_Std*Dt
  End If

! Imposes a limit on Dt due to the type of met - currently only UKV set (to 60s)
  Dt = Min(Dt,Flow%MaxDT)

End Subroutine CalcRandomWalkDt

!-------------------------------------------------------------------------------------------------------------

Subroutine EvolveParticle(                                                                  &
             Flow, Dt, H1, H2, Turbulence, MesoscaleMotions, Damping, VerticalVelocity, Zs, &
             HCoordParticle,  WMeanAndTurb, Particle, Extra, iP                             &
           )
! Time-step Particles.

  Implicit None
  ! Argument list:
  Type(Flow_),     Intent(In)    :: Flow             !
  Real(Std),       Intent(In)    :: Dt               !
  Real(Std),       Intent(In)    :: H1               ! Use for XY
  Real(Std),       Intent(In)    :: H2               ! use for xy
  Real(Std),       Intent(In)    :: Zs               ! Height below which dry deposit
  Real(Std),       Intent(Out)   :: WMeanAndTurb     ! Vertical velocity excluding sedimentation
  Logical,         Intent(In)    :: Turbulence       !
  Logical,         Intent(In)    :: MesoscaleMotions !
  Logical,         Intent(In)    :: Damping          !
  Logical,         Intent(In)    :: VerticalVelocity !
  Type(Particle_), Intent(InOut) :: Particle         !
  Type(Extra_),    Intent(InOut) :: Extra            !
  Type(HCoord_),   Intent(In)    :: HCoordParticle 
  Integer,         Intent(In)    :: iP               ! Particle index
  ! Locals:
  Real(Std) :: A(3)      !
  Real(Std) :: B(3)      !
  Real(Std) :: U(3)      ! Velocity due to mean flow and mean and random plume rise.
  Real(Std) :: R         !
  Real(Std) :: UP        ! U VELOCITY OF THE PARTICLE/PLUME
  Real(Std) :: VP        ! V VELOCITY OF THE PARTICLE/PLUME
  Real(Std) :: WP        ! W VELOCITY OF THE PARTICLE/PLUME
  Real(Std) :: RANDPLUMX ! RANDOM PLUME INDUCED TURBULENT COMPONENT
  Real(Std) :: RANDPLUMY ! RANDOM PLUME INDUCED TURBULENT COMPONENT
  Real(Std) :: RANDPLUMZ ! RANDOM PLUME INDUCED TURBULENT COMPONENT
  Real(Std) :: RANDPLUMSIG ! Standard deviation of random displacement in
                           ! each coord direction due to plume induced
                           ! turbulence.
  Real(Std) :: IPartAge
  Real(Std) :: KFactor(3)
  Real(Std) :: KFactorM
  Real(Std) :: TravelTime
  Real(Std) :: ZTemp
  Real(Std) :: DtLeft
  Real(Std) :: dX(2)

  Real ZGradientDamping

  dX(:) = 0.0
  Particle%XOld       = Particle%X
  Particle%iHCoordOld = Particle%iHCoord
  Particle%iZCoordOld = Particle%iZCoord

  If (Extra%FMass > 0.0) Then ! move out to Case.f90 and remove use stmt?
    IPartAge = ShortTime2RealTime(Particle%T)
    If (IPartAge /= 0.0) Then
      If (Extra%Inhomog) Then
        Call PlumeRise(                                                   &
             Flow%U(1), Flow%U(2), Flow%U(3), Flow%T, Flow%Theta, Flow%P, &
             Extra%UAOLD(1), Extra%UAOLD(2), Extra%UAOLD(3),              &
             Extra%TAOLD, Extra%ThetaAOLD, Extra%PAOLD,                   &
             Flow%TauUU(3), Sqrt(Flow%SigUU(3)), Flow%Eps,                &
             IPartAge, dT, .false.,                                       &
             UP, VP, WP, RANDPLUMSig,                                     &
             Extra%FMass, Extra%FM(1), Extra%FM(2), Extra%FM(3),          &
             Extra%FH, Extra%FMASS0, Extra%B0OLD                          &
           )
      Else
        Call PlumeRise(                                                   &
             Flow%U(1), Flow%U(2), Flow%U(3), Flow%T, Flow%Theta, Flow%P, &
             Extra%UAOLD(1), Extra%UAOLD(2), Extra%UAOLD(3),              &
             Extra%TAOLD, Extra%ThetaAOLD, Extra%PAOLD,                   &
             Flow%HTauUU(3), Sqrt(Flow%HSigUU(3)), Flow%HEps,             &
             IPartAge, dT, .false.,                                       &
             UP, VP, WP, RANDPLUMSig,                                     &
             Extra%FMass, Extra%FM(1), Extra%FM(2), Extra%FM(3),          &
             Extra%FH, Extra%FMASS0, Extra%B0OLD                          &
           )
    End If
  Else
      UP=0.0
      VP=0.0
      WP = Flow%U(3) + Extra%FM(3)/Extra%FMass
    End If
    U(1) = UP + RandPlumSig*Gauss(iP)/Dt
    U(2) = VP + RandPlumSig*Gauss(iP)/Dt
    U(3) = WP + RandPlumSig*Gauss(iP)/Dt
  Else
    U(1) = Flow%U(1)
    U(2) = Flow%U(2)
    U(3) = Flow%U(3)
  End If

  If (Turbulence) Then

    If (Extra%VelMem) Then

  !  If (Particle%X(3) * Particle%X(3) < 9.0 * 10000.0) Then ! Note 9 hardwired $$
  !    ZGradientDamping = Erf(Particle%X(3) / Sqrt(2.0 * 10000.0))
  !  Else
      ZGradientDamping = 1.0
  !  End If
      Call XUSdeTerms(Particle, Extra, Flow, Dt, U(3), A, B, Zs, ZGradientDamping, 1.0_Std)

      Extra%U(1) = Extra%U(1) + A(1) + B(1) * Gauss(iP)
      Extra%U(2) = Extra%U(2) + A(2) + B(2) * Gauss(iP)
      Extra%U(3) = Extra%U(3) + A(3) + B(3) * Gauss(iP)
!     Extra%U(3) = Extra%U(3) * sqrt(Flow%SigUU(3)) / sqrt(Extra%SigUU(3)) + A(3) + B(3) * Gauss(iP)
!     Extra%U(3) = A(3) + B(3) * Gauss(iP)

      dX(1) = dX(1) + (U(1) + Extra%U(1)) * Dt 
      dX(2) = dX(2) + (U(2) + Extra%U(2)) * Dt 
      If (VerticalVelocity) Then
        WMeanAndTurb = U(3) + Extra%U(3)
        If (Particle%WSed > 0.0 ) Then
          If (Particle%X(3) < Zs) Then
            If (Particle%WSed * Dt / Zs < 0.01) Then
              ZTemp = Particle%X(3) + (WMeanAndTurb - Particle%WSed * Particle%X(3) / Zs) * Dt
            Else
              ZTemp = WMeanAndTurb * Zs / Particle%WSed +            &
                (Particle%X(3) - WMeanAndTurb * Zs / Particle%WSed) * Exp(-Particle%WSed * Dt / Zs)
            End If
            If (ZTemp > Zs) Then
              If (Particle%WSed *(Particle%X(3) - Zs) / (Particle%X(3) * Particle%WSed &
                  - WMeanAndTurb * Zs) < 0.01) Then
                DtLeft = Dt + Zs * (Zs - Particle%X(3)) / (Particle%X(3) * Particle%WSed - WMeanAndTurb * Zs)
              Else
                DtLeft = Dt + Zs / Particle%WSed * Log(Zs * (Particle%WSed - WMeanAndTurb) /    &
                  (Particle%X(3) * Particle%WSed - WMeanAndTurb * Zs))
              End If
              ZTemp = Zs + (WMeanAndTurb - Particle%WSed) * DtLeft
            End If
          Else
            ZTemp = Particle%X(3) + (WMeanAndTurb - Particle%WSed) * Dt
            If (ZTemp < Zs) Then
              DtLeft = Dt - (Zs - Particle%X(3)) / (WMeanAndTurb - Particle%WSed)
              If (Particle%WSed * DtLeft / Zs > 0.01) Then
                ZTemp = WMeanAndTurb * Zs / Particle%WSed +                                      &
                  (Zs  - WMeanAndTurb * Zs / Particle%WSed) * Exp(-Particle%WSed * DtLeft / Zs)
              End If
            End If
          End If
        Else
          ZTemp = Particle%X(3) + WMeanAndTurb * Dt
        End If

        Particle%X(3) = ZTemp

      End If

      Extra%SigUU(1) = Flow%SigUU(1)
      Extra%SigUU(2) = Flow%SigUU(2)
      Extra%SigUU(3) = Flow%SigUU(3)

    Else

      If (Extra%Inhomog) Then
        Call XInhomogSdeTerms(Particle, Flow, Dt, A, B, 1.0_Std)
      Else
        Call XHomogSdeTerms(Particle, Flow, Dt, A, B)
      End If

      ! Near source correction. Not applied to 3rd component because can upset
      ! profiles (probably because it changes dK/dZ at constant travel time if tau
      ! varies with height). $$
      ! Note the algebraic expression is significantly faster than the exp. Change
      ! unresolved mesoscale motions too? $$
      If (Damping) Then
        TravelTime = ShortTime2RealTime(Particle%T) + Dt/2.0
        If (Extra%Inhomog) Then
          KFactor(1:2) = TravelTime / Sqrt(TravelTime**2 + Flow%TauUU(1:2)**2)
 !        KFactor(1:2) = 1.0 - Exp(- TravelTime / Flow%TauUU(1:2))
        Else
          KFactor(1:2) = TravelTime / Sqrt(TravelTime**2 + Flow%HTauUU(1:2)**2)
 !        KFactor(1:2) = 1.0 - Exp(- TravelTime / Flow%HTauUU(1:2))
        End If
        A(1:2) = A(1:2) * KFactor(1:2)
        B(1:2) = B(1:2) * Sqrt(KFactor(1:2))
      End If

      dX(1) = dX(1) + (U(1) * Dt + A(1) + B(1) * Gauss(iP))
      dX(2) = dX(2) + (U(2) * Dt + A(2) + B(2) * Gauss(iP))
      If (VerticalVelocity) Then
        WMeanAndTurb = U(3) + (A(3) + B(3) * Gauss(iP)) / Dt
        If (Particle%WSed > 0.0 ) Then
          If (Particle%X(3) < Zs) Then
            If (Particle%WSed * Dt / Zs < 0.01) Then
              ZTemp = Particle%X(3) + (WMeanAndTurb - Particle%WSed * Particle%X(3) / Zs) * Dt
            Else
              ZTemp = WMeanAndTurb * Zs / Particle%WSed +            &
                (Particle%X(3) - WMeanAndTurb * Zs / Particle%WSed) * Exp(-Particle%WSed * Dt / Zs)
            End If
            If (ZTemp > Zs) Then
              If (Particle%WSed *(Particle%X(3) - Zs) / (Particle%X(3) * Particle%WSed &
                  - WMeanAndTurb * Zs) < 0.01) Then
                DtLeft = Dt + Zs * (Zs - Particle%X(3)) / (Particle%X(3) * Particle%WSed - WMeanAndTurb * Zs)
              Else
                DtLeft = Dt + Zs / Particle%WSed * Log(Zs * (Particle%WSed - WMeanAndTurb) /    &
                  (Particle%X(3) * Particle%WSed - WMeanAndTurb * Zs))
              End If
              ZTemp = Zs + (WMeanAndTurb - Particle%WSed) * DtLeft
            End If
          Else
            ZTemp = Particle%X(3) + (WMeanAndTurb - Particle%WSed) * Dt
            If (ZTemp < Zs) Then
              DtLeft = Dt - (Zs - Particle%X(3)) / (WMeanAndTurb - Particle%WSed)
              If (Particle%WSed * DtLeft / Zs > 0.01) Then
                ZTemp = WMeanAndTurb * Zs / Particle%WSed +                                      &
                  (Zs  - WMeanAndTurb * Zs / Particle%WSed) * Exp(-Particle%WSed * DtLeft / Zs)
              End If
            End If
          End If
        Else
          ZTemp = Particle%X(3) + WMeanAndTurb * Dt
        End If

        Particle%X(3) = ZTemp

      End If

    End If

  Else

    dX(1) = dX(1) + U(1) * Dt 
    dX(2) = dX(2) + U(2) * Dt 
    If (VerticalVelocity) Then
      WMeanAndTurb = U(3)
      If (Particle%WSed > 0.0 ) Then
        If (Particle%X(3) < Zs) Then
          If (Particle%WSed * Dt / Zs < 0.01) Then
            ZTemp = Particle%X(3) + (WMeanAndTurb - Particle%WSed * Particle%X(3) / Zs) * Dt
          Else
            ZTemp = WMeanAndTurb * Zs / Particle%WSed +            &
              (Particle%X(3) - WMeanAndTurb * Zs / Particle%WSed) * Exp(-Particle%WSed * Dt / Zs)
          End If
          If (ZTemp > Zs) Then
            If (Particle%WSed *(Particle%X(3) - Zs) / (Particle%X(3) * Particle%WSed &
               - WMeanAndTurb * Zs) < 0.01) Then
              DtLeft = Dt + Zs * (Zs - Particle%X(3)) / (Particle%X(3) * Particle%WSed - WMeanAndTurb * Zs)
            Else
              DtLeft = Dt + Zs / Particle%WSed * Log(Zs * (Particle%WSed - WMeanAndTurb) /    &
                (Particle%X(3) * Particle%WSed - WMeanAndTurb * Zs))
            End If
            ZTemp = Zs + (WMeanAndTurb - Particle%WSed) * DtLeft
          End If
        Else
          ZTemp = Particle%X(3) + (WMeanAndTurb - Particle%WSed) * Dt
          If (ZTemp < Zs) Then
            DtLeft = Dt - (Zs - Particle%X(3)) / (WMeanAndTurb - Particle%WSed)
            If (Particle%WSed * DtLeft / Zs > 0.01) Then
              ZTemp = WMeanAndTurb * Zs / Particle%WSed +                                      &
                (Zs  - WMeanAndTurb * Zs / Particle%WSed) * Exp(-Particle%WSed * DtLeft / Zs)
            End If
          End If
        End If
      Else
        ZTemp = Particle%X(3) + WMeanAndTurb * Dt
      End If

      Particle%X(3) = ZTemp

    End If

  End If

  If (MesoscaleMotions) Then

    If (Extra%MVelMem) Then

      Extra%UM(1) = Extra%UM(1) * (1.0 - Dt/Flow%TauUUM) +  &
                    Sqrt(2.0 * Dt * Flow%SigUUM / Flow%TauUUM) * Gauss(iP)
      Extra%UM(2) = Extra%UM(2) * (1.0 - Dt/Flow%TauUUM) +  &
                    Sqrt(2.0 * Dt * Flow%SigUUM / Flow%TauUUM) * Gauss(iP)
      dX(1) = dX(1) + Extra%UM(1) * Dt
      dX(2) = dX(2) + Extra%UM(2) * Dt

    Else

      If (Damping) Then
        KFactorM =                                                                 &
          Sqrt(                                                                    &
            2.0 * Flow%SigUUM * Flow%TauUUM * Dt *                                 &
            (1.0 - Exp(- (ShortTime2RealTime(Particle%T) + Dt/2.0) / Flow%TauUUM)) &
          )
      Else
        KFactorM = Sqrt(2.0 * Flow%SigUUM * Flow%TauUUM * Dt)
      End If
      dX(1) = dX(1) + KFactorM * Gauss(iP)
      dX(2) = dX(2) + KFactorM * Gauss(iP)

    End If

  End If

  If (HCoordParticle%CoordType == H_LatLong) Then
    Call TrajectoryCalc( HCoordParticle, Particle, dX)
  Else
    Particle%X(1) = Particle%X(1) + dX(1) / H1
    Particle%X(2) = Particle%X(2) + dX(2) / H2
  End If

  Extra%UAOLD(1) = Flow%U(1)
  Extra%UAOLD(2) = Flow%U(2)
  Extra%UAOLD(3) = Flow%U(3)
  Extra%TAOLD    = Flow%T
  Extra%ThetaAOLD= Flow%Theta
  Extra%PAOLD    = Flow%P

  Particle%H = Flow%H

End Subroutine EvolveParticle

!-------------------------------------------------------------------------------------------------------------

Subroutine TrajectoryCalc(HCoordParticle, Particle, dX)
! Computes trajectory more accurately near the pole when in a lat-long coordinate system.

  Implicit None
  ! Argument list:
  Type(HCoord_),   Intent(In)      :: HCoordParticle !
  Type(Particle_), Intent(InOut)   :: Particle       !
  Real(Std),       Intent(In)      :: dX(2)          !
  ! Locals
  Real(Std)  :: d                  ! Total distance moved by particle
  Real(Std)  :: theta              ! Direction of movement by particle
  Real(Std)  :: x, y, z            !
  Real(Std)  :: x2, y2, z2         !

  ! Put position in radians with latitude measured from equator
  Particle%X(1) = Particle%X(1) * HCoordParticle%Unit(1) + HCoordParticle%Origin(1)
  Particle%X(2) = Particle%X(2) * HCoordParticle%Unit(2) + HCoordParticle%Origin(2)

  ! Convert dX into a magnitude and angle (in radians)
  d = Sqrt(dX(1) ** 2 + dX(2) ** 2) / EarthRadius
  theta = ATan2ZeroTest( dX(1), dX(2) )

  ! Final position in rotated x, y, z coords
  x = - Sin(d) * Cos(theta)
  y =   Sin(d) * Sin(theta)
  z =   Cos(d)

  ! Final position in normal x, y, z coords
  x2 = x * Sin(Particle%X(2)) + z * Cos(Particle%X(2))
  y2 = y
  z2 = - x * Cos(Particle%X(2)) + z * Sin(Particle%X(2))

  ! Final position in standard lat-long coordinates
  Particle%X(2) = ATan2(z2, sqrt(x2 ** 2 + y2 ** 2))
  Particle%X(1) = Particle%X(1) + ATan2ZeroTest(y2, x2)

  ! Convert back to previous coordinate system
  Particle%X(1) = (Particle%X(1) - HCoordParticle%Origin(1)) / HCoordParticle%Unit(1)
  Particle%X(2) = (Particle%X(2) - HCoordParticle%Origin(2)) / HCoordParticle%Unit(2)

End Subroutine TrajectoryCalc

!-------------------------------------------------------------------------------------------------------------

Subroutine XHomogSdeTerms(Particle, Flow, Dt, A, B)
! Code to compute terms in the random walk stochastic differential equations.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In)  :: Particle !
  Type(Flow_),     Intent(In)  :: Flow     !
  Real(Std),       Intent(In)  :: Dt       !
  Real(Std),       Intent(Out) :: A(3)     !
  Real(Std),       Intent(Out) :: B(3)     !

  A(1:3) = 0.0

  B(1) = Sqrt(2.0 * Dt * Flow%HK(1))
  B(2) = Sqrt(2.0 * Dt * Flow%HK(2))
  B(3) = Sqrt(2.0 * Dt * Flow%HK(3))

End Subroutine XHomogSdeTerms

!-------------------------------------------------------------------------------------------------------------

Subroutine XInhomogSdeTerms(Particle, Flow, Dt, A, B, ZGradientDamping)
! Code to compute terms in the random walk stochastic differential equations.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In)  :: Particle         !
  Type(Flow_),     Intent(In)  :: Flow             !
  Real(Std),       Intent(In)  :: Dt               !
  Real(Std),       Intent(In)  :: ZGradientDamping !
  Real(Std),       Intent(Out) :: A(3)             !
  Real(Std),       Intent(Out) :: B(3)             !

  A(1) = 0.0
  A(2) = 0.0
  A(3) = Flow%dKdX(3) * Dt * ZGradientDamping

  B(1) = Sqrt(2.0 * Dt * Flow%K(1))
  B(2) = Sqrt(2.0 * Dt * Flow%K(2))
  B(3) = Sqrt(2.0 * Dt * Flow%K(3))

End Subroutine XInhomogSdeTerms

!-------------------------------------------------------------------------------------------------------------

Subroutine XUSdeTerms(Particle, Extra, Flow, Dt, Umean, A, B, Zs, ZGradientDamping, WVarianceDamping)
! Code to compute terms in the random walk stochastic differential equations.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In)  :: Particle         !
  Type(Extra_),    Intent(In)  :: Extra            !
  Type(Flow_),     Intent(In)  :: Flow             !
  Real(Std),       Intent(In)  :: Dt               !
  Real(Std),       Intent(In)  :: Umean            ! Mean vertical velocity
  Real(Std),       Intent(In)  :: ZGradientDamping !} Indicates results are to be calculated with vertical
  Real(Std),       Intent(In)  :: WVarianceDamping !} velocity variance and vertical gradients damped by these
                                                   !  factors.
  Real(Std),       Intent(In)  :: Zs               ! Height below which dry deposit
  Real(Std),       Intent(Out) :: A(3)             !
  Real(Std),       Intent(Out) :: B(3)             !
  ! Locals:
  Real(Std) :: SigWWDamped
  Real(Std) :: SkDamped
  Real(Std) :: dSigWWdZDamped
  Real(Std) :: dSkdZDamped
  Real(Std) :: Alpha
  Real(Std) :: Beta
  Real(Std) :: Gamma
  Real(Std) :: SigWA
  Real(Std) :: SigWB
  Real(Std) :: WMeanA
  Real(Std) :: WMeanB
  Real(Std) :: Az
  Real(Std) :: Bz
  Real(Std) :: Va
  Real(Std) :: Vb
  Real(Std) :: ErfVaM1       ! Erf(Va / Sqrt(2)) - 1.
  Real(Std) :: ErfVbM1       ! Erf(Vb / Sqrt(2)) - 1.
  Real(Std) :: ErfVaP1       ! Erf(Va / Sqrt(2)) + 1.
  Real(Std) :: ErfVbP1       ! Erf(Vb / Sqrt(2)) + 1.
  Real(Std) :: ErfVa         ! Erf(Va / Sqrt(2)).
  Real(Std) :: ErfVb         ! Erf(Vb / Sqrt(2)).
  Real(Std) :: PdfA
  Real(Std) :: PdfB
  Real(Std) :: TotalPdf
  Real(Std) :: dAlphadZ
  Real(Std) :: dBetadZ
  Real(Std) :: dGammadZ
  Real(Std) :: dSigWAdZ
  Real(Std) :: dSigWBdZ
  Real(Std) :: dWMeanAdZ
  Real(Std) :: dWMeanBdZ
  Real(Std) :: dAzSigWAdZ
  Real(Std) :: dAzdZ
  Real(Std) :: Phi
  Real(Std) :: Q
  Real(Std) :: mu            ! mean vertical velocity

  A(1) = - Dt * Extra%U(1) / Flow%TauUU(1)
  A(2) = - Dt * Extra%U(2) / Flow%TauUU(2)

  ! Gaussian turbulence formula. We use this unless skew turbulence requested, the skewness is non-zero and
  ! WVarianceDamping is non-zero (the skew formula will not work if skewness or WVarianceDamping is zero).
  If (.not.Extra%Skew .or. Flow%Sk == 0.0 .or. WVarianceDamping == 0.0) Then

!    A(3) = Dt * ( - Extra%U(3) / Flow%TauUU(3) + 0.5 * Flow%dSigUUdX(3))
!    A(3) = A(3) + 0.5 * (Flow%SigUU(3) - Extra%SigUU(3)) * Extra%U(3) / Flow%SigUU(3)
    A(3) = Extra%U(3) * sqrt(Flow%SigUU(3)) / sqrt(Extra%SigUU(3)) * (1.0 - Dt / Flow%TauUU(3))
    A(3) = A(3) + Dt * 0.5 * Flow%dSigUUdX(3) * WVarianceDamping * ZGradientDamping - Extra%U(3)

  ! Skew turbulence formula.
  Else

    SigWWDamped    = Flow%SigUU(3) * WVarianceDamping
    SkDamped       = Flow%Sk / WVarianceDamping**1.5
    dSigWWdZDamped = Flow%dSigUUdX(3) * WVarianceDamping
    dSkdZDamped    = Flow%dSkdZ / WVarianceDamping**1.5

    Alpha = SkDamped**(1.0/3.0)
    Beta  = SigWWDamped / (1.0 + Alpha * Alpha)
    Gamma = SkDamped * SigWWDamped**1.5 / (3.0 * Alpha + Alpha**3.0)

    SigWB = 0.5 * (Sqrt(Gamma * Gamma / (Beta * Beta) + 4.0 * Beta) - Gamma / Beta)
    SigWA = SigWB + Gamma / Beta

    WMeanA = Alpha * SigWA
    WMeanB = Alpha * SigWB

    Az = SigWB / (SigWA + SigWB)
    Bz = 1.0 - Az

    Va       = (Extra%U(3) - WMeanA)/SigWA
    Vb       = (Extra%U(3) + WMeanB)/SigWB

    dAlphadZ   = dSkdZDamped / (3.0 * SkDamped**(2.0/3.0))
    dBetadZ    = (dSigWWdZDamped - 2.0 * Alpha * SigWWDamped * dAlphadZ / (1.0 +  Alpha * Alpha)) / &
                 (1.0 +  Alpha * Alpha)
    dGammadZ   = (0.5 * dSigWWdZDamped + 2.0 * SigWWDamped * dAlphadZ / (Alpha * (3.0 + Alpha * Alpha))) * &
                 3.0 * Alpha * Alpha * Sqrt(SigWWDamped) / (3.0 + Alpha * Alpha)
    dSigWBdZ   = (0.5 * Gamma * dGammadZ / (Beta * Beta) + (1.0 - 0.5 * Gamma * Gamma / Beta**3.0) &
               * dBetadZ) / &
                 Sqrt(Gamma * Gamma / (Beta * Beta) + 4.0 * Beta) - 0.5 * dGammadZ / Beta + &
                 0.5 * Gamma * dBetadZ / (Beta * Beta)
    dSigWAdZ   = dSigWBdZ + dGammadZ / Beta - Gamma * dBetadZ / (Beta * Beta)
    dWMeanAdZ  = SigWA * dAlphadZ + Alpha * dSigWAdZ
    dWMeanBdZ  = SigWB * dAlphadZ + Alpha * dSigWBdZ
    dAzSigWAdZ = Az * Az * dSigWAdZ + Bz * Bz * dSigWBdZ
    dAzdZ      = (SigWA * dSigWBdZ - SigWB * dSigWAdZ) / (SigWA + SigWB)**2.0

    If (Particle%XOld(3) > Zs) Then
      mu = Umean - Particle%WSed
    Else
      mu = Umean - Particle%WSed * Particle%XOld(3) / Zs
    End If

    ! For extreme positive velocities we approximate A to avoid numerical inaccuracies and/or 0/0.
    If (Min(Va, Vb) > 10.0) Then

      ! See else section for full expression.
      ! Here we approximate 1 - Erf (V / Sqrt(2)) = exp(-V**2/2) [1 - 1 / V**2] / [Sqrt(Pi) V / Sqrt(2)].
      ! This gives
      ! ErfVaM1 = - Exp(- 0.5 * Va * Va) * (1.0 - 1.0 / Va**2) / (Sqrt(Pi/2.0) * Va)
      ! ErfVbM1 = - Exp(- 0.5 * Vb * Vb) * (1.0 - 1.0 / Vb**2) / (Sqrt(Pi/2.0) * Vb)
      ! We also use ErfVaM1, ErfVbM1, PdfA, PdfB, TotalPdf to indicate the actual values divided by PdfA or
      ! PdfB, depending on which of PdfA and PdfB is largest (as judged from the exponential term).
      If (Va <= Vb) Then ! Here we divide by PdfA.
        ErfVaM1  = - (1.0 - 1.0 / Va**2) * 2.0 * SigWA / Va
        ErfVbM1  = - (1.0 - 1.0 / Vb**2) * 2.0 * SigWA / Vb * Exp(- 0.5 * Vb * Vb + 0.5 * Va * Va)
        PdfA     = 1.0
        PdfB     = (SigWA / SigWB) * Exp(- 0.5 * Vb * Vb + 0.5 * Va * Va)
        TotalPdf = Az * PdfA + Bz * PdfB
      Else ! Here we divide by PdfB.
        ErfVaM1  = (1.0 - 1.0 / Va**2) * 2.0 * SigWB / Va * Exp(- 0.5 * Va * Va + 0.5 * Vb * Vb)
        ErfVbM1  = (1.0 - 1.0 / Vb**2) * 2.0 * SigWB / Vb
        PdfA     = (SigWB / SigWA) * Exp(- 0.5 * Va * Va + 0.5 * Vb * Vb)
        PdfB     = 1.0
        TotalPdf = Az * PdfA + Bz * PdfB
      End If

      Phi = -0.5 * ErfVaM1 * (Az * SigWA * dAlphadZ + Alpha * dAzSigWAdZ + mu * dAzdZ) + &
             0.5 * ErfVbM1 * (Bz * SigWB * dAlphadZ + Alpha * dAzSigWAdZ + mu * dAzdZ) + &
             (dAzSigWAdZ + (Alpha + mu / SigWA) * Az * dWMeanAdZ +                       &
             (Az * dWMeanAdZ + (Alpha + mu / SigWA) * Az * dSigWAdZ) * Va +              &
             Az * dSigWAdZ * Va * Va) * PdfA * SigWA +                                   &
             (dAzSigWAdZ + (Alpha - mu / SigWB) * Bz * dWMeanBdZ -                       &
             (Bz * dWMeanBdZ + (Alpha - mu / SigWB) * Bz * dSigWBdZ) * Vb +              &
             Bz * dSigWBdZ * Vb * Vb) * PdfB * SigWB
      Phi = Phi * ZGradientDamping

      Q = Az * PdfA * Va / SigWA + Bz * PdfB * Vb / SigWB

      A(3) = Dt * (- SigWWDamped * Q / Flow%TauUU(3) + Phi) / TotalPdf

    ! For extreme negative velocities we approximate A to avoid numerical inaccuracies and/or 0/0.
    Else If (Max(Va, Vb) < - 10.0) Then

      ! See else section for full expression.
      ! Here we approximate 1 + Erf (V / Sqrt(2)) = - exp(-V**2/2) [1 - 1 / V**2] / [Sqrt(Pi) V / Sqrt(2)].
      ! This gives
      ! ErfVaP1 = - Exp(- 0.5 * Va * Va) * (1.0 - 1.0 / Va**2) / (Sqrt(Pi/2.0) * Va)
      ! ErfVbP1 = - Exp(- 0.5 * Vb * Vb) * (1.0 - 1.0 / Vb**2) / (Sqrt(Pi/2.0) * Vb)
      ! We also use ErfVaM1, ErfVbM1, PdfA, PdfB, TotalPdf to indicate the actual values divided by PdfA or
      ! PdfB, depending on which of PdfA and PdfB is largest (as judged from the exponential term).
      If (Vb >= Va) Then ! Here we divide by PdfB.
        ErfVaP1  = - (1.0 - 1.0 / Va**2) * 2.0 * SigWB / Va * Exp(- 0.5 * Va * Va + 0.5 * Vb * Vb)
        ErfVbP1  = - (1.0 - 1.0 / Vb**2) * 2.0 * SigWB / Vb
        PdfA     = (SigWB / SigWA) * Exp(- 0.5 * Va * Va + 0.5 * Vb * Vb)
        PdfB     = 1.0
        TotalPdf = Az * PdfA + Bz * PdfB
      Else ! Here we divide by PdfA.
        ErfVaP1  = - (1.0 - 1.0 / Va**2) * 2.0 * SigWA / Va
        ErfVbP1  = - (1.0 - 1.0 / Vb**2) * 2.0 * SigWA / Vb * Exp(- 0.5 * Vb * Vb + 0.5 * Va * Va)
        PdfA     = 1.0
        PdfB     = (SigWA / SigWB) * Exp(- 0.5 * Vb * Vb + 0.5 * Va * Va)
        TotalPdf = Az * PdfA + Bz * PdfB
      End If

      Phi = -0.5 * ErfVaP1 * (Az * SigWA * dAlphadZ + Alpha * dAzSigWAdZ + mu * dAzdZ) + &
             0.5 * ErfVbP1 * (Bz * SigWB * dAlphadZ + Alpha * dAzSigWAdZ + mu * dAzdZ) + &
             (dAzSigWAdZ + (Alpha + mu / SigWA) * Az * dWMeanAdZ +                       &
             (Az * dWMeanAdZ + (Alpha + mu / SigWA) * Az * dSigWAdZ) * Va +              &
             Az * dSigWAdZ * Va * Va) * PdfA * SigWA +                                   &
             (dAzSigWAdZ + (Alpha - mu / SigWB) * Bz * dWMeanBdZ -                       &
             (Bz * dWMeanBdZ + (Alpha - mu / SigWB) * Bz * dSigWBdZ) * Vb +              &
             Bz * dSigWBdZ * Vb * Vb) * PdfB * SigWB
      Phi = Phi * ZGradientDamping

      Q = Az * PdfA * Va / SigWA + Bz * PdfB * Vb / SigWB

      A(3) = Dt * (- SigWWDamped * Q / Flow%TauUU(3) + Phi) / TotalPdf

    ! For non-extreme velocities we evaluate A exactly.
    Else

      ErfVa    = Erf(Va / Sqrt(2.0))
      ErfVb    = Erf(Vb / Sqrt(2.0))
      PdfA     = Exp(- 0.5 * Va * Va) / (Sqrt(2.0 * Pi) * SigWA)
      PdfB     = Exp(- 0.5 * Vb * Vb) / (Sqrt(2.0 * Pi) * SigWB)
      TotalPdf = Az * PdfA + Bz * PdfB

      Phi = -0.5 * ErfVa * (Az * SigWA * dAlphadZ + Alpha * dAzSigWAdZ + mu * dAzdZ) + &
             0.5 * ErfVb * (Bz * SigWB * dAlphadZ + Alpha * dAzSigWAdZ + mu * dAzdZ) + &
             (dAzSigWAdZ + (Alpha + mu / SigWA) * Az * dWMeanAdZ +                     &
             (Az * dWMeanAdZ + (Alpha + mu / SigWA) * Az * dSigWAdZ) * Va +            &
             Az * dSigWAdZ * Va * Va) * PdfA * SigWA +                                 &
             (dAzSigWAdZ + (Alpha - mu / SigWB) * Bz * dWMeanBdZ -                     &
             (Bz * dWMeanBdZ + (Alpha - mu / SigWB) * Bz * dSigWBdZ) * Vb +            &
             Bz * dSigWBdZ * Vb * Vb) * PdfB * SigWB
      Phi = Phi * ZGradientDamping

      Q = Az * PdfA * Va / SigWA + Bz * PdfB * Vb / SigWB

      A(3) = Dt * (- SigWWDamped * Q / Flow%TauUU(3) + Phi) / TotalPdf

    End If

  End If

  B(1) = Sqrt(2.0 * Dt * Flow%SigUU(1) / Flow%TauUU(1))
  B(2) = Sqrt(2.0 * Dt * Flow%SigUU(2) / Flow%TauUU(2))
  B(3) = Sqrt(2.0 * Dt * Flow%SigUU(3) * WVarianceDamping / Flow%TauUU(3))

End Subroutine XUSdeTerms

!-------------------------------------------------------------------------------------------------------------

Function InitW(Flow, WVarianceDamping, Skew, iP) Result(W)
! Generates a random vertical velocity.

  Implicit None
  ! Argument list:
  Type(Flow_), Intent(In) :: Flow             !
  Real(Std),   Intent(In) :: WVarianceDamping !
  Logical,     Intent(In) :: Skew
  Integer,     Intent(In) :: iP               ! Particle or puff index
  ! Function result:
  Real(Std) :: W
  ! Locals:
  Real(Std) :: SigWWDamped
  Real(Std) :: SkDamped
  Real(Std) :: Alpha
  Real(Std) :: Beta
  Real(Std) :: Gamma
  Real(Std) :: SigWA
  Real(Std) :: SigWB
  Real(Std) :: WMeanA
  Real(Std) :: WMeanB
  Real(Std) :: Az
  Real(Std) :: Bz

  If (.not.Skew .or. Flow%Sk == 0.0 .or. WVarianceDamping == 0.0) Then

    SigWWDamped = Flow%SigUU(3) * WVarianceDamping

    W = Sqrt(SigWWDamped) * Gauss(iP)

  Else

    SigWWDamped = Flow%SigUU(3) * WVarianceDamping
    SkDamped    = Flow%Sk / WVarianceDamping**1.5

    Alpha = SkDamped**(1.0/3.0)
    Beta  = SigWWDamped / (1.0 + Alpha * Alpha)
    Gamma = SkDamped * SigWWDamped**1.5 / (3.0 * Alpha + Alpha**3.0)

    SigWB = 0.5 * (Sqrt(Gamma * Gamma / (Beta * Beta) + 4.0 * Beta) - Gamma / Beta)
    SigWA = SigWB + Gamma / Beta

    WMeanA = Alpha * SigWA
    WMeanB = Alpha * SigWB

    Az = SigWB / (SigWA + SigWB)
    Bz = 1.0 - Az

!    Va       = (Extra%U(3) - WMeanA) / SigWA
!    Vb       = (Extra%U(3) + WMeanB) / SigWB
!    PdfA     = Exp(- 0.5 * Va * Va) / (Sqrt(2.0 * Pi) * SigWA)
!    PdfB     = Exp(- 0.5 * Vb * Vb) / (Sqrt(2.0 * Pi) * SigWB)
!    TotalPdf = Az * PdfA + Bz * PdfB

    Call GetRandomNumber(W, iP)
    If (W < Az) Then
      W = WMeanA + SigWA * Gauss(iP)
    Else
      W = - WMeanB + SigWB * Gauss(iP)
    End If

  End If

End Function InitW

!-------------------------------------------------------------------------------------------------------------

Function ReflectW(Flow, WVarianceDamping, Skew, WIn) Result(WOut)
! Generates a random vertical velocity.

  Implicit None
  ! Argument list:
  Type(Flow_), Intent(In) :: Flow             !
  Real(Std),   Intent(In) :: WVarianceDamping !
  Logical,     Intent(In) :: Skew
  Real(Std),   Intent(In) :: WIn              !
  ! Function result:
  Real(Std) :: WOut
  ! Locals:
  Real(Std) :: SigWWDamped
  Real(Std) :: SkDamped
  Real(Std) :: Alpha
  Real(Std) :: Beta
  Real(Std) :: Gamma
  Real(Std) :: SigWA
  Real(Std) :: SigWB
  Real(Std) :: WMeanA
  Real(Std) :: WMeanB
  Real(Std) :: Az
  Real(Std) :: Bz
  Real(Std) :: Va
  Real(Std) :: Vb
  Real(Std) :: Vr ! reference value
  Real(Std) :: WOutL
  Real(Std) :: WOutU
  Real(Std) :: IntegralWIn
  Real(Std) :: IntegralWOut
  Real(Std) :: IntegralWOutL
  Real(Std) :: IntegralWOutU
  Integer   :: i

  If (WIn > 0.0) Then

    WOut = WIn

  Else If (.not.Skew .or. Flow%Sk == 0.0) Then

    WOut = - WIn

  Else

    SigWWDamped = Flow%SigUU(3) * WVarianceDamping
    SkDamped    = Flow%Sk / WVarianceDamping**1.5

    Alpha = SkDamped**(1.0/3.0)
    Beta  = SigWWDamped / (1.0 + Alpha * Alpha)
    Gamma = SkDamped * SigWWDamped**1.5 / (3.0 * Alpha + Alpha**3.0)

    SigWB = 0.5 * (Sqrt(Gamma * Gamma / (Beta * Beta) + 4.0 * Beta) - Gamma / Beta)
    SigWA = SigWB + Gamma / Beta

    WMeanA =   Alpha * SigWA
    WMeanB = - Alpha * SigWB ! Note WMeanB < 0 - different convention to elsewhere.

    Az = SigWB / (SigWA + SigWB)
    Bz = 1.0 - Az

    Va = (WIn - WMeanA) / SigWA
    Vb = (WIn - WMeanB) / SigWB

    ! $$ extreme value (8.0) should be given some thought.

    ! For extreme positive velocities we approximate IntegralWIn to avoid numerical inaccuracies.
    If (Min(Va, Vb) > 8.0) Then

      ! See else section for full expression.
      ! Here we approximate 1 - Erf (V / Sqrt(2)) = exp(-V**2/2) [1 - 1 / V**2] / [Sqrt(Pi) V / Sqrt(2)].
      ! This gives
      ! Erf(Va / Sqrt(2)) - 1 = - Exp(- 0.5 * Va * Va) * (1 - 1 / Va**2) / (Sqrt(Pi/2) * Va)
      ! Erf(Vb / Sqrt(2)) - 1 = - Exp(- 0.5 * Vb * Vb) * (1 - 1 / Vb**2) / (Sqrt(Pi/2) * Vb)
      ! We divide IntegralWIn by Exp(- 0.5 * Vr * Vr) / Sqrt(2 * Pi).
      If (Va <= Vb) Then
        Vr = Va
      Else
        Vr = Vb
      End If
      IntegralWIn = + Az * SigWA                               * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                    + Az * WMeanA * ((1.0 - 1.0 / Va**2) / Va) * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                    + Bz * SigWB                               * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr) &
                    + Bz * WMeanB * ((1.0 - 1.0 / Vb**2) / Vb) * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr)

    ! For extreme negative velocities we approximate IntegralWIn to avoid numerical inaccuracies.
    Else If (Max(Va, Vb) < - 8.0) Then

      ! See else section for full expression.
      ! Here we approximate 1 + Erf (V / Sqrt(2)) = - exp(-V**2/2) [1 - 1 / V**2] / [Sqrt(Pi) V / Sqrt(2)].
      ! This gives
      ! Erf(Va / Sqrt(2)) + 1 = - Exp(- 0.5 * Va * Va) * (1 - 1 / Va**2) / (Sqrt(Pi/2) * Va)
      ! Erf(Va / Sqrt(2)) + 1 = - Exp(- 0.5 * Vb * Vb) * (1 - 1 / Vb**2) / (Sqrt(Pi/2) * Vb)
      ! We divide IntegralWIn by Exp(- 0.5 * Vr * Vr) / Sqrt(2 * Pi).
      If (Vb >= Va) Then
        Vr = Vb
      Else
        Vr = Va
      End If
      IntegralWIn = + Az * SigWA                               * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                    + Az * WMeanA * ((1.0 - 1.0 / Va**2) / Va) * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                    + Bz * SigWB                               * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr) &
                    + Bz * WMeanB * ((1.0 - 1.0 / Vb**2) / Vb) * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr)

    ! For non-extreme velocities we evaluate IntegralWIn exactly.
    ! Here we divide by Exp(- 0.5 * Vr * Vr) / Sqrt(2 * Pi).
    Else

      Vr = 0.0
      IntegralWIn = + Az * SigWA * Exp(- 0.5 * Va * Va)                  &
                    - Az * WMeanA * Sqrt(Pi / 2.0) * Erf(Va / Sqrt(2.0)) &
                    + Bz * SigWB * Exp(- 0.5 * Vb * Vb)                  &
                    - Bz * WMeanB * Sqrt(Pi / 2.0) * Erf(Vb / Sqrt(2.0))

    End If

    WOutU = 0.0
    Do i = 1, 20
      WOutU = WOutU - WIn
      Va   = (WOutU - WMeanA) / SigWA
      Vb   = (WOutU - WMeanB) / SigWB
      ! Evaluate IntegralWOut using the same approximations as for IntegralWIn above. Again we divide
      ! IntegralWOut by Exp(- 0.5 * Vr * Vr) / Sqrt(2 * Pi).
      If (Min(Va, Vb) > 8.0) Then
        IntegralWOutU = + Az * SigWA                               * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                        + Az * WMeanA * ((1.0 - 1.0 / Va**2) / Va) * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                        + Bz * SigWB                               * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr) &
                        + Bz * WMeanB * ((1.0 - 1.0 / Vb**2) / Vb) * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr)
      Else If (Max(Va, Vb) < - 8.0) Then
        IntegralWOutU = + Az * SigWA                               * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                        + Az * WMeanA * ((1.0 - 1.0 / Va**2) / Va) * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                        + Bz * SigWB                               * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr) &
                        + Bz * WMeanB * ((1.0 - 1.0 / Vb**2) / Vb) * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr)
      Else
        IntegralWOutU = + Az * SigWA * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr)                       &
                        - Az * WMeanA * Sqrt(Pi / 2.0) * Erf(Va / Sqrt(2.0)) * Exp(0.5 * Vr * Vr) &
                        + Bz * SigWB * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr)                       &
                        - Bz * WMeanB * Sqrt(Pi / 2.0) * Erf(Vb / Sqrt(2.0)) * Exp(0.5 * Vr * Vr)
      End If
      If (IntegralWOutU < IntegralWIn) Exit
      If (i == 20) Then
        write(6,*) 'loop 1 exit' ! $$ Convert to proper warning
        Exit
      End If
    End Do

    WOutL = WOutU
    Do i = 1, 20
      WOutL = WOutL / 2.0
      Va    = (WOutL - WMeanA) / SigWA
      Vb    = (WOutL - WMeanB) / SigWB
      ! Evaluate IntegralWOut using the same approximations as for IntegralWIn above. Again we divide
      ! IntegralWOut by Exp(- 0.5 * Vr * Vr) / Sqrt(2 * Pi).
      If (Min(Va, Vb) > 8.0) Then
        IntegralWOutL = + Az * SigWA                               * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                        + Az * WMeanA * ((1.0 - 1.0 / Va**2) / Va) * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                        + Bz * SigWB                               * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr) &
                        + Bz * WMeanB * ((1.0 - 1.0 / Vb**2) / Vb) * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr)
      Else If (Max(Va, Vb) < - 8.0) Then
        IntegralWOutL = + Az * SigWA                               * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                        + Az * WMeanA * ((1.0 - 1.0 / Va**2) / Va) * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                        + Bz * SigWB                               * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr) &
                        + Bz * WMeanB * ((1.0 - 1.0 / Vb**2) / Vb) * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr)
      Else
        IntegralWOutL = + Az * SigWA * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr)                       &
                        - Az * WMeanA * Sqrt(Pi / 2.0) * Erf(Va / Sqrt(2.0)) * Exp(0.5 * Vr * Vr) &
                        + Bz * SigWB * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr)                       &
                        - Bz * WMeanB * Sqrt(Pi / 2.0) * Erf(Vb / Sqrt(2.0)) * Exp(0.5 * Vr * Vr)
      End If
      If (IntegralWOutL > IntegralWIn) Exit
      If (i == 20) Then
        write(6,*) 'loop 2 exit' ! $$ Convert to proper warning
        Exit
      End If
    End Do

    Do i = 1, 20
      WOut = (WOutL + WOutU) / 2.0
      Va   = (WOut - WMeanA) / SigWA
      Vb   = (WOut - WMeanB) / SigWB
      ! Evaluate IntegralWOut using the same approximations as for IntegralWIn above. Again we divide
      ! IntegralWOut by Exp(- 0.5 * Vr * Vr) / Sqrt(2 * Pi).
      If (Min(Va, Vb) > 8.0) Then
        IntegralWOut = + Az * SigWA                               * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                       + Az * WMeanA * ((1.0 - 1.0 / Va**2) / Va) * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                       + Bz * SigWB                               * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr) &
                       + Bz * WMeanB * ((1.0 - 1.0 / Vb**2) / Vb) * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr)
      Else If (Max(Va, Vb) < - 8.0) Then
        IntegralWOut = + Az * SigWA                               * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                       + Az * WMeanA * ((1.0 - 1.0 / Va**2) / Va) * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr) &
                       + Bz * SigWB                               * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr) &
                       + Bz * WMeanB * ((1.0 - 1.0 / Vb**2) / Vb) * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr)
      Else
        IntegralWOut = + Az * SigWA * Exp(- 0.5 * Va * Va + 0.5 * Vr * Vr)                       &
                       - Az * WMeanA * Sqrt(Pi / 2.0) * Erf(Va / Sqrt(2.0)) * Exp(0.5 * Vr * Vr) &
                       + Bz * SigWB * Exp(- 0.5 * Vb * Vb + 0.5 * Vr * Vr)                       &
                       - Bz * WMeanB * Sqrt(Pi / 2.0) * Erf(Vb / Sqrt(2.0)) * Exp(0.5 * Vr * Vr)
      End If
      If (IntegralWOut < IntegralWIn) Then
        WOutU         = WOut
        IntegralWOutU = IntegralWOut
      Else
        WOutL         = WOut
        IntegralWOutL = IntegralWOut
      End If
      If (Abs(WOutU - WOutL) < 0.001 * Sqrt(Flow%SigUU(3))) Exit
      If (i == 20) Then
        write(6,*) 'loop 3 exit' ! $$ Convert to proper warning
        Exit
      End If
    End Do

  End If

End Function ReflectW

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateTime(SkewTime, VelMemTime, InhomogTime, MVelMemTime, Dt, Particle, Extra)
! .

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In)    :: SkewTime
  Type(ShortTime_), Intent(In)    :: VelMemTime
  Type(ShortTime_), Intent(In)    :: InhomogTime
  Type(ShortTime_), Intent(In)    :: MVelMemTime
  Type(ShortTime_), Intent(In)    :: Dt
  Type(Particle_),  Intent(InOut) :: Particle
  Type(Extra_),     Intent(InOut) :: Extra

  Particle%TOld = Particle%T
  Particle%T    = Particle%T + Dt

  If (Particle%T >= SkewTime   ) Extra%Skew    = .false.
  If (Particle%T >= VelMemTime ) Extra%VelMem  = .false.
  If (Particle%T >= InhomogTime) Extra%Inhomog = .false.
  If (Particle%T >= MVelMemTime) Extra%MVelMem = .false.

End Subroutine UpdateTime

!-------------------------------------------------------------------------------------------------------------

Subroutine DiffusionProcessReflect(Particle, Flow, RDt, Zs, WMeanAndTurb, DrydepFrac)

! Note this routine assumes Zs < Flow%ZInterface(1) (the boundary layer depth).

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(InOut) :: Particle
  Type(Flow_),     Intent(In)    :: Flow
  Real(Std),       Intent(In)    :: RDt
  Real(Std),       Intent(In)    :: Zs
  Real(Std),       Intent(InOut) :: WMeanAndTurb
  Real(Std),       Intent(Out)   :: DrydepFrac
  ! Locals:
  Logical   :: Crossed   ! True if particle has crossed an interface between turbulent layers
  Real(Std) :: W         ! Vertical velocity
  Real(Std) :: dTE       ! Time elapsed before crossing boundary
  Real(Std) :: Dt        ! Time remaining
  Real(Std) :: ZOld      !
  Real(Std) :: WUpper    ! Upper W in bisection method
  Real(Std) :: WLower    ! Lower W in bisection method
  Real(Std) :: WMid      ! midpoint W in bisection method
  Real(Std) :: WTemp1    ! Temporary value for interval bisection
  Real(Std) :: WTemp2    ! Temporary value for interval bisection
  Real(Std) :: WTempMid  ! Temporary value for interval bisection
  Real(Std) :: Tempwt2   ! Temporary value for transmitted velocity
  Real(Std) :: K3TurbLayers(MaxTurbLayers)  ! K values in each layer
  Integer   :: Layer     ! Turbulent layer
  Integer   :: BisLoop   ! Bisection Loop Index

  Dt     = RDt
  ZOld   = Particle%XOld(3)

  ! Determine which layer the particle is in at the start of the timestep
  Layer = 1
  Do
    If ( ZOld < Flow%ZInterface(Layer) ) Then
      Exit
    Else
      Layer = Layer + 1
      If ( Layer == Flow%nTurbLayers ) Exit
    End If
  End Do

  ! Number of Interfaces = Number of turbulent layers - 1
  ! Flow%nTurbLayers = number of turbulent layers
  ! Layer = 1 denotes boundary layer
  ! Flow%ZInterface(Layer) = interface at top of Layer
  ! Flow%ZInterface(1) = boundary layer height

  DrydepFrac = 0.0

  Do

    If (ZOld < Zs) Then
      W = WMeanAndTurb - Particle%WSed * ZOld / Zs
    Else
      W = WMeanAndTurb - Particle%WSed
    End If


    If ( Layer < Flow%nTurbLayers .and. Layer > 1 ) Then
      ! Test for particle crossing interface at top of layer
      If ( Particle%X(3) > Flow%ZInterface(Layer) ) Then
        Crossed = .true.
      ! Test for particle crossing interface at bottom of layer
      Else If ( Particle%X(3) < Flow%ZInterface(Layer - 1) ) Then
        Crossed = .true.
      Else
        Crossed = .false.
      End If
    ! For top Layer test only for particle crossing interface at bottom of layer
    Else If ( Layer == Flow%nTurbLayers) Then
      If ( Particle%X(3) < Flow%ZInterface(Layer - 1) ) Then
        Crossed = .true.
      Else
        Crossed = .false.
      End If
    ! For bottom layer test only for particle crossing interface at top of layer
    Else
      If ( Particle%X(3) > Flow%ZInterface(Layer) ) Then
        Crossed = .true.
      Else
        Crossed = .false.
      End If
    End If

    If (Crossed) Then

      ! For upward moving particles
      If ( W > 0.0 ) Then

        If (ZOld < Zs) Then
          If (Particle%WSed > 0.0) Then
            If (Particle%WSed * (ZOld - Zs) / (Particle%WSed * ZOld - WMeanAndTurb * Zs) < 0.01) Then
              dTE = Zs * (ZOld - Zs) / (Particle%WSed * ZOld - WMeanAndTurb * Zs)
            Else
              dTE = - Zs / Particle%WSed * Log(Zs * (Particle%WSed - WMeanAndTurb) /  &
                (Particle%WSed * ZOld - WMeanAndTurb * Zs))
            End If
          Else
            dTE = (Zs - ZOld) / WMeanAndTurb
          End If
          DrydepFrac = DrydepFrac + dTE / RDt
          dTE = dTE + (Flow%ZInterface(Layer) - Zs) / (WMeanAndTurb - Particle%WSed)
        Else
          dTE = (Flow%ZInterface(Layer) - ZOld) / (WMeanAndTurb - Particle%WSed)
        End If

        If (Dt - dTE < 0.0) Then
          dTE = Dt
          Dt = 0.0
        Else
          Dt  = Dt - dTE
        End If

        If (Particle%WSed > 0.0) Then
! Solve for transmitted velocity

! Interval bisection

          WLower = 0.0

          WUpper = 100.0

          Call TransmitInterfaceW(                                                            &
                 WMeanAndTurb, WLower, RDt, Flow, Layer, Particle%WSed, .true., WTemp1        &
               )

          Call TransmitInterfaceW(                                                            &
                 WMeanAndTurb, WUpper, RDt, Flow, Layer, Particle%WSed, .true., WTemp2        &
               )


! Particle transmitted
          If ((WTemp1 > 0.0 .and. WTemp2 < 0.0) .or. (WTemp1 < 0.0 .and. WTemp2 > 0.0)) Then

! Solve by iteration for transmitted velocity
            Do BisLoop=1,50
              WMid = (WLower + WUpper) / 2.0

              Call TransmitInterfaceW(                                                        &
                     WMeanAndTurb, WMid, RDt, Flow, Layer, Particle%WSed, .true., WTempMid    &
                   )

              If ((WTemp1 > 0.0 .and. WTempMid > 0.0) .or. (WTemp1 < 0.0 .and. WTempMid < 0.0)) Then
                WLower = WMid
              Else
                WUpper = WMid
              End If
            End Do

            WMeanAndTurb = (WLower + WUpper) / 2.0 + Particle%WSed
            Particle%X(3) = Flow%ZInterface(Layer) + (WMeanAndTurb - Particle%WSed) * Dt
            ZOld = Flow%ZInterface(Layer)
            Layer = Layer + 1

! else Particle reflected
          Else
            WUpper = -100.0

            Call ReflectInterfaceW(                                                           &
                   WMeanAndTurb, WLower, RDt, Flow, Layer, Particle%WSed, .true., WTemp1      &
                 )

            Call ReflectInterfaceW(                                                           &
                   WMeanAndTurb, WUpper, RDt, Flow, Layer, Particle%WSed, .true., WTemp2      &
                 )


            If ((WTemp1 > 0.0 .and. WTemp2 > 0.0) .or. (WTemp1 < 0.0 .and. WTemp2 < 0.0)) Then

! Precision problems possible when argument of exp() small.
! Set total reflected velocity equal to zero - should be as good as any approximation

              WMeanAndTurb = Particle%WSed
              Particle%X(3) = Flow%ZInterface(Layer)
              ZOld = Flow%ZInterface(Layer)
              Crossed = .false.

!              write (6,*) 'Error: bisection method failed'
!              write (6,*) 'reflected according to transmission'
!              write (6,*) 'but not reflected according to reflection'
!              stop
            Else
! Solve by iteration for refelected velocity
              Do BisLoop=1,50
                WMid = (WLower + WUpper) / 2.0

                Call ReflectInterfaceW(                                                       &
                       WMeanAndTurb, WMid, RDt, Flow, Layer, Particle%Wsed, .true., WTempMid  &
                     )

                If ((WTemp1 > 0.0 .and. WTempMid > 0.0) .or. (WTemp1 < 0.0 .and. WTempMid < 0.0)) Then
                  WLower = WMid
                Else
                  WUpper = WMid
                End If
              End Do

              WMeanAndTurb = (WLower + WUpper) / 2.0 + Particle%WSed
              Particle%X(3) = Flow%ZInterface(Layer) + (WMeanAndTurb - Particle%WSed) * Dt
              ZOld = Flow%ZInterface(Layer)
              If (Particle%X(3) < Zs) Then
                dTE = (Zs - Flow%ZInterface(Layer)) / (WMeanAndTurb - Particle%WSed)
                If (Particle%WSed / Zs * (Dt - dTE) > 0.01) Then
                  Particle%X(3) = WMeanAndTurb / Particle%WSed * Zs + (Zs - WMeanAndTurb * Zs / &
                    Particle%WSed) * Exp(-Particle%WSed / Zs * (Dt - dTE))
                End If
              End If

            End If

          End If

! no sedimetation - can solve explicitly
        Else

          K3TurbLayers(Layer)     = Flow%SigW2TurbLayers(Layer) * Flow%TauWTurbLayers(Layer)
          K3TurbLayers(Layer + 1) = Flow%SigW2TurbLayers(Layer + 1) * Flow%TauWTurbLayers(Layer + 1)

          Tempwt2 = 2.0 * K3TurbLayers(Layer + 1) / RDt * (WMeanAndTurb**2 * RDt / &
            (2.0 * K3TurbLayers(Layer)) + Log(K3TurbLayers(Layer + 1) / K3TurbLayers(Layer)))
          If (Tempwt2 > 0.0) Then
            WMeanAndTurb = Sqrt(Tempwt2)
            Particle%X(3) = Flow%ZInterface(Layer) + WMeanAndTurb * Dt
            ZOld = Flow%ZInterface(Layer)
            Layer = Layer + 1

          Else
            WMeanAndTurb = -WMeanAndTurb
            Particle%X(3) = Flow%ZInterface(Layer) + WMeanAndTurb * Dt
            ZOld = Flow%ZInterface(Layer)

          End If
        End If

! For particles moving downwards
      Else

        dTE = (Flow%ZInterface(Layer - 1) - ZOld) / (WMeanAndTurb - Particle%WSed)
        If (Dt - dTE < 0.0) Then
          dTE = Dt
          Dt = 0.0
        Else
          Dt  = Dt - dTE
        End If
        If (Particle%WSed > 0.0) Then
! Solve for transmitted velocity

! Interval bisection
! test transmitted velocity =0.0

          WLower = 0.0
          WUpper = -100.0

          Call TransmitInterfaceW(                                                             &
                 WMeanAndTurb, WLower, RDt, Flow, Layer, Particle%WSed, .false., WTemp1        &
               )

          Call TransmitInterfaceW(                                                             &
                 WMeanAndTurb, WUpper, RDt, Flow, Layer, Particle%WSed, .false., WTemp2        &
               )

! Particle transmitted
          If ((WTemp1 > 0.0 .and. WTemp2 < 0.0) .or. (WTemp1 < 0.0 .and. WTemp2 > 0.0)) Then

! Solve by iteration for transmitted velocity
            Do BisLoop=1,50
              WMid = (WLower + WUpper) / 2.0

              Call TransmitInterfaceW(                                                         &
                     WMeanAndTurb, WMid, RDt, Flow, Layer, Particle%WSed, .false., WTempMid    &
                   )

              If ((WTemp1 > 0.0 .and. WTempMid > 0.0) .or. (WTemp1 < 0.0 .and. WTempMid < 0.0)) Then
                WLower = WMid
              Else
                WUpper = WMid
              End If
            End Do

            WMeanAndTurb = (WLower + WUpper) / 2.0 + Particle%WSed
            Particle%X(3) = Flow%ZInterface(Layer - 1) + (WMeanAndTurb - Particle%WSed) * Dt
            ZOld = Flow%ZInterface(Layer - 1)
            Layer = Layer - 1
            If (Particle%X(3) < Zs) Then
              dTE = (Zs - Flow%ZInterface(Layer)) / (WMeanAndTurb - Particle%WSed)
              If (Particle%WSed / Zs * (Dt - dTE) > 0.01) Then
                Particle%X(3) = WMeanAndTurb * Zs / Particle%WSed + (Zs - WMeanAndTurb * Zs / &
                  Particle%WSed) * Exp(-Particle%WSed / Zs * (Dt - dTE))
              End If
            End If

          Else
            WUpper = 100.0
! Particle reflected

            Call ReflectInterfaceW(                                                            &
                   WMeanAndTurb, WLower, RDt, Flow, Layer, Particle%WSed, .false., WTemp1      &
                 )

            Call ReflectInterfaceW(                                                            &
                   WMeanAndTurb, WUpper, RDt, Flow, Layer, Particle%WSed, .false., WTemp2      &
                 )

            If ((WTemp1 > 0.0 .and. WTemp2 > 0.0) .or. (WTemp1 < 0.0 .and. WTemp2 < 0.0)) Then
! Precision problems possible if argument of exp() small.
! Set total reflected velocity equal to zero - should be as good as any approximation

              WMeanAndTurb = Particle%WSed
              Particle%X(3) = Flow%ZInterface(Layer - 1)
              ZOld = Flow%ZInterface(Layer - 1)
              Crossed = .false.

!              write (6,*) 'Error: bisection method failed'
!              write (6,*) 'particle reflected according to transmission formula'
!              write (6,*) 'but particle not reflected according to reflection formula'
!              stop
            Else
! Solve by iteration for refelected velocity
              Do BisLoop=1,50
                WMid = (WLower + WUpper) / 2.0

                Call ReflectInterfaceW(                                                        &
                       WMeanAndTurb, WMid, RDt, Flow, Layer, Particle%WSed, .false., WTempMid  &
                     )

                If ((WTemp1 > 0.0 .and. WTempMid > 0.0) .or. (WTemp1 < 0.0 .and. WTempMid < 0.0)) Then
                  WLower = WMid
                Else
                  WUpper = WMid
                End If
              End Do
              WMeanAndTurb = (WLower + WUpper) / 2.0 + Particle%WSed
              Particle%X(3) = Flow%ZInterface(Layer - 1) + (WMeanAndTurb - Particle%WSed) * Dt
              ZOld = Flow%ZInterface(Layer - 1)

            End If

          End If

! no sedimentation - solve explicitly
        Else

          K3TurbLayers(Layer - 1) = Flow%SigW2TurbLayers(Layer - 1) * Flow%TauWTurbLayers(Layer - 1)
          K3TurbLayers(Layer)     = Flow%SigW2TurbLayers(Layer) * Flow%TauWTurbLayers(Layer)

          Tempwt2 = 2.0 * K3TurbLayers(Layer - 1) / RDt * (WMeanAndTurb**2 * RDt / &
            (2.0 * K3TurbLayers(Layer)) + Log(K3TurbLayers(Layer - 1) / K3TurbLayers(Layer)))
          If (Tempwt2 > 0.0) Then
            WMeanAndTurb = -Sqrt(Tempwt2)
            Particle%X(3) = Flow%ZInterface(Layer - 1) + WMeanAndTurb * Dt
            ZOld = Flow%ZInterface(Layer - 1)
            Layer = Layer - 1
          Else
            WMeanAndTurb = - WMeanAndTurb
            Particle%X(3) = Flow%ZInterface(Layer - 1) + WMeanAndTurb * Dt
            ZOld = Flow%ZInterface(Layer - 1)
          End If
        End If

      End If

    End If

    ! $$ Is there a need to trap things going wrong due to precision - e.g. Dt ending up < 0?

    ! $$ Could simplify by having extra 0th layer below ground with K3 = 0.0? Efficiency implications?

    ! Mark advection as complete
    If ( Particle%X(3) >= 0.0 ) Then     ! $$ check - is ground level = 0.0?
      If ( .not.Crossed ) Then
        If (ZOld < Zs) Then
          If (Particle%X(3) > Zs) Then
            If (Particle%WSed > 0.0) Then
              If (Particle%WSed * (ZOld - Zs) / (Particle%WSed * ZOld - WMeanAndTurb * Zs) < 0.01) Then
                dTE = Zs * (ZOld - Zs) / (Particle%WSed * ZOld  - WMeanAndTurb * Zs)
              Else
                dTE = -Zs / Particle%WSed * Log(Zs * (Particle%WSed - WMeanAndTurb) /      &
                  (ZOld * Particle%WSed - WMeanAndTurb * Zs))
              End If
            Else
              dTE = (Zs - ZOld) / WMeanAndTurb
            End If
            DrydepFrac = DrydepFrac + dTE / RDt
          Else
            DrydepFrac = DrydepFrac + Dt / RDt
          End If
        Else
          If (Particle%X(3) < Zs) Then
            dTE = (Zs - ZOld) / (WMeanAndTurb - Particle%WSed)
            DrydepFrac = DrydepFrac + Max((Dt - dTE),0.0) / RDt
          End If
        End If
        Exit
      End If
    ! Reflect only if in BL
    Else If ( Layer == 1 ) Then
      If (ZOld > Zs) Then
        dTE = (Zs - ZOld) / (WMeanAndTurb - Particle%WSed)
        If (Dt - dTE < 0.0) Then
          dTE = Dt
          Dt = 0.0
        Else
          Dt = Dt - dTE
        End If
      End If
      If (Particle%WSed > 0.0) Then
        If (Particle%WSed * Min(Zs, ZOld) / (Particle%WSed * Min(Zs, ZOld) - WMeanAndTurb * Zs) < 0.01) Then
          dTE = Zs * Min(Zs, ZOld) / (Particle%WSed * Min(Zs, ZOld) - WMeanAndTurb * Zs)
        Else
          dTE = -Zs / Particle%WSed * Log(-WMeanAndTurb * Zs /                             &
            (Particle%WSed * Min(ZOld, Zs) - WMeanAndTurb * Zs))
        End If
      Else
        dTE = -Min(Zs, ZOld) / WMeanAndTurb
      End If                 ! $$ check ground level = 0.0
      DrydepFrac = DrydepFrac + Min(dTE , Dt) / RDt
      If (Dt - dTE < 0.0) Then
        dTE = Dt
        Dt = 0.0
      Else
        Dt = Dt - dTE
      End If
      WMeanAndTurb = -WMeanAndTurb
      If (Particle%WSed > 0.0) Then
        If (Particle%WSed * Dt / Zs < 0.01) Then
          Particle%X(3) = WMeanAndTurb * Dt
        Else
          Particle%X(3) = WMeanAndTurb * Zs / Particle%WSed - WMeanAndTurb * Zs /           &
            Particle%WSed * Exp(-Particle%WSed / Zs * Dt)
        End If
        If (Particle%X(3) > Zs) Then
          If (Particle%WSed / WMeanAndTurb < 0.01) Then
            dTE = Zs / WMeanAndTurb
          Else
            dTE = -Zs / Particle%WSed * Log((WMeanAndTurb - Particle%WSed) / WMeanAndTurb)
          End If
          If (Dt - dTE > 0.0) Then
            Particle%X(3) = Zs + (WMeanAndTurb - Particle%WSed) * (Dt - dTE)
          End If
        End If
      Else
        Particle%X(3) = WMeanAndTurb * Dt
      End If
      ZOld = 0.0             ! $$ check ground level = 0.0
    End If

  End Do

End Subroutine DiffusionProcessReflect

!-------------------------------------------------------------------------------------------------------------

Subroutine TransmitInterfaceW(WMeanAndTurb, WTest, RDt, Flow, Layer, WSed, Up, WReturn)
! Evaluates function f where f(x) = 0 is the equation to solve to get the transmitted velocity.

  Implicit None
  ! Argument list:
  Real(Std),       Intent(In)    :: WSed
  Type(Flow_),     Intent(In)    :: Flow
  Real(Std),       Intent(In)    :: RDt
  Real(Std),       Intent(In)    :: WMeanAndTurb
  Real(Std),       Intent(In)    :: WTest            ! transmitted W tested
  Real(Std),       Intent(Out)   :: WReturn          ! value of transmition equation for WTest
  Integer,         Intent(In)    :: Layer
  Logical,         Intent(In)    :: Up               ! flag = true if particle moving upwards

  ! Locals:
  Real(Std) :: Tempwi    ! Temporary value related to WMeanAndTurb
  Real(Std) :: Tempwrt   ! Temporary value related to WTest
  Real(Std) :: K3TurbLayers(MaxTurbLayers)  ! K values in each layer
  Integer   :: i                            ! Loop variable


! For upward moving particles

  If (Up) Then

! Calculate K values in each layer
    Do i = Layer, Layer + 1
      K3TurbLayers(i) = Flow%SigW2TurbLayers(i) * Flow%TauWTurbLayers(i)
    End Do

    Tempwi = WMeanAndTurb * Sqrt(RDt / K3TurbLayers(Layer)) / 2.0
    Tempwrt = (WTest + WSed) * Sqrt(RDt / K3TurbLayers(Layer + 1)) / 2.0

    If (Tempwi > 4.0) Then
      If (Tempwrt > 4.0) Then
        If ((Tempwrt**2 - Tempwi**2) < 0.0) Then
          WReturn = -Exp(Tempwrt**2 - Tempwi**2) * (WMeanAndTurb - WSed) / WMeanAndTurb *    &
                      Sqrt(2.0 * K3TurbLayers(Layer) / RDt) + WTest / (WTest + WSed) *       &
                      Sqrt(2.0 * K3TurbLayers(Layer + 1) / RDt)
        Else
          WReturn = -(WMeanAndTurb - WSed) * Sqrt(2.0 * K3TurbLayers(Layer) / RDt) /         &
                       WMeanAndTurb + Exp(Tempwi**2 - Tempwrt**2) * WTest / (WTest + WSed) * &
                       Sqrt(2.0 * K3TurbLayers(Layer + 1) / RDt)
        End If
      Else
        WReturn = -Exp(-Tempwi**2) * (WMeanAndTurb - WSed) / WMeanAndTurb *                  &
                    Sqrt(2.0 * K3TurbLayers(Layer) / RDt) + Exp(-Tempwrt**2) *               &
                    Sqrt(2.0 * K3TurbLayers(Layer + 1) / RDt) + Sqrt(Pi / 2.0) * WSed *      &
                    (Erf(Tempwrt) - 1.0)
      End If

    Else
      If (Tempwrt > 4.0) Then
        WReturn = -Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * Exp(-Tempwi**2) + WTest /         &
                    (WTest + WSed) * Sqrt(2.0 * K3TurbLayers(Layer + 1) / RDt) *             &
                    Exp(-Tempwrt**2) + Sqrt(Pi / 2.0) * WSed * (1.0 - Erf(Tempwi))
      Else
        WReturn = -Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * Exp(-Tempwi**2) +                 &
                    Sqrt(2.0 * K3TurbLayers(Layer + 1) / RDt) * Exp(-Tempwrt**2) -           &
                    Sqrt(Pi / 2.0) * WSed * (Erf(Tempwi) - Erf(Tempwrt))
      End If
    End If

 ! For particles moving downwards
  Else

! Calculate K values in each layer
    Do i = Layer - 1, Layer
      K3TurbLayers(i) = Flow%SigW2TurbLayers(i) * Flow%TauWTurbLayers(i)
    End Do

    Tempwi = WMeanAndTurb * Sqrt(RDt / K3TurbLayers(Layer)) / 2.0
    Tempwrt = (WTest + WSed) * Sqrt(RDt / K3TurbLayers(Layer - 1)) / 2.0

    If (Tempwi < -4.0) Then
      If (Tempwrt < -4.0) Then
        If ((Tempwrt**2-Tempwi**2) < 0.0) Then
          WReturn = -Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * (WMeanAndTurb - WSed) *             &
                      Exp(Tempwrt**2-Tempwi**2) / WMeanAndTurb + WTest *                         &
                      Sqrt(2.0 * K3TurbLayers(Layer - 1) / RDt) / (WTest + WSed)
        Else
          WReturn = -Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * (WMeanAndTurb - WSed) /             &
                      WMeanAndTurb + WTest * Exp(Tempwi**2-Tempwrt**2) *                         &
                      Sqrt(2.0 * K3TurbLayers(Layer - 1) / RDt) / (WTest + WSed)
        End If
      Else
        WReturn = -Exp(-Tempwi**2) * (WMeanAndTurb - WSed) / WMeanAndTurb *                   &
                    Sqrt(2.0 * K3TurbLayers(Layer) / RDt) + Exp(-Tempwrt**2) *                &
                    Sqrt(2.0 * K3TurbLayers(Layer - 1) / RDt) + Sqrt(Pi / 2.0) * WSed *       &
                    (Erf(Tempwrt) + 1.0)
      End If
    Else
      If (Tempwrt < -4.0) Then
        WReturn = -Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * Exp(-Tempwi**2) + WTest *             &
                    Sqrt(2.0 * K3TurbLayers(Layer - 1) / RDt) * Exp(-Tempwrt**2) /               &
                    (WTest + WSed) - Sqrt(Pi / 2.0) * WSed * (Erf(Tempwi) + 1.0)
      Else
        WReturn = -Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * Exp(-Tempwi**2) +                     &
                    Sqrt(2.0 * K3TurbLayers(Layer - 1) / RDt) * Exp(-Tempwrt**2) -               &
                    Sqrt(Pi / 2.0) * WSed * (Erf(Tempwi) - Erf(Tempwrt))
      End If
    End If
  End If



End Subroutine TransmitInterfaceW

!-------------------------------------------------------------------------------------------------------------

Subroutine ReflectInterfaceW(WMeanAndTurb, WTest, RDt, Flow, Layer, WSed, Up, WReturn)
! Evaluates function f where f(x) = 0 is the equation to solve to get the reflected velocity.

  Implicit None
  ! Argument list:
  Real(Std),       Intent(In)    :: WSed
  Type(Flow_),     Intent(In)    :: Flow
  Real(Std),       Intent(In)    :: RDt
  Real(Std),       Intent(In)    :: WMeanAndTurb
  Real(Std),       Intent(In)    :: WTest            ! reflected W tested
  Real(Std),       Intent(Out)   :: WReturn          ! value of reflection equation for WTest
  Integer,         Intent(In)    :: Layer
  Logical,         Intent(In)    :: Up               ! flag = true if particle moving upwards
  ! Locals:
  Real(Std) :: Tempwi    ! Temporary value related to WMeanAndTurb
  Real(Std) :: Tempwrt   ! Temporary value related to WTest
  Real(Std) :: K3TurbLayers(MaxTurbLayers)  ! K values in each layer


! Calculate K value in Layer
  K3TurbLayers(Layer) = Flow%SigW2TurbLayers(Layer) * Flow%TauWTurbLayers(Layer)


! For upward moving particles

  If (Up) Then
    Tempwi = WMeanAndTurb * Sqrt(RDt / K3TurbLayers(Layer)) / 2.0
    Tempwrt = (WTest + WSed) * Sqrt(RDt / K3TurbLayers(Layer)) / 2.0

    If (Tempwi > 4.0) Then
      If (Tempwrt < -4.0) Then
        WReturn = -(WMeanAndTurb - WSed) * Sqrt(2.0 * K3TurbLayers(Layer) / RDt) *            &
                    Exp(-Tempwi**2 ) / WMeanAndTurb + WTest * Exp (-Tempwrt**2) *             &
                    Sqrt(2.0 * K3TurbLayers(Layer) / RDt) / (WTest + WSed) -                  &
                    Sqrt(2.0 * Pi) * WSed
      Else
        WReturn = -(WMeanAndTurb - WSed) * Sqrt(2.0 * K3TurbLayers(Layer) / RDt) *            &
                    Exp(-Tempwi**2) / WMeanAndTurb + Sqrt(2.0 * K3TurbLayers(Layer) / RDt) *  &
                    Exp(-Tempwrt**2) + Sqrt(Pi / 2.0) * WSed * (Erf(Tempwrt) - 1.0)
      End If

    Else
      If (Tempwrt < -4.0) Then
        WReturn = -Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * Exp(-Tempwi**2) + WTest *          &
                    Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * Exp(-Tempwrt**2) /                &
                    (WTest + WSed) - Sqrt(Pi / 2.0) * WSed * (1.0 + Erf(Tempwi))
      Else
        WReturn = -Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * Exp(-Tempwi**2) +                  &
                    Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * Exp(-Tempwrt**2) -                &
                    Sqrt(Pi / 2.0) * WSed * (Erf(Tempwi) - Erf(Tempwrt))
      End If
    End If

 ! For particles moving downwards
  Else
    Tempwi = WMeanAndTurb * Sqrt(RDt / K3TurbLayers(Layer)) / 2.0
    Tempwrt = (WTest + WSed) * Sqrt(RDt / K3TurbLayers(Layer)) / 2.0

    If (Tempwi < -4.0) Then
      If (Tempwrt > 4.0) Then
        WReturn = -Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * (WMeanAndTurb - WSed) *            &
                    Exp(-Tempwi**2) / WMeanAndTurb + Sqrt(2.0 * K3TurbLayers(Layer) / RDt) *  &
                    WTest * Exp(-Tempwrt**2) / (WTest + WSed) + Sqrt(2.0 * Pi) * WSed
      Else
        WReturn = -Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * (WMeanAndTurb - WSed) *            &
                    Exp(-Tempwi**2) / WMeanAndTurb + Sqrt(2.0 * K3TurbLayers(Layer) / RDt) *  &
                    Exp(-Tempwrt**2) + Sqrt(Pi / 2.0) * WSed * (Erf(Tempwrt) + 1.0)
      End If
    Else
      If (Tempwrt > 4.0) Then
        WReturn = -Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * Exp(-Tempwi**2) + WTest *          &
                    Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * Exp(-Tempwrt**2) /                &
                    (WTest + WSed) + Sqrt(Pi / 2.0) * WSed * (1.0 - Erf(Tempwi))
      Else
        WReturn = -Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * Exp(-Tempwi**2) +                  &
                    Sqrt(2.0 * K3TurbLayers(Layer) / RDt) * Exp(-Tempwrt**2) -                &
                    Sqrt(Pi / 2.0) * WSed * (Erf(Tempwi) - Erf(Tempwrt))
      End If
    End If
  End If

End Subroutine ReflectInterfaceW

!-------------------------------------------------------------------------------------------------------------

Subroutine RadioactiveDecay(Mass, Specieses, RDt, Field)
! Applies radioactive decay to a particle.

! Code copied from decplume.f (NAME v6.9) and decpower.f (NAME v8.11)

  Implicit None
  ! Argument list:
  Real(Std),        Intent(InOut)           :: Mass(:)   ! Mass associated with particle.
  Type(Specieses_), Intent(In),    Target   :: Specieses ! A collection of species.
  Real(Std),        Intent(In)              :: RDt       ! Time step.
  Logical,          Intent(In)              :: Field     ! Flag to field or particle
  ! Locals:
  Real(Std) :: NewMass(Size(Mass))   ! Local variable of Mass to ensure the Mass is updated appropriately.
  Integer   :: nTracers              ! Number of tracers
  Integer   :: iTracer               ! Loop variable over species on particles or fields 
  Integer   :: iDaughter             ! Loop variable over all products in decay chain.
  Integer   :: iSize                 ! Field particle size index 
  Integer   :: j, k                  ! For Do Loops
  Integer   :: p, n                  ! For calculating factorials
  Integer   :: j1, j2, k1            ! Indices in Specieses of species in decay chain.
  Integer   :: i                     ! Index in Specieses of daughter product.
  Real(P64) :: L                     ! Variable used to give more precision to Log(2.0) and
                                     ! hence the calculation of Mass(:)
  Real(P64) :: D, E, A               ! For defining subsets of the radioactive decay calculation
                                     ! D = PRODUCT(decayconstant_k - decayconstant_j) from k = 1 to i
                                     ! but not for k = j.
                                     ! E = SUM((EXP^-decayconstant_j*RDt)/D) from j=1 to i.
                                     ! A = PRODUCT(decayconstant_j+1*BranchingRatio_j,j+1) from
                                     ! j = 1 to i - 1.
                                     ! Note that decayconstant = Log(2.0)*InvHalfLife
  Type(Species_), Pointer :: Species

  ! Test size of Mass? Other routines in Particle.F90 of Puff.F90 too (and Masses too)? $$

  NewMass(:) = 0.0
  If ( Field ) Then
    nTracers = Specieses%nFields
  Else
    nTracers = Specieses%nParticleSpecieses
  End If

  ! Loop over all species
  Do iTracer = 1, nTracers
   
    If ( Field ) Then
      Species => Specieses%Specieses( Specieses%iField2Species(iTracer) )
    Else
      Species => Specieses%Specieses( Specieses%iParticle2Species(iTracer) )
    End If

    If (Mass(iTracer) > 0.0) Then

      ! Calculate radioactive decay from parent and all subsequent daughter products.
      If (Species%InvHalfLife > 0.0) Then

        ! Loop over all products in decay chain A1 -> A2 -> A3 -> ... -> An
        Do iDaughter = 1, Species%DecayChainLength

          L = LOG(2.0)
          E = 0.0

          Do j = 1, iDaughter
            j1 = Species%DecayChain(j)
            D = 1.0

            Do k = 1, iDaughter
              If (k == j) CYCLE
              k1 = Species%DecayChain(k)
              D = D * (L * Specieses%Specieses(k1)%InvHalfLife &
                  - L * Specieses%Specieses(j1)%InvHalfLife)
            End Do

            E = E + Exp(-RDt * L * Specieses%Specieses(j1)%InvHalfLife) / D
          End Do

          A = 1.0

          Do j = 1, iDaughter - 1
            j1 = Species%DecayChain(j + 1)
            j2 = Species%DecayChain(j)
            A = A * (L * Specieses%Specieses(j1)%InvHalfLife) &
                * Specieses%Specieses(j2)%BranchingRatio
          End Do

          i = Species%DecayChain(iDaughter)

          ! Code to constrain the value of E which may diverge erroneously
          ! for decay chains including radionuclides with very long half lives.
          If (iDaughter > 1) Then
            If (E < 0.0) Then
              E = 0.0
            End If

            p = 1
            Do n = 1, iDaughter - 1
              p = p * n
            End Do

            If (E > RDt**(iDaughter - 1) / p) Then
              E = RDt**(iDaughter - 1) / p
            End If
          End If

          ! Update the 'new' mass.
          If ( Field ) Then
            iSize = Specieses%iField2Size(iTracer)
            NewMass(Specieses%iSpeciesAndSize2Field(i, iSize)) =                          &
              NewMass(Specieses%iSpeciesAndSize2Field(i, iSize)) + Mass(iTracer) * A * E
          Else
            NewMass(Specieses%iSpecies2Particle(i)) =                                     &
              NewMass(Specieses%iSpecies2Particle(i)) + Mass(iTracer) * A * E
          End If

        End Do

      End If

    End If
  End Do

  ! Updating mass, ensuring no modification of the mass for 'stable' radionuclides
  Do iTracer = 1, nTracers
    If ( Field ) Then
      Species => Specieses%Specieses( Specieses%iField2Species(iTracer) )
    Else 
      Species => Specieses%Specieses( Specieses%iParticle2Species(iTracer) )
    End If
    If (Species%InvHalfLife > 0.0) Then
      Mass(iTracer) = NewMass(iTracer)
    End If
  End Do

End Subroutine RadioactiveDecay

!-------------------------------------------------------------------------------------------------------------

Subroutine AgentDecay(Time, Mass, Specieses, RDt, Rain, RelHumidity, ZenAngle, Temperature, Field, TravelTime)
! Calculate the mass loss of a particle due to UV decay of biological and chemical
! agents, and deactivation of viable foot-and-mouth virus by low relative humidity and 
! midges by high rainfall.

! Code adapted from decayuv.f and  virus.f (NAME v8.11)
! $$ NAME code checks that UVLossRate lies between 0 and 100 in this routine
! $$ (i.e. at each time step and for each particle) - better to check input
! $$ values at the input stage.

  Implicit None
  ! Argument list:
  Type(Time_),      Intent(In)           :: Time
  Real(Std),        Intent(InOut)        :: Mass(:)
  Type(Specieses_), Intent(In),   Target :: Specieses
  Real(Std),        Intent(In)           :: RDt          ! Time step (seconds).
  Type(Rain_),      Intent(In)           :: Rain
  Real(Std),        Intent(In)           :: RelHumidity
  Real(Std),        Intent(In)           :: ZenAngle
  Real(Std),        Intent(In)           :: Temperature  ! Temperature at particle or field location (K)
  Logical,          Intent(In)           :: Field        ! Flag to indicate field or particle
  Type(ShortTime_), Intent(In), Optional :: TravelTime   !
  ! Locals:
  Type(Species_), Pointer   :: Species
  Integer                   :: nTracers             ! Number of tracers
  Integer                   :: iTracer              ! Species (or field) index 
  Real(Std),      Parameter :: RHCrit = 50.0        ! Threshold value for relative humidity.
  Real(Std),      Parameter :: RainCrit  = 1.0      ! Threshold value for rain deactivation of midges (in mm/hr).
  Real(Std)                 :: DecayRate            ! Actual UV decay rate per hour.
  Real(Std)                 :: k_H2                 ! Reaction Rate H2+OH
  Real(Std),      Parameter :: OHConcPerMonth(12) =                                                  &
       (/ 2.0e5, 3.0e5, 5.0e5, 1.1e6, 1.5e6, 1.8e6, 2.6e6, 2.3e6, 1.1e6, 3.0e5, 2.0e5, 1.0e5 /)
!      Estimated OH concentration (molecules cm-3) over NW Europe as provided by Dudley Shallcross
!      This could be improved by providing 3D field of OH via STOCHEM $$
  Real(Std) :: ParticleAge           ! Age of particle (in seconds) at start of timestep.
  Real(Std) :: EffectiveTime1        ! Effective time for decay at start.
  Real(Std) :: EffectiveTime2        ! Effective time for decay at end.

  If ( Field ) Then
    nTracers = Specieses%nFields
  Else
    nTracers = Specieses%nParticleSpecieses
  End If

  Do iTracer = 1, nTracers

    If ( Field ) Then
      Species => Specieses%Specieses( Specieses%iField2Species(iTracer) )
    Else
      Species => Specieses%Specieses( Specieses%iParticle2Species(iTracer) )
    End If

    ! Consider only those species with positive mass.
    If (Mass(iTracer) <= 0.0) Cycle

    ! Apply mass loss for each fmd species (assuming virus loss rate of 50% per hour).
    If ((Species%Name .CIEq. 'FOOT-AND-MOUTH') .or. (Species%Name .CIEq. 'TCID50')) Then
      If (RelHumidity < RHCrit) Then
        Mass(iTracer) = Mass(iTracer) * Exp(RDt * Log(0.5)/3600.0)
      End If
    End If

    ! Apply mass loss for midges at later time steps due to rain 
    If (Species%Name .CIEq. 'MIDGE') Then
      If (Rain%ConPpt + Rain%DynPpt > RainCrit) Then
        Mass(iTracer) = 0.0 
      End If
    End If

    ! Apply mass loss due to UV decay for each species.
    If (Species%UVLossRate > 0.0) Then
      DecayRate = 1.0 - Species%UVLossRate * Cos(ZenAngle)/100.0
      Mass(iTracer) = Mass(iTracer) * Exp(RDt * Log(DecayRate)/3600.0)
    End If

    ! Power law decay.
    If (.not. Field .and. Species%UsePowerLawDecay) Then
      If (.not.Present(TravelTime)) Then
        Call Message('UNEXPECTED FATAL ERROR in AgentDecay: travel time is not available', 4)
      End If
      ParticleAge    = ShortTime2RealTime(TravelTime)
      EffectiveTime1 = Max(Species%PowerLawDecayDelay, ParticleAge)
      EffectiveTime2 = Max(Species%PowerLawDecayDelay, ParticleAge + RDt)
      Mass(iTracer) = Mass(iTracer) * (EffectiveTime1 / EffectiveTime2) ** Species%PowerLawDecayExponent
    End If

    ! Removal of H2 by OH oxidation: for particles only 
    ! OH concentration varies by month (strictly should be a 3D field) $$
    ! k_H2 (rate constant for H2 OH reaction) has units cm3 molecule-1 s-1
    ! Scheme used in Price et al 2007 JGR Vol 112
    If ( Species%Name(1:8) .CIEq. 'HYDROGEN' ) Then
      k_H2 = 5.5e-12 * Exp(-2000.0/Temperature)
      Mass(iTracer) = Mass(iTracer) * Exp(- RDt * k_H2 * OHConcPerMonth(Time%Month))
    End If

  End Do

End Subroutine AgentDecay

!-------------------------------------------------------------------------------------------------------------

Subroutine DryDeposition(Specieses, Flow, Cloud, RDt, YLatLong, Surface,  &
                          Plant, DryDep, Particle, Mass, Zs, DrydepFrac)
!.

! Code copied from drydep.f (NAME v6.9)

  Implicit None
  ! Argument list:
  Type(Specieses_), Intent(In),   Target :: Specieses
  Type(Flow_),      Intent(In)           :: Flow
  Type(Cloud_),     Intent(In)           :: Cloud
  Real(Std),        Intent(In)           :: RDt
  Real(Std),        Intent(In)           :: Zs
  Real(Std),        Intent(In)           :: DrydepFrac
  Real(Std),        Intent(In)           :: YLatLong
  Type(Plant_),     Intent(In)           :: Plant
  Real(Std),        Intent(Out)          :: DryDep(MaxSpecieses)
  Type(Surface_),   Intent(In)           :: Surface
  Type(Particle_),  Intent(In)           :: Particle
  Real(Std),        Intent(InOut)        :: Mass(:)
  ! Locals:
  Type(Time_)          :: Time
  Real(Std)            :: Vd
  Integer              :: iParticleSpecies ! Loop variable over species on particles.
  Type(Species_), Pointer :: Species


! Calculate deposition for each species.

  If (DrydepFrac > 0.0) Then
    Do iParticleSpecies = 1, Specieses%nParticleSpecieses
     
      Species => Specieses%Specieses( Specieses%iParticle2Species(iParticleSpecies) )
      If (Mass(iParticleSpecies) > 0.0) Then

        Time = ShortTime2Time(ParticleTime(Particle))
        Vd = CalcVd(                                            &
               Species, Flow, Surface, Plant, Time,             &
               Cloud%Cloud3d, Particle%WSed, Particle%Diameter, &
               Particle%X(3), Zs, YLatLong                      &
             )

        If (Vd == 0.0) Then
          DryDep(iParticleSpecies) = 0.0
        Else
          DryDep(iParticleSpecies) = Mass(iParticleSpecies) * (1.0 - Exp( -Vd * DrydepFrac * RDt / Zs ))
        End If

        If (DryDep(iParticleSpecies) >= Mass(iParticleSpecies)) Then
          DryDep(iParticleSpecies) = 0.9 * Mass(iParticleSpecies)
        End If

        If (DryDep(iParticleSpecies) < 0.0) Then
          DryDep(iParticleSpecies) = 0.0
        End If

        Mass(iParticleSpecies) = Mass(iParticleSpecies) - DryDep(iParticleSpecies)

      Else

        DryDep(iParticleSpecies) = 0.0

      End If
    End Do

    DryDep(:) = DryDep(:)/RDt ! Dry dep rate

  Else

    DryDep(:) = 0.0

  End If

End Subroutine DryDeposition

!-------------------------------------------------------------------------------------------------------------

Subroutine WetDeposition(Specieses, Flow, Cloud, Rain, RDt, WetDep, Particle, Mass)
!.

! Code copied from wetdep.f (NAME v6.9)

  Implicit None
  ! Argument list:
  Type(Specieses_), Intent(In),   Target :: Specieses
  Type(Flow_),      Intent(In)           :: Flow
  Type(Cloud_),     Intent(In)           :: Cloud
  Type(Rain_),      Intent(In)           :: Rain
  Real(Std),        Intent(In)           :: RDt
  Real(Std),        Intent(Out)          :: WetDep(MaxSpecieses)
  Type(Particle_),  Intent(InOut)        :: Particle
  Real(Std),        Intent(InOut)        :: Mass(:)
  ! Locals:
  Real(Std) :: Lambda              ! Scavenging coefficient
  Integer   :: iParticleSpecies    ! Loop variable over species on particles.
  Type(Species_), Pointer :: Species



! Calculate deposition for each species.

  Do iParticleSpecies = 1, Specieses%nParticleSpecieses

    Species => Specieses%Specieses( Specieses%iParticle2Species(iParticleSpecies) )
    
    If (Mass(iParticleSpecies) > 0.0) Then
    
      ! Calculate wet scavenging coefficient
      Lambda = CalcWetScavCoeff(Species, Cloud, Rain, Particle%X(3), Flow%T)
                   
      If (Lambda > 0.0) Then
        WetDep(iParticleSpecies) = Mass(iParticleSpecies) * (1.0 - Exp(-Lambda * RDt))
        If (WetDep(iParticleSpecies) >= Mass(iParticleSpecies)) Then
          WetDep(iParticleSpecies) = Mass(iParticleSpecies)
        End If
        Mass(iParticleSpecies) = Mass(iParticleSpecies) - WetDep(iParticleSpecies)
      Else
        WetDep(iParticleSpecies) = 0.0
      End If

    Else
      WetDep(iParticleSpecies) = 0.0
    End If
  End Do

  WetDep(:) = WetDep(:)/RDt ! Wet dep rate

End Subroutine WetDeposition

!-------------------------------------------------------------------------------------------------------------

Subroutine DeepConvectionNew(Cloud, RDt, iP, Rain, Flow, XPa)
! Updated convection scheme based on mass fluxes.
! Note the vertical coordinate for the scheme is pressure (Pa).
! This subroutine calls the subroutines: Updraughts, Subsidence and Polint.

  Implicit None
  ! Argument list:
  Type(Cloud_),    Intent(In)    :: Cloud
  Real(Std),       Intent(In)    :: RDt         ! time step (Delta t)
  Integer,         Intent(In)    :: iP          ! Particle index
  Type(Rain_),     Intent(In)    :: Rain        
  Type(Flow_),     Intent(In)    :: Flow
  Real(Std),       Intent(InOut) :: XPa(3)      ! Particle position in pressure (Pa)
  ! Local variables:
  Integer              :: iLevel                    ! Counter.
  Integer              :: iIndex                    !
  Integer              :: MaxLevelCloud             ! max # of pressure levels.
  Integer, parameter   :: iLevelMax = 51            ! Max number of vertical pressure levels.
  Real(Std)            :: ConCloudBasePa            ! Local copy of convective cloud base in Pa.
  Real(Std)            :: ConCloudTopPa             ! Local copy of convective cloud top in Pa.
  Real(Std)            :: RNum                      ! Random number.
  Real(Std)            :: MMax, Beta1, Beta2        ! Variables for estimation of mass fluxes.
  Real(Std)            :: IntMassFlux               ! Integral of temporary mass fluxes over the cloud depth.
  Real(Std)            :: MassFluxCloudBase         ! Cloud base mass flux (Pa/s).
  Real(Std)            :: MFInterp                  !
  Real(Std)            :: LevelPa(iLevelMax)        ! Array of pressure corresponding to vertical levels.
  Real(Std)            :: MassFluxTemp(iLevelMax)   ! Array of temporary and non-dimensional mass fluxes.
  Real(Std)            :: MassFlux(iLevelMax)       ! Array of final mass fluxes (Pa/s).
  Real(Std)            :: DeltaPa                   ! Delta Pressure.
  Real(Std)            :: EntRate(iLevelMax-1)      ! Entrainment rate.
  Real(Std)            :: EntFlux(iLevelMax-1)      ! Entrainment flux.
  Real(Std)            :: DetFlux(iLevelMax-1)      ! Detrainment flux.
  Real(Std)            :: ProbEn(iLevelMax-1)       ! Probability of entrainment.
  Real(Std)            :: ProbUp(iLevelMax-1)       ! Probability of updraughts.
  Real(Std)            :: ProbDe(iLevelMax-1)       ! Probability of detrainment.
  Real(Std), Parameter :: PMin       = 10000.0      ! Minimum pressure (Pa).
  Real(Std), Parameter :: PRef       = 60000.0      ! Reference pressure (Pa).
  Real(Std), Parameter :: A1         = 0.5          ! Parameters for estimation of mass fluxe profiles.
  Real(Std), Parameter :: A2         = 0.2          !
  Real(Std), Parameter :: RegrFactor = 16556.057897 ! Regression factor for obtaining clous base mass flux
                                                    ! from total amount of convective precip. at the ground.
                                                    ! Includes unit convers. from mm/hr to mm/s=1kg/(m^2 s).
  Real(Std), Parameter :: EntFacDp   = 0.9          ! Coefficients for entrainment rate. See UM version.
  Real(Std), Parameter :: Ae2        = 1.5          !
  Real(Std), Parameter :: EntCoef    = 3.0          !
  Logical              :: PlumeUp                   ! Logical flag for updraughts.
  
  If (Flow%Backwards) Then
    Call Message("ERROR: New convection scheme not yet allowed to be used in backwards mode", 3)
  End If
  
  ! Counters initialisation to zero.
  IntMassFlux = 0
  
  ! Set local copies of convective cloud base and convective cloud top in Pa.
  ! This is necessary as they might be modified locally in this routine.
  ConCloudBasePa = Cloud%ConCloudBasePa
  ConCloudTopPa  = Cloud%ConCloudTopPa
  
  ! Set up of 'ad hoc' vertical pressure layers (min=3, max=50) based on the vertical extent of the cloud.
  ! Levels are inside the cloud + 1 below + 1 above.
  MaxLevelCloud = ((ConCloudBasePa - ConCloudTopPa)/2000)
  If (MaxLevelCloud > (iLevelMax-1) ) MaxLevelCloud = iLevelMax - 1
  If (MaxLevelCloud < 2) MaxLevelCloud = 2
  DeltaPa = min( ((ConCloudBasePa - ConCloudTopPa)/(MaxLevelCloud-1)), &
               ((Flow%PS - ConCloudTopPa)/(MaxLevelCloud-0.5)) )
  If ( ((ConCloudBasePa - ConCloudTopPa)/(MaxLevelCloud-1)) < &
     ((Flow%PS - ConCloudTopPa)/(MaxLevelCloud-0.5)) ) Then
     LevelPa(1) = ConCloudBasePa + DeltaPa*0.5
  Else
     LevelPa(1) = Flow%PS
  End If
  MaxLevelCloud = MaxLevelCloud+1
  Do iLevel=1, MaxLevelCloud
      LevelPa(iLevel) = LevelPa(1) - (iLevel-1) * DeltaPa
  Enddo
  
  ! Check whether cloud base is above the surface and adjust if not.
  If (ConCloudBasePa > Flow%PS) Then
      ConCloudBasePa = LevelPa(1) - DeltaPa*0.5
  End If
  
  ! Parametrization of mass fluxes using empirical formulas from CRM.
  ! Profiles depend on the height of the freezing level.
  ! Non-dimensional mass fluxes, defined on pressure levels, first and last are zero.
  Do iLevel=1, MaxLevelCloud
     If  ( (iLevel>1).And.(iLevel<MaxLevelCloud) ) Then
       If (                                                &
             ((Flow%FLPa - ConCloudTopPa) >= PMin)   .And. &
             ((ConCloudBasePa - Flow%FLPa) >= PMin)        &
          ) Then
               If (Flow%FLPa >= PRef) Then
                  MMax = 1.0 + A1 * (Flow%FLPa - ConCloudBasePa + PMin) / &
                                    (PRef - ConCloudBasePa + PMin)
               Else
                  MMax = A1*3
               End If
               Beta1 = Log(MMax)
               Beta2 = Log(MMax/A2)
               If (LevelPa(iLevel) >= Flow%FLPa) Then
                   MassFluxTemp(iLevel) = MMax * Exp(-Beta1 * ( (LevelPa(iLevel) - Flow%FLPa) / &
                                                              (ConCloudBasePa - Flow%FLPa) )**2)
               Else
                   MassFluxTemp(iLevel) = MMax * Exp(-Beta2 * ( (LevelPa(iLevel) - Flow%FLPa) / &
                                                              (ConCloudTopPa - Flow%FLPa) )**2)
               End If
          Else
            MMax = A1*2
            Beta2 = Log(MMax/A2)
            MassFluxTemp(iLevel) = MMax * Exp(-Beta2 * ( (LevelPa(iLevel) - ConCloudBasePa) / &
                                                         (ConCloudTopPa - ConCloudBasePa) )**2)
          End If
     Else
         MassFluxTemp(iLevel)=0.0
     End If
  End Do
  
  ! Calculation of integral of mass fluxes over the cloud depth.
  Do iLevel=1, MaxLevelCloud-1
      IntMassFlux = IntMassFlux + ( (MassFluxTemp(iLevel) + MassFluxTemp(iLevel+1)) / 2 ) * DeltaPa
  End Do

  ! The cloud base mass flux can be estimated knowing the total amount of precipitation at the ground.
  ! and the integral of the temporary mass fluxes
   MassFluxCloudBase = (RegrFactor * Rain%ConPpt) / IntMassFlux
  
  ! Calculation of final Mass Fluxes (Pa/s).
  Do iLevel=1, MaxLevelCloud
       MassFlux(iLevel) = MassFluxCloudBase * MassFluxTemp(iLevel)
  End Do

  ! Estimation of entrainment and detrainment fluxes and associated probabilities for particles.
  ! defined on number of pressure levels -1.
  Do iLevel=1, MaxLevelCloud-1
     EntRate(iLevel) = EntFacDp * EntCoef * Ae2 * (LevelPa(iLevel)/(Flow%PS**2))
     ! Fix what happens at the first level (just below the cloud bottom).
     If (iLevel == 1) Then
         EntFlux(iLevel) = MassFlux(iLevel+1)
         DetFlux(iLevel) = 0.0
         ProbUp(iLevel) = 0.0
         ProbDe(iLevel)= Abs( DetFlux(iLevel) ) / ( MassFlux(iLevel) + EntFlux(iLevel) )
         ProbEn(iLevel) = RDt * EntFlux(iLevel) / DeltaPa
     ! Intermediate levels.
     Else if ( (iLevel>1).And.(iLevel<(MaxLevelCloud-1)) ) Then
         EntFlux(iLevel) = EntRate(iLevel) * ( MassFlux(iLevel) * DeltaPa)
         DetFlux(iLevel) = MassFlux(iLevel) - MassFlux(iLevel+1) + EntFlux(iLevel)
         If (DetFlux(iLevel) < 0.0) Then
            EntFlux(iLevel) = EntFlux(iLevel) - DetFlux(iLevel)
            DetFlux(iLevel) = MassFlux(iLevel) - MassFlux(iLevel+1) + EntFlux(iLevel)
         EndIf
         ProbEn(iLevel) = RDt * EntFlux(iLevel) / DeltaPa
         ProbDe(iLevel) = Abs( DetFlux(iLevel) ) / ( MassFlux(iLevel) + EntFlux(iLevel) )
         ProbUp(iLevel) = MassFlux(iLevel+1) / ( MassFlux(iLevel) + EntFlux(iLevel) )
     ! Top Level.
     Else
         EntFlux = 0.0
         DetFlux = MassFlux(iLevel) - MassFlux(iLevel+1) + EntFlux(iLevel)
         ProbUp(iLevel) = 0.0
         ProbDe(iLevel) = 1.0
         ProbEn(iLevel) = 0.0
     End If
  End Do
  
! Vertical displacement.
! Find level at which the particle is.
  Do iIndex=1, MaxLevelCloud-1
     If (                                              &
         (XPa(3) > LevelPa(iIndex+1))            .And. &
         (XPa(3) <= LevelPa(iIndex))                   &
     ) Then
          iLevel = iIndex
          ! Check whether particle entrains.
          Call GetRandomNumber(Rnum, iP)
          ! Move parcels into the right plume column (for the time being only to the updraught one).
          If (Rnum <= ProbEn(iLevel)) Then
              PlumeUp = .true.
              Call Updraughts(iLevel, MaxLevelCloud, iP, XPa, PlumeUp, LevelPa, DeltaPa, ProbDe)
          End If
     End If
  End Do
  ! Apply subsidence to all the particles.
  Call Subsidence(RDt, iP, XPa, LevelPa, MassFlux, MaxLevelCloud, PlumeUp, Flow)

End Subroutine DeepConvectionNew

!--------------------------------------------------------------------------------------------------------------------

Subroutine Subsidence(RDt, iP, XPa, LevelPa, MassFlux, MaxLevelCloud, PlumeUp, Flow)
! This subroutine applies a subsidence flux equal to the updraught flux to all particles
! to compensate for upward motion, so that the resulting particle pressure is larger than the
! one in place before applying subsidence.
  
  Implicit None
  ! Argument list:
  Type(Flow_), Intent(In)           :: Flow
  Real(Std), Intent(In)             :: RDt
  Integer, Intent(In)               :: iP
  Integer, Intent(In)               :: MaxLevelCloud
  Real(Std), Intent(In)             :: LevelPa(MaxLevelCloud)
  Real(Std), Intent(In)             :: MassFlux(MaxLevelCloud)
  Real(Std), Intent(InOut)          :: XPa(3)
  Logical, Intent(InOut)            :: PlumeUp
  ! Locals
  Real(Std) :: MFInterp
  Real(Std) :: Pres(2), MassF(2)
  Integer   :: iLevel
  Real(Std) :: RNum

  !At the end of each timestep, parcels are brought back to the environment.
  Do iLevel=1, MaxLevelCloud-1
      If (                                                     &
          (XPa(3) > LevelPa(iLevel+1))                  .And.  &
          (XPa(3) <= LevelPa(iLevel))                          &
      ) Then
           Pres(1) = LevelPa(iLevel)
           MassF(1) = MassFlux(iLevel)
           Pres(2) = LevelPa(iLevel+1)
           MassF(2) = MassFlux(iLevel+1)
           Call Polint(Pres, MassF, XPa(3), MFInterp)
           XPa(3) = XPa(3) + MFInterp * RDt
           If (XPa(3) > LevelPa(1)) Then
              Call GetRandomNumber(RNum, iP)
              XPa(3) = LevelPa(1) - RNum * (LevelPa(1) - LevelPa(2))
           End If
      End If
  End Do
  PlumeUp = .false.

End Subroutine Subsidence

!--------------------------------------------------------------------------------------------------------------------------

Subroutine Updraughts(iLevel, MaxLevelCloud, iP, XPa, PlumeUp, LevelPa, DeltaPa, ProbDe)
! This subroutine determines whether, after entraining in the cloud, a particle detrains or is caught in an updraught
! and, if the second condition is verified, updates the particle pressure.

  Implicit None
  ! Argument list:
  Integer, Intent(In)            :: iLevel, iP
  Integer, Intent(In)            :: MaxLevelCloud
  Real(Std), Intent(In)          :: LevelPa(MaxLevelCloud), DeltaPa, ProbDe(MaxLevelCloud)
  Real(Std), Intent(InOut)       :: XPa(3)
  Logical, Intent(InOut)         :: PlumeUp
  ! Local
  Real(Std)                :: RNum
  Integer                  :: iIndex
  
  iIndex = iLevel
  Do While (iIndex <= MaxLevelCloud-1)
    If (PlumeUp) Then
      Call GetRandomNumber(Rnum, iP)
      If (RNum <= ProbDe(iIndex)) Then
          PlumeUp = .false.
      Else
          Call GetRandomNumber(Rnum, iP)
          XPa(3) = LevelPa(iIndex+1) - Rnum * DeltaPa
      End If
      If ((XPa(3) >= 110000.0) .Or. (XPa(3) <= 1000.0)) Then ! Check this arbitrary values?
          Call Message("ERROR: Updraughts have been moving particles below the surface or above the tropopause.",2)
      End If
    End If
   iIndex = iIndex + 1
  End Do

End Subroutine Updraughts

!-----------------------------------------------------------------------------------------------------------------------

Subroutine Polint(press,massf,given_p,desired_massflux)
! This subroutine is used inside subroutine Subsidence to calculate the mass flux pressure at the particle position
! given the mass fluxes at the pressure levels below and above the particle.
! Note: from Fortran recipes.

  Implicit None
  ! Argument list:
  Integer, Parameter  :: n=2, NMAX=10
  Real, Intent(In)    :: given_p, press(n),massf(n)
  Real, Intent(Out)   :: desired_massflux
  ! Locals
  Integer             :: i,m,ns
  Real                :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX),dy

  ns = 1
  dif = Abs(given_p-press(1))
  Do i = 1,n
     dift = Abs(given_p-press(i))
     If (dift < dif) Then
        ns = i
        dif = dift
     End If
     c(i) = massf(i)
     d(i) = massf(i)
  End Do
  desired_massflux = massf(ns)
  ns = ns-1
  Do m=1, n-1
     Do i=1, n-m
        ho = press(i)-given_p
        hp = press(i+m)-given_p
        w = c(i+1)-d(i)
        den = ho-hp
        If (den == 0) Call Message("UNEXPECTED FATAL ERROR: Failure in Subroutine Polint",4)
           den = w/den
           d(i) = hp*den
           c(i) = ho*den
     End Do
     If (2*ns < (n-m)) Then
        dy = c(ns+1)
     Else
        dy = d(ns)
        ns = ns-1
     End If
     desired_massflux = desired_massflux+dy
  End Do

End Subroutine Polint

!-------------------------------------------------------------------------------------------------------

Subroutine DeepConvectionOld(Cloud, RDt, Particle, iP)

  Implicit None
  ! Argument list:
  Type(Cloud_),    Intent(In)    :: Cloud
  Real(Std),       Intent(In)    :: RDt         ! time step
  Type(Particle_), Intent(InOut) :: Particle
  Integer,         Intent(In)    :: iP          ! Particle index
  ! Locals:
  real(STD)   :: temp                           ! temp random number
  real(STD)   :: CloudFracAdjusted
  real(STD)   :: CloudDepth
  ! note cloud bases/tops assumed in same units as particle%X(3) and
  ! increasing upwards $$
  ! $$ Really want uniform dist in p not z

! cloud depth and amount

  CloudDepth = Cloud%ConCloudTop - Cloud%ConCloudBase

  CloudFracAdjusted = (Cloud%ConCloud/25.0) * (RDT / 3600.0)


! apply convective mixing if particle lies between cloud bottom and top

  If (                                          &
    (Particle%x(3) <= Cloud%ConCloudTop)  .and. &
    (Particle%x(3) >= Cloud%ConCloudBase) .and. &
    (Cloud%ConCloud > 0.0)                .and. &
    (CloudDepth > 500.0)                        &
  ) then

    call GetRandomNumber(temp, iP)

    if (temp < CloudFracAdjusted) then

      call GetRandomNumber(temp, iP)
      Particle%x(3) = Cloud%ConCloudBase + (temp * CloudDepth)

! set flag to re-initialise particle velocities ($$ but this isn't currently done)

      Particle%NeedToSetVel = .true.

    endif
  endif

End Subroutine DeepConvectionOld

!-------------------------------------------------------------------------------------------------------------

Subroutine ParticleEulerianFields(                         &
             Particle, Particle2, Mass, H3, Coords, Grids, &
             iT, iHGrid, iZGrid, Concentration             &
           )
! Calculate contribution to Eulerian fields from a particle.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In)           :: Particle  ! Particle (X etc in wrong coord sys).
  Type(Particle_), Intent(In)           :: Particle2 ! Particle with only X, iHCoord,
                                                     ! iZCoord defined (in right coord sys).
  Real(Std),       Intent(In)           :: Mass(:)   ! mass or dry, wet or total dep rate.
  Real(Std),       Intent(In)           :: H3        ! 
  Type(Coords_),   Intent(In)           :: Coords    ! 
  Type(Grids_),    Intent(In),   Target :: Grids     ! 
  Integer,         Intent(In)           :: iT        ! 
  Integer,         Intent(In)           :: iHGrid    ! 
  Integer,         Intent(In)           :: iZGrid    ! 
  Real(Con),       Intent(InOut)        :: Concentration(-1:,-1:,:,:,:) !$$ Halo for EulerianField. Should be better way.
  ! Locals:
  Integer   :: iX       !
  Integer   :: iY       !
  Integer   :: iZ       !
  Integer   :: jX       !
  Integer   :: jY       !
  Real(Std) :: X(2)
  Real(Std) :: Area
  Logical   :: VariableZGrid
  Type(HGrid_), Pointer :: HGrid
  Type(ZGrid_), Pointer :: ZGrid

  HGrid => Grids%HGrids(iHGrid)
  ZGrid => Grids%ZGrids(iZGrid)

  X(:) = Particle2%X(1:2)
  If (HGrid%XCycle > 0) Then
    Do While (X(1) < HGrid%XCentre - 0.5 * HGrid%XCycle)
      X(1) = X(1) + HGrid%XCycle
    End Do
    Do While (X(1) > HGrid%XCentre + 0.5 * HGrid%XCycle)
      X(1) = X(1) - HGrid%XCycle
    End Do
  End If
  If (HGrid%YCycle > 0) Then
    Do While (X(2) < HGrid%YCentre - 0.5 * HGrid%YCycle)
      X(2) = X(2) + HGrid%YCycle
    End Do
    Do While (X(2) > HGrid%YCentre + 0.5 * HGrid%YCycle)
      X(2) = X(2) - HGrid%YCycle
    End Do
  End If

  iX = NInt((X(1) - HGrid%X0) / HGrid%dX + 1.0)
  iY = NInt((X(2) - HGrid%Y0) / HGrid%dY + 1.0)
  If (iX > HGrid%nX .or. iX < 1 .or. iY > HGrid%nY .or. iY < 1) Return

  Call CalcArea(Coords%HCoords(Particle2%iHCoord), HGrid, iX, iY, Area)

  If (ZGrid%Variable) Then
    Do iZ = 1, ZGrid%nZ
      If ((Particle2%X(3) >= ZGrid%Z(iZ) - 0.5 * ZGrid%AvZ(iZ)) .and. &
          (Particle2%X(3) <= ZGrid%Z(iZ) + 0.5 * ZGrid%AvZ(iZ))) Then
!! !$OMP ATOMIC
        Concentration(iX, iY, iZ, :, iT) = Concentration(iX, iY, iZ, :, iT) + &
          Mass(:) / Abs(Area * ZGrid%AvZ(iZ) * H3)
      End If
    End Do
  Else
    iZ = NInt((Particle2%X(3) - ZGrid%Z0) / ZGrid%dZ + 1.0)
    If (iZ > ZGrid%nZ .or. iZ < 1) Return
!! !$OMP ATOMIC
    Concentration(iX, iY, iZ, :, iT) = Concentration(iX, iY, iZ, :, iT) + &
      Mass(:) / Abs(Area * ZGrid%dZ * H3)
  End If

End Subroutine ParticleEulerianFields

!-------------------------------------------------------------------------------------------------------------

Subroutine ParticleConc(                    &
             Particle, Particle2, Mass, H3, &
             Coords, Grids, Reqs,           &
             iField, iT, iLast,             &
             Results                        &
           )
! Calculate contribution to concentration from particle.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In)           :: Particle  ! Particle (X etc in wrong coord sys).
  Type(Particle_), Intent(In)           :: Particle2 ! Particle with only X, iHCoord,
                                                     ! iZCoord defined (in right coord sys).
  Real(Std),       Intent(In)           :: Mass(:)   ! mass or
                                                        ! dry, wet or total dep rate
  Real(Std),       Intent(In)           :: H3
  Type(Coords_),   Intent(In)           :: Coords    !
  Type(Grids_),    Intent(In),   Target :: Grids     !
  Type(Reqs_),     Intent(In)           :: Reqs      !
  Integer,         Intent(In)           :: iField    !
  Integer,         Intent(In)           :: iT        ! $$ rename iTA
  Integer,         Intent(In)           :: iLast     !
  Type(Results_),  Intent(InOut)        :: Results   !
  ! Locals:
  Integer   :: iX               !
  Integer   :: iY               !
  Integer   :: iZ               !
  Integer   :: jX               !
  Integer   :: jY               !
  Integer   :: iParticleSpecies !
  Integer   :: iHGrid           !
  Integer   :: iZGrid           !
  Real(Std) :: Volume           !
  Real(Std) :: X(2)
  Real(Std) :: Area 
  Logical   :: Unstructured
  Logical   :: VariableZGrid
  Type(HGrid_), Pointer :: HGrid
  Type(ZGrid_), Pointer :: ZGrid

  iParticleSpecies = Reqs%FieldReqs(iField)%iParticleSpecies
  iHGrid           = Reqs%FieldReqs(iField)%iHGrid
  iZGrid           = Reqs%FieldReqs(iField)%iZGrid

  Volume = 1.0

  If (iHGrid /= 0) Then
    HGrid => Grids%HGrids(iHGrid)
    If (HGrid%Unstructured) Then
      Unstructured = .true.
      iY = 1
    Else

      X(:) = Particle2%X(1:2)
      If (HGrid%XCycle > 0) Then
        Do While (X(1) < HGrid%XCentre - 0.5 * HGrid%XCycle)
          X(1) = X(1) + HGrid%XCycle
        End Do
        Do While (X(1) > HGrid%XCentre + 0.5 * HGrid%XCycle)
          X(1) = X(1) - HGrid%XCycle
        End Do
      End If
      If (HGrid%YCycle > 0) Then
        Do While (X(2) < HGrid%YCentre - 0.5 * HGrid%YCycle)
          X(2) = X(2) + HGrid%YCycle
        End Do
        Do While (X(2) > HGrid%YCentre + 0.5 * HGrid%YCycle)
          X(2) = X(2) - HGrid%YCycle
        End Do
      End If

      Unstructured = .false.
      iX = NInt((X(1) - HGrid%X0) / HGrid%dX + 1.0)
      iY = NInt((X(2) - HGrid%Y0) / HGrid%dY + 1.0)
      If (iX > HGrid%nX .or. iX < 1 .or. iY > HGrid%nY .or. iY < 1) Return
    End If

    Call CalcArea(Coords%HCoords(Particle2%iHCoord), HGrid, iX, iY, Area)

    Volume = Volume * Area

  Else
    Unstructured = .false.
    iX = 1
    iY = 1
  End If

  VariableZGrid = .false.
  If (iZGrid /= 0) Then
    ZGrid => Grids%ZGrids(iZGrid)
    If (ZGrid%Variable) Then
      VariableZGrid = .true.
    Else
      iZ = NInt((Particle2%X(3) - ZGrid%Z0) / ZGrid%dZ + 1.0)
      If (iZ > ZGrid%nZ .or. iZ < 1) Return
      Volume = Volume * ZGrid%dZ * H3
    End If
  Else If (Reqs%FieldReqs(iField)%AvBL) Then
    If (Particle2%X(3) <= Particle%H) Then
      iZ = 1
      Volume = Volume * Particle%H
    Else
      Return
    End If
  Else
    iZ = 1
  End If

  If (Unstructured) Then
    Do iX = 1, HGrid%nX

      X(:) = Particle2%X(1:2)
      If (HGrid%XCycle > 0) Then
        Do While (X(1) < HGrid%X(iX) - 0.5 * HGrid%XCycle)
          X(1) = X(1) + HGrid%XCycle
        End Do
        Do While (X(1) > HGrid%X(iX) + 0.5 * HGrid%XCycle)
          X(1) = X(1) - HGrid%XCycle
        End Do
      End If
      If (HGrid%YCycle > 0) Then
        Do While (X(2) < HGrid%Y(iX) - 0.5 * HGrid%YCycle)
          X(2) = X(2) + HGrid%YCycle
        End Do
        Do While (X(2) > HGrid%Y(iX) + 0.5 * HGrid%YCycle)
          X(2) = X(2) - HGrid%YCycle
        End Do
      End If

      jX = NInt((X(1) - HGrid%X(iX)) / HGrid%dX + 1.0)
      jY = NInt((X(2) - HGrid%Y(iX)) / HGrid%dY + 1.0)
      If (jX == 1 .and. jY == 1) Then
        If (VariableZGrid) Then
          Do iZ = 1, ZGrid%nZ
            If ((Particle2%X(3) >= ZGrid%Z(iZ) - 0.5 * ZGrid%AvZ(iZ)) .and. &
                (Particle2%X(3) <= ZGrid%Z(iZ) + 0.5 * ZGrid%AvZ(iZ))) Then
!$OMP ATOMIC
              Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iT, iLast) =       &
                Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iT, iLast) +     &
                Mass(iParticleSpecies) / Abs(Volume * ZGrid%AvZ(iZ) * H3)
            End If
          End Do
        Else
!$OMP ATOMIC
          Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iT, iLast) =   &
            Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iT, iLast) + &
            Mass(iParticleSpecies) / Abs(Volume)
        End If
      End If
    End Do
  Else
    If (VariableZGrid) Then
      Do iZ = 1, ZGrid%nZ
        If ((Particle2%X(3) >= ZGrid%Z(iZ) - 0.5 * ZGrid%AvZ(iZ)) .and. &
            (Particle2%X(3) <= ZGrid%Z(iZ) + 0.5 * ZGrid%AvZ(iZ))) Then
!$OMP ATOMIC
          Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iT, iLast) =       &
            Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iT, iLast) +     &
            Mass(iParticleSpecies) / Abs(Volume * ZGrid%AvZ(iZ) * H3)
        End If
      End Do
    Else
!$OMP ATOMIC
      Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iT, iLast) =   &
        Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iT, iLast) + &
        Mass(iParticleSpecies) / Abs(Volume)
    End If
  End If

End Subroutine ParticleConc

!-------------------------------------------------------------------------------------------------------------

Subroutine ParticleCloudGamma(               &
             Particle, Particle2, Mass,      &
             Specieses,                      &
             CloudGammaParamses,             &
             Coords, Grids, Reqs,            &
             iField, iT,                     &
             Results                         &
           )
! Sums photon flux over photon energy, defines the grids upon which the photon fluxes are calculated and
! determines the particle and receptor positions for the calculation of the distance between the two.
! $$ Assumes simplifications concerning coord systems: ZCoord is m agl; HGrid and ZGrid must
!    both be specified; HGrid is not wrapped

  Implicit None
  ! Argument list:
  Type(Particle_),           Intent(In)           :: Particle     ! Particle (X etc in wrong coord sys).
  Type(Particle_),           Intent(In)           :: Particle2    ! Particle with only X, iHCoord,
                                                                  ! iZCoord defined (in right coord sys).
  Real(Std),                 Intent(In)           :: Mass(:)      ! Mass (activity in Bq).
  Type(Coords_),             Intent(In)           :: Coords       ! Collection of co-ordinates.
  Type(Grids_),              Intent(In),   Target :: Grids        ! Collection of grids.
  Type(Reqs_),               Intent(In)           :: Reqs         ! Collection of requirements.
  Integer,                   Intent(In)           :: iField       ! Index of Field.
  Integer,                   Intent(In)           :: iT           ! $$ rename iTA
  Type(Results_),            Intent(InOut)        :: Results      ! Collection of results.
  Type(Specieses_),          Intent(In),   Target :: Specieses          ! A collection of species.
  Type(CloudGammaParamses_), Intent(In),   Target :: CloudGammaParamses ! A collection of sets of cloud
                                                                        ! gamma parameters.

  ! Locals:
  Integer                          :: iX               ! Index of x co-ordinate.
  Integer                          :: iY               ! Index of y co-ordinate.
  Integer                          :: iZ               ! Index of z co-ordinate.
  Integer                          :: iSpecies         ! Index of species in set of all species.
  Integer                          :: iParticleSpecies ! Index of species in set of species on particles.
  Integer                          :: iHGrid           ! Index of horizontal grid.
  Integer                          :: iZGrid           ! Index of vertical grid.
  Integer                          :: iHCoord          ! Index of horizontal co-ordinate system.
  Integer                          :: iZCoord          ! Index of vertical co-ordinate system.
  Integer                          :: iE               ! Photon Energy index for looping over nEnergies
                                                       ! for a single species.
  Integer                          :: nEnergies        ! Number of photon energies associated with a
                                                       ! single species.
  Real(Std)                        :: X_A(3)           ! Position of particle.
  Real(Std)                        :: X_B(3)           ! Position of receptor.
  Real(Std)                        :: DistanceAB       ! Distance between particle and receptor.
  Real(Std)                        :: PF(MaxEnergies)  ! Photon Flux.
  Type(HCoord_)                    :: HCoord           !} Horizontal and vertical coord systems (in
  Type(ZCoord_)                    :: ZCoord           !} correct systems for distance calculation).
  Type(HGrid_),            Pointer :: HGrid            ! Horizontal grid for field.
  Type(ZGrid_),            Pointer :: ZGrid            ! Vertical grid for field.
  Type(Species_),          Pointer :: Species          ! A species and associated parameters.
  Type(CloudGammaParams_), Pointer :: CloudGammaParams ! A set of cloud gamma parameters.

  iSpecies         = Reqs%FieldReqs(iField)%iSpecies
  iParticleSpecies = Reqs%FieldReqs(iField)%iParticleSpecies
  iHGrid           = Reqs%FieldReqs(iField)%iHGrid
  iZGrid           = Reqs%FieldReqs(iField)%iZGrid

  Species => Specieses%Specieses(iSpecies)
  If (Species%iCloudGammaParams == 0) Return
  CloudGammaParams => CloudGammaParamses%CloudGammaParamses(Species%iCloudGammaParams)

  ! Determine particle position.
  X_A(:)  = Particle2%X(1:3)
  iHCoord = Particle2%iHCoord
  iZCoord = Particle2%iZCoord

  ! Set coord systems used by the two locations
  HCoord = Coords%HCoords(iHCoord)
  ZCoord = Coords%ZCoords(iZCoord)

  ! Check that vertical coords are 'm agl'
  If (ZCoord%Name /= 'm agl') Then
    Call Message('UNEXPECTED FATAL ERROR in ParticleCloudGamma: the vertical coord system' // &
                 ' is not "m agl"', 4)
  End If

  ! Check that horizontal coord is lat-long or a cartesian-based system
  If (.not.(                               &
    HCoord%CoordType == H_LatLong     .or. &
    HCoord%CoordType == H_PSCartesian .or. &
    HCoord%CoordType == H_TMCartesian      &
  )) Then
    Call Message('UNEXPECTED FATAL ERROR in ParticleCloudGamma: unsupported horizontal coord system', 4)
  End If

  ! Determine horizontal and vertical grid.
  If (iHGrid /= 0) Then
    HGrid => Grids%HGrids(iHGrid)
  Else
    Call Message('FATAL ERROR in ParticleCloudGamma: horizontal grid required for ' // &
                 'cloud gamma dose calculation', 3)
  End If

  If (iZGrid /= 0) Then
    ZGrid => Grids%ZGrids(iZGrid)
  Else
    Call Message('FATAL ERROR in ParticleCloudGamma: vertical grid required for ' // &
                 'cloud gamma dose calculation', 3)
  End If

  ! Loop over vertical grid levels
  Do iZ = 1, ZGrid%nZ

    If (HGrid%Unstructured) Then
    ! Handle unstructured grids
      iY = 1
      Do iX = 1, HGrid%nX
        X_B(:) = (/ HGrid%X(iX), HGrid%Y(iX), ZGrid%Z(iZ) /)
        
        ! Call subroutine calculating the distance between the particle and the receptor.
        Call CalcDistanceBetweenTwoPoints(  &
               HCoord,                      &
               X_A(:), X_B(:),              &
               DistanceAB                   &
             )

        ! Call subroutine calculating the photon flux.
        Call ParticleFluence(CloudGammaParams, Mass(iParticleSpecies), DistanceAB, PF)

        ! Sum the photon fluxes over photon energy
        Do iE = 1, CloudGammaParams%nEnergies
!$OMP ATOMIC
          Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iT, iE) =   &
            Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iT, iE) + PF(iE)
        End Do
      End Do
    Else

    ! Handle structured grids
      Do iY = 1, HGrid%nY
        Do iX = 1, HGrid%nX
          X_B(:) = (/ HGrid%X(iX), HGrid%Y(iY), ZGrid%Z(iZ) /)
          
          ! Call subroutine calculating the distance between the particle and the receptor.
          Call CalcDistanceBetweenTwoPoints(  &
                 HCoord,                      &
                 X_A(:), X_B(:),              &
                 DistanceAB                   &
               )

          ! Call subroutine calculating the photon flux.
          Call ParticleFluence(CloudGammaParams, Mass(iParticleSpecies), DistanceAB, PF)

          ! Sum the photon fluxes over photon energy
          Do iE = 1, CloudGammaParams%nEnergies
!$OMP ATOMIC
            Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iT, iE) =   &
              Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iT, iE) + PF(iE)
          End Do
        End Do
      End Do
    End If

  End Do

End Subroutine ParticleCloudGamma

!-------------------------------------------------------------------------------------------------------------

Subroutine ParticleFluence(                  &
             CloudGammaParams,               &
             Mass, r, PF                     &
           )

! Performs calculations of photon flux (also known as the fluence) for individual photon energies (and
! for all species/radionuclides), for all particle -> receptor distances; for use in the cloud
! gamma dose calculation.

  Implicit None

  ! Argument list:
  Type(CloudGammaParams_),   Intent(In)    :: CloudGammaParams     ! A collection of cloud gamma parameters.
  Real(Std),                 Intent(In)    :: Mass                 ! Mass (Activity in Bq).
  Real(Std),                 Intent(InOut) :: r                    ! Distance between particle & receptor
                                                                   ! (metres).
  Real(Std),                 Intent(InOut) :: PF(:)                ! Photon flux - number of photons
                                                                   ! incident on
                                                                   ! the cross sectional area of a
                                                                   ! sphere (m^-2 s^-1).
  ! Locals:
  Integer   :: j                    ! Photon energy index.
  Integer   :: nEnergies            ! Number of photon energies for a single species.
  Real(Std) :: lac                  ! Linear attenuation coefficient (unitless) -
                                    ! parameterises the amount of photon attenuation.
  Real(Std) :: a                    ! Berger Build-up factor coeffs } Defines the amount
  Real(Std) :: b                    ! Berger Build-up factor coeffs } of scattering back
  Real(Std) :: BU                   ! Berger Build-up factor } to the recpetor.
  Real(Std) :: f                    ! Photon Intensity - defines the frequency of photon
                                    ! emissions.
  Real(Std) :: rmax                 ! rmax is the distance between the receptor and the particle,
                                    ! for Exp^(-lac * r * (1-b)) = Exp^(-100). rmax is required
                                    ! to enable very small insignificant doses as a result of very large
                                    ! mean free paths to be ignored.

  PF(:) = 0.0

  nEnergies = CloudGammaParams%nEnergies

  ! Loop over all photon energies of the species of interest
  Do j = 1, nEnergies

    ! Defining the build-up factor coefficients
    a = 0.0
    a = CloudGammaParams%BBuildUpFactora(j)
    b = 0.0
    b = CloudGammaParams%BBuildUpFactorb(j)

    ! Determining the linear attenuation coefficient
    lac = 0.0
    lac = CloudGammaParams%LinearAttenCoeff(j)

    ! Calculating the Berger approx of the build-up factor
    !BU = 0.0
    !BU = 1.0 + (a * lac * r * EXP(b * lac * r))

    ! Determining the Photon Intensity
    f = 0.0
    f = CloudGammaParams%PhotonIntensity(j)

    ! Setting lower boundary condition on the distance between the particle and receptor, r
    r = MAX(r,1.0E-03)

    ! Setting upper boundary condition on the distance between the particle and receptor, r,
    ! accounting for the energy of the photon emission
    rmax = 50 / (lac * (1 - b))

    If (r > rmax) Then
      PF(j) = 0.0
    Else
      ! Calculating the photon flux as a function of photon energy (m^-2 s^-1)
      PF(j) = (f * (Exp(-(lac * r)) + (a * lac * r * EXP(lac * r * (b - 1)))) * Mass)/(4.0 * Pi * r**2.0)
    EndIf

  End Do

End Subroutine ParticleFluence

!-------------------------------------------------------------------------------------------------------------

Subroutine ParticleNumbers(        &
             Particle, Particle2,  &
             Coords, Grids, Reqs,  &
             iField, iT, iLast,    & ! rename iT iTA $$
             Results               &
           )
! Calculate contribution to particle numbers from an individual particle.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In)           :: Particle  ! Particle (in wrong coord sys).
  Type(Particle_), Intent(In)           :: Particle2 ! Particle with only X, iHCoord,
                                                     ! iZCoord defined (in right coord sys).
  Type(Coords_),   Intent(In)           :: Coords    ! Collection of coord systems.
  Type(Grids_),    Intent(In),   Target :: Grids     ! Collection of grids.
  Type(Reqs_),     Intent(In)           :: Reqs      ! Collection of requirements.
  Integer,         Intent(In)           :: iField    ! Index of field.
  Integer,         Intent(In)           :: iT        ! Index of output time.
  Integer,         Intent(In)           :: iLast     ! Single index for 7th dimension
  Type(Results_),  Intent(InOut)        :: Results   ! Collection of results.
  ! Locals:
  Integer :: iX            !} Indices for particle numbers array.
  Integer :: iY            !}
  Integer :: iZ            !}
  Integer :: jX            !] Variables used with unstructured grids.
  Integer :: jY            !]
  Integer :: iHGrid        ! Index of horizontal grid (in Grids)
  Integer :: iZGrid        ! Index of vertical grid (in Grids)
  Real(Std) :: X(2)
  Logical :: Unstructured  ! Indicates that horizontal grid is unstructured.
  Logical :: VariableZGrid ! Indicates that a variable vertical grid is used.
  Type(HGrid_), Pointer :: HGrid ! Horizontal grid of particle numbers field.
  Type(ZGrid_), Pointer :: ZGrid ! Vertical grid of particle numbers field.

  iHGrid   = Reqs%FieldReqs(iField)%iHGrid
  iZGrid   = Reqs%FieldReqs(iField)%iZGrid

  If (iHGrid /= 0) Then
    HGrid => Grids%HGrids(iHGrid)
    If (HGrid%Unstructured) Then
      Unstructured = .true.
      iY = 1
    Else

      X(:) = Particle2%X(1:2)
      If (HGrid%XCycle > 0) Then
        Do While (X(1) < HGrid%XCentre - 0.5 * HGrid%XCycle)
          X(1) = X(1) + HGrid%XCycle
        End Do
        Do While (X(1) > HGrid%XCentre + 0.5 * HGrid%XCycle)
          X(1) = X(1) - HGrid%XCycle
        End Do
      End If
      If (HGrid%YCycle > 0) Then
        Do While (X(2) < HGrid%YCentre - 0.5 * HGrid%YCycle)
          X(2) = X(2) + HGrid%YCycle
        End Do
        Do While (X(2) > HGrid%YCentre + 0.5 * HGrid%YCycle)
          X(2) = X(2) - HGrid%YCycle
        End Do
      End If

      Unstructured = .false.
      iX = NInt((X(1) - HGrid%X0) / HGrid%dX + 1.0)
      iY = NInt((X(2) - HGrid%Y0) / HGrid%dY + 1.0)
      If (iX > HGrid%nX .or. iX < 1 .or. iY > HGrid%nY .or. iY < 1) Return
    End If
  Else
    Unstructured = .false.
    iX = 1
    iY = 1
  End If

  VariableZGrid = .false.
  If (iZGrid /= 0) Then
    ZGrid => Grids%ZGrids(iZGrid)
    If (ZGrid%Variable) Then
      VariableZGrid = .true.
    Else
      iZ = NInt((Particle2%X(3) - ZGrid%Z0) / ZGrid%dZ + 1.0)
      If (iZ > ZGrid%nZ .or. iZ < 1) Return
    End If
  Else If (Reqs%FieldReqs(iField)%AvBL) Then
    If (Particle2%X(3) <= Particle%H) Then
      iZ = 1
    Else
      Return
    End If
  Else
    iZ = 1
  End If

  If (Unstructured) Then
    Do iX = 1, HGrid%nX

      X(:) = Particle2%X(1:2)
      If (HGrid%XCycle > 0) Then
        Do While (X(1) < HGrid%X(iX) - 0.5 * HGrid%XCycle)
          X(1) = X(1) + HGrid%XCycle
        End Do
        Do While (X(1) > HGrid%X(iX) + 0.5 * HGrid%XCycle)
          X(1) = X(1) - HGrid%XCycle
        End Do
      End If
      If (HGrid%YCycle > 0) Then
        Do While (X(2) < HGrid%Y(iX) - 0.5 * HGrid%YCycle)
          X(2) = X(2) + HGrid%YCycle
        End Do
        Do While (X(2) > HGrid%Y(iX) + 0.5 * HGrid%YCycle)
          X(2) = X(2) - HGrid%YCycle
        End Do
      End If

      jX = NInt((X(1) - HGrid%X(iX)) / HGrid%dX + 1.0)
      jY = NInt((X(2) - HGrid%Y(iX)) / HGrid%dY + 1.0)
      If (jX == 1 .and. jY == 1) Then
        If (VariableZGrid) Then
          Do iZ = 1, ZGrid%nZ
            If ((Particle2%X(3) >= ZGrid%Z(iZ) - 0.5 * ZGrid%AvZ(iZ)) .and. &
                (Particle2%X(3) <= ZGrid%Z(iZ) + 0.5 * ZGrid%AvZ(iZ))) Then
!$OMP ATOMIC
              Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) =       &
                Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) + 1.0
            End If
          End Do
        Else
!$OMP ATOMIC
          Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) =       &
            Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) + 1.0
        End If
      End If
    End Do
  Else
    If (VariableZGrid) Then
      Do iZ = 1, ZGrid%nZ
        If ((Particle2%X(3) >= ZGrid%Z(iZ) - 0.5 * ZGrid%AvZ(iZ)) .and. &
            (Particle2%X(3) <= ZGrid%Z(iZ) + 0.5 * ZGrid%AvZ(iZ))) Then
!$OMP ATOMIC
          Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) =       &
            Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) + 1.0
        End If
      End Do
    Else
!$OMP ATOMIC
      Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) =       &
        Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) + 1.0
    End If
  End If

End Subroutine ParticleNumbers

!-------------------------------------------------------------------------------------------------------------

Subroutine ParticleNumbersBySpecies(    &
             Particle, Particle2, Mass, &
             Coords, Grids, Reqs,       &
             iField, iT, iLast,         & ! rename iT iTA $$
             Results                    &
           )
! Calculate contribution to particle numbers from an individual particle (by species).

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In)           :: Particle  ! Particle (in wrong coord sys).
  Type(Particle_), Intent(In)           :: Particle2 ! Particle with only X, iHCoord,
                                                     ! iZCoord defined (in right coord sys).
  Real(Std),       Intent(In)           :: Mass(:)
  Type(Coords_),   Intent(In)           :: Coords    ! Collection of coord systems.
  Type(Grids_),    Intent(In),   Target :: Grids     ! Collection of grids.
  Type(Reqs_),     Intent(In)           :: Reqs      ! Collection of requirements.
  Integer,         Intent(In)           :: iField    ! Index of field.
  Integer,         Intent(In)           :: iT        ! Index of output time.
  Integer,         Intent(In)           :: iLast     ! Single index for 7th dimension
  Type(Results_),  Intent(InOut)        :: Results   ! Collection of results.
  ! Locals:
  Integer :: iX               !} Indices for particle numbers array.
  Integer :: iY               !}
  Integer :: iZ               !}
  Integer :: jX               !] Variables used with unstructured grids.
  Integer :: jY               !]
  Integer :: iParticleSpecies ! Index of species for particle numbers field.
  Integer :: iHGrid           ! Index of horizontal grid (in Grids)
  Integer :: iZGrid           ! Index of vertical grid (in Grids)
  Real(Std) :: X(2)
  Logical :: Unstructured  ! Indicates that horizontal grid is unstructured.
  Logical :: VariableZGrid ! Indicates that a variable vertical grid is used.
  Type(HGrid_), Pointer :: HGrid ! Horizontal grid of particle numbers field.
  Type(ZGrid_), Pointer :: ZGrid ! Vertical grid of particle numbers field.

  iParticleSpecies = Reqs%FieldReqs(iField)%iParticleSpecies
  iHGrid           = Reqs%FieldReqs(iField)%iHGrid
  iZGrid           = Reqs%FieldReqs(iField)%iZGrid

  ! No contribution if particle does not carry any mass for this species.
  If (Mass(iParticleSpecies) <= 0.0) Return

  If (iHGrid /= 0) Then
    HGrid => Grids%HGrids(iHGrid)
    If (HGrid%Unstructured) Then
      Unstructured = .true.
      iY = 1
    Else

      X(:) = Particle2%X(1:2)
      If (HGrid%XCycle > 0) Then
        Do While (X(1) < HGrid%XCentre - 0.5 * HGrid%XCycle)
          X(1) = X(1) + HGrid%XCycle
        End Do
        Do While (X(1) > HGrid%XCentre + 0.5 * HGrid%XCycle)
          X(1) = X(1) - HGrid%XCycle
        End Do
      End If
      If (HGrid%YCycle > 0) Then
        Do While (X(2) < HGrid%YCentre - 0.5 * HGrid%YCycle)
          X(2) = X(2) + HGrid%YCycle
        End Do
        Do While (X(2) > HGrid%YCentre + 0.5 * HGrid%YCycle)
          X(2) = X(2) - HGrid%YCycle
        End Do
      End If

      Unstructured = .false.
      iX = NInt((X(1) - HGrid%X0) / HGrid%dX + 1.0)
      iY = NInt((X(2) - HGrid%Y0) / HGrid%dY + 1.0)
      If (iX > HGrid%nX .or. iX < 1 .or. iY > HGrid%nY .or. iY < 1) Return
    End If
  Else
    Unstructured = .false.
    iX = 1
    iY = 1
  End If

  VariableZGrid = .false.
  If (iZGrid /= 0) Then
    ZGrid => Grids%ZGrids(iZGrid)
    If (ZGrid%Variable) Then
      VariableZGrid = .true.
    Else
      iZ = NInt((Particle2%X(3) - ZGrid%Z0) / ZGrid%dZ + 1.0)
      If (iZ > ZGrid%nZ .or. iZ < 1) Return
    End If
  Else If (Reqs%FieldReqs(iField)%AvBL) Then
    If (Particle2%X(3) <= Particle%H) Then
      iZ = 1
    Else
      Return
    End If
  Else
    iZ = 1
  End If

  If (Unstructured) Then
    Do iX = 1, HGrid%nX

      X(:) = Particle2%X(1:2)
      If (HGrid%XCycle > 0) Then
        Do While (X(1) < HGrid%X(iX) - 0.5 * HGrid%XCycle)
          X(1) = X(1) + HGrid%XCycle
        End Do
        Do While (X(1) > HGrid%X(iX) + 0.5 * HGrid%XCycle)
          X(1) = X(1) - HGrid%XCycle
        End Do
      End If
      If (HGrid%YCycle > 0) Then
        Do While (X(2) < HGrid%Y(iX) - 0.5 * HGrid%YCycle)
          X(2) = X(2) + HGrid%YCycle
        End Do
        Do While (X(2) > HGrid%Y(iX) + 0.5 * HGrid%YCycle)
          X(2) = X(2) - HGrid%YCycle
        End Do
      End If

      jX = NInt((X(1) - HGrid%X(iX)) / HGrid%dX + 1.0)
      jY = NInt((X(2) - HGrid%Y(iX)) / HGrid%dY + 1.0)
      If (jX == 1 .and. jY == 1) Then
        If (VariableZGrid) Then
          Do iZ = 1, ZGrid%nZ
            If ((Particle2%X(3) >= ZGrid%Z(iZ) - 0.5 * ZGrid%AvZ(iZ)) .and. &
                (Particle2%X(3) <= ZGrid%Z(iZ) + 0.5 * ZGrid%AvZ(iZ))) Then
!$OMP ATOMIC
              Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) =       &
                Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) + 1.0
            End If
          End Do
        Else
!$OMP ATOMIC
          Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) =       &
            Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) + 1.0
        End If
      End If
    End Do
  Else
    If (VariableZGrid) Then
      Do iZ = 1, ZGrid%nZ
        If ((Particle2%X(3) >= ZGrid%Z(iZ) - 0.5 * ZGrid%AvZ(iZ)) .and. &
            (Particle2%X(3) <= ZGrid%Z(iZ) + 0.5 * ZGrid%AvZ(iZ))) Then
!$OMP ATOMIC
          Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) =       &
            Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) + 1.0
        End If
      End Do
    Else
!$OMP ATOMIC
      Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) =       &
        Results%Fields(iField)%P64(1, iX, iY, iZ, 1, iT, iLast) + 1.0
    End If
  End If

End Subroutine ParticleNumbersBySpecies

!-------------------------------------------------------------------------------------------------------------

!## Updated to include the travel time
!## discrepancy between the output grid time and the actual particle time.

Subroutine ParticleMeanTravelTime(Particle, Particle2, Mass, Coords, Grids, Reqs, &
                                  iField, iT, iTA, iLast, Results)
! Calculate contribution to (mass-weighted) mean travel-time field from particle.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In)           :: Particle  !
  Type(Particle_), Intent(In)           :: Particle2 !
  Real(Std),       Intent(In)           :: Mass(:)
  Type(Coords_),   Intent(In)           :: Coords    !
  Type(Grids_),    Intent(In),   Target :: Grids     !
  Type(Reqs_),     Intent(In)           :: Reqs      !
  Integer,         Intent(In)           :: iField    ! Index of field.
  Integer,         Intent(In)           :: iT        ! Time index of field.
  Integer,         Intent(In)           :: iTA       ! Array index at which results for a
                                                     ! given time are stored.
  Integer,         Intent(In)           :: iLast     ! Single index for 7th dimension
  Type(Results_),  Intent(InOut)        :: Results   !
  ! Locals:
  Integer          :: iParticleSpecies ! Index of species used for mass-weighting.
  Integer          :: iHGrid           ! Index of horizontal grid.
  Integer          :: iZGrid           ! Index of vertical grid.
  Integer          :: iTGrid           ! Index of temporal grid.
  Integer          :: iX               ! X index of field.
  Integer          :: iY               ! Y index of field.
  Integer          :: iZ               ! Z index of field.
  Real(Std) :: X(2)
  Type(ShortTime_) :: OffsetTime  ! Offset time of particle relative to the output time.
  Real(Std)        :: TravelTime  ! Travel time of particle (as a real).
  Logical          :: VariableZGrid ! Indicates that a variable vertical grid is used.
  Type(HGrid_), Pointer :: HGrid ! Horizontal grid of particle numbers field.
  Type(ZGrid_), Pointer :: ZGrid ! Vertical grid of particle numbers field.

  iParticleSpecies = Reqs%FieldReqs(iField)%iParticleSpecies
  iHGrid           = Reqs%FieldReqs(iField)%iHGrid
  iZGrid           = Reqs%FieldReqs(iField)%iZGrid
  iTGrid           = Reqs%FieldReqs(iField)%iTGrid

  If (iHGrid /= 0) Then
    HGrid => Grids%HGrids(iHGrid)

    X(:) = Particle2%X(1:2)
    If (HGrid%XCycle > 0) Then ! $$ support for unstructured grids?
      Do While (X(1) < HGrid%XCentre - 0.5 * HGrid%XCycle)
        X(1) = X(1) + HGrid%XCycle
      End Do
      Do While (X(1) > HGrid%XCentre + 0.5 * HGrid%XCycle)
        X(1) = X(1) - HGrid%XCycle
      End Do
    End If
    If (HGrid%YCycle > 0) Then
      Do While (X(2) < HGrid%YCentre - 0.5 * HGrid%YCycle)
        X(2) = X(2) + HGrid%YCycle
      End Do
      Do While (X(2) > HGrid%YCentre + 0.5 * HGrid%YCycle)
        X(2) = X(2) - HGrid%YCycle
      End Do
    End If

    iX = NInt((X(1) - HGrid%X0) / HGrid%dX + 1.0)
    iY = NInt((X(2) - HGrid%Y0) / HGrid%dY + 1.0)
    If (iX > HGrid%nX .or. iX < 1 .or. iY > HGrid%nY .or. iY < 1) Return
  Else
    iX = 1
    iY = 1
  End If

  VariableZGrid = .false.
  If (iZGrid /= 0) Then
    ZGrid => Grids%ZGrids(iZGrid)
    If (ZGrid%Variable) Then
      VariableZGrid = .true.
    Else
      iZ = NInt((Particle2%X(3) - ZGrid%Z0) / ZGrid%dZ + 1.0)
      If (iZ > ZGrid%nZ .or. iZ < 1) Return
    End If
  Else If (Reqs%FieldReqs(iField)%AvBL) Then
    If (Particle2%X(3) <= Particle%H) Then
      iZ = 1
    Else
      Return
    End If
  Else
    iZ = 1
  End If

  OffsetTime = ParticleTime(Particle) - TInTGrid(Grids%TGrids(iTGrid), iT)
  TravelTime = ShortTime2RealTime(Particle%T - OffsetTime)

  If (VariableZGrid) Then
    Do iZ = 1, ZGrid%nZ
      If ((Particle2%X(3) >= ZGrid%Z(iZ) - 0.5 * ZGrid%AvZ(iZ)) .and. &
          (Particle2%X(3) <= ZGrid%Z(iZ) + 0.5 * ZGrid%AvZ(iZ))) Then
!$OMP ATOMIC
        Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iTA, iLast) =    &
          Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iTA, iLast) +  &
          Mass(iParticleSpecies) * TravelTime
      End If
    End Do
  Else
!$OMP ATOMIC
    Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iTA, iLast) =    &
      Results%Fields(iField)%Std(1, iX, iY, iZ, 1, iTA, iLast) +  &
      Mass(iParticleSpecies) * TravelTime
  End If

End Subroutine ParticleMeanTravelTime

!-------------------------------------------------------------------------------------------------------------

Subroutine ParticleMass(Particle, Particle2, Mass, Reqs, iField, iT, iLast, Results)
! Compute mean z.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In)    :: Particle  !
  Type(Particle_), Intent(In)    :: Particle2 !
  Real(Std),       Intent(In)    :: Mass(:)
  Type(Reqs_),     Intent(In)    :: Reqs      !
  Integer,         Intent(In)    :: iField    !
  Integer,         Intent(In)    :: iT        ! rename iT iTA $$
  Integer,         Intent(In)    :: iLast     ! Single index for 7th dimension
  Type(Results_),  Intent(InOut) :: Results   !
  ! Locals:
  Integer :: iParticleSpecies  ! Index of species used for mass-weighting.

  iParticleSpecies = Reqs%FieldReqs(iField)%iParticleSpecies

!$OMP ATOMIC
  Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, iLast) = Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, iLast) + &
                                                     Mass(iParticleSpecies)

End Subroutine ParticleMass

!-------------------------------------------------------------------------------------------------------------

Subroutine ParticleMeanZ(Particle, Particle2, Mass, Reqs, iField, iT, iLast, Results)
! Compute mean z.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In)    :: Particle  !
  Type(Particle_), Intent(In)    :: Particle2 !
  Real(Std),       Intent(In)    :: Mass(:)
  Type(Reqs_),     Intent(In)    :: Reqs      !
  Integer,         Intent(In)    :: iField    !
  Integer,         Intent(In)    :: iT        ! rename iT iTA $$
  Integer,         Intent(In)    :: iLast     ! Single index for 7th dimension
  Type(Results_),  Intent(InOut) :: Results   !
  ! Locals:
  Integer :: iParticleSpecies  ! Index of species used for mass-weighting.

  iParticleSpecies = Reqs%FieldReqs(iField)%iParticleSpecies

!$OMP ATOMIC
  Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, iLast) = Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, iLast) + &
                                                     Mass(iParticleSpecies) * Particle2%X(3)

End Subroutine ParticleMeanZ

!-------------------------------------------------------------------------------------------------------------

Subroutine ParticleSigZ2(Particle, Particle2, Mass, Reqs, iField, iT, iLast, Results)
! Compute sigma_z^2.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In)    :: Particle  !
  Type(Particle_), Intent(In)    :: Particle2 !
  Real(Std),       Intent(In)    :: Mass(:)
  Type(Reqs_),     Intent(In)    :: Reqs      !
  Integer,         Intent(In)    :: iField    !
  Integer,         Intent(In)    :: iT        ! rename iT iTA $$
  Integer,         Intent(In)    :: iLast     ! Single index for 7th dimension
  Type(Results_),  Intent(InOut) :: Results   !
  ! Locals:
  Integer :: iParticleSpecies  ! Index of species used for mass-weighting.

  iParticleSpecies = Reqs%FieldReqs(iField)%iParticleSpecies

!$OMP ATOMIC
  Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, iLast) = Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, iLast) + &
                                                     Mass(iParticleSpecies) * Particle2%X(3)**2

End Subroutine ParticleSigZ2

!-------------------------------------------------------------------------------------------------------------

Subroutine ParticleXStats(Particle, Particle2, Mass, Coords, Grids, Reqs, iField, iS, Results)
! Calculate contribution to X Stats from particle.

  Implicit None
  ! Argument list:
  Type(Particle_), Intent(In)    :: Particle  ! Single particle under consideration.
  Type(Particle_), Intent(In)    :: Particle2 ! Single particle under consideration.
  Real(Std),       Intent(In)    :: Mass(:)
  Type(Coords_),   Intent(In)    :: Coords    !
  Type(Grids_),    Intent(In)    :: Grids     !
  Type(Reqs_),     Intent(In)    :: Reqs      !
  Integer,         Intent(In)    :: iField    ! Index of field.
  Integer,         Intent(In)    :: iS        ! Travel time index of X Stats field.
  Type(Results_),  Intent(InOut) :: Results   !
  ! Locals:
  Type(ShortTime_) :: ParticleT        ! Particle time.
  Integer          :: i                ! Loop index.
  Integer          :: iTGrid           ! Index of temporal grid.
  Integer          :: iT               ! Time index of X Stats field.
  Integer          :: iTA              ! Array index at which results for a given time are
                                       ! stored.
  Integer          :: iParticleSpecies ! Index of species used for mass-weighting.
  Real(Std)        :: XStats(10)       ! X Stats contribution from particle.

  ! Calculate time index of X Stats field nearest to the current particle time.
  !
  ! Now changed to next grid time (DJT, 5/4/05). Needs further thought though - e.g.
  ! shouldn't contribute to first grid time if particle time <= T(1) - dT. $$
  ! If symmetric time av wanted, could treat handle in processing fields. $$
  ParticleT = ParticleTime(Particle)
  iTGrid    = Reqs%FieldReqs(iField)%iTGrid
  iT        = 1
  If (iTGrid /= 0) Then
    Do
      If (iT >= Grids%TGrids(iTGrid)%nT) Exit
!      If (                                                                    &
!        Abs(ShortTime2RealTime(ParticleT - Grids%TGrids(iTGrid)%T(iT    ))) < &
!        Abs(ShortTime2RealTime(ParticleT - Grids%TGrids(iTGrid)%T(iT + 1)))   &
!      ) Exit
      If (ParticleT <= TInTGrid(Grids%TGrids(iTGrid), iT)) Exit
      iT = iT + 1
    End Do
  End If

  iParticleSpecies = Reqs%FieldReqs(iField)%iParticleSpecies

  ! Compute the contribution to X Stats.
  XStats(:) = (/                                                            & ! $$ Problem with wrapped coords?
                 Mass(iParticleSpecies),                                    &
                 Mass(iParticleSpecies) * Particle2%X(1),                   &
                 Mass(iParticleSpecies) * Particle2%X(2),                   &
                 Mass(iParticleSpecies) * Particle2%X(3),                   &
                 Mass(iParticleSpecies) * Particle2%X(1)**2,                &
                 Mass(iParticleSpecies) * Particle2%X(2)**2,                &
                 Mass(iParticleSpecies) * Particle2%X(3)**2,                &
                 Mass(iParticleSpecies) * Particle2%X(1) * Particle2%X(2),  &
                 Mass(iParticleSpecies) * Particle2%X(2) * Particle2%X(3),  &
                 Mass(iParticleSpecies) * Particle2%X(3) * Particle2%X(1)   &
              /)

  iTA = CalciTA(Results%Fields(iField), iT) ! $$ If this fails, may need to store more field-times.

!$OMP CRITICAL (PTCLXSTATS)
  Results%Fields(iField)%Std(1, 1, 1, 1, iS, iTA, :) = Results%Fields(iField)%Std(1, 1, 1, 1, iS, iTA, :) + &
                                                       XStats(:)
!$OMP END CRITICAL (PTCLXSTATS)
End Subroutine ParticleXStats

!-------------------------------------------------------------------------------------------------------------

Subroutine ParticleInfo(                                     &
             iCase,                                          &
             Particle, Particle2, Extra, Mass,               &
             Flow, Cloud, Rain,                              &
             OutputOpts,                                     &
             Coords, Grids, Flows, Specieses, Sources, Reqs, &
             iPPInfo, iT,                                    &
             Units, Results                                  &
           )
! Outputs particle info and infomation on flow properties etc which have been used
! to time-step the particle.

  Implicit None
  ! Argument list:
  Integer,           Intent(In)    :: iCase      ! Number of case.
  Type(Particle_),   Intent(In)    :: Particle   !
  Type(Particle_),   Intent(In)    :: Particle2  !
  Type(Extra_),      Intent(In)    :: Extra      !
  Real(Std),         Intent(In)    :: Mass(:)
  Type(Flow_),       Intent(In)    :: Flow       !
  Type(Cloud_),      Intent(In)    :: Cloud      !
  Type(Rain_),       Intent(In)    :: Rain       !
  Type(OutputOpts_), Intent(In)    :: OutputOpts ! Output options.
  Type(Coords_),     Intent(In)    :: Coords
  Type(Grids_),      Intent(In)    :: Grids      !
  Type(Flows_),      Intent(In)    :: Flows
  Type(Specieses_),  Intent(In)    :: Specieses  !
  Type(Sources_),    Intent(In)    :: Sources    !
  Type(Reqs_),       Intent(In)    :: Reqs       !
  Integer,           Intent(In)    :: iPPInfo    !
  Integer,           Intent(In)    :: iT         !
  Type(Units_),      Intent(InOut) :: Units      ! Collection of information on input/output unit numbers.
  Type(Results_),    Intent(InOut) :: Results    !
  ! Locals:
  Integer                        :: iTGrid
  Type(Time_)                    :: T         ! Current time.
  Type(Time_)                    :: T0        ! Release time.
  Character(MaxOutputLineLength) :: Line
  Character(8)                   :: Skew
  Character(8)                   :: VelMem
  Character(8)                   :: Inhomog
  Integer                        :: i
  Real(Std) :: U(3), SigUU(3) ! $$ for temp fix up $$
  Real(Std) :: RH             ! Relative humidity (%).
  Real(Std) :: CloudOktas     ! Cloud Cover in Oktas.
  Real(Std) :: WindSpeed      ! Wind Speed.
  Real(Std) :: WindDirection  ! Wind Direction (degrees).

  ! Index of TGrid.
  iTGrid = Reqs%PPInfoReqs(iPPInfo)%iTGrid

  ! Numerical screen and disk output.
  If (Reqs%PPInfoReqs(iPPInfo)%Screen .or. &
      Reqs%PPInfoReqs(iPPInfo)%Disk) Then

    If (iTGrid == 0) Then
      T = ShortTime2Time(Particle%T0 + Particle%T)
    Else
      T = ShortTime2Time(TInTGrid(Grids%TGrids(iTGrid), iT))
    End If
    T0 = ShortTime2Time(Particle%T0)

    ! Construct output line.

    ! $$ Temp fix up for Extra%U, Extra%SigUU undefined (note the equiv puff code seems OK)
    If (Extra%VelMem) Then
      U(:)     = Flow%U(:)
      SigUU(:) = Flow%SigUU(:)
    Else
      U(:)     = 0.0
      SigUU(:) = 0.0
    End If

    ! 1. Basic information. (Output travel time as a time? $$)
    Line =                                                                           &
        Trim(FormatChar('N',                                    15, .true., 'R')) // &
        Trim(  Int2Char(Particle%iUP,                           15, .true., 'R')) // &
        Trim(FormatChar(Sources%Sources(Particle%iSource)%Name, 15, .true., 'R')) // &
        Trim(FormatChar(Time2Char(T0, .true., 0, .false.),      25, .true., 'R')) // &
        Trim(FormatChar(Time2Char(T,  .true., 0, .false.),      25, .true., 'R')) // &
        Trim(  Std2Char(ShortTime2RealTime(Particle%T),         15, .true., 'R')) // &
        Trim(  Std2Char(Particle2%X(1),                         15, .true., 'R')) // &
        Trim(  Std2Char(Particle2%X(2),                         15, .true., 'R')) // &
        Trim(  Std2Char(Particle2%X(3),                         15, .true., 'R')) // &
        Trim(  Std2Char(U(1),                                   15, .true., 'R')) // &
        Trim(  Std2Char(U(2),                                   15, .true., 'R')) // &
        Trim(  Std2Char(U(3),                                   15, .true., 'R')) // &
        Trim(FormatChar('-',                                    15, .true., 'R')) // &
        Trim(FormatChar('-',                                    15, .true., 'R')) // &
        Trim(FormatChar('-',                                    15, .true., 'R'))

    ! 2. Met information. $$ Currently from particle (i.e. previous time step) - use Flow instead?
    If (Reqs%PPInfoReqs(iPPInfo)%Met) Then

      ! Limit relative humidity.
      RH = CalcRH(Flow%Q, Flow%T, Flow%P)
      If (RH <   0.0) RH =   0.0
      If (RH > 100.0) RH = 100.0

      CloudOktas = Cloud%Cloud * 8.0

      ! Calculate wind speed and direction at particle location
      WindSpeed = Sqrt(Flow%U(1)**2 + Flow%U(2)**2)

      WindDirection = ATan2ZeroTest(-Flow%U(1), -Flow%U(2)) * 180.0 / Pi
      If (WindDirection < 0.0) WindDirection = WindDirection + 360.0

      ! $$ Issues:
      ! RH < 0, > 100 tests shouldn't be needed if everything defined. (applies elsewhere in code too).
      ! Implement wind direction convention for calms/northerlies = 0/360.

      Line = Trim(Line)                                       // & ! $$ check old values
             Trim(Std2Char(Flow%U(1),       15, .true., 'R')) // & ! always defined.
             Trim(Std2Char(Flow%U(2),       15, .true., 'R')) // &
             Trim(Std2Char(Flow%U(3),       15, .true., 'R')) // &
             Trim(Std2Char(SigUU(1),        15, .true., 'R')) // &
             Trim(Std2Char(SigUU(2),        15, .true., 'R')) // &
             Trim(Std2Char(SigUU(3),        15, .true., 'R')) // &
             Trim(Std2Char(Flow%T,          15, .true., 'R')) // &
             Trim(Std2Char(Flow%P,          15, .true., 'R')) // &
             Trim(Std2Char(Flow%Theta,      15, .true., 'R')) // &
             Trim(Std2Char(Particle%H,      15, .true., 'R')) // &
             Trim(Std2Char(CloudOktas,      15, .true., 'R')) // &
             Trim(Std2Char(RH,              15, .true., 'R')) // &
             Trim(Std2Char(WindSpeed,       15, .true., 'R')) // &
             Trim(Std2Char(WindDirection,   15, .true., 'R'))

    End If

    ! 3. Mass information
    If (Reqs%PPInfoReqs(iPPInfo)%Mass) Then
      Do i = 1, Specieses%nParticleSpecieses
        Line = Trim(Line) // Std2Char(Mass(i), 15, .true., 'R')
      End Do
    End If

    ! 4. Plume rise information. $$ Need to check what these values are actually?
    If (Reqs%PPInfoReqs(iPPInfo)%PlumeRise) Then
      Line = Trim(Line)                                      // &
             Trim(Std2Char(Extra%FMass,    15, .true., 'R')) // &
             Trim(Std2Char(Extra%FMass0,   15, .true., 'R')) // &
             Trim(Std2Char(Extra%FM(1),    15, .true., 'R')) // &
             Trim(Std2Char(Extra%FM(2),    15, .true., 'R')) // &
             Trim(Std2Char(Extra%FM(3),    15, .true., 'R')) // &
             Trim(Std2Char(Extra%FH,       15, .true., 'R'))
    End If

    ! 5. Dispersion scheme information.
    If (Reqs%PPInfoReqs(iPPInfo)%DispersionScheme) Then
      Write (Skew,    *) Extra%Skew
      Write (VelMem,  *) Extra%VelMem
      Write (Inhomog, *) Extra%Inhomog
      Line = Trim(Line)                                 // &
             Trim(FormatChar(Skew,    15, .true., 'R')) // &
             Trim(FormatChar(VelMem,  15, .true., 'R')) // &
             Trim(FormatChar(Inhomog, 15, .true., 'R'))
    End If

    ! $$ need something for puff family information (e.g. ' -  , -  , -  , -') ?

! $$ Worth including:
!   If (Particle%FMASS /= 0.0) Then
!     Line =
!       Trim(Line)                                                                 // &
!       Trim(Std2Char((Particle%FM(1)/Particle%FMASS+Flow%U(1)), 15, .true., 'R')) // &
!       Trim(Std2Char((Particle%FM(2)/Particle%FMASS+Flow%U(2)), 15, .true., 'R')) // &
!       Trim(Std2Char((Particle%FM(3)/Particle%FMASS+Flow%U(3)), 15, .true., 'R')) // &
!       Trim(Std2Char((Particle%FH/(1004.0*Particle%FMASS)+Flow%Theta), 15, .true., 'R'))
!   End If
! Headers would be something like:
! Header =                                              &
!     Trim(FormatChar('UP',       15, .true.,  'R')) // &
!     Trim(FormatChar('VP',       15, .true.,  'R')) // &
!     Trim(FormatChar('WP',       15, .true.,  'R')) // &
!     Trim(FormatChar('THETAP',   15, .true.,  'R'))

  !$OMP CRITICAL(OutputPPInfosCritical)
    Call OutputPPInfos(                                           &
           Line,                                                  &
           Particle2%X, Particle%XOld, Particle%T, Particle%TOld, & ! $$ coords for XOld may be suspect
           iCase, iPPInfo, iT, Particle%iUP, .false.,             &
           OutputOpts,                                            &
           Coords, Grids, Flows, Specieses, Sources, Reqs,        &
           Units, Results                                         &
         )
  !$OMP END CRITICAL(OutputPPInfosCritical)

  End If

End Subroutine ParticleInfo

!-------------------------------------------------------------------------------------------------------------

End Module ParticleModule

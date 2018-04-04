! Module:  Source Module

Module SourceModule

! This module provides code for handling sources.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowsModule, Only: Flow_, Surface_, Soil_, Rain_
Use SizeDistModule
Use SpeciesModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: SourceTypeNames                 ! Names of met dependent source types.
Public  :: S_Generic                       ! Code for generic source type (no met dependence).
Public  :: S_Dust                          ! Code for dust source type.
Public  :: S_SeaSalt                       ! Code for sea-salt source type.
Public  :: S_IterativePlumeModel           ! Code for iterative plume-rise model
Public  :: Source_                         ! Information defining a source.
Public  :: SourceGroup_                    ! Information defining a source group.
Public  :: Sources_                        ! A collection of sources.
Public  :: TimeDep_                        ! A source time dependency.
Public  :: TimeDeps_                       ! A collection of source time dependencies.
Public  :: InitSources                     ! Initialises a collection of sources.
Public  :: InitSource                      ! Initialises a source.
Public  :: AddSource                       ! Adds a source to a collection of sources.
Public  :: FindSourceIndex                 ! Finds the index of a source in a collection of sources.
Public  :: FindSourceGroupIndex            ! Finds the index of a source group in a collection of sources.
Public  :: InitTimeDeps                    ! Initialises a collection of source time dependencies.
Public  :: InitTimeDep                     ! Initialises a source time dependency.
Public  :: AddTimeDep                      ! Adds a source time dependency to a collection of source time
                                           ! dependencies.
Public  :: FindTimeDepIndex                ! Finds the index of a source time dependency in a collection of
                                           ! source time dependencies.
Public  :: SetUpSources                    ! Sets up Sources.
Public  :: SetUpiSources                   ! Sets up indices in Sources.
Public  :: FirstLastReleaseTimes           ! Determines the first source emission time and last possible
                                           ! source emission time for all sources, and the first release time
                                           ! and the last possible release time for each source and for all
                                           ! sources (accounting for the source information, the computational
                                           ! domain's time structure, and whether the met is fixed). Note that
                                           ! here we distinguish between source emission times and release
                                           ! times. A release time is the actual time of release of a
                                           ! particle/puff in the model and may correspond to a range of
                                           ! source emission times if the particle/puff has a time spread.
Public  :: ParticleReleaseInfo             ! Calculates, for a particle release, the mass to release (for each
                                           ! species), the number of particles to release, and the next
                                           ! release time.
Public  :: PuffReleaseInfo                 ! Calculates, for a puff release, the mass to release (for each
                                           ! species).
Public  :: SourcePositionInStandardLatLong ! Returns the location of the named source in the standard lat-long
                                           ! coord system.

!-------------------------------------------------------------------------------------------------------------

Character(MaxCharLength), Parameter :: SourceTypeNames(3) =           &
                                       (/                             &
                                         'Dust    ',                  & ! 1
                                         'Sea Salt',                  & ! 2
                                         'Iterative Plume Rise Model' & ! 3
                                       /)
! SourceTypeNames :: Names of met dependent source types.

! Codes for met dependent source types.
Integer, Parameter :: S_Generic             = 0 ! Code for generic source type (no met dependence).
Integer, Parameter :: S_Dust                = 1 ! Code for dust source type.
Integer, Parameter :: S_SeaSalt             = 2 ! Code for sea-salt source type.
Integer, Parameter :: S_IterativePlumeModel = 3 ! Code for iterative plume-rise model.

!-------------------------------------------------------------------------------------------------------------

Type :: Source_ ! A source.
  Character(MaxCharLength) :: Name
  Integer                  :: nSourceGroups
  Character(MaxCharLength) :: SourceGroupNames(MaxSourceGroupsPerSource)
  Integer                  :: iSourceGroups(MaxSourceGroupsPerSource)
  Integer                  :: SourceType
  Character(MaxCharLength) :: HGridName
  Character(MaxCharLength) :: ZGridName
  Integer                  :: iHGrid
  Integer                  :: iZGrid
  Character(MaxCharLength) :: HCoordName
  Character(MaxCharLength) :: ZCoordName
  Character(MaxCharLength) :: dHCoordName ! $$ currently dHCoordName and dZCoordName
  Character(MaxCharLength) :: dZCoordName ! $$ must equal HCoordName and ZCoordName
  Logical                  :: dHMetres
  Logical                  :: dZMetres
  Integer                  :: iHCoord
  Integer                  :: iZCoord
  Integer                  :: idHCoord
  Integer                  :: idZCoord
  Real(Std)                :: X(3)
  Real(Std)                :: dX(3)
  Real(Std)                :: Angle
  Character(MaxCharLength) :: Shape
  Logical                  :: TopHat
  Logical                  :: UniformArea
  Logical                  :: NoReflect
  Type(Time_)              :: StartTime
  Type(Time_)              :: StopTime
  Type(Time_)              :: MaxTravelTime
  Type(ShortTime_)         :: sStartTime
  Type(ShortTime_)         :: sStopTime
  Type(ShortTime_)         :: sMaxTravelTime
  Real(Std)                :: nParticles
  Logical                  :: ParticleRate
  Integer                  :: nSpecieses
  Character(MaxCharLength) :: SpeciesNames(MaxSpecieses)
  Integer                  :: iSpecieses(MaxSpecieses)
  Integer                  :: iParticleSpecieses(MaxSpecieses)
  Real(Std)                :: SourceStrengths(MaxSpecieses)
  Logical                  :: SourceRate(MaxSpecieses)
  Character(MaxCharLength) :: MaterialUnitName(MaxSpecieses)
  Integer                  :: iMaterialUnit(MaxSpecieses)
  Integer                  :: iCritSpecies
  Logical                  :: PlumeRise
  Real(Std)                :: VolumeFlowRate ! $$ Issues of coord system and area calculation.
  Real(Std)                :: FlowVelocity   ! Vector for jets? $$ Do we need both variables?
  Real(Std)                :: Temperature  ! Do we need molecular weight?
  Logical                  :: Particulate
  Real(Std)                :: Density
  Real(Std)                :: Diameter
  Real(Std)                :: ParticleShape
  Logical                  :: DensityPresent
  Logical                  :: ParticleShapePresent
  Character(MaxCharLength) :: ShapeScheme
  Integer                  :: ShapeSchemeCode
  Character(MaxCharLength) :: SizeDistName
  Integer                  :: iSizeDist
  Character(MaxCharLength) :: TimeDepName
  Integer                  :: iTimeDep
  Logical                  :: FixedMet
  Integer                  :: iBioIsoprene
  Integer                  :: iBioAPinene
  Integer                  :: MetDependent
  Type(ShortTime_)         :: EulerianTime
  Character(MaxCharLength) :: ZGridIterativePlumeRise
  Integer                  :: iZGridIterativePlumeRise
  Real(Std)                :: WaterVapourFraction0
  Real(Std)                :: GasFraction0
  Real(Std)                :: PlumeRiseHeight
  Real(Std)                :: SummitElevation
  Real(Std)                :: DistalFineAshFraction
  ! Name               :: Name of source.
  ! SourceGroupNames   :: Names of source groups to which the source belongs.
  ! iSourceGroups      :: Indices of source groups to which the source belongs.
  ! SourceType         :: Code for met dependent source type.
  ! HGridName          :} Names of horizontal and vertical grids (for horizontally and vertically gridded
  ! ZGridName          :} sources).
  ! iHGrid             :] Indices of horizontal and vertical grids.
  ! iZGrid             :]
  ! HCoordName         :} Name of horizontal and vertical coord systems used for defining
  ! ZCoordName         :} the source centre X.
  ! dHCoordName        :] Name of horizontal and vertical coord systems used for defining
  ! dZCoordName        :] source shape and extent - see Shape, TopHat, UniformArea,
  ! dHMetres           :] NoReflect, dX and Angle. If dHMetres and/or dZMetres are true
  ! dZMetres           :] then the coord systems are rescaled so as to be in metres at
  !                       the source centre.
  ! iHCoord            :} Indices of coord systems.
  ! iZCoord            :}
  ! idHCoord           :}
  ! idZCoord           :}
  ! X                  :: Source centre in the coord system defined by HCoordName and
  !                       ZCoordName.
  ! dX                 :: dX is the extent of the source in the directions of the
  !                       principle axes (in the coord system defined by dHCoordName,
  !                       dZCoordName, dHMetres and dZMetres), with the third principle
  !                       axis assumed vertical.
  ! Angle              :: Angle of the first principle axis (dX(1)) to the X axis
  !                       (measured anticlockwise in radians).
  ! Shape              :: Source shape, which can be 'cuboid' or 'ellipsoid'. The source
  !                       is assumed uniformly distributed within its extent (unless
  !                       TopHat is false). If shape is ellipsoid and one or two of dX,
  !                       dY and dZ are zero then source is uniformly distributed over a
  !                       disk or line segment. Note 'cuboid', 'ellipsoid' and 'uniform'
  !                       mean cuboid, ellipsoid and uniform in the coord system defined
  !                       by dHCoordName, dZCoordName, dHMetres and dZMetres. This will
  !                       not be truly cuboid, ellipsoid or uniform to the extent that
  !                       the coord system is not cartesian (but see uniform area).
  ! TopHat             :: If false, the source is approximated by a Gaussian distribution
  !                       with the same covariance matrix in the coord system defined by
  !                       dHCoordName, dZCoordName, dHMetres and dZMetres.
  ! UniformArea        :: For TopHat true, indicates that the source is uniform across
  !                       its horizontal area (as opposed to being uniform in terms of
  !                       the area dXdY defined by the coord system dHCoordName). For
  !                       TopHat false, indicates that the Gaussian distribution is
  !                       multiplied by the ratio of the local true area to the
  !                       coordinate area (and then normalised to integrate to 1).
  ! NoReflect          :: Indicates no reflection of the source at the ground, with any
  !                       part of the source below the ground being lost.
  ! StartTime          :: Start time of the source.
  ! StopTime           :: End time of the source.
  ! MaxTravelTime      :: Maximum travel time for which particles are to be followed.
  ! sStartTime         :} Short time versions of StartTime, StopTime and MaxTravelTime.
  ! sStopTime          :}
  ! sMaxTravelTime     :}
  ! nParticles         :: Number of particles released or particle release rate.
  ! ParticleRate       :: Indicates that nParticles is the particle release rate, not the
  !                    :: number of particles released. Is false for fixed met and
  !                       instantaneous sources and true otherwise.
  ! nSpecieses         :: Number of species.
  ! SpeciesNames       :: Species names.
  ! iSpecieses         :: Indices of species in set of all species.
  ! iParticleSpecieses :: Indices of species in set of species on particles.
  ! SourceStrengths    :: Total mass released or mass release rate for each species.
  ! SourceRate         :: Indicates that SourceStrength is the mass release rate, not the
  !                       total mass released. Is false for instantaneous sources and
  !                       true otherwise.
  ! iMaterialUnit      :: Unit in which the quantity of material is measured.
  ! iCritSpecies       :: Index (in the arrays SpeciesNames etc) of the species which,
  !                       through the MassLimit of the species, will limit the mass of
  !                       the particles. 0 indicates no limit due to the MassLimits of
  !                       the specieses.
  ! PlumeRise          :: Indicates plume rise is to be modelled.
  ! VolumeFlowRate     :: Volume flow rate of the source. } Defined only if PlumeRise is
  ! FlowVelocity       :: Flow velocity.                  } true. The plume rise
  ! Temperature        :: Temperature of emission (K).    } parameters  refer to the
  !                                                         total emission, not just the
  !                                                         species of interest.
  ! Particulate          :: Indicates that the species are emitted as particulates,
  !                         particulates being anything with significant fall speed and/or
  !                         inertia. If some species are emitted as particulates and some
  !                         are not, or if some species have different particle density,
  !                         particle diameter or particle size distribution, then they need
  !                         to be treated in different sources.
  ! Diameter             :: Particle diameter. Defined only if          }
  !                         SizeDistName is blank.                      } 
  ! Density              :: Particle density.                           } Defined only if
  ! ParticleShape        :: Particle shape                              } Particulate is
  ! DensityPresent       :: Indicates whether particle density is given } true
  ! ParticleShapePresent :: Indicates whether ParticleShape is given    } 
  ! ShapeScheme          :: Shape Scheme                                } 
  ! ShapeSchemeCode      :: Shape Scheme code                           } 
  ! SizeDistName         :: Name of particle size distribution. Blank   }
  !                         indicates no particle size distribution.    }
  ! iSizeDist            :: Index of particle size distribution. 0      }
  !                         indicates no particle size distribution.    }
  ! TimeDepName          :: Name of source time dependency. Blank indicates no source time
  !                         dependency.
  ! iTimeDep             :: Index of source time dependency. 0 indicates no source time
  !                         dependency.
  ! FixedMet             :: Indicates that the met is fixed.
  ! iBioIsoprene         :: Index of biogenic isoprene in list of source species.
  ! iBioAPinene          :: Index of biogenic a-pinene in list of source species.
  ! MetDependent         :: 1 = source with a met dependency, 0 = no met dependency.
  ! EulerianTime         :: Time at which particles' mass is transferred to Eulerian field.
  ! ZGridIterativePlumeRise  :: Name of Z-grid for iterative plume-rise model
  ! iZGridIterativePlumeRise :: Index of  Z-grid for iterative plume-rise model
  ! WaterVapourFraction0     :: Initial mass fraction of water vapour for iterative plume-rise model
  ! GasFraction0             :: Total initial mass fraction of gas for iterative plume-rise model
  ! PlumeRiseHeight          :: Observed rise height. This is required when applying the iterative
  !                             plume-rise model to gas-only sources. For volcanic ash it should
  !                             ordinarily equal the source height
  ! SummitElevation          :: Height of source (above ground/sea level)
  ! DistalFineAshFraction    :: Fraction of ash that survives into the far field

  ! Note StartTime, StopTime, MaxTravelTime, nParticles, ParticleRate,
  ! SourceStrengths, SourceRate and the source time dependency do not need to ensure
  ! that anything is released. This provides a variety of ways to turn off a source
  ! without deleting it completely.

  ! ParticleRate and SourceRate are false for instantaneous sources, but are true
  ! otherwise, with the exception of fixed met cases where ParticleRate is always
  ! false.

  ! $$ Note potential problem here with interacting species e.g. condensing plumes
  ! with water vapour and particulates. Solution may depend on modelling approach
  ! adopted, but could allow extra plume rise parameters (water vapour etc) rather
  ! than putting water vapour in as a species. If want water vapour and particulate
  ! concentrations then need to treat two sources, both with water vapour affecting
  ! plume rise, but only one having water vapour as a species. Phase changes are an
  ! extra complication, best treated as a species change.

End Type Source_

!-------------------------------------------------------------------------------------------------------------

Type :: SourceGroup_ ! A source group.
  Character(MaxCharLength) :: Name ! Name of source group.
  ! $$ This has been made a derived type to make it easier in the future for each source group to store the
  !    source names. This would make it easy to add the ability to have a source group input block.
End Type SourceGroup_

!-------------------------------------------------------------------------------------------------------------

Type :: Sources_ ! A collection of sources.
  Integer                    :: nSources                      ! Number of sources.
  Type(Source_),     Pointer :: Sources(:)                    ! Collection of sources.
  Integer                    :: nSourceGroups                 ! Number of source groups.
  Type(SourceGroup_)         :: SourceGroups(MaxSourceGroups) ! Collection of source groups.
End Type Sources_

!-------------------------------------------------------------------------------------------------------------

Type :: TimeDep_ ! A source time dependency.
  Character(MaxCharLength) :: Name
  Integer                  :: nFactors
  Real(Std)                :: Factors(MaxTimeDepFactors)
  Integer                  :: nWildTimePairs(MaxTimeDepFactors)
  Type(WildTimePair_)      :: WildTimePairs(MaxTimeDepPairs, MaxTimeDepFactors)
  ! Name           :: Name of source time dependency.
  ! nFactors       :: Number of source time dependency factors in source time
  !                   dependency.
  ! Factors        :: Source time dependency factors.
  ! nWildTimePairs :: Number of wild-card time pairs used with each source time
  !                   dependency factor.
  ! WildTimePairs  :: Wild-card time pairs.

  ! If the time lies, for each j in (1, nWildTimePairs(i)), in the time intervals
  ! defined by the wild-card time pair WildTimePairs(j, i), then the source strength
  ! is multiplied by Factors(i). This is done for each i, so several factors could be
  ! applied at the same time.
End Type TimeDep_

!-------------------------------------------------------------------------------------------------------------

Type :: TimeDeps_ ! A collection of source time dependencies.
  Integer        :: nTimeDeps             ! Number of source time dependencies.
  Type(TimeDep_) :: TimeDeps(MaxTimeDeps) ! Collection of source time dependencies.
End Type TimeDeps_

!-------------------------------------------------------------------------------------------------------------

Interface Operator(==) ! Equality of source time dependencies.
  Module Procedure TimeDepEq
End Interface

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitSources(MaxSources) Result (Sources)
! Initialises a collection of sources.

  Implicit None
  ! Argument list:
  Integer, Intent(In) :: MaxSources ! Maximum number of sources.
  ! Function result:
  Type(Sources_) :: Sources ! Initialised collection of sources.
  ! Locals:
  Integer :: ErrorCode ! Error code.

  Allocate(                      &
    Sources%Sources(MaxSources), &
    Stat = ErrorCode             &
  )
  If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate array for sources', 3)

  Sources%nSources      = 0
  Sources%nSourceGroups = 0

End Function InitSources

!-------------------------------------------------------------------------------------------------------------

Function InitSource(                                                             &
           Name,                                                                 &
           SourceGroupNames,                                                     &
           SourceType, HGridName, ZGridName,                                     &
           HCoordName, ZCoordName, dHCoordName, dZCoordName, dHMetres, dZMetres, &
           X, dX, Angle, Shape, TopHat, UniformArea, NoReflect,                  &
           StartTime, StopTime, MaxTravelTime,                                   &
           CharNParticles, CharSourceStrengths,                                  &
           PlumeRise, VolumeFlowRate, FlowVelocity, Temperature,                 &
           Particulate, Diameter, Density, ParticleShape,                        &
           DensityPresent, ParticleShapePresent, ShapeScheme, SizeDistName,      &
           TimeDepName,                                                          &
           FixedMet,                                                             &
           EulerianTime,                                                         &
           ZGridIterativePlumeRise,                                              &
           DryGasFraction0,                                                      &
           WaterVapourFraction0,                                                 &
           PlumeRiseHeight,                                                      &
           SummitElevation,                                                      &
           DistalFineAshFraction                                                 &
         )                                                                       &
Result (Source)
! Initialises a source.

! Note the checks on StartTime, StopTime, MaxTravelTime, nParticles, ParticleRate,
! SourceStrengths, SourceRate and the source time dependency do not ensure that
! anything is released. This provides a variety of ways to turn off a source without
! deleting it completely.

! $$ should VolumeFlowRate and FlowVelocity both be needed?

  Implicit None
  ! Argument list:
  Character(*),             Intent(In) :: Name
  Character(*),             Intent(In) :: SourceGroupNames(:)
  Integer,                  Intent(In) :: SourceType
  Character(*),             Intent(In) :: HGridName
  Character(*),             Intent(In) :: ZGridName
  Character(*),             Intent(In) :: HCoordName
  Character(*),             Intent(In) :: ZCoordName
  Character(*),             Intent(In) :: dHCoordName
  Character(*),             Intent(In) :: dZCoordName
  Logical,                  Intent(In) :: dHMetres    ! $$ add option of dHCoord being set up
  Logical,                  Intent(In) :: dZMetres    ! internally as stereographic projection at
                                                      ! source location
  Real(Std),                Intent(In) :: X(3)
  Real(Std),                Intent(In) :: dX(3)
  Real(Std),                Intent(In) :: Angle
  Character(*),             Intent(In) :: Shape
  Logical,                  Intent(In) :: TopHat
  Logical,                  Intent(In) :: UniformArea
  Logical,                  Intent(In) :: NoReflect
  Type(Time_),              Intent(In) :: StartTime
  Type(Time_),              Intent(In) :: StopTime
  Type(Time_),              Intent(In) :: MaxTravelTime
  Character(*),             Intent(In) :: CharNParticles
  Character(*),             Intent(In) :: CharSourceStrengths(:)
  Logical,                  Intent(In) :: PlumeRise
  Real(Std),                Intent(In) :: VolumeFlowRate
  Real(Std),                Intent(In) :: FlowVelocity
  Real(Std),                Intent(In) :: Temperature
  Logical,                  Intent(In) :: Particulate
  Real(Std),                Intent(In) :: Diameter
  Real(Std),                Intent(In) :: Density
  Real(Std),                Intent(In) :: ParticleShape
  Logical,                  Intent(In) :: DensityPresent
  Logical,                  Intent(In) :: ParticleShapePresent
  Character(MaxCharLength), Intent(In) :: ShapeScheme
  Character(*),             Intent(In) :: SizeDistName
  Character(*),             Intent(In) :: TimeDepName
  Logical,                  Intent(In) :: FixedMet
  Type(ShortTime_),         Intent(In) :: EulerianTime
  Character(*),             Intent(In) :: ZGridIterativePlumeRise
  Real(Std),                Intent(In) :: DryGasFraction0
  Real(Std),                Intent(In) :: WaterVapourFraction0
  Real(Std),                Intent(In) :: PlumeRiseHeight
  Real(Std),                Intent(In) :: SummitElevation
  Real(Std),                Intent(In) :: DistalFineAshFraction
  ! Name                :: Name of source.
  ! SourceGroupNames    :: Names of source groups to which the source belongs.
  ! SourceType          :: Code for met dependent source type.
  ! HGridName           :: Names of horizontal and vertical grids (for horizontally and vertically gridded
  ! ZGridName           :: sources).
  ! HCoordName          :} Name of horizontal and vertical coord systems used for
  ! ZCoordName          :} defining the source centre X.
  ! dHCoordName         :] Name of horizontal and vertical coord systems used for
  ! dZCoordName         :] defining source shape and extent - see Shape, TopHat,
  ! dHMetres            :] UniformArea, NoReflect, dX and Angle. If dHMetres and/or
  ! dZMetres            :] dZMetres are true then the coord systems are rescaled so as
  !                        to be in metres at the source centre.
  ! X                   :: Source centre in the coord system defined by HCoordName and
  !                        ZCoordName.
  ! dX                  :: dX is the extent of the source in the directions of the
  !                        principle axes (in the coord system defined by dHCoordName,
  !                        dZCoordName, dHMetres and dZMetres), with the third
  !                        principle axis assumed vertical.
  ! Angle               :: Angle of the first principle axis (dX(1)) to the X axis
  !                        (measured anticlockwise in degrees).
  ! Shape               :: Source shape, which can be 'cuboid' or 'ellipsoid'. The
  !                        source is assumed uniformly distributed within its extent
  !                        (unless TopHat is false). If shape is ellipsoid and one or
  !                        two of dX, dY and dZ are zero then source is uniformly
  !                        distributed over a disk or line segment. Note 'cuboid',
  !                        'ellipsoid' and 'uniform' mean cuboid, ellipsoid and
  !                        uniform in the coord system defined by dHCoordName,
  !                        dZCoordName, dHMetres and dZMetres. This will not be truly
  !                        cuboid, ellipsoid or uniform to the extent that the coord
  !                        system is not cartesian (but see uniform area).
  ! TopHat              :: If false, the source is approximated by a Gaussian
  !                        distribution with the same covariance matrix in the coord
  !                        system defined by dHCoordName, dZCoordName, dHMetres and
  !                        dZMetres.
  ! UniformArea         :: For TopHat true, indicates that the source is uniform
  !                        across its horizontal area (as opposed to being uniform in
  !                        terms of the area dXdY defined by the coord system
  !                        dHCoordName). For TopHat false, indicates that the Gaussian
  !                        distribution is multiplied by the ratio of the local true
  !                        area to the coordinate area (and then normalised to
  !                        integrate to 1).
  ! NoReflect           :: Indicates no reflection of the source at the ground, with
  !                        any part of the source below the ground being lost.
  ! StartTime           :: Start time of the source.
  ! StopTime            :: End time of the source.
  ! MaxTravelTime       :: Maximum travel time for which particles are to be followed.
  ! CharNParticles      :: Character string indicating number of particles released. Character format for
  !                        number of particles to be released is
  !                                        #Particles{/[s|sec|min|hr|day]}
  !                        Here {...} indicates optional components, and [ | ] indicates alternative options.
  !                        Spaces are not allowed within #Particles or within the words 'sec', 'min', 'hr' or
  !                        'day'. Spaces are permitted either side of / but are not required. 's', 'sec',
  !                        'min', 'hr' and 'day' are case insensitive.
  ! CharSourceStrengths :: Character string indicating the source strength for each species released.
  !                        Character format for source strength is
  !                                SpeciesName SourceStrength MaterialUnit{/[s|sec|min|hr|day]}
  !                        Here {...} indicates optional components, and [ | ] indicates alternative options.
  !                        Spaces are not generally allowed within SpeciesName, SourceStrength or MaterialUnit
  !                        or within the words 'sec', 'min', 'hr' or 'day'; however they are permitted in
  !                        SpeciesName if it is enclosed in quotes. Quotes can be '...' or "...", but quotes
  !                        of the same type mustn't be used within SpeciesName. Spaces between SpeciesName,
  !                        SourceStrength and MaterialUnit are compulsory. Spaces are permitted either side of
  !                        / but are not required. SpeciesName, MaterialUnit, 's', 'sec', 'min', 'hr' and
  !                        'day' are case insensitive.
  ! PlumeRise               :: Indicates plume rise is to be modelled.
  ! VolumeFlowRate          :: Volume flow rate of the source. } Defined only if PlumeRise
  ! FlowVelocity            :: Flow velocity.                  } is true. The plume rise
  ! Temperature             :: Temperature of emission (K).    } parameters  refer to the
  !                                                              total emission, not just
  !                                                              the species of interest. 
  !                                                              For iterative plume-rise 
  !                                                              model temperature and flow 
  !                                                              (exit) velocity used to  
  !                                                              initialise model. 
  ! Particulate          :: Indicates that the species are emitted as particulates.
  ! Diameter             :: Particle diameter. Defined only if                } Defined only if
  !                         SizeDistName is blank.                            } Particulate is 
  ! Density              :: Particle density.                                 } true.  
  ! ParticleShape        :: Particle shape.                                   }
  ! DensityPresent       :: Indicates whether a particle shape has been given }
  ! ParticleShapePresent :: Indicates whether a particle shape has been given }
  ! ShapeScheme          :: Identifies which scheme to calculate CD           }
  ! SizeDistName         :: Name of particle size distribution. Blank         }
  !                         indicates no particle size distribution.          }
  ! TimeDepName          :: Name of source time dependency. Blank indicates no source
  !                         time dependency.
  ! FixedMet             :: Indicates that the met is fixed.
  ! EulerianTime         :: Time after which particles' mass is transferred to Eulerian field.
  ! ZGridIterativePlumeRise :: Name of Z-grid for iterative plume-rise model 
  ! DryGasFraction0         :: Mass fraction of dry gas for iterative plume-rise model 
  ! WaterVapourFraction0    :: Mass fraction of water vapour for iterative plume-rise model 
  ! PlumeRiseHeight         :: Maximum rise height of plume
  ! SummitElevation         :: Height of source
  ! DistalFineAshFraction   :: Fraction of ash that survives into the far field
  ! Function result:
  Type (Source_) :: Source ! Initialised source.
  ! Locals:
  Real(Std) :: SourceDuration             ! Duration of source.
  Logical   :: QuestionMark(MaxSpecieses) ! Indicates a question mark has been used for the source strength in
                                          ! CharSourceStrengths.
  Integer   :: i                          ! Loop index.
  Character(MaxCharLength) :: MaterialUnitName

  ! $$ check species and source group names not duplicated.

  ! Name.
  Call TokenLengthTest(             &
         C         = Name,          &
         Length    = MaxCharLength, &
         Zero      = .true.,        &
         BlockKey  = 'Sources',     &
         Item      = ' ',           &
         ColumnKey = 'Name'         &
       )
  Source%Name = Name

  ! Source groups.
  If (Size(SourceGroupNames) > MaxSourceGroupsPerSource) Then
    Call Message(                                    &
           'FATAL ERROR: Source '                 // &
           Trim(Source%Name)                      // &
           ' belongs to too many source groups.',    &
           3                                         &
         )
  End If
  Do i = 1, Size(SourceGroupNames)
    Call TokenLengthTest(                   &
           C         = SourceGroupNames(i), &
           Length    = MaxCharLength,       &
           Zero      = .true.,              &
           BlockKey  = 'Sources',           &
           Item      = Source%Name,         &
           ColumnKey = 'Source Groups',     &
           List      = .true.               &
         )
  End Do
  Source%nSourceGroups                            = Size(SourceGroupNames)
  Source%SourceGroupNames(1:Source%nSourceGroups) = SourceGroupNames

  ! Source type.
  If (SourceType < 0 .or. SourceType > Size(SourceTypeNames)) Then
    Call Message('FATAL ERROR: Source type wrong', 3) ! $$ better message.
                                                      ! $$ input to init source as char and store char and int
  End If
  Source%SourceType = SourceType

  ! Grids and coordinate systems.
  Call TokenLengthTest(             &
         C         = HGridName,     &
         Length    = MaxCharLength, &
         Zero      = .false.,       &
         BlockKey  = 'Sources',     &
         Item      = Name,          &
         ColumnKey = 'H-Grid'       &
       )
  Call TokenLengthTest(             &
         C         = ZGridName,     &
         Length    = MaxCharLength, &
         Zero      = .false.,       &
         BlockKey  = 'Sources',     &
         Item      = Name,          &
         ColumnKey = 'Z-Grid'       &
       )
  Call TokenLengthTest(             &
         C         = HCoordName,    &
         Length    = MaxCharLength, &
         Zero      = .false.,       &
         BlockKey  = 'Sources',     &
         Item      = Name,          &
         ColumnKey = 'H-Coord'      &
       )
  Call TokenLengthTest(             &
         C         = dHCoordName,   &
         Length    = MaxCharLength, &
         Zero      = .false.,       &
         BlockKey  = 'Sources',     &
         Item      = Name,          &
         ColumnKey = 'dH-Coord'     &
       )
  Call TokenLengthTest(             &
         C         = ZCoordName,    &
         Length    = MaxCharLength, &
         Zero      = .false.,       &
         BlockKey  = 'Sources',     &
         Item      = Name,          &
         ColumnKey = 'Z-Coord'      &
       )
  Call TokenLengthTest(             &
         C         = dZCoordName,   &
         Length    = MaxCharLength, &
         Zero      = .false.,       &
         BlockKey  = 'Sources',     &
         Item      = Name,          &
         ColumnKey = 'dZ-Coord'     &
       )
  Call TokenLengthTest(                         &
         C         = ZGridIterativePlumeRise,   &
         Length    = MaxCharLength,             &
         Zero      = .false.,                   &
         BlockKey  = 'Sources',                 &
         Item      = Name,                      &
         ColumnKey = 'Z-GridIterativePlumeRise' &
       )

  If (HGridName == ' ') Then

    If (HCoordName == ' ') Then
      Call Message('FATAL ERROR: H-Coord is blank', 3)
    End If
    If (dHCoordName == ' ') Then
      Call Message('FATAL ERROR: dH-Coord is blank', 3)
    End If
    Source%HGridName   = HGridName
    Source%HCoordName  = HCoordName
    Source%dHCoordName = dHCoordName
    Source%dHMetres    = dHMetres

  Else

    If (HCoordName /= ' ') Then
      Call Message('FATAL ERROR: H-Coord should be blank', 3)
    End If
    If (dHCoordName /= ' ') Then
      Call Message('FATAL ERROR: dH-Coord should be blank', 3)
    End If
    If (dHMetres) Then
      Call Message('FATAL ERROR: dH-Metres should be false', 3)
    End If
    Source%HGridName = HGridName
    Source%dHMetres  = dHMetres

  End If

  If (ZGridName == ' ') Then

    If (ZCoordName == ' ') Then
      Call Message('FATAL ERROR: Z-Coord is blank', 3)
    End If
    If (dZCoordName == ' ') Then
      Call Message('FATAL ERROR: dZ-Coord is blank', 3)
    End If
    Source%ZGridName   = ZGridName
    Source%ZCoordName  = ZCoordName
    Source%dZCoordName = dZCoordName
    Source%dZMetres    = dZMetres

  Else

    If (ZCoordName /= ' ') Then
      Call Message('FATAL ERROR: Z-Coord should be blank', 3)
    End If
    If (dZCoordName /= ' ') Then
      Call Message('FATAL ERROR: dZ-Coord should be blank', 3)
    End If
    If (dZMetres) Then
      Call Message('FATAL ERROR: dZ-Metres should be false', 3)
    End If
    Source%ZGridName = ZGridName
    Source%dZMetres  = dZMetres

  End If

  ! Source spatial distribution.
  If (.not.(Shape .CIEq. 'Ellipsoid' ) .and. &
      .not.(Shape .CIEq. 'Cuboid'    ) .and. &
      .not.(Shape .CIEq. 'Cylindroid')       &
     ) Then
    Call Message('FATAL ERROR: Source shape is not "Ellipsoid", "Cylindroid" or "Cuboid"', 3)
  End If

  If (HGridName /= ' ' .or. ZGridName /= ' ') Then
    If (.not. (Shape .CIEq. 'Cuboid')) Then
      Call Message('FATAL ERROR: Source shape is not "Cuboid"', 3)
    End If
    If (.not. TopHat) Then
      Call Message('FATAL ERROR: Source must have Top Hat = Yes', 3)
    End If
  End If

  If (HGridName /= ' ') Then
    Source%X(1:2)  = 0.0
    Source%dX(1:2) = 0.0
    Source%Angle   = 0.0
  Else
    Source%X(1:2)  = X(1:2)
    Source%dX(1:2) = dX(1:2)
    Source%Angle = Angle * Pi / 180.0
  End If
  If (ZGridName /= ' ') Then
    Source%X(3)  = 0.0
    Source%dX(3) = 0.0
  Else
    Source%X(3)  = X(3)
    Source%dX(3) = dX(3)
  End If
  Source%Shape       = Shape
  Source%TopHat      = TopHat
  Source%UniformArea = UniformArea
  Source%NoReflect   = NoReflect

  ! Source times.
  If (IsTimeInterval(StartTime) .or. IsTimeInterval(StopTime) .or. .not.IsTimeInterval(MaxTravelTime)) Then
    Call Message('UNEXPECTED FATAL ERROR in InitSource', 4)
  End If
  If (MaxTravelTime < ZeroTime()) Then
    Call Message(                                            &
           'FATAL ERROR: The max travel time of source "' // &
           Trim(Name)                                     // &
           '" is < 0',                                       &
           3                                                 &
         )
  End If
  Source%StartTime      = StartTime
  Source%StopTime       = StopTime
  Source%MaxTravelTime  = MaxTravelTime
  Source%sStartTime     = Time2ShortTime(StartTime)
  Source%sStopTime      = Time2ShortTime(StopTime)
  Source%sMaxTravelTime = Time2ShortTime(MaxTravelTime)

  ! # particles.
  Call ParseNParticles(CharNParticles, Source%nParticles, Source%ParticleRate)
  If (Source%nParticles < 0) Then
    Call Message(                                     &
           'FATAL ERROR: # Particles for source "' // &
           Trim(Name)                              // &
           '" is < 0',                                &
           3                                          &
         )
  End If

  ! Source strengths and species.
  If (Size(CharSourceStrengths) == 0) Then
    Call Message(                                                           &
           'FATAL ERROR: No source strengths are specified for source "' // &
           Trim(Name)                                                    // &
           '"',                                                             &
           3                                                                &
         )
  End If
  If (Size(CharSourceStrengths) > MaxSpecieses) Then
    Call Message(                                                             &
           'FATAL ERROR: Too many source strengths specified for source "' // &
           Trim(Name)                                                      // &
           '"',                                                               &
           3                                                                  &
         )
  End If

  Source%nSpecieses = Size(CharSourceStrengths)

  Do i = 1, Size(CharSourceStrengths)
    Call ParseSourceStrength(          &
           CharSourceStrengths(i),     &
           Source%SpeciesNames(i),     &
           Source%SourceStrengths(i),  &
           Source%MaterialUnitName(i), &
           Source%SourceRate(i),       &
           QuestionMark(i)             &
         )
    If (Source%SourceStrengths(i) < 0) Then
      Call Message(                                                          &
             'FATAL ERROR: '                                              // &
             'A negative source strength has been specified for source "' // &
             Trim(Name)                                                   // &
             '"',                                                            &
             3                                                               &
           )
    End If
  End Do

  Source%iBioIsoprene = 0
  Source%iBioAPinene  = 0
  Source%MetDependent = 0

  Do i = 1, Size(CharSourceStrengths)
    If (Source%SpeciesNames(i) .CIEq. 'C10H16bio') Then
      Source%MetDependent = 1
      Source%iBioAPinene  = i
    Else If (Source%SpeciesNames(i) .CIEq. 'C5H8bio') Then
! $$ set c5h8bio to zero for now - emissions are in the source file but
! a light dependent release for isoprene is needed before c5h8bio can be used
       Source%SourceStrengths(i) = 0.0
!      Source%MetDependent       = 1
!      Source%iBioIsoprene       = i
    Else If (Source%SpeciesNames(i) .CIEq. 'RESUSPENDED_ASH') Then
       Source%MetDependent       = 1  
    Else If (Source%SpeciesNames(i) .CIEq. 'MIDGE') Then
       Source%MetDependent       = 1     
    End If
  End Do

  ! For instantaneous sources, convert to not use rates (except for particle rate with
  ! fixed met, which is treated separately below).
  If (Source%sStartTime == Source%sStopTime) Then

    If (Source%ParticleRate .and. .not. FixedMet) Then
      Source%nParticles   = 0
      Source%ParticleRate = .false.
    End If
    Do i = 1, Source%nSpecieses
      If (Source%SourceRate(i)) Then
        Source%SourceStrengths(i) = 0
        Source%SourceRate     (i) = .false.
      End If
    End Do

  ! For infinite-duration sources, convert to rates (except for particle rate with
  ! fixed met, which is treated separately below).
  Else If (IsInfFuture(Source%sStopTime - Source%sStartTime)) Then

    If (.not. Source%ParticleRate .and. .not. FixedMet) Then
      Source%nParticles   = 0
      Source%ParticleRate = .true.
    End If
    Do i = 1, Source%nSpecieses
      If (.not. Source%SourceRate(i)) Then
        Source%SourceStrengths(i) = 0
        Source%SourceRate     (i) = .true.
      End If
    End Do

  ! For non-instantaneous non-infinite-duration sources, convert to rates (except for
  ! particle rate with fixed met, which is treated separately below). If StopTime <
  ! StartTime we do nothing as source won't release.
  Else If (Source%sStopTime > Source%sStartTime) Then

    If (.not. Source%ParticleRate .and. .not. FixedMet) Then
      SourceDuration      = ShortTime2RealTime(Source%sStopTime - Source%sStartTime)
      Source%nParticles   = Source%nParticles / SourceDuration
      Source%ParticleRate = .true.
    End If
    Do i = 1, Source%nSpecieses
      If (.not. Source%SourceRate(i)) Then
        SourceDuration            = ShortTime2RealTime(Source%sStopTime - Source%sStartTime)
        Source%SourceStrengths(i) = Source%SourceStrengths(i) / SourceDuration
        Source%SourceRate     (i) = .true.
      End If
    End Do

  End If

  ! For fixed met, must not give nParticles as a rate. (For fixed met, particles and
  ! puffs for a source are all released at the same time, which can be chosen
  ! arbitrarily and is representative of all release times).
  If (Source%ParticleRate .and. FixedMet) Then
    Call Message(                                         &
           'FATAL ERROR: For fixed met "# Particles" ' // &
           'must not be given as a rate.',                &
           3                                              &
         )
  End If

  ! Plume rise.
  Source%PlumeRise = PlumeRise
  If (PlumeRise) Then                          !! $$ may want to change this 
    Source%Temperature    = Temperature
    Source%VolumeFlowRate = VolumeFlowRate
    Source%FlowVelocity   = FlowVelocity
  End If

  ! Iterative plume rise model
  If (Source%SourceType == S_IterativePlumeModel) Then
    Source%MetDependent = 2 !$$ Would be better to use SourceType to flag special met dependency.
    If ( Source%nSpecieses > 1 ) Then
      Call Message(                                                                &
             'FATAL ERROR: iterative plume-rise model only works for one species', &
             3                                                                     &
           )
    End If
    If ( Source%StartTime == Source%StopTime ) Then
      Call Message(                                                                             &
             'FATAL ERROR: iterative plume-rise model does not work for instantaneous sources', &
             3                                                                                  &
           )
    End If
    Source%WaterVapourFraction0    = WaterVapourFraction0
    Source%GasFraction0            = DryGasFraction0 + WaterVapourFraction0
    Source%Temperature             = Temperature
    Source%FlowVelocity            = FlowVelocity
    Source%PlumeRiseHeight         = PlumeRiseHeight
    Source%SummitElevation         = SummitElevation
    Source%DistalFineAshFraction   = DistalFineAshFraction
    Source%ZGridIterativePlumeRise = ZGridIterativePlumeRise
  Else
    Source%ZGridIterativePlumeRise = ' '
    Source%WaterVapourFraction0  = 0.0
    Source%GasFraction0          = 0.0
    Source%PlumeRiseHeight       = 0.0
    Source%SummitElevation       = 0.0
    Source%DistalFineAshFraction = 0.0
  End If

  ! Particulates.
  Source%Particulate = Particulate
  If (Particulate) Then
    If (SizeDistName == ' ') Then
      If (Diameter <= 0.0) Then
        Call Message(                                    &
               'FATAL ERROR: Particle diameter is <= 0', &
               3                                         &
             )
      End If
    Else
      Call TokenLengthTest(                           &
             C         = SizeDistName,                &
             Length    = MaxCharLength,               &
             Zero      = .true.,                      &
             BlockKey  = 'Sources',                   &
             Item      = Name,                        &
             ColumnKey = 'Particle Size Distribution' &
           )
   
    End If

    Source%Diameter             = Diameter
    Source%Density              = Density
    Source%DensityPresent       = DensityPresent 
    Source%ParticleShape        = ParticleShape
    Source%ParticleShapePresent = ParticleShapePresent 
    Source%ShapeScheme          = ShapeScheme 
    Source%SizeDistName         = SizeDistName

    If (ShapeScheme .CIEq. 'Ganser') Then
      Source%ShapeSchemeCode = ShapeScheme_Ganser
    Else If (ShapeScheme .CIEq. 'WH') Then
      Source%ShapeSchemeCode = ShapeScheme_WH
    Else If (ShapeScheme .CIEq. 'White') Then
      Source%ShapeSchemeCode = ShapeScheme_White
    Else If (ShapeScheme .CIEq. '  ') Then
      Source%ShapeSchemeCode = ShapeScheme_White  
    End If

  End If

  ! Source time dependency.
  Call TokenLengthTest(                &
         C         = TimeDepName,      &
         Length    = MaxCharLength,    &
         Zero      = .false.,          &
         BlockKey  = 'Sources',        &
         Item      = Name,             &
         ColumnKey = 'Time Dependency' &
       )
  Source%TimeDepName = TimeDepName

  ! Fixed met.
  Source%FixedMet = FixedMet

  ! Lagrangian-Eulerian time
  Source%EulerianTime = EulerianTime

  ! Restrictions for specific source types.

  ! Generic sources.
  If (SourceType == S_Generic .Or. SourceType == S_IterativePlumeModel) Then

    ! Coord system and grid restrictions.
    If (HGridName /= ' ') Then
      Call Message('FATAL ERROR: HGridName is not blank for a non met dependent source', 3)
    End If
    If (ZGridName /= ' ') Then
      Call Message('FATAL ERROR: ZGridName is not blank for a non met dependent source', 3)
    End If

    ! Source strength restrictions.
    Do i = 1, Source%nSpecieses
      If (QuestionMark(i)) Then
        Call Message('FATAL ERROR: The source strength is "?" for a non met dependent source', 3)
      End If
    End Do

  ! Dust and sea-salt sources.
  Else If (SourceType == S_Dust .or. SourceType == S_SeaSalt) Then

    ! Coord system and grid restrictions.
    If (HGridName == ' ') Then
      Call Message('FATAL ERROR: H-Grid is blank for a dust or sea-salt source', 3)
    End If
    If (ZGridName /= ' ') Then
      Call Message('FATAL ERROR: Z-Grid is not blank for a dust or sea-salt source', 3)
    End If
    If (.not.(ZCoordName .CIEq. 'm agl')) Then
      Call Message('FATAL ERROR: Z-Coord is not "m agl" for a dust or sea-salt source', 3)
    End If

    ! Source spatial distribution restrictions.
    If (X(3) /= 0.0) Then
      Call Message('FATAL ERROR: Z is not 0 for a dust or sea-salt source', 3)
    End If
    If (dX(3) /= 0.0) Then
      Call Message('FATAL ERROR: dZ is not 0 for a dust or sea-salt source', 3)
    End If

    ! Source strength restrictions.
    If (Source%nSpecieses > 1) Then
      Call Message('FATAL ERROR: Number of species is > 1 for a dust or sea-salt source', 3)
    End If
    Do i = 1, Source%nSpecieses
      If (Trim(Source%MaterialUnitName(i))  /= 'g') Then
        Call Message('FATAL ERROR: Material unit is not "g" for a dust or sea-salt source', 3)
      End If
      If (.not. QuestionMark(i)) Then
        Call Message('FATAL ERROR: The source strength is not "?" for a dust or sea-salt source', 3)
      End If
    End Do

    ! Plume rise restrictions.
    If (PlumeRise) Then
      Call Message('FATAL ERROR: plume rise is not possible for a dust or sea-salt source', 3)
    End If

    ! Particulates.
    If (.not. Particulate) Then
      Call Message('FATAL ERROR: only particulates are possible for a dust or sea-salt source', 3)
    End If
    If (SizeDistName == ' ') Then
      Call Message('FATAL ERROR: a particle size distribution must be given for a dust or sea-salt source', 3)
    End If

    ! Source time dependency.
    If (TimeDepName /= ' ') Then
      Call Message('FATAL ERROR: a time dependency must not be given for a dust or sea-salt source', 3)
    End If

  End If

End Function InitSource

!-------------------------------------------------------------------------------------------------------------

Subroutine AddSource(Source, Sources)
! Adds a source to a collection of sources.

  Implicit None
  ! Argument list:
  Type(Source_),  Intent(In)    :: Source  ! Source to be added.
  Type(Sources_), Intent(InOut) :: Sources ! Collection of sources.
  ! Locals:
  Integer :: i ! Loop index.

  ! Note the test for identical source names is done later on in SetupSources.

  If (Sources%nSources >= Size(Sources%Sources)) Then
    Call Message(                                                                                        &
           'FATAL ERROR in adding the source "'                                                       // &
           Trim(Source%Name)                                                                          // &
           '": there are too many sources. The maximum number of sources is set to '                  // &
           Trim(Int2Char(Size(Sources%Sources)))                                                      // &
           ' and can be altered via the variable "Max # Sources" in the "Main Options" input block.',    &
           3                                                                                             &
         )
  End If

  Sources%nSources                  = Sources%nSources + 1
  Sources%Sources(Sources%nSources) = Source

  ! $$ Check consistency of fixed met between sources (as in Flows).

End Subroutine AddSource

!-------------------------------------------------------------------------------------------------------------

Function FindSourceIndex(Name, Sources)
! Finds the index of a source in a collection of sources.

  Implicit None
  ! Argument list:
  Character(*),   Intent(In) :: Name    ! Name of the source.
  Type(Sources_), Intent(In) :: Sources ! Collection of sources.
  ! Function result:
  Integer :: FindSourceIndex ! Index of the source.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Sources%nSources
    If (Name .CIEq. Sources%Sources(i)%Name) Then
      FindSourceIndex = i
      Return
    End If
  End Do

  Call Message(                     &
         'FATAL ERROR: source "' // &
         Trim(Name)              // &
         '" not found',             &
         3                          &
       )

End Function FindSourceIndex

!-------------------------------------------------------------------------------------------------------------

Function FindSourceGroupIndex(Name, Sources)
! Finds the index of a source group in a collection of sources.

  Implicit None
  ! Argument list:
  Character(*),   Intent(In) :: Name    ! Name of the source group.
  Type(Sources_), Intent(In) :: Sources ! Collection of sources.
  ! Function result:
  Integer :: FindSourceGroupIndex ! Index of the source group.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Sources%nSourceGroups
    If (Name .CIEq. Sources%SourceGroups(i)%Name) Then
      FindSourceGroupIndex = i
      Return
    End If
  End Do

  Call Message(                           &
         'FATAL ERROR: source group "' // &
         Trim(Name)                    // &
         '" not found',                   &
         3                                &
       )

End Function FindSourceGroupIndex

!-------------------------------------------------------------------------------------------------------------

Function InitTimeDeps() Result (TimeDeps)
! Initialises a collection of source time dependencies.

  Implicit None
  ! Function result:
  Type(TimeDeps_) :: TimeDeps ! Initialised collection of source time dependencies.

  TimeDeps%nTimeDeps = 0

End Function InitTimeDeps

!-------------------------------------------------------------------------------------------------------------

Function InitTimeDep(Name, Factors, nWildTimePairs, WildTimePairs) Result (TimeDep)
! Initialises a source time dependency.

  Implicit None
  ! Argument list:
  Character(*),        Intent(In) :: Name
  Real(Std),           Intent(In) :: Factors(:)
  Integer,             Intent(In) :: nWildTimePairs(:)
  Type(WildTimePair_), Intent(In) :: WildTimePairs(:, :)
  ! Name           :: Name of source time dependency.
  ! Factors        :: The source time dependency factors.
  ! nWildTimePairs :: Number of wild-card time pairs used with each source time
  !                   dependency factor.
  ! WildTimePairs  :: Wild-card time pairs. The (j, i)th element is the jth wild-card
  !                   time pair for the ith factor. The second dimension must match
  !                   Factors in size. The first dimension must not be less than, but
  !                   can be larger than, the number of time pairs for each factor.
  ! Function result:
  Type(TimeDep_) :: TimeDep ! Initialised source time dependency.
  ! Locals:
  Integer :: i !} Loop indices.
  Integer :: j !}

  If (Name == ' ') Then
    Call Message(                                                                 &
           'FATAL ERROR in InitTimeDep: name of source time dependency is blank', &
           3                                                                      &
         )
  End If
  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message(                                            &
           'FATAL ERROR in InitTimeDep: '                 // &
           'name of source time dependency is given as "' // &
           Trim(Name)                                     // &
           '" and is too long',                              &
           3                                                 &
         )
  End If
  TimeDep%Name = Name

  If (Size(Factors) > MaxTimeDepFactors) Then
    Call Message(                                            &
           'FATAL ERROR in InitTimeDep: '                 // &
           'too many factors in source time dependency "' // &
           Trim(Name)                                     // &
           '"',                                              &
           3                                                 &
         )
  End If
  TimeDep%nFactors = Size(Factors)

  TimeDep%Factors(1:TimeDep%nFactors) = Factors(1:TimeDep%nFactors)

  If (Size(WildTimePairs, 2) /= TimeDep%nFactors) Then
    Call Message('UNEXPECTED FATAL ERROR in InitTimeDep', 4)
  End If

  Do i = 1, TimeDep%nFactors
    If (Size(WildTimePairs, 1) < nWildTimePairs(i)) Then
      Call Message('UNEXPECTED FATAL ERROR in InitTimeDep', 4)
    End If
    TimeDep%nWildTimePairs(i) = nWildTimePairs(i)
    Do j = 1, nWildTimePairs(i)
      TimeDep%WildTimePairs(j, i) = WildTimePairs(j, i)
    End Do
  End Do

End Function InitTimeDep

!-------------------------------------------------------------------------------------------------------------

Subroutine AddTimeDep(TimeDep, TimeDeps)
! Adds a source time dependency to a collection of source time dependencies.

  Implicit None
  ! Argument list:
  Type(TimeDep_),  Intent(In)    :: TimeDep  ! Source time dependency to be added.
  Type(TimeDeps_), Intent(InOut) :: TimeDeps ! Collection of source time dependencies.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, TimeDeps%nTimeDeps
    If (TimeDep%Name .CIEq. TimeDeps%TimeDeps(i)%Name) Then
      If (TimeDep == TimeDeps%TimeDeps(i)) Then
        Return
      Else
        Call Message(                                                         &
               'FATAL ERROR in adding the source time dependency "'        // &
               Trim(TimeDep%Name)                                          // &
               '": a different source time dependency with the same name ' // &
               'already exists',                                              &
               3                                                              &
             )
      End If
    End If
  End Do

  If (TimeDeps%nTimeDeps >= MaxTimeDeps) Then
    Call Message(                                                  &
           'FATAL ERROR in adding the source time dependency "' // &
           Trim(TimeDep%Name)                                   // &
           '": there are too many source time dependencies',       &
           3                                                       &
         )
  End If

  TimeDeps%nTimeDeps                    = TimeDeps%nTimeDeps + 1
  TimeDeps%TimeDeps(TimeDeps%nTimeDeps) = TimeDep

End Subroutine AddTimeDep

!-------------------------------------------------------------------------------------------------------------

Function FindTimeDepIndex(Name, TimeDeps)
! Finds the index of a source time dependency in a collection of source time
! dependencies.

  Implicit None
  ! Argument list:
  Character(*),    Intent(In) :: Name     ! Name of the source time dependency.
  Type(TimeDeps_), Intent(In) :: TimeDeps ! Collection of source time dependencies.
  ! Function result:
  Integer :: FindTimeDepIndex ! Index of the source time dependency.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, TimeDeps%nTimeDeps
    If (Name .CIEq. TimeDeps%TimeDeps(i)%Name) Then
      FindTimeDepIndex = i
      Return
    End If
  End Do

  Call Message(                                     &
         'FATAL ERROR: source time dependency "' // &
         Trim(Name)                              // &
         '" not found',                             &
         3                                          &
       )

End Function FindTimeDepIndex

!-------------------------------------------------------------------------------------------------------------

Function TimeDepEq(TimeDep1, TimeDep2)
! Tests for equality of source time dependencies.

  Implicit None
  ! Argument list:
  Type(TimeDep_), Intent(In) :: TimeDep1 !} The source time dependencies.
  Type(TimeDep_), Intent(In) :: TimeDep2 !}
  ! Function result:
  Logical :: TimeDepEq ! Indicates if source time dependencies are equal.
  ! Locals:
  Integer :: i !} Loop indices.
  Integer :: j !}

  TimeDepEq =                                         &
    (TimeDep1%Name     .CIEq. TimeDep2%Name)    .and. &
     TimeDep1%nFactors   ==   TimeDep2%nFactors
  Do i = 1, Min(TimeDep1%nFactors, TimeDep2%nFactors)
    TimeDepEq =                                                        &
      TimeDepEq                                                  .and. &
      (TimeDep1%nWildTimePairs(i) == TimeDep2%nWildTimePairs(i)) .and. &
      (TimeDep1%Factors       (i) == TimeDep2%Factors       (i))
    Do j = 1, Min(TimeDep1%nWildTimePairs(i), TimeDep2%nWildTimePairs(i))
      TimeDepEq =                                                            &
        TimeDepEq                                                      .and. &
        (TimeDep1%WildTimePairs(j, i) == TimeDep2%WildTimePairs(j, i))
    End Do
  End Do

End Function TimeDepEq

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpSources(Grids, SizeDists, Specieses, MaterialUnits, Sources)
! Sets up Sources.

  Implicit None
  ! Argument list:
  Type(Grids_),         Intent(In),    Target :: Grids         ! Collection of grids.
  Type(SizeDists_),     Intent(In),    Target :: SizeDists     ! Collection of size distributions.
  Type(Specieses_),     Intent(In),    Target :: Specieses     ! Collection of species.
  Type(MaterialUnits_), Intent(In)            :: MaterialUnits ! Collection of material units.
  Type(Sources_),       Intent(InOut), Target :: Sources       ! Collection of sources.
  ! Locals:
  Type(HGrid_),            Pointer :: HGrid
  Type(ZGrid_),            Pointer :: ZGrid
  Type(SizeDist_),         Pointer :: SizeDist 
  Type(Species_),          Pointer :: Species
  Type(Source_),           Pointer :: Source
  Real(Std)                        :: ParticleRate
  Integer                          :: iHGrid
  Integer                          :: iZGrid
  Integer                          :: iSpecies
  Integer                          :: i
  Integer                          :: j
  Integer                          :: k
  Character(MaxCharLength)         :: Names(Sources%nSources*MaxSourceGroupsPerSource)
  Type(MaterialUnit_)              :: SourceMaterialUnit
  Type(MaterialUnit_)              :: SpeciesMaterialUnit
  ! HGrid        :} Abbreviations for grids, size distributions, species and sources.
  ! ZGrid        :}
  ! SizeDist     :}
  ! Species      :}
  ! Source       :}
  ! ParticleRate :: Particle release rate (for continuous sources) or source mass as a fraction of particle
  !                 mass (for instantaneous sources) if release rate limited by the MassLimit of a particular
  !                 species.
  ! iHGrid       :} Indices for grids and species.
  ! iZGrid       :}
  ! iSpecies     :}
  ! i            :] Loop indices.
  ! j            :]
  ! k            :]
  ! Names        :: Copy of source or source group names.

  ! Check for identical source names.
  Do i = 1, Sources%nSources
    Names(i) = Sources%Sources(i)%Name
  End Do
  Call SortInPlace(OrderingType = 'CI', C = Names(1:Sources%nSources))
  Do i = 1, Sources%nSources - 1
    If (Names(i) .CIEq. Names(i + 1)) Then
      Do j = 1, Sources%nSources - 1
        If (Names(i) .CIEq. Sources%Sources(j)%Name) Then
          Call Message(                                                                &
                 'FATAL ERROR: Two sources have the same (case insensitive) name "' // &
                 Trim(Sources%Sources(j)%Name)                                      // &
                 '".',                                                                 &
                 3                                                                     &
               )
        End If
      End Do
      Call Message('UNEXPECTED FATAL ERROR in SetUpSources.', 4)
    End If
  End Do

  ! Assemble list of source groups.
  ! $$ this could sometimes, but not always, be more efficient if sorting used.
  ! The If false block below does this.
  If (.true.) Then
    iLoop: Do i = 1, Sources%nSources
      Source => Sources%Sources(i)
      jLoop: Do j = 1, Source%nSourceGroups
        kLoop: Do k = 1, Sources%nSourceGroups
          If (Source%SourceGroupNames(j) .CIEq. Sources%SourceGroups(k)%Name) Cycle jLoop
        End Do kLoop
        Sources%nSourceGroups                            = Sources%nSourceGroups + 1
        Sources%SourceGroups(Sources%nSourceGroups)%Name = Source%SourceGroupNames(j)
      End Do jLoop
    End Do iLoop
  Else
    ! Assemble list of source group names.
    k = 0
    Do i = 1, Sources%nSources
      Source => Sources%Sources(i)
      Do j = 1, Source%nSourceGroups
        k = k + 1
        Names(k) = Source%SourceGroupNames(j)
      End Do
    End Do
    ! Sort list.
    Call SortInPlace(OrderingType = 'CI', C = Names(1:k))
    ! Store source groups in Sources, eliminating duplicates.
    Do i = 1, k
      If (i == 1) Then
        Sources%nSourceGroups                            = Sources%nSourceGroups + 1
        Sources%SourceGroups(Sources%nSourceGroups)%Name = Names(i)
      Else If (.not. (Names(i) .CIEq. Names(i - 1))) Then
        Sources%nSourceGroups                            = Sources%nSourceGroups + 1
        Sources%SourceGroups(Sources%nSourceGroups)%Name = Names(i)
      End If
    End Do
  End If

  Do i = 1, Sources%nSources

    Source => Sources%Sources(i)

    ! Gridded sources - get coord system names.
    If (Source%HGridName /= ' ') Then

      iHGrid             =  FindHGridIndex(Source%HGridName, Grids)
      HGrid              => Grids%HGrids(iHGrid)
      Source%HCoordName  =  HGrid%HCoordName
      Source%dHCoordName =  HGrid%HCoordName

    End If

    If (Source%ZGridName /= ' ') Then

      iZGrid             =  FindZGridIndex(Source%ZGridName, Grids)
      ZGrid              => Grids%ZGrids(iZGrid)
      Source%ZCoordName  =  ZGrid%ZCoordName
      Source%dZCoordName =  ZGrid%ZCoordName

    End If

    Source%iCritSpecies = 0
    ParticleRate        = 0.0

    Do j = 1, Source%nSpecieses

      iSpecies = FindSpeciesIndex(Source%SpeciesNames(j), Specieses)
      Species => Specieses%Specieses(iSpecies)

      ! $$ In future allow conversions.
      SourceMaterialUnit = MaterialUnits%MaterialUnits(                                      &
                             FindMaterialUnitIndex(Source%MaterialUnitName(j),MaterialUnits) &
                           )
      SpeciesMaterialUnit = MaterialUnits%MaterialUnits(                                     &
                              FindMaterialUnitIndex(Species%MaterialUnitName,MaterialUnits)  &
                            )
      If (.not.(SourceMaterialUnit == SpeciesMaterialUnit)) Then
        Call Message(                                                              &
               'FATAL ERROR: The material unit given for source strength must ' // &
               'match that given for the species',                                 &
               3                                                                   &
             )
      End If

      ! Calculate which species if any will give lower limit on particle release rate
      ! due to Species%MassLimit. Note Species%MassLimit only applies for non-fixed
      ! met.
      If (.not.Source%FixedMet) Then
        If (Species%UseMassLimit) Then
          If (Source%SourceStrengths(j) / Species%MassLimit > ParticleRate) Then
            Source%iCritSpecies = j
            ParticleRate = Source%SourceStrengths(j) / Species%MassLimit
          End If
        End If
      End If

    End Do

    ! Check density available for particulates, either via the source or the particle size distribution.
    If (Source%Particulate)  Then
      If (Source%SizeDistName /= ' ') Then
        SizeDist => SizeDists%SizeDists(FindSizeDistIndex(Source%SizeDistName, SizeDists))
        If (.not. (Source%DensityPresent .or. SizeDist%DensityPresent)) Then ! $$ here an below could check both not present?
                                                                             ! If allow both present, which do we want to use?
                                                                             ! At present the PSD value is used if both present.
          Call Message(                                                                       &
                 'FATAL ERROR: A density needs to be given for particulates, '             // &
                 'either via the source or via an associated particle size distribution.',    &
                 3                                                                            &
               )  
        End If
      Else
        If (.not. Source%DensityPresent) Then
          Call Message(                                                                       &
                 'FATAL ERROR: A density needs to be given for particulates, '             // &
                 'either via the source or via an associated particle size distribution.',    &
                 3                                                                            &
               )  
        End If
      End If
    End If

    ! Check shape scheme specified if particle shape available.
    ! Note Source%ShapeScheme, unlike Source%ShapeSchemeCode, distinguishes between ' ' and 'White'.
    If (Source%Particulate)  Then
      If (Source%SizeDistName /= ' ')  Then
        SizeDist => SizeDists%SizeDists(FindSizeDistIndex(Source%SizeDistName, SizeDists))
        If (SizeDist%ParticleShapePresent) Then
          If (Source%ShapeScheme == ' ') Then
            Call Message(                                              &
                   'FATAL ERROR: When a particle shape value is   ' // &
                   'defined, an associated "Sedimentation Scheme" ' // &
                   'must also be given.',                              &
                   3                                                   &
                 )
          End If
        End If
      End If
      If (Source%ParticleShapePresent) Then
        If (Source%ShapeScheme == ' ') Then
          Call Message(                                              &
                 'FATAL ERROR: When a particle shape value is   ' // &
                 'defined, an associated "Sedimentation Scheme" ' // &
                 'must also be given.',                              &
                 3                                                   &
               )
        End If
      End If
    End If 

    ! Check particle shape available for particulates if required, either via the source or the particle size distribution.
    ! Note Source%ShapeScheme, unlike Source%ShapeSchemeCode, distinguishes between ' ' and 'White'.
    If (Source%Particulate)  Then
      If (Source%ShapeScheme /= ' ')  Then
        If (Source%SizeDistName /= ' ') Then
          SizeDist => SizeDists%SizeDists(FindSizeDistIndex(Source%SizeDistName, SizeDists))
          If (.not. (Source%ParticleShapePresent .or. SizeDist%ParticleShapePresent)) Then
            Call Message(                                                     &
                   'FATAL ERROR: When a sedimentation scheme is defined, ' // &
                   'an associated "Particle Shape" must also be given.',      &
                   3                                                          &
                 )
          End If
        Else
          If (.not. Source%ParticleShapePresent) Then
            Call Message(                                                     &
                   'FATAL ERROR: When a sedimentation scheme is defined, ' // &
                   'an associated "Particle Shape" must also be given.',      &
                   3                                                          &
                 )
          End If
        End If
      Else
        Source%ParticleShape        = 1.0
        Source%ParticleShapePresent = .true.
      End If
    End If

  End Do

End Subroutine SetUpSources

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpiSources(Coords, Grids, Specieses, MaterialUnits, SizeDists, TimeDeps, Sources)
! Sets up indices in Sources.

  Implicit None
  ! Argument list:
  Type(Coords_),        Intent(In)            :: Coords        ! Collection of coord systems.
  Type(Grids_),         Intent(In)            :: Grids         ! Collection of grids.
  Type(Specieses_),     Intent(In)            :: Specieses     ! Collection of species.
  Type(MaterialUnits_), Intent(In)            :: MaterialUnits ! Collection of material units
  Type(SizeDists_),     Intent(In)            :: SizeDists     ! Collection of particle size distributions.
  Type(TimeDeps_),      Intent(In)            :: TimeDeps      ! Collection of source time dependencies.
  Type(Sources_),       Intent(InOut), Target :: Sources       ! Collection of sources.
  ! Locals:
  Type(Source_),  Pointer :: Source ! Abbreviation for sources.
  Integer                 :: i      !} Loop indices.
  Integer                 :: j      !}

  Do i = 1, Sources%nSources

    Source => Sources%Sources(i)

    ! Coords systems.
    If (Source%HCoordName  /= ' ') Source%iHCoord  = FindHCoordIndex(Source%HCoordName,  Coords)
    If (Source%ZCoordName  /= ' ') Source%iZCoord  = FindZCoordIndex(Source%ZCoordName,  Coords)
    If (Source%dHCoordName /= ' ') Source%idHCoord = FindHCoordIndex(Source%dHCoordName, Coords)
    If (Source%dZCoordName /= ' ') Source%idZCoord = FindZCoordIndex(Source%dZCoordName, Coords)

    ! Grids.
    If (Source%HGridName /= ' ') Source%iHGrid = FindHGridIndex(Source%HGridName, Grids)
    If (Source%ZGridName /= ' ') Source%iZGrid = FindZGridIndex(Source%ZGridName, Grids)
    If (Source%ZGridIterativePlumeRise /= ' ') Then
      Source%iZGridIterativePlumeRise = FindZGridIndex(Source%ZGridIterativePlumeRise, Grids)
    End If

    ! Species.
    Do j = 1, Source%nSpecieses
      Source%iSpecieses(j) = FindSpeciesIndex(Source%SpeciesNames(j), Specieses)
      Source%iParticleSpecieses(j) = Specieses%iSpecies2Particle( Source%iSpecieses(j) )
      If (Source%iParticleSpecieses(j) == 0) Then
        Call Message(                               &
               'FATAL ERROR: Species "'          // &
               Trim(Source%SpeciesNames(j))      // &
               '" in source "'                   // &
               Trim(Source%Name)                 // &
               '" is not carried on particles.',    &
               3                                    &
             )
      End If
    End Do

    ! Particle size distributions.
    If (Source%Particulate) Then
      If (Source%SizeDistName == ' ') Then
        Source%iSizeDist = 0
      Else
        Source%iSizeDist = FindSizeDistIndex(Source%SizeDistName, SizeDists)
      End If
    End If

    ! Source time dependencies.
    If (Source%TimeDepName == ' ') Then
      Source%iTimeDep = 0
    Else
      Source%iTimeDep = FindTimeDepIndex(Source%TimeDepName, TimeDeps)
    End If

    ! Source groups
    Do j = 1, Source%nSourceGroups
      Source%iSourceGroups(j) = FindSourceGroupIndex(Source%SourceGroupNames(j), Sources)
    End Do

    ! Material Unit names
    Do j = 1, Source%nSpecieses
      Source%iMaterialUnit(j) = FindMaterialUnitIndex(Source%MaterialUnitName(j),MaterialUnits)
    End Do
  End Do

End Subroutine SetUpiSources

!-------------------------------------------------------------------------------------------------------------

Subroutine FirstLastReleaseTimes(                &
             ReleaseTimeLimit, CompDomain,       &
             Sources,                            &
             FirstSourceTime,   LastSourceTime,  &
             FirstReleaseTime,  LastReleaseTime, &
             FirstReleaseTimes, LastReleaseTimes &
           )
! Determines the first source emission time and last possible source emission time for all sources, and the
! first release time and the last possible release time for each source and for all sources (accounting for
! the source information, the computational domain's time structure, and whether the met is fixed). Note that
! here we distinguish between source emission times and release times. A release time is the actual time of
! release of a particle/puff in the model and may correspond to a range of source emission times if the
! particle/puff has a time spread.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In)         :: ReleaseTimeLimit
  Type(Domain_),    Intent(In)         :: CompDomain
  Type(Sources_),   Intent(In), Target :: Sources
  Type(ShortTime_), Intent(Out)        :: FirstSourceTime
  Type(ShortTime_), Intent(Out)        :: LastSourceTime
  Type(ShortTime_), Intent(Out)        :: FirstReleaseTime
  Type(ShortTime_), Intent(Out)        :: LastReleaseTime
  Type(ShortTime_), Intent(Out)        :: FirstReleaseTimes(:)
  Type(ShortTime_), Intent(Out)        :: LastReleaseTimes(:)
  ! ReleaseTimeLimit  :: Upper limit for release times used for source emissions in the infinite past.
  !                      Infinite past or future values have no effect.
  ! CompDomain        :: The computational domain.
  ! Sources           :: Collection of sources.
  ! FirstSourceTime   :: First source emission time for all sources.
  ! LastSourceTime    :: Last possible source emission time for all sources.
  ! FirstReleaseTime  :: First release time for all sources.
  ! LastReleaseTime   :: Last possible release time for all sources.
  ! FirstReleaseTimes :: First release times for the sources.
  ! LastReleaseTimes  :: Last possible release time for the sources.
  ! Locals:
  Type(Source_),   Pointer :: Source      ! Abbreviation for source.
  Type(ShortTime_)         :: ReleaseTime ! Release time for sources with emissions in the infinite past.
  Integer                  :: i           ! Loop index.

  ! 1) Check array sizes.

  If (Size(FirstReleaseTimes) < Sources%nSources .or. Size(LastReleaseTimes) < Sources%nSources) Then
    Call Message('UNEXPECTED FATAL ERROR in FirstLastReleaseTimes', 4)
  End If

  ! 2) Find first and last possible source emission time for each source, accounting for the start and end
  ! times of the source and the computational domain and the fact that the the source might not actually emit
  ! (due to the start time, stop time, max travel time or source strength of the source or due to the
  ! computational domain).

  Do i = 1, Sources%nSources

    Source => Sources%Sources(i)

    If (IsInfFuture(Source%sStartTime)) Then

      Call Message(                                           &
             'Warning: Start time for source "'            // &
             Trim(Source%Name)                             // &
             '" is infinity. The source will be ignored.',    &
             1                                                &
           )
      FirstReleaseTimes(i) = InfFutureShortTime()
      LastReleaseTimes (i) = InfPastShortTime  ()

    Else If (IsInfPast(Source%sStopTime)) Then

      Call Message(                                            &
             'Warning: Stop time for source "'              // &
             Trim(Source%Name)                              // &
             '" is -infinity. The source will be ignored.',    &
             1                                                 &
           )
      FirstReleaseTimes(i) = InfFutureShortTime()
      LastReleaseTimes (i) = InfPastShortTime  ()

    Else If (Source%sStartTime > Source%sStopTime) Then

      Call Message(                                                                  &
             'Warning: The start time is greater than the stop time for source "' // &
             Trim(Source%Name)                                                    // &
             '". The source will be ignored.',                                       &
             1                                                                       &
           )
      FirstReleaseTimes(i) = InfFutureShortTime()
      LastReleaseTimes (i) = InfPastShortTime  ()

    Else If (Source%sMaxTravelTime == ZeroShortTime()) Then

      Call Message(                                       &
             'Warning: Max travel time for source "'   // &
             Trim(Source%Name)                         // &
             '" is zero. The source will be ignored.',    &
             1                                            &
           )
      FirstReleaseTimes(i) = InfFutureShortTime()
      LastReleaseTimes (i) = InfPastShortTime  ()

    Else If (All(Source%SourceStrengths(1:Source%nSpecieses) == 0)) Then

      Call ControlledMessage(                                                                              &
             'All source strengths for source "'                                                        // &
             Trim(Source%Name)                                                                          // &
             '" are zero (possibly because its specified as a rate or as met dependent for an '         // &
             'instantaneous source, or as total amount released for an infinite duration source). The ' // &
             'source will be ignored.',                                                                    &
             MessageControls    = GlobalMessageControls,                                                   &
             MessageControlName = 'Zero source strength',                                                  &
             ErrorCode          = 1                                                                        &
           )
      FirstReleaseTimes(i) = InfFutureShortTime()
      LastReleaseTimes (i) = InfPastShortTime  ()

    Else If (                                                                                       &
      Source%sStartTime >= DomainEnd  (CompDomain)                                             .or. &
      Source%sStopTime  <  DomainStart(CompDomain)                                             .or. &
      (Source%sStopTime == DomainStart(CompDomain) .and. Source%sStartTime < Source%sStopTime) .or. &
      .not.SInDomain(ZeroShortTime(), CompDomain)                                                   &
    ) Then

      FirstReleaseTimes(i) = InfFutureShortTime()
      LastReleaseTimes (i) = InfPastShortTime  ()

    Else

      FirstReleaseTimes(i) = TMax(Source%sStartTime, DomainStart(CompDomain))
      LastReleaseTimes (i) = TMin(Source%sStopTime,  DomainEnd  (CompDomain))

    End If

  End Do

  ! 3) Find first and last possible source emission time for all sources.

  FirstSourceTime = InfFutureShortTime()
  LastSourceTime  = InfPastShortTime()
  Do i = 1, Sources%nSources
    FirstSourceTime = TMin(FirstSourceTime, FirstReleaseTimes(i))
    LastSourceTime  = TMax(LastSourceTime,  LastReleaseTimes (i))
  End Do

  ! 4) Find first and last possible release time for each source and deal with sources with start time in the
  ! infinite past. For fixed met, the first release time for sources with start time in the infinite past is
  ! set to the earliest of the finite values of FirstReleaseTimes, LastReleaseTimes and ReleaseTimeLimit (if
  ! any) or to a default value. For non-fixed met we have a fatal error. For fixed met, the last possible
  ! release time is set to the first release time. (For fixed met, the particles/puffs for a source are all
  ! released at the same time, which can be chosen arbitrarily within the source period and is representative
  ! of all release times).

  If (Sources%nSources > 0) Then

    If (Source%FixedMet) Then

      ReleaseTime = InfFutureShortTime()
      Do i = 1, Sources%nSources
        If (.not.IsInfPast(FirstReleaseTimes(i))) ReleaseTime = TMin(ReleaseTime, FirstReleaseTimes(i))
        If (.not.IsInfPast(LastReleaseTimes (i))) ReleaseTime = TMin(ReleaseTime, LastReleaseTimes (i))
      End Do
      If (.not.IsInfPast(ReleaseTimeLimit)) ReleaseTime = TMin(ReleaseTime, ReleaseTimeLimit)
      If (IsInfFuture(ReleaseTime)) ReleaseTime = ReferenceShortTime()

      Do i = 1, Sources%nSources
        If (.not.IsInfFuture(FirstReleaseTimes(i))) Then
          If (IsInfPast(FirstReleaseTimes(i))) FirstReleaseTimes(i) = ReleaseTime
          LastReleaseTimes(i) = FirstReleaseTimes(i)
        End If
      End Do

    Else

      Do i = 1, Sources%nSources
        If (IsInfPast(FirstReleaseTimes(i))) Then
          Call Message(                                                       &
                 'FATAL ERROR: For non-fixed met, either the '             // &
                 'start time of the computational domain or the start '    // &
                 'time of every source must not be in the infinite past.',    &
                 3                                                            &
               )
        End If
      End Do

    End If

  End If

  ! 5) Find first and last possible release time for all sources.

  FirstReleaseTime = InfFutureShortTime()
  LastReleaseTime  = InfPastShortTime()
  Do i = 1, Sources%nSources
    FirstReleaseTime = TMin(FirstReleaseTime, FirstReleaseTimes(i))
    LastReleaseTime  = TMax(LastReleaseTime,  LastReleaseTimes (i))
  End Do

End Subroutine FirstLastReleaseTimes

!-------------------------------------------------------------------------------------------------------------

Subroutine ParticleReleaseInfo(                                                &
             Time, Source,                                                     &
             ParticleCeiling, ParticleFactor, MaxParticles, nOldParticles,     &
             Specieses, SizeDists, TimeDeps,                                   &
             Flow, Rain, Surface, Soil, Area, iSizeRange,                      &
             Mass, Diameter, Density, ParticleShape, Rate,                     &
             nNewParticles, NextTime, ParticleFactorType, Qm0                  &
           )
! Calculates, for a particle release, the mass to release (for each species), the
! number of particles to release, and the next release time.

! Note source duration is [StartTime, StopTime] although, for particles with varying
! met, each particle is released at the start of the source period relevant to its
! mass, and so, for non-instantaneous sources, no particle can be released at
! StopTime.

! The following method is adopted for adjusting the number of particles released and
! the mass per particle (but not the total mass released) if nearing the maximum
! number of particles allowed, to try to prevent the model running out of particles.
!
! Let c = ParticleCeiling, pf = ParticleFactor and n = number of particles after the
! next particle is released (= number of existing particles + 1). Then mass per
! particle is increaed by a factor
!     f = { 1   for n <= c
!         { pf  otherwise.
! The cumulative mass-increase-factor per particle y(N) = sum_1^N f is given by
!     y = { N               for N <= c
!         { c + pf (N - c)  otherwise
! with inverse
!     N = { y                 for y <= c
!         { (y - c) / pf + c  otherwise.
!
! If several particles are released at once, the number of particles is determined as
! above, but the mass is distributed evenly across the particles.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In)  :: Time
  Type(Source_),    Intent(In)  :: Source
  Integer,          Intent(In)  :: ParticleCeiling
  Real(Std),        Intent(In)  :: ParticleFactor
  Integer,          Intent(In)  :: MaxParticles
  Integer,          Intent(In)  :: nOldParticles
  Type(Specieses_), Intent(In)  :: Specieses
  Type(SizeDists_), Intent(In)  :: SizeDists
  Type(TimeDeps_),  Intent(In)  :: TimeDeps
  Type(Flow_),      Intent(In)  :: Flow
  Type(Rain_),      Intent(In)  :: Rain
  Type(Surface_),   Intent(In)  :: Surface
  Type(Soil_),      Intent(In)  :: Soil
  Real(Std),        Intent(In)  :: Area
  Integer,          Intent(In)  :: iSizeRange
  Real(Std),        Intent(Out) :: Mass(MaxSpecieses)
  Real(Std),        Intent(Out) :: Diameter
  Real(Std),        Intent(Out) :: Density
  Real(Std),        Intent(Out) :: ParticleShape
  Logical,          Intent(Out) :: Rate
  Integer,          Intent(Out) :: nNewParticles
  Type(ShortTime_), Intent(Out) :: NextTime
  Integer,          Intent(Out) :: ParticleFactorType
  Real(Std),        Intent(In)  :: Qm0
  ! Time               :: Time of release. On first call (for a given source) this
  !                       should be the first release time returned by
  !                       FirstLastReleaseTimes. On subsequent calls NextTime from the
  !                       previous call should be used. If Time calculated this way is
  !                       infinity (indicating no release) this routine should not be
  !                       called).
  ! Source             :: Source.
  ! ParticleCeiling    :} Particle number ceiling above which the number of particles
  !                    :} released and the mass per particle (but not the total mass
  !                       released) is adjusted to try to prevent the model running
  !                       out of particles.
  ! ParticleFactor     :: Factor by which the mass per particle is increased above
  !                       ParticleCeiling.
  ! MaxParticles       :: Maximum number of particles allowed.
  ! nOldParticles      :: Number of particles currently used.
  ! Specieses          :: Collection of specieses.
  ! SizeDists          :: Collection of particle size distributions.
  ! TimeDeps           :: Collection of source time dependencies.
  ! Flow               :: Flow information.
  ! Rain               :: Rain information.
  ! Surface            :: Surface information.
  ! Soil               :: Soil information.
  ! Area               :: Area of grid cell.
  ! iSizeRange         :: Index of range within the particle size distribution. Used when a release from a
  !                       specified size range is required.
  ! Mass               :: Mass to release or mass release rate (for each species).
  ! Diameter           :} Diameter, density and shape of particle to release.
  ! Density            :}
  ! ParticleShape      :}
  ! Rate               :: Indicates Mass is a release rate. Is true for
  !                       non-instantaneous sources with fixed met and false
  !                       otherwise.
  ! nNewParticles      :: Number of particles to release (can be zero due to source
  !                       time dependency, due to Source%iCritSpecies =
  !                       Source%nParticles = 0, or due to MaxParticles).
  ! ParticleFactorType :: Indicates type of adjustment to number of particles released
  !                       to try to prevent the model running out of particles:
  !                           0 - no adjustment,
  !                           1 - ParticleCeiling reached, numbers reduced and masses
  !                               increased,
  !                           2 - MaxParticles reached and release stopped.
  !                       Note 0 may mean ParticleCeiling or MaxParticles reached but
  !                       there would be nothing to release anyway and so no
  !                       adjustment has been made.
  ! NextTime           :: Time of following release or infinity if no following
  !                       release.
  ! Qm0                :: Source strength (mass flux) computed from iterative plume model
  !
  ! Note Flow, Surface, Soil, Area and iSizeRange are only used for dust and sea-salt sources.
  ! $$ could treat iSizeRange similarly for generic sources (possibly as an option).
  ! Local parameters:
  Real(Std), Parameter :: Beta = 0.09  ! Unit K-1 empirical coefficient Guenther et al 1993.
  Real(Std), Parameter :: Ts   = 303.0 ! Unit K standard temperature Guenther et al 1993.
  ! Locals:
  Real(Std)        :: Factor                     ! Source time dependency factor.
  Real(Std)        :: PFactor                    ! Factor to adjust particle release rate to try to prevent
                                                 ! the model running out of particles. Values > 1 reduce
                                                 ! particle numbers and increase mass per particle.
  Integer          :: P1                         !} Temporary variables related to particle numbers.
  Integer          :: P2                         !}
  Real(Std)        :: ReleaseRate                ! Particle release rate.
  Type(ShortTime_) :: ReleaseInterval            ! Particle release interval.
  Real(Std)        :: SourceStrength             ! Source strength for met-dependent sources.
  Integer          :: iSpecies                   ! Species index in set of all species.
  Integer          :: iParticleSpecies           ! Species index in set of species on particles.
  Integer          :: i                          ! Loop index.
  Real(Std)        :: TFactor(Source%nSpecieses) ! Met dependent source factor.
  Type(ShortTime_) :: StartOfYear                ! Start time of calendar year.
  Real(Std)        :: JulDay                     ! Julian day   
  Real(Std)        :: phi                        ! phase of year
  Real(Std)        :: logMidge                   ! log of midge number
  Real(Std), Parameter :: b0        = -0.89 ! Intercept parameter from midge dynamic equation
  Real(Std), Parameter :: b11       = -1.25 ! Sin, 12 month parameter from midge dynamics equation
  Real(Std), Parameter :: b21       = -3.21 ! Cos, 12 month parameter from midge dynamics equation
  Real(Std), Parameter :: b12       = -1.20 ! Sin, 6 month parameter from midge dynamics equation
  Real(Std), Parameter :: b22       = -0.62 ! Sin, 6 month parameter from midge dynamics equation
  Real(Std), Parameter :: c1        =  0.10 ! Temperature parameter from midge dynamics equation
  Real(Std), Parameter :: c2        = -0.14 ! Wind speed parameter from midge dynamics equation
  Real(Std), Parameter :: c3        =  0.04 ! Precipitation parameter from midge dynamics equation
  Real(Std), Parameter :: UStarCrit =  0.4  ! Threshold friction velocity for ash resuspension.

# ifdef ExtraChecks
    If (                                                                         &
      Time < Source%sStartTime                                              .or. &
      Time > Source%sStopTime                                               .or. &
      (Time == Source%sStopTime .and. Source%sStartTime < Source%sStopTime) .or. &
      IsInfFuture(Time)                                                          &
    ) Then
      Call Message('UNEXPECTED FATAL ERROR in ParticleReleaseInfo', 4)
    End If
# endif

  ! Default for ParticleFactorType.
  ParticleFactorType = 0

  !------------------------------------!
  ! Particle diameter, density, shape. !
  !------------------------------------!

  If (Source%Particulate) Then
    If (Source%iSizeDist > 0) Then
      Call CalculateParticleDiameterEtc(Source, SizeDists, iSizeRange, Diameter, Density, ParticleShape)
    Else
      Diameter      = Source%Diameter
      Density       = Source%Density
      ParticleShape = Source%ParticleShape
    End If
  Else
    Diameter = 0.0
  End If

  !------------------------------------------------------------------------!
  ! Nothing to release due to Source%iCritSpecies = Source%nParticles = 0. !
  !------------------------------------------------------------------------!

  If (Source%iCritSpecies == 0 .and. Source%nParticles == 0.0) Then

    If (Source%FixedMet) Then
      Call Message(                                                         &
             'Warning: No particles are due to be released for source "' // &
             Trim(Source%Name)                                           // &
             '". The source will be ignored.',                              &
             1                                                              &
           )
    Else
      Call Message(                                                                &
             'Warning: No particles are due to be released for source "'        // &
             Trim(Source%Name)                                                  // &
             '" (possibly because "# particles" is specified as a rate for an ' // &
             'instantaneous source or as a total for an infinite duration '     // &
             'source). The source will be ignored.',                               &
             1                                                                     &
           )
    End If

    NextTime      = InfFutureShortTime()
    nNewParticles = 0

  !------------------------------------------------------------------------!
  ! Fixed Met. Here non-instantaneous sources have square-shaped t-spread. !
  !------------------------------------------------------------------------!

  Else If (Source%FixedMet) Then

    ! $$ Dust/sea-salt sources are not currently handled for fixed met cases
    If (Source%SourceType == S_Dust .or. Source%SourceType == S_SeaSalt) Then
      Call Message(                                                                         &
             'Error in ParticleReleaseInfo: a dust or sea-salt source is not supported ' // &
             'for the fixed met option',                                                    &
             3                                                                              &
           )
    End If

    ! SourceTimeDependencyFactor - needs treating elsewhere $$

    NextTime = InfFutureShortTime()

    nNewParticles = Ceiling(Source%nParticles)

    ! Adjust nNewParticles to try to prevent the model running out of particles.
    P1 = nOldParticles                 - ParticleCeiling
    P2 = nOldParticles + nNewParticles - ParticleCeiling
    If (P1 >= 0) Then
      nNewParticles = Ceiling(Real(nNewParticles, Std) / ParticleFactor)
      ParticleFactorType = 1
    Else If (P2 > 0) Then
      nNewParticles = Ceiling(Real(P2, Std) / ParticleFactor - Real(P1, Std))
      ParticleFactorType = 1
    End If

    Mass(:) = 0.0
    Do i = 1, Source%nSpecieses
      iParticleSpecies       = Source%iParticleSpecieses(i)
      Mass(iParticleSpecies) = Source%SourceStrengths(i) / nNewParticles
    End Do

    Rate = Source%SourceRate(i)

    ! If exceed MaxParticles, reduce nNewParticles and set ParticleFactorType to 2.
    If (nOldParticles + nNewParticles > MaxParticles) Then
      nNewParticles      = MaxParticles - nOldParticles
      ParticleFactorType = 2
    End If

  !---------------------------------------------------------------!
  ! Non-fixed met, instantaneous sources. These have no t-spread. !
  !---------------------------------------------------------------!

  Else If (Source%sStopTime == Source%sStartTime) Then

    If (Source%iTimeDep == 0) Then
      Factor = 1.0
    Else
      Call SourceTimeDependencyFactor(                                         &
             ShortTime2Time(Time), Source, TimeDeps%TimeDeps(Source%iTimeDep), &
             Factor                                                            &
           )
    End If

    If (Factor == 0.0) Then

      NextTime = InfFutureShortTime()

      nNewParticles = 0

    Else

      NextTime = InfFutureShortTime()

      If (Source%iCritSpecies == 0) Then
        nNewParticles = Ceiling(Source%nParticles)
      Else
        iSpecies   = Source%iSpecieses(Source%iCritSpecies)
        nNewParticles = Ceiling(                                                   &
                          Max(                                                     &
                            Source%nParticles,                                     &
                            Factor * Source%SourceStrengths(Source%iCritSpecies) / &
                            Specieses%Specieses(iSpecies)%MassLimit                &
                          )                                                        &
                        )
      End If

      ! Adjust nNewParticles to try to prevent the model running out of particles.
      P1 = nOldParticles                 - ParticleCeiling
      P2 = nOldParticles + nNewParticles - ParticleCeiling
      If (P1 >= 0) Then
        nNewParticles = Ceiling(Real(nNewParticles, Std) / ParticleFactor)
        ParticleFactorType = 1
      Else If (P2 > 0) Then
        nNewParticles = Ceiling(Real(P2, Std) / ParticleFactor - Real(P1, Std))
        ParticleFactorType = 1
      End If

      Mass(:) = 0.0
      Do i = 1, Source%nSpecieses
        iParticleSpecies       = Source%iParticleSpecieses(i)
        Mass(iParticleSpecies) = Factor * Source%SourceStrengths(i) / nNewParticles
      End Do

      Rate = .false.

      ! If exceed MaxParticles, reduce nNewParticles and set ParticleFactorType to 2.
      If (nOldParticles + nNewParticles > MaxParticles) Then
        nNewParticles      = MaxParticles - nOldParticles
        ParticleFactorType = 2
      End If

    End If

  !----------------------------------------------------------------------------------!
  ! Non-fixed met, non-instantaneous dust/sea-salt sources. These have no t-spread.  !
  ! Note the mass released is applicable to the meteorology and surface conditions   !
  ! at the actual release time (rather than averaged over the period that the        !
  ! particle represents).                                                            !
  ! Note: Dust/sea-salt sources are not presently considered for fixed met cases.    !
  !----------------------------------------------------------------------------------!

  Else If (Source%SourceType == S_Dust .or. Source%SourceType == S_SeaSalt) Then

    If (Source%SourceType == S_Dust) Then
      ! Calculate dust source strength (in g).
      ! This is the source strength for a particular grid cell and size range.
      Call CalcDustSourceStrength(                &
             Flow, Surface, Soil, Area, Diameter, &
             iSizeRange, SourceStrength           &
           )
    Else
      ! Calculate sea-salt source strength (in g).
      ! This is the source strength for a particular grid cell and size range.
      Call CalcSeaSaltSourceStrength(                   &
             Flow, Surface, Area, Diameter, iSizeRange, &
             SourceStrength                             &
           )
    End If

    If (SourceStrength <= 0.0) Then

      nNewParticles = 0

      NextTime = Time + Time2ShortTime(Char2Time('00:05', Interval = .true.))
      If (NextTime >= Source%sStopTime) Then
        NextTime = InfFutureShortTime()
      End If

    Else

      ! Factor to adjust particle release rate to try to prevent the model running out
      ! of particles.
      If (nOldParticles < ParticleCeiling) Then
        PFactor = 1.0
      Else
        PFactor            = ParticleFactor
        ParticleFactorType = 1
      End If

      If (Source%iCritSpecies == 0) Then
        ReleaseRate = Source%nParticles
      Else
        Call Message(                                                                                        &
               'Error in ParticleReleaseInfo: the particle mass limit option is currently not supported ' // &
               'for a dust or sea-salt source',                                                              &
               3                                                                                             &
             )
        ! iSpecies    = Source%iSpecieses(Source%iCritSpecies)
        ! ReleaseRate = Max(                                                     &
        !                 Source%nParticles,                                     &
        !                 Factor * Source%SourceStrengths(Source%iCritSpecies) / &
        !                 Specieses%Specieses(iSpecies)%MassLimit                &
        !               )
      End If

      nNewParticles = 1

      ! (a) The model hasn't run out of particles.
      If (nOldParticles + nNewParticles <= MaxParticles) Then

        ReleaseRate     = ReleaseRate / PFactor
        ReleaseInterval = RealTime2ShortTime(1.0 / ReleaseRate)

        ! $$ if ReleaseInterval = 0 (possible for very short sources) set to Tiny and
        ! re-evaluate ReleaseRate. This imposes a max
        ! particle release rate. (Check if this is possible earlier and give warning?) $$

        NextTime = Time + ReleaseInterval

        If (NextTime >= Source%sStopTime) Then
          ReleaseInterval = Source%sStopTime - Time
          ReleaseRate     = 1.0 / ShortTime2RealTime(ReleaseInterval)
          NextTime        = InfFutureShortTime()
        End If

        Mass(:) = 0.0
        Do i = 1, Source%nSpecieses
          iParticleSpecies       = Source%iParticleSpecieses(i)
          Mass(iParticleSpecies) = SourceStrength / ReleaseRate
        End Do

        Rate = .false.

      ! (b) The model has run out of particles. Reduce nNewParticles and set
      ! ParticleFactorType to 2. Also evaluate NextTime so that the next chance to
      ! release isn't delayed by PFactor.
      Else

        nNewParticles      = 0
        ParticleFactorType = 2

        ReleaseInterval = RealTime2ShortTime(1.0 / ReleaseRate)

        ! $$ if ReleaseInterval = 0 (possible for very short sources) set to Tiny and
        ! re-evaluate ReleaseRate. This imposes a max
        ! particle release rate. (Check if this is possible earlier and give warning?) $$

        NextTime = Time + ReleaseInterval

        If (NextTime >= Source%sStopTime) Then
          NextTime = InfFutureShortTime()
        End If

      End If

    End If

  !----------------------------------------------------------------------------------!
  ! Non-fixed met, non-instantaneous sources. These have no t-spread. Note the       !
  ! source time dependency factor is applied only at the actual release time (rather !
  ! than averaged over the period that the particle represents).                     !
  !----------------------------------------------------------------------------------!

  Else

    If (Source%iTimeDep == 0) Then
      Factor = 1.0
    Else
      Call SourceTimeDependencyFactor(                                         &
             ShortTime2Time(Time), Source, TimeDeps%TimeDeps(Source%iTimeDep), &
             Factor                                                            &
           )
    End If

    If (Factor == 0.0) Then

      Call SourceTimeDependencyFactor(                                         &
             ShortTime2Time(Time), Source, TimeDeps%TimeDeps(Source%iTimeDep), &
             Factor, NextTime                                                  &
           )

      nNewParticles = 0

      If (NextTime >= Source%sStopTime) Then
        NextTime = InfFutureShortTime()
      End If

    Else

      ! Factor to adjust particle release rate to try to prevent the model running out
      ! of particles.
      If (nOldParticles < ParticleCeiling) Then
        PFactor = 1.0
      Else
        PFactor            = ParticleFactor
        ParticleFactorType = 1
      End If

      If (Source%iCritSpecies == 0) Then
        ReleaseRate = Source%nParticles
      Else
        iSpecies    = Source%iSpecieses(Source%iCritSpecies)

! $$ nb: iCritSpecies is calculated before we adjust the biogenic emission as the met
! is not available at the time 'setupsources_speciesgrids' subroutine is run. This may result in
! weighting too many particles on the biogenic emission as the actual emission is likely to be less
! than in the sources file. Sources file must be set to release at least 1 particle per timestep.
!        print*,'species limiting particle release=',ispecies,Source%SourceStrengths(Source%iCritSpecies)

        ReleaseRate = Max(                                                     &
                        Source%nParticles,                                     &
                        Factor * Source%SourceStrengths(Source%iCritSpecies) / &
                        Specieses%Specieses(iSpecies)%MassLimit                &
                      )
      End If

      nNewParticles = 1

      ! (a) The model hasn't run out of particles.
      If (nOldParticles + nNewParticles <= MaxParticles) Then

        ReleaseRate     = ReleaseRate / PFactor
        ReleaseInterval = RealTime2ShortTime(1.0 / ReleaseRate)

        ! $$ if ReleaseInterval = 0 (possible for very short sources) set to Tiny and
        ! re-evaluate ReleaseRate. This imposes a max
        ! particle release rate. (Check if this is possible earlier and give warning?) $$

        NextTime = Time + ReleaseInterval

        If (NextTime >= Source%sStopTime) Then
          ReleaseInterval = Source%sStopTime - Time
          ReleaseRate     = 1.0 / ShortTime2RealTime(ReleaseInterval)
          NextTime        = InfFutureShortTime()
        End If

        ! TFactor (= 1 except for met dependent sources).
        Do i = 1, Source%nSpecieses
          TFactor(i) = 1.0
        End Do

        ! Met dependent sources.
        If (Source%MetDependent == 1) Then
          
          Do i = 1, Source%nSpecieses
            iSpecies = Source%iSpecieses(i)

            If (Specieses%Specieses(iSpecies)%Name .CIEq. 'MIDGE') Then

              ! Calculating Julian Day (as a real number) from time
              StartOfYear = RoundToYear(Time, .false., .false.)
              JulDay      = ShortTime2RealTime(Time - StartOfYear) / 86400.0 

              ! Applying Simon Gubbins/Chris Sanders Midge Dynamics Equation
              ! to determine 'mass' of midge population becoming airborne
              phi = (2.0 * Pi * JulDay) / 366.0
              logMidge = b0 + (b11*sin(phi))     + (b21*cos(phi))                 &
                            + (b12*sin(2.0*phi)) + (b22*cos(2.0*phi))             &
                            + (c1*(Flow%T - 273.15)) + (c2*((Flow%U(1)**2 + Flow%U(2)**2)**0.5))
              If (Rain%ConPpt+Rain%DynPpt > 0.0) logMidge = logMidge + c3

              TFactor(i) = Exp(logMidge)

            Else If (Specieses%Specieses(iSpecies)%Name .CIEq. 'C10H16bio') Then

              If (Source%iBioAPinene > 0) Then
                TFactor(Source%iBioAPinene) = Exp(Beta * (Flow%T - Ts))
              Else
                print *, 'source%iBioAPinene less then 0:', Source%iBioAPinene ! $$ del in due course?
              End If

            Else If (Specieses%Specieses(iSpecies)%Name .CIEq. 'RESUSPENDED_ASH') Then

              If (Flow%UStar > UStarCrit .and. Rain%ConPpt+Rain%DynPpt < 0.01) Then
                TFactor(i) = (Flow%UStar - UStarCrit)**3
              Else
                TFactor(i) = 0.0
              End If

            End If

          End Do

        End If

! This bit of code allows the ammonia emission cycle to only be applied to ammonia and
! not all the other species that are on the same particle - traffic cycle is applied to 
! all species on the source line

        Mass(:) = 0.0
        Do i = 1, Source%nSpecieses
          iParticleSpecies = Source%iParticleSpecieses(i)
          If (Source%TimeDepName .CIEq. 'AMMONIA') Then
            If (Source%SpeciesNames(i) .CIEq. 'AMMONIA') Then
              Mass(iParticleSpecies) = Factor * TFactor(i) * Source%SourceStrengths(i) / ReleaseRate
            Else
              Mass(iParticleSpecies) = TFactor(i) * Source%SourceStrengths(i) / ReleaseRate
            Endif
          Else If ( Source%MetDependent == 2 ) Then
            Mass(iParticleSpecies) = Factor * TFactor(i) * Qm0 / ReleaseRate   
          Else
            Mass(iParticleSpecies) = Factor * TFactor(i) * Source%SourceStrengths(i) / ReleaseRate
          Endif
        End Do

        Rate = .false.

      ! (b) The model has run out of particles. Reduce nNewParticles and set
      ! ParticleFactorType to 2. Also evaluate NextTime so that the next chance to
      ! release isn't delayed by PFactor.
      Else

        nNewParticles      = 0
        ParticleFactorType = 2

        ReleaseInterval = RealTime2ShortTime(1.0 / ReleaseRate)

        ! $$ if ReleaseInterval = 0 (possible for very short sources) set to Tiny and
        ! re-evaluate ReleaseRate. This imposes a max
        ! particle release rate. (Check if this is possible earlier and give warning?) $$

        NextTime = Time + ReleaseInterval

        If (NextTime >= Source%sStopTime) Then
          NextTime = InfFutureShortTime()
        End If

      End If

    End If

  End If

End Subroutine ParticleReleaseInfo

!-------------------------------------------------------------------------------------------------------------

Subroutine PuffReleaseInfo(              &
             TimeM, Time, TimeP, Source, &
             Specieses, TimeDeps,        &
             Mass, Rate                  &
           )
! Calculates, for a puff release, the mass to release (for each species).

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In)  :: TimeM
  Type(ShortTime_), Intent(In)  :: Time
  Type(ShortTime_), Intent(In)  :: TimeP
  Type(Source_),    Intent(In)  :: Source
  Type(Specieses_), Intent(In)  :: Specieses
  Type(TimeDeps_),  Intent(In)  :: TimeDeps
  Real(Std),        Intent(Out) :: Mass(MaxSpecieses)
  Logical,          Intent(Out) :: Rate
  ! TimeM     :: Trailing time edge of puff.
  ! Time      :: Time of release.
  ! TimeP     :: Leading time edge of puff.
  ! Source    :: Source.
  ! Specieses :: Collection of specieses.
  ! TimeDeps  :: Collection of source time dependencies.
  ! Mass      :: Mass to release or mass release rate (for each species).
  ! Rate      :: Indicates Mass is a release rate. Is true for non-instantaneous
  !              sources with fixed met and false otherwise.
  ! Locals:
  Real(Std) :: Factor           ! Source time dependency factor.
  Integer   :: iSpecies         ! Species index in set of all species.
  Integer   :: iParticleSpecies ! Species index in set of species on particles. 
  Integer   :: i                ! Loop index.

# ifdef ExtraChecks
    If (                                                               &
      Source%sStartTime > TimeM                                   .or. &
      TimeM             > Time                                    .or. &
      Time              > TimeP                                   .or. &
      TimeP             > Source%sStopTime                        .or. &
      (TimeM == TimeP .and. Source%sStartTime < Source%sStopTime) .or. &
      IsInfFuture(Time)                                                &
    ) Then
      Call Message('UNEXPECTED ERROR in PuffReleaseInfo', 4)
    End If
# endif

  ! Fixed Met. Here non-instantaneous sources have square-shaped t-spread.
  If (Source%FixedMet) Then

    ! SourceTimeDependencyFactor - needs treating elsewhere $$

    Mass(:) = 0.0
    Do i = 1, Source%nSpecieses
      iParticleSpecies       = Source%iParticleSpecieses(i)
      Mass(iParticleSpecies) = Source%SourceStrengths(i)
    End Do

    Rate = Source%SourceRate(i)

  ! Instantaneous sources.
  Else If (Source%sStopTime == Source%sStartTime) Then

    If (Source%iTimeDep == 0) Then
      Factor = 1.0
    Else
      Call SourceTimeDependencyFactor(                                         &
             ShortTime2Time(Time), Source, TimeDeps%TimeDeps(Source%iTimeDep), &
             Factor                                                            &
           )
    End If

    Mass(:) = 0.0
    Do i = 1, Source%nSpecieses
      iParticleSpecies       = Source%iParticleSpecieses(i)
      Mass(iParticleSpecies) = Factor * Source%SourceStrengths(i)
    End Do

    Rate = .false.

  ! Non-instantaneous sources. These have triangle-shaped t-spread. Note source time
  ! dependency factor is applied only at the actual release time (rather than averaged
  ! over the period that the puff represents).
  Else

    If (Source%iTimeDep == 0) Then
      Factor = 1.0
    Else
      Call SourceTimeDependencyFactor(                                         &
             ShortTime2Time(Time), Source, TimeDeps%TimeDeps(Source%iTimeDep), &
             Factor                                                            &
           )
    End If

    Mass(:) = 0.0
    Do i = 1, Source%nSpecieses
      iParticleSpecies       = Source%iParticleSpecieses(i)
      Mass(iParticleSpecies) = Factor * Source%SourceStrengths(i) * &
                       0.5 * ShortTime2RealTime(TimeP - TimeM)
    End Do

    Rate = .false.

  End If

End Subroutine PuffReleaseInfo

!-------------------------------------------------------------------------------------------------------------

Subroutine CalculateParticleDiameterEtc(Source, SizeDists, iSizeRange, Diameter, Density, ParticleShape)
! Calculates the diameter, density and particle shape parameter of particulates to release.

! For generic sources a random diameter (log-uniform distribution) from a random size range is returned.
! For dust sources the mid point of a specified size range is returned.
! For sea-salt sources a random diameter (polynomial distribution) from a specified size range is returned.

! For the sea-salt case there are only two size ranges. The lower is solved by a bisection method and the
! upper by a Newton-Raphson iteration.

! $$ The generic and dust cases could be combined/generalised and made switchable from the input file.
! The sea-salt case requires exactly two ranges and specified end points (although in fact the input values
! are ignored). The number of ranges should be checked and a message issues about range boundaries being
! ignored. For dust a message about range quantities (but not boundaries) being ignored would be useful).

! For particle density and particle shape, the values are taken from the relevant particle size range bin
! (if these characteristics have been defined as part of the particle size distribution) and otherwise
! are taken from the (size-independent) characteristics defined in the source term.

  Implicit None
  ! Argument list:
  Type(Source_),    Intent(In)         :: Source          ! Source.
  Type(SizeDists_), Intent(In), Target :: SizeDists       ! Collection of particle size distributions.
  Integer,          Intent(In)         :: iSizeRange      ! Index of size range within the particle size
                                                          ! distribution (as specified by the calling routine).
  Real(Std),        Intent(Out)        :: Diameter        ! Diameter of particle to release.
  Real(Std),        Intent(Out)        :: Density         !} Particle Shape and density
  Real(Std),        Intent(Out)        :: ParticleShape   !}
  ! Locals:
  Type(SizeDist_), Pointer :: SizeDist              ! Abbreviation for particle size distribution.
  Integer                  :: jSizeRange            ! Index of size range within the particle size
                                                    ! distribution (as calculated within this routine).
  Real(Std)                :: LowerProbabilityBound !} Lower and upper bounds on the probability and diameter
  Real(Std)                :: UpperProbabilityBound !} for a particle size distribution.
  Real(Std)                :: LowerDiameterBound    !}
  Real(Std)                :: UpperDiameterBound    !}
  Real(Std)                :: U                     ! Uniform random number.
  Real(Std)                :: A0                    !} Used in iteratively solving the size of the particle
  Real(Std)                :: Top                   !} for the sea-salt distribution.
  Real(Std)                :: Bottom                !}
  Real(Std)                :: Radius                !}
  Real(Std)                :: Value                 !}
  Real(Std)                :: Derivative            !}
  Integer                  :: i                     ! Loop index.

  SizeDist => SizeDists%SizeDists(Source%iSizeDist)

  If (Source%SourceType == S_Generic .Or. Source%SourceType == S_IterativePlumeModel) Then

    ! Determine size range.
    Call GetRandomNumber(U, 0)
    jSizeRange = 1 ! Probably unnecessary, but avoids any chance of jSizeRange not being set due to rounding.
    Do i = 1, SizeDist%nSizeRanges
      LowerProbabilityBound = SizeDist%CumulativeFrac(i)
      UpperProbabilityBound = SizeDist%CumulativeFrac(i + 1)
      If (LowerProbabilityBound <= U .and. U <= UpperProbabilityBound) Then
        jSizeRange = i
        Exit
      End If
    End Do

    ! Determine diameter bounds for the size range.
    LowerDiameterBound = SizeDist%DiameterRangeBoundary(jSizeRange)
    UpperDiameterBound = SizeDist%DiameterRangeBoundary(jSizeRange + 1)

    ! Determine diameter such that log(diameter) is uniformly distributed in the range.
    Call GetRandomNumber(U, 0)
    Diameter = (Log(UpperDiameterBound) - Log(LowerDiameterBound)) * U + Log(LowerDiameterBound)
    Diameter = Exp(Diameter)
    
    If (SizeDist%DensityPresent) Then
      Density = SizeDist%Density(jSizeRange)
    Else
      Density = Source%Density
    End If

    If (SizeDist%ParticleShapePresent) Then
      ParticleShape = SizeDist%ParticleShape(jSizeRange)
    Else
      ParticleShape = Source%ParticleShape
    End If

  Else If (Source%SourceType == S_Dust) Then

    ! Determine diameter bounds for the specified size range.
    LowerDiameterBound = SizeDist%DiameterRangeBoundary(iSizeRange)
    UpperDiameterBound = SizeDist%DiameterRangeBoundary(iSizeRange + 1)

    ! Determine particle diameter as mid-point in the range on a log scale.
    Diameter = (Log(UpperDiameterBound) - Log(LowerDiameterBound)) * 0.5 + Log(LowerDiameterBound)
    Diameter = Exp(Diameter)

    If (SizeDist%DensityPresent) Then
      Density = SizeDist%Density(iSizeRange)
    Else
      Density = Source%Density
    End If
        
    If (SizeDist%ParticleShapePresent) Then
      ParticleShape = SizeDist%ParticleShape(iSizeRange)
    Else
      ParticleShape = Source%ParticleShape
    End If

  Else If (Source%SourceType == S_SeaSalt) Then

    Call GetRandomNumber(U, 0)

    ! Smaller bin.
    If (iSizeRange == 1) Then

      Top    = 3.2
      Bottom = 0.1

      A0 = U * 2.4396E-10 - 2.5188E-16

      Do i = 1, 10

        Radius = (Top + Bottom) / 2.0
        Value  =   1.1922E-12 * Radius**6 - 1.3038E-11 * Radius**5 + 4.9943E-11 * Radius**4      &
                 - 7.5987E-11 * Radius**3 + 5.942E-11  * Radius**2 - 5.2333E-12 * Radius    - A0

        If (Value < 0.0) Bottom = Radius
        If (Value > 0.0) Top    = Radius
        If (Value == 0.0) Exit

      End Do

      Diameter = Radius * 2.0

    ! Larger bin.
    Else If (iSizeRange == 2)  Then

      Radius = 3.5

      A0 = U * 2.1935E-10 + 6.324E-10

      Do i = 1, 4
        Value      = -5.9313E-13 * Radius**3 - 1.029E-11  * Radius**2 + 2.3663E-10 * Radius - A0
        Derivative = -1.7794E-12 * Radius**2 - 2.0579E-11 * Radius    + 2.3663E-10
        Radius     = Radius - Value / Derivative
      End Do

      Diameter = Radius * 2.0

    Else

      Diameter = 0.0

    End If

    If (SizeDist%DensityPresent) Then
      Density = SizeDist%Density(iSizeRange)
    Else
      Density = Source%Density
    End If
        
    If (SizeDist%ParticleShapePresent) Then
      ParticleShape = SizeDist%ParticleShape(iSizeRange)
    Else
      ParticleShape = Source%ParticleShape
    End If

  End If

End Subroutine CalculateParticleDiameterEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine SourceTimeDependencyFactor(Time, Source, TimeDep, Factor, NextTime)
! Calculates the source time dependency factor due to a source's source time
! dependency.

  Implicit None
  ! Argument list:
  Type(Time_),      Intent(In)            :: Time
  Type(Source_),    Intent(In)            :: Source  ! $$ not used
  Type(TimeDep_),   Intent(In)            :: TimeDep
  Real(Std),        Intent(Out)           :: Factor
  Type(ShortTime_), Intent(Out), Optional :: NextTime
  ! Time     :: Time.
  ! Source   :: Source.
  ! TimeDep  :: Source time dependency.
  ! Factor   :: Source time dependency factor.
  ! NextTime :: Next time that a change in the factor might occur.
  ! Locals:
  Logical     :: InAll     ! Indicates Time lies in all time intervals associated with
                           ! a given factor in the source time dependency.
  Logical     :: In1       ! Indicates Time lies in a time interval associated with a
                           ! given factor in the source time dependency.
  Integer     :: i         !} Loop indices.
  Integer     :: j         !}
  Type(Time_) :: NextTime1 ! Next time that a change might occur in whether the time
                           ! lies in a time interval associated with a given factor in
                           ! the source time dependency.

  Factor = 1.0

  If (Present(NextTime)) Then

    NextTime = InfFutureShortTime()

    Do i = 1, TimeDep%nFactors
      InAll = .true.
      Do j = 1, TimeDep%nWildTimePairs(i)
        Call InWildTimeInterval(Time, TimeDep%WildTimePairs(j, i), In1, NextTime1)
        If (Present(NextTime)) NextTime = TMin(NextTime, Time2ShortTime(NextTime1))
        If (.not.In1) Then
          InAll = .false.
          Exit
        End If
      End Do
      If (InAll) Factor = Factor * TimeDep%Factors(i)
    End Do

  Else

    Do i = 1, TimeDep%nFactors
      InAll = .true.
      Do j = 1, TimeDep%nWildTimePairs(i)
        Call InWildTimeInterval(Time, TimeDep%WildTimePairs(j, i), In1)
        If (.not.In1) Then
          InAll = .false.
          Exit
        End If
      End Do
      If (InAll) Factor = Factor * TimeDep%Factors(i)
    End Do

  End If

End Subroutine SourceTimeDependencyFactor

!-------------------------------------------------------------------------------------------------------------

Function SourcePositionInStandardLatLong(Name, Sources, Coords)
! $$ move elsewhere or replace by more general routine?
! Returns the location of the named source in the standard lat-long coord system.

  Implicit None
  ! Argument list:
  Character(*),   Intent(In) :: Name     ! Name of the source to be located.
  Type(Sources_), Intent(In) :: Sources
  Type(Coords_),  Intent(In) :: Coords
  ! Function result:
  Real(Std) :: SourcePositionInStandardLatLong(2)  ! Long/lat of named source.
  ! Locals:
  Integer                  :: i            ! Source index.
  Integer                  :: iHCoord      ! Index of HCoord system used by source.
  Character(MaxCharLength) :: HCoordName   ! Name of HCoord system used by source.
  Type(HCoord_)            :: HCoord1      ! Horizontal coord system used by source.
  Type(HCoord_)            :: HCoord2      ! Standard lat-long horizontal coord system.
  Real(Std)                :: SourcePos(2) ! Horizontal coords of named source.

  ! Determine the source index.
  i = FindSourceIndex(Name, Sources)

  ! Determine the name and index of the horizontal coord system used by this source.
  HCoordName = Sources%Sources(i)%HCoordName
  iHCoord    = FindHCoordIndex(HCoordName, Coords)

  ! Determine the horizontal coord systems for coordinate conversion.
  HCoord1 = Coords%HCoords(iHCoord)
  HCoord2 = HCoord_LatLong()

  ! Calculate the source position.
  SourcePos = Sources%Sources(i)%X(1:2)
  SourcePositionInStandardLatLong = ConvertH(HCoord1, HCoord2, SourcePos)

End Function SourcePositionInStandardLatLong

!-------------------------------------------------------------------------------------------------------------

Subroutine ParseNParticles(CharNParticles, nParticles, Rate)
! Parses a character string to obtain information on the number of particles to be
! released.

! Character format for number of particles to be released is
!                 #Particles{/[s|sec|min|hr|day]}
! Here {...} indicates optional components, and [ | ] indicates alternative options.
! Spaces are not allowed within #Particles or within the words 'sec', 'min', 'hr' or 'day'.
! Spaces are permitted either side of / but are not required.
! 's', 'sec', 'min', 'hr' and 'day' are case insensitive.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)  :: CharNParticles
  Real(Std),    Intent(Out) :: nParticles
  Logical,      Intent(Out) :: Rate
  ! CharNParticles :: Character string to be parsed.
  ! nParticles     :: Number of particles.
  ! Rate           :: If true, indicates that nParticles is the release rate (per
  !                   second), while, if false, indicates that nParticles is the total
  !                   number released.
  ! Locals:
  Character(MaxCharLength + 1) :: CharTemp ! Temporary copy of CharNParticles for
                                           ! processing. Note the length is declared
                                           ! as MaxCharLength + 1 to ensure there is a
                                           ! trailing space.
  Integer                      :: s1       !} Positions of delimiters for bits of
  Integer                      :: s2       !} interest.
  Integer                      :: s3       !}
  Integer                      :: Error    ! Non-zero values indicates an error has
                                           ! occurred in a subroutine or function
                                           ! called from here.

  If (Len_Trim(CharNParticles) > MaxCharLength) Then
    Call Message('Error in parsing # particles: string too long', 3)
  End If

  CharTemp = CharNParticles

  s1 = Verify(CharTemp(1:), ' ')
  If (s1 == 0) Go To 1
  s1 = s1 - 1

  s2 = Scan(CharTemp(s1 + 1:), ' /')
  If (s2 == 0) Go To 1
  s2 = s2 + s1

  nParticles = Char2Std(CharTemp(s1 + 1:s2 - 1), Error)
  If (Error /= 0) Go To 1

  s3 = Verify(CharTemp(s2:), ' ')
  If (s3 == 0) Then
    Rate = .false.
  Else
    s3 = s3 + s2 - 1
    If (CharTemp(s3:s3) == '/') Then
      If ((AdjustL(CharTemp(s3 + 1:)) .CIEq. 's') .or. (AdjustL(CharTemp(s3 + 1:)) .CIEq. 'sec')) Then
        Rate = .true.
      Else If (AdjustL(CharTemp(s3 + 1:)) .CIEq. 'min') Then
        Rate       = .true.
        nParticles = nParticles / 60.0
      Else If (AdjustL(CharTemp(s3 + 1:)) .CIEq. 'hr') Then
        Rate       = .true.
        nParticles = nParticles / 3600.0
      Else If (AdjustL(CharTemp(s3 + 1:)) .CIEq. 'day') Then
        Rate       = .true.
        nParticles = nParticles / 86400.0
      Else
        Go To 1
      End If
    Else
      Go To 1
    End If
  End If

  Return

1 Continue
  Call Message('Error in parsing # particles', 3)

End Subroutine ParseNParticles

!-------------------------------------------------------------------------------------------------------------

Subroutine ParseSourceStrength(CharSourceStrength, SpeciesName, Q, MaterialUnitName, Rate, QuestionMark)
! Parses a character string to obtain information on the source strength for a given
! species.

! Character format for source strength is
!                 SpeciesName SourceStrength MaterialUnit{/[s|sec|min|hr|day]}
! Here {...} indicates optional components, and [ | ] indicates alternative options.
! Spaces are not generally allowed within SpeciesName, SourceStrength or MaterialUnit or within the words
! 'sec', 'min', 'hr' or 'day'; however they are permitted in SpeciesName if it is enclosed in quotes. Quotes
! can be '...' or "...", but quotes of the same type mustn't be used within SpeciesName.
! Spaces between SpeciesName, SourceStrength and MaterialUnit are compulsory.
! Spaces are permitted either side of / but are not required.
! SpeciesName, MaterialUnit, 's', 'sec', 'min', 'hr' and 'day' are case insensitive.
! When the source strength is met dependent, a question mark should be used in place of SourceStrength. In
! this case the source strength is set to 1.

  Implicit None
  ! Argument list:
  Character(*),             Intent(In)  :: CharSourceStrength
  Character(MaxCharLength), Intent(Out) :: SpeciesName
  Real(Std),                Intent(Out) :: Q
  Character(MaxCharLength), Intent(Out) :: MaterialUnitName
  Logical,                  Intent(Out) :: Rate
  Logical,                  Intent(Out) :: QuestionMark
  ! CharSourceStrength :: Character string to be parsed.
  ! SpeciesName        :: Species name.
  ! Q                  :: Source strength for the species.
  ! MaterialUnitName   :: Name of Material Unit used for the species.
  ! Rate               :: If true, indicates that Q is the release rate (per second),
  !                       while, if false, indicates that Q is the total amount
  !                       released.
  ! QuestionMark       :: Indicates a question mark has been used for the source strength in
  !                       CharSourceStrength.
  ! Locals:
  Character(MaxCharLength + 1) :: CharTemp ! Temporary copy of CharSourceStrength for
                                           ! processing. Note the length is declared
                                           ! as MaxCharLength + 1 to ensure there is a
                                           ! trailing space.
  Integer                      :: s1       !} Positions of delimiters for bits of
  Integer                      :: s2       !} interest.
  Integer                      :: s3       !}
  Integer                      :: s4       !}
  Integer                      :: s5       !}
  Integer                      :: s6       !}
  Integer                      :: s7       !}
  Integer                      :: Error    ! Non-zero values indicates an error has
                                           ! occurred in a subroutine or function
                                           ! called from here.

  If (Len_Trim(CharSourceStrength) > MaxCharLength) Then
    Call Message('Error in parsing source strength: string too long', 3)
  End If

  CharTemp = CharSourceStrength

  s1 = Verify(CharTemp(1:), ' ''"')
  If (s1 == 0) Go To 1
  s1 = s1 - 1

  If (s1 == 0) Then
    s2 = Scan(CharTemp(s1 + 1:), ' ')
  Else
    s2 = Scan(CharTemp(s1 + 1:), CharTemp(s1:s1))
  End If
  If (s2 == 0) Go To 1
  s2 = s2 + s1

  SpeciesName = CharTemp(s1 + 1:s2 - 1)

  s3 = Verify(CharTemp(s2 + 1:), ' ')
  If (s3 == 0) Go To 1
  s3 = s3 + s2 - 1

  s4 = Scan(CharTemp(s3 + 1:), ' ')
  If (s4 == 0) Go To 1
  s4 = s4 + s3

  If (CharTemp(s3 + 1:s4 - 1) == '?') Then
    QuestionMark = .true.
    Q            = 1.0
  Else
    QuestionMark = .false.
    Q            = Char2Std(CharTemp(s3 + 1:s4 - 1), Error)
    If (Error /= 0) Go To 1
  End If

  s5 = Verify(CharTemp(s4 + 1:), ' ')
  If (s5 == 0) Go To 1
  s5 = s5 + s4 - 1

  s6 = Scan(CharTemp(s5 + 1:), ' /')
  If (s6 == 0) Go To 1
  s6 = s6 + s5

  MaterialUnitName = CharTemp(s5 + 1:s6 - 1)

  s7 = Verify(CharTemp(s6:), ' ')
  If (s7 == 0) Then
    Rate = .false.
  Else
    s7 = s7 + s6 - 1
    If (CharTemp(s7:s7) == '/') Then
      If ((AdjustL(CharTemp(s7 + 1:)) .CIEq. 's') .or. (AdjustL(CharTemp(s7 + 1:)) .CIEq. 'sec')) Then
        Rate = .true.
      Else If (AdjustL(CharTemp(s7 + 1:)) .CIEq. 'min') Then
        Rate = .true.
        Q    = Q / 60.0
      Else If (AdjustL(CharTemp(s7 + 1:)) .CIEq. 'hr') Then
        Rate = .true.
        Q    = Q / 3600.0
      Else If (AdjustL(CharTemp(s7 + 1:)) .CIEq. 'day') Then
        Rate = .true.
        Q    = Q / 86400.0
      Else
        Go To 1
      End If
    Else
      Go To 1
    End If
  End If

  If (QuestionMark) Q = 1.0

  Return

1 Continue
  Call Message('Error in parsing source strength', 3)

End Subroutine ParseSourceStrength

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcDustSourceStrength(       &
             Flow, Surface, Soil, Area,  &
             Diameter, iSizeRange,       &
             SourceStrength              &
           )
! Calculates an emission rate for PM10 resulting from dust.

  Implicit None
  ! Argument list:
  Type(Flow_),    Intent(In)  :: Flow
  Type(Surface_), Intent(In)  :: Surface
  Type(Soil_),    Intent(In)  :: Soil
  Real(Std),      Intent(In)  :: Area
  Real(Std),      Intent(In)  :: Diameter
  Integer,        Intent(In)  :: iSizeRange
  Real(Std),      Intent(Out) :: SourceStrength
  ! Locals:
  Real(8) :: TFV ! Threshold fraction velocity.
  Real(8) :: H   ! Horizontal flux.

  SourceStrength = 0.0

  ! Test for significant amount of bare soil. Note all elements of Surface are negative over the sea - we test
  ! for not being over the sea in ClayFrac and SoilMoisture in case the sea boundary is slightly different
  ! from that in LandUseFracs(8).
  ! $$ SoilFracs not tested. Perhaps better to make consistent in AncillaryMet/NWPFlow.
  If (Surface%LandUseFracs(8) > 0.3 .and. Soil%ClayFrac >= 0.0 .and. Surface%SoilMoisture >= 0.0) Then

    If ((Soil%SoilFracs(iSizeRange) > 0.0) .and. (Flow%UStar >= 1.0E-5)) Then

      ! Calculate the threshold flux velocity.
      TFV = -0.2 * Log10(Diameter/1.0E6) + 0.05 * Surface%SoilMoisture - 1.2

      ! If TFV less than 0 then set to 0.
      If (TFV < 0.0) Then
        TFV = 0.0
      End If

      If (TFV >= Flow%UStar) Then
        H = 0.0
      Else

        H = 2.61 * Flow%Rho * Surface%LandUseFracs(8) * Flow%UStar * Flow%UStar * Flow%UStar * &
            (1.0 + TFV / Flow%UStar) *                                                       &
            (1.0 - ((TFV / Flow%UStar)*(TFV / Flow%UStar))) *                                &
            Soil%SoilFracs(iSizeRange) / Gravity

      End If

      ! This is in terms of g emitted over the whole area per sec.
      SourceStrength = 1000.0 * Area * H * 10.0**(13.4 * Min(Soil%ClayFrac, 0.2) - 6.0)

    End If

  End If

End Subroutine CalcDustSourceStrength

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcSeaSaltSourceStrength( &
             Flow, Surface, Area,     &
             Diameter, iSizeRange,    &
             SourceStrength           &
           )
! Calculates an emission rate for PM10 resulting from sea salt.

  Implicit None
  ! Argument list:
  Type(Flow_),    Intent(In)  :: Flow
  Type(Surface_), Intent(In)  :: Surface
  Real(Std),      Intent(In)  :: Area
  Real(Std),      Intent(In)  :: Diameter
  Integer,        Intent(In)  :: iSizeRange
  Real(Std),      Intent(Out) :: SourceStrength
  ! Locals:
  Real(8) :: Speed                   ! Wind speed at 10m.
  Real(8) :: Calc_Emission_Rate_Salt ! These have been hard-wired as the bin sizes have been decided upon.

  SourceStrength = 0.0

  ! If over sea, do sea-salt calculation.
  ! Surface%LandFrac is the land fraction (0 = totally sea, 1 = totally land)
  If (Surface%LandFrac < 1.0) Then

    ! 10m wind speed is the root of the U squared + V squared.
    Speed = Sqrt(Flow%U(1) * Flow%U(1) + Flow%U(2) * Flow%U(2))

    ! Multiplying by the area to calculate the total emission over the grid box
    ! gives an emission rate in grams / second.

    If (iSizeRange == 1) Then

      Calc_Emission_Rate_Salt = 2.4396E-10

    Else If (iSizeRange == 2) Then

      Calc_Emission_Rate_Salt = 2.1935E-10

    Else

      Call Message('FATAL ERROR: with sea-salt routine', 3)

    End If

    SourceStrength = Area * (1.0 - Surface%LandFrac) * (Speed**3.41) *    &
                       Calc_Emission_Rate_Salt / 9.0

  End If

End Subroutine CalcSeaSaltSourceStrength

!-------------------------------------------------------------------------------------------------------------

End Module SourceModule

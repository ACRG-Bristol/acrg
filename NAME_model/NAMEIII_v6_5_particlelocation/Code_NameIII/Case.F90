! Module:  Case Module

Module CaseModule

! This module provides code to compute one case.

!-------------------------------------------------------------------------------------------------------------

Use TimerModule
Use ServiceModule
Use FlowAndFlowProfileModule, Only: ModifyTurbByTerminalVelocity
Use FlowsModule, Only: Flows_, FlowMemory_,                                       &
                       MetsFlowsOverallTValid, UpdateMetsFlows,                   &
                       ConvertToZ, GetAttrib, ResetFlowMemory, Reflect, CalcdZdZ, &
                       XInTheDomain,                                              &
                       A_Flow, A_Cloud, A_Rain, A_Surface, A_Soil, A_Plant,       &
                       Mets_,                                                     &
                       Flow_, Cloud_, Rain_, Surface_, Soil_, Plant_,             &
                       ConvertFlow, ReverseFlow
Use SizeDistModule
Use SpeciesModule
Use SourceModule
Use OutputModule, Only: OutputOpts_, FieldReq_, PPInfoReq_, Reqs_, Field_, Results_, ReqInfo_, CalciTA,     &
                        Q_AirConc, Q_DryDep, Q_WetDep, Q_PhotonFlux,                                        &
                        ! Q_Dep,                                                                            &
                        Q_nParticles, Q_nParticlesBySpec, Q_nPuffs, Q_nParticleSteps, Q_nPuffSteps, Q_Mass, &
                        Q_MeanZ, Q_SigmaZ, Q_XStats, Q_MeanS, Q_PuffCentres,                                &
                        !Q_SigmaC,                                                                          &
                        Q_ChemistryField,                                                                   &
                        Q_SigmaWW, Q_TauWW, Q_MeanFlowU, Q_MeanFlowV, Q_MeanFlowW, Q_TemperatureK,          &
                        Q_PotentialTempK, Q_SpecificHumidity, Q_PressurePa, Q_Density, Q_Topography,        &
                        Q_UStar, Q_HeatFlux, Q_BLDepth, Q_WindSpeed, Q_WindDirectionDeg, Q_PptRateMMHR,     &
                        Q_TemperatureC, Q_CloudOktas, Q_RHPercent, Q_Pasquill, Q_ProgressPercent,           &
                        Q_ClockTime, Q_X, Q_Y, Q_SigmaVV, Q_MesoscaleSigmaVV,                               &
                        Q_TotalOrDynCloudWater, Q_TotalOrDynCloudIce, Q_Cloud3d, Q_RoughnessLength,         &
                        Q_PSeaLevelPa, Q_LandUseFracs, Q_CanopyWater, Q_LAI, Q_CanopyHeight,                &
                        Q_StomataConduct, Q_SoilMoisture, Q_LandFrac, Q_EulerianConcentration,              & 
                        Q_OriginalSourceStrength, Q_RevisedSourceStrength,                                  &
                        Q_ConCloudBase, Q_ConCloudTop,                                                      &
                        Q_EulerianTotalDep, Q_EulerianDryDep, Q_EulerianWetDep,                             &
                        FirstLastReqTimes, PrepareResults, UpdateResults,                                   &
                        CalciLast,                                                                          &
                        CalcReqInfoNoT, CalcReqInfoT, CalcReqInfoS, CalcReqInfoEveryT,                      &
                        ProcessAndOutputResults
Use ParticleModule
Use PuffModule
Use ChemistryModule
Use OpenMPModule, Only : OpenMPOpts_
Use EulerianInterfaceModule
Use SemiLagrangianModule
Use IterativePlumeRiseModel

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: Time3D_                   ! 3D array of times. $$ needs to be public on sun workstations, probably
                                     ! because components of DispState are accessible in MainNameIII.F90.
                                     ! Would perhaps be better if design allows them to be private.
Public  :: MainOpts_                 ! Main options.
Public  :: MultiCaseOpts_            ! Multiple case options.
Public  :: DispOpts_                 ! A set of dispersion model options.
Public  :: DispOptses_               ! A collection of sets of dispersion model options.
Public  :: DispState_                ! Information on the state of the dispersion calculation.
Public  :: PreInitMainOpts           ! Pre-initialises the main options.
Public  :: InitMainOpts              ! Initialises the main options.
Public  :: CheckMainOpts             ! Checks the main options.
Public  :: PreInitMultiCaseOpts      ! Pre-initialises the multiple case options.
Public  :: InitMultiCaseOpts         ! Initialises the multiple case options.
Public  :: CheckMultiCaseOpts        ! Checks the multiple case options.
Public  :: InitDispOptses            ! Initialises a collection of sets of dispersion model options.
Public  :: InitDispOpts              ! Initialises a set of dispersion model options.
Public  :: AddDispOpts               ! Add a set of dispersion model options to a collection of sets of
                                     ! dispersion model options.
Public  :: WriteRestartFileDispState
Public  :: ReadRestartFileDispState
Public  :: InitDispState
Public  :: RestartAdjustmentForInputChanges ! Adjusts the dispersion model state to account for
                                            ! permitted changes to the input file on restart.
Public  :: PrepareForNextCase        ! Prepares the dispersion model state for running the next case.
Public  :: StartCase                 !
Public  :: RunToRestartDumpOrSuspendOrEndOfCase ! This routine carries out the loop over successive
                                     ! updates of the met and flow module instances,
                                     ! calling LoopOverSyncTimes once between each update
                                     ! time.
Public  :: CaseModuleTimerSummary
Public  :: CaseModuleTimerInitialise

!-------------------------------------------------------------------------------------------------------------

Type :: MainOpts_ ! Main options.
  Logical                      :: Initialised
  Character(MaxCharLength)     :: RunName
  Character(MaxCharLength)     :: TimeFrameName
  Logical                      :: FlatEarth
  Logical                      :: Backwards
  Logical                      :: FixedMet
  Character(1)                 :: RandomMode
  Character(MaxFileNameLength) :: RunToFile
  Integer                      :: MaxSources
  Logical                      :: SameResultsWithUpdateOnDemand
  Integer                      :: MaxFieldReqs
  Integer                      :: MaxFieldOutputGroups
  
  ! Initialised     :: Indicates the main options have been initialised.
  ! RunName         :: Name of run.
  ! TimeFrameName   :: Name of time frame.
  ! FlatEarth       :: Indicates a flat Earth calculation.
  ! Backwards       ::
  ! FixedMet        :: Indicates a calculation with fixed met/flow.
  ! RandomMode      :: Indicates initial random number seed is fixed/random/given
  ! RunToFile       :: File containing time to suspend run.
  ! MaxSources      :: Maximum number of sources.
  ! SameResultsWithUpdateOnDemand :: Indicates results should be made the same whether the met and flow module
  !                                  instances use update-on-demand or update-at-once (so it should really be
  !                                  called SameResultsWithOrWithoutUpdateOnDemand).
  ! MaxFieldReqs                  :: Maximum number of Field Output Requirements.
  ! MaxFieldOutputGroups          :: Maximum number of Field Output Groups.
  
End Type MainOpts_

!-------------------------------------------------------------------------------------------------------------

Type :: MultiCaseOpts_ ! Multiple case options.
  Logical :: Initialised            ! Indicates the multiple case options have been initialised.
  Integer :: DispOptsesEnsembleSize ! Number of sets of dispersion options to use.
  Integer :: MetEnsembleSize        ! Size of the met ensemble (i.e. number of met realisations).
End Type MultiCaseOpts_

!-------------------------------------------------------------------------------------------------------------

Type :: DispOpts_ ! A set of dispersion model options.
  Type(Time_)              :: FixedMetTime
  Integer                  :: MaxParticles
  Integer                  :: MaxFullParticles
  Integer                  :: MaxPuffs
  Integer                  :: MaxOriginalPuffs
  Integer                  :: ParticleCeiling
  Real(Std)                :: ParticleFactor
  Real(Std)                :: Zs
  Type(Time_)              :: SkewTime
  Type(Time_)              :: VelMemTime
  Type(Time_)              :: InhomogTime
  Type(Time_)              :: MVelMemTime
  Type(Time_)              :: PuffTime
  Type(ShortTime_)         :: sSkewTime
  Type(ShortTime_)         :: sVelMemTime
  Type(ShortTime_)         :: sInhomogTime
  Type(ShortTime_)         :: sMVelMemTime
  Type(ShortTime_)         :: sPuffTime
  Type(Time_)              :: SyncdT
  Character(MaxCharLength) :: DomainName
  Type(Time_)              :: PuffInterval
  Integer                  :: DeltaOpt
  Integer                  :: DeepConvectionCode
  Logical                  :: RadioactiveDecay
  Logical                  :: AgentDecay
  Logical                  :: DryDep
  Logical                  :: WetDep
  Logical                  :: Turbulence
  Logical                  :: MesoscaleMotions
  Logical                  :: VerticalVelocity
  Logical                  :: Damping
  Logical                  :: Chemistry
  Type(PuffOpts_)          :: PuffOpts
  ! FixedMetTime     :: Time of fixed met data (for FixedMet = .true. only).
  ! MaxParticles     :: Maximum number of particles.
  ! MaxFullParticles :: Maximum number of full particles (i.e. particles for which a
  !                     set of extra information is provided).
  ! MaxPuffs         :: Maximum number of puffs.
  ! MaxOriginalPuffs :: Maximum number of original puffs (i.e. puffs released at the
  !                     source and not created by puff splitting).
  ! Zs               :: Height below which particles are subject to dry deposition
  ! SkewTime         :] Travel time over which the model allows for (i) skewness, (ii)
  ! VelMemTime       :] velocity memory, and (iii) inhomogeneity. These travel times
  ! InhomogTime      :] must satisfy SkewTime <= VelMemTime <= InhomogTime.
  ! MVelMemTime      :: unresolved mesoscale motions vel mem time
  ! PuffTime         :: Travel time over which puffs are used.
  ! sSkewTime        :} Short time versions of SkewTime, VelMemTime, InhomogTime and
  ! sVelMemTime      :} PuffTime. These are meaningful only if the time is not
  ! sInhomogTime     :} infinite.
  ! sPuffTime        :}
  ! SyncdT           :: Time interval at which particles and puffs are synchronised.
  ! DomainName       :: Name of computational domain.
  ! PuffInterval     :: Time interval between release of puffs (for non-fixed met
  !                     cases).
  ! DeltaOpt         :: 0 indicates Delta = infinity; 1 indicates Delta limited by
  !                  :: inhomogeneity; 2 indicates Delta limited by inhomogeneity and
  !                  :: by sigma_p.
  ! DeepConvectionCode :: Code for deep convection scheme used. 
  ! RadioactiveDecay :} Indicates that radioactive decay,
  ! AgentDecay       :} agent decay, dry deposition, wet deposition, turbulence
  ! DryDep           :} and unresolved mesoscale motions are modelled.
  ! WetDep           :}
  ! Turbulence       :}
  ! MesoscaleMotions :}
  ! VerticalVelocity :: Indicates that particles are advected in 3-dimensions.
  !                     If false, 2-d trajectories are generated. Defaults to true.
  ! Damping - near source damping of eddy diffusive growth.
  ! Chemistry        :: Indicates that chemistry scheme is run.
End Type DispOpts_

!-------------------------------------------------------------------------------------------------------------

Type :: DispOptses_ ! A collection of sets of dispersion model options.
  Integer         :: nDispOptses               ! Number of sets of dispersion model options.
  Type(DispOpts_) :: DispOptses(MaxDispOptses) ! Sets of dispersion model options.
End Type DispOptses_

!-------------------------------------------------------------------------------------------------------------

Type :: Time3D_ ! 3D array of times.
  Type(ShortTime_), Pointer :: T(:, :, :)
End Type Time3D_

!-------------------------------------------------------------------------------------------------------------

Type :: DispState_ ! A set of information on the state of the dispersion calculation.
  ! If changing this type, remember restart read and write routines.
  Logical          :: LastCase
  Integer          :: iDispOpts
  Integer          :: iMetCase
  Integer          :: iDomain
  Type(ShortTime_) :: sLastDispTime
  Type(ShortTime_) :: sFirstReqTime
  Type(ShortTime_) :: sLastReqTime
  Type(ShortTime_) :: sFirstReleaseTime
  Type(ShortTime_) :: sLastReleaseTime ! Latest possible release time,
                                       ! accounting for sources, time part of comp
                                       ! domain and fixed met flag.
  Type(ShortTime_) :: sFirstMetTime !} First and last times when met might be needed.
  Type(ShortTime_) :: sLastMetTime  !}
  Type(Time_)      :: Time
  Type(Time_)      :: SyncTime
  Type(Time_)      :: MetTime
  Logical          :: PPsFinished

  Logical          :: SpaceAllocated
  Integer          :: MaxParticles
  Integer          :: MaxFullParticles
  Integer          :: MaxPuffs
  Integer          :: MaxOriginalPuffs
  Integer          :: nParticleSpecieses

  Integer                   :: nParticles
  Integer                   :: LastParticle
  Integer,          Pointer :: FreeParticleStack(:)
  Integer                   :: nParticleExtras
  Integer                   :: LastParticleExtra
  Integer,          Pointer :: FreeParticleExtraStack(:)
  Integer,          Pointer :: iParticleExtras(:) ! Indices of the Extra corresponding to a given Particle
  Type(Particle_),  Pointer :: Particles(:)
  Real(Std),        Pointer :: ParticleMasses(:, :) ! 1st index species, 2nd index particles/puffs
  Type(Extra_),     Pointer :: ParticleExtras(:)    ! Extra(0) is a default extra for
                                                    ! cheap particles, rest are
                                                    ! useable for expensive particles.

  Integer                   :: nOriginalPuffs
  Integer                   :: LastOriginalPuff
  Integer,          Pointer :: FreeOriginalPuffStack(:)
  Integer                   :: nPuffs
  Integer                   :: LastPuff
  Integer,          Pointer :: FreePuffStack(:)
  Integer,          Pointer :: iPuffs(:) ! Indices of the Puff corresponding to a given original puff
  Type(Puff_),      Pointer :: Puffs(:)
  Real(Std),        Pointer :: PuffMasses(:, :)
  Type(Extra_),     Pointer :: PuffExtras(:)

  Type(ShortTime_), Pointer :: TLast(:)
  Type(ShortTime_), Pointer :: T(:)
  Real(Std),        Pointer :: SigA2Last(:)
  Real(Std),        Pointer :: SigA2(:)
  Integer,          Pointer :: SigA2Defined(:)

  Integer                   :: iULastParticle     ! Unique index of last particle introduced. $$ Integer(8)
  Integer                   :: iULastPuff         ! Unique index of last puff introduced. $$ Integer(8)
  Integer                   :: iULastOriginalPuff ! Unique index of last original puff introduced.
                                                  ! $$ Integer(8)

  Integer          :: ParticleFactorType
  Integer(I64)     :: nParticleTimeSteps
  Integer(I64)     :: nPuffTimeSteps
  Type(ShortTime_) :: MaxTSpread
  Type(Time3D_),    Pointer :: NextReleaseTime(:)
  Type(ShortTime_), Pointer :: NextReleaseTimeM(:)
  Integer          :: iZCoordMagl
  Integer          :: iHCoordLatLong

  Logical          :: PPsLostDueToFlow
  Type(EulerianField_) :: EulerianField

  ! LastCase           :: Indicates this is the last case.
  ! iDispOpts          :: Index of the set of dispersion model options being used.
  ! iMetCase           :: Number of the met realisation in the met ensemble.
  ! iDomain            :: Index of computational domain.
  ! Time               :: Current time.
  ! SyncTime           :: Time particles and puffs were last synchronised.
  ! MetTime            :: Time of met.

  ! ParticleMasses     :: Mass of each species carried on particle. The elements in this mass array are
  !                       mapped to the collection of species in Specieses using Specieses%iParticle2Species
  !                       and Specieses%iSpecies2Particle.

  ! LastParticle       :} Index of last particle and puff in the collections of
  ! LastPuff           :} particles and puffs.
  ! nParticles         :: Number of particles which are currently active.
  ! FreeParticleStack  :: Indices of free particles available for release.
  ! nOriginalPuffs     :: Number of original puffs (i.e. puffs released at the source
  !                       and not created by puff splitting).
  ! Particles          :] Collections of particles and puffs.
  ! Puffs              :]
  ! ParticleFactorType :: 0 - normal
  !                       1 - particle release rate damped
  !                       2 - no particles being released - see source.f90 for details
  ! nParticleTimeSteps :} Number of particle and puff time steps made.
  ! nPuffTimeSteps     :}
  ! MaxTSpread         :: Maximum time spread of particle/puffs about their nominal
  !                       time.
  ! TLast              :: !
  ! T                  :: ! $$ comments needed here $$
  ! SigA2Last          :: ! $$ also needs reviewing/moving to puff?
  ! SigA2              :: !
  ! SigA2Defined       :: 0 = no values defined, 1 = SigA2 and T defined, 2 = SigA2,
  !                    :: T, SigA2Last and TLast defined.
  ! NextReleaseTime
  ! NextReleaseTimeM   :: T- time for puffs

  ! PPsLostDueToFlow   :: Particles/puffs lost due to absence of convert/flow/cloud etc info
End Type DispState_

! Used for checking whether Timers need to be initialised.
Logical :: TimersInitialised = .false.
! Timers (see Timer module)
Type(Timer_), Save :: LPPATsTimer,    &       ! Main particle loop
                      CPRSTimer,      &       ! Particle update loop
                      OutputTimer             ! ProcessAndOutputresults

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

! Initialise timers

Subroutine CaseModuleTimerInitialise(TimerOpts)

  Implicit None

  Type(TimerOpts_) :: TimerOpts

  ! Create Timers
  If (.not. TimersInitialised) Then
    Call TimerCreate(LPPATsTimer,   "WorkerThread_LPPATs"                  ,TimerOpts)
    Call TimerCreate(CPRSTimer,     "WorkerThread_CPRS"                    ,TimerOpts)
    Call TimerCreate(OutputTimer,   "WorkerThread_ProcessAndOutputResults" ,TimerOpts)
    TimersInitialised = .true.
  End If

End Subroutine CaseModuleTimerInitialise

!-------------------------------------------------------------------------------------------------------------

! Output timer summary information

Subroutine CaseModuleTimerSummary()

  Implicit None

  If (TimersInitialised) Then
    Call TimerWriteSummary(LPPATsTimer)
    Call TimerWriteSummary(CPRSTimer)
    Call TimerWriteSummary(OutputTimer)
  End If

End Subroutine CaseModuleTimerSummary

!-------------------------------------------------------------------------------------------------------------

Function PreInitMainOpts()
! Pre-initialises the main options.

  Implicit None
  ! Function result:
  Type(MainOpts_) :: PreInitMainOpts ! Main options.
  ! Locals:
  Type(MainOpts_) :: MainOpts ! Local copy of function result.

  MainOpts%Initialised = .false.

  PreInitMainOpts = MainOpts

End Function PreInitMainOpts

!-------------------------------------------------------------------------------------------------------------

Subroutine InitMainOpts(                    &
             RunName,                       &
             TimeFrameName,                 &
             FlatEarth,                     &
             Backwards,                     &
             FixedMet,                      &
             RandomMode,                    &
             RunToFile,                     &
             MaxSources,                    &
             SameResultsWithUpdateOnDemand, &
             MaxFieldReqs,                  &
             MaxFieldOutputGroups,          &
             MainOpts                       &
           )
! Initialises the main options.

  Implicit None
  ! Argument list:
  Character(*),    Intent(In)    :: RunName
  Character(*),    Intent(In)    :: TimeFrameName
  Logical,         Intent(In)    :: FlatEarth
  Logical,         Intent(In)    :: Backwards
  Logical,         Intent(In)    :: FixedMet
  Character(1),    Intent(In)    :: RandomMode
  Character(*),    Intent(In)    :: RunToFile
  Integer,         Intent(In)    :: MaxSources
  Logical,         Intent(In)    :: SameResultsWithUpdateOnDemand
  Integer,         Intent(In)    :: MaxFieldReqs
  Integer,         Intent(In)    :: MaxFieldOutputGroups
  Type(MainOpts_), Intent(InOut) :: MainOpts
  ! RunName         ::
  ! TimeFrameName   :: Name of time frame.
  ! FlatEarth       :: Indicates a flat Earth calculation.
  ! Backwards       ::
  ! FixedMet        :: Indicates a calculation with fixed met/flow.
  ! RandomMode      ::
  ! MaxSources      :: Maximum number of sources.
  ! MainOpts        :: Main options.
  ! SameResultsWithUpdateOnDemand :: Indicates results should be made the same whether the met and flow module
  !                                  instances use update-on-demand or update-at-once (so it should really be
  !                                  called SameResultsWithOrWithoutUpdateOnDemand).
  ! MaxFieldReqs                  :: Maximum number of field requirements.
  ! MaxFieldOutputGroups          :: Maximum number of field output groups.

  If (MainOpts%Initialised) Then
    Call Message(                                                          &
           'FATAL ERROR in InitMainOpts: an attempt has been made to '  // &
           're-initialise the main options - check that there is only ' // &
           'one set of main options specified in the input',               &
           3                                                               &
         )
  End If

  Call TokenLengthTest(TimeFrameName, MaxCharLength, .true., 'Main Options', ' ', &
                       'Absolute or Relative Time?')

  MainOpts%TimeFrameName   = TimeFrameName
  MainOpts%FixedMet        = FixedMet
  MainOpts%FlatEarth       = FlatEarth
  MainOpts%Backwards       = Backwards
  MainOpts%RunName         = RunName     ! $$ test length

  MainOpts%RandomMode      = RandomMode

  MainOpts%RunToFile       = RunToFile

  MainOpts%MaxSources      = MaxSources

  MainOpts%SameResultsWithUpdateOnDemand = SameResultsWithUpdateOnDemand

  MainOpts%MaxFieldReqs                  = MaxFieldReqs
  
  MainOpts%MaxFieldOutputGroups          = MaxFieldOutputGroups
  
  MainOpts%Initialised                   = .true.

End Subroutine InitMainOpts

!-------------------------------------------------------------------------------------------------------------

Subroutine CheckMainOpts(MainOpts)
! Checks the main options.

  Implicit None
  ! Argument list:
  Type(MainOpts_), Intent(In) :: MainOpts ! Main options.

  If (.not.MainOpts%Initialised) Then
    Call Message('FATAL ERROR in CheckMainOpts', 3)
  End If

End Subroutine CheckMainOpts

!-------------------------------------------------------------------------------------------------------------

Function PreInitMultiCaseOpts()
! Pre-initialises the multiple case options.

  Implicit None
  ! Function result:
  Type(MultiCaseOpts_) :: PreInitMultiCaseOpts ! Multiple case options.
  ! Locals:
  Type(MultiCaseOpts_) :: MultiCaseOpts ! Local copy of function result.

  MultiCaseOpts%Initialised = .false.

  PreInitMultiCaseOpts = MultiCaseOpts

End Function PreInitMultiCaseOpts

!-------------------------------------------------------------------------------------------------------------

Subroutine InitMultiCaseOpts(DispOptsesEnsembleSize, MetEnsembleSize, MultiCaseOpts)
! Initialises the multiple case options.

  Implicit None
  ! Argument list:
  Integer,              Intent(In)    :: DispOptsesEnsembleSize
  Integer,              Intent(In)    :: MetEnsembleSize
  Type(MultiCaseOpts_), Intent(InOut) :: MultiCaseOpts
  ! DispOptsesEnsembleSize :: Number of dispersion options scenarios to be considered.
  ! MetEnsembleSize        :: Number of met scenarios to be considered.
  ! MultiCaseOpts          :: Multiple case options.

  If (MultiCaseOpts%Initialised) Then
    Call Message(                                                                            &
           'FATAL ERROR in InitMultiCaseOpts: an attempt has been made to re-initialise ' // &
           'the multiple case options - check that there is only one set of '             // &
           'multiple case options specified in the input',                                   &
           3                                                                                 &
         )
  End If

  If (DispOptsesEnsembleSize <= 0) Then
    Call Message('ERROR in reading multiple case options: Dispersion Options Ensemble Size should be >= 1', 3)
  End If
  MultiCaseOpts%DispOptsesEnsembleSize = DispOptsesEnsembleSize
  If (MetEnsembleSize <= 0) Then
    Call Message('ERROR in reading multiple case options: Met Ensemble Size should be >= 1', 3)
  End If
  MultiCaseOpts%MetEnsembleSize = MetEnsembleSize

  MultiCaseOpts%Initialised = .true.

End Subroutine InitMultiCaseOpts

!-------------------------------------------------------------------------------------------------------------

Subroutine CheckMultiCaseOpts(MultiCaseOpts, DispOptses)
! Checks the multiple case options.

  Implicit None
  ! Argument list:
  Type(MultiCaseOpts_), Intent(In) :: MultiCaseOpts ! Multiple case options.
  Type(DispOptses_),    Intent(In) :: DispOptses    ! Collection of sets of dispersion model options.

  If (.not.MultiCaseOpts%Initialised) Then
    Call Message('FATAL ERROR in CheckMultiCaseOpts', 3)
  End If

  ! Note met ensemble size checked in individual met modules.

  If (DispOptses%nDispOptses < MultiCaseOpts%DispOptsesEnsembleSize) Then
    Call Message(                                                                                &
           'FATAL ERROR in checking multiple case options: the "Dispersion Options Ensemble ' // &
           'Size" is larger than the number of sets of dispersion options',                      &
           3                                                                                     &
         )
  End If

End Subroutine CheckMultiCaseOpts

!-------------------------------------------------------------------------------------------------------------

Function InitDispOptses()
! Initialises a collection of sets of dispersion model options.

  Implicit None
  ! Function result:
  Type(DispOptses_) :: InitDispOptses ! Collection of sets of dispersion model options.
  ! Locals:
  Type(DispOptses_) :: DispOptses ! Local copy of function result.

  DispOptses%nDispOptses = 0

  InitDispOptses = DispOptses

End Function InitDispOptses

!-------------------------------------------------------------------------------------------------------------

Subroutine InitDispOpts(                                                     &
             FixedMetTime,                                                   &
             MaxParticles, MaxFullParticles, MaxPuffs, MaxOriginalPuffs,     &
             ParticleCeiling, ParticleFactor,                                &
             SkewTime, VelMemTime, InhomogTime, MVelMemTime, PuffTime,       &
             SyncdT, DomainName,                                             &
             PuffInterval, DeltaOpt,                                         &
             DeepConvection,                                                 &
             RadioactiveDecay, AgentDecay,                                   &
             DryDep, WetDep, Turbulence, MesoscaleMotions, Damping,          &
             VerticalVelocity, Chemistry,                                    &
             A1, A5, A7, Zs,                                                 &
             DispOpts                                                        &
           )
! Initialises a set of dispersion model options.

  Implicit None
  ! Argument list:
  Type(Time_),     Intent(In)    :: FixedMetTime
  Integer,         Intent(In)    :: MaxParticles
  Integer,         Intent(In)    :: MaxFullParticles
  Integer,         Intent(In)    :: MaxPuffs
  Integer,         Intent(In)    :: MaxOriginalPuffs
  Integer,         Intent(In)    :: ParticleCeiling
  Real(Std),       Intent(In)    :: ParticleFactor
  Type(Time_),     Intent(In)    :: SkewTime
  Type(Time_),     Intent(In)    :: VelMemTime
  Type(Time_),     Intent(In)    :: InhomogTime
  Type(Time_),     Intent(In)    :: MVelMemTime
  Type(Time_),     Intent(In)    :: PuffTime
  Type(Time_),     Intent(In)    :: SyncdT
  Character(*),    Intent(In)    :: DomainName
  Type(Time_),     Intent(In)    :: PuffInterval
  Integer,         Intent(In)    :: DeltaOpt
  Character(*),    Intent(In)    :: DeepConvection
  Logical,         Intent(In)    :: RadioactiveDecay
  Logical,         Intent(In)    :: AgentDecay
  Logical,         Intent(In)    :: DryDep
  Logical,         Intent(In)    :: WetDep
  Logical,         Intent(In)    :: Turbulence
  Logical,         Intent(In)    :: MesoscaleMotions
  Logical,         Intent(In)    :: VerticalVelocity
  Logical,         Intent(In)    :: Damping
  Logical,         Intent(In)    :: Chemistry
  Real(Std),       Intent(In)    :: A1
  Real(Std),       Intent(In)    :: A5
  Real(Std),       Intent(In)    :: A7
  Real(Std),       Intent(In)    :: Zs
  Type(DispOpts_), Intent(InOut) :: DispOpts
  ! FixedMetTime     :: Time of fixed met data (for FixedMet = .true. only).
  ! MaxParticles     :: Maximum number of particles.
  ! MaxFullParticles :: Maximum number of full particles (i.e. particles for which a
  !                     set of extra information is provided).
  ! MaxPuffs         :: Maximum number of puffs.
  ! MaxOriginalPuffs :: Maximum number of original puffs (i.e. puffs released at the
  !                     source and not created by puff splitting).
  ! ParticleCeiling  ::
  ! ParticleFactor   ::
  ! Zs               :: Height below which particles are subject to dry deposition
  ! SkewTime         :} Travel time over which the model allows for (i) skewness, (ii)
  ! VelMemTime       :} velocity memory, and (iii) inhomogeneity. These travel times
  ! InhomogTime      :} must satisfy SkewTime <= VelMemTime <= InhomogTime.
  ! MVelMemTime      :: unresolved mesoscale motions vel mem time
  ! PuffTime         :: Travel time over which puffs are used.
  ! SyncdT           :: Time interval at which particles and puffs are synchronised.
  ! DomainName       :: Name of computational domain.
  ! PuffInterval     :: Time interval between release of puffs (for non-fixed met cases).
  ! DeltaOpt         :: 0 indicates Delta = infinity; 1 indicates Delta limited by inhomogeneity; 2 indicates
  !                     Delta limited by inhomogeneity and by sigma_p.
  ! DeepConvection   :} Indicates that deep convection, radioactive decay,
  ! RadioactiveDecay :} agent decay, dry deposition, wet deposition, turbulence
  ! AgentDecay       :} and unresolved mesoscale motions are modelled.
  ! DryDep           :}
  ! WetDep           :}
  ! Turbulence       :}
  ! MesoscaleMotions :}
  ! VerticalVelocity :: Indicates that particles are advected in 3-dimensions.
  !                     If false, 2-d trajectories are generated. Defaults to true.
  ! Chemistry        :: Indicates that chemistry scheme is run.

  DispOpts%FixedMetTime = FixedMetTime

  If (MaxFullParticles > MaxParticles) Then
    Call Message('ERROR in reading dispersion options: Max # Full Particles should ' // &
                 'not exceed Max # Particles', 3)
  End If
  If (MaxOriginalPuffs > MaxPuffs) Then
    Call Message('ERROR in reading dispersion options: Max # Original Puffs should not exceed Max # Puffs', 3)
  End If
  DispOpts%MaxParticles     = MaxParticles
  DispOpts%MaxFullParticles = MaxFullParticles
  DispOpts%MaxPuffs         = MaxPuffs
  DispOpts%MaxOriginalPuffs = MaxOriginalPuffs

  DispOpts%ParticleCeiling  = ParticleCeiling ! check >= 1 <= MaxParticles $$
  DispOpts%ParticleFactor   = ParticleFactor  ! check >= 1.0 $$

  If (VelMemTime < SkewTime) Then
    Call Message('Error in InitDispOpts', 3)
  End If
  If (InhomogTime < VelMemTime) Then
    Call Message('Error in InitDispOpts', 3)
  End If
  DispOpts%SkewTime     = SkewTime
  DispOpts%sSkewTime    = Time2ShortTime(SkewTime)
  DispOpts%VelMemTime   = VelMemTime
  DispOpts%sVelMemTime  = Time2ShortTime(VelMemTime)
  DispOpts%InhomogTime  = InhomogTime
  DispOpts%sInhomogTime = Time2ShortTime(InhomogTime)
  DispOpts%MVelMemTime  = MVelMemTime
  DispOpts%sMVelMemTime = Time2ShortTime(MVelMemTime)

  DispOpts%PuffTime  = PuffTime
  DispOpts%sPuffTime = Time2ShortTime(PuffTime)

  If (SyncdT <= ZeroTime()) Then
    Call Message('FATAL ERROR: Sync Time must be > 0', 3)
  End If
  DispOpts%SyncdT = SyncdT

  If (Len_Trim(DomainName) > MaxCharLength) Then
    Call Message(                                                        &
           'FATAL ERROR in InitDispOpts: computational domain name (' // &
           Trim(DomainName)                                           // &
           ') is too long',                                              &
           3                                                             &
         )
  End If
  If (Len_Trim(DomainName) <= 0) Then
    Call Message('FATAL ERROR in InitDispOpts: computational domain name is blank', 3)
  End If
  DispOpts%DomainName = DomainName

  DispOpts%PuffInterval = PuffInterval
  DispOpts%DeltaOpt     = DeltaOpt

  ! If test to ensure sensible input values to DeepConvection
  If (DeepConvection .CIEq. 'no') Then
    DispOpts%DeepConvectionCode = Convection_None
  Else If (DeepConvection .CIEq. 'old') Then
    DispOpts%DeepConvectionCode = Convection_Old
  Else If (DeepConvection .CIEq. 'new') Then
    DispOpts%DeepConvectionCode = Convection_New
  Else
    Call Message('FATAL ERROR: Unsuitable input choice for "Deep Convection?"', 3)
  End If
      
  DispOpts%RadioactiveDecay = RadioactiveDecay
  DispOpts%AgentDecay       = AgentDecay
  DispOpts%DryDep           = DryDep
  DispOpts%WetDep           = WetDep
  DispOpts%Turbulence       = Turbulence
  DispOpts%MesoscaleMotions = MesoscaleMotions
  DispOpts%VerticalVelocity = VerticalVelocity
  DispOpts%Damping          = Damping
  DispOpts%Chemistry        = Chemistry
  DispOpts%PuffOpts%A1      = A1
  DispOpts%PuffOpts%A5      = A5
  DispOpts%PuffOpts%A7      = A7

  DispOpts%Zs               = Zs

End Subroutine InitDispOpts

!-------------------------------------------------------------------------------------------------------------

Subroutine AddDispOpts(DispOpts, DispOptses)
! Add a set of dispersion model options to a collection of sets of dispersion model
! options.

  Implicit None
  ! Argument list:
  Type(DispOpts_),   Intent(In)    :: DispOpts   ! A set of dispersion model options.
  Type(DispOptses_), Intent(InOut) :: DispOptses ! A collection of sets of dispersion model options.

  If (DispOptses%nDispOptses == MaxDispOptses) Then
    Call Message('FATAL ERROR in AddDispOpts', 3)
  End If

  DispOptses%nDispOptses                        = DispOptses%nDispOptses + 1
  DispOptses%DispOptses(DispOptses%nDispOptses) = DispOpts

End Subroutine AddDispOpts

!-------------------------------------------------------------------------------------------------------------

Subroutine WriteRestartFileDispState(Unit, DispState)
!

  Implicit None

  ! Argument list:
  Integer,          Intent(In)         :: Unit
  Type(DispState_), Intent(In), Target :: DispState
  ! Locals:
  Integer                       :: i
  Type(EulerianField_), Pointer :: EulerianField

  Write (Unit) DispState%LastCase
  Write (Unit) DispState%iDispOpts
  Write (Unit) DispState%iMetCase
  Write (Unit) DispState%iDomain
  Write (Unit) DispState%sLastDispTime
  Write (Unit) DispState%sFirstReqTime,     DispState%sLastReqTime
  Write (Unit) DispState%sFirstReleaseTime, DispState%sLastReleaseTime
  Write (Unit) DispState%sFirstMetTime,     DispState%sLastMetTime
  Write (Unit) DispState%Time
  Write (Unit) DispState%SyncTime
  Write (Unit) DispState%MetTime
  Write (Unit) DispState%PPsFinished

  Write (Unit) DispState%SpaceAllocated
  Write (Unit) DispState%MaxParticles
  Write (Unit) DispState%MaxFullParticles
  Write (Unit) DispState%MaxPuffs
  Write (Unit) DispState%MaxOriginalPuffs
  Write (Unit) DispState%nParticleSpecieses

  Write (Unit) DispState%nParticles, DispState%LastParticle
  Write (Unit) DispState%FreeParticleStack(1:DispState%LastParticle)
  Write (Unit) DispState%nParticleExtras, DispState%LastParticleExtra
  Write (Unit) DispState%FreeParticleExtraStack(1:DispState%LastParticleExtra)
  Write (Unit) DispState%iParticleExtras(1:DispState%LastParticle)
  Write (Unit) DispState%Particles(1:DispState%LastParticle) ! May need to subdivide to
                                                             ! prevent stack overflow $$
  Write (Unit) DispState%ParticleMasses(:, 1:DispState%LastParticle)
  Write (Unit) DispState%ParticleExtras(0:DispState%LastParticleExtra)

  Write (Unit) DispState%nOriginalPuffs, DispState%LastOriginalPuff
  Write (Unit) DispState%FreeOriginalPuffStack(1:DispState%LastOriginalPuff)
  Write (Unit) DispState%nPuffs, DispState%LastPuff
  Write (Unit) DispState%FreePuffStack(1:DispState%LastPuff)
  Write (Unit) DispState%iPuffs(1:DispState%LastOriginalPuff)
  Write (Unit) DispState%Puffs(1:DispState%LastPuff)
  Write (Unit) DispState%PuffMasses(:, 1:DispState%LastPuff)
  Write (Unit) DispState%PuffExtras(1:DispState%LastPuff)
  Write (Unit) DispState%TLast(:)
  Write (Unit) DispState%T(:)
  Write (Unit) DispState%SigA2Last(:)
  Write (Unit) DispState%SigA2(:)
  Write (Unit) DispState%SigA2Defined(:)

  Write (Unit) DispState%iULastParticle, DispState%iULastPuff, DispState%iULastOriginalPuff

  Write (Unit) DispState%ParticleFactorType
  Write (Unit) DispState%nParticleTimeSteps
  Write (Unit) DispState%nPuffTimeSteps
  Write (Unit) DispState%MaxTSpread

  Write (Unit) Size(DispState%NextReleaseTime)
  Do i = 1, Size(DispState%NextReleaseTime)
    Write (Unit) Shape(DispState%NextReleaseTime(i)%T)
    Write (Unit) DispState%NextReleaseTime(i)%T(:, :, :)
  End Do

  Write (Unit) Size(DispState%NextReleaseTimeM)
  Write (Unit) DispState%NextReleaseTimeM(:)

  Write (Unit) DispState%iZCoordMagl, DispState%iHCoordLatLong

  Write (Unit) DispState%PPsLostDueToFlow

  ! Eulerian model
  EulerianField => DispState%EulerianField

  Write (Unit) EulerianField%iOld
  Write (Unit) EulerianField%iNew
  Write (Unit) EulerianField%InitialiseFields

  Write (Unit) Associated(EulerianField%ZAboveGround)
  If (Associated(EulerianField%ZAboveGround)) Then
    Write (Unit) Shape(EulerianField%ZAboveGround)
    Write (Unit) EulerianField%ZAboveGround(:,:,:)
  End If

  Write (Unit) Associated(EulerianField%U)
  If (Associated(EulerianField%U)) Then
    Write (Unit) Shape(EulerianField%U)
    Write (Unit) EulerianField%U(:,:,:,:)
  End If

  Write (Unit) Associated(EulerianField%V)
  If (Associated(EulerianField%V)) Then
    Write (Unit) Shape(EulerianField%V)
    Write (Unit) EulerianField%V(:,:,:,:)
  End If

  Write (Unit) Associated(EulerianField%W)
  If (Associated(EulerianField%W)) Then
    Write (Unit) Shape(EulerianField%W)
    Write (Unit) EulerianField%W(:,:,:,:)
  End If

  Write (Unit) Associated(EulerianField%ModifiedW)
  If (Associated(EulerianField%ModifiedW)) Then
    Write (Unit) Shape(EulerianField%ModifiedW)
    Write (Unit) EulerianField%ModifiedW(:,:,:,:)
  End If

  Write (Unit) Associated(EulerianField%Density)
  If (Associated(EulerianField%Density)) Then
    Write (Unit) Shape(EulerianField%Density)
    Write (Unit) EulerianField%Density(:,:,:,:)
  End If

  Write (Unit) Associated(EulerianField%NewSource)
  If (Associated(EulerianField%NewSource)) Then
    Write (Unit) Shape(EulerianField%NewSource)
    !Write (Unit) EulerianField%NewSource(:,:,:,:)  $$ Currently zero 
  End If

  Write (Unit) Associated(EulerianField%Concentration)
  If (Associated(EulerianField%Concentration)) Then
    Write (Unit) Shape(EulerianField%Concentration)
    Write (Unit) EulerianField%Concentration(:,:,:,:,:)
  End If

  Write (Unit) Associated(EulerianField%Diffusivity)
  If (Associated(EulerianField%Diffusivity)) Then
    Write (Unit) Shape(EulerianField%Diffusivity)
    Write (Unit) EulerianField%Diffusivity(:,:,:,:)
  End If

  Write (Unit) Associated(EulerianField%Temperature)
  If (Associated(EulerianField%Temperature)) Then
    Write (Unit) Shape(EulerianField%Temperature)
    Write (Unit) EulerianField%Temperature(:,:,:,:)
  End If

  Write (Unit) Associated(EulerianField%Humidity)
  If (Associated(EulerianField%Humidity)) Then
    Write (Unit) Shape(EulerianField%Humidity)
    Write (Unit) EulerianField%Humidity(:,:,:,:)
  End If

  Write (Unit) Associated(EulerianField%dZdZ)
  If (Associated(EulerianField%dZdZ)) Then
    Write (Unit) Shape(EulerianField%dZdZ)
    Write (Unit) EulerianField%dZdZ(:,:,:,:)
  End If

  Write (Unit) Associated(EulerianField%Vd)
  If (Associated(EulerianField%Vd)) Then
    Write (Unit) Shape(EulerianField%Vd)
    Write (Unit) EulerianField%Vd(:,:,:,:)
  End If

  Write (Unit) Associated(EulerianField%WSedAtGround)
  If (Associated(EulerianField%WSedAtGround)) Then
    Write (Unit) Shape(EulerianField%WSedAtGround)
    Write (Unit) EulerianField%WSedAtGround(:,:,:,:)
  End If

  Write (Unit) Associated(EulerianField%Lambda)
  If (Associated(EulerianField%Lambda)) Then
    Write (Unit) Shape(EulerianField%Lambda)
    Write (Unit) EulerianField%Lambda(:,:,:,:)
  End If

  Write (Unit) Associated(EulerianField%TotalDepositionRate)
  If (Associated(EulerianField%TotalDepositionRate)) Then
    Write (Unit) Shape(EulerianField%TotalDepositionRate)
    Write (Unit) EulerianField%TotalDepositionRate(:,:,:)
  End If

  Write (Unit) Associated(EulerianField%DryDepositionRate)
  If (Associated(EulerianField%DryDepositionRate)) Then
    Write (Unit) Shape(EulerianField%DryDepositionRate)
    Write (Unit) EulerianField%DryDepositionRate(:,:,:)
  End If

  Write (Unit) Associated(EulerianField%WetDepositionRate)
  If (Associated(EulerianField%WetDepositionRate)) Then
    Write (Unit) Shape(EulerianField%WetDepositionRate)
    Write (Unit) EulerianField%WetDepositionRate(:,:,:)
  End If

End Subroutine WriteRestartFileDispState

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadRestartFileDispState(nSources, iFirstNewSource, NewSourcesDetected, Unit, DispState)
!

  Implicit None

  ! Argument list:
  Integer,          Intent(In)             :: nSources
  Integer,          Intent(InOut)         :: iFirstNewSource
  Logical,          Intent(InOut)         :: NewSourcesDetected
  Integer,          Intent(InOut)         :: Unit
  Type(DispState_), Intent(InOut), Target :: DispState
  ! Locals:
  Integer                       :: i
  Integer                       :: n1, n2, n3, n4, n11 ! array dimensions
  Type(EulerianField_), Pointer :: EulerianField
  Logical                       :: IsAssociated     ! Indicates the allocatable array should be allocated.
  Integer                       :: ArrayShape(5)    ! Shape of array to be allocated.
  Integer                       :: ErrorCode        ! Error code.

  Read (Unit) DispState%LastCase
  Read (Unit) DispState%iDispOpts
  Read (Unit) DispState%iMetCase
  Read (Unit) DispState%iDomain
  Read (Unit) DispState%sLastDispTime
  Read (Unit) DispState%sFirstReqTime,     DispState%sLastReqTime
  Read (Unit) DispState%sFirstReleaseTime, DispState%sLastReleaseTime
  Read (Unit) DispState%sFirstMetTime,     DispState%sLastMetTime
  Read (Unit) DispState%Time
  Read (Unit) DispState%SyncTime
  Read (Unit) DispState%MetTime
  Read (Unit) DispState%PPsFinished

  Read (Unit) DispState%SpaceAllocated
  Read (Unit) DispState%MaxParticles
  Read (Unit) DispState%MaxFullParticles
  Read (Unit) DispState%MaxPuffs
  Read (Unit) DispState%MaxOriginalPuffs
  Read (Unit) DispState%nParticleSpecieses

  Allocate(DispState%FreeParticleStack     (                                DispState%MaxParticles    ))
  Allocate(DispState%FreeParticleExtraStack(                                DispState%MaxFullParticles))
  Allocate(DispState%iParticleExtras       (                                DispState%MaxParticles    ))
  Allocate(DispState%Particles             (                                DispState%MaxParticles    ))
  Allocate(DispState%ParticleMasses        (DispState%nParticleSpecieses,   DispState%MaxParticles    ))
  Allocate(DispState%ParticleExtras        (                              0:DispState%MaxFullParticles))
  Allocate(DispState%FreeOriginalPuffStack (                                DispState%MaxOriginalPuffs))
  Allocate(DispState%FreePuffStack         (                                DispState%MaxPuffs        ))
  Allocate(DispState%iPuffs                (                                DispState%MaxOriginalPuffs))
  Allocate(DispState%Puffs                 (                                DispState%MaxPuffs        ))
  Allocate(DispState%PuffMasses            (DispState%nParticleSpecieses,   DispState%MaxPuffs        ))
  Allocate(DispState%PuffExtras            (                                DispState%MaxPuffs        ))
  Allocate(DispState%TLast                 (                                DispState%MaxOriginalPuffs))
  Allocate(DispState%T                     (                                DispState%MaxOriginalPuffs))
  Allocate(DispState%SigA2Last             (                                DispState%MaxOriginalPuffs))
  Allocate(DispState%SigA2                 (                                DispState%MaxOriginalPuffs))
  Allocate(DispState%SigA2Defined          (                                DispState%MaxOriginalPuffs))

  Read (Unit) DispState%nParticles, DispState%LastParticle
  Do i = 1, Size(DispState%FreeParticleStack) ! This saves space in restart file.
    DispState%FreeParticleStack(i) = i
  End Do
  Read (Unit) DispState%FreeParticleStack(1:DispState%LastParticle)
  Read (Unit) DispState%nParticleExtras, DispState%LastParticleExtra
  Do i = 1, Size(DispState%FreeParticleExtraStack) ! This saves space in restart file.
    DispState%FreeParticleExtraStack(i) = i
  End Do
  Read (Unit) DispState%FreeParticleExtraStack(1:DispState%LastParticleExtra)
  Read (Unit) DispState%iParticleExtras(1:DispState%LastParticle)
  Read (Unit) DispState%Particles(1:DispState%LastParticle) ! May need to subdivide to prevent
                                                            ! stack overflow $$
  Read (Unit) DispState%ParticleMasses(:, 1:DispState%LastParticle)
  Read (Unit) DispState%ParticleExtras(0:DispState%LastParticleExtra)

  Read (Unit) DispState%nOriginalPuffs, DispState%LastOriginalPuff
  Do i = 1, Size(DispState%FreeOriginalPuffStack) ! This saves space in restart file.
    DispState%FreeOriginalPuffStack(i) = i
  End Do
  Read (Unit) DispState%FreeOriginalPuffStack(1:DispState%LastOriginalPuff)
  Read (Unit) DispState%nPuffs, DispState%LastPuff
  Do i = 1, Size(DispState%FreePuffStack) ! This saves space in restart file.
    DispState%FreePuffStack(i) = i
  End Do
  Read (Unit) DispState%FreePuffStack(1:DispState%LastPuff)
  Read (Unit) DispState%iPuffs(1:DispState%LastOriginalPuff)
  Read (Unit) DispState%Puffs(1:DispState%LastPuff)         ! May need to subdivide
                                                            ! to prevent stack overflow $$
  Read (Unit) DispState%PuffMasses(:, 1:DispState%LastPuff)
  Read (Unit) DispState%PuffExtras(1:DispState%LastPuff)
  Read (Unit) DispState%TLast(:)
  Read (Unit) DispState%T(:)
  Read (Unit) DispState%SigA2Last(:)
  Read (Unit) DispState%SigA2(:)
  Read (Unit) DispState%SigA2Defined(:)

  Read (Unit) DispState%iULastParticle, DispState%iULastPuff, DispState%iULastOriginalPuff

  Read (Unit) DispState%ParticleFactorType
  Read (Unit) DispState%nParticleTimeSteps
  Read (Unit) DispState%nPuffTimeSteps
  Read (Unit) DispState%MaxTSpread

  ! Check and adjust number of sources
  Read (Unit) n1
  If (n1 > nSources) Then
    Call Message(                                                                                        &
           'FATAL ERROR: the number of sources in the input file(s) has been reduced on restart'      // &
           ' - this can have unexpected consequences and it is recommended to start the run from new',   &
           3                                                                                             &
         )
  Else If (n1 < nSources) Then
    Call Message(                                                                                        &
           'WARNING: the number of sources in the input file(s) has been increased on restart'        // &
           ' - this is allowed provided that the new sources are listed after all existing sources',     &
           1                                                                                             &
         )
    NewSourcesDetected = .true.
    iFirstNewSource    = n1 + 1
    Allocate(DispState%NextReleaseTime(nSources))
    Allocate(DispState%NextReleaseTimeM(nSources))
  Else
    Allocate(DispState%NextReleaseTime(n1))
    Allocate(DispState%NextReleaseTimeM(n1))
  End If

  Do i = 1, n1
    Read (Unit) n2, n3, n4
    Allocate(DispState%NextReleaseTime(i)%T(n2, n3, n4))
    Read (Unit) DispState%NextReleaseTime(i)%T(:, :, :)
  End Do

  Read (Unit) n11
  If (n11 /= n1) Then
    Call Message(                                                                        &
           'UNEXPECTED FATAL ERROR: the arrays NextReleaseTime and NextReleaseTimeM ' // &
           'should have the same size',                                                  &
           4                                                                             &
         )
  End If
  Read (Unit) DispState%NextReleaseTimeM(1:n1)

  Read (Unit) DispState%iZCoordMagl, DispState%iHCoordLatLong

  Read (Unit) DispState%PPsLostDueToFlow

  ! Eulerian model
  EulerianField => DispState%EulerianField

  Read (Unit) EulerianField%iOld
  Read (Unit) EulerianField%iNew
  Read (Unit) EulerianField%InitialiseFields

  ! ZAboveGround
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:3)
    Allocate(EulerianField%ZAboveGround(ArrayShape(1), ArrayShape(2), ArrayShape(3)), &
             Stat = ErrorCode                                                         &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%ZAboveGround(:,:,:)
  End If

  ! U
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(EulerianField%U(-1:ArrayShape(1) - 2, -1:ArrayShape(2) - 2, ArrayShape(3), ArrayShape(4)), &
             Stat = ErrorCode                                                                           &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%U(:,:,:,:)
  End If

  ! V
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(EulerianField%V(-1:ArrayShape(1) - 2, -1:ArrayShape(2) - 2, ArrayShape(3), ArrayShape(4)), &
             Stat = ErrorCode                                                                           &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%V(:,:,:,:)
  End If

  ! W
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(EulerianField%W(-1:ArrayShape(1) - 2, -1:ArrayShape(2) - 2, ArrayShape(3), ArrayShape(4)), &
             Stat = ErrorCode                                                                           &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%W(:,:,:,:)
  End If

  ! ModifiedW
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(EulerianField%ModifiedW(-1:ArrayShape(1) - 2, -1:ArrayShape(2) - 2, ArrayShape(3), ArrayShape(4)), &
             Stat = ErrorCode                                                                                   &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%ModifiedW(:,:,:,:)
  End If

  ! Density
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(EulerianField%Density(-1:ArrayShape(1) - 2, -1:ArrayShape(2) - 2, ArrayShape(3), ArrayShape(4)), &
             Stat = ErrorCode                                                                                 &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%Density(:,:,:,:)
  End If

  ! NewSource
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(EulerianField%NewSource(-1:ArrayShape(1) - 2, -1:ArrayShape(2) - 2, ArrayShape(3), ArrayShape(4)), &
             Stat = ErrorCode                                                                                   &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    !Read (Unit, End = 1, Err = 2) EulerianField%NewSource(:,:,:,:) $$ Currently zero
    EulerianField%NewSource(:,:,:,:) = 0.0
  End If

  ! Concentration
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:5)
    Allocate(EulerianField%Concentration(-1:ArrayShape(1) - 2, &
                                         -1:ArrayShape(2) - 2, &
                                         ArrayShape(3),        &
                                         ArrayShape(4),        &
                                         ArrayShape(5)),       &
             Stat = ErrorCode                                  &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%Concentration(:,:,:,:,:)
  End If

  ! Diffusivity
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(EulerianField%Diffusivity(ArrayShape(1), ArrayShape(2), ArrayShape(3), ArrayShape(4)), &
             Stat = ErrorCode                                                                       &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%Diffusivity(:,:,:,:)
  End If

  ! Temperature
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(EulerianField%Temperature(ArrayShape(1), ArrayShape(2), ArrayShape(3), ArrayShape(4)), &
             Stat = ErrorCode                                                                       &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%Temperature(:,:,:,:)
  End If

  ! Humidity
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(EulerianField%Humidity(ArrayShape(1), ArrayShape(2), ArrayShape(3), ArrayShape(4)), &
             Stat = ErrorCode                                                                       &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%Humidity(:,:,:,:)
  End If


  ! dZdZ
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(EulerianField%dZdZ(ArrayShape(1), ArrayShape(2), ArrayShape(3), ArrayShape(4)), &
             Stat = ErrorCode                                                                &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%dZdZ(:,:,:,:)
  End If

  ! Vd
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(EulerianField%Vd(ArrayShape(1), ArrayShape(2), ArrayShape(3), ArrayShape(4)), &
             Stat = ErrorCode                                                              &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%Vd(:,:,:,:)
  End If

  ! WSedAtGround
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(EulerianField%WSedAtGround(ArrayShape(1), ArrayShape(2), ArrayShape(3), ArrayShape(4)), &
             Stat = ErrorCode                                                                        &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%WSedAtGround(:,:,:,:)
  End If

  ! Lambda
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(EulerianField%Lambda(ArrayShape(1), ArrayShape(2), ArrayShape(3), ArrayShape(4)), &
             Stat = ErrorCode                                                                  &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%Lambda(:,:,:,:)
  End If

  ! Total deposition rate
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:3)
    Allocate(EulerianField%TotalDepositionRate(ArrayShape(1), ArrayShape(2), ArrayShape(3)), &
             Stat = ErrorCode                                                                &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%TotalDepositionRate(:,:,:)
  End If

  ! Dry deposition rate
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:3)
    Allocate(EulerianField%DryDepositionRate(ArrayShape(1), ArrayShape(2), ArrayShape(3)), &
             Stat = ErrorCode                                                              &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%DryDepositionRate(:,:,:)
  End If

  ! Wet deposition rate
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:3)
    Allocate(EulerianField%WetDepositionRate(ArrayShape(1), ArrayShape(2), ArrayShape(3)), &
             Stat = ErrorCode                                                              &
            )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for Eulerian fields', 3)
    Read (Unit, End = 1, Err = 2) EulerianField%WetDepositionRate(:,:,:)
  End If


  Return

1 Continue
  Call Message('FATAL ERROR: The restart file is shorter than expected', 3)

2 Continue
  Call Message('FATAL ERROR: An error occurred when reading the restart file', 3)

End Subroutine ReadRestartFileDispState

!-------------------------------------------------------------------------------------------------------------

Subroutine InitDispState(DispState)

  Implicit None
  ! Argument list:
  Type(DispState_), Intent(Out) :: DispState ! Dispersion model state.

  DispState%SpaceAllocated = .false.

End Subroutine InitDispState

!-------------------------------------------------------------------------------------------------------------

Subroutine RestartAdjustmentForInputChanges(      &
             Grids, Domains,                      &
             Sources, SizeDists,                  &
             Reqs,                                &
             NewSourcesDetected, iFirstNewSource, &
             DispState                            &
           )
! Adjusts the dispersion model state to account for permitted changes to the input file on restart.
! The allowed changes are i) stopping an existing non-instantaneous source, ii) addition of extra sources
! (but only after all pre-existing sources), and iii) addition of extra output requests (again only after
! all existing requests). $$ Note that other sections of the code do not support iii) yet, so it is
! currently not advisable to add extra output requests at restart.

  Implicit None
  ! Argument list:
  Type(Grids_),     Intent(In), Target :: Grids              ! Collection of grids.
  Type(Domains_),   Intent(In)         :: Domains            ! Collection of domains.
  Type(Sources_),   Intent(In), Target :: Sources            ! Collection of sources.
  Type(SizeDists_), Intent(In), Target :: SizeDists          ! Collection of particle size distributions.
  Type(Reqs_),      Intent(In)         :: Reqs               ! Collection of requirements.
  Logical,          Intent(In)         :: NewSourcesDetected ! Indicates that new sources have been detected.
  Integer,          Intent(In)         :: iFirstNewSource    ! Index of first new source term in input file.
  Type(DispState_), Intent(InOut)      :: DispState          ! Dispersion model state.
  ! Locals:
  Integer                  :: i, j, k, l
  Integer                  :: ArrayDim(3)
  Integer                  :: nExistingSources
  Integer                  :: iNewSource
  Integer                  :: PrevParticleIndex
  Real(Std)                :: TimeFrac
  Type(ShortTime_)         :: FirstSourceTime
  Type(ShortTime_)         :: LastSourceTime
  Type(ShortTime_)         :: NextReleaseTime(Sources%nSources)
  Type(ShortTime_)         :: LastReleaseTime(Sources%nSources)
  Type(ShortTime_)         :: NextReleaseTimeForSource
  Type(HGrid_),    Pointer :: HGrid
  Type(Source_),   Pointer :: Source
  Type(SizeDist_), Pointer :: SizeDist
  Type(ShortTime_)         :: PrevParticleReleaseTime
  

  ! Update the variables in DispState controlling the model run to account for new sources or output requests:
  !   sFirstReleaseTime, sLastReleaseTime
  !   sFirstReqTime, sLastReqTime, sFirstMetTime, sLastMetTime, sLastDispTime
  ! Note that these controls may not be fully implemented during the present case (e.g. if DispState%Time
  ! is later than DispState%sFirstReleaseTime or DispState%sFirstReqTime then any changes in these latter
  ! two times will have no effect). Furthermore, the model will not rerun over a time period that it has
  ! previously covered, even if new requests/sources have been added within that period. However all future
  ! events starting from the restart value of DispState%Time should be correctly handled.

  ! 1) Compute interim first and last release times for all sources (including the existing sources).
  !    The interim times may not be the most appropriate for fixed met cases and will be revised below in 3).

  Call FirstLastReleaseTimes(                                       &
         InfFutureShortTime(), Domains%Domains(DispState%iDomain),  &
         Sources,                                                   &
         FirstSourceTime,             LastSourceTime,               &
         DispState%sFirstReleaseTime, DispState%sLastReleaseTime,   &
         NextReleaseTime,             LastReleaseTime               &
       )

  ! 2) Compute first and last times that need consideration for output (including the existing output).

  Call FirstLastReqTimes(                                 &
         FirstSourceTime, DispState%sFirstReleaseTime,    &
         Domains%Domains(DispState%iDomain),              &
         Grids, Reqs,                                     &
         DispState%sFirstReqTime, DispState%sLastReqTime, &
         DispState%sFirstMetTime, DispState%sLastMetTime, &
         DispState%sLastDispTime                          &
       )

  ! 3) Compute first and last release times for all sources (including the existing sources). This second call
  !    provides a more appropriate release time for (fixed met) sources with emissions in the infinite past.
  !    Note that these times will not necessarily be applied, but will be used as constraints below.

  Call FirstLastReleaseTimes(                                         &
         DispState%sFirstReqTime, Domains%Domains(DispState%iDomain), &
         Sources,                                                     &
         FirstSourceTime,             LastSourceTime,                 &
         DispState%sFirstReleaseTime, DispState%sLastReleaseTime,     &
         NextReleaseTime,             LastReleaseTime                 &
       )

  ! 4) Update FirstMetTime to take account of the requirement for met data from the first release time.
  
  DispState%sFirstMetTime = TMin(DispState%sFirstMetTime, DispState%sFirstReleaseTime)
  
  If (NewSourcesDetected) Then
    nExistingSources = iFirstNewSource - 1
    iNewSource       = iFirstNewSource
  Else
    nExistingSources = Sources%nSources
    iNewSource       = Sources%nSources + 1
  End If
  
  ! For existing sources, check that the source has not been stopped.
  
  Do i = 1, nExistingSources
    
    Source => Sources%Sources(i)
    
    If (Source%SourceType == S_Generic .Or. Source%SourceType == S_IterativePlumeModel) Then
      NextReleaseTimeForSource = DispState%NextReleaseTime(i)%T(1, 1, 1)
    Else If (Source%SourceType == S_Dust .or. Source%SourceType == S_SeaSalt) Then
      NextReleaseTimeForSource = InfFutureShortTime()
      ArrayDim(:) = Shape(DispState%NextReleaseTime(i)%T)
      Do j = 1, ArrayDim(1)
      Do k = 1, ArrayDim(2)
      Do l = 1, ArrayDim(3)
        NextReleaseTimeForSource = TMin(NextReleaseTimeForSource, DispState%NextReleaseTime(i)%T(j, k, l))
      End Do
      End Do
      End Do
    End If
    
    If (.not.IsInfFuture(NextReleaseTimeForSource)) Then
      
      ! Only treat non-instantaneous sources here
      ! In the event that an existing non-instantaneous source is modified such that the start time and
      ! stop time become the same (which would appear to look like an instantaneous release and so the
      ! following code block would not be executed). However this situation would be picked up when the
      ! source information is parsed (and a suitable message reported). In particular no new particles
      ! would be released.
      
      If (Source%StartTime < Source%StopTime) Then
        
        If (Source%sStopTime <= NextReleaseTimeForSource) Then
          ! Deactivate source by setting next release time to the infinite future
          DispState%NextReleaseTimeM(i)           = InfFutureShortTime()
          DispState%NextReleaseTime(i)%T(:, :, :) = InfFutureShortTime()
        End If
        
        ! Produce a warning if the new stop time is before the restart time
        If (Source%StopTime < DispState%Time) Then
          Call Message(                                                                                       &
             'WARNING: the stop time of an existing source with the name "'                                // &
             Trim(Source%Name)                                                                             // &
             '" has been changed in a restarted run to an earlier time than the restart and any particles' // &
             ' that have already been released by that source prior to the restart may still be active.'   // &
             ' Consider restarting from an earlier restart file.',                                            &
             1                                                                                                &
           )
        ! Also produce a warning if existing particles already contain mass representative of the killed time
        ! period. The checking here will not capture all eventualities and could be regarded as somewhat
        ! unnecessarily complicated, but it should help to reduce warnings where the effects are trivial.
        Else If (Source%sStopTime < NextReleaseTimeForSource) Then
          ! Although it is difficult to establish categorically whether mass has already been released
          ! for the period after which the source should have stopped, it is possible to constrain when
          ! a warning message might be issued (e.g. when the time period is very small relative to the
          ! interval that the previously released particle represents). This helps to avoid warnings
          ! which might be caused by small rounding errors in times, say. A warning will be suppressed
          ! if the possible fraction of mass after the revised stop time is less than 5% of the total
          ! mass on the last particle that was released. 
          ! First determine the release time of the last particle released by this source
          PrevParticleReleaseTime = InfPastShortTime()
          PrevParticleIndex       = -1
          Do k = 1, DispState%LastParticle
            If (DispState%Particles(k)%iSource /= i) Cycle
            If (DispState%Particles(k)%T0 > PrevParticleReleaseTime) Then
              PrevParticleReleaseTime = DispState%Particles(k)%T0
              PrevParticleIndex       = k
            End If
          End Do
          ! Then produce warning if particle is still active and there is a significant fraction of the time
          ! interval after the source has stopped
          If (PrevParticleIndex /= -1) Then
            TimeFrac = ShortTime2RealTime(NextReleaseTimeForSource - Source%sStopTime)        /  &
                       ShortTime2RealTime(NextReleaseTimeForSource - PrevParticleReleaseTime)
            If (ParticleActive(DispState%Particles(PrevParticleIndex)) .and. TimeFrac > 0.05) Then
              Call Message(                                                                               &
                 'WARNING: the stop time of an existing source with the name "'                        // &
                 Trim(Source%Name)                                                                     // &
                 '" has been changed in a restarted run, but some mass may have already been released' // &
                 ' representative of the time period after the new stop time',                            &
                 1                                                                                        &
              )
            End If
          End If
        End If
        
      End If
      
    End If
    
  End Do

  ! Update the variables DispState%NextReleaseTime and DispState%NextReleaseTimeM for any new sources found.
  ! For each new source, allocate the required 3-d time array and set the next release time to be the later
  ! of the calculated next release time from step 3) and DispState%Time.
  
  Do i = iNewSource, Sources%nSources
    
    Source => Sources%Sources(i)
    
    If (Source%SourceType == S_Generic .Or. Source%SourceType == S_IterativePlumeModel) Then
      Allocate(DispState%NextReleaseTime(i)%T(1, 1, 1))
    Else If (Source%SourceType == S_Dust .or. Source%SourceType == S_SeaSalt) Then
      HGrid    => Grids%HGrids(Source%iHGrid)
      SizeDist => SizeDists%SizeDists(Source%iSizeDist)
      Allocate(DispState%NextReleaseTime(i)%T(HGrid%nX, HGrid%nY, SizeDist%nSizeRanges + 1))
                                                                  !$$ +1 is probably unnecessary - check
    End If
    
    DispState%NextReleaseTimeM(i)           = TMax(NextReleaseTime(i), Time2ShortTime(DispState%Time))
    DispState%NextReleaseTime(i)%T(:, :, :) = DispState%NextReleaseTimeM(i)
    
    ! When the new source has a stop time before the restart time
    ! (issue a warning message and reset the next release time to the infinite future)
    If (LastReleaseTime(i) < Time2ShortTime(DispState%Time)) Then
      Call Message(                                                                             &
         'WARNING: a new source with the name "'                                             // &
         Trim(Source%Name)                                                                   // &
         '" has been added in a restarted run that has a stop time earlier than the restart' // &
         ' - this source will have no effect. '                                              // &
         'Consider restarting from an earlier restart file to capture this source.',            &
         1                                                                                      &
       )
       
      DispState%NextReleaseTimeM(i)           = InfFutureShortTime()
      DispState%NextReleaseTime(i)%T(:, :, :) = DispState%NextReleaseTimeM(i)
      
    ! Produce a warning if new source is set to start before restart time
    Else If (NextReleaseTime(i) < Time2ShortTime(DispState%Time)) Then
      Call Message(                                                                                          &
         'WARNING: a new source with the name "'                                                          // &
         Trim(Source%Name)                                                                                // &
         '" has been added in a restarted run that has a start time earlier than the restart'             // &
         ' - the source will only be active from the restart time. '                                      // &
         'Consider restarting from an earlier restart file to capture earlier releases from this source.',   &
         1                                                                                                   &
       )
    End If
    
  End Do

  ! Particles and puffs have finished if all previous particles and puffs have finished and no further
  ! ones are due to be released.
  DispState%PPsFinished = DispState%PPsFinished .and.                                   &
                          (Time2ShortTime(DispState%Time) > DispState%sLastReleaseTime)

End Subroutine RestartAdjustmentForInputChanges

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForNextCase(                                            &
             iCase,                                                       &
             Coords, Grids, Domains, Specieses, Sources, SizeDists, Reqs, &
             MainOpts, MultiCaseOpts, DispOptses,                         &
             Results, DispState,                                          &
             EndOfCase, EndOfCases                                        &
           )
! Prepares the dispersion model state for running the next case.

  Implicit None
  ! Argument list:
  Integer,              Intent(In)         :: iCase         ! Number of case.
  Type(Coords_),        Intent(In)         :: Coords        ! Collection of coord systems.
  Type(Grids_),         Intent(In), Target :: Grids         ! Collection of grids.
  Type(Domains_),       Intent(In)         :: Domains       ! Collection of domains.
  Type(Specieses_),     Intent(In)         :: Specieses     ! Collection of specieses.
  Type(Sources_),       Intent(In), Target :: Sources       ! Collection of sources.
  Type(SizeDists_),     Intent(In), Target :: SizeDists     ! Collection of particle size distributions.
  Type(Reqs_),          Intent(In)         :: Reqs          ! Collection of requirements.
  Type(MainOpts_),      Intent(In)         :: MainOpts      ! Main options.
  Type(MultiCaseOpts_), Intent(In)         :: MultiCaseOpts ! Multiple case options.
  Type(DispOptses_),    Intent(In), Target :: DispOptses    ! Collection of sets of dispersion model
                                                            ! options.
  Type(Results_),       Intent(InOut)      :: Results
  Type(DispState_),     Intent(InOut)      :: DispState     ! Dispersion model state.
  Logical,              Intent(InOut)      :: EndOfCase     ! Indicates that there is no more to do for
                                                            ! this case.
  Logical,              Intent(InOut)      :: EndOfCases    ! Indicates that there are no more cases to be
                                                            ! run.
  ! Locals:
  Type(Time_)              :: MaxPPTimeSpread ! Max time for which fields need to be remembered.
  Type(DispOpts_), Pointer :: DispOpts ! Abbreviation.
  Integer                  :: i        ! Loop index.
  Logical                  :: Error    ! Error flag.
  Type(ShortTime_)         :: Dummy
  Type(ShortTime_)         :: TDummy(Sources%nSources)
  Type(ShortTime_)         :: FirstSourceTime
  Type(ShortTime_)         :: LastSourceTime
  Type(ShortTime_)         :: NextReleaseTime(Sources%nSources)
  Type(HGrid_),    Pointer :: HGrid
  Type(Source_),   Pointer :: Source
  Type(SizeDist_), Pointer :: SizeDist

! $$ in following better to talk of 'field reqs' and 'fields' (the values computed to meet the req). eg
! the fields are processed, not the field reqs.
! $$ any mods needed for chemistry?
!
! The computational domain
! ------------------------
!
! The computational domain (CD) limits the domain of physical interest by space, time and particle/puff travel
! time. CDEndTime < CDStartTime, or CDStartTime = infinite future or CDEndTime = infinite past: fatal error.
! (true for all domains).
!
! Note fixed met option implies the flow module domains are unbounded in time.
! Particle/puff reuse option implies fixed met.
!
! Time spread of particles/puffs
! ------------------------------
!
! Dispersion from material released at nearby times is similar, and, for fixed met, identical. To take
! advantage of this, particles/puffs can have a spread in time which means that they represent material
! released at a range of times. If a particle/puff has a time spread, the particle/puff at a given nominal
! time will contribute to output for a range of times.
!
! MaxPPTimeLag = max interval between trailing time-edge and nominal time of particles/puffs.
!
! If MaxPPTimeLag > 0 particles/puffs may need to be followed beyond the end time of the computational domain
! in order to treat all physics occuring in the computational domain. Particles/puffs are not released before
! the start of the computational domain (see below) and so there is not a similar issue at the start of the
! computational domain.
!
! Particle/puff reuse option
! --------------------------
!
! As noted above, for fixed met, a release at any one time is representative of any other and so releases made
! at one time can be 'reused' for another. The particle/puff reuse option (which requires fixed met)
! implements this and causes (i) any given source to release all its material together at the same time, and
! (ii) each particle or puff to be given a time spread equal to the release duration.
!
! How time spread of particles/puffs is used with and without the particle/puff reuse option
! ------------------------------------------------------------------------------------------
!
! (a) Without the  particle/puff reuse option:
!     Time spread occurs for puffs from non-instantaneous sources. The time spread covers the interval between
!     puff releases which is specified by the user in order to resolve the changing met. The nominal time
!     associated with a puff and the time for which they get met data is the trailing time-edge of the puff.
!     This avoids non-causal influence of met and the need for met beyond the times of interest. The time
!     spread is triangular in shape to provide a smooth match between successively released puffs. When the
!     material from puffs released at different times overlaps significantly, the time spread may be turned
!     off.
! (b) With the particle/puff reuse option:
!     Time spread occurs for particles/puffs from non-instantaneous sources. The time spread covers the entire
!     source duration (restricted by the computational domain) with all material for the source being released
!     at the same time. The nominal time associated with a particle/puff (but not the time for which they get
!     met data because the met data is fixed) is the trailing time-edge of the particle/puff unless the
!     release starts at -infinity (after restricting by the start time of the computational domain) when this
!     is impossible. In that case a point within the particle/puff is selected to be the nominal time. The
!     choice determines the model time at which the particle/puff is released because, at release, the
!     particle's/puff's nominal time equals the model time at which the particle/puff is released. The choice
!     is made so that the model time of the release occurs at or before that of all other particles/puffs (see
!     below for details). The time spread is rectangular in shape to cover the entire source duration
!     (restricted by the computational domain) with a single release time.
!
! Particle/puff inactivation
! --------------------------
!
! Particles/puffs are inactivated
! (i)   when leaving the computation domain (CD) in space
! (ii)  when travel time crosses the CD travel time limit
! (iii) when travel time crosses the source travel time limit
! (iv)  when the trailing time-edge of the particle/puff crosses the CD end time
! (v)   if certain errors associated with flow modules occur
! (vi)  when they are inactivated as part of the release mechanism
! (vii) when no more contributions to the output will be made.
!
! LastPPTime = first time with no particles/puffs and no releases due.
!
! Release times
! -------------
!
! First and last release times for each source are calculated taking account of any restrictions caused by the
! computational domain and the need to use instantaneous sources with the particle/puff reuse option. The
! actual range of release times may be further restricted by the overall stop time for the model if not all
! releases are necessary (see below).
!
! FirstReleaseForSource:
!     if no release possible: +infinity; else: max(CDStart, SourceStart).
!     if result = -infinity and if reuse: set to first of the finite max(CDStart, SourceStart) or
!                                         min(CDEnd, SourceStop) values or, if no such values, an arbitrary
!                                         value (the same for all such cases). This will put
!                                         FirstReleaseForSource within the true source period and at or before
!                                         all other FirstReleaseForSource values. When there are
!                                         particles/puffs present from any sources that fall in this if-block,
!                                         MaxPPTimeLag = +infinity; otherwise MaxPPTimeLag = 0.
!     if result = -infinity: fatal error.
!
! LastReleaseForSource:
!     if reuse: min(CDEnd, FirstReleaseForSource); else: min(CDEnd, SourceStop).
!
! FirstRelease = min(FirstReleaseForSource) over sources with releases possible.
!
! LastRelease  = max(LastReleaseForSource)  over sources with releases possible.
!
! Req times
! ---------
!
! The first and last times required for calculating and outputting a req are calculated in order to fix the
! run start and stop times, taking account of any restrictions caused by the computational domain, the need to
! obtain appropriate contributions to any averaging period, and the need to output the req at the right time.
!
! FirstReqTime
! LastReqTime
! FirstFlowTime
! LastFlowTime
! LastDispTime
!
! Run start and stop times
! ------------------------
!
! The run start and stop times are controlled (assuming no fatal errors) by the computational domain, the
! release times, and the output requirements.
!
! StartTime = min(FirstReqTime, FirstRelease).
!
! StopTime occurs when the following is true:
!                 { LastReqTime + MaxPPTimeLag  if before demise of all tracers
!         Time >= {
!                 { LastReqTime                 if at or after demise of all tracers
! (-infinity + anything = -infinity here and causes the run to stop at the start).
!
! StopTime will be infinite if and only if LastReqTime /= -infinity and one of the following occurs:
! (i)  LastReqTime = +infinity,
! (ii) MaxPPTimeLag = +infinity and CDomain is unbounded in space and travel time and not all sources
!      have a finite travel time limit.
!
! If StopTime will be infinite and there is no means to suspend the model (i.e. through restart options or arg
! list): fatal error
!
! Met and flow times
! ------------------
! Dispersion times
! ----------------

  ! 1) Set DispState%iDispOpts, DispState%iMetCase, LastCase and EndOfCases.

  ! Prepare first case.
  If (iCase == 1) Then

    DispState%iDispOpts = 1
    DispState%iMetCase  = 1

  ! Prepare other cases if a multiple case run.
  Else

    If (DispState%iMetCase < MultiCaseOpts%MetEnsembleSize) Then
      DispState%iMetCase = DispState%iMetCase + 1
    Else If (DispState%iDispOpts < MultiCaseOpts%DispOptsesEnsembleSize) Then
      DispState%iDispOpts = DispState%iDispOpts + 1
      DispState%iMetCase  = 1
    Else
      EndOfCases = .true.
      Return
    End If

  End If

  DispState%LastCase = DispState%iDispOpts == MultiCaseOpts%DispOptsesEnsembleSize .and. &
                       DispState%iMetCase  == MultiCaseOpts%MetEnsembleSize

  DispOpts => DispOptses%DispOptses(DispState%iDispOpts)

  ! 2) Computational domain.

  DispState%iDomain = FindDomainIndex(DispOpts%DomainName, Domains, Error)
  If (Error) Then
    If (MultiCaseOpts%DispOptsesEnsembleSize > 1) Then ! $$ any value to this?
      EndOfCase = .true.
      Call Message('Error in PrepareForNextCase: Domain ' // Trim(DispOpts%DomainName) // ' not found.', 2)
      Return
    Else
      Call Message('Error in PrepareForNextCase: Domain ' // Trim(DispOpts%DomainName) // ' not found.', 3)
    End If
  End If

  ! $$ Check and note that this routine will work with no sources, requirements
  ! or dispersion options - may need to trap start time of case = infinity

  ! 3-) Allocate NextReleaseTimeM.
  Allocate(DispState%NextReleaseTimeM(Sources%nSources))

  ! 3) First and last release times for sources.

  Call FirstLastReleaseTimes(                                       &
         InfFutureShortTime(), Domains%Domains(DispState%iDomain),  &
         Sources,                                                   &
         FirstSourceTime,             LastSourceTime,               &
         DispState%sFirstReleaseTime, DispState%sLastReleaseTime,   &
         NextReleaseTime,             TDummy                        &
       )

  Do i = 1, Sources%nSources
    DispState%NextReleaseTimeM(i) = NextReleaseTime(i)
  End Do

  DispState%PPsFinished = IsInfFuture(DispState%sFirstReleaseTime)

  ! 4) First and last times that need consideration for output.

  Call FirstLastReqTimes(                                 &
         FirstSourceTime, DispState%sFirstReleaseTime,    &
         Domains%Domains(DispState%iDomain),              &
         Grids, Reqs,                                     &
         DispState%sFirstReqTime, DispState%sLastReqTime, &
         DispState%sFirstMetTime, DispState%sLastMetTime, &
         DispState%sLastDispTime                          &
       )
  ! $$ Error trap for no output required.

  ! 3a) First and last release times for sources. Note this second call isn't necessary, but it will provide a
  ! more appropriate release time for sources with emissions in the infinite past.

  Call FirstLastReleaseTimes(                                         &
         DispState%sFirstReqTime, Domains%Domains(DispState%iDomain), &
         Sources,                                                     &
         FirstSourceTime,             LastSourceTime,                 &
         DispState%sFirstReleaseTime, DispState%sLastReleaseTime,     &
         NextReleaseTime,             TDummy                          &
       )

  Do i = 1, Sources%nSources
    DispState%NextReleaseTimeM(i) = NextReleaseTime(i)
  End Do

  DispState%PPsFinished = IsInfFuture(DispState%sFirstReleaseTime)

  ! 4+) Allocate NextReleaseTime and duplicate NextRelaseTime for sources with grids.

  ! For multiple cases, we should deallocate or avoid reallocation $$

  Allocate(DispState%NextReleaseTime(Sources%nSources))
  Do i = 1, Sources%nSources
    Source => Sources%Sources(i)
    If (Source%SourceType == S_Generic .Or. Source%SourceType == S_IterativePlumeModel) Then
      Allocate(DispState%NextReleaseTime(i)%T(1, 1, 1))
    Else If (Source%SourceType == S_Dust .or. Source%SourceType == S_SeaSalt) Then
      HGrid    => Grids%HGrids(Source%iHGrid)
      SizeDist => SizeDists%SizeDists(Source%iSizeDist)
      Allocate(DispState%NextReleaseTime(i)%T(HGrid%nX, HGrid%nY, SizeDist%nSizeRanges + 1))
                                                                  !$$ +1 is probably unnecessary - check
    End If
  End Do

  Do i = 1, Sources%nSources
    DispState%NextReleaseTime(i)%T(:, :, :) = NextReleaseTime(i)
  End Do

  ! 5) Start time of case.

  DispState%Time = ShortTime2Time(TMin(DispState%sFirstReqTime, DispState%sFirstReleaseTime))

  ! 6) SyncTime and MetTime.

  ! FirstMetTime.
  DispState%sFirstMetTime = TMin(DispState%sFirstMetTime, DispState%sFirstReleaseTime)

  DispState%SyncTime    = DispState%Time
  If (.not.MainOpts%FixedMet) Then
    DispState%MetTime = DispState%Time
  Else
    DispState%MetTime = DispOpts%FixedMetTime
  End If

  ! 7) Allocate arrays.

  DispState%MaxParticles       = DispOpts%MaxParticles
  DispState%MaxFullParticles   = DispOpts%MaxFullParticles
  DispState%MaxPuffs           = DispOpts%MaxPuffs
  DispState%MaxOriginalPuffs   = DispOpts%MaxOriginalPuffs
  DispState%nParticleSpecieses = Specieses%nParticleSpecieses

  If (.not.DispState%SpaceAllocated) Then

    Allocate(DispState%FreeParticleStack     (                                DispState%MaxParticles    ))
    Allocate(DispState%FreeParticleExtraStack(                                DispState%MaxFullParticles))
    Allocate(DispState%iParticleExtras       (                                DispState%MaxParticles    ))
    Allocate(DispState%Particles             (                                DispState%MaxParticles    ))
    Allocate(DispState%ParticleMasses        (DispState%nParticleSpecieses,   DispState%MaxParticles    ))
    Allocate(DispState%ParticleExtras        (                              0:DispState%MaxFullParticles))
    Allocate(DispState%FreeOriginalPuffStack (                                DispState%MaxOriginalPuffs))
    Allocate(DispState%FreePuffStack         (                                DispState%MaxPuffs        ))
    Allocate(DispState%iPuffs                (                                DispState%MaxOriginalPuffs))
    Allocate(DispState%Puffs                 (                                DispState%MaxPuffs        ))
    Allocate(DispState%PuffMasses            (DispState%nParticleSpecieses,   DispState%MaxPuffs        ))
    Allocate(DispState%PuffExtras            (                                DispState%MaxPuffs        ))
    Allocate(DispState%TLast                 (                                DispState%MaxOriginalPuffs))
    Allocate(DispState%T                     (                                DispState%MaxOriginalPuffs))
    Allocate(DispState%SigA2Last             (                                DispState%MaxOriginalPuffs))
    Allocate(DispState%SigA2                 (                                DispState%MaxOriginalPuffs))
    Allocate(DispState%SigA2Defined          (                                DispState%MaxOriginalPuffs))
    DispState%SpaceAllocated = .true.

  Else

    If (DispState%MaxParticles /= Size(DispState%Particles)) Then

      DeAllocate(DispState%FreeParticleStack)
      DeAllocate(DispState%iParticleExtras  )
      DeAllocate(DispState%Particles        )
      DeAllocate(DispState%ParticleMasses   )
      Allocate(DispState%FreeParticleStack(                              DispState%MaxParticles))
      Allocate(DispState%iParticleExtras  (                              DispState%MaxParticles))
      Allocate(DispState%Particles        (                              DispState%MaxParticles))
      Allocate(DispState%ParticleMasses   (DispState%nParticleSpecieses, DispState%MaxParticles))

    End If

    If (DispState%MaxFullParticles /= Size(DispState%FreeParticleExtraStack)) Then

      DeAllocate(DispState%FreeParticleExtraStack)
      DeAllocate(DispState%ParticleExtras)
      Allocate(DispState%FreeParticleExtraStack(  DispState%MaxFullParticles))
      Allocate(DispState%ParticleExtras        (0:DispState%MaxFullParticles))

    End If

    If (DispState%MaxPuffs /= Size(DispState%Puffs)) Then

      DeAllocate(DispState%FreePuffStack)
      DeAllocate(DispState%Puffs        )
      DeAllocate(DispState%PuffMasses   )
      DeAllocate(DispState%PuffExtras   )
      Allocate(DispState%FreePuffStack(                              DispState%MaxPuffs))
      Allocate(DispState%Puffs        (                              DispState%MaxPuffs))
      Allocate(DispState%PuffMasses   (DispState%nParticleSpecieses, DispState%MaxPuffs))
      Allocate(DispState%PuffExtras   (                              DispState%MaxPuffs))

    End If

    If (DispState%MaxOriginalPuffs /= Size(DispState%SigA2)) Then

      DeAllocate(DispState%FreeOriginalPuffStack)
      DeAllocate(DispState%iPuffs               )
      DeAllocate(DispState%TLast                )
      DeAllocate(DispState%T                    )
      DeAllocate(DispState%SigA2Last            )
      DeAllocate(DispState%SigA2                )
      DeAllocate(DispState%SigA2Defined)
      Allocate(DispState%FreeOriginalPuffStack(DispState%MaxOriginalPuffs))
      Allocate(DispState%iPuffs               (DispState%MaxOriginalPuffs))
      Allocate(DispState%TLast                (DispState%MaxOriginalPuffs))
      Allocate(DispState%T                    (DispState%MaxOriginalPuffs))
      Allocate(DispState%SigA2Last            (DispState%MaxOriginalPuffs))
      Allocate(DispState%SigA2                (DispState%MaxOriginalPuffs))
      Allocate(DispState%SigA2Defined         (DispState%MaxOriginalPuffs))

    End If

  End If

  ! 8) Reset various values in DispState.

  DispState%nParticles   = 0
  DispState%LastParticle = 0
  ! Implied do-loop slow to compile when MaxParticles is large.
  !DispState%FreeParticleStack  = (/ (i, i = 1, MaxParticles) /)
  Do i = 1, DispState%MaxParticles
    DispState%FreeParticleStack(i) = i
  End Do

  DispState%nParticleExtras   = 0
  DispState%LastParticleExtra = 0
  Do i = 1, DispState%MaxFullParticles
    DispState%FreeParticleExtraStack(i) = i
  End Do
  DispState%ParticleExtras(0)%FMass   = 0.0     ! Set up Extra(0) for cheap particles.
  DispState%ParticleExtras(0)%VelMem  = .false.
  DispState%ParticleExtras(0)%Inhomog = .false.
  DispState%ParticleExtras(0)%MVelMem = .false.
  DispState%ParticleExtras(0)%TPlus   = ZeroShortTime()
  DispState%ParticleExtras(0)%TMinus  = ZeroShortTime()
  DispState%ParticleExtras(0)%TInst   = InfPastShortTime(Interval = .true.)

  DispState%ParticleExtras(0)%Skew    = .false.
  DispState%ParticleExtras(0)%U(:)    = 0.0
  DispState%ParticleExtras(0)%UM(:)   = 0.0
  DispState%ParticleExtras(0)%FM(:)   = 0.0
  DispState%ParticleExtras(0)%FH      = 0.0
  DispState%ParticleExtras(0)%FMass0  = 0.0
  DispState%ParticleExtras(0)%B0Old   = 0.0

  DispState%nOriginalPuffs   = 0
  DispState%LastOriginalPuff = 0
  Do i = 1, Size(DispState%FreeOriginalPuffStack)
    DispState%FreeOriginalPuffStack(i) = i
    DispState%iPuffs               (i) = 0
  End Do
  DispState%nPuffs           = 0
  DispState%LastPuff         = 0
  Do i = 1, Size(DispState%FreePuffStack)
    DispState%FreePuffStack(i) = i
  End Do

  DispState%nParticleTimeSteps = 0
  DispState%nPuffTimeSteps     = 0
  DispState%MaxTSpread         = ZeroShortTime()
  DispState%SigA2Last(:)       = 0
  DispState%SigA2(:)           = 0
  DispState%SigA2Defined(:)    = 0

  DispState%iULastParticle     = 0
  DispState%iULastPuff         = 0
  DispState%iULastOriginalPuff = 0

  ! Note m agl always available cos added by flows.
  DispState%iZCoordMagl = FindZCoordIndex('m agl', Coords)

  ! Note Lat_Long always available
  DispState%iHCoordLatLong = FindHCoordIndex('Lat-Long', Coords)

  ! MaxPPTimeSpread. - replace by lead and lag

  If (MainOpts%FixedMet) Then
    MaxPPTimeSpread = InfFutureTime(Interval = .true.)
  Else
    MaxPPTimeSpread = DispOpts%PuffInterval * 2
  End If ! $$ Could store in Dispstate and check still valid throughout run
         ! $$ Need to maximise over cases.

  DispState%ParticleFactorType = 0

  DispState%PPsLostDueToFlow = .false.

  ! 9) Check that plume rise is not requested in combination with 2-d trajectories
  ! $$ Is this a sensible place for this test?
  If (.not.DispOpts%VerticalVelocity) Then
    Do i = 1, Sources%nSources
      If (Sources%Sources(i)%PlumeRise) Then
        Call Message('FATAL ERROR: plume rise cannot be used when Vertical Velocity? is false', 3)
      End If
    End Do
  End If

  ! Prepare Results to store the results of the case.
  Call PrepareResults(                                                                          &
         Time2ShortTime(DispState%Time), TMax(DispState%sLastReqTime, DispState%sLastDispTime), &
         Time2ShortTime(DispOpts%SyncdT),                                                       &
         Time2ShortTime(DispOpts%PuffInterval), Time2ShortTime(DispOpts%PuffInterval), & ! $$ correct
                                                                                         ! for fixed met
         Domains%Domains(DispState%iDomain),                                           &
         Grids, Reqs,                                                                  &
         Results                                                                       &
       )

  If (IsInfFuture(DispState%Time)) Then ! $$ move to where EndOfCase is set? (currently Main, but should
                                        ! probably be moved). Currently this message appears
                                        ! before 'Case n started'
    Call Message(                                                                                       &
           'WARNING: There is nothing to compute for this case, possibly because sources or output ' // &
           'requirements lie outside the computational domain',                                         &
           1                                                                                            &
         )
  End If

! Set up Eulerian field
  If ( Specieses%nFields > 0 ) Call SetUpEulerianArrays(DispState%EulerianField, &
                                                        Specieses,               &
                                                        DispOpts%Turbulence,     &
                                                        DispOpts%WetDep,         & 
                                                        DispOpts%AgentDecay)

End Subroutine PrepareForNextCase

!-------------------------------------------------------------------------------------------------------------

Subroutine StartCase(                &
             iCase,                  &
             Coords, Grids, Domains, &
             Specieses,              &
             CloudGammaParamses,     &
             Sources,                &
             Reqs,                   &
             MaterialUnits,          &
             OutputOpts,             &
             SizeDists,              &
             OpenMPOpts,             &
             Units,                  &
             Mets, Flows,            &
             DispState,              &
             Results                 &
           )
! This routine carries out the loop over successive updates of the met and flow
! module instances, calling LoopOverSyncTimes once between each update time.

  Implicit None
  ! Argument list:
  Integer,                   Intent(In)    :: iCase
  Type(Coords_),             Intent(In)    :: Coords
  Type(Grids_),              Intent(In)    :: Grids
  Type(Domains_),            Intent(In)    :: Domains
  Type(Specieses_),          Intent(In)    :: Specieses
  Type(CloudGammaParamses_), Intent(In)    :: CloudGammaParamses
  Type(Sources_),            Intent(In)    :: Sources
  Type(Reqs_),               Intent(In)    :: Reqs
  Type(MaterialUnits_),      Intent(In)    :: MaterialUnits
  Type(OutputOpts_),         Intent(In)    :: OutputOpts
  Type(SizeDists_),          Intent(In)    :: SizeDists
  Type(OpenMPOpts_),         Intent(In)    :: OpenMPOpts
  Type(Units_),              Intent(InOut) :: Units
  Type(Mets_),               Intent(InOut) :: Mets
  Type(Flows_),              Intent(InOut) :: Flows
  Type(DispState_),          Intent(InOut) :: DispState
  Type(Results_),            Intent(InOut) :: Results
  ! iCase              :: Number of case.
  ! Coords             :: Collection of coord systems.
  ! Grids              :: Collection of grids.
  ! Domains            :: Collection of domains.
  ! Specieses          :: Collection of species.
  ! CloudGammaParamses :: Collection of sets of cloud gamma parameters.
  ! Sources            :: Collection of sources.
  ! Reqs               :: Collection of requirements.
  ! MaterialUnits      :: Collection of material units
  ! OutputOpts         :: Output options.
  ! OpenMPOpts         :: OpenMP options.
  ! Units              :: Collection of information on input/output unit numbers.
  ! Mets               :: Collection of met module instance states.
  ! Flows              :: Set of flow module instance states.
  ! DispState          :: Dispersion model state.
  ! Results            :: Results.
  ! Locals:
  Type(ReqInfo_) :: ReqInfo(1)

  Call CalcReqInfoNoT(Reqs, ReqInfo(1))

  ! Calculate type O (other) fields. $$ calc runtotime and pass to this routine
  ! to use instead of inffuturetime
  Call CalcOtherResults(ReqInfo, InfFutureTime(), Coords, Grids, Domains, Reqs, DispState, Results)

  Call TimerOn(OutputTimer)
  ! Process and output results. Note this call is not necessary - results calculated to date would be
  ! processed and output anyway at the next call to ProcessAndOutputResults. However this call will make
  ! results available sooner.
  Call ProcessAndOutputResults(                                                                   &
         iCase, .false., .false.,                                                                 &
         InfPastShortTime(), .true., .false.,                                                     &
            ! $$ review
         ZeroShortTime(),                                                                         &
         .false., .false.,                                                                        &
         TMin(DispState%sFirstReqTime, DispState%sFirstReleaseTime),                              &
         TMax(TMax(DispState%sLastDispTime, DispState%sLastReqTime), DispState%sLastReleaseTime), &
         OutputOpts, SizeDists, OpenMPOpts,                                                       &
         Coords, Grids, Domains,                                                                  &
         Specieses, CloudGammaParamses, Sources, Reqs, MaterialUnits,                             &
         Units, Mets, Flows, Results                                                              &
       )
  Call TimerOff(OutputTimer)
  Call TimerWriteLast(OutputTimer)

End Subroutine StartCase

!-------------------------------------------------------------------------------------------------------------

Subroutine RunToRestartDumpOrSuspendOrEndOfCase(      &
             NextDumpTime, RunToTime,                 &
             iCase,                                   &
             Coords, Grids, Domains,                  &
             Specieses, CloudGammaParamses,           &
             SizeDists, Sources, TimeDeps,            &
             Reqs, MaterialUnits,                     &
             MainOpts, OutputOpts, DispOpts,          &
             Units,                                   &
             Mets, Flows,                             &
             Results,                                 &
             DispState,                               &
             ChemOpts, ChemistryDefn, ChemistryState, &
             OpenMPOpts,                              &
             Suspend, EndOfCase, EndOfCases           &
           )
! This routine carries out the loop over successive updates of the met and flow
! module instances, calling LoopOverSyncTimes once between each update time.

  Implicit None
  ! Argument list:
  Type(Time_),               Intent(In)    :: NextDumpTime ! Exit at next break point
                                                           ! (sync time or met update) after this time.
                                                           ! (rename variable to reflect suspend options too)
  Type(Time_),               Intent(In)    :: RunToTime    ! Suspend run at this point.
  Integer,                   Intent(In)    :: iCase
  Type(Coords_),             Intent(In)    :: Coords
  Type(Grids_),              Intent(In)    :: Grids
  Type(Domains_),            Intent(In)    :: Domains
  Type(Specieses_),          Intent(In)    :: Specieses
  Type(CloudGammaParamses_), Intent(In)    :: CloudGammaParamses
  Type(SizeDists_),          Intent(In)    :: SizeDists
  Type(Sources_),            Intent(In)    :: Sources
  Type(TimeDeps_),           Intent(In)    :: TimeDeps
  Type(Reqs_),               Intent(In)    :: Reqs
  Type(MaterialUnits_),      Intent(In)    :: MaterialUnits
  Type(MainOpts_),           Intent(In)    :: MainOpts
  Type(OutputOpts_),         Intent(In)    :: OutputOpts
  Type(DispOpts_),           Intent(In)    :: DispOpts
  Type(Units_),              Intent(InOut) :: Units
  Type(Mets_),               Intent(InOut) :: Mets
  Type(Flows_),              Intent(InOut) :: Flows
  Type(Results_),            Intent(InOut) :: Results
  Type(DispState_),          Intent(InOut) :: DispState
  Type(ChemOpts_),           Intent(In)    :: ChemOpts
  Type(ChemistryDefn_),      Intent(In)    :: ChemistryDefn
  Type(ChemistryState_),     Intent(InOut) :: ChemistryState
  Type(OpenMPOpts_),         Intent(In)    :: OpenMPOpts
  Logical,                   Intent(InOut) :: Suspend
  Logical,                   Intent(InOut) :: EndOfCase
  Logical,                   Intent(InOut) :: EndOfCases
  ! iCase              :: Number of case.
  ! Coords             :: Collection of coord systems.
  ! Grids              :: Collection of grids.
  ! Domains            :: Collection of domains.
  ! Specieses          :: Collection of species.
  ! CloudGammaParamses :: Collection of sets of cloud gamma parameters.
  ! SizeDists          :: Collection of particle size distributions.
  ! Sources            :: Collection of sources.
  ! TimeDeps           :: Collection of source time dependencies.
  ! Reqs               :: Collection of requirements.
  ! MaterialUnits      :: Collection of material units
  ! MainOpts           :: Main options.
  ! DispOpts           :: Set of dispersion model options.
  ! Units              :: Collection of information on input/output unit numbers.
  ! Mets               :: Set of met module instance states.
  ! Flows              :: Set of flow module instance states.
  ! Results            :: Results.
  ! DispState          :: Dispersion model state.
  ! ChemOpts           :: Chemistry options.
  ! ChemistryDefn      :: Information defining a chemistry scheme.
  ! ChemistryState     :: State of chemistry calculation.
  ! OpenMPOpts         :: OpenMP options
  ! TimerOpts          :: Timer module options
  ! EndOfCase          :: Indicates that there is no more to do for this case.
  ! EndOfCases         :: Indicates that there are no more cases to be run.
  ! Locals:
  Type(Time_) :: TargetTime             ! Target time to which LoopOverSyncTimes should progress.
  Type(Time_) :: NextMetFlowsUpdateTime ! Target time to which LoopOverSyncTimes should progress.

  NextMetFlowsUpdateTime = ShortTime2Time(DispState%sFirstMetTime)

  ! Loop over successive updates of the met and flow module instances.
  Do While (.not.(Suspend .or. EndOfCase .or. EndOfCases .or. DispState%Time >= NextDumpTime))

    ! Update MetTime.
    If (.not.MainOpts%FixedMet) DispState%MetTime = DispState%Time

    ! Update met and flow module instances and calculate target time.
    If (NextMetFlowsUpdateTime <= DispState%Time) Then
    !  $$ This if test has been commented out and should probably be removed eventually.
    !     If restored, will need modifying to cope with Eulerian advection scheme. 
    !     The effect of removing it is that met is always read (or prepared for reading if update-on-demand
    !     is used) for as long as the run continues. This is an efficiency loss for runs which have potential
    !     to be long but lose all their particles early, but using update-on-demand would probably be a 
    !     satisfactory solution in such cases. The reduced complexity benefits (this code has caused problems 
    !     in the past) probably out-weigh the efficiency loss. Some consequential simplifications may well be 
    !     possible and should be considered - search for EndOfCase and PPsFinished for ideas. In particular
    !     the use of PPsFinished in the if test around the call to ProcessAndOutputResults seems complex and
    !     probably unnecessary (??) and PPsFinsihed could just be a local variable in that routine.
    !  Note the third "or" clause is for runs with no particles and with requests for met output at a
    !  single time only. Without the clause the met would not be read in for this case.
    !  If (                                                           &
    !    DispState%Time < ShortTime2Time(DispState%sLastMetTime) .or. &
    !    .not.DispState%PPsFinished                              .or. &
    !    (                                                            &
    !      TMin(DispState%sFirstReqTime, DispState%sFirstReleaseTime) == Time2ShortTime(DispState%SyncTime) &
    !      .and. DispState%Time == ShortTime2Time(DispState%sLastMetTime)                                   &
    !    )                                                                                                  &
    !  ) Then
        Call UpdateMetsFlows(             &
               Coords, Grids, Domains,    &
               iCase, DispState%iMetCase, &
               DispState%MetTime,         &
               Mets, Flows,               &
               Units                      &
             )
        NextMetFlowsUpdateTime = MetsFlowsOverallTValid(Flows)
    !  Else
    !    NextMetFlowsUpdateTime = InfFutureTime()
    !  End If
    End If

    ! Calculate target time for LoopParticlesPuffsAndTimeSteps.
    TargetTime = TMin(NextMetFlowsUpdateTime, DispState%SyncTime + DispOpts%SyncdT)

    ! Set up chemistry state after the first met and flow update step (for chemistry runs only).
    If (DispOpts%Chemistry .and. .not.ChemistryState%Initialised) Then
      Call InitChemistryState(                       &
             DispState%MaxParticles,                 &
             Coords, Grids, Domains, DispState%Time, &
             ChemOpts, ChemistryDefn,                &
             Units, Mets, Flows,                     &
             DispState%EulerianField,                &
             ChemistryState                          &
           )
    End If

    ! Loop over successive times when the particle/puffs are synchronised.
    Call RunToSyncTimeOrMetFlowUpdateOrEndOfCase(   &
           iCase,                                   &
           TargetTime, RunToTime,                   &
           Coords, Grids, Domains,                  &
           Specieses, CloudGammaParamses,           &
           SizeDists, Sources, TimeDeps,            &
           Mets, Flows, Reqs, MaterialUnits,        &
           MainOpts, OutputOpts,                    &
           DispOpts,                                &
           Units, Results, DispState,               &
           ChemOpts, ChemistryDefn, ChemistryState, &
           OpenMPOpts,                              &
           Suspend, EndOfCase, EndOfCases           &
         )

  End Do

End Subroutine RunToRestartDumpOrSuspendOrEndOfCase

!-------------------------------------------------------------------------------------------------------------

Subroutine RunToSyncTimeOrMetFlowUpdateOrEndOfCase(   &
             iCase,                                   &
             TargetTime, RunToTime,                   &
             Coords, Grids, Domains,                  &
             Specieses, CloudGammaParamses,           &
             SizeDists, Sources, TimeDeps,            &
             Mets, Flows, Reqs, MaterialUnits,        &
             MainOpts, OutputOpts, DispOpts,          &
             Units, Results, DispState,               &
             ChemOpts, ChemistryDefn, ChemistryState, &
             OpenMPOpts,                              &
             Suspend, EndOfCase, EndOfCases           &
           )
! This routine carries out the loop over successive times when the particle/puffs are synchronised, calling
! LoopParticlesPuffsAndTimeSteps once between each synchronisation time.

  Implicit None

  ! Argument list:
  Integer,                   Intent(In)            :: iCase
  Type(Time_),               Intent(In)            :: TargetTime
  Type(Time_),               Intent(In)            :: RunToTime
  Type(Coords_),             Intent(In)            :: Coords
  Type(Grids_),              Intent(In)            :: Grids
  Type(Domains_),            Intent(In)            :: Domains
  Type(Specieses_),          Intent(In)            :: Specieses
  Type(CloudGammaParamses_), Intent(In)            :: CloudGammaParamses
  Type(SizeDists_),          Intent(In),    Target :: SizeDists
  Type(Sources_),            Intent(In),    Target :: Sources
  Type(TimeDeps_),           Intent(In)            :: TimeDeps
  Type(Mets_),               Intent(InOut)         :: Mets
  Type(Flows_),              Intent(InOut)         :: Flows
  Type(Reqs_),               Intent(In)            :: Reqs
  Type(MaterialUnits_),      Intent(In)            :: MaterialUnits
  Type(MainOpts_),           Intent(In)            :: MainOpts
  Type(OutputOpts_),         Intent(In)            :: OutputOpts
  Type(DispOpts_),           Intent(In)            :: DispOpts
  Type(Units_),              Intent(InOut)         :: Units
  Type(Results_),            Intent(InOut)         :: Results
  Type(DispState_),          Intent(InOut), Target :: DispState
  Type(ChemOpts_),           Intent(In)            :: ChemOpts
  Type(ChemistryDefn_),      Intent(In)            :: ChemistryDefn
  Type(ChemistryState_),     Intent(InOut)         :: ChemistryState
  Type(OpenMPOpts_),         Intent(In)            :: OpenMPOpts
  Logical,                   Intent(InOut)         :: Suspend
  Logical,                   Intent(InOut)         :: EndOfCase
  Logical,                   Intent(InOut)         :: EndOfCases
  ! iCase              :: Number of case.
  ! TargetTime         :: Target time to which LoopOverSyncTimes should progress.
  ! Coords             :: Collection of coord systems.
  ! Grids              :: Collection of grids.
  ! Domains            :: Collection of domains.
  ! Specieses          :: Collection of species.
  ! CloudGammaParamses :: Collection of sets of cloud gamma parameters.
  ! SizeDists          :: Collection of particle size distributions.
  ! Sources            :: Collection of sources.
  ! TimeDeps           :: Collection of source time dependencies.
  ! Mets               :: Set of met module instance states.
  ! Flows              :: Set of flow module instance states.
  ! Reqs               :: Collection of requirements.
  ! MainOpts           :: Main options.
  ! DispOpts           :: Set of dispersion model options.
  ! Units              :: Collection of information on input/output unit numbers.
  ! Results            :: Results.
  ! DispState          :: Dispersion model state.
  ! ChemOpts           :: Chemistry options.
  ! ChemistryDefn      :: Information defining a chemistry scheme.
  ! ChemistryState     :: State of chemistry calculation.
  ! OpenMPOpts         :: OpenMP options
  ! EndOfCase          :: Indicates that there is no more to do for this case.
  ! EndOfCases         :: Indicates that there are no more cases to be run.
  ! Locals:
  ! Type(Time_)                      :: TargetTime1
  Type(ReqInfo_)                   :: ReqInfoL(MaxReqInfoTimes)
  Type(ReqInfo_)                   :: ReqInfoE(MaxReqInfoTimes)
  Type(ReqInfo_)                   :: ReqInfoF(MaxReqInfoTimes)
  Type(ReqInfo_)                   :: ReqInfoO(MaxReqInfoTimes)
  Type(ReqInfo_)                   :: ReqInfoR(MaxReqInfoTimes)
  Integer                          :: i
  Integer                          :: iReqInfo
  Type(ShortTime_)                 :: OldSyncTime
  Type(Flow_)                      :: Flow
  Type(Cloud_)                     :: Cloud
  Type(Rain_)                      :: Rain
  Logical                          :: PuffCentre !
  Integer                          :: j
  Type(Puff_),             Pointer :: Puff
  Type(Particle_),         Pointer :: Particle
  Type(Extra_),            Pointer :: Extra
  Type(FlowMemory_)                :: FlowMemory
  Type(Position_)                  :: Position
  Type(Source_),           Pointer :: Source
  Integer                          :: iP
  Integer                          :: iE
  Integer                          :: iUP
  Logical                          :: PPsFinished ! Indicates all particles/puffs are finished.
  Logical                          :: Exist
  Integer                          :: RunToUnit
  Character(MaxCharLength)         :: CharRunToTimeFromFile
  Type(Time_)                      :: RunToTimeFromFile
  Integer                          :: Char2TimeErrorCode
  Integer                          :: IOStat
  Integer                          :: nParticles
  Real(Std)                        :: HMax, H(3)
  Real(Std)                        :: TimeStep
  Integer                          :: iOld
  Integer                          :: iNew
  Type(EulerianField_),    Pointer :: EulerianField
  Integer                          :: iMetMod
  Real(Std)                        :: Time0
  Real(Std)                        :: Time1
  Real(Std)                        :: TimeOfInterest
  Integer                          :: iSpecies
  Integer                          :: iSpeciesUse
  Integer                          :: iParticleSize
  Integer                          :: iTracer
  Type(SizeDist_),         Pointer :: SizeDist
  Integer                          :: nX
  Integer                          :: nY
  Integer                          :: ix
  Integer                          :: iy
  Integer                          :: iz
  Real(Std)                        :: DummyConcentration(Specieses%nFields)
  Integer                          :: iFirst
  Integer                          :: iLast

  ! TargetTime1 :: Target time to which LoopParticlesPuffsAndTimeSteps should progress.
  ! ReqInfo     :: Information on what requirements are required at the next sync time.
  ! i           :: Loop index.
  ! iReqInfo    ::
  
  ! iFirst      ::} Shorthand references for the first and last advected non-sedimenting fields
  ! iLast       ::}

  ! DryDep and WetDep output not supported at sync times. Need to add support or check
  ! to ensure not requested. $$ Also PPInfo output involving flow, cloud or rain

  ! Ensure message about particle/puff loss generated every sync step
  DispState%PPsLostDueToFlow = .false.

  ! Determine output requirements for source strength (if any) 
  If ( Any(Sources%Sources(:)%MetDependent == 2 ) ) Then  !! necessary?
    Call CalcReqInfoT(                                           &
           .true., .false., 'A', .true., .false.,                & ! $$ use IncTime1 at start of run to
           Time2ShortTime(DispState%SyncTime) ,                  & ! include Time1 in range $$
           Time2ShortTime(DispState%SyncTime + DispOpts%SyncdT), & ! $$ should this be TargetTime ?
           ZeroShortTime(), ZeroShortTime(),                     &
           Grids, Reqs,                                          &
           ReqInfoR                                              &
         )
  End If

  ! Initialise particle/puffs.
  Call Release(                                                                        &
         Time2ShortTime(DispState%Time), Time2ShortTime(TargetTime),                   &
         Coords, Grids, Domains, Mets, Flows, Specieses, Sources, SizeDists, TimeDeps, &
         DispOpts,                                                                     &
         Units,                                                                        &
         MainOpts%FixedMet,                                                            &
         DispState,                                                                    &
         ReqInfoR,                                                                     &
         Reqs,                                                                         &
         Results                                                                       &
       )

  ! Evolve particle puffs.
  Call TimerOn(LPPATsTimer)
  !$OMP PARALLEL NUM_THREADS(OpenMPOpts%nParticleThreads)
  Call LoopParticlesPuffsAndTimeSteps(                                       &
         iCase,                                                              &
         Time2ShortTime(DispState%Time),                                     &
         Time2ShortTime(TargetTime),                                         &
         OutputOpts,                                                         &
         Coords, Grids, Domains, Mets, Flows,                                &
         Specieses, CloudGammaParamses, Sources, SizeDists, Reqs,            &
         DispOpts, Units, Results, DispState                                 &
       )
  !$OMP END PARALLEL
  Call TimerOff(LPPATsTimer)
  Call TimerWriteLast(LPPATsTimer)

  ! Call Eulerian advection model
  If (Specieses%nFields > 0) Then

    EulerianField => DispState%EulerianField

    ! Read flow field 
    ! $$ note this will fill halos for static fields - OK because halo's of such fields not used.
    !    Could treat halos in separate routine

    If ( EulerianField%InitialiseFields .And. Specieses%AdvectedFields ) Then
      Call EulerianFlowField(EulerianField,       &
                             Specieses,           &
                             SizeDists,           &
                             DispState%Time,      &
                             Coords,              &
                             Grids,               &
                             Domains,             &
                             Mets,                &
                             Flows,               &
                             FlowMemory,          &
                             Units,               &
                             DispOpts%Turbulence, &
                             DispOpts%DryDep,     &
                             DispOpts%WetDep,     &
                             DispOpts%AgentDecay)
    End If
    EulerianField%iOld = EulerianField%iNew
    EulerianField%iNew = 3 - EulerianField%iNew
    If (Specieses%AdvectedFields) Then
      Call EulerianFlowField(EulerianField,       &
                             Specieses,           &
                             SizeDists,           &
                             TargetTime,          &
                             Coords,              &
                             Grids,               &
                             Domains,             &
                             Mets,                &
                             Flows,               &
                             FlowMemory,          &
                             Units,               &
                             DispOpts%Turbulence, &
                             DispOpts%DryDep,     &
                             DispOpts%WetDep,     &
                             DispOpts%AgentDecay)
    End If

    ! Advection
    ! Call semi-Lagrangian solver: advection scheme
    ! Call sl_advect(fld_n, fld_np1, rho_n, rho_np1, source_np1,  &
    !                u_n, v_n, w_n, u_np1, v_np1, w_np1, timestep )

    iOld = EulerianField%iOld
    iNew = EulerianField%iNew
    nX = EulerianField%HGrid%nX
    nY = EulerianField%HGrid%nY

    TimeStep = ShortTime2RealTime(Time2ShortTime(TargetTime) - Time2ShortTime(DispState%SyncTime))

    ! Non-sedimenting species
    iFirst = Specieses%iFirstAdvectedNonSedimentingField
    iLast  = Specieses%iLastAdvectedNonSedimentingField
    If (iLast >= iFirst) Then
      Call sl_advect(EulerianField%Concentration(:, :, :, iFirst:iLast, iOld), &
                     EulerianField%Concentration(:, :, :, iFirst:iLast, iNew), &
                     EulerianField%Density(:, :, :, iOld),                     &
                     EulerianField%Density(:, :, :, iNew),                     &
                     EulerianField%NewSource(:, :, :, iFirst:iLast),           &
                     EulerianField%U(:, :, :, iOld),                           &
                     EulerianField%V(:, :, :, iOld),                           &
                     EulerianField%W(:, :, :, iOld),                           &
                     EulerianField%U(:, :, :, iNew),                           &
                     EulerianField%V(:, :, :, iNew),                           &
                     EulerianField%W(:, :, :, iNew),                           &
                     TimeStep,                                                 &
                     iLast - iFirst + 1)
    End If

    ! Sedimenting particles
    If (Specieses%AdvectedSedimentingFields) Then

      Do iTracer = 1, Specieses%nFields
        iSpecies = Specieses%iField2Species(iTracer)
        iSpeciesUse = Specieses%iSpecies2SpeciesUses(iSpecies)
        If (.not. Specieses%SpeciesUseses(iSpeciesUse)%AdvectField) Cycle
        If (iFirst <= iTracer .and. iTracer <= iLast) Cycle
        iParticleSize = Specieses%iField2Size(iTracer)
        SizeDist => SizeDists%SizeDists(Specieses%SpeciesUseses(iSpeciesUse)%iSizeDist)
        Call EulerianSedimentationVelocity(                    &
               iTracer,                                        &
               SizeDist%Density(iParticleSize),                &
               SizeDist%RepresentativeDiameter(iParticleSize), &
               EulerianField                                   &
             )
        Call sl_advect(EulerianField%Concentration(:, :, :, iTracer:iTracer, iOld), &
                       EulerianField%Concentration(:, :, :, iTracer:iTracer, iNew), &
                       EulerianField%Density(:, :, :, iOld),                        &
                       EulerianField%Density(:, :, : ,iNew),                        &
                       EulerianField%NewSource(:, :, :, iTracer:iTracer),           &
                       EulerianField%U(:, :, :, iOld),                              &
                       EulerianField%V(:, :, :, iOld),                              &
                       EulerianField%ModifiedW(:, :, :, iOld),                      &
                       EulerianField%U(:, :, :, iNew),                              &
                       EulerianField%V(:, :, :, iNew),                              &
                       EulerianField%ModifiedW(:, :, :, iNew),                      &
                       TimeStep,                                                    &
                       1)
        EulerianField%TotalDepositionRate(:, :, iTracer) =  0.5 * (      &
            EulerianField%Concentration(1:nX, 1:nY, 1, iTracer, iOld) *  &
              EulerianField%WSedAtGround(:, :, iTracer, iOld) +          &
            EulerianField%Concentration(1:nX, 1:nY, 1, iTracer, iNew) *  &
              EulerianField%WSedAtGround(:, :, iTracer, iNew)            &
          )
      End Do
    End If

    ! Non-advected fields.
    ! $$ could be more efficient with an array of iOld's and iNew's - one for each field.
    Do iTracer = 1, Specieses%nFields
      iSpecies = Specieses%iField2Species(iTracer) 
      iSpeciesUse = Specieses%iSpecies2SpeciesUses(iSpecies)
      If (.not. Specieses%SpeciesUseses(iSpeciesUse)%AdvectField) Then
        EulerianField%Concentration(:, :, :, iTracer, iNew) =  &
          EulerianField%Concentration(:, :, :, iTracer, iOld)
      End If
    End Do
    
    Do iTracer = 1, Specieses%nFields
      iSpecies = Specieses%iField2Species(iTracer) 
      iSpeciesUse = Specieses%iSpecies2SpeciesUses(iSpecies)
      If (.not. Specieses%SpeciesUseses(iSpeciesUse)%AdvectField) Cycle

      ! Vertical diffusion
      If ( DispOpts%Turbulence ) Then

        ! Update concentration 
        EulerianField%Concentration(:, :, :, iTracer, iOld) =  &
          EulerianField%Concentration(:, :, :, iTracer, iNew)

        ! Calculate diffusion
        Call EulerianVerticalDiffusion(EulerianField%Concentration(:, :, :, iTracer:iTracer, iOld), &
                                       EulerianField%Concentration(:, :, :, iTracer:iTracer, iNew), &
                                       EulerianField%Diffusivity(:, :, :, iOld),                    &
                                       EulerianField%Diffusivity(:, :, :, iNew),                    &
                                       EulerianField%Density(:, :, :, iOld),                        &
                                       EulerianField%Density(:, :, :, iNew),                        &
                                       TimeStep,                                                    &
                                       1,                                                           &
                                       EulerianField%HGrid%nX,                                      &
                                       EulerianField%HGrid%nY,                                      &
                                       EulerianField%ZGrid%nZ,                                      &
                                       EulerianField%ZAboveGround(:, :, :),                         &
                                       EulerianField%Vd(:, :, iTracer:iTracer, iOld),               &
                                       EulerianField%Vd(:, :, iTracer:iTracer, iNew)                &
             )

      End If

      ! Calculate wet and dry deposition 
      Call EulerianDeposition(TimeStep, iTracer, DispOpts%DryDep, DispOpts%WetDep, DispOpts%Turbulence, &
                              EulerianField)

    End Do

    ! Calculate radioactive decay of fields 
    If (DispOpts%RadioactiveDecay) Then
      Do ix = 1, nX
        Do iy = 1, nY
          Do iz = 1, EulerianField%ZGrid%nZ
            DummyConcentration(:) = Real(EulerianField%Concentration(ix, iy, iz, :, iNew), Std)
            Call RadioactiveDecay(DummyConcentration, Specieses, TimeStep, Field = .True.)
            EulerianField%Concentration(ix, iy, iz, :, iNew) = Real(DummyConcentration(:), Con)
          End Do
        End Do
      End Do
    End If

    ! Agent (virus) decay of fields 
    If (DispOpts%AgentDecay) Then
      Call EulerianAgentDecay(TargetTime, TimeStep, EulerianField, Specieses, Coords, Grids, Domains, &
                              Mets, Flows, FlowMemory, Units)
    End If

    ! Calculate contribution to Eulerian concentration fields from each particle
    Call Particles2Fields(Coords, Grids, Domains, Sources, Mets, Flows, FlowMemory, Units, Results, &
                          DispState, SizeDists, Specieses)

  End If

  !--------------------------------!
  ! Any active puffs or particles? !
  !--------------------------------!

  PPsFinished = .true.

  OneTimeLoop: Do

    ! Loop over puffs.
    Do iP = 1, DispState%LastPuff
      If (ParticleReallyActive(DispState%Puffs(iP)%P)) Then
        PPsFinished = .false.
        Exit OneTimeLoop
      End If
    End Do

    ! Loop over particles.
    Do iP = 1, DispState%LastParticle
      If (ParticleReallyActive(DispState%Particles(iP))) Then
        PPsFinished = .false.
        Exit OneTimeLoop
      End If
    End Do

    Exit OneTimeLoop

  End Do OneTimeLoop

  ! Ensure doesn't stop if all particles/puffs killed, but more to release.
  If (                                                                                &
    PPsFinished                                                                 .and. &
    Time2ShortTime(TargetTime) <= DispState%sLastReleaseTime                    .and. &
    Time2ShortTime(TargetTime) - DispState%MaxTSpread < DispState%sLastDispTime       &
  ) Then
    PPsFinished = .false.
  End If

  DispState%PPsFinished = PPsFinished

  ! EndOfCase = PPsFinished ! Now replaced by the following:
  ! If PPsFinished, MaxTSpread will equal zero at end of this routine. Hence can take it as zero here
  ! in determining if run will end.
  EndOfCase = Time2ShortTime(TargetTime) >= DispState%sLastReqTime
  If (PPsFinished) Then
  Else
    EndOfCase = EndOfCase .and. &
                Time2ShortTime(TargetTime) >= DispState%sLastDispTime + DispState%MaxTSpread
  End If

  ! Perform chemistry calculation for the current synchronisation interval
  ! (for chemistry runs only).
  ! NOTE: any output needing chemistry must be requested as Sync output
  ! (otherwise output occurs before the chemistry update).
  !$$ Need to ensure that time step is ok here.
  If (DispOpts%Chemistry) Then
    Call UpdateChemistry(                                                         &
           Coords, Grids, Domains, Units, Mets, Flows, Specieses,                 &
           DispState%Time, TargetTime,                                            &
           ShortTime2RealTime(Time2ShortTime(DispOpts%SyncdT)),                   &
           DispState%LastParticle, DispState%Particles, DispState%ParticleMasses, &
           ChemOpts, ChemistryDefn, ChemistryState, DispState%EulerianField,      &
           OpenMPOpts                                                             &
         )
  End If
  If (DispOpts%Chemistry .or. Specieses%nFields > 0) Then
    ! Calculate requirements for any contributions to chemistry outputs within the synchronisation time step.
    ! $$ this needs tidying up, e.g. moving these routines to more appropriate places below. Also consider
    !    what happens here if the met update time is not a sync time?
    !    Also might want to output zeros if no fields.
    Call CalcReqInfoT(                                           &
           .true., .true., 'E', .true., .false.,                 & ! $$ use IncTime1 at start of run to
           Time2ShortTime(DispState%SyncTime) ,                  & ! include Time1 in range $$
           Time2ShortTime(DispState%SyncTime + DispOpts%SyncdT), & ! $$ should this be TargetTime ?
           ZeroShortTime(), ZeroShortTime(),                     &
           Grids, Reqs,                                          &
           ReqInfoE                                              &
         )
    ! Process output requests for chemistry fields.
    Call CalcEulerianResults(                                  &
           TargetTime, ReqInfoE, Coords, Grids, Domains, Reqs, &
           Units, Mets, Flows, Results,                        &
           ChemOpts, ChemistryDefn, ChemistryState, DispState, &
           SizeDists, Specieses                                &
         )
  End If

  ! Perform splitting of massive particles (when particle splitting is turned on).
  ! $$ Note this subroutine is currently located in the Chemistry module, but could
  ! $$ be invoked in more general applications - consider moving it elsewhere?
  ! If (DispOpts%ParticleSplitting) Then  $$ need to set user input via a DispOpts flag
  If (DispOpts%Chemistry) Then
    Call SplitMassiveParticles(                                                                  &
           Specieses,                                                                            &
           DispState%LastParticle, DispState%nParticles,                                         &
           DispOpts%ParticleCeiling, DispOpts%ParticleFactor, DispOpts%MaxParticles,             &
           DispState%FreeParticleStack, DispState%Particles, DispState%ParticleMasses,           &
           DispState%ParticleFactorType,                                                         &
           DispState%nParticleExtras, DispState%LastParticleExtra,                               &
           DispState%FreeParticleExtraStack, DispState%ParticleExtras, DispState%iParticleExtras &
         )
  End If

  ! Update current time (value may be wrong if end of case, but that doesn't
  ! matter).
  DispState%Time = TargetTime

  ! If synchronisation time reached or end of case, calculate output.
  If ((DispState%Time >= DispState%SyncTime + DispOpts%SyncdT) .or. EndOfCase) Then

    ! Calculate requirements which the particles and puffs within the
    ! synchronisation time step might contribute to.
    Call CalcReqInfoT(                                                                               &
           .true., .true., 'L', .true.,                                                              &
           .false.,                                                                                  &
             ! $$ use IncTime1 at start of run ??
           Time2ShortTime(DispState%SyncTime), Time2ShortTime(DispState%SyncTime + DispOpts%SyncdT), &
           DispState%MaxTSpread, DispState%MaxTSpread,                                               &
           Grids, Reqs,                                                                              &
           ReqInfoL                                                                                  &
         )
    Call CalcReqInfoT(                                                                                       &
           .true., .true., 'F', .true.,                                                                      &
           TMin(DispState%sFirstReqTime, DispState%sFirstReleaseTime) == Time2ShortTime(DispState%SyncTime), &
           ! IncTime1 used at start of run to include Time1 in range
           Time2ShortTime(DispState%SyncTime), Time2ShortTime(DispState%SyncTime + DispOpts%SyncdT),         &
           DispState%MaxTSpread, DispState%MaxTSpread,                                                       &
           Grids, Reqs,                                                                                      &
           ReqInfoF                                                                                          &
         )
    Call CalcReqInfoT(                                                                               &
           .true., .true., 'O', .true.,                                                              &
           .false.,                                                                                  &
           ! $$ use IncTime1 at start of run to include Time1 in range $$
           Time2ShortTime(DispState%SyncTime), Time2ShortTime(DispState%SyncTime + DispOpts%SyncdT), &
           DispState%MaxTSpread, DispState%MaxTSpread,                                               &
           Grids, Reqs,                                                                              &
           ReqInfoO                                                                                  &
         )

    ! Update synchronisation time (value may be wrong if end of case, but that
    ! doesn't matter).
    OldSyncTime = Time2ShortTime(DispState%SyncTime)
    DispState%SyncTime = DispState%SyncTime + DispOpts%SyncdT

    ! Calculate contribution to output. $$ need to test T > TOld to avoid double
    !                                      counting for particles killed in
    !                                      OneParticle(or Puff)TimeStep ??
    Do i = 1, DispState%LastPuff

      Puff     => DispState%Puffs(i)
      Particle => DispState%Puffs(i)%P
      Extra    => DispState%PuffExtras(i)

      If (.not.ParticleActive(Puff%P)) Cycle

      If (UseParticleForOutput(Particle)) Then

        iReqInfo = 1
        Do
          If (ReqInfoL(iReqInfo)%NoReqs) Exit
          If (ParticleInstTime(Puff%P, Extra) <= ReqInfoL(iReqInfo)%Time) Then
            If (OldSyncTime < ReqInfoL(iReqInfo)%Time) Exit
          Else
            If (OldSyncTime + Extra%TMinus < ReqInfoL(iReqInfo)%Time) Exit
          End If
          iReqInfo = iReqInfo + 1
        End Do

        Call ResetFlowMemory(Flows, FlowMemory)
        Position = X2Position(Coords, Puff%P%X, Puff%P%iHCoord, Puff%P%iZCoord)

        Do
          If (ReqInfoL(iReqInfo)%NoReqs) Exit
          If (ParticleInstTime(Puff%P, Extra) <= ReqInfoL(iReqInfo)%Time) Then
            If (Time2ShortTime(TargetTime) < ReqInfoL(iReqInfo)%Time) Exit
          Else
            If (                                                                 &
              Time2ShortTime(TargetTime) + Extra%TPlus < ReqInfoL(iReqInfo)%Time &
            ) Exit
          End If
          PuffCentre =                                            &
            OldSyncTime < ReqInfoL(iReqInfo)%Time .and.           &
            Time2ShortTime(TargetTime) >= ReqInfoL(iReqInfo)%Time
          Call CalcPuffResults(                                                                      &
                 iCase,                                                                              &
                 Puff, DispState%PuffMasses(:, i), Extra,                                            &
                 (/ (0.0, j = 1, MaxSpecieses) /),                                                   & ! temp
                                                                                                       ! fix
                 (/ (0.0, j = 1, MaxSpecieses) /),                                                   & ! temp
                                                                                                       ! fix
                 ShortTime2RealTime(Time2ShortTime(DispOpts%SyncdT)),                                &
                 PuffCentre,                                                                         &
                 ReqInfoL(iReqInfo),                                                                 &
                 OutputOpts,                                                                         &
                 Coords, Grids, Domains, Mets, Flows, Specieses, Sources, SizeDists, Reqs, DispOpts, &
                 Position, FlowMemory,                                                               &
                 Units, Results, DispState                                                           &
               )
          iReqInfo = iReqInfo + 1
        End Do

      End If

      ! Kill puff.
      If (.not.ParticleReallyActive(Particle)) Then
        Call InactivatePuff(                                                                          &
               i,                                                                                     &
               DispState%Puffs,                                                                       &
               DispState%nPuffs, DispState%LastPuff, DispState%FreePuffStack,                         &
               DispState%nOriginalPuffs, DispState%LastOriginalPuff, DispState%FreeOriginalPuffStack, &
               DispState%iPuffs                                                                       &
             )
      End If

    End Do

    Call TimerOn(CPRSTimer)
    !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(OpenMPOpts%nParticleUpdateThreads)
    Call CalcParticleResultsSub(DispState,ReqInfoL,OldSyncTime,Mets,Flows,&
                                Coords,TargetTime,iCase,OutputOpts,&
                                Grids,Domains,Specieses,CloudGammaParamses,Sources,SizeDists,Reqs,&
                                Units,Results)
    !$OMP END PARALLEL
    Call TimerOff(CPRSTimer)
    Call TimerWriteLast(CPRSTimer)

    ! Convert puffs to particles.
    Do i = 1, DispState%LastPuff

      Puff     => DispState%Puffs(i)
      Particle => DispState%Puffs(i)%P
      Extra    => DispState%PuffExtras(i)

      If (.not. ParticleActive(Puff%P)) Cycle

      If (TravelTime(Particle) < DispOpts%sPuffTime) Cycle

      ! Turn off time spread of puff if this has not already happened by setting TInst - the output travel
      ! time beyond which no time spread is assumed. Note this will not take effect immediately.
      If (Extra%TInst == InfFutureShortTime(Interval = .true.)) Then
        Extra%TInst = Puff%P%T + Extra%TPlus
      End If

      If (ParticleMinTime(Particle, Extra) >= ParticleInstTime(Particle, Extra)) Then

        ! Calculate number of particles. !$$ need to think about fixed met. Also Part mass restrictions.
        Source => Sources%Sources(Particle%iSource)
        If (Source%ParticleRate) Then
          nParticles = ShortTime2RealTime(ParticleMaxTime(Particle, Extra) &
                     - ParticleMinTime(Particle, Extra)) &
                       * 0.5                &
                       * Source%nParticles  &
                       * 0.5 ** Puff%N
        Else
          nParticles = Source%nParticles * 0.5 ** Puff%N
        End If

        Call MetricCoeffs(Coords%HCoords(Particle%iHCoord), Particle%X(:), HMax, H(1), H(2))
        H(3) = 1.0 ! $$ OK if height above ground but note more sophisticated treatment in PuffConc.

        ! Create particles.
        Do j = 1, nParticles

          ! Replace by next block once limiting particle releases done smoothly as in routine Release.
          If (DispState%nParticles + 1 > DispState%MaxParticles) Then
            Call Message('WARNING: particles lost on puff-to-particle conversion - too many particles', 1)
          End If
          If (DispState%nParticleExtras + 1 > DispState%MaxFullParticles) Then
            Call Message('WARNING: particles lost on puff-to-particle conversion - ' // &
                         'too many "full" particles', 1)
          End If

    !      If (DispState%nParticles + 1 > DispState%MaxParticles) Then
    !        Call Message('UNEXPECTED FATAL ERROR in Case.Release', 4)
    !      End If
    !      If (DispState%nParticleExtras + 1 > DispState%MaxFullParticles) Then
    !        Call Message('FATAL ERROR: too many "full" particles', 3)
    !      End If

          DispState%nParticles = DispState%nParticles + 1
          iP = DispState%FreeParticleStack(DispState%nParticles)
          If (iP > DispState%LastParticle) DispState%LastParticle = iP
          DispState%nParticleExtras = DispState%nParticleExtras + 1
          iE = DispState%FreeParticleExtraStack(DispState%nParticleExtras)
          If (iE > DispState%LastParticleExtra) DispState%LastParticleExtra = iE
          DispState%iParticleExtras(iP) = iE
          DispState%iULastParticle = DispState%iULastParticle + 1
          iUP = DispState%iULastParticle

          DispState%Particles(iP)              = Particle
          DispState%Particles(iP)%iUP          = iUP
          DispState%Particles(iP)%NeedToSetVel = .true. ! $$ could use puff random velocity
          DispState%ParticleExtras(iP)         = Extra
          DispState%ParticleExtras(iP)%TInst   = InfPastShortTime(Interval = .true.) ! $$ Needed at present
                                                                                     ! because of nature
                                                                                     ! of test at
                                                                                     ! end of
                                                                                     ! oneparticletimestep
                                                                                     ! to transform
                                                                                     ! to a cheap particle.
          DispState%ParticleMasses(:, iP)      = DispState%PuffMasses(:, i) / nParticles

          ! Perturb position
          DispState%Particles(iP)%X(1) = DispState%Particles(iP)%X(1) + Sqrt(Puff%XXp(1)) * Gauss(iP) / H(1)
          DispState%Particles(iP)%X(2) = DispState%Particles(iP)%X(2) + Sqrt(Puff%XXp(2)) * Gauss(iP) / H(2)
          DispState%Particles(iP)%X(3) = DispState%Particles(iP)%X(3) + Sqrt(Puff%XXp(3)) * Gauss(iP) / H(3)
          If (DispState%Particles(iP)%X(3) < 0.0) DispState%Particles(iP)%X(3) = &
                                                - DispState%Particles(iP)%X(3)

        End Do

        ! Kill puff.
        Call MarkParticle(UseForOutput = .false., Particle = Particle)
        If (.not.ParticleReallyActive(Particle)) Then
          Call InactivatePuff(                                                                          &
                 i,                                                                                     &
                 DispState%Puffs,                                                                       &
                 DispState%nPuffs, DispState%LastPuff, DispState%FreePuffStack,                         &
                 DispState%nOriginalPuffs, DispState%LastOriginalPuff, DispState%FreeOriginalPuffStack, &
                 DispState%iPuffs                                                                       &
               )
        End If

      End If

    End Do






    ! Calculate type F (depending on flow information) fields.
    Call CalcFlowResults(ReqInfoF, Coords, Grids, Domains, Reqs, Units, Mets, Flows, Results)

    ! Calculate type O (other) fields.
    Call CalcOtherResults(ReqInfoO, RunToTime, Coords, Grids, Domains, Reqs, DispState, Results)

    Call TimerOn(OutputTimer)
    ! Output results.
    If (EndOfCase .or. DispState%PPsFinished) Then ! could attempt to avoid calling output if next output time
                                                   ! and end of case not reached
!      If (DispState%LastCase) EndOfCases = .true. $$ this causes metextraction runs to stop prematurely
      Call ProcessAndOutputResults(                                                                   &
             iCase, EndOfCase, EndOfCases .or. DispState%LastCase,                                    &
             Time2ShortTime(TargetTime), .true., .true.,                                              &
               ! $$ review
             ZeroShortTime(),                                                                         &
             DispOpts%RadioActiveDecay, DispOpts%Chemistry,                                           &
               ! Note need max run end time
             TMin(DispState%sFirstReqTime, DispState%sFirstReleaseTime),                              &
               ! over all cases for TimeForEField
             TMax(TMax(DispState%sLastDispTime, DispState%sLastReqTime), DispState%sLastReleaseTime), &
             OutputOpts, SizeDists, OpenMPOpts, Coords, Grids, Domains,                               &
             Specieses, CloudGammaParamses, Sources, Reqs, MaterialUnits,                             &
             Units, Mets, Flows, Results                                                              &
           )
    Else
      Call ProcessAndOutputResults(                                                                   &
             iCase, EndOfCase, EndOfCases,                                                            &
             Time2ShortTime(TargetTime), .true., .false.,                                             &
                ! $$ review
             DispState%MaxTSpread,                                                                    &
             DispOpts%RadioActiveDecay, DispOpts%Chemistry,                                           &
             TMin(DispState%sFirstReqTime, DispState%sFirstReleaseTime),                              &
             TMax(TMax(DispState%sLastDispTime, DispState%sLastReqTime), DispState%sLastReleaseTime), &
             OutputOpts, SizeDists, OpenMPOpts, Coords, Grids, Domains,                               &
             Specieses, CloudGammaParamses, Sources, Reqs, MaterialUnits,                             &
             Units, Mets, Flows, Results                                                              &
           )
    End If
    Call TimerOff(OutputTimer)
    Call TimerWriteLast(OutputTimer)

  End If

  If (DispState%Time == DispState%SyncTime) Then
    Call UpdateResults(                                          &
           Time2ShortTime(DispState%Time), DispState%MaxTSpread, &
           Grids, Reqs,                                          &
           Results                                               &
         )
  End If

  If (MainOpts%RunToFile /= ' ') Then
    Inquire(File = Trim(MainOpts%RunToFile), Exist = Exist)
    If (Exist) Then
      RunToUnit = OpenFile(MainOpts%RunToFile, Units, Status = 'Old', Action = 'Read') ! $$ May need
      ! error traps here, for editing file at same time. Or pause and try again.
      Read (RunToUnit, '(A)', IOStat = IOStat) CharRunToTimeFromFile
      Call CloseUnit(RunToUnit, Units)
      If (IOStat /= 0) Then
        Suspend = .true.
        ! Call Message $$
      Else If (CharRuntoTimeFromFile == ' ') Then
        Suspend = .true.
      Else
        RunToTimeFromFile = Char2Time(CharRunToTimeFromFile, ErrorCode = Char2TimeErrorCode)
        If (Char2TimeErrorCode > 0) Then
          ! Call Message $$
          Suspend = .true.
        Else If (DispState%Time >= RunToTimeFromFile) Then
          Suspend = .true.
        End If
      End If
    End If
  End If
  ! Give message to say run suspended by run-to file at ... $$
  ! Need to extend this to specify case to stop at (like command line runTo).
  ! Put this if block in a subroutine.

End Subroutine RunToSyncTimeOrMetFlowUpdateOrEndOfCase

!-------------------------------------------------------------------------------------------------------------

subroutine CalcParticleResultsSub(DispState,ReqInfoL,OldSyncTime,Mets,Flows,&
                                  Coords,TargetTime,iCase,&
                                  OutputOpts,Grids,Domains,&
                                  Specieses,CloudGammaParamses,Sources,SizeDists,Reqs,Units,&
                                  Results)

    implicit none

    Type(DispState_),      Intent(InOut), Target :: DispState
    Type(ReqInfo_),        Intent(In)            :: ReqInfoL(MaxReqInfoTimes)
    Type(ShortTime_),      Intent(In)            :: OldSyncTime
    Type(Mets_),           Intent(InOut)         :: Mets
    Type(Flows_),          Intent(InOut)         :: Flows
    Type(Coords_),         Intent(In)            :: Coords
    Integer,               Intent(In)            :: iCase
    Type(Time_),           Intent(In)            :: TargetTime
    Type(OutputOpts_),     Intent(In)            :: OutputOpts
    Type(Grids_),          Intent(In)            :: Grids
    Type(Domains_),        Intent(In)            :: Domains
    Type(Specieses_),      Intent(In)            :: Specieses
    Type(SizeDists_),      Intent(In)            :: SizeDists
    Type(CloudGammaParamses_), Intent(In)        :: CloudGammaParamses
    Type(Sources_),        Intent(In)            :: Sources
    Type(Reqs_),           Intent(In)            :: Reqs
    Type(Units_),          Intent(InOut)         :: Units
    Type(Results_),        Intent(InOut)         :: Results

!   locals
    Type(FlowMemory_)                :: FlowMemory
    Type(Particle_),         Pointer :: Particle
    Type(Extra_),            Pointer :: Extra
    Type(Position_)                  :: Position
    integer :: iReqInfo
    integer :: i,j

    logical :: Reprod=.false.

    !$OMP DO                                                    &
    !$OMP SCHEDULE(dynamic,16)
    Do i = 1, DispState%LastParticle

      Particle => DispState%Particles(i)
      Extra    => DispState%ParticleExtras(DispState%iParticleExtras(i))

      If (.not.ParticleActive(Particle)) Cycle

      If (UseParticleForOutput(Particle)) Then

        iReqInfo = 1
        Do
          If (ReqInfoL(iReqInfo)%NoReqs) Exit
          If (OldSyncTime < ReqInfoL(iReqInfo)%Time) Exit
          iReqInfo = iReqInfo + 1
        End Do

        Call ResetFlowMemory(Flows, FlowMemory)
        Position = X2Position(Coords, Particle%X, Particle%iHCoord, Particle%iZCoord)

        Do

          If (ReqInfoL(iReqInfo)%NoReqs) Exit

          If (Time2ShortTime(TargetTime) < ReqInfoL(iReqInfo)%Time) Exit

          Call CalcParticleResults(                                                  &
                 iCase,                                                              &
                 Particle, DispState%ParticleMasses(:, i), Extra,                    &
                 (/ (0.0, j = 1, MaxSpecieses) /),                                   & ! temp fix
                 (/ (0.0, j = 1, MaxSpecieses) /),                                   & ! temp fix
                 0.0, .true.,                                                        & ! temp fix
                 ReqInfoL(iReqInfo),                                                 &
                 OutputOpts,                                                         &
                 Coords, Grids, Domains, Mets, Flows,                                &
                 Specieses, CloudGammaParamses, Sources, SizeDists, Reqs,            &
                 Position, FlowMemory,                                               &
                 Units, Results, DispState                                           &
               )
          iReqInfo = iReqInfo + 1
        End Do

      End If

      If (.not.Reprod) Then
      ! Kill particle.
        If (.not.ParticleReallyActive(Particle)) Then
          !$OMP CRITICAL
          Call InactivateParticle(                                            &
               i,                                                           &
!!! AJM
               Coords,DispState%iHCoordLatLong,                             &
!!! AJM
               DispState%Particles,                                         &
               DispState%nParticles, DispState%FreeParticleStack,           &
               DispState%nParticleExtras, DispState%FreeParticleExtraStack, &
               DispState%iParticleExtras                                    &
             )
          !$OMP END CRITICAL
        End If
      End If

    End Do
    !$OMP END DO

  if (Reprod) THEN
    !$OMP SINGLE
    Do i = 1, DispState%LastParticle
      Particle => DispState%Particles(i)

      If (.not.ParticleActive(Particle)) Cycle

      ! Kill particle.
      If (.not.ParticleReallyActive(Particle)) Then
        Call InactivateParticle(                                            &
               i,                                                           &
!!! AJM
               Coords,DispState%iHCoordLatLong,                             &
!!! AJM
               DispState%Particles,                                         &
               DispState%nParticles, DispState%FreeParticleStack,           &
               DispState%nParticleExtras, DispState%FreeParticleExtraStack, &
               DispState%iParticleExtras                                    &
             )
      End If
    end do
    !$OMP END SINGLE
  end if

end subroutine CalcParticleResultsSub


!-------------------------------------------------------------------------------------------------------------

Subroutine LoopParticlesPuffsAndTimeSteps(                                       &
             iCase, Time, TargetTime,                                            &
             OutputOpts,                                                         &
             Coords, Grids, Domains, Mets, Flows,                                &
             Specieses, CloudGammaParamses, Sources, SizeDists, Reqs,            &
             DispOpts, Units, Results, DispState                                 &
           )
! Evolve particle/puffs, looping over particle/puffs and time steps.

  Implicit None
  ! Argument list:
  Integer,                   Intent(In)            :: iCase
  Type(ShortTime_),          Intent(In)            :: Time       !
  Type(ShortTime_),          Intent(In)            :: TargetTime !
  Type(OutputOpts_),         Intent(In)            :: OutputOpts
  Type(Coords_),             Intent(In)            :: Coords     !
  Type(Grids_),              Intent(In)            :: Grids      !
  Type(Domains_),            Intent(In)            :: Domains    !
  Type(Mets_),               Intent(InOut)         :: Mets       !
  Type(Flows_),              Intent(InOut)         :: Flows      !
  Type(Specieses_),          Intent(In)            :: Specieses  !
  Type(CloudGammaParamses_), Intent(In)            :: CloudGammaParamses
  Type(Sources_),            Intent(In)            :: Sources    !
  Type(SizeDists_),          Intent(In)            :: SizeDists
  Type(Reqs_),               Intent(In)            :: Reqs       !
  Type(DispOpts_),           Intent(In)            :: DispOpts   !
  Type(Units_),              Intent(InOut)         :: Units
  Type(Results_),            Intent(InOut)         :: Results    !
  Type(DispState_),          Intent(InOut), Target :: DispState  !
  ! iCase              :: Number of case.
  ! Time               :: Current time.
  ! TargetTime         :: Target time to which LoopOverSyncTimes should progress.
  ! OutputOpts         :: Output options.
  ! Coords             :: Collection of coord systems.
  ! Grids              :: Collection of grids.
  ! Domains            :: Collection of domains.
  ! Mets               :: Set of met module instance states.
  ! Flows              :: Set of flow module instance states.
  ! Specieses          :: Collection of species.
  ! CloudGammaParamses :: Collection of sets of cloud gamma parameters.
  ! Sources            :: Collection of sources.
  ! SizeDists          :: Collection of particle size distributions.
  ! Reqs               :: Collection of requirements.
  ! DispOpts           :: Set of dispersion model options.
  ! Units              :: Collection of information on input/output unit numbers.
  ! Results            :: Results.
  ! DispState          :: Dispersion model state.
  ! Locals:
  Type(Particle_),  Pointer :: Particle
  Type(Puff_),      Pointer :: Puff
  Type(Extra_),     Pointer :: Extra
  Integer                   :: iP
  Integer                   :: iOP
  Integer                   :: j
  Integer                   :: iReqInfoT
  Integer                   :: iReqInfoT1
  Integer                   :: iReqInfoS
  Type(ReqInfo_)            :: ReqInfoT(MaxReqInfoTimes)
  Type(ReqInfo_)            :: ReqInfoTPPInfos(MaxReqInfoTimes)
  Type(ReqInfo_)            :: ReqInfoS(MaxReqInfoTravelTimes)
  Type(ReqInfo_)            :: ReqInfoE
  Type(Flow_)               :: Flow
  Type(Cloud_)              :: Cloud
  Type(Rain_)               :: Rain
  Type(Surface_)            :: Surface
  Real(Std)                 :: DryDep(MaxSpecieses)
  Real(Std)                 :: WetDep(MaxSpecieses)
  Real(Std)                 :: RDt
  Integer                   :: iHCoord
  Integer                   :: iZCoord
  Real(Std)                 :: X(3)
  Type(Position_)           :: Position
  Type(FlowMemory_)         :: FlowMemory
  Logical                   :: PuffCentre
  Integer                   :: ErrorCode
  ! Particle        :: An abbreviation for the particle being considered or the
  !                    particle part of the puff being considered (in order to
  !                    enable a more unified treatment of particles and puffs).
  ! Puff            ::
  ! iP              ::
  ! j               ::
  ! iReqInfoT       ::
  ! iReqInfoT1      ::
  ! iReqInfoS       ::
  ! ReqInfoT        ::
  ! ReqInfoTPPInfos :: Info on requirements - sets of particle/puff information only.
  ! ReqInfoS        ::
  ! ReqInfoE        :: Every time output info
  ! Flow            ::
  ! Cloud           ::
  ! Rain            ::
  ! Surface         ::
  ! DryDep          :: Dry deposition losses (by species) over time-step RDt.
  ! WetDep          :: Wet deposition losses (by species) over time-step RDt.
  ! RDt             :: Time-step of particle/puff evolution.
  ! iHCoord         ::} Index of coord systems used when checking particles
  ! iZCoord         ::} are in the dispersion domain.
  ! X               :: Particle position in native coord systems of the domain.
  ! Position        :: Particle position.
  ! FlowMemory      ::
  ! ErrorCode       ::
  Logical          :: In
  Type(ShortTime_) :: ParticleTargetTime ! Target time for a particular particle/puff.
  Integer          :: i
  Logical          :: NoMorePuffs

  ! Calculate all t-based requirements over the current synchronisation step.
  ! Store these requirements in an array of ReqInfo_ variables (that is,
  ! the ReqInfo specifications at successive output times).
  Call CalcReqInfoT(                                 &
         .true., .true., 'L', .false.,               &
         .false., Time, TargetTime,                  & ! $$ use IncTime1 at start of run ??
         DispState%MaxTSpread, DispState%MaxTSpread, &
         Grids, Reqs,                                &
         ReqInfoT                                    &
       )

  ! Calculate PP Info t-based requirements
  Call CalcReqInfoT(                                 &
         .false., .true., 'L', .false.,              &
         .true., Time, TargetTime,                   &
         DispState%MaxTSpread, DispState%MaxTSpread, &
         Grids, Reqs,                                &
         ReqInfoTPPInfos                             &
       )

  ! Calculate all s-based requirements over the current synchronisation step.
  ! Store these requirements in an array of ReqInfo_ variables (that is,
  ! the ReqInfo specifications at successive output travel times).
  Call CalcReqInfoS(                                                                  &
         .false.,                                                                     &
         .false., Time, TargetTime,                                                   &
            ! $$ use IncTime1 at start of run ??
         .false., ZeroShortTime(), Domains%Domains(DispState%iDomain)%sMaxTravelTime, &
         DispState%MaxTSpread, DispState%MaxTSpread,                                  &
         Grids, Reqs,                                                                 &
         ReqInfoS                                                                     &
       )

  Call CalcReqInfoEveryT(                            &
         .false.,                                    &
         Time, TargetTime,                           &
         DispState%MaxTSpread, DispState%MaxTSpread, &
         Grids, Reqs,                                &
         ReqInfoE                                    &
       )

  ! Treat puffs when j = 1, particles when j = 2 (the following code
  ! is for the puffs.

    ! Do the aspects of the puff evolution which need to be carried out when the puffs
    ! are synchronised.
!$OMP SINGLE
      DispState%TLast(:)     = DispState%T(:)
      DispState%SigA2Last(:) = DispState%SigA2(:)

      ! Compute sigma_a.
      Call Step2(                                          &
             DispState%LastPuff, DispState%nOriginalPuffs, &
             DispState%Puffs, DispState%PuffMasses,        &
             DispState%SigA2, DispState%T                  &
           )

      ! Re-combine puffs.
      Call Step3(                                                                                   &
             DispOpts%PuffOpts,                                                                     &
             DispState%SigA2,                                                                       &
             DispState%Puffs, DispState%PuffExtras, DispState%PuffMasses,                           &
             DispState%nPuffs, DispState%LastPuff, DispState%FreePuffStack,                         &
             DispState%nOriginalPuffs, DispState%LastOriginalPuff, DispState%FreeOriginalPuffStack, &
             DispState%iPuffs                                                                       &
           )

      Do i = 1, DispState%nOriginalPuffs
        If (DispState%SigA2Defined(i) < 2) Then
          DispState%SigA2Defined(i) = DispState%SigA2Defined(i) + 1
        End If
      End Do

    ! Loop over puffs when j = 1 and over particles when j = 2.
    iOP = 0
    iP  = 0
    j=1 ! Puffs
    PuffLoop: Do

      ! Exit loop when all particle/puffs treated.
      Call NextPuff(DispState%LastOriginalPuff, DispState%Puffs, DispState%iPuffs, iP, iOP, NoMorePuffs)
      If (NoMorePuffs) Exit PuffLoop

      Call LoopParticlesPuffsAndTimeStepsSub(                                &
         j,iP,iOP,                                                           &
         iCase,                                                              &
         Time,                                                               &
         TargetTime,                                                         &
         OutputOpts,                                                         &
         Coords, Grids, Domains, Mets, Flows,                                &
         Specieses, CloudGammaParamses, Sources, SizeDists, Reqs,            &
         DispOpts, Units, Results, DispState,                                &
         ReqInfoT, ReqInfoTPPInfos, ReqInfoS, ReqInfoE                       &
       )

    End Do PuffLoop
!$OMP END SINGLE

! Process particles
    iOP = 0
    j=2 ! Particles
!$OMP DO SCHEDULE(static,16)
    ParticleLoop: Do iP = 1, DispState%LastParticle
      Call LoopParticlesPuffsAndTimeStepsSub(                                &
         j,iP,iOP,                                                           &
         iCase,                                                              &
         Time,                                                               &
         TargetTime,                                                         &
         OutputOpts,                                                         &
         Coords, Grids, Domains, Mets, Flows,                                &
         Specieses, CloudGammaParamses, Sources, SizeDists, Reqs,            &
         DispOpts, Units, Results, DispState,                                &
         ReqInfoT, ReqInfoTPPInfos, ReqInfoS, ReqInfoE                       &
       )

    End Do ParticleLoop
  ! $$ Kill particles with mass = 0?

End Subroutine LoopParticlesPuffsAndTimeSteps

!-------------------------------------------------------------------------------------------------------------

Subroutine LoopParticlesPuffsAndTimeStepsSub(                                    &
             j, iP, iOP,                                                         &
             iCase, Time, TargetTime,                                            &
             OutputOpts,                                                         &
             Coords, Grids, Domains, Mets, Flows,                                &
             Specieses, CloudGammaParamses, Sources, SizeDists, Reqs,            &
             DispOpts, Units, Results, DispState,                                &
             ReqInfoT, ReqInfoTPPInfos, ReqInfoS, ReqInfoE                       &
           )
! Evolve particle/puffs, looping over particle/puffs and time steps.

  Implicit None
  ! Argument list:
  Integer,                   Intent(In)            :: j
  Integer,                   Intent(In)            :: iP
  Integer,                   Intent(In)            :: iOP
  Integer,                   Intent(In)            :: iCase
  Type(ShortTime_),          Intent(In)            :: Time       !
  Type(ShortTime_),          Intent(In)            :: TargetTime !
  Type(OutputOpts_),         Intent(In)            :: OutputOpts
  Type(Coords_),             Intent(In)            :: Coords     !
  Type(Grids_),              Intent(In)            :: Grids      !
  Type(Domains_),            Intent(In)            :: Domains    !
  Type(Mets_),               Intent(InOut)         :: Mets       !
  Type(Flows_),              Intent(InOut)         :: Flows      !
  Type(Specieses_),          Intent(In)            :: Specieses  !
  Type(CloudGammaParamses_), Intent(In)            :: CloudGammaParamses
  Type(Sources_),            Intent(In)            :: Sources    !
  Type(SizeDists_),          Intent(In)            :: SizeDists
  Type(Reqs_),               Intent(In)            :: Reqs       !
  Type(DispOpts_),           Intent(In)            :: DispOpts   !
  Type(Units_),              Intent(InOut)         :: Units
  Type(Results_),            Intent(InOut)         :: Results    !
  Type(DispState_),          Intent(InOut), Target :: DispState  !
  Type(ReqInfo_),            Intent(In)      :: ReqInfoT(MaxReqInfoTimes)
  Type(ReqInfo_),            Intent(In)      :: ReqInfoTPPInfos(MaxReqInfoTimes)
  Type(ReqInfo_),            Intent(In)      :: ReqInfoS(MaxReqInfoTravelTimes)
  Type(ReqInfo_),            Intent(In)      :: ReqInfoE
  ! iCase              :: Number of case.
  ! Time               :: Current time.
  ! TargetTime         :: Target time to which LoopOverSyncTimes should progress.
  ! OutputOpts         :: Output options.
  ! Coords             :: Collection of coord systems.
  ! Grids              :: Collection of grids.
  ! Domains            :: Collection of domains.
  ! Mets               ::
  ! Flows              :: Set of flow module instance states.
  ! Specieses          :: Collection of species.
  ! CloudGammaParamses :: Collection of sets of cloud gamma parameters.
  ! Sources            :: Collection of sources.
  ! SizeDists          :: Collection of particle size distributions.
  ! Reqs               :: Collection of requirements.
  ! DispOpts           :: Set of dispersion model options.
  ! Units              :: Collection of information on input/output unit numbers.
  ! Results            :: Results.
  ! DispState          :: Dispersion model state.
  ! Locals:
  Type(Particle_),  Pointer :: Particle
  Type(Puff_),      Pointer :: Puff
  Type(Extra_),     Pointer :: Extra
  Integer                   :: iReqInfoT
  Integer                   :: iReqInfoT1
  Integer                   :: iReqInfoS
  Type(Flow_)               :: Flow
  Type(Cloud_)              :: Cloud
  Type(Rain_)               :: Rain
  Type(Surface_)            :: Surface
  Type(Soil_)               :: Soil
  Type(Plant_)              :: Plant
  Real(Std)                 :: DryDep(MaxSpecieses)
  Real(Std)                 :: WetDep(MaxSpecieses)
  Real(Std)                 :: RDt
  Integer                   :: iHCoord
  Integer                   :: iZCoord
  Real(Std)                 :: X(3)
  Type(Position_)           :: Position
  Type(FlowMemory_)         :: FlowMemory
  Logical                   :: PuffCentre
  Integer                   :: ErrorCode
  ! Particle        :: An abbreviation for the particle being considered or the
  !                    particle part of the puff being considered (in order to
  !                    enable a more unified treatment of particles and puffs).
  ! Puff            ::
  ! iP              ::
  ! j               ::
  ! iReqInfoT       ::
  ! iReqInfoT1      ::
  ! iReqInfoS       ::
  ! ReqInfoT        ::
  ! ReqInfoTPPInfos :: Info on requirements - sets of particle/puff information only.
  ! ReqInfoS        ::
  ! ReqInfoE        :: Every time output info
  ! Flow            ::
  ! Cloud           ::
  ! Rain            ::
  ! Surface         ::
  ! Soil            ::
  ! Plant           ::
  ! DryDep          :: Dry deposition losses (by species) over time-step RDt.
  ! WetDep          :: Wet deposition losses (by species) over time-step RDt.
  ! RDt             :: Time-step of particle/puff evolution.
  ! iHCoord         ::} Index of coord systems used when checking particles
  ! iZCoord         ::} are in the dispersion domain.
  ! X               :: Particle position in native coord systems of the domain.
  ! Position        :: Particle position.
  ! FlowMemory      ::
  ! ErrorCode       ::
  Logical          :: In
  Type(ShortTime_) :: ParticleTargetTime ! Target time for a particular particle/puff.
  Integer          :: i
  Logical          :: NoMorePuffs

      ! Set pointer to particle or to particle part of puff.
      If (j == 1) Then
        Particle => DispState%Puffs(iP)%P
        Puff     => DispState%Puffs(iP)
        Extra    => DispState%PuffExtras(iP)
      Else
        Particle => DispState%Particles(iP)
        Extra    => DispState%ParticleExtras(DispState%iParticleExtras(iP))
      End If

      ! If particle/puff not active, treat next particle/puff.
!GDR      If (.not.ParticleActive(Particle)) Cycle ParticlePuffLoop2
!GDR better before subroutine call in calling routine?!
      If (ParticleActive(Particle)) Then

      ! Output of sets of particle/puff information for particles/puffs that have just been released.
      If (TravelTime(Particle) == ZeroShortTime()) Then

        ! Reset index of ReqInfo arrays to beginning of each array.
        iReqInfoT = 1
        iReqInfoS = 1

        ! Calculate contribution to t-based output.

        ! Find an output time which the particle/puff can contribute to.
        Do
          If (ReqInfoTPPInfos(iReqInfoT)%NoReqs) Exit
          If (ParticleInstTime(Particle, Extra) <= ReqInfoTPPInfos(iReqInfoT)%Time) Then
            If (OldParticleTime(Particle) <= ReqInfoTPPInfos(iReqInfoT)%Time) Exit
              ! Note <= is different from similar block below
          Else
            If (OldParticleMinTime(Particle, Extra) <= ReqInfoTPPInfos(iReqInfoT)%Time) Exit
              ! Note <= is different from similar block below
          End If
          iReqInfoT = iReqInfoT + 1
        End Do

        iReqInfoT1 = iReqInfoT

        ! Calculate contributions for output times which the particle/puff can
        ! contribute to.
        Do
          If (ReqInfoTPPInfos(iReqInfoT1)%NoReqs) Exit
          If (ParticleInstTime(Particle, Extra) <= ReqInfoTPPInfos(iReqInfoT1)%Time) Then
            If (ParticleTime(Particle) < ReqInfoTPPInfos(iReqInfoT1)%Time) Exit
          Else
            If (ParticleMaxTime(Particle, Extra) < ReqInfoTPPInfos(iReqInfoT1)%Time) Exit
          End If
          
          Call ResetFlowMemory(Flows, FlowMemory)
          Position = X2Position(Coords, Particle%X, Particle%iHCoord, Particle%iZCoord)
   
          If (j == 1) Then
            PuffCentre =                                                          &
              OldParticleTime(Particle) <  ReqInfoTPPInfos(iReqInfoT1)%Time .and. &
              ParticleTime(Particle)    >= ReqInfoTPPInfos(iReqInfoT1)%Time
            Call CalcPuffResults(                                                                      &
                   iCase,                                                                              &
                   Puff, DispState%PuffMasses(:, iP), Extra, DryDep, WetDep, RDt,                      &
                   PuffCentre,                                                                         &
                   ReqInfoTPPInfos(iReqInfoT1),                                                        &
                   OutputOpts,                                                                         &
                   Coords, Grids, Domains, Mets, Flows, Specieses, Sources, SizeDists, Reqs, DispOpts, &
                   Position, FlowMemory,                                                               &
                   Units, Results, DispState                                                           &
                 )
          Else
            Call CalcParticleResults(                                                  &
                   iCase,                                                              &
                   Particle, DispState%ParticleMasses(:, iP), Extra, DryDep, WetDep, RDt, .true., &
                   ReqInfoTPPInfos(iReqInfoT1),                                        &
                   OutputOpts,                                                         &
                   Coords, Grids, Domains, Mets, Flows,                                &
                   Specieses, CloudGammaParamses, Sources, SizeDists, Reqs,            &
                   Position, FlowMemory,                                               &
                   Units, Results, DispState                                           &
                 )
          End If
          iReqInfoT1 = iReqInfoT1 + 1
        End Do
      End If

      ! Reset index of ReqInfo arrays to beginning of each array.
      iReqInfoT = 1
      iReqInfoS = 1

      ! Target time for this particle.
      ParticleTargetTime = TargetTime
      If (ParticleTargetTime > Particle%T0 + Domains%Domains(DispState%iDomain)%sMaxTravelTime) Then
        ParticleTargetTime = Particle%T0 + Domains%Domains(DispState%iDomain)%sMaxTravelTime
      End If
      If (                                                             &
        ParticleTargetTime >                                           &
        Particle%T0 + Sources%Sources(Particle%iSource)%sMaxTravelTime &
      ) Then
        ParticleTargetTime = Particle%T0 +                                    &
                             Sources%Sources(Particle%iSource)%sMaxTravelTime
      End If
      If (ParticleTargetTime > Domains%Domains(DispState%iDomain)%sEndTime) Then
        ParticleTargetTime = Domains%Domains(DispState%iDomain)%sEndTime
      End If

      ! $$ The above can result in ParticleTime > ParticleTargetTime (e.g. backwards run
      ! with domain end = infinity (incorrectly)). The following loop is not then
      ! executed with no opportunity to kill the particle. Should check end > start in
      ! setting up domain and either kill particles here or ensure situation cannot
      ! arise (e.g. test at release).

      ! Loop over time steps.
      TimeStepLoop: Do While (.not.(ParticleTime(Particle) >= ParticleTargetTime))

        ! Time-step particle/puff.
        If (j == 1) Then
          Call OnePuffTimeStep(          &
                 iP,                     &
                 ParticleTargetTime,     &
                 Coords, Grids, Domains, &
                 Units, Mets, Flows,     &
                 Specieses,              &
                 DispOpts,               &
                 Flow, Cloud, Rain,      &
                 Surface, Soil, Plant,   &
                 DryDep, WetDep, RDt,    &
                 DispState               &
               )
        Else
          Call OneParticleTimeStep(      &
                 iP,                     &
                 ParticleTargetTime,     &
                 Coords, Grids, Domains, &
                 Units, Mets, Flows,     &
                 Specieses,              &
                 DispOpts,               &
                 Flow, Cloud, Rain,      &
                 Surface, Soil, Plant,   &
                 DryDep, WetDep, RDt,    &
                 DispState               &
               )
        End If

        ! If particle/puff inactivated in OneParticleTimeStep/OnePuffTimeStep, exit
        ! TimeStepLoop. Note that, if a particle/puff is inactivated here, it
        ! doesn't have to be considered in calculating the (non-sync) output.
        If (.not.ParticleReallyActive(Particle)) Exit TimeStepLoop

        !---------------------------------------!
        ! Mark particle/puffs for inactivation. !
        !---------------------------------------!

        Call ResetFlowMemory(Flows, FlowMemory)

        Position = X2Position(Coords, Particle%X, Particle%iHCoord, Particle%iZCoord)

        ! Kill particles outside the computational domain in space.
        Call XInTheDomain(                                    &
               Coords, Grids, Domains, DispState%iDomain,     &
               Time, .false., TravelTime(Particle), Position, &
               Units, Mets, Flows,                            &
               FlowMemory,                                    &
               In, ErrorCode                                  &
             )
        If (.not.In) Then
          ! $$ If (ErrorCode == ..
          Call MarkParticle(UseForOutput = .false., Particle = Particle)
          Exit TimeStepLoop
        End If

        !------------------------------------!
        ! Calculate contributions to output. !
        !------------------------------------!

        ! Calculate contribution to every-time output.
        If (.not.ReqInfoE%NoReqs) Then
          If (j == 1) Then
            Call CalcPuffResults(                                                                      &
                   iCase,                                                                              &
                   Puff, DispState%PuffMasses(:, iP), Extra, DryDep, WetDep, RDt,                      &
                   .true.,                                                                             &
                   ReqInfoE,                                                                           &
                   OutputOpts,                                                                         &
                   Coords, Grids, Domains, Mets, Flows, Specieses, Sources, SizeDists, Reqs, DispOpts, &
                   Position, FlowMemory,                                                               &
                   Units, Results, DispState                                                           &
                 )
          Else
            Call CalcParticleResults(                                                  &
                   iCase,                                                              &
                   Particle, DispState%ParticleMasses(:, iP), Extra, DryDep, WetDep, RDt, .true., &
                   ReqInfoE,                                                           &
                   OutputOpts,                                                         &
                   Coords, Grids, Domains, Mets, Flows,                                &
                   Specieses, CloudGammaParamses, Sources, SizeDists, Reqs,            &
                   Position, FlowMemory,                                               &
                   Units, Results, DispState                                           &
                 )
          End If
        End If

        ! Calculate contribution to t-based output.

        ! Find an output time which the particle/puff can contribute to.
        Do
          If (ReqInfoT(iReqInfoT)%NoReqs) Exit
          If (ParticleInstTime(Particle, Extra) <= ReqInfoT(iReqInfoT)%Time) Then
            If (OldParticleTime(Particle) < ReqInfoT(iReqInfoT)%Time) Exit
          Else
            If (OldParticleMinTime(Particle, Extra) < ReqInfoT(iReqInfoT)%Time) Exit
          End If
          iReqInfoT = iReqInfoT + 1
        End Do

        iReqInfoT1 = iReqInfoT

        ! Calculate contributions for output times which the particle/puff can
        ! contribute to.
        Do
          If (ReqInfoT(iReqInfoT1)%NoReqs) Exit
          If (ParticleInstTime(Particle, Extra) <= ReqInfoT(iReqInfoT1)%Time) Then
            If (ParticleTime(Particle) < ReqInfoT(iReqInfoT1)%Time) Exit
          Else
            If (ParticleMaxTime(Particle, Extra) < ReqInfoT(iReqInfoT1)%Time) Exit
          End If
          If (j == 1) Then
            PuffCentre =                                                   &
              OldParticleTime(Particle) <  ReqInfoT(iReqInfoT1)%Time .and. &
              ParticleTime(Particle)    >= ReqInfoT(iReqInfoT1)%Time
            Call CalcPuffResults(                                                                      &
                   iCase,                                                                              &
                   Puff, DispState%PuffMasses(:, iP), Extra, DryDep, WetDep, RDt,                      &
                   PuffCentre,                                                                         &
                   ReqInfoT(iReqInfoT1),                                                               &
                   OutputOpts,                                                                         &
                   Coords, Grids, Domains, Mets, Flows, Specieses, Sources, SizeDists, Reqs, DispOpts, &
                   Position, FlowMemory,                                                               &
                   Units, Results, DispState                                                           &
                 )
          Else
           Call CalcParticleResults(                                                  &
                   iCase,                                                              &
                   Particle, DispState%ParticleMasses(:, iP), Extra, DryDep, WetDep, RDt, .true., &
                   ReqInfoT(iReqInfoT1),                                               &
                   OutputOpts,                                                         &
                   Coords, Grids, Domains, Mets, Flows,                                &
                   Specieses, CloudGammaParamses, Sources, SizeDists, Reqs,            &
                   Position, FlowMemory,                                               &
                   Units, Results, DispState                                           &
                 )
          End If
          iReqInfoT1 = iReqInfoT1 + 1
        End Do

        ! Calculate contribution to s-based output.

        ! Find an output time which the particle/puff can contribute to.
        Do
          If (ReqInfoS(iReqInfoS)%NoReqs) Exit
          If (OldTravelTime(Particle) < ReqInfoS(iReqInfoS)%Time) Exit
          iReqInfoS = iReqInfoS + 1
        End Do

        ! Calculate contributions for output times which the particle/puff can
        ! contribute to.
        Do
          If (ReqInfoS(iReqInfoS)%NoReqs) Exit
          If (TravelTime(Particle) < ReqInfoS(iReqInfoS)%Time) Exit
          If (j == 1) Then
            ! $$ Travel time output currently not supported for puffs.
          Else
            Call CalcParticleResults(                                                             &
                   iCase,                                                                         &
                   Particle, DispState%ParticleMasses(:, iP), Extra, DryDep, WetDep, RDt, .true., &
                   ReqInfoS(iReqInfoS),                                                           &
                   OutputOpts,                                                                    &
                   Coords, Grids, Domains, Mets, Flows,                                           &
                   Specieses, CloudGammaParamses, Sources, SizeDists, Reqs,                       &
                   Position, FlowMemory,                                                          &
                   Units, Results, DispState                                                      &
                 )
          End If
          iReqInfoS = iReqInfoS + 1
        End Do

        !---------------------------------------!
        ! Mark particle/puffs for inactivation. !
        !---------------------------------------!

        ! Kill particles/puffs which won't now contribute to the requirements.
        If (ParticleMinTime(Particle, Extra) >= DispState%sLastDispTime) Then
          Call MarkParticle(UseForOutput = .true., Particle = Particle)
          Exit TimeStepLoop
        End If

        ! Kill particle/puffs with travel time >= max travel time for the
        ! computational domain. ! $$ Use SInDomain function
        If (TravelTime(Particle) >= Domains%Domains(DispState%iDomain)%sMaxTravelTime) Then
          Call MarkParticle(UseForOutput = .true., Particle = Particle)
          Exit TimeStepLoop
        End If

        ! Kill particle/puffs with travel time >= max travel time for the source.
        If (                                                                       &
          TravelTime(Particle) >= Sources%Sources(Particle%iSource)%sMaxTravelTime &
        ) Then
          Call MarkParticle(UseForOutput = .true., Particle = Particle)
          Exit TimeStepLoop
        End If

        ! Kill particle/puffs outside the computational domain in time. $$ Not for fixed met.
        If (ParticleTime(Particle) >= Domains%Domains(DispState%iDomain)%sEndTime) Then
          Call MarkParticle(UseForOutput = .true., Particle = Particle)
          Exit TimeStepLoop
        End If

      End Do TimeStepLoop
      End if ! ParticleActive

End Subroutine LoopParticlesPuffsAndTimeStepsSub

!-------------------------------------------------------------------------------------------------------------

Subroutine Release(                                                                        &
             Time1, Time2,                                                                 &
             Coords, Grids, Domains, Mets, Flows, Specieses, Sources, SizeDists, TimeDeps, &
             DispOpts,                                                                     &
             Units,                                                                        &
             FixedMet,                                                                     &
             DispState,                                                                    &
             ReqInfoR,                                                                     &
             Reqs,                                                                         &
             Results                                                                       &
           )
! Releases particles and puffs.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In)            :: Time1
  Type(ShortTime_), Intent(In)            :: Time2
  Type(Coords_),    Intent(In)            :: Coords
  Type(Grids_),     Intent(In),    Target :: Grids
  Type(Domains_),   Intent(In)            :: Domains
  Type(Mets_),      Intent(InOut)         :: Mets
  Type(Flows_),     Intent(InOut)         :: Flows
  Type(Specieses_), Intent(In)            :: Specieses
  Type(Sources_),   Intent(In),    Target :: Sources
  Type(SizeDists_), Intent(In),    Target :: SizeDists
  Type(TimeDeps_),  Intent(In)            :: TimeDeps
  Type(DispOpts_),  Intent(In)            :: DispOpts
  Type(Units_),     Intent(InOut)         :: Units
  Logical,          Intent(In)            :: FixedMet
  Type(DispState_), Intent(InOut), Target :: DispState
  Type(ReqInfo_),   Intent(In)            :: ReqInfoR(:)
  Type(Reqs_),      Intent(In)            :: Reqs
  Type(Results_),   Intent(InOut)         :: Results
  ! Time1     :} Particles and puffs due for release
  ! Time2     :} in [Time1, Time2) are released.
  ! Coords    :: Collection of coord systems.
  ! Grids     :: Collection of grids.
  ! Domains   :: Collection of domains.
  ! Mets      :: Set of met module instance states.
  ! Flows     :: Set of flow module instance states.
  ! Specieses :: Collection of species.
  ! Sources   :: Collection of sources.
  ! SizeDists :: Collection of particle size distributions.
  ! TimeDeps  :: Collection of source time dependencies.
  ! DispOpts  :: Collection of sets of dispersion model options.
  ! Units     :: Collection of information on input/output unit numbers.
  ! FixedMet  :: Indicates a calculation with fixed met/flow.
  ! DispState :: Dispersion model state.
  ! ReqInfoR  :: Ouput requirements for source information

  ! Locals:
  Integer                       :: i                       ! Loop index.
  Integer                       :: iUP                     ! Unique index of particle/puff (i.e. an index unique to the
                                                           ! particle/puff which is not reused despite particle/puff
                                                           ! recycling). $$ Integer(8)
  Integer                       :: iP                      ! Index of particle/puff.
  Integer                       :: iUOP                    ! Unique index of original puff (i.e. an index unique to the
                                                           ! original puff which is not reused despite
                                                           ! puff recycling). $$ Integer(8)
  Integer                       :: iOP
  Integer                       :: iE                      ! Index of extra.
  Integer                       :: iS                      ! Index of source.
  Integer                       :: iX                      !
  Integer                       :: iY                      !
  Integer                       :: iSizeRange              !
  Integer                       :: nX                      !
  Integer                       :: nY                      !
  Integer                       :: nSizeRanges             !
  Type(ShortTime_)              :: PTimeM                  !} PTime is the particle/puff time, and PTimeP and PTimeM
  Type(ShortTime_)              :: PTime                   !} are the leading and trailing puff times.
  Type(ShortTime_)              :: PTimeP                  !}
  Integer                       :: nParticles              ! Number of particles to release.
  Real(Std)                     :: Area                    !
  Real(Std)                     :: Mass(MaxSpecieses)      ! Mass to release or mass release rate (for each species).
  Real(Std)                     :: Diameter                !
  Real(Std)                     :: ParticleShape           !
  Real(Std)                     :: Density                 !
  Logical                       :: Rate                    ! Indicates Mass is a release rate. Is true for
                                                           ! non-instantaneous sources with fixed met and false
                                                           ! otherwise.
  Type(HGrid_),    Pointer      :: HGrid                   !
  Type(SizeDist_), Pointer      :: SizeDist                !
  Type(Source_),   Pointer      :: Source                  !
  Type(Position_)               :: Position
  Type(Flow_)                   :: Flow
  Type(Cloud_)                  :: Cloud
  Type(Rain_)                   :: Rain
  Type(Surface_)                :: Surface
  Type(Soil_)                   :: Soil
  Type(Plant_)                  :: Plant
  Type(FlowMemory_)             :: FlowMemory
  Integer                       :: ParticleFactorType      !
  Integer                       :: ErrorCode               !
  Real(Std)                     :: HMax
  Real(Std)                     :: H1
  Real(Std)                     :: H2
  Real(Std)                     :: X(3)                    !
  Real(Std)                     :: dX(3)                   !
  Type(IterativePlumeRiseVars_) :: IterativePlumeRiseVars  !
  Real(Std)                     :: QmFinal                 !
  Real(Std)                     :: QmOriginal              !

  ! Loop over sources.
  SourcesLoop: Do iS = 1, Sources%nSources

    Source => Sources%Sources(iS)

    ! Test for source duration intersecting [Time1, Time2). Note source duration is
    ! [StartTime, StopTime] although, for particles with varying met, each particle
    ! is released at the start of the source period relevant to its mass, and so, for
    ! non-instantaneous sources, no particle can be released at StopTime.
    If (Source%sStopTime < Time1 .or. Source%sStartTime >= Time2) Cycle

    ! Varying  met.
    If (.not.FixedMet) Then

      ! Release particles.
      If (ZeroTime() >= DispOpts%PuffTime) Then

        If (Source%SourceType == S_Generic .Or. Source%SourceType == S_IterativePlumeModel) Then
          nX          = 1
          nY          = 1
          nSizeRanges = 1
        Else If (Source%SourceType == S_Dust .or. Source%SourceType == S_SeaSalt) Then
          HGrid       => Grids%HGrids(Source%iHGrid)
          SizeDist    => SizeDists%SizeDists(Source%iSizeDist)
          nX          =  HGrid%nX
          nY          =  HGrid%nY
          nSizeRanges =  SizeDist%nSizeRanges
        End If

        YLoop: Do iY = 1, nY
        XLoop: Do iX = 1, nX

          ! Get flow and surface information for dust & sea-salt sources.
          If (Source%SourceType == S_Dust .or. Source%SourceType == S_SeaSalt) Then

            ! Set X, dX.
            X(1)  = HGrid%X(iX)
            X(2)  = HGrid%Y(iY)
            X(3)  = Source%X(3)
            dX(1) = HGrid%dX
            dX(2) = HGrid%dY
            dX(3) = Source%dX(3)

            ! Calculate area.
            Call MetricCoeffs(Coords%HCoords(HGrid%iHCoord), (/ HGrid%X(iX), HGrid%Y(iY) /), HMax, H1, H2)
            Area = HGrid%dX * HGrid%dY * H1 * H2 ! $$ Assumes regular grid.

            Call ResetFlowMemory(Flows, FlowMemory)
! $$ Note height hard wired to 10m to ensure correct wind speed for sea salt releases
            Position = X2Position(                                &
                         Coords,                                  &
                         (/ HGrid%X(iX), HGrid%Y(iY), 10.0_Std /), &
                         HGrid%iHCoord,                           &
                         DispState%iZCoordMagl                    &
                       )
            Call GetAttrib(                                                       &
                   iAttribParam  = A_Flow,                                        &
                   Coords        = Coords,                                        &
                   Grids         = Grids,                                         &
                   Domains       = Domains,                                       &
                   Moisture      = .true.,                                        & ! $$
                   Inhomog       = DispOpts%sInhomogTime > ZeroShortTime(),       &
                   Homog         = .not. DispOpts%sInhomogTime > ZeroShortTime(), &
                   Time          = Time1,                                         &
       !                $$  DispState%NextReleaseTime(iS)%T(iX, iY, iSizeRange) would be better, but
       !                    would have to be inside size range loop, and
       !                    inside particle loop for best accuracy.
                   AnyTravelTime = .false.,                                       &
                   TravelTime    = ZeroShortTime(),                               &
                   Position      = Position,                                      &
                   Units         = Units,                                         &
                   Mets          = Mets,                                          &
                   Flows         = Flows,                                         &
                   FlowMemory    = FlowMemory,                                    &
                   Flow          = Flow,                                          &
                   Cloud         = Cloud,                                         &
                   Rain          = Rain,                                          &
                   Surface       = Surface,                                       &
                   Soil          = Soil,                                          &
                   Plant         = Plant,                                         &
                   ErrorCode     = ErrorCode                                      &
                 )
            If (ErrorCode > 0) Then
              If (ErrorCode == 3) Then
                Call Message(                                                         &
                       'FATAL ERROR: no flow attribute defined in the collection ' // &
                       'of flow module instance attributes',                          &
                       3                                                              &
                     )
              Else If (.not. DispState%PPsLostDueToFlow) Then
                DispState%PPsLostDueToFlow = .true.
                Call ControlledMessage(                                                             &
                       'Particles and/or puffs are being lost due to lack of a valid flow module.', &
                       MessageControls    = GlobalMessageControls,                                  &
                       MessageControlName = 'No flow for particle/puff',                            &
                       ErrorCode          = 2                                                       &
                     )
              End If
      !       ! If (ErrorCode == 1) Then
      !       !   Call Message('Error: no flow module instance suitable for ' // &
      !       !                'vertical coordinate conversion', 2)
      !       ! Else
      !       !   Call Message('Error: no flow module instance suitable for ' // &
      !       !                'supplying the flow information', 2)
      !       ! End If
              Cycle XLoop
            End If

            Call GetAttrib(                                                       &
                   iAttribParam  = A_Surface,                                     &
                   Coords        = Coords,                                        &
                   Grids         = Grids,                                         &
                   Domains       = Domains,                                       &
                   Moisture      = .true.,                                        & ! $$
                   Inhomog       = DispOpts%sInhomogTime > ZeroShortTime(),       &
                   Homog         = .not. DispOpts%sInhomogTime > ZeroShortTime(), &
                   Time          = Time1,                                         & ! $$ see comment
                                                                                    ! above on previous call
                   AnyTravelTime = .false.,                                       &
                   TravelTime    = ZeroShortTime(),                               &
                   Position      = Position,                                      &
                   Units         = Units,                                         &
                   Mets          = Mets,                                          &
                   Flows         = Flows,                                         &
                   FlowMemory    = FlowMemory,                                    &
                   Flow          = Flow,                                          &
                   Cloud         = Cloud,                                         &
                   Rain          = Rain,                                          &
                   Surface       = Surface,                                       &
                   Soil          = Soil,                                          &
                   Plant         = Plant,                                         &
                   ErrorCode     = ErrorCode                                      &
                 )
            If (ErrorCode > 0) Then
              If (ErrorCode == 3) Then
                Call Message(                                                            &
                       'FATAL ERROR: no surface attribute defined in the collection ' // &
                       'of flow module instance attributes',                             &
                       3                                                                 &
                     )
              Else If (.not. DispState%PPsLostDueToFlow) Then
                DispState%PPsLostDueToFlow = .true. ! $$ this message could be more helpful and
                                                    ! mention "surface"
                Call ControlledMessage(                                                             &
                       'Particles and/or puffs are being lost due to lack of a valid flow module.', &
                       MessageControls    = GlobalMessageControls,                                  &
                       MessageControlName = 'No flow for particle/puff',                            &
                       ErrorCode          = 2                                                       &
                     )
              End If
      !       ! If (ErrorCode == 1) Then
      !       !   Call Message('Error: no flow module instance suitable for ' // &
      !       !                'vertical coordinate conversion', 2)
      !       ! Else
      !       !   Call Message('Error: no flow module instance suitable for ' // &
      !       !                'supplying the flow information', 2)
      !       ! End If
              Cycle XLoop
            End If


            Call GetAttrib(                                                       &
                   iAttribParam  = A_Soil,                                        &
                   Coords        = Coords,                                        &
                   Grids         = Grids,                                         &
                   Domains       = Domains,                                       &
                   Moisture      = .true.,                                        & ! $$
                   Inhomog       = DispOpts%sInhomogTime > ZeroShortTime(),       &
                   Homog         = .not. DispOpts%sInhomogTime > ZeroShortTime(), &
                   Time          = Time1,                                         & ! $$ see comment
                                                                                    ! above on previous call
                   AnyTravelTime = .false.,                                       &
                   TravelTime    = ZeroShortTime(),                               &
                   Position      = Position,                                      &
                   Units         = Units,                                         &
                   Mets          = Mets,                                          &
                   Flows         = Flows,                                         &
                   FlowMemory    = FlowMemory,                                    &
                   Flow          = Flow,                                          &
                   Cloud         = Cloud,                                         &
                   Rain          = Rain,                                          &
                   Surface       = Surface,                                       &
                   Soil          = Soil,                                          &
                   Plant         = Plant,                                         &
                   ErrorCode     = ErrorCode                                      &
                 )
            If (ErrorCode > 0) Then
              If (ErrorCode == 3) Then
                Call Message(                                                            &
                       'FATAL ERROR: no soil attribute defined in the collection ' // &
                       'of flow module instance attributes',                             &
                       3                                                                 &
                     )
              Else If (.not. DispState%PPsLostDueToFlow) Then
                DispState%PPsLostDueToFlow = .true. ! $$ this message could be more helpful and mention "soil"
                Call ControlledMessage(                                                             &
                       'Particles and/or puffs are being lost due to lack of a valid flow module.', &
                       MessageControls    = GlobalMessageControls,                                  &
                       MessageControlName = 'No flow for particle/puff',                            &
                       ErrorCode          = 2                                                       &
                     )
              End If
              Cycle XLoop
            End If

            If (IsBackwards()) Call ReverseFlow(Flow)

          ! Get flow information for met dependent chemistry sources.
          Else If (Source%MetDependent == 1) Then

            ! Set X, dX.
            X(:)  = Source%X(:)
            dX(:) = Source%dX(:)
            !force metdependent sources to take surface temperature, rather
            ! than T at particle release ht
            X(3) = 0.0

            Call ResetFlowMemory(Flows, FlowMemory)
            Position = X2Position(Coords, X, Source%iHCoord, Source%iZcoord)
            Call GetAttrib(                                                       &
                   iAttribParam  = A_Flow,                                        &
                   Coords        = Coords,                                        &
                   Grids         = Grids,                                         &
                   Domains       = Domains,                                       &
                   Moisture      = .true.,                                        & ! $$
                   Inhomog       = DispOpts%sInhomogTime > ZeroShortTime(),       &
                   Homog         = .not. DispOpts%sInhomogTime > ZeroShortTime(), &
                   Time          = Time1,                                         &
       !                $$  DispState%NextReleaseTime(iS)%T(iX, iY, iSizeRange) would be better, but
       !                    would have to be inside size range loop, and inside
       !                    particle loop for best accuracy.
                   AnyTravelTime = .false.,                                       &
                   TravelTime    = ZeroShortTime(),                               &
                   Position      = Position,                                      &
                   Units         = Units,                                         &
                   Mets          = Mets,                                          &
                   Flows         = Flows,                                         &
                   FlowMemory    = FlowMemory,                                    &
                   Flow          = Flow,                                          &
                   Cloud         = Cloud,                                         &
                   Rain          = Rain,                                          &
                   Surface       = Surface,                                       &
                   Soil          = Soil,                                          &
                   Plant         = Plant,                                         &
                   ErrorCode     = ErrorCode                                      &
                 )
            If (ErrorCode > 0) Then
              If (ErrorCode == 3) Then
                Call Message(                                                         &
                       'FATAL ERROR: no flow attribute defined in the collection ' // &
                       'of flow module instance attributes',                          &
                       3                                                              &
                     )
              Else If (.not. DispState%PPsLostDueToFlow) Then
                DispState%PPsLostDueToFlow = .true.
                Call ControlledMessage(                                                             &
                       'Particles and/or puffs are being lost due to lack of a valid flow module.', &
                       MessageControls    = GlobalMessageControls,                                  &
                       MessageControlName = 'No flow for particle/puff',                            &
                       ErrorCode          = 2                                                       &
                     )
              End If
      !       ! If (ErrorCode == 1) Then
      !       !   Call Message('Error: no flow module instance suitable for ' // &
      !       !                'vertical coordinate conversion', 2)
      !       ! Else
      !       !   Call Message('Error: no flow module instance suitable for ' // &
      !       !                'supplying the flow information', 2)
      !       ! End If
              Cycle XLoop
            End If
            !reset X, dX (from forcing metdependent sources to take surface T above)
            X(:)  = Source%X(:)
            dX(:) = Source%dX(:)
            Call ResetFlowMemory(Flows, FlowMemory)
            Position = X2Position(Coords, X, Source%iHCoord, Source%iZcoord)

            If (IsBackwards()) Call ReverseFlow(Flow)

            Call GetAttrib(                                                       &
                   iAttribParam  = A_Rain,                                        &
                   Coords        = Coords,                                        &
                   Grids         = Grids,                                         &
                   Domains       = Domains,                                       &
                   Moisture      = .false.,                                       & ! $$
                   Inhomog       = .false.,                                       &
                   Homog         = .false.,                                       &
                   Time          = Time1,                                         &
       !                $$  DispState%NextReleaseTime(iS)%T(iX, iY, iSizeRange) would be better, but would
       !                    have to be inside size range loop, and inside particle loop for best accuracy.
                   AnyTravelTime = .false.,                                       &
                   TravelTime    = ZeroShortTime(),                               &
                   Position      = Position,                                      &
                   Units         = Units,                                         &
                   Mets          = Mets,                                          &
                   Flows         = Flows,                                         &
                   FlowMemory    = FlowMemory,                                    &
                   Flow          = Flow,                                          &
                   Cloud         = Cloud,                                         &
                   Rain          = Rain,                                          &
                   Surface       = Surface,                                       &
                   Soil          = Soil,                                          &
                   Plant         = Plant,                                         &
                   ErrorCode     = ErrorCode                                      &
                 )
            If (ErrorCode > 0) Then
              If (ErrorCode == 3) Then
                Call Message(                                                         &
                       'FATAL ERROR: no rain attribute defined in the collection ' // &
                       'of flow module instance attributes',                          &
                       3                                                              &
                     )
              Else If (.not. DispState%PPsLostDueToFlow) Then
                DispState%PPsLostDueToFlow = .true.
                Call ControlledMessage(                                                             &
                       'Particles and/or puffs are being lost due to lack of a valid flow module.', &
                       MessageControls    = GlobalMessageControls,                                  &
                       MessageControlName = 'No flow for particle/puff',                            &
                       ErrorCode          = 2                                                       &
                     ) ! $$ Should use 'No rain for particle/puff' control
              End If
      !       ! If (ErrorCode == 1) Then
      !       !   Call Message('Error: no flow module instance suitable for ' // &
      !       !                'vertical coordinate conversion', 2)
      !       ! Else
      !       !   Call Message('Error: no flow module instance suitable for ' // &
      !       !                'supplying the rain information', 2)
      !       ! End If
              Cycle XLoop
            End If

          Else If (Source%MetDependent == 2) Then                               ! Iterative plume model 
 
            Call FlowProfilesAtSource(Time1, Source, Coords, Grids, Domains, &
                                      Mets, Flows, FlowMemory, Units,        &
                                      IterativePlumeRiseVars) 
            Call IteratePlumeRiseModel(Source, IterativePlumeRiseVars,       &
                                      QmOriginal, QmFinal) 

            ! Convert Qm from kg/s to g/s (species material unit - assumed to be 'g')    
            If ( Trim(AdjustL(Specieses%Specieses(1)%MaterialUnitName)) /= 'g' ) Then
              Call Message('FATAL ERROR: material unit of species must be in g when using iterative plume rise model', 3)
            End If
            QmOriginal = QmOriginal * 1000.0
            QmFinal    = QmFinal * 1000.0

            ! Pass mass fluxes to output routine
            Call CalcSourceResults(ReqInfoR, Reqs, Results, (/QmOriginal, QmFinal/))

            ! For volcanic ash, scale QmFinal by the input distal fine ash fraction
            If ( Source%SpeciesNames(1) .CIEq. 'VOLCANIC_ASH' ) QmFinal = QmFinal * Source%DistalFineAshFraction

            X(:)  = Source%X(:)
            dX(:) = Source%dX(:)

          Else

            ! Set X, dX.
            X(:)  = Source%X(:)
            dX(:) = Source%dX(:)

          End If

          SizeRangeLoop: Do iSizeRange = 1, nSizeRanges

            ParticleLoop: Do

              PTime = DispState%NextReleaseTime(iS)%T(iX, iY, iSizeRange)

              ! Exit if next release not due in this call to this routine.
              If (PTime >= Time2) Exit

              Call ParticleReleaseInfo(                                                                     &
                     PTime, Source,                                                                         &
                     DispOpts%ParticleCeiling, DispOpts%ParticleFactor, DispOpts%MaxParticles,              &
                     DispState%nParticles,                                                                  &
                     Specieses, SizeDists, TimeDeps,                                                        &
                     Flow, Rain, Surface, Soil, Area, iSizeRange,                                           &
                     Mass, Diameter, Density, ParticleShape,                                                &
                     Rate, nParticles,                                                                      &
                     DispState%NextReleaseTime(iS)%T(iX, iY, iSizeRange),                                   &
                     ParticleFactorType, QmFinal                                                            &
                   )

              ! Warning messages for reductions in number of particles released.
              If (                                                       &
                DispState%ParticleFactorType /= ParticleFactorType .and. &
                .not.(nParticles == 0 .and. ParticleFactorType == 0)     &
              ) Then
                If (ParticleFactorType == 0) Then
                  Call Message(                                                                             &
                         'Particle releases returned to normal (this will apply to particle splitting too)' &
                       )
                Else If (ParticleFactorType == 1 .and. DispState%ParticleFactorType == 0) Then
                  Call Message(                                                             &
                         'WARNING: Particle releases reduced as nearing particle limit ' // &
                         '(this will apply to particle splitting too)',                     &
                         1                                                                  &
                       )
                Else If (ParticleFactorType == 1 .and. DispState%ParticleFactorType == 2) Then
                  Call Message(                                                                            &
                         'Particle releases restarted but at a reduced rate as nearing particle limit ' // &
                         '(this will apply to particle splitting too)'                                     &
                       )
                Else If (ParticleFactorType == 2) Then
                  Call Message(                                                             &
                         'WARNING: Particle releases stopped as particle limit reached ' // &
                         '(this will apply to particle splitting too)',                     &
                         1                                                                  &
                       )
                End If
                DispState%ParticleFactorType = ParticleFactorType
              End If

              Do i = 1, nParticles

                If (DispState%nParticles + 1 > DispState%MaxParticles) Then
                  Call Message('UNEXPECTED FATAL ERROR in Case.Release', 4)
                End If
                If (DispState%nParticleExtras + 1 > DispState%MaxFullParticles) Then
                  Call Message('FATAL ERROR: too many "full" particles', 3)
                End If

                DispState%nParticles = DispState%nParticles + 1
                iP = DispState%FreeParticleStack(DispState%nParticles)
                If (iP > DispState%LastParticle) DispState%LastParticle = iP
                DispState%nParticleExtras = DispState%nParticleExtras + 1
                iE = DispState%FreeParticleExtraStack(DispState%nParticleExtras)
                If (iE > DispState%LastParticleExtra) DispState%LastParticleExtra = iE
                DispState%iParticleExtras(iP) = iE
                DispState%iULastParticle = DispState%iULastParticle + 1
                iUP = DispState%iULastParticle

                Call ReleaseParticle(                                  &
                       PTime, PTime, PTime,                            &
                       DispOpts, Units,                                &
                       Coords, Grids, Domains, Mets, Flows,            & 
                       Specieses, iS, Source, SizeDists,               &
                       Mass, Diameter, Density, ParticleShape,         &
                       X, dX,                                          &
                       iUP, iP, iE, DispState                          &
                     )

              End Do

              ! Exit if last release has been made for this source (subsequent entries to
              ! this routine won't make further releases because NextReleaseTime has been
              ! set to infinity by ParticleReleaseInfo).
              If (DispState%NextReleaseTime(iS)%T(iX, iY, iSizeRange) >= Source%sStopTime) Exit

            End Do ParticleLoop

          End Do SizeRangeLoop

        End Do XLoop
        End Do YLoop

      ! Release puffs.
      Else

        Do

          PTimeM = DispState%NextReleaseTimeM(iS)
          PTime  = DispState%NextReleaseTime (iS)%T(1, 1, 1)
          PTimeP = TMin(PTime + Time2ShortTime(DispOpts%PuffInterval), Source%sStopTime)

          ! Exit if next release not due in this call to this routine.
          If (PTime >= Time2) Exit

          Call PuffReleaseInfo(                 &
                 PTimeM, PTime, PTimeP, Source, &
                 Specieses, TimeDeps,           &
                 Mass, Rate                     &
               )
          DispState%NextReleaseTime(iS)%T(1, 1, 1) = PTimeP
          DispState%NextReleaseTimeM(iS)           = PTime

          DispState%MaxTSpread = TMax(DispState%MaxTSpread, PTimeP - PTime)

          Call GetNewOriginalPuffIndex(                                                                 &
                 iUOP, iUP, iOP, iP,                                                                    &
                 DispState%nOriginalPuffs, DispState%LastOriginalPuff, DispState%FreeOriginalPuffStack, &
                 DispState%iULastOriginalPuff,                                                          &
                 DispState%nPuffs, DispState%LastPuff, DispState%FreePuffStack, DispState%iULastPuff,   &
                 DispState%iPuffs                                                                       &
               )

          Call ReleasePuff(                                      &
                 PTimeM, PTime, PTimeP,                          &
                 DispOpts, Units,                                &
                 Coords, Grids, Domains, Mets, Flows, Specieses, &
                 iS, Source, SizeDists, Mass,                    &
                 iUP, iP, iP, iUOP, iOP, DispState               &
               )

          ! Exit if last release has been made for this source (subsequent entries to
          ! this routine won't make further releases).
          If (PTime == PTimeP) Exit

        End Do

      End If

    ! Fixed met.
    Else

      ! Release particles.
      If (ZeroTime() >= DispOpts%PuffTime) Then

        PTimeM = Source%sStartTime
        PTime  = DispState%NextReleaseTime(iS)%T(1, 1, 1)
        PTimeP = Source%sStopTime

        ! Cycle if release not due in this call to this routine.
        If (PTime >= Time2) Cycle

        Call ParticleReleaseInfo(                                                          &
               PTime, Source,                                                              &
               DispOpts%ParticleCeiling, DispOpts%ParticleFactor, DispOpts%MaxParticles,   &
               DispState%nParticles,                                                       &
               Specieses, SizeDists, TimeDeps,                                             &
               Flow, Rain, Surface, Soil, Area, iSizeRange,                                &
               Mass, Diameter, Density, ParticleShape,                                     &
               Rate, nParticles,                                                           &
               DispState%NextReleaseTime(iS)%T(1, 1, 1),                                   &
               ParticleFactorType, QmFinal                                                 &
             )

        ! Warning messages for reductions in number of particles released.
        If (                                                       &
          DispState%ParticleFactorType /= ParticleFactorType .and. &
          .not.(nParticles == 0 .and. ParticleFactorType == 0)     &
        ) Then
          If (ParticleFactorType == 0) Then
            Call Message(                                                                             &
                   'Particle releases returned to normal (this will apply to particle splitting too)' &
                 )
          Else If (ParticleFactorType == 1 .and. DispState%ParticleFactorType == 0) Then
            Call Message(                                                             &
                   'WARNING: Particle releases reduced as nearing particle limit ' // &
                   '(this will apply to particle splitting too)',                     &
                   1                                                                  &
                 )
          Else If (ParticleFactorType == 1 .and. DispState%ParticleFactorType == 2) Then
            Call Message(                                                                            &
                   'Particle releases restarted but at a reduced rate as nearing particle limit ' // &
                   '(this will apply to particle splitting too)'                                     &
                 )
          Else If (ParticleFactorType == 2) Then
            Call Message(                                                             &
                   'WARNING: Particle releases stopped as particle limit reached ' // &
                   '(this will apply to particle splitting too)',                     &
                   1                                                                  &
                 )
          End If
          DispState%ParticleFactorType = ParticleFactorType
        End If

        Do i = 1, nParticles

          If (DispState%nParticles + 1 > DispState%MaxParticles) Then
            Call Message('UNEXPECTED FATAL ERROR in Case.Release', 4)
          End If
          If (DispState%nParticleExtras + 1 > DispState%MaxFullParticles) Then
            Call Message('FATAL ERROR: too many "full" particles', 3)
          End If

          DispState%nParticles = DispState%nParticles + 1
          iP = DispState%FreeParticleStack(DispState%nParticles)
          If (iP > DispState%LastParticle) DispState%LastParticle = iP
          DispState%nParticleExtras = DispState%nParticleExtras + 1
          iE = DispState%FreeParticleExtraStack(DispState%nParticleExtras)
          If (iE > DispState%LastParticleExtra) DispState%LastParticleExtra = iE
          DispState%iParticleExtras(iP) = iE
          DispState%iULastParticle = DispState%iULastParticle + 1
          iUP = DispState%iULastParticle

          Call ReleaseParticle(                                              &
                 PTimeM, PTime, PTimeP,                                      & ! $$ Triangle v square t-shape
                 DispOpts, Units,                                            &
                 Coords, Grids, Domains, Mets, Flows,                        &
                 Specieses, iS, Source, SizeDists,                           &
                 Mass, Diameter, Density, ParticleShape,                     &
                 Source%X, Source%dX,                                        &
                 iUP, iP, iE, DispState                                      &
               )

        End Do

      ! Release puffs.
      Else

        PTimeM = Source%sStartTime
        PTime  = DispState%NextReleaseTime(iS)%T(1, 1, 1)
        PTimeP = Source%sStopTime

        ! Cycle if release not due in this call to this routine.
        If (PTime >= Time2) Cycle

        Call PuffReleaseInfo(                 &
               PTimeM, PTime, PTimeP, Source, &
               Specieses, TimeDeps,           &
               Mass, Rate                     &
             )
        DispState%NextReleaseTime(iS)%T(1, 1, 1) = InfFutureShortTime()

        Call GetNewOriginalPuffIndex(                                                                 &
               iUOP, iUP, iOP, iP,                                                                    &
               DispState%nOriginalPuffs, DispState%LastOriginalPuff, DispState%FreeOriginalPuffStack, &
               DispState%iULastOriginalPuff,                                                          &
               DispState%nPuffs, DispState%LastPuff, DispState%FreePuffStack, DispState%iULastPuff,   &
               DispState%iPuffs                                                                       &
             )

        Call ReleasePuff(                                      &
               PTimeM, PTime, PTimeP,                          & ! $$ Triangle v square t-shape
               DispOpts, Units,                                &
               Coords, Grids, Domains, Mets, Flows, Specieses, &
               iS, Source, SizeDists, Mass,                    &
               iUP, iP, iP, iUOP, iOP, DispState               &
             )

      End If

    End If

  End Do SourcesLoop

End Subroutine Release

!-------------------------------------------------------------------------------------------------------------

Subroutine ReleaseParticle(                                  &
             TimeM, Time, TimeP,                             &
             DispOpts,                                       &
             Units,                                          &
             Coords, Grids, Domains, Mets, Flows, Specieses, &
             iS, Source, SizeDists,                          &
             Mass, Diameter, Density, ParticleShape,         &
             X, dX,                                          &
             iUP, iP, iE, DispState                          &
           )
! Releases a particle.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In)            :: TimeM
  Type(ShortTime_), Intent(In)            :: Time
  Type(ShortTime_), Intent(In)            :: TimeP
  Type(DispOpts_),  Intent(In)            :: DispOpts
  Type(Units_),     Intent(InOut)         :: Units
  Type(Coords_),    Intent(In)            :: Coords
  Type(Grids_),     Intent(In)            :: Grids
  Type(Domains_),   Intent(In)            :: Domains
  Type(Mets_),      Intent(InOut)         :: Mets
  Type(Flows_),     Intent(InOut)         :: Flows
  Type(Specieses_), Intent(In)            :: Specieses
  Integer,          Intent(In)            :: iS
  Type(Source_),    Intent(In)            :: Source
  Type(SizeDists_), Intent(In)            :: SizeDists
  Real(Std),        Intent(In)            :: Mass(MaxSpecieses)
  Real(Std),        Intent(In)            :: Diameter
  Real(Std),        Intent(In)            :: Density
  Real(Std),        Intent(In)            :: ParticleShape
  Real(Std),        Intent(In)            :: X(3)
  Real(Std),        Intent(In)            :: dX(3)
  Integer,          Intent(In)            :: iUP
  Integer,          Intent(In)            :: iP
  Integer,          Intent(In)            :: iE
  Type(DispState_), Intent(InOut), Target :: DispState
  ! Time      ::
  ! DispOpts  :: Collection of sets of dispersion model options.
  ! Units     :: Collection of information on input/output unit numbers.
  ! Coords    :: Collection of coord systems.
  ! Grids     :: Collection of grids.
  ! Domains   :: Collection of domains.
  ! Mets      :: Set of met module instance states.
  ! Flows     :: Set of flow module instance states.
  ! Specieses :: Collection of species.
  ! iS        ::
  ! Source    :: source.
  ! SizeDists :: Collection of particle size distributions.
  ! Mass      ::
  ! iP        ::
  ! DispState :: Dispersion model state.
  ! Locals:
  Type(FlowMemory_)        :: FlowMemory
  Type(Flow_)              :: Flow
  Type(Position_)          :: OldPosition
  Type(Position_)          :: Position
  Real(Std)                :: H3
  Integer                  :: ErrorCode
  Type(Particle_), Pointer :: Particle
  Type(Extra_),    Pointer :: Extra
  Logical                  :: Reflected
  Type(Cloud_)             :: Cloud   !} Dummy variables for arguments in GetAttrib.
  Type(Rain_)              :: Rain    !}
  Type(Surface_)           :: Surface !}
  Type(Soil_)              :: Soil    !}
  Type(Plant_)             :: Plant   !}

  Particle => DispState%Particles(iP)
  Extra    => DispState%ParticleExtras(iE)

  If (Source%dZMetres) Then
    Call ResetFlowMemory(Flows, FlowMemory)
    Position = X2Position(Coords, Source%X, Source%iHCoord, Source%iZCoord)
    Call CalcdZdZ(                                                        &
           Coords, Grids, Domains, Source%iZCoord, DispState%iZCoordMagl, &
           Time, .false., ZeroShortTime(), Position,                      & ! $$
           Units, Mets, Flows,                                            &
           FlowMemory,                                                    &
           H3,                                                            &
           ErrorCode                                                      &
         )
    H3 = 1.0/H3
  Else
    H3 = 1.0
  End If

  Call InitParticle(                                                      &
         iUP,                                                             &
         Time,                                                            &
         DispOpts%sSkewTime, DispOpts%sVelMemTime, DispOpts%sInhomogTime, &
         DispOpts%sMVelMemTime,                                           &
         Coords, Specieses, SizeDists,                                    &
         iS, Source, H3, Diameter, Density, ParticleShape, X, dX,         &
         Particle, Extra, iP                                              &
       )
  DispState%ParticleMasses(1:DispState%nParticleSpecieses, iP) = Mass(1:DispState%nParticleSpecieses)
  ! Move to release? $$

  Call ResetFlowMemory(Flows, FlowMemory)

  ! Reflect.
  OldPosition = X2Position(Coords, Particle%XOld, Particle%iHCoord, Particle%iZCoord)
  Position    = X2Position(Coords, Particle%X,    Particle%iHCoord, Particle%iZCoord)
  Call Reflect(                                        &
         Coords, Grids, Domains,                       &
         Time, .false., TravelTime(Particle),          & ! $$
         Extra%VelMem, OldPosition, Position, Extra%U, &
         Units, Mets, Flows,                           &
         FlowMemory,                                   & ! $$ should reset
         Reflected,                                    &
         ErrorCode                                     & ! $$ test error
       )
  ! Note reflect routine must work even if Xold = X (see building esp)

  If (Source%NoReflect .and. Reflected) Then
    ! These first two calls are needed because the particle will attempt a single time step in the
    ! current code structure (the test for particle being 'really active' is after the first step).
    ! Perhaps should change this, but this should work for now. Possible precision problems?
    ! Perhaps should also ensure DiffusionProcessReflect works for X below ground (could fail with attempt
    ! to access Flow%ZInterface(0) due to precision) or kill particles properly before oneparticletimestep. $$
    Call ResetFlowMemory(Flows, FlowMemory)
    Call Position2XUnknownCoord(                          &
           Coords, Position,                              &
           Particle%iHCoord, Particle%iZCoord, Particle%X &
         )

    ! Kill Extra and reassign particle to Extra(0).
    If (DispState%iParticleExtras(iP) > 0) Then
      If (                                                        &
        Extra%FMass == 0.0                                  .and. &
        .not.Extra%VelMem                                   .and. &
        .not.Extra%Inhomog                                  .and. &
        .not.Extra%MVelMem                                  .and. &
        Extra%TPlus  == ZeroShortTime()                     .and. &
        Extra%TMinus == ZeroShortTime()                     .and. &
        Extra%TInst  == InfPastShortTime(Interval = .true.)       &
      ) Then
        ! Return extra to the stack of free extras and reset number of active extras.
        DispState%FreeParticleExtraStack(DispState%nParticleExtras) = DispState%iParticleExtras(iP)
        DispState%nParticleExtras = DispState%nParticleExtras - 1
        ! Assigne particle to Extra(0).
        DispState%iParticleExtras(iP) = 0
      End If
    End If

    Call MarkParticle(UseForOutput = .false., Particle = Particle)
    Return
  End If

  If (Reflected) Call ResetFlowMemory(Flows, FlowMemory)

  Call Position2XUnknownCoord(                          &
         Coords, Position,                              &
         Particle%iHCoord, Particle%iZCoord, Particle%X &
       )

!  Position = X2Position(Coords, Particle%X, Particle%iHCoord, Particle%iZCoord)
!  (Might be needed for bit reproducible results)

  If (Extra%VelMem .or. Extra%MVelMem .or. Extra%FMass /= 0.0) Then ! test on real dangerous $$

    ! Get flow info.
    Call GetAttrib(                              &
           iAttribParam  = A_Flow,               &
           Coords        = Coords,               &
           Grids         = Grids,                &
           Domains       = Domains,              &
           Moisture      = .true.,               & ! $$
           Inhomog       = Extra%Inhomog,        &
           Homog         = .not.Extra%Inhomog,   &
           Time          = Time,                 &
           AnyTravelTime = .false.,              &
           TravelTime    = TravelTime(Particle), &
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
      If (.not.DispState%PPsLostDueToFlow) Then
        DispState%PPsLostDueToFlow = .true.
        Call ControlledMessage(                                                             &
               'Particles and/or puffs are being lost due to lack of a valid flow module.', &
               MessageControls    = GlobalMessageControls,                                  &
               MessageControlName = 'No flow for particle/puff',                            &
               ErrorCode          = 2                                                       &
             )
      End If
     ! If (ErrorCode == 1) Then
     !   Call Message('Error: no flow module instance suitable for ' // &
     !                'vertical coordinate conversion', 2)
     ! Else
     !   Call Message('Error: no flow module instance suitable for ' // &
     !                'supplying the flow information', 2)
     ! End If

      ! Kill Extra and reassign particle to Extra(0).
      If (DispState%iParticleExtras(iP) > 0) Then
        If (                                                        &
          Extra%FMass == 0.0                                  .and. &
          .not.Extra%VelMem                                   .and. &
          .not.Extra%Inhomog                                  .and. &
          .not.Extra%MVelMem                                  .and. &
          Extra%TPlus  == ZeroShortTime()                     .and. &
          Extra%TMinus == ZeroShortTime()                     .and. &
          Extra%TInst  == InfPastShortTime(Interval = .true.)       &
        ) Then
          ! Return extra to the stack of free extras and reset number of active extras.
          DispState%FreeParticleExtraStack(DispState%nParticleExtras) = DispState%iParticleExtras(iP)
          DispState%nParticleExtras = DispState%nParticleExtras - 1
          ! Assigne particle to Extra(0).
          DispState%iParticleExtras(iP) = 0
        End If
      End If

      Call MarkParticle(UseForOutput = .false., Particle = Particle)
      Return
    End If
    If (IsBackwards()) Call ReverseFlow(Flow)

  End If

  ! Initialise particle velocity.
  Call InitParticle2(Flow, Particle, Extra, iP)

  ! Kill Extra and reassign particle to Extra(0).
  If (DispState%iParticleExtras(iP) > 0) Then
    If (                                                        &
      Extra%FMass == 0.0                                  .and. &
      .not.Extra%VelMem                                   .and. &
      .not.Extra%Inhomog                                  .and. &
      .not.Extra%MVelMem                                  .and. &
      Extra%TPlus  == ZeroShortTime()                     .and. &
      Extra%TMinus == ZeroShortTime()                     .and. &
      Extra%TInst  == InfPastShortTime(Interval = .true.)       &
    ) Then
      ! Return extra to the stack of free extras and reset number of active extras.
      DispState%FreeParticleExtraStack(DispState%nParticleExtras) = DispState%iParticleExtras(iP)
      DispState%nParticleExtras = DispState%nParticleExtras - 1
      ! Assigne particle to Extra(0).
      DispState%iParticleExtras(iP) = 0
    End If
  End If

End Subroutine ReleaseParticle

!-------------------------------------------------------------------------------------------------------------

Subroutine ReleasePuff(                                      &
             TimeM, Time, TimeP,                             &
             DispOpts, Units,                                &
             Coords, Grids, Domains, Mets, Flows, Specieses, &
             iS, Source, SizeDists, Mass,                    &
             iUP, iP, iE, iUOP, iOP, DispState               &
           )
! Releases a puff.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In)            :: TimeM
  Type(ShortTime_), Intent(In)            :: Time
  Type(ShortTime_), Intent(In)            :: TimeP
  Type(DispOpts_),  Intent(In)            :: DispOpts
  Type(Units_),     Intent(InOut)         :: Units
  Type(Coords_),    Intent(In)            :: Coords
  Type(Grids_),     Intent(In)            :: Grids
  Type(Domains_),   Intent(In)            :: Domains
  Type(Mets_),      Intent(InOut)         :: Mets
  Type(Flows_),     Intent(InOut)         :: Flows
  Type(Specieses_), Intent(In)            :: Specieses
  Integer,          Intent(In)            :: iS
  Type(Source_),    Intent(In)            :: Source
  Type(SizeDists_), Intent(In)            :: SizeDists
  Real(Std),        Intent(In)            :: Mass(MaxSpecieses)
  Integer,          Intent(InOut)         :: iUP  ! $$ Why inout - change to In after checking
  Integer,          Intent(InOut)         :: iP
  Integer,          Intent(InOut)         :: iE
  Integer,          Intent(InOut)         :: iUOP
  Integer,          Intent(InOut)         :: iOP
  Type(DispState_), Intent(InOut), Target :: DispState
  ! TimeM     ::
  ! Time      ::
  ! TimeP     ::
  ! DispOpts  :: Collection of sets of dispersion model options.
  ! Units     :: Collection of information on input/output unit numbers.
  ! Coords    :: Collection of coord systems.
  ! Grids     :: Collection of grids.
  ! Domains   :: Collection of domains.
  ! Mets      :: Set of met module instance states.
  ! Flows     :: Set of flow module instance states.
  ! Specieses :: Collection of species.
  ! iS        ::
  ! Source    :: source.
  ! SizeDists :: Collection of particle size distributions.
  ! Mass      ::
  ! iUP       :: Unique index of puff
  ! iP        ::
  ! DispState :: Dispersion model state.
  ! Locals:
  Type(FlowMemory_)         :: FlowMemory
  Type(Flow_)               :: Flow
  Type(Position_)           :: Position
  Real(Std)                 :: H3
  Integer                   :: ErrorCode
  Type(Puff_),      Pointer :: Puff
  Type(Particle_),  Pointer :: Particle
  Type(Extra_),     Pointer :: Extra
  Real(Std)                 :: Temp
  Type(Cloud_)              :: Cloud   !} Dummy variables for arguments in GetAttrib.
  Type(Rain_)               :: Rain    !}
  Type(Surface_)            :: Surface !}
  Type(Soil_)               :: Soil    !}
  Type(Plant_)              :: Plant   !}

  Puff     => DispState%Puffs(iP)
  Particle => DispState%Puffs(iP)%P
  Extra    => DispState%PuffExtras(iE)

  If (.not.Source%dZMetres) Then
    Call ResetFlowMemory(Flows, FlowMemory)
    Position = X2Position(Coords, Source%X, Source%iHCoord, Source%iZCoord)
    Call CalcdZdZ(                                                        &
           Coords, Grids, Domains, Source%iZCoord, DispState%iZCoordMagl, &
           Time, .false., ZeroShortTime(), Position,                      & ! $$
           Units, Mets, Flows,                                            &
           FlowMemory,                                                    &
           H3,                                                            &
           ErrorCode                                                      &
         )
    H3 = 1.0/H3
  Else
    H3 = 1.0
  End If

  Call InitPuff(                                                          &
         iUP,                                                             &
         TimeM, Time, TimeP, TimeM == TimeP,                              &
         DispOpts%sSkewTime, DispOpts%sVelMemTime, DispOpts%sInhomogTime, &
         DispOpts%sMVelMemTime,                                           &
         Coords, Specieses,                                               &
!         SizeDists,                                                       &
         iS, Source, H3,                                                  &
         iUOP, iOP,                                                       &
         Puff, Extra                                                      &
       )
  DispState%PuffMasses(1:DispState%nParticleSpecieses, iP) = Mass(1:DispState%nParticleSpecieses)
  ! Move to release? $$

  Call ResetFlowMemory(Flows, FlowMemory)

  ! True centroid location, after reflection.
  If (Puff%XXp(3) == 0.0) Then
    Temp = Puff%P%X(3)
  Else
    Temp = Puff%P%X(3) / Sqrt(2.0*Puff%XXp(3))
    Temp = Puff%P%X(3)*Erf(Temp) + Sqrt(2.0*Puff%XXp(3)/Pi)*Exp(-Temp**2)
  End If
  Position = X2Position(Coords, (/ Particle%X(1:2), Temp /), Particle%iHCoord, Particle%iZCoord)

  ! Get flow info.
  Call GetAttrib(                              &
         iAttribParam  = A_Flow,               &
         Coords        = Coords,               &
         Grids         = Grids,                &
         Domains       = Domains,              &
         Moisture      = .true.,               & ! $$
         Inhomog       = Extra%Inhomog,        &
         Homog         = .not.Extra%Inhomog,   &
         Time          = Time,                 &
         AnyTravelTime = .false.,              &
         TravelTime    = TravelTime(Particle), &
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
    If (.not.DispState%PPsLostDueToFlow) Then
      DispState%PPsLostDueToFlow = .true.
      Call ControlledMessage(                                                             &
             'Particles and/or puffs are being lost due to lack of a valid flow module.', &
             MessageControls    = GlobalMessageControls,                                  &
             MessageControlName = 'No flow for particle/puff',                            &
             ErrorCode          = 2                                                       &
           )
    End If
   ! If (ErrorCode == 1) Then
   !   Call Message('Error: no flow module instance suitable for ' // &
   !                'vertical coordinate conversion', 2)
   ! Else
   !   Call Message('Error: no flow module instance suitable for ' // &
   !                'supplying the flow information', 2)
   ! End If
    Call MarkParticle(UseForOutput = .false., Particle = Particle)
    Return
  End If
  If (IsBackwards()) Call ReverseFlow(Flow)

  ! Do second initialisation of puff.
  Call InitPuff2(DispOpts%DeltaOpt, Flow, Puff, Extra, iP)

End Subroutine ReleasePuff

!-------------------------------------------------------------------------------------------------------------

Subroutine OneParticleTimeStep(      &
             iParticle,              &
             TargetTime,             &
             Coords, Grids, Domains, &
             Units, Mets, Flows,     &
             Specieses,              &
             DispOpts,               &
             Flow, Cloud, Rain,      &
             Surface, Soil, Plant,   &
             DryDep, WetDep, RDt,    &
             DispState               &
           )
! Time-step particles.

  Implicit None
  ! Argument list:
  Integer,          Intent(In)            :: iParticle            !
  Type(ShortTime_), Intent(In)            :: TargetTime           !
  Type(Coords_),    Intent(In)            :: Coords               !
  Type(Grids_),     Intent(In)            :: Grids                !
  Type(Domains_),   Intent(In)            :: Domains              !
  Type(Units_),     Intent(InOut)         :: Units
  Type(Mets_),      Intent(InOut)         :: Mets
  Type(Flows_),     Intent(InOut)         :: Flows                !
  Type(Specieses_), Intent(In)            :: Specieses            !
  Type(DispOpts_),  Intent(In)            :: DispOpts             !
  Type(Flow_),      Intent(Out)           :: Flow                 !
  Type(Cloud_),     Intent(Out)           :: Cloud                !
  Type(Rain_),      Intent(Out)           :: Rain                 !
  Type(Surface_),   Intent(Out)           :: Surface              !
  Type(Soil_),      Intent(Out)           :: Soil                 !
  Type(Plant_),     Intent(Out)           :: Plant                !
  Real(Std),        Intent(Out)           :: DryDep(MaxSpecieses) !
  Real(Std),        Intent(Out)           :: WetDep(MaxSpecieses) !
  Real(Std),        Intent(Out)           :: RDt                  !
  Type(DispState_), Intent(InOut), Target :: DispState            !
  ! Locals:
  Type(Particle_),  Pointer :: Particle    ! An abbreviation for the particle being considered.
  Type(Extra_),     Pointer :: Extra
  Type(Position_)           :: OldPosition !
  Type(Position_)           :: Position    !
  Type(FlowMemory_)         :: FlowMemory  !
  Type(ShortTime_)          :: SDt         !
  Type(ShortTime_)          :: SDtLim      !
  Real(Std)                 :: RDt2
  Real(Std)                 :: HMax            ! Max metric coefficient.
  Real(Std)                 :: H1          !
  Real(Std)                 :: H2          !
  Real(Std)                 :: Zs          ! Height below which dry deposit
  Real(Std)                 :: DrydepFrac  ! Fraction of timestep below Zs
  Real(Std)                 :: DtUsed      !
  Real(Std)                 :: WMeanAndTurb    !
  Type(ShortTime_)          :: Time        !
  Real(Std)                 :: ZenithAngle ! Solar zenith angle at particle location.
  Real(Std)                 :: RelHumidity ! Relative humidity at particle location.
  Real(Std)                 :: XLatLong(3) ! Position in Lat Long coord system
                                           ! (land use dry deposition scheme)
  Real(Std)                 :: XPa(3)      ! Particle position in Pa for DeepConvNew
  Integer                   :: ErrorCode   !
  Logical                   :: In
  Logical                   :: Reflected

  Type(Position_)           :: PositionTemp    !
  Real Temp
  Logical                   :: UseLandUseDryDep
  Integer                   :: iParticleSpecies
  Integer                   :: iZCoordPa   ! Index of the vertical coord system 'Pa' in
                                           ! the collection of all vertical coord systems

  !$OMP ATOMIC
  DispState%nParticleTimeSteps = DispState%nParticleTimeSteps + 1

  Particle => DispState%Particles(iParticle)
  Extra    => DispState%ParticleExtras(DispState%iParticleExtras(iParticle))

  Time = ParticleTime(Particle)

  Call ResetFlowMemory(Flows, FlowMemory)

  Position = X2Position(Coords, Particle%X, Particle%iHCoord, Particle%iZCoord)

    ! True centroid location, after reflection.
!    Temp = Particle%X(3) / Sqrt(2.0*10000.0)
!    Temp = Particle%X(3)*Erf(Temp) + Sqrt(2.0*10000.0/Pi)*Exp(-Temp**2)
!    PositionTemp = X2Position(Coords, (/ Particle%X(1:2), Temp /), Particle%iHCoord, Particle%iZCoord)
  ! Get flow info.
  Call GetAttrib(                              &
         iAttribParam  = A_Flow,               &
         Coords        = Coords,               &
         Grids         = Grids,                &
         Domains       = Domains,              &
         Moisture      = .true.,               & ! $$
         Inhomog       = Extra%Inhomog,        &
         Homog         = .not.Extra%Inhomog,   &
         Time          = Time,                 &
         AnyTravelTime = .false.,              &
         TravelTime    = TravelTime(Particle), &
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
 ! Flow%Sk = Flow%Sk * Particle%X(3) / Temp
 ! Flow%dSkdz = Flow%dSkdz * Particle%X(3) / Temp
  If (ErrorCode > 0) Then
    If (.not.DispState%PPsLostDueToFlow) Then
      DispState%PPsLostDueToFlow = .true.
      Call ControlledMessage(                                                             &
             'Particles and/or puffs are being lost due to lack of a valid flow module.', &
             MessageControls    = GlobalMessageControls,                                  &
             MessageControlName = 'No flow for particle/puff',                            &
             ErrorCode          = 2                                                       &
           )
    End If
   ! If (ErrorCode == 1) Then
   !   Call Message('Error: no flow module instance suitable for ' // &
   !                'vertical coordinate conversion', 2)
   ! Else
   !   Call Message('Error: no flow module instance suitable for ' // &
   !                'supplying the flow information', 2)
   ! End If
    Call MarkParticle(UseForOutput = .false., Particle = Particle)
    Return
  End If
  If (IsBackwards()) Call ReverseFlow(Flow)

  ! Get cloud info if necessary.
  If ((DispOpts%DeepConvectionCode /= Convection_None) .or. DispOpts%DryDep .or. DispOpts%WetDep) Then
    Call GetAttrib(                              &
           iAttribParam  = A_Cloud,              &
           Coords        = Coords,               &
           Grids         = Grids,                &
           Domains       = Domains,              &
           Moisture      = .false.,              &
           Inhomog       = .false.,              &
           Homog         = .false.,              &
           Time          = Time,                 &
           AnyTravelTime = .false.,              &
           TravelTime    = TravelTime(Particle), &
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
      If (ErrorCode == 3) Then
        Call Message(                                                          &
               'FATAL ERROR: no cloud attribute defined in the collection ' // &
               'of flow module instance attributes',                           &
               3                                                               &
             )
      Else If (.not.DispState%PPsLostDueToFlow) Then
        DispState%PPsLostDueToFlow = .true.
        Call ControlledMessage(                                                             &
               'Particles and/or puffs are being lost due to lack of a valid flow module.', &
               MessageControls    = GlobalMessageControls,                                  &
               MessageControlName = 'No flow for particle/puff',                            &
               ErrorCode          = 2                                                       &
             ) ! $$ Should use 'No cloud for particle/puff' control
     ! Else If (ErrorCode == 2)
     !   Call Message('Error: no flow module instance suitable for ' // &
     !                'supplying the cloud information', 2)
     ! Else
     !   Call Message('Error: no flow module instance suitable for ' // &
     !                'vertical coordinate conversion', 2)
      End If
      Call MarkParticle(UseForOutput = .false., Particle = Particle)
      Return
    End If
  End If

  ! Get rain info if necessary.
  If (DispOpts%WetDep .or. DispOpts%AgentDecay .or. (DispOpts%DeepConvectionCode /= Convection_None)) Then
    Call GetAttrib(                              &
           iAttribParam  = A_Rain,               &
           Coords        = Coords,               &
           Grids         = Grids,                &
           Domains       = Domains,              &
           Moisture      = .false.,              &
           Inhomog       = .false.,              &
           Homog         = .false.,              &
           Time          = Time,                 &
           AnyTravelTime = .false.,              &
           TravelTime    = TravelTime(Particle), &
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
      If (ErrorCode == 3) Then
        Call Message(                                                         &
               'FATAL ERROR: no rain attribute defined in the collection ' // &
               'of flow module instance attributes',                          &
               3                                                              &
             )
      Else If (.not.DispState%PPsLostDueToFlow) Then
        DispState%PPsLostDueToFlow = .true.
        Call ControlledMessage(                                                             &
               'Particles and/or puffs are being lost due to lack of a valid flow module.', &
               MessageControls    = GlobalMessageControls,                                  &
               MessageControlName = 'No flow for particle/puff',                            &
               ErrorCode          = 2                                                       &
             ) ! $$ Should use 'No rain for particle/puff' control
     ! Else If (ErrorCode == 2)
     !   Call Message('Error: no flow module instance suitable for ' // &
     !                'supplying the rain information', 2)
     ! Else
     !   Call Message('Error: no flow module instance suitable for ' // &
     !                'vertical coordinate conversion', 2)
      End If
      Call MarkParticle(UseForOutput = .false., Particle = Particle)
      Return
    End If
  End If

  ! Calculate zenith angle and relative humidity (used in depletion calculations).
  !$$ Zenith angle calculation needs a calendar time frame.
  If (DispOpts%AgentDecay) Then
    ZenithAngle = CalcZenithAngle(        &
                    Coords,               &
                    ShortTime2Time(Time), &
                    Position              &
                  )
    RelHumidity = CalcRH(Flow%Q, Flow%T, Flow%P)
  End If

  ! Calculate time step for random walk.
  Call CalcRandomWalkDt(Flow, Particle, Extra, RDt)

!!! May Want to switch to this if very slow
!!! IF (ShortTime2RealTime(TravelTime(Particle)).GT.(2.24.3600.)) RDt=Huge(1.)

  SDt    = RealTime2ShortTime(RDt)

  SDtLim = TargetTime - Time
  If (SDt > SDtLim) Then
    SDt = SDtLim
    RDt = ShortTime2RealTime(SDt)
  End If

!  WRITE(6,*)'Particle time step (s):',RDt

  ! Convert to right coord system for turbulent advection.
  ! $$ Velocities, sigmas too using Particle%iZCoord to Flow%iZCoord
  Call ConvertToH(Coords, Flow%iHCoord, Position)
  !  If (iConvert == 0) Then         ! $$ This is needed on NAG HP Compiler without
  !    Write (6, *) 're1', iConvert  ! gline option. Not sure why. Without, code stops
  !  End If                          ! in ConvertToZ with 'Error in ConvertZ: flow
                                    ! module not found'
  Call ConvertToZ(                                      &
         Coords, Grids, Domains,                        &
         Flow%iZCoord,                                  &
         Time, .false., TravelTime(Particle), Position, &
         Units, Mets, Flows,                            &
         FlowMemory, ErrorCode                          &
       ) ! $$ probably can't fail here but should probably check error code
  Particle%X = Position2X(Coords, Position, Flow%iHCoord, Flow%iZCoord)
  Particle%iHCoord = Flow%iHCoord
  Particle%iZCoord = Flow%iZCoord

  ! Get surface info if necessary (for new land use dependent dry deposition scheme).
  If (DispOpts%DryDep) Then

    UseLandUseDryDep = .false. ! $$ could precalculate this
    Do iParticleSpecies = 1, Specieses%nParticleSpecieses
      If (Specieses%Specieses( Specieses%iParticle2Species(iParticleSpecies) )%LandUseDryDep) Then
        UseLandUseDryDep = .true.
        Exit
      End If
    End Do
   
    If (UseLandUseDryDep) Then

      Call GetAttrib(                              &
             iAttribParam  = A_Surface,            &
             Coords        = Coords,               &
             Grids         = Grids,                &
             Domains       = Domains,              &
             Moisture      = .false.,              &
             Inhomog       = .false.,              &
             Homog         = .false.,              &
             Time          = Time,                 &
             AnyTravelTime = .false.,              &
             TravelTime    = TravelTime(Particle), &
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
        If (ErrorCode == 3) Then
          Call Message(                                                            &
                 'FATAL ERROR: no surface attribute defined in the collection ' // &
                 'of flow module instance attributes',                             &
                 3                                                                 &
               )
        Else If (.not.DispState%PPsLostDueToFlow) Then
          DispState%PPsLostDueToFlow = .true.
          Call Message(                                                            &
                 'FATAL ERROR: no surface information available',                  &
                 3                                                                 &
               )
!          Call ControlledMessage(                                                             &
!                 'Particles and/or puffs are being lost due to lack of a valid flow module.', &
!                 MessageControls    = GlobalMessageControls,                                  &
!                 MessageControlName = 'No flow for particle/puff',                            &
!                 ErrorCode          = 2                                                       &
!               ) ! $$ Should use 'No surface for particle/puff' control
       ! Else If (ErrorCode == 2)
       !   Call Message('Error: no flow module instance suitable for ' // &
       !                'supplying the rain information', 2)
       ! Else
       !   Call Message('Error: no flow module instance suitable for ' // &
       !                'vertical coordinate conversion', 2)
        End If
        Call MarkParticle(UseForOutput = .false., Particle = Particle)
        Return
      End If

      Call GetAttrib(                              &
             iAttribParam  = A_Plant,              &
             Coords        = Coords,               &
             Grids         = Grids,                &
             Domains       = Domains,              &
             Moisture      = .false.,              &
             Inhomog       = .false.,              &
             Homog         = .false.,              &
             Time          = Time,                 &
             AnyTravelTime = .false.,              &
             TravelTime    = TravelTime(Particle), &
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
        If (ErrorCode == 3) Then
          Call Message(                                                          &
                 'FATAL ERROR: no plant attribute defined in the collection ' // &
                 'of flow module instance attributes',                           &
                 3                                                               &
               )
        Else If (.not.DispState%PPsLostDueToFlow) Then
          DispState%PPsLostDueToFlow = .true.
          Call Message(                                                            &
                 'FATAL ERROR: no plant information available',                    &
                 3                                                                 &
               )
!          Call ControlledMessage(                                                             &
!                 'Particles and/or puffs are being lost due to lack of a valid flow module.', &
!                 MessageControls    = GlobalMessageControls,                                  &
!                 MessageControlName = 'No flow for particle/puff',                            &
!                 ErrorCode          = 2                                                       &
!               ) ! $$ Should use 'No surface for particle/puff' control
       ! Else If (ErrorCode == 2)
       !   Call Message('Error: no flow module instance suitable for ' // &
       !                'supplying the rain information', 2)
       ! Else
       !   Call Message('Error: no flow module instance suitable for ' // &
       !                'vertical coordinate conversion', 2)
        End If
        Call MarkParticle(UseForOutput = .false., Particle = Particle)
        Return
      End If
    End If
  End If

  ! Set Zs
  Zs = Min(DispOpts%Zs, Flow%H)

  ! Calculate sedimentation velocity
  Particle%WSed = TerminalVelocity(                                                                        &
                    Particle%Diameter, Particle%Density, Particle%ParticleShape, Particle%ShapeSchemeCode, &
                    Flow%T, Flow%P, Flow%Rho                                                               &
                  )

  ! Modify turbulence parameters for trajectory-crossing effect and inertial effect
  If (Particle%WSed > 0.0) Then

    Call ModifyTurbByTerminalVelocity(Particle%WSed, Particle%X(3), Extra%Inhomog, Flow)

    ! Extra restrictions on time-step $$ Note that Flow%MaxDZ and previous RDt / SDt are calculated
    !                                 $$ based on unmodified sigmas, taus, Ks etc.
    RDt2 = Huge(RDt2)
    If (Extra%VelMem) Then
      If (Flow%MaxDZ /= Huge(Flow%MaxDZ)) Then
        RDt2 = Min(RDt2, 0.1_Std * Flow%MaxDZ / Max(Abs(Flow%U(3) + Extra%U(3) - Particle%WSed), &
               Sqrt(Flow%SigUU(3))))
      End If
      RDt2 = Min(RDt2, 0.1_Std * Flow%TauUU(3))
    Else
      If (Extra%Inhomog) Then
        If (Flow%MaxDZ /= Huge(Flow%MaxDZ)) RDt2 = Min(RDt2, 0.03_Std * Flow%MaxDZ**2 / (2.0 * Flow%K(3)))
      Else
        If (Flow%MaxDZ /= Huge(Flow%MaxDZ)) RDt2 = Min(RDt2, 0.03_Std * Flow%MaxDZ**2 / (2.0 * Flow%HK(3)))
      End If
    End If


    If (RDt2 < RDt) Then
      SDt = RealTime2ShortTime(RDt2)
      RDt = RDt2
    End If

  End If

  ! Evolve particle.
  Call MetricCoeffs(Coords%HCoords(Particle%iHCoord), Particle%X(1:2), HMax, H1, H2)
  Call EvolveParticle(                                                                                &
         Flow, RDt, H1, H2,                                                                           &
         DispOpts%Turbulence, DispOpts%MesoscaleMotions, DispOpts%Damping, DispOpts%VerticalVelocity, &
         Zs, Coords%HCoords(Particle%iHCoord), WMeanAndTurb, Particle, Extra, iParticle               &
       )

  ! Correct position for any particles beyond coordinate system pole.
  Call RationaliseHPosition(Particle%X(1:2),Coords%HCoords(Particle%iHCoord))

  ! Reflect.
  If (Extra%Inhomog) Then
  !  OldPosition = X2Position(Coords, Particle%XOld, Particle%iHCoord, Particle%iZCoord)
  !  Position    = X2Position(Coords, Particle%X,    Particle%iHCoord, Particle%iZCoord)
  !  Call Reflect(                                        &
  !         Coords, Grids, Domains,                       &
  !         Time, .false., TravelTime(Particle),          & ! $$
  !         Extra%VelMem, OldPosition, Position, Extra%U, &
  !         Units, Mets, Flows,                           &
  !         FlowMemory,                                   & ! $$ should reset
  !         Reflected,                                    &
  !         ErrorCode                                     & ! $$ test error
  !       )
  !  Particle%X = Position2X(Coords, Position, Particle%iHCoord, Particle%iZCoord)

    ! New dry deposition changes
    If (Particle%WSed > 0.0) Then

      If (Particle%XOld(3) > Zs) Then

        If (Particle%X(3) > Zs) Then
          DrydepFrac =0.0
        Else
          DtUsed = (Zs - Particle%XOld(3))/(WMeanAndTurb - Particle%WSed)
          If (Particle%X(3) > 0.0) Then
            DrydepFrac = (RDt - DtUsed) / RDt
          Else
            If (Abs(Particle%WSed / (Particle%WSed - WMeanAndTurb)) > 0.01) Then
              DrydepFrac = - Zs / (RDt * Particle%WSed) *                                       &
                             ALog(- WMeanAndTurb / (Particle%WSed - WMeanAndTurb))
              DtUsed = DtUsed - Zs / Particle%WSed *                                            &
                         ALog(- WMeanAndTurb / (Particle%WSed - WMeanAndTurb))
            Else
              DrydepFrac = Zs / (RDt * (Particle%WSed - WMeanAndTurb))
              DtUsed = DtUsed + Zs / (Particle%WSed - WMeanAndTurb)
            End If
           If (Abs(Particle%WSed * (RDt - DtUsed) / Zs) > 0.01) Then
             Particle%X(3) = -  WMeanAndTurb * Zs / Particle%WSed *                             &
                               (1.0 - Exp(-Particle%WSed * (RDt - DtUsed) / Zs))
            Else
              Particle%X(3) = WMeanAndTurb * (DtUsed - RDt)
            End If
            If (Extra%VelMem) Extra%U(3) = ReflectW(Flow, 1.0, Extra%Skew, Extra%U(3))
            If (Particle%X(3) < Zs) Then
              DrydepFrac = DrydepFrac + (RDt - DtUsed) / RDt
            Else
              If (Abs(Particle%WSed / WMeanAndTurb) > 0.01) Then
                DrydepFrac = DrydepFrac - Zs / (RDt * Particle%WSed) *                          &
                               ALog((Particle%WSed + WMeanAndTurb) / WMeanAndTurb)
                DtUsed = DtUsed - Zs / Particle%WSed *                                          &
                           ALog((Particle%WSed + WMeanAndTurb) / WMeanAndTurb)
              Else
                DrydepFrac = DrydepFrac - Zs / (RDt * WMeanAndTurb)
                DtUsed = DtUsed - Zs / WMeanAndTurb
              End If
              Particle%X(3) = Zs - (Particle%WSed + WMeanAndTurb) * (RDt - DtUsed)
            End If
          End If
        End If

      Else

        If (Particle%X(3) > Zs) Then
          If (Abs((Zs - Particle%XOld(3)) * Particle%WSed /                                     &
                (Particle%WSed * Particle%XOld(3) - Zs * WMeanAndTurb)) > 0.01) Then
            DtUsed = - Zs / Particle%WSed * ALog(Zs * (Particle%WSed - WMeanAndTurb) /          &
                       (Particle%XOld(3) * Particle%WSed - Zs * WMeanAndTurb))
          Else
            DtUsed = Zs *(Particle%XOld(3) - Zs) /                                              &
                     (Particle%XOld(3) * Particle%WSed - Zs * WMeanAndTurb)
          End If
          DrydepFrac = DtUsed / RDt
        Else
          If (Particle%X(3) > 0.0) Then
            DrydepFrac = 1.0
          Else
            If (Abs(Particle%XOld(3) * Particle%WSed /                                          &
                  (Particle%XOld(3) * Particle%WSed - Zs * WMeanAndTurb)) > 0.01) Then
              DtUsed = - Zs / Particle%WSed * ALog(- Zs * WMeanAndTurb /                        &
                         (Particle%XOld(3) * Particle%WSed - Zs * WMeanAndTurb))
            Else
              DtUsed = Zs * Particle%XOld(3) /                                                  &
                         (Particle%XOld(3) * Particle%WSed - Zs * WMeanAndTurb)
            End If
            DrydepFrac = DtUsed / RDt
            If (Abs(Particle%WSed * (RDt - DtUsed) / Zs) > 0.01) Then
              Particle%X(3) = - WMeanAndTurb * Zs / Particle%WSed *                             &
                                (1.0 - Exp (- Particle%WSed * (RDt - DtUsed) / Zs))
            Else
              Particle%X(3) = WMeanAndTurb * (DtUsed - RDt)
            End If
            If (Extra%VelMem) Extra%U(3) = ReflectW(Flow, 1.0, Extra%Skew, Extra%U(3))
            If (Particle%X(3) < Zs) Then
              DrydepFrac = 1.0
            Else
              If (Abs(Particle%WSed / WMeanAndTurb) > 0.01) Then
                DrydepFrac = DrydepFrac - Zs / (RDt * Particle%WSed) *                          &
                               ALog((Particle%WSed + WMeanAndTurb) / WMeanAndTurb)
                DtUsed = DtUsed - Zs / Particle%WSed *                                          &
                           ALog((Particle%WSed + WMeanAndTurb) / WMeanAndTurb)
              Else
                DrydepFrac = DrydepFrac - Zs / (RDt * WMeanAndTurb)
                DtUsed = DtUsed - Zs /  WMeanAndTurb
              End If
              Particle%X(3) = Zs - (Particle%WSed + WMeanAndTurb) * (RDt - DtUsed)
            End If
          End If
        End If

      End If

    Else

      If (Particle%X(3) > 0.0) Then

        If ( Abs( Particle%X(3) - Particle%XOld(3) ) /= 0.0 ) Then
          DrydepFrac = Abs( Min(Zs, Particle%XOld(3)) - Min(Zs, Particle%X(3)) )
          DrydepFrac = DrydepFrac / Abs( Particle%X(3) - Particle%XOld(3) )
        Else If ( Particle%XOld(3) < Zs ) Then
          DrydepFrac = 1.0
        Else
          DrydepFrac = 0.0
        End If

      Else

        If ( Abs( Particle%X(3) - Particle%XOld(3) ) /= 0.0 ) Then
          DrydepFrac = Min(Zs, Particle%XOld(3)) + Min(Zs, -Particle%X(3))
          DrydepFrac = DrydepFrac / Abs( Particle%X(3) - Particle%XOld(3) )
        Else If ( Particle%XOld(3) < Zs ) Then
          DrydepFrac = 1.0
        Else
          DrydepFrac = 0.0
        End If
        Particle%X(3) = - Particle%X(3)
        If (Extra%VelMem) Extra%U(3) = ReflectW(Flow, 1.0, Extra%Skew, Extra%U(3))

      End If

    End If
    ! End of New Dry Deposition changes

    ! Might fail if reflect changes coord - use Position2XUnknownCoord instead $$
  Else
    Call DiffusionProcessReflect(Particle, Flow, RDt, Zs, WMeanAndTurb, DrydepFrac)
    ! $$ follow diffusion reflect by proper reflect, to avoid particle inside building.
    ! or call reflect from inside diffusion reflect.
  End If
  ! (i) If particle now inside building (because it entered a building domain during
  ! the step) there is a problem.
  ! (ii) Also problem of multiple building effects regions - particle could be
  ! in one b effects region (outside the building) but in another building!
  ! Best to post process output to remove concs in buildings, and ensure (i) doesn't
  ! happen often or at all (e.g. big enough b effects region, small time steps (could
  ! maintain list of building locations in various coord systems and use to limit dt,
  ! independently of which flow module being used), use b effects region only for
  ! sources within region).
  ! If can happen sometimes, still need fixup - set flow to zero in next step and then
  ! do crude reflection, randomising within building effects region minus building?
  ! (doing this next step will reduce cost)

  ! Radioactive decay, agent decay (UV, FMD virus, Midges) dry and wet deposition.
  ! Note that the call to Deep convection has been moved after dry and wet deposition to avoid 
  ! reconverting particle position from pressure to metres.   
  If (DispOpts%RadioactiveDecay) Then
    Call RadioactiveDecay(DispState%ParticleMasses(:, iParticle), Specieses, RDt, Field = .false.)
  End If
  If (DispOpts%AgentDecay) Then
    Call AgentDecay(                                                                        &
           ShortTime2Time(Time), DispState%ParticleMasses(:, iParticle), Specieses, RDt,    &
           Rain, RelHumidity, ZenithAngle, Flow%T, Field = .False., TravelTime = Particle%T &
         )
  End If
  If (DispOpts%DryDep) Then
    If (DrydepFrac > 1.1 .or. Drydepfrac < -0.1) Then ! $$
      Call Message(                                                          &
        'ERROR: DrydepFrac out of range ' // Trim(Std2Char(DrydepFrac)), 2)
    End If
    DrydepFrac = Min(DryDepFrac, 1.0)
    DrydepFrac = Max(DryDepFrac, 0.0)

    ! Calculate particle position in lat long so can test if y>60 deg N for tundra
    ! in land use dependent dry deposition scheme
    Position = X2Position(Coords, Particle%X, Particle%iHCoord, Particle%iZCoord)

    Call ConvertToH(Coords, DispState%iHCoordLatLong, Position)

    XLatLong = Position2X(Coords, Position, DispState%iHCoordLatLong, Particle%iZCoord)

    Call DryDeposition(                                              &
           Specieses, Flow, Cloud, RDt, XLatLong(2),                 &
           Surface, Plant, DryDep, Particle,                         &
           DispState%ParticleMasses(:, iParticle), Zs, DrydepFrac    &
         )
  Else
    DryDep = 0.0
  End If
  If (DispOpts%WetDep) Then
    Call WetDeposition(Specieses, Flow, Cloud, Rain, RDt, WetDep, Particle, &
                       DispState%ParticleMasses(:, iParticle))
  Else
    WetDep = 0.0
  End If

! Deep Convection
  If (DispOpts%DeepConvectionCode == Convection_Old) Then
    Call DeepConvectionOld(Cloud, RDt, Particle, iParticle)
  Else If ( (DispOpts%DeepConvectionCode == Convection_New) .And. (Rain%ConPpt > 0.0) ) Then
    If ( (Cloud%ConCloudBasePa > -1).And.(Cloud%ConCloudTopPa > -1) ) Then
      ! To find particle position in Pa for DeepConvectionNew subroutine.  $$ move iZCoordPa to DispState
      iZCoordPa = FindZCoordIndex('Pa', Coords)
      Call ResetFlowMemory(Flows, FlowMemory)
      Position = X2Position(Coords, Particle%X, Particle%iHCoord, Particle%iZCoord)
      Call ConvertToZ(                                   &
          Coords, Grids, Domains,                        &
          iZCoordPa,                                     &
          Time, .false., TravelTime(Particle), Position, &
          Units, Mets, Flows,                            &
          FlowMemory, ErrorCode                          &
        )
      If (ErrorCode == 0) Then ! $$ an error can probably only occur if particle leaves domain in earlier
                               ! advection step, in which case ignoring the particle is fine.
        XPa = Position2X(Coords, Position, Flow%iHCoord, iZCoordPa)
        Call DeepConvectionNew(Cloud, RDt, iParticle, Rain, Flow, XPa)
        ! Now return the updated particle position to Particle%X(3)
        Particle%X(3) = XPa(3)
        Particle%iZCoord = iZCoordPa
      End If
    End If
  End If

  ! Update particle time and flag indicating whether we use x- or xu-Markov model.
  Call UpdateTime(                                                        &
         DispOpts%sSkewTime, DispOpts%sVelMemTime, DispOpts%sInhomogTime, &
         DispOpts%sMVelMemTime,                                           &
         SDt,                                                             &
         Particle, Extra                                                  &
       )

  ! Kill Extra and reassign particle to Extra(0).
  !$OMP CRITICAL (EXTRAFREE)
  If (DispState%iParticleExtras(iParticle) > 0) Then
    If (                                                        &
      Extra%FMass == 0.0                                  .and. &
      .not.Extra%VelMem                                   .and. &
      .not.Extra%Inhomog                                  .and. &
      .not.Extra%MVelMem                                  .and. &
      Extra%TPlus  == ZeroShortTime()                     .and. &
      Extra%TMinus == ZeroShortTime()                     .and. &
      Extra%TInst  == InfPastShortTime(Interval = .true.)       & ! $$ could check ParticleMinT > TInst
                                                                  ! Would help in creating particles
                                                                  ! from puffs
    ) Then
      ! Return extra to the stack of free extras and reset number of active extras.
      DispState%FreeParticleExtraStack(DispState%nParticleExtras) = DispState%iParticleExtras(iParticle)
      DispState%nParticleExtras = DispState%nParticleExtras - 1
      ! Assigne particle to Extra(0).
      DispState%iParticleExtras(iParticle) = 0
    End If
  End If
  !$OMP END CRITICAL (EXTRAFREE)

End Subroutine OneParticleTimeStep

!-------------------------------------------------------------------------------------------------------------

Subroutine OnePuffTimeStep(          &
             iPuff,                  &
             TargetTime,             &
             Coords, Grids, Domains, &
             Units, Mets, Flows,     &
             Specieses,              &
             DispOpts,               &
             Flow, Cloud, Rain,      &
             Surface, Soil, Plant,   &
             DryDep, WetDep, RDt,    &
             DispState               &
           )
! Time-step puffs.

  Implicit None
  ! Argument list:
  Integer,          Intent(In)            :: iPuff                !
  Type(ShortTime_), Intent(In)            :: TargetTime           !
  Type(Coords_),    Intent(In)            :: Coords               !
  Type(Grids_),     Intent(In)            :: Grids                !
  Type(Domains_),   Intent(In)            :: Domains              !
  Type(Units_),     Intent(InOut)         :: Units
  Type(Mets_),      Intent(InOut)         :: Mets
  Type(Flows_),     Intent(InOut)         :: Flows                !
  Type(Specieses_), Intent(In)            :: Specieses            !
  Type(DispOpts_),  Intent(In)            :: DispOpts             !
  Type(Flow_),      Intent(Out)           :: Flow                 !
  Type(Cloud_),     Intent(Out)           :: Cloud                !
  Type(Rain_),      Intent(Out)           :: Rain                 !
  Type(Surface_),   Intent(Out)           :: Surface              !
  Type(Soil_),      Intent(Out)           :: Soil                 !
  Type(Plant_),     Intent(Out)           :: Plant                !
  Real(Std),        Intent(Out)           :: DryDep(MaxSpecieses) !
  Real(Std),        Intent(Out)           :: WetDep(MaxSpecieses) !
  Real(Std),        Intent(Out)           :: RDt                  !
  Type(DispState_), Intent(InOut), Target :: DispState            !
  ! Locals:
  Type(Puff_),      Pointer :: Puff        ! An abbreviation for the puff being considered.
  Type(Particle_),  Pointer :: Particle    ! An abbreviation for the particle being considered.
  Type(Extra_),     Pointer :: Extra
  Type(Position_)           :: OldPosition !
  Type(Position_)           :: Position    !
  Type(FlowMemory_)         :: FlowMemory  !
  Type(ShortTime_)          :: SDt         !
  Type(ShortTime_)          :: SDtLim      !
  Real(Std)                 :: RDt2
  Real(Std)                 :: HMax            ! Max metric coefficient.
  Real(Std)                 :: H1          !
  Real(Std)                 :: H2          !
  Real(Std)                 :: Zs          ! Height below which dry deposit
  Real(Std)                 :: WMeanAndTurb ! Mean and Turbulent vertical velocity
  Real(Std)                 :: DrydepFrac  ! Fraction of timestep below Zs
  Type(ShortTime_)          :: Time        !
  Real(Std)                 :: ZenithAngle ! Solar zenith angle at puff centre.
  Real(Std)                 :: RelHumidity ! Relative humidity at puff centre.
  Real(Std)                 :: XLatLong(3)    ! position in Lat Long coord system
                                              ! (land use dry deposition scheme)
  Real(Std)                 :: Delta       ! Limit on puff size
  Logical                   :: SigPGtDelta !
  Real(Std)                 :: SigA2       !
  Logical                   :: PuffMoved   !
  Integer                   :: ErrorCode
  Logical                   :: In
  Integer                   :: iOP
  Logical                   :: Reflected
  Real(Std)                 :: Temp         !
  Real(Std)                 :: DtUsed       !
  Type(Position_)           :: PositionTemp ! Position at puff centroid height

  DryDep = 0.0 ! $$
  WetDep = 0.0 ! $$

  DispState%nPuffTimeSteps = DispState%nPuffTimeSteps + 1

  Puff     => DispState%Puffs(iPuff)
  Particle => DispState%Puffs(iPuff)%P
  Extra    => DispState%PuffExtras(iPuff)

  Time = ParticleTime(Particle)

  Call ResetFlowMemory(Flows, FlowMemory)

  Position = X2Position(Coords, Particle%X, Particle%iHCoord, Particle%iZCoord)

  iOP   = Puff%OriginalPuff

  If (DispState%SigA2Defined(iOP) < 2) Then
    SigA2 = Puff%XXh(3)
  Else
    SigA2 = DispState%SigA2(iOP) +                                        &
            (DispState%SigA2(iOP) - DispState%SigA2Last(iOP)) *           &
            ShortTime2RealTime(ParticleTime(Puff%P) - DispState%T(iOP)) / &
            ShortTime2RealTime(DispState%T(iOP) - DispState%TLast(iOP))
  End If

  ! Approximate by Gaussian.
  Call Step4(Puff)

  ! True centroid location, after reflection.
  If (Puff%XXp(3) == 0.0) Then
    Temp = Puff%P%X(3)
  Else
    Temp = Puff%P%X(3) / Sqrt(2.0*Puff%XXp(3))
    Temp = Puff%P%X(3)*Erf(Temp) + Sqrt(2.0*Puff%XXp(3)/Pi)*Exp(-Temp**2)
  End If
  PositionTemp = X2Position(Coords, (/ Particle%X(1:2), Temp /), Particle%iHCoord, Particle%iZCoord)

  ! Get flow info.
  Call GetAttrib(                              &
         iAttribParam  = A_Flow,               &
         Coords        = Coords,               &
         Grids         = Grids,                &
         Domains       = Domains,              &
         Moisture      = .true.,               & ! $$
         Inhomog       = Extra%Inhomog,        &
         Homog         = .not.Extra%Inhomog,   &
         Time          = Time,                 &
         AnyTravelTime = .false.,              &
         TravelTime    = TravelTime(Particle), &
         Position      = PositionTemp,         &
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
    If (.not.DispState%PPsLostDueToFlow) Then
      DispState%PPsLostDueToFlow = .true.
      Call ControlledMessage(                                                             &
             'Particles and/or puffs are being lost due to lack of a valid flow module.', &
             MessageControls    = GlobalMessageControls,                                  &
             MessageControlName = 'No flow for particle/puff',                            &
             ErrorCode          = 2                                                       &
           )
    End If
   ! If (ErrorCode == 1) Then
   !   Call Message('Error: no flow module instance suitable for ' // &
   !                'vertical coordinate conversion', 2)
   ! Else
   !   Call Message('Error: no flow module instance suitable for ' // &
   !                'supplying the flow information', 2)
   ! End If
    Call MarkParticle(UseForOutput = .false., Particle = Particle)
    Return
  End If
  If (IsBackwards()) Call ReverseFlow(Flow)

  ! Get cloud info if necessary.
  If ((DispOpts%DeepConvectionCode /= Convection_None) .or. DispOpts%DryDep .or. DispOpts%WetDep) Then
    Call GetAttrib(                              &
           iAttribParam  = A_Cloud,              &
           Coords        = Coords,               &
           Grids         = Grids,                &
           Domains       = Domains,              &
           Moisture      = .false.,              &
           Inhomog       = .false.,              &
           Homog         = .false.,              &
           Time          = Time,                 &
           AnyTravelTime = .false.,              &
           TravelTime    = TravelTime(Particle), &
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
      If (ErrorCode == 3) Then
        Call Message(                                                          &
               'FATAL ERROR: no cloud attribute defined in the collection ' // &
               'of flow module instance attributes',                           &
               3                                                               &
             )
      Else If (.not.DispState%PPsLostDueToFlow) Then
        DispState%PPsLostDueToFlow = .true.
        Call ControlledMessage(                                                             &
               'Particles and/or puffs are being lost due to lack of a valid flow module.', &
               MessageControls    = GlobalMessageControls,                                  &
               MessageControlName = 'No flow for particle/puff',                            &
               ErrorCode          = 2                                                       &
             ) ! $$ Should use 'No cloud for particle/puff' control
     ! Else If (ErrorCode == 2)
     !   Call Message('Error: no flow module instance suitable for ' // &
     !                'supplying the cloud information', 2)
     ! Else
     !   Call Message('Error: no flow module instance suitable for ' // &
     !                'vertical coordinate conversion', 2)
      End If
      Call MarkParticle(UseForOutput = .false., Particle = Particle)
      Return
    End If
  End If

  ! Get rain info if necessary.
  If (DispOpts%WetDep .or. (DispOpts%DeepConvectionCode /= Convection_None)) Then
    Call GetAttrib(                              &
           iAttribParam  = A_Rain,               &
           Coords        = Coords,               &
           Grids         = Grids,                &
           Domains       = Domains,              &
           Moisture      = .false.,              &
           Inhomog       = .false.,              &
           Homog         = .false.,              &
           Time          = Time,                 &
           AnyTravelTime = .false.,              &
           TravelTime    = TravelTime(Particle), &
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
      If (ErrorCode == 3) Then
        Call Message(                                                         &
               'FATAL ERROR: no rain attribute defined in the collection ' // &
               'of flow module instance attributes',                          &
               3                                                              &
             )
      Else If (.not.DispState%PPsLostDueToFlow) Then
        DispState%PPsLostDueToFlow = .true.
        Call ControlledMessage(                                                             &
               'Particles and/or puffs are being lost due to lack of a valid flow module.', &
               MessageControls    = GlobalMessageControls,                                  &
               MessageControlName = 'No flow for particle/puff',                            &
               ErrorCode          = 2                                                       &
             ) ! $$ Should use 'No rain for particle/puff' control
     ! Else If (ErrorCode == 2)
     !   Call Message('Error: no flow module instance suitable for ' // &
     !                'supplying the rain information', 2)
     ! Else
     !   Call Message('Error: no flow module instance suitable for ' // &
     !                'vertical coordinate conversion', 2)
      End If
      Call MarkParticle(UseForOutput = .false., Particle = Particle)
      Return
    End If
  End If

  ! Calculate zenith angle and relative humidity (used in depletion calculations).
  !$$ Zenith angle calculation needs a calendar time frame.
  If (DispOpts%AgentDecay) Then
    ZenithAngle = CalcZenithAngle(        &
                    Coords,               &
                    ShortTime2Time(Time), &
                    Position              &
                  )
    RelHumidity = CalcRH(Flow%Q, Flow%T, Flow%P)
  End If

  If (DispOpts%DeltaOpt == 0) Then
    Delta = Huge(Delta)
  Else
    Delta = Flow%DeltaI
  End If
!  If (DispOpts%DeltaOpt == 2) Then
!    Delta = Min(DeltaI, A6*Sqrt(Puff%XXh(3)))
!  Else
!    Delta = DeltaI
!  End If

  ! Reduce sigma_p.
  Call Step5(DispOpts%PuffOpts, DispOpts%DeltaOpt, Delta, Flow, SigA2, SigPGtDelta,&
  PuffMoved, Puff, Extra, iPuff)

    ! $$ temp-fixup. If puff has moved, might be below ground, so need to apply a
    ! reflection. This should use appropriate coords and a flow module reflect routine
    ! (e.g. in a building region the following wont reflect out of the building).
    ! However currently vertical coords are always height above ground and the
    ! building case doesn't seem to cause problems. May need to update Puff%P%XOld and
    ! Puff%P%SigUU too??
    !
    ! Possibly a better long term solution is to ensure step5 can't move the puff a
    ! long way in a single step, and don't bother to call getattrib again (if getattrib
    ! not called again there's no problem in not doing the reflection till later).
    If (Puff%P%X(3) < 0.0) Then
      If (Extra%VelMem) Extra%U(3) = - Extra%U(3)
      Puff%P%X(3) = - Puff%P%X(3)
    End If

  ! Split puffs.
  Call Step6(                                                                                &
         DispOpts%PuffOpts,                                                                  &
         iPuff, DispOpts%DeltaOpt, Delta, SigA2, SigPGtDelta,                                &
         DispState%Puffs, DispState%PuffExtras, DispState%PuffMasses,                        &
         DispState%nPuffs, DispState%LastPuff, DispState%FreePuffStack, DispState%iULastPuff &
       )

  ! $$ Alternative scheme to avoid reducing puffs in size. Comment out step5 to step6 inc
  ! Split puffs.
!  PuffMoved = .false.
!  SigPGtDelta = .false.
!  If (Puff%XXp(3) > Delta**2) Puff%R = .true.
!  Call Step6(                                                          &
!         iPuff, Delta, SigA2, SigPGtDelta,                             &
!         DispState%Puffs, DispState%PuffExtras, DispState%PuffMasses,  &
!         DispState%nPuffs, DispState%LastPuff, DispState%FreePuffStack &
!       )

  If (.not.ParticleActive(Puff%P)) Return

  If (PuffMoved) Then

    ! $$ temp-fixup. If puff has moved, might be below ground, so need to apply a
    ! reflection. This should use appropriate coords and a flow module reflect routine
    ! (e.g. in a building region the following wont reflect out of the building).
    ! However currently vertical coords are always height above ground and the
    ! building case doesn't seem to cause problems. May need to update Puff%P%XOld and
    ! Puff%P%SigUU too??
    !
    ! Possibly a better long term solution is to ensure step5 can't move the puff a
    ! long way in a single step, and don't bother to call getattrib again (if getattrib
    ! not called again there's no problem in not doing the reflection till later).
!    If (Puff%P%X(3) < 0.0) Then
!      Puff%P%U(3) = - Puff%P%U(3)
!      Puff%P%X(3) = - Puff%P%X(3)
!    End If

    Call ResetFlowMemory(Flows, FlowMemory)

    If (Particle%X(3) < 0.0) Then
      Particle%X(3) = - Particle%X(3)
      If (Extra%VelMem) Extra%U(3) = ReflectW(Flow, 1.0 - A6**2, Extra%Skew, Extra%U(3))
    End If

    Position = X2Position(Coords, Particle%X, Particle%iHCoord, Particle%iZCoord)

    ! True centroid location, after reflection.
    If (Puff%XXp(3) == 0.0) Then
      Temp = Puff%P%X(3)
    Else
      Temp = Puff%P%X(3) / Sqrt(2.0*Puff%XXp(3))
      Temp = Puff%P%X(3)*Erf(Temp) + Sqrt(2.0*Puff%XXp(3)/Pi)*Exp(-Temp**2)
    End If
    PositionTemp = X2Position(Coords, (/ Particle%X(1:2), Temp /), Particle%iHCoord, Particle%iZCoord)

    ! Get flow info.
    Call GetAttrib(                              &
           iAttribParam  = A_Flow,               &
           Coords        = Coords,               &
           Grids         = Grids,                &
           Domains       = Domains,              &
           Moisture      = .true.,               & ! $$
           Inhomog       = Extra%Inhomog,        &
           Homog         = .not.Extra%Inhomog,   &
           Time          = Time,                 &
           AnyTravelTime = .false.,              &
           TravelTime    = TravelTime(Particle), &
           Position      = PositionTemp,         &
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
      If (.not.DispState%PPsLostDueToFlow) Then
        DispState%PPsLostDueToFlow = .true.
        Call ControlledMessage(                                                             &
               'Particles and/or puffs are being lost due to lack of a valid flow module.', &
               MessageControls    = GlobalMessageControls,                                  &
               MessageControlName = 'No flow for particle/puff',                            &
               ErrorCode          = 2                                                       &
             )
      End If
     ! If (ErrorCode == 1) Then
     !   Call Message('Error: no flow module instance suitable for ' // &
     !                'vertical coordinate conversion', 2)
     ! Else
     !   Call Message('Error: no flow module instance suitable for ' // &
     !                'supplying the flow information', 2)
     ! End If
      Call MarkParticle(UseForOutput = .false., Particle = Particle)
      Return
    End If
    If (IsBackwards()) Call ReverseFlow(Flow)

  End If

  ! Convert to right coord system for turbulent advection.
  ! $$ Velocities, sigmas too using Particle%iZCoord to Flow%iZCoord
  Call ConvertToH(Coords, Flow%iHCoord, Position)
  !  If (iConvert == 0) Then         ! $$ This is needed on NAG HP Compiler without
  !    Write (6, *) 're1', iConvert  ! gline option. Not sure why. Without, code stops
  !  End If                          ! in ConvertToZ with 'Error in ConvertZ: flow
                                    ! module not found'
  Call ConvertToZ(                                      &
         Coords, Grids, Domains,                        &
         Flow%iZCoord,                                  &
         Time, .false., TravelTime(Particle), Position, &
         Units, Mets, Flows,                            &
         FlowMemory, ErrorCode                          &
       ) ! $$ probably can't fail here but should probably check error code
  Particle%X = Position2X(Coords, Position, Flow%iHCoord, Flow%iZCoord)
  Particle%iHCoord = Flow%iHCoord
  Particle%iZCoord = Flow%iZCoord

  ! Calculate time step and time-step puffs.
  SDtLim = TargetTime - ParticleTime(Puff%P)
  ! $$ Note dt is calculated at a slightly different point to that for particles. This is because it was
  ! at one time thought that it was necessary to impose an additional limit on dt based on max(sigw,w)/a
  ! (this wasn't necessary once A was calculated for skew case using asymptotics where appropriate).
  ! This could be done by moving the next 5 lines into Step7. Ultimately we should consider rationalising dt
  ! strategy - e.g. min step for advection might be different from that for deposition.
  Call CalcRandomWalkDt(Flow, Puff%P, Extra, RDt)
  SDt = RealTime2ShortTime(RDt)
  If (SDt > SDtLim) Then
    SDt = SDtLim
    RDt  = ShortTime2RealTime(SDt)
  End If

  ! Set Zs
  Zs = Flow%H       ! $$ Need to investigate varying Zs for puff scheme using DispOpts%Zs (as for particles)

  ! Calculate sedimentation velocity
  Particle%WSed = TerminalVelocity(                                                                        &
                    Particle%Diameter, Particle%Density, Particle%ParticleShape, Particle%ShapeSchemeCode, &
                    Flow%T, Flow%P, Flow%Rho                                                               &
                  )

  ! Modify turbulence parameters for trajectory-crossing effect and inertial effect
  If (Particle%WSed > 0.0) Then

    Call ModifyTurbByTerminalVelocity(Particle%WSed, Particle%X(3), Extra%Inhomog, Flow)

    ! Extra restrictions on time-step $$ Note that Flow%MaxDZ and previous RDt / SDt are calculated
    !                                 $$ based on unmodified sigmas, taus, Ks etc.
    RDt2 = Huge(RDt2)
    If (Extra%VelMem) Then
      If (Flow%MaxDZ /= Huge(Flow%MaxDZ)) Then
        RDt2 = Min(RDt2, 0.1_Std * Flow%MaxDZ / Max(Abs(Flow%U(3) + Extra%U(3) - Particle%WSed), &
               Sqrt(Flow%SigUU(3))))
      End If
      RDt2 = Min(RDt2, 0.1_Std * Flow%TauUU(3))
    Else
      If (Extra%Inhomog) Then
        If (Flow%MaxDZ /= Huge(Flow%MaxDZ)) RDt2 = Min(RDt2, 0.03_Std * Flow%MaxDZ**2 / (2.0 * Flow%K(3)))
      Else
        If (Flow%MaxDZ /= Huge(Flow%MaxDZ)) RDt2 = Min(RDt2, 0.03_Std * Flow%MaxDZ**2 / (2.0 * Flow%HK(3)))
      End If
    End If

    If (RDt2 < RDt) Then
      SDt = RealTime2ShortTime(RDt2)
      RDt = RDt2
    End If

  End If

  Call MetricCoeffs(Coords%HCoords(Particle%iHCoord), Particle%X(1:2), HMax, H1, H2)
  Call Step7(                                            &
         DispOpts%DeltaOpt, Delta, Flow, RDt,            &
         H1, H2, Zs, WMeanAndTurb,                       &
         DispOpts%Turbulence, DispOpts%MesoscaleMotions, &
         Puff, Extra, iPuff                              &
       )

  ! Correct position for any particles beyond coordinate system pole.
  Call RationaliseHPosition(Particle%X(1:2),Coords%HCoords(Particle%iHCoord))

  ! Reflect.
  If (Extra%Inhomog) Then
 !   OldPosition = X2Position(Coords, Particle%XOld, Particle%iHCoord, Particle%iZCoord)
 !   Position    = X2Position(Coords, Particle%X,    Particle%iHCoord, Particle%iZCoord)
 !   Call Reflect(                                        &
 !          Coords, Grids, Domains,                       &
 !          Time, .false., TravelTime(Particle),          & ! $$
 !          Extra%VelMem, OldPosition, Position, Extra%U, &
 !          Units, Mets, Flows,                           &
 !          FlowMemory,                                   & ! $$ should reset
 !          Reflected,                                    &
 !          ErrorCode                                     & ! $$ test error
 !        )
 !   Particle%X = Position2X(Coords, Position, Particle%iHCoord, Particle%iZCoord)
    ! Might fail if reflect changes coord - use Position2XUnknownCoord instead $$

    !New dry deposition changes
    If (Particle%WSed > 0.0) Then

      If (Particle%XOld(3) > Zs) Then

        If (Particle%X(3) > Zs) Then
          DrydepFrac =0.0
        Else
          DtUsed = (Zs - Particle%XOld(3))/(WMeanAndTurb - Particle%WSed)
          If (Particle%X(3) > 0.0) Then
            DrydepFrac = (RDt - DtUsed) / RDt
          Else
            If (Abs(Particle%WSed / (Particle%WSed - WMeanAndTurb)) > 0.01) Then
              DrydepFrac = - Zs / (RDt * Particle%WSed) *                                       &
                             ALog(- WMeanAndTurb / (Particle%WSed - WMeanAndTurb))
              DtUsed = DtUsed - Zs / Particle%WSed *                                            &
                         ALog(- WMeanAndTurb / (Particle%WSed - WMeanAndTurb))
            Else
              DrydepFrac = Zs / (RDt * (Particle%WSed - WMeanAndTurb))
              DtUsed = DtUsed + Zs / (Particle%WSed - WMeanAndTurb)
            End If
            If (Abs(Particle%WSed * (RDt - DtUsed) / Zs) > 0.01) Then
              Particle%X(3) = - WMeanAndTurb * Zs / Particle%WSed *                             &
                                (1.0 - Exp(-Particle%WSed * (RDt - DtUsed) / Zs))
            Else
              Particle%X(3) = WMeanAndTurb * (DtUsed - RDt)
            End If
            If (Extra%VelMem) Extra%U(3) = ReflectW(Flow, 1.0 - A6**2, Extra%Skew, Extra%U(3))
            If (Particle%X(3) < Zs) Then
              DrydepFrac = DrydepFrac + (RDt - DtUsed) / RDt
            Else
              If (Abs(Particle%WSed / WMeanAndTurb) > 0.01) Then
                DrydepFrac = DrydepFrac - Zs / (RDt * Particle%WSed) *                          &
                               ALog((Particle%WSed + WMeanAndTurb) / WMeanAndTurb)
                DtUsed = DtUsed - Zs / Particle%WSed *                                          &
                           ALog((Particle%WSed + WMeanAndTurb) / WMeanAndTurb)
              Else
                DrydepFrac = DrydepFrac - Zs / (RDt * WMeanAndTurb)
                DtUsed = DtUsed - Zs / WMeanAndTurb
              End If
              Particle%X(3) = Zs - (Particle%WSed + WMeanAndTurb) * (RDt - DtUsed)
            End If
          End If
        End If

      Else

        If (Particle%X(3) > Zs) Then
          If (Abs((Zs - Particle%XOld(3)) * Particle%WSed /                                     &
                (Particle%WSed * Particle%XOld(3) - Zs * WMeanAndTurb)) > 0.01) Then
            DtUsed = - Zs / Particle%WSed * ALog(Zs * (Particle%WSed - WMeanAndTurb) /          &
                       (Particle%XOld(3) * Particle%WSed - Zs * WMeanAndTurb))
          Else
            DtUsed = Zs *(Particle%XOld(3) - Zs) /                                              &
                     (Particle%XOld(3) * Particle%WSed - Zs * WMeanAndTurb)
          End If
          DrydepFrac = DtUsed / RDt
        Else
          If (Particle%X(3) > 0.0) Then
            DrydepFrac = 1.0
          Else
            If (Abs(Particle%XOld(3) * Particle%WSed /                                          &
                  (Particle%XOld(3) * Particle%WSed - Zs * WMeanAndTurb)) > 0.01) Then
              DtUsed = - Zs / Particle%WSed * ALog(- Zs * WMeanAndTurb /                        &
                         (Particle%XOld(3) * Particle%WSed - Zs * WMeanAndTurb))
            Else
              DtUsed = Zs * Particle%XOld(3) /                                                  &
                         (Particle%XOld(3) * Particle%WSed - Zs * WMeanAndTurb)
            End If
            DrydepFrac = DtUsed / RDt
            If (Abs(Particle%WSed * (RDt - DtUsed) / Zs) > 0.01) Then
              Particle%X(3) = - WMeanAndTurb * Zs / Particle%WSed *                             &
                                (1.0 - Exp (- Particle%WSed * (RDt - DtUsed) / Zs))
            Else
              Particle%X(3) = WMeanAndTurb * (DtUsed - RDt)
            End If
            If (Extra%VelMem) Extra%U(3) = ReflectW(Flow, 1.0 - A6**2, Extra%Skew, Extra%U(3))
            If (Particle%X(3) < Zs) Then
              DrydepFrac = 1.0
            Else
              If (Abs(Particle%WSed / WMeanAndTurb) > 0.01) Then
                DrydepFrac = DrydepFrac - Zs / (RDt * Particle%WSed) *                          &
                               ALog((Particle%WSed + WMeanAndTurb) / WMeanAndTurb)
                DtUsed = DtUsed - Zs / Particle%WSed *                                          &
                           ALog((Particle%WSed + WMeanAndTurb) / WMeanAndTurb)
              Else
                DrydepFrac = DrydepFrac - Zs / (RDt * WMeanAndTurb)
                DtUsed = DtUsed - Zs /  WMeanAndTurb
              End If
              Particle%X(3) = Zs - (Particle%WSed + WMeanAndTurb) * (RDt - DtUsed)
            End If
          End If
        End If

      End If

    Else

      If (Particle%X(3) > 0.0) Then

        If ( Abs( Particle%X(3) - Particle%XOld(3) ) /= 0.0 ) Then
          DrydepFrac = Abs( Min(Zs, Particle%XOld(3)) - Min(Zs, Particle%X(3)) )
          DrydepFrac = DrydepFrac / Abs( Particle%X(3) - Particle%XOld(3) )
        Else If ( Particle%XOld(3) < Zs ) Then
          DrydepFrac = 1.0
        Else
          DrydepFrac = 0.0
        End If

      Else

        If ( Abs( Particle%X(3) - Particle%XOld(3) ) /= 0.0 ) Then
          DrydepFrac = Min(Zs, Particle%XOld(3)) + Min(Zs, -Particle%X(3))
          DrydepFrac = DrydepFrac / Abs( Particle%X(3) - Particle%XOld(3) )
        Else If ( Particle%XOld(3) < Zs ) Then
          DrydepFrac = 1.0
        Else
          DrydepFrac = 0.0
        End If
        Particle%X(3) = - Particle%X(3)
        If (Extra%VelMem) Extra%U(3) = ReflectW(Flow, 1.0 - A6**2, Extra%Skew, Extra%U(3))

      End If

    End If
    ! End of New Dry Deposition changes

  Else
    Call DiffusionProcessReflect(Particle, Flow, RDt, Zs, WMeanAndTurb, DrydepFrac)
    ! $$ follow diffusion reflect by proper reflect, to avoid particle inside building.
    ! or call reflect from inside diffusion reflect.
  End If

  ! Deep convection, radioactive decay, agent decay (UV, FMD virus) and deposition.
  If (DispOpts%DeepConvectionCode /= Convection_None) Then
!    Call DeepConvectionNew(Cloud, RDt, iPuff, Rain, Flow, XPa)
    Call DeepConvectionOld(Cloud, RDt, Particle, iPuff)
  End If
  If (DispOpts%RadioactiveDecay) Then
    Call RadioactiveDecay(DispState%PuffMasses(:, iPuff), Specieses, RDt, Field = .false.)
  End If
  If (DispOpts%AgentDecay) Then
    Call AgentDecay(                                                                        &
           ShortTime2Time(Time), DispState%PuffMasses(:, iPuff), Specieses, RDt,            &
           Rain, RelHumidity, ZenithAngle, Flow%T, Field = .False., TravelTime = Particle%T &
         )
  End If
  If (DispOpts%DryDep) Then
!    If (Particle%X(3) < Flow%H) Then ! $$ Better way to calculate DryDepFrac (as for particles)?
!      DryDepFrac = 1.0
!    Else
!      DryDepFrac = 0.0
!    End If
    If (DrydepFrac > 1.1 .or. Drydepfrac < -0.1) Then ! $$
      Call Message(                                                          &
        'ERROR: DrydepFrac out of range ' // Trim(Std2Char(DrydepFrac)), 2)
    End If
    DrydepFrac = Min(DryDepFrac, 1.0)
    DrydepFrac = Max(DryDepFrac, 0.0)


    ! Calculate particle position in lat long so can test if y>60 deg N for tundra
    ! in land use dependent dry deposition scheme
    Position = X2Position(Coords, Particle%X, Particle%iHCoord, Particle%iZCoord)

    Call ConvertToH(Coords, DispState%iHCoordLatLong, Position)

    XLatLong = Position2X(Coords, Position, DispState%iHCoordLatLong, Particle%iZCoord)

    If (Any(Specieses%Specieses(1:Specieses%nSpecieses)%LandUseDryDep)) Then
      Call Message(                                                       &
             'FATAL ERROR: the land use dry deposition scheme '        // &
             'cannot be used with puffs',                                 &
             3                                                            &
           )
    End If

    Call DryDeposition(                                   &
           Specieses, Flow, Cloud, RDt, XLatLong(2),      & ! $$How does DrydepFrac work for the puff scheme?
           Surface, Plant, DryDep, Particle,              &
           DispState%PuffMasses(:, iPuff), Zs, DrydepFrac &  ! $$How does DrydepFrac work for the puff scheme?
         )
  Else
    DryDep = 0.0
  End If
  If (DispOpts%WetDep) Then
    Call WetDeposition(Specieses, Flow, Cloud, Rain, RDt, WetDep, Particle, DispState%PuffMasses(:, iPuff))
  Else
    WetDep = 0.0
  End If

  Call UpdateTime(                                                        &
         DispOpts%sSkewTime, DispOpts%sVelMemTime, DispOpts%sInhomogTime, &
         DispOpts%sMVelMemTime,                                           &
         SDt,                                                             &
         Particle, Extra                                                  &
       )

End Subroutine OnePuffTimeStep

!-------------------------------------------------------------------------------------------------------------

Subroutine Particles2Fields(                                             &
             Coords, Grids, Domains, Sources, Mets, Flows,               &
             FlowMemory, Units, Results, DispState, SizeDists, Specieses &
           )

! Coordinates calculation of Eulerian field concentrations from particles (as source term for Eulerian model)

  Implicit None
  ! Argument list:
  Type(Coords_),             Intent(In)            :: Coords
  Type(Grids_),              Intent(In)            :: Grids
  Type(Domains_),            Intent(In)            :: Domains
  Type(Sources_),            Intent(In)            :: Sources
  Type(Mets_),               Intent(InOut)         :: Mets
  Type(Flows_),              Intent(InOut)         :: Flows
  Type(FlowMemory_),         Intent(InOut)         :: FlowMemory
  Type(Units_),              Intent(InOut)         :: Units
  Type(Results_),            Intent(InOut)         :: Results
  Type(DispState_),          Intent(InOut), Target :: DispState ! $$ for this and some other
                                                                ! variables add comments below
  Type(Specieses_),          Intent(In)            :: Specieses
  Type(SizeDists_),          Intent(In),    Target :: SizeDists

  ! Coords             :: Collection of coord systems.
  ! Grids              :: Collection of grids.
  ! Domains            :: Collection of domains.
  ! Mets               :: Collection of met module instance states.
  ! Flows              :: Collection of flow module instance states.
  ! DispState          :: Information on the state of the dispersion calculation.
  ! Position           :: Coords of the Particle%X in various coord systems in Coords, with flags to indicate
  !                       whether the values are valid.
  ! FlowMemory         :: Flow memory.
  ! Units              :: Collection of information on input/output unit numbers.

  ! Locals:
  Integer                       :: iHCoord    !} Indices of required coord systems.
  Integer                       :: iZCoord    !}
  Integer                       :: iP
  Type(Particle_)               :: Particle2  ! Version of Particle with only certain elements defined,
                                              ! and these being defined in the coord systems indicated by
                                              ! iHCoord and iZCoord.
  Real(Std)                     :: dZdZ       ! Rate of change of coordinate with respect to m agl.
  Integer                       :: ErrorCode  ! Error code for coord conversion.
  Type(Flow_)                   :: Flow       !} Flow, cloud, rain and surface information at the particle
  Type(Particle_),      Pointer :: Particle
  Type(EulerianField_), Pointer :: EulerianField
  Type(Position_)               :: Position
  Integer                       :: iSpecies
  Integer                       :: iParticleSpecies
  Integer                       :: iTracer
  Integer                       :: iParticleSize
  Real(Std)                     :: ParticleMass(Specieses%nFields)
  Integer                       :: iSpeciesUse
  Type(SizeDist_),      Pointer :: SizeDist

  EulerianField => DispState%EulerianField

  Do iP = 1, DispState%LastParticle
    Particle => DispState%Particles(iP)
    If (ParticleReallyActive(Particle)) Then
      If ( TravelTime(Particle) >= Sources%Sources(Particle%iSource)%EulerianTime ) Then

        ParticleMass(:) = 0.0
        Do iParticleSpecies = 1, Specieses%nParticleSpecieses
          iSpecies = Specieses%iParticle2Species(iParticleSpecies)
          iSpeciesUse = Specieses%iSpecies2SpeciesUses(iSpecies)    
          If ( Specieses%SpeciesUseses(iSpeciesUse)%iSizeDist == 0) Then 
            ! Species is a non-sedimenting particle
            iTracer = Specieses%iSpeciesAndSize2Field(iSpecies, 1)
          Else
            ! Species is a sedimenting particle with a range of bin sizes 
            SizeDist => SizeDists%SizeDists(Specieses%SpeciesUseses(iSpeciesUse)%iSizeDist)
            ! Assign particle to bin.
            ! Note this is coded to avoid particles slipping between size ranges due to rounding errors with
            ! iTracer left unassigned.
            ! The ith size range is [d(i), d(i+1)) except for the largest range which is [d(i), d(i+1)],
            ! where d(i) = DiameterRangeBoundary(i), although rounding errors may affect this.
            ! Particles outside the size ranges have material retained on the particle.
            ! $$ Although rounding errors are unlikely to be significant, could, in order to avoid leaving
            ! material on particles when this is not intended, allow a tolerance, or advise that the field size
            ! boundaries be chosen a bit wider, or adopt convention that outer boundaries are interpreted
            ! as 0 and infinity. 
            iTracer = 0
            If (Particle%Diameter <= SizeDist%DiameterRangeBoundary(SizeDist%nSizeRanges + 1)) Then
              Do iParticleSize = SizeDist%nSizeRanges, 1, -1 
                If (Particle%Diameter >= SizeDist%DiameterRangeBoundary(iParticleSize)) Then
                  iTracer = Specieses%iSpeciesAndSize2Field(iSpecies, iParticleSize)
                  Exit
                End If
              End Do
            End If
          End If
          ! Cycle if iTracer = 0. This can be due to species not being held on fields or due to particle size
          ! being outside the range of sizes held on fields.
          If (iTracer == 0) Cycle
          ParticleMass(iTracer) = DispState%ParticleMasses(iParticleSpecies, iP)
          ! Set mass associated with species on iP-th particle to zero
          DispState%ParticleMasses(iParticleSpecies, iP) = 0.0
        End Do

        Call ResetFlowMemory(Flows, FlowMemory)
        Position = X2Position(Coords, Particle%X, Particle%iHCoord, Particle%iZCoord)

        ! Find indices of required coord systems
        iHCoord = Grids%HGrids(EulerianField%iHGrid)%iHCoord
        iZCoord = Grids%ZGrids(EulerianField%iZGridBoundary)%iZCoord

        If (iHCoord /= Particle%iHCoord .or. iZCoord /= Particle%iZCoord) Then
          Call ConvertToH(Coords, iHCoord, Position)
          Call ConvertToZ(                                                        &
                 Coords, Grids, Domains,                                          &
                 iZCoord,                                                         &
                 ParticleTime(Particle), .false., TravelTime(Particle), Position, &
                 Units, Mets, Flows,                                              &
                 FlowMemory,                                                      &
                 ErrorCode                                                        &
               )
          If (ErrorCode > 0) Then
            If (.not.DispState%PPsLostDueToFlow) Then
              DispState%PPsLostDueToFlow = .true.
              Call ControlledMessage(                                                             &
                     'Particles and/or puffs are being lost due to lack of a valid flow module.', &
                     MessageControls    = GlobalMessageControls,                                  &
                     MessageControlName = 'No flow for particle/puff',                            &
                     ErrorCode          = 2                                                       &
                   ) ! $$ Should use 'No conversion for particle/puff' control
            End If ! $$ move message to subroutine where can deal with the variety of error codes and give
                   ! better message.
            Call MarkParticle(UseForOutput = .false., Particle = Particle)
            Return
          End If
          Particle2%X = Position2X(Coords, Position, iHCoord, iZCoord)
        Else
          Particle2%X(:) = Particle%X(:)
        End If
        Particle2%iHCoord = iHCoord
        Particle2%iZCoord = iZCoord

        Call CalcdZdZ(                                                          &
               Coords, Grids, Domains,                                          &
               iZCoord, DispState%iZCoordMagl,                                  &
               ParticleTime(Particle), .false., TravelTime(Particle), Position, &
               Units, Mets, Flows,                                              &
               FlowMemory,                                                      &
               dZdZ,                                                            &
               ErrorCode                                                        &
             )

        Call ParticleEulerianFields(                                                   &
               Particle, Particle2, ParticleMass, 1.0/dZdZ,                            &
               Coords, Grids,                                                          &
               EulerianField%iNew, EulerianField%iHGrid, EulerianField%iZGridBoundary, &
               EulerianField%Concentration                                             &
             )

        ! Kill particle if no longer any mass associatd with iP-th particle
        If ( Sum(DispState%ParticleMasses(:, iP)) == 0.0 ) Then
          Call MarkParticle(UseForOutput = .false., Particle = Particle)
        End If

      End If
    End If
  End Do

End Subroutine Particles2Fields

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcParticleResults(                                                  &
             iCase,                                                              &
             Particle, Mass, Extra, DryDep, WetDep, RDt, TimeCentre,             &
             ReqInfos,                                                           &
             OutputOpts,                                                         &
             Coords, Grids, Domains, Mets, Flows,                                &
             Specieses, CloudGammaParamses, Sources, SizeDists, Reqs,            &
             Position, FlowMemory,                                               &
             Units, Results, DispState                                           &
           )
! Coordinates calculation of type L (depending on Lagrangian information) fields and sets of particle/puff
! information from particles.

  Implicit None
  ! Argument list:
  Integer,                   Intent(In)           :: iCase
  Type(Particle_),           Intent(InOut)        :: Particle
  Real(Std),                 Intent(In)           :: Mass(:)
  Type(Extra_),              Intent(In)           :: Extra
  Real(Std),                 Intent(In)           :: DryDep(:)
  Real(Std),                 Intent(In)           :: WetDep(:)
  Real(Std),                 Intent(In)           :: RDt
  Logical,                   Intent(In)           :: TimeCentre
  Type(ReqInfo_),            Intent(In)           :: ReqInfos
  Type(OutputOpts_),         Intent(In)           :: OutputOpts
  Type(Coords_),             Intent(In)           :: Coords
  Type(Grids_),              Intent(In)           :: Grids
  Type(Domains_),            Intent(In)           :: Domains
  Type(Mets_),               Intent(InOut)        :: Mets
  Type(Flows_),              Intent(InOut)        :: Flows
  Type(Specieses_),          Intent(In)           :: Specieses
  Type(CloudGammaParamses_), Intent(In)           :: CloudGammaParamses
  Type(Sources_),            Intent(In),   Target :: Sources
  Type(SizeDists_),          Intent(In)           :: SizeDists
  Type(Reqs_),               Intent(In),   Target :: Reqs
  Type(Position_),           Intent(InOut)        :: Position
  Type(FlowMemory_),         Intent(InOut)        :: FlowMemory
  Type(Units_),              Intent(InOut)        :: Units
  Type(Results_),            Intent(InOut)        :: Results
  Type(DispState_),          Intent(InOut)        :: DispState ! $$ for this and some other
                                                               ! variables add comments below
  ! iCase              :: Number of case.
  ! Particle           :: Particle.
  ! Extra              :: Set of extra information for the particle.
  ! DryDep             :} Dry and wet deposition loss from the particle (for each species) over the last time
  ! WetDep             :} step.
  ! RDt                :: Time step.
  ! TimeCentre         :: Indicates time centre of the particle has crossed the requirement time.
  ! ReqInfos           :: Information on what requirements are required now.
  ! OutputOpts         :: Output options.
  ! Coords             :: Collection of coord systems.
  ! Grids              :: Collection of grids.
  ! Domains            :: Collection of domains.
  ! Mets               :: Collection of met module instance states.
  ! Flows              :: Collection of flow module instance states.
  ! Specieses          :: Collection of specieses.
  ! CloudGammaParamses :: Collection of sets of cloud gamma parameters.
  ! Sources            :: Collection of sources.
  ! SizeDists          :: Collection of particle size distributions.
  ! Reqs               :: Collection of requirements.
  ! DispState          :: Information on the state of the dispersion calculation.
  ! Position           :: Coords of the Particle%X in various coord systems in Coords, with flags to indicate
  !                       whether the values are valid.
  ! FlowMemory         :: Flow memory.
  ! Units              :: Collection of information on input/output unit numbers.
  ! Results            :: Collection of results.
  ! Locals:
  Type(Source_),    Pointer :: Source     !} Abbreviations.
  Type(FieldReq_),  Pointer :: FieldReq   !}
  Type(PPInfoReq_), Pointer :: PPInfoReq  !}
  Integer                   :: iField     ! Field index.
  Integer                   :: iPPInfo    ! Particle/puff info index.
  Integer                   :: iT         ! Time index.
  Integer                   :: iTA        ! Array index at which results for a given
                                          ! time are stored.
  Integer                   :: iS         ! Travel time index.
  Integer                   :: iHCoord    !} Indices of required coord systems.
  Integer                   :: iZCoord    !}
  Integer                   :: iComponent ! Component index.
  Integer                   :: iRange     ! Size range index.
  Integer                   :: iLast      ! Last index in field array.
  Logical                   :: NeedCalciLast ! Indicates need to call routine CalciLast (could always call
                                             ! CalciLast, but it is more efficient to avoid if not necessary).
  Integer                   :: i          !} Loop indices.
  Integer                   :: j          !}
  Type(Particle_)           :: Particle2  ! Version of Particle with only certain elements defined, and these
                                          ! being defined in the coord systems indicated by iHCoord and
                                          ! iZCoord.
  Real(Std)                 :: dZdZ       ! Rate of change of coordinate with respect to m agl.
  Integer                   :: ErrorCode  ! Error code for coord conversion.
  Type(Flow_)               :: Flow       !} Flow, cloud, rain and surface information at the particle
  Type(Cloud_)              :: Cloud      !} location.
  Type(Rain_)               :: Rain       !}
  Type(Surface_)            :: Surface    !}
  Type(Soil_)               :: Soil       !}
  Type(Plant_)              :: Plant      !}

  ! Fields.
  Do i = 1, ReqInfos%nFields

    iField   =  ReqInfos%iField(i)
    FieldReq => Reqs%FieldReqs(iField)

    If (FieldReq%iSource /= 0 .and. FieldReq%iSource /= Particle%iSource) Cycle

    If (FieldReq%iSourceGroup /= 0) Then
      Source => Sources%Sources(Particle%iSource)
      If (All(FieldReq%iSourceGroup /= Source%iSourceGroups(1:Source%nSourceGroups))) Cycle
    End If ! $$ This could be more efficient for large numbers of source groups per source.

    ! If FieldReq depends on a species which is not carried on particles, ignore FieldReq.
    ! Contribution will be zero.
    If (FieldReq%iSpecies /= 0 .and. FieldReq%iParticleSpecies == 0) Cycle

    If (FieldReq%iParticleSpecies /= 0) Then
      If (Mass(FieldReq%iParticleSpecies) == 0) Cycle
    End If

    NeedCalciLast = .false.
    iComponent = 1  
    If (FieldReq%iSizeDist /= 0) Then
      iRange = 0
      If (Particle%Diameter >= SizeDists%SizeDists(FieldReq%iSizeDist)%DiameterRangeBoundary(1)) Then
        Do j = 2, SizeDists%SizeDists(FieldReq%iSizeDist)%nSizeRanges + 1
          If (Particle%Diameter <= SizeDists%SizeDists(FieldReq%iSizeDist)%DiameterRangeBoundary(j)) Then
            iRange = j - 1
            Exit
          End If
        End Do
      End If
      If (iRange == 0) Cycle
      NeedCalciLast = .true.
    Else
      iRange = 1
    End If
    If (NeedCalciLast) Then
      iLast = CalciLast(FieldReq, iComponent, iRange)
    Else
      iLast = 1
    End If

    ! Find indices of required coord systems (set to current coord system if no coord system needed).
    ! $$ Note this uses the grid's hcoord if /= 0, then hcoord if grid absent. This
    ! copes with Xstats - no grid but require particle in particular coord system.
    ! However in future might need both coord systems.
    If (FieldReq%iHGrid /= 0) Then
      iHCoord = FieldReq%iHGridCoord
    Else If (FieldReq%iHCoord /= 0) Then
      iHCoord = FieldReq%iHCoord
    Else
      iHCoord = Particle%iHCoord
    End If
    If (FieldReq%iZGrid /= 0) Then
      iZCoord = FieldReq%iZGridCoord
    Else If (FieldReq%iZCoord /= 0) Then
      iZCoord = FieldReq%iZCoord
    Else
      iZCoord = Particle%iZCoord
    End If

    If (iHCoord /= Particle%iHCoord .or. iZCoord /= Particle%iZCoord) Then
      Call ConvertToH(Coords, iHCoord, Position)
      Call ConvertToZ(                                                        &
             Coords, Grids, Domains,                                          &
             iZCoord,                                                         &
             ParticleTime(Particle), .false., TravelTime(Particle), Position, &
             Units, Mets, Flows,                                              &
             FlowMemory,                                                      &
             ErrorCode                                                        &
           )
      If (ErrorCode > 0) Then
        If (.not.DispState%PPsLostDueToFlow) Then
          DispState%PPsLostDueToFlow = .true.
          Call ControlledMessage(                                                             &
                 'Particles and/or puffs are being lost due to lack of a valid flow module.', &
                 MessageControls    = GlobalMessageControls,                                  &
                 MessageControlName = 'No flow for particle/puff',                            &
                 ErrorCode          = 2                                                       &
               ) ! $$ Should use 'No conversion for particle/puff' control
        End If ! $$ move message to subroutine where can deal with the variety of error codes and give
               ! better message.
        Call MarkParticle(UseForOutput = .false., Particle = Particle)
        Return
      End If
      Particle2%X = Position2X(Coords, Position, iHCoord, iZCoord)
    Else
      Particle2%X(:) = Particle%X(:)
    End If
    Particle2%iHCoord = iHCoord
    Particle2%iZCoord = iZCoord

    ! May need old position in future, as for puffs. $$

    ! Process t-based requirements.
    If (FieldReq%iSGrid == 0) Then

      iT  = ReqInfos%iT(i)
      iTA = CalciTA(Results%Fields(iField), iT)

      If (FieldReq%iQuantity == Q_AirConc) Then

        Call CalcdZdZ(                                                          &
               Coords, Grids, Domains,                                          &
               iZCoord, DispState%iZCoordMagl,                                  &
               ParticleTime(Particle), .false., TravelTime(Particle), Position, &
               Units, Mets, Flows,                                              &
               FlowMemory,                                                      &
               dZdZ,                                                            &
               ErrorCode                                                        &
             )
        Call ParticleConc(                          &
               Particle, Particle2, Mass, 1.0/dZdZ, &
               Coords, Grids, Reqs,                 &
               iField, iTA, iLast,                  &
               Results                              &
             )

      Else If (FieldReq%iQuantity == Q_DryDep) Then

        Call ParticleConc(                       &
               Particle, Particle2, DryDep, 1.0, &
               Coords, Grids, Reqs,              &
               iField, iTA, iLast,               &
               Results                           &
             )

      Else If (FieldReq%iQuantity == Q_WetDep) Then

        Call ParticleConc(                       &
               Particle, Particle2, WetDep, 1.0, &
               Coords, Grids, Reqs,              &
               iField, iTA, iLast,               &
               Results                           &
             )

      Else If (FieldReq%iQuantity == Q_PhotonFlux) Then

        Call ParticleCloudGamma(                 &
               Particle, Particle2, Mass,        &
               Specieses,                        &
               CloudGammaParamses,               &
               Coords, Grids, Reqs,              &
               iField, iTA,                      &
               Results                           &
             )

      Else If (FieldReq%iQuantity == Q_MeanS) Then

        Call ParticleMeanTravelTime(      &
               Particle, Particle2, Mass, &
               Coords, Grids, Reqs,       &
               iField, iT, iTA, iLast,    &
               Results                    &
             )

      Else If (FieldReq%iQuantity == Q_SigmaZ) Then

        Call ParticleSigZ2(               &
               Particle, Particle2, Mass, &
               Reqs,                      &
               iField, iTA, iLast,        &
               Results                    &
             )

      Else If (FieldReq%iQuantity == Q_MeanZ) Then

        Call ParticleMeanZ(               &
               Particle, Particle2, Mass, &
               Reqs,                      &
               iField, iTA, iLast,        &
               Results                    &
             )

      Else If (FieldReq%iQuantity == Q_Mass) Then

        Call ParticleMass(                &
               Particle, Particle2, Mass, &
               Reqs,                      &
               iField, iTA, iLast,        &
               Results                    &
             )

      Else If (FieldReq%iQuantity == Q_nParticles) Then

        Call ParticleNumbers(       &
               Particle, Particle2, &
               Coords, Grids, Reqs, &
               iField, iTA, iLast,  &
               Results              &
             )

      Else If (FieldReq%iQuantity == Q_nParticlesBySpec) Then

        Call ParticleNumbersBySpecies(    &
               Particle, Particle2, Mass, &
               Coords, Grids, Reqs,       &
               iField, iTA, iLast,        &
               Results                    &
             )

      Else If (FieldReq%iQuantity == Q_nParticleSteps) Then

        Results%Fields(iField)%P64(1, 1, 1, 1, 1, iTA, 1) = DispState%nParticleTimeSteps
        ! $$ no good for source restricted requirements.

      End If

    ! Process s-based requirements.
    Else

      iS = ReqInfos%iS(i)

      If (FieldReq%iQuantity == Q_XStats) Then

        Call ParticleXStats(                    &
               Particle, Particle2, Mass,       &
               Coords, Grids, Reqs, iField, iS, &
               Results                          &
             )

      End If

    End If

  End Do

  ! Sets of particle/puff information.
  If (TimeCentre) Then

    Do i = 1, ReqInfos%nPPInfos

      iPPInfo   =  ReqInfos%iPPInfo(i)
      PPInfoReq => Reqs%PPInfoReqs(iPPInfo)

      If (PPInfoReq%Particles == 'N') Cycle

      If (PPInfoReq%iSource /= 0) Then
        If (Particle%iSource /= PPInfoReq%iSource) Cycle
      End If

      If (PPInfoReq%Particles == 'S') Then
        If (Particle%iUP < PPInfoReq%FirstParticle .or. Particle%iUP > PPInfoReq%LastParticle) Cycle
      End If

      iT = ReqInfos%iPPInfoT(i)

      ! Find indices of required coord systems (set to current coord system if no
      ! coord system needed).
      If (PPInfoReq%iHCoord /= 0) Then
        iHCoord = PPInfoReq%iHCoord
      Else
        iHCoord = Particle%iHCoord
      End If
      If (PPInfoReq%iZCoord /= 0) Then
        iZCoord = PPInfoReq%iZCoord
      Else
        iZCoord = Particle%iZCoord
      End If

      If (iHCoord /= Particle%iHCoord .or. iZCoord /= Particle%iZCoord) Then
        Call ConvertToH(Coords, iHCoord, Position)
        Call ConvertToZ(                                                        &
               Coords, Grids, Domains,                                          &
               iZCoord,                                                         &
               ParticleTime(Particle), .false., TravelTime(Particle), Position, &
               Units, Mets, Flows,                                              &
               FlowMemory,                                                      &
               ErrorCode                                                        &
             )
        If (ErrorCode > 0) Then
          If (.not.DispState%PPsLostDueToFlow) Then
            DispState%PPsLostDueToFlow = .true.
            Call ControlledMessage(                                                             &
                   'Particles and/or puffs are being lost due to lack of a valid flow module.', &
                   MessageControls    = GlobalMessageControls,                                  &
                   MessageControlName = 'No flow for particle/puff',                            &
                   ErrorCode          = 2                                                       &
                 ) ! $$ Should use 'No conversion for particle/puff' control
                   ! $$ move message to subroutine where can deal with the variety of error codes and give
                   ! better message.
                   ! $$ need to actually mark particles
          End If
          Return
        End If
        Particle2%X       = Position2X(Coords, Position, iHCoord, iZCoord)
        Particle2%iHCoord = iHCoord
        Particle2%iZCoord = iZCoord
      Else
        Particle2%X(:)    = Particle%X(:)
        Particle2%iHCoord = Particle%iHCoord
        Particle2%iZCoord = Particle%iZCoord
      End If

      ! Get flow/cloud/rain information if required.
      If (PPInfoReq%Met) Then
        Call GetAttrib(                                &
               iAttribParam  = A_Flow,                 &
               Coords        = Coords,                 &
               Grids         = Grids,                  &
               Domains       = Domains,                &
               Moisture      = .true.,                 & ! $$
               Inhomog       = .true.,                 & ! $$
               Homog         = .true.,                 & ! $$
               Time          = ParticleTime(Particle), &
               AnyTravelTime = .false.,                &
               TravelTime    = TravelTime(Particle),   &
               Position      = Position,               &
               Units         = Units,                  &
               Mets          = Mets,                   &
               Flows         = Flows,                  &
               FlowMemory    = FlowMemory,             &
               Flow          = Flow,                   &
               Cloud         = Cloud,                  &
               Rain          = Rain,                   &
               Surface       = Surface,                &
               Soil          = Soil,                   &
               Plant         = Plant,                  &
               ErrorCode     = ErrorCode               &
             )
        If (ErrorCode > 0) Cycle ! $$ i.e. consider next output request. This treatment is different
                                 ! from above. Both need reviewing.
        If (IsBackwards()) Call ReverseFlow(Flow)
        Call GetAttrib(                                &
               iAttribParam  = A_Cloud,                &
               Coords        = Coords,                 &
               Grids         = Grids,                  &
               Domains       = Domains,                &
               Moisture      = .false.,                &
               Inhomog       = .false.,                &
               Homog         = .false.,                &
               Time          = ParticleTime(Particle), &
               AnyTravelTime = .false.,                &
               TravelTime    = TravelTime(Particle),   &
               Position      = Position,               &
               Units         = Units,                  &
               Mets          = Mets,                   &
               Flows         = Flows,                  &
               FlowMemory    = FlowMemory,             &
               Flow          = Flow,                   &
               Cloud         = Cloud,                  &
               Rain          = Rain,                   &
               Surface       = Surface,                &
               Soil          = Soil,                   &
               Plant         = Plant,                  &
               ErrorCode     = ErrorCode               &
             )
        ! If (ErrorCode > 0) Cycle ! $$ no error check because cloud/rain may not be available
        Call GetAttrib(                                &
               iAttribParam  = A_Rain,                 &
               Coords        = Coords,                 &
               Grids         = Grids,                  &
               Domains       = Domains,                &
               Moisture      = .false.,                &
               Inhomog       = .false.,                &
               Homog         = .false.,                &
               Time          = ParticleTime(Particle), &
               AnyTravelTime = .false.,                &
               TravelTime    = TravelTime(Particle),   &
               Position      = Position,               &
               Units         = Units,                  &
               Mets          = Mets,                   &
               Flows         = Flows,                  &
               FlowMemory    = FlowMemory,             &
               Flow          = Flow,                   &
               Cloud         = Cloud,                  &
               Rain          = Rain,                   &
               Surface       = Surface,                &
               Soil          = Soil,                   &
               Plant         = Plant,                  &
               ErrorCode     = ErrorCode               &
             )
        ! If (ErrorCode > 0) Cycle ! $$ no error check because cloud/rain may not be available
      End If

      Call ParticleInfo(                                     &
             iCase,                                          &
             Particle, Particle2, Extra, Mass,               &
             Flow, Cloud, Rain,                              &
             OutputOpts,                                     &
             Coords, Grids, Flows, Specieses, Sources, Reqs, &
             iPPInfo, iT,                                    &
             Units, Results                                  &
           )

    End Do

  End If

End Subroutine CalcParticleResults

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcPuffResults(                                                                      &
             iCase,                                                                              &
             Puff, Mass, Extra, DryDep, WetDep, RDt, TimeCentre,                                 &
             ReqInfos,                                                                           &
             OutputOpts,                                                                         &
             Coords, Grids, Domains, Mets, Flows, Specieses, Sources, SizeDists, Reqs, DispOpts, &
             Position, FlowMemory,                                                               &
             Units, Results, DispState                                                           &
           )
! Coordinates calculation of type L (depending on Lagrangian information) fields and sets of particle/puff
! information from puffs.

  Implicit None
  ! Argument list:
  Integer,           Intent(In)            :: iCase
  Type(Puff_),       Intent(InOut), Target :: Puff
  Real(Std),         Intent(In)            :: Mass(:)
  Type(Extra_),      Intent(In)            :: Extra
  Real(Std),         Intent(In)            :: DryDep(:)
  Real(Std),         Intent(In)            :: WetDep(:)
  Real(Std),         Intent(In)            :: RDt
  Logical,           Intent(In)            :: TimeCentre
  Type(ReqInfo_),    Intent(In)            :: ReqInfos
  Type(OutputOpts_), Intent(In)            :: OutputOpts
  Type(Coords_),     Intent(In)            :: Coords
  Type(Grids_),      Intent(In)            :: Grids
  Type(Domains_),    Intent(In)            :: Domains
  Type(Mets_),       Intent(InOut)         :: Mets
  Type(Flows_),      Intent(InOut)         :: Flows
  Type(Specieses_),  Intent(In)            :: Specieses
  Type(Sources_),    Intent(In),    Target :: Sources
  Type(SizeDists_),  Intent(In)            :: SizeDists
  Type(Reqs_),       Intent(In),    Target :: Reqs
  Type(DispOpts_),   Intent(In)            :: DispOpts ! A set of dispersion model options.
  Type(Position_),   Intent(InOut)         :: Position
  Type(FlowMemory_), Intent(InOut)         :: FlowMemory
  Type(Units_),      Intent(InOut)         :: Units
  Type(Results_),    Intent(InOut)         :: Results
  Type(DispState_),  Intent(InOut)         :: DispState ! $$ for this and some other variables
                                                        ! add comments below
  ! iCase      :: Number of case.
  ! Puff       :: Puff.
  ! Extra      :: Set of extra information for the puff.
  ! DryDep     :} Dry and wet deposition loss from the puff (for each species) over the last time step.
  ! WetDep     :}
  ! RDt        :: Time step.
  ! TimeCentre :: Indicates time centre of the puff has crossed the requirement time.
  ! ReqInfos   :: Information on what requirements are required now.
  ! OutputOpts :: Output options.
  ! Coords     :: Collection of coord systems.
  ! Grids      :: Collection of grids.
  ! Domains    :: Collection of domains.
  ! Mets       :: Collection of met module instance states.
  ! Flows      :: Collection of flow module instance states.
  ! Specieses  :: Collection of specieses.
  ! Sources    :: Collection of sources.
  ! SizeDists  :: Collection of particle size distributions.
  ! Reqs       :: Collection of requirements.
  ! DispOpts   :: A set of dispersion model options.
  ! DispState  :: Information on the state of the dispersion calculation.
  ! Position   :: Coords of the Puff%P%X in various coord systems in Coords, with flags to indicate whether
  !               the values are valid.
  ! FlowMemory :: Flow memory.
  ! Units      :: Collection of information on input/output unit numbers.
  ! Results    :: Collection of results.
  ! Locals:
  Type(Source_),    Pointer :: Source        !} Abbreviations.
  Type(FieldReq_),  Pointer :: FieldReq      !}
  Type(PPInfoReq_), Pointer :: PPInfoReq     !}
  Type(Particle_),  Pointer :: Particle      ! Abbreviation for particle part of puff.
  Integer                   :: iField        ! Field index.
  Integer                   :: iPPInfo       ! Particle/puff info index.
  Integer                   :: iT            ! Time index.
  Integer                   :: iTA           ! Array index at which results for a given time are stored.
  Integer                   :: iS            ! Travel time index.
  Integer                   :: iHCoord       !} Indices of required coord systems.
  Integer                   :: iZCoord       !}
  Integer                   :: i             ! Loop index.
  Type(FlowMemory_)         :: FlowMemoryOld !} Flow memory and position associated with Puff%P%XOld.
  Type(Position_)           :: PositionOld   !}
  Type(Puff_)               :: Puff2         ! Version of Puff with only certain elements defined, and these
                                             ! being defined in the coord systems indicated by iHCoord and
                                             ! iZCoord.
  Real(Std)                 :: Frac          ! Fraction of puff mass contributing.
  Real(Std)                 :: XOffset(3)    ! Puff position offset due to time spread of puff.
  Real(Std)                 :: dZdZ          ! Rate of change of coordinate with respect to m agl.
  Integer                   :: ErrorCode     ! Error code for coord conversion.
  Type(Flow_)               :: Flow          !} Flow, cloud, rain and surface information at the puff
  Type(Cloud_)              :: Cloud         !} location.
  Type(Rain_)               :: Rain          !}
  Type(Surface_)            :: Surface       !}
  Type(Soil_)               :: Soil          !}
  Type(Plant_)              :: Plant         !}

  Particle => Puff%P

  Call ResetFlowMemory(Flows, FlowMemoryOld)
  PositionOld = X2Position(Coords, Particle%XOld, Particle%iHCoordOld, Particle%iZCoordOld)

  ! Fields.
  Do i = 1, ReqInfos%nFields

    iField   =  ReqInfos%iField(i)
    FieldReq => Reqs%FieldReqs(iField)

    If (FieldReq%iSource /= 0 .and. FieldReq%iSource /= Particle%iSource) Cycle

    If (FieldReq%iSourceGroup /= 0) Then
      Source => Sources%Sources(Particle%iSource)
      If (All(FieldReq%iSourceGroup /= Source%iSourceGroups(1:Source%nSourceGroups))) Cycle
    End If ! $$ This could be more efficient for large numbers of source groups per source.

    ! If FieldReq depends on a species which is not carried on particles, ignore FieldReq.
    ! Contribution will be zero.
    If (FieldReq%iSpecies /= 0 .and. FieldReq%iParticleSpecies == 0) Cycle

    If (FieldReq%iParticleSpecies /= 0) Then
      If (Mass(FieldReq%iParticleSpecies) == 0) Cycle
    End If

    ! Find indices of required coord systems (set to current coord system if no coord system needed).
    ! $$ Note this uses the grid's hcoord if /= 0, then hcoord if grid absent. This
    ! copes with Xstats - no grid but require particle in particular coord system.
    ! However in future might need both coord systems.
    If (FieldReq%iHGrid /= 0) Then
      iHCoord = FieldReq%iHGridCoord
    Else If (FieldReq%iHCoord /= 0) Then
      iHCoord = FieldReq%iHCoord
    Else
      iHCoord = Particle%iHCoord
    End If
    If (FieldReq%iZGrid /= 0) Then
      iZCoord = FieldReq%iZGridCoord
    Else If (FieldReq%iZCoord /= 0) Then
      iZCoord = FieldReq%iZCoord
    Else
      iZCoord = Particle%iZCoord
    End If

    If (iHCoord /= Particle%iHCoord .or. iZCoord /= Particle%iZCoord) Then
      Call ConvertToH(Coords, iHCoord, Position)
      Call ConvertToZ(                                                        &
             Coords, Grids, Domains,                                          &
             iZCoord,                                                         &
             ParticleTime(Particle), .false., TravelTime(Particle), Position, &
             Units, Mets, Flows,                                              &
             FlowMemory,                                                      &
             ErrorCode                                                        &
           )
      If (ErrorCode > 0) Then
        If (.not.DispState%PPsLostDueToFlow) Then
          DispState%PPsLostDueToFlow = .true.
          Call ControlledMessage(                                                             &
                 'Particles and/or puffs are being lost due to lack of a valid flow module.', &
                 MessageControls    = GlobalMessageControls,                                  &
                 MessageControlName = 'No flow for particle/puff',                            &
                 ErrorCode          = 2                                                       &
               ) ! $$ Should use 'No conversion for particle/puff' control
                 ! $$ move message to subroutine where can deal with the variety of error codes and give
                 ! better message.
        End If
        Call MarkParticle(UseForOutput = .false., Particle = Particle)
        Return
      End If
      Puff2%P%X = Position2X(Coords, Position, iHCoord, iZCoord)
    Else
      Puff2%P%X(:) = Puff%P%X(:)
    End If
    Puff2%P%iHCoord = iHCoord
    Puff2%P%iZCoord = iZCoord

    If (iHCoord /= Particle%iHCoordOld .or. iZCoord /= Particle%iZCoordOld) Then
      Call ConvertToH(Coords, iHCoord, PositionOld)
      Call ConvertToZ(                                                    &
             Coords, Grids, Domains,                                      &
             iZCoord,                                                     &
             OldParticleTime(Particle), .false., OldTravelTime(Particle), &
             PositionOld,                                                 &
             Units, Mets, Flows,                                          &
             FlowMemoryOld,                                               &
             ErrorCode                                                    &
           )
      If (ErrorCode > 0) Then
        If (.not.DispState%PPsLostDueToFlow) Then
          DispState%PPsLostDueToFlow = .true.
            Call ControlledMessage(                                                             &
                   'Particles and/or puffs are being lost due to lack of a valid flow module.', &
                   MessageControls    = GlobalMessageControls,                                  &
                   MessageControlName = 'No flow for particle/puff',                            &
                   ErrorCode          = 2                                                       &
                 ) ! $$ Should use 'No conversion for particle/puff' control
                 ! $$ move message to subroutine where can deal with the variety of error codes and give
                 ! better message.
        End If
        Call MarkParticle(UseForOutput = .false., Particle = Particle)
        Return
      End If
      Puff2%P%XOld = Position2X(Coords, PositionOld, iHCoord, iZCoord)
    Else
      Puff2%P%XOld(:) = Puff%P%XOld(:)
    End If
    Puff2%P%iHCoordOld = iHCoord
    Puff2%P%iZCoordOld = iZCoord

    Call PuffFraction(Puff, Extra, ReqInfos%Time, RdT, Frac, XOffset)

    ! Process t-based requirements.
    If (FieldReq%iSGrid == 0) Then

      iT  = ReqInfos%iT(i)
      iTA = CalciTA(Results%Fields(iField), iT)

      If (FieldReq%iQuantity == Q_AirConc) Then

        Call CalcdZdZ(                                                          &
               Coords, Grids, Domains,                                          &
               iZCoord, DispState%iZCoordMagl,                                  &
               ParticleTime(Particle), .false., TravelTime(Particle), Position, &
               Units, Mets, Flows,                                              &
               FlowMemory,                                                      &
               dZdZ,                                                            &
               ErrorCode                                                        &
             )
        Call PuffConc(                                     &
               DispOpts%PuffOpts,                          &
               Puff, Puff2, Mass, Frac, XOffset, 1.0/dZdZ, &
               Coords, Grids, Reqs,                        &
               iField, iTA,                                &
               Results                                     &
             )

      Else If (FieldReq%iQuantity == Q_PuffCentres .and. TimeCentre) Then

        Call PuffCentres(           &
               Puff, Puff2, Mass,   &
               Coords, Grids, Reqs, &
               iField, iTA,         &
               Results              &
             )

      Else If (FieldReq%iQuantity == Q_DryDep) Then

        Call PuffConc(                                  &
               DispOpts%PuffOpts,                       &
               Puff, Puff2, DryDep, Frac, XOffset, 1.0, &
               Coords, Grids, Reqs,                     &
               iField, iTA,                             &
               Results                                  &
             )

      Else If (FieldReq%iQuantity == Q_WetDep) Then

        Call PuffConc(                                  &
               DispOpts%PuffOpts,                       &
               Puff, Puff2, WetDep, Frac, XOffset, 1.0, &
               Coords, Grids, Reqs,                     &
               iField, iTA,                             &
               Results                                  &
             )

      Else If (FieldReq%iQuantity == Q_SigmaZ) Then

        Call PuffSigZ2(                 &
               Puff, Puff2, Mass, Frac, &
               Reqs,                    &
               iField, iTA,             &
               Results                  &
             )

      Else If (FieldReq%iQuantity == Q_MeanZ) Then

        Call PuffMeanZ(                 &
               Puff, Puff2, Mass, Frac, &
               Reqs,                    &
               iField, iTA,             &
               Results                  &
             )

      Else If (FieldReq%iQuantity == Q_Mass) Then

        Call PuffMass(                  &
               Puff, Puff2, Mass, Frac, &
               Reqs,                    &
               iField, iTA,             &
               Results                  &
             )

      Else If (FieldReq%iQuantity == Q_nPuffs .and. TimeCentre) Then

        Results%Fields(iField)%Std(1, 1, 1, 1, 1, iTA, 1) =       &
          Results%Fields(iField)%Std(1, 1, 1, 1, 1, iTA, 1) + 1.0

      Else If (FieldReq%iQuantity == Q_nPuffSteps) Then

        Results%Fields(iField)%P64(1, 1, 1, 1, 1, iTA, 1) = DispState%nPuffTimeSteps
        ! $$ no good for source restricted requirements.

      End If

    ! Process s-based requirements.
    Else

      iS = ReqInfos%iS(i)

      ! $$

    End If

  End Do

  ! Sets of particle/puff information.
  If (TimeCentre) Then

    Do i = 1, ReqInfos%nPPInfos

      iPPInfo   =  ReqInfos%iPPInfo(i)
      PPInfoReq => Reqs%PPInfoReqs(iPPInfo)

      If (PPInfoReq%Puffs == 'N') Cycle

      If (PPInfoReq%iSource /= 0) Then
        If (Puff%P%iSource /= PPInfoReq%iSource) Cycle
      End If

      If (PPInfoReq%Puffs == 'S') Then
        If (Puff%P%iUP < PPInfoReq%FirstPuff .or. Puff%P%iUP > PPInfoReq%LastPuff) Cycle
      End If

      iT = ReqInfos%iPPInfoT(i)

      ! Find indices of required coord systems (set to current coord system if no
      ! coord system needed).
      If (PPInfoReq%iHCoord /= 0) Then
        iHCoord = PPInfoReq%iHCoord
      Else
        iHCoord = Particle%iHCoord
      End If
      If (PPInfoReq%iZCoord /= 0) Then
        iZCoord = PPInfoReq%iZCoord
      Else
        iZCoord = Particle%iZCoord
      End If

      If (iHCoord /= Particle%iHCoord .or. iZCoord /= Particle%iZCoord) Then
        Call ConvertToH(Coords, iHCoord, Position)
        Call ConvertToZ(                                                        &
               Coords, Grids, Domains,                                          &
               iZCoord,                                                         &
               ParticleTime(Particle), .false., TravelTime(Particle), Position, &
               Units, Mets, Flows,                                              &
               FlowMemory,                                                      &
               ErrorCode                                                        &
             )
        If (ErrorCode > 0) Then
          If (.not.DispState%PPsLostDueToFlow) Then
            DispState%PPsLostDueToFlow = .true.
            Call ControlledMessage(                                                             &
                   'Particles and/or puffs are being lost due to lack of a valid flow module.', &
                   MessageControls    = GlobalMessageControls,                                  &
                   MessageControlName = 'No flow for particle/puff',                            &
                   ErrorCode          = 2                                                       &
                 ) ! $$ Should use 'No conversion for particle/puff' control
                   ! $$ move message to subroutine where can deal with the variety of error codes and give
                   ! better message.
                   ! $$ need to actually mark particles
          End If
          Return
        End If
        Puff2%P%X       = Position2X(Coords, Position, iHCoord, iZCoord)
        Puff2%P%iHCoord = iHCoord
        Puff2%P%iZCoord = iZCoord
      Else
        Puff2%P%X(:)    = Puff%P%X(:)
        Puff2%P%iHCoord = Puff%P%iHCoord
        Puff2%P%iZCoord = Puff%P%iZCoord
      End If

      ! Get flow/cloud/rain information if required.
      If (PPInfoReq%Met) Then ! $$ should this be at the reflected puff centroid?
        Call GetAttrib(                                &
               iAttribParam  = A_Flow,                 &
               Coords        = Coords,                 &
               Grids         = Grids,                  &
               Domains       = Domains,                &
               Moisture      = .true.,                 & ! $$
               Inhomog       = .true.,                 & ! $$
               Homog         = .true.,                 & ! $$
               Time          = ParticleTime(Particle), &
               AnyTravelTime = .false.,                &
               TravelTime    = TravelTime(Particle),   &
               Position      = Position,               &
               Units         = Units,                  &
               Mets          = Mets,                   &
               Flows         = Flows,                  &
               FlowMemory    = FlowMemory,             &
               Flow          = Flow,                   &
               Cloud         = Cloud,                  &
               Rain          = Rain,                   &
               Surface       = Surface,                &
               Soil          = Soil,                   &
               Plant         = Plant,                  &
               ErrorCode     = ErrorCode               &
             )
        If (ErrorCode > 0) Cycle ! $$ i.e. consider next output request. This treatment is different
                                 ! from above. Both need reviewing.
        If (IsBackwards()) Call ReverseFlow(Flow)
        Call GetAttrib(                                &
               iAttribParam  = A_Cloud,                &
               Coords        = Coords,                 &
               Grids         = Grids,                  &
               Domains       = Domains,                &
               Moisture      = .false.,                &
               Inhomog       = .false.,                &
               Homog         = .false.,                &
               Time          = ParticleTime(Particle), &
               AnyTravelTime = .false.,                &
               TravelTime    = TravelTime(Particle),   &
               Position      = Position,               &
               Units         = Units,                  &
               Mets          = Mets,                   &
               Flows         = Flows,                  &
               FlowMemory    = FlowMemory,             &
               Flow          = Flow,                   &
               Cloud         = Cloud,                  &
               Rain          = Rain,                   &
               Surface       = Surface,                &
               Soil          = Soil,                   &
               Plant         = Plant,                  &
               ErrorCode     = ErrorCode               &
             )
        ! If (ErrorCode > 0) Cycle ! $$ no error check because cloud/rain may not be available
        Call GetAttrib(                                &
               iAttribParam  = A_Rain,                 &
               Coords        = Coords,                 &
               Grids         = Grids,                  &
               Domains       = Domains,                &
               Moisture      = .false.,                &
               Inhomog       = .false.,                &
               Homog         = .false.,                &
               Time          = ParticleTime(Particle), &
               AnyTravelTime = .false.,                &
               TravelTime    = TravelTime(Particle),   &
               Position      = Position,               &
               Units         = Units,                  &
               Mets          = Mets,                   &
               Flows         = Flows,                  &
               FlowMemory    = FlowMemory,             &
               Flow          = Flow,                   &
               Cloud         = Cloud,                  &
               Rain          = Rain,                   &
               Surface       = Surface,                &
               Soil          = Soil,                   &
               Plant         = Plant,                  &
               ErrorCode     = ErrorCode               &
             )
        ! If (ErrorCode > 0) Cycle ! $$ no error check because cloud/rain may not be available
      End If

      Call PuffInfo(                                         &
             iCase,                                          &
             Puff, Puff2, Extra, Mass,                       &
             Flow, Cloud, Rain,                              &
             OutputOpts,                                     &
             Coords, Grids, Flows, Specieses, Sources, Reqs, &
             iPPInfo, iT,                                    &
             Units, Results                                  &
           )

    End Do

  End If

End Subroutine CalcPuffResults

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcEulerianResults(                                  &
             Time, ReqInfos, Coords, Grids, Domains, Reqs,       &
             Units, Mets, Flows, Results,                        &
             ChemOpts, ChemistryDefn, ChemistryState, DispState, &
             SizeDists, Specieses                                &
           )
! Coordinates calculation of type E (depending on Eulerian information) fields from chemistry fields.
!$$ Need to add background STOCHEM fields here??

  Implicit None
  ! Argument list:
  Type(Time_),           Intent(In)            :: Time
  Type(ReqInfo_),        Intent(In)            :: ReqInfos(:)
  Type(Coords_),         Intent(In)            :: Coords
  Type(Grids_),          Intent(In),    Target :: Grids
  Type(Domains_),        Intent(In)            :: Domains
  Type(Reqs_),           Intent(In),    Target :: Reqs
  Type(Units_),          Intent(InOut)         :: Units
  Type(Mets_),           Intent(InOut)         :: Mets
  Type(Flows_),          Intent(InOut)         :: Flows
  Type(Results_),        Intent(InOut), Target :: Results
  Type(ChemOpts_),       Intent(In)            :: ChemOpts
  Type(ChemistryDefn_),  Intent(In)            :: ChemistryDefn
  Type(ChemistryState_), Intent(In)            :: ChemistryState
  Type(DispState_),      Intent(In)            :: DispState
  Type(SizeDists_),      Intent(In),    Target :: SizeDists
  Type(Specieses_),      Intent(In)            :: Specieses
  ! Time           :: Current synchronisation time.
  ! ReqInfos       :: Requirements information for Eulerian results.
  ! Coords         :: Collection of coord systems.
  ! Grids          :: Collection of grids.
  ! Domains        :: Collection of domains.
  ! Flows          :: Collection of flows.
  ! Reqs           :: Collection of requirements.
  ! Units          :: Collection of information on input/output unit numbers.
  ! Mets           :: Collection of met module instance states.
  ! Flows          :: Collection of flow module instance states.
  ! Results        :: Collection of results.
  ! ChemOpts       :: Chemistry options.
  ! ChemistryDefn  :: A chemistry definition.
  ! ChemistryState :: State of a chemistry calculation.
  ! DispState      :: Information on the state of the dispersion calculation.
  ! Locals:
  Type(HGrid_),     Pointer :: HGridChem     !} Horizontal and vertical grids and interface
  Type(ZGrid_),     Pointer :: ZGridChem     !} levels used by the chemistry scheme.
  Type(ZGrid_),     Pointer :: ZLevelsBdy    !}
  Type(HGrid_),     Pointer :: HGrid         !] Horizontal, vertical and temporal grids
  Type(ZGrid_),     Pointer :: ZGrid         !] of the field requirement.
  Type(TGrid_),     Pointer :: TGrid         !]
  Type(ShortTime_)          :: SyncTime      ! Current synchronisation time (as short time).
  Type(FieldReq_),  Pointer :: FieldReq      !} Field requirements.
  Type(FieldReq_),  Pointer :: FieldReq1     !}
  Type(Field_),     Pointer :: Field         ! Local pointer to a field.
  Real(Std)         :: X(3)               ! Coordinates of output grid point.
  Type(Position_)   :: Position           ! Position of output grid point.
  Type(FlowMemory_) :: FlowMemory         ! Flow memory information.
  Type(ShortTime_)  :: TravelTime         ! Dummy variable for TravelTime argument in ConvertToZ.
  Type(HCoeffs_)    :: HCoeffs            ! Horizontal interpolation coefficients.
  Type(ZCoeffs_)    :: ZCoeffs            ! Vertical interpolation coefficients.
  Logical :: Processed(Reqs%MaxFieldReqs) ! Array of flags to indicate that fields have been processed.
  Integer :: iReqInfo                     ! Index of set of requirements in ReqInfos.
  Integer :: i                            !} Loop indices over field requirements.
  Integer :: j                            !}
  Integer :: iHGrid                       ! Index of horizontal grid for field requirement.
  Integer :: iZGrid                       ! Index of vertical grid for field requirement.
  Integer :: iTGrid                       ! Index of time grid for field requirement.
  Integer :: nX                           !} Size of horizontal grid.
  Integer :: nY                           !}
  Integer :: nZ                           ! Size of vertical grid.
  Integer :: nT                           ! Size of time grid.
  Integer :: nFields                      ! Number of field requirements having these same grids.
  Integer :: iFields(Reqs%MaxFieldReqs)   ! Indices of field requirements having these same grids.
                                          ! (field requirements with the same grids are processed together)
  Integer :: iT                           !} Indices of the field.
  Integer :: iX                           !}
  Integer :: iY                           !}
  Integer :: iZ                           !}
  Integer :: Err                          ! Error code from height conversion.
  Integer :: iXChem                       !} Indices of coincident chemistry gridbox.
  Integer :: iYChem                       !}
  Integer :: iZChem                       !}
  Integer :: iTA                          ! Array index at which results for a given time are stored.
  Integer :: iChemField                   ! Loop index over chemistry field species.
  Integer :: iChem                        ! Index of chemistry field in ChemistryDefn.
  Integer :: iEulerianTracer              ! Index of species on Eulerian field including particle size
  Integer :: nParticleSizes               ! Number of particle sizes - 1 for sedimenting particles; 
                                          ! is equal to 0 for non-sedimenting particles 
  Integer :: kMax, kMin, k
  Real(Std) :: ZBoundary
  Type(SizeDist_), Pointer :: SizeDist

  ! Calculate current synchronisation time as a short time.
  SyncTime = Time2ShortTime(Time)

  ! Loop over each set of requirements in ReqInfos.
  Do iReqInfo = 1, Size(ReqInfos)

    If (ReqInfos(iReqInfo)%NoReqs) Exit

    Processed(:) = .false.

    ! Loop over each field requirement to be processed.
    Do j = 1, ReqInfos(iReqInfo)%nFields

      If (Processed(j)) Cycle

      FieldReq => Reqs%FieldReqs( ReqInfos(iReqInfo)%iField(j) )

      If ( Specieses%iSpeciesAndSize2Field(FieldReq%iSpecies, 1) == 0 ) Cycle

      ! Identify horizontal and vertical grids of the Eulerian solver.
      HGridChem => Grids%HGrids(DispState%EulerianField%iHGrid)
      ZGridChem => Grids%ZGrids(DispState%EulerianField%iZGrid)
      ! Identify vertical boundary levels of the Eulerian solver.
      ZLevelsBdy => Grids%ZGrids(DispState%EulerianField%iZGridBoundary)

      ! Index and size of grids for this field requirement. Note that these grids need not
      ! necessarily be the same as the computational grids used by the chemistry scheme.
      iHGrid = FieldReq%iHGrid
      If (iHGrid /= 0) Then
        HGrid => Grids%HGrids(iHGrid)
      Else
        Call Message('ERROR in CalcEulerianResults: requested field must have a H-Grid', 3)
      End If
      nX = HGrid%nX
      nY = HGrid%nY

      iZGrid = FieldReq%iZGrid
      If (iZGrid /= 0) Then
        ZGrid => Grids%ZGrids(iZGrid)
        nZ = ZGrid%nZ
      Else
        nZ = 1
        If ( FieldReq%iQuantity == Q_EulerianConcentration .or. FieldReq%iQuantity == Q_ChemistryField ) Then
          Call Message('ERROR in CalcEulerianResults: requested field must have a Z-Grid', 3)
        End If
      End If

      iTGrid = FieldReq%iTGrid
      If (iTGrid /= 0) Then
        TGrid => Grids%TGrids(iTGrid)
      Else
        Call Message('ERROR in CalcEulerianResults: requested field must have a T-Grid', 3)
      End If
      nT = TGrid%nT

      ! Determine any fields on the same grids (all such fields will be processed together).
      nFields    = 1
      iFields(1) = j
      Do i = j + 1, ReqInfos(iReqInfo)%nFields

        If (Processed(i)) Cycle

        FieldReq1 => Reqs%FieldReqs( ReqInfos(iReqInfo)%iField(i) )

        If ( Specieses%iSpeciesAndSize2Field(FieldReq1%iSpecies, 1) == 0 ) Cycle

        If (FieldReq1%iHGrid /= iHGrid .or.   &
            FieldReq1%iZGrid /= iZGrid .or.   &
            FieldReq1%iTGrid /= iTGrid) Cycle

        nFields          = nFields + 1
        iFields(nFields) = i

      End Do

      ! Time index of the field.
      iT = ReqInfos(iReqInfo)%iT(j)

      ! Loop over space, extracting chemistry information at each output grid point
      ! (taking the value of each chemistry field from the coincident chemistry gridbox).
      Do iX = 1, nX
      Do iY = 1, nY
      Do iZ = 1, nZ
        ! Calculate coordinates of output grid point in the chemistry scheme coord systems
        ! (converting coords if necessary).
        If ( iZGrid == 0 ) Then
          
          If (HGrid%Unstructured) Then
            X = (/ HGrid%X(iX), HGrid%Y(iX), 0.0 /)
          Else
            X = (/ HGrid%X(iX), HGrid%Y(iY), 0.0 /)
          End If
          
          If ( FieldReq%iHGridCoord /= HGridChem%iHCoord ) Then
            Position = X2Position(Coords, X, FieldReq%iHGridCoord, ZGridChem%iZCoord)
            Call ResetFlowMemory(Flows, FlowMemory)
            Call ConvertToH(Coords, HGridChem%iHCoord, Position)
            X = Position2X(Coords, Position, HGridChem%iHCoord, ZGridChem%iZCoord)
          End If
          
        Else
          
          If (HGrid%Unstructured) Then
            X = (/ HGrid%X(iX), HGrid%Y(iX), ZGrid%Z(iZ) /)
          Else
            X = (/ HGrid%X(iX), HGrid%Y(iY), ZGrid%Z(iZ) /)
          End If
          
          If ( FieldReq%iHGridCoord /= HGridChem%iHCoord .or. FieldReq%iZGridCoord /= ZGridChem%iZCoord ) Then
            Position = X2Position(Coords, X, FieldReq%iHGridCoord, FieldReq%iZGridCoord)
            Call ResetFlowMemory(Flows, FlowMemory)
            Call ConvertToH(Coords, HGridChem%iHCoord, Position)
            Call ConvertToZ(                               &
                   Coords, Grids, Domains,                 &
                   ZGridChem%iZCoord,                      &
                   SyncTime, .true., TravelTime, Position, &
                   Units, Mets, Flows,                     &
                   FlowMemory, Err                         &
                 )
            If (Err /= 0) Then
              Call Message('ERROR in CalcEulerianResults: height conversion error', 3)
            End If
            X = Position2X(Coords, Position, HGridChem%iHCoord, ZGridChem%iZCoord)
          End If
          
        End If

        ! Find coincident chemistry gridbox.
        Call GetHCoeffs(X(1), X(2), HGridChem, HCoeffs)

        iXChem = HCoeffs%iXNear
        iYChem = HCoeffs%iYNear

        ! We now want to do:
        !  Call GetZCoeffs(X(3), ZLevelsBdy, ZCoeffs)
        !  iZChem = ZCoeffs%iZNear
        ! with GetZCoeffs using the specified boundary levels ("z on boundaries") to determine
        ! the nearest grid point, rather than taking the nearest "layer-centre grid point".
        ! But GetZCoeffs doesn't currently do this so we've coded it here in line. 
        ! $$ Longer term should support linked pairs of grids (grid points and boundaries) better
        !    in GridAndDomain.F90. 

        ! Find levels surrounding point (if relevant)
        If ( iZGrid /= 0 ) Then
          kMax = ZLevelsBdy%nZ + 2 ! Note there are ZLevelsBdy%nZ + 1 grid cell boundaries.
          kMin = 0
          Do
            k = kMax - kMin
            If (k <= 1) Exit
            k = kMin + k/2
            If (k == ZLevelsBdy%nZ + 1) Then
              ZBoundary = ZLevelsBdy%Z(k-1) + 0.5 * ZLevelsBdy%AvZ(k-1) 
            Else
              ZBoundary = ZLevelsBdy%Z(k) - 0.5 * ZLevelsBdy%AvZ(k) 
            End If
            If ((X(3) - ZBoundary) * ZLevelsBdy%IncreasingValue > 0.0) Then
              kMin = k
            Else
              kMax = k
            End If
          End Do
          ! Set iZChem.
          iZChem = Max(kMin, 1)
        End If
   
        ! Calculate all chemistry requests at this output grid point.
        Do i = 1, nFields

          FieldReq => Reqs%FieldReqs( ReqInfos(iReqInfo)%iField(iFields(i)) )
          Field    => Results%Fields( ReqInfos(iReqInfo)%iField(iFields(i)) )

          iTA = CalciTA(Field, iT)

          ! Calculate index of species in the list of chemistry fields.
          iEulerianTracer = Specieses%iSpeciesAndSize2Field(FieldReq%iSpecies, 1)
          If (Specieses%SpeciesUseses(Specieses%iSpecies2SpeciesUses(FieldReq%iSpecies))%iSizeDist == 0) Then
            nParticleSizes = 0
          Else
            SizeDist => SizeDists%SizeDists(                                  &
                          Specieses%SpeciesUseses(                            &
                            Specieses%iSpecies2SpeciesUses(FieldReq%iSpecies) &
                          )%iSizeDist                                         &
                        )
            nParticleSizes = SizeDist%nSizeRanges - 1
          End If

          If (FieldReq%iQuantity == Q_EulerianConcentration .or. FieldReq%iQuantity == Q_ChemistryField) Then
            Field%Std(1, iX, iY, iZ, 1, iTA, 1) =                                                         &
              Sum(DispState%EulerianField%Concentration(iXChem, iYChem, iZChem,                           &
                                                        iEulerianTracer:iEulerianTracer + nParticleSizes, &
                                                        DispState%EulerianField%iNew))
                                                        !$$ review when allow particle size segregated output
          Else If (FieldReq%iQuantity == Q_EulerianTotalDep) Then 
            Field%Std(1, iX, iY, 1, 1, iTA, 1) =                                                                 &
              Sum(DispState%EulerianField%TotalDepositionRate(iXChem, iYChem,                                    &
                                                              iEulerianTracer:iEulerianTracer + nParticleSizes))
          Else If (FieldReq%iQuantity == Q_EulerianDryDep) Then
            If (Associated(DispState%EulerianField%DryDepositionRate)) Then 
              Field%Std(1, iX, iY, 1, 1, iTA, 1) =                                                               &
                Sum(DispState%EulerianField%DryDepositionRate(iXChem, iYChem,                                    &
                                                              iEulerianTracer:iEulerianTracer + nParticleSizes))
            End If
          Else If (FieldReq%iQuantity == Q_EulerianWetDep) Then
            If (Associated(DispState%EulerianField%WetDepositionRate)) Then 
              Field%Std(1, iX, iY, 1, 1, iTA, 1) =                                                               &
                Sum(DispState%EulerianField%WetDepositionRate(iXChem, iYChem,                                    &
                                                              iEulerianTracer:iEulerianTracer + nParticleSizes))
            End If
          End If

        End Do

      End Do
      End Do
      End Do

      ! Update the information on which fields have been processed.
      Do i = 1, nFields
        Processed(iFields(i)) = .true.
      End Do

    End Do

  End Do

End Subroutine CalcEulerianResults

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcFlowResults(ReqInfos, Coords, Grids, Domains, Reqs, Units, Mets, Flows, Results)
! Coordinates calculation of type F (depending on flow information) fields.

  Implicit None
  ! Argument list:
  Type(ReqInfo_), Intent(In)            :: ReqInfos(:) ! Information on what requirements are required next.
  Type(Coords_),  Intent(In)            :: Coords      ! Collection of coord systems.
  Type(Grids_),   Intent(In),    Target :: Grids       ! Collection of grids.
  Type(Domains_), Intent(In)            :: Domains     ! Collection of domains.
  Type(Reqs_),    Intent(In),    Target :: Reqs        ! Collection of requirements.
  Type(Units_),   Intent(InOut)         :: Units       ! Collection of information on input/output unit
                                                       ! numbers.
  Type(Mets_),    Intent(InOut)         :: Mets        ! Collection of met module instance states.
  Type(Flows_),   Intent(InOut)         :: Flows       ! Collection of flow module instance states.
  Type(Results_), Intent(InOut), Target :: Results     ! Collection of results.
  ! Locals:
  Integer                    :: nFields                      !} Number and indices of fields to be processed
  Integer                    :: iFields(Reqs%MaxFieldReqs)   !} together.
  Logical                    :: Processed(Reqs%MaxFieldReqs) ! Indicates fields processed.
                                                             ! $$ could use MaxFieldReqsPerReqInfo?
  Real(Std)                  :: Value(9)                     ! Value of field.
  Type(TGrid_),      Pointer :: TGrid                        ! Abbreviation for time grid.
  Type(HGrid_),      Pointer :: HGrid                        !
  Type(ZGrid_),      Pointer :: ZGrid                        !
  Integer                    :: iHCoord                      !} Indices of coord systems.
  Integer                    :: iZCoord                      !}
  Integer                    :: iHGrid                       !] Indices of grids.
  Integer                    :: iZGrid                       !]
  Integer                    :: iTGrid                       !]
  Integer                    :: nX                           !} Size of grids.
  Integer                    :: nY                           !}
  Integer                    :: nZ                           !}
  Integer                    :: nT                           !}
  Integer                    :: iT                           !} Loop indices for looping over T/X/Y/Z.
  Integer                    :: iX                           !}
  Integer                    :: iY                           !}
  Integer                    :: iZ                           !}
  Integer                    :: iTA                          ! Array index at which results for a given time
                                                             ! are stored.
  Integer                    :: i                            !} Loop indices.
  Integer                    :: j                            !}
  Integer                    :: k                            !}
  Type(Flow_)                :: Flow                         !] Flow, cloud, rain and surface information.
  Type(Cloud_)               :: Cloud                        !]
  Type(Rain_)                :: Rain                         !]
  Type(Surface_)             :: Surface                      !]
  Type(Soil_)                :: Soil                         !]
  Type(Plant_)               :: Plant                        !]
  Type(FlowMemory_)          :: FlowMemory                   ! Flow memory.
  Type(Position_)            :: Position                     ! Coords of location in various coord systems
                                                             ! in Coords, with flags to indicate whether
                                                             ! the values are valid.
  Integer                    :: FlowErrorCode                !} Error codes for calls to GetAttrib.
  Integer                    :: CloudErrorCode               !}
  Integer                    :: RainErrorCode                !}
  Integer                    :: SurfaceErrorCode             !}
  Integer                    :: PlantErrorCode               !}
  Type(ShortTime_)           :: TravelTime                   ! Dummy variable for TravelTime argument
                                                             ! in GetAttrib.

  Type(FieldReq_),   Pointer :: FieldReq
  Type(FieldReq_),   Pointer :: FieldReq1
  Type(Field_),      Pointer :: Field
  Integer                    :: iReqInfo

  Do iReqInfo = 1, Size(ReqInfos)

    If (ReqInfos(iReqInfo)%NoReqs) Exit

    Processed(:) = .false.

    Do j = 1, ReqInfos(iReqInfo)%nFields

      If (Processed(j)) Cycle

      FieldReq => Reqs%FieldReqs(ReqInfos(iReqInfo)%iField(j))

      ! Indices and size of grids.
      iHGrid  = FieldReq%iHGrid
      iZGrid  = FieldReq%iZGrid
      iTGrid  = FieldReq%iTGrid
      iHCoord = FieldReq%iHCoord
      iZCoord = FieldReq%iZCoord
      nX      = Grids%HGrids(iHGrid)%nX
      nY      = Grids%HGrids(iHGrid)%nY
      nZ      = Grids%ZGrids(iZGrid)%nZ
      nT      = Grids%TGrids(iTGrid)%nT

      ! Abbreviation for time grid.
      TGrid => Grids%TGrids(iTGrid)
      HGrid => Grids%HGrids(iHGrid)
      ZGrid => Grids%ZGrids(iZGrid)

      ! Determine the fields on same grid with the same values of iHCoord and iZCoord if these are non-zero.
      ! These fields will be processed together.
      nFields    = 1
      iFields(1) = j
      Do i = j + 1, ReqInfos(iReqInfo)%nFields
        FieldReq1 => Reqs%FieldReqs(ReqInfos(iReqInfo)%iField(i))
        If (Processed(i)) Cycle
        If (FieldReq1%iHGrid /= iHGrid .or. FieldReq1%iZGrid /= iZGrid .or. FieldReq1%iTGrid /= iTGrid) Cycle
        If (FieldReq1%iHCoord /= iHCoord .and. FieldReq1%iHCoord /= 0 .and. iHCoord /= 0) Cycle
        If (FieldReq1%iZCoord /= iZCoord .and. FieldReq1%iZCoord /= 0 .and. iZCoord /= 0 ) Cycle
        nFields          = nFields + 1
        iFields(nFields) = i
        iHCoord = Max(FieldReq1%iHCoord, iHCoord)
        iZCoord = Max(FieldReq1%iZCoord, iZCoord)
      End Do

      iT = ReqInfos(iReqInfo)%iT(j)

      ! Loop over space.
      Do iX = 1, nX
      Do iY = 1, nY
      Do iZ = 1, nZ

        ! Get information from the flow module instances.
        Call ResetFlowMemory(Flows, FlowMemory)
        If (HGrid%Unstructured) Then
          Position = X2Position(                                    &
                       Coords,                                      &
                       (/ HGrid%X(iX), HGrid%Y(iX), ZGrid%Z(iZ) /), &
                       HGrid%iHCoord,                               &
                       ZGrid%iZCoord                                &
                     )
        Else
          Position = X2Position(                                    &
                       Coords,                                      &
                       (/ HGrid%X(iX), HGrid%Y(iY), ZGrid%Z(iZ) /), &
                       HGrid%iHCoord,                               &
                       ZGrid%iZCoord                                &
                     )
        End If
        Call GetAttrib(                             &
               iAttribParam  = A_Flow,              &
               Coords        = Coords,              &
               Grids         = Grids,               &
               Domains       = Domains,             &
               Moisture      = .true.,              &
               Inhomog       = .true.,              &
               Homog         = .false.,             &
               Time          = TInTGrid(TGrid, iT), &
               AnyTravelTime = .true.,              &
               TravelTime    = TravelTime,          &
               Position      = Position,            &
               Units         = Units,               &
               Mets          = Mets,                &
               Flows         = Flows,               &
               FlowMemory    = FlowMemory,          &
               Flow          = Flow,                &
               Cloud         = Cloud,               &
               Rain          = Rain,                &
               Surface       = Surface,             &
               Soil          = Soil,                &
               Plant         = Plant,               &
               ErrorCode     = FlowErrorCode        &
             )
        Call GetAttrib(                             &
               iAttribParam  = A_Cloud,             &
               Coords        = Coords,              &
               Grids         = Grids,               &
               Domains       = Domains,             &
               Moisture      = .false.,             &
               Inhomog       = .false.,             &
               Homog         = .false.,             &
               Time          = TInTGrid(TGrid, iT), &
               AnyTravelTime = .true.,              &
               TravelTime    = TravelTime,          &
               Position      = Position,            &
               Units         = Units,               &
               Mets          = Mets,                &
               Flows         = Flows,               &
               FlowMemory    = FlowMemory,          &
               Flow          = Flow,                &
               Cloud         = Cloud,               &
               Rain          = Rain,                &
               Surface       = Surface,             &
               Soil          = Soil,                &
               Plant         = Plant,               &
               ErrorCode     = CloudErrorCode       &
             )
        Call GetAttrib(                             &
               iAttribParam  = A_Rain,              &
               Coords        = Coords,              &
               Grids         = Grids,               &
               Domains       = Domains,             &
               Moisture      = .false.,             &
               Inhomog       = .false.,             &
               Homog         = .false.,             &
               Time          = TInTGrid(TGrid, iT), &
               AnyTravelTime = .true.,              &
               TravelTime    = TravelTime,          &
               Position      = Position,            &
               Units         = Units,               &
               Mets          = Mets,                &
               Flows         = Flows,               &
               FlowMemory    = FlowMemory,          &
               Flow          = Flow,                &
               Cloud         = Cloud,               &
               Rain          = Rain,                &
               Surface       = Surface,             &
               Soil          = Soil,                &
               Plant         = Plant,               &
               ErrorCode     = RainErrorCode        &
             )
        Call GetAttrib(                             &
               iAttribParam  = A_Surface,           &
               Coords        = Coords,              &
               Grids         = Grids,               &
               Domains       = Domains,             &
               Moisture      = .false.,             &
               Inhomog       = .false.,             &
               Homog         = .false.,             &
               Time          = TInTGrid(TGrid, iT), &
               AnyTravelTime = .true.,              &
               TravelTime    = TravelTime,          &
               Position      = Position,            &
               Units         = Units,               &
               Mets          = Mets,                &
               Flows         = Flows,               &
               FlowMemory    = FlowMemory,          &
               Flow          = Flow,                &
               Cloud         = Cloud,               &
               Rain          = Rain,                &
               Surface       = Surface,             &
               Soil          = Soil,                &
               Plant         = Plant,               &
               ErrorCode     = SurfaceErrorCode     &
             )
        Call GetAttrib(                             &
               iAttribParam  = A_Plant,             &
               Coords        = Coords,              &
               Grids         = Grids,               &
               Domains       = Domains,             &
               Moisture      = .false.,             &
               Inhomog       = .false.,             &
               Homog         = .false.,             &
               Time          = TInTGrid(TGrid, iT), &
               AnyTravelTime = .true.,              &
               TravelTime    = TravelTime,          &
               Position      = Position,            &
               Units         = Units,               &
               Mets          = Mets,                &
               Flows         = Flows,               &
               FlowMemory    = FlowMemory,          &
               Flow          = Flow,                &
               Cloud         = Cloud,               &
               Rain          = Rain,                &
               Surface       = Surface,             &
               Soil          = Soil,                &
               Plant         = Plant,               &
               ErrorCode     = PlantErrorCode       &
             )


        ! Convert information from the flow module instances to the desired coord system.
        If (FlowErrorCode == 0 .and. (iHCoord /= 0 .or. iZCoord /= 0)) Then
          If (iHCoord == 0) iHCoord = Flow%iHCoord
          If (iZCoord == 0) iZCoord = Flow%iZCoord
          Call ConvertFlow(Coords, iHCoord, iZCoord, Position, Flow)
          ! $$ other attributes?
        End If

        ! Calculate the field values.
        Do i = 1, nFields

          FieldReq1 => Reqs%FieldReqs( ReqInfos(iReqInfo)%iField(iFields(i)) )
          Field     => Results%Fields( ReqInfos(iReqInfo)%iField(iFields(i)) )

          iTA = CalciTA(Field, iT)

          ! $$ write error message to warn of missing output?
          ! If so check flow results needed - i.e. PPs exist (i.e. between FirstReleaseTime and time of
          ! demise of tracers) or between FirstFlowTimeForReqs & LastFlowTimeForReqs

          ! No flow data available.
          If (FlowErrorCode > 0) Then
            Select Case (FieldReq1%iQuantity)
              Case (                &
                Q_SigmaWW,          &
                Q_TauWW,            &
                Q_MeanFlowU,        &
                Q_MeanFlowV,        &
                Q_MeanFlowW,        &
                Q_TemperatureK,     &
                Q_TemperatureC,     &
                Q_PotentialTempK,   &
                Q_SpecificHumidity, &
                Q_PressurePa,       &
                Q_Density,          &
                Q_Topography,       &
                Q_UStar,            &
                Q_HeatFlux,         &
                Q_BLDepth,          &
                Q_WindSpeed,        &
                Q_WindDirectionDeg, &
                Q_RHPercent,        &
                Q_Pasquill,         &
                Q_SigmaVV,          &
                Q_MesoscaleSigmaVV  &
              )
                Field%Std(1, iX, iY, iZ, 1, iTA, 1) = -Huge(Value(1))
                Cycle
            End Select
          End If

          ! No cloud data available.
          If (CloudErrorCode > 0) Then
            Select Case (FieldReq1%iQuantity)
              Case (                    &
                Q_CloudOktas,           &
                Q_TotalOrDynCloudWater, &
                Q_TotalOrDynCloudIce,   &
                Q_Cloud3d,              &
                Q_ConCloudBase,         &
                Q_ConCloudTop,          &
                Q_Pasquill              &
              )
                Field%Std(1, iX, iY, iZ, 1, iTA, 1) = -Huge(Value(1))
                Cycle
            End Select
          End If

          ! No rain data available.
          If (RainErrorCode > 0) Then
            Select Case (FieldReq1%iQuantity)
              Case (          &
                Q_PptRateMMHR &
              )
                Field%Std(1, iX, iY, iZ, 1, iTA, 1) = -Huge(Value(1))
                Cycle
            End Select
          End If

          ! No surface data available.
          If (SurfaceErrorCode > 0) Then
            Select Case (FieldReq1%iQuantity)
              Case (               &
                Q_LandUseFracs,    &
                Q_SoilMoisture,    &
                Q_LandFrac         &
              )
                Do k = 1, FieldReq1%nComponents
                  Field%Std(1, iX, iY, iZ, 1, iTA, k) = -Huge(Value(1))
                End Do
                Cycle
            End Select
          End If

          ! No plant data available.
          If (PlantErrorCode > 0) Then
            Select Case (FieldReq1%iQuantity)
              Case (               &
                Q_CanopyWater,     &
                Q_LAI,             &
                Q_CanopyHeight,    &
                Q_StomataConduct   &
              )
                Do k = 1, FieldReq1%nComponents
                  Field%Std(1, iX, iY, iZ, 1, iTA, k) = -Huge(Value(1))
                End Do
                Cycle
            End Select
          End If

          ! All required data available.
          Do k = 1, FieldReq1%nComponents

            Select Case (FieldReq1%iQuantity)

              Case (Q_SigmaWW) ! $$ put these in a sensible order, above too

                Value(k) = Flow%SigUU(3)

              Case (Q_TauWW)

                Value(k) = Flow%TauUU(3)

              Case (Q_MeanFlowU)

                Value(k) = Flow%U(1)

              Case (Q_MeanFlowV)

                Value(k) = Flow%U(2)

              Case (Q_MeanFlowW)

                Value(k) = Flow%U(3)

              Case (Q_TemperatureK)

                Value(k) = Flow%T

              Case (Q_TemperatureC)

                Value(k) = Flow%T - TKAtTCEq0

              Case (Q_PotentialTempK)

                Value(k) = Flow%Theta

              Case (Q_SpecificHumidity)

                Value(k) = Flow%Q

              Case (Q_PressurePa)

                Value(k) = Flow%P

              Case (Q_Density)

                Value(k) = Flow%Rho

              Case (Q_Topography)

                Value(k) = Flow%Topog

              Case (Q_UStar)

                Value(k) = Flow%UStar

              Case (Q_HeatFlux)

                Value(k) = Flow%WT * CP * Flow%PS / (Flow%T0 * GasConstant)

              Case (Q_BLDepth)

                Value(k) = Flow%H

              Case (Q_WindSpeed)

                Value(k) = Sqrt(Flow%U(1)**2 + Flow%U(2)**2)

              Case (Q_WindDirectionDeg)

                Value(k) = ATan2ZeroTest(-Flow%U(1), -Flow%U(2)) * 180.0 / Pi
                If (Value(k) < 0.0) Value(k) = Value(k) + 360.0

              Case (Q_PptRateMMHR)

                Value(k) = Rain%ConPpt + Rain%DynPpt

              Case (Q_CloudOktas)

                Value(k) = Cloud%Cloud * 8.0

              Case (Q_RHPercent)

                Value(k) = CalcRH(Flow%Q, Flow%T, Flow%P)

              Case (Q_Pasquill)

                Value(k) = CalcPasquill(                                       &
                             Sqrt(Flow%U(1)**2 + Flow%U(2)**2),                &
                             Flow%WT * CP * Flow%PS / (Flow%T0 * GasConstant), &
                             Cloud%Cloud                                       &
                           )

              Case (Q_SigmaVV)

                Value(k) = Flow%SigUU(2)

              Case (Q_MesoscaleSigmaVV)

                Value(k) = Flow%SigUUM

              Case (Q_TotalOrDynCloudWater)

                Value(k) = Cloud%TotalOrDynCloudWater

              Case (Q_TotalOrDynCloudIce)

                Value(k) = Cloud%TotalOrDynCloudIce

              Case (Q_Cloud3d)

                Value(k) = Cloud%Cloud3d

              Case (Q_RoughnessLength)

                Value(k) = Flow%Z0

              Case (Q_PSeaLevelPa)

                Value(k) = Flow%PSeaLevel

              Case (Q_LandUseFracs)

                Value(k) = Surface%LandUseFracs(k)

              Case (Q_CanopyWater)

                Value(k) = Plant%CanopyWater(k)

              Case (Q_LAI)

                Value(k) = Plant%LAI(k)

              Case (Q_CanopyHeight)

                Value(k) = Plant%CanopyHeight(k)

              Case (Q_StomataConduct)

                Value(k) = Plant%StomataConduct(k)

              Case (Q_SoilMoisture)

                Value(k) = Surface%SoilMoisture

              Case (Q_LandFrac)

                Value(k) = Surface%LandFrac
              
              Case (Q_ConCloudBase)
              
                Value(k) = Cloud%ConCloudBase
                
              Case (Q_ConCloudTop)
              
                Value(k) = Cloud%ConCloudTop 

            End Select

            Field%Std(1, iX, iY, iZ, 1, iTA, k) = Value(k)

          End Do

        End Do

      End Do
      End Do
      End Do

      ! Update the information on which fields processed.
      Do i = 1, nFields
        Processed(iFields(i)) = .true.
      End Do

    End Do

  End Do

End Subroutine CalcFlowResults

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcSourceResults(ReqInfos, Reqs, Results, Q)
! Coordinates calculation of type A fields

  Implicit None
  ! Argument list:
  Type(ReqInfo_), Intent(In)            :: ReqInfos(:)
  Type(Reqs_),    Intent(In),    Target :: Reqs
  Type(Results_), Intent(InOut), Target :: Results
  Real(Std),      Intent(In)            :: Q(:)
  ! ReqInfos :: Requirements information for source.
  ! Results  :: Collection of results.

  ! Locals:
  Type(FieldReq_), Pointer :: FieldReq      ! Local pointer to a field request.
  Type(Field_),    Pointer :: Field         ! Local pointer to a field.
  Integer                  :: iReqInfo      ! Index of set of requirements in ReqInfos.
  Integer                  :: i             ! Loop indices over field requirements.
  Integer                  :: iT            ! Index of the field.
  Integer                  :: iTA           ! Array index at which results for a given time are stored.

  ! Loop over each set of requirements in ReqInfos.  
  Do iReqInfo = 1, Size(ReqInfos)

    If (ReqInfos(iReqInfo)%NoReqs) Exit

    Do i = 1, ReqInfos(iReqInfo)%nFields

      FieldReq => Reqs%FieldReqs( ReqInfos(iReqInfo)%iField(i) )
      Field => Results%Fields( ReqInfos(iReqInfo)%iField(i) )

      ! Time index of the field.
      iT = ReqInfos(iReqInfo)%iT(i)
      iTA = CalciTA(Field, iT)

      If (FieldReq%iQuantity == Q_OriginalSourceStrength) Then
        Field%Std(1, 1, 1, 1, 1, iTA, 1) = Q(1)
      Else
        Field%Std(1, 1, 1, 1, 1, iTA, 1) = Q(2)
      End If

    End Do
  End Do

End Subroutine CalcSourceResults

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcOtherResults(ReqInfos, RunToTime, Coords, Grids, Domains, Reqs, DispState, Results)
! Coordinates calculation of type O (other) fields.

  Implicit None
  ! Argument list:
  Type(ReqInfo_),   Intent(In)           :: ReqInfos(:) ! Information on what requirements are required next.
  Type(Time_),      Intent(In)           :: RunToTime   ! Time at which run will be suspended.
  Type(Coords_),    Intent(In),   Target :: Coords      ! Collection of coord systems.
  Type(Grids_),     Intent(In),   Target :: Grids       ! Collection of grids.
  Type(Domains_),   Intent(In)           :: Domains     ! Collection of domains.
  Type(Reqs_),      Intent(In),   Target :: Reqs        ! Collection of requirements.
  Type(DispState_), Intent(In),   Target :: DispState   ! Information on the state of the dispersion
                                                        ! calculation.
  Type(Results_),   Intent(InOut)        :: Results     ! Collection of results.
  ! Locals:
  Integer                    :: nFields                      !} Number and indices of fields to be processed
  Integer                    :: iFields(Reqs%MaxFieldReqs)   !} together.
  Logical                    :: Processed(Reqs%MaxFieldReqs) ! Indicates fields processed.
  Real(Std)                  :: Value                        ! Value of field.
  Type(TGrid_),      Pointer :: TGrid                        ! Abbreviation for time grid.
  Type(FieldReq_),   Pointer :: FieldReq                     ! Abbreviation for field requirement.
  Integer                    :: iHCoord                      !} Indices of coord systems.
  Integer                    :: iZCoord                      !}
  Integer                    :: iHGrid                       !] Indices of grids.
  Integer                    :: iZGrid                       !]
  Integer                    :: iTGrid                       !]
  Integer                    :: nX                           !} Size of grids.
  Integer                    :: nY                           !}
  Integer                    :: nZ                           !}
  Integer                    :: nT                           !}
  Integer                    :: iTU                          !] Upper and lower limits for time loop.
  Integer                    :: iTL                          !]
  Integer                    :: iT                           !} Loop indices for looping over T/X/Y/Z.
  Integer                    :: iX                           !}
  Integer                    :: iY                           !}
  Integer                    :: iZ                           !}
  Integer                    :: iTA                          ! Array index at which results for a given time
                                                             ! are stored.
  Integer                    :: i                            !} Loop indices.
  Integer                    :: j                            !}
  Integer                    :: iReqInfo                     !}
  Integer                    :: iField                       !}
  Type(Position_)            :: Position                     ! Coords of location in various coord systems
                                                             ! in Coords, with flags to indicate whether
                                                             ! the values are valid.
  Type(ShortTime_)         :: StartTime
  Type(ShortTime_)         :: EndTime
  Type(ShortTime_)         :: Time
  Type(HCoord_),   Pointer :: HCoordIn
  Type(HCoord_),   Pointer :: HCoordOut
  Type(HGrid_),    Pointer :: HGrid
  Real(Pos)                :: Point(2)

  Do iReqInfo = 1, Size(ReqInfos)

    If (ReqInfos(iReqInfo)%NoReqs) Exit

    Do i = 1, ReqInfos(iReqInfo)%nFields

      iField   =  ReqInfos(iReqInfo)%iField(i)
      FieldReq => Reqs%FieldReqs(iField)

      ! Process t-based requirements.
      If (FieldReq%iTGrid /= 0) Then

        iT  = ReqInfos(iReqInfo)%iT(i)
        iTA = CalciTA(Results%Fields(iField), iT)

        If (FieldReq%iQuantity == Q_ProgressPercent) Then

          StartTime = TMin(DispState%sFirstReqTime, DispState%sFirstReleaseTime)
          EndTime   = TMax(TMax(DispState%sLastDispTime, DispState%sLastReqTime), DispState%sLastReleaseTime)
          EndTIme   = TMin(EndTime, Time2ShortTime(RunToTime))
          Time      = Time2ShortTime(DispState%Time)

          Results%Fields(iField)%Std(1, 1, 1, 1, 1, iTA, 1) =                                      &
            100.0 * ShortTime2RealTime(Time - StartTime) / ShortTime2RealTime(EndTime - StartTime)

        Else If (FieldReq%iQuantity == Q_ClockTime) Then

          Results%Fields(iField)%T(1, 1, 1, 1, 1, iTA, 1) = CurrentClockTime()

        End If

      ! Process t-independent requirements.
      Else

        If (FieldReq%iQuantity == Q_X) Then

          HCoordIn  => Coords%HCoords(FieldReq%iHGridCoord)
          HCoordOut => Coords%HCoords(FieldReq%iHCoord)
          HGrid     => Grids%HGrids(FieldReq%iHGrid)
          Do iX = 1, HGrid%nX
          Do iY = 1, HGrid%nY
            If (HGrid%Unstructured) Then
              Point = ConvertH(HCoordIn, HCoordOut, (/ HGrid%X(iX), HGrid%Y(iX) /))
            Else
              Point = ConvertH(HCoordIn, HCoordOut, (/ HGrid%X(iX), HGrid%Y(iY) /))
            End If
            Results%Fields(iField)%Std(1, iX, iY, 1, 1, 1, 1) = Point(1)
          End Do
          End Do

        Else If (FieldReq%iQuantity == Q_Y) Then

          HCoordIn  => Coords%HCoords(FieldReq%iHGridCoord)
          HCoordOut => Coords%HCoords(FieldReq%iHCoord)
          HGrid     => Grids%HGrids(FieldReq%iHGrid)
          Do iX = 1, HGrid%nX
          Do iY = 1, HGrid%nY
            If (HGrid%Unstructured) Then
              Point = ConvertH(HCoordIn, HCoordOut, (/ HGrid%X(iX), HGrid%Y(iX) /))
            Else
              Point = ConvertH(HCoordIn, HCoordOut, (/ HGrid%X(iX), HGrid%Y(iY) /))
            End If
            Results%Fields(iField)%Std(1, iX, iY, 1, 1, 1, 1) = Point(2)
          End Do
          End Do

        End If

      End If

    End Do

  End Do

End Subroutine CalcOtherResults

!-------------------------------------------------------------------------------------------------------------

End Module CaseModule

! Module:  Output Module

Module OutputModule

! This module provides code for storing and outputting results.

! Module overview
! ---------------

! The items to be calculated are specified in a number of 'requirements'. Requirements can be for a 'field' or
! a 'set of particle/puff information'.
!
! (a) Fields
!
! Fields are grided quantities which may depend on T, S (travel time), H (horizontal position), Z, and, for
! exceedence probabilities or percentiles, on a value or percentage (denoted D for 'data value'). In addition
! they may sometimes have multiple components.
!
! They may be items required by the user or internally generated requirements which are needed as an
! intermediate step in calculating the user's requirements.
!
! Raw fields can be of one of four types:
!   type L - quantities derived from (Lagrangian) particle/puff data,
!   type E - quantities derived from (Eulerian) grided concentration data,
!   type F - quantities derived from flow data,
!   type O - quantities not depending on any of the above.
!   type A - quantities associated with source properties
! Raw fields can be combined to produce derived fields which can be of one of three types:
!   type P - quantities which, when concentrations = 0, are zero,
!   type Q - quantities which, when concentrations = 0, don't depend on type F quantities,
!   type R - other quantities.
! Raw and derived fields can be 'processed' by a sequence of processing steps drawn from the following list:
!   (1) calculation of averages or integrals,
!   (2) calculation of exceedence probabilities,
!   (3) calculation of percentiles (can include 100 and 0 percentiles),
!   (4) calculation of max or min (more efficient than (3) if only 100 and/or 0 percentiles are needed),
!   (5) calculation of moments about zero,
!   (6) calculation of central moments.
!
! A dependence on S is only possible for type L quantities and for type P, Q and R quantities which depend on
! type L quantities which depend on S.
!
! Type F quantities and type P, Q and R quantities which depend on type F quantities must depend on T, H and Z
! (at least prior to any processing steps) in order to determine which flow module instance to use.
!
! Type L, E, F, P, Q, and R quantities must depend on time (at least prior to any processing steps).
!
! Note the difference between the dependence of a field (after any processing steps) on T/S/H/Z and the
! dependence of the underlying quantity (prior to any processing steps) on T/S/H/Z. The field may have fewer
! dependencies following processing.
!
! Except for averaging/integrating processing steps with only a single input grid value contributing to each
! output grid value, each derived field calculation or processing step calculates a new field from the
! previously available ones. ! $$ not sure the exception is quite correct.
!
! Contributions to a field are always restricted to the computational domain. Zero contributions are assumed
! outside the domain. Warning/error messages are given for type F, O, Q or R quantities if contributions from
! outside the computational domain in space/time are requested. Contributions to type F quantities, or via
! type F quantities to type Q or R quantities, from locations which are below ground or inside buildings etc
! but which are inside the computational domain will depend on what information the flow modules return.
! However flow modules should be designed to give plausible, if not accurate values in such cases.
!
! $$ do stuff in this para. Note buildings needs to give zero wind and sensible temp etc inside buildings
!
! Processing steps:
!
! Any processing step (3) must be preceeded by a processing step (2), and it doesn't calculate the percentiles
! of the probabilities but calculates the percentiles of the thing the probabilities were calculated from (by
! inverting the probability distribution). Processing steps (2) to (6) use a 'data grid' to specify the
! probabilities, percentiles, max/min or moments required. Averages/integrals, probabilities, percentiles,
! max/min and moments can be over space, time and the ensemble of cases. A percentile calculation is always
! referred to as being over one point in space and time and as not being over the ensemble, as the preceeding
! probability processing step will give the appropriate region. Of course 'the ensemble of cases' cannot be
! invoked in more than one step.
!
! The processing steps can add an extra dimension to the resulting field, namely if the data grid specified
! for calculating probabilities, percentiles, max/min or moments has more than one point. The processing steps
! must not produce results which have more than one such extra dimension. This means that the probability
! processing steps without a following percentile step, the percentile processing steps, the max/min
! processing steps, and the moment processing steps which are used in producing a given output must not have
! more than one step with a data grid with more than one point. Also, after such a processing step, not even
! probability processing steps with a following percentile step can have a data grid with more than one point.
! In addition, probability processing steps with a following percentile processing step must have a data grid
! with more than one point. There would be little value in this situation in allowing data grids with only one
! point and, by preventing this, we can be sure there is no data grid with more than one point which exists
! before the probability processing step and needs to be 'carried through' the percentile processing step. In
! imposing these restrictions, 'a data grid with more than one point' is interpreted as including floating
! data grids (a floating data grid might in principle have only one point; however it is hard to relax the
! restrictions for this case). Finally, we note that floating data grids can only be used for probability
! processing steps with a following percentile processing step or for probability processing steps with no
! following processing step. This ensures floating data grids do not need to be 'carried though' processing
! steps.
!
! $$ could also have a second type of distribution - consisting of the collection of values, or the
! highest/lowest so many values.
!
! If the processing is complex it is difficult to label the output to reflect all the processing steps. The
! output column headers describe the processing only if it is 'simple'. Otherwise the processing is described
! simply as 'complex processing'. The processing is simple if it consists only of
!   (1) a single averaging/integrating processing step, followed by
!   (2) a probability processing step without a following percentile processing step, a probability processing
!       step with a following percentile processing step, a max/min processing step, or a moment processing
!       step.
! The column headers may not identify the output uniquely (especially for complex processing). However the
! field requirement name, which also appears in the column header, can be used to identify the output.
!
! For simple or complex processing, the data grid values and units given in the output are determined as
! follows. If there is only a single processing step with a data grid or there is a processing step with a
! data grid with more than one point (discounting probability steps with following percentile steps), then
! this grid is given. Otherwise the output is presented as if there was no data grid (and labelled 'complex
! processing'). The units given correspond to the result of all the processing steps, although, where
! probabilities or moments are plotted, the unit is that of the random variable whose probability or moment is
! plotted. ! $$ Here > 1 point includes any floating grid (as above)?
!
! Averaging/integrating:
!
! This is always labelled with the centre of the averaging/integrating region in space and the end of the
! averaging/integrating interval in time (relative to the direction of the run). Travel time
! averaging/integrating is not possible.
!
! Where the first processing step involves averaging/integrating, any type L quantities used can (as an
! option) be 'intrinsically averaged/integrated' before storing on a grid using the fact that the raw
! information is not intrinsically grid based. This will involve contributions from the actual particle/puff
! positions and every particle/puff evolution time step. Of course for type L quantities calculated from
! particles, there is always some spatial intrinsic averaging whether this is requested or not (equivalent to
! specifying averaging with 'H/Z window defined by grid').
!
! Note an averaging/integrating processing step with no ensemble averaging and only a single input grid value
! contributing to each output grid value may be used
!   (1) to give the region to use for intrinsic averaging/integrating (first processing step only),
!   (2) to give the averaging/integrating region to take account of in fluctuation predictions (first
!       processing step only),
!   (3) to multiple the value by the size of the integrating region (but not with intrinsic integrating which
!       is covered in (1) - if there is intrinsic integrating the quantity is computed over the integrating
!       region and not multiplied afterwards in the processing), or
!   (4) to do nothing at all.
! Except for use (3), this leads to no processing step (in the sense of calculating a new field from the
! previously available ones) actually being required.
!
! The way in which averaging/integrating processing steps occuring before any other types of processing are
! organised requires some comments (averaging/integrating processing steps occuring after any other types of
! processing are organised straightforwardly).
!   (1) For type E, F and O fields this is straightforward.
!   (2) For type L fields, the first processing step may involve intrinsic and grid averaging/integrating. If
!       this is the case, an extra field requirement is set up for the intrinsically averaged/integrated
!       field, the result of which is used in the grid averaging/integrating. (This is no different from if
!       there is only grid averaging/integrating - an extra field requirement is set up in this case too.)
!   (3) For type P, Q and R fields, the calculation of the derived field may occur at any point in these
!       averaging/integrating processing steps. This is determined by the underlying quantity. For example a
!       quantity A may be derived from B and A may be such that the first processing step (if it exists and is
!       an averaging/integrating step) is applied before calculating A from B. Suppose a field requirement for
!       A has two processing steps. Then a field requirement for B will be created with the first processing
!       step (or A's reinterpretation of this - see below) specified. The results will be used to calculate A
!       and then the second processing step will be applied. Equally the quantity A may be such that all
!       processing steps are applied to A after calculating A from B. Note however that any intrinsic
!       averaging/integrating is always passed through to any required type L fields. The above provides the
!       only good reason why one might want to specify a processing step which does nothing and one of the
!       only two good reasons why one might want to specify two or more consecutive averaging/integrating
!       steps.
!
! $$ document choices made in (3) for different quantities and derived fields needed.
!
! In creating field requirements for contributions to a derived field it is possible to reinterpret the
! processing steps. However this mustn't require contributions from outside the time extent of the original
! processing step (and quantities which don't depend on time are outside any time extent).
!
! $$ check this last sentence - cf $-comments at start of block defining 'extra requirements'. Also comment
! in Type(Process_):
!            (c) Averaging/integrating over infinite regions of space or time is not allowed for quantities
!                of type F, O, Q or R. It can however be applied to quantities of other types used to
!                derived the type F/O/Q/R quantity.
!                $$ for example integrated deposition rate treated as an instantaneous quantity (total dep)
!                for further processing.
!
! Note there is an option to invoke the decay of deposited material contributing to a time average or integral
! (so that an integral is the total material on the ground at the end of the integration period, not the
! integral of the deposition rate, and an average is the equivalent integral divided by the averaging period).
! This is implemented only for the first processing step. This provides the second good reason why one might
! want to specify two or more consecutive averaging/integrating steps.
!
! The precise interpretation of averaging/integrating for averaging/integrating processing steps occuring
! before any other types of processing, is as follows (averaging/integrating processing steps occuring after
! any other types of processing are interpreted straightforwardly):
!   (1) for quantities which can be defined at a point, e.g. 'air concentration' or 'deposition rate',
!       averaging/integrating over space and/or time is interpreted straightforwardly,
!   (2) for extensive quantities, e.g. 'mass' or '# particles', it is interpreted as summing over space and
!       averaging/integrating over time,
!   (3) for mass-weighted statistics of quantities evaluated at contaminant positions, e.g. 'Mean Travel Time'
!       or 'X Stats', it is interpreted as a mass-weighted average/integral.
!
! replace third option with derived field calc (remember to remove M as a valid character when
! checking QInfo) $$
! or add mass component to mean travel time (like X Stats).
!
! For quantities which depend on both travel-time and time, a particle/puff contributes to the requirement
! only at the given travel-time. Hence results must have some intrinsic averaging/integrating over time, since
! there may not be any particle/puffs with exactly the specified values for both travel-time and time.
!
! $$ not quite true for puffs
!
! Fluctuations:
!
! In some cases a probability or moment calculation step can take account of unresolved fluctuations, provided
! there is at most a single previous processing step and this step is an averaging/integrating step not
! involving an ensemble average.
!
! Fluctuation parameters (quantities determining the effect of fluctuations on the probabilities or moments)
! are calculated as 'derived fields'. Of course if such quantities are output then this is a second type of
! output which takes account of fluctuations. At most a single processing step is allowed prior to calculating
! fluctuation parameters and this step must be an averaging/integrating step not involving an ensemble
! average.
!
! The size of the averaging/integrating region in each of the T, X, Y and Z directions to be used in
! estimating fluctuations is taken from the single averaging/integrating step before calculating the
! fluctuation parameters or calculating the probabilities/moments or is zero if there is no such step, unless
! separate fluctuation averaging/integrating scales are given.
!
! Results obtained using separate fluctuation averaging/integrating scales require careful interpretation.
! This option can be used to allow consideration of small scale stochastic fluctuations while smoothing the
! resolved scales, e.g. in order to reduce sampling noise due to the finite number of particles/puffs. For a
! good fluctuations model results should be relatively insensitive to this smoothing (other than in showing an
! improvement due to noise reduction and replacing a deterministic treatment of some scales with a stochastic
! treatment). It can also be used to avoid the need to smooth the resolved scales if these are in any case
! slowly varying relative to the fluctuation averaging/integrating scales. It doesn't really make sense to
! make the fluctuation averaging/integrating scales larger than the scale of variation of the resolved scales,
! but since the latter isn't known we don't test for this. It also doesn't make sense to make the fluctuation
! averaging/integrating scales differ from the resolved averaging/integrating scales when integrating rather
! than averaging is being considered, and this is prohibited. It is also not possible to set the fluctuation
! averaging/integrating scales to infinity using separate fluctuation averaging/integrating scales.
!
! Note fluctuation parameters must allow one processing step before calculating the derived field that is the
! fluctuation parameter, in order to specify the size of the averaging/integrating region and to use this
! information to calculate any field requirements needed to derive the fluctuation parameter.
!
! (b) Sets of particle/puff information
!
! The second type of output which can be produced is detailed information about the particles/puffs. Unlike
! field information, this information is not stored in arrays prior to output but is output as it is computed,
! and output files are not produced unless there is at least a single piece of particle/puff information to
! report.
!
! (c) Computation/output of requirements and the run duration
!
! The question of what output is produced and the run duration is guided by the following principles:
!   (a) Output is limited by the computational domain in time.
!   (b) The run (for each case) starts early enough to produce all output within the computational domain in
!       time and to release all tracers due for release within the computational domain in time.
!   (c) The run (for each case) keeps going for long enough to produce all output within the computational 
!       domain in time.
!
! Not all of this of course is controlled from within this module.
!
! We will now describe these matters and their implementation in more detail, including the issue of when the
! output is computed and produced.
!
! At the first time when there is no (and there is no possibility of there being in the future) tracer which
! might influence the values of the output which is due to be produced according to (a)-(d) above, any tracer
! which is present is killed. This time is called the 'time of demise of all tracers'.
! $$ Tracers here means Lagrangian tracers. Perhaps revisit with the advent of Eulerian tracers.
!
! Fields depending on time are not output if they are outside the computational domain in time. Fields which 
! don't depend on time due to the underlying quantity not depending on time or due to the field being 
! processed over all time are always output.
!
! Computation of fields with nominal time ReqTime is completed at the first opportunity after time
!            { ReqTime + MaxPPTimeLag  if before demise of all tracers
!            { ReqTime                 if at or after demise of all tracers
! is reached (where MaxPPTimeLag is the maximum possible interval which can occur at any time between the
! trailing time-edge and the nominal time of the particles/puffs), except for fields which don't depend on
! time. Computation of fields where the underlying quantity doesn't depend on time is completed at the first
! opportunity after the start of the case, and computation of fields which are processed over all time is
! completed at the end of the case. Computation of fields processed over all cases is completed at the end of
! all cases.
!
! Fields are output at the first opportunity after their computation is completed with the exception of cases
! where different times are output to the same file with times varying across (rather than down) the page -
! here it is not possible to output the fields until all results for the file are ready for output.
!
! Particle/puff information is output as it is produced and so is only produced when particles/puffs are
! present.
!
! Note the run may need to carry on past the end of the computational domain if MaxPPTimeLag > 0.

! Module use
! ----------

! $$

! Module call tree
! ----------------

! PreInitOutputOpts
!
! InitOutputOpts
!
! CheckOutputOpts
!
! InitReqs
!
! InitProcess
!
! AddProcess---ProcessEq(==)
!
! InitFieldReq
!
! AddFieldReq---FieldReqEq(==)
!
! InitPdfReq $$
!
! AddPdfReq $$
!
! InitPPInfoReq
!
! AddPPInfoReq---PPInfoReqEq(==)
!
!                (InitProcess
!                (
!                (AddProcess---ProcessEq(==)
!                (
!                (InitFieldReq
! $$ SetUpReqs---(
!                (AddFieldReq---FieldReqEq(==)
!                (
!                (FindFieldReqIndex---FieldReqEq(==)
!                (
!                (FindEquivFieldReqIndex---FieldReqEquiv(.Equiv.)
!
!              (SpeciesModule.FindSpeciesIndex
!              (
! SetUpiReqs---(SourceModule.FindSourceIndex
!              (
!              (FindProcessIndex
!
! FirstLastReqTimes
!
! InitResults
!
! WriteRestartFileResults
!
! ReadRestartFileResults
!
! RestartAdjustmentForResults
!
! PrepareResults
!
! UpdateResults
!
! CalciTA
!
! CalcReqInfoNoT
!
! CalcReqInfoT---CalcNextReqInfoT
!
! CalcReqInfoS---CalcNextReqInfoS
!
! CalcReqInfoEveryT
!
!                              (ProcessFields---(Species.CalcDecayFactor
!                              (                (
! $$ ProcessAndOutputResults---(DerivePdfs---
!                              (
!                              (OutputFields---(FieldNumHeaders
!                              (               (FieldGraphHeaders
!                              (OutputPdfs---
!                              (
!                              (ClosePPInfos
!
! OutputPPInfos---(PPInfoNumHeaders
!                 (PPInfoColumnHeaders
!
! Calls to routines in the service modules are not included.
! Calls to routines in other modules are included (indicated by ModuleName.RoutineName), but subsequent calls
! from the called routine are not.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowsModule, Only: Flows_, FlowMemory_,                              &
                       FlowList, ResetFlowMemory, ConvertToZ, GetAttrib, &
                       A_Flow,                                           &
                       Mets_,                                            &
                       Flow_, Cloud_, Rain_, Surface_, Soil_, Plant_
Use SizeDistModule
Use SpeciesModule
Use SourceModule
Use FluctuationModule
Use OpenMPModule, Only: OpenMPOpts_
Use TimerModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: Quantities                  ! Names of quantities useable in field requirements.
                                       ! $$ private if tokens2fieldreqs
Public  :: QInfo                       ! Table of information on the quantities.
                                       ! $$ moved to this module
Public  :: Q_AirConc                   !} Codes for the quantities in fields requirements.
Public  :: Q_DryDep                    !}
Public  :: Q_WetDep                    !}
Public  :: Q_Dep                       !}
Public  :: Q_nParticles                !}
Public  :: Q_nParticlesBySpec          !}
Public  :: Q_nPuffs                    !}
Public  :: Q_nParticleSteps            !}
Public  :: Q_nPuffSteps                !}
Public  :: Q_Mass                      !}
Public  :: Q_MeanZ                     !}
Public  :: Q_SigmaZ                    !}
Public  :: Q_XStats                    !}
Public  :: Q_MeanS                     !}
Public  :: Q_PuffCentres               !}
Public  :: Q_SigmaC                    !}
Public  :: Q_ChemistryField            !}
Public  :: Q_SigmaWW                   !}
Public  :: Q_TauWW                     !}
Public  :: Q_MeanFlowU                 !}
Public  :: Q_MeanFlowV                 !}
Public  :: Q_MeanFlowW                 !}
Public  :: Q_TemperatureK              !}
Public  :: Q_PotentialTempK            !}
Public  :: Q_SpecificHumidity          !}
Public  :: Q_PressurePa                !}
Public  :: Q_Density                   !}
Public  :: Q_Topography                !}
Public  :: Q_UStar                     !}
Public  :: Q_HeatFlux                  !}
Public  :: Q_BLDepth                   !}
Public  :: Q_WindSpeed                 !}
Public  :: Q_WindDirectionDeg          !}
Public  :: Q_PptRateMMHR               !}
Public  :: Q_TemperatureC              !}
Public  :: Q_CloudOktas                !}
Public  :: Q_RHPercent                 !}
Public  :: Q_Pasquill                  !}
Public  :: Q_ProgressPercent           !}
Public  :: Q_ClockTime                 !}
Public  :: Q_X                         !}
Public  :: Q_Y                         !}
Public  :: Q_SigmaVV                   !}
Public  :: Q_MesoscaleSigmaVV          !}
Public  :: Q_TotalOrDynCloudWater      !}
Public  :: Q_TotalOrDynCloudIce        !}
Public  :: Q_Cloud3d                   !}
Public  :: Q_RoughnessLength           !}
Public  :: Q_PSeaLevelPa               !}
Public  :: Q_PhotonFlux                !}
Public  :: Q_AduEffCloudGammaDose      !}
Public  :: Q_AduLunCloudGammaDose      !}
Public  :: Q_AduThyCloudGammaDose      !}
Public  :: Q_AduBoSCloudGammaDose      !}
Public  :: Q_AreaAtRisk                !}
Public  :: Q_LandUseFracs              !}
Public  :: Q_CanopyWater               !}
Public  :: Q_LAI                       !}
Public  :: Q_CanopyHeight              !}
Public  :: Q_StomataConduct            !}
Public  :: Q_SoilMoisture              !}
Public  :: Q_LandFrac                  !}
Public  :: Q_MixingRatio               !}
Public  :: Q_EulerianConcentration     !}
Public  :: Q_EMixingRatio              !}
Public  :: Q_ConCloudBase              !}
Public  :: Q_ConCloudTop               !}
Public  :: Q_EulerianTotalDep          !}
Public  :: Q_EulerianDryDep            !}
Public  :: Q_EulerianWetDep            !}
Public  :: Q_OriginalSourceStrength    !}
Public  :: Q_RevisedSourceStrength     !}
Public  :: P_AvInt                     ! Code for calculating averages or integrals.
Public  :: P_Prob                      ! Code for calculating exceedence probabilities.
Public  :: P_Percent                   ! Code for calculating percentiles.          $$ P/PP-codes private if
Public  :: P_MaxMin                    ! Code for calculating max or min.           $$ tokens2fieldreqs moved
Public  :: P_Moments                   ! Code for calculating moments about zero.   $$ to this module
Public  :: P_CentralMoments            ! Code for calculating central moments.
Public  :: P_Av                        ! Code for averaging.
Public  :: P_Int                       ! Code for integrating.
Public  :: OutputOpts_                 ! Output options.
Public  :: Process_                    ! A processing step.
Public  :: FieldReq_                   ! A field requirement.
Public  :: PdfReq_                     ! A pdf requirement.
Public  :: PPInfoReq_                  ! A requirement for a set of particle/puff information.
Public  :: Reqs_                       ! A collection of requirements.
Public  :: Field_                      ! A field.
Public  :: Pdf_                        ! A pdf.
Public  :: PPInfo_                     ! A set of particle/puff information.
Public  :: Results_                    ! A collection of results.
Public  :: ReqInfo_                    ! Information on what requirements are required.
Public  :: PreInitOutputOpts           ! Pre-initialises the output options.
Public  :: InitOutputOpts              ! Initialises the output options.
Public  :: CheckOutputOpts             ! Checks the output options.
Public  :: InitReqs                    ! Initialises a collection of requirements.
Public  :: InitProcess                 ! Initialises a processing step.
Public  :: AddProcess                  ! Adds a processing step to a collection of requirements.
Public  :: InitFieldReq                ! Initialises a field requirement.
Public  :: AddFieldReq                 ! Adds a field requirement to the collection of requirements.
Public  :: InitPdfReq                  ! Initialises a pdf requirement.
Public  :: AddPdfReq                   ! Adds a pdf requirement to the collection of requirements.
Public  :: InitPPInfoReq               ! Initialises a requirement for a set of particle/puff information.
Public  :: AddPPInfoReq                ! Adds a requirement for a set of particle/puff information to the
                                       ! collection of requirements.
Public  :: SetUpReqs                   ! Sets up collection of requirements, including adding any coord
                                       ! systems or grids required.
Public  :: SetUpiReqs                  ! Sets up indices in Reqs.
Public  :: FirstLastReqTimes           ! Calculates the first/last times which the run, flow modules, and
                                       ! dispersion calculation must cover as a result of the requirements.
Public  :: InitResults                 ! Initialises a collection of results.
Public  :: WriteRestartFileResults     ! Writes the collection of results to the restart file.
Public  :: ReadRestartFileResults      ! Reads the collection of results from the restart file.
Public  :: RestartAdjustmentForResults ! Adjusts results read from the restart file.
Public  :: PrepareResults              ! Prepares the collection of results for a new case.
Public  :: UpdateResults               ! Updates the collection of results by deciding which times will be
                                       ! stored from now on.
Public  :: CalciTA                     ! Calculates the time index in a field array at which the results for a
                                       ! given time are stored.
Public  :: CalciLast                   ! Calculates last index in a field array corresponding to a given component and size range.
Public  :: CalcReqInfoNoT              ! Calculates requirements for results which don't depend on time due to
                                       ! the underlying quantity not depending on time.
Public  :: CalcReqInfoT                ! Calculates requirements for results in a given time range which
                                       ! require contributions at discrete times.
Public  :: CalcReqInfoS                ! Calculates requirements for results in a given time and travel-time
                                       ! range which require contributions at discrete travel-times.
Public  :: CalcReqInfoEveryT           ! Calculates requirements for results in a given time range which
                                       ! require contributions at every time in some time interval.
Public  :: ProcessAndOutputResults     ! Processes and outputs results.
Public  :: OutputPPInfos               ! Outputs particle/puff information.
Public  :: OutputModuleTimerInitialise
Public  :: OutputModuleTimerSummary

!-------------------------------------------------------------------------------------------------------------

! Information on quantities in field requirements.
Integer,                  Parameter :: nQuantities = 73
Character(MaxCharLength), Parameter :: Quantities(nQuantities) =                &
                                       (/                                       &
                                         'Air Concentration                  ', & !  1
                                         'Dry Deposition Rate                ', & !  2
                                         'Wet Deposition Rate                ', & !  3
                                         'Deposition Rate                    ', & !  4
                                         '# Particles                        ', & !  5
                                         '# Puffs                            ', & !  6
                                         '# Particle Steps                   ', & !  7
                                         '# Puff Steps                       ', & !  8
                                         'Mass                               ', & !  9
                                         'Mean Z                             ', & ! 10
                                         'Sigma Z                            ', & ! 11
                                         'X Stats                            ', & ! 12
                                         'Mean Travel Time                   ', & ! 13
                                         'Puff Centres                       ', & ! 14
                                         'Sigma C                            ', & ! 15
                                         'Chemistry Field                    ', & ! 16
                                         'Sigma WW                           ', & ! 17
                                         'Tau WW                             ', & ! 18
                                         'Mean Flow U                        ', & ! 19
                                         'Mean Flow V                        ', & ! 20
                                         'Mean Flow W                        ', & ! 21
                                         'Temperature (K)                    ', & ! 22
                                         'Potential Temperature (K)          ', & ! 23
                                         'Specific Humidity                  ', & ! 24
                                         'Pressure (Pa)                      ', & ! 25
                                         'Density                            ', & ! 26
                                         'Topography                         ', & ! 27
                                         'u-star                             ', & ! 28
                                         'Sensible Heat Flux                 ', & ! 29
                                         'Boundary Layer Depth               ', & ! 30
                                         'Wind Speed                         ', & ! 31
                                         'Wind Direction (degrees)           ', & ! 32
                                         'Precipitation Rate (mm/hr)         ', & ! 33
                                         'Temperature (C)                    ', & ! 34
                                         'Cloud Amount (oktas)               ', & ! 35
                                         'Relative Humidity (%)              ', & ! 36
                                         'Pasquill Stability                 ', & ! 37
                                         '# Particles By Species             ', & ! 38
                                         'Progress (%)                       ', & ! 39
                                         'Clock Time                         ', & ! 40
                                         'X                                  ', & ! 41
                                         'Y                                  ', & ! 42
                                         'Sigma VV                           ', & ! 43
                                         'Mesoscale Sigma VV                 ', & ! 44
                                         'Cloud Water (kg/kg)                ', & ! 45 $$ Could be total
                                         'Cloud Ice (kg/kg)                  ', & ! 46 $$ or dynamic
                                         '3d Cloud (Fraction)                ', & ! 47
                                         'Roughness Length                   ', & ! 48
                                         'Sea Level Pressure (Pa)            ', & ! 49 $$ Currently only
                                                                                  !    $$ works for NWP met
                                         'Photon Flux                        ', & ! 50
                                         'Adult Effective Cloud Gamma Dose   ', & ! 51
                                         'Adult Lung Cloud Gamma Dose        ', & ! 52
                                         'Adult Thyroid Cloud Gamma Dose     ', & ! 53
                                         'Adult Bone Surface Cloud Gamma Dose', & ! 54
                                         'Area at risk                       ', & ! 55
                                         'Land Use Fractions                 ', & ! 56
                                         'Canopy Water                       ', & ! 57
                                         'Leaf Area Index                    ', & ! 58
                                         'Canopy Height                      ', & ! 59
                                         'Stomatal Conductance               ', & ! 60
                                         'Soil Moisture                      ', & ! 61
                                         'Land Fraction                      ', & ! 62
                                         'Mixing Ratio                       ', & ! 63
                                         'Eulerian Concentration             ', & ! 64
                                         'E Mixing Ratio                     ', & ! 65
                                         'Concentration                      ', & ! 66
                                         'Convective Cloud Base              ', & ! 67
                                         'Convective Cloud Top               ', & ! 68
                                         'Eulerian Total Deposition Rate     ', & ! 69
                                         'Eulerian Dry Deposition Rate       ', & ! 70
                                         'Eulerian Wet Deposition Rate       ', & ! 71
                                         'Original Source Strength           ', & ! 72 
                                         'Revised Source Strength            '  & ! 73
                                       /)
Character(MaxCharLength), Parameter :: QInfo(nQuantities) =        &
                                       (/                          &
                                         ' 1 L       f T HZ  srd', & !  1 Air Concentration
                                         ' 1 L         T H   srd', & !  2 Dry Deposition Rate
                                         ' 1 L         T H   srd', & !  3 Wet Deposition Rate
                                         ' 1 P 0       T H   srd', & !  4 Deposition Rate
                                         ' 1 L   8 U   T HZ   rd', & !  5 # Particles
                                         ' 1 L     U   T      r ', & !  6 # Puffs
                                         ' 1 L   8     T      r ', & !  7 # Particle Steps
                                         ' 1 L   8     T      r ', & !  8 # Puff Steps
                                         ' 1 L         T     srd', & !  9 Mass
                                         ' 1 L     M   T    zsrd', & ! 10 Mean Z
                                         ' 1 L     M   T    zsrd', & ! 11 Sigma Z
                                         '10 L     M   TS  hzsr ', & ! 12 X Stats
                                         ' 1 L     M   T HZ  srd', & ! 13 Mean Travel Time
                                         ' 1 L         T HZ  sr ', & ! 14 Puff Centres
                                         ' 1 P 1     p T HZ  sr ', & ! 15 Sigma C
                                         ' 1 E         T HZ  s  ', & ! 16 Chemistry Field
                                         ' 1 F         T HZ     ', & ! 17 Sigma WW
                                         ' 1 F         T HZ     ', & ! 18 Tau WW
                                         ' 1 F         T HZh    ', & ! 19 Mean Flow U
                                         ' 1 F         T HZh    ', & ! 20 Mean Flow V
                                         ' 1 F         T HZ     ', & ! 21 Mean Flow W
                                         ' 1 F         T HZ     ', & ! 22 Temperature (K)
                                         ' 1 F         T HZ     ', & ! 23 Potential Temperature (K)
                                         ' 1 F         T HZ     ', & ! 24 Specific Humidity
                                         ' 1 F         T HZ     ', & ! 25 Pressure (Pa)
                                         ' 1 F         T HZ     ', & ! 26 Density
                                         ' 1 F         T HZ     ', & ! 27 Topography
                                         ' 1 F         T HZ     ', & ! 28 u-star
                                         ' 1 F         T HZ     ', & ! 29 Sensible Heat Flux
                                         ' 1 F         T HZ     ', & ! 30 Boundary Layer Depth
                                         ' 1 F         T HZ     ', & ! 31 Wind Speed
                                         ' 1 F         T HZh    ', & ! 32 Wind Direction (degrees)
                                         ' 1 F         T HZ     ', & ! 33 Precipitation Rate (mm/hr)
                                         ' 1 F         T HZ     ', & ! 34 Temperature (C)
                                         ' 1 F         T HZ     ', & ! 35 Cloud Amount (oktas)
                                         ' 1 F         T HZ     ', & ! 36 Relative Humidity (%)
                                         ' 1 F         T HZ     ', & ! 37 Pasquill Stability
                                         ' 1 L   8 U   T HZ  srd', & ! 38 # Particles By Species
                                         ' 1 O         T        ', & ! 39 Progress (%)
                                         ' 1 O   :     T        ', & ! 40 Clock Time
                                         ' 1 O           H h    ', & ! 41 X
                                         ' 1 O           H h    ', & ! 42 Y
                                         ' 1 F         T HZh    ', & ! 43 Sigma VV
                                         ' 1 F         T HZh    ', & ! 44 Mesoscale Sigma VV
                                         ' 1 F         T HZ     ', & ! 45 Cloud Water (kg/kg)
                                         ' 1 F         T HZ     ', & ! 46 Cloud Ice (kg/kg)
                                         ' 1 F         T HZ     ', & ! 47 3d Cloud (Fraction)
                                         ' 1 F         T HZ     ', & ! 48 Roughness Length
                                         ' 1 F         T HZ     ', & ! 49 Sea Level Pressure (Pa)
                                         ' 1 L         T HZ  sr ', & ! 50 Photon Flux
                                         ' 1 P 0       T HZ  sr ', & ! 51 Adult Effective Cloud Gamma Dose
                                         ' 1 P 0       T HZ  sr ', & ! 52 Adult Lung Cloud Gamma Dose
                                         ' 1 P 0       T HZ  sr ', & ! 53 Adult Thyroid Cloud Gamma Dose
                                         ' 1 P 0       T HZ  sr ', & ! 54 Adult Bone Surface Cloud Gamma Dose
                                         ' 1 P 1       T HZ  sr ', & ! 55 Area at risk
                                         ' 9 F         T HZ     ', & ! 56 Land Use Fractions
                                         ' 5 F         T HZ     ', & ! 57 Canopy Water
                                         ' 5 F         T HZ     ', & ! 58 Leaf Area Index
                                         ' 5 F         T HZ     ', & ! 59 Canopy Height
                                         ' 5 F         T HZ     ', & ! 60 Stomatal Conductance
                                         ' 1 F         T HZ     ', & ! 61 Soil Moisture
                                         ' 1 F         T HZ     ', & ! 62 Land Fraction
                                         ' 1 P 0       T HZ  srd', & ! 63 Mixing Ratio
                                         ' 1 E         T HZ  s  ', & ! 64 Eulerian Concentration
                                         ' 1 P 0       T HZ  s  ', & ! 65 E Mixing Ratio
                                         ' 1 P 0       T HZ  s  ', & ! 66 Concentration
                                         ' 1 F         T HZ     ', & ! 67 Convective Cloud Base
                                         ' 1 F         T HZ     ', & ! 68 Convective Cloud Top
                                         ' 1 E         T H   s  ', & ! 69 Eulerian Total Deposition Rate
                                         ' 1 E         T H   s  ', & ! 70 Eulerian Dry Deposition Rate
                                         ' 1 E         T H   s  ', & ! 71 Eulerian Wet Deposition Rate
                                         ' 1 A         T     sr ', & ! 72 Original Source Strength
                                         ' 1 A         T     sr '  & ! 73 Revised Source Strength
                                       /)
Character(MaxCharLength), Parameter :: QUnits(nQuantities) = &
                                       (/                    &
                                         'Variable ',        & !  1 Air Concentration
                                         'Variable ',        & !  2 Dry Deposition Rate
                                         'Variable ',        & !  3 Wet Deposition Rate
                                         'Variable ',        & !  4 Deposition Rate
                                         '         ',        & !  5 # Particles
                                         '         ',        & !  6 # Puffs
                                         '         ',        & !  7 # Particle Steps
                                         '         ',        & !  8 # Puff Steps
                                         'Variable ',        & !  9 Mass
                                         'Variable ',        & ! 10 Mean Z
                                         'Variable ',        & ! 11 Sigma Z
                                         'Variable ',        & ! 12 X Stats
                                         's        ',        & ! 13 Mean Travel Time
                                         'Variable ',        & ! 14 Puff Centres
                                         'Variable ',        & ! 15 Sigma C
                                         'Variable ',        & ! 16 Chemistry Field
                                         'm^2 / s^2',        & ! 17 Sigma WW
                                         's        ',        & ! 18 Tau WW
                                         'm / s    ',        & ! 19 Mean Flow U
                                         'm / s    ',        & ! 20 Mean Flow V
                                         'm / s    ',        & ! 21 Mean Flow W
                                         'K        ',        & ! 22 Temperature (K)
                                         'K        ',        & ! 23 Potential Temperature (K)
                                         'kg / kg  ',        & ! 24 Specific Humidity
                                         'Pa       ',        & ! 25 Pressure (Pa)
                                         'kg / m^3 ',        & ! 26 Density
                                         'm        ',        & ! 27 Topography
                                         'm / s    ',        & ! 28 u-star
                                         'W / m^2  ',        & ! 29 Sensible Heat Flux
                                         'm        ',        & ! 30 Boundary Layer Depth
                                         'm / s    ',        & ! 31 Wind Speed
                                         'deg      ',        & ! 32 Wind Direction (degrees)
                                         'mm / hr  ',        & ! 33 Precipitation Rate (mm/hr)
                                         'deg C    ',        & ! 34 Temperature (C)
                                         'oktas    ',        & ! 35 Cloud Amount (oktas)
                                         '%        ',        & ! 36 Relative Humidity (%)
                                         '         ',        & ! 37 Pasquill Stability
                                         '         ',        & ! 38 # Particles By Species
                                         '%        ',        & ! 39 Progress (%)
                                         '         ',        & ! 40 Clock Time
                                         'Variable ',        & ! 41 X
                                         'Variable ',        & ! 42 Y
                                         'm^2 / s^2',        & ! 43 Sigma VV
                                         'm^2 / s^2',        & ! 44 Mesoscale Sigma VV
                                         'kg / kg  ',        & ! 45 Cloud Water (kg/kg)
                                         'kg / kg  ',        & ! 46 Cloud Ice (kg/kg)
                                         'fraction ',        & ! 47 3d Cloud (Fraction)
                                         'm        ',        & ! 48 Roughness Length
                                         'Pa       ',        & ! 49 Sea Level Pressure (Pa)
                                         'm^-2 s^-1',        & ! 50 Photon Flux
                                         'Sv s^-1  ',        & ! 51 Adult Effective Cloud Gamma Dose
                                         'Gy s^-1  ',        & ! 52 Adult Lung Cloud Gamma Dose
                                         'Gy s^-1  ',        & ! 53 Adult Thyroid Cloud Gamma Dose
                                         'Gy s^-1  ',        & ! 54 Adult Bone Surface Cloud Gamma Dose
                                         'N/A      ',        & ! 55 Area at risk
                                         '         ',        & ! 56 Land Use Fractions
                                         'kg / m^2 ',        & ! 57 Canopy Water
                                         '         ',        & ! 58 Leaf Area Index
                                         'm        ',        & ! 59 Canopy Height
                                         'm / s    ',        & ! 60 Stomatal Conductance
                                         'kg / m^2 ',        & ! 61 Soil Moisture
                                         '         ',        & ! 62 Land Fraction
                                         'Variable ',        & ! 63 Mixing Ratio
                                         'Variable ',        & ! 64 Eulerian Concentration
                                         'Variable ',        & ! 65 E Mixing Ratio
                                         'Variable ',        & ! 66 Concentration
                                         'm agl    ',        & ! 67 Convective Cloud Base in m agl
                                         'm agl    ',        & ! 68 Convective Cloud Top in m agl
                                         'Variable ',        & ! 69 Eulerian Total Deposition Rate
                                         'Variable ',        & ! 70 Eulerian Dry Deposition Rate
                                         'Variable ',        & ! 71 Eulerian Wet Deposition Rate
                                         'Variable ',        & ! 72 Original Source Strength
                                         'Variable '         & ! 73 Revised Source Strength
                                       /)
Character(MaxCharLength), Parameter :: XStatsNames(10) = &
                                       (/                &
                                         "Mass        ", & !  1
                                         "<X>         ", & !  2
                                         "<Y>         ", & !  3
                                         "<Z>         ", & !  4
                                         "<X'X'>      ", & !  5
                                         "<Y'Y'>      ", & !  6
                                         "<Z'Z'>      ", & !  7
                                         "<X'Y'>      ", & !  8
                                         "<Y'Z'>      ", & !  9
                                         "<X'Z'>      "  & ! 10
                                       /)
! nQuantities :: Number of quantities useable in field requirements.
! Quantities  :: Names of quantities useable in field requirements. These are mostly self explanatory, but the
!                following notes should cover cases where this is not so.
!                  X Stats: The components of this are
!                             (mass,<x>,<y>,<z>,<x'x'>,<y'y'>,<z'z'>,<x'y'>,<y'z'>,<z'x'>).
!                           It is initially calculated as
!                             mass*(1,<x>,<y>,<z>,<xx>,<yy>,<zz>,<xy>,<yz>,<zx>)
!                           and then processed to get final values.
!                  Puff Centres: Air concentration calculated as if the puffs are particles located at the
!                                puff centres.
!                  $$ more notes on other fields?
!                  $$ table of processing chains (not 'processing steps' but calc of type PQR fields) with
!                     comments on any extra processes introduced and a comment that:
!                     Although extra processes can be introduced at this stage, they mustn't extend before the
!                     start of the original time range for F, O, Q, R quantities or after the end of the
!                     original time range for any quantities.
! QInfo       :: Table of information on the quantities.
!                  1st 2 characters - number of components
!                  4th character    - LEFOPQR-type (this takes values L, E, F, O, P, Q and R as defined in
!                                     comments at start of module)
!                  6th character    - for P, Q and R quantities, the number of averaging/integrating steps to
!                                     be applied (if present) before calculating the derived quantity
!                  8,:              - 8 means the type of the type used to store the results is Real(P64)
!                                     : means the type of the type used to store the results is Type(Time_)
!                                     otherwise the type of the type used to store the results is Real(Std)
!                  U,M              - U means spatial averaging/integrating is interpreted as summing
!                                     M means averaging/integrating is interpreted as mass weighted $$ remove
!                  f                - quantity for which fluctuations can be allowed for (but not a
!                                     fluctuation
!                                     parameter)
!                  p                - fluctuation parameter
!                  T                - quantity varies with T                       } at least before averaging
!                  S                - quantity varies with S (travel-time)         } or integrating
!                  H                - quantity varies with H (horizontal position) }
!                  Z                - quantity varies with Z                       }
!                  h                - quantity depends on choice of H-Coord
!                  z                - quantity depends on choice of Z-Coord
!                  s                - quantity varies with species (or species group)
!                  r                - quantity can be restricted by source (or source group).
!                  d                - quantity can be restricted by particle size.
! QUnits      :: Units in which the quantities are expressed. In some cases the unit is not given but is
!                determined in OutputFields. These cases are indicated by 'Variable'.
! XStatsNames :: Names of the components of XStats.

! Codes for the quantities in field requirements.
Integer, Parameter :: Q_AirConc                =  1 !} Codes for the quantities in field requirements.
Integer, Parameter :: Q_DryDep                 =  2 !}
Integer, Parameter :: Q_WetDep                 =  3 !}
Integer, Parameter :: Q_Dep                    =  4 !}
Integer, Parameter :: Q_nParticles             =  5 !}
Integer, Parameter :: Q_nParticlesBySpec       = 38 !}
Integer, Parameter :: Q_nPuffs                 =  6 !}
Integer, Parameter :: Q_nParticleSteps         =  7 !}
Integer, Parameter :: Q_nPuffSteps             =  8 !}
Integer, Parameter :: Q_Mass                   =  9 !}
Integer, Parameter :: Q_MeanZ                  = 10 !}
Integer, Parameter :: Q_SigmaZ                 = 11 !}
Integer, Parameter :: Q_XStats                 = 12 !}
Integer, Parameter :: Q_MeanS                  = 13 !}
Integer, Parameter :: Q_PuffCentres            = 14 !}
Integer, Parameter :: Q_SigmaC                 = 15 !}
Integer, Parameter :: Q_ChemistryField         = 16 !}
Integer, Parameter :: Q_SigmaWW                = 17 !}
Integer, Parameter :: Q_TauWW                  = 18 !}
Integer, Parameter :: Q_MeanFlowU              = 19 !}
Integer, Parameter :: Q_MeanFlowV              = 20 !}
Integer, Parameter :: Q_MeanFlowW              = 21 !}
Integer, Parameter :: Q_TemperatureK           = 22 !}
Integer, Parameter :: Q_PotentialTempK         = 23 !}
Integer, Parameter :: Q_SpecificHumidity       = 24 !}
Integer, Parameter :: Q_PressurePa             = 25 !}
Integer, Parameter :: Q_Density                = 26 !}
Integer, Parameter :: Q_Topography             = 27 !}
Integer, Parameter :: Q_UStar                  = 28 !}
Integer, Parameter :: Q_HeatFlux               = 29 !}
Integer, Parameter :: Q_BLDepth                = 30 !}
Integer, Parameter :: Q_WindSpeed              = 31 !}
Integer, Parameter :: Q_WindDirectionDeg       = 32 !}
Integer, Parameter :: Q_PptRateMMHR            = 33 !}
Integer, Parameter :: Q_TemperatureC           = 34 !}
Integer, Parameter :: Q_CloudOktas             = 35 !}
Integer, Parameter :: Q_RHPercent              = 36 !}
Integer, Parameter :: Q_Pasquill               = 37 !}
Integer, Parameter :: Q_ProgressPercent        = 39 !}
Integer, Parameter :: Q_ClockTime              = 40 !}
Integer, Parameter :: Q_X                      = 41 !}
Integer, Parameter :: Q_Y                      = 42 !}
Integer, Parameter :: Q_SigmaVV                = 43 !}
Integer, Parameter :: Q_MesoscaleSigmaVV       = 44 !}
Integer, Parameter :: Q_TotalOrDynCloudWater   = 45 !}
Integer, Parameter :: Q_TotalOrDynCloudIce     = 46 !}
Integer, Parameter :: Q_Cloud3d                = 47 !}
Integer, Parameter :: Q_RoughnessLength        = 48 !}
Integer, Parameter :: Q_PSeaLevelPa            = 49 !}
Integer, Parameter :: Q_PhotonFlux             = 50 !}
Integer, Parameter :: Q_AduEffCloudGammaDose   = 51 !}
Integer, Parameter :: Q_AduLunCloudGammaDose   = 52 !}
Integer, Parameter :: Q_AduThyCloudGammaDose   = 53 !}
Integer, Parameter :: Q_AduBoSCloudGammaDose   = 54 !}
Integer, Parameter :: Q_AreaAtRisk             = 55 !}
Integer, Parameter :: Q_LandUseFracs           = 56 !}
Integer, Parameter :: Q_CanopyWater            = 57 !}
Integer, Parameter :: Q_LAI                    = 58 !}
Integer, Parameter :: Q_CanopyHeight           = 59 !}
Integer, Parameter :: Q_StomataConduct         = 60 !}
Integer, Parameter :: Q_SoilMoisture           = 61 !}
Integer, Parameter :: Q_LandFrac               = 62 !}
Integer, Parameter :: Q_MixingRatio            = 63 !}
Integer, Parameter :: Q_EulerianConcentration  = 64 !}
Integer, Parameter :: Q_EMixingRatio           = 65 !}
Integer, Parameter :: Q_Concentration          = 66 !}
Integer, Parameter :: Q_ConCloudBase           = 67 !}
Integer, Parameter :: Q_ConCloudTop            = 68 !}
Integer, Parameter :: Q_EulerianTotalDep       = 69 !}
Integer, Parameter :: Q_EulerianDryDep         = 70 !}
Integer, Parameter :: Q_EulerianWetDep         = 71 !}
Integer, Parameter :: Q_OriginalSourceStrength = 72 !}
Integer, Parameter :: Q_RevisedSourceStrength  = 73 !}

! Codes for processing results.
Integer, Parameter :: P_AvInt          = 1 ! Code for calculating averages or integrals.
Integer, Parameter :: P_Prob           = 2 ! Code for calculating exceedence probabilities.
Integer, Parameter :: P_Percent        = 3 ! Code for calculating percentiles.
Integer, Parameter :: P_MaxMin         = 4 ! Code for calculating max or min.
Integer, Parameter :: P_Moments        = 5 ! Code for calculating moments about zero.
Integer, Parameter :: P_CentralMoments = 6 ! Code for calculating central moments.
Integer, Parameter :: P_Av             = 7 ! Code for averaging.
Integer, Parameter :: P_Int            = 8 ! Code for integrating.

! Codes for output information on concentration pdfs.
Integer, Parameter :: Pdf_ExceedenceProbs = 1 ! Probabilities of exceeding specified
                                              ! concentration thresholds.
Integer, Parameter :: Pdf_Percentiles     = 2 ! Percentiles of concentration.

! Codes for the mode of pdf calculation.
Integer, Parameter :: Pdf_AutoMode = 1 ! Default threshold values.
Integer, Parameter :: Pdf_UserMode = 2 ! User-specified threshold values.

!-------------------------------------------------------------------------------------------------------------

Type :: OutputOpts_ ! Output options.
  Logical                      :: Initialised   ! Indicates output options have been initialised.
  Character(MaxFileNameLength) :: Folder        ! Folder used for output (full or relative path). Ends in '\'
                                                ! or '/'.
  Logical                      :: Seconds       ! Indicates seconds are to be included in time output.
  Integer                      :: DecimalPlaces ! Number of decimal places to be included for seconds in time
                                                ! output.
  Logical                      :: Pre65Format   ! $$ remove this when not needed
End Type OutputOpts_

!-------------------------------------------------------------------------------------------------------------

Type :: Process_ ! A processing step.

  Character(MaxCharLength) :: Name      ! Name of processing step.
  Integer                  :: Code      ! Code to indicate averaging/integrating, calculating exceedence
                                        ! probabilities, calculating percentiles, calculating max/min,
                                        ! calculating moments about zero, or calculating central moments.

  Logical                  :: Ensemble  ! Indicates whether the processing is over the ensemble of cases.

  Integer                  :: TAvInt    ! Code to indicate, for averaging/integrating, whether the time
                                        ! processing is averaging, integrating or nothing.
  Integer                  :: nT        !} Number of times contributing to the processing, total duration
  Type(Time_)              :: T         !} contributing to the processing, and time interval between times
  Type(Time_)              :: dT        !} contributing to the processing.
  Logical                  :: UseTGrid  ! Indicates T and dT are determined by the T-grid and nT, with T = the
                                        ! (finite) processing time implied by the T grid and dT = T / nT (but
                                        ! these values are not stored in this structure).
  Type(Time_)              :: T0        ! Time origin for contributions when time processing is over all T.
  Type(ShortTime_)         :: sT        !] Short time versions of T, dT and T0.
  Type(ShortTime_)         :: sdT       !]
  Type(ShortTime_)         :: sT0       !]

  Integer                  :: HAvInt    !} Quantities analogous to TAvInt, nT, T, dT and UseTGrid for
  Integer                  :: nX        !} horizontal processing.
  Integer                  :: nY        !}
  Real(Pos)                :: X         !}
  Real(Pos)                :: Y         !}
  Real(Pos)                :: dX        !}
  Real(Pos)                :: dY        !}
  Logical                  :: UseHGrid  !}

  Integer                  :: ZAvInt    !] Quantities analogous to TAvInt, nT, T, dT and UseTGrid for vertical
  Integer                  :: nZ        !] processing.
  Real(Pos)                :: Z         !]
  Real(Pos)                :: dZ        !]
  Logical                  :: UseZGrid  !]
  Logical                  :: BL        ! Indicates that the vertical processing is over the boundary layer.

  Character(MaxCharLength) :: DGridName !} Name and index of D-Grid used for calculating exceedence
  Integer                  :: iDGrid    !} probabilities, calculating percentiles, calculating max/min,
                                        !} calculating moments about zero, or calculating central moments.
                                        !} Set to blank/zero if averaging integrating.

  ! For time processing the possibilities are as follows:
  !                              TGrid?    TAvInt for Code =      nT     T         dT     UseTGrid     T0
  !                                       P_AvInt   Otherwise
  ! ----------------------------------------------------------------------------------------------------------
  ! Depending on T:
  ! T window defined by grid       Yes   P_Av/P_Int     0*       > 0     0*        0*        Yes    ref-time*
  ! Rolling T window               Yes   P_Av/P_Int     0*     { > 0  dT * nT     > 0 }      No     ref-time*
  !                                                            {  1      0         0  }
  ! T window from T = -infinity    Yes   P_Av/P_Int     0*        0*  infinity    > 0        No     ref-time*
  ! T window from T = -infinity    Yes   P_Av/P_Int     -         1   infinity  infinity     No     ref-time*
  ! T window over all T            No    P_Av/P_Int     0*        0*  infinity    > 0        No*       any
  ! T window over all T            No    P_Av/P_Int     -         1   infinity  infinity     No*    ref-time*
  !
  ! Not depending on T:            No       P_Av        0*        1      0         0         No     ref-time*
  !
  ! (For no T averaging use:     Yes/No     P_Av        -         1      0         0         No     ref-time*)
  ! ----------------------------------------------------------------------------------------------------------
  !
  ! For horizontal processing the possibilities are as follows:
  !                           HGrid?    HAvInt for Code =     nX/Y      X/Y        dX/Y    UseHGrid
  !                                    P_AvInt   Otherwise
  ! -----------------------------------------------------------------------------------------------
  ! Depending on H:
  ! H window defined by grid    Yes   P_Av/P_Int     0*       > 0        0*         0*        Yes
  ! Rolling H window            Yes   P_Av/P_Int     0*     { > 0   dX/Y * nX/Y    > 0 }      No
  !                                                         {  1         0          0  }
  ! H window over all H         No    P_Av/P_Int     -         1     infinity    infinity     No*
  !
  ! Not depending on H:         No       P_Av        0*        1         0          0         No
  !
  ! (For no H averaging use:  Yes/No     P_Av        -         1         0          0         No)
  ! -----------------------------------------------------------------------------------------------
  !
  ! For vertical processing the possibilities are as follows:
  !                              ZGrid?    ZAvInt for Code =      nZ      Z        dZ     UseZGrid  BL
  !                                       P_AvInt   Otherwise
  ! ---------------------------------------------------------------------------------------------------
  ! Depending on Z:
  ! Z window defined by grid       Yes   P_Av/P_Int     0*       > 0      0*       0*        Yes    No
  ! Rolling Z window               Yes   P_Av/P_Int     0*     { > 0   dZ * nZ    > 0 }      No     No
  !                                                            {  1       0        0  }
  ! Z window defined by b-layer    No    P_Av/P_Int     0*       > 0      0*       0*        No     Yes
  ! Z window over all Z            No    P_Av/P_Int     -         1   infinity  infinity     No*    No
  !
  ! Not depending on Z:            No       P_Av        0*        1       0        0         No     No
  !
  ! (For no Z averaging use:     Yes/No     P_Av        -         1       0        0         No     No)
  ! ---------------------------------------------------------------------------------------------------
  !
  ! Note (1) * implies value chosen is not physically meaningful but it is computationally convenient to set a
  !          particular value.
  !      (2) > 0 implies finite.
  !      (3) Infinite values of X/Y/Z/dX/dY/dZ are stored as -1.
  !      (4) TGrid?, HGrid?, ZGrid? refers to the presence of the grid in the field requirement (i.e. it
  !          refers to the output of the processing step).
  !      (5) 'Depending on T/H/Z' refers to the input to the processing step. This input might not depend on
  !          T/H/Z due to a previous processing step or due to the underlying quantity. In this case an
  !          averaging/integrating processing step has no effect, but other types of process still apply.
  !      (6) Some cases are not allowed.
  !            (a) nT/X/Y/Z = 1, T/X/Y/Z = dT/X/Y/Z = infinity is not allowed except for intrinsic
  !                averaging/integrating (the other possible uses - to give a fluctuations
  !                averaging/integrating region, to multiply by the integrating region, and to do nothing -
  !                are of no value here, although the first of these is still possible in combination with
  !                intrinsic averaging/integrating).
  !            (b) nT/X/Y/Z = 1, T/X/Y/Z = dT/X/Y/Z > 0 is not allowed except for averaging/integrating.
  !            (c) Averaging/integrating over infinite regions of space or time is not allowed for quantities
  !                of type F, O, Q or R. It can however be applied to quantities of other types used to
  !                derived the type F/O/Q/R quantity.
  !                $$ for example integrated deposition rate treated as an instantaneous quantity (total dep)
  !                for further processing.
  !      (7) X/Y/Z and dX/dY/dZ are defined in the H/ZGrids coord system.
  !      (8) For dT > 0, the times used are the nominal time - i * dT for i = 0, nT - 1.
  !      (9) For dX > 0, the x-values used are the nominal X + i * dX for i = -nX/2 + 1/2, nX/2 - 1/2.
  !          Similarly for Y/Z. If BL is true, the processing is over the boundary layer in the vertical,
  !          and the z-values used are i * H/nZ for i = 1/2, nZ - 1/2 where H is the boundary layer depth and
  !          the m agl coord system is used.
  ! See also comments at start of module.

End Type Process_

!-------------------------------------------------------------------------------------------------------------

Type :: FieldReq_ ! A field requirement.

  ! Quantities defining the field.
  Character(MaxCharLength) :: Name
  Character(MaxCharLength) :: Quantity
  Integer                  :: iQuantity
  Character(1)             :: LEFOPQRType
  Character(1)             :: TypeType
  Integer                  :: nComponents
  Integer                  :: nRanges
  Character(MaxCharLength) :: SpeciesName
  Integer                  :: iSpecies
  Integer                  :: iParticleSpecies
  Character(MaxCharLength) :: MaterialUnitName
  Integer                  :: iMaterialUnit
  Character(MaxCharLength) :: SourceName
  Character(MaxCharLength) :: SourceGroupName
  Integer                  :: iSource
  Integer                  :: iSourceGroup
  Character(MaxCharLength) :: SizeDistName
  Integer                  :: iSizeDist
  Logical                  :: DecayDep
  Logical                  :: SemiInfApprox
  Character(MaxCharLength) :: TGridName
  Character(MaxCharLength) :: SGridName
  Character(MaxCharLength) :: HGridName
  Character(MaxCharLength) :: ZGridName
  Character(MaxCharLength) :: HCoordName
  Character(MaxCharLength) :: ZCoordName
  Integer                  :: iTGrid
  Integer                  :: iSGrid
  Integer                  :: iHGrid
  Integer                  :: iZGrid
  Integer                  :: iHCoord
  Integer                  :: iZCoord
  Integer                  :: iHGridCoord
  Integer                  :: iZGridCoord
  Integer                  :: nProcesses
  Character(MaxCharLength) :: ProcessNames(MaxProcessesPerFieldReq)
  Integer                  :: iProcesses(MaxProcessesPerFieldReq)
  Logical                  :: IntrinsicTAv
  Logical                  :: IntrinsicHAv
  Logical                  :: IntrinsicZAv
  Logical                  :: Fluctuations
  Type(ShortTime_)         :: FlAvT
  Real(Std)                :: FlAvX ! $$ ought to add metres/grid coord choices
  Real(Std)                :: FlAvY ! $$ ought to add metres/grid coord choices
  Real(Std)                :: FlAvZ ! $$ ought to add metres/grid coord choices
  Logical                  :: UseFlAvT
  Logical                  :: UseFlAvH
  Logical                  :: UseFlAvZ
  Logical                  :: Sync
  ! Name             :: Name of field requirement.
  ! Quantity         :: Quantity required.
  ! iQuantity        :: Index of quantity.
  ! LEFOPQRType      :: LEFOPQR-type (as defined in comments at start of module). Equals L, E, F, O, P, Q or R.
  ! TypeType         :: Indicates the type of type used to store the field. Values are as follows:
  !                       blank = Real(Std)
  !                       8     = Real(P64)
  !                       :     = Type(Time_).
  ! nComponents      :: Number of components of the field.
  ! nRanges          :: Number of particle size ranges in the field.
  ! SpeciesName      :: Name of species (blank if field not species dependent).
  ! iSpecies         :: Index of species in set of all species (0 if field not species dependent).
  ! iParticleSpecies :: Index of species in set of species on particles (0 if field not species dependent or species not on 
  !                     particles).
  ! MaterialUnitName :: Name of material unit associated with field request.
  ! iMaterialUnit    :: Index of material unit associated with field request.
  ! SourceName       :} Results are due to the specified source or source-group only. Blank indicates no
  ! SourceGroupName  :} restriction by source or source-group. Cannot both be set.
  ! iSource          :] Indices of the source and source-group named in SourceName and SourceGroupName. 0
  ! iSourceGroup     :] indicates no restriction by source or source-group.
  ! SizeDistName     :: Name of particle size distribution giving the particle size ranges to be used in
  !                     calculating the field. Blank indicates no restriction by particle size.
  ! iSizeDist        :: Index of particle size distribution. 0 indicates no restriction by particle size.
  ! DecayDep         :: Indicates that the deposits made during an averaging or integrating time interval are
  !                     to be decayed to account for the time till the end of the interval. Must be false
  !                     unless Quantity is 'Deposition Rate', 'Dry Deposition Rate' or 'Wet Deposition Rate'.
  ! SemiInfApprox    :: Indicates that cloud gamma doses are to be calculated using the semi-infinite cloud
  !                     approximation. Must be false unless Quantity is 'Adult Effective Cloud Gamma Dose',
  !                     'Adult Lung Cloud Gamma Dose', 'Adult Thyroid Cloud Gamma Dose' or 'Adult Bone Surface
  !                     Cloud Gamma Dose'.
  ! TGridName        :} Name of time/travel-time/horizontal/vertical grids if used and blank otherwise.
  ! SGridName        :}
  ! HGridName        :}
  ! ZGridName        :}
  ! HCoordName       :] Name of horizontal/vertical coord systems if used, and blank otherwise. If a coord name
  ! ZCoordName       :] is blank, the field does not depend on the choice of coord system.
  ! iTGrid           :} Indices of grids in Grids (external to this module). 0 indicates no grid used.
  ! iSGrid           :}
  ! iHGrid           :}
  ! iZGrid           :}
  ! iHCoord          :] Indices of coord systems in Coords (external to this module). 0 indicates no coord
  ! iZCoord          :] systems used.
  ! iHGridCoord      :} Indices in Coords of the coord systems associated with the horizontal and vertical
  ! iZGridCoord      :} grids. If AvBL is true, iZGridCoord is used to indicate the index of the m agl coord
  !                     system. 0 indicates no such coord system.
  ! nProcesses       :: Number of processing steps.
  ! ProcessNames     :} Names and indices of the processing steps.
  ! iProcesses       :}
  ! IntrinsicTAv     :] Indicates intrinsic averaging/integrating over T/H/Z in first processing step. Must be
  ! IntrinsicHAv     :] false for type E, F or O quantities.
  ! IntrinsicZAv     :]
  ! Fluctuations     :: Indicates fluctuations are accounted for in the first probability or moment calculation
  !                     step and/or in the calculation of fluctuation parameters.
  ! FlAvT            :} Averaging scales for fluctuations (set to zero if UseFlAvT/H/Z is false).
  ! FlAvX            :}
  ! FlAvY            :}
  ! FlAvZ            :}
  ! UseFlAvT         :] Indicates FlAvT/X/Y/Z are the values to be used (set to false if Fluctuations is
  ! UseFlAvH         :] false). If false and Fluctuations true, the scales are determined by the processing
  ! UseFlAvZ         :] steps.
  ! Sync             :: Results are calculated when the particles/puffs are synchronised in time. Sync must be
  !                     true for type E, F or O quantities and false for deposition quantities of type L or if
  !                     IntrinsicTAv is true.

  ! Quantities defining how the field is output.
  Logical                  :: Disk
  Logical                  :: Screen
  Character(MaxCharLength) :: OutputGroup
  Integer                  :: iOutputGroup
  Logical                  :: TAcross
  Logical                  :: SAcross
  Logical                  :: XAcross
  Logical                  :: YAcross
  Logical                  :: ZAcross
  Logical                  :: DAcross
  Logical                  :: TSeparateFile
  Logical                  :: SSeparateFile
  Logical                  :: XSeparateFile
  Logical                  :: YSeparateFile
  Logical                  :: ZSeparateFile
  Logical                  :: DSeparateFile
  Logical                  :: TNewFile
  Logical                  :: GridIndices
  Logical                  :: AlignCols
  Logical                  :: ZeroLines
  Logical                  :: Flush
  Logical                  :: NameII
  Logical                  :: Graph
  Real(Std)                :: XScale
  Real(Std)                :: YScale
  ! Disk          :} Indicates results are output numerically to disk and/or screen.
  ! Screen        :}
  ! OutputGroup   :: Name of output group for grouping numerical screen and disk output. Blank if no numerical
  !                  output.
  ! iOutputGroup  :: Index of output group. 0 if no numerical output.
  ! TAcross       :} Controls the splitting of the t,s,x,y,z,d dependence of output fields into rows, cols and
  ! SAcross       :} separate files. The 'across' settings are relevant even if separate files are selected
  ! XAcross       :} since they influence whether the t,s,x,y,z,d appears to the left of a row or above a
  ! YAcross       :} column. These settings can be used even if there is no t,s,x,y,z,d dependence (and this
  ! ZAcross       :} may be necessary to meet restrictions on variations between field requirements in the
  ! DAcross       :} same output group).
  ! TSeparateFile :}
  ! SSeparateFile :}
  ! XSeparateFile :}
  ! YSeparateFile :}
  ! ZSeparateFile :}
  ! DSeparateFile :}
  ! TNewFile      :: TNewFile causes new output files to be used after a restart (possibly overwriting any
  !                  previous files).
  ! GridIndices   :: Includes columns of grid indices for the non-across coords.
  ! AlignCols     :: Adds spaces to align columns.
  ! ZeroLines     :: Includes lines with all results zero.
  ! Flush         :: Flushes buffer after each output time.
  ! NameII        :: Makes output look more like Name II fields output if TAcross is true and like Name II
  !                  time series output if TAcross is false. The effects of this are as follows.
  !                    (1) Different file headers.
  !                    (2) Different number and order of the column-header lines. The lines used for TAcross
  !                        true are
  !                          Species category
  !                          Species
  !                          Time averaging/integrating
  !                          Quantity
  !                          Units
  !                          Z plus Z averaging information when Z is 'across'
  !                          T when this is 'across'
  !                          Blank line
  !                        and the lines used for TAcross false are
  !                          Y when this is 'across'
  !                          X when this is 'across'
  !                          X-Y location name when this is 'across'
  !                          Species category
  !                          Species
  !                          Quantity
  !                          Z plus Z averaging information when Z is 'across'
  !                          Units
  !                          Blank line.
  !                        Also YAcross must be true if TAcross is false in order to ensure the first column
  !                        header is actually present (if the first column header isn't present OutputFields
  !                        doesn't work).
  !                    (3) No case number information in file name.
  !                    (4) Time information in file name (if present) is e.g. 200101190600 instead of
  !                          T6_200101190600.
  !                    (5) Horizontal grid origin (for a structured grid) is output as (X0 - dX/2, Y0 - dY/2)
  !                        and grid indices are increased by 1/2.
  !                    (6) Altered column widths.
  !                    (7) Some differences in the content of column-header lines. Many of these differences
  !                        are inconsequential. However
  !                          (a) Time averaging/integrating times output to the nearest hour.
  !                          (b) 'Boundary layer' written for things which don't depend on Z.
  !                          (c) Heights output as 'From Z-dZ/2 - Z+dZ/2 unit' even if no vertical averaging.
  !                          (d) Time format changed (in particular UTC always used and no seconds given).
  !                          (e) 'Air concentration' written as 'Dosage' when time integrated.
  !                          (f) 'Max' is used as a prefix for the quantity when 100th percentiles are output.
  !                    (8) Different format for the output values in the preliminary and field columns.
  !                    (9) Probabilities and percentiles other than 100th percentiles and complex processing
  !                        are disabled.
  !                  The implications of this are as follows.
  !                    (2)  implies that
  !                           source restricted output can't be distinguished in the output files,
  !                           ensemble averaging information is not given,
  !                           time averaging information is not given (for T 'down' only),
  !                           horizontal averaging information is not given,
  !                           horizontal locations can't be distinguished unless X and Y are 'down' (for T
  !                             'across' only),
  !                           the period over which 100th percentiles are calculated is not given,
  !                           information on whether 100th percentiles are calculated over the ensemble is not
  !                             given.
  !                         It also implies, in conjunction with (7c), that vertical averaging information may
  !                         be incomplete or absent. Also YAcross must be true if TAcross is false in order to
  !                         ensure the first column header is actually present (if the first column header
  !                         isn't present OutputFields doesn't work).
  !                    (3)  can lead to output for different cases and output processed over the ensemble
  !                         overwriting each other.
  !                    (4)  can lead to output for different times overwriting each other if the times are
  !                         closer than a minute.
  !                    (7a) can lead to different averaging times being indistinguishable if the times are
  !                         closer than an hour.
  !                    (7d) can lead to different times being indistinguishable in the files if the times are
  !                         closer than a minute.
  !                    (9)  implies probabilities and percentiles other than 100th percentiles and complex
  !                         processing are disabled.
  !                  Hence the NameII flag is best avoided for:
  !                    (1) multi-case runs,
  !                    (2) > 1 choice of source restriction on the output (including none),
  !                    (3) output time grids with times closer than a minute,
  !                    (4) X or Y 'across' if > 1 horizontal grid point and T 'across',
  !                    (5) two choices of averaging/integrating time which are closer than an hour,
  !                    (6) > 1 choice of time averaging (including none) and T 'down',
  !                    (7) > 1 choice of horizontal averaging (including none),
  !                    (8) Choices of vertical averaging which can't be distinguished,
  !                    (9) > 1 choice of 100th percentile calculation time.
  !                  To make output look most like NameII 'fields output' or 'time series output' (and in
  !                  particular to work with the PV-Wave code), it is necessary to have the following for
  !                  fields output
  !                    (1) Separate File = 'T'
  !                    (2) Across = 'TZ' (D also needed for quantities depending on a data grid - must not be
  !                        a 'floating' D-Grid; can also include D if no D dependence)
  !                    (3) Output Format options to include I, A and 2
  !                    (4) Output Group beginning 'Fields_' (case sensitive) followed by a string which needs
  !                        to be specified as the selgrid variable when running PV-Wave
  !                    (5) A T-grid, a structured regular H-grid and no S-grid
  !                  and the following for time series output
  !                    (1) Separate File = 'XY' or blank
  !                    (2) Across = 'XYZ' (D also needed for quantities depending on a data grid - must not be
  !                        a 'floating' D-Grid; can also include D if no D dependence)
  !                    (3) Output Format options to include A, Z and 2 but not I
  !                    (4) Output Group beginning 'Time_series_' (case sensitive) followed by a string which
  !                        needs to be specified as the selgrid variable when running PV-Wave
  !                    (5) A T-grid, an unstructured H-grid of named points and no S-grid.
  !                  Where there are outputs which the PV-Wave code cannot distinguish, the PV-Wave code will
  !                  ignore all instances after the first (for both fields and time series). In addition to
  !                  applying to outputs which are indistinguishable because of the output format, the PV-Wave
  !                  code does not distinguish fields by the time averaging information or time series by the
  !                  X and Y values (but does use the X-Y location name).
  !                  $$ could correct for time averaging info
  !                  $$ keep these comments aligned with those in PVWave Readme.txt
  ! Graph         :: Indicates results are output graphically to the screen.
  ! XScale        :} Scales used for plotting the results graphically. 0 if no graphical output.
  ! YScale        :}

  ! Note all field requirements in the same output group must have the same values for T/S/X/Y/Z/DAcross,
  ! T/S/X/Y/Z/DSeparateFile, TNewFile, GridIndices, AlignCols, ZeroLines, Flush, NameII, Ensemble and, unless
  ! T/S/X/Y/Z/DAcross is true and T/S/X/Y/Z/DSeparateFile is false, T/S/H/Z/DGrid. Here 'the same values for
  ! T/S/H/Z/DGrid' can mean they all have no grid, although, for TGrid, they must all fail to have a TGrid in
  ! the same way (i.e. because the underlying quantity doesn't depend on time or because the fields are
  ! processed over all time). If a field requirement has an unstructured HGrid, then the field requirement
  ! must have XAcross = YAcross and XSeparateFile = YSeparateFile.

  ! Quantities derived from the processing steps used.
  Logical                  :: ComplexProc
  Logical                  :: Ensemble
  Logical                  :: AvEnsemble
  Integer                  :: AvTAvInt
  Type(Time_)              :: AvT
  Type(ShortTime_)         :: sAvT
  Logical                  :: AvUseTGrid
  Integer                  :: AvHAvInt
  Real(Std)                :: AvX
  Real(Std)                :: AvY
  Logical                  :: AvUseHGrid
  Integer                  :: AvZAvInt
  Real(Std)                :: AvZ
  Logical                  :: AvUseZGrid
  Logical                  :: AvBL
  Integer                  :: PCode
  Logical                  :: PEnsemble
  Type(Time_)              :: PT
  Type(ShortTime_)         :: sPT
  Logical                  :: PUseTGrid
  Real(Std)                :: PX
  Real(Std)                :: PY
  Logical                  :: PUseHGrid
  Real(Std)                :: PZ
  Logical                  :: PUseZGrid
  Logical                  :: PBL
  Character(MaxCharLength) :: DGridName
  Integer                  :: iDGrid
  Logical                  :: DGridFloating
  Integer                  :: iPQRProc
  Integer                  :: iFOQRProc
  Integer                  :: iFluctuationsProc
  Logical                  :: NoProc(MaxProcessesPerFieldReq)
  Integer                  :: iLastProcess
  ! ComplexProc       :: Indicates the processing is complex (as defined in comments at start of module).
  ! Ensemble          :: Indicates there is a processing step with Ensemble true. In other words the results
  !                      depend on the ensemble of cases rather than on an individual case (because they
  !                      involve either taking an ensemble average or calculating
  !                      probabilities/percentiles/moments across the ensemble of cases).
  ! AvEnsemble        :} The values of Ensemble, TAvInt, T, sT, UseTGrid, HAvInt, X, Y, UseHGrid, ZAvInt, Z,
  ! AvTAvInt          :} UseZGrid and BL from the first processing step if this is an averaging/integrating
  ! AvT               :} processing step, with the exception that T/H/ZAvInt is zero (instead of being P_Av
  ! sAvT              :} with a zero size to the averaging region) if there is no averaging/integrating in the
  ! AvUseTGrid        :} corresponding direction. For simple processing (as defined in comments at start of
  ! AvHAvInt          :} module) they are used to label the output (see comments at start of module). They are
  ! AvX               :} also used to control any intrinsic averaging/integrating and in calculating
  ! AvY               :} fluctuation averaging/integrating scales. Set to false, 0, 0.0 or 00:00 if no such
  ! AvUseHGrid        :} step.
  ! AvZAvInt          :}
  ! AvZ               :}
  ! AvUseZGrid        :}
  ! AvBL              :}
  ! PCode             :] For simple processing (as defined in comments at start of module), copy of Code,
  ! PEnsemble         :] Ensemble, T, sT, UseTGrid, X, Y, UseHGrid, Z, UseZGrid and BL from the probability
  ! PT                :] step without a following percentile step, the probability step with a following
  ! sPT               :] percentile step, the max/min step or the moment step to be used to label the output
  ! PUseTGrid         :] (see comments at start of module). For the case of a probability step with a
  ! PX                :] following percentile step, Code comes from the percentile step and the other
  ! PY                :] variables from the probability step. Set to false, 0, 0.0, 00:00 or blank if no such
  ! PUseHGrid         :] step.
  ! PZ                :]
  ! PUseZGrid         :]
  ! PBL               :]
  ! DGridName         :} Name and index of data grid used to label output (see comments at start of module).
  ! iDGrid            :} Blank and zero if no such grid.
  ! DGridFloating     :: Indicates whether the data grid DGridName is floating.
  ! iPQRProc          :: Index of the processing step (in ProcessNames and iProcesses) after which the derived
  !                      quantity is calculated for type P, Q and R quantities. 0 for type L, E, F or O
  !                      quantities and when the derived quantity is calculated before any processing steps.
  ! iFOQRProc         :: Index of the processing step after which quantities of type F, O, Q or R appear in
  !                      the processing chain. 0 for type L, E or P quantities and when quantities of type F,
  !                      O, Q or R appear before any steps in the processing chain.
  ! iFluctuationsProc :: Index of processing step (in ProcessNames and iProcesses) in which fluctuations are
  !                      accounted for. 0 if no such processing step.
  ! NoProc            :: Indicates which processing steps don't actually require processing (see comments at
  !                      start of module) or, for type P, Q and R fields, which processing steps only apply to
  !                      the fields from which the derived field is derived.
  ! iLastProcess      :: Index (in ProcessNames and iProcesses) of last processing step with NoProc false. 0
  !                      if no such processing step.

  ! Indices of other fields needed to process this field.
  Integer :: iExt         ! Index of a field which extends the field and which can be used instead. 0 if no
                          ! such field.
  Integer :: iProc        !} Index of field required for last processing step which actually requires 
                          !} processing
                          ! (i.e. has NoProc false). 0 if no such field.
  Integer :: iDryDep      !} For derived fields, indices of fields required to produce the derived field. 0 if
  Integer :: iWetDep      !} no such fields.
  Integer :: iAirConc     !}
  Integer :: iChemField   !}
  Integer :: iEulConc     !}
  Integer :: iDensity     !}
  Integer :: iMeanS       !}
  Integer :: iXStats      !}
  Integer :: iSigmaC      !}
  Integer :: iMass        !}
  Integer :: iMeanZ       !}
  Integer :: iPhotonFlux  !}
End Type FieldReq_

!-------------------------------------------------------------------------------------------------------------

Type :: PdfReq_ ! A pdf requirement.
  Character(MaxCharLength) :: Name
  Integer                  :: iQuantity
  Integer                  :: PdfType
  Character(MaxCharLength) :: SpeciesName
  Integer                  :: iSpecies
  Character(MaxCharLength) :: MaterialUnitName
  Integer                  :: iMaterialUnit
  Integer                  :: nSources
  Character(MaxCharLength) :: SourceName
  Character(MaxCharLength) :: SourceGroupName
  Integer                  :: iSource
  Integer                  :: iSourceGroup
  Character(MaxCharLength) :: HGridName
  Character(MaxCharLength) :: ZGridName
  Character(MaxCharLength) :: TGridName
  Character(MaxCharLength) :: HCoordName
  Character(MaxCharLength) :: ZCoordName
  Integer                  :: iHGrid
  Integer                  :: iZGrid
  Integer                  :: iTGrid
  Integer                  :: iHCoord
  Integer                  :: iZCoord
  Integer                  :: CalcType
  Integer                  :: TAvOrInt
  Real(Std)                :: AvTime
  Integer                  :: PdfMode
  Integer                  :: PdfSize
  Real(Std)                :: PdfThresholds(MaxPdfSize)
  Logical                  :: Screen
  Logical                  :: Disk
  Logical                  :: AvEnsemble
  Integer                  :: iOutputGroup
  Character(5)             :: SeparateFile
  Character(MaxCharLength) :: OutputGroup
  Integer                  :: iAirConc
  Integer                  :: iMeanS
  Integer                  :: iXStats
  Integer                  :: iMass
  Integer                  :: iMeanZ
  Integer                  :: iSigmaC
  ! Name         :: Name of pdf.
  ! Code         :: Code for pdf.
  ! PdfType      :: Type of pdf (1 = exceedence probabilities, 2 = percentiles).
  ! SpeciesName  :: Name of species.
  ! iSpecies     :: Index of species.
  ! SpeciesName  :: Name of material unit associated with pdf request
  ! iSpecies     :: Index of material unit associated with pdf request
  ! nSources     :: Number of sources in fluctuations calculation. 0 indicates all
  !                 sources are considered.
  ! SourceName      :} Results are due to the specified source or source-group (defined
  ! SourceGroupName :} only if nSources is not 0).
  ! iSource      :] Indices of the sources and source-group named in SourceName and
  ! iSourceGroup :] SourceGroupName (defined only if nSources is not 0).
  ! HGridName    :} Name of horizontal/vertical/temporal grids if used, blank
  ! ZGridName    :} otherwise. If blank, either the pdf has no variation in that
  ! TGridName    :} variable (e.g. continuous-release plume calculations have no
  !                 temporal variation) or it refers to a time-integrated calculation.
  ! HCoordName   :] Name of horizontal/vertical coord systems used.
  ! ZCoordName   :] (ZCoordName defined only if ZGridName is defined).
  !                 Note these are needed because, unlike in Field_, we need the names
  !                 to set up extra requirements.
  ! iHGrid       :} Index of grids and their coord systems in Grids and Coords
  ! iZGrid       :} (external to this module). 0 indicates no grid or coord system
  ! iTGrid       :} is used.
  ! iHCoord      :}
  ! iZCoord      :}
  ! CalcType     :: Indicates the type of fluctuations calculation to be carried out:
  !                 1 - continuous releases (with possibility of time-averaging),
  !                 2 - time-integrated quantities for finite-duration releases,
  !                 3 - instantaneous quantities for finite-duration releases.
  ! TAvOrInt     :: Indicates time-averaged concentrations are considered
  !                 (for type 1 calculation only).
  ! AvTime       :: Averaging time (in seconds) for the time-averaged concentrations
  !                 (for type 1 calculation only).
  ! PdfMode      :: Mode of specifying thresholds for the pdf calculation
  !                 (1 = auto mode using default values, 2 = user-specified values).
  ! PdfSize      :: Number of threshold values in user-specified input.
  ! PdfThresholds:: User-specified threshold values for the pdf calculation.
  ! Screen       :: Indicates results are output numerically to the screen.
  ! Disk         :: Indicates results are output numerically to disk.
  ! AvEnsemble   :: Indicates results are statistically processed.
  ! iOutputGroup :: Index of output group. 0 = unassigned. All pdfs in the same
  !                 output group are output to the same unit for numerical screen or disk
  !                 output.
  ! SeparateFile :: SeparateFile = 'T' will put different times in different files.
  ! OutputGroup  :: Name of output group - used for file names, titles, grouping fields.
  !                 $$ Need to read these in via input.

  ! For all indices 0 indicates that the index value is unknown.

  ! Note that SourceName and SourceGroupName (equally iSource and iSourceGroup) are
  ! mutually exclusive options. If one is specified, the other must be kept blank.

  ! The output group is not user controllable - instead all pdfs which use
  ! the same grid are assigned to the same output group.
  ! In 'processing' pdfs, only other pdfs from the same group can be used. $$ ?? still true?

  ! The PdfReq_ type is primarily designed to represent requirements for probability
  ! distributions of concentration fluctuations, although, in principle, it could be
  ! applicable to other distributions (e.g. those of deposition fields).

  ! Averaging and integrating over space is not supported in the fluctuations scheme.

  ! Results are always calculated when the particle/puffs are synchronised in time.
End Type PdfReq_

!-------------------------------------------------------------------------------------------------------------

Type :: PPInfoReq_ ! A requirement for a set of particle/puff information.
  Character(MaxCharLength) :: Name
  Character(1)             :: Particles
  Character(1)             :: Puffs
  Integer                  :: FirstParticle
  Integer                  :: LastParticle
  Integer                  :: FirstPuff
  Integer                  :: LastPuff
  Character(MaxCharLength) :: SourceName
  Integer                  :: iSource
  Logical                  :: Met
  Logical                  :: Mass
  Logical                  :: PlumeRise
  Logical                  :: DispersionScheme
  Logical                  :: PuffFamily
  Character(MaxCharLength) :: HCoordName
  Character(MaxCharLength) :: ZCoordName
  Character(MaxCharLength) :: TGridName
  Integer                  :: iHCoord
  Integer                  :: iZCoord
  Integer                  :: iTGrid
  Logical                  :: Sync
  Logical                  :: Disk
  Logical                  :: Screen
  Logical                  :: TSeparateFile
  Logical                  :: PSeparateFile
  Logical                  :: Flush
  Logical                  :: Graph
  ! Name             :: Name of requirement for a set of particle/puff information.
  ! Particles        :: A, S and N indicate that All, Some and No particles are included.
  ! Puffs            :: A, S and N indicate that All, Some and No puffs are included.
  ! FirstParticle    :} Output restricted to particle/puff numbers in this range. Zero indicates no
  ! LastParticle     :} restriction or particles/puffs not included.
  ! FirstPuff        :}
  ! LastPuff         :}
  ! SourceName       :: Results are due to the specified source only. Blank indicates no restriction by
  !                     source.
  ! iSource          :: Index of the source. 0 indicates no restriction by source.
  ! Met              :} Indicates information is included on met, mass, plume rise, the dispersion scheme, and
  ! Mass             :} (for puffs only) the puff family.
  ! PlumeRise        :}
  ! DispersionScheme :}
  ! PuffFamily       :}
  ! HCoordName       :] Name of horizontal/vertical coord systems used.
  ! ZCoordName       :]
  ! TGridName        :: Name of time grid if used, blank otherwise. If blank, results are output every time
  !                     step.
  ! iHCoord          :} Indices of coord systems in Coords (external to this module).
  ! iZCoord          :}
  ! iTGrid           :: Index of the time grid in Grids (external to this module). 0 indicates no grid used.
  ! Sync             :: Results are calculated when the particles/puffs are synchronised in time.
  ! Disk             :} Indicates results are output numerically to disk and/or screen.
  ! Screen           :}
  ! TSeparateFile    :] Controls the splitting of times and particles/puffs into separate files.
  ! PSeparateFile    :]
  ! Flush            :: Flushes buffer after each output time.
  ! Graph            :: Indicates results are output graphically to the screen.

  ! TSeparateFile is only possible if there is a time grid.
  ! PSeparateFile is only possible if TSeparateFile is false. $$ and any other SeparateFile,
  ! if these are added.

  ! Add S-grid options / S-Separate file. Need to add to appropriate CalcReqInfo routine & exclude from others
  ! and ensure calculated. Can't have T and S grid. $$
  ! Add domain restrictions. $$
  ! Restricting by iUOP as well as iUP would be useful. $$

End Type PPInfoReq_

!-------------------------------------------------------------------------------------------------------------

Type :: Reqs_ ! A collection of requirements.

  ! Processing steps used for processing fields.
  Integer        :: nProcesses              ! Number of processing steps.
  Type(Process_) :: Processes(MaxProcesses) ! Processing steps.

  ! Fields.
  Integer                  :: nFieldReqs      ! Number of field requirements.
  Integer                  :: nFieldGroups    ! Number of output groups of fields.
  Type(FieldReq_), Pointer :: FieldReqs(:)    ! Field requirements.
  Type(FieldReq_), Pointer :: FieldGroups(:)  ! Output field groups.
  
  Integer :: MaxFieldReqs              ! Maximum number of field requests.
  Integer :: MaxFieldOutputGroups      ! Maximum number of field output groups.

  ! Pdfs.
  Integer       :: nPdfs               ! Number of pdfs.
  Integer       :: nPdfGroups          ! Number of groups of pdfs.
  Type(PdfReq_) :: PdfReqs(MaxPdfReqs) ! Pdf requirements.

  ! Sets of particle/puff information.
  Integer          :: nPPInfoReqs               ! Number of requirements for sets of particle/puff
                                                ! information.
  Type(PPInfoReq_) :: PPInfoReqs(MaxPPInfoReqs) ! Requirements for sets of particle/puff information.

End Type Reqs_

!-------------------------------------------------------------------------------------------------------------

Type :: Field_ ! A field.

  ! If changing this type, remember WriteRestartFileResults and ReadRestartFileResults routines.

  Character(1)         :: TypeType
  Integer              :: nT
  Integer              :: iTL
  Real(Std),   Pointer :: Std(:, :, :, :, :, :, :)
  Real(P64),   Pointer :: P64(:, :, :, :, :, :, :)
  Type(Time_), Pointer :: T  (:, :, :, :, :, :, :)
  Real(Std),   Pointer :: MaxStd(:, :, :, :, :, :)
  Real(Std),   Pointer :: MinStd(:, :, :, :, :, :)
  Integer,     Pointer :: S     (:, :, :, :, :, :)
  Integer              :: LastProcThisCase
  Integer              :: LastProc
  Integer              :: LastNumOutput
  Integer              :: LastGraphOutput
  Integer              :: GraphUnit
  ! TypeType         :: Indicates the type of type used to store the field. Values are as follows:
  !                       blank = Real(Std)
  !                       8     = Real(P64)
  !                       :     = Type(Time_).
  ! nT               :: Number of time levels allocated to store the field.
  ! iTL              :: Index of lowest time level is stored.
  ! Std              :} The field values. The indices correspond to D, X, Y, Z, S, T and components/ranges of
  ! P64              :} the field. Usually Std will be used, but some fields use P64 or T (see QInfo table and
  ! T                :} TypeType).
  ! MaxStd           :] If the field is a probability distribution using a D-grid with more than one point or
  ! MinStd           :] a floating D-grid, this gives the max and min values of the underlying field.
  ! S                :: If the field is a probability distribution and a floating D-grid is used, this gives
  !                     the 'shift' or 'offset' associated with it. The indices correspond to X, Y, Z, S, T
  !                     and components of the field.
  ! LastProcThisCase :: Index of the last time to have been processed for the current case.
  ! LastProc         :: Index of the last time to have been processed completely. For fields processed over
  !                     the ensemble of cases, processing is only completed once contributions from all cases
  !                     have been made.
  ! LastNumOutput    :: Index of the last time to have been output numerically.
  ! LastGraphOutput  :: Index of the last time to have been output graphically.
  ! GraphUnit        :: Unit for the window which is used for graphical output. The window only counts if it
  !                     is open and OutputFields is executing, or if it has been opened by OutputFields and
  !                     may need to be written to again (note 'may' does not mean 'will'). Equivalently the
  !                     windows which count are the open ones minus windows that won't need to be written to
  !                     again but are left open. (Windows that might have needed to be written to again but
  !                     are currently closed because of a restart are not included. Any further output is
  !                     regarded as being to a new window, so the window doesn't need to be written to again).
  !                     If the window doesn't count, GraphUnit is set to zero.

  ! If time level iT is stored, it is stored with time index iTA in the array P/D/T where iTA = iT mod nT. The
  ! time levels stored are those time levels which lie in the range iTL to iTL + nT - 1. iTL and nT are chosen
  ! so that these always lie in the relevant time grid.

  ! For fields not depending on time, nT = iTL = 1 and LastProcThisCase, LastProc, LastNumOutput and
  ! LastGraphOutput are 0 or 1 depending on whether the field has been processed etc.

  ! When there is a field which extends the field and which can be used instead, the elements Std, P64, T,
  ! MaxStd, MinStd and S are not valid. However the other elements are valid and take the values they
  ! otherwise would.

End Type Field_

!-------------------------------------------------------------------------------------------------------------

!## Includes other information defining the pdf within this type (although
!## to some extent this is just replicating existing data in the requirements).
Type :: Pdf_ ! A pdf.
  Integer        :: PdfType                ! Type of pdf
                                           ! (1 = exceedence probs, 2 = percentiles).
  Integer        :: PdfMode                ! Mode of specifying thresholds
                                           ! (1 = auto-mode, 2 = user-specified).
  Logical        :: FixedThresholds        ! Indicates that concentration thresholds are
                                           ! fixed in the calculation of exceedence probs.
  Integer        :: PdfSize                ! Number of points calculated in each pdf.
  Real(Std)      :: Thresholds(MaxPdfSize) ! Fixed threshold values in pdf calculation
                                           ! (except exceedence probs in auto-mode).
  Integer, Pointer   :: Scale(:,:,:,:)     ! Index of lowest concentration threshold in
                                           ! auto-mode calculation of exceedence probs
                                           ! for each individual pdf (this is essentially
                                           ! a scale factor, see documentation).
  Real(Std), Pointer :: Data(:,:,:,:,:)    ! Pdf data.
End Type Pdf_

!-------------------------------------------------------------------------------------------------------------

Type :: PPInfo_ ! A set of particle/puff information. More correctly this should be called 'information about
                ! a set of particle/puff information' as the particle/puff information itself is output as
                ! soon as its available and hence, unlike in the analogous Field_ type, it is not stored here.

  ! If changing this type, remember WriteRestartFileResults and ReadRestartFileResults routines.

  Integer :: nLines       (MaxPPInfoOutputFiles)
  Integer :: DiskUnits    (MaxPPInfoOutputFiles)
  Integer :: ScreenUnits  (MaxPPInfoOutputFiles)
  Integer :: ParticleFiles(MaxPPInfoOutputFiles)
  Integer :: PuffFiles    (MaxPPInfoOutputFiles)
  Integer :: LastFile
  Integer :: nT
  Integer :: iTL
  Integer :: GraphUnit
  ! nLines        :} The number of lines written to file, the disk unit and the screen unit for a number of
  ! DiskUnits     :} the files/windows which are used for numerical output. Files/windows are only included if
  ! ScreenUnits   :} they have been opened by OutputPPInfos and may need to be written to again (note 'may'
  !                  does not mean 'will'). Equivalently the included files/windows are the open ones plus
  !                  files that have been opened by OutputPPInfos and may need to be written to again but are
  !                  currently closed due to a restart (for these files DiskUnits = 0, but nLines /= 0) minus
  !                  windows that won't need to be written to again but are left open. (Windows that might
  !                  have needed to be written to again but are currently closed because of a restart are not
  !                  included. Any further output is regarded as being to a new window, so the window doesn't
  !                  need to be written to again). The array index ranges over relevant files/windows (for
  !                  details of the index see the calculation of iFile in OutputPPInfos). Unused slots in the
  !                  arrays are filled with zeros.
  ! ParticleFiles :] When separate particles/puffs are output in separate files, these arrays give the indices
  ! PuffFiles     :] of the files used for the various particles/puffs. The array index is the unique index of
  !                  the particle/puff (+ 1 - PPInfoReq%FirstParticle/Puff if not all particles/puffs are to
  !                  be output).
  ! LastFile      :: Index of the file with the highest index which has been used so far.
  ! nT            :: Number of time levels which may currently be output. !} Significant for separate times in
  ! iTL           :: Index of lowest time which may currently be output.  !} separate files only.
  ! GraphUnit     :: Unit for the window which is used for graphical output. The window only counts if it has
  !                  been opened by OutputPPInfos and may need to be written to again (note 'may' does not
  !                  mean 'will'). Equivalently the windows which count are the open ones minus windows that
  !                  won't need to be written to again but are left open. (Windows that might have needed to
  !                  be written to again but are currently closed because of a restart are not included. Any
  !                  further output is regarded as being to a new window, so the window doesn't need to be
  !                  written to again). If the window doesn't count, GraphUnit is set to zero.

  ! A time level iT which may currently be output can be referred to by an adjusted index iTA in the range
  ! [1, nT] where iTA = iT mod nT. The time levels which may currently be output are those time levels which
  ! lie in the range iTL to iTL + nT - 1. iTL and nT are chosen so that these always lie in the relevant time
  ! grid. The adjusted index iTA (and hence nT and iTL) is only used for the case of separate times in
  ! separate files, where it is used to index the files/windows.

End Type PPInfo_

!-------------------------------------------------------------------------------------------------------------

Type :: Results_ ! A collection of results.

  ! If changing this type, remember WriteRestartFileResults and ReadRestartFileResults routines.

  ! Run information.
  Character(MaxCharLength) :: RunName        ! Name of run.
  Type(Time_)              :: StartClockTime ! Clock time of the start of the run.

  ! Fields.
  Integer                :: nFields
  Type(Field_), Pointer  :: Fields          (:)
  Integer, Pointer        :: FieldnLines     (:, :)
  Integer, Pointer        :: FieldDiskUnits  (:, :)
  Integer, Pointer        :: FieldScreenUnits(:, :)
  Type(ShortTime_)       :: MaxTimeField
  ! nFields          :: Number of fields.
  ! Fields           :: Fields (ordered to correspond to FieldReqs in Reqs).
  ! FieldnLines      :} The number of lines written to file, the disk unit and the screen unit for a number of
  ! FieldDiskUnits   :} the files/windows which are used for numerical output. Files/windows are only included
  ! FieldScreenUnits :} if they are open and OutputFields is executing, or if they have been opened by
  !                     OutputFields and may need to be written to again (note 'may' does not mean 'will').
  !                     Equivalently the included files/windows are the open ones plus files that have been
  !                     opened by OutputFields and may need to be written to again but are currently closed
  !                     due to a restart (for these files FieldDiskUnits = 0, but FieldnLines /= 0) minus
  !                     windows that won't need to be written to again but are left open. (Windows that might
  !                     have needed to be written to again but are currently closed because of a restart are
  !                     not included. Any further output is regarded as being to a new window, so the window
  !                     doesn't need to be written to again). The first array index ranges over output groups
  !                     and the second over relevant files/windows within the output group (for details of the
  !                     second index see the calculation of iFile in OutputFields). Unused slots in the arrays
  !                     are filled with zeros. These variables are stored here rather than in the Fields_ type
  !                     because the information is for an output group not a field.
  ! MaxTimeField     :: Maximum of the times for which processing and output of fields has been requested.
  !                     This is a maximum over all cases and is used to control the processing and output of
  !                     fields processed over the ensemble of cases.

  ! Pdfs.
  Integer                      :: nPdfs
  Type(Pdf_)                   :: Pdfs(MaxPdfReqs)
  Logical                      :: PdfsProc(9000, MaxPdfReqs)
                                  ! $$ 9000 hardwired for MaxTPoints now that MaxTPoints
                                  ! removed. In due course won't be needed.
  Integer                      :: PdfsNextOutput(MaxPdfReqs)
  Integer                      :: PdfScreenUnits(MaxPdfReqs)
  Integer                      :: PdfDiskUnits(MaxPdfReqs)
  Character(MaxFileNameLength) :: PdfFiles(MaxPdfReqs)
  ! nPdfs           :} See fields above, with 'fields' replaced by 'pdfs'.
  ! Pdfs            :}
  ! PdfsProc        :}
  ! PdfsNextOutput  :}
  ! PdfScreenUnits  :}
  ! PdfDiskUnits    :}
  ! PdfFiles        :}

  ! Sets of particle/puff information.
  Integer       :: nPPInfos               ! Number of sets of particle/puff information.
  Type(PPInfo_) :: PPInfos(MaxPPInfoReqs) ! Sets of particle/puff information (ordered to correspond to
                                          ! PPInfoReqs in Reqs).

End Type Results_

!-------------------------------------------------------------------------------------------------------------

Type :: ReqInfo_ ! Information on what requirements are required.

  ! ReqInfo_ is used in 4 ways to give 4 different types of information:
  ! (1) A list of requirements which don't depend on time due to the underlying quantity not depending on
  !     time.
  ! (2) A list of requirements which require contributions at discrete times, with the list giving all the
  !     requirements requiring a contribution at a given time.
  ! (3) A list of requirements which require contributions at discrete travel-times, with the list giving all
  !     the requirements requiring a contribution at a given travel-time.
  ! (4) A list of all requirements which require contributions at every time in some time interval, with the
  !     list giving all the requirements for which the time interval overlaps with a given time interval. Note
  !     that if a requirement has two or more time intervals which require contributions and which overlap
  !     with the given time interval, then the requirement is listed once for each such interval.

  Logical          :: NoReqs
  Type(ShortTime_) :: Time
  Integer          :: nFields
  Integer          :: iField        (MaxFieldReqsPerReqInfo)
  Integer          :: iT            (MaxFieldReqsPerReqInfo)
  Integer          :: iS            (MaxFieldReqsPerReqInfo)
  Type(ShortTime_) :: FieldStartTime(MaxFieldReqsPerReqInfo)
  Type(ShortTime_) :: FieldEndTime  (MaxFieldReqsPerReqInfo)
  Integer          :: nPPInfos
  Integer          :: iPPInfo (MaxPPInfoReqsPerReqInfo)
  Integer          :: iPPInfoT(MaxPPInfoReqsPerReqInfo)
  ! NoReqs         :: Indicates that there are no or no more requirements or times or travel-times for which
  !                   contributions are needed.
  ! Time           :: Depending on whether (1), (2), (3) or (4) above applies, this is
  !                     (1) undefined,
  !                     (2) time of contributions,
  !                     (3) travel-time of contributions,
  !                     (4) end of the given time interval.
  ! nFields        :} Number of fields, their indices, time indices (for type 2 and 4 use), travel-time
  ! iField         :} indices (for type 3 use), and start and end of the time interval for which contributions
  ! iT             :} are required (for type 4 use).
  ! iS             :}
  ! FieldStartTime :}
  ! FieldEndTime   :}
  ! nPPInfos       :] Number of sets of particle/puff information, their indices and time indices.
  ! iPPInfo        :]
  ! iPPInfoT       :]

End Type ReqInfo_

! Used for checking whether Timers need to be initialised.
Logical :: TimersInitialised = .false.
! Timers (see Timer module)
Type(Timer_), Save :: ProcessFieldsTimer,     &
                      OutputFieldsTimer


!-------------------------------------------------------------------------------------------------------------

Interface Operator(==) ! Equality of processing steps, field requirements and requirements for a set of
                       ! particle/puff information.
  Module Procedure ProcessEq
  Module Procedure FieldReqEq
  Module Procedure PPInfoReqEq
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface Operator(.Equiv.) ! Equivalence of field requirements (equivalence means the requirements are for
                            ! numerically equal fields, but output options may differ).
  Module Procedure FieldReqEquiv
End Interface

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

! Initialise timers

Subroutine OutputModuleTimerInitialise(TimerOpts)

  Implicit None

  Type(TimerOpts_) :: TimerOpts

  ! Create Timers
  If (.not. TimersInitialised) Then
    Call TimerCreate(ProcessFieldsTimer, "ProcessFields" ,TimerOpts)
    Call TimerCreate(OutputFieldsTimer,  "OutputFields"  ,TimerOpts)
    TimersInitialised = .true.
  End If

End Subroutine OutputModuleTimerInitialise

!-------------------------------------------------------------------------------------------------------------

! Output timer summary information

Subroutine OutputModuleTimerSummary()

  Implicit None

  If (TimersInitialised) Then
    Call TimerWriteSummary(ProcessFieldsTimer)
    Call TimerWriteSummary(OutputFieldsTimer)
  End If

End Subroutine OutputModuleTimerSummary

!-------------------------------------------------------------------------------------------------------------

Function PreInitOutputOpts() Result(OutputOpts)
! Pre-initialises the output options.

  Implicit None
  ! Function result:
  Type(OutputOpts_) :: OutputOpts ! Pre-initialised output options.

  OutputOpts%Initialised = .false.

End Function PreInitOutputOpts

!-------------------------------------------------------------------------------------------------------------

Subroutine InitOutputOpts(Folder, Seconds, DecimalPlaces, Pre65Format, OutputOpts)
! Initialises the output options.

  Implicit None
  ! Argument list:
  Character(*),      Intent(In)    :: Folder
  Logical,           Intent(In)    :: Seconds
  Integer,           Intent(In)    :: DecimalPlaces
  Logical,           Intent(In)    :: Pre65Format ! $$ remove this when not needed
  Type(OutputOpts_), Intent(InOut) :: OutputOpts
  ! Folder        :: Folder used for output (full or relative path). Relative paths should start with '.' or
  !                  should be blank to indicate the current folder. The path can optionally end in '\' or
  !                  '/'.
  ! Seconds       :: Indicates seconds are to be included in time output.
  ! DecimalPlaces :: Number of decimal places to be included for seconds in time output.
  ! OutputOpts    :: Initialised output options.

  If (OutputOpts%Initialised) Then
    Call Message('FATAL ERROR: The output options have been given more than once in the input file(s)', 3)
  End If

  Call TokenLengthTest(Folder, MaxFileNameLength, .false., 'Output Options', ' ', 'Folder')

  If (Folder == ' ') Then
    OutputOpts%Folder = '.\'
  Else If (                                               &
    Folder(Len_Trim(Folder):Len_Trim(Folder)) == '\' .or. &
    Folder(Len_Trim(Folder):Len_Trim(Folder)) == '/'      &
  ) Then
    OutputOpts%Folder = Trim(Folder)
  Else
    Call TokenLengthTest(Folder, MaxFileNameLength - 1, .false., 'Output Options', ' ', 'Folder')
    OutputOpts%Folder = Trim(Folder) // '\'
  End If

  OutputOpts%Seconds       = Seconds
  OutputOpts%DecimalPlaces = DecimalPlaces
  OutputOpts%Pre65Format   = Pre65Format

  OutputOpts%Initialised = .true.

End Subroutine InitOutputOpts

!-------------------------------------------------------------------------------------------------------------

Subroutine CheckOutputOpts(OutputOpts)
! Checks the output options.

  Implicit None
  ! Argument list:
  Type(OutputOpts_), Intent(In) :: OutputOpts ! Output options.

  If (.not.OutputOpts%Initialised) Then
    Call Message('FATAL ERROR: The output options have not been given in the input file(s)', 3)
  End If

End Subroutine CheckOutputOpts

!-------------------------------------------------------------------------------------------------------------

Function InitReqs(MaxFieldReqs, MaxFieldOutputGroups) Result (Reqs)
! Initialises a collection of requirements.

  Implicit None
  
  ! Argument list:
  Integer, Intent(In) :: MaxFieldReqs
  Integer, Intent(In) :: MaxFieldOutputGroups
  
  ! Function result:
  Type(Reqs_) :: Reqs ! Initialised collection of requirements.
  ! Locals:
  Integer :: i ! Loop index.

  Reqs%nProcesses   = 0

  Reqs%nFieldReqs   = 0
  Reqs%nFieldGroups = 0

  Reqs%nPdfs        = 0
  Reqs%nPdfGroups   = 0

  Reqs%nPPInfoReqs  = 0
  
  Reqs%MaxFieldReqs = MaxFieldReqs
  
  Reqs%MaxFieldOutputGroups = MaxFieldOutputGroups
  
  Allocate(Reqs%FieldReqs(MaxFieldReqs))
  Allocate(Reqs%FieldGroups(MaxFieldOutputGroups)) 

  ! Check QInfo table.
  Do i = 1, nQuantities
    If (                                                                    &
      Verify(QInfo(i),      ' 1234567890LEFOPQRA8:UMfpTSHZhzsrd') /= 0 .or. &
      QInfo(i)(1:2) == ' '                                             .or. &
      Verify(QInfo(i)(1:2), ' 1234567890'                   ) /= 0     .or. &
      Verify(QInfo(i)(4:4), 'LEFOPQRA'                      ) /= 0     .or. &
      (QInfo(i)(6:6) /= '1' .and. Scan(QInfo(i),      'p'   ) /= 0)    .or. &
      (QInfo(i)(6:6) /= ' ' .and. Scan(QInfo(i)(4:4), 'LEFOA') /= 0)   .or. &
      (QInfo(i)(6:6) == ' ' .and. Scan(QInfo(i)(4:4), 'PQR' ) /= 0)    .or. &
      Verify(QInfo(i)(6:6), ' 1234567890'                   ) /= 0     .or. &
      (Scan(QInfo(i), 'f'     ) /= 0 .and. Scan(QInfo(i), 'p') /= 0)   .or. &
      (Scan(QInfo(i), 'EFO'   ) /= 0 .and. Scan(QInfo(i), 'S') /= 0)   .or. &
      (Scan(QInfo(i), 'F'     ) /= 0 .and. Scan(QInfo(i), 'T') == 0)   .or. &
      (Scan(QInfo(i), 'F'     ) /= 0 .and. Scan(QInfo(i), 'H') == 0)   .or. &
      (Scan(QInfo(i), 'F'     ) /= 0 .and. Scan(QInfo(i), 'Z') == 0)   .or. &
      (Scan(QInfo(i), 'LEFPQR') /= 0 .and. Scan(QInfo(i), 'T') == 0)        &
    ) Then
      Call Message('UNEXPECTED FATAL ERROR in Init Reqs', 4)
    End If
  End Do

! $$ We can't check that P, Q and R quantities which depend on S depend on
! type L quantities which depend on S.

End Function InitReqs

!-------------------------------------------------------------------------------------------------------------

Function InitProcess(                                  &
           Name, Code, Ensemble,                       &
           TAvInt, nT,     T,    dT,     UseTGrid, T0, &
           HAvInt, nX, nY, X, Y, dX, dY, UseHGrid,     &
           ZAvInt, nZ,     Z,    dZ,     UseZGrid, BL, &
           DGridName,                                  &
           BlockKey, Item                              &
         )                                             &
Result (Process)
! Initialises a processing step.

! Note there are three types of error messages determined by BlockKey.
! BlockKey = Processing steps             : processing step assumed input in such a block
! BlockKey = Output Requirements - Fields : processing step assumed input in such a block
! BlockKey = blank                        : processing step assumed internally generated.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name
  Integer,      Intent(In) :: Code
  Logical,      Intent(In) :: Ensemble
  Integer,      Intent(In) :: TAvInt
  Integer,      Intent(In) :: nT
  Type(Time_),  Intent(In) :: T
  Type(Time_),  Intent(In) :: dT
  Logical,      Intent(In) :: UseTGrid
  Type(Time_),  Intent(In) :: T0
  Integer,      Intent(In) :: HAvInt
  Integer,      Intent(In) :: nX
  Integer,      Intent(In) :: nY
  Real(Std),    Intent(In) :: X
  Real(Std),    Intent(In) :: Y
  Real(Std),    Intent(In) :: dX
  Real(Std),    Intent(In) :: dY
  Logical,      Intent(In) :: UseHGrid
  Integer,      Intent(In) :: ZAvInt
  Integer,      Intent(In) :: nZ
  Real(Std),    Intent(In) :: Z
  Real(Std),    Intent(In) :: dZ
  Logical,      Intent(In) :: UseZGrid
  Logical,      Intent(In) :: BL
  Character(*), Intent(In) :: DGridName
  Character(*), Intent(In) :: BlockKey
  Character(*), Intent(In) :: Item
  ! Name      :: Name of processing step. Must be blank unless processing step input through a 'processing
  !              steps' block.
  ! Code      :: Code to indicate averaging/integrating, calculating exceedence probabilities, calculating
  !              percentiles, calculating max/min, calculating moments about zero, or calculating central
  !              moments.
  ! Ensemble  :: Indicates whether the processing is over the ensemble of cases.
  ! TAvInt    :: Code to indicate, for averaging/integrating, whether the time processing is averaging,
  !              integrating or nothing.
  ! nT        :} Number of times contributing to the processing, total duration contributing to the
  ! T         :} processing, and time interval between times contributing to the processing.
  ! dT        :}
  ! UseTGrid  :: Indicates T and dT are determined by the T-grid and nT, with T = the (finite) processing time
  !              implied by the T grid and dT = T / nT (as opposed to by the input values).
  ! T0        :: Time origin for contributions when time processing is over all T.
  ! HAvInt    :} Quantities analogous to TAvInt, nT, T, dT and UseTGrid for horizontal processing.
  ! nX        :}
  ! nY        :}
  ! X         :}
  ! Y         :}
  ! dX        :}
  ! dY        :}
  ! UseHGrid  :}
  ! ZAvInt    :] Quantities analogous to TAvInt, nT, T, dT and UseTGrid for vertical processing.
  ! nZ        :]
  ! Z         :]
  ! dZ        :]
  ! UseZGrid  :]
  ! BL        :: Indicates that the vertical processing is over the boundary layer.
  ! DGridName :: Name of D-Grid used for calculating exceedence probabilities, calculating percentiles,
  !              calculating max/min, calculating moments about zero, or calculating central moments. Blank if
  !              averaging integrating.
  ! BlockKey  :: Input block keyword. Blank if internally generated processing step.
  ! Item      :: Input item giving rise to this call. Blank if BlockKey = 'processing steps' or if internally
  !              generated processing step.
  !
  ! For time processing the possibilities are as follows:
  !                              TGrid?    TAvInt for Code =      nT     T         dT     UseTGrid     T0
  !                                       P_AvInt   Otherwise
  ! ----------------------------------------------------------------------------------------------------------
  ! Depending on T:
  ! T window defined by grid       Yes   P_Av/P_Int     0*       > 0     0*        0*        Yes    ref-time*
  ! Rolling T window               Yes   P_Av/P_Int     0*     { > 0  dT * nT     > 0 }      No     ref-time*
  !                                                            {  1      0         0  }
  ! T window from T = -infinity    Yes   P_Av/P_Int     0*        0*  infinity    > 0        No     ref-time*
  ! T window from T = -infinity    Yes   P_Av/P_Int     -         1   infinity  infinity     No     ref-time*
  ! T window over all T            No    P_Av/P_Int     0*        0*  infinity    > 0        No*       any
  ! T window over all T            No    P_Av/P_Int     -         1   infinity  infinity     No*    ref-time*
  !
  ! Not depending on T:            No       P_Av        0*        1      0         0         No     ref-time*
  !
  ! (For no T averaging use:     Yes/No     P_Av        -         1      0         0         No     ref-time*)
  ! ----------------------------------------------------------------------------------------------------------
  !
  ! For horizontal processing the possibilities are as follows:
  !                           HGrid?    HAvInt for Code =     nX/Y      X/Y        dX/Y    UseHGrid
  !                                    P_AvInt   Otherwise
  ! -----------------------------------------------------------------------------------------------
  ! Depending on H:
  ! H window defined by grid    Yes   P_Av/P_Int     0*       > 0        0*         0*        Yes
  ! Rolling H window            Yes   P_Av/P_Int     0*     { > 0   dX/Y * nX/Y    > 0 }      No
  !                                                         {  1         0          0  }
  ! H window over all H         No    P_Av/P_Int     -         1     infinity    infinity     No*
  !
  ! Not depending on H:         No       P_Av        0*        1         0          0         No
  !
  ! (For no H averaging use:  Yes/No     P_Av        -         1         0          0         No)
  ! -----------------------------------------------------------------------------------------------
  !
  ! For vertical processing the possibilities are as follows:
  !                              ZGrid?    ZAvInt for Code =      nZ      Z        dZ     UseZGrid  BL
  !                                       P_AvInt   Otherwise
  ! ---------------------------------------------------------------------------------------------------
  ! Depending on Z:
  ! Z window defined by grid       Yes   P_Av/P_Int     0*       > 0      0*       0*        Yes    No
  ! Rolling Z window               Yes   P_Av/P_Int     0*     { > 0   dZ * nZ    > 0 }      No     No
  !                                                            {  1       0        0  }
  ! Z window defined by b-layer    No    P_Av/P_Int     0*       > 0      0*       0*        No     Yes
  ! Z window over all Z            No    P_Av/P_Int     -         1   infinity  infinity     No*    No
  !
  ! Not depending on Z:            No       P_Av        0*        1       0        0         No     No
  !
  ! (For no Z averaging use:     Yes/No     P_Av        -         1       0        0         No     No)
  ! ---------------------------------------------------------------------------------------------------
  !
  ! Note (1) * implies value chosen is not physically meaningful but it is computationally convenient to set a
  !          particular value.
  !      (2) > 0 implies finite.
  !      (3) Infinite values of X/Y/Z/dX/dY/dZ are stored as -1.
  !      (4) TGrid?, HGrid?, ZGrid? refers to the presence of the grid in the field requirement (i.e. it
  !          refers to the output of the processing step).
  !      (5) 'Depending on T/H/Z' refers to the input to the processing step. This input might not depend on
  !          T/H/Z due to a previous processing step or due to the underlying quantity. In this case an
  !          averaging/integrating processing step has no effect, but other types of process still apply.
  !      (6) Some cases are not allowed.
  !            (a) nT/X/Y/Z = 1, T/X/Y/Z = dT/X/Y/Z = infinity is not allowed except for intrinsic
  !                averaging/integrating (the other possible uses - to give a fluctuations
  !                averaging/integrating region, to multiply by the integrating region, and to do nothing -
  !                are of no value here, although the first of these is still possible in combination with
  !                intrinsic averaging/integrating).
  !            (b) nT/X/Y/Z = 1, T/X/Y/Z = dT/X/Y/Z > 0 is not allowed except for averaging/integrating.
  !            (c) Averaging/integrating over infinite regions of space or time is not allowed for quantities
  !                of type F, O, Q or R. It can however be applied to quantities of other types used to
  !                derived the type F/O/Q/R quantity.
  !                $$ for example integrated deposition rate treated as an instantaneous quantity (total dep)
  !                for further processing.
  !      (7) X/Y/Z and dX/dY/dZ are defined in the H/ZGrids coord system.
  !      (8) For dT > 0, the times used are the nominal time - i * dT for i = 0, nT - 1.
  !      (9) For dX > 0, the x-values used are the nominal X + i * dX for i = -nX/2 + 1/2, nX/2 - 1/2.
  !          Similarly for Y/Z. If BL is true, the processing is over the boundary layer in the vertical,
  !          and the z-values used are i * H/nZ for i = 1/2, nZ - 1/2 where H is the boundary layer depth and
  !          the m agl coord system is used.
  ! See also comments at start of module.
  !
  ! Function result:
  Type (Process_) :: Process ! Initialised processing step.

  ! Name.
  If (BlockKey .CIEq. 'Processing Steps') Then
    Call TokenLengthTest(Name, MaxCharLength, .true., BlockKey, ' ', 'Name')
    If (Name(1:17) .CIEq. 'Internal Process ') Then
      Call Message(                                                         &
             'FATAL ERROR in reading item "'                             // &
             Trim(Name)                                                  // &
             '" from block "'                                            // &
             Trim(BlockKey)                                              // &
             '": this item name is of a form reserved for internal use',    &
             3                                                              &
           )
    End If
  Else
    If (Name /= ' ') Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)
  End If
  Process%Name = Name

  ! Code.
  If (                             &
    Code /= P_AvInt          .and. &
    Code /= P_Prob           .and. &
    Code /= P_Percent        .and. &
    Code /= P_MaxMin         .and. &
    Code /= P_Moments        .and. &
    Code /= P_CentralMoments       &
  ) Then
    Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)
  End If
  Process%Code = Code

  ! Ensemble.
  Process%Ensemble = Ensemble

  ! T.

  ! Check TAvInt a possible value.
  If (TAvInt /= 0 .and. TAvInt /= P_Av .and. TAvInt /= P_Int) Then
    Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)

  ! 'T window defined by grid' case.
  Else If (                                    &
    (Code == P_AvInt .neqv. TAvInt == 0) .and. &
    nT  > 0                              .and. &
    T  == ZeroTime()                     .and. &
    dT == ZeroTime()                     .and. &
    UseTGrid                             .and. &
    T0 == ReferenceTime()                      &
  ) Then

  ! 'Rolling T window' & 'not depending on T' cases.
  Else If (                                                           &
    (Code == P_AvInt .neqv. TAvInt == 0)                        .and. &
    nT  > 0                                                     .and. &
    T  >= ZeroTime() .and. .not.IsInfFuture(T)                  .and. &
    dT >= ZeroTime() .and. .not.IsInfFuture(dT)                 .and. &
    .not. (nT  > 1 .and. T == ZeroTime())                       .and. &
    .not. (nT == 1 .and. T /= ZeroTime() .and. Code /= P_AvInt) .and. &
    .not.UseTGrid                                               .and. &
    T0 == ReferenceTime()                                             &
  ) Then
    If (T /= dT * nT) Then
      If (BlockKey .CIEq. 'Processing Steps') Then
        Call Message(                             &
               'FATAL ERROR in reading item "' // &
               Trim(Name)                      // &
               '" from block "'                // &
               Trim(BlockKey)                  // &
               '": T is not a multiple of dT',    &
               3                                  &
             )
      Else If (BlockKey .CIEq. 'Output Requirements - Fields') Then
        If (Code == P_AvInt) Then
          Call Message(                                   &
                 'FATAL ERROR in reading item "'       // &
                 Trim(Item)                            // &
                 '" from block "'                      // &
                 Trim(BlockKey)                        // &
                 '": Av T is not a multiple of Av dT',    &
                 3                                        &
               )
        Else
          Call Message(                                 &
                 'FATAL ERROR in reading item "'     // &
                 Trim(Item)                          // &
                 '" from block "'                    // &
                 Trim(BlockKey)                      // &
                 '": P T is not a multiple of P dT',    &
                 3                                      &
               )
        End If
      Else
        Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)
      End If
    End If

  ! 'T window from T = -infinity' & 'T window over all T' cases with infinity > dT > 0.
  Else If (                                          &
    (Code == P_AvInt .neqv. TAvInt == 0)       .and. &
    nT == 0                                    .and. &
    IsInfFuture(T)                             .and. &
    dT > ZeroTime() .and. .not.IsInfFuture(dT) .and. &
    .not.UseTGrid                                    &
  ) Then

  ! 'T window from T = -infinity' & 'T window over all T' cases with dT = infinity.
  Else If (                                    &
    (Code == P_AvInt .neqv. TAvInt == 0) .and. &
    nT == 1                              .and. &
    IsInfFuture(T)                       .and. &
    IsInfFuture(dT)                      .and. &
    Code == P_AvInt                      .and. &
    .not.UseTGrid                        .and. &
    T0 == ReferenceTime()                      &
  ) Then

  ! Error.
  Else
    If (BlockKey .CIEq. 'Processing Steps') Then
      Call Message(                                                              &
             'FATAL ERROR in reading item "'                                  // &
             Trim(Name)                                                       // &
             '" from block "'                                                 // &
             Trim(BlockKey)                                                   // &
             '": The time processing options are not an allowed combination',    &
             3                                                                   &
           )
    Else If (BlockKey .CIEq. 'Output Requirements - Fields') Then
      Call Message(                                                              &
             'FATAL ERROR in reading item "'                                  // &
             Trim(Item)                                                       // &
             '" from block "'                                                 // &
             Trim(BlockKey)                                                   // &
             '": The time processing options are not an allowed combination',    &
             3                                                                   &
           )
    Else
      Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)
    End If
  End If

  Process%TAvInt   = TAvInt
  Process%nT       = nT
  Process%T        = T
  Process%dT       = dT
  Process%UseTGrid = UseTGrid
  Process%T0       = T0
  Process%sT       = Time2ShortTime(Process%T)
  Process%sdT      = Time2ShortTime(Process%dT)
  Process%sT0      = Time2ShortTime(Process%T0)

  ! X and Y.

  ! Check HAvInt a possible value.
  If (HAvInt /= 0 .and. HAvInt /= P_Av .and. HAvInt /= P_Int) Then
    Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)

  ! 'H window defined by grid' case.
  Else If (                                    &
    (Code == P_AvInt .neqv. HAvInt == 0) .and. &
    nX  > 0   .and. nY  > 0              .and. &
    X  == 0.0 .and. Y  == 0.0            .and. &
    dX == 0.0 .and. dY == 0.0            .and. &
    UseHGrid                                   &
  ) Then

  ! 'Rolling H window' & 'not depending on H' cases.
  Else If (                                                    &
    (Code == P_AvInt .neqv. HAvInt == 0)                 .and. &
    nX  > 0   .and. nY  > 0                              .and. &
    X  >= 0.0 .and. Y  >= 0.0                            .and. &
    dX >= 0.0 .and. dY >= 0.0                            .and. &
    .not. (nX  > 1 .and. X == 0.0)                       .and. &
    .not. (nY  > 1 .and. Y == 0.0)                       .and. &
    .not. (nX == 1 .and. X /= 0.0 .and. Code /= P_AvInt) .and. &
    .not. (nY == 1 .and. Y /= 0.0 .and. Code /= P_AvInt) .and. &
    .not.UseHGrid                                              &
  ) Then
    If (Abs(X - dX * nX) > X * 10.0 * Epsilon(1.0)) Then ! $$ could, for T,X,Y,Z, allow
                                                         ! non-exact T=nT*dT and adjust
      If (BlockKey .CIEq. 'Processing Steps') Then       ! T appropriately?
        Call Message(                             &
               'FATAL ERROR in reading item "' // &
               Trim(Name)                      // &
               '" from block "'                // &
               Trim(BlockKey)                  // &
               '": X is not a multiple of dX',    &
               3                                  &
             )
      Else If (BlockKey .CIEq. 'Output Requirements - Fields') Then
        If (Code == P_AvInt) Then
          Call Message(                                   &
                 'FATAL ERROR in reading item "'       // &
                 Trim(Item)                            // &
                 '" from block "'                      // &
                 Trim(BlockKey)                        // &
                 '": Av X is not a multiple of Av dX',    &
                 3                                        &
               )
        Else
          Call Message(                                 &
                 'FATAL ERROR in reading item "'     // &
                 Trim(Item)                          // &
                 '" from block "'                    // &
                 Trim(BlockKey)                      // &
                 '": P X is not a multiple of P dX',    &
                 3                                      &
               )
        End If
      Else
        Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)
      End If
    End If
    If (Abs(Y - dY * nY) > Y * 10.0 * Epsilon(1.0)) Then
      If (BlockKey .CIEq. 'Processing Steps') Then
        Call Message(                             &
               'FATAL ERROR in reading item "' // &
               Trim(Name)                      // &
               '" from block "'                // &
               Trim(BlockKey)                  // &
               '": Y is not a multiple of dY',    &
               3                                  &
             )
      Else If (BlockKey .CIEq. 'Output Requirements - Fields') Then
        If (Code == P_AvInt) Then
          Call Message(                                   &
                 'FATAL ERROR in reading item "'       // &
                 Trim(Item)                            // &
                 '" from block "'                      // &
                 Trim(BlockKey)                        // &
                 '": Av Y is not a multiple of Av dY',    &
                 3                                        &
               )
        Else
          Call Message(                                 &
                 'FATAL ERROR in reading item "'     // &
                 Trim(Item)                          // &
                 '" from block "'                    // &
                 Trim(BlockKey)                      // &
                 '": P Y is not a multiple of P dY',    &
                 3                                      &
               )
        End If
      Else
        Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)
      End If
    End If

  ! 'H window over all H' case.
  Else If (                                    &
    (Code == P_AvInt .neqv. HAvInt == 0) .and. &
    nX   == 1    .and. nY == 1           .and. &
    X    == -1.0 .and. Y  == -1.0        .and. &
    dX   == -1.0 .and. dY == -1.0        .and. &
    Code == P_AvInt                      .and. &
    .not.UseHGrid                              &
  ) Then

  ! Error.
  Else
    If (BlockKey .CIEq. 'Processing Steps') Then
      Call Message(                                                                    &
             'FATAL ERROR in reading item "'                                        // &
             Trim(Name)                                                             // &
             '" from block "'                                                       // &
             Trim(BlockKey)                                                         // &
             '": The horizontal processing options are not an allowed combination',    &
             3                                                                         &
           )
    Else If (BlockKey .CIEq. 'Output Requirements - Fields') Then
      Call Message(                                                                    &
             'FATAL ERROR in reading item "'                                        // &
             Trim(Item)                                                             // &
             '" from block "'                                                       // &
             Trim(BlockKey)                                                         // &
             '": The horizontal processing options are not an allowed combination',    &
             3                                                                         &
           )
    Else
      Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)
    End If
  End If

  Process%HAvInt   = HAvInt
  Process%nX       = nX
  Process%nY       = nY
  Process%X        = X
  Process%Y        = Y
  Process%dX       = dX
  Process%dY       = dY
  Process%UseHGrid = UseHGrid

  ! Z.

  ! Check ZAvInt a possible value.
  If (ZAvInt /= 0 .and. ZAvInt /= P_Av .and. ZAvInt /= P_Int) Then
    Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)

  ! 'Z window defined by grid' case.
  Else If (                                    &
    (Code == P_AvInt .neqv. ZAvInt == 0) .and. &
    nZ  > 0                              .and. &
    Z  == 0.0                            .and. &
    dZ == 0.0                            .and. &
    UseZGrid                             .and. &
    .not.BL                                    &
  ) Then

  ! 'Rolling Z window' & 'not depending on Z' cases.
  Else If (                                                    &
    (Code == P_AvInt .neqv. ZAvInt == 0)                 .and. &
    nZ  > 0                                              .and. &
    Z  >= 0.0                                            .and. &
    dZ >= 0.0                                            .and. &
    .not. (nZ  > 1 .and. Z == 0.0)                       .and. &
    .not. (nZ == 1 .and. Z /= 0.0 .and. Code /= P_AvInt) .and. &
    .not.UseZGrid                                        .and. &
    .not.BL                                                    &
  ) Then
    If (Abs(Z - dZ * nZ) > Z * 10.0 * Epsilon(1.0)) Then
      If (BlockKey .CIEq. 'Processing Steps') Then
        Call Message(                             &
               'FATAL ERROR in reading item "' // &
               Trim(Name)                      // &
               '" from block "'                // &
               Trim(BlockKey)                  // &
               '": Z is not a multiple of dZ',    &
               3                                  &
             )
      Else If (BlockKey .CIEq. 'Output Requirements - Fields') Then
        If (Code == P_AvInt) Then
          Call Message(                                   &
                 'FATAL ERROR in reading item "'       // &
                 Trim(Item)                            // &
                 '" from block "'                      // &
                 Trim(BlockKey)                        // &
                 '": Av Z is not a multiple of Av dZ',    &
                 3                                        &
               )
        Else
          Call Message(                                 &
                 'FATAL ERROR in reading item "'     // &
                 Trim(Item)                          // &
                 '" from block "'                    // &
                 Trim(BlockKey)                      // &
                 '": P Z is not a multiple of P dZ',    &
                 3                                      &
               )
        End If
      Else
        Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)
      End If
    End If

  ! 'Z window defined by b-layer' case.
  Else If (                                    &
    (Code == P_AvInt .neqv. ZAvInt == 0) .and. &
    nZ  > 0                              .and. &
    Z  == 0.0                            .and. &
    dZ == 0.0                            .and. &
    .not.UseZGrid                        .and. &
    BL                                         &
  ) Then

  ! 'Z window over all Z' case.
  Else If (                                    &
    (Code == P_AvInt .neqv. ZAvInt == 0) .and. &
    nZ   == 1                            .and. &
    Z    == -1.0                         .and. &
    dZ   == -1.0                         .and. &
    Code == P_AvInt                      .and. &
    .not.UseZGrid                        .and. &
    .not.BL                                    &
  ) Then

  ! Error.
  Else
    If (BlockKey .CIEq. 'Processing Steps') Then
      Call Message(                                                                  &
             'FATAL ERROR in reading item "'                                      // &
             Trim(Name)                                                           // &
             '" from block "'                                                     // &
             Trim(BlockKey)                                                       // &
             '": The vertical processing options are not an allowed combination',    &
             3                                                                       &
           )
    Else If (BlockKey .CIEq. 'Output Requirements - Fields') Then
      Call Message(                                                                  &
             'FATAL ERROR in reading item "'                                      // &
             Trim(Item)                                                           // &
             '" from block "'                                                     // &
             Trim(BlockKey)                                                       // &
             '": The vertical processing options are not an allowed combination',    &
             3                                                                       &
           )
    Else
      Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)
    End If
  End If

  Process%ZAvInt   = ZAvInt
  Process%nZ       = nZ
  Process%Z        = Z
  Process%dZ       = dZ
  Process%UseZGrid = UseZGrid
  Process%BL       = BL

  ! Percentile processing steps must always be over one point in space and time and must not be over the
  ! ensemble.
  If (Process%Code == P_Percent) Then
    If (                                             &
      (Process%Ensemble .neqv. .false.        ) .or. &
      (Process%nT         /=   1              ) .or. &
      (Process%T          /=   ZeroTime()     ) .or. &
      (Process%dT         /=   ZeroTime()     ) .or. &
      (Process%UseTGrid .neqv. .false.        ) .or. &
      (Process%T0         /=   ReferenceTime()) .or. &
      (Process%nX         /=   1              ) .or. &
      (Process%nY         /=   1              ) .or. &
      (Process%X          /=   0.0            ) .or. &
      (Process%Y          /=   0.0            ) .or. &
      (Process%dX         /=   0.0            ) .or. &
      (Process%dY         /=   0.0            ) .or. &
      (Process%UseHGrid .neqv. .false.        ) .or. &
      (Process%nZ         /=   1              ) .or. &
      (Process%Z          /=   0.0            ) .or. &
      (Process%dZ         /=   0.0            ) .or. &
      (Process%UseZGrid .neqv. .false.        ) .or. &
      (Process%BL       .neqv. .false.        )      &
    ) Then
      If (BlockKey .CIEq. 'Processing Steps') Then
        Call Message(                                                                          &
               'FATAL ERROR in reading item "'                                              // &
               Trim(Name)                                                                   // &
               '" from block "'                                                             // &
               Trim(BlockKey)                                                               // &
               '": Percentile processing steps must not involve processing over a spatial ' // &
               'or temporal extent or over the ensemble of cases.',                            &
               3                                                                               &
             )
      Else
        Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)
      End If
    End If
  End If

  ! DGridName.
  If (Code /= P_AvInt) Then
    If (BlockKey .CIEq. 'Processing Steps') Then
      Call TokenLengthTest(DGridName, MaxCharLength, .true., BlockKey, Name, 'D-Grid')
    Else
      Call TokenLengthTest(DGridName, MaxCharLength, .true., ' ', ' ', ' ')
    End If
  Else
    If (DGridName /= ' ') Then
      If (BlockKey .CIEq. 'Processing Steps') Then
        Call Message(                                                                           &
               'FATAL ERROR in reading item "'                                               // &
               Trim(Name)                                                                    // &
               '" from block "'                                                              // &
               Trim(BlockKey)                                                                // &
               '": The D-Grid name should be blank for an averaging or integrating process',    &
               3                                                                                &
             )
      Else
        Call Message('UNEXPECTED FATAL ERROR in InitProcess', 4)
      End If
    End If
  End If
  Process%DGridName = DGridName

End Function InitProcess

!-------------------------------------------------------------------------------------------------------------

Subroutine AddProcess(Process, Reqs, Name)
! Adds a processing step to a collection of requirements.

  Implicit None
  ! Argument list:
  Type(Process_),           Intent(In)             :: Process ! Processing step to be added.
  Type(Reqs_),              Intent(InOut)          :: Reqs    ! Collection of requirements.
  Character(MaxCharLength), Intent(Out),  Optional :: Name    ! If present, generates and returns a name for
                                                              ! the processing step.
  ! Locals:
  Type(Process_) :: ProcessL ! Local copy of processing step to be added.
  Integer        :: i        ! Loop index.

  If ((Present(Name) .and. Process%Name /= ' ') .or. (.not.Present(Name) .and. Process%Name == ' ')) Then
    Call Message('UNEXPECTED FATAL ERROR in AddProcess', 4)
  End If

  If (Present(Name)) Then
    Do i = 1, Reqs%nProcesses
      ProcessL      = Process
      ProcessL%Name = Reqs%Processes(i)%Name
      If (ProcessL == Reqs%Processes(i)) Then
        Name = Reqs%Processes(i)%Name
        Return
      End If
    End Do
  Else
    Do i = 1, Reqs%nProcesses
      If (Process%Name .CIEq. Reqs%Processes(i)%Name) Then
        If (Process == Reqs%Processes(i)) Then
          Return
        Else
          Call Message(                                                                 &
                 'FATAL ERROR in adding the processing step "'                       // &
                 Trim(Process%Name)                                                  // &
                 '": a different processing step with the same name already exists',    &
                 3                                                                      &
               )
        End If
      Endif
    End Do
  End If

  If (Reqs%nProcesses >= MaxProcesses) Then
    If (Present(Name)) Then
      Call Message(                                        &
             'FATAL ERROR in adding a processing step ' // &
             ': there are too many processing steps',      &
             3                                             &
           )
    Else
      Call Message(                                           &
             'FATAL ERROR in adding the processing step "' // &
             Trim(Process%Name)                            // &
             '": there are too many processing steps',        &
             3                                                &
           )
    End If
  End If

  Reqs%nProcesses                 = Reqs%nProcesses + 1
  Reqs%Processes(Reqs%nProcesses) = Process

  If (Present(Name)) Then
    Name = 'Internal Process ' // Int2Char(Reqs%nProcesses)
    Reqs%Processes(Reqs%nProcesses)%Name = Name
  End If

End Subroutine AddProcess

!-------------------------------------------------------------------------------------------------------------

Function FindProcessIndex(Name, Reqs)
! Finds the index of a processing step.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name ! Name of the processing step.
  Type(Reqs_),  Intent(In) :: Reqs ! Collection of requirements.
  ! Function result:
  Integer :: FindProcessIndex ! Index of the processing step.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Reqs%nProcesses
    If (Name .CIEq. Reqs%Processes(i)%Name) Then
      FindProcessIndex = i
      Return
    End If
  End Do

  Call Message(                              &
         'FATAL ERROR: processing step "' // &
         Trim(Name)                       // &
         '" not found',                      &
         3                                   &
       )

End Function FindProcessIndex

!-------------------------------------------------------------------------------------------------------------

Function ProcessEq(Process1, Process2)
! Tests for equality of processing steps.

  Implicit None
  ! Argument list:
  Type(Process_), Intent(In) :: Process1 !} The two processing steps.
  Type(Process_), Intent(In) :: Process2 !}
  ! Function result:
  Logical :: ProcessEq ! Indicates if processing steps are equal.

  ProcessEq = (Process1%Name      .CIEq. Process2%Name     ) .and. &
              (Process1%Code        ==   Process2%Code     ) .and. &
              (Process1%Ensemble  .eqv.  Process2%Ensemble ) .and. &
              (Process1%TAvInt      ==   Process2%TAvInt   ) .and. &
              (Process1%nT          ==   Process2%nT       ) .and. &
              (Process1%T           ==   Process2%T        ) .and. &
              (Process1%dT          ==   Process2%dT       ) .and. &
              (Process1%UseTGrid  .eqv.  Process2%UseTGrid ) .and. &
              (Process1%T0          ==   Process2%T0       ) .and. &
              (Process1%HAvInt      ==   Process2%HAvInt   ) .and. &
              (Process1%nX          ==   Process2%nX       ) .and. &
              (Process1%nY          ==   Process2%nY       ) .and. &
              (Process1%X           ==   Process2%X        ) .and. &
              (Process1%Y           ==   Process2%Y        ) .and. &
              (Process1%dX          ==   Process2%dX       ) .and. &
              (Process1%dY          ==   Process2%dY       ) .and. &
              (Process1%UseHGrid  .eqv.  Process2%UseHGrid ) .and. &
              (Process1%ZAvInt      ==   Process2%ZAvInt   ) .and. &
              (Process1%nZ          ==   Process2%nZ       ) .and. &
              (Process1%Z           ==   Process2%Z        ) .and. &
              (Process1%dZ          ==   Process2%dZ       ) .and. &
              (Process1%UseZGrid  .eqv.  Process2%UseZGrid ) .and. &
              (Process1%BL        .eqv.  Process2%BL       ) .and. &
              (Process1%DGridName .CIEq. Process2%DGridName)

End Function ProcessEq

!-------------------------------------------------------------------------------------------------------------

Function InitFieldReq(                                                             &
           Name,                                                                   &
           Quantity,                                                               &
           SpeciesName, SourceName, SourceGroupName, SizeDistName,                 &
           DecayDep, SemiInfApprox,                                                &
           TGridName, SGridName, HGridName, ZGridName, HCoordName, ZCoordName,     &
           ProcessNames,                                                           &
           IntrinsicTAv, IntrinsicHAv, IntrinsicZAv,                               &
           Fluctuations, FlAvT, FlAvX, FlAvY, FlAvZ, UseFlAvT, UseFlAvH, UseFlAvZ, &
           Sync, MaterialUnits, MaterialUnitName,                                  &
           BlockKey,                                                               &
           OutputRoute,                                                            &
           OutputGroup, Across, SeparateFile, OutputFormat,                        &
           XScale, YScale                                                          &
         )                                                                         &
Result (FieldReq)
! Initialises a field requirement.

! Note there are two types of error messages determined by BlockKey.
! BlockKey = Output Requirements - Fields : field requirement assumed input in such a block
! BlockKey = blank                        : field requirement assumed internally generated.

  Implicit None
  ! Argument list:
  Character(*),             Intent(In)           :: Name
  Character(*),             Intent(In)           :: Quantity
  Character(*),             Intent(In)           :: SpeciesName
  Character(*),             Intent(In)           :: SourceName
  Character(*),             Intent(In)           :: SourceGroupName
  Character(*),             Intent(In)           :: SizeDistName
  Logical,                  Intent(In)           :: DecayDep
  Logical,                  Intent(In)           :: SemiInfApprox
  Character(*),             Intent(In)           :: TGridName
  Character(*),             Intent(In)           :: SGridName
  Character(*),             Intent(In)           :: HGridName
  Character(*),             Intent(In)           :: ZGridName
  Character(*),             Intent(In)           :: HCoordName
  Character(*),             Intent(In)           :: ZCoordName
  Character(MaxCharLength), Intent(In)           :: ProcessNames(:) ! $$ changed from * to support Intel compiler on Cray
  Logical,                  Intent(In)           :: IntrinsicTAv
  Logical,                  Intent(In)           :: IntrinsicHAv
  Logical,                  Intent(In)           :: IntrinsicZAv
  Logical,                  Intent(In)           :: Fluctuations
  Type(ShortTime_),         Intent(In)           :: FlAvT
  Real(Std),                Intent(In)           :: FlAvX
  Real(Std),                Intent(In)           :: FlAvY
  Real(Std),                Intent(In)           :: FlAvZ
  Logical,                  Intent(In)           :: UseFlAvT
  Logical,                  Intent(In)           :: UseFlAvH
  Logical,                  Intent(In)           :: UseFlAvZ
  Logical,                  Intent(In)           :: Sync
  Type(MaterialUnits_),     Intent(InOut)        :: MaterialUnits
  Character(*),             Intent(In)           :: MaterialUnitName
  Character(*),             Intent(In)           :: BlockKey
  Character(*),             Intent(In), Optional :: OutputRoute
  Character(*),             Intent(In), Optional :: OutputGroup
  Character(*),             Intent(In), Optional :: Across
  Character(*),             Intent(In), Optional :: SeparateFile
  Character(*),             Intent(In), Optional :: OutputFormat
  Real(Std),                Intent(In), Optional :: XScale
  Real(Std),                Intent(In), Optional :: YScale
  ! Name            :: Name of field requirement.
  ! Quantity        :: Quantity required.
  ! SpeciesName     :: Name of species (blank if field not species dependent).
  ! SourceName      :} Results are due to the specified source or source-group only. Blank indicates no
  ! SourceGroupName :} restriction by source or source-group. Cannot both be set.
  ! SizeDistName    :: Name of particle size distribution giving the particle size ranges to be used in
  !                    calculating the field. Blank indicates no restriction by particle size.
  ! DecayDep        :: Indicates that the deposits made during an averaging or integrating time interval are
  !                    to be decayed to account for the time till the end of the interval. Must be false
  !                    unless Quantity is 'Deposition Rate', 'Dry Deposition Rate' or 'Wet Deposition Rate'.
  ! SemiInfApprox   :: Indicates that cloud gamma doses are to be calculated using the semi-infinite cloud
  !                    approximation. Must be false unless Quantity is 'Adult Effective Cloud Gamma Dose',
  !                    'Adult Lung Cloud Gamma Dose', 'Adult Thyroid Cloud Gamma Dose' or 'Adult Bone Surface
  !                    Cloud Gamma Dose'.
  ! TGridName       :] Name of time/travel-time/horizontal/vertical grids if used and blank otherwise.
  ! SGridName       :]
  ! HGridName       :]
  ! ZGridName       :]
  ! HCoordName      :} Name of horizontal/vertical coord systems if used, and blank otherwise. If a coord name
  ! ZCoordName      :} is blank, the field does not depend on the choice of coord system.
  ! ProcessNames    :: Names of the processing steps.
  ! IntrinsicTAv    :} Indicates intrinsic averaging/integrating over T/H/Z in first processing step. Must be
  ! IntrinsicHAv    :} false for type E, F or O quantities.
  ! IntrinsicZAv    :}
  ! Fluctuations    :: Indicates fluctuations are accounted for in the first probability or moment calculation
  !                    step and/or in the calculation of fluctuation parameters.

  ! FlAvT           :} Averaging scales for fluctuations (must be zero if UseFlAvT/H/Z is false).
  ! FlAvX           :}
  ! FlAvY           :}
  ! FlAvZ           :}
  ! UseFlAvT        :] Indicates FlAvT/X/Y/Z are the values to be used (must be false if Fluctuations is
  ! UseFlAvH        :] false). If false and Fluctuations true, the scales are determined by the processing
  ! UseFlAvZ        :] steps.
  ! Sync            :: Results are calculated when the particles/puffs are synchronised in time. Sync must be
  !                    true for type E, F or O quantities and false for deposition quantities of type L or if
  !                    IntrinsicTAv is true.
  ! MaterialUnitName :: Name of material unit for this request
  ! BlockKey        :: Input block keyword. Blank if internally generated field requirement.
  ! OutputRoute     :: The presence of various characters (case insensitively) have the following meanings:
  !                      D - numerical output to disk
  !                      S - numerical output to screen
  !                      G - graphical output to screen.
  ! OutputGroup     :: Name of output group for grouping numerical disk and screen output. Must be blank if no
  !                    numerical output.
  ! Across          :} The presence of T,S,X,Y,Z (case insensitively) implies the corresponding coord has the
  ! SeparateFile    :} 'Across' or 'Separate File' attribute. This controls the splitting of the t,s,x,y,z,d
  !                    dependence of output fields into rows, cols, and separate files. The 'across' settings
  !                    are relevant even if separate files are selected since they influence whether the
  !                    t,s,x,y,z,d appears to the left of a row or above a column. The presence of N (case
  !                    insensitively) in SeparateFile causes new output files to be used after a restart
  !                    (possibly overwriting any previous files). These settings can be used even if there is
  !                    no t,s,x,y,z,d dependence (and this may be necessary to meet restrictions on variations
  !                    between field requirements in the same output group).
  ! OutputFormat    :: The presence of various characters (case insensitively) have the following meanings:
  !                      I - include columns of grid Indices for the non-across coords
  !                      A - add spaces to Align columns
  !                      Z - include lines with all results Zero
  !                      F - Flushes buffer after each output time
  !                      2 - make output look more like Name II fields output.
  ! XScale          :} Scales used for plotting the results graphically. Values not significant if no
  ! YScale          :} graphical output.
  !
  ! Note the optional arguments must either all be present (for fields to be output) or all absent (for fields
  ! used only to calculate other fields).
  !
  ! Function result:
  Type(FieldReq_) :: FieldReq ! Initialised field requirement.
  ! Locals:
  Integer :: i ! Loop index.

  ! $$ Note error messages may be misleading when this is called from within the module
  ! rather than from input (they are worded assuming the error is an input error) - correct using BlockKey

  !--------------------------------!
  ! Quantities defining the field. !
  !--------------------------------!

  ! Name.
  If (BlockKey .CIEq. 'Output Requirements - Fields') Then
    Call TokenLengthTest(Name, MaxCharLength, .false., BlockKey, ' ', 'Name')
    If ((Name(1:18) .CIEQ. 'Unnamed Field Req ') .or. (Name(1:19) .CIEQ. 'Internal Field Req ')) Then
      Call Message(                                                         &
             'FATAL ERROR in reading item "'                             // &
             Trim(Name)                                                  // &
             '" from block "'                                            // &
             Trim(BlockKey)                                              // &
             '": this item name is of a form reserved for internal use',    &
             3                                                              &
           )
    End If
    If (Name == ' ') Then
      FieldReq%Name = 'Unnamed Field Req'
    Else
      FieldReq%Name = Name
    End If
  Else
    If (Name /= ' ') Call Message('UNEXPECTED FATAL ERROR in InitFieldReq', 4)
    FieldReq%Name = Name
  End If

  ! Quantity.
  Call TokenLengthTest(Quantity, MaxCharLength, .true., BlockKey, Name, 'Quantity')
  FieldReq%Quantity = Quantity

  ! iQuantity.
  FieldReq%iQuantity = 0
  Do i = 1, Size(Quantities)
    If (Quantity .CIEq. Quantities(i)) Then
      FieldReq%iQuantity = i
      Exit
    End If
  End Do
  If (FieldReq%iQuantity == 0) Then
    Call Message(                                                     &
           'FATAL ERROR in reading Output Requirements - Fields: ' // &
           'Quantity "'                                            // &
           Trim(Quantity)                                          // &
           '" not recognised',                                        &
           3                                                          &
         )
  End If

  ! LEFOPQRType.
  FieldReq%LEFOPQRType = QInfo(FieldReq%iQuantity)(4:4)

  ! TypeType.
  If      (Scan(QInfo(FieldReq%iQuantity), '8') /= 0) Then
    FieldReq%TypeType = '8'
  Else If (Scan(QInfo(FieldReq%iQuantity), ':') /= 0) Then
    FieldReq%TypeType = ':'
  Else
    FieldReq%TypeType = ' '
  End If

  ! nComponents.
  FieldReq%nComponents = Char2Int(QInfo(FieldReq%iQuantity)(1:2))

  ! Species.
  If (Scan(QInfo(FieldReq%iQuantity), 's') /= 0) Then
    If (SpeciesName == ' ') Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'Species should be specified for '                      // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  Else
    If (SpeciesName /= ' ') Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'Species should not be specified for '                  // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  End If
  Call TokenLengthTest(SpeciesName, MaxCharLength, .false., BlockKey, Name, 'Species')
  FieldReq%SpeciesName = SpeciesName

  ! Material Unit
  If ((FindMaterialUnitIndex(MaterialUnitName, MaterialUnits) .eq. -1) .and. &
      (Trim(MaterialUnitName) .ne. '' )) Then
    Call AddMaterialUnit(MaterialUnitName,'unknown unit',UnknownUnitType,1.0,MaterialUnits)
  End If

  If (Scan(QInfo(FieldReq%iQuantity), 's') == 0) Then
    If (MaterialUnitName /= ' ') Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'Material Unit should not be specified for '            // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  End If
  
  Call TokenLengthTest(MaterialUnitName, MaxCharLength, .false., BlockKey, Name, 'Material Unit')

  FieldReq%MaterialUnitName = MaterialUnitName
    
  ! Source and SourceGroup.
  Call TokenLengthTest(SourceName,      MaxCharLength, .false., BlockKey, Name, 'Source'      )
  Call TokenLengthTest(SourceGroupName, MaxCharLength, .false., BlockKey, Name, 'Source Group')
  If (Scan(QInfo(FieldReq%iQuantity), 'r') == 0) Then
    If (SourceName /= ' ' .or. SourceGroupName /= ' ') Then
      Call Message(                                                    &
             'FATAL ERROR in reading Output Requirements - Fields:' // &
             'values must not be given for Source or Source Group ' // &
             'for the quantity '                                    // &
             Trim(Quantity),                                           &
             3                                                         &
           )
    End If
  End If
  If (SourceName /= ' ' .and. SourceGroupName /= ' ') Then
    Call Message(                                                          &
           'FATAL ERROR in reading Output Requirements - Fields:'       // &
           'values must not be given for both Source and Source Group',    &
           3                                                               &
         )
  End If
  FieldReq%SourceName      = SourceName
  FieldReq%SourceGroupName = SourceGroupName

  ! SizeDist.
  Call TokenLengthTest(SizeDistName, MaxCharLength, .false., BlockKey, Name, 'Particle Size Distribution')
  If (Scan(QInfo(FieldReq%iQuantity), 'd') == 0) Then
    If (SizeDistName /= ' ') Then
      Call Message(                                                        &
             'FATAL ERROR in reading Output Requirements - Fields:'     // &
             'values must not be given for Particle Size Distribution ' // &
             'for the quantity '                                        // &
             Trim(Quantity),                                               &
             3                                                             &
           )
    End If
  End If
  FieldReq%SizeDistName = SizeDistName

  ! DecayDep.
  If (DecayDep) Then
    If (                                         &
      .not. FieldReq%iQuantity == Q_DryDep .and. &
      .not. FieldReq%iQuantity == Q_WetDep .and. &
      .not. FieldReq%iQuantity == Q_Dep          &
    ) Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             '"Decay deposition?" should be false for '              // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  End If
  FieldReq%DecayDep = DecayDep

  ! SemiInfApprox.
  If (SemiInfApprox) Then
    If (                                                       &
      .not. FieldReq%iQuantity == Q_AduEffCloudGammaDose .and. &
      .not. FieldReq%iQuantity == Q_AduLunCloudGammaDose .and. &
      .not. FieldReq%iQuantity == Q_AduThyCloudGammaDose .and. &
      .not. FieldReq%iQuantity == Q_AduBoSCloudGammaDose       &
    ) Then
      Call Message(                                                             &
             'FATAL ERROR in reading Output Requirements - Fields: '         // &
             'Cloud Gamma Dose must be specified for '                       // &
             '"Semi-infinite approx?" set to "Yes"',                            &
             3                                                                  &
           )
    End If
  End If
  FieldReq%SemiInfApprox = SemiInfApprox

  ! Grids and coords.
  If (Scan(QInfo(FieldReq%iQuantity), 'T') == 0) Then
    If (TGridName /= ' ') Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'T-Grid must not be specified for '                     // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  Else
    If (TGridName == ' ' .and. Size(ProcessNames) == 0) Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'T-Grid must be specified for field requirement "'      // &
             Trim(Name)                                              // &
             '"',                                                       &
             3                                                          &
           )
    End If
  End If
  If (Scan(QInfo(FieldReq%iQuantity), 'S') == 0) Then
    If (SGridName /= ' ') Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'S-Grid must not be specified for '                     // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  Else
    If (SGridName == ' ') Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'S-Grid must be specified for '                         // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  End If
  If (Scan(QInfo(FieldReq%iQuantity), 'H') == 0) Then
    If (HGridName /= ' ') Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'H-Grid must not be specified for '                     // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  Else
    If (HGridName == ' ' .and. Size(ProcessNames) == 0) Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'H-Grid must be specified for field requirement "'      // &
             Trim(Name)                                              // &
             '"',                                                       &
             3                                                          &
           )
    End If
  End If
  If (Scan(QInfo(FieldReq%iQuantity), 'Z') == 0) Then
    If (ZGridName /= ' ') Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'Z-Grid must not be specified for '                     // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  Else
    If (ZGridName == ' ' .and. Size(ProcessNames) == 0) Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'Z-Grid must be specified for field requirement "'      // &
             Trim(Name)                                              // &
             '"',                                                       &
             3                                                          &
           )
    End If
  End If
  If (Scan(QInfo(FieldReq%iQuantity), 'h') == 0) Then
    If (HCoordName /= ' ') Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'H-Coord must not be specified for '                    // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  Else
    If (HCoordName == ' ') Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'H-Coord must be specified for '                        // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  End If
  If (Scan(QInfo(FieldReq%iQuantity), 'z') == 0) Then
    If (ZCoordName /= ' ') Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'Z-Coord must not be specified for '                    // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  Else
    If (ZCoordName == ' ') Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'Z-Coord must be specified for '                        // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  End If
  Call TokenLengthTest(TGridName,  MaxCharLength, .false., BlockKey, ' ', 'T-Grid' )
  Call TokenLengthTest(SGridName,  MaxCharLength, .false., BlockKey, ' ', 'S-Grid' )
  Call TokenLengthTest(HGridName,  MaxCharLength, .false., BlockKey, ' ', 'H-Grid' )
  Call TokenLengthTest(ZGridName,  MaxCharLength, .false., BlockKey, ' ', 'Z-Grid' )
  Call TokenLengthTest(HCoordName, MaxCharLength, .false., BlockKey, ' ', 'H-Coord')
  Call TokenLengthTest(ZCoordName, MaxCharLength, .false., BlockKey, ' ', 'Z-Coord')
  FieldReq%TGridName  = TGridName
  FieldReq%SGridName  = SGridName
  FieldReq%HGridName  = HGridName
  FieldReq%ZGridName  = ZGridName
  FieldReq%HCoordName = HCoordName
  FieldReq%ZCoordName = ZCoordName

  ! Processing steps.
  If (Size(ProcessNames) > 0 .and. (FieldReq%TypeType /= ' ' .and. FieldReq%TypeType /= '8')) Then
  ! $$ Prevents processing of fields other than std or p64 precision - add support for other cases in future?
    Call Message(                                                                                     &
           'FATAL ERROR in reading Output Requirements - Fields: '                                 // &
           'no processing, averaging, probabilities, percentiles, max/min or moments can be used ' // &
           'with the quantity "'                                                                   // &
           Trim(FieldReq%Quantity)                                                                 // &
           '"',                                                                                       &
           3                                                                                          &
         )
  End If
  FieldReq%nProcesses                          = Size(ProcessNames)
  FieldReq%ProcessNames(1:FieldReq%nProcesses) = ProcessNames(:)

  ! Intrinsic averaging/integrating.
  If (                                                                                              &
    (IntrinsicTAv .or. IntrinsicHAv .or. IntrinsicZAv)                                              &
    .and.                                                                                           &
    (FieldReq%LEFOPQRType == 'E' .or. FieldReq%LEFOPQRType == 'F' .or. FieldReq%LEFOPQRType == 'O') &
  ) Then
    Call Message(                                                         &
           'FATAL ERROR in reading Output Requirements - Fields: '     // &
           'Intrinsic averaging/integrating can not be specified for ' // &
           Trim(Quantity),                                                &
           3                                                              &
         )
  End If
  If (IntrinsicTAv .and. Scan(QInfo(FieldReq%iQuantity), 'T') == 0) Then
    Call Message(                                                              &
           'FATAL ERROR in reading Output Requirements - Fields: '          // &
           'Intrinsic time averaging/integrating can not be specified for ' // &
           Trim(Quantity),                                                     &
           3                                                                   &
         )
  End If
  If (IntrinsicHAv .and. Scan(QInfo(FieldReq%iQuantity), 'H') == 0 .and. FieldReq%LEFOPQRType == 'L') Then
    Call Message(                                                                    &
           'FATAL ERROR in reading Output Requirements - Fields: '                // &
           'Intrinsic horizontal averaging/integrating can not be specified for ' // &
           Trim(Quantity),                                                           &
           3                                                                         &
         )
  End If ! $$ this can arise from type PQR and not be previously trapped (for Z too). Word message carefully
         ! for BlockKey = ' '
  If (IntrinsicZAv .and. Scan(QInfo(FieldReq%iQuantity), 'Z') == 0 .and. FieldReq%LEFOPQRType == 'L') Then
    Call Message(                                                                  &
           'FATAL ERROR in reading Output Requirements - Fields: '              // &
           'Intrinsic vertical averaging/integrating can not be specified for ' // &
           Trim(Quantity),                                                         &
           3                                                                       &
         )
  End If
  If (.not.IntrinsicTAv .and. Scan(QInfo(FieldReq%iQuantity), 'S') /= 0) Then
    Call Message(                                                                                          &
           'FATAL ERROR in reading Output Requirements - Fields: Intrinsic time averaging/integrating ' // &
           'must be specified for '                                                                     // &
           Trim(FieldReq%Quantity),                                                                        &
           3                                                                                               &
         )
  End If
  FieldReq%IntrinsicTAv = IntrinsicTAv
  FieldReq%IntrinsicHAv = IntrinsicHAv
  FieldReq%IntrinsicZAv = IntrinsicZAv

  ! Fluctuations.
  If (Fluctuations .and. Scan(QInfo(FieldReq%iQuantity), 'fp') == 0) Then
    Call Message(                                                     &
           'FATAL ERROR in reading Output Requirements - Fields: ' // &
           '"Fluctuations?" can not be specified for '             // &
           Trim(Quantity),                                            &
           3                                                          &
         )
  End If
  If (.not.Fluctuations .and. Scan(QInfo(FieldReq%iQuantity), 'p') /= 0) Then
    Call Message(                                                     &
           'FATAL ERROR in reading Output Requirements - Fields: ' // &
           '"Fluctuations?" must be specified for '                // &
           Trim(Quantity),                                            &
           3                                                          &
         )
  End If
  FieldReq%Fluctuations = Fluctuations

  ! FlAvT/X/Y/Z.
  If (.not.Fluctuations .and. (UseFlAvT .or. UseFlAvH .or. UseFlAvZ)) Then
    Call Message(                                                                      &
           'FATAL ERROR in reading Output Requirements - Fields: Fluctuation '      // &
           'averaging scales must not be given without the "Fluctuations?" switch',    &
           3                                                                           &
         )
  End If
  If (UseFlAvT) Then
    If (FlAvT < ZeroShortTime()) Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'Fluctuation averaging scales must be >= 0',               &
             3                                                          &
           )
    End If
  Else
    If (FlAvT /= ZeroShortTime()) Call Message('UNEXPECTED FATAL ERROR in InitFieldReq', 4)
  End If
  If (UseFlAvH) Then
    If (FlAvX < 0.0 .or. FlAvY < 0.0) Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'Fluctuation averaging scales must be >= 0',               &
             3                                                          &
           )
    End If
  Else
    If (FlAvX /= 0.0 .or. FlAvY /= 0.0) Call Message('UNEXPECTED FATAL ERROR in InitFieldReq', 4)
  End If
  If (UseFlAvZ) Then
    If (FlAvZ < 0.0) Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             'Fluctuation averaging scales must be >= 0',               &
             3                                                          &
           )
    End If
  Else
    If (FlAvZ /= 0.0) Call Message('UNEXPECTED FATAL ERROR in InitFieldReq', 4)
  End If
  FieldReq%FlAvT    = FlAvT
  FieldReq%FlAvX    = FlAvX
  FieldReq%FlAvY    = FlAvY
  FieldReq%FlAvZ    = FlAvZ
  FieldReq%UseFlAvT = UseFlAvT
  FieldReq%UseFlAvH = UseFlAvH
  FieldReq%UseFlAvZ = UseFlAvZ

  ! Sync.
  If (Sync) Then
    If (FieldReq%IntrinsicTAv) Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             '"Sync?" and "Intrinsic T Av?" cannot both be true',       &
             3                                                          &
           )
    End If
    If (                                  &
      FieldReq%iQuantity == Q_DryDep .or. &
      FieldReq%iQuantity == Q_WetDep .or. &
      FieldReq%iQuantity == Q_Dep         &
    ) Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             '"Sync?" must be false for '                            // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  Else
    If (FieldReq%LEFOPQRType == 'E' .or. FieldReq%LEFOPQRType == 'F' .or. FieldReq%LEFOPQRType == 'O') Then
      Call Message(                                                     &
             'FATAL ERROR in reading Output Requirements - Fields: ' // &
             '"Sync?" must be true for '                             // &
             Trim(Quantity),                                            &
             3                                                          &
           )
    End If
  End If
  FieldReq%Sync = Sync

  !----------------------------------------------!
  ! Quantities defining how the field is output. !
  !----------------------------------------------!

  ! All output options present.
  If (                          &
    Present(OutputRoute)  .and. &
    Present(OutputGroup)  .and. &
    Present(Across)       .and. &
    Present(SeparateFile) .and. &
    Present(OutputFormat) .and. &
    Present(XScale)       .and. &
    Present(YScale)             &
  ) Then

    If (Verify(OutputFormat, 'IiAaZzFf2 ') /= 0) Then
      Call Message(                                                                                   &
             'WARNING reading item "'                                                              // &
             Trim(Name)                                                                            // &
             '" from "Output Requirements - Fields": '                                             // &
             'the variable "Output Format" has an inappropriate character which will be ignored. ' // &
             'Allowed characters are I, i, A, a, Z, z, F, f and 2.',                                  &
             2                                                                                        &
           )
    End If   
   
    FieldReq%Disk   = Scan(OutputRoute, 'Dd') /= 0
    FieldReq%Screen = Scan(OutputRoute, 'Ss') /= 0

    If (OutputGroup == ' ' .and. (FieldReq%Disk .or. FieldReq%Screen)) Then
      Call Message(                                                                          &
             'FATAL ERROR: Output group cannot be blank for field output to screen or disk', &
             3                                                                               &
           )
    Else If (OutputGroup /= ' ' .and. .not.(FieldReq%Disk .or. FieldReq%Screen)) Then
      Call Message(                                                                            &
             'FATAL ERROR: Output group must be blank for field not output to screen or disk', &
             3                                                                                 &
           )
    Else
      FieldReq%OutputGroup = OutputGroup
    End If

    FieldReq%TAcross = Scan(Across, 'Tt') /= 0
    FieldReq%SAcross = Scan(Across, 'Ss') /= 0
    FieldReq%XAcross = Scan(Across, 'Xx') /= 0
    FieldReq%YAcross = Scan(Across, 'Yy') /= 0
    FieldReq%ZAcross = Scan(Across, 'Zz') /= 0
    FieldReq%DAcross = Scan(Across, 'Dd') /= 0

    FieldReq%TSeparateFile = Scan(SeparateFile, 'Tt') /= 0
    FieldReq%SSeparateFile = Scan(SeparateFile, 'Ss') /= 0
    FieldReq%XSeparateFile = Scan(SeparateFile, 'Xx') /= 0
    FieldReq%YSeparateFile = Scan(SeparateFile, 'Yy') /= 0
    FieldReq%ZSeparateFile = Scan(SeparateFile, 'Zz') /= 0
    FieldReq%DSeparateFile = Scan(SeparateFile, 'Dd') /= 0

    FieldReq%TNewFile = Scan(SeparateFile, 'Nn') /= 0

    FieldReq%GridIndices = Scan(OutputFormat, 'Ii') /= 0
    FieldReq%AlignCols   = Scan(OutputFormat, 'Aa') /= 0
    FieldReq%ZeroLines   = Scan(OutputFormat, 'Zz') /= 0
    FieldReq%Flush       = Scan(OutputFormat, 'Ff') /= 0
    FieldReq%NameII      = Scan(OutputFormat, '2' ) /= 0

    ! $$ NameII checks and warnings - if Name II flag set need to check other
    ! conditions needed for name II output are satisfied - possibly test elsewhere.
    ! Fatal error for:
    !   probs or percentiles other than 100th
    !   complex processing,
    !   YAcross false with TAcross false.
    ! Warnings for:
    !   multicase runs,
    !   source restricted output,
    !   T-grid with dt < 1 min.
    !   X or Y across with T across,
    !   Possibly also runs with two or more choices of averaging/integrating or 100th percentile
    !   calculation time
    !   which cannot be distinguished (but this is hard to test for).
    !   Possibly also options which don't match Name II field or time series output formats.
    ! See definition of NameII flag above for details.

    FieldReq%Graph = Scan(OutputRoute, 'Gg') /= 0

    If (FieldReq%Graph) Then
      FieldReq%XScale = XScale
      FieldReq%YScale = YScale
    Else
      FieldReq%XScale = 0.0
      FieldReq%YScale = 0.0
    End If

    ! Note graphical output only supports the cases
    !   (1) nZ > 1, nX = nY = nS = 0, FieldReq%nComponents = 1 (nT can be 1 or > 1)
    !   (2) nS > 1, nX = nY = nZ = 0, FieldReq%nComponents = 1 (nT can be 1 or > 1).
    ! Add checks for this somewhere (and notes on it) $$

    If (.not. FieldReq%Disk .and. .not. FieldReq%Screen .and. .not. FieldReq%Graph) Then
      Call Message('FATAL ERROR: No output route is specified', 3) ! reword $$
    End If

  ! No output options present.
  Else If (                          &
    .not.Present(OutputRoute)  .and. &
    .not.Present(OutputGroup)  .and. &
    .not.Present(Across)       .and. &
    .not.Present(SeparateFile) .and. &
    .not.Present(OutputFormat) .and. &
    .not.Present(XScale)       .and. &
    .not.Present(YScale)             &
  ) Then

    FieldReq%Disk   = .false.
    FieldReq%Screen = .false.

    FieldReq%OutputGroup = ' '

    FieldReq%TAcross = .false.
    FieldReq%SAcross = .false.
    FieldReq%XAcross = .false.
    FieldReq%YAcross = .false.
    FieldReq%ZAcross = .false.
    FieldReq%DAcross = .false.

    FieldReq%TSeparateFile = .false.
    FieldReq%SSeparateFile = .false.
    FieldReq%XSeparateFile = .false.
    FieldReq%YSeparateFile = .false.
    FieldReq%ZSeparateFile = .false.
    FieldReq%DSeparateFile = .false.

    FieldReq%TNewFile      = .false.

    FieldReq%GridIndices = .false.
    FieldReq%AlignCols   = .false.
    FieldReq%ZeroLines   = .false.
    FieldReq%Flush       = .false.
    FieldReq%NameII      = .false.

    FieldReq%Graph = .false.

    FieldReq%XScale = 0.0
    FieldReq%YScale = 0.0

  ! Some but not all output options present.
  Else

    Call Message('UNEXPECTED FATAL ERROR in InitFieldReq', 4)

  End If

  !-------------------------------------------------------!
  ! Indices of other fields needed to process this field. ! $$ needed?
  !-------------------------------------------------------!

  FieldReq%iProc        = 0
  FieldReq%iDryDep      = 0
  FieldReq%iWetDep      = 0
  FieldReq%iAirConc     = 0
  FieldReq%iEulConc     = 0
  FieldReq%iChemField   = 0
  FieldReq%iDensity     = 0
  FieldReq%iMeanS       = 0
  FieldReq%iXStats      = 0
  FieldReq%iSigmaC      = 0
  FieldReq%iMass        = 0
  FieldReq%iMeanZ       = 0
  FieldReq%iPhotonFlux  = 0

End Function InitFieldReq

!-------------------------------------------------------------------------------------------------------------

Subroutine AddFieldReq(FieldReq, Reqs, Name)
! Adds a field requirement to the collection of requirements.

  Implicit None
  ! Argument list:
  Type(FieldReq_),          Intent(In)             :: FieldReq ! Field requirement to be added.
  Type(Reqs_),              Intent(InOut)          :: Reqs     ! Collection of requirements.
  Character(MaxCharLength), Intent(Out),  Optional :: Name     ! If present, generates and returns a name for
                                                               ! the field requirement.
  ! Locals:
  Type(FieldReq_) :: FieldReqL ! Local copy of field requirement to be added.
  Integer         :: i         ! Loop index.

  If ((Present(Name) .and. FieldReq%Name /= ' ') .or. (.not.Present(Name) .and. FieldReq%Name == ' ')) Then
    Call Message('UNEXPECTED FATAL ERROR in AddProcess', 4)
  End If

  If (Present(Name)) Then
    Do i = 1, Reqs%nFieldReqs
      FieldReqL      = FieldReq
      FieldReqL%Name = Reqs%FieldReqs(i)%Name
      If (FieldReqL == Reqs%FieldReqs(i)) Then
        Name = Reqs%FieldReqs(i)%Name
        Return
      End If
    End Do
  Else
    Do i = 1, Reqs%nFieldReqs
      If (FieldReq%Name .CIEq. Reqs%FieldReqs(i)%Name) Then
        If (FieldReq == Reqs%FieldReqs(i)) Then
          Return
        Else
          Call Message(                                                                   &
                 'FATAL ERROR in adding the field requirement "'                       // &
                 Trim(FieldReq%Name)                                                   // &
                 '": a different field requirement with the same name already exists',    &
                 3                                                                        &
               )
        End If
      Endif
    End Do
  End If

  ! Check room for another field requirement.
  If (Reqs%nFieldReqs >= Reqs%MaxFieldReqs) Then
    If (Present(Name)) Then
      Call Message(                                                                                                        &
             'FATAL ERROR in adding a field requirement '                                                               // &
             ': there are too many field requirements. '                                                                // & 
             'The maximum number of field requirements is set to: '                                                     // &
              Trim(Int2Char(Reqs%MaxFieldReqs))                                                                         // &
             ' and can be altered via the variable "Max # Field Reqs" in the "Main Options" input block. '              // &
             'Note the number of field requirements may be greater than the number explicitly requested due to fields ' // &
             'required for e.g. time averaging or calculating derived quantities such as mixing ratios.',                  &
             3                                                                                                             &
           )
    Else
      Call Message(                                                                                                        &
             'FATAL ERROR in adding the field requirement "'                                                            // &
             Trim(FieldReq%Name)                                                                                        // &
             '": there are too many field requirements. '                                                               // &
             'The maximum number of field requirements is set to: '                                                     // &
              Trim(Int2Char(Reqs%MaxFieldReqs))                                                                         // &
             ' and can be altered via the variable "Max # Field Reqs" in the "Main Options" input block. '              // &
             'Note the number of field requirements may be greater than the number explicitly requested due to fields ' // &
             'required for e.g. time averaging or calculating derived quantities such as mixing ratios.',                  &            
             3                                                                                                             &
           )
    End If
  End If

  Reqs%nFieldReqs                 = Reqs%nFieldReqs + 1
  Reqs%FieldReqs(Reqs%nFieldReqs) = FieldReq

  If (FieldReq%Name .CIEq. 'Unnamed Field Req') Then
    Reqs%FieldReqs(Reqs%nFieldReqs)%Name = 'Unnamed Field Req ' // Int2Char(Reqs%nFieldReqs)
  End If

  If (Present(Name)) Then
    Name = 'Internal Field Req ' // Int2Char(Reqs%nFieldReqs)
    Reqs%FieldReqs(Reqs%nFieldReqs)%Name = Name
  End If

End Subroutine AddFieldReq

!-------------------------------------------------------------------------------------------------------------

Function FindFieldReqIndex(Name, Reqs, Error)
! Finds the index of a field requirement.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)            :: Name  ! Name of the field requirement.
  Type(Reqs_),  Intent(In)            :: Reqs  ! Collection of requirements.
  Logical,      Intent(Out), Optional :: Error ! Error flag for field requirement not found.
  ! Function result:
  Integer :: FindFieldReqIndex ! Index of the field requirement.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Reqs%nFieldReqs
    If (Name .CIEq. Reqs%FieldReqs(i)%Name) Then
      FindFieldReqIndex = i
      If (Present(Error)) Error = .false.
      Return
    End If
  End Do

  If (Present(Error)) Then
    FindFieldReqIndex = 0
    Error = .true.
  Else
    Call Message(                                &
           'FATAL ERROR: field requirement "' // &
           Trim(Name)                         // &
           '" not found',                        &
           3                                     &
         )
  End If

End Function FindFieldReqIndex

!-------------------------------------------------------------------------------------------------------------

Function FieldReqEq(FieldReq1, FieldReq2)
! Tests for equality of two field requirements.

  Implicit None
  ! Argument list:
  Type(FieldReq_), Intent(In) :: FieldReq1 !} The two field requirements.
  Type(FieldReq_), Intent(In) :: FieldReq2 !}
  ! Function result:
  Logical :: FieldReqEq ! Indicates if field requirements are equal.
  ! Locals:
  Integer :: i ! Loop index.

  FieldReqEq = (FieldReq1%Name            .CIEq. FieldReq2%Name            ) .and. &
               (FieldReq1%Quantity        .CIEq. FieldReq2%Quantity        ) .and. &
               (FieldReq1%SpeciesName     .CIEq. FieldReq2%SpeciesName     ) .and. &
               (FieldReq1%SourceName      .CIEq. FieldReq2%SourceName      ) .and. &
               (FieldReq1%SourceGroupName .CIEq. FieldReq2%SourceGroupName ) .and. &
               (FieldReq1%SizeDistName    .CIEq. FieldReq2%SizeDistName    ) .and. &
               (FieldReq1%DecayDep        .eqv.  FieldReq2%DecayDep        ) .and. &
               (FieldReq1%SemiInfApprox   .eqv.  FieldReq2%SemiInfApprox   ) .and. &
               (FieldReq1%TGridName       .CIEq. FieldReq2%TGridName       ) .and. &
               (FieldReq1%SGridName       .CIEq. FieldReq2%SGridName       ) .and. &
               (FieldReq1%HGridName       .CIEq. FieldReq2%HGridName       ) .and. &
               (FieldReq1%ZGridName       .CIEq. FieldReq2%ZGridName       ) .and. &
               (FieldReq1%HCoordName      .CIEq. FieldReq2%HCoordName      ) .and. &
               (FieldReq1%ZCoordName      .CIEq. FieldReq2%ZCoordName      ) .and. &
               (FieldReq1%nProcesses        ==   FieldReq2%nProcesses      ) .and. &
               (FieldReq1%IntrinsicTAv    .eqv.  FieldReq2%IntrinsicTAv    ) .and. &
               (FieldReq1%IntrinsicHAv    .eqv.  FieldReq2%IntrinsicHAv    ) .and. &
               (FieldReq1%IntrinsicZAv    .eqv.  FieldReq2%IntrinsicZAv    ) .and. &
               (FieldReq1%Fluctuations    .eqv.  FieldReq2%Fluctuations    ) .and. &
               (FieldReq1%FlAvT             ==   FieldReq2%FlAvT           ) .and. &
               (FieldReq1%FlAvX             ==   FieldReq2%FlAvX           ) .and. &
               (FieldReq1%FlAvY             ==   FieldReq2%FlAvY           ) .and. &
               (FieldReq1%FlAvZ             ==   FieldReq2%FlAvZ           ) .and. &
               (FieldReq1%UseFlAvT        .eqv.  FieldReq2%UseFlAvT        ) .and. &
               (FieldReq1%UseFlAvH        .eqv.  FieldReq2%UseFlAvH        ) .and. &
               (FieldReq1%UseFlAvZ        .eqv.  FieldReq2%UseFlAvZ        ) .and. &
               (FieldReq1%Sync            .eqv.  FieldReq2%Sync            ) .and. &
               (FieldReq1%Disk            .eqv.  FieldReq2%Disk            ) .and. &
               (FieldReq1%Screen          .eqv.  FieldReq2%Screen          ) .and. &
               (FieldReq1%OutputGroup     .CIEq. FieldReq2%OutputGroup     ) .and. &
               (FieldReq1%TAcross         .eqv.  FieldReq2%TAcross         ) .and. &
               (FieldReq1%SAcross         .eqv.  FieldReq2%SAcross         ) .and. &
               (FieldReq1%XAcross         .eqv.  FieldReq2%XAcross         ) .and. &
               (FieldReq1%YAcross         .eqv.  FieldReq2%YAcross         ) .and. &
               (FieldReq1%ZAcross         .eqv.  FieldReq2%ZAcross         ) .and. &
               (FieldReq1%DAcross         .eqv.  FieldReq2%DAcross         ) .and. &
               (FieldReq1%TSeparateFile   .eqv.  FieldReq2%TSeparateFile   ) .and. &
               (FieldReq1%SSeparateFile   .eqv.  FieldReq2%SSeparateFile   ) .and. &
               (FieldReq1%XSeparateFile   .eqv.  FieldReq2%XSeparateFile   ) .and. &
               (FieldReq1%YSeparateFile   .eqv.  FieldReq2%YSeparateFile   ) .and. &
               (FieldReq1%ZSeparateFile   .eqv.  FieldReq2%ZSeparateFile   ) .and. &
               (FieldReq1%DSeparateFile   .eqv.  FieldReq2%DSeparateFile   ) .and. &
               (FieldReq1%TNewFile        .eqv.  FieldReq2%TNewFile        ) .and. &
               (FieldReq1%GridIndices     .eqv.  FieldReq2%GridIndices     ) .and. &
               (FieldReq1%AlignCols       .eqv.  FieldReq2%AlignCols       ) .and. &
               (FieldReq1%ZeroLines       .eqv.  FieldReq2%ZeroLines       ) .and. &
               (FieldReq1%Flush           .eqv.  FieldReq2%Flush           ) .and. &
               (FieldReq1%NameII          .eqv.  FieldReq2%NameII          ) .and. &
               (FieldReq1%Graph           .eqv.  FieldReq2%Graph           ) .and. &
               (FieldReq1%XScale            ==   FieldReq2%XScale          ) .and. &
               (FieldReq1%YScale            ==   FieldReq2%YScale          ) .and. &
               (FieldReq1%MaterialUnitName  ==   FieldReq2%MaterialUnitName)

  Do i = 1, FieldReq1%nProcesses
    FieldReqEq = FieldReqEq .and. (FieldReq1%ProcessNames(i) .CIEq. FieldReq2%ProcessNames(i))
  End Do

End Function FieldReqEq

!-------------------------------------------------------------------------------------------------------------

Function FindEquivFieldReqIndex(FieldReq, Reqs, Error)
! Finds the index of an equivalent field requirement. Equivalence means the requirements are for numerically
! equal fields.

! Note, in using this routine, Error should be present if there might not be an equivalent field requirement.

  Implicit None
  ! Argument list:
  Type(FieldReq_), Intent(In)            :: FieldReq ! Field requirement.
  Type(Reqs_),     Intent(In)            :: Reqs     ! Collection of requirements.
  Logical,         Intent(Out), Optional :: Error    ! Error flag for field requirement not found.
  ! Function result:
  Integer :: FindEquivFieldReqIndex ! Index of the field requirement.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Reqs%nFieldReqs
    If (FieldReq .Equiv. Reqs%FieldReqs(i)) Then
      FindEquivFieldReqIndex = i
      If (Present(Error)) Error = .false.
      Return
    End If
  End Do

  If (Present(Error)) Then
    FindEquivFieldReqIndex = 0
    Error = .true.
  Else
    Call Message('UNEXPECTED FATAL ERROR in FindEquivFieldReqIndex: field requirement not found', 4)
  End If

End Function FindEquivFieldReqIndex

!-------------------------------------------------------------------------------------------------------------

Function FieldReqEquiv(FieldReq1, FieldReq2)
! Tests for equivalence of two field requirements. Equivalence means the requirements are for numerically
! equal fields.
!
! $$ If two field requirements only differ in their material unit they are treated as different.
! $$ This is potentially wasteful as the two corresponding fields will often only differ by a
! $$ simple rescaling factor.

  Implicit None
  ! Argument list:
  Type(FieldReq_), Intent(In) :: FieldReq1 !} The two field requirements.
  Type(FieldReq_), Intent(In) :: FieldReq2 !}
  ! Function result:
  Logical :: FieldReqEquiv ! Indicates if field requirements are equal.
  ! Locals:
  Integer :: i ! Loop index.

  FieldReqEquiv = (FieldReq1%Quantity        .CIEq. FieldReq2%Quantity        ) .and. &
                  (FieldReq1%SpeciesName     .CIEq. FieldReq2%SpeciesName     ) .and. &
                  (FieldReq1%SourceName      .CIEq. FieldReq2%SourceName      ) .and. &
                  (FieldReq1%SourceGroupName .CIEq. FieldReq2%SourceGroupName ) .and. &
                  (FieldReq1%SizeDistName    .CIEq. FieldReq2%SizeDistName    ) .and. &
                  (FieldReq1%DecayDep        .eqv.  FieldReq2%DecayDep        ) .and. &
                  (FieldReq1%SemiInfApprox   .eqv.  FieldReq2%SemiInfApprox   ) .and. &
                  (FieldReq1%TGridName       .CIEq. FieldReq2%TGridName       ) .and. &
                  (FieldReq1%SGridName       .CIEq. FieldReq2%SGridName       ) .and. &
                  (FieldReq1%HGridName       .CIEq. FieldReq2%HGridName       ) .and. &
                  (FieldReq1%ZGridName       .CIEq. FieldReq2%ZGridName       ) .and. &
                  (FieldReq1%HCoordName      .CIEq. FieldReq2%HCoordName      ) .and. &
                  (FieldReq1%ZCoordName      .CIEq. FieldReq2%ZCoordName      ) .and. &
                  (FieldReq1%nProcesses        ==   FieldReq2%nProcesses      ) .and. &
                  (FieldReq1%IntrinsicTAv    .eqv.  FieldReq2%IntrinsicTAv    ) .and. &
                  (FieldReq1%IntrinsicHAv    .eqv.  FieldReq2%IntrinsicHAv    ) .and. &
                  (FieldReq1%IntrinsicZAv    .eqv.  FieldReq2%IntrinsicZAv    ) .and. &
                  (FieldReq1%Fluctuations    .eqv.  FieldReq2%Fluctuations    ) .and. &
                  (FieldReq1%FlAvT             ==   FieldReq2%FlAvT           ) .and. &
                  (FieldReq1%FlAvX             ==   FieldReq2%FlAvX           ) .and. &
                  (FieldReq1%FlAvY             ==   FieldReq2%FlAvY           ) .and. &
                  (FieldReq1%FlAvZ             ==   FieldReq2%FlAvZ           ) .and. &
                  (FieldReq1%UseFlAvT        .eqv.  FieldReq2%UseFlAvT        ) .and. &
                  (FieldReq1%UseFlAvH        .eqv.  FieldReq2%UseFlAvH        ) .and. &
                  (FieldReq1%UseFlAvZ        .eqv.  FieldReq2%UseFlAvZ        ) .and. &
                  (FieldReq1%Sync            .eqv.  FieldReq2%Sync            ) .and. &
                  (FieldReq1%MaterialUnitName  ==   FieldReq2%MaterialUnitName)

  Do i = 1, FieldReq1%nProcesses
    FieldReqEquiv = FieldReqEquiv .and. (FieldReq1%ProcessNames(i) .CIEq. FieldReq2%ProcessNames(i))
  End Do

End Function FieldReqEquiv

!-------------------------------------------------------------------------------------------------------------

Function InitPdfReq(                                       &
           Name, SpeciesName, SourceName, SourceGroupName, &
           PdfType,                                        &
           HGridName, ZGridName, TGridName,                &
           HCoordName, ZCoordName,                         &
           CalcType, TAvOrInt, AvTime,                     &
           PdfMode, PdfSize, PdfThresholds,                &
           Screen, Disk, AvEnsemble                        &
         )
! Initialises a pdf requirement.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name        ! Name of pdf.
  Character(*), Intent(In) :: SpeciesName ! Name of species.
  Character(*), Intent(In) :: SourceName  !} Results are due to the specified source
  Character(*), Intent(In) :: SourceGroupName   !} or source-group only. Blank indicates no
                                          !} restriction by source or source-group.
  Integer,      Intent(In) :: PdfType
                       ! Type of pdf
                       ! (1 = exceedence probabilities, 2 = percentiles).
  Character(*), Intent(In) :: HGridName
  Character(*), Intent(In) :: ZGridName
  Character(*), Intent(In) :: TGridName
                       !} Name of horizontal/vertical/temporal grids if used, blank
                       !} otherwise. If blank, either the pdf has no variation
                       !} in that variable (e.g. continuous-release plume
                       !} calculations have no temporal variation) or
                       !} it refers to a time-integrated calculation.
  Character(*), Intent(In) :: HCoordName
  Character(*), Intent(In) :: ZCoordName
  Integer,      Intent(In) :: CalcType
                       ! Type of fluctuations calculation to be carried out:
                       !  1 - continuous releases (with possibility of time-averaging),
                       !  2 - time-integrated quantities for finite-duration releases,
                       !  3 - instantaneous quantities for finite-duration releases.
  Integer,      Intent(In) :: TAvOrInt
                       ! Indicates that time-averaged concentrations are considered
                       ! (valid for type 1 calculation only). TAvOrInt should be false
                       ! if CalcType is not equal to 1.
  Real(Std),    Intent(In) :: AvTime
                       ! Averaging time used for the time-averaged concentrations
                       ! when TAvOrInt is true (valid for type 1 calculation only).
                       ! The averaging time is given in seconds. AvTime should
                       ! be zero if TAvOrInt is false.
  Integer,      Intent(In) :: PdfMode     ! Mode of specifying pdf threshold values
                                          ! (1 = auto-mode, 2 = user-specified).
  Integer,      Intent(In) :: PdfSize     ! Number of user-specified threshold values
                                          ! to be used in the calculation of the pdf.
                                          ! PdfSize should be set to -1 in auto-mode.
                                          ! $$ replace with zero?
  Real(Std),    Intent(In) :: PdfThresholds(MaxPdfSize)
                                          ! Array of user-specified threshold values
                                          ! to be used in the calculation of the pdf.
                                          !$$ What about auto-mode here ???
                                          ! $$ -1 is coded as a special value - replace with 0?
  Logical,      Intent(In) :: Screen      ! Indicates results are output numerically
                                          ! to the screen.
  Logical,      Intent(In) :: Disk        ! Indicates results are output numerically
                                          ! to disk.
  Logical,      Intent(In) :: AvEnsemble  ! Indicates results are statistically
                                          ! processed.
  ! Function result:
  Type(PdfReq_) :: InitPdfReq ! Initialised pdf requirement.
  ! Locals:
  Type(PdfReq_) :: PdfReq ! Local copy of function result.
  Integer       :: i      ! Loop index.

  ! Check inputs OK $$

  ! iQuantity.
  PdfReq%iQuantity = 0
  Do i = 1, Size(Quantities)
    If (Name .CIEq. Quantities(i)) Then
      PdfReq%iQuantity = i
      Exit
    End If
  End Do
  If (PdfReq%iQuantity == 0) Then
    Call Message(                                                   &
           'ERROR: Quantity "' // Trim(Name) // '" not recognised', &
           3                                                        &
         )
  End If

  If (Len_Trim(Name) > MaxCharLength .or. Len_Trim(Name) == 0) Then
    Call Message('Error in InitPdfReq: invalid name', 3)
  End If
  PdfReq%Name         = Name

  If (PdfType < 1 .or. PdfType > 2) Then
    Call Message('Error in InitPdfReq: invalid pdf type', 3)
  End If
  PdfReq%PdfType  = PdfType

  If (Len_Trim(SpeciesName) > MaxCharLength .or. Len_Trim(SpeciesName) == 0) Then
    Call Message('Error in InitPdfReq: invalid species name', 3)
  End If
  PdfReq%SpeciesName = SpeciesName

  If (                                        &
    Len_Trim(SourceName) > MaxCharLength .or. &
    Len_Trim(SourceGroupName) > MaxCharLength &
  ) Then
    Call Message('Error in InitPdfReq: invalid source/group name', 3)
  End If
  PdfReq%SourceName      = SourceName
  PdfReq%SourceGroupName = SourceGroupName

  If (Len_Trim(SourceName) > 0) Then
    If (Len_Trim(SourceGroupName) > 0) Then
      ! Case not permitted.
      Call Message('Error in InitPdfReq: source name and group name both specified', 3)
    Else
      ! Single source specified.
      PdfReq%nSources = 1
    End If
  Else
    If (Len_Trim(SourceGroupName) > 0) Then
      ! Source group specified.
      PdfReq%nSources = 1 !$$ More work needed on this !
    Else
      ! All sources considered.
      PdfReq%nSources = 0
    End If
  End If

  If (Len_Trim(HGridName) > MaxCharLength .or. Len_Trim(HGridName) == 0) Then
    Call Message('Error in InitPdfReq: invalid name for horizontal grid', 3)
  End If
  PdfReq%HGridName = HGridName

  If (Len_Trim(ZGridName) > MaxCharLength) Then
    Call Message('Error in InitPdfReq: invalid name for vertical grid', 3)
  End If
  PdfReq%ZGridName = ZGridName

  If (Len_Trim(TGridName) > MaxCharLength) Then
    Call Message('Error in InitPdfReq: invalid name for temporal grid', 3)
  End If
  PdfReq%TGridName = TGridName

  PdfReq%HCoordName = HCoordName ! $$ some checks may be needed here.
  PdfReq%ZCoordName = ZCoordName

  If (CalcType < 1 .or. CalcType > 3) Then
    Call Message('Error in InitPdfReq: invalid calculation type', 3)
  End If
  PdfReq%CalcType    = CalcType

  If (TAvOrInt == P_Av .and. CalcType /= 1) Then
    Call Message('Error in InitPdfReq: time-averaging option is not valid here', 3)
  End If
  If (TAvOrInt == P_Av) Then
    ! Averaging time should be greater than zero but less than limit of MaxAvTime.
    If (AvTime <= 0.0) Then
      Call Message('Error in InitPdfReq: averaging time should be positive', 3)
    Else If (AvTime > MaxAvTime) Then
      Call Message(                                                             &
             'Error in InitPdfReq: averaging time exceeds maximum limit of ' // &
             Trim(Std2Char(MaxAvTime)),                                         &
             3                                                                  &
           )
    End If
  Else
    ! Averaging time should be zero.
    If (AvTime /= 0.0) Then
      Call Message('Error in InitPdfReq: time-averaging option not requested', 3)
    End If
  End If
  PdfReq%TAvOrInt = TAvOrInt
  PdfReq%AvTime   = AvTime

  If (PdfMode == 1) Then
    ! Auto-mode for pdf threshold values.
    PdfReq%PdfSize       = -1
    PdfReq%PdfThresholds = -1 !$$ Need to consider this case ?
  Else If (PdfMode == 2) Then
    ! User-specified threshold values for the pdf.
    PdfReq%PdfSize       = PdfSize
    PdfReq%PdfThresholds = PdfThresholds
  Else
    Call Message('Error in InitPdfReq: invalid value of pdf mode', 3)
  End If
  PdfReq%PdfMode      = PdfMode

  PdfReq%Screen       = Screen
  PdfReq%Disk         = Disk
  PdfReq%AvEnsemble   = AvEnsemble
  PdfReq%SeparateFile = ' '
  PdfReq%OutputGroup  = ' '

  InitPdfReq = PdfReq

End Function InitPdfReq

!-------------------------------------------------------------------------------------------------------------

Subroutine AddPdfReq(PdfReq, Reqs)
! Adds a pdf requirement to the collection of requirements.

  Implicit None
  ! Argument list:
  Type(PdfReq_), Intent(In)    :: PdfReq ! Pdf requirement.
  Type(Reqs_),   Intent(InOut) :: Reqs   ! Description of requirements.
  ! Locals:
  Integer         :: i             ! Index.
  Integer         :: iG            ! Index of output group.
  Type(FieldReq_) :: ExtraFieldReq !

  ! Check that there is room for another pdf requirement.
  If (Reqs%nPdfs == MaxPdfReqs) Then
    Call Message('Error in AddPdfReq: no space to store new requirement', 3)
  End If

  ! Check that the request is for a concentration pdf.
  If (.not.(PdfReq%Name .CIEq. 'Air Concentration')) Then
    Call Message('Error in AddPdfReq: pdfs other than for concentration are not supported', 3)
  End If

  ! Calculate output group index of pdf requirement.
  ! [Pdfs with same grids are assigned to the same output group.]
  iG = 0
  Do i = 1, Reqs%nPdfs
    If (PdfReq%HGridName .CIEq. Reqs%PdfReqs(i)%HGridName) Then
    If (PdfReq%ZGridName .CIEq. Reqs%PdfReqs(i)%ZGridName) Then
    If (PdfReq%TGridName .CIEq. Reqs%PdfReqs(i)%TGridName) Then
      iG = Reqs%PdfReqs(i)%iOutputGroup
      Exit
    End If
    End If
    End If
  End Do

  ! Add request for pdf to Reqs.
  Reqs%nPdfs               = Reqs%nPdfs + 1
  Reqs%PdfReqs(Reqs%nPdfs) = PdfReq

  ! Update number of output groups and assign output group to latest pdf requirement.
  If (iG == 0) Then
    iG              = Reqs%nPdfGroups + 1
    Reqs%nPdfGroups = iG
  End If
  Reqs%PdfReqs(Reqs%nPdfs)%iOutputGroup = iG

End Subroutine AddPdfReq

!-------------------------------------------------------------------------------------------------------------

Function InitPPInfoReq(                                                        &
           Name,                                                               &
           Particles, Puffs, FirstParticle, LastParticle, FirstPuff, LastPuff, &
           SourceName,                                                         &
           Met, Mass, PlumeRise, DispersionScheme, PuffFamily,                 &
           HCoordName, ZCoordName, TGridName,                                  &
           Sync,                                                               &
           OutputRoute, OutputFormat                                           &
         )                                                                     &
Result (PPInfoReq)
! Initialises a requirement for a set of particle/puff information.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name             ! Name of requirement for a set of particle/puff information.
  Character(1), Intent(In) :: Particles        !} A, S and N indicate that All, Some and No particles/puffs
  Character(1), Intent(In) :: Puffs            !} are included.
  Integer,      Intent(In) :: FirstParticle    !] Output restricted to particle/puff numbers in this range.
  Integer,      Intent(In) :: LastParticle     !] Not significant unless Particles/Puffs = 'S'.
  Integer,      Intent(In) :: FirstPuff        !]
  Integer,      Intent(In) :: LastPuff         !]
  Character(*), Intent(In) :: SourceName       ! Results are due to the specified source only. Blank indicates
                                               ! no restriction by source.
  Logical,      Intent(In) :: Met              !} Indicates information is included on met, mass, plume rise,
  Logical,      Intent(In) :: Mass             !} the dispersion scheme, and (for puffs only) the puff family.
  Logical,      Intent(In) :: PlumeRise        !}
  Logical,      Intent(In) :: DispersionScheme !}
  Logical,      Intent(In) :: PuffFamily       !}
  Character(*), Intent(In) :: HCoordName       !] Name of horizontal/vertical coord systems used.
  Character(*), Intent(In) :: ZCoordName       !]
  Character(*), Intent(In) :: TGridName        ! Name of time grid if used, blank otherwise. If blank, results
                                               ! are output every time step.
  Logical,      Intent(In) :: Sync             ! Results are calculated when the particles/puffs are
                                               ! synchronised in time.
  Character(*), Intent(In) :: OutputRoute      ! The presence of various characters (case insensitively) have
                                               ! the following meanings:
                                               !   D - numerical output to disk
                                               !   S - numerical output to screen.
  Character(*), Intent(In) :: OutputFormat     ! The presence of various characters (case insensitively) have
                                               ! the following meanings:
                                               !   P - separate Particles/Puffs in separate files
                                               !   T - separate Times in separate files
                                               !   F - Flushes buffer after each output time.
  ! Function result:
  Type(PPInfoReq_) :: PPInfoReq ! Initialised requirement for a set of particle/puff information.

  ! Name.
  Call TokenLengthTest(                                             &
         Name, MaxCharLength, .true.,                               &
         'Output Requirements - Sets of Particle/Puff Information', &
         ' ', 'Name'                                                &
       )
  PPInfoReq%Name = Name

  ! Particles, Puffs.
  If (Scan(Particles, 'ASN') == 0 .or. Scan(Puffs, 'ASN') == 0) Then
    Call Message('UNEXPECTED FATAL ERROR in InitPPInfoReq', 4)
  End If
  PPInfoReq%Particles = Particles
  PPInfoReq%Puffs     = Puffs

  ! FirstParticle, LastParticle, FirstPuff, LastPuff.
  If (Particles == 'S') Then
    If (FirstParticle < 1 .or. LastParticle < FirstParticle) Then
      Call Message(                                                                                   &
             'FATAL ERROR reading item "'                                                          // &
             Trim(Name)                                                                            // &
             '" from "Output Requirements - Sets of Particle/Puff Information": '                  // &
             'the range of particles requested must satisfy 1 <= first particle <= last particle',    &
             3                                                                                        &
           )
    End If
    PPInfoReq%FirstParticle = FirstParticle
    PPInfoReq%LastParticle  = LastParticle
  Else
    PPInfoReq%FirstParticle = 0
    PPInfoReq%LastParticle  = 0
  End If
  If (Puffs == 'S') Then
    If (FirstPuff < 1 .or. LastPuff < FirstPuff) Then
      Call Message(                                                                       &
             'FATAL ERROR reading item "'                                              // &
             Trim(Name)                                                                // &
             '" from "Output Requirements - Sets of Particle/Puff Information": '      // &
             'the range of puffs requested must satisfy 1 <= first puff <= last puff',    &
             3                                                                            &
           )
    End If
    PPInfoReq%FirstPuff = FirstPuff
    PPInfoReq%LastPuff  = LastPuff
  Else
    PPInfoReq%FirstPuff = 0
    PPInfoReq%LastPuff  = 0
  End If

  ! SourceName.
  Call TokenLengthTest(                                             &
         SourceName, MaxCharLength, .false.,                        &
         'Output Requirements - Sets of Particle/Puff Information', &
         Name, 'Source'                                             &
       )
  PPInfoReq%SourceName = SourceName

  ! Optional parts of output.
  PPInfoReq%Met              = Met
  PPInfoReq%Mass             = Mass
  PPInfoReq%PlumeRise        = PlumeRise
  PPInfoReq%DispersionScheme = DispersionScheme
  PPInfoReq%PuffFamily       = PuffFamily

  ! Coords.
  Call TokenLengthTest(                                             &
         HCoordName, MaxCharLength, .true.,                         &
         'Output Requirements - Sets of Particle/Puff Information', &
         Name, 'H-Coord'                                            &
       )
  PPInfoReq%HCoordName = HCoordName
  Call TokenLengthTest(                                             &
         ZCoordName, MaxCharLength, .true.,                         &
         'Output Requirements - Sets of Particle/Puff Information', &
         Name, 'Z-Coord'                                            &
       )
  PPInfoReq%ZCoordName = ZCoordName

  ! T-Grid.
  Call TokenLengthTest(                                             &
         TGridName, MaxCharLength, .false.,                         &
         'Output Requirements - Sets of Particle/Puff Information', &
         Name, 'T-Grid'                                             &
       )
  PPInfoReq%TGridName = TGridName

  ! Sync.
  PPInfoReq%Sync = Sync

  ! OutputRoute.
  If (Verify(OutputRoute, 'DdSsGg ') /= 0) Then ! $$ similar tests for other char string inputs
                                                ! - e.g. in fieldreq
    Call Message(                                                                  &
           'FATAL ERROR reading item "'                                         // &
           Trim(Name)                                                           // &
           '" from "Output Requirements - Sets of Particle/Puff Information": ' // &
           'the variable "Output Route" has an inappropriate value',               &
           3                                                                       &
         )
  End If
  PPInfoReq%Disk   = Scan(OutputRoute, 'Dd') /= 0
  PPInfoReq%Screen = Scan(OutputRoute, 'Ss') /= 0
  PPInfoReq%Graph  = Scan(OutputRoute, 'Gg') /= 0
  If (.not.PPInfoReq%Disk .and. .not.PPInfoReq%Screen .and. .not.PPInfoReq%Graph) Then
    Call Message(                                                                     &
            ! $$ similar test for fieldreqs?
           'FATAL ERROR reading item "'                                            // &
           Trim(Name)                                                              // &
           '" from "Output Requirements - Sets of Particle/Puff Information": '    // &
           'the value of "Output Route" means the requirement will not be output',    &
           3                                                                          &
         )
  End If

  ! OutputFormat.
  If (Verify(OutputFormat, 'TtPpFf ') /= 0) Then
    Call Message(                                                                                   &
           'WARNING reading item "'                                                              // &
           Trim(Name)                                                                            // &
           '" from "Output Requirements - Sets of Particle/Puff Information": '                  // &
           'the variable "Output Format" has an inappropriate character which will be ignored. ' // &
           'Allowed characters are T, t, P, p, F and f.',                                           &
           2                                                                                        &
         )
  End If
  PPInfoReq%TSeparateFile = Scan(OutputFormat, 'Tt') /= 0
  PPInfoReq%PSeparateFile = Scan(OutputFormat, 'Pp') /= 0
  PPInfoReq%Flush         = Scan(OutputFormat, 'Ff') /= 0
  If (PPInfoReq%TSeparateFile .and. TGridName == ' ') Then
    Call Message('FATAL ERROR in InitPPInfoReq', 3)
  End If
  If (PPInfoReq%PSeparateFile .and. PPInfoReq%TSeparateFile) Then
    Call Message('FATAL ERROR in InitPPInfoReq', 3)
  End If

End Function InitPPInfoReq

!-------------------------------------------------------------------------------------------------------------

Subroutine AddPPInfoReq(PPInfoReq, Reqs)
! Adds a requirement for a set of particle/puff information to the collection of requirements.

  Implicit None
  ! Argument list:
  Type(PPInfoReq_), Intent(In)    :: PPInfoReq ! Requirement for a set of particle/puff information.
  Type(Reqs_),      Intent(InOut) :: Reqs      ! Collection of requirements.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Reqs%nPPInfoReqs
    If (PPInfoReq%Name .CIEq. Reqs%PPInfoReqs(i)%Name) Then
      Call Message(                                                                              &
             'FATAL ERROR in adding the requirement for a set of particle/puff information "' // &
             Trim(PPInfoReq%Name)                                                             // &
             '": a requirement with the same name already exists',                               &
             3                                                                                   &
           )
    Endif
  End Do

  ! Check room for another requirement for a set of particle/puff information.
  If (Reqs%nPPInfoReqs >= MaxPPInfoReqs) Then
    Call Message(                                                                              &
           'FATAL ERROR in adding the requirement for a set of particle/puff information "' // &
           Trim(PPInfoReq%Name)                                                             // &
           '": there are too many such requirements',                                          &
           3                                                                                   &
         )
  End If

  ! Add requirement for a set of particle/puff information to Reqs.
  Reqs%nPPInfoReqs                  = Reqs%nPPInfoReqs + 1
  Reqs%PPInfoReqs(Reqs%nPPInfoReqs) = PPInfoReq

End Subroutine AddPPInfoReq

!-------------------------------------------------------------------------------------------------------------

Function PPInfoReqEq(PPInfoReq1, PPInfoReq2)
! Tests for equality of two requirement for a set of particle/puff information.

  Implicit None
  ! Argument list:
  Type(PPInfoReq_), Intent(In) :: PPInfoReq1 !} The two requirement for a set of particle/puff information.
  Type(PPInfoReq_), Intent(In) :: PPInfoReq2 !}
  ! Function result:
  Logical :: PPInfoReqEq ! Indicates if requirements for a set of particle/puff information are equal.

  PPInfoReqEq = (PPInfoReq1%Name             .CIEq. PPInfoReq2%Name            ) .and. &
                (PPInfoReq1%Particles          ==   PPInfoReq2%Particles       ) .and. &
                (PPInfoReq1%Puffs              ==   PPInfoReq2%Puffs           ) .and. &
                (PPInfoReq1%FirstParticle      ==   PPInfoReq2%FirstParticle   ) .and. &
                (PPInfoReq1%LastParticle       ==   PPInfoReq2%LastParticle    ) .and. &
                (PPInfoReq1%FirstPuff          ==   PPInfoReq2%FirstPuff       ) .and. &
                (PPInfoReq1%LastPuff           ==   PPInfoReq2%LastPuff        ) .and. &
                (PPInfoReq1%SourceName       .CIEq. PPInfoReq2%SourceName      ) .and. &
                (PPInfoReq1%Met              .eqv.  PPInfoReq2%Met             ) .and. &
                (PPInfoReq1%Mass             .eqv.  PPInfoReq2%Mass            ) .and. &
                (PPInfoReq1%PlumeRise        .eqv.  PPInfoReq2%PlumeRise       ) .and. &
                (PPInfoReq1%DispersionScheme .eqv.  PPInfoReq2%DispersionScheme) .and. &
                (PPInfoReq1%PuffFamily       .eqv.  PPInfoReq2%PuffFamily      ) .and. &
                (PPInfoReq1%HCoordName       .CIEq. PPInfoReq2%HCoordName      ) .and. &
                (PPInfoReq1%ZCoordName       .CIEq. PPInfoReq2%ZCoordName      ) .and. &
                (PPInfoReq1%TGridName        .CIEq. PPInfoReq2%TGridName       ) .and. &
                (PPInfoReq1%Sync             .eqv.  PPInfoReq2%Sync            ) .and. &
                (PPInfoReq1%Disk             .eqv.  PPInfoReq2%Disk            ) .and. &
                (PPInfoReq1%Screen           .eqv.  PPInfoReq2%Screen          ) .and. &
                (PPInfoReq1%TSeparateFile    .eqv.  PPInfoReq2%TSeparateFile   ) .and. &
                (PPInfoReq1%PSeparateFile    .eqv.  PPInfoReq2%PSeparateFile   ) .and. &
                (PPInfoReq1%Flush            .eqv.  PPInfoReq2%Flush           )

End Function PPInfoReqEq

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpReqs(Sources, Coords, Grids, Specieses, Reqs, MaterialUnits)
! Sets up collection of requirements, including adding any coord systems or grids required.

  Implicit None
  ! Argument list:
  Type(Sources_),        Intent(In),    Target :: Sources        ! Collection of sources.
  Type(Coords_),         Intent(InOut)         :: Coords         ! Collection of coord systems.
  Type(Grids_),          Intent(InOut), Target :: Grids          ! Collection of grids.
  Type(Specieses_),      Intent(InOut), Target :: Specieses      ! Collection of species.
  Type(Reqs_),           Intent(InOut), Target :: Reqs           ! Collection of requirements.
  Type(MaterialUnits_),  Intent(InOut)         :: MaterialUnits  ! Collection of material units.
  ! Locals:
  Type(TGrid_),            Pointer :: TGrid
  Type(DGrid_),            Pointer :: DGrid
  Type(Species_),          Pointer :: Species
  Type(Source_),           Pointer :: Source
  Type(Process_),          Pointer :: Process
  Type(Process_),          Pointer :: Process1
  Type(FieldReq_),         Pointer :: FieldReq
  Type(FieldReq_),         Pointer :: FieldReq1
  Integer                          :: iTGrid
  Integer                          :: iDGrids(MaxDGrids)
  Integer                          :: iProcesses(MaxProcessesPerFieldReq)
  Type(Process_)                   :: ExtraProcess
  Type(FieldReq_)                  :: ExtraFieldReq
  Character(MaxCharLength)         :: TGridName
  Character(MaxCharLength)         :: ProcessName
  Character(MaxCharLength)         :: FieldReqName
  Logical                          :: IntrinsicAvRegionsOverlap
  Logical                          :: Changes
  Logical                          :: Error
  Integer                          :: i
  Integer                          :: j
  Integer                          :: k
  Real(Std)                        :: Temp_1(2)
  Real(Std)                        :: Temp_2(2)
  Type(ShortTime_)                 :: T0 ! $$ remove in due course
  Integer                          :: iSpecies
  Type(MaterialUnit_)              :: SpeciesMaterialUnit
  Character(MaxCharLength)         :: SpeciesMaterialUnitName
  Type(MaterialUnit_)              :: FieldReqMaterialUnit
  ! TGrid                     :} Abbreviations for grids, processing steps and field requirements.
  ! DGrid                     :}
  ! Process                   :}
  ! Process1                  :}
  ! FieldReq                  :}
  ! FieldReq1                 :}
  ! iTGrid                    :] Indices of grids and processes
  ! iDGrids                   :]
  ! iProcesses                :]
  ! ExtraProcess              :} Extra processing steps and field requirements to be added.
  ! ExtraFieldReq             :}
  ! TGridName                 :] Names of grids, processing steps and field requirements.
  ! ProcessName               :]
  ! FieldReqName              :]
  ! IntrinsicAvRegionsOverlap :: Indicates only intrinsic averaging/integrating to do but intrinsic
  !                              averaging/integrating regions overlap.
  ! Changes                   :: Indicates changes have occurred and so the process should be iterated.
  ! Error                     :: Error flag.
  ! i                         :} Loop indices.
  ! j                         :}
  ! k                         :}
  ! Temp_1                    :] Temporary arrays needed to avoid array-copying warnings on Linux Intel
  ! Temp_2                    :] compiler with -C option. $$ review whether still needed.
  ! iSpecies                  : Index of species associated with a particular request
  ! SpeciesMaterialUnit       : Material unit of species associated with a field request
  ! SpeciesMaterialUnitName   : Name of material unit of species associated with a field request
  ! FieldReqMaterialUnit      : Material unit of a particular field request

  ! Pdfs.
  i = 1

  Do While (i <= Reqs%nPdfs)

    ! Concentration pdf requires Mean Concentration and Sigma C.
    If (Reqs%PdfReqs(i)%iQuantity == Q_AirConc) Then
      iSpecies = FindSpeciesIndex(Reqs%PdfReqs(i)%SpeciesName, Specieses)
      SpeciesMaterialUnitName = Specieses%Specieses(iSpecies)%MaterialUnitName
      ExtraFieldReq = InitFieldReq(                                        &
                        Name            = ' ',                             &
                        Quantity        = Reqs%PdfReqs(i)%Name,            &
                        SpeciesName     = Reqs%PdfReqs(i)%SpeciesName,     &
                        SourceName      = Reqs%PdfReqs(i)%SourceName,      &
                        SourceGroupName = Reqs%PdfReqs(i)%SourceGroupName, &
                        SizeDistName    = ' ',                             &
                        DecayDep        = .false.,                         &
                        SemiInfApprox   = .false.,                         &
                        TGridName       = Reqs%PdfReqs(i)%TGridName,       &
                        SGridName       = ' ',                             &
                        HGridName       = Reqs%PdfReqs(i)%HGridName,       &
                        ZGridName       = Reqs%PdfReqs(i)%ZGridName,       &
                        HCoordName      = ' ',                             &
                        ZCoordName      = ' ',                             &
                  !      TAvOrInt        = Reqs%PdfReqs(i)%TAvOrInt,        &
                  !      AvTime          = ZeroTime(),                      &
                  !      AvdT            = ZeroTime(),                      &
                        ProcessNames    = Reqs%FieldReqs(1)%ProcessNames(1:0), &
                        IntrinsicTAv    = .false.,                         &
                        IntrinsicHAv    = .false.,                         &
                        IntrinsicZAv    = .false.,                         &
                        Fluctuations    = .false.,                         & ! $$ (below too)
                        FlAvT           = ZeroShortTime(), & ! $$
                        FlAvX           = 0.0,             & ! $$
                        FlAvY           = 0.0,             & ! $$
                        FlAvZ           = 0.0,             & ! $$
                        UseFlAvT        = .false.,         & ! $$
                        UseFlAvH        = .false.,         & ! $$
                        UseFlAvZ        = .false.,         & ! $$
                        Sync            = .false.,                         &
                        MaterialUnits   = MaterialUnits,                   &
                        MaterialUnitName = SpeciesMaterialUnitName,        &
                        BlockKey        = ' '                              &
                      )
      ! $$ meaning of PdfReq%TAvOrInt need a bit more thought - probably
      ! reflects averaging of inputs or averaging over met changes,
      ! not averaging of fluctuations (below too). FieldReq%TAvOrInt too.
      ! Currently Av time for fluctuations set to zero. - need to set up processnames
      Reqs%PdfReqs(i)%iAirConc = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        Reqs%PdfReqs(i)%iAirConc = FindEquivFieldReqIndex(ExtraFieldReq, Reqs)
      End If

      ExtraFieldReq = InitFieldReq(                                        &
                        Name            = ' ',                             &
                        Quantity        = 'Sigma C',                       &
                        SpeciesName     = Reqs%PdfReqs(i)%SpeciesName,     &
                        SourceName      = Reqs%PdfReqs(i)%SourceName,      &
                        SourceGroupName = Reqs%PdfReqs(i)%SourceGroupName, &
                        SizeDistName    = ' ',                             &
                        DecayDep        = .false.,                         &
                        SemiInfApprox   = .false.,                         &
                        TGridName       = Reqs%PdfReqs(i)%TGridName,       &
                        SGridName       = ' ',                             &
                        HGridName       = Reqs%PdfReqs(i)%HGridName,       &
                        ZGridName       = Reqs%PdfReqs(i)%ZGridName,       &
                        HCoordName      = ' ',                             &
                        ZCoordName      = ' ',                             &
                   !     TAvOrInt        = Reqs%PdfReqs(i)%TAvOrInt,        &
                   !     AvTime          = ZeroTime(),                      &
                   !     AvdT            = ZeroTime(),                      &
                        ProcessNames    = Reqs%FieldReqs(1)%ProcessNames(1:0), &
                        IntrinsicTAv    = .false.,                         &
                        IntrinsicHAv    = .false.,                         &
                        IntrinsicZAv    = .false.,                         &
                        Fluctuations    = .true.,                          &
                        FlAvT           = ZeroShortTime(), & ! $$
                        FlAvX           = 0.0,             & ! $$
                        FlAvY           = 0.0,             & ! $$
                        FlAvZ           = 0.0,             & ! $$
                        UseFlAvT        = .false.,         & ! $$
                        UseFlAvH        = .false.,         & ! $$
                        UseFlAvZ        = .false.,         & ! $$
                        Sync            = .false.,                         &
                        MaterialUnits   = MaterialUnits,                   &
                        MaterialUnitName = SpeciesMaterialUnitName,        &
                        BlockKey        = ' '                              &
                      )
      Reqs%PdfReqs(i)%iSigmaC = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        Reqs%PdfReqs(i)%iSigmaC = FindEquivFieldReqIndex(ExtraFieldReq, Reqs)
      End If

    End If

    i = i + 1

  End Do

  !-------------------------------!
  ! Loop over field requirements. !
  !-------------------------------!

  i = 1

  Do While (i <= Reqs%nFieldReqs)

    FieldReq => Reqs%FieldReqs(i)

    !----------------------------------------!
    ! Set up iProcesses, iDGrids and iTGrid. !
    !----------------------------------------!

    Do j = 1, FieldReq%nProcesses
      iProcesses(j) = FindProcessIndex(FieldReq%ProcessNames(j), Reqs)
      Process => Reqs%Processes(iProcesses(j))
      If (Process%Code /= P_AvInt) Then
        iDGrids(j) = FindDGridIndex(Process%DGridName, Grids)
      Else
        iDGrids(j) = 0
      End If
    End Do

    If (FieldReq%TGridName /= ' ') Then
      iTGrid = FindTGridIndex(FieldReq%TGridName, Grids)
    Else
      iTGrid = 0
    End If

    !--------------------------------------------------------!
    ! Check processing steps are valid when using P64 reals. !
    ! Note: only averages/integrals are permitted here       !
    ! $$ Remove this block once all P64 output is treated in !
    !    a consistent way (extended to pdfs, etc.)           !
    !--------------------------------------------------------!

    Do j = 1, FieldReq%nProcesses
      Process => Reqs%Processes(iProcesses(j))
      If (FieldReq%TypeType == '8' .and. Process%Code /= P_AvInt) Then
        Call Message(                                                             &
               'FATAL ERROR: Processing step "'                                // &
               Trim(Process%Name)                                              // &
               '" is not permitted for a quantity stored with P64 precision',     &
               3                                                                  &
             )
      End If
    End Do

    !---------------------------------------------------!
    ! Check processing steps and data grids consistent. !
    !---------------------------------------------------!

    ! Check data grids OK for processing steps.
    Do j = 1, FieldReq%nProcesses
      Process => Reqs%Processes(iProcesses(j))
      If (Process%Code /= P_AvInt .and. Process%Code /= P_Prob) Then
        DGrid => Grids%DGrids(iDGrids(j))
        If (FloatingDGrid(DGrid)) Then
          Call Message(                                                              &
                 'FATAL ERROR: Processing step "'                                 // &
                 Trim(Process%Name)                                               // &
                 '" uses a floating data grid which is not allowed except for a ' // &
                 'probability processing step',                                      &
                 3                                                                   &
               )
        End If
      End If
      If (Process%Code == P_Percent) Then
        Do k = 1, nDInDGrid(DGrid)
          If (DInDGrid(DGrid, k, 0) < 0.0 .or. DInDGrid(DGrid, k, 0) > 100.0) Then
            Call Message(                                                      &
                   'FATAL ERROR: Processing step "'                         // &
                   Trim(Process%Name)                                       // &
                   '" uses a data grid with values outside [0, 100] '       // &
                   'which is not allowed for a percentile processing step',    &
                   3                                                           &
                 )
          End If
        End Do
      End If
      If (Process%Code == P_MaxMin) Then
        Do k = 1, nDInDGrid(DGrid)
          If (DInDGrid(DGrid, k, 0) /= 0.0 .and. DInDGrid(DGrid, k, 0) /= 100.0) Then
            Call Message(                                                    &
                   'FATAL ERROR: Processing step "'                       // &
                   Trim(Process%Name)                                     // &
                   '" uses a data grid with values other than 0 and 100 ' // &
                   'which is not allowed for a max/min processing step',     &
                   3                                                         &
                 )
          End If
        End Do
      End If
    End Do

    ! Check percentile processing steps are preceeded by probability steps with a D-grid with more than one
    ! point or a floating D-grid.
    Do j = 1, FieldReq%nProcesses
      Process => Reqs%Processes(iProcesses(j))
      If (Process%Code == P_Percent) Then
        If (j == 1) Then
          Call Message(                                                         &
                 'FATAL ERROR: Processing step "'                            // &
                 Trim(Process%Name)                                          // & ! $$ these errors should
                 '" is a percentile processing step which is not preceeded ' // & ! give field req name.
                 'by a probability processing step.',                           &
                 3                                                              &
               )
        End If
        Process => Reqs%Processes(iProcesses(j - 1))
        If (Process%Code /= P_Prob) Then
          Call Message(                                                         &
                 'FATAL ERROR: Processing step "'                            // &
                 Trim(Process%Name)                                          // &
                 '" is a percentile processing step which is not preceeded ' // &
                 'by a probability processing step.',                           &
                 3                                                              &
               )
        End If
        DGrid => Grids%DGrids(iDGrids(j - 1))
        If (.not.(nDInDGrid(DGrid) > 1 .or. FloatingDGrid(DGrid))) Then
          Call Message(                                                                     &
                 'FATAL ERROR: Processing step "'                                        // &
                 Trim(Process%Name)                                                      // &
                 '" is a percentile processing step which is not preceeded '             // &
                 'by a probability processing step with a data grid with more than one ' // &
                 'point or a floating D-grid.',                                             &
                 3                                                                          &
               )
        End If
      End If
    End Do

    ! Allow only one processing step across ensemble of cases.
    Do j = 1, FieldReq%nProcesses
      Process => Reqs%Processes(iProcesses(j))
      If (Process%Ensemble) Then
        Do k = j + 1, FieldReq%nProcesses
          Process => Reqs%Processes(iProcesses(k))
          If (Process%Ensemble) Then
            Call Message(                                                                 &
                   'FATAL ERROR: Field requirement "'                                  // &
                   Trim(FieldReq%Name)                                                 // &
                   '" has more than one processing step across the ensemble of cases',    &
                   3                                                                      &
                 )
          End If
        End Do
        Exit
      End If
    End Do

    ! Allow only one processing step with D-grids with more than one point or with floating D-grids,
    ! discounting probability processing steps which are followed immediately by a percentile processing step.
    ! In addition, after such a processing step, not even probability processing steps with a following
    ! percentile step are allowed to have a D-grid with more than one point or a floating D-grid.
    Do j = 1, FieldReq%nProcesses
      Process => Reqs%Processes(iProcesses(j))
      If (Process%Code == P_AvInt) Cycle
      If (Process%Code == P_Prob .and. j + 1 <= FieldReq%nProcesses) Then
        Process1 => Reqs%Processes(iProcesses(j + 1))
        If (Process1%Code == P_Percent) Cycle
      End If
      DGrid => Grids%DGrids(iDGrids(j))
      If (nDInDGrid(DGrid) > 1 .or. FloatingDGrid(DGrid)) Then
        Do k = j + 1, FieldReq%nProcesses
          DGrid => Grids%DGrids(iDGrids(k))
          If (nDInDGrid(DGrid) > 1 .or. FloatingDGrid(DGrid)) Then
            Call Message(                                                                                    &
                   'FATAL ERROR: Field requirement "'                                                     // &
                   Trim(FieldReq%Name)                                                                    // &
                   '" has more than one processing step which requires an extra dimension in the output',    &
                   3                                                                                         &
                 )
          End If
        End Do
        Exit
      End If
    End Do

    ! Check floating data grids are only used for probability processing steps with a following percentile
    ! processing step or for probability processing steps with no following processing step. Note we already
    ! have checked floating data grids are only used for probability processing steps.
    Do j = 1, FieldReq%nProcesses
      Process => Reqs%Processes(iProcesses(j))
      If (Process%Code == P_Prob) Then
        DGrid => Grids%DGrids(iDGrids(j))
        If (FloatingDGrid(DGrid)) Then
          If (j + 1 <= FieldReq%nProcesses) Then
            Process1 => Reqs%Processes(iProcesses(j + 1))
            If (Process1%Code /= P_Percent) Then
              Call Message(                                                                              &
                     'FATAL ERROR: Field requirement "'                                               // &
                     Trim(FieldReq%Name)                                                              // &
                     '" uses a floating data grid and is followed by a processing step other than a ' // &
                     'percentile processing step',                                                       &
                     3                                                                                   &
                   )
            End If
          End If
        End If
      End If
    End Do

    !--------------------------------------------------------!
    ! Check processing steps, grids and quantity consistent. !
    !--------------------------------------------------------!

    If (FieldReq%nProcesses > 0) Then

      Process => Reqs%Processes(iProcesses(FieldReq%nProcesses))

      ! Exclude non-allowed combinations - T dependencies.

      ! No T dependence of quantity and no T-grid.
      If (Scan(QInfo(FieldReq%iQuantity), 'T') == 0) Then
        ! Exclude the following processing steps:
        !   T window defined by grid
        !   Rolling T window (except where indistinguishable from "Not depending on T")
        !   T window from T = -infinity
        !   T window over all T.
        If (                                &
          Process%UseTGrid             .or. &
          Process%T      /= ZeroTime() .or. &
          Process%TAvInt == P_Int           &
        ) Then
          Call Message(                                                               &
                 'FATAL ERROR in reading Output Requirements - Fields: '           // &
                 'The quantity "'                                                  // &
                 Trim(FieldReq%Quantity)                                           // &
                 '" used in field requirement "'                                   // &
                 Trim(FieldReq%Name)                                               // &
                 '" does not depend on T and so cannot be processed as specified',    &
                 3                                                                    &
               )
        End If
      End If

      ! T dependence of quantity and no T-grid.
      If (Scan(QInfo(FieldReq%iQuantity), 'T') /= 0 .and. FieldReq%TGridName == ' ') Then
        ! Exclude the following processing steps:
        !   T window defined by grid
        !   Rolling T window (except where indistinguishable from "Not depending on T").
        !   However for FieldReq%nProcesses = 1 we exclude all of "Rolling T window"/"Not depending on T".
        If (                                                               &
          Process%UseTGrid                                            .or. &
          (Process%T  > ZeroTime() .and. .not.IsInfFuture(Process%T)) .or. &
          (Process%T == ZeroTime() .and. Process%TAvInt == P_Int    ) .or. &
          (Process%T == ZeroTime() .and. FieldReq%nProcesses == 1   )      &
        ) Then
          Call Message(                                                                               &
                 'FATAL ERROR in reading Output Requirements - Fields: '                           // &
                 'The field requirement "'                                                         // &
                 Trim(FieldReq%Name)                                                               // &
                 '" does not have a T-grid and so cannot be produced by the processing specified',    &
                 3                                                                                    &
               )
        End If
      End If

      ! T dependence of quantity and T-grid.
      If (FieldReq%TGridName /= ' ') Then
        ! Exclude the following processing steps:
        !   T window over all T, dT > 0 (except where indistinguishable from "T window from T = -infinity").
        If (Process%T0 /= ReferenceTime()) Then
          Call Message(                                                                     &
                 'FATAL ERROR in reading Output Requirements - Fields: '                 // &
                 'The field requirement "'                                               // &
                 Trim(FieldReq%Name)                                                     // &
                 '" has a T-grid and so cannot be produced by the processing specified',    &
                 3                                                                          &
               )
        End If
      End If

      ! Exclude non-allowed combinations - H dependencies.

      ! No H dependence of quantity and no H-grid.
      If (Scan(QInfo(FieldReq%iQuantity), 'H') == 0) Then
        ! Exclude the following processing steps:
        !   H window defined by grid
        !   Rolling H window (except where indistinguishable from "Not depending on H")
        !   H window over all H.
        If (                           &
          Process%UseHGrid        .or. &
          Process%X      /= 0.0   .or. &
          Process%Y      /= 0.0   .or. &
          Process%HAvInt == P_Int      &
        ) Then
          Call Message(                                                                   &
                 'FATAL ERROR in reading Output Requirements - Fields: '               // &
                 'The quantity "'                                                      // &
                 Trim(FieldReq%Quantity)                                               // &
                 '" used in field requirement "'                                       // &
                 Trim(FieldReq%Name)                                                   // &
                 '" does not depend on X & Y and so cannot be processed as specified',    &
                 3                                                                        &
               )
        End If
      End If

      ! H dependence of quantity and no H-grid.
      If (Scan(QInfo(FieldReq%iQuantity), 'H') /= 0 .and. FieldReq%HGridName == ' ') Then
        ! Exclude the following processing steps:
        !   H window defined by grid
        !   Rolling H window (except where indistinguishable from "Not depending on H").
        !   However for FieldReq%nProcesses = 1 we exclude all of "Rolling H window"/"Not depending on H".
        If (                                                                            &
          Process%UseHGrid                                                         .or. &
          Process%X > 0.0                                                          .or. &
          Process%Y > 0.0                                                          .or. &
          (Process%X == 0.0 .and. Process%Y == 0.0 .and. Process%HAvInt == P_Int ) .or. &
          (Process%X == 0.0 .and. Process%Y == 0.0 .and. FieldReq%nProcesses == 1)      &
        ) Then
          Call Message(                                                                                &
                 'FATAL ERROR in reading Output Requirements - Fields: '                            // &
                 'The field requirement "'                                                          // &
                 Trim(FieldReq%Name)                                                                // &
                 '" does not have an H-grid and so cannot be produced by the processing specified',    &
                 3                                                                                     &
               )
        End If
      End If

      ! H dependence of quantity and H-grid.
      If (FieldReq%HGridName /= ' ') Then
        ! Exclude the following processing steps:
        !   H window over all H.
        If (                   &
          Process%X < 0.0 .or. &
          Process%Y < 0.0      &
        ) Then
          Call Message(                                                                      &
                 'FATAL ERROR in reading Output Requirements - Fields: '                  // &
                 'The field requirement "'                                                // &
                 Trim(FieldReq%Name)                                                      // &
                 '" has an H-grid and so cannot be produced by the processing specified',    &
                 3                                                                           &
               )
        End If
      End If

      ! Exclude non-allowed combinations - Z dependencies.

      ! No Z dependence of quantity and no Z-grid.
      If (Scan(QInfo(FieldReq%iQuantity), 'Z') == 0) Then
        ! Exclude the following processing steps:
        !   Z window defined by grid
        !   Rolling Z window (except where indistinguishable from "Not depending on Z")
        !   Z window defined by b-layer
        !   Z window over all Z.
        If (                           &
          Process%UseZGrid        .or. &
          Process%Z      /= 0.0   .or. &
          Process%ZAvInt == P_Int .or. &
          Process%BL                   &
        ) Then
          Call Message(                                                               &
                 'FATAL ERROR in reading Output Requirements - Fields: '           // &
                 'The quantity "'                                                  // &
                 Trim(FieldReq%Quantity)                                           // &
                 '" used in field requirement "'                                   // &
                 Trim(FieldReq%Name)                                               // &
                 '" does not depend on Z and so cannot be processed as specified',    &
                 3                                                                    &
               )
        End If
      End If

      ! Z dependence of quantity and no Z-grid.
      If (Scan(QInfo(FieldReq%iQuantity), 'Z') /= 0 .and. FieldReq%ZGridName == ' ') Then
        ! Exclude the following processing steps:
        !   Z window defined by grid
        !   Rolling Z window (except where indistinguishable from "Not depending on Z").
        !   However for FieldReq%nProcesses = 1 we exclude all of "Rolling Z window"/"Not depending on Z".
        If (                                                                           &
          Process%UseZGrid                                                        .or. &
          Process%Z > 0.0                                                         .or. &
          (Process%Z == 0.0 .and. Process%ZAvInt == P_Int  .and. .not.Process%BL) .or. &
          (Process%Z == 0.0 .and. FieldReq%nProcesses == 1 .and. .not.Process%BL)      &
        ) Then
          Call Message(                                                                               &
                 'FATAL ERROR in reading Output Requirements - Fields: '                           // &
                 'The field requirement "'                                                         // &
                 Trim(FieldReq%Name)                                                               // &
                 '" does not have a Z-grid and so cannot be produced by the processing specified',    &
                 3                                                                                    &
               )
        End If
      End If

      ! Z dependence of quantity and Z-grid.
      If (FieldReq%ZGridName /= ' ') Then
        ! Exclude the following processing steps:
        !   Z window defined by b-layer
        !   Z window over all Z.
        If (                   &
          Process%Z < 0.0 .or. &
          Process%BL           &
        ) Then
          Call Message(                                                                     &
                 'FATAL ERROR in reading Output Requirements - Fields: '                 // &
                 'The field requirement "'                                               // &
                 Trim(FieldReq%Name)                                                     // &
                 '" has a Z-grid and so cannot be produced by the processing specified',    &
                 3                                                                          &
               )
        End If
      End If

    End If

    ! Check Process%dT | TGrid%dT (including checking that Process%dT is finite).
    If (FieldReq%TGridName /= ' ') Then

      TGrid => Grids%TGrids(iTGrid)

      If (FieldReq%nProcesses > 0) Then

        Process => Reqs%Processes(iProcesses(FieldReq%nProcesses))

        If (Process%dT > ZeroTime()) Then
          If (IsInfFuture(Process%dT)) Then
            Call Message(                                                                                  &
                   'FATAL ERROR processing field requirement '                                          // &
                   Trim(FieldReq%Name)                                                                  // &
                   ': The time interval at which contributions to a processing step '                   // &
                   '(i.e. an averaging/integrating/probability/percentile/max/min/moment calculation) ' // &
                   'are required must divide the interval at which the results from the processing '    // &
                   'step are required (note infinity does not divide a finite interval)',                  &
                   3                                                                                       &
                 )
          Else
            If (TGrid%dT /= Process%dT * (TGrid%dT / Process%dT)) Then
              Call Message(                                                                                  &
                     'FATAL ERROR processing field requirement '                                          // &
                     Trim(FieldReq%Name)                                                                  // &
                     ': The time interval at which contributions to a processing step '                   // &
                     '(i.e. an averaging/integrating/probability/percentile/max/min/moment calculation) ' // &
                     'are required must divide the interval at which the results from the processing '    // &
                     'step are required (note infinity does not divide a finite interval)',                  &
                     3                                                                                       &
                   )
            End If
          End If
        End If

      End If

    End If

    ! $$ Could relax this test if allow more complex grids.
    ! More complex grids could also reduce output calculation costs (e.g. if only want one hour av
    ! every 24 hours)

    ! $$ XYZ too - but need tolerance for division testing.

    !----------------------------------------------------------------------------------!
    ! Set up field requirements with information from processing steps and data grids. !
    !----------------------------------------------------------------------------------!

    FieldReq%ComplexProc   = .false.
    FieldReq%Ensemble      = .false.
    FieldReq%AvEnsemble    = .false.
    FieldReq%AvTAvInt      = 0
    FieldReq%AvT           = ZeroTime()
    FieldReq%sAvT          = ZeroShortTime()
    FieldReq%AvUseTGrid    = .false.
    FieldReq%AvHAvInt      = 0
    FieldReq%AvX           = 0.0
    FieldReq%AvY           = 0.0
    FieldReq%AvUseHGrid    = .false.
    FieldReq%AvZAvInt      = 0
    FieldReq%AvZ           = 0.0
    FieldReq%AvUseZGrid    = .false.
    FieldReq%AvBL          = .false.
    FieldReq%PCode         = 0
    FieldReq%PEnsemble     = .false.
    FieldReq%PT            = ZeroTime()
    FieldReq%sPT           = ZeroShortTime()
    FieldReq%PUseTGrid     = .false.
    FieldReq%PX            = 0.0
    FieldReq%PY            = 0.0
    FieldReq%PUseHGrid     = .false.
    FieldReq%PZ            = 0.0
    FieldReq%PUseZGrid     = .false.
    FieldReq%PBL           = .false.
    FieldReq%DGridName     = ' '
    FieldReq%DGridFloating = .false.

    ! ComplexProc.
    Do j = 1, FieldReq%nProcesses
      Process => Reqs%Processes(iProcesses(j))
      If (j == 1                   .and. Process%Code == P_AvInt) Cycle
      If (j == FieldReq%nProcesses .and. Process%Code /= P_AvInt) Exit
      Process1 => Reqs%Processes(iProcesses(FieldReq%nProcesses))
      If (j == FieldReq%nProcesses - 1 .and. Process%Code == P_Prob .and. Process1%Code == P_Percent) Exit
      FieldReq%ComplexProc = .true.
    End Do

    ! Ensemble.
    Do j = 1, FieldReq%nProcesses
      Process => Reqs%Processes(iProcesses(j))
      If (Process%Ensemble) Then
        FieldReq%Ensemble = .true.
      End If
    End Do

    ! AvEnsemble, AvTAvInt, AvT, sAvT, AvUseTGrid, AvHAvInt, AvX, AvY, AvUseHGrid, AvZAvInt, AvZ, AvUseZGrid
    ! and AvBL.
    If (FieldReq%nProcesses > 0) Then
      Process => Reqs%Processes(iProcesses(1))
      If (Process%Code == P_AvInt) Then
        FieldReq%AvEnsemble = Process%Ensemble
        FieldReq%AvTAvInt   = Process%TAvInt
        FieldReq%AvT        = Process%T
        FieldReq%sAvT       = Process%sT
        FieldReq%AvUseTGrid = Process%UseTGrid
        FieldReq%AvHAvInt   = Process%HAvInt
        FieldReq%AvX        = Process%X
        FieldReq%AvY        = Process%Y
        FieldReq%AvUseHGrid = Process%UseHGrid
        FieldReq%AvZAvInt   = Process%ZAvInt
        FieldReq%AvZ        = Process%Z
        FieldReq%AvUseZGrid = Process%UseZGrid
        FieldReq%AvBL       = Process%BL
        If (                                 &
          Process%TAvInt == P_Av       .and. &
          Process%T      == ZeroTime() .and. &
          .not.Process%UseTGrid              &
        ) FieldReq%AvTAvInt = 0
        If (                           &
          Process%HAvInt == P_Av .and. &
          Process%X      == 0.0  .and. &
          Process%Y      == 0.0  .and. &
          .not.Process%UseHGrid        &
        ) FieldReq%AvHAvInt = 0
        If (                           &
          Process%ZAvInt == P_Av .and. &
          Process%Z      == 0.0  .and. &
          .not.Process%UseZGrid  .and. &
          .not.Process%BL              &
        ) FieldReq%AvZAvInt = 0
      End If
    End If

    ! PCode, PEnsemble, PT, sPT, PUseTGrid, PX, PY, PUseHGrid, PZ, PUseZGrid and PBL.
    If (.not. FieldReq%ComplexProc) Then
      If (FieldReq%nProcesses > 0) Then
        Process => Reqs%Processes(iProcesses(FieldReq%nProcesses))
        If (Process%Code /= P_AvInt) Then
          FieldReq%PCode = Process%Code
          If (Process%Code == P_Percent) Then
            Process => Reqs%Processes(iProcesses(FieldReq%nProcesses - 1))
          End If
          FieldReq%PEnsemble = Process%Ensemble
          FieldReq%PT        = Process%T
          FieldReq%sPT       = Process%sT
          FieldReq%PUseTGrid = Process%UseTGrid
          FieldReq%PX        = Process%X
          FieldReq%PY        = Process%Y
          FieldReq%PUseHGrid = Process%UseHGrid
          FieldReq%PZ        = Process%Z
          FieldReq%PUseZGrid = Process%UseZGrid
          FieldReq%PBL       = Process%BL
        End If
      End If
    End If

    ! DGridName and DGridFloating.
    Do j = 1, FieldReq%nProcesses
      Process => Reqs%Processes(iProcesses(j))
      If (Process%Code == P_AvInt) Cycle
      If (Process%Code == P_Prob .and. j + 1 <= FieldReq%nProcesses) Then
        Process1 => Reqs%Processes(iProcesses(j + 1))
        If (Process1%Code == P_Percent) Cycle
      End If
      DGrid => Grids%DGrids(iDGrids(j))
      If (nDInDGrid(DGrid) > 1 .or. FloatingDGrid(DGrid)) Then
        FieldReq%DGridName     = Process%DGridName
        FieldReq%DGridFloating = FloatingDGrid(DGrid)
        Exit
      End If
    End Do
    If (FieldReq%DGridName == ' ') Then
      Do j = 1, FieldReq%nProcesses
        Process => Reqs%Processes(iProcesses(j))
        If (Process%Code == P_AvInt) Cycle
        If (Process%Code == P_Prob .and. j + 1 <= FieldReq%nProcesses) Then
          Process1 => Reqs%Processes(iProcesses(j + 1))
          If (Process1%Code == P_Percent) Cycle
        End If
        DGrid => Grids%DGrids(iDGrids(j))
        If (FieldReq%DGridName /= ' ') Then
          FieldReq%DGridName     = ' '
          FieldReq%DGridFloating = .false.
          Exit
        End If
        FieldReq%DGridName     = Process%DGridName
        FieldReq%DGridFloating = FloatingDGrid(DGrid)
      End Do
    End If

    ! iPQRProc (index of the processing step after which the derived quantity is calculated).
    FieldReq%iPQRProc = 0
    If (FieldReq%LEFOPQRType == 'P' .or. FieldReq%LEFOPQRType == 'Q' .or. FieldReq%LEFOPQRType == 'R') Then
      FieldReq%iPQRProc = Char2Int(QInfo(FieldReq%iQuantity)(6:6))
      FieldReq%iPQRProc = Min(FieldReq%iPQRProc, FieldReq%nProcesses)
      Do j = 1, FieldReq%iPQRProc
        Process => Reqs%Processes(iProcesses(j))
        If (Process%Code /= P_AvInt) Then
          FieldReq%iPQRProc = j - 1
          Exit
        End If
      End Do
    End If

    ! iFOQRProc (index of the processing step after which quantities of type F, O, Q or R appear in the
    ! processing chain).
    FieldReq%iFOQRProc = 0
    If (FieldReq%LEFOPQRType == 'Q' .or. FieldReq%LEFOPQRType == 'R') Then
      ! $$ currently there are no Q or R quantities. When there are, need to either calculate iFOQRProc from
      ! QInfo table or add iFOQRProc to QInfo table. Good to check table if possible in the code.
    End If

    ! iFluctuationsProc (index of processing step in which fluctuations are accounted for). If Fluctuations is
    ! true and the quantity is one for which fluctuations can be allowed for (but is not a fluctuation
    ! parameter), set to first probability/moment processing step and check there is such a processing step.
    FieldReq%iFluctuationsProc = 0
    If (FieldReq%Fluctuations .and. Scan(QInfo(FieldReq%iQuantity), 'f') /= 0) Then

      Do j = 1, FieldReq%nProcesses
        Process => Reqs%Processes(iProcesses(j))
        If (                                    &
          Process%Code == P_Prob           .or. &
          Process%Code == P_Moments        .or. &
          Process%Code == P_CentralMoments      &
        ) Then
          FieldReq%iFluctuationsProc = j
          Exit
        End If
      End Do
      If (FieldReq%iFluctuationsProc == 0) Then
        Call Message(                                      & ! $$ reword
               'FATAL ERROR: Field requirement "'       // &
               Trim(FieldReq%Name)                      // &
               '" cannot take account of fluctuations',    &
               3                                           &
             )
      End If
    End If

    ! NoProc. Set NoProc true if (i) processing step is an averaging/integrating step with no ensemble
    ! averaging, only a single input grid value contributing to each output grid value, and no integrating
    ! other than intrinsic integrating, or if (ii) processing step will be applied before calculating the
    ! derived field.
    FieldReq%NoProc(1:FieldReq%nProcesses) = .false.
    Do j = 1, FieldReq%nProcesses
      Process => Reqs%Processes(iProcesses(j))
      If (                                                                       &
        Process%Code == P_AvInt                                            .and. &
        .not.Process%Ensemble                                              .and. &
        Process%nT == 1                                                    .and. &
        Process%nX == 1                                                    .and. &
        Process%nY == 1                                                    .and. &
        Process%nZ == 1                                                    .and. &
        (Process%TAvInt == P_Av .or. (FieldReq%IntrinsicTAv .and. j == 1)) .and. &
        (Process%HAvInt == P_Av .or. (FieldReq%IntrinsicHAv .and. j == 1)) .and. &
        (Process%ZAvInt == P_Av .or. (FieldReq%IntrinsicZAv .and. j == 1))       &
      ) Then
        FieldReq%NoProc(j) = .true.
      End If
    End Do
    FieldReq%NoProc(1:FieldReq%iPQRProc) = .true.

    ! iLastProcess (index of last processing step with NoProc false).
    FieldReq%iLastProcess = 0
    Do j = 1, FieldReq%nProcesses
      If (.not.FieldReq%NoProc(j)) FieldReq%iLastProcess = j
    End Do

    !-----------------------------------------------------------------------------------------!
    ! Check that for intrinsic averaging/integrating there is an appropriate processing step. !
    !-----------------------------------------------------------------------------------------!

    ! Intrinsic T-averaging/integrating.
    If (FieldReq%IntrinsicTAv) Then

      If (FieldReq%nProcesses > 0) Then
        Process => Reqs%Processes(iProcesses(1))
        If (                                                         &
          Process%Code /= P_AvInt                               .or. &
          (Process%T == ZeroTime() .and. .not.Process%UseTGrid)      &
        ) Then
          Call Message(                                                                         &
                 'FATAL ERROR in reading Output Requirements - Fields: Intrinsic T '         // &
                 'averaging/integrating is specified for '                                   // &
                 Trim(FieldReq%Quantity)                                                     // &
                 ', but the first processing step (if any) is not an averaging/integrating ' // &
                 'step over a non-zero time interval',                                          &
                 3                                                                              &
               )
        End If
      Else
        Call Message(                                                                         &
               'FATAL ERROR in reading Output Requirements - Fields: Intrinsic T '         // &
               'averaging/integrating is specified for '                                   // &
               Trim(FieldReq%Quantity)                                                     // &
               ', but the first processing step (if any) is not an averaging/integrating ' // &
               'step over a non-zero time interval',                                          &
               3                                                                              &
             )
      End If

    End If

    ! Intrinsic H-averaging/integrating.
    If (FieldReq%IntrinsicHAv) Then

      If (FieldReq%nProcesses > 0) Then
        Process => Reqs%Processes(iProcesses(1))
        If (                                                                         &
          Process%Code /= P_AvInt                                               .or. &
          (Process%X == 0.0 .and. Process%Y == 0.0 .and. .not.Process%UseHGrid)      &
        ) Then
          Call Message(                                                                         &
                 'FATAL ERROR in reading Output Requirements - Fields: Intrinsic H '         // &
                 'averaging/integrating is specified for '                                   // &
                 Trim(FieldReq%Quantity)                                                     // &
                 ', but the first processing step (if any) is not an averaging/integrating ' // &
                 'step over a non-zero horizontal region',                                      &
                 3                                                                              &
               )
        End If
      Else
        Call Message(                                                                         &
               'FATAL ERROR in reading Output Requirements - Fields: Intrinsic H '         // &
               'averaging/integrating is specified for '                                   // &
               Trim(FieldReq%Quantity)                                                     // &
               ', but the first processing step (if any) is not an averaging/integrating ' // &
               'step over a non-zero horizontal region',                                      &
               3                                                                              &
             )
      End If

    End If

    ! Intrinsic Z-averaging/integrating.
    If (FieldReq%IntrinsicZAv) Then

      If (FieldReq%nProcesses > 0) Then
        Process => Reqs%Processes(iProcesses(1))
        If (                                                                        &
          Process%Code /= P_AvInt                                              .or. &
          (Process%Z == 0.0 .and. .not.Process%UseZGrid .and. .not.Process%BL)      &
        ) Then
          Call Message(                                                                         &
                 'FATAL ERROR in reading Output Requirements - Fields: Intrinsic Z '         // &
                 'averaging/integrating is specified for '                                   // &
                 Trim(FieldReq%Quantity)                                                     // &
                 ', but the first processing step (if any) is not an averaging/integrating ' // &
                 'step over a non-zero vertical region',                                        &
                 3                                                                              &
               )
        End If
      Else
        Call Message(                                                                         &
               'FATAL ERROR in reading Output Requirements - Fields: Intrinsic Z '         // &
               'averaging/integrating is specified for '                                   // &
               Trim(FieldReq%Quantity)                                                     // &
               ', but the first processing step (if any) is not an averaging/integrating ' // &
               'step over a non-zero vertical region',                                        &
               3                                                                              &
             )
      End If

    End If

    ! Ideally check no intrinsic av if type P Q R and doesn't depend on type L field. $$

    !---------------------------------------------------------------------------!
    ! Check dT/X/Y/Z = infinity only used with intrinsic averaging/integrating. !
    !---------------------------------------------------------------------------!

    If (FieldReq%nProcesses > 0) Then

      Do j = 1, FieldReq%nProcesses

        Process => Reqs%Processes(iProcesses(j))

        If (                                                                        &
          IsInfFuture(Process%dT)                                                   &
                                                                              .and. &
          (                                                                         &
            j /= 1                                                       .or.       &
            .not. FieldReq%IntrinsicTAv                                  .or.       &
            (FieldReq%LEFOPQRType /= 'L' .and. FieldReq%nProcesses == 1)            &
          )                                                                         &
        ) Then
          Call Message('FATAL ERROR in SetUpReqs', 3) ! $$ reword
        End If

        If (                                                                        &
          Process%dX == -1.0                                                        &
                                                                              .and. &
          (                                                                         &
            j /= 1                                                       .or.       &
            .not. FieldReq%IntrinsicHAv                                  .or.       &
            (FieldReq%LEFOPQRType /= 'L' .and. FieldReq%nProcesses == 1)            &
          )                                                                         &
        ) Then
          Call Message('FATAL ERROR in SetUpReqs', 3) ! $$ reword
        End If

        If (                                                                        &
          Process%dY == -1.0                                                        &
                                                                              .and. &
          (                                                                         &
            j /= 1                                                       .or.       &
            .not. FieldReq%IntrinsicHAv                                  .or.       &
            (FieldReq%LEFOPQRType /= 'L' .and. FieldReq%nProcesses == 1)            &
          )                                                                         &
        ) Then
          Call Message('FATAL ERROR in SetUpReqs', 3) ! $$ reword
        End If

        If (                                                                        &
          Process%dZ == -1.0                                                        &
                                                                              .and. &
          (                                                                         &
            j /= 1                                                       .or.       &
            .not. FieldReq%IntrinsicZAv                                  .or.       &
            (FieldReq%LEFOPQRType /= 'L' .and. FieldReq%nProcesses == 1)            &
          )                                                                         &
        ) Then
          Call Message('FATAL ERROR in SetUpReqs', 3) ! $$ reword
        End If

      End Do

    End If

    !-------------------------------------------------------------------------!
    ! Check T/X/Y/Z = infinity not used for quantities of type F, O, Q and R. !
    !-------------------------------------------------------------------------!

    If (                               &
      FieldReq%LEFOPQRType == 'F' .or. &
      FieldReq%LEFOPQRType == 'O' .or. &
      FieldReq%LEFOPQRType == 'Q' .or. &
      FieldReq%LEFOPQRType == 'R'      &
    ) Then

      If (FieldReq%nProcesses > 0) Then

        Do j = FieldReq%iPQRProc + 1, FieldReq%nProcesses

          Process => Reqs%Processes(iProcesses(j))

          If (IsInfFuture(Process%T)) Then
            Call Message('FATAL ERROR in SetUpReqs', 3) ! $$ reword
          End If

          If (Process%X == -1.0) Then
            Call Message('FATAL ERROR in SetUpReqs', 3) ! $$ reword
          End If

          If (Process%Y == -1.0) Then
            Call Message('FATAL ERROR in SetUpReqs', 3) ! $$ reword
          End If

          If (Process%Z == -1.0) Then
            Call Message('FATAL ERROR in SetUpReqs', 3) ! $$ reword
          End If

        End Do

      End If

    End If

    !-------------------------------------------------------------!
    ! Check fluctuations options and processing steps consistent. !
    !-------------------------------------------------------------!

    If (FieldReq%Fluctuations) Then

      ! Check fluctuations are on only if (1) the quantity is one for which fluctuation calculations are
      ! possible and there is a probability or moment calculation step (the latter has already been checked
      ! above), or (2) the quantity is a fluctuation parameter. In both cases the only processing steps
      ! allowed before the probability/moment calculation step or the calculation of the fluctuation parameter
      ! is a single averaging/integrating step not involving an ensemble average.
      If (Scan(QInfo(FieldReq%iQuantity), 'f') /= 0) Then

        Do j = 1, FieldReq%iFluctuationsProc - 1
          Process => Reqs%Processes(iProcesses(j))
          If (Process%Code /= P_AvInt .or. Process%Ensemble .or. j /= 1) Then
            Call Message(                                      &
                   'FATAL ERROR: Field requirement "'       // &
                   Trim(FieldReq%Name)                      // &
                   '" cannot take account of fluctuations',    &
                   3                                           &
                 )
          End If
        End Do

      Else If (Scan(QInfo(FieldReq%iQuantity), 'p') /= 0) Then

        Do j = 1, FieldReq%iPQRProc
          Process => Reqs%Processes(iProcesses(j))
          If (Process%Code /= P_AvInt .or. Process%Ensemble .or. j /= 1) Then
            Call Message(                                      &
                   'FATAL ERROR: Field requirement "'       // &
                   Trim(FieldReq%Name)                      // &
                   '" cannot take account of fluctuations',    &
                   3                                           &
                 )
          End If
        End Do

      End If

      ! Check that separate fluctuation scales are not used when integrating rather than averaging is being
      ! considered.
      If (Scan(QInfo(FieldReq%iQuantity), 'f') /= 0) Then

        Do j = 1, FieldReq%iFluctuationsProc - 1
          Process => Reqs%Processes(iProcesses(j))
          If (                                                     &
            (FieldReq%UseFlAvT .and. Process%TAvInt == P_Int) .or. &
            (FieldReq%UseFlAvH .and. Process%HAvInt == P_Int) .or. &
            (FieldReq%UseFlAvH .and. Process%HAvInt == P_Int) .or. &
            (FieldReq%UseFlAvZ .and. Process%ZAvInt == P_Int)      &
          ) Then
            Call Message(                                      & ! $$ reword
                   'FATAL ERROR: Field requirement "'       // &
                   Trim(FieldReq%Name)                      // &
                   '" cannot take account of fluctuations',    &
                   3                                           &
                 )
          End If
        End Do

      Else If (Scan(QInfo(FieldReq%iQuantity), 'p') /= 0) Then

        Do j = 1, FieldReq%iPQRProc
          Process => Reqs%Processes(iProcesses(j))
          If (                                                     &
            (FieldReq%UseFlAvT .and. Process%TAvInt == P_Int) .or. &
            (FieldReq%UseFlAvH .and. Process%HAvInt == P_Int) .or. &
            (FieldReq%UseFlAvH .and. Process%HAvInt == P_Int) .or. &
            (FieldReq%UseFlAvZ .and. Process%ZAvInt == P_Int)      &
          ) Then
            Call Message(                                      & ! $$ reword
                   'FATAL ERROR: Field requirement "'       // &
                   Trim(FieldReq%Name)                      // &
                   '" cannot take account of fluctuations',    &
                   3                                           &
                 )
          End If
        End Do

      End If

    End If

    ! $$ Other checks? currently flucs only available in very
    ! restricted conditions. Perhaps at point of asking for derived field.
    ! $$ check flucs are on for a fluc parameter?

    !------------------------------------------------------!
    ! Checks for decay of deposition with power law decay. !
    !------------------------------------------------------!

    If (FieldReq%DecayDep) Then

      Species => Specieses%Specieses(FindSpeciesIndex(FieldReq%SpeciesName, Specieses))

      If (Species%UsePowerLawDecay) Then

        If (FieldReq%SourceName == ' ') Then

          Do j = 1, Sources%nSources

            Source => Sources%Sources(j)
            If (Sources%Sources(1)%StartTime /= Source%StartTime .or. Source%StartTime /= Source%StopTime) Then
              Call Message(                                                                             &
                     'FATAL ERROR: Field Requirement "'                                              // &
                     Trim(FieldReq%Name)                                                             // &
                     '" involves (i) decay of deposited material, (ii) species "'                    // &
                     Trim(FieldReq%SpeciesName)                                                      // &
                     '" which undergoes power law decay, and (iii) source(s) which do not emit all ' // &
                     'their material at the same time. This is not an allowed combination',             &
                     3                                                                                  &
                   )
            End If

          End Do

        Else

          Source => Sources%Sources(FindSourceIndex(FieldReq%SourceName, Sources))
          If (Source%StartTime /= Source%StopTime) Then
            Call Message(                                                                             &
                   'FATAL ERROR: Field Requirement "'                                              // &
                   Trim(FieldReq%Name)                                                             // &
                   '" involves (i) decay of deposited material, (ii) species "'                    // &
                   Trim(FieldReq%SpeciesName)                                                      // &
                   '" which undergoes power law decay, and (iii) source(s) which do not emit all ' // &
                   'their material at the same time. This is not an allowed combination',             &
                   3                                                                                  &
                 )
          End If

        End If

      End If

    End If
    ! $$ add check (in species) that power law decay and chemistry can't be used together.
    ! $$ if sources individually instantaneous, could produce derived requests for results
    !    from individual sources, or, in all cases, could use intrinsic averaging (perhaps an
    !    'allow overlapping regions' switch.

    !--------------------------------------------------------------------------------------!
    ! Check output options consistent for field requirements within the same output group. !
    !--------------------------------------------------------------------------------------!

    If (FieldReq%OutputGroup /= ' ') Then

      Do j = 1, i - 1

        FieldReq1 => Reqs%FieldReqs(j)

        If (FieldReq1%OutputGroup .CIEq. FieldReq%OutputGroup) Then

          If (                                                           &
            (FieldReq1%TAcross       .neqv. FieldReq%TAcross      ) .or. &
            (FieldReq1%SAcross       .neqv. FieldReq%SAcross      ) .or. &
            (FieldReq1%XAcross       .neqv. FieldReq%XAcross      ) .or. &
            (FieldReq1%YAcross       .neqv. FieldReq%YAcross      ) .or. &
            (FieldReq1%ZAcross       .neqv. FieldReq%ZAcross      ) .or. &
            (FieldReq1%DAcross       .neqv. FieldReq%DAcross      ) .or. &
            (FieldReq1%TSeparateFile .neqv. FieldReq%TSeparateFile) .or. &
            (FieldReq1%SSeparateFile .neqv. FieldReq%SSeparateFile) .or. &
            (FieldReq1%XSeparateFile .neqv. FieldReq%XSeparateFile) .or. &
            (FieldReq1%YSeparateFile .neqv. FieldReq%YSeparateFile) .or. &
            (FieldReq1%ZSeparateFile .neqv. FieldReq%ZSeparateFile) .or. &
            (FieldReq1%DSeparateFile .neqv. FieldReq%DSeparateFile) .or. &
            (FieldReq1%TNewFile      .neqv. FieldReq%TNewFile     ) .or. &
            (FieldReq1%GridIndices   .neqv. FieldReq%GridIndices  ) .or. &
            (FieldReq1%AlignCols     .neqv. FieldReq%AlignCols    ) .or. &
            (FieldReq1%ZeroLines     .neqv. FieldReq%ZeroLines    ) .or. &
            (FieldReq1%Flush         .neqv. FieldReq%Flush        ) .or. &
            (FieldReq1%NameII        .neqv. FieldReq%NameII       )      &
          ) Then

            Call Message(                                                             &
                   'FATAL ERROR in adding a field requirement with output group "' // &
                   Trim(FieldReq%OutputGroup)                                      // &
                   '": all field requirements in the same output group must have ' // &
                   'the same values '                                              // &
                   'for "Across", "Separate File" and "Output Format" '            // &
                   '(except possibly for the End of output option)',                  &
                   3                                                                  &
                 )

          Else If (FieldReq1%Ensemble .neqv. FieldReq%Ensemble) Then

            Call Message(                                                             &
                   'FATAL ERROR in adding a field requirement with output group "' // &
                   Trim(FieldReq%OutputGroup)                                      // &
                   '": all field requirements in the same output group must have ' // &
                   'the same value for "Ensemble Av?"',                               & ! reword
                   3                                                                  &
                 )

          Else If (                                                       &
            (                                                             &
              .not. (FieldReq%TGridName .CIEq. FieldReq1%TGridName)       &
              .and.                                                       &
              .not. (FieldReq%TAcross .and. .not. FieldReq%TSeparateFile) &
            )                                                             &
            .or.                                                          &
            (                                                             &
              .not. (FieldReq%SGridName .CIEq. FieldReq1%SGridName)       &
              .and.                                                       &
              .not. (FieldReq%SAcross .and. .not. FieldReq%SSeparateFile) &
            )                                                             &
            .or.                                                          &
            (                                                             &
              .not. (FieldReq%HGridName .CIEq. FieldReq1%HGridName)       &
              .and.                                                       &
              .not. (FieldReq%XAcross .and. .not. FieldReq%XSeparateFile) &
            )                                                             &
            .or.                                                          &
            (                                                             &
              .not. (FieldReq%HGridName .CIEq. FieldReq1%HGridName)       &
              .and.                                                       &
              .not. (FieldReq%YAcross .and. .not. FieldReq%YSeparateFile) &
            )                                                             &
            .or.                                                          &
            (                                                             &
              .not. (FieldReq%ZGridName .CIEq. FieldReq1%ZGridName)       &
              .and.                                                       &
              .not. (FieldReq%ZAcross .and. .not. FieldReq%ZSeparateFile) &
            )                                                             &
            .or.                                                          &
            (                                                             &
              .not. (FieldReq%DGridName .CIEq. FieldReq1%DGridName)       &
              .and.                                                       &
              .not. (FieldReq%DAcross .and. .not. FieldReq%DSeparateFile) &
            )                                                             &
          ) Then

            Call Message(                                                                            &
                   'FATAL ERROR in adding a field requirement with output group "'                // &
                   Trim(FieldReq%OutputGroup)                                                     // &
                   '": all field requirements in the same output group must have '                // &
                   'the same grids (including any data grid), except in directions or times for ' // &
                   'which "Across" is set and "Separate File" isn''t',                               &
                   3                                                                                 &
                 )

          ! Fields in same output group with no TGrid must have no TGrid for same reason.
          Else If (                                                        &
            (                                                              &
              (FieldReq%TGridName == ' ' .and. FieldReq1%TGridName == ' ') &
              .and.                                                        &
              .not. (                                                      &
                Scan(QInfo(FieldReq%iQuantity), 'T') == 0                  &
                .eqv.                                                      &
                Scan(QInfo(FieldReq1%iQuantity), 'T') == 0                 &
               )                                                           &
              .and.                                                        &
              .not. (FieldReq%TAcross .and. .not. FieldReq%TSeparateFile)  &
            )                                                              &
          ) Then

            Call Message(                                                             & ! $$ reword
                   'FATAL ERROR in adding a field requirement with output group "' // &
                   Trim(FieldReq%OutputGroup)                                      // &
                   '": all field requirements in the same output group must have ' // &
                   'the same grids, except in directions or times for '            // &
                   'which "Across" is set and "Separate File" isn''t',                &
                   3                                                                  &
                 )

          End If

          Exit

        End If

      End Do

    End If

    !--------------------------------------------------------------------!
    ! Check X/Y Across/SeparateFile are the same for unstructured grids. !
    !--------------------------------------------------------------------!

    If (FieldReq%HGridName /= ' ') Then
      If (Grids%HGrids(FindHGridIndex(FieldReq%HGridName, Grids))%Unstructured) Then
        If (                                                          &
          (FieldReq%XAcross       .neqv. FieldReq%YAcross      ) .or. &
          (FieldReq%XSeparateFile .neqv. FieldReq%YSeparateFile)      &
        ) Then
          Call Message(                                                              &
                 'FATAL ERROR: For field requirements using unstructured grids, ' // &
                 'the X and Y settings for Across and Separate File must be '     // &
                 'the same',                                                         &
                 3                                                                   &
               )
        End If
      End If
    End If

    !------------------------------------------------------------------!
    ! Check that Material Unit of species is compatible with           !
    ! Material Unit of output requirement                              !
    !------------------------------------------------------------------!

    If (.not. (FieldReq%SpeciesName .CIEq. ' ') ) Then

      ! Get species material unit.
      iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
      SpeciesMaterialUnit  = MaterialUnits%MaterialUnits(                        &
                               FindMaterialUnitIndex(                            &
                                 Specieses%Specieses(iSpecies)%MaterialUnitName, &
                                 MaterialUnits                                   &
                               )                                                 &
                             )

      ! Get field req material unit and check consistency of units.
      ! $$ "material" is not appropriate for all cases - should rename

      ! 1. Mass based output quantities. Currently, only allow conversion for the following quantities:
      If ( ( FieldReq%iQuantity == Q_AirConc                ) .or. &
           ( FieldReq%iQuantity == Q_ChemistryField         ) .or. &
           ( FieldReq%iQuantity == Q_EulerianConcentration  ) .or. &
           ( FieldReq%iQuantity == Q_DryDep                 ) .or. &
           ( FieldReq%iQuantity == Q_WetDep                 ) .or. &
           ( FieldReq%iQuantity == Q_Dep                    ) .or. &
           ( FieldReq%iQuantity == Q_EulerianTotalDep       ) .or. &
           ( FieldReq%iQuantity == Q_EulerianDryDep         ) .or. &
           ( FieldReq%iQuantity == Q_EulerianWetDep         ) .or. &
           ( FieldReq%iQuantity == Q_OriginalSourceStrength ) .or. &
           ( FieldReq%iQuantity == Q_RevisedSourceStrength  )      &
         ) Then
        ! Get field req material unit and set default if appropriate.
        If ( Trim(FieldReq%MaterialUnitName) == '') Then
          FieldReqMaterialUnit = SpeciesMaterialUnit        
          FieldReq%MaterialUnitName = FieldReqMaterialUnit%Name
        Else
          FieldReqMaterialUnit = MaterialUnits%MaterialUnits(                       &
                                   FindMaterialUnitIndex(FieldReq%MaterialUnitName, &
                                     MaterialUnits                                  &
                                   )                                                &
                                 )
        End If
        ! Check whether units are allowed.
        ! Check whether unit type is allowed
        If ( (FieldReqMaterialUnit%UnitType .ne. MassUnitType     ) .and. &
             (FieldReqMaterialUnit%UnitType .ne. ActivityUnitType ) .and. &
             (FieldReqMaterialUnit%UnitType .ne. DobsonUnitType   ) .and. &
             (FieldReqMaterialUnit%UnitType .ne. UnknownUnitType  )       &
           ) Then
          Call Message( 'FATAL ERROR: Unit "'             // &
                        Trim(FieldReqMaterialUnit%Name)   // &
                        '" in Field Request "'            // &
                        Trim(FieldReq%Name)               // &
                        '" is not allowed for quantity "' // &
                        Trim(FieldReq%Quantity)           // &
                        '".', 3                              &
               )
        End If
        ! Check whether units are compatible or, for unknown units, identical
        If ( ( FieldReqMaterialUnit%UnitType .eq. UnknownUnitType     ) .and.  &
             ( FieldReqMaterialUnit          .ne. SpeciesMaterialUnit )        &
           ) Then
          Call Message( 'FATAL ERROR: Unit "'               // &
                        Trim(FieldReqMaterialUnit%Name)     // &
                        '" in Field Request "'              // &
                        Trim(FieldReq%Name)                 // &
                        '" does not match unit "'           // &
                        Trim(SpeciesMaterialUnit%Name)      // &
                        '" of species "'                    // &
                        Trim(FieldReq%SpeciesName)          // &
                        '".', 3                                &
               )
        End If
        If ( (SpeciesMaterialUnit%UnitType  .ne. FieldReqMaterialUnit%UnitType) .and. &
             (FieldReqMaterialUnit%UnitType .ne. DobsonUnitType               )       &
           ) Then
          Call Message( 'FATAL ERROR: Unit "'               // &
                        Trim(FieldReqMaterialUnit%Name)     // &
                        '" in Field Request "'              // &
                        Trim(FieldReq%Name)                 // &
                        '" is not compatible with unit "'   // &
                        Trim(SpeciesMaterialUnit%Name)      // &
                        '" of species "'                    // &
                        Trim(FieldReq%SpeciesName)          // &
                        '".', 3                                &
               )
        End If
        ! Only allow the use of Dobson units for air concentrations 
        If ( (FieldReq%iQuantity            .ne. Q_AirConc      ) .and.     &
             (FieldReqMaterialUnit%UnitType .eq. DobsonUnitType )           &
           ) Then
          Call Message( 'FATAL ERROR in Field Request "'      // &
                        Trim(FieldReq%Name)                   // &
                        '". Dobson units can only be used  '  // &
                        'for air concentrations.', 3             &
               )         
        End If
        ! If output in Dobson units is requested, the species must be given
        ! in mass units.
        If ( (FieldReqMaterialUnit%UnitType .eq. DobsonUnitType ) .and. &
             (SpeciesMaterialUnit%UnitType  .ne. MassUnitType   )       &            
           ) Then
          Call Message( 'FATAL ERROR in Field Request "'      // &
                        Trim(FieldReq%Name)                   // &
                        '". Dobson units can only be used  '  // &
                        'if species unit is mass unit.', 3       &
               )         
        End If
        ! Dobson units can only be used if 
        ! (i) the horizontal grid is specified
        ! (ii) the vertical grid is unspecified and 
        ! (iii) we do not request boundary layer averages
        If (FieldReqMaterialUnit%UnitType .eq. DobsonUnitType ) Then
          If (FieldReq%HGridName == ' ') Then
            Call Message('ERROR in setting up field request '    // &
                         Trim(FieldReq%Name)                     // &
                         ' horizontal grid has to be specified'  // &
                         ' if output is requested in Dobson units', &
                         3                                          &
                 )
          End If
          If (FieldReq%ZGridName .ne. ' ') Then
            Call Message('ERROR in setting up field request '    // &
                         Trim(FieldReq%Name)                     // &
                         ' vertical grid to be unspecified'      // &
                         ' if output is requested in Dobson units', &
                         3                                          &
                 )
          End If
          If (FieldReq%AvBL) Then
            Call Message('ERROR in setting up field request '    // &
                         Trim(FieldReq%Name)                     // &
                         ' BL average can not be used'           // &
                         ' if output is requested in Dobson units', &
                         3                                          &
                 )
          End If
        End If
      ! 2. Mixing ratio based output.
      Else If ( FieldReq%iQuantity == Q_MixingRatio .or. FieldReq%iQuantity == Q_EMixingRatio ) Then
        ! Get field req material unit and set default if appropriate.
        If ( Trim(FieldReq%MaterialUnitName) == '') Then
          FieldReqMaterialUnit = MaterialUnits%MaterialUnits(                  &
                                   FindMaterialUnitIndex('ppm', MaterialUnits) &
                                 )
          FieldReq%MaterialUnitName = FieldReqMaterialUnit%Name
        Else
          FieldReqMaterialUnit = MaterialUnits%MaterialUnits(                       &
                                   FindMaterialUnitIndex(FieldReq%MaterialUnitName, &
                                     MaterialUnits                                  &
                                   )                                                &
                                 )
        End If
        ! Check whether units are allowed.
        If ( FieldReqMaterialUnit%UnitType .ne. VolumetricUnitType ) Then
          Call Message( 'FATAL ERROR: Unit "'             // &
                        Trim(FieldReqMaterialUnit%Name)   // &
                        '" in Field Request "'            // &
                        Trim(FieldReq%Name)               // &
                        '" is not allowed for quantity "' // &
                        Trim(FieldReq%Quantity)           // &
                        '".', 3                              &
               )
        End If
        If ( SpeciesMaterialUnit%UnitType .ne. MassUnitType ) Then
          Call Message( 'FATAL ERROR: Unit "'             // &
                        Trim(FieldReqMaterialUnit%Name)   // &
                        '" of Species "'                  // &
                        Trim(FieldReq%SpeciesName)        // &
                        '" used in Field Request "'       // &
                        Trim(FieldReq%Name)               // &
                        '" is not allowed for quantity "' // &
                        Trim(FieldReq%Quantity)           // &
                        '".', 3                              &
               )
        End If
      ! 3. Other output.  
      Else
        ! Get field req material unit and set default if appropriate.
        If ( Trim(FieldReq%MaterialUnitName) == '') Then
          FieldReqMaterialUnit = SpeciesMaterialUnit        
          FieldReq%MaterialUnitName = FieldReqMaterialUnit%Name
        Else
          FieldReqMaterialUnit = MaterialUnits%MaterialUnits(                       &
                                   FindMaterialUnitIndex(FieldReq%MaterialUnitName, &
                                     MaterialUnits                                  &
                                   )                                                &
                                 )
        End If
        ! Check whether units are identical
        If ( SpeciesMaterialUnit .ne. FieldReqMaterialUnit ) Then
          Call Message( 'FATAL ERROR: Unit "'               // &
                        Trim(FieldReqMaterialUnit%Name)     // &
                        '" in Field Request "'              // &
                        Trim(FieldReq%Name)                 // &
                        '" does not match unit "'           // &
                        Trim(SpeciesMaterialUnit%Name)      // &
                        '" of species "'                    // &
                        Trim(FieldReq%SpeciesName)          // &
                        '".', 3                                &
               )
        End If
      End If
    End If

    
    !------------------------------------------------------------------!
    ! Create m agl coord system for boundary layer average quantities. !
    !------------------------------------------------------------------!

    If (FieldReq%AvBL) Then
      Call AddZCoord(ZCoord_m_agl(), Coords)
    End If

    ! $$ Prevent probs and %iles < 100 for Name II format.

    !-------------------------------------------------------------------------------------------------------!
    ! Add extra field requirements needed by the existing field requirements, together with any extra coord !
    ! systems, grids or processes needed.                                                                   !
    !-------------------------------------------------------------------------------------------------------!

    ! Convention. In the calls to InitFieldReq, if the value needed for an argument is known and equal to the
    ! corresponding element of FieldReq, the argument is given as FieldReq%... instead of giving the value
    ! explicitly. Hence the values that are given explicitly indicate arguments that may need to be different
    ! from the corresponding element of FieldReq.

    ! Restriction. Although extra processes can be introduced at this stage, they mustn't extend before the
    ! start of the original time range for F, O, Q, R quantities or after the end of the original time range
    ! for any quantities.

    ! $$ Could check this restriction ? Is it right? (e.g. think it should read '...extra processes and
    ! reqs ...
    ! if the req requesting the introduction is a type F, O, Q or R quantity ...). Check also v comments near
    ! start of module (also marked with $).
    ! $$ need to think also about what time levels are stored. This implies (?) mustn't go outside range
    ! except
    ! by intrinsic averaging.

    ! Process. Note Process is always assigned to allow use of Process%Code in if tests.
    If (FieldReq%iLastProcess > 0) Then
      Process => Reqs%Processes(iProcesses(FieldReq%iLastProcess))
    Else
      Process => Reqs%Processes(1)
    End If

    ! If only intrinsic averaging/integrating to do, work out if intrinsic averaging/integrating regions
    ! overlap.
    IntrinsicAvRegionsOverlap = .false.
    If (                                                                                  &
      FieldReq%iLastProcess == 0                                                    .and. &
      FieldReq%iPQRProc     == 0                                                    .and. &
      (FieldReq%IntrinsicTAv .or. FieldReq%IntrinsicHAv .or. FieldReq%IntrinsicZAv)       &
    ) Then
      ! $$ need to consider setting to true only if remove Avdt | Grid dT restriction
      ! or for semi-infinite intrinsic av up to T. Later case needs care - especially for fixed met.
      ! could have two reqs - integral to first time and t-grid integral.Is first time always finite?
    End If

    ! Percentiles. $$ combine with next if block??
    If (FieldReq%iLastProcess > 0 .and. Process%Code == P_Percent) Then

      ! $$ note SpeciesMaterialUnitName isn't really the right name here for the mixing ratio cases.
      If (FieldReq%SpeciesName .CIEq. ' ') Then 
        SpeciesMaterialUnitName  = ' '
      Else If (FieldReq%Quantity .CIEq. 'Mixing Ratio') Then
        SpeciesMaterialUnitName = FieldReq%MaterialUnitName
      Else If (FieldReq%Quantity .CIEq. 'E Mixing Ratio') Then
        SpeciesMaterialUnitName = FieldReq%MaterialUnitName
      Else
        iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
        SpeciesMaterialUnitName = Specieses%Specieses(iSpecies)%MaterialUnitName
      End If

      ExtraFieldReq = InitFieldReq(                                                           &
                        Name            = ' ',                                                &
                        Quantity        = FieldReq%Quantity,                                  &
                        SpeciesName     = FieldReq%SpeciesName,                               &
                        SourceName      = FieldReq%SourceName,                                &
                        SourceGroupName = FieldReq%SourceGroupName,                           &
                        SizeDistName    = FieldReq%SizeDistName,                              &
                        DecayDep        = FieldReq%DecayDep,                                  &
                        SemiInfApprox   = FieldReq%SemiInfApprox,                             &
                        TGridName       = FieldReq%TGridName,                                 &
                        SGridName       = FieldReq%SGridName,                                 &
                        HGridName       = FieldReq%HGridName,                                 &
                        ZGridName       = FieldReq%ZGridName,                                 &
                        HCoordName      = FieldReq%HCoordName,                                &
                        ZCoordName      = FieldReq%ZCoordName,                                &
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%iLastProcess - 1), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                              &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                              &
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                              &
                        Fluctuations    = FieldReq%Fluctuations,                              &
                        FlAvT           = FieldReq%FlAvT,                                     &
                        FlAvX           = FieldReq%FlAvX,                                     &
                        FlAvY           = FieldReq%FlAvY,                                     &
                        FlAvZ           = FieldReq%FlAvZ,                                     &
                        UseFlAvT        = FieldReq%UseFlAvT,                                  &
                        UseFlAvH        = FieldReq%UseFlAvH,                                  &
                        UseFlAvZ        = FieldReq%UseFlAvZ,                                  &
                        Sync            = FieldReq%Sync,                                      &
                        MaterialUnits   = MaterialUnits,                                      &
                        MaterialUnitName = SpeciesMaterialUnitName,                           &
                        BlockKey        = ' '                                                 &
                      )

      FieldReq%iProc = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iProc = FindFieldReqIndex(FieldReqName, Reqs)
      End If

    ! Averaging/integrating, probabilities, max/min and moments.
    Else If (FieldReq%iLastProcess > 0 .and. Process%Code /= P_Percent) Then

      If (Process%nT /= 1) Then

        If (FieldReq%TGridName /= ' ') Then

          TGrid => Grids%TGrids(iTGrid)

          T0 = TInTGrid(TGrid, 1) - (Process%sdT * (Process%nT - 1))
          TGridName = Trim(TGrid%Name)                                // &
                      Trim(Std2Char(ShortTime2RealTime(Process%sdT))) // &
                      Trim(Int2Char(Process%nT))

          Call AddTGrid(                                                      &
                 InitTGrid(                                                   &
                   Name = TGridName,                                          &
                   nT   = (TInTGrid(TGrid, TGrid%nT) - T0) / Process%sdT + 1, &
                   dT   = Process%dT,                                         &
                   T0   = ShortTime2Time(T0)                                  &
                 ),                                                           &
                 Grids                                                        &
               ) ! $$ better to be able to ask AddTGrid to assign and return name
               ! $$ Note will need to set correct grid averaging region if useT/H/ZGrid or BL (and av int)
               ! $$ For Bl (nZ > 1), need derived grid defined in terms of bl depth
               ! $$ Note CalcFlowResults etc should be able to cope with BL av with nZ = 1

        Else If (.not.IsInfFuture(Process%dT)) Then

        Else

          TGridName = ' '

        End If

      Else

        TGridName = FieldReq%TGridName

      End If

      ! $$ X, Y and Z grids

      ! If last processing step but still intrinsic averaging/integrating to do, create new processing step to
      ! do just the intrinsic averaging/integrating.
      If (                                                                                  &
        FieldReq%iLastProcess == 1                                                    .and. &
        (FieldReq%IntrinsicTAv .or. FieldReq%IntrinsicHAv .or. FieldReq%IntrinsicZAv)       &
      ) Then

        If (FieldReq%IntrinsicTAv) Then
          ExtraProcess%TAvInt   = Process%TAvInt
          ExtraProcess%dT       = Process%dT
          ExtraProcess%UseTGrid = Process%UseTGrid
        Else
          ExtraProcess%TAvInt   = P_Av
          ExtraProcess%dT       = ZeroTime()
          ExtraProcess%UseTGrid = .false.
        End If
        If (FieldReq%IntrinsicHAv) Then
          ExtraProcess%HAvInt   = Process%HAvInt
          ExtraProcess%dX       = Process%dX
          ExtraProcess%dY       = Process%dY
          ExtraProcess%UseHGrid = Process%UseHGrid
        Else
          ExtraProcess%HAvInt   = P_Av
          ExtraProcess%dX       = 0.0
          ExtraProcess%dY       = 0.0
          ExtraProcess%UseHGrid = .false.
        End If
        If (FieldReq%IntrinsicZAv) Then
          ExtraProcess%ZAvInt   = Process%ZAvInt
          ExtraProcess%dZ       = Process%dZ
          ExtraProcess%UseZGrid = Process%UseZGrid .or. (Process%BL .and. Process%nZ /= 1)
          ExtraProcess%BL       = Process%BL .and. Process%nZ == 1
        Else
          ExtraProcess%ZAvInt   = P_Av
          ExtraProcess%dZ       = 0.0
          ExtraProcess%UseZGrid = .false.
          ExtraProcess%BL       = .false.
        End If

        Call AddProcess(                            &
               InitProcess(                         &
                 Name      = ' ',                   &
                 Code      = P_AvInt,               &
                 Ensemble  = .false.,               &
                 TAvInt    = ExtraProcess%TAvInt,   &
                 nT        = 1,                     &
                 T         = ExtraProcess%dT,       &
                 dT        = ExtraProcess%dT,       &
                 UseTGrid  = ExtraProcess%UseTGrid, &
                 T0        = ReferenceTime(),       &
                 HAvInt    = ExtraProcess%HAvInt,   &
                 nX        = 1,                     &
                 nY        = 1,                     &
                 X         = ExtraProcess%dX,       &
                 Y         = ExtraProcess%dY,       &
                 dX        = ExtraProcess%dX,       &
                 dY        = ExtraProcess%dY,       &
                 UseHGrid  = ExtraProcess%UseHGrid, &
                 ZAvInt    = ExtraProcess%ZAvInt,   &
                 nZ        = 1,                     &
                 Z         = ExtraProcess%dZ,       &
                 dZ        = ExtraProcess%dZ,       &
                 UseZGrid  = ExtraProcess%UseZGrid, &
                 BL        = ExtraProcess%BL,       &
                 DGridName = ' ',                   &
                 BlockKey  = ' ',                   &
                 Item      = ' '                    &
               ),                                   &
               Reqs,                                &
               Name = ProcessName                   &
             )

        ! $$ note SpeciesMaterialUnitName isn't really the right name here for the mixing ratio cases.
        If (FieldReq%SpeciesName .CIEq. ' ') Then 
          SpeciesMaterialUnitName  = ' '
        Else If (FieldReq%Quantity .CIEq. 'Mixing Ratio') Then
          SpeciesMaterialUnitName = FieldReq%MaterialUnitName
        Else If (FieldReq%Quantity .CIEq. 'E Mixing Ratio') Then
          SpeciesMaterialUnitName = FieldReq%MaterialUnitName
        Else
          iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
          SpeciesMaterialUnitName = Specieses%Specieses(iSpecies)%MaterialUnitName
        End If

        ExtraFieldReq = InitFieldReq(                                 &
                          Name            = ' ',                      &
                          Quantity        = FieldReq%Quantity,        &
                          SpeciesName     = FieldReq%SpeciesName,     &
                          SourceName      = FieldReq%SourceName,      &
                          SourceGroupName = FieldReq%SourceGroupName, &
                          SizeDistName    = FieldReq%SizeDistName,    &
                          DecayDep        = FieldReq%DecayDep,        &
                          SemiInfApprox   = FieldReq%SemiInfApprox,   &
                          TGridName       = TGridName,                &
                          SGridName       = FieldReq%SGridName,       &
                          HGridName       = FieldReq%HGridName,       &
                          ZGridName       = FieldReq%ZGridName,       &
                          HCoordName      = FieldReq%HCoordName,      &
                          ZCoordName      = FieldReq%ZCoordName,      &
                          ProcessNames    = (/ ProcessName /),        &
                          IntrinsicTAv    = FieldReq%IntrinsicTAv,    &
                          IntrinsicHAv    = FieldReq%IntrinsicHAv,    &
                          IntrinsicZAv    = FieldReq%IntrinsicZAv,    &
                          Fluctuations    = FieldReq%Fluctuations,    &
                          FlAvT           = FieldReq%FlAvT,           &
                          FlAvX           = FieldReq%FlAvX,           &
                          FlAvY           = FieldReq%FlAvY,           &
                          FlAvZ           = FieldReq%FlAvZ,           &
                          UseFlAvT        = FieldReq%UseFlAvT,        &
                          UseFlAvH        = FieldReq%UseFlAvH,        &
                          UseFlAvZ        = FieldReq%UseFlAvZ,        &
                          Sync            = FieldReq%Sync,            &
                          MaterialUnits   = MaterialUnits,            &
                          MaterialUnitName = SpeciesMaterialUnitName, &
                          BlockKey        = ' '                       &
                        )

        FieldReq%iProc = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
        If (Error) Then
          Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
          FieldReq%iProc = FindFieldReqIndex(FieldReqName, Reqs)
        End If

      ! Otherwise.
      Else

        ! $$ note SpeciesMaterialUnitName isn't really the right name here for the mixing ratio cases.
        If (FieldReq%SpeciesName .CIEq. ' ') Then 
          SpeciesMaterialUnitName  = ' '
        Else If (FieldReq%Quantity .CIEq. 'Mixing Ratio') Then
          SpeciesMaterialUnitName = FieldReq%MaterialUnitName
        Else If (FieldReq%Quantity .CIEq. 'E Mixing Ratio') Then
          SpeciesMaterialUnitName = FieldReq%MaterialUnitName
        Else
          iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
          SpeciesMaterialUnitName = Specieses%Specieses(iSpecies)%MaterialUnitName
        End If

        ExtraFieldReq = InitFieldReq(                                                           &
                          Name            = ' ',                                                &
                          Quantity        = FieldReq%Quantity,                                  &
                          SpeciesName     = FieldReq%SpeciesName,                               &
                          SourceName      = FieldReq%SourceName,                                &
                          SourceGroupName = FieldReq%SourceGroupName,                           &
                          SizeDistName    = FieldReq%SizeDistName,                              &
                          DecayDep        = FieldReq%DecayDep,                                  &
                          SemiInfApprox   = FieldReq%SemiInfApprox,                             &
                          TGridName       = TGridName,                                          &
                          SGridName       = FieldReq%SGridName,                                 &
                          HGridName       = FieldReq%HGridName,                                 &
                          ZGridName       = FieldReq%ZGridName,                                 &
                          HCoordName      = FieldReq%HCoordName,                                &
                          ZCoordName      = FieldReq%ZCoordName,                                &
                          ProcessNames    = FieldReq%ProcessNames(1:FieldReq%iLastProcess - 1), &
                          IntrinsicTAv    = FieldReq%IntrinsicTAv,                              &
                          IntrinsicHAv    = FieldReq%IntrinsicHAv,                              &
                          IntrinsicZAv    = FieldReq%IntrinsicZAv,                              &
                          Fluctuations    = FieldReq%Fluctuations,                              &
                          FlAvT           = FieldReq%FlAvT,                                     &
                          FlAvX           = FieldReq%FlAvX,                                     &
                          FlAvY           = FieldReq%FlAvY,                                     &
                          FlAvZ           = FieldReq%FlAvZ,                                     &
                          UseFlAvT        = FieldReq%UseFlAvT,                                  &
                          UseFlAvH        = FieldReq%UseFlAvH,                                  &
                          UseFlAvZ        = FieldReq%UseFlAvZ,                                  &
                          Sync            = FieldReq%Sync,                                      &
                          MaterialUnits   = MaterialUnits,                                      &
                          MaterialUnitName = SpeciesMaterialUnitName,                           &
                          BlockKey        = ' '                                                 &
                        )

        FieldReq%iProc = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
        If (Error) Then
          Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
          FieldReq%iProc = FindFieldReqIndex(FieldReqName, Reqs)
        End If

      End If

      ! $$ If FieldReq%iLastProcess = iFluctuationProc { add fluc params
      !                                                { set fluctuation false except in fluc params

    ! Intrinsic averaging/integrating regions overlap.
    Else If (IntrinsicAvRegionsOverlap) Then

  ! $$ Here there will be a only single process with only intrinsic av/int
  ! $$ Only relevant if relax Avdt | Grid dT restriction or for semi-infinite intrinsic av up to T
  ! $$ Set up grids and process needed.
  ! $$ Set NoProc of FieldReq to false since this will now be processed from another field. Also may be
  !    useful to add an IntrinsicAvRegionsOverlap flag to fieldreq to indicate nature of processing needed.
  !    ExtraFieldReq = InitFieldReq(                                 &
  !                      Name            = ' ',                      &
  !                      Quantity        = FieldReq%Quantity,        &
  !                      SpeciesName     = FieldReq%SpeciesName,     &
  !                      SourceName      = FieldReq%SourceName,      &
  !                      SourceGroupName = FieldReq%SourceGroupName, &
  !                      SizeDistName    = FieldReq%SizeDistName,    &
  !                      DecayDep        = FieldReq%DecayDep,        &
  !                      SemiInfApprox   = FieldReq%SemiInfApprox,   &
  !                      TGridName       = TGridName,                &
  !                      SGridName       = FieldReq%SGridName,       &
  !                      HGridName       = FieldReq%HGridName,       &
  !                      ZGridName       = FieldReq%ZGridName,       &
  !                      HCoordName      = FieldReq%HCoordName,      &
  !                      ZCoordName      = FieldReq%ZCoordName,      &
  !                      ProcessNames    = (/ ProcessName /),        &
  !                      IntrinsicTAv    = FieldReq%IntrinsicTAv,    &
  !                      IntrinsicHAv    = FieldReq%IntrinsicHAv,    &
  !                      IntrinsicZAv    = FieldReq%IntrinsicZAv,    &
  !                      Fluctuations    = FieldReq%Fluctuations,    &
  !                      FlAvT           = FieldReq%FlAvT,           &
  !                      FlAvX           = FieldReq%FlAvX,           &
  !                      FlAvY           = FieldReq%FlAvY,           &
  !                      FlAvZ           = FieldReq%FlAvZ,           &
  !                      UseFlAvT        = FieldReq%UseFlAvT,        &
  !                      UseFlAvH        = FieldReq%UseFlAvH,        &
  !                      UseFlAvZ        = FieldReq%UseFlAvZ,        &
  !                      Sync            = FieldReq%Sync,            &
  !                      BlockKey        = ' '                       &
  !                    )
  !    FieldReq%iProc = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
  !    If (Error) Then
  !      Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
  !      FieldReq%iProc = FindFieldReqIndex(FieldReqName, Reqs) ! $$ note fieldreq will say no explicit
  !                                                             ! processing
  !                                                             ! so setting iProc is a bit inconsistent.
  !                                                             ! Use another
  !                                                             ! variable?
  !    End If

    ! Type P, Q, R fields - fluctuation parameters. Note Fluctuations turned off in any field requirements
    ! needed to calculate the fluctuation parameters. Intrinsic averaging is turned off for type E, F or O
    ! requirements.

    ! Sigma C - requires Air Concentration, Mean Travel Time and X Stats.
    Else If (FieldReq%iQuantity == Q_SigmaC) Then

      iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
      SpeciesMaterialUnitName  = Specieses%Specieses(iSpecies)%MaterialUnitName
      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Air Concentration',                          &
                        SpeciesName     = FieldReq%SpeciesName,                         &
                        SourceName      = FieldReq%SourceName,                          &
                        SourceGroupName = FieldReq%SourceGroupName,                     &
                        SizeDistName    = FieldReq%SizeDistName,                        &
                        DecayDep        = FieldReq%DecayDep,                            & ! = false
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       & ! = false
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           & ! = ' '
                        HGridName       = FieldReq%HGridName,                           &
                        ZGridName       = FieldReq%ZGridName,                           &
                        HCoordName      = FieldReq%HCoordName,                          & ! = ' '
                        ZCoordName      = FieldReq%ZCoordName,                          & ! = ' '
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        &
                        Fluctuations    = .false.,                                      &
                        FlAvT           = ZeroShortTime(),                              &
                        FlAvX           = 0.0,                                          &
                        FlAvY           = 0.0,                                          &
                        FlAvZ           = 0.0,                                          &
                        UseFlAvT        = .false.,                                      &
                        UseFlAvH        = .false.,                                      &
                        UseFlAvZ        = .false.,                                      &
                        Sync            = FieldReq%Sync,                                &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = SpeciesMaterialUnitName,                     &
                        BlockKey        = ' '                                           &
                      )
      FieldReq%iAirConc = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iAirConc = FindFieldReqIndex(FieldReqName, Reqs)
      End If

      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Mean Travel Time',                           &
                        SpeciesName     = FieldReq%SpeciesName,                         &
                        SourceName      = FieldReq%SourceName,                          &
                        SourceGroupName = FieldReq%SourceGroupName,                     &
                        SizeDistName    = FieldReq%SizeDistName,                        &
                        DecayDep        = FieldReq%DecayDep,                            & ! = false
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       & ! = false
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           & ! = ' '
                        HGridName       = FieldReq%HGridName,                           &
                        ZGridName       = FieldReq%ZGridName,                           &
                        HCoordName      = FieldReq%HCoordName,                          & ! = ' '
                        ZCoordName      = FieldReq%ZCoordName,                          & ! = ' '
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        &
                        Fluctuations    = .false.,                                      &
                        FlAvT           = ZeroShortTime(),                              &
                        FlAvX           = 0.0,                                          &
                        FlAvY           = 0.0,                                          &
                        FlAvZ           = 0.0,                                          &
                        UseFlAvT        = .false.,                                      &
                        UseFlAvH        = .false.,                                      &
                        UseFlAvZ        = .false.,                                      &
                        Sync            = FieldReq%Sync,                                &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = SpeciesMaterialUnitName,                     &
                        BlockKey        = ' '                                           &
                      )
      FieldReq%iMeanS = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iMeanS = FindFieldReqIndex(FieldReqName, Reqs)
      End If

      ! Add stereographic cartesian coord system at source location for X Stats.
      ! $$ Treatment of Multiple source case needs clarifying here.
      ! $$ get AddHCoord to return name
      Temp_1 = (/ 0.0_Std, 0.0_Std /)
      Temp_2 = (/ 1.0_Std, 1.0_Std /)
      Call AddHCoord(                                                                               &
             InitHCoord(                                                                            &
               Name        = 'Cartesian System at Source ' // FieldReq%SourceName,                  &
               CoordType   = H_PSCartesian,                                                         &
               Pole        = SourcePositionInStandardLatLong(FieldReq%SourceName, Sources, Coords), &
               Angle       = 0.0_Std,                                                               &
               Origin      = Temp_1,                                                                &
               Unit        = Temp_2,                                                                &
               ThetaOrigin = 0.0_Std,                                                               &
               ScaleFactor = 1.0_Std                                                                &
             ),                                                                                     &
             Coords                                                                                 &
           )
      ! Add m agl coord system for X Stats.
      Call AddZCoord(ZCoord_m_agl(), Coords)
      ! Add travel-time grid (S-grid) for X Stats.
      ! $$ Currently set up to use the same S-grid for all X Stats.
      ! $$ Also the choice of grid could be improved (e.g. variable spacing).
      ! $$ get AddTGrid to return name
      Call AddTGrid(                                            &
             InitTGrid(                                         &
               Name = 'S Grid for X Stats',                     &
               nT   = 40,                                       &
               dT   = Char2Time('00:00:05', Interval = .true.), &
               T0   = ZeroTime()                                &
             ),                                                 &
             Grids                                              &
           )
      ! Add process for T-averaging of X Stats.
      Call AddProcess(                      &
             InitProcess(                   &
               Name      = ' ',             &
               Code      = P_AvInt,         &
               Ensemble  = .false.,         &
               TAvInt    = P_Av,            &
               nT        = 1,               &
               T         = ZeroTime(),      & !}
               dT        = ZeroTime(),      & !}
               UseTGrid  = .true.,          & !} $$ control over av T? possibly relate to any
                                              !} $$ t-averaging in FieldReq?
               T0        = ReferenceTime(), &
               HAvInt    = P_Av,            &
               nX        = 1,               &
               nY        = 1,               &
               X         = 0.0,             &
               Y         = 0.0,             &
               dX        = 0.0,             &
               dY        = 0.0,             &
               UseHGrid  = .false.,         &
               ZAvInt    = P_Av,            &
               nZ        = 1,               &
               Z         = 0.0,             &
               dZ        = 0.0,             &
               UseZGrid  = .false.,         &
               BL        = .false.,         &
               DGridName = ' ',             &
               BlockKey  = ' ',             &
               Item      = ' '              &
             ),                             &
             Reqs,                          &
             Name = ProcessName             &
           )
      ExtraFieldReq = InitFieldReq(                                         &
                        Name            = ' ',                              &
                        Quantity        = 'X Stats',                        &
                        SpeciesName     = FieldReq%SpeciesName,             &
                        SourceName      = FieldReq%SourceName,              &
                        SourceGroupName = FieldReq%SourceGroupName,         &
                        SizeDistName    = FieldReq%SizeDistName,            &
                        DecayDep        = FieldReq%DecayDep,                & ! = false
                        SemiInfApprox   = FieldReq%SemiInfApprox,           & ! = false
                        TGridName       = FieldReq%TGridName,               &
                        SGridName       = 'S Grid for X Stats',             & ! $$
                        HGridName       = ' ',                              &
                        ZGridName       = ' ',                              &
                        HCoordName      = 'Cartesian System at Source ' //  & ! $$
                                          FieldReq%SourceName,              &
                        ZCoordName      = 'm agl',                          & ! $$
                        ProcessNames    = (/ ProcessName /),                &
                        IntrinsicTAv    = .true.,                           &
                        IntrinsicHAv    = .false.,                          &
                        IntrinsicZAv    = .false.,                          &
                        Fluctuations    = .false.,                          &
                        FlAvT           = ZeroShortTime(),                  &
                        FlAvX           = 0.0,                              &
                        FlAvY           = 0.0,                              &
                        FlAvZ           = 0.0,                              &
                        UseFlAvT        = .false.,                          &
                        UseFlAvH        = .false.,                          &
                        UseFlAvZ        = .false.,                          &
                        Sync            = FieldReq%Sync,                    &
                        MaterialUnits   = MaterialUnits,                    &
                        MaterialUnitName = SpeciesMaterialUnitName,         &
                        BlockKey        = ' '                               &
                      )
      FieldReq%iXStats = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iXStats = FindFieldReqIndex(FieldReqName, Reqs)
      End If

    ! Type P, Q, R fields - other than fluctuation parameters. Intrinsic averaging is turned off for type E, F
    ! or O requirements.

    ! Mixing Ratio - requires Air Concentration, Temperature and pressure.
    Else If (FieldReq%iQuantity == Q_MixingRatio) Then

      ! Check that both horizontal and vertical grids are specified
      If (FieldReq%HGridName == ' ') Then
        Call Message('ERROR: in setting up field request '  // &
                     Trim(FieldReq%Name)                    // &
                     ' horizontal grid has to be specified' // &
                     ' for quantity MixingRatio ',             &
                     3                                         &
             )
      End If
      If (FieldReq%ZGridName == ' ') Then
        Call Message('ERROR: in setting up field request ' // &
                     Trim(FieldReq%Name)                   // &
                     ' vertical grid has to be specified'  // &
                     ' for quantity MixingRatio ',            &
                     3                                        &
             )
      End If
    
      iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
      SpeciesMaterialUnitName  = Specieses%Specieses(iSpecies)%MaterialUnitName
      ! Check that Molecular Weight is specified for species and is a positive number.
      If (Specieses%Specieses(iSpecies)%MolecularWeight <= 0.0) Then
        Call Message ("ERROR in setting up Field Requirement '"       // &
                      Trim(FieldReq%Name) // "':"                     // &
                      " Molecular weight of species '"                // &
                      Trim(Specieses%Specieses(iSpecies)%Name)        // &
                      "' is not specified or not a positive number.",    &
                      3                                                  &
             )
      End If
      
      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Air Concentration',                          &
                        SpeciesName     = FieldReq%SpeciesName,                         &
                        SourceName      = FieldReq%SourceName,                          &
                        SourceGroupName = FieldReq%SourceGroupName,                     &
                        SizeDistName    = FieldReq%SizeDistName,                        &
                        DecayDep        = FieldReq%DecayDep,                            &
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       &
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           &
                        HGridName       = FieldReq%HGridName,                           &
                        ZGridName       = FieldReq%ZGridName,                           &
                        HCoordName      = FieldReq%HCoordName,                          &
                        ZCoordName      = FieldReq%ZCoordName,                          &
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        &
                        Fluctuations    = FieldReq%Fluctuations,                        &
                        FlAvT           = FieldReq%FlAvT,                               &
                        FlAvX           = FieldReq%FlAvX,                               &
                        FlAvY           = FieldReq%FlAvY,                               &
                        FlAvZ           = FieldReq%FlAvZ,                               &
                        UseFlAvT        = FieldReq%UseFlAvT,                            &
                        UseFlAvH        = FieldReq%UseFlAvH,                            &
                        UseFlAvZ        = FieldReq%UseFlAvZ,                            &
                        Sync            = FieldReq%Sync,                                &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = SpeciesMaterialUnitName,                     &
                        BlockKey        = ' '                                           &
                      )
      FieldReq%iAirConc = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iAirConc = FindFieldReqIndex(FieldReqName, Reqs)
      End If

      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Density',                                    &
                        SpeciesName     = ' ',                                          &
                        SourceName      = ' ',                                          &
                        SourceGroupName = ' ',                                          &
                        SizeDistName    = ' ',                                          &
                        DecayDep        = FieldReq%DecayDep,                            &
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       &
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           &
                        HGridName       = FieldReq%HGridName,                           &
                        ZGridName       = FieldReq%ZGridName,                           &
                        HCoordName      = FieldReq%HCoordName,                          &
                        ZCoordName      = FieldReq%ZCoordName,                          &
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        &
                        Fluctuations    = FieldReq%Fluctuations,                        &
                        FlAvT           = FieldReq%FlAvT,                               &
                        FlAvX           = FieldReq%FlAvX,                               &
                        FlAvY           = FieldReq%FlAvY,                               &
                        FlAvZ           = FieldReq%FlAvZ,                               &
                        UseFlAvT        = FieldReq%UseFlAvT,                            &
                        UseFlAvH        = FieldReq%UseFlAvH,                            &
                        UseFlAvZ        = FieldReq%UseFlAvZ,                            &
                        Sync            = .True.,                                       &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = '',                                          &
                        BlockKey        = ' '                                           &
                      )
                      
      FieldReq%iDensity = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iDensity = FindFieldReqIndex(FieldReqName, Reqs)
      End If

    ! Eulerian Mixing Ratio - requires Chemistry Field, Temperature and pressure.
    Else If (FieldReq%iQuantity == Q_EMixingRatio) Then

      ! Check that both horizontal and vertical grids are specified
      If (FieldReq%HGridName == ' ') Then
        Call Message('ERROR: in setting up field request '  // &
                     Trim(FieldReq%Name)                    // &
                     ' horizontal grid has to be specified' // &
                     ' for quantity '                       // &
                     Trim(FieldReq%Quantity),                  &
                     3                                         &
             )
      End If
      If (FieldReq%ZGridName == ' ') Then
        Call Message('ERROR: in setting up field request ' // &
                     Trim(FieldReq%Name)                   // &
                     ' vertical grid has to be specified'  // &
                     ' for quantity '                      // &
                     Trim(FieldReq%Quantity),                 &
                     3                                        &
             )
      End If
    
      iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
      SpeciesMaterialUnitName  = Specieses%Specieses(iSpecies)%MaterialUnitName
      ! Check that Molecular Weight is specified for species and is a positive number.
      If (Specieses%Specieses(iSpecies)%MolecularWeight <= 0.0) Then
        Call Message ("ERROR in setting up Field Requirement '"       // &
                      Trim(FieldReq%Name) // "':"                     // &
                      " Molecular weight of species '"                // &
                      Trim(Specieses%Specieses(iSpecies)%Name)        // &
                      "' is not specified or not a positive number.",    &
                      3                                                  &
             )
      End If
     
      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Chemistry Field',                            &
                        SpeciesName     = FieldReq%SpeciesName,                         &
                        SourceName      = FieldReq%SourceName,                          &
                        SourceGroupName = FieldReq%SourceGroupName,                     &
                        SizeDistName    = FieldReq%SizeDistName,                        &
                        DecayDep        = FieldReq%DecayDep,                            &
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       &
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           &
                        HGridName       = FieldReq%HGridName,                           &
                        ZGridName       = FieldReq%ZGridName,                           &
                        HCoordName      = FieldReq%HCoordName,                          &
                        ZCoordName      = FieldReq%ZCoordName,                          &
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        &
                        Fluctuations    = FieldReq%Fluctuations,                        &
                        FlAvT           = FieldReq%FlAvT,                               &
                        FlAvX           = FieldReq%FlAvX,                               &
                        FlAvY           = FieldReq%FlAvY,                               &
                        FlAvZ           = FieldReq%FlAvZ,                               &
                        UseFlAvT        = FieldReq%UseFlAvT,                            &
                        UseFlAvH        = FieldReq%UseFlAvH,                            &
                        UseFlAvZ        = FieldReq%UseFlAvZ,                            &
                        Sync            = FieldReq%Sync,                                &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = SpeciesMaterialUnitName,                     &
                        BlockKey        = ' '                                           &
                      )
      FieldReq%iChemField = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iChemField = FindFieldReqIndex(FieldReqName, Reqs)
      End If

      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Density',                                    &
                        SpeciesName     = ' ',                                          &
                        SourceName      = ' ',                                          &
                        SourceGroupName = ' ',                                          &
                        SizeDistName    = ' ',                                          &
                        DecayDep        = FieldReq%DecayDep,                            &
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       &
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           &
                        HGridName       = FieldReq%HGridName,                           &
                        ZGridName       = FieldReq%ZGridName,                           &
                        HCoordName      = FieldReq%HCoordName,                          &
                        ZCoordName      = FieldReq%ZCoordName,                          &
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        &
                        Fluctuations    = FieldReq%Fluctuations,                        &
                        FlAvT           = FieldReq%FlAvT,                               &
                        FlAvX           = FieldReq%FlAvX,                               &
                        FlAvY           = FieldReq%FlAvY,                               &
                        FlAvZ           = FieldReq%FlAvZ,                               &
                        UseFlAvT        = FieldReq%UseFlAvT,                            &
                        UseFlAvH        = FieldReq%UseFlAvH,                            &
                        UseFlAvZ        = FieldReq%UseFlAvZ,                            &
                        Sync            = .True.,                                       &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = '',                                          &
                        BlockKey        = ' '                                           &
                      )
                      
      FieldReq%iDensity = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iDensity = FindFieldReqIndex(FieldReqName, Reqs)
      End If

    ! Adult Cloud Gamma Effective or Organ Dose - requires Photon Flux and/or Air Conc.
    Else If (                                           &
      FieldReq%iQuantity == Q_AduEffCloudGammaDose .or. &
      FieldReq%iQuantity == Q_AduLunCloudGammaDose .or. &
      FieldReq%iQuantity == Q_AduThyCloudGammaDose .or. &
      FieldReq%iQuantity == Q_AduBoSCloudGammaDose      &
    ) Then

      If (FieldReq%SemiInfApprox) Then

        iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
        SpeciesMaterialUnitName  = Specieses%Specieses(iSpecies)%MaterialUnitName
        ExtraFieldReq = InitFieldReq(                                                     &
                          Name            = ' ',                                          &
                          Quantity        = 'Air Concentration',                          &
                          SpeciesName     = FieldReq%SpeciesName,                         &
                          SourceName      = FieldReq%SourceName,                          &
                          SourceGroupName = FieldReq%SourceGroupName,                     &
                          SizeDistName    = FieldReq%SizeDistName,                        &
                          DecayDep        = FieldReq%DecayDep,                            & ! = false
                          SemiInfApprox   = .false.,                                      &
                          TGridName       = FieldReq%TGridName,                           &
                          SGridName       = FieldReq%SGridName,                           & ! = ' '
                          HGridName       = FieldReq%HGridName,                           &
                          ZGridName       = FieldReq%ZGridName,                           &
                          HCoordName      = FieldReq%HCoordName,                          & ! = ' '
                          ZCoordName      = FieldReq%ZCoordName,                          & ! = ' '
                          ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                          IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                          IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                          IntrinsicZAv    = FieldReq%IntrinsicZAv,                        &
                          Fluctuations    = FieldReq%Fluctuations,                        & ! = false
                          FlAvT           = FieldReq%FlAvT,                               &
                                                                         ! = ZeroShortTime()
                          FlAvX           = FieldReq%FlAvX,                               & ! = 0
                          FlAvY           = FieldReq%FlAvY,                               & ! = 0
                          FlAvZ           = FieldReq%FlAvZ,                               & ! = 0
                          UseFlAvT        = FieldReq%UseFlAvT,                            & ! = false
                          UseFlAvH        = FieldReq%UseFlAvH,                            & ! = false
                          UseFlAvZ        = FieldReq%UseFlAvZ,                            & ! = false
                          Sync            = FieldReq%Sync,                                &
                          MaterialUnits   = MaterialUnits,                                &
                          MaterialUnitName = SpeciesMaterialUnitName,                     &
                          BlockKey        = ' '                                           &
                        )
        FieldReq%iAirConc = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
        If (Error) Then
          Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
          FieldReq%iAirConc = FindFieldReqIndex(FieldReqName, Reqs)
        End If

      Else

        iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
        SpeciesMaterialUnitName = Specieses%Specieses(iSpecies)%MaterialUnitName
        ExtraFieldReq = InitFieldReq(                                                     &
                          Name            = ' ',                                          &
                          Quantity        = 'Photon Flux',                                &
                          SpeciesName     = FieldReq%SpeciesName,                         &
                          SourceName      = FieldReq%SourceName,                          &
                          SourceGroupName = FieldReq%SourceGroupName,                     &
                          SizeDistName    = FieldReq%SizeDistName,                        &
                          DecayDep        = FieldReq%DecayDep,                            & ! = false
                          SemiInfApprox   = FieldReq%SemiInfApprox,                       & ! = false
                          TGridName       = FieldReq%TGridName,                           &
                          SGridName       = FieldReq%SGridName,                           & ! = ' '
                          HGridName       = FieldReq%HGridName,                           &
                          ZGridName       = FieldReq%ZGridName,                           &
                          HCoordName      = FieldReq%HCoordName,                          & ! = ' '
                          ZCoordName      = FieldReq%ZCoordName,                          & ! = ' '
                          ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                          IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                          IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                          IntrinsicZAv    = FieldReq%IntrinsicZAv,                        &
                          Fluctuations    = FieldReq%Fluctuations,                        & ! = false
                          FlAvT           = FieldReq%FlAvT,                               &
                                                                                     ! = ZeroShortTime()
                          FlAvX           = FieldReq%FlAvX,                               & ! = 0
                          FlAvY           = FieldReq%FlAvY,                               & ! = 0
                          FlAvZ           = FieldReq%FlAvZ,                               & ! = 0
                          UseFlAvT        = FieldReq%UseFlAvT,                            & ! = false
                          UseFlAvH        = FieldReq%UseFlAvH,                            & ! = false
                          UseFlAvZ        = FieldReq%UseFlAvZ,                            & ! = false
                          Sync            = FieldReq%Sync,                                &
                          MaterialUnits   = MaterialUnits,                                &
                          MaterialUnitName = SpeciesMaterialUnitName,                     &
                          BlockKey        = ' '                                           &
                         )
        FieldReq%iPhotonFlux = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
        If (Error) Then
          Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
          FieldReq%iPhotonFlux = FindFieldReqIndex(FieldReqName, Reqs)
        End If

      End If

    ! Concentration - requires Air Concentration and Eulerian Concentration.
    Else If (FieldReq%iQuantity == Q_Concentration) Then
     
      iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
      SpeciesMaterialUnitName = Specieses%Specieses(iSpecies)%MaterialUnitName
      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Air Concentration',                          &
                        SpeciesName     = FieldReq%SpeciesName,                         &
                        SourceName      = FieldReq%SourceName,                          & ! = ' '
                        SourceGroupName = FieldReq%SourceGroupName,                     & ! = ' '
                        SizeDistName    = FieldReq%SizeDistName,                        &
                        DecayDep        = FieldReq%DecayDep,                            & ! = false
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       & ! = false
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           & ! = ' '
                        HGridName       = FieldReq%HGridName,                           &
                        ZGridName       = FieldReq%ZGridName,                           &
                        HCoordName      = FieldReq%HCoordName,                          & ! = ' '
                        ZCoordName      = FieldReq%ZCoordName,                          & ! = ' '
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        &
                        Fluctuations    = FieldReq%Fluctuations,                        & ! = false
                        FlAvT           = FieldReq%FlAvT,                               & ! = ZeroShortTime()
                        FlAvX           = FieldReq%FlAvX,                               & ! = 0
                        FlAvY           = FieldReq%FlAvY,                               & ! = 0
                        FlAvZ           = FieldReq%FlAvZ,                               & ! = 0
                        UseFlAvT        = FieldReq%UseFlAvT,                            & ! = false
                        UseFlAvH        = FieldReq%UseFlAvH,                            & ! = false
                        UseFlAvZ        = FieldReq%UseFlAvZ,                            & ! = false
                        Sync            = FieldReq%Sync,                                &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = SpeciesMaterialUnitName,                     &
                        BlockKey        = ' '                                           &
                      )
      FieldReq%iAirConc = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iAirConc = FindFieldReqIndex(FieldReqName, Reqs)
      End If

      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Eulerian Concentration',                     &
                        SpeciesName     = FieldReq%SpeciesName,                         &
                        SourceName      = FieldReq%SourceName,                          & ! = ' '
                        SourceGroupName = FieldReq%SourceGroupName,                     & ! = ' '
                        SizeDistName    = FieldReq%SizeDistName,                        &
                        DecayDep        = FieldReq%DecayDep,                            & ! = false
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       & ! = false
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           & ! = ' '
                        HGridName       = FieldReq%HGridName,                           &
                        ZGridName       = FieldReq%ZGridName,                           &
                        HCoordName      = FieldReq%HCoordName,                          & ! = ' '
                        ZCoordName      = FieldReq%ZCoordName,                          & ! = ' '
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        &
                        Fluctuations    = FieldReq%Fluctuations,                        & ! = false
                        FlAvT           = FieldReq%FlAvT,                               & ! = ZeroShortTime()
                        FlAvX           = FieldReq%FlAvX,                               & ! = 0
                        FlAvY           = FieldReq%FlAvY,                               & ! = 0
                        FlAvZ           = FieldReq%FlAvZ,                               & ! = 0
                        UseFlAvT        = FieldReq%UseFlAvT,                            & ! = false
                        UseFlAvH        = FieldReq%UseFlAvH,                            & ! = false
                        UseFlAvZ        = FieldReq%UseFlAvZ,                            & ! = false
                        Sync            = .true.,                                       & ! $$ logic of sync needs improving
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = SpeciesMaterialUnitName,                     &
                        BlockKey        = ' '                                           &
                      )
      FieldReq%iEulConc = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iEulConc = FindFieldReqIndex(FieldReqName, Reqs)
      End If
   
    ! Area At Risk - requires Air Concentration.
    Else If (FieldReq%iQuantity == Q_AreaAtRisk) Then

      iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
      SpeciesMaterialUnitName = Specieses%Specieses(iSpecies)%MaterialUnitName
      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Air Concentration',                          &
                        SpeciesName     = FieldReq%SpeciesName,                         &
                        SourceName      = FieldReq%SourceName,                          &
                        SourceGroupName = FieldReq%SourceGroupName,                     &
                        SizeDistName    = FieldReq%SizeDistName,                        &
                        DecayDep        = FieldReq%DecayDep,                            & ! = false
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       & ! = false
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           & ! = ' '
                        HGridName       = FieldReq%HGridName,                           &
                        ZGridName       = FieldReq%ZGridName,                           &
                        HCoordName      = FieldReq%HCoordName,                          & ! = ' '
                        ZCoordName      = FieldReq%ZCoordName,                          & ! = ' '
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        &
                        Fluctuations    = FieldReq%Fluctuations,                        & ! = false
                        FlAvT           = FieldReq%FlAvT,                               & ! = ZeroShortTime()
                        FlAvX           = FieldReq%FlAvX,                               & ! = 0
                        FlAvY           = FieldReq%FlAvY,                               & ! = 0
                        FlAvZ           = FieldReq%FlAvZ,                               & ! = 0
                        UseFlAvT        = FieldReq%UseFlAvT,                            & ! = false
                        UseFlAvH        = FieldReq%UseFlAvH,                            & ! = false
                        UseFlAvZ        = FieldReq%UseFlAvZ,                            & ! = false
                        Sync            = FieldReq%Sync,                                &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = SpeciesMaterialUnitName,                     &
                        BlockKey        = ' '                                           &
                      )
      FieldReq%iAirConc = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iAirConc = FindFieldReqIndex(FieldReqName, Reqs)
      End If

    ! Deposition Rate - requires Dry Deposition Rate and Wet Deposition Rate.
    Else If (FieldReq%iQuantity == Q_Dep) Then

      iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
      SpeciesMaterialUnitName = Specieses%Specieses(iSpecies)%MaterialUnitName
      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Dry Deposition Rate',                        &
                        SpeciesName     = FieldReq%SpeciesName,                         &
                        SourceName      = FieldReq%SourceName,                          &
                        SourceGroupName = FieldReq%SourceGroupName,                     &
                        SizeDistName    = FieldReq%SizeDistName,                        &
                        DecayDep        = FieldReq%DecayDep,                            &
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       & ! = false
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           & ! = ' '
                        HGridName       = FieldReq%HGridName,                           &
                        ZGridName       = FieldReq%ZGridName,                           & ! = ' '
                        HCoordName      = FieldReq%HCoordName,                          & ! = ' '
                        ZCoordName      = FieldReq%ZCoordName,                          & ! = ' '
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        & ! = ' '
                        Fluctuations    = FieldReq%Fluctuations,                        & ! = false
                        FlAvT           = FieldReq%FlAvT,                               & ! = ZeroShortTime()
                        FlAvX           = FieldReq%FlAvX,                               & ! = 0
                        FlAvY           = FieldReq%FlAvY,                               & ! = 0
                        FlAvZ           = FieldReq%FlAvZ,                               & ! = 0
                        UseFlAvT        = FieldReq%UseFlAvT,                            & ! = false
                        UseFlAvH        = FieldReq%UseFlAvH,                            & ! = false
                        UseFlAvZ        = FieldReq%UseFlAvZ,                            & ! = false
                        Sync            = FieldReq%Sync,                                &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = SpeciesMaterialUnitName,                     &
                        BlockKey        = ' '                                           &
                      )
      FieldReq%iDryDep = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iDryDep = FindFieldReqIndex(FieldReqName, Reqs)
      End If

      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Wet Deposition Rate',                        &
                        SpeciesName     = FieldReq%SpeciesName,                         &
                        SourceName      = FieldReq%SourceName,                          &
                        SourceGroupName = FieldReq%SourceGroupName,                     &
                        SizeDistName    = FieldReq%SizeDistName,                        &
                        DecayDep        = FieldReq%DecayDep,                            &
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       & ! = false
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           & ! = ' '
                        HGridName       = FieldReq%HGridName,                           &
                        ZGridName       = FieldReq%ZGridName,                           & ! = ' '
                        HCoordName      = FieldReq%HCoordName,                          & ! = ' '
                        ZCoordName      = FieldReq%ZCoordName,                          & ! = ' '
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        & ! = ' '
                        Fluctuations    = FieldReq%Fluctuations,                        & ! = false
                        FlAvT           = FieldReq%FlAvT,                               & ! = ZeroShortTime()
                        FlAvX           = FieldReq%FlAvX,                               & ! = 0
                        FlAvY           = FieldReq%FlAvY,                               & ! = 0
                        FlAvZ           = FieldReq%FlAvZ,                               & ! = 0
                        UseFlAvT        = FieldReq%UseFlAvT,                            & ! = false
                        UseFlAvH        = FieldReq%UseFlAvH,                            & ! = false
                        UseFlAvZ        = FieldReq%UseFlAvZ,                            & ! = false
                        Sync            = FieldReq%Sync,                                &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = SpeciesMaterialUnitName,                     &
                        BlockKey        = ' '                                           &
                      )
      FieldReq%iWetDep = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iWetDep = FindFieldReqIndex(FieldReqName, Reqs)
      End If

    ! Mean travel time - requires mean concentration.
    Else If (FieldReq%iQuantity == Q_MeanS) Then

      iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
      SpeciesMaterialUnitName = Specieses%Specieses(iSpecies)%MaterialUnitName
      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Air Concentration',                          &
                        SpeciesName     = FieldReq%SpeciesName,                         &
                        SourceName      = FieldReq%SourceName,                          &
                        SourceGroupName = FieldReq%SourceGroupName,                     &
                        SizeDistName    = FieldReq%SizeDistName,                        &
                        DecayDep        = FieldReq%DecayDep,                            & ! = false
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       & ! = false
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           & ! = ' '
                        HGridName       = FieldReq%HGridName,                           &
                        ZGridName       = FieldReq%ZGridName,                           &
                        HCoordName      = FieldReq%HCoordName,                          & ! = ' '
                        ZCoordName      = FieldReq%ZCoordName,                          & ! = ' '
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        &
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        &
                        Fluctuations    = FieldReq%Fluctuations,                        & ! = false
                        FlAvT           = FieldReq%FlAvT,                               & ! = ZeroShortTime()
                        FlAvX           = FieldReq%FlAvX,                               & ! = 0
                        FlAvY           = FieldReq%FlAvY,                               & ! = 0
                        FlAvZ           = FieldReq%FlAvZ,                               & ! = 0
                        UseFlAvT        = FieldReq%UseFlAvT,                            & ! = false
                        UseFlAvH        = FieldReq%UseFlAvH,                            & ! = false
                        UseFlAvZ        = FieldReq%UseFlAvZ,                            & ! = false
                        Sync            = FieldReq%Sync,                                &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = SpeciesMaterialUnitName,                     &
                        BlockKey        = ' '                                           &
                      )
      ! $$ need to ensure that same set of particles are considered as in travel time
      ! calculation (i.e. same mass). Best to allow grid box mass calc.
      FieldReq%iAirConc = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iAirConc = FindFieldReqIndex(FieldReqName, Reqs)
      End If

    ! Sigma Z - requires Mean Z and Mass.
    Else If (FieldReq%iQuantity == Q_SigmaZ) Then

      iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
      SpeciesMaterialUnitName = Specieses%Specieses(iSpecies)%MaterialUnitName
      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Mean Z',                                     &
                        SpeciesName     = FieldReq%SpeciesName,                         &
                        SourceName      = FieldReq%SourceName,                          &
                        SourceGroupName = FieldReq%SourceGroupName,                     &
                        SizeDistName    = FieldReq%SizeDistName,                        &
                        DecayDep        = FieldReq%DecayDep,                            & ! = false
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       & ! = false
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           & ! = ' '
                        HGridName       = FieldReq%HGridName,                           & ! = ' '
                        ZGridName       = FieldReq%ZGridName,                           & ! = ' '
                        HCoordName      = FieldReq%HCoordName,                          & ! = ' '
                        ZCoordName      = FieldReq%ZCoordName,                          &
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        & ! = false
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        & ! = false
                        Fluctuations    = FieldReq%Fluctuations,                        & ! = false
                        FlAvT           = FieldReq%FlAvT,                               & ! = ZeroShortTime()
                        FlAvX           = FieldReq%FlAvX,                               & ! = 0
                        FlAvY           = FieldReq%FlAvY,                               & ! = 0
                        FlAvZ           = FieldReq%FlAvZ,                               & ! = 0
                        UseFlAvT        = FieldReq%UseFlAvT,                            & ! = false
                        UseFlAvH        = FieldReq%UseFlAvH,                            & ! = false
                        UseFlAvZ        = FieldReq%UseFlAvZ,                            & ! = false
                        Sync            = FieldReq%Sync,                                &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = SpeciesMaterialUnitName,                     &
                        BlockKey        = ' '                                           &
                      )
      FieldReq%iMeanZ = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iMeanZ = FindFieldReqIndex(FieldReqName, Reqs)
      End If

      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Mass',                                       &
                        SpeciesName     = FieldReq%SpeciesName,                         &
                        SourceName      = FieldReq%SourceName,                          &
                        SourceGroupName = FieldReq%SourceGroupName,                     &
                        SizeDistName    = FieldReq%SizeDistName,                        &
                        DecayDep        = FieldReq%DecayDep,                            & ! = false
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       & ! = false
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           & ! = ' '
                        HGridName       = FieldReq%HGridName,                           & ! = ' '
                        ZGridName       = FieldReq%ZGridName,                           & ! = ' '
                        HCoordName      = FieldReq%HCoordName,                          & ! = ' '
                        ZCoordName      = ' ',                                          &
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        & ! = false
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        & ! = false
                        Fluctuations    = FieldReq%Fluctuations,                        & ! = false
                        FlAvT           = FieldReq%FlAvT,                               & ! = ZeroShortTime()
                        FlAvX           = FieldReq%FlAvX,                               & ! = 0
                        FlAvY           = FieldReq%FlAvY,                               & ! = 0
                        FlAvZ           = FieldReq%FlAvZ,                               & ! = 0
                        UseFlAvT        = FieldReq%UseFlAvT,                            & ! = false
                        UseFlAvH        = FieldReq%UseFlAvH,                            & ! = false
                        UseFlAvZ        = FieldReq%UseFlAvZ,                            & ! = false
                        Sync            = FieldReq%Sync,                                &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = SpeciesMaterialUnitName,                     &
                        BlockKey        = ' '                                           &
                      )
      FieldReq%iMass = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iMass = FindFieldReqIndex(FieldReqName, Reqs)
      End If

    ! Mean Z - requires Mass.
    Else If (FieldReq%iQuantity == Q_MeanZ) Then

      iSpecies = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
      SpeciesMaterialUnitName = Specieses%Specieses(iSpecies)%MaterialUnitName
      ExtraFieldReq = InitFieldReq(                                                     &
                        Name            = ' ',                                          &
                        Quantity        = 'Mass',                                       &
                        SpeciesName     = FieldReq%SpeciesName,                         &
                        SourceName      = FieldReq%SourceName,                          &
                        SourceGroupName = FieldReq%SourceGroupName,                     &
                        SizeDistName    = FieldReq%SizeDistName,                        &
                        DecayDep        = FieldReq%DecayDep,                            & ! = false
                        SemiInfApprox   = FieldReq%SemiInfApprox,                       & ! = false
                        TGridName       = FieldReq%TGridName,                           &
                        SGridName       = FieldReq%SGridName,                           & ! = ' '
                        HGridName       = FieldReq%HGridName,                           & ! = ' '
                        ZGridName       = FieldReq%ZGridName,                           & ! = ' '
                        HCoordName      = FieldReq%HCoordName,                          & ! = ' '
                        ZCoordName      = ' ',                                          &
                        ProcessNames    = FieldReq%ProcessNames(1:FieldReq%nProcesses), &
                        IntrinsicTAv    = FieldReq%IntrinsicTAv,                        &
                        IntrinsicHAv    = FieldReq%IntrinsicHAv,                        & ! = false
                        IntrinsicZAv    = FieldReq%IntrinsicZAv,                        & ! = false
                        Fluctuations    = FieldReq%Fluctuations,                        & ! = false
                        FlAvT           = FieldReq%FlAvT,                               & ! = ZeroShortTime()
                        FlAvX           = FieldReq%FlAvX,                               & ! = 0
                        FlAvY           = FieldReq%FlAvY,                               & ! = 0
                        FlAvZ           = FieldReq%FlAvZ,                               & ! = 0
                        UseFlAvT        = FieldReq%UseFlAvT,                            & ! = false
                        UseFlAvH        = FieldReq%UseFlAvH,                            & ! = false
                        UseFlAvZ        = FieldReq%UseFlAvZ,                            & ! = false
                        Sync            = FieldReq%Sync,                                &
                        MaterialUnits   = MaterialUnits,                                &
                        MaterialUnitName = SpeciesMaterialUnitName,                     &
                        BlockKey        = ' '                                           &
                      )
      FieldReq%iMass = FindEquivFieldReqIndex(ExtraFieldReq, Reqs, Error)
      If (Error) Then
        Call AddFieldReq(ExtraFieldReq, Reqs, FieldReqName)
        FieldReq%iMass = FindFieldReqIndex(FieldReqName, Reqs)
      End If

    End If

    i = i + 1

  End Do

  !-----------------------------------------------------------------------------------------------!
  ! Check for one field extending another, set iExt, and avoid chains of references between field !
  ! requirements.                                                                                 !
  !-----------------------------------------------------------------------------------------------!

  ! Set iExt = 0.
  Do i = 1, Reqs%nFieldReqs
    FieldReq => Reqs%FieldReqs(i)
    FieldReq%iExt = 0
  End Do

  ! Find field extensions.
  Do i = 1, Reqs%nFieldReqs

    FieldReq => Reqs%FieldReqs(i)

    Do j = 1, Reqs%nFieldReqs

      If (j == i) Cycle
      FieldReq1 => Reqs%FieldReqs(j)
      If (FieldReq1%iExt /= 0) Cycle

      If (FieldReq1 .equiv. FieldReq) Then ! $$ extend definition of extension to include grids?
                                           ! $$ define function IsAnExtension?
                                           ! $$ If A extends B in some ways and B extends A in other
                                           !    ways, define C which extends both. May want to mark req with
                                           !  'AllInMemory flag' as cant say .not.TSeparateFile .and. TAcross
                                           !    for non-output reqs.
        If (                                                               &
               (.not.FieldReq1%TSeparateFile .and. FieldReq1%TAcross) .or. &
          .not.(.not.FieldReq%TSeparateFile  .and. FieldReq%TAcross )      &
        ) Then
          FieldReq%iExt = j
          Exit
        End If
      End If

    End Do

  End Do

  ! Avoid chains of references between field requirements.
  Do

    Changes = .false.

    Do i = 1, Reqs%nFieldReqs

      FieldReq => Reqs%FieldReqs(i)

      If (FieldReq%iExt /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iExt)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iExt = FieldReq1%iExt
          Changes = .true.
        End If
      End If

      If (FieldReq%iProc /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iProc)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iProc = FieldReq1%iExt
          Changes = .true.
        End If
      End If

      If (FieldReq%iDryDep /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iDryDep)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iDryDep = FieldReq1%iExt
          Changes = .true.
        End If
      End If

      If (FieldReq%iWetDep /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iWetDep)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iWetDep = FieldReq1%iExt
          Changes = .true.
        End If
      End If

      If (FieldReq%iAirConc /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iAirConc)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iAirConc = FieldReq1%iExt
          Changes = .true.
        End If
      End If

      If (FieldReq%iEulConc /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iEulConc)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iEulConc = FieldReq1%iExt
          Changes = .true.
        End If
      End If

      If (FieldReq%iChemField /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iChemField)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iChemField = FieldReq1%iExt
          Changes = .true.
        End If
      End If

      If (FieldReq%iDensity /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iDensity)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iDensity = FieldReq1%iExt
          Changes = .true.
        End If
      End If

      If (FieldReq%iMeanS /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iMeanS)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iMeanS = FieldReq1%iExt
          Changes = .true.
        End If
      End If

      If (FieldReq%iXStats /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iXStats)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iXStats = FieldReq1%iExt
          Changes = .true.
        End If
      End If

      If (FieldReq%iSigmaC /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iSigmaC)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iSigmaC = FieldReq1%iExt
          Changes = .true.
        End If
      End If

      If (FieldReq%iMass /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iMass)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iMass = FieldReq1%iExt
          Changes = .true.
        End If
      End If

      If (FieldReq%iMeanZ /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iMeanZ)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iMeanZ = FieldReq1%iExt
          Changes = .true.
        End If
      End If

      If (FieldReq%iPhotonFlux /= 0) Then
        FieldReq1 => Reqs%FieldReqs(FieldReq%iPhotonFlux)
        If (FieldReq1%iExt /= 0) Then
          FieldReq%iPhotonFlux = FieldReq1%iExt
          Changes = .true.
        End If
      End If

    End Do

    If (.not.Changes) Exit

  End Do

End Subroutine SetUpReqs

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpiReqs(Coords, Grids, Specieses, Sources, SizeDists, Reqs, MaterialUnits)
! Sets up indices in Reqs.

  Implicit None
  ! Argument list:
  Type(Coords_),         Intent(In)            :: Coords         ! Collection of coord systems.
  Type(Grids_),          Intent(In)            :: Grids          ! Collection of grids.
  Type(Specieses_),      Intent(In)            :: Specieses      ! Collection of specieses.
  Type(Sources_),        Intent(In)            :: Sources        ! Collection of sources.
  Type(SizeDists_),      Intent(In)            :: SizeDists      ! Collection of particle size distributions.
  Type(Reqs_),           Intent(InOut), Target :: Reqs           ! Collection of requirements.
  Type(MaterialUnits_),  Intent(In)            :: MaterialUnits  ! Collection of material units.
  ! Locals:
  Type(Process_),   Pointer :: Process   !} Abbreviations for processing steps and requirements.
  Type(FieldReq_),  Pointer :: FieldReq  !}
  Type(PPInfoReq_), Pointer :: PPInfoReq !}
  Integer                   :: i         !] Loop indices.
  Integer                   :: j         !]

  ! Processes.
  Do i = 1, Reqs%nProcesses

    Process => Reqs%Processes(i)

    ! iDGrid.
    If (Process%DGridName == ' ') Then
      Process%iDGrid = 0
    Else
      Process%iDGrid = FindDGridIndex(Process%DGridName, Grids)
    End If

  End Do

  ! Fields.
  Do i = 1, Reqs%nFieldReqs

    FieldReq => Reqs%FieldReqs(i)

    ! iSpecies.
    If (FieldReq%SpeciesName == ' ' ) Then 
      FieldReq%iSpecies         = 0
      FieldReq%iParticleSpecies = 0
    Else
      FieldReq%iSpecies         = FindSpeciesIndex(FieldReq%SpeciesName, Specieses)
      FieldReq%iParticleSpecies = Specieses%iSpecies2Particle(FieldReq%iSpecies)
      If (FieldReq%LEFOPQRType == 'L' .and. FieldReq%iParticleSpecies == 0) Then
        Call Message(                                               &
               'WARNING: The output field requirement '          // &
               Trim(FieldReq%Name)                               // &
               ' involving the quantity '                        // &
               Trim(FieldReq%Quantity)                           // &
               ' involves a species '                            // &
               Trim(Specieses%Specieses(FieldReq%iSpecies)%Name) // &
               ' that is not carried on particles. '             // &
               'Hence the output will be zero.',                    &
               1                                                    &
             )
      End If
    End If
    Select Case (FieldReq%iQuantity)
      Case (Q_PhotonFlux, Q_AduEffCloudGammaDose, Q_AduLunCloudGammaDose, Q_AduThyCloudGammaDose, Q_AduBoSCloudGammaDose)
        If (Specieses%Specieses(FieldReq%iSpecies)%iCloudGammaParams == 0) Then
          Call Message(                                                        &
                 'WARNING: The output field requirement '                   // &
                 Trim(FieldReq%Name)                                        // &
                 ' involving the quantity '                                 // &
                 Trim(FieldReq%Quantity)                                    // &
                 ' involves a species '                                     // &
                 Trim(Specieses%Specieses(FieldReq%iSpecies)%Name)          // &
                 ' that does not have any cloud gamma parameters defined. ' // &
                 'Hence the output will be zero.',                             &
                 1                                                             &
               )
        End If
    End Select

    ! iMaterialUnit.
    If (FieldReq%MaterialUnitName .CIEq. ' ') Then
      FieldReq%iMaterialUnit = 0
    Else
      FieldReq%iMaterialUnit = FindMaterialUnitIndex(FieldReq%MaterialUnitName, MaterialUnits)
    End If

    ! iSource.
    If (FieldReq%SourceName == ' ') Then
      FieldReq%iSource = 0
    Else
      FieldReq%iSource = FindSourceIndex(FieldReq%SourceName, Sources)
    End If

    ! iSourceGroup.
    If (FieldReq%SourceGroupName == ' ') Then
      FieldReq%iSourceGroup = 0
    Else
      FieldReq%iSourceGroup = FindSourceGroupIndex(FieldReq%SourceGroupName, Sources)
    End If

    ! iSizeDist.
    If (FieldReq%SizeDistName == ' ') Then
      FieldReq%iSizeDist = 0
      FieldReq%nRanges   = 1
    Else
      FieldReq%iSizeDist = FindSizeDistIndex(FieldReq%SizeDistName, SizeDists)
      FieldReq%nRanges   = SizeDists%SizeDists(FieldReq%iSizeDist)%nSizeRanges
    End If

    ! iTGrid.
    If (FieldReq%TGridName == ' ') Then
      FieldReq%iTGrid = 0
    Else
      FieldReq%iTGrid = FindTGridIndex(FieldReq%TGridName, Grids)
    End If

    ! iSGrid.
    If (FieldReq%SGridName == ' ') Then
      FieldReq%iSGrid = 0
    Else
      FieldReq%iSGrid = FindTGridIndex(FieldReq%SGridName, Grids)
    End If

    ! iHCoord, iHGrid and iHGridCoord.
    If (FieldReq%HCoordName == ' ') Then
      FieldReq%iHCoord = 0
    Else
      FieldReq%iHCoord = FindHCoordIndex(FieldReq%HCoordName, Coords)
    End If
    If (FieldReq%HGridName == ' ') Then
      FieldReq%iHGrid      = 0
      FieldReq%iHGridCoord = 0
    Else
      FieldReq%iHGrid      = FindHGridIndex(FieldReq%HGridName, Grids)
      FieldReq%iHGridCoord = FindHCoordIndex(Grids%HGrids(FieldReq%iHGrid)%HCoordName, Coords)
    End If

    ! iZCoord, iZGrid and iZGridCoord.
    If (FieldReq%ZCoordName == ' ') Then
      FieldReq%iZCoord = 0
    Else
      FieldReq%iZCoord = FindZCoordIndex(FieldReq%ZCoordName, Coords)
    End If
    If (FieldReq%ZGridName == ' ') Then
      FieldReq%iZGrid = 0
      If (FieldReq%AvBL) Then
        FieldReq%iZGridCoord = FindZCoordIndex('m agl', Coords)
      Else
        FieldReq%iZGridCoord = 0
      End If
    Else
      FieldReq%iZGrid      = FindZGridIndex(FieldReq%ZGridName, Grids)
      FieldReq%iZGridCoord = FindZCoordIndex(Grids%ZGrids(FieldReq%iZGrid)%ZCoordName, Coords)
    End If

    ! iDGrid.
    If (FieldReq%DGridName == ' ') Then
      FieldReq%iDGrid = 0
    Else
      FieldReq%iDGrid = FindDGridIndex(FieldReq%DGridName, Grids)
    End If

    ! iProcesses.
    Do j = 1, FieldReq%nProcesses
      FieldReq%iProcesses(j) = FindProcessIndex(FieldReq%ProcessNames(j), Reqs)
    End Do

    ! iOutputGroup.
    FieldReq%iOutputGroup = 0
    If (FieldReq%OutputGroup /= ' ') Then
      Do j = 1, i - 1
        If (FieldReq%OutputGroup .CIEq. Reqs%FieldReqs(j)%OutputGroup) Then
          FieldReq%iOutputGroup = Reqs%FieldReqs(j)%iOutputGroup
          Exit
        End If
      End Do
      If (FieldReq%iOutputGroup == 0) Then
        If (Reqs%nFieldGroups >= Reqs%MaxFieldOutputGroups) Then
          Call Message(                                                                                                       &
                 'FATAL ERROR: too many output groups requested. '                  // &
                 'The maximum number of field output groups is set to: '            // &
                  Trim(Int2Char(Reqs%MaxFieldOutputGroups))                         // &
                 ' and can be altered via the variable "Max # Field Output Groups"' // &
                 ' in the "Main Options" input block',                                 &
                 3                                                                                                           &
               )
        End If
        Reqs%nFieldGroups     = Reqs%nFieldGroups + 1
        FieldReq%iOutputGroup = Reqs%nFieldGroups
      End If
    End If

  End Do

  ! Sets of particle/puff information.
  Do i = 1, Reqs%nPPInfoReqs

    PPInfoReq => Reqs%PPInfoReqs(i)

    ! iSource.
    If (PPInfoReq%SourceName == ' ') Then
      PPInfoReq%iSource = 0
    Else
      PPInfoReq%iSource = FindSourceIndex(PPInfoReq%SourceName, Sources)
    End If

    ! iHCoord.
    PPInfoReq%iHCoord = FindHCoordIndex(PPInfoReq%HCoordName, Coords)

    ! iZCoord.
    PPInfoReq%iZCoord = FindZCoordIndex(PPInfoReq%ZCoordName, Coords)

    ! iTGrid.
    If (PPInfoReq%TGridName == ' ') Then
      PPInfoReq%iTGrid = 0
    Else
      PPInfoReq%iTGrid = FindTGridIndex(PPInfoReq%TGridName, Grids)
    End If

  End Do

  ! Pdfs.
  Do i = 1, Reqs%nPdfs
    Reqs%PdfReqs(i)%iHGrid = FindHGridIndex(Reqs%PdfReqs(i)%HGridName, Grids)
    Reqs%PdfReqs(i)%iHCoord = FindHCoordIndex(Grids%HGrids(Reqs%PdfReqs(i)%iHGrid)%HCoordName, Coords)
    If (Reqs%PdfReqs(i)%ZGridName == ' ') Then
      Reqs%PdfReqs(i)%iZGrid  = 0
      Reqs%PdfReqs(i)%iZCoord = 0
    Else
      Reqs%PdfReqs(i)%iZGrid = FindZGridIndex(Reqs%PdfReqs(i)%ZGridName, Grids)
      Reqs%PdfReqs(i)%iZCoord = FindZCoordIndex(Grids%ZGrids(Reqs%PdfReqs(i)%iZGrid)%ZCoordName, Coords)
    End If
    If (Reqs%PdfReqs(i)%TGridName == ' ') Then
      Reqs%PdfReqs(i)%iTGrid = 0
    Else
      Reqs%PdfReqs(i)%iTGrid = FindTGridIndex(Reqs%PdfReqs(i)%TGridName, Grids)
    End If
    If (Reqs%PdfReqs(i)%SpeciesName == ' ') Then
      Reqs%PdfReqs(i)%iSpecies = 0
    Else
      Reqs%PdfReqs(i)%iSpecies = FindSpeciesIndex(Reqs%PdfReqs(i)%SpeciesName, Specieses)
    End If
  End Do

End Subroutine SetUpiReqs

!-------------------------------------------------------------------------------------------------------------

Subroutine FirstLastReqTimes(                           &
             FirstSourceTime, FirstReleaseTime,         &
             CompDomain,                                &
             Grids, Reqs,                               &
             FirstRunTimeForReqs,  LastRunTimeForReqs,  &
             FirstFlowTimeForReqs, LastFlowTimeForReqs, &
             LastDispTimeForReqs                        &
           )
! Calculates the first/last times which the run, flow modules, and dispersion calculation must cover as a
! result of the requirements.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In)         :: FirstSourceTime
  Type(ShortTime_), Intent(In)         :: FirstReleaseTime
  Type(Domain_),    Intent(In)         :: CompDomain
  Type(Grids_),     Intent(In), Target :: Grids
  Type(Reqs_),      Intent(In), Target :: Reqs
  Type(ShortTime_), Intent(Out)        :: FirstRunTimeForReqs
  Type(ShortTime_), Intent(Out)        :: LastRunTimeForReqs
  Type(ShortTime_), Intent(Out)        :: FirstFlowTimeForReqs
  Type(ShortTime_), Intent(Out)        :: LastFlowTimeForReqs
  Type(ShortTime_), Intent(Out)        :: LastDispTimeForReqs
  ! FirstSourceTime      :: First source emission time for all sources.
  ! FirstReleaseTime     :: First release time for all sources.
  ! CompDomain           :: Computational domain.
  ! Grids                :: Collection of grids.
  ! Reqs                 :: Collection of requirements.
  ! FirstRunTimeForReqs  :} First and last times which run must cover as a result of the requirements.
  ! LastRunTimeForReqs   :}
  ! FirstFlowTimeForReqs :] First and last times which flow modules must cover as a result of the
  ! LastFlowTimeForReqs  :] requirements.
  ! LastDispTimeForReqs  :: Last time which dispersion calculation must cover as a result of the requirements.
  !
  ! Note that in FirstSourceTime and FirstReleaseTime we distinguish between source emission times and release
  ! times. A release time is the actual time of release of a particle/puff in the model and may correspond to
  ! a range of source emission times if the particle/puff has a time spread.
  !
  ! The precise interpretation of FirstRunTimeForReqs, LastRunTimeForReqs, FirstFlowTimeForReqs,
  ! LastFlowTimeForReqs and LastDispTimeForReqs is quite subtle.
  !
  ! [FirstRunTimeForReqs, LastRunTimeForReqs] is a lower limit on the period which the run must cover
  ! resulting from the requirements. The run must also cover the dispersion calculation period, from the first
  ! release till the demise of all tracers.
  !
  ! [FirstFlowTimeForReqs, LastFlowTimeForReqs] is a lower limit on the period within the run which the flow
  ! modules must cover resulting from the requirements (but note there is no need to extend the run to cover
  ! this period if the run would otherwise be shorter). The flow must also cover the dispersion calculation
  ! period, from the first release till the demise of all tracers.
  !
  ! (-infinity, LastDispTimeForReqs] is an upper limit on the dispersion calculation period resulting from the
  ! requirements, although it does not account for the time (if any) between the trailing time-edge and the
  ! nominal time of the particles/puffs. The dispersion calculation period is also restricted by the period
  ! from the first release of tracer to the demise of all tracers.
  !
  ! Locals:
  Type(TGrid_),     Pointer :: TGrid     !} Abbreviations for grids, processes and requirements.
  Type(Process_),   Pointer :: Process   !}
  Type(FieldReq_),  Pointer :: FieldReq  !}
  Type(PPInfoReq_), Pointer :: PPInfoReq !}
  Integer                   :: iField    ! Index of field.
  Integer                   :: iPPInfo   ! Index of set of particle/puff information.
  Type(ShortTime_)          :: T         ! Various times which limit FirstRunTimeForReqs, LastRunTimeForReqs,
                                         ! FirstFlowTimeForReqs, LastFlowTimeForReqs and/or
                                         ! LastDispTimeForReqs.
  Integer                   :: i         ! Loop index.

  !--------------------------------------------!
  ! FirstRunTimeForReqs, FirstFlowTimeForReqs. !
  !--------------------------------------------!

  ! Note: in the following 'output field reqs' means field reqs that are output and not just used in
  ! calculating other reqs.
  !
  ! FirstRunTimeForReqs is calculated only from output field reqs where the underlying quantity depends on
  ! time, with, for type F, O, Q and R reqs only, account taken of the time of earlier contributions to these
  ! reqs due to processing steps (for other types, contributions are zero if concentrations are zero).
  !
  ! FirstFlowTimeForReqs is calculated only from output field reqs of type F and R where the underlying
  ! quantity depends on time, with account taken of the time of earlier contributions to these reqs due to
  ! processing steps (for other types, contributions do not depend on the flow if concentrations are zero).
  !
  ! If FirstSourceTime = FirstReleaseTime we know the run/flow must start at least as early as FirstSourceTime
  ! and so there is no need to account for any more reqs. However if FirstSourceTime < FirstReleaseTime we
  ! only know the run must start at least as early as FirstReleaseTime and so we need to account for some more
  ! reqs. The extra reqs that we need to account for are as follows:
  !
  ! FirstRunTimeForReqs needs to account for field reqs at or after FirstSourceTime where the underlying
  ! quantity depends on time.
  !
  ! FirstFlowTimeForReqs needs to account for field reqs of type F at or after FirstSourceTime where the
  ! underlying quantity depends on time.

  ! For FirstSourceTime < FirstReleaseTime, the extra constraints will work completely generally although they
  ! may not be optimal for the general case and may be more complex than needed for the way the module is
  ! actually used (the use of FirstSourceTime < FirstReleaseTime is restricted to cases where the
  ! particle/puff reuse option is used and FirstSourceTime is in the infinite past - this implies also that
  ! the met is fixed and all flow module domains are unbounded in time; however these restrictions are not
  ! assumed in this module).
  !
  ! Non-optimality:
  ! Firstly the limiting requirement may be one which isn't used. This can arise if, in order to keep the time
  ! grid regular, the grid contains more times than required. However this is not an issue for the way the
  ! module is actually used because FirstSourceTime will be in the infinite past. Secondly the limiting req
  ! may be a type 'L' req (such reqs are not calculated till after FirstReleaseTime). However here it is
  ! convenient to start earlier to provide a uniform way of handling Field%iTL, Field%LastProcThisCase etc. It
  ! would also seem strange to the user if, e.g., one had a continuous source over all time with output at a
  ! single time T and the run started later than T.
  !
  ! Unnecessary complexity:
  ! For the way the module is actually used the met will be fixed and all flow module domains will be
  ! unbounded in time. Hence we could account for field reqs of all types in FirstFlowTimeForReqs without
  ! altering things.

  FirstRunTimeForReqs  = InfFutureShortTime()
  FirstFlowTimeForReqs = InfFutureShortTime()

  ! Loop over fields.
  Do iField = 1, Reqs%nFieldReqs

    FieldReq => Reqs%FieldReqs(iField)

    If (.not.FieldReq%Disk .and. .not.FieldReq%Screen .and. .not.FieldReq%Graph) Cycle

    If (Scan(QInfo(FieldReq%iQuantity), 'T') == 0) Cycle

    ! For this requirement, calculate first time at or after the start of the computational domain.
    If (FieldReq%iTGrid == 0) Then
      T = DomainEnd(CompDomain)
    Else
      TGrid => Grids%TGrids(FieldReq%iTGrid)
      Call FirstTAfterT(TGrid, DomainStart(CompDomain), .false., T)
      If (T > DomainEnd(CompDomain) .or. IsInfFuture(T)) Cycle
    End If

    ! Adjust T for processing for type F, O, Q and R reqs.
    If (                               &
      FieldReq%LEFOPQRType == 'F' .or. &
      FieldReq%LEFOPQRType == 'O' .or. &
      FieldReq%LEFOPQRType == 'Q' .or. &
      FieldReq%LEFOPQRType == 'R'      &
    ) Then

      Do i = FieldReq%nProcesses, FieldReq%iFOQRProc + 1, -1

        ! If T is -infinity, no need to adjust further.
        If (IsInfPast(T)) Exit

        Process => Reqs%Processes(FieldReq%iProcesses(i))

        ! 'T window defined by grid' case.
        If (Process%UseTGrid) Then

          ! $$ needs more support from GridAndDomain to do this

        ! 'Rolling T window' case.
        Else If (Process%nT > 1) Then

          T = T - Process%sdT * (Process%nT - 1)
          If (T < DomainStart(CompDomain)) Then ! $$ move error message to setupreqs?
            Call Message(                                                                             &
                   'FATAL ERROR: The processing specified for field requirement "'                 // &
                   Trim(FieldReq%Name)                                                             // &
                   '" means that contributions from before the start of the computational domain ' // &
                   'are needed for results at times lying in the computational domain.',              &
                   3                                                                                  &
                 )
          End If

        ! Other processing cases have no effect here or are prohibited.
        Else

        End If

      End Do

    End If

    ! FirstRunTimeForReqs.
    FirstRunTimeForReqs = TMin(FirstRunTimeForReqs, T)

    ! FirstFlowTimeForReqs.
    If (FieldReq%LEFOPQRType == 'F' .or. FieldReq%LEFOPQRType == 'R') Then
      FirstFlowTimeForReqs = TMin(FirstFlowTimeForReqs, T)
    End If

  End Do

  ! Error for FirstRunTimeForReqs = -infinity.
  If (IsInfPast(FirstRunTimeForReqs)) Then
    Call Message(                                                                            &
           'FATAL ERROR: Either the start time of the computational domain or the first ' // &
           'output time of every requirement must not be in the infinite past.',             &
           3                                                                                 &
         )
  End If

  If (FirstSourceTime < FirstReleaseTime) Then

    ! Loop over fields.
    Do iField = 1, Reqs%nFieldReqs

      FieldReq => Reqs%FieldReqs(iField)

      If (Scan(QInfo(FieldReq%iQuantity), 'T') == 0) Cycle

      ! For this requirement, calculate first time at or after the start of the computational domain.
      If (FieldReq%iTGrid == 0) Then
        T = DomainEnd(CompDomain)
      Else
        TGrid => Grids%TGrids(FieldReq%iTGrid)
        Call FirstTAfterT(TGrid, DomainStart(CompDomain), .false., T)
        If (T > DomainEnd(CompDomain) .or. IsInfFuture(T)) Cycle
      End If

      If (T < FirstSourceTime) Cycle

      ! FirstRunTimeForReqs.
      FirstRunTimeForReqs = TMin(FirstRunTimeForReqs, T)

      ! FirstFlowTimeForReqs.
      If (FieldReq%LEFOPQRType == 'F') Then
        FirstFlowTimeForReqs = TMin(FirstFlowTimeForReqs, T)
      End If

    End Do

    ! Error for FirstRunTimeForReqs = -infinity.
    If (IsInfPast(FirstRunTimeForReqs)) Then
      Call Message(                                                                                       &
             'FATAL ERROR: The start time of the computational domain and at least one source are in ' // &
             'the infinite past and there is an output requirement that requires contributions '       // &
             'in the infinite past (other than through intrinsic averaging/integrating).',                &
             3                                                                                            &
           )
    End If

  End If

  !---------------------------------------------------------------!
  ! LastRunTimeForReqs, LastFlowTimeForReqs, LastDispTimeForReqs. !
  !---------------------------------------------------------------!

  ! LastRunTimeForReqs is calculated from t-grided output field reqs.
  !
  ! LastFlowTimeForReqs is calculated from t-grided output field reqs of type F and R only.
  !
  ! LastDispTimeForReqs is calculated from field reqs of type L and E (here the underlying quantity always
  ! depends on time) and from reqs for sets of particle/puff information.

  LastRunTimeForReqs  = InfPastShortTime()
  LastFlowTimeForReqs = InfPastShortTime()
  LastDispTimeForReqs = InfPastShortTime()

  ! Loop over fields.
  Do iField = 1, Reqs%nFieldReqs

    FieldReq => Reqs%FieldReqs(iField)

    If (Scan(QInfo(FieldReq%iQuantity), 'T') == 0) Cycle

    ! For this requirement, calculate last time at or before the end of the computational domain.
    If (FieldReq%iTGrid == 0) Then
      T = DomainEnd(CompDomain)
    Else
      TGrid => Grids%TGrids(FieldReq%iTGrid)
      Call LastTBeforeT(TGrid, DomainEnd(CompDomain), .false., T)
      If (T < DomainStart(CompDomain) .or. IsInfPast(T)) Cycle
    End If

    ! LastRunTimeForReqs.
    If (                                                             &
      FieldReq%iTGrid /= 0                                     .and. &
      (FieldReq%Disk .or. FieldReq%Screen .or. FieldReq%Graph)       &
    ) Then
      LastRunTimeForReqs = TMax(LastRunTimeForReqs, T)
    End If

    ! LastFlowTimeForReqs.
    If (                                                                   &
      FieldReq%iTGrid /= 0                                           .and. &
      (FieldReq%Disk .or. FieldReq%Screen .or. FieldReq%Graph)       .and. &
      (FieldReq%LEFOPQRType == 'F' .or. FieldReq%LEFOPQRType == 'R')       &
    ) Then
      LastFlowTimeForReqs = TMax(LastFlowTimeForReqs, T)
    End If

    ! LastDispTimeForReqs.
    If (FieldReq%LEFOPQRType == 'L' .or. FieldReq%LEFOPQRType == 'E') Then
      LastDispTimeForReqs = TMax(LastDispTimeForReqs, T)
    End If

  End Do

  ! Loop over sets of particle/puff information.
  Do iPPInfo = 1, Reqs%nPPInfoReqs

    PPInfoReq => Reqs%PPInfoReqs(iPPInfo)

    ! For this requirement, calculate last time at or before the end of the computational domain.
    If (PPInfoReq%iTGrid == 0) Then
      T = DomainEnd(CompDomain)
    Else
      TGrid => Grids%TGrids(PPInfoReq%iTGrid)
      Call LastTBeforeT(TGrid, DomainEnd(CompDomain), .false., T)
      If (T < DomainStart(CompDomain) .or. IsInfPast(T)) Cycle
    End If

    ! LastDispTimeForReqs.
    LastDispTimeForReqs = TMax(LastDispTimeForReqs, T)

  End Do

End Subroutine FirstLastReqTimes

!-------------------------------------------------------------------------------------------------------------

Function InitResults(RunName, StartClockTime, Reqs) Result(Results)
! Initialises a collection of results.

  Implicit None
  ! Argument list:
  Character(MaxCharLength), Intent(In)    :: RunName        ! Name of run.
  Type(Time_),              Intent(In)    :: StartClockTime ! Clock time of the start of the run.
  Type(Reqs_),              Intent(In)    :: Reqs           ! Collection of requirements.
  
  ! Function result:
  Type(Results_) :: Results ! Collection of  results.
  ! Locals:
  Integer :: i ! Loop index.

  ! Name and clock time of the start time of the run.
  Results%RunName        = RunName
  Results%StartClockTime = StartClockTime
  
  ! Allocation of memory for fields
  Allocate(Results%Fields(Reqs%MaxFieldReqs))
  Allocate(Results%FieldnLines(Reqs%MaxFieldOutputGroups, MaxFieldOutputFiles))
  Allocate(Results%FieldDiskUnits(Reqs%MaxFieldOutputGroups, MaxFieldOutputFiles))
  Allocate(Results%FieldScreenUnits(Reqs%MaxFieldOutputGroups, MaxFieldOutputFiles))

  ! Set pointers to null.
  Do i = 1, Reqs%MaxFieldReqs
    Results%Fields(i)%Std    => Null()
    Results%Fields(i)%P64    => Null()
    Results%Fields(i)%T      => Null()
    Results%Fields(i)%MaxStd => Null()
    Results%Fields(i)%MinStd => Null()
    Results%Fields(i)%S      => Null()
  End Do
  Do i = 1, MaxPdfReqs
    Results%Pdfs(i)%Scale => Null()
    Results%Pdfs(i)%Data  => Null()
  End Do

  ! Set MaxTimeField to the infinite past.
  Results%MaxTimeField = InfPastShortTime()

End Function InitResults

!-------------------------------------------------------------------------------------------------------------

Subroutine WriteRestartFileResults(Unit, Results)
! Writes the collection of results to the restart file.

  Implicit None
  ! Argument list:
  Integer,        Intent(In)         :: Unit    ! Input/output unit number.
  Type(Results_), Intent(In), Target :: Results ! Collection of results.
  ! Locals:
  Type(Field_), Pointer :: Field ! Abbreviations for results.
  Integer               :: i     ! Loop index.

  ! Note the I/O units stored in Results (i.e. FieldDiskUnits, FieldScreenUnits, Field%GraphUnit,
  ! PdfDiskUnits, PdfScreenUnits, PPInfo%DiskUnits, PPInfo%ScreenUnits and PPInfo%GraphUnit) are not needed in
  ! the restart file and so are only written if it is simpler to do so than not.

  ! Run information.
  Write (Unit) Results%RunName, Results%StartClockTime

  ! Fields.
  Write (Unit) Results%nFields

  Do i = 1, Results%nFields

    Field => Results%Fields(i)

    Write (Unit) Field%TypeType, Field%nT, Field%iTL

    Write (Unit) Associated(Field%Std)
    If (Associated(Field%Std)) Then
      Write (Unit) Shape(Field%Std)
      Write (Unit) Field%Std(:, :, :, :, :, :, :)
    End If

    Write (Unit) Associated(Field%P64)
    If (Associated(Field%P64)) Then
      Write (Unit) Shape(Field%P64)
      Write (Unit) Field%P64(:, :, :, :, :, :, :)
    End If

    Write (Unit) Associated(Field%T)
    If (Associated(Field%T)) Then
      Write (Unit) Shape(Field%T)
      Write (Unit) Field%T(:, :, :, :, :, :, :)
    End If

    Write (Unit) Associated(Field%MaxStd)
    If (Associated(Field%MaxStd)) Then
      Write (Unit) Shape(Field%MaxStd)
      Write (Unit) Field%MaxStd(:, :, :, :, :, :)
    End If

    Write (Unit) Associated(Field%MinStd)
    If (Associated(Field%MinStd)) Then
      Write (Unit) Shape(Field%MinStd)
      Write (Unit) Field%MinStd(:, :, :, :, :, :)
    End If

    Write (Unit) Associated(Field%S)
    If (Associated(Field%S)) Then
      Write (Unit) Shape(Field%S)
      Write (Unit) Field%S(:, :, :, :, :, :)
    End If

    Write (Unit) Field%LastProcThisCase, Field%LastProc, Field%LastNumOutput, Field%LastGraphOutput

  End Do

  Write (Unit) Results%FieldnLines(:, :)
  Write (Unit) Results%MaxTimeField

  ! Pdfs.
  Write (Unit) Results%nPdfs
  Do i = 1, Results%nPdfs
    Write (Unit) Shape(Results%Pdfs(i)%Scale)
    Write (Unit) Shape(Results%Pdfs(i)%Data)
    Write (Unit) Results%Pdfs(i)%PdfType
    Write (Unit) Results%Pdfs(i)%PdfMode
    Write (Unit) Results%Pdfs(i)%FixedThresholds
    Write (Unit) Results%Pdfs(i)%PdfSize
    Write (Unit) Results%Pdfs(i)%Thresholds(:)
    Write (Unit) Results%Pdfs(i)%Scale(:, :, :, :)
    Write (Unit) Results%Pdfs(i)%Data(:, :, :, :, :)
  End Do
  Write (Unit) Results%PdfsProc(:, :)
  Write (Unit) Results%PdfsNextOutput(:)
  Write (Unit) Results%PdfFiles(:)

  ! Sets of particle/puff information.
  Write (Unit) Results%nPPInfos
  Write (Unit) Results%PPInfos(1:Results%nPPInfos)

End Subroutine WriteRestartFileResults

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadRestartFileResults(Unit, Results)
! Reads the collection of results from the restart file.

  Implicit None
  ! Argument list:
  Integer,        Intent(In)          :: Unit    ! Input/output unit number.
  Type(Results_), Intent(Out), Target :: Results ! Collection of results.
  ! Locals:
  Type(Field_), Pointer :: Field         ! Abbreviation for field result.
  Logical               :: IsAssociated  ! Indicates the allocatable array should be allocated.
  Integer               :: ArrayShape(7) ! Shape of array to be allocated.
  Integer               :: i             ! Loop index.
  Integer               :: ErrorCode     ! Error code.

  ! Note the I/O units stored in Results (i.e. FieldDiskUnits, FieldScreenUnits, Field%GraphUnit,
  ! PdfDiskUnits, PdfScreenUnits, PPInfo%DiskUnits, PPInfo%ScreenUnits and PPInfo%GraphUnit) are set to zero
  ! to indicate the files need opening.

  ! Run information.
  Read (Unit, End = 1, Err = 2) Results%RunName, Results%StartClockTime

  ! Fields.
  Read (Unit, End = 1, Err = 2) Results%nFields

  Do i = 1, Results%nFields

    Field => Results%Fields(i)

    Read (Unit, End = 1, Err = 2) Field%TypeType, Field%nT, Field%iTL

    Read (Unit, End = 1, Err = 2) IsAssociated
    If (IsAssociated) Then
      Read (Unit, End = 1, Err = 2) ArrayShape
      Allocate(          &
        Field%Std(       &
          ArrayShape(1), &
          ArrayShape(2), &
          ArrayShape(3), &
          ArrayShape(4), &
          ArrayShape(5), &
          ArrayShape(6), &
          ArrayShape(7)  &
        ),               &
        Stat = ErrorCode &
      )
      If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for results', 3)
      Read (Unit, End = 1, Err = 2) Field%Std(:, :, :, :, :, :, :)
    End If

    Read (Unit, End = 1, Err = 2) IsAssociated
    If (IsAssociated) Then
      Read (Unit, End = 1, Err = 2) ArrayShape
      Allocate(          &
        Field%P64(       &
          ArrayShape(1), &
          ArrayShape(2), &
          ArrayShape(3), &
          ArrayShape(4), &
          ArrayShape(5), &
          ArrayShape(6), &
          ArrayShape(7)  &
        ),               &
        Stat = ErrorCode &
      )
      If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for results', 3)
      Read (Unit, End = 1, Err = 2) Field%P64(:, :, :, :, :, :, :)
    End If

    Read (Unit, End = 1, Err = 2) IsAssociated
    If (IsAssociated) Then
      Read (Unit, End = 1, Err = 2) ArrayShape
      Allocate(          &
        Field%T(         &
          ArrayShape(1), &
          ArrayShape(2), &
          ArrayShape(3), &
          ArrayShape(4), &
          ArrayShape(5), &
          ArrayShape(6), &
          ArrayShape(7)  &
        ),               &
        Stat = ErrorCode &
      )
      If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for results', 3)
      Read (Unit, End = 1, Err = 2) Field%T(:, :, :, :, :, :, :)
    End If

    Read (Unit, End = 1, Err = 2) IsAssociated
    If (IsAssociated) Then
      Read (Unit, End = 1, Err = 2) ArrayShape(1:6)
      Allocate(          &
        Field%MaxStd(    &
          ArrayShape(1), &
          ArrayShape(2), &
          ArrayShape(3), &
          ArrayShape(4), &
          ArrayShape(5), &
          ArrayShape(6)  &
        ),               &
        Stat = ErrorCode &
      )
      If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for results', 3)
      Read (Unit, End = 1, Err = 2) Field%MaxStd(:, :, :, :, :, :)
    End If

    Read (Unit, End = 1, Err = 2) IsAssociated
    If (IsAssociated) Then
      Read (Unit, End = 1, Err = 2) ArrayShape(1:6)
      Allocate(          &
        Field%MinStd(    &
          ArrayShape(1), &
          ArrayShape(2), &
          ArrayShape(3), &
          ArrayShape(4), &
          ArrayShape(5), &
          ArrayShape(6)  &
        ),               &
        Stat = ErrorCode &
      )
      If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for results', 3)
      Read (Unit, End = 1, Err = 2) Field%MinStd(:, :, :, :, :, :)
    End If

    Read (Unit, End = 1, Err = 2) IsAssociated
    If (IsAssociated) Then
      Read (Unit, End = 1, Err = 2) ArrayShape(1:6)
      Allocate(          &
        Field%S(         &
          ArrayShape(1), &
          ArrayShape(2), &
          ArrayShape(3), &
          ArrayShape(4), &
          ArrayShape(5), &
          ArrayShape(6)  &
        ),               &
        Stat = ErrorCode &
      )
      If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for results', 3)
      Read (Unit, End = 1, Err = 2) Field%S(:, :, :, :, :, :)
    End If

    Read (Unit, End = 1, Err = 2) Field%LastProcThisCase, Field%LastProc, Field%LastNumOutput, &
                                  Field%LastGraphOutput
    Field%GraphUnit = 0

  End Do

  Read (Unit, End = 1, Err = 2) Results%FieldnLines(:, :)
  Results%FieldDiskUnits  (:, :) = 0
  Results%FieldScreenUnits(:, :) = 0
  Read (Unit, End = 1, Err = 2) Results%MaxTimeField

  ! Pdfs.
  Read (Unit, End = 1, Err = 2) Results%nPdfs
  Do i = 1, Results%nPdfs
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(                &
      Results%Pdfs(i)%Scale( &
        ArrayShape(1),       &
        ArrayShape(2),       &
        ArrayShape(3),       &
        ArrayShape(4)        &
      ),                     &
      Stat = ErrorCode       &
    )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for results', 3)
    Read (Unit, End = 1, Err = 2) ArrayShape(1:5)
    Allocate(               &
      Results%Pdfs(i)%Data( &
        ArrayShape(1),      &
        ArrayShape(2),      &
        ArrayShape(3),      &
        ArrayShape(4),      &
        ArrayShape(5)       &
      ),                    &
      Stat = ErrorCode      &
    )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for results', 3)
    Read (Unit, End = 1, Err = 2) Results%Pdfs(i)%PdfType
    Read (Unit, End = 1, Err = 2) Results%Pdfs(i)%PdfMode
    Read (Unit, End = 1, Err = 2) Results%Pdfs(i)%FixedThresholds
    Read (Unit, End = 1, Err = 2) Results%Pdfs(i)%PdfSize
    Read (Unit, End = 1, Err = 2) Results%Pdfs(i)%Thresholds(:)
    Read (Unit, End = 1, Err = 2) Results%Pdfs(i)%Scale(:, :, :, :)
    Read (Unit, End = 1, Err = 2) Results%Pdfs(i)%Data(:, :, :, :, :)
  End Do
  Read (Unit, End = 1, Err = 2) Results%PdfsProc(:, :)
  Read (Unit, End = 1, Err = 2) Results%PdfsNextOutput(:)
  Read (Unit, End = 1, Err = 2) Results%PdfFiles(:)
  Results%PdfDiskUnits  (:) = 0
  Results%PdfScreenUnits(:) = 0

  ! Sets of particle/puff information.
  Read (Unit, End = 1, Err = 2) Results%nPPInfos
  Read (Unit, End = 1, Err = 2) Results%PPInfos(1:Results%nPPInfos)

  Do i = 1, Results%nPPInfos
    Results%PPInfos(i)%DiskUnits  (:) = 0
    Results%PPInfos(i)%ScreenUnits(:) = 0
    Results%PPInfos(i)%GraphUnit      = 0
  End Do

  Return

1 Continue
  Call Message('FATAL ERROR: The restart file is shorter than expected', 3)

2 Continue
  Call Message('FATAL ERROR: An error occurred when reading the restart file', 3)

End Subroutine ReadRestartFileResults

!-------------------------------------------------------------------------------------------------------------

Subroutine RestartAdjustmentForResults(Reqs, Results)
! Adjusts results read from the restart file.

  Implicit None
  ! Argument list:
  Type(Reqs_),    Intent(In)    :: Reqs    ! Collection of requirements.
  Type(Results_), Intent(InOut) :: Results ! Collection of results.
  ! Locals:
  Integer :: i ! Loop index.

  ! If TNewFile is set, ensure output files start afresh (possibly overwriting any previous files).
  Do i = 1, Reqs%nFieldReqs
    If (Reqs%FieldReqs(i)%TNewFile) Then
      Results%FieldnLines(Reqs%FieldReqs(i)%iOutputGroup, :) = 0
    End If
  End Do

End Subroutine RestartAdjustmentForResults

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareResults(                                                           &
             StartTime, MaxEndTime, Syncdt, MaxPPTimeLead, MaxPPTimeLag, CompDomain, &
             Grids, Reqs,                                                            &
             Results                                                                 &
           )
! Prepares the collection of results for a new case.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In)            :: StartTime
  Type(ShortTime_), Intent(In)            :: MaxEndTime
  Type(ShortTime_), Intent(In)            :: SyncdT
  Type(ShortTime_), Intent(In)            :: MaxPPTimeLead
  Type(ShortTime_), Intent(In)            :: MaxPPTimeLag
  Type(Domain_),    Intent(In)            :: CompDomain
  Type(Grids_),     Intent(In),    Target :: Grids
  Type(Reqs_),      Intent(In),    Target :: Reqs
  Type(Results_),   Intent(InOut), Target :: Results
  ! StartTime     :: Start time of run. This must be the earlier of FirstRunTimeForReqs as computed by
  !                  FirstLastReqTimes and FirstReleaseTime as computed by SourceModule.FirstLastReleaseTimes.
  ! MaxEndTime    :: Latest possible end time of run. This must be the later of LastRunTimeForReqs and
  !                  LastDispTimeForReqs as computed by FirstLastReqTimes. However, like LastDispTimeForReqs,
  !                  it does not account for the time (if any) between the trailing time-edge and the nominal
  !                  time of the particles/puffs.
  ! SyncdT        :: The time interval at which particles/puffs are synchronised.
  ! MaxPPTimeLead :: The maximum possible interval which can occur at any time between the leading time-edge
  !                  and the nominal time of the particles/puffs.
  ! MaxPPTimeLag  :: The maximum possible interval which can occur at any time between the trailing time-edge
  !                  and the nominal time of the particles/puffs.
  ! CompDomain    :: Computational domain.
  ! Grids         :: Collection of grids.
  ! Reqs          :: Collection of requirements.
  ! Results       :: Collection of results.
  ! Locals:
  Type(TGrid_),     Pointer :: TGrid         !} Abbreviations for grids, requirements and results.
  Type(DGrid_),     Pointer :: DGrid         !}
  Type(Process_),   Pointer :: Process       !}
  Type(FieldReq_),  Pointer :: FieldReq      !}
  Type(PPInfoReq_), Pointer :: PPInfoReq     !}
  Type(Field_),     Pointer :: Field         !}
  Type(PPInfo_),    Pointer :: PPInfo        !}
  Integer                   :: iTGrid        !] Indices of grids.
  Integer                   :: iSGrid        !]
  Integer                   :: iHGrid        !]
  Integer                   :: iZGrid        !]
  Integer                   :: iDGrid        !]
  Integer                   :: nT            !} Size of grids.
  Integer                   :: nS            !}
  Integer                   :: nX            !}
  Integer                   :: nY            !}
  Integer                   :: nZ            !}
  Integer                   :: nD            !}
  Integer                   :: nLast         ! Size of last (7-th) dimension of field arrays
                                             ! (storing the components and/or size ranges).
  Integer                   :: iTL           ! Index of next time in time grid which is >= StartTime.
  Type(ShortTime_)          :: ReqMemoryTime ! Duration for which results may need to be remembered, excluding
                                             ! any memory needed for time processing purposes.
  Type(ShortTime_)          :: ProcTime      ! Duration for which fields may need to be remembered for time
                                             ! processing purposes.
  Integer                   :: iT1           !} Smallest and largest grid indices which may be used.
  Integer                   :: iT2           !}
  Integer                   :: i             ! Loop index.
  Integer                   :: Error         ! Error code.
  Type(ShortTime_)          :: DummyT        ! Dummy variable needed for an argument list.

  ! ReqMemoryTime for fields.
  ReqMemoryTime = SyncdT + MaxPPTimeLead + MaxPPTimeLag

  ! Fields.
  Do i = 1, Reqs%nFieldReqs

    ! Abbreviations for requirement and result.
    FieldReq => Reqs%FieldReqs(i)
    Field    => Results%Fields(i)

    iTGrid = FieldReq%iTGrid
    If (iTGrid /= 0) Then
      TGrid => Grids%TGrids(iTGrid)
    End If

    If (.not.Associated(Field%Std) .and. .not.Associated(Field%P64) .and. .not.Associated(Field%T)) Then

      iSGrid = FieldReq%iSGrid
      iHGrid = FieldReq%iHGrid
      iZGrid = FieldReq%iZGrid
      iDGrid = FieldReq%iDGrid

      ! ProcTime.

      ProcTime = ZeroShortTime()

      If (FieldReq%nProcesses > FieldReq%iPQRProc) Then

        Process => Reqs%Processes(FieldReq%iProcesses(FieldReq%nProcesses))
        ! $$ should use iLastProcess? - needs careful consideration.
        ! Don't want to exclude intrnisic av cases which will contribute to proc time.
        ! But if proc step does nothing at all, need to check proc time of previous step is accounted for
        ! (might not be if fieldreq chain bypasses the unnecessary fieldreq - need to check).

        ! 'T window defined by grid' case.
        If (Process%UseTGrid) Then

          ! $$ needs more support from GridAndDomain to do this

        ! 'Rolling T window' case.
        Else If (Process%nT > 1) Then

          If (FieldReq%nProcesses == 1 .and. FieldReq%IntrinsicTAv) Then
            ProcTime = Process%sT
          Else
            ProcTime = Process%sdT * (Process%nT - 1)
          End If

        ! 'T window from T = -infinity' case with infinity > dT > 0.
        Else If (IsInfFuture(Process%sT) .and. .not. IsInfFuture(Process%sdT)) Then

          ProcTime = Process%sdT ! $$ review when writing rest of code for this case.

        ! 'T window from T = -infinity' case with dT = infinity.
        Else If (IsInfFuture(Process%sT) .and. IsInfFuture(Process%sdT)) Then

          ! $$ review when writing rest of code for this case.
          If (FieldReq%nProcesses == 1 .and. FieldReq%IntrinsicTAv) Then
            ProcTime = Process%sT
          End If

        End If

      End If

      ! nT.
      If (iTGrid == 0) Then
        nT = 1
      Else
        Call FirstTAfterT(TGrid, StartTime,                             .false., DummyT, iT1)
        Call LastTBeforeT(TGrid, MaxEndTime + ReqMemoryTime + ProcTime, .false., DummyT, iT2)
        If (iT2 < iT1) Then
          nT = 0
        Else If (iT1 == -Huge(iT1)/3 .or. iT2 == Huge(iT2)/3) Then
          nT = Huge(nT)/3
        Else
          nT = iT2 - iT1 + 1
        End If
        If (IsInfFuture(ReqMemoryTime)) Then
        Else If (FieldReq%Ensemble) Then
        Else If (.not.FieldReq%TSeparateFile .and. FieldReq%TAcross) Then
        Else
          nT = Min((ReqMemoryTime + ProcTime) / TGrid%sdT + 1, nT)
        End If
        If (nT == Huge(nT)/3) Then
          ! $$ trap too large nT once have infinite grids.
        End If
      End If

      ! nS.
      If (iSGrid == 0) Then
        nS = 1
      Else
        nS = Grids%TGrids(iSGrid)%nT
      End If

      ! nX, nY.
      If (iHGrid == 0) Then
        nX = 1
        nY = 1
      Else
        nX = Grids%HGrids(iHGrid)%nX
        nY = Grids%HGrids(iHGrid)%nY
      End If

      ! nZ.
      If (iZGrid == 0) Then
        nZ = 1
      Else
        nZ = Grids%ZGrids(iZGrid)%nZ
      End If

      ! nD.
      If (iDGrid == 0) Then
        nD = 1
      Else
        nD = nDInDGrid(Grids%DGrids(iDGrid))
      End If

      ! nLast
      nLast = FieldReq%nComponents * FieldReq%nRanges

      If (FieldReq%iExt == 0 .and. nT /= 0) Then

        If (FieldReq%TypeType == ' ') Then
          If (FieldReq%iQuantity == Q_PhotonFlux) Then ! $$ better not to hard code different quantities here.
            Allocate(Field%Std(nD, nX, nY, nZ, nS, nT, MaxEnergies), Stat = Error)
          Else
            Allocate(Field%Std(nD, nX, nY, nZ, nS, nT, nLast), Stat = Error)
          End If
        Else If (FieldReq%TypeType == '8') Then
          Allocate(Field%P64(nD, nX, nY, nZ, nS, nT, nLast), Stat = Error)
        Else If (FieldReq%TypeType == ':') Then
          Allocate(Field%T  (nD, nX, nY, nZ, nS, nT, nLast), Stat = Error)
        End If
        If (Error /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for results', 3)

        If (iDGrid /= 0) Then
          If (FloatingDGrid(Grids%DGrids(iDGrid))) Then
            Allocate(Field%S(nX, nY, nZ, nS, nT, nLast), Stat = Error)
          End If
        End If
        If (Error /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for results', 3)

        If (FieldReq%iLastProcess > 0) Then ! $$ if use nProcesses would need to associate Process
                                            ! Not completely clear if nProcesses or iLastProcess is correct.
          If (Process%Code == P_Prob) Then
            DGrid => Grids%DGrids(Process%iDGrid)
            If (nDInDGrid(DGrid) > 1 .or. FloatingDGrid(DGrid)) Then
              Allocate(Field%MaxStd(nX, nY, nZ, nS, nT, nLast), Stat = Error)
              If (Error /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for results', 3)
              Allocate(Field%MinStd(nX, nY, nZ, nS, nT, nLast), Stat = Error)
              If (Error /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for results', 3)
            End If
          End If
        End If

        If (Associated(Field%Std   )) Field%Std(:, :, :, :, :, :, :) = 0.0
        If (Associated(Field%P64   )) Field%P64(:, :, :, :, :, :, :) = 0.0
        If (Associated(Field%T     )) Field%T  (:, :, :, :, :, :, :) = ZeroTime()
        If (Associated(Field%MaxStd)) Field%MaxStd(:, :, :, :, :, :) = -Huge(Field%MaxStd)
        If (Associated(Field%MinStd)) Field%MinStd(:, :, :, :, :, :) =  Huge(Field%MinStd)
        If (Associated(Field%S     )) Field%S     (:, :, :, :, :, :) = -Huge(Field%S)

      End If

      Field%TypeType = FieldReq%TypeType

      Field%nT = nT

    End If

    If (.not. FieldReq%Ensemble) Then
      If (Associated(Field%Std   )) Field%Std(:, :, :, :, :, :, :) = 0.0
      If (Associated(Field%P64   )) Field%P64(:, :, :, :, :, :, :) = 0.0
      If (Associated(Field%T     )) Field%T  (:, :, :, :, :, :, :) = ZeroTime()
      If (Associated(Field%MaxStd)) Field%MaxStd(:, :, :, :, :, :) = -Huge(Field%MaxStd)
      If (Associated(Field%MinStd)) Field%MinStd(:, :, :, :, :, :) =  Huge(Field%MinStd)
      If (Associated(Field%S     )) Field%S     (:, :, :, :, :, :) = -Huge(Field%S)
    End If

    If (iTGrid == 0) Then

      Field%iTL              = 1
      Field%LastProcThisCase = 0
      Field%LastProc         = 0
      Field%LastNumOutput    = 0
      Field%LastGraphOutput  = 0

    Else

      ! Calculate index of next time in time grid which is >= StartTime and adjust to ensure the range iTL to
      ! iTL + nT - 1 lies in the time grid. This will be the earliest time stored from now on.
      Call FirstTAfterT(TGrid, StartTime, .false., DummyT, iTL)
      Field%iTL = iTL
      Field%iTL = Min(Field%iTL, TGrid%nT - Field%nT + 1)
      Field%iTL = Max(Field%iTL, 1)

      ! Ensure data which won't be computed isn't processed or output.
      Field%LastProcThisCase = iTL - 1
      Field%LastProc         = iTL - 1
      Field%LastNumOutput    = iTL - 1
      Field%LastGraphOutput  = iTL - 1

    End If

  End Do

  Results%Fields(:)%GraphUnit = 0

  Results%FieldnLines     (:, :) = 0
  Results%FieldDiskUnits  (:, :) = 0
  Results%FieldScreenUnits(:, :) = 0

  Results%nFields = Reqs%nFieldReqs

  ! Pdfs.
  Do i = 1, Reqs%nPdfs

    iTGrid = Reqs%PdfReqs(i)%iTGrid
    iHGrid = Reqs%PdfReqs(i)%iHGrid
    iZGrid = Reqs%PdfReqs(i)%iZGrid

    Results%Pdfs(i)%PdfType = Reqs%PdfReqs(i)%PdfType
    Results%Pdfs(i)%PdfMode = Reqs%PdfReqs(i)%PdfMode
    Results%Pdfs(i)%FixedThresholds = .true.

    nX = Grids%HGrids(iHGrid)%nX
    nY = Grids%HGrids(iHGrid)%nY
    If (iZGrid == 0) Then
      nZ = 1
    Else
      nZ = Grids%ZGrids(iZGrid)%nZ
    End If
    If (iTGrid == 0) Then
      nT = 1
    Else
      nT = Grids%TGrids(iTGrid)%nT
    End If

    If (Reqs%PdfReqs(i)%PdfMode == 1) Then
      If (Reqs%PdfReqs(i)%PdfType == 1) Then
        nD = AutoPdfNumProbabilities
        Results%Pdfs(i)%FixedThresholds = .false.
        Allocate(Results%Pdfs(i)%Scale(nX, nY, nZ, nT), Stat = Error)
        If (Error /= 0) Then
          Call Message('Error in PrepareResults', 3)
        End If
      Else
        nD = AutoPdfNumPercentiles
        Results%Pdfs(i)%Thresholds(1:AutoPdfNumPercentiles) = &
          AutoPdfPercentiles(1:AutoPdfNumPercentiles)
      End If
    Else
      nD = Reqs%PdfReqs(i)%PdfSize
      Results%Pdfs(i)%Thresholds = Reqs%PdfReqs(i)%PdfThresholds
    End If
    Results%Pdfs(i)%PdfSize = nD
    Allocate(Results%Pdfs(i)%Data(nX, nY, nZ, nT, nD), Stat = Error)
    If (Error /= 0) Then
      Call Message('Error in PrepareResults', 3)
    End If
    If (.not.Results%Pdfs(i)%FixedThresholds) Then
      Results%Pdfs(i)%Scale(:, :, :, :) = 0.0
    End If
    Results%Pdfs(i)%Data(:, :, :, :, :) = 0.0
    Results%PdfsProc(:,i)               = .false.

  End Do

  Results%nPdfs = Reqs%nPdfs
  Results%PdfsNextOutput(1:Reqs%nPdfGroups) = 1

  ! ReqMemoryTime for sets of particle/puff information.
  ReqMemoryTime = SyncdT

  ! Sets of particle/puff information.
  Do i = 1, Reqs%nPPInfoReqs

    ! Abbreviations for requirement and result.
    PPInfoReq => Reqs%PPInfoReqs(i)
    PPInfo    => Results%PPInfos(i)

    ! nLines, DiskUnits, ScreenUnits, ParticleFiles, PuffFiles, LastFile and GraphUnit.
    PPInfo%nLines       (:) = 0
    PPInfo%DiskUnits    (:) = 0
    PPInfo%ScreenUnits  (:) = 0
    PPInfo%ParticleFiles(:) = 0
    PPInfo%PuffFiles    (:) = 0
    PPInfo%LastFile         = 0
    PPInfo%GraphUnit        = 0

    ! If not separate times in separate files, cycle.
    If (.not.PPInfoReq%TSeparateFile) Cycle

    ! Abbreviation for T-grid.
    TGrid => Grids%TGrids(PPInfoReq%iTGrid)

    ! nT.
    Call FirstTAfterT(TGrid, StartTime,  .false., DummyT, iT1)
    Call LastTBeforeT(TGrid, MaxEndTime, .false., DummyT, iT2)
    If (iT2 < iT1) Then
      nT = 1
    Else If (iT1 == -Huge(iT1)/3 .or. iT2 == Huge(iT2)/3) Then
      nT = Huge(nT)/3
    Else
      nT = iT2 - iT1 + 1
    End If
    If (IsInfFuture(ReqMemoryTime)) Then
      PPInfo%nT = nT
    Else
      PPInfo%nT = Min(ReqMemoryTime / TGrid%sdT + 1, nT)
    End If
    If (PPInfo%nT > MaxPPInfoOutputFiles) Then
      Call Message(                                                                     &
             'FATAL ERROR: too many output files/windows required at the same time ' // &
             'for the set of particle/puff information "'                            // &
             Trim(PPInfoReq%Name)                                                    // &
             '"',                                                                       &
             3                                                                          &
           )
    End If

    ! Calculate index of next time in time grid which is >= StartTime and adjust to ensure the range iTL to
    ! iTL + nT - 1 lies in the time grid. This will be the earliest time stored from now on.
    Call FirstTAfterT(TGrid, StartTime, .false., DummyT, PPInfo%iTL)
    PPInfo%iTL = Min(PPInfo%iTL, TGrid%nT - PPInfo%nT + 1)
    PPInfo%iTL = Max(PPInfo%iTL, 1)

  End Do

  Results%nPPInfos = Reqs%nPPInfoReqs

End Subroutine PrepareResults

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateResults(Time, MaxPPTimeLag, Grids, Reqs, Results)
! Updates the collection of results by deciding which times will be stored from now on.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In)            :: Time         ! Time.
  Type(ShortTime_), Intent(In)            :: MaxPPTimeLag ! The maximum possible interval which can occur at
                                                          ! any time between the trailing time-edge and the
                                                          ! nominal time of the particles/puffs.
  Type(Grids_),     Intent(In),    Target :: Grids        ! Collection of grids.
  Type(Reqs_),      Intent(In),    Target :: Reqs         ! Collection of requirements.
  Type(Results_),   Intent(InOut), Target :: Results      ! Collection of results.
  ! Locals:
  Type(TGrid_),     Pointer :: TGrid     !} Abbreviations for grids, requirements and results.
  Type(FieldReq_),  Pointer :: FieldReq  !}
  Type(PPInfoReq_), Pointer :: PPInfoReq !}
  Type(Field_),     Pointer :: Field     !}
  Type(PPInfo_),    Pointer :: PPInfo    !}
  Integer                   :: iTL       ! Updated index of lowest time level stored.
  Integer                   :: iTLL      !} Time indices delimiting any overlap between times stored up to now
  Integer                   :: iTUU      !} and times stored from now on.
  Integer                   :: iTAL      !] Time indices in a field array delimiting any overlap between times
  Integer                   :: iTAU      !] stored up to now and times stored from now on.
  Integer                   :: i         !} Loop indices.
  Integer                   :: j         !}
  Type(ShortTime_)          :: DummyT    ! Dummy variable needed for an argument list.

  ! Fields.
  Do i = 1, Reqs%nFieldReqs

    ! Abbreviations for requirement and result.
    FieldReq => Reqs%FieldReqs(i)
    Field    => Results%Fields(i)

    ! If no t-grid, cycle.
    If (FieldReq%iTGrid == 0) Cycle

    ! Abbreviation for time grid.
    TGrid => Grids%TGrids(FieldReq%iTGrid)

    ! Calculate index of next time in time grid which is >= Time - MaxPPTimeLag and adjust to ensure iTL does
    ! not decrease and the range iTL to iTL + nT - 1 lies in the time grid. This will be the earliest time
    ! stored from now on. The fields processed over the ensemble of cases or with TAcross and not
    ! TSeparateFile are an exception because here all data times need to be retained in memory.
    If (FieldReq%Ensemble) Then
      iTL = Field%iTL
    Else If (.not.FieldReq%TSeparateFile .and. FieldReq%TAcross) Then
      iTL = Field%iTL
    Else
      Call FirstTAfterT(TGrid, Time - MaxPPTimeLag, .false., DummyT, iTL)
      iTL = Max(iTL, Field%iTL)
      iTL = Min(iTL, TGrid%nT - Field%nT + 1)
      iTL = Max(iTL, 1)
    End If

    ! Time indices delimiting any overlap between times stored up to now and times stored from now on.
    iTLL = Max(iTL, Field%iTL)
    iTUU = Min(iTL, Field%iTL) + Field%nT - 1

    ! There is an overlap.
    If (iTLL <= iTUU) Then

      ! Time indices in the field array delimiting any overlap between times stored up to now and times stored
      ! from now on.
      iTAL = CalciTA(Field, iTLL)
      iTAU = CalciTA(Field, iTUU)

      ! Set Field to zero at non-overlap times.
      If (iTAL <= iTAU) Then
        Do j = 1, iTAL - 1
          If (Associated(Field%Std   )) Field%Std(:, :, :, :, :, j, :) = 0.0
          If (Associated(Field%P64   )) Field%P64(:, :, :, :, :, j, :) = 0.0
          If (Associated(Field%T     )) Field%T  (:, :, :, :, :, j, :) = ZeroTime()
          If (Associated(Field%MaxStd)) Field%MaxStd(:, :, :, :, j, :) = -Huge(Field%MaxStd)
          If (Associated(Field%MinStd)) Field%MinStd(:, :, :, :, j, :) =  Huge(Field%MinStd)
          If (Associated(Field%S     )) Field%S     (:, :, :, :, j, :) = -Huge(Field%S)
        End Do
        Do j = iTAU + 1, Field%nT
          If (Associated(Field%Std   )) Field%Std(:, :, :, :, :, j, :) = 0.0
          If (Associated(Field%P64   )) Field%P64(:, :, :, :, :, j, :) = 0.0
          If (Associated(Field%T     )) Field%T  (:, :, :, :, :, j, :) = ZeroTime()
          If (Associated(Field%MaxStd)) Field%MaxStd(:, :, :, :, j, :) = -Huge(Field%MaxStd)
          If (Associated(Field%MinStd)) Field%MinStd(:, :, :, :, j, :) =  Huge(Field%MinStd)
          If (Associated(Field%S     )) Field%S     (:, :, :, :, j, :) = -Huge(Field%S)
        End Do
      Else
        Do j = iTAU + 1, iTAL - 1
          If (Associated(Field%Std   )) Field%Std(:, :, :, :, :, j, :) = 0.0
          If (Associated(Field%P64   )) Field%P64(:, :, :, :, :, j, :) = 0.0
          If (Associated(Field%T     )) Field%T  (:, :, :, :, :, j, :) = ZeroTime()
          If (Associated(Field%MaxStd)) Field%MaxStd(:, :, :, :, j, :) = -Huge(Field%MaxStd)
          If (Associated(Field%MinStd)) Field%MinStd(:, :, :, :, j, :) =  Huge(Field%MinStd)
          If (Associated(Field%S     )) Field%S     (:, :, :, :, j, :) = -Huge(Field%S)
        End Do
      End If

    ! There is no overlap.
    Else

      ! Set Field to zero at all times.
      Do j = 1, Field%nT
        If (Associated(Field%Std   )) Field%Std(:, :, :, :, :, j, :) = 0.0
        If (Associated(Field%P64   )) Field%P64(:, :, :, :, :, j, :) = 0.0
        If (Associated(Field%T     )) Field%T  (:, :, :, :, :, j, :) = ZeroTime()
        If (Associated(Field%MaxStd)) Field%MaxStd(:, :, :, :, j, :) = -Huge(Field%MaxStd)
        If (Associated(Field%MinStd)) Field%MinStd(:, :, :, :, j, :) =  Huge(Field%MinStd)
        If (Associated(Field%S     )) Field%S     (:, :, :, :, j, :) = -Huge(Field%S)
      End Do

    End If

    ! Update Field%iTL.
    Field%iTL = iTL

  End Do

  ! Sets of particle/puff information.
  Do i = 1, Reqs%nPPInfoReqs

    ! Abbreviations for requirement and result.
    PPInfoReq => Reqs%PPInfoReqs(i)
    PPInfo    => Results%PPInfos(i)

    ! If not separate times in separate files, cycle.
    If (.not.PPInfoReq%TSeparateFile) Cycle

    ! Abbreviation for T-grid.
    TGrid => Grids%TGrids(PPInfoReq%iTGrid)

    ! Calculate index of next time in time grid which is >= Time and adjust to ensure the range iTL to iTL +
    ! nT - 1 lies in the time grid. This will be the earliest time stored from now on.
    Call FirstTAfterT(TGrid, Time, .false., DummyT, PPInfo%iTL)
    PPInfo%iTL = Min(PPInfo%iTL, TGrid%nT - PPInfo%nT + 1)
    PPInfo%iTL = Max(PPInfo%iTL, 1)

  End Do

End Subroutine UpdateResults

!-------------------------------------------------------------------------------------------------------------

Function CalciTA(Field, iT) Result(iTA)
! Calculates the time index in a field array at which the results for a given time are stored.

  Implicit None
  ! Argument list:
  Type(Field_), Intent(In) :: Field ! The field of interest.
  Integer,      Intent(In) :: iT    ! The time level of interest.
  ! Function result:
  Integer :: iTA ! The time index in the field array at which the results are stored.

  If (iT < Field%iTL .or. iT > Field%iTL + Field%nT - 1) Then
    Call Message('UNEXPECTED FATAL ERROR in CalciTA', 4)
  End If

  iTA = Mod(iT, Field%nT)
  If (iTA <= 0) iTA = iTA + Field%nT

End Function CalciTA

!-------------------------------------------------------------------------------------------------------------

Function CalciLast(FieldReq, iComponent, iRange) Result(iLast)
! Calculates last index in a field array corresponding to a given component and size range.

  Implicit None
  ! Argument list:
  Type(FieldReq_), Intent(In) :: FieldReq   ! The field requirement corresponding to the field of interest.
  Integer,         Intent(In) :: iComponent ! Index of component.
  Integer,         Intent(In) :: iRange     ! Index of size range.
  ! Function result:
  Integer :: iLast ! The last index in the field array.

# ifdef ExtraChecks
    If (                                                         &
      iComponent < 1 .or. iComponent > FieldReq%nComponents .or. &
      iRange     < 1 .or. iRange     > FieldReq%nRanges          &
    ) Then
      Call Message('UNEXPECTED FATAL ERROR in CalciLast', 4)
    End If
# endif

  iLast = FieldReq%nComponents * (iRange - 1) + iComponent

End Function CalciLast

!-------------------------------------------------------------------------------------------------------------

Function CalciComponent(FieldReq, iLast) Result(iComponent)
! Calculates index of component corresponding to a given last index in a field array.

  Implicit None
  ! Argument list:
  Type(FieldReq_), Intent(In) :: FieldReq ! The field requirement corresponding to the field of interest.
  Integer,         Intent(In) :: iLast    ! The last index in the field array.
  ! Function result:
  Integer :: iComponent ! Index of component.

# ifdef ExtraChecks
    If (iLast < 1 .or. iLast > FieldReq%nComponents * FieldReq%nRanges) Then
      Call Message('UNEXPECTED FATAL ERROR in CalciComponent', 4)
    End If
# endif

  ! This is written this way to look similar to generalisations for more indices.
  iComponent = Mod(                   &
                 iLast - 1,           &
                 FieldReq%nComponents &
               ) + 1

End Function CalciComponent

!-------------------------------------------------------------------------------------------------------------

Function CalciRange(FieldReq, iLast) Result(iRange)
! Calculates index of size range corresponding to a given last index in a field array.

  Implicit None
  ! Argument list:
  Type(FieldReq_), Intent(In) :: FieldReq ! The field requirement corresponding to the field of interest.
  Integer,         Intent(In) :: iLast    ! The last index in the field array.
  ! Function result:
  Integer :: iRange ! Index of size range.

# ifdef ExtraChecks
    If (iLast < 1 .or. iLast > FieldReq%nComponents * FieldReq%nRanges) Then
      Call Message('UNEXPECTED FATAL ERROR in CalciRange', 4)
    End If
# endif

  ! This is written this way to look similar to generalisations for more indices.
  iRange = Mod(                                  &
             (iLast - 1) / FieldReq%nComponents, &
             FieldReq%nRanges                    &
           ) + 1

End Function CalciRange

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcReqInfoNoT(Reqs, ReqInfo)
! Calculates requirements for results which don't depend on time due to the underlying quantity not depending
! on time.

! Included here are:
! fields where the underlying quantity doesn't depend on time and with no explict processing step.

! Note that avoiding contributions being made more than once is the responsibility of the calling routine.

  Implicit None
  ! Argument list:
  Type(Reqs_),    Intent(In), Target :: Reqs    ! Collection of requirements.
  Type(ReqInfo_), Intent(Out)        :: ReqInfo ! Information on what requirements are required.
  ! Locals:
  Type(FieldReq_), Pointer :: FieldReq ! Abbreviation for field requirement.
  Integer                  :: iField   ! Field requirement index.

  ReqInfo%NoReqs = .true.

  Do iField = 1, Reqs%nFieldReqs

    FieldReq => Reqs%FieldReqs(iField)

    If (FieldReq%iExt /= 0) Cycle

    If (FieldReq%iLastProcess                >  0) Cycle
    If (Scan(QInfo(FieldReq%iQuantity), 'T') /= 0) Cycle

    If (ReqInfo%NoReqs) Then

      If (ReqInfo%nFields == MaxFieldReqsPerReqInfo) Then
        Call Message(                                                                              &
               'FATAL ERROR: The requirements are too numerous. In this regard only '           // &
               'field requirements for quantities with no time dependence count. The limit is ' // &
               'controlled by MaxFieldReqsPerReqInfo in the global parameters module.',            &
               3                                                                                   &
             )
      End If

      ReqInfo%NoReqs                  = .false.
      ReqInfo%nPPInfos                = 0
      ReqInfo%nFields                 = 1
      ReqInfo%iField(ReqInfo%nFields) = iField
      ReqInfo%iT    (ReqInfo%nFields) = 1

    Else If (.not.ReqInfo%NoReqs) Then

      If (ReqInfo%nFields == MaxFieldReqsPerReqInfo) Then
        Call Message(                                                                              &
               'FATAL ERROR: The requirements are too numerous. In this regard only '           // &
               'field requirements for quantities with no time dependence count. The limit is ' // &
               'controlled by MaxFieldReqsPerReqInfo in the global parameters module.',            &
               3                                                                                   &
             )
      End If

      ReqInfo%nFields                 = ReqInfo%nFields + 1
      If (ReqInfo%nFields > MaxFieldReqsPerReqInfo) Then
        Call Message('FATAL ERROR: MaxFieldReqsPerReqInfo is not large enough', 3)
      End If
      ReqInfo%iField(ReqInfo%nFields) = iField
      ReqInfo%iT    (ReqInfo%nFields) = 1

    End If

  End Do

End Subroutine CalcReqInfoNoT

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcReqInfoT(                  &
             IncFields, IncPPInfos,       &
             LEFOType, Sync,              &
             IncTime1, Time1, Time2,      &
             MaxPPTimeLead, MaxPPTimeLag, &
             Grids, Reqs,                 &
             ReqInfo                      &
           )
! Calculates requirements for results in a given time range which require contributions at discrete times.

! Included here are:
! (i)  fields of type LEFOType, with FieldReq%Sync = Sync, with a T-Grid, with no S-grid, with no intrinsic
!      T-averaging/integrating, and with no explict processing step,
! (ii) sets of particle/puff information (if LEFOType = 'L' only) with FieldReq%Sync = Sync and with a T-grid.

! Note that avoiding contributions being made more than once is the responsibility of the calling routine.

  Implicit None
  ! Argument list:
  Logical,          Intent(In)  :: IncFields
  Logical,          Intent(In)  :: IncPPInfos
  Character(1),     Intent(In)  :: LEFOType
  Logical,          Intent(In)  :: Sync
  Logical,          Intent(In)  :: IncTime1
  Type(ShortTime_), Intent(In)  :: Time1
  Type(ShortTime_), Intent(In)  :: Time2
  Type(ShortTime_), Intent(In)  :: MaxPPTimeLead
  Type(ShortTime_), Intent(In)  :: MaxPPTimeLag
  Type(Grids_),     Intent(In)  :: Grids
  Type(Reqs_),      Intent(In)  :: Reqs
  Type(ReqInfo_),   Intent(Out) :: ReqInfo(:)
  ! IncFields     :} Indicates fields and sets of particle/puff information are to be included.
  ! IncPPInfos    :}
  ! LEFOType      :: Type of field requirements (L, E, F or O) sought.
  ! Sync          :: Indicates whether to calculate sync time or non sync time requirements.
  ! IncTime1      :: Indicates the time range is [Time1, Time2]. Otherwise its (Time1, Time2].
  ! Time1         :} Times defining time range.
  ! Time2         :}
  ! MaxPPTimeLead :: The maximum possible interval which can occur at any time between the leading time-edge
  !                  and the nominal time of the particles/puffs.
  ! MaxPPTimeLag  :: The maximum possible interval which can occur at any time between the trailing time-edge
  !                  and the nominal time of the particles/puffs.
  ! Grids         :: Collection of grids.
  ! Reqs          :: Collection of requirements.
  ! ReqInfo       :: Information on what requirements are required.
  ! Locals:
  Logical          :: IncTime1L ! Local copy of IncTime1.
  Type(ShortTime_) :: Time1L    ! Local copy of Time1.
  Type(ShortTime_) :: Time2L    ! Local copy of Time2.
  Type(ShortTime_) :: LastTime  ! Last time found so far at which contributions required.
  Integer          :: i         ! Loop index.

  ! $$ improve efficiency by not accounting for MaxPPTimeLag/Lead for PPInfos
  !    - possibly other calcreqinfo routines too.

  IncTime1L = IncTime1
  If (LEFOType == 'L') Then
    Time1L = Time1 - MaxPPTimeLag
    Time2L = Time2 + MaxPPTimeLead
  Else
    Time1L = Time1
    Time2L = Time2
  End If

  LastTime = Time1L

  ! Find all requirements within the time range.
  Do i = 1, Size(ReqInfo) + 1

    If (i > Size(ReqInfo)) Then
      Call Message(                                                                                         &
             'FATAL ERROR: The times at which requirements need to be computed are too numerous or too ' // &
             'densely spaced. In this regard only requirements with a time grid, with no travel-time '   // &
             'grid and with no intrinsic time averaging/integrating count. The limit is controlled '     // &
             'by MaxReqInfoTimes in the global parameters module.',                                         &
             3                                                                                              &
           )
    End If

    Call CalcNextReqInfoT(IncFields, IncPPInfos, LEFOType, Sync, IncTime1L, LastTime, Grids, Reqs, ReqInfo(i))

    If (ReqInfo(i)%NoReqs) Exit

    If (Time2L < ReqInfo(i)%Time) Then
      ReqInfo(i)%NoReqs = .true.
      Exit
    End If

    LastTime  = ReqInfo(i)%Time
    IncTime1L = .false.

  End Do

End Subroutine CalcReqInfoT

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcReqInfoS(                               &
             Sync,                                     &
             IncTime1, Time1, Time2,                   &
             IncTravelTime1, TravelTime1, TravelTime2, &
             MaxPPTimeLead, MaxPPTimeLag,              &
             Grids, Reqs,                              &
             ReqInfo                                   &
           )
! Calculates requirements for results in a given time and travel-time range which require contributions at
! discrete travel-times.

! Included here are:
! fields of type L, with FieldReq%Sync = Sync, with an S-grid, and with no explict processing step.

! Note that avoiding contributions being made more than once is the responsibility of the calling routine.

  Implicit None
  ! Argument list:
  Logical,          Intent(In)  :: Sync
  Logical,          Intent(In)  :: IncTime1
  Type(ShortTime_), Intent(In)  :: Time1
  Type(ShortTime_), Intent(In)  :: Time2
  Logical,          Intent(In)  :: IncTravelTime1
  Type(ShortTime_), Intent(In)  :: TravelTime1
  Type(ShortTime_), Intent(In)  :: TravelTime2
  Type(ShortTime_), Intent(In)  :: MaxPPTimeLead
  Type(ShortTime_), Intent(In)  :: MaxPPTimeLag
  Type(Grids_),     Intent(In)  :: Grids
  Type(Reqs_),      Intent(In)  :: Reqs
  Type(ReqInfo_),   Intent(Out) :: ReqInfo(:)
  ! Sync           :: Indicates whether to calculate sync time or non sync time requirements.
  ! IncTime1       :: Indicates the time range is [Time1, Time2]. Otherwise its (Time1, Time2].
  ! Time1          :} Times defining time range.
  ! Time2          :}
  ! IncTravelTime1 :: Indicates the travel-time range is [TravelTime1, TravelTime2]. Otherwise its
  !                   (TravelTime1, TravelTime2].
  ! TravelTime1    :] Travel-times defining travel-time range.
  ! TravelTime2    :]
  ! MaxPPTimeLead  :: The maximum possible interval which can occur at any time between the leading time-edge
  !                   and the nominal time of the particles/puffs.
  ! MaxPPTimeLag   :: The maximum possible interval which can occur at any time between the trailing time-edge
  !                   and the nominal time of the particles/puffs.
  ! Grids          :: Collection of grids.
  ! Reqs           :: Collection of requirements.
  ! ReqInfo        :: Information on what requirements are required.
  ! Locals:
  Logical          :: IncTravelTime1L ! Local copy of IncTravelTime1.
  Type(ShortTime_) :: Time1L          ! Local copy of Time1.
  Type(ShortTime_) :: Time2L          ! Local copy of Time2.
  Type(ShortTime_) :: LastTravelTime  ! Last travel-time found so far at which contributions required.
  Integer          :: i               ! Loop index.

  IncTravelTime1L = IncTravelTime1
  Time1L          = Time1 - MaxPPTimeLag
  Time2L          = Time2 + MaxPPTimeLead

  LastTravelTime = TravelTime1

  ! Find all requirements within the travel-time range.
  Do i = 1, Size(ReqInfo) + 1

    If (i > Size(ReqInfo)) Then
      Call Message(                                                                                         &
             'FATAL ERROR: The travel-times at which requirements need to be computed are too numerous ' // &
             'or too densely spaced. In this regard only requirements with a travel-time grid '          // &
             'count. The limit is controlled by MaxReqInfoTravelTimes in the global parameters module.',    &
             3                                                                                              &
           )
    End If

    Call CalcNextReqInfoS(                                            &
           Sync,                                                      &
           IncTime1, Time1L, Time2L, IncTravelTime1L, LastTravelTime, &
           Grids, Reqs,                                               &
           ReqInfo(i)                                                 &
         )

    If (ReqInfo(i)%NoReqs) Exit

    If (TravelTime2 < ReqInfo(i)%Time) Then
      ReqInfo(i)%NoReqs = .true.
      Exit
    End If

    LastTravelTime  = ReqInfo(i)%Time
    IncTravelTime1L = .false.

  End Do

End Subroutine CalcReqInfoS

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcReqInfoEveryT(             &
             Sync,                        &
             Time1, Time2,                &
             MaxPPTimeLead, MaxPPTimeLag, &
             Grids, Reqs,                 &
             ReqInfo                      &
           )
! Calculates requirements for results in a given time range which require contributions at every time in some
! time interval.

! Included here are:
! (i)   fields of type L, with FieldReq%Sync = Sync, with a T-grid, with no S-grid, with intrinsic
!       T-averaging/integrating, and with no explict processing step,
! (ii)  fields of type L, with FieldReq%Sync = Sync, with no T-grid, with no S-grid, with intrinsic
!       T-averaging/integrating over all T, and with no explict processing step,
! (iii) sets of particle/puff information with no T-grid.
! Note that fields with intrinsic T-averaging/integrating and with no explict processing step must have
! non-overlapping averaging/integrating periods and Process%nT = 1 with Process%T = Process%dT.

! Note that avoiding contributions being made more than once is the responsibility of the calling routine.

  Implicit None
  ! Argument list:
  Logical,                  Intent(In)  :: Sync
  Type(ShortTime_),         Intent(In)  :: Time1
  Type(ShortTime_),         Intent(In)  :: Time2
  Type(ShortTime_),         Intent(In)  :: MaxPPTimeLead
  Type(ShortTime_),         Intent(In)  :: MaxPPTimeLag
  Type(Grids_),     Target, Intent(In)  :: Grids
  Type(Reqs_),      Target, Intent(In)  :: Reqs
  Type(ReqInfo_),           Intent(Out) :: ReqInfo
  ! Sync          :: Indicates whether to calculate sync time or non sync time requirements.
  ! Time1         :} Times defining time range.
  ! Time2         :}
  ! MaxPPTimeLead :: The maximum possible interval which can occur at any time between the leading time-edge
  !                  and the nominal time of the particles/puffs.
  ! MaxPPTimeLag  :: The maximum possible interval which can occur at any time between the trailing time-edge
  !                  and the nominal time of the particles/puffs.
  ! Grids         :: Collection of grids.
  ! Reqs          :: Collection of requirements.
  ! ReqInfo       :: Information on what requirements are required.
  ! Locals:
  Type(TGrid_),     Pointer :: TGrid     !} Abbreviations for grids and requirements.
  Type(FieldReq_),  Pointer :: FieldReq  !}
  Type(PPInfoReq_), Pointer :: PPInfoReq !}
  Integer                   :: iField    !] Indices of requirements.
  Integer                   :: iPPInfo   !]
  Integer                   :: iT        ! Index of time in time grid.
  Type(ShortTime_)          :: Time1L    ! Local copy of Time1.
  Type(ShortTime_)          :: Time2L    ! Local copy of Time2.
  Type(ShortTime_)          :: T         ! Time of possible future requirements.
  Type(ShortTime_)          :: LastT     ! Last value of T.

# ifdef ExtraChecks
    If (Time2 <= Time1) Call Message('UNEXPECTED FATAL ERROR in CalcReqInfoEveryT', 4)
# endif

  Time1L = Time1 - MaxPPTimeLag
  Time2L = Time2 + MaxPPTimeLead

  ReqInfo%NoReqs = .true.
  ReqInfo%Time   = Time2L

  ReqInfo%nFields  = 0
  ReqInfo%nPPInfos = 0

  ! Fields with a T-grid.
  Do iField = 1, Reqs%nFieldReqs

    FieldReq => Reqs%FieldReqs(iField)

    If (FieldReq%iExt /= 0) Cycle

    If (FieldReq%iLastProcess  >    0   ) Cycle
    If (FieldReq%LEFOPQRType   /=   'L' ) Cycle
    If (FieldReq%Sync        .neqv. Sync) Cycle
    If (FieldReq%iTGrid        ==   0   ) Cycle
    If (FieldReq%iSGrid        /=   0   ) Cycle
    If (.not. FieldReq%IntrinsicTAv     ) Cycle

    TGrid => Grids%TGrids(FieldReq%iTGrid)

    LastT = Time1L

    Do

      Call FirstTAfterT(TGrid, LastT, .true., T, iT)
      If (iT == Huge(iT)/3) Cycle

      If (T - FieldReq%sAvT >= Time2L) Cycle

      If (ReqInfo%nFields == MaxFieldReqsPerReqInfo) Then
        Call Message(                                                                                     &
               'FATAL ERROR: The times at which requirements need to be computed are too numerous or ' // &
               'too densely spaced. In this regard only field requirements with no travel-time grid '  // &
               'and with intrinsic time averaging/integrating count, but the same time counts more '   // &
               'than once if needed for more than one field requirement and a field averaged or '      // &
               'integrated over all time counts as a time at infinity. The limit is controlled by '    // &
               'MaxFieldReqsPerReqInfo in the global parameters module.',                                 &
               3                                                                                          &
             )
      End If

      ReqInfo%NoReqs                          = .false.
      ReqInfo%nFields                         = ReqInfo%nFields + 1
      If (ReqInfo%nFields > MaxFieldReqsPerReqInfo) Then
        Call Message('FATAL ERROR: MaxFieldReqsPerReqInfo is not large enough', 3)
      End If
      ReqInfo%iField        (ReqInfo%nFields) = iField
      ReqInfo%iT            (ReqInfo%nFields) = iT
      ReqInfo%FieldStartTime(ReqInfo%nFields) = T - FieldReq%sAvT
      ReqInfo%FieldEndTime  (ReqInfo%nFields) = T

      LastT = T

    End Do

  End Do

  ! Fields with no T-grid.
  Do iField = 1, Reqs%nFieldReqs

    FieldReq => Reqs%FieldReqs(iField)

    If (FieldReq%iExt /= 0) Cycle

    If (FieldReq%iLastProcess  >    0                                 ) Cycle
    If (FieldReq%LEFOPQRType   /=   'L'                               ) Cycle
    If (FieldReq%Sync        .neqv. Sync                              ) Cycle
    If (FieldReq%iTGrid        /=   0                                 ) Cycle
    If (FieldReq%iSGrid        /=   0                                 ) Cycle
    If (.not. (FieldReq%IntrinsicTAv .and. IsInfFuture(FieldReq%sAvT))) Cycle

    If (ReqInfo%nFields == MaxFieldReqsPerReqInfo) Then
      Call Message(                                                                                     &
             'FATAL ERROR: The times at which requirements need to be computed are too numerous or ' // &
             'too densely spaced. In this regard only field requirements with no travel-time grid '  // &
             'and with intrinsic time averaging/integrating count, but the same time counts more '   // &
             'than once if needed for more than one field requirement and a field averaged or '      // &
             'integrated over all time counts as a time at infinity. The limit is controlled by '    // &
             'MaxFieldReqsPerReqInfo in the global parameters module.',                                 &
             3                                                                                          &
           )
    End If

    ReqInfo%NoReqs                          = .false.
    ReqInfo%nFields                         = ReqInfo%nFields + 1
    If (ReqInfo%nFields > MaxFieldReqsPerReqInfo) Then
      Call Message('FATAL ERROR: MaxFieldReqsPerReqInfo is not large enough', 3)
    End If
    ReqInfo%iField        (ReqInfo%nFields) = iField
    ReqInfo%iT            (ReqInfo%nFields) = 1
    ReqInfo%FieldStartTime(ReqInfo%nFields) = InfPastShortTime()
    ReqInfo%FieldEndTIme  (ReqInfo%nFields) = InfFutureShortTime()

  End Do

  ! Sets of particle/puff information.
  Do iPPInfo = 1, Reqs%nPPInfoReqs

    PPInfoReq => Reqs%PPInfoReqs(iPPInfo)

    If (PPInfoReq%iTGrid == 0) Then

      ReqInfo%NoReqs                    = .false.
      ReqInfo%nPPInfos                  = ReqInfo%nPPInfos + 1
      If (ReqInfo%nPPInfos > MaxPPInfoReqsPerReqInfo) Then
        Call Message('FATAL ERROR: MaxPPInfoReqsPerReqInfo is not large enough', 3)
      End If
      ReqInfo%iPPInfo(ReqInfo%nPPInfos) = iPPInfo

    End If

  End Do

End Subroutine CalcReqInfoEveryT

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcNextReqInfoT(IncFields, IncPPInfos, LEFOType, Sync, IncTime, Time, Grids, Reqs, ReqInfo)
! Calculates the next time for which there are requirements which require contributions at discrete times, and
! calculates what requirements require contributions at that time.

! Included here are
! (i)  fields of type LEFOType, with FieldReq%Sync = Sync, with a T-Grid, with no S-grid, with no intrinsic
!      T-averaging/integrating, and with no explict processing step,
! (ii) sets of particle/puff information (if LEFOType = 'L' only) with FieldReq%Sync = Sync and with a T-grid.

  Implicit None
  ! Argument list:
  Logical,          Intent(In)         :: IncFields  !} Indicates fields and sets of particle/puff information
  Logical,          Intent(In)         :: IncPPInfos !} are to be included.
  Character(1),     Intent(In)         :: LEFOType   ! Type of field requirements (L, E, F or O) sought.
  Logical,          Intent(In)         :: Sync       ! Indicates whether to calculate sync time or non sync
                                                     ! time requirements.
  Logical,          Intent(In)         :: IncTime    ! Indicates the times of interest include Time.
  Type(ShortTime_), Intent(In)         :: Time       ! Earliest time of interest.
  Type(Grids_),     Intent(In), Target :: Grids      ! Collection of grids.
  Type(Reqs_),      Intent(In), Target :: Reqs       ! Collection of requirements.
  Type(ReqInfo_),   Intent(Out)        :: ReqInfo    ! Information on what requirements are required next.
  ! Locals:
  Type(TGrid_),     Pointer :: TGrid     !} Abbreviations for grids and requirements.
  Type(FieldReq_),  Pointer :: FieldReq  !}
  Type(PPInfoReq_), Pointer :: PPInfoReq !}
  Integer                   :: iField    !] Indices of requirements.
  Integer                   :: iPPInfo   !]
  Integer                   :: iT        ! Index of time in time grid.
  Type(ShortTime_)          :: T         ! Time of possible future requirements.

  ReqInfo%NoReqs = .true.
  ReqInfo%Time   = InfFutureShortTime()

  ! Fields.
  If (IncFields) Then

    Do iField = 1, Reqs%nFieldReqs

      FieldReq => Reqs%FieldReqs(iField)

      If (FieldReq%iExt /= 0) Cycle

      If (FieldReq%iLastProcess   >    0       ) Cycle
      If (FieldReq%LEFOPQRType    /=   LEFOType) Cycle
      If (FieldReq%Sync         .neqv. Sync    ) Cycle
      If (FieldReq%iTGrid         ==   0       ) Cycle
      If (FieldReq%iSGrid         /=   0       ) Cycle
      If (FieldReq%IntrinsicTAv                ) Cycle

      TGrid => Grids%TGrids(FieldReq%iTGrid)

      Call FirstTAfterT(TGrid, Time, .not.IncTime, T, iT)
      If (iT == Huge(iT)/3) Cycle

      If (T < ReqInfo%Time .or. ReqInfo%NoReqs) Then

        ReqInfo%NoReqs                  = .false.
        ReqInfo%Time                    = T
        ReqInfo%nPPInfos                = 0
        ReqInfo%nFields                 = 1
        ReqInfo%iField(ReqInfo%nFields) = iField
        ReqInfo%iT    (ReqInfo%nFields) = iT

      Else If (ReqInfo%Time == T .and. .not.ReqInfo%NoReqs) Then

        ReqInfo%nFields                 = ReqInfo%nFields + 1
        If (ReqInfo%nFields > MaxFieldReqsPerReqInfo) Then
          Call Message('FATAL ERROR: MaxFieldReqsPerReqInfo is not large enough', 3)
        End If
        ReqInfo%iField(ReqInfo%nFields) = iField
        ReqInfo%iT    (ReqInfo%nFields) = iT

      End If

    End Do

  End If

  If (LEFOType /= 'L') Return

  ! Sets of particle/puff information.
  If (IncPPInfos) Then

    Do iPPInfo = 1, Reqs%nPPInfoReqs

      PPInfoReq => Reqs%PPInfoReqs(iPPInfo)

      If (PPInfoReq%Sync   .neqv. Sync) Cycle
      If (PPInfoReq%iTGrid   ==   0   ) Cycle

      TGrid => Grids%TGrids(PPInfoReq%iTGrid)

      Call FirstTAfterT(TGrid, Time, .not.IncTime, T, iT)
      If (iT == Huge(iT)/3) Cycle

      If (T < ReqInfo%Time .or. ReqInfo%NoReqs) Then

        ReqInfo%NoReqs                     = .false.
        ReqInfo%Time                       = T
        ReqInfo%nFields                    = 0
        ReqInfo%nPPInfos                   = 1
        ReqInfo%iPPInfo (ReqInfo%nPPInfos) = iPPInfo
        ReqInfo%iPPInfoT(ReqInfo%nPPInfos) = iT

      Else If (ReqInfo%Time == T .and. .not.ReqInfo%NoReqs) Then

        ReqInfo%nPPInfos                   = ReqInfo%nPPInfos + 1
        If (ReqInfo%nPPInfos > MaxPPInfoReqsPerReqInfo) Then
          Call Message('FATAL ERROR: MaxPPInfoReqsPerReqInfo is not large enough', 3)
        End If
        ReqInfo%iPPInfo (ReqInfo%nPPInfos) = iPPInfo
        ReqInfo%iPPInfoT(ReqInfo%nPPInfos) = iT

      End If

    End Do

  End If

End Subroutine CalcNextReqInfoT

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcNextReqInfoS(Sync, IncTime1, Time1, Time2, IncTravelTime, TravelTime, Grids, Reqs, ReqInfo)
! Calculates the next travel-time for which there are requirements which require contributions at discrete
! travel-times in a given time range, and calculates what requirements require contributions at that
! travel-time.

! Included here are:
! fields of type L, with FieldReq%Sync = Sync, with an S-grid, and with no explict processing step.

  Implicit None
  ! Argument list:
  Logical,                  Intent(In)  :: Sync          ! Indicates whether to calculate sync time or non
                                                         ! sync time requirements.
  Logical,                  Intent(In)  :: IncTime1      ! Indicates the time range is [Time1, Time2].
                                                         ! Otherwise its (Time1, Time2].
  Type(ShortTime_),         Intent(In)  :: Time1         !} Times defining time range.
  Type(ShortTime_),         Intent(In)  :: Time2         !}
  Logical,                  Intent(In)  :: IncTravelTime ! Indicates the travel-times of interest include
                                                         ! TravelTime.
  Type(ShortTime_),         Intent(In)  :: TravelTime    ! Earliest travel-time of interest.
  Type(Grids_),     Target, Intent(In)  :: Grids         ! Collection of grids.
  Type(Reqs_),      Target, Intent(In)  :: Reqs          ! Collection of requirements.
  Type(ReqInfo_),           Intent(Out) :: ReqInfo       ! Information on what requirements are required next.
  ! Locals:
  Type(TGrid_),    Pointer :: SGrid    !} Abbreviations for grids and requirements.
  Type(FieldReq_), Pointer :: FieldReq !}
  Integer                  :: iField   ! Index of field requirement.
  Integer                  :: iS       ! Index of travel-time in travel-time grid.
  Type(ShortTime_)         :: S        ! Travel-time of possible future requirements.

  ReqInfo%NoReqs = .true.
  ReqInfo%Time   = InfFutureShortTime(Interval = .true.)

  ! Fields.
  Do iField = 1, Reqs%nFieldReqs

    FieldReq => Reqs%FieldReqs(iField)

    If (FieldReq%iExt /= 0) Cycle

    If (FieldReq%iLastProcess  >    0   ) Cycle
    If (FieldReq%LEFOPQRType   /=   'L' ) Cycle
    If (FieldReq%Sync        .neqv. Sync) Cycle
    If (FieldReq%iSGrid        ==   0   ) Cycle

    SGrid => Grids%TGrids(FieldReq%iSGrid)

    ! Use time1, time2, inctime1 restriction $$

    Call FirstTAfterT(SGrid, TravelTime, .not.IncTravelTime, S, iS)
    If (iS == Huge(iS)/3) Cycle

    If (S < ReqInfo%Time .or. ReqInfo%NoReqs) Then

      ReqInfo%NoReqs                  = .false.
      ReqInfo%Time                    = S
      ReqInfo%nPPInfos                = 0
      ReqInfo%nFields                 = 1
      ReqInfo%iField(ReqInfo%nFields) = iField
      ReqInfo%iS    (ReqInfo%nFields) = iS

    Else If (ReqInfo%Time == S .and. .not.ReqInfo%NoReqs) Then

      ReqInfo%nFields                 = ReqInfo%nFields + 1
      If (ReqInfo%nFields > MaxFieldReqsPerReqInfo) Then
        Call Message('FATAL ERROR: MaxFieldReqsPerReqInfo is not large enough', 3)
      End If
      ReqInfo%iField(ReqInfo%nFields) = iField
      ReqInfo%iS    (ReqInfo%nFields) = iS

    End If

  End Do

End Subroutine CalcNextReqInfoS

!-------------------------------------------------------------------------------------------------------------

Subroutine ProcessAndOutputResults(                            &
             iCase, EndOfCase, EndOfCases,                     &
             Time, IncTimeField, IncPPInfoTime,                &
             MaxPPTimeLag,                                     &
             RadioactiveDecay, Chemistry, StartTime, EndTime,  &
             OutputOpts, SizeDists, OpenMPOpts,                &
             Coords, Grids, Domains, Specieses,                &
             CloudGammaParamses, Sources, Reqs, MaterialUnits, &
             Units, Mets, Flows, Results                       &
           )
! Processes and outputs results.

! This routine must be called at every sync time and at the end of each case, but can be called more often.

  Implicit None
  ! Argument list:
  Integer,                   Intent(In)    :: iCase
  Logical,                   Intent(In)    :: EndOfCase
  Logical,                   Intent(In)    :: EndOfCases
  Type(ShortTime_),          Intent(In)    :: Time
  Logical,                   Intent(In)    :: IncTimeField
  Logical,                   Intent(In)    :: IncPPInfoTime
  Type(ShortTime_),          Intent(In)    :: MaxPPTimeLag
  Logical,                   Intent(In)    :: RadioactiveDecay
  Logical,                   Intent(In)    :: Chemistry
  Type(ShortTime_),          Intent(In)    :: StartTime
  Type(ShortTime_),          Intent(In)    :: EndTime
  Type(OutputOpts_),         Intent(In)    :: OutputOpts
  Type(SizeDists_),          Intent(In)    :: SizeDists
  Type(OpenMPOpts_),         Intent(In)    :: OpenMPOpts
  Type(Coords_),             Intent(In)    :: Coords
  Type(Grids_),              Intent(In)    :: Grids
  Type(Domains_),            Intent(In)    :: Domains
  Type(Specieses_),          Intent(In)    :: Specieses
  Type(CloudGammaParamses_), Intent(In)    :: CloudGammaParamses
  Type(Sources_),            Intent(In)    :: Sources
  Type(Reqs_),               Intent(In)    :: Reqs
  Type(MaterialUnits_),      Intent(In)    :: MaterialUnits
  Type(Units_),              Intent(InOut) :: Units
  Type(Mets_),               Intent(InOut) :: Mets
  Type(Flows_),              Intent(InOut) :: Flows
  Type(Results_),            Intent(InOut) :: Results
  ! iCase              :: Number of case.
  ! EndOfCase          :: Indicates the end of the case. True on last call in a case, once all computation
  !                       completed for the case.
  ! EndOfCases         :: Indicates the end of the cases. True on last call in last case, once all computation
  !                       completed for the last case.
  ! Time               :: Time up to which output can be completed (subject to adjustments for MaxPPTimeLag).
  ! IncTimeField       :} Indicates that the period for which field or particle/puff information output can be
  ! IncPPInfoTime      :} completed includes Time. If false, these flags postpone processing and output of
  !                       output due at time Time (subject to adjustments for MaxPPTimeLag) and this option
  !                       must be used if the output might receive contributions at the start of the next sync
  !                       time (i.e. if any of the CalcReqInfo routines are subsequently called for this case
  !                       with Time1 set to Time and IncTime1 set to true). Has no effect if EndOfCase is
  !                       true.
  !                       $$ could check that, if true, subsequent calls to calcreqinfo are appropriate
  !                       - but need
  !                       to distinguish IncTime1Fields and IncTime1PPInfos in calls to calcreqinfo
  ! MaxPPTimeLag       :: The maximum possible interval which can occur at any time between the trailing
  !                       time-edge and the nominal time of the particles/puffs.
  ! RadioactiveDecay   :: Indicates that radioactive decay is modelled.
  ! Chemistry          :: Indicates that chemistry is modelled.
  ! StartTime          :} Start time and expected end time of the run.
  ! EndTime            :}
  ! OutputOpts         :: Output options.
  ! OpenMPOpts         :: OpenMP options.
  ! Coords             :: Collection of coord systems.
  ! Grids              :: Collection of grids.
  ! Domains            :: Collection of domains.
  ! Specieses          :: Collection of specieses.
  ! CloudGammaParamses :: Collection of sets of cloud gamma parameters.
  ! Sources            :: Collection of sources.
  ! Reqs               :: Collection of requirements.
  ! MaterialUnits      :: Collection of material units.
  ! Units              :: Collection of information on input/output unit numbers.
  ! Mets               :: Collection of met module instance states.
  ! Flows              :: Collection of flow module instance states.
  ! Results            :: Collection of results.
  ! Locals:
  Type(ShortTime_) :: TimeField ! Time up to which output can be completed for fields, adjusted for
                                ! MaxPPTimeLag.

  ! TimeField.
  If (EndOfCase) Then
    TimeField = Time
  Else
    TimeField = Time - MaxPPTimeLag
  End If

  ! Process fields.
  Call TimerOn(ProcessFieldsTimer)
  Call ProcessFields(                                      &
         iCase, EndOfCase, EndOfCases,                     &
         TimeField, IncTimeField .or. EndOfCase,           &
         RadioactiveDecay,                                 &
         OpenMPOpts,                                       &
         Coords, Grids, Domains, Specieses,                &
         CloudGammaParamses, Sources, Reqs, MaterialUnits, &
         Units, Mets, Flows, Results                       &
       )
  Call TimerOff(ProcessFieldsTimer)
  Call TimerWriteLast(ProcessFieldsTimer)

  ! Process pdfs.
  If (.not.IsInfPast(TimeField)) &
    Call DerivePdfs(ShortTime2Time(TimeField), Coords, Grids, Flows, Reqs, Results)

  ! Output fields.
  Call TimerOn(OutputFieldsTimer)
  Call OutputFields(                                      &
         iCase, EndOfCase, EndOfCases,                    &
         RadioactiveDecay, Chemistry, StartTime, EndTime, &
         OutputOpts, SizeDists, OpenMPOpts,               &
         Coords, Grids, Flows, Specieses, Sources, Reqs,  &
         MaterialUnits,                                   &
         Units, Results                                   &
       )
  Call TimerOff(OutputFieldsTimer)
  Call TimerWriteLast(OutputFieldsTimer)

  ! Output pdfs.
  If (.not.IsInfPast(TimeField)) Call OutputPdfs(                        &
                                        ShortTime2Time(TimeField),       &
                                        Coords, Grids, Reqs, OutputOpts, &
                                        Results, Units, iCase, Sources   &
                                      )

  ! Close and/or flush files used for output of particle/puff information.
  Call ClosePPInfos(                         &
         EndOfCase,                          &
         Time, IncPPInfoTime .or. EndOfCase, &
         Grids, Reqs,                        &
         Units, Results                      &
       )

End Subroutine ProcessAndOutputResults

!-------------------------------------------------------------------------------------------------------------

Subroutine ProcessFields(                                      &
             iCase, EndOfCase, EndOfCases,                     &
             Time, IncTime,                                    &
             RadioactiveDecay,                                 &
             OpenMPOpts,                                       &
             Coords, Grids, Domains, Specieses,                &
             CloudGammaParamses, Sources, Reqs, MaterialUnits, &
             Units, Mets, Flows, Results                       &
           )
! Processes fields.

  Implicit None
  ! Argument list:
  Integer,                   Intent(In)            :: iCase
  Logical,                   Intent(In)            :: EndOfCase
  Logical,                   Intent(In)            :: EndOfCases
  Type(ShortTime_),          Intent(In)            :: Time
  Logical,                   Intent(In)            :: IncTime
  Logical,                   Intent(In)            :: RadioactiveDecay
  Type(OpenMPOpts_),         Intent(In)            :: OpenMPOpts
  Type(Coords_),             Intent(In)            :: Coords
  Type(Grids_),              Intent(In),    Target :: Grids
  Type(Domains_),            Intent(In)            :: Domains
  Type(Specieses_),          Intent(In)            :: Specieses
  Type(CloudGammaParamses_), Intent(In)            :: CloudGammaParamses
  Type(Sources_),            Intent(In)            :: Sources
  Type(Reqs_),               Intent(In),    Target :: Reqs
  Type(MaterialUnits_),      Intent(In)            :: MaterialUnits
  Type(Units_),              Intent(InOut)         :: Units
  Type(Mets_),               Intent(InOut)         :: Mets
  Type(Flows_),              Intent(InOut)         :: Flows
  Type(Results_),            Intent(InOut), Target :: Results
  ! iCase              :: Number of case.
  ! EndOfCase          :: Indicates the end of the case.
  ! EndOfCases         :: Indicates the end of the cases.
  ! Time               :: Time up to which fields can be processed.
  ! IncTime            :: Indicates period for which fields can be processed includes Time.
  ! RadioactiveDecay   :: Indicates that radioactive decay is modelled.
  ! OpenMPOpts         :: OpenMP options.
  ! Coords             :: Collection of coord systems.
  ! Grids              :: Collection of grids.
  ! Domains            :: Collection of domains.
  ! Specieses          :: Collection of specieses.
  ! CloudGammaParamses :: Collection of sets of cloud gamma parameters.
  ! Sources            :: Collection of sources.
  ! Reqs               :: Collection of requirements.
  ! MaterialUnits      :: Collection of material units.
  ! Units              :: Collection of information on input/output unit numbers.
  ! Mets               :: Collection of met module instance states.
  ! Flows              :: Collection of flow module instance states.
  ! Results            :: Collection of results.
  ! Locals:
  Type(TGrid_),    Pointer :: TGrid                          !} Abbreviations for grids, processes, field
  Type(TGrid_),    Pointer :: TGridC                         !} requirements and fields. 'C' denotes items
  Type(DGrid_),    Pointer :: DGrid                          !} related to a field contributing to a
  Type(DGrid_),    Pointer :: DGridC                         !} processing step.
  Type(Process_),  Pointer :: Process                        !}
  Type(Process_),  Pointer :: ProcessC                       !}
  Type(FieldReq_), Pointer :: FieldReq                       !}
  Type(FieldReq_), Pointer :: FieldReqC                      !}
  Type(Field_),    Pointer :: Field                          !}
  Type(Field_),    Pointer :: FieldC                         !}
  Type(Field_),    Pointer :: FieldAirConc                   !}
  Type(Field_),    Pointer :: FieldEulConc                   !}
  Type(Field_),    Pointer :: FieldChemField                 !}
  Type(Field_),    Pointer :: FieldDensity                   !}
  Type(Field_),    Pointer :: FieldMeanS                     !}
  Type(Field_),    Pointer :: FieldXStats                    !}
  Type(Field_),    Pointer :: FieldMeanZ                     !}
  Type(Field_),    Pointer :: FieldMass                      !}
  Type(Field_),    Pointer :: FieldDryDep                    !}
  Type(Field_),    Pointer :: FieldWetDep                    !}
  Type(Field_),    Pointer :: FieldPhotonFlux                !}
  Integer                  :: iTGrid                         !] Indices of grids and fields. 'C' denotes items
  Integer                  :: iTGridC                        !] related to a field contributing to a
  Integer                  :: iField                         !] processing step.
  Integer                  :: iFieldC                        !]
  Integer                  :: iAirConc                       !]
  Integer                  :: iEulConc                       !]
  Integer                  :: iChemField                     !]
  Integer                  :: iMeanS                         !]
  Integer                  :: iXStats                        !]
  Integer                  :: iMeanZ                         !]
  Integer                  :: iMass                          !]
  Integer                  :: iDryDep                        !]
  Integer                  :: iWetDep                        !]
  Integer                  :: iDensity                       !]
  Integer                  :: iPhotonFlux                    !]
  Integer                  :: nT                             ! Size of time grid (1 if no time grid).
  Type(ShortTime_)         :: T                              !} Times, indices of times and corresponding
  Type(ShortTime_)         :: TC                             !} indices in field arrays. 'C' denotes items
  Integer                  :: iT                             !} related to a field contributing to a
  Integer                  :: iTL                            !}
  Integer                  :: iTU                            !} 
  Integer                  :: iTC                            !} processing step.
  Integer                  :: iTA                            !}
  Integer                  :: iTAC                           !}
  Integer                  :: LastProcThisCase(Reqs%MaxFieldReqs) !] Values of Results%Fields(:)%LastProcThisCase
  Integer                  :: LastProc(Reqs%MaxFieldReqs)         !] and Results%Fields(:)%LastProc on entry to
                                                                  !] this routine.
  Logical                  :: FieldProcessed(Reqs%MaxFieldReqs)   ! Indicates which fields have been processed.
  Type(ShortTime_)         :: DummyT                              ! Dummy variable needed for an argument list.
  Type(MaterialUnit_)      :: SpeciesMaterialUnit                 ! Material unit associated with species
  Type(MaterialUnit_)      :: FieldReqMaterialUnit                ! Material unit associated with field request
  Real                     :: ThresholdConversionFactor           !} Rescaling factor for probability thresholds, 
                                                                  !} if material unit of species differs from
                                                                  !} material unit in field request
  Real                     :: UnitConversionFactor                !} Conversion factor between units of 
                                                                  !} field request and associated species
  Real                     :: Mr_Ozone = 47.998                   ! molar mass of ozone [g mol^{-1}]

  Results%MaxTimeField = TMax(Results%MaxTimeField, Time)

  FieldProcessed(1:Reqs%nFieldReqs) = .false.

  Do iField = 1, Reqs%nFieldReqs
    LastProcThisCase(iField) = Results%Fields(iField)%LastProcThisCase
    LastProc(iField)         = Results%Fields(iField)%LastProc
  End Do

  ! Loop to repeat processing until all fields processed.
  Do While (.not. All(FieldProcessed(1:Reqs%nFieldReqs)))

    ! Loop over fields. This loop runs backwards for improved efficiency - in general fields needed to process
    ! a given field will be later in the collection of fields.
    !
    ! For each field we determine whether the other fields needed for processing are available and, if they
    ! are, we process the field. When a field is processed, the processing does as much as possible consistent
    ! with Time, IncTime, Results%MaxTimeField, EndOfCase and EndOfCases. Hence, in any one call to this
    ! routine, there is no need to process a field more than once.
    !
    ! Note all fields are processed soon as possible and, if used by other fields, are used as soon as
    ! possible after processing (and before leaving this routine). This means that Field%LastProcThisCase and
    ! Field%LastProc at the end of this routine are calculable from Time, IncTime, Results%MaxTimeField,
    ! EndOfCase, EndOfCases and FieldReq (including its T-grid).
    FieldsLoop: Do iField = Reqs%nFieldReqs, 1, -1

      ! Cycle if field already processed.
      If (FieldProcessed(iField)) Cycle

      ! Abbreviations for field requirement and field.
      FieldReq => Reqs%FieldReqs(iField)
      Field    => Results%Fields(iField)

      ! Set up time grid and its size.
      iTGrid = FieldReq%iTGrid
      If (iTGrid == 0) Then
        nT = 1
      Else
        TGrid => Grids%TGrids(iTGrid)
        nT    =  TGrid%nT
      End If

      ! Fields with an explicit processing step. Note "FieldReq%iProc /= 0" and "FieldReq%iLastProcess > 0"
      ! are equivalent and both conditions mean there is an explicit processing step.
      If (FieldReq%iExt == 0 .and. FieldReq%iLastProcess > 0) Then

        ! Locate other field required.
        iFieldC   =  FieldReq%iProc
        FieldReqC => Reqs%FieldReqs(iFieldC)
        FieldC    => Results%Fields(iFieldC)

        ! Cycle if required fields not processed or if required fields need processing over the ensemble of
        ! cases and EndOfCases not set.
        If (.not.FieldProcessed(iFieldC)) Cycle FieldsLoop
        If (.not.EndOfCases .and. FieldReqC%Ensemble) Then
          FieldProcessed(iField) = .true.
          Cycle FieldsLoop
        End If

        ! Set up Process, DGrid, ProcessC and DGridC.
        Process => Reqs%Processes(FieldReq%iProcesses(FieldReq%iLastProcess))
        If (Process%iDGrid /= 0) DGrid => Grids%DGrids(Process%iDGrid)
        If (FieldReqC%iLastProcess > 0) Then
          ProcessC => Reqs%Processes(FieldReqC%iProcesses(FieldReqC%iLastProcess))
          If (ProcessC%iDGrid /= 0) DGridC => Grids%DGrids(ProcessC%iDGrid)
        End If

        ! Set up iTGridC and TGridC.
        iTGridC = FieldReqC%iTGrid
        If (iTGridC /= 0) TGridC => Grids%TGrids(FieldReqC%iTGrid)

        ! Loop over time indices of field.
        Do iT = Field%LastProcThisCase + 1, nT
          ! $$ Amend if allow infinite grids - use Field%iTL + Field%nT - 1?

          ! Calculate T. Note T used in Field2AvInt to calculate decay of deposited material.
          If (Scan(QInfo(FieldReq%iQuantity), 'T') == 0) Then
            T = InfPastShortTime()
          Else If (iTGrid == 0) Then
            T = InfFutureShortTime()
          Else
            T = TInTGrid(TGrid, iT)
          End If

          ! Loop over time indices of contributing field.
          Do iTC = LastProc(iFieldC) + 1, FieldC%LastProc

            ! Calculate TC and do tests to cycle and exit loop. Note TC used in Field2AvInt to calculate decay
            ! of deposited material.

            ! T-gridded field.
            If (iTGrid /= 0) Then

              TC = TInTGrid(TGridC, iTC)
              If (TC > T) Exit
              If (TC < T - Process%sdT * (Process%nT - 1)) Cycle
              ! $$ for irreg grids, also test if TC /= round(TC, T, Process%dt)
              ! $$ other processing types

            ! Non-T-gridded field with T-gridded contributing field.
            Else If (iTGridC /= 0) Then

              TC = TInTGrid(TGridC, iTC)
              ! $$ check contrib time correct

            ! Non-T-gridded contributing field due to underlying quantity not depending on time.
            Else If (Scan(QInfo(FieldReq%iQuantity), 'T') == 0) Then

              TC = InfPastShortTime()

            ! Non-T-gridded contributing field due to field processed over all time.
            Else

              TC = InfFutureShortTime()

            End If

            ! Indices in field arrays corresponding to time indices.
            iTA  = CalciTA(Field,  iT )
            iTAC = CalciTA(FieldC, iTC)

            ! Process field.
            If (Process%Code == P_AvInt) Then
              Call Field2AvInt(                                         &
                     Specieses, Sources, RadioactiveDecay,              &
                     OpenMPOpts,                                        &
                     Process, FieldReq, TC, T, iTAC, iTA, FieldC, Field &
                   )
            Else If (Process%Code == P_Prob) Then
              ! If the Material unit of the Field Request differs from that of the 
              ! associated species, work out how the probability thresholds need to be adjusted
              If ( (FieldReq%iQuantity == Q_AirConc              ) .or. &
                   (FieldReq%iQuantity == Q_ChemistryField       ) .or. &
                   (FieldReq%iQuantity == Q_EulerianConcentration) .or. &
                   (FieldReq%iQuantity == Q_Concentration        ) .or. &
                   (FieldReq%iQuantity == Q_DryDep               ) .or. &
                   (FieldReq%iQuantity == Q_WetDep               ) .or. &
                   (FieldReq%iQuantity == Q_Dep                  ) .or. &
                   (FieldReq%iQuantity == Q_EulerianTotalDep     ) .or. &
                   (FieldReq%iQuantity == Q_EulerianDryDep       ) .or. &
                   (FieldReq%iQuantity == Q_EulerianWetDep       )      &
              ) Then
                SpeciesMaterialUnit = MaterialUnits%MaterialUnits(                           &
                                        Specieses%Specieses(FieldReq%iSpecies)%iMaterialUnit &
                                      )
                FieldReqMaterialUnit = MaterialUnits%MaterialUnits( &
                                          FieldReq%iMaterialUnit    &
                                       )
                If (FieldReqMaterialUnit .ne. SpeciesMaterialUnit) Then 
                  ThresholdConversionFactor = SpeciesMaterialUnit%ConversionFactor  &
                                            / FieldReqMaterialUnit%ConversionFactor
                  ! For Dobson units, convert to correct molar mass
                  If (FieldReqMaterialUnit%UnitType .eq. DobsonUnitType) Then
                    ThresholdConversionFactor = ThresholdConversionFactor             &
                                         * Mr_Ozone / Specieses%Specieses(FieldReq%iSpecies)%MolecularWeight
                  End If
                  Call Field2Prob(Process, iTAC, iTA, DGrid, FieldC, Field, ThresholdConversionFactor)
                Else
                  Call Field2Prob(Process, iTAC, iTA, DGrid, FieldC, Field)
                End If
              Else
                Call Field2Prob(Process, iTAC, iTA, DGrid, FieldC, Field)
              End If
            Else If (Process%Code == P_Percent) Then
              Call Prob2Percentile(Process, iTAC, iTA, DGridC, DGrid, FieldC, Field)
            Else If (Process%Code == P_MaxMin) Then
              Call Field2MaxMin(Process, iTAC, iTA, DGrid, FieldC, Field)
            End If

          End Do

        End Do

        ! Update Field%LastProcThisCase.
        If (Scan(QInfo(FieldReq%iQuantity), 'T') == 0) Then
          Field%LastProcThisCase = 1
        Else If (iTGrid == 0) Then
          If (EndOfCase) Then
            Field%LastProcThisCase = 1
          Else
            Field%LastProcThisCase = 0
          End If
        Else
          If (FieldReq%Ensemble .and. EndOfCases) Then
            Call LastTBeforeT(TGrid, Results%MaxTimeField, .false.,      DummyT, Field%LastProcThisCase)
          Else
            Call LastTBeforeT(TGrid, Time,                 .not.IncTime, DummyT, Field%LastProcThisCase)
          End If
          Field%LastProcThisCase = Max(Field%LastProcThisCase, LastProcThisCase(iField))
          Field%LastProcThisCase = Min(Field%LastProcThisCase, nT) ! $$ Amend if allow infinite grids
        End If

        ! Finalise field and update Field%LastProc.
        If ((FieldReq%Ensemble .and. EndOfCases) .or. .not. FieldReq%Ensemble) Then
          Do iT = Field%LastProc + 1, Field%LastProcThisCase
            If (Process%Code == P_AvInt) Then
              iTA = CalciTA(Field, iT)
              If (Process%TAvInt == P_Av) Then
                If (Field%TypeType == ' ') Then
                  Field%Std(:, :, :, :, 1, iTA, :) = Field%Std(:, :, :, :, 1, iTA, :) / Process%nT
                Else If (Field%TypeType == '8') Then
                  Field%P64(:, :, :, :, 1, iTA, :) = Field%P64(:, :, :, :, 1, iTA, :) / Process%nT
                End If
              Else If (Process%TAvInt == P_Int) Then
                If (Field%TypeType == ' ') Then
                  Field%Std(:, :, :, :, 1, iTA, :) = Field%Std(:, :, :, :, 1, iTA, :) * &
                                                     ShortTime2RealTime(Process%sdT)
                Else If (Field%TypeType == '8') Then
                  Field%P64(:, :, :, :, 1, iTA, :) = Field%P64(:, :, :, :, 1, iTA, :) * &
                                                     ShortTime2RealTime(Process%sdT)
                End If
              End If
              If (Process%Ensemble) Then
                If (Field%TypeType == ' ') Then
                  Field%Std(:, :, :, :, 1, iTA, :) = Field%Std(:, :, :, :, 1, iTA, :) / Real(iCase, Std)
                Else If (Field%TypeType == '8') Then
                  Field%P64(:, :, :, :, 1, iTA, :) = Field%P64(:, :, :, :, 1, iTA, :) / Real(iCase, Std)
                End If
              End If
            Else If (Process%Code == P_Prob) Then
              iTA = CalciTA(Field, iT)
              Field%Std(:, :, :, :, 1, iTA, :) = Field%Std(:, :, :, :, 1, iTA, :) / Process%nT
              If (Process%Ensemble) Then
                Field%Std(:, :, :, :, 1, iTA, :) = Field%Std(:, :, :, :, 1, iTA, :) / Real(iCase, Std)
              End If
            End If
          End Do
          Field%LastProc = Field%LastProcThisCase
        End If

        FieldProcessed(iField) = .true.

      ! $$ For fields with overlap of intrinsic averaging regions the above should work fine provided
      ! $$ the right processes are used.
      ! $$ Note decay of deposition needs to be taken account of in Case.F90 in generating intrinsic
      ! averaged fields.

      ! Fields of type P, Q or R.
      ! $$ Currently includes Deposition rate, Sigma C, X Stats, Mean Travel Time, Mean Z, Sigma Z and
      ! $$ adult effective and organ cloud gamma dose.
      !    Make these all P/Q/R, make if test on pqr attribute
      Else If (                                                   &
        FieldReq%iExt == 0                                  .and. &
        (                                                         &
          FieldReq%iQuantity == Q_Dep                  .or.       &
          FieldReq%iQuantity == Q_Concentration        .or.       &
          FieldReq%iQuantity == Q_MixingRatio          .or.       &
          FieldReq%iQuantity == Q_EMixingRatio         .or.       &
          FieldReq%iQuantity == Q_SigmaC               .or.       &
          FieldReq%iQuantity == Q_XStats               .or.       &
          FieldReq%iQuantity == Q_MeanS                .or.       &
          FieldReq%iQuantity == Q_MeanZ                .or.       &
          FieldReq%iQuantity == Q_SigmaZ               .or.       &
          FieldReq%iQuantity == Q_AreaAtRisk           .or.       &
          FieldReq%iQuantity == Q_AduEffCloudGammaDose .or.       &
          FieldReq%iQuantity == Q_AduLunCloudGammaDose .or.       &
          FieldReq%iQuantity == Q_AduThyCloudGammaDose .or.       &
          FieldReq%iQuantity == Q_AduBoSCloudGammaDose            &
        )                                                         &
      ) Then

        ! Locate other fields required.
        iDryDep      = FieldReq%iDryDep
        iWetDep      = FieldReq%iWetDep
        iAirConc     = FieldReq%iAirConc
        iEulConc     = FieldReq%iEulConc
        iChemField   = FieldReq%iChemField
        iDensity     = FieldReq%iDensity  
        iMeanS       = FieldReq%iMeanS
        iXStats      = FieldReq%iXStats
        iMass        = FieldReq%iMass
        iMeanZ       = FieldReq%iMeanZ
        iPhotonFlux  = FieldReq%iPhotonFlux

        ! Cycle if required fields not processed or if required fields need processing over the ensemble of
        ! cases and EndOfCases not set.
        If (iDryDep /= 0) Then
          If (.not.FieldProcessed(iDryDep)) Cycle FieldsLoop
          If (.not.EndOfCases .and. Reqs%FieldReqs(iDryDep)%Ensemble) Then
            FieldProcessed(iField) = .true.
            Cycle FieldsLoop
          End If
          FieldDryDep => Results%Fields(iDryDep)
        End If

        If (iWetDep /= 0) Then
          If (.not.FieldProcessed(iWetDep)) Cycle FieldsLoop
          If (.not.EndOfCases .and. Reqs%FieldReqs(iWetDep)%Ensemble) Then
            FieldProcessed(iField) = .true.
            Cycle FieldsLoop
          End If
          FieldWetDep => Results%Fields(iWetDep)
        End If

        If (iAirConc /= 0) Then
          If (.not.FieldProcessed(iAirConc)) Cycle FieldsLoop
          If (.not.EndOfCases .and. Reqs%FieldReqs(iAirConc)%Ensemble) Then
            FieldProcessed(iField) = .true.
            Cycle FieldsLoop
          End If
          FieldAirConc => Results%Fields(iAirConc)
        End If

        If (iEulConc /= 0) Then
          If (.not.FieldProcessed(iEulConc)) Cycle FieldsLoop
          If (.not.EndOfCases .and. Reqs%FieldReqs(iEulConc)%Ensemble) Then
            FieldProcessed(iField) = .true.
            Cycle FieldsLoop
          End If
          FieldEulConc => Results%Fields(iEulConc)
        End If

        If (iChemField /= 0) Then
          If (.not.FieldProcessed(iChemField)) Cycle FieldsLoop
          If (.not.EndOfCases .and. Reqs%FieldReqs(iChemField)%Ensemble) Then
            FieldProcessed(iField) = .true.
            Cycle FieldsLoop
          End If
          FieldChemField => Results%Fields(iChemField)
        End If

        If (iDensity /= 0) Then
          If (.not.FieldProcessed(iDensity)) Cycle FieldsLoop
          If (.not.EndOfCases .and. Reqs%FieldReqs(iDensity)%Ensemble) Then
            FieldProcessed(iField) = .true.
            Cycle FieldsLoop
          End If
          FieldDensity => Results%Fields(iDensity)
        End If

        If (iMeanS /= 0) Then
          If (.not.FieldProcessed(iMeanS)) Cycle FieldsLoop
          If (.not.EndOfCases .and. Reqs%FieldReqs(iMeanS)%Ensemble) Then
            FieldProcessed(iField) = .true.
            Cycle FieldsLoop
          End If
          FieldMeanS => Results%Fields(iMeanS)
        End If

        If (iXStats /= 0) Then
          If (.not.FieldProcessed(iXStats)) Cycle FieldsLoop
          If (.not.EndOfCases .and. Reqs%FieldReqs(iXStats)%Ensemble) Then
            FieldProcessed(iField) = .true.
            Cycle FieldsLoop
          End If
          FieldXStats => Results%Fields(iXStats)
        End If

        If (iMass /= 0) Then
          If (.not.FieldProcessed(iMass)) Cycle FieldsLoop
          If (.not.EndOfCases .and. Reqs%FieldReqs(iMass)%Ensemble) Then
            FieldProcessed(iField) = .true.
            Cycle FieldsLoop
          End If
          FieldMass => Results%Fields(iMass)
        End If

        If (iMeanZ /= 0) Then
          If (.not.FieldProcessed(iMeanZ)) Cycle FieldsLoop
          If (.not.EndOfCases .and. Reqs%FieldReqs(iMeanZ)%Ensemble) Then
            FieldProcessed(iField) = .true.
            Cycle FieldsLoop
          End If
          FieldMeanZ => Results%Fields(iMeanZ)
        End If

        If (iPhotonFlux /= 0) Then
          If (.not.FieldProcessed(iPhotonFlux)) Cycle FieldsLoop
          If (.not.EndOfCases .and. Reqs%FieldReqs(iPhotonFlux)%Ensemble) Then
            FieldProcessed(iField) = .true.
            Cycle FieldsLoop
          End If
          FieldPhotonFlux => Results%Fields(iPhotonFlux)
        End If

        ! Update Field%LastProcThisCase. Note that, for convenience, it is updated before the processing is
        ! done.
        If (Scan(QInfo(FieldReq%iQuantity), 'T') == 0) Then
          Field%LastProcThisCase = 1
        Else If (iTGrid == 0) Then
          If (EndOfCase) Then
            Field%LastProcThisCase = 1
          Else
            Field%LastProcThisCase = 0
          End If
        Else
          If (FieldReq%Ensemble .and. EndOfCases) Then
            Call LastTBeforeT(TGrid, Results%MaxTimeField, .false.,      DummyT, Field%LastProcThisCase)
          Else
            Call LastTBeforeT(TGrid, Time,                 .not.IncTime, DummyT, Field%LastProcThisCase)
          End If
          Field%LastProcThisCase = Max(Field%LastProcThisCase, LastProcThisCase(iField))
          Field%LastProcThisCase = Min(Field%LastProcThisCase, nT) ! $$ Amend if allow infinite grids
        End If

        ! Loop over time indices of field.
        Do iT = LastProcThisCase(iField) + 1, Field%LastProcThisCase

          ! Index in field array corresponding to time index.
          iTA = CalciTA(Field, iT)
          ! $$ used fields might not have same t-grid once field extensions implemented

          ! Process field. $$ put processing in separate routines
          If (FieldReq%iQuantity == Q_Dep) Then

            Field%Std(1, :, :, :, 1, iTA, :) = FieldDryDep%Std(1, :, :, :, 1, iTA, :) + &
                                               FieldWetDep%Std(1, :, :, :, 1, iTA, :)

          Else If (FieldReq%iQuantity == Q_Concentration) Then

            Field%Std(1, :, :, :, 1, iTA, :) = FieldAirConc%Std(1, :, :, :, 1, iTA, :) + &
                                               FieldEulConc%Std(1, :, :, :, 1, iTA, :)

          Else If (FieldReq%iQuantity == Q_MixingRatio) Then

            Call DeriveMixingRatio(                            &
                   iField, iAirConc, iDensity,                 &
                   FieldReq%iSpecies, FieldReq%iMaterialUnit,  &
                   iTA, Specieses, Results, MaterialUnits      &
                 )

           Else If (FieldReq%iQuantity == Q_EMixingRatio) Then

            Call DeriveEMixingRatio(                         &
                   iField, iChemField, iDensity,             &
                   FieldReq%iSpecies,FieldReq%iMaterialUnit, &
                   iTA, Specieses, Results, MaterialUnits    &
                 )
           
          Else If (FieldReq%iQuantity == Q_SigmaC) Then

            Call DeriveSigmaCField(                            &
                   iField, iAirConc, iMeanS, iXStats, iT, iTA, &
                   Coords, Grids, Domains, Sources, Reqs,      &
                   Units, Mets, Flows, Results                 &
                 )


          Else If (                                           &
            FieldReq%iQuantity == Q_AduEffCloudGammaDose .or. &
            FieldReq%iQuantity == Q_AduLunCloudGammaDose .or. &
            FieldReq%iQuantity == Q_AduThyCloudGammaDose .or. &
            FieldReq%iQuantity == Q_AduBoSCloudGammaDose      &
          ) Then

            If (FieldReq%SemiInfApprox) Then
              Call Derive_CloudGammaDose_via_SemiInfApprox(              &
                     iField, iAirConc, iT, iTA,                          &
                     Grids, Reqs, Results, Specieses, CloudGammaParamses &
                   )
            Else
              Call DeriveCloudGammaDose(                                 &
                     iField, iPhotonFlux, iT, iTA,                       &
                     Grids, Reqs, Results, Specieses, CloudGammaParamses &
                   )
            End If

          Else If (FieldReq%iQuantity == Q_AreaAtRisk) Then

            iTAC = CalciTA(FieldAirConc, iT)
            ! $$ Should do this for other derived quantities, but only relevant if avaerging before deriving.
            ! $$ Also used fields might not have same t-grid once field extensions implemented
            Call DeriveAreaAtRisk(                                       &
                   iField, iAirConc,                                     &
                   iTA, iTAC,                                            &
                   Coords, Grids, Domains, Flows, Sources, Reqs, Results &
                 )

          Else If (FieldReq%iQuantity == Q_XStats) Then

            Call DeriveXStats(Grids, FieldReq%iSGrid, iTA, Field)

          Else If (FieldReq%iQuantity == Q_MeanS) Then

            Call DeriveMeanTravelTimeField(iField, iAirConc, iTA, Coords, Grids, Reqs, Results)

          Else If (FieldReq%iQuantity == Q_MeanZ) Then

            If (FieldMass%Std(1, 1, 1, 1, 1, iTA, 1) == 0.0) Then
              Field%Std(1, 1, 1, 1, 1, iTA, 1) = 0.0
            Else
              Field%Std(1, 1, 1, 1, 1, iTA, 1) = Field%Std(1, 1, 1, 1, 1, iTA, 1) /   &
                                                 FieldMass%Std(1, 1, 1, 1, 1, iTA, 1)
            End If

          Else If (FieldReq%iQuantity == Q_SigmaZ) Then

            If (FieldMass%Std(1, 1, 1, 1, 1, iTA, 1) == 0.0) Then
              Field%Std(1, 1, 1, 1, 1, iTA, 1) = 0.0
            Else
              Field%Std(1, 1, 1, 1, 1, iTA, 1) = Field%Std(1, 1, 1, 1, 1, iTA, 1) /   &
                                                 FieldMass%Std(1, 1, 1, 1, 1, iTA, 1)
              Field%Std(1, 1, 1, 1, 1, iTA, 1) = Sqrt(                                         &
                                                   Max(                                        &
                                                     Field%Std(1, 1, 1, 1, 1, iTA, 1) -        &
                                                     FieldMeanZ%Std(1, 1, 1, 1, 1, iTA, 1)**2, &
                                                     0.0                                       &
                                                   )                                           &
                                                 )
            End If

          Else

            Call Message('UNEXPECTED ERROR in ProcessFields', 4)

          End If

        End Do

        If ((FieldReq%Ensemble .and. EndOfCases) .or. .not. FieldReq%Ensemble) Then
          Field%LastProc = Field%LastProcThisCase
        End If

        FieldProcessed(iField) = .true.

      ! Other fields (null processing).
      Else

        ! Update Field%LastProcThisCase.
        If (Scan(QInfo(FieldReq%iQuantity), 'T') == 0) Then
          Field%LastProcThisCase = 1
        Else If (iTGrid == 0) Then
          If (EndOfCase) Then
            Field%LastProcThisCase = 1
          Else
            Field%LastProcThisCase = 0
          End If
        Else
          If (FieldReq%Ensemble .and. EndOfCases) Then
            Call LastTBeforeT(TGrid, Results%MaxTimeField, .false.,      DummyT, Field%LastProcThisCase)
          Else
            Call LastTBeforeT(TGrid, Time,                 .not.IncTime, DummyT, Field%LastProcThisCase)
          End If
          Field%LastProcThisCase = Max(Field%LastProcThisCase, LastProcThisCase(iField))
          Field%LastProcThisCase = Min(Field%LastProcThisCase, nT) ! $$ Amend if allow infinite grids
        End If

        If ((FieldReq%Ensemble .and. EndOfCases) .or. .not. FieldReq%Ensemble) Then
          Field%LastProc = Field%LastProcThisCase
        End If

        FieldProcessed(iField) = .true.

      End If

    End Do FieldsLoop

  End Do

! Convert units
  
  ConvertUnitsLoop : Do iField = 1, Reqs%nFieldReqs
     Field    => Results%Fields(iField)
     FieldReq => Reqs%FieldReqs(iField)

! Only convert fields that are output to disk or screen, only convert ensemble requests at the end of the run
     If (FieldReq%Ensemble .and. .not.EndOfCases) Cycle
     If ((.not. FieldReq%Disk) .and. (.not. FieldReq%Screen)) Cycle
! Do not scale probabilities
     If (FieldReq%PCode == P_Prob) Cycle
! Currently only convert the following quantities:
     If ( (.not. FieldReq%iQuantity == Q_AirConc               ) .and. &
          (.not. FieldReq%iQuantity == Q_ChemistryField        ) .and. &
          (.not. FieldReq%iQuantity == Q_EulerianConcentration ) .and. &
          (.not. FieldReq%iQuantity == Q_Concentration         ) .and. &
          (.not. FieldReq%iQuantity == Q_DryDep                ) .and. &
          (.not. FieldReq%iQuantity == Q_WetDep                ) .and. &
          (.not. FieldReq%iQuantity == Q_Dep                   ) .and. &
          (.not. FieldReq%iQuantity == Q_OriginalSourceStrength) .and. &
          (.not. FieldReq%iQuantity == Q_RevisedSourceStrength)        &
        ) Cycle
     SpeciesMaterialUnit = MaterialUnits%MaterialUnits(                           &
                             Specieses%Specieses(FieldReq%iSpecies)%iMaterialUnit &
                           )
     FieldReqMaterialUnit = MaterialUnits%MaterialUnits(FieldReq%iMaterialUnit)
     
! No need to convert if same units are used in Species and Field request
     If (FieldReqMaterialUnit .eq. SpeciesMaterialUnit) Cycle
     UnitConversionFactor = SpeciesMaterialUnit%ConversionFactor / FieldReqMaterialUnit%ConversionFactor
! For output in Dobson units, correct for different molar mass of species
     If (FieldReqMaterialUnit%UnitType .eq. DobsonUnitType) Then
       UnitConversionFactor = UnitConversionfactor                                              &
                            * Mr_Ozone / Specieses%Specieses(FieldReq%iSpecies)%MolecularWeight
     End If
     If (Reqs%FieldReqs(iField)%Ensemble) Then
       iTL = LastProc(iField) + 1
       iTU = Field%LastProc
     Else
       iTL = LastProcThisCase(iField) + 1
       iTU = Field%LastProcThisCase
     End If

     Do iT = iTL, iTU
       ! Index in field array corresponding to time index.
       iTA = CalciTA(Field, iT)
       ! $$ used fields might not have same t-grid once field extensions implemented

       If (Field%TypeType .eq. ' ') Then
         Field%Std(:, :, :, :, :, iTA, :) = Field%Std(:, :, :, :, :, iTA, :)*UnitConversionFactor
       Else If (Field%TypeType .eq. '8') Then
         Field%P64(:, :, :, :, :, iTA, :) = Field%P64(:, :, :, :, :, iTA, :)*UnitConversionFactor
       Else
         Call Message('ERROR: Cannot convert units for field of type ' // Trim(Field%TypeType), &
                      3                                                                         &
              )
       End If
     End Do
    
  End Do ConvertUnitsLoop

End Subroutine ProcessFields

!-------------------------------------------------------------------------------------------------------------

Subroutine Field2AvInt(                                                          &
             Specieses, Sources, RadioactiveDecay,                               &
             OpenMPOpts,                                                         &
             Process, AvIntReq, TField, TAvInt, iTAField, iTAAvInt, Field, AvInt &
           )
! Calculate contribution to an averaged/integrated field from one time level of a field.

! Note the contributions are simply added to aid precision. Appropriate scaling must be done later.

! Note any pre-existing D-grid (which cannot be a floating D-grid) is carried forward.

  Implicit None
  ! Argument list:
  Type(Specieses_), Intent(In)    :: Specieses        ! Collection of specieses.
  Type(Sources_),   Intent(In)    :: Sources          ! Collection of sources.
  Logical,          Intent(In)    :: RadioactiveDecay ! Indicates that radioactive decay is modelled.
  Type(OpenMPOpts_),Intent(In)    :: OpenMPOpts       ! OpenMP options
  Type(Process_),   Intent(In)    :: Process          ! Process to be used in generating averaged/integrated
                                                      ! field.
  Type(FieldReq_),  Intent(In)    :: AvIntReq         ! Field requirement for averaged/integrated field.
  Type(ShortTime_), Intent(In)    :: TField           ! Time of interest in Field.
  Type(ShortTime_), Intent(In)    :: TAvInt           ! Time of interest in AvInt.
  Integer,          Intent(In)    :: iTAField         ! Index in Field corresponding to time of interest.
  Integer,          Intent(In)    :: iTAAvInt         ! Index in AvInt corresponding to time of interest.
  Type(Field_),     Intent(In)    :: Field            ! Contributing field.
  Type(Field_),     Intent(InOut) :: AvInt            ! Averaged/integrated field.
  ! Locals:
  Real(Std) :: DecayFactor             ! Decay factor for deposition.
  Integer   :: iD, iX, iY, iZ, iS, iC  ! Loop variables
  Integer   :: lBounds(7), uBounds(7)  ! Lower and upper bounds of field arrays

  ! Calculate decay factor for deposition.
  ! This is applied only on the first processing step.
  ! TAvInt and/or TField could be in the infinite future (with TAvInt >= TField).
  ! TField can be infinite even if AvIntReq%nProcesses = 1 (conside a single processing step with intrinsic
  ! averaging over all time and some spatial averaging).
  ! TAvInt == TField is trapped here which will ensure TField is finite and prevent the time subtractions
  ! failing.
  If (                                             &
    RadioactiveDecay .and. AvIntReq%DecayDep .and. &
    AvIntReq%nProcesses == 1                 .and. &
    TAvInt /= TField                               &
  ) Then
    Call CalcDecayFactor(                                                                        &
           Specieses%Specieses(AvIntReq%iSpecies),                                               &
           TMax(TField - Sources%Sources(Max(AvIntReq%iSource, 1))%sStartTime, ZeroShortTime()), &
           TAvInt - TField,                                                                      &
           DecayFactor                                                                           &
         )
  Else
    DecayFactor = 1.0
  End If

  ! $$ need to treat same physics for implicitly averaged dep

  ! $$ if do probs/percentiles these will be in std even if quantity is p64. (Need to allocate as std.)
  ! For av/int less clear. Could use p64 (as here) or revert to std.

  If (Field%TypeType == ' ') Then
    lBounds = lBound(AvInt%Std)
    uBounds = uBound(AvInt%Std)
    Do iC = lBounds(7), uBounds(7)
      Do iS = lBounds(5), uBounds(5)
        Do iZ = lBounds(4), uBounds(4)
          !$OMP PARALLEL DO                                                                           &
#ifdef    UseOpenMP3_0
          !$OMP COLLAPSE(2)                                                                           &
#endif
          !$OMP NUM_THREADS(OpenMPOpts%nOutputProcessThreads)                                         &
          !$OMP DEFAULT(None)                                                                         &
          !$OMP PRIVATE(iX,iY,iD)                                                                     &
          !$OMP SHARED(iZ,iS,iC,lBounds,uBounds,AvInt,Field,iTAAvInt,iTAField,DecayFactor,OpenMPOpts)
          Do iY = lBounds(3), uBounds(3)
            Do iX = lBounds(2), uBounds(2)
              Do iD = lBounds(1), uBounds(1)
                AvInt%Std(iD, iX, iY, iZ, iS, iTAAvInt, iC) = &
                AvInt%Std(iD, iX, iY, iZ, iS, iTAAvInt, iC) + &
                Field%Std(iD, iX, iY, iZ, iS, iTAField, iC) * DecayFactor
              End Do
            End Do
          End Do
          !$OMP END PARALLEL DO
        End Do
      End Do
    End Do
    
  Else If (Field%TypeType == '8') Then
    lBounds = lBound(AvInt%P64)
    uBounds = uBound(AvInt%P64)
    Do iC = lBounds(7), uBounds(7)
      Do iS = lBounds(5), uBounds(5)
        Do iZ = lBounds(4), uBounds(4)
          !$OMP PARALLEL DO                                                                           &
#ifdef    UseOpenMP3_0
          !$OMP COLLAPSE(2)                                                                           &
#endif
          !$OMP NUM_THREADS(OpenMPOpts%nOutputProcessThreads)                                         &
          !$OMP DEFAULT(None)                                                                         &
          !$OMP PRIVATE(iX,iY,iD)                                                                     &
          !$OMP SHARED(iZ,iS,iC,lBounds,uBounds,AvInt,Field,iTAAvInt,iTAField,DecayFactor,OpenMPOpts)
          Do iY = lBounds(3), uBounds(3)
            Do iX = lBounds(2), uBounds(2)
              Do iD = lBounds(1), uBounds(1)
                AvInt%P64(iD, iX, iY, iZ, iS, iTAAvInt, iC) = &
                AvInt%P64(iD, iX, iY, iZ, iS, iTAAvInt, iC) + &
                Field%P64(iD, iX, iY, iZ, iS, iTAField, iC) * DecayFactor
              End Do
            End Do
          End Do
          !$OMP END PARALLEL DO
        End Do
      End Do
    End Do
  End If
End Subroutine Field2AvInt

!-------------------------------------------------------------------------------------------------------------

Subroutine Field2Prob(Process, iTAField, iTAProb, DGridProb, Field, Prob, ThresholdConversionFactor)
! Calculate contribution to a probability distribution from one time level of a field.

! Note the contributions are simply added to aid precision. Appropriate scaling must be done later.

! Note any pre-existing D-grid (which cannot be a floating D-grid) is carried forward.

  Implicit None
  ! Argument list:
  Type(Process_), Intent(In)    :: Process   ! Process to be used in generating probability distribution.
  Integer,        Intent(In)    :: iTAField  ! Index in Field corresponding to time of interest.
  Integer,        Intent(In)    :: iTAProb   ! Index in Prob corresponding to time of interest.
  Type(DGrid_),   Intent(In)    :: DGridProb ! D-grid for probability distribution.
  Type(Field_),   Intent(In)    :: Field     ! Contributing field.
  Type(Field_),   Intent(InOut) :: Prob      ! Probability distribution.
  Real,           Optional      :: ThresholdConversionFactor !} Rescaling factor for probability thresholds
                                                             !} if material unit of species differs from that
                                                             !} of the field request.
  ! Locals:
  Integer :: nD    ! Number of points in D-Grid.
  Integer :: nPED  ! Number of points in pre-existing D-Grid (1 if no pre-existing D-Grid).
  Integer :: s     ! Shift (offset) for floating D-Grid.
  Integer :: ds    ! Change in shift (offset) for floating D-Grid.
  Integer :: iD    !} Loop indices.
  Integer :: iX    !}
  Integer :: iY    !}
  Integer :: iZ    !}
  Integer :: iLast !}
  Real    :: TCF   ! Threshold conversion factor, see variable ThresholdConversionFactor

  If (Present(ThresholdConversionFactor)) Then
    TCF = ThresholdConversionFactor
  Else
    TCF = 1.0
  End If
  
  nD = nDInDGrid(DGridProb)

  nPED = Size(Field%Std, 1)

  ! Check don't have pre-existing D-grid with more than one point and a new D-grid with more than one point or
  ! a new floating D-grid.
  If (nPED > 1 .and. (nD > 1 .or. FloatingDGrid(DGridProb))) Then
    Call Message('UNEXPECTED FATAL ERROR in Field2Prob', 4)
  End If

  ! If floating D-Grid, adjust shift (offset).
  If (FloatingDGrid(DGridProb)) Then

    Do iLast = 1, Size(Prob%Std, 7)
    Do iZ    = 1, Size(Prob%Std, 4)
    Do iY    = 1, Size(Prob%Std, 3)
    Do iX    = 1, Size(Prob%Std, 2)

      If (Field%Std(1, iX, iY, iZ, 1, iTAField, 1) > 0) Then

        s = ShiftIndexInDGrid(DGridProb, TCF*Field%Std(1, iX, iY, iZ, 1, iTAField, iLast))

        If (s > Prob%S(iX, iY, iZ, 1, iTAProb, iLast)) Then

          ! Note Prob%S can be -Huge(Prob%S), so next statement phrased to avoid integer overflow.
          ds = Min(s, nD + Prob%S(iX, iY, iZ, 1, iTAProb, iLast)) - Prob%S(iX, iY, iZ, 1, iTAProb, iLast)

          Prob%S(iX, iY, iZ, 1, iTAProb, iLast) = s

          Do iD = 1, nD - ds
            Prob%Std(iD, iX, iY, iZ, 1, iTAProb, iLast) = Prob%Std(iD + ds, iX, iY, iZ, 1, iTAProb, iLast)
          End Do

          Do iD = nD - ds + 1, nD
            Prob%Std(iD, iX, iY, iZ, 1, iTAProb, iLast) = 0.0
          End Do

        End If

      End If

    End Do
    End Do
    End Do
    End Do

  End If

  ! Update Prob with contribution from Field.
  Do iLast = 1, Size(Prob%Std, 7)
  Do iZ    = 1, Size(Prob%Std, 4)
  Do iY    = 1, Size(Prob%Std, 3)
  Do iX    = 1, Size(Prob%Std, 2)

    If (FloatingDGrid(DGridProb)) Then
      s = Prob%S(iX, iY, iZ, 1, iTAProb, iLast)
    Else
      s = 0
    End If

    ! No pre-existing D-Grid (with more than one point).
    If (nPED == 1) Then

      Do iD = 1, nD
        If (TCF*Field%Std(1, iX, iY, iZ, 1, iTAField, 1) > DInDGrid(DGridProb, iD, s)) Then
          Prob%Std(iD, iX, iY, iZ, 1, iTAProb, iLast) = Prob%Std(iD, iX, iY, iZ, 1, iTAProb, iLast) + 1.0
        Else
          Exit
        End If
      End Do

      If (nD > 1 .or. FloatingDGrid(DGridProb)) Then
        Prob%MaxStd(iX, iY, iZ, 1, iTAProb, iLast) = Max(                                            &
                                                   Prob%MaxStd( iX, iY, iZ, 1, iTAProb,  iLast),     &
                                                   TCF*Field%Std(1, iX, iY, iZ, 1, iTAField, iLast)  &
                                                 )
        Prob%MinStd(iX, iY, iZ, 1, iTAProb, iLast) = Min(                                            &
                                                   Prob%MinStd( iX, iY, iZ, 1, iTAProb,  iLast),     &
                                                   TCF*Field%Std(1, iX, iY, iZ, 1, iTAField, iLast)  &
                                                 )
      End If

    ! Pre-existing D-Grid. Here the D-Grid is carried forward.
    Else

      Do iD = 1, nPED
        If (TCF*Field%Std(iD, iX, iY, iZ, 1, iTAField, iLast) > DInDGrid(DGridProb, 1, s)) Then
          Prob%Std(iD, iX, iY, iZ, 1, iTAProb, iLast) = Prob%Std(iD, iX, iY, iZ, 1, iTAProb, iLast) + 1.0
        Else
          Exit
        End If
      End Do

    End If

  End Do
  End Do
  End Do
  End Do

End Subroutine Field2Prob

!-------------------------------------------------------------------------------------------------------------

Subroutine Prob2Percentile(Process, iTAProb, iTAPercent, DGridProb, DGridPercent, Prob, Percent)
! Calculate percentiles from a probability distribution.

! The percentiles are calculated by interpolating the exceedence probability data.

! Note DGridProb must have more than one point or a floating D-Grid. This ensures that there is no D-Grid with
! more than one point or floating D-Grid which exists before the probability processing step.

  Implicit None
  ! Argument list:
  Type(Process_), Intent(In)    :: Process      ! Process to be used in generating percentiles.
  Integer,        Intent(In)    :: iTAProb      ! Index in Prob corresponding to time of interest.
  Integer,        Intent(In)    :: iTAPercent   ! Index in Percent corresponding to time of interest.
  Type(DGrid_),   Intent(In)    :: DGridProb    ! D-grid for probability distribution.
  Type(DGrid_),   Intent(In)    :: DGridPercent ! D-grid for percentiles.
  Type(Field_),   Intent(In)    :: Prob         ! Probability distribution.
  Type(Field_),   Intent(InOut) :: Percent      ! Collection of percentiles.
  ! Locals:
  Integer   :: nDProb    !} Number of points in D-Grids.
  Integer   :: nDPercent !}
  Integer   :: s         ! Shift (offset) for floating D-Grid.
  Integer   :: k         !} Indices, exceedence probabilities and values (percentiles) used in interpolation.
  Integer   :: kMax      !}
  Integer   :: kMin      !}
  Real(Std) :: P         !}
  Real(Std) :: PMax      !}
  Real(Std) :: PMin      !}
  Real(Std) :: Value     !}
  Real(Std) :: MaxValue  !}
  Real(Std) :: MinValue  !}
  Integer   :: iD        !] Loop indices.
  Integer   :: iX        !]
  Integer   :: iY        !]
  Integer   :: iZ        !]
  Integer   :: iLast     !]

  nDProb = nDInDGrid(DGridProb)

  nDPercent = nDInDGrid(DGridPercent)

  ! Check DGridProb has more than one point or is a floating D-Grid.
  If (.not.(nDProb > 1 .or. FloatingDGrid(DGridProb))) Then
    Call Message('UNEXPECTED FATAL ERROR in Prob2Percentile', 4)
  End If

  Do iLast = 1, Size(Percent%Std, 7)
  Do iZ    = 1, Size(Percent%Std, 4)
  Do iY    = 1, Size(Percent%Std, 3)
  Do iX    = 1, Size(Percent%Std, 2)

    Do iD = 1, nDPercent

      P = 1.0 - DInDGrid(DGridPercent, iD, 0) / 100.0
      s = Prob%S(iX, iY, iZ, 1, iTAProb, iLast)

      ! Find indices surrounding point. $$ Need dgrid increasing.
      kMax = nDProb + 1
      kMin = 0
      Do
        k = kMax - kMin
        If (k <= 1) Exit
        k = kMin + k/2
        If (Prob%Std(k, iX, iY, iZ, 1, iTAProb, iLast) < P) Then
          kMax = k
        Else
          kMin = k
        End If
      End Do

      ! Interpolate.
      If (kMin == nDProb) Then
        PMax     = 0.0
        PMin     = Prob%Std(kMin, iX, iY, iZ, 1, iTAProb, iLast)
        MaxValue = Prob%MaxStd(   iX, iY, iZ, 1, iTAProb, iLast)
        MinValue = DInDGrid(DGridProb, kMin, s)
      Else If (kMax == 1) Then
        PMax     = Prob%Std(kMax, iX, iY, iZ, 1, iTAProb, iLast)
        PMin     = 1.0
        MaxValue = DInDGrid(DGridProb, kMax, s)
        MinValue = Prob%MinStd(iX, iY, iZ, 1, iTAProb, iLast)
      Else
        PMax     = Prob%Std(kMax, iX, iY, iZ, 1, iTAProb, iLast)
        PMin     = Prob%Std(kMin, iX, iY, iZ, 1, iTAProb, iLast)
        MaxValue = DInDGrid(DGridProb, kMax, s)
        MinValue = DInDGrid(DGridProb, kMin, s)
      End If
      If (PMax == 0.0) MaxValue = Prob%MaxStd(iX, iY, iZ, 1, iTAProb, iLast)
      If (PMin == 0.0) MinValue = Prob%MaxStd(iX, iY, iZ, 1, iTAProb, iLast)
      If (PMax == 1.0) MaxValue = Prob%MinStd(iX, iY, iZ, 1, iTAProb, iLast)
      If (PMin == 1.0) MinValue = Prob%MinStd(iX, iY, iZ, 1, iTAProb, iLast)
      If (PMax == PMin) Then ! $$ should make similar test in other interpolations (in other routines).
        Value = MaxValue
      Else
        Value = MinValue + (MaxValue - MinValue) * (P - PMin) / (PMax - PMin)
      End If

      Percent%Std(iD, iX, iY, iZ, 1, iTAPercent, iLast) = Value

    End Do

  End Do
  End Do
  End Do
  End Do

End Subroutine Prob2Percentile

!-------------------------------------------------------------------------------------------------------------

Subroutine Field2MaxMin(Process, iTAField, iTAMaxMin, DGridMaxMin, Field, MaxMin)
! Calculate contribution to a max/min from one time level of a field.

! Note any pre-existing D-grid (which cannot be a floating D-grid) is carried forward.

  Implicit None
  ! Argument list:
  Type(Process_), Intent(In)    :: Process     ! Process to be used in generating max/mins.
  Integer,        Intent(In)    :: iTAField    ! Index in Field corresponding to time of interest.
  Integer,        Intent(In)    :: iTAMaxMin   ! Index in MaxMin corresponding to time of interest.
  Type(DGrid_),   Intent(In)    :: DGridMaxMin ! D-grid for max/mins.
  Type(Field_),   Intent(In)    :: Field       ! Contributing field.
  Type(Field_),   Intent(InOut) :: MaxMin      ! Collection of max/mins.
  ! Locals:
  Integer :: nD   ! Number of points in D-Grid.
  Integer :: nPED ! Number of points in pre-existing D-Grid (1 if no pre-existing D-Grid).
  Integer :: iD   ! Loop index.

  nD = nDInDGrid(DGridMaxMin)

  nPED = Size(Field%Std, 1)

  ! Check don't have pre-existing D-grid with more than one point and a new D-grid with more than one point or
  ! a new floating D-grid.
  If (nPED > 1 .and. (nD > 1 .or. FloatingDGrid(DGridMaxMin))) Then
    Call Message('UNEXPECTED FATAL ERROR in Field2MaxMin', 4)
  End If

  ! Check new D-grid isn't floating.
  If (FloatingDGrid(DGridMaxMin)) Then
    Call Message('UNEXPECTED FATAL ERROR in Field2MaxMin', 4)
  End If

  ! No pre-existing D-Grid (with more than one point).
  If (nPED == 1) Then

    Do iD = 1, nD
      If (DInDGrid(DGridMaxMin, iD, 0) == 100.0) Then
        MaxMin%Std(iD, :, :, :, 1, iTAMaxMin, :) = Max(                                        &
                                                     MaxMin%Std(iD, :, :, :, 1, iTAMaxMin, :), &
                                                     Field%Std (1,  :, :, :, 1, iTAField,  :)  &
                                                   )
      Else If (DInDGrid(DGridMaxMin, iD, 0) == 0.0) Then
        MaxMin%Std(iD, :, :, :, 1, iTAMaxMin, :) = Min(                                        &
                                                     MaxMin%Std(iD, :, :, :, 1, iTAMaxMin, :), &
                                                     Field%Std (1,  :, :, :, 1, iTAField,  :)  &
                                                   )
      Else
        Call Message('UNEXPECTED FATAL ERROR in Field2MaxMin', 4)
      End If
    End Do

  ! Pre-existing D-Grid. Here the D-Grid is carried forward.
  Else

    Do iD = 1, nPED
      If (DInDGrid(DGridMaxMin, 1, 0) == 100.0) Then
        MaxMin%Std(iD, :, :, :, 1, iTAMaxMin, :) = Max(                                        &
                                                     MaxMin%Std(iD, :, :, :, 1, iTAMaxMin, :), &
                                                     Field%Std (iD, :, :, :, 1, iTAField,  :)  &
                                                   )
      Else If (DInDGrid(DGridMaxMin, 1, 0) == 0.0) Then
        MaxMin%Std(iD, :, :, :, 1, iTAMaxMin, :) = Min(                                        &
                                                     MaxMin%Std(iD, :, :, :, 1, iTAMaxMin, :), &
                                                     Field%Std (iD, :, :, :, 1, iTAField,  :)  &
                                                   )
      Else
        Call Message('UNEXPECTED FATAL ERROR in Field2MaxMin', 4)
      End If
    End Do

  End If

End Subroutine Field2MaxMin

!-------------------------------------------------------------------------------------------------------------

Subroutine DeriveMixingRatio(             &
             iField, iAirConc, iDensity,  &
             iSpecies, iMaterialUnit,     &
             iTA, Specieses, Results,     &
             MaterialUnits                &
           )
! Calculates volumetric Mixing Ratio field at a given time.

  Implicit None
  ! Argument list:
  Integer,                  Intent(In)    :: iField           !} Index of field request that is processed,
                                                              !} i.e. mixing ratio field request
  Integer,                  Intent(In)    :: iAirConc         ! Index of air concentration field request
  Integer,                  Intent(In)    :: iDensity         ! Index of density field request
  Integer,                  Intent(In)    :: iSpecies         ! index of species
  Integer,                  Intent(In)    :: iMaterialUnit    ! index of Material unit of field request
  Integer,                  Intent(In)    :: iTA              ! Timestep for which conversion is carried out
  Type(Specieses_), Target, Intent(In)    :: Specieses        ! Collection of species
  Type(Results_), Target,   Intent(InOut) :: Results          ! Collection of results
  Type(MaterialUnits_),     Intent(In)    :: MaterialUnits    ! Collection of material units.

  ! Local variables:
  Type(Field_), Pointer     :: FieldMixingRatio
  Type(Field_), Pointer     :: FieldAirConc
  Type(Field_), Pointer     :: FieldDensity
  Type(Species_), Pointer   :: Species
  
  Integer                   :: iLast

  Type(MaterialUnit_)       :: SpeciesMaterialUnit
  Type(MaterialUnit_)       :: mgMaterialUnit
  Type(MaterialUnit_)       :: FieldReqMaterialUnit
  Type(MaterialUnit_)       :: ppmMaterialUnit

  FieldMixingRatio => Results%Fields(iField   )
  FieldAirConc     => Results%Fields(iAirConc )
  FieldDensity     => Results%Fields(iDensity )

  Species => Specieses%Specieses(iSpecies)

! Work out conversion factors
  SpeciesMaterialUnit = MaterialUnits%MaterialUnits(Species%iMaterialUnit)
  mgMaterialUnit = MaterialUnits%MaterialUnits(FindMaterialUnitIndex("mg",MaterialUnits))

  FieldReqMaterialUnit = MaterialUnits%MaterialUnits(iMaterialUnit)
  ppmMaterialUnit = MaterialUnits%MaterialUnits(FindMaterialUnitIndex("ppm",MaterialUnits))
    
  Do iLast = 1, Size(FieldAirConc%Std  (:, :, :, :, :, :, :),7)
    FieldMixingRatio%Std(1, :, :, :, 1, iTA, iLast) = FieldAirConc%Std  (1, :, :, :, 1, iTA, iLast) &
                                                    / FieldDensity%Std  (1, :, :, :, 1, iTA, 1)     &
                                                    * MoleMassAir / Species%MolecularWeight         &
                                                    * SpeciesMaterialUnit%ConversionFactor          &
                                                    / mgMaterialUnit%ConversionFactor               &
                                                    * FieldReqMaterialUnit%ConversionFactor         &
                                                    / ppmMaterialUnit%ConversionFactor
  End Do

End Subroutine DeriveMixingRatio

!-------------------------------------------------------------------------------------------------------------

Subroutine DeriveEMixingRatio(             &
             iField, iChemField, iDensity, &
             iSpecies, iMaterialUnit,      &
             iTA, Specieses, Results,      &
             MaterialUnits                 &
           )
! Calculates volumetric Mixing Ratio field at a given time for chemistry field species on Eulerian grid.

  Implicit None
  ! Argument list:
  Integer,                      Intent(In)    :: iField        !} Index of field request that is processed,
                                                               !} i.e. mixing ratio field request
  Integer,                      Intent(In)    :: iChemField    ! Index of air concentration field request
  Integer,                      Intent(In)    :: iDensity      ! Index of density field request
  Integer,                      Intent(In)    :: iSpecies      ! index of species
  Integer,                      Intent(In)    :: iMaterialUnit ! index of Material unit of field request
  Integer,                      Intent(In)    :: iTA           ! Timestep for which conversion is carried out
  Type(Specieses_),     Target, Intent(In)    :: Specieses     ! Collection of species
  Type(Results_),       Target, Intent(InOut) :: Results       ! Collection of results
  Type(MaterialUnits_),         Intent(In)    :: MaterialUnits ! Collection of material units.


  ! Local variables:
  Type(Field_),   Pointer   :: FieldMixingRatio
  Type(Field_),   Pointer   :: FieldAirConc
  Type(Field_),   Pointer   :: FieldDensity
  Type(Species_), Pointer   :: Species
  
  Integer                   :: iLast

  Type(MaterialUnit_)       :: SpeciesMaterialUnit
  Type(MaterialUnit_)       :: mgMaterialUnit
  Type(MaterialUnit_)       :: FieldReqMaterialUnit
  Type(MaterialUnit_)       :: ppmMaterialUnit

  FieldMixingRatio => Results%Fields(iField   )
  FieldAirConc     => Results%Fields(iChemField )
  FieldDensity     => Results%Fields(iDensity )

  Species => Specieses%Specieses(iSpecies)

! Work out conversion factors
  SpeciesMaterialUnit = MaterialUnits%MaterialUnits(Species%iMaterialUnit)
  mgMaterialUnit = MaterialUnits%MaterialUnits(FindMaterialUnitIndex("mg",MaterialUnits))

  FieldReqMaterialUnit = MaterialUnits%MaterialUnits(iMaterialUnit)
  ppmMaterialUnit = MaterialUnits%MaterialUnits(FindMaterialUnitIndex("ppm",MaterialUnits))

  Do iLast = 1, Size(FieldAirConc%Std  (:, :, :, :, :, :, :),7)
    FieldMixingRatio%Std(1, :, :, :, 1, iTA, iLast) = FieldAirConc%Std  (1, :, :, :, 1, iTA, iLast) &
                                                    / FieldDensity%Std  (1, :, :, :, 1, iTA, 1)     &
                                                    * MoleMassAir / Species%MolecularWeight         &
                                                    * SpeciesMaterialUnit%ConversionFactor          &
                                                    / mgMaterialUnit%ConversionFactor               &
                                                    * FieldReqMaterialUnit%ConversionFactor         &
                                                    / ppmMaterialUnit%ConversionFactor
  End Do

End Subroutine DeriveEMixingRatio

!-------------------------------------------------------------------------------------------------------------

Subroutine DeriveAreaAtRisk(                                       &
             iField, iAirConc,                                     &
             iTA, iTAAirConc,                                      &
             Coords, Grids, Domains, Flows, Sources, Reqs, Results &
           )
! Calculates Area At Risk field at one time.

  Implicit None
  ! Argument list:
  Integer,        Intent(In)            :: iField      ! Index of area at risk field.
  Integer,        Intent(In)            :: iAirConc    ! Index of air concentration field.
  Integer,        Intent(In)            :: iTA         ! Index for T in area at risk field.
  Integer,        Intent(In)            :: iTAAirConc  ! Index for T in air concentration field.
  Type(Coords_),  Intent(In),    Target :: Coords      ! Collection of coord systems.
  Type(Grids_),   Intent(In),    Target :: Grids       ! Collection of grids.
  Type(Domains_), Intent(In)            :: Domains     ! Collection of domains.
  Type(Flows_),   Intent(In)            :: Flows       ! Collection of flow module instance states.
  Type(Sources_), Intent(In),    Target :: Sources     ! Collection of sources.
  Type(Reqs_),    Intent(In),    Target :: Reqs        ! Collection of requirements.
  Type(Results_), Intent(InOut), Target :: Results     ! Collection of results.
  ! Locals:
  Integer                      :: iHGrid          ! Index for horizontal grid.
  Integer                      :: iZGrid          ! Index for vertical grid.
  Integer                      :: nX              !} Size of horizontal grid.
  Integer                      :: nY              !}
  Integer                      :: nZ              ! Size of vertical grid.
  Integer                      :: iX              ! Index for X.
  Integer                      :: iY              ! Index for Y.
  Integer                      :: iZ              ! Index for Z.
  Real(Std)                    :: Value           !
  Real(Std),       Allocatable :: MaxAirConc(:)   !
  Real(Std)                    :: HMax            !
  Real(Std)                    :: H1              !
  Real(Std)                    :: H2              !
  Real(Std)                    :: SourceX(2)      !
  Real(Std)                    :: X(2)            !
  Real(Std)                    :: dR              !
  Integer                      :: nR              !
  Integer                      :: iR              !
  Integer                      :: iTheta          !
  Type(HCoeffs_)               :: HCoeffs         !
  Type(HCoord_),   Pointer     :: HCoord          !} Abbreviations for coord systems, grids, sources,
  Type(HCoord_),   Pointer     :: HCoordSource    !} requirements and results.
  Type(HGrid_),    Pointer     :: HGrid           !}
  Type(Source_),   Pointer     :: Source          !}
  Type(FieldReq_), Pointer     :: FieldReq        !}
  Type(FieldReq_), Pointer     :: FieldReqAirConc !}
  Type(Field_),    Pointer     :: Field           !}
  Type(Field_),    Pointer     :: FieldAirConc    !}
  ! Local parameters:
  Integer, Parameter :: nTheta = 3600

  FieldReq        => Reqs%FieldReqs(iField)
  FieldReqAirConc => Reqs%FieldReqs(iAirConc)
  Field           => Results%Fields(iField)
  FieldAirConc    => Results%Fields(iAirConc)

  iHGrid = FieldReq%iHGrid
  iZGrid = FieldReq%iZGrid

  If (iHGrid /= 0) Then
    nX = Grids%HGrids(iHGrid)%nX
    nY = Grids%HGrids(iHGrid)%nY
    HGrid  => Grids%HGrids(iHGrid)
    HCoord => Coords%HCoords(HGrid%iHCoord)
  Else
    nX = 1
    nY = 1
  End If
  If (iZGrid /= 0) Then
    nZ = Grids%ZGrids(iZGrid)%nZ
  Else
    nZ = 1
  End If

  ! Identify source. Note that, if no sources, result is zero (and will have been set to zero earlier).
  If (FieldReq%iSource > 0) Then
    Source => Sources%Sources(FieldReq%iSource)
  Else If (FieldReq%iSourceGroup > 0) Then
    ! $$ not yet supported.
    Call Message ('FATAL ERROR: Area at risk only works for single sources.', 3)
  Else If (Sources%nSources == 1) Then
    Source => Sources%Sources(1)
  Else If (Sources%nSources == 0) Then
    Return
  Else
    Call Message ('FATAL ERROR: Area at risk only works for single sources.', 3)
    ! $$ check earlier (and ideally remove restriction by setting up multiple airconc requests)
  End If

  ! No horizontal grid.
  If (iHGrid == 0) Then

    Do iZ = 1, nZ
      If (FieldAirConc%Std(1, 1, 1, iZ, 1, iTAAirConc, 1) == 0.0) Then
        Value = 0.0
      Else
        Value = 1.0
      End If
      ! $$ Note this has same effect as Field%Std(1, iX, iY, iZ, 1, iTA, 1) = Value, but anticipates
      ! combining results for multiple sources.
      Field%Std(1, 1, 1, iZ, 1, iTA, 1) = Max(Value, Field%Std(1, 1, 1, iZ, 1, iTA, 1))
    End Do

  ! Horizontal grid.
  Else

    ! $$ check earlier. Ideally remove this restriction by setting up extra regular grid (difficult).
    If (HGrid%Unstructured .or. HGrid%Variable) Then
      Call Message ('FATAL ERROR: Area at risk only works for regular grids.', 3)
    End If

    HCoordSource => Coords%HCoords(Source%iHCoord)

    SourceX(:) = ConvertH(HCoordSource, HCoord, Source%X(1:2))

    ! $$ Could do better than using h1, h2 at long range.
    Call MetricCoeffs(HCoord, SourceX, HMax, H1, H2)

    dR = Min(H1 * HGrid%dX, H2 * HGrid%dY)

    nR = Ceiling(                                                              &
           Sqrt(                                                               &
             H1**2 * Max(                                                      &
                       (HGrid%X0 + (HGrid%nX - 1) * HGrid%dX - SourceX(1))**2, &
                       (SourceX(1) - HGrid%X0)**2                              &
                     )                                                         &
             +                                                                 &
             H2**2 * Max(                                                      &
                       (HGrid%Y0 + (HGrid%nY - 1) * HGrid%dY - SourceX(2))**2, &
                       (SourceX(2) - HGrid%Y0)**2                              &
                     )                                                         &
           )                                                                   &
           / dR                                                                &
         )

    Allocate(MaxAirConc(nR))

    Do iZ = 1, nZ

      MaxAirConc(:) = 0.0
      Do iR = 1, nR
      Do iTheta = 1, nTheta
        X(1) = SourceX(1) + dR * Real(iR) * Cos(2.0 * Pi * Real(iTheta) / Real(nTheta)) / H1
        X(2) = SourceX(2) + dR * Real(iR) * Sin(2.0 * Pi * Real(iTheta) / Real(nTheta)) / H2
        Call GetHCoeffs(X(1), X(2), HGrid, HCoeffs)
        MaxAirConc(iR) = Max(                                                                 &
                           MaxAirConc(iR),                                                    &
                           InterpXY(HCoeffs, FieldAirConc%Std(1, :, :, iZ, 1, iTAAirConc, 1)) &
                         )
      End Do
      End Do

      ! Loop over grid points.
      Do iX = 1, nX
      Do iY = 1, nY

        X(1) = HGrid%X0 + Real(iX - 1) * HGrid%dX - SourceX(1)
        X(2) = HGrid%Y0 + Real(iY - 1) * HGrid%dY - SourceX(2)
        iR = NInt(Sqrt((H1 * X(1))**2 + (H2 * X(2))**2) / dR)
        iR = Min(Max(iR, 1), nR)

        ! Calculate "Area at risk" (= sigma_y / y for Gaussian plume).
        Value = FieldAirConc%Std(1, iX, iY, iZ, 1, iTAAirConc, 1)
        If (Value == 0.0 .or. MaxAirConc(iR) == 0.0) Then
          Value = 0.0
        Else
          Value = Value / MaxAirConc(iR)
          Value = Min(Value, 1.0)
          Value = 1.0 / Max(Sqrt(- 2.0 * Log(Value)), 1.0)
        End If

        ! $$ Note this has same effect as Field%Std(1, iX, iY, iZ, 1, iTA, 1) = Value, but anticipates
        ! combining results for multiple sources.
        Field%Std(1, iX, iY, iZ, 1, iTA, 1) = Max(Value, Field%Std(1, iX, iY, iZ, 1, iTA, 1))

      End Do
      End Do

    End Do

    DeAllocate(MaxAirConc)

  End If

End Subroutine DeriveAreaAtRisk

!-------------------------------------------------------------------------------------------------------------

Subroutine DeriveCloudGammaDose(                     &
             iField, iPhotonFlux,                    &
             iT, iTA, Grids,                         &
             Reqs, Results,                          &
             Specieses, CloudGammaParamses           &
           )
! Calculates Cloud Gamma Dose field at one time.

  Implicit None
  ! Argument list:
  Integer,                   Intent(In)            :: iField             ! Index of cloud gamma dose field.
  Integer,                   Intent(In)            :: iPhotonFlux        ! Index of photon flux field.
  Integer,                   Intent(In)            :: iT                 ! Index for T.
  Integer,                   Intent(In)            :: iTA                ! Index for T in a field array.
  Type(Grids_),              Intent(In),    Target :: Grids              ! Collection of grids.
  Type(Reqs_),               Intent(In),    Target :: Reqs                ! Collection of requirements.
  Type(Results_),            Intent(InOut), Target :: Results            ! Collection of results.
  Type(Specieses_),          Intent(In),    Target :: Specieses          ! Collection of species.
  Type(CloudGammaParamses_), Intent(In),    Target :: CloudGammaParamses ! Collection of sets of cloud
                                                                         ! gamma parameters.
  ! Locals:
  Integer                          :: iHGrid              ! Index for horizontal grid.
  Integer                          :: iZGrid              ! Index for vertical grid.
  Integer                          :: iX                  ! Index for X.
  Integer                          :: iY                  ! Index for Y.
  Integer                          :: iZ                  ! Index for Z.
  Integer                          :: iSpecies            ! Index for Species.
  Integer                          :: nEnergies           ! The number of photon energies associated
                                                          ! to a particular species.
  Integer                          :: iE                  ! Index for looping over nEnergies.
  Real(Std)                        :: a                   ! Local value of AirKermapuFluence.
  Real(Std)                        :: b                   ! Local value of DosepuAirKerma.
  Real(Std)                        :: PhotonFlux          ! Photon flux.
  Real(Std)                        :: CloudGammaDose      ! Effective or Organ Cloud Gamma Dose.
  Type(HGrid_),            Pointer :: HGrid               !} Abbreviations for grids, requirements, results,
  Type(ZGrid_),            Pointer :: ZGrid               !} species and collections of cloud gamma
                                                          !} parameters.
  Type(FieldReq_),         Pointer :: FieldReq            !}
  Type(FieldReq_),         Pointer :: FieldReqPhotonFlux  !}
  Type(Field_),            Pointer :: Field               !}
  Type(Field_),            Pointer :: FieldPhotonFlux     !}
  Type(Species_),          Pointer :: Species             !}
  Type(CloudGammaParams_), Pointer :: CloudGammaParams    !}

  FieldReq           => Reqs%FieldReqs(iField)
  ! $$ Don't think FieldReqPhotonFlux is used.
  FieldReqPhotonFlux => Reqs%FieldReqs(iPhotonFlux)
  Field              => Results%Fields(iField)
  FieldPhotonFlux    => Results%Fields(iPhotonFlux)

  iSpecies = FieldReq%iSpecies
  iHGrid = FieldReq%iHGrid
  iZGrid = FieldReq%iZGrid

  Species => Specieses%Specieses(iSpecies)
  If (Species%iCloudGammaParams == 0) Return
  CloudGammaParams => CloudGammaParamses%CloudGammaParamses(Species%iCloudGammaParams)

  If (iHGrid /= 0) Then
    HGrid => Grids%HGrids(iHGrid)
  Else
    Call Message('FATAL ERROR in DeriveCloudGammaDose: horizontal grid required for '   // &
                 'cloud gamma dose calculation', 3)
  End If

  If (iZGrid /= 0) Then
    ZGrid => Grids%ZGrids(iZGrid)
  Else
    Call Message('FATAL ERROR in DeriveCloudGammaDose: vertical grid required for '     // &
                 'cloud gamma dose calculation', 3)
  End If

  ! Loop over vertical grid levels
  Do iZ = 1, ZGrid%nZ

    If (HGrid%Unstructured) Then

    ! Handle unstructured grids
      iY = 1
      Do iX = 1, HGrid%nX

        CloudGammaDose = 0.0
        Field%Std(1, iX, iY, iZ, 1, iTA, 1) = 0.0

        nEnergies = CloudGammaParams%nEnergies
        ! Loop over all photon energies of the species of interest
        Do iE = 1, nEnergies

          ! Defining the Air Kerma per unit fleunce (photon fluxes)
          a = 0.0
          a = CloudGammaParams%AirKermapuFluence(iE)


          ! Defining the Adult Effective/Organ Cloud Gamma Dose per unit Air Kerma
          b = 0.0
          If (FieldReq%iQuantity == Q_AduEffCloudGammaDose) Then
            b = CloudGammaParams%AdEffDosepuAirKerma(iE)
          Else If (FieldReq%iQuantity == Q_AduLunCloudGammaDose) Then
            b = CloudGammaParams%AdLunDosepuAirKerma(iE)
          Else If (FieldReq%iQuantity == Q_AduThyCloudGammaDose) Then
            b = CloudGammaParams%AdThyDosepuAirKerma(iE)
          Else If (FieldReq%iQuantity == Q_AduBoSCloudGammaDose) Then
            b = CloudGammaParams%AdBoSDosepuAirKerma(iE)
          End If

          PhotonFlux = FieldPhotonFlux%Std(1, iX, iY, iZ, 1, iTA, iE)

          ! Calculating Adult Effective Cloud Gamma Dose i.e. the whole body dose (Sv s^-1)
          ! or Adult Organ Cloud Gamma Dose (Gy s^-1)
          CloudGammaDose = PhotonFlux * a * b
          Field%Std(1, iX, iY, iZ, 1, iTA, 1) = Field%Std(1, iX, iY, iZ, 1, iTA, 1) + CloudGammaDose
        End Do
      End Do

    Else

      ! Handle structured grids
      Do iY = 1, HGrid%nY
        Do iX = 1, HGrid%nX

          CloudGammaDose = 0.0
          Field%Std(1, iX, iY, iZ, 1, iTA, 1) = 0.0

          nEnergies = CloudGammaParams%nEnergies
          ! Loop over all photon energies of the species of interest
          Do iE = 1, nEnergies

            ! Defining the Air Kerma per unit fleunce (photon fluxes)
            a = 0.0
            a = CloudGammaParams%AirKermapuFluence(iE)

            ! Defining the Adult Effective/Organ Cloud Gamma Dose per unit Air Kerma
            b = 0.0
            If (FieldReq%iQuantity == Q_AduEffCloudGammaDose) Then
              b = CloudGammaParams%AdEffDosepuAirKerma(iE)
            Else If (FieldReq%iQuantity == Q_AduLunCloudGammaDose) Then
              b = CloudGammaParams%AdLunDosepuAirKerma(iE)
            Else If (FieldReq%iQuantity == Q_AduThyCloudGammaDose) Then
              b = CloudGammaParams%AdThyDosepuAirKerma(iE)
            Else If (FieldReq%iQuantity == Q_AduBoSCloudGammaDose) Then
              b = CloudGammaParams%AdBoSDosepuAirKerma(iE)
            End If

            PhotonFlux = FieldPhotonFlux%Std(1, iX, iY, iZ, 1, iTA, iE)

            ! Calculating Adult Effective Cloud Gamma Dose i.e. the whole body dose (Sv s^-1)
            ! or Adult Organ Cloud Gamma Dose (Gy s^-1)
            CloudGammaDose = PhotonFlux * a * b
            Field%Std(1, iX, iY, iZ, 1, iTA, 1) = Field%Std(1, iX, iY, iZ, 1, iTA, 1) + CloudGammaDose
          End Do
        End Do
      End Do

    End If

  End Do

End Subroutine DeriveCloudGammaDose

!-------------------------------------------------------------------------------------------------------------

Subroutine Derive_CloudGammaDose_via_SemiInfApprox( &
             iField, iAirConc,                      &
             iT, iTA, Grids,                        &
             Reqs, Results,                         &
             Specieses, CloudGammaParamses          &
           )
! Calculates Cloud Gamma Dose field using the semi-infinite cloud approximation, at one time.

  Implicit None
  ! Argument list:
  Integer,                   Intent(In)            :: iField             ! Index of cloud gamma dose field.
  Integer,                   Intent(In)            :: iAirConc           ! Index of air conc field.
  Integer,                   Intent(In)            :: iT                 ! Index for T.
  Integer,                   Intent(In)            :: iTA                ! Index for T in a field array.
  Type(Grids_),              Intent(In),    Target :: Grids              ! Collection of grids.
  Type(Reqs_),               Intent(In),    Target :: Reqs               ! Collection of requirements.
  Type(Results_),            Intent(InOut), Target :: Results            ! Collection of results.
  Type(Specieses_),          Intent(In),    Target :: Specieses          ! Collection of species.
  Type(CloudGammaParamses_), Intent(In),    Target :: CloudGammaParamses ! Collection of sets of cloud
                                                                         ! gamma parameters.
  ! Locals:
  Integer                          :: iHGrid              ! Index for horizontal grid.
  Integer                          :: iZGrid              ! Index for vertical grid.
  Integer                          :: iX                  ! Index for X.
  Integer                          :: iY                  ! Index for Y.
  Integer                          :: iZ                  ! Index for Z.
  Integer                          :: iSpecies            ! Index for Species.
  Integer                          :: nEnergies           ! The number of photon energies associated
                                                          ! to a particular species.
  Integer                          :: iE                  ! Index for looping over nEnergies.
  Real(Std)                        :: k                   ! Conversion factor (Gy s^-1 per MeV m^-3 s^-1).
  Real(Std)                        :: E                   ! Photon Energy
  Real(Std)                        :: I                   ! Photon Intensity
  Real(Std)                        :: a                   ! Sum of the product of the Photon Energy & Photon
                                                          ! Intensity
  Real(Std)                        :: b                   ! Local value of DosepuAirKerma.
  Real(Std)                        :: AirConc             ! Air concentration.
  Real(Std)                        :: CloudGammaDose      ! Effective or Organ Cloud Gamma Dose.
  Type(HGrid_),            Pointer :: HGrid               !} Abbreviations for grids, requirements, results,
  Type(ZGrid_),            Pointer :: ZGrid               !} species and collections of cloud gamma
                                                          !}parameters.
  Type(FieldReq_),         Pointer :: FieldReq            !}
  Type(FieldReq_),         Pointer :: FieldReqAirConc     !}
  Type(Field_),            Pointer :: Field               !}
  Type(Field_),            Pointer :: FieldAirConc        !}
  Type(Species_),          Pointer :: Species             !}
  Type(CloudGammaParams_), Pointer :: CloudGammaParams    !}

  ! Defining Conversion factor
  k = 6.3376E-14

  FieldReq        => Reqs%FieldReqs(iField)
  ! $$ Don't think FieldReqAirConc is used.
  FieldReqAirConc => Reqs%FieldReqs(iAirConc)
  Field           => Results%Fields(iField)
  FieldAirConc    => Results%Fields(iAirConc)

  iSpecies = FieldReq%iSpecies
  iHGrid = FieldReq%iHGrid
  iZGrid = FieldReq%iZGrid

  Species => Specieses%Specieses(iSpecies)
  If (Species%iCloudGammaParams == 0) Return
  CloudGammaParams => CloudGammaParamses%CloudGammaParamses(Species%iCloudGammaParams)

  If (iHGrid /= 0) Then
    HGrid => Grids%HGrids(iHGrid)
  Else
    Call Message('FATAL ERROR in Derive_CloudGammaDose_via_SemiInfApprox: '   // &
                 'horizontal grid required for cloud gamma dose calculation', 3)
  End If

  If (iZGrid /= 0) Then
    ZGrid => Grids%ZGrids(iZGrid)
  Else
    Call Message('FATAL ERROR in Derive_CloudGammaDose_via_SemiInfApprox: '     // &
                 'vertical grid required for cloud gamma dose calculation', 3)
  End If

  ! Loop over vertical grid levels
  Do iZ = 1, ZGrid%nZ

    If (HGrid%Unstructured) Then
    ! Handle unstructured grids
      iY = 1
      Do iX = 1, HGrid%nX

        CloudGammaDose = 0.0
        Field%Std(1, iX, iY, iZ, 1, iTA, 1) = 0.0
        AirConc = 0.0
        AirConc = FieldAirConc%Std(1, iX, iY, iZ, 1, iTA, 1)
        a = 0.0

        nEnergies = CloudGammaParams%nEnergies
        ! Loop over all photon energies of the species of interest
        Do iE = 1, nEnergies

          ! Defining the photon energy
          E = 0.0
          E = CloudGammaParams%PhotonEnergy(iE)

          ! Defining the photon intensity
          I = 0.0
          I = CloudGammaParams%PhotonIntensity(iE)

          ! Defining the Adult Effective/Organ Cloud Gamma Dose per unit Air Kerma
          b = 0.0
          If (FieldReq%iQuantity == Q_AduEffCloudGammaDose) Then
            b = CloudGammaParams%AdEffDosepuAirKerma(iE)
          Else If (FieldReq%iQuantity == Q_AduLunCloudGammaDose) Then
            b = CloudGammaParams%AdLunDosepuAirKerma(iE)
          Else If (FieldReq%iQuantity == Q_AduThyCloudGammaDose) Then
            b = CloudGammaParams%AdThyDosepuAirKerma(iE)
          Else If (FieldReq%iQuantity == Q_AduBoSCloudGammaDose) Then
            b = CloudGammaParams%AdBoSDosepuAirKerma(iE)
          End If

          ! Calculating the sum of the product of the photon energy and photon intensity
          a = a + (E * I * b)

        End Do
        ! Calculating Adult Effective Cloud Gamma Dose Rate i.e. the whole body dose (Sv s^-1)
        ! or Adult Organ Cloud Gamma Dose Rate (Gy s^-1)
        CloudGammaDose = k * AirConc * a
        Field%Std(1, iX, iY, iZ, 1, iTA, 1) = CloudGammaDose

      End Do
    Else

    ! Handle structured grids
      Do iY = 1, HGrid%nY
        Do iX = 1, HGrid%nX

          CloudGammaDose = 0.0
          Field%Std(1, iX, iY, iZ, 1, iTA, 1) = 0.0
          AirConc = 0.0
          AirConc = FieldAirConc%Std(1, iX, iY, iZ, 1, iTA, 1)
          a = 0.0

          nEnergies = CloudGammaParams%nEnergies
          ! Loop over all photon energies of the species of interest
          Do iE = 1, nEnergies

            ! Defining the photon energy
            E = 0.0
            E = CloudGammaParams%PhotonEnergy(iE)

            ! Defining the photon intensity
            I = 0.0
            I = CloudGammaParams%PhotonIntensity(iE)

            ! Defining the Adult Effective/Organ Cloud Gamma Dose per unit Air Kerma
            b = 0.0
            If (FieldReq%iQuantity == Q_AduEffCloudGammaDose) Then
              b = CloudGammaParams%AdEffDosepuAirKerma(iE)
            Else If (FieldReq%iQuantity == Q_AduLunCloudGammaDose) Then
              b = CloudGammaParams%AdLunDosepuAirKerma(iE)
            Else If (FieldReq%iQuantity == Q_AduThyCloudGammaDose) Then
              b = CloudGammaParams%AdThyDosepuAirKerma(iE)
            Else If (FieldReq%iQuantity == Q_AduBoSCloudGammaDose) Then
              b = CloudGammaParams%AdBoSDosepuAirKerma(iE)
            End If

            ! Calculating the sum of the product of the photon energy and photon intensity
            a = a + (E * I * b)

          End Do
          ! Calculating Adult Effective Cloud Gamma Dose Rate i.e. the whole body dose (Sv s^-1)
          ! or Adult Organ Cloud Gamma Dose Rate (Gy s^-1)
          CloudGammaDose = k * AirConc * a
          Field%Std(1, iX, iY, iZ, 1, iTA, 1) = CloudGammaDose
        End Do
      End Do
    End If

  End Do

End Subroutine Derive_CloudGammaDose_via_SemiInfApprox

!-------------------------------------------------------------------------------------------------------------

Subroutine DeriveXStats(Grids, iSGrid, iTA, Field)
!

  Implicit None
  ! Argument list:
  Type(Grids_), Intent(In)    :: Grids  !
  Integer,      Intent(In)    :: iSGrid !
  Integer,      Intent(In)    :: iTA    !
  Type(Field_), Intent(InOut) :: Field  !
  ! Locals:
  Integer :: iS ! Loop index.

  Do iS = 1, Grids%TGrids(iSGrid)%nT
    If (Field%Std(1, 1, 1, 1, iS, iTA, 1) == 0.0) Then
      Field%Std(1, 1, 1, 1, iS, iTA, :) = 0.0  ! this should be true already.
    Else
      ! Scale displacement stats by total mass of particles.
      Field%Std(1, 1, 1, 1, iS, iTA, 2:) = Field%Std(1, 1, 1, 1, iS, iTA, 2:) / &
                                           Field%Std(1, 1, 1, 1, iS, iTA, 1)
      ! Remove mean values from displacement stats.
      Field%Std(1, 1, 1, 1, iS, iTA, 5)  = Max(                                   &
                                             0.0,                                 &
                                             Field%Std(1, 1, 1, 1, iS, iTA, 5) -  &
                                             Field%Std(1, 1, 1, 1, iS, iTA, 2)**2 &
                                           )
      Field%Std(1, 1, 1, 1, iS, iTA, 6)  = Max(                                   &
                                             0.0,                                 &
                                             Field%Std(1, 1, 1, 1, iS, iTA, 6) -  &
                                             Field%Std(1, 1, 1, 1, iS, iTA, 3)**2 &
                                           )
      Field%Std(1, 1, 1, 1, iS, iTA, 7)  = Max(                                   &
                                             0.0,                                 &
                                             Field%Std(1, 1, 1, 1, iS, iTA, 7) -  &
                                             Field%Std(1, 1, 1, 1, iS, iTA, 4)**2 &
                                           )
      Field%Std(1, 1, 1, 1, iS, iTA, 8)  = Field%Std(1, 1, 1, 1, iS, iTA,  8) - &
                                           Field%Std(1, 1, 1, 1, iS, iTA,  2) * &
                                           Field%Std(1, 1, 1, 1, iS, iTA,  3)
      Field%Std(1, 1, 1, 1, iS, iTA, 9)  = Field%Std(1, 1, 1, 1, iS, iTA,  9) - &
                                           Field%Std(1, 1, 1, 1, iS, iTA,  3) * &
                                           Field%Std(1, 1, 1, 1, iS, iTA,  4)
      Field%Std(1, 1, 1, 1, iS, iTA, 10) = Field%Std(1, 1, 1, 1, iS, iTA, 10) - &
                                           Field%Std(1, 1, 1, 1, iS, iTA,  4) * &
                                           Field%Std(1, 1, 1, 1, iS, iTA,  2)
    End If
  End Do

End Subroutine DeriveXStats

!-------------------------------------------------------------------------------------------------------------

Subroutine DeriveSigmaCField(                             &
             i, iMeanConc, iTravelTime, iXStats, iT, iTA, &
             Coords, Grids, Domains, Sources, Reqs,       &
             Units, Mets, Flows, Results                  &
           )
! Calculates Sigma C field at any particular time.

  Implicit None
  ! Argument list:
  Integer,        Intent(In)           :: i           ! Index of Sigma C field.
  Integer,        Intent(In)           :: iMeanConc   ! Index of mean concentration field.
  Integer,        Intent(In)           :: iTravelTime ! Index of mean travel time field.
  Integer,        Intent(In)           :: iXStats     ! Index of X Stats field.
  Integer,        Intent(In)           :: iT          ! Index for T.
  Integer,        Intent(In)           :: iTA         ! Index for T in a field array.
  Type(Coords_),  Intent(In)           :: Coords      !
  Type(Grids_),   Intent(In)           :: Grids       !
  Type(Domains_), Intent(In)           :: Domains     !
  Type(Sources_), Intent(In)           :: Sources
  Type(Reqs_),    Intent(In),   Target :: Reqs        !
  Type(Units_),   Intent(InOut)        :: Units       !
  Type(Mets_),    Intent(InOut)        :: Mets        !
  Type(Flows_),   Intent(InOut)        :: Flows       !
  Type(Results_), Intent(InOut)        :: Results     !
  ! Locals:
  Integer           :: iHGrid         ! Index for horizontal grid.
  Integer           :: iZGrid         ! Index for vertical grid.
  Integer           :: iSGrid         ! Index for travel-time grid of X-Stats field.
  Integer           :: nX             !} Size of horizontal grid.
  Integer           :: nY             !}
  Integer           :: nZ             ! Size of vertical grid.
  Integer           :: nT             ! Size of travel-time grid.
  Integer           :: iX             ! Index for X.
  Integer           :: iY             ! Index for Y.
  Integer           :: iZ             ! Index for Z.
  Integer           :: iS             ! Index for S.
  Real(Std)         :: GridPoint(3)   ! Coords of 'receptor' grid-point.
  Real(Std)         :: MeanConc       ! Mean concentration at grid-point.
  Real(Std)         :: TravelTime     ! Mean travel time at grid-point (as real).
  Type(ShortTime_)  :: TravelTimeShort ! Mean travel time at grid-point (as short time).
  Logical           :: DualTimes      ! Plume interpolated from X-Stats at two times.
  Real(Std)         :: S1             !} Nearest travel-times on the travel-time grid
  Real(Std)         :: S2             !} (S1 is not defined if DualTimes is false).
  Real(Std)         :: InterpCoeff    ! Interpolation coefficient if two travel-times.
  Real(Std)         :: XStats1(10)    !} Displacement statistics of the plume
  Real(Std)         :: XStats2(10)    !} at travel-times S1 and S2, respectively.
  Real(Std)         :: PlumeCentre(3) ! Coords of centre-of-mass of the plume.
  Real(Std)         :: YOffset        !} Cross-wind and vertical offsets of grid-point
  Real(Std)         :: ZOffset        !} relative to the plume centre-of-mass.
  Real(Std)         :: HorzWind       !] Met data evaluated at the plume centre-of-mass.
  Real(Std)         :: SigmaVel2      !]
  Real(Std)         :: Epsilon        !]
  Real(Std)         :: Z0             !]
  Real(Std)         :: Lc             ! Length scale of the in-plume fluctuations.
  Real(Std)         :: Zb             ! Estimate of advection=diffusion height.
  Real(Std)         :: EffHorzWindAtReceptor ! Effective advection speed at grid-point.
  Real(Std)         :: SigmaY         !} Estimates of the lateral and vertical spread
  Real(Std)         :: SigmaZ         !} of the plume.
  Type(FlowMemory_) :: FlowMemory     !
  Type(Position_)   :: Position       !
  Type(Flow_)       :: FlowPlume      ! Flow information at the plume centre-of-mass.
  Type(Flow_)       :: FlowReceptor   ! Flow information at the receptor grid-point.
  Integer           :: ErrorCode

  Type(FieldReq_), Pointer :: FieldReq !}
  Type(Process_),  Pointer :: Process  !}
  Type(Time_)              :: ProcAvTime ! time needed for processing $$ replace by shorttime? elsewhere too?
  Integer           :: j

  Type(Cloud_)   :: Cloud   !} Dummy variables for arguments in GetAttrib.
  Type(Rain_)    :: Rain    !}
  Type(Surface_) :: Surface !}
  Type(Soil_)    :: Soil    !}
  Type(Plant_)   :: Plant   !}

  ! Parameters of the fluctuations scheme.
  Real(Std), Parameter :: A2 = 1.0    !
  Real(Std), Parameter :: A3 = 1.0    !

  FieldReq => Reqs%FieldReqs(i)

  ! Set up index of horizontal/vertical grids for this field.
  ! Note that these spatial grids must be defined for any fluctuations calculation.
  iHGrid = Reqs%FieldReqs(i)%iHGrid
  iZGrid = Reqs%FieldReqs(i)%iZGrid
  ! Set up size of horizontal/vertical grids.
  If (iHGrid /= 0) Then
    nX = Grids%HGrids(iHGrid)%nX
    nY = Grids%HGrids(iHGrid)%nY
  Else
    Call Message('Error in DeriveSigmaCField: horizontal grid not specified', 3)
  End If
  If (iZGrid /= 0) Then
    nZ = Grids%ZGrids(iZGrid)%nZ
  Else
    Call Message('Error in DeriveSigmaCField: vertical grid not specified', 3)
  End If
  ! Set up index and size of travel-time grid for the relevant X-Stats field.
  iSGrid = Reqs%FieldReqs(iXStats)%iSGrid
  nT     = Grids%TGrids(iSGrid)%nT

  ! Loop over grid points.
  Do iX = 1, nX
  Do iY = 1, nY
  Do iZ = 1, nZ

    ! Location of grid point and values of other fields at this grid point.
    GridPoint  = (/ Grids%HGrids(iHGrid)%X(iX),  &
                    Grids%HGrids(iHGrid)%Y(iY),  &
                    Grids%ZGrids(iZGrid)%Z(iZ)   &
                 /)
    MeanConc   = Results%Fields(iMeanConc)%Std(1, iX, iY, iZ, 1, iTA, 1)
    TravelTime = Results%Fields(iTravelTime)%Std(1, iX, iY, iZ, 1, iTA, 1)

    ! SigmaC field is zero wherever mean concentration or travel-time is zero.
    If (MeanConc == 0.0 .or. TravelTime == 0.0) Cycle

    ! Determine X-Stats of the plume based on the mean travel time for grid box.
    ! $$ Currently using nearest travel-time on S-grid which is greater than (or equal)
    ! $$ to the grid box travel-time. Need to consider interpolation instead here.
    iS = 1
    DualTimes = .true.
    TravelTimeShort = RealTime2ShortTime(TravelTime)
    Do While ((TInTGrid(Grids%TGrids(iSGrid), iS) < TravelTimeShort) .and. iS < nT)
      iS = iS + 1
    End Do
    If (TInTGrid(Grids%TGrids(iSGrid), nT) < TravelTimeShort) Then
      Call Message('Warning: grid-box travel-time exceeds maximum travel-time available in X-Stats', 1)
    End If
    If ((iS == 1) .or. (TravelTimeShort >= TInTGrid(Grids%TGrids(iSGrid), nT))) Then
      ! If the travel-time is less than the shortest travel-time in the S-grid
      ! or greater than the largest travel-time in the S-grid, then use X-Stats
      ! at the respective extreme point (with no interpolation).
      DualTimes = .false.
    Else If (Results%Fields(iXStats)%Std(1, 1, 1, 1, iS, iTA, 1) == 0.0  .or. &
             Results%Fields(iXStats)%Std(1, 1, 1, 1, iS, iTA, 5) <= 0.1  .or. &
             Results%Fields(iXStats)%Std(1, 1, 1, 1, iS, iTA, 6) <= 0.1  .or. &
             Results%Fields(iXStats)%Std(1, 1, 1, 1, iS, iTA, 7) <= 0.1) Then
      ! $$ Fudge to avert problem of particles not reaching the greater travel-time
      ! $$ of the next S-grid point - use previous grid time instead here.
      ! $$ (more of an issue with short-interval releases than continuous ones).
      iS = iS - 1
      DualTimes = .false.
    End If
    If (DualTimes) Then
      S1 = ShortTime2RealTime(TInTGrid(Grids%TGrids(iSGrid), iS - 1))
      S2 = ShortTime2RealTime(TInTGrid(Grids%TGrids(iSGrid), iS))
      InterpCoeff = (TravelTime - S1) / (S2 - S1)
      XStats1 = Results%Fields(iXStats)%Std(1, 1, 1, 1, iS - 1, iTA, :)
      XStats2 = Results%Fields(iXStats)%Std(1, 1, 1, 1, iS,     iTA, :)
      ! Check that both X-Stats exist at these time/travel-time combinations.
      ! $$ Current fix-up allows for the case where S1 is zero (X-Stats default to zero)
      ! $$ OK for negligible source size but need to consider inital contribution otherwise.
      If ((XStats1(1) <= 0.0 .and. S1 > 0.0) .or. XStats2(1) <= 0.0) Then
        Call Message('Error in DeriveSigmaCField: X-Stats information not available', 3)
      End If
      ! Location of centre-of-mass of the plume (interpolated from X-Stats).
      PlumeCentre = (1.0 - InterpCoeff)*XStats1(2:4) + InterpCoeff*XStats2(2:4)
    Else
      XStats1 = Results%Fields(iXStats)%Std(1, 1, 1, 1, iS, iTA, :)
      ! Check that X-Stats exist at this time/travel-time combination.
      If (XStats1(1) <= 0.0) Then
        Call Message('Error in DeriveSigmaCField: X-Stats information not available', 3)
      End If
      ! Location of centre-of-mass of the plume (obtained from X-Stats).
      PlumeCentre = XStats1(2:4)
    End If

    ! Location of grid-point (in same coordinate systems as used by X-Stats).
    Call ResetFlowMemory(Flows, FlowMemory)
    Position = X2Position(                                                    &
                 Coords,                                                      &
                 GridPoint,                                                   &
                 Reqs%FieldReqs(i)%iHGridCoord, Reqs%FieldReqs(i)%iZGridCoord &
               )
    Call ConvertToH(Coords, Reqs%FieldReqs(iXStats)%iHCoord, Position)
    Call ConvertToZ(                                             &
           Coords, Grids, Domains,                               &
           Reqs%FieldReqs(iXStats)%iZCoord,                      &
           TInTGrid(Grids%TGrids(Reqs%FieldReqs(i)%iTGrid), iT), &
           .true., TravelTimeShort,                              &
           Position,                                             &
           Units, Mets, Flows,                                   &
           FlowMemory, ErrorCode                                 &
         )
    GridPoint = Position2X(                                                        &
                  Coords,                                                          &
                  Position,                                                        &
                  Reqs%FieldReqs(iXStats)%iHCoord, Reqs%FieldReqs(iXStats)%iZCoord &
                )

    ! Get flow information at the plume's centre-of-mass.
    Call ResetFlowMemory(Flows, FlowMemory)
    Position = X2Position(                                                        &
                 Coords,                                                          &
                 PlumeCentre,                                                     &
                 Reqs%FieldReqs(iXStats)%iHCoord, Reqs%FieldReqs(iXStats)%iZCoord &
               )
    Call GetAttrib(                                                              &
           A_Flow,                                                               &
           Coords, Grids, Domains,                                               &
           Moisture      = .false.,                                              &
           Inhomog       = .true.,                                               &
           Homog         = .false.,                                              & ! $$
           Time          = TInTGrid(Grids%TGrids(Reqs%FieldReqs(i)%iTGrid), iT), &
           AnyTravelTime = .true.,                                               &
           TravelTime    = TravelTimeShort,                                      &
           Position      = Position,                                             &
           Units         = Units,                                                &
           Mets          = Mets,                                                 &
           Flows         = Flows,                                                &
           FlowMemory    = FlowMemory,                                           &
           Flow          = FlowPlume,                                            &
           Cloud         = Cloud,                                                &
           Rain          = Rain,                                                 &
           Surface       = Surface,                                              &
           Soil          = Soil,                                                 &
           Plant         = Plant,                                                &
           ErrorCode     = ErrorCode                                             &
         )

         ! $$ need to deal with errors

    ! Process met data for the fluctuations scheme.
    HorzWind  = Sqrt(FlowPlume%U(1)**2 + FlowPlume%U(2)**2)
    SigmaVel2 = (FlowPlume%SigUU(1) + FlowPlume%SigUU(2) + FlowPlume%SigUU(3)) / 3.0
    Epsilon   = FlowPlume%Eps
    Z0        = FlowPlume%Z0
    ! Check that Epsilon is non-zero (otherwise a division by zero results).
    If (Epsilon == 0.0) Then
      Call Message('Error in DeriveSigmaCField: energy dissipation rate is zero', 3)
    End If

    ! Cross-wind and vertical offsets of grid-point relative to the centre-of-mass.
    YOffset = Abs( (GridPoint(1) - PlumeCentre(1)) * FlowPlume%U(2) -            &
                   (GridPoint(2) - PlumeCentre(2)) * FlowPlume%U(1) ) / HorzWind
    ZOffset = Abs(GridPoint(3) - PlumeCentre(3))

    ! Calculate length-scale Lc of the in-plume fluctuations.
    Lc = 1.0 / (A2 * Sqrt(Epsilon * TravelTime**3)) +       &
         1.0 / (A3 * SigmaVel2 * Sqrt(SigmaVel2) / Epsilon)
    Lc = 1.0 / Lc
    ! Calculate estimate of advection-diffusion equivalence height Zb.
    Zb = 0.16 * Lc / Log((0.16*Lc+Z0)/Z0)
    ! Get flow information at the effective height of the grid-point.
    GridPoint(3) = Max(GridPoint(3), Zb)
    Call ResetFlowMemory(Flows, FlowMemory)
    Position = X2Position(                                                        &
                 Coords,                                                          &
                 GridPoint,                                                       &
                 Reqs%FieldReqs(iXStats)%iHCoord, Reqs%FieldReqs(iXStats)%iZCoord &
               )
    Call GetAttrib(                                                              &
           A_Flow,                                                               &
           Coords, Grids, Domains,                                               &
           Moisture      = .false.,                                              &
           Inhomog       = .false.,                                              &
           Homog         = .false.,                                              & ! $$
           Time          = TInTGrid(Grids%TGrids(Reqs%FieldReqs(i)%iTGrid), iT), &
           AnyTravelTime = .true.,                                               &
           TravelTime    = TravelTimeShort,                                      &
           Position      = Position,                                             &
           Units         = Units,                                                &
           Mets          = Mets,                                                 &
           Flows         = Flows,                                                &
           FlowMemory    = FlowMemory,                                           &
           Flow          = FlowReceptor,                                         &
           Cloud         = Cloud,                                                &
           Rain          = Rain,                                                 &
           Surface       = Surface,                                              &
           Soil          = Soil,                                                 &
           Plant         = Plant,                                                &
           ErrorCode     = ErrorCode                                             &
         )

        ! $$ need to deal with errors

    EffHorzWindAtReceptor = Sqrt(FlowReceptor%U(1)**2 + FlowReceptor%U(2)**2)

    ! Calculate SigmaY and SigmaZ of the plume.
    ! $$ Currently using flow information (wind direction) at plume centre
    ! $$ - consider getting separate flow details at the two end points?
    If (DualTimes) Then
      Call CalculatePlumeSigmas(.true., InterpCoeff,                             &
                                XStats1(5), XStats1(8), XStats1(6), XStats1(7),  &
                                XStats2(5), XStats2(8), XStats2(6), XStats2(7),  &
                                FlowPlume%U(1), FlowPlume%U(2),                  &
                                FlowPlume%U(1), FlowPlume%U(2),                  &
                                SigmaY, SigmaZ)
    Else
      Call CalculatePlumeSigmas(.false., 0.0,                                    &
                                XStats1(5), XStats1(8), XStats1(6), XStats1(7),  &
                                0.0, 0.0, 0.0, 0.0,                              &
                                FlowPlume%U(1), FlowPlume%U(2), 0.0, 0.0,        &
                                SigmaY, SigmaZ)
    End If

    ProcAvTime = ZeroTime()
    Do j = 1, FieldReq%nProcesses ! $$ this needs reviewing (e.g. treat Int?)
      Process => Reqs%Processes(FieldReq%iProcesses(j))
      If (Process%TAvInt /= P_Av) Cycle
      If (Process%Code /= P_AvInt) Exit
      ProcAvTime = ProcAvTime + Process%T
    End Do

    ! Calculate an estimate of the concentration standard deviation at the grid point.
    Call CalculateSigmaC(CalcType        = 1,                                     &
                         TravelTime      = TravelTime,                            &
                         YOffset         = YOffset,                               &
                         ZOffset         = ZOffset,                               &
                         MeanConc        = MeanConc,                              &
                         SigmaY          = SigmaY,                                &
                         SigmaZ          = SigmaZ,                                &
                         U               = HorzWind,                              &
                         SigmaVel2       = SigmaVel2,                             &
                         Epsilon         = Epsilon,                               &
                         SourceDiameter  = 0.0,                                   & !$$
                         ReleaseRate     = Sources%Sources(1)%SourceStrengths(1), & !$$ Assumes (1),
                         SigmaPlumeRise2 = 0.0,                                   & !$$ (1) and source _Rate_
                         TAverage        = ProcAvTime > ZeroTime(),               &   ! $$ is this needed?
                         AveragingTime   = ShortTime2RealTime(Time2ShortTime(     &
                                             ProcAvTime                           &
                                           )),                                    & !$$ Needs to be set up.
                         ReleaseTime     = 0.0,                                   & !$$
                         FluctsAdvectVel = EffHorzWindAtReceptor,                 &
                         SigmaC          = Results%Fields(i)%Std(1, iX, iY, iZ, 1, iTA, 1))
  End Do
  End Do
  End Do

End Subroutine DeriveSigmaCField

!-------------------------------------------------------------------------------------------------------------

!## Added subroutine to process a mean travel time field (scaling by total mass
!## in each grid-box).
!$$ In the longer term, it might be more efficient to calculate mass fields rather
!$$ than conc fields - and to then derive concentrations by dividing by volume.
!$$ That is, the opposite way around to the current set up - but needs some work.

Subroutine DeriveMeanTravelTimeField(i, iMeanConc, iTA, Coords, Grids, Reqs, Results)
! Processes mean travel time field at any particular time in its temporal grid.

  Implicit None
  ! Argument list:
  Integer,        Intent(In)    :: i         ! Index of mean travel time field.
  Integer,        Intent(In)    :: iMeanConc ! Index of mean concentration field.
  Integer,        Intent(In)    :: iTA       ! Index for T in a field array.
  Type(Coords_),  Intent(In)    :: Coords    !
  Type(Grids_),   Intent(In)    :: Grids     !
  Type(Reqs_),    Intent(In)    :: Reqs      !
  Type(Results_), Intent(InOut) :: Results   !
  ! Locals:
  Integer   :: iHGrid       ! Index for horizontal grid.
  Integer   :: iZGrid       ! Index for vertical grid.
  Integer   :: nX           !} Size of horizontal grid.
  Integer   :: nY           !}
  Integer   :: nZ           ! Size of vertical grid.
  Integer   :: iX           ! Index for X.
  Integer   :: iY           ! Index for Y.
  Integer   :: iZ           ! Index for Z.
  Real(Std) :: MeanConc     ! Mean concentration at grid-point.
  Real(Std) :: GridPoint(3) ! Coords of grid-point.
  Real(Std) :: Volume       ! Volume of grid-box.
  Real(Std) :: HMax         ! Max metric coefficient.
  Real(Std) :: H1           !} Metric terms for horizontal distances
  Real(Std) :: H2           !} evaluated at the grid-point.

  ! Set up index of horizontal/vertical grids for this field.
  iHGrid = Reqs%FieldReqs(i)%iHGrid
  iZGrid = Reqs%FieldReqs(i)%iZGrid
  ! Set up size of horizontal/vertical grids.
  If (iHGrid /= 0) Then
    nX = Grids%HGrids(iHGrid)%nX
    nY = Grids%HGrids(iHGrid)%nY
  Else
    nX = 1
    nY = 1
  End If
  If (iZGrid /= 0) Then
    nZ = Grids%ZGrids(iZGrid)%nZ
  Else
    nZ = 1
  End If

  ! Loop over grid points.
  Do iX = 1, nX
  Do iY = 1, nY
  Do iZ = 1, nZ

    ! Mean concentration at grid-point.
    MeanConc = Results%Fields(iMeanConc)%Std(1, iX, iY, iZ, 1, iTA, 1)

    ! Set travel-time to zero if mean concentration (~ mass in grid box) is zero.
    If (MeanConc == 0.0) Then
      Results%Fields(i)%Std(1, iX, iY, iZ, 1, iTA, 1) = 0.0
      Cycle
    End If

    ! Check that grids are defined. $$ in principle it should be possible to set up
    ! routine for horz/vert integrated values, but extra work will be needed.
    If (iHGrid == 0 .or. iZGrid == 0) Then
      Call Message(                                                        &
             'Error in DeriveMeanTravelTimeField: horizontal/vertical ' // &
             'integrated values are currently not supported',              &
             3                                                             &
           )
    End If

    ! Location of grid point.
    GridPoint  = (/                              &
                    Grids%HGrids(iHGrid)%X(iX),  &
                    Grids%HGrids(iHGrid)%Y(iY),  &
                    Grids%ZGrids(iZGrid)%Z(iZ)   &
                 /)

    ! Calculate volume of grid-box for mass calculation.
    Volume = 1.0
    Call MetricCoeffs(Coords%HCoords(Grids%HGrids(iHGrid)%iHCoord),  &
                      GridPoint(1:2), HMax, H1, H2)
    Volume = Volume * Grids%HGrids(iHGrid)%dX * Grids%HGrids(iHGrid)%dY * H1 * H2
    Volume = Volume * Min(                                              &
                          Grids%ZGrids(iZGrid)%dZ,                      &
                          GridPoint(3) + 0.5 * Grids%ZGrids(iZGrid)%dZ  &
                         ) ! $$ assumes ground at z = 0.0

    ! Scale value in mean travel time array by mass in grid-box.
    Results%Fields(i)%Std(1, iX, iY, iZ, 1, iTA, 1) = Results%Fields(i)%Std(1, iX, iY, iZ, 1, iTA, 1) /  &
                                                      (MeanConc * Volume)

  End Do
  End Do
  End Do

End Subroutine DeriveMeanTravelTimeField

!-------------------------------------------------------------------------------------------------------------

!## Revised the subroutine for processing pdfs.
!## Pdfs are now processed sequentially in a single block (rather than by groups).
!## There is a single run through the processing (unlike the iterative approach
!## used for fields), so that all calculations must be achievable on the first run.
!## This routine must be called *after* the ProcessFields routine.
!$$ Could include fields and pdfs together in a single iterative processing scheme
!$$ if some fields were to require pdf information (e.g. single-percentile fields).
!$$ Consider for future!

Subroutine DerivePdfs(Time, Coords, Grids, Flows, Reqs, Results)
! Processes pdfs.

  Implicit None
  ! Argument list:
  Type(Time_),    Intent(In)    :: Time    ! Current synchronisation time.
  Type(Coords_),  Intent(In)    :: Coords  !
  Type(Grids_),   Intent(In)    :: Grids   !
  Type(Flows_),   Intent(In)    :: Flows   !
  Type(Reqs_),    Intent(In)    :: Reqs    !
  Type(Results_), Intent(InOut) :: Results !
  ! Locals:
  Integer :: i       ! Index for pdf.
  Integer :: iHGrid  ! Index for horizontal grid.
  Integer :: iZGrid  ! Index for vertical grid.
  Integer :: iTGrid  ! Index for time grid.
  Integer :: nX      !} Size of horizontal grid.
  Integer :: nY      !}
  Integer :: nZ      ! Size of vertical grid.
  Integer :: nT      ! Size of time grid.
  Integer :: nP      ! Size of pdf at each grid point.
  Integer :: iStart  ! Start index of processing times for pdf.
  Integer :: iT      ! Index for processing time T.
  Integer :: iX      ! Index for X.
  Integer :: iY      ! Index for Y.
  Integer :: iZ      ! Index for Z.
  Integer :: iTA
  Integer :: iMeanC  !} Indices for mean concentration field
  Integer :: iSigmaC !} and sigma_c field on required grid.

  Do i = 1, Reqs%nPdfs

    ! Determine whether the pdf needs processing.
    If (Reqs%PdfReqs(i)%AvEnsemble) Cycle

    ! Index and size of grids for this pdf.
    iHGrid = Reqs%PdfReqs(i)%iHGrid
    iZGrid = Reqs%PdfReqs(i)%iZGrid
    iTGrid = Reqs%PdfReqs(i)%iTGrid
    If (iHGrid == 0) Then
      nX = 1
      nY = 1
    Else
      nX = Grids%HGrids(iHGrid)%nX
      nY = Grids%HGrids(iHGrid)%nY
    End If
    If (iZGrid == 0) Then
      nZ = 1
    Else
      nZ = Grids%ZGrids(iZGrid)%nZ
    End If
    If (iTGrid == 0) Then     ! $$ in fact, TGrid must exist at present.
      nT = 1
    Else
      nT = Grids%TGrids(iTGrid)%nT
    End If
    ! Pdf size.
    nP = Results%Pdfs(i)%PdfSize

    ! Set up start index of processing times.
    ! Note that a starting index is only updated for pdfs requested as output.
    iStart = Results%PdfsNextOutput(Reqs%PdfReqs(i)%iOutputGroup)

    ! Locate fields required in the calculation of the pdf.
    iMeanC  = Reqs%PdfReqs(i)%iAirConc
    iSigmaC = Reqs%PdfReqs(i)%iSigmaC

    ! Loop over times in temporal grid starting at the next required processing time.
    Do iT = iStart, nT
      If (Results%PdfsProc(iT, i)) Cycle
      If (Time2ShortTime(Time) < TInTGrid(Grids%TGrids(iTGrid), iT)) Exit
      iTA = CalciTA(Results%Fields(iMeanC), iT)
      ! Check that pre-requisite fields are available for this time.
      ! Since all fields are processed before pdfs, there should be no problems.
      If (.not. Results%Fields(iMeanC) %LastProcThisCase >= iT      .or.  &
          .not. Results%Fields(iSigmaC)%LastProcThisCase >= iT) Then
        Call Message('Error in DerivePdfs: fields not available here', 3)
      End If
      ! Calculate pdfs for this time.
      If (Results%Pdfs(i)%FixedThresholds) Then
        Do iX = 1, nX
        Do iY = 1, nY
        Do iZ = 1, nZ
          If (Results%Fields(iMeanC)%Std(1, iX, iY, iZ, 1, iTA, 1)  == 0.0  .or. &
              Results%Fields(iSigmaC)%Std(1, iX, iY, iZ, 1, iTA, 1) == 0.0) Cycle ! $$
          Call CalculatePdf(                                                              &
                 PdfType         = Results%Pdfs(i)%PdfType,                               &
                 FixedThresholds = Results%Pdfs(i)%FixedThresholds,                       &
                 PdfSize         = Results%Pdfs(i)%PdfSize,                               &
                 PdfThresholds   = Results%Pdfs(i)%Thresholds(1:MaxPdfSize),              &
                 MeanC           = Results%Fields(iMeanC)%Std(1, iX, iY, iZ, 1, iTA, 1),  &
                 SigmaC          = Results%Fields(iSigmaC)%Std(1, iX, iY, iZ, 1, iTA, 1), &
                 PdfValues       = Results%Pdfs(i)%Data(iX, iY, iZ, iT, 1:nP)             &
               )
        End Do
        End Do
        End Do
      Else
        Do iX = 1, nX
        Do iY = 1, nY
        Do iZ = 1, nZ
          If (Results%Fields(iMeanC)%Std(1, iX, iY, iZ, 1, iTA, 1)  == 0.0  .or. &
              Results%Fields(iSigmaC)%Std(1, iX, iY, iZ, 1, iTA, 1) == 0.0) Cycle ! $$
          Call CalculatePdf(                                                              &
                 PdfType         = Results%Pdfs(i)%PdfType,                               &
                 FixedThresholds = Results%Pdfs(i)%FixedThresholds,                       &
                 PdfSize         = Results%Pdfs(i)%PdfSize,                               &
                 MeanC           = Results%Fields(iMeanC)%Std(1, iX, iY, iZ, 1, iTA, 1),  &
                 SigmaC          = Results%Fields(iSigmaC)%Std(1, iX, iY, iZ, 1, iTA, 1), &
                 PdfScale        = Results%Pdfs(i)%Scale(iX, iY, iZ, iT),                 &
                 PdfValues       = Results%Pdfs(i)%Data(iX, iY, iZ, iT, 1:nP)             &
               )
        End Do
        End Do
        End Do
      End If
      Results%PdfsProc(iT, i) = .true.
    End Do
  End Do

End Subroutine DerivePdfs

!-------------------------------------------------------------------------------------------------------------

Subroutine OutputFields(                                      &
             iCase, EndOfCase, EndOfCases,                    &
             RadioactiveDecay, Chemistry, StartTime, EndTime, &
             OutputOpts, SizeDists,                           &
             OpenMPOpts,                                      &
             Coords, Grids, Flows, Specieses, Sources, Reqs,  &
             MaterialUnits,                                   &
             Units, Results                                   &
           )
! Outputs fields.

! Before editing this routine think about:
! (1) Unstructured and named horizontal grids
! (2) FieldReq%iTGrid = 0 (and the two possible causes of this - 'underlying quantity doesn't depend on T' and
!     'fields are processed over all T')
! (3) FieldReq%TAcross = true (and the possibility of different time grids or none for field in the same
!     output group)
! (4) Output of fields where the underlying quantity doesn't depend on T
! (5) EndOfCase and output of fields which are processed over all T
! (6) EndOfCases and output processed over cases
! (7) Restartability
! (8) Definitions of components of Results.
!
! Note that DiskUnit and ScreenUnit determine whether a file/window needs to be opened (a zero value indicates
! it needs to be opened) while nLines determines whether, after opening a file, headers need to be written or
! whether the file needs to be positioned after the first nLines lines. Normally opening is followed by
! writing the headers - however this may not be the case after a restart. On some operating systems all screen
! output goes to a single terminal window which is associated with unit 6 (in which case the question of
! opening the window doesn't arise).

  Implicit None
  ! Argument list:
  Integer,              Intent(In)            :: iCase
  Logical,              Intent(In)            :: EndOfCase
  Logical,              Intent(In)            :: EndOfCases
  Logical,              Intent(In)            :: RadioactiveDecay
  Logical,              Intent(In)            :: Chemistry
  Type(ShortTime_),     Intent(In)            :: StartTime
  Type(ShortTime_),     Intent(In)            :: EndTime
  Type(OutputOpts_),    Intent(In)            :: OutputOpts
  Type(SizeDists_),     Intent(In)            :: SizeDists
  Type(OpenMPOpts_),    Intent(In)            :: OpenMPOpts
  Type(Coords_),        Intent(In)            :: Coords
  Type(Grids_),         Intent(In),    Target :: Grids
  Type(Flows_),         Intent(In)            :: Flows
  Type(Specieses_),     Intent(In)            :: Specieses
  Type(Sources_),       Intent(In)            :: Sources
  Type(Reqs_),          Intent(In),    Target :: Reqs
  Type(MaterialUnits_), Intent(In)            :: MaterialUnits
  Type(Units_),         Intent(InOut)         :: Units
  Type(Results_),       Intent(InOut), Target :: Results
  ! iCase            :: Number of case.
  ! EndOfCase        :: Indicates the end of the case.
  ! EndOfCases       :: Indicates the end of the cases.
  ! RadioactiveDecay :: Indicates that radioactive decay is modelled.
  ! Chemistry        :: Indicates that chemistry is modelled.
  ! StartTime        :} Start time and expected end time of the run.
  ! EndTime          :}
  ! OutputOpts       :: Output options.
  ! OpenMPOpts       :: OpenMP options.
  ! Coords           :: Collection of coord systems.
  ! Grids            :: Collection of grids.
  ! Flows            :: Collection of flow module instance states.
  ! Specieses        :: Collection of specieses.
  ! Sources          :: Collection of sources.
  ! Reqs             :: Collection of requirements.
  ! MaterialUnits    :: Collection of material units.
  ! Units            :: Collection of information on input/output unit numbers.
  ! Results          :: Collection of results.
  ! Locals:
  Integer                                :: nNum(2)
  Integer                                :: NumList(Reqs%MaxFieldReqs, 2)
  Integer                                :: LastNumOutput(Reqs%MaxFieldReqs)
  Integer                                :: LastGraphOutput(Reqs%MaxFieldReqs)
  Integer                                :: iFile
  Integer,                       Pointer :: nLines
  Integer,                       Pointer :: DiskUnit
  Integer,                       Pointer :: ScreenUnit
  Integer,                       Pointer :: GraphUnit
  Character(MaxFileNameLength)           :: File
  Character(MaxCharLength3)              :: Title
  Integer                                :: nPrelimCols
  Integer                                :: nFieldCols
  Integer                                :: PCW
  Integer                                :: FCW
  Character(MaxCharLength)               :: T
  Logical                                :: NonZeroLine
  Character(MaxOutputLineLength)         :: Line
  Integer                                :: CharPos
  Character(MaxCharLength)               :: Header
  Character(MaxCharLength)               :: Unit
  Type(TGrid_),                  Pointer :: TGrid
  Type(HGrid_),                  Pointer :: HGrid
  Type(ZGrid_),                  Pointer :: ZGrid
  Type(DGrid_),                  Pointer :: DGrid
  Type(FieldReq_),               Pointer :: FieldReq
  Type(Field_),                  Pointer :: Field
  Integer                                :: iTGrid
  Integer                                :: iHGrid
  Integer                                :: iZGrid
  Integer                                :: iDGrid
  Integer                                :: iField
  Integer                                :: iOutputGroup
  Integer                                :: nT
  Integer                                :: nX
  Integer                                :: nY
  Integer                                :: nZ
  Integer                                :: nD
  Integer                                :: iTU
  Integer                                :: iXU
  Integer                                :: iYU
  Integer                                :: iZU
  Integer                                :: iDU
  Integer                                :: iTL
  Integer                                :: iXL
  Integer                                :: iYL
  Integer                                :: iZL
  Integer                                :: iDL
  Integer                                :: iT
  Integer                                :: iX
  Integer                                :: iY
  Integer                                :: iZ
  Integer                                :: iD
  Integer                                :: jTU
  Integer                                :: jXU
  Integer                                :: jYU
  Integer                                :: jZU
  Integer                                :: jDU
  Integer                                :: jTL
  Integer                                :: jXL
  Integer                                :: jYL
  Integer                                :: jZL
  Integer                                :: jDL
  Integer                                :: jT
  Integer                                :: jX
  Integer                                :: jY
  Integer                                :: jZ
  Integer                                :: jD
  Integer                                :: kTU
  Integer                                :: kXU
  Integer                                :: kYU
  Integer                                :: kZU
  Integer                                :: kDU
  Integer                                :: kTL
  Integer                                :: kXL
  Integer                                :: kYL
  Integer                                :: kZL
  Integer                                :: kDL
  Integer                                :: kT
  Integer                                :: kX
  Integer                                :: kY
  Integer                                :: kZ
  Integer                                :: kD
  Integer                                :: iTA
  Integer                                :: iComponent
  Integer                                :: iRange
  Integer                                :: iDS
  Integer                                :: i
  Integer                                :: j
  Integer                                :: k
  Integer                                :: k1
  Integer                                :: k2
  Integer                                :: IOStat
  ! nNum            :: Number of fields from an output group to output numerically to disk and screen.
  ! NumList         :: Indices of the fields from an output group to output numerically to disk and screen.
  ! LastNumOutput   :: Index (for each field) of the last time to have been output numerically.
  ! LastGraphOutput :: Index (for each field) of the last time to have been output graphically.
  ! iFile           :: Index of file or window for numerical output, corresponding to the second index of the
  !                    arrays Results%FieldnLines, Results%FieldDiskUnits and Results%FieldScreenUnits.
  ! nLines          :} Abbreviations for number of lines written to file, disk unit and screen unit for
  ! DiskUnit        :} numerical output, and unit for graphical output.
  ! ScreenUnit      :}
  ! GraphUnit       :}
  ! File            :: File name for numerical output.
  ! Title           :: Window title.
  ! nPrelimCols     :} Number of preliminary and field columns (columns of comma separated data not columns of
  ! nFieldCols      :} characters) in numerical output.
  ! PCW             :: Preliminary column width (excluding separating commas) for numerical output.
  ! FCW             :: Field column width (excluding separating commas) for numerical output.
  ! T               :: Character version of the time corresponding to the data for numerical output.
  ! NonZeroLine     :: Indicates that some of the field values in Line are non-zero.
  ! Line            :: Character string for writing output.
  ! CharPos         :: Pointer to position in line where next entry is to be inserted
  ! Header          :: Character string for column header.
  ! Unit            :: Character string for the unit in which a quantity is expressed.
  ! TGrid           :} Abbreviations for grids, requirements and results.
  ! HGrid           :}
  ! ZGrid           :}
  ! DGrid           :}
  ! FieldReq        :}
  ! Field           :}
  ! iTGrid          :] Indices of grids.
  ! iHGrid          :]
  ! iZGrid          :]
  ! iDGrid          :]
  ! iField          :: Index of field requirement and/or field.
  ! iOutputGroup    :: Index of output group.
  ! nT              :} Size of grids.
  ! nX              :}
  ! nY              :}
  ! nZ              :}
  ! nD              :}
  ! iTU             :] Upper and lower indices and loop indices for looping over output files in connection
  ! iXU             :] with T/X/Y/Z for numerical output. iTU, iTL and iT are also used for looping over times
  ! iYU             :] for graphical output.
  ! iZU             :]
  ! iDU             :]
  ! iTL             :]
  ! iXL             :]
  ! iYL             :]
  ! iZL             :]
  ! iDL             :]
  ! iT              :]
  ! iX              :]
  ! iY              :]
  ! iZ              :]
  ! iD              :]
  ! jTU             :} Upper and lower indices and loop indices for looping over output rows in connection
  ! jXU             :} with T/X/Y/Z for numerical output.
  ! jYU             :}
  ! jZU             :}
  ! jDU             :}
  ! jTL             :}
  ! jXL             :}
  ! jYL             :}
  ! jZL             :}
  ! jDL             :}
  ! jT              :}
  ! jX              :}
  ! jY              :}
  ! jZ              :}
  ! jD              :}
  ! kTU             :] Upper and lower indices and loop indices for looping over output columns in connection
  ! kXU             :] with T/X/Y/Z for numerical output.
  ! kYU             :]
  ! kZU             :]
  ! kDU             :]
  ! kTL             :]
  ! kXL             :]
  ! kYL             :]
  ! kZL             :]
  ! kDL             :]
  ! kT              :]
  ! kX              :]
  ! kY              :]
  ! kZ              :]
  ! kD              :]
  ! iTA             :: Time index in a field array at which the results for a given time are stored.
  ! iDS             :: Loop index for looping over the two cases of disk numerical output and screen numerical
  !                    output.
  ! i               :} Loop indices.
  ! j               :}
  ! k               :: Loop index for looping over column header lines for numerical output.
  ! k1              :} Functions of k for controlling column header lines for numerical output. k1 and k2
  ! k2              :} control the preliminary and field columns respectively.
  ! IOStat          :: Error code for open and read statements.

  ! $$ Need to treat travel-time fields.

  !-----------------------------------------------!
  ! Initialise LastNumOutput and LastGraphOutput. !
  !-----------------------------------------------!

  Do iField = 1, Reqs%nFieldReqs
    LastNumOutput  (iField) = 0
    LastGraphOutput(iField) = 0
  End Do

  !-----------------------------------!
  ! Numerical disk and screen output. !
  !-----------------------------------!

  ! Loop over output groups.
  !$OMP PARALLEL DO                                 &
  !$OMP SCHEDULE(DYNAMIC,1)                         &
  !$OMP NUM_THREADS(OpenMPOpts%nOutputGroupThreads) &
  !$OMP DEFAULT(None)                               &
  !$OMP SHARED(iCase,                               &
  !$OMP   EndOfCase,                                &
  !$OMP   EndOfCases,                               &
  !$OMP   RadioactiveDecay,                         &
  !$OMP   Chemistry,                                &
  !$OMP   StartTime,                                &
  !$OMP   EndTime,                                  &
  !$OMP   OpenMPOpts,                               &
  !$OMP   OutputOpts,                               &
  !$OMP   Coords,                                   &
  !$OMP   Grids,                                    &
  !$OMP   Flows,                                    &
  !$OMP   Specieses,                                &
  !$OMP   SizeDists,                                &
  !$OMP   Sources,                                  &
  !$OMP   Reqs,                                     &
  !$OMP   MaterialUnits,                            &
  !$OMP   Units,                                    &
  !$OMP   Results,                                  &
  !$OMP   LastNumOutput,                            &
  !$OMP   LastGraphOutput)                          &
  !$OMP PRIVATE(nNum,                               &
  !$OMP   NumList,                                  &
  !$OMP   iFile,                                    &
  !$OMP   nLines,                                   &
  !$OMP   DiskUnit,                                 &
  !$OMP   ScreenUnit,                               &
  !$OMP   GraphUnit,                                &
  !$OMP   File,                                     &
  !$OMP   Title,                                    &
  !$OMP   nPrelimCols,                              &
  !$OMP   nFieldCols,                               &
  !$OMP   PCW,                                      &
  !$OMP   FCW,                                      &
  !$OMP   T,                                        &
  !$OMP   NonZeroLine,                              &
  !$OMP   Line,                                     &
  !$OMP   CharPos,                                  &
  !$OMP   Header,                                   &
  !$OMP   Unit,                                     &
  !$OMP   TGrid, HGrid, ZGrid, DGrid,               &
  !$OMP   FieldReq,                                 &
  !$OMP   Field,                                    &
  !$OMP   iTGrid, iHGrid, iZGrid, iDGrid,           &
  !$OMP   iField,                                   &
  !$OMP   iOutputGroup,                             &
  !$OMP   nT, nX, nY, nZ, nD,                       &
  !$OMP   iTU, iXU, iYU, iZU, iDU,                  &
  !$OMP   iTL, iXL, iYL, iZL, iDL,                  &
  !$OMP   iT, iX, iY, iZ, iD,                       &
  !$OMP   jTU, jXU, jYU, jZU, jDU,                  &
  !$OMP   jTL, jXL, jYL, jZL, jDL,                  &
  !$OMP   jT, jX, jY, jZ, jD,                       &
  !$OMP   kTU, kXU, kYU, kZU, kDU,                  &
  !$OMP   kTL, kXL, kYL, kZL, kDL,                  &
  !$OMP   kT, kX, kY, kZ, kD,                       &
  !$OMP   iTA,                                      &
  !$OMP   iRange,                                   &
  !$OMP   iComponent,                               &
  !$OMP   iDS,                                      &
  !$OMP   i, j, k, k1, k2,                          &
  !$OMP   IOStat                                    &
  !$OMP )
  OutputGroupLoop: Do iOutputGroup = 1, Reqs%nFieldGroups
    ! Determine which fields in the output group are to be output to disk and to screen.
    nNum(1) = 0
    nNum(2) = 0
    Do i = 1, Reqs%nFieldReqs
      If (Reqs%FieldReqs(i)%Ensemble .and. .not.EndOfCases) Cycle
      If (Reqs%FieldReqs(i)%iOutputGroup /= iOutputGroup) Cycle
      If (Reqs%FieldReqs(i)%Disk) Then
        nNum(1)             = nNum(1) + 1
        NumList(nNum(1), 1) = i
        iField              = i
      End If
      If (Reqs%FieldReqs(i)%Screen) Then
        nNum(2)             = nNum(2) + 1
        NumList(nNum(2), 2) = i
        iField              = i
      End If
    End Do
    If (nNum(1) == 0 .and. nNum(2) == 0) Cycle

    ! Abbreviations for a representative field requirement and field for the output group, for grid indices,
    ! for grids and for grid sizes. Grid sizes are set to 1 if the grid is missing.
    FieldReq => Reqs%FieldReqs(iField)
    Field    => Results%Fields(iField)
    iTGrid = FieldReq%iTGrid
    iHGrid = FieldReq%iHGrid
    iZGrid = FieldReq%iZGrid
    iDGrid = FieldReq%iDGrid
    If (iTGrid == 0) Then
      nT = 1
    Else
      TGrid => Grids%TGrids(iTGrid)
      nT = TGrid%nT
    End If
    If (iHGrid == 0) Then
      nX = 1
      nY = 1
    Else
      HGrid => Grids%HGrids(iHGrid)
      nX = HGrid%nX
      nY = HGrid%nY
    End If
    If (iZGrid == 0) Then
      nZ = 1
    Else
      ZGrid => Grids%ZGrids(iZGrid)
      nZ = ZGrid%nZ
    End If
    If (iDGrid == 0) Then
      nD = 1
    Else
      DGrid => Grids%DGrids(iDGrid)
      nD = nDInDGrid(DGrid)
    End If

    ! Column width.
    If (FieldReq%AlignCols) Then
      If (FieldReq%NameII .and. FieldReq%TAcross) Then
        PCW = 19
        FCW = 24
      Else if (FieldReq%NameII .and. .not. FieldReq%TAcross) Then
        PCW = 29
        FCW = 24
      Else
        PCW = OutputPrelimColumnWidth
        FCW = OutputFieldColumnWidth
      End If
    Else
      PCW = -1
      FCW = -1
    End If

    ! Determine range of T-grid indices for looping over output files/windows.
    !
    ! Note iTL and iTU are limited by what data is available for output, by what data has already been output,
    ! and by how the data are arranged between and within different files. Data are available for output if
    ! they have been processed.
    !
    ! For the case with TSeparateFile false and TAcross true, this limiting is hard to do. This is because (i)
    ! all the data to be output needs to be available before any can be output, and (ii) the time grids (and
    ! their presence/absence) may differ for the various fields in the output group. As a result we set iTU
    ! equal to iTL (this will let us enter the loop over output files/windows) if and only if EndOfCase is
    ! set, even though there may turn out to be no output. If it turns out that there is no output, we simply
    ! make sure we don't open files/windows or write any output (see 'If no output, cycle iDS Loop'). The code
    ! cannot be easily extended to output before EndOfCase is set since one would need to check that all the
    ! data to be output is available.
    If (FieldReq%TSeparateFile) Then
      iTL = Field%LastNumOutput + 1
      iTU = Field%LastProc
    Else If (.not. FieldReq%TAcross) Then
      iTL = Field%LastNumOutput + 1
      iTU = Field%LastProc
      iTU = Min(iTU, iTL)
    Else
      iTL = 1
      If (EndOfCase) Then
        iTU = 1
      Else
        iTU = 0
      End If
    End If

    ! Determine range of X-grid indices for looping over output files/windows.
    If (FieldReq%XSeparateFile) Then
      iXL = 1
      iXU = nX
    Else
      iXL = 1
      iXU = 1
    End If

    ! Determine range of Y-grid indices for looping over output files/windows.
    If (FieldReq%YSeparateFile) Then
      iYL = 1
      iYU = nY
    Else
      iYL = 1
      iYU = 1
    End If

    ! Determine range of Z-grid indices for looping over output files/windows.
    If (FieldReq%ZSeparateFile) Then
      iZL = 1
      iZU = nZ
    Else
      iZL = 1
      iZU = 1
    End If

    ! Determine range of D-grid indices for looping over output files/windows.
    If (FieldReq%DSeparateFile) Then
      iDL = 1
      iDU = nD
    Else
      iDL = 1
      iDU = 1
    End If

    ! Loop over output files/windows.
    iTLoop: Do iT = iTL, iTU
    iXLoop: Do iX = iXL, iXU
    iYLoop: Do iY = iYL, iYU
    iZLoop: Do iZ = iZL, iZU
    iDLoop: Do iD = iDL, iDU

      ! Index for storing number of lines written to file, disk unit and screen unit in Results. Note that,
      ! after exiting this routine, Results only needs to store details on files and windows that may need to
      ! be written to again. The index is chosen as small as possible consistent with the need to distinguish
      ! such files/windows after exiting this routine (there is no need to distinguish within the routine as
      ! we are inside the loop over files/windows).
      If (FieldReq%TSeparateFile .or. FieldReq%TAcross) Then
        iFile = 1
      Else
        iFile = iD + (iZ - 1)*iDU + (iY - 1)*iZU*iDU + (iX - 1)*iYU*iZU*iDU
      End If
      If (iFile > MaxFieldOutputFiles) Then
        If (nNum(1) > 0) Then
          Call Message(                                                             &
                 'FATAL ERROR: too many output files required at the same time ' // &
                 'for output group "'                                            // &
                 Trim(FieldReq%OutputGroup)                                      // &
                 '"',                                                               &
                 3                                                                  &
               )
        Else
          Call Message(                                              &
                 'FATAL ERROR: too many output windows required ' // &
                 'for output group "'                             // &
                 Trim(FieldReq%OutputGroup)                       // &
                 '"',                                                &
                 3                                                   &
               )
        End If
      End If

      ! Abbreviations for number of lines written to file, disk unit and screen unit.
      nLines     => Results%FieldnLines     (iOutputGroup, iFile)
      DiskUnit   => Results%FieldDiskUnits  (iOutputGroup, iFile)
      ScreenUnit => Results%FieldScreenUnits(iOutputGroup, iFile)

      ! Loop over disk output and screen output. In this loop we open files/windows, write headers and write
      ! column headers. Note the files/windows are opened and the headers written only after computing the
      ! first line of the column headers. This enables the number of field columns and the record length to be
      ! computed.
      iDSLoopHeaders: Do iDS = 1, 2

        ! Cycle if file/window already open (tested by DiskUnit /= 0 or ScreenUnit /= 0 as appropriate) or if
        ! the disk or screen case is not needed (tested by nNum(iDS) == 0).
        If (                                    &
          (iDS == 1 .and. DiskUnit   /= 0) .or. &
          (iDS == 2 .and. ScreenUnit /= 0) .or. &
          nNum(iDS) == 0                        &
        ) Cycle

        ! Loop over column headers.
        kLoop: Do k = 1, 21

          ! Options for order of column headers. Note the first pass (k = 1) must be included with k1 = 1 and
          ! with k2 /= 9, 22, 10, 11, 12 or 14 unless it is certain that the loop will not be 'cycled'. Other
          ! passes can be cycled to reduce the number of headers. k1 = 3 and k2 = 15 must only be used
          ! together.
          ! $$ Could arrange code to allow k1 = 2 on first pass and possibly k1 = 3 too.
          ! $$ Make better use of 'Header'. Currently only used for 2, 6-8.
          ! $$ could give location name when X, Y are down for unstructured grids.
          ! Preliminary column headers:
          ! k1 = 1: Commas.
          ! k1 = 2: T, X, Y, Z and D when these are 'down'.
          ! k1 = 3: Blank line.
          ! Field column headers: $$ reorder (see comments above too)
          ! k2 =  1: Species category.
          ! k2 = 21: Name.
          ! k2 =  2: Quantity.
          ! k2 =  3: Species.
          ! k2 =  4: Units.
          ! k2 =  5: Source/source group.
          ! k2 = 23: Particle size range.
          ! k2 = 20: Ensemble averaging information.
          ! k2 =  6: Time averaging/integrating information.
          ! k2 =  7: Horizontal averaging/integrating information.
          ! k2 =  8: Vertical averaging/integrating information.
          ! k2 = 16: Probabilities and percentiles.
          ! k2 = 19: Probabilities and percentiles - over ensemble.
          ! k2 = 17: Probabilities and percentiles - over time.
          ! k2 =  9: T when this is 'across' (otherwise line omitted).                 } If would otherwise be blank, these give
          ! k2 = 22: X-Y location name when this is 'across' (otherwise line omitted). } averaging/integrating information.
          ! k2 = 10: X when this is 'across' (otherwise line omitted).                 } $$ could change this convention when
          ! k2 = 11: Y when this is 'across' (otherwise line omitted).                 } other ideas (see "Plan for headers"
          ! k2 = 12: Z when this is 'across' (otherwise line omitted).                 } below) implemented.
          ! k2 = 18: D when this is 'across' (otherwise line omitted).
          ! k2 = 13: Commas.
          ! k2 = 14: Z plus Z averaging information when Z is 'across'.
          ! k2 = 15: Blank line.

          ! Convention for reading column headers (other than last row of column headers):
          ! For each row, should combine the (last) preliminary column header entry with the field column header entry.
          ! E.g. "  Quantity:,   AirConcentration," becomes "Quantity: AirConcentration".
          ! If field column header is blank, ignore preliminary column header.
          ! Ignore text in square brackets in preliminary column header.
          ! If preliminary column header begins with "...", concatenate with previous line.
          ! E.g. "    Prob & %ile info:,     Dth percentile," 
          !      "  ... over [ensemble],  ensemble of cases,"
          ! becomes
          ! "Prob & %ile info: Dth percentile over ensemble of cases".

          ! 1) NAME II fields format.
          If (FieldReq%NameII .and. FieldReq%TAcross) Then

            Select Case (k)
              Case (1)
                k1 = 1
                k2 = 1
              Case (2)
                k1 = 1
                k2 = 3
              Case (3)
                k1 = 1
                k2 = 6
              Case (4)
                k1 = 1
                k2 = 2
              Case (5)
                k1 = 1
                k2 = 4
              Case (6)
                k1 = 1
                k2 = 14
              Case (7)
                k1 = 2
                k2 = 9
              Case (8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21)
                Cycle
              Case (20)
                k1 = 3
                k2 = 15
            End Select

          ! 2) NAME II time series format.
          Else If (FieldReq%NameII .and. .not. FieldReq%TAcross) Then

            Select Case (k)
              Case (1)
                k1 = 1
                k2 = 11
              Case (2)
                k1 = 1
                k2 = 10
              Case (3)
                k1 = 1
                k2 = 22
              Case (4)
                k1 = 1
                k2 = 1
              Case (5)
                k1 = 1
                k2 = 3
              Case (6)
                k1 = 1
                k2 = 2
              Case (7)
                k1 = 1
                k2 = 14
              Case (8)
                k1 = 2
                k2 = 4
              Case (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21)
                Cycle
              Case (20)
                k1 = 3
                k2 = 15
            End Select

          ! 3) Standard format.
          Else

            Select Case (k)
              Case (21)
                k1 = 2
              Case Default
                k1 = 1
            End Select
            k2 = k
            Select Case (k) ! $$ reorder headers to avoid need for this.
              Case (1)
                k2 = 1
              Case (2)
                k2 = 21
              Case (3)
                k2 = 2
              Case (4)
                k2 = 3
              Case (5)
                k2 = 4
              Case (6)
                k2 = 5
              Case (7)
                k2 = 23
                If (OutputOpts%Pre65Format) Cycle
              Case (8)
                k2 = 20
              Case (9)
                k2 = 6
              Case (10)
                k2 = 7
              Case (11)
                k2 = 8
              Case (12)
                k2 = 16
              Case (13)
                k2 = 19
              Case (14)
                k2 = 17
              Case (15)
                k2 = 9
              Case (16)
                k2 = 22
              Case (17)
                k2 = 10
              Case (18)
                k2 = 11
              Case (19)
                k2 = 12
              Case (20)
                k2 = 18
              Case (21)
                k2 = 13
            End Select

          End If

          If (                                           &
            (k2 ==  9 .and. .not. FieldReq%TAcross) .or. &
            (k2 == 22 .and. .not. FieldReq%XAcross .and. .not. FieldReq%YAcross) .or. &
            (k2 == 10 .and. .not. FieldReq%XAcross) .or. &
            (k2 == 11 .and. .not. FieldReq%YAcross) .or. &
            (k2 == 12 .and. .not. FieldReq%ZAcross) .or. &
            (k2 == 14 .and. .not. FieldReq%ZAcross)      &
          ) Then
            If (k == 1) Call Message('UNEXPECTED FATAL ERROR in OutputFields', 4)
            Cycle
          End If

          If (k2 == 18 .and. .not. FieldReq%DAcross) Then
            If (OutputOpts%Pre65Format) Then
              k2 = 13
            Else
              If (k == 1) Call Message('UNEXPECTED FATAL ERROR in OutputFields', 4)
              Cycle
            End If
          End If

          nPrelimCols = 0
          nFieldCols  = 0

          Line = ' '

          ! Preliminary column headers.

          ! Commas.
          If (k1 == 1) Then

            If (FieldReq%GridIndices) Then
              If (.not. FieldReq%TAcross .and. iTGrid /= 0) Then
                Line = Trim(Line) // FormatChar(' ', PCW, .true., 'R')
                nPrelimCols = nPrelimCols + 1
              End If
              If (.not. FieldReq%XAcross .and. iHGrid /= 0) Then
                Line = Trim(Line) // FormatChar(' ', PCW, .true., 'R')
                nPrelimCols = nPrelimCols + 1
              End If
              If (.not. FieldReq%YAcross .and. iHGrid /= 0) Then
                If (.not. HGrid%Unstructured) Then
                  Line = Trim(Line) // FormatChar(' ', PCW, .true., 'R')
                  nPrelimCols = nPrelimCols + 1
                End If
              End If
              If (.not. FieldReq%ZAcross .and. iZGrid /= 0) Then
                Line = Trim(Line) // FormatChar(' ', PCW, .true., 'R')
                nPrelimCols = nPrelimCols + 1
              End If
              If (.not. FieldReq%DAcross .and. iDGrid /= 0) Then
                Line = Trim(Line) // FormatChar(' ', PCW, .true., 'R')
                nPrelimCols = nPrelimCols + 1
              End If
            End If

            If (.not. FieldReq%TAcross .and. iTGrid /= 0) Then
              Line = Trim(Line) // FormatChar(' ', PCW, .true., 'R')
              nPrelimCols = nPrelimCols + 1
            End If
            If (.not. FieldReq%XAcross .and. iHGrid /= 0) Then
              Line = Trim(Line) // FormatChar(' ', PCW, .true., 'R')
              nPrelimCols = nPrelimCols + 1
            End If
            If (.not. FieldReq%YAcross .and. iHGrid /= 0) Then
              If (.not. HGrid%Named) Then
                Line = Trim(Line) // FormatChar(' ', PCW, .true., 'R')
                nPrelimCols = nPrelimCols + 1
              End If
            End If
            If (.not. FieldReq%ZAcross .and. iZGrid /= 0) Then
              Line = Trim(Line) // FormatChar(' ', PCW, .true., 'R')
              nPrelimCols = nPrelimCols + 1
            End If
            If (.not. FieldReq%DAcross .and. iDGrid /= 0) Then
              Line = Trim(Line) // FormatChar(' ', PCW, .true., 'R')
              nPrelimCols = nPrelimCols + 1
            End If

            ! $$ Plan for headers (see "convention for reading headers" above):
            ! add "... from Z ="
            !     "... to Z =" (for av/int)
            !     "... over Z from Z ='
            !     "... to Z =' (for prob/%ile).
            If (.not. FieldReq%NameII .and. .not. OutputOpts%Pre65Format) Then
              Select Case (k2) 
                Case (1)
                  Line = Trim(Line) // FormatChar('Species category:',       PCW, .true., 'R')
                Case (21)
                  Line = Trim(Line) // FormatChar('Field Name:',             PCW, .true., 'R')
                Case (2)
                  Line = Trim(Line) // FormatChar('Quantity:',               PCW, .true., 'R')
                Case (3)
                  Line = Trim(Line) // FormatChar('Species:',                PCW, .true., 'R')
                Case (4)
                  Line = Trim(Line) // FormatChar('Units:',                  PCW, .true., 'R')
                Case (5)
                  Line = Trim(Line) // FormatChar('Source/source group:',    PCW, .true., 'R')
                Case (23)
                  Line = Trim(Line) // FormatChar('Particle size range:',    PCW, .true., 'R')
                Case (20)
                  Line = Trim(Line) // FormatChar('Ensemble av info:',       PCW, .true., 'R')
                Case (6)
                  Line = Trim(Line) // FormatChar('Time av/int info:',       PCW, .true., 'R')
                Case (7)
                  Line = Trim(Line) // FormatChar('Horizontal av/int info:', PCW, .true., 'R')
                Case (8)
                  Line = Trim(Line) // FormatChar('Vertical av/int info:',   PCW, .true., 'R')
                Case (16)
                  Line = Trim(Line) // FormatChar('Prob & %-ile info:',      PCW, .true., 'R')
                Case (19)
                  Line = Trim(Line) // FormatChar('... over [ensemble]',     PCW, .true., 'R')
                Case (17)
                  Line = Trim(Line) // FormatChar('... over [T]',            PCW, .true., 'R')
                Case (9)
                  Line = Trim(Line) // FormatChar('T:',                      PCW, .true., 'R')
                Case (22)
                  Line = Trim(Line) // FormatChar('X-Y location:',           PCW, .true., 'R')
                Case (10)
                  Line = Trim(Line) // FormatChar('X:',                      PCW, .true., 'R')
                Case (11)
                  Line = Trim(Line) // FormatChar('Y:',                      PCW, .true., 'R')
                Case (12)
                  Line = Trim(Line) // FormatChar('Z:',                      PCW, .true., 'R')
                Case (18)
                  Line = Trim(Line) // FormatChar('D:',                      PCW, .true., 'R')
                Case (14)
                  Line = Trim(Line) // FormatChar('Z plus Z av info:',       PCW, .true., 'R')
                Case Default
                  Line = Trim(Line) // FormatChar(' ',                       PCW, .true., 'R')
              End Select
              nPrelimCols = nPrelimCols + 1
            End If
            
          ! T, X, Y, Z and D when these are 'down'.
          Else If (k1 == 2) Then

            If (FieldReq%GridIndices) Then
              If (.not. FieldReq%TAcross .and. iTGrid /= 0) Then
                Line = Trim(Line) // FormatChar('T Index', PCW, .true., 'R')
              End If
              If (.not. FieldReq%XAcross .and. iHGrid /= 0) Then
                If (HGrid%Unstructured) Then
                  Line = Trim(Line) // FormatChar('Location Index', PCW, .true., 'R')
                Else
                  If (FieldReq%NameII) Then
                    Line = Trim(Line) // FormatChar('X grid', PCW, .true., 'R')
                  Else
                    Line = Trim(Line) // FormatChar('X Index', PCW, .true., 'R')
                  End If
                End If
              End If
              If (.not. FieldReq%YAcross .and. iHGrid /= 0) Then
                If (.not. HGrid%Unstructured) Then
                  If (FieldReq%NameII) Then
                    Line = Trim(Line) // FormatChar('Y grid', PCW, .true., 'R')
                  Else
                    Line = Trim(Line) // FormatChar('Y Index', PCW, .true., 'R')
                  End If
                End If
              End If
              If (.not. FieldReq%ZAcross .and. iZGrid /= 0) Then
                Line = Trim(Line) // FormatChar('Z Index', PCW, .true., 'R')
              End If
              If (.not. FieldReq%DAcross .and. iDGrid /= 0) Then
                Line = Trim(Line) // FormatChar('D Index', PCW, .true., 'R')
              End If
            End If

            If (.not. FieldReq%TAcross .and. iTGrid /= 0) Then
              Line = Trim(Line) // FormatChar('T', PCW, .true., 'R')
            End If
            If (.not. FieldReq%XAcross .and. iHGrid /= 0) Then
              If (HGrid%Named) Then
                Line = Trim(Line) // FormatChar('Location', PCW, .true., 'R')
              Else
                If (FieldReq%NameII .and. (Coords%HCoords(HGrid%iHCoord) .equiv. HCoord_LatLong())) Then
                  Line = Trim(Line) // FormatChar('Longitude', PCW, .true., 'R')
                Else
                  Line = Trim(Line) // FormatChar('X (' // Trim(HGrid%HCoordName) // ')', PCW, .true., 'R')
                End If
              End If
            End If
            If (.not. FieldReq%YAcross .and. iHGrid /= 0) Then
              If (.not. HGrid%Named) Then
                If (FieldReq%NameII .and. (Coords%HCoords(HGrid%iHCoord) .equiv. HCoord_LatLong())) Then
                  Line = Trim(Line) // FormatChar('Latitude', PCW, .true., 'R')
                Else
                  Line = Trim(Line) // FormatChar('Y (' // Trim(HGrid%HCoordName) // ')', PCW, .true., 'R')
                End If
              End If
            End If
            If (.not. FieldReq%ZAcross .and. iZGrid /= 0) Then
              Line = Trim(Line) // FormatChar('Z (' // Trim(ZGrid%ZCoordName) // ')', PCW, .true., 'R')
            End If
            If (.not. FieldReq%DAcross .and. iDGrid /= 0) Then
              Line = Trim(Line) // FormatChar('D', PCW, .true., 'R')
            End If

            If (.not. FieldReq%NameII .and. .not. OutputOpts%Pre65Format) Then
              Line = Trim(Line) // FormatChar(' ', PCW, .true., 'R')
            End If

          ! Blank line.
          Else If (k1 == 3) Then

          End If

          ! Field column headers.

          ! Loop over fields in the output group.
          iLoopHeaders: Do i = 1, nNum(iDS)

            ! Abbreviations for field requirement, field, grid indices, grids and grid sizes. Grid sizes are
            ! set to 1 if the grid is missing.
            FieldReq => Reqs%FieldReqs(NumList(i, iDS))
            Field    => Results%Fields(NumList(i, iDS))
            iTGrid = FieldReq%iTGrid
            iHGrid = FieldReq%iHGrid
            iZGrid = FieldReq%iZGrid
            iDGrid = FieldReq%iDGrid
            If (iTGrid == 0) Then
              nT = 1
            Else
              TGrid => Grids%TGrids(iTGrid)
              nT = TGrid%nT
            End If
            If (iHGrid == 0) Then
              nX = 1
              nY = 1
            Else
              HGrid => Grids%HGrids(iHGrid)
              nX = HGrid%nX
              nY = HGrid%nY
            End If
            If (iZGrid == 0) Then
              nZ = 1
            Else
              ZGrid => Grids%ZGrids(iZGrid)
              nZ = ZGrid%nZ
            End If
            If (iDGrid == 0) Then
              nD = 1
            Else
              DGrid => Grids%DGrids(iDGrid)
              nD = nDInDGrid(DGrid)
            End If

            ! Determine range of T/X/Y/Z/D-grid indices for row. Note that for the 'not across' case all that
            ! matters is that the lower and upper indices are the same. This bit of code could be simplified
            ! with 'Else' blocks but this has not been done to maintain similarity with code below.
            If (FieldReq%TSeparateFile .or. .not. FieldReq%TAcross) Then
              kTL = iT
              kTU = iT
            End If
            If (FieldReq%XSeparateFile .or. .not. FieldReq%XAcross) Then
              kXL = iX
              kXU = iX
            End If
            If (FieldReq%YSeparateFile .or. .not. FieldReq%YAcross) Then
              kYL = iY
              kYU = iY
            End If
            If (FieldReq%ZSeparateFile .or. .not. FieldReq%ZAcross) Then
              kZL = iZ
              kZU = iZ
            End If
            If (FieldReq%DSeparateFile .or. .not. FieldReq%DAcross) Then
              kDL = iD
              kDU = iD
            End If
            If (.not. FieldReq%TSeparateFile .and. FieldReq%TAcross) Then
              kTL = Field%LastNumOutput + 1
              kTU = Field%LastProc
            End If
            If (.not. FieldReq%XSeparateFile .and. FieldReq%XAcross) Then
              kXL = 1
              kXU = nX
            End If
            If (.not. FieldReq%YSeparateFile .and. FieldReq%YAcross) Then
              kYL = 1
              kYU = nY
            End If
            If (.not. FieldReq%ZSeparateFile .and. FieldReq%ZAcross) Then
              kZL = 1
              kZU = nZ
            End If
            If (.not. FieldReq%DSeparateFile .and. FieldReq%DAcross) Then
              kDL = 1
              kDU = nD
            End If

            ! Loop over output columns.
             jLoopHeaders: Do j  = 1,   FieldReq%nComponents * FieldReq%nRanges
            iComponent = CalciComponent(FieldReq, j)
            iRange     = CalciRange    (FieldReq, j)
            kTLoopHeaders: Do kT = kTL, kTU
            If (iTGrid /= 0) Then
              If (FieldReq%NameII .and. FieldReq%TAcross) Then
                T = Time2Char(                             &
                      ShortTime2Time(TInTGrid(TGrid, kT)), &
                      OutputOpts%Seconds,                  &
                      OutputOpts%DecimalPlaces,            &
                      .true.,                              &
                      .true.                               &
                    )
              Else If (FieldReq%NameII .and. .not. FieldReq%TAcross) Then
                T = Time2Char(                             &
                      ShortTime2Time(TInTGrid(TGrid, kT)), &
                      OutputOpts%Seconds,                  &
                      OutputOpts%DecimalPlaces,            &
                      .true.,                              &
                      .false.,                             &
                      .true.                               &
                    )
              Else
                T = Time2Char(                             &
                      ShortTime2Time(TInTGrid(TGrid, kT)), &
                      OutputOpts%Seconds,                  &
                      OutputOpts%DecimalPlaces,            &
                      .true.                               &
                    )
              End If
            End If
            kXLoopHeaders: Do kX = kXL, kXU
            kYLoopHeaders: Do kY = kYL, kYU
            kZLoopHeaders: Do kZ = kZL, kZU
            kDLoopHeaders: Do kD = kDL, kDU

              nFieldCols = nFieldCols + 1

              ! For floating D-grid, add extra column for D values.
              If (FieldReq%DGridFloating) Then
                If (k2 == 16 .or. k2 == 17 .or. k2 == 18) Then ! $$ add 19, or have only one k2 number?
                  Line = Trim(Line) // FormatChar('D value for next column', FCW, .true., 'R')
                Else
                  Line = Trim(Line) // FormatChar(' ', FCW, .true., 'R')
                End If
                nFieldCols = nFieldCols + 1
              End If

              ! Time averaging/integrating information. Also used for T when this is 'across' if this would
              ! otherwise be blank.
              Select Case (k2)

                Case (6, 9)

                  If (FieldReq%AvTAvInt == P_Av ) Header = 'average'
                  If (FieldReq%AvTAvInt == P_Int) Header = 'integral'
                  If (Scan(QInfo(FieldReq%iQuantity), 'T') == 0) Then
                    Header = FormatChar(' ', FCW, .true., 'R')
                  Else If (FieldReq%ComplexProc) Then
                    Header = FormatChar('Complex processing', FCW, .true., 'R')
                  Else If (FieldReq%AvTAvInt == 0) Then
                    Header = FormatChar('No time averaging', FCW, .true., 'R')
                  Else If (FieldReq%AvUseTGrid) Then
                    Header = FormatChar('Grid interval T-' // Trim(Header), FCW, .true., 'R')
                      ! $$ give range when T across?
                  Else If (.not. IsInfFuture(FieldReq%sAvT)) Then
                    If (FieldReq%NameII) Then
                      If (FieldReq%AvTAvInt == P_Av ) Header = 'averaged'
                      If (FieldReq%AvTAvInt == P_Int) Header = 'integrated'
                      Header =                                                                                 &
                        FormatChar(                                                                            &
                          Trim(                                                                                &
                            Int2Char(NInt(ShortTime2RealTime(FieldReq%sAvT)/3600.0), FormatString = 'I3.3')    &
                          )                                                                                 // &
                          ' hr time '                                                                       // &
                          Trim(Header),                                                                        &
                          FCW, .true., 'R'                                                                     &
                        )
                    Else
                      Header = FormatChar(                        &
                                 Trim(                            &
                                   Time2Char(FieldReq%AvT,        &
                                     OutputOpts%Seconds,          &
                                     OutputOpts%DecimalPlaces,    &
                                     .true.                       &
                                   )                              &
                                 )                             // &
                                 ' '                           // &
                                 Trim(Header),                    &
                                 FCW, .true., 'R'                 &
                               )
                    End If
                  Else If (FieldReq%iTGrid /= 0) Then
                    Header = FormatChar('Time ' // Trim(Header) // ' from start', FCW, .true., 'R')
                  Else
                    Header = FormatChar('Time ' // Trim(Header), FCW, .true., 'R')
                  End If

                ! Horizontal averaging/integrating information. Also used for X-Y location name when this is
                ! 'across', X when this is 'across', and Y when this is 'across', if these would otherwise be
                ! blank.
                Case (7, 22, 10, 11)

                  If (FieldReq%AvHAvInt == P_Av ) Header = 'average'
                  If (FieldReq%AvHAvInt == P_Int) Header = 'integral'
                  If (Scan(QInfo(FieldReq%iQuantity), 'H') == 0) Then
                    Header = FormatChar(' ', FCW, .true., 'R')
                  Else If (FieldReq%ComplexProc) Then
                    Header = FormatChar('Complex processing', FCW, .true., 'R')
                  Else If (FieldReq%AvHAvInt == 0) Then
                    Header = FormatChar('No horizontal averaging', FCW, .true., 'R')
                  Else If (FieldReq%AvUseHGrid) Then
                    Header = FormatChar('Grid box H-' // Trim(Header), FCW, .true., 'R')
                      ! $$ give range when H across?
                  Else If (FieldReq%AvX >= 0) Then
                    Header = FormatChar('Finite H-' // Trim(Header), FCW, .true., 'R') ! $$ give dX, dY?
                  Else
                    Header = FormatChar('Horizontal ' // Trim(Header), FCW, .true., 'R')
                  End If

                ! Vertical averaging/integrating information. Also used for Z when this is 'across' if this
                ! would otherwise be blank.
                Case (8, 12)

                  If (FieldReq%AvZAvInt == P_Av ) Header = 'average'
                  If (FieldReq%AvZAvInt == P_Int) Header = 'integral'
                  If (Scan(QInfo(FieldReq%iQuantity), 'Z') == 0) Then
                    Header = FormatChar(' ', FCW, .true., 'R')
                  Else If (FieldReq%ComplexProc) Then
                    Header = FormatChar('Complex processing', FCW, .true., 'R')
                  Else If (FieldReq%AvZAvInt == 0) Then
                    Header = FormatChar('No vertical averaging', FCW, .true., 'R')
                  Else If (FieldReq%AvUseZGrid) Then
                    Header = FormatChar('Grid box Z-' // Trim(Header), FCW, .true., 'R')
                      ! $$ give range when Z across?
                  Else If (FieldReq%AvBL) Then
                    Header = FormatChar('Boundary layer ' // Trim(Header), FCW, .true., 'R')
                  Else If (FieldReq%AvZ >= 0) Then
                    Header = FormatChar('Finite Z-' // Trim(Header), FCW, .true., 'R') ! $$ give dZ?
                  Else
                    Header = FormatChar('Vertical ' // Trim(Header), FCW, .true., 'R')
                  End If

              End Select

              Select Case (k2)

                ! Species category.
                Case (1)

                  If (FieldReq%iSpecies == 0) Then
                    Line = Trim(Line) // FormatChar(' ', FCW, .true., 'R')
                  Else
                    Line = Trim(Line)                                                                    // &
                           FormatChar(Specieses%Specieses(FieldReq%iSpecies)%Category, FCW, .true., 'R')
                  End If

                ! Name.
                Case (21)

                  Line = Trim(Line) // FormatChar(FieldReq%Name, FCW, .true., 'R')

                ! Quantity.
                Case (2)

                  ! Note that if dry deposition rate, wet deposition rate or deposition rate is time integrated,
                  ! 'rate' is omitted from the output.
                  If (FieldReq%NameII) Then
                    If (     FieldReq%iQuantity == Q_AirConc .and. FieldReq%AvTAvInt == P_Int) Then
                      Header = FormatChar('Dosage',           FCW, .true., 'R')
                    Else If (FieldReq%iQuantity == Q_DryDep  .and. FieldReq%AvTAvInt == P_Int) Then
                      Header = FormatChar('Dry deposition',   FCW, .true., 'R')
                    Else If (FieldReq%iQuantity == Q_WetDep  .and. FieldReq%AvTAvInt == P_Int) Then
                      Header = FormatChar('Wet deposition',   FCW, .true., 'R')
                    Else If (FieldReq%iQuantity == Q_Dep     .and. FieldReq%AvTAvInt == P_Int) Then
                      Header = FormatChar('Total deposition', FCW, .true., 'R')
                    Else
                      Header = FormatChar(FieldReq%Quantity,  FCW, .true., 'R', .true.)
                    End If
                    If (FieldReq%PCode == P_Percent .or. FieldReq%PCode == P_MaxMin) Then
                      ! $$ limit name II to 100%
                      Line = Trim(Line) // FormatChar('Max ' // Trim(AdjustL(Header)), FCW + 1, .false., 'R')
                    Else
                      Line = Trim(Line) // Header
                    End If
                  Else
                    If (     FieldReq%iQuantity == Q_DryDep .and. FieldReq%AvTAvInt == P_Int) Then
                      Header = FormatChar(FieldReq%Quantity(1:14), FCW, .true., 'R')
                    Else If (FieldReq%iQuantity == Q_WetDep .and. FieldReq%AvTAvInt == P_Int) Then
                      Header = FormatChar(FieldReq%Quantity(1:14), FCW, .true., 'R')
                    Else If (FieldReq%iQuantity == Q_Dep    .and. FieldReq%AvTAvInt == P_Int) Then
                      Header = FormatChar(FieldReq%Quantity(1:10), FCW, .true., 'R')
                    Else If (FieldReq%iQuantity == Q_AduEffCloudGammaDose) Then
                      Header = FormatChar('Ad Eff Cld Gamma Dose Rate', FCW, .true., 'R')
                      If (FieldReq%AvTAvInt == P_Int) Then
                        Header = FormatChar('Ad Eff Cld Gamma Dose', FCW, .true., 'R')
                      End If
                    Else If (FieldReq%iQuantity == Q_AduLunCloudGammaDose) Then
                      Header = FormatChar('Ad Lun Cld Gamma Dose Rate', FCW, .true., 'R')
                      If (FieldReq%AvTAvInt == P_Int) Then
                        Header = FormatChar('Ad Lun Cld Gamma Dose', FCW, .true., 'R')
                      End If
                    Else If (FieldReq%iQuantity == Q_AduThyCloudGammaDose) Then
                      Header = FormatChar('Ad Thy Cld Gamma Dose Rate', FCW, .true., 'R')
                      If (FieldReq%AvTAvInt == P_Int) Then
                        Header = FormatChar('Ad Thy Cld Gamma Dose', FCW, .true., 'R')
                      End If
                    Else If (FieldReq%iQuantity == Q_AduBoSCloudGammaDose) Then
                      Header = FormatChar('Ad BoS Cld Gamma Dose Rate', FCW, .true., 'R')
                      If (FieldReq%AvTAvInt == P_Int) Then
                        Header = FormatChar('Ad BoS Cld Gamma Dose', FCW, .true., 'R')
                      End If
                    Else
                      Header = FormatChar(FieldReq%Quantity, FCW, .true., 'R', .true.)
                    End If
                    Line = Trim(Line) // Header
                  End If

                ! Species.
                Case (3)

                  Line = Trim(Line) // FormatChar(FieldReq%SpeciesName, FCW, .true., 'R')

                ! Units.
                Case (4)

                  If (                                                 &
                    FieldReq%iQuantity == Q_AirConc               .or. &
                    FieldReq%iQuantity == Q_Mass                  .or. &
                    FieldReq%iQuantity == Q_PuffCentres           .or. &
                    FieldReq%iQuantity == Q_SigmaC                .or. &
                    FieldReq%iQuantity == Q_ChemistryField        .or. &
                    FieldReq%iQuantity == Q_EulerianConcentration .or. &
                    FieldReq%iQuantity == Q_Concentration              &
                  ) Then
                    Unit = MaterialUnits%MaterialUnits(FieldReq%iMaterialUnit)%Name
                    If (MaterialUnits%MaterialUnits(FieldReq%iMaterialUnit)%UnitType == DobsonUnitType) Then
                      If (FieldReq%NameII) Then
                        If (FieldReq%AvTAvInt == P_Int) Then
                          Unit = Trim(Unit) // 's'
                        End If
                        If (iHGrid /= 0 .and. (iZGrid == 0 .and. .not. FieldReq%AvBL)) Then
                          Unit = Trim(Unit)
                        Else 
                          Call Message('ERROR: Vertical grid has to be unspecified '           // & 
                                       'and boundary layer average can not be used if output ' // &
                                       'is requested in Dobson units.',                           &
                                      3                                                           &
                               )
                        End If
                      Else
                        If (FieldReq%AvTAvInt == P_Int) Then
                          Unit = Trim(Unit) // ' s'
                        End If
                        If (iHGrid /= 0 .and. (iZGrid == 0 .and. .not. FieldReq%AvBL)) Then
                          Unit = Trim(Unit)
                        Else 
                          Call Message('ERROR: Vertical grid has to be unspecified '           // & 
                                       'and boundary layer average can not be used if output ' // &
                                       'is requested in Dobson units.',                           &
                                      3                                                           &
                               )
                        End If
                      End If
                    Else
                      If (FieldReq%NameII) Then
                        If (FieldReq%AvTAvInt == P_Int) Then
                          Unit = Trim(Unit) // 's'
                        End If
                        If (iHGrid == 0 .and. (iZGrid /= 0 .or. FieldReq%AvBL)) Then
                          Unit = Trim(Unit) // '/m'
                        Else If (iHGrid /= 0 .and. (iZGrid == 0 .and. .not. FieldReq%AvBL)) Then
                          Unit = Trim(Unit) // '/m2'
                        Else If (iHGrid /= 0 .and. (iZGrid /= 0 .or. FieldReq%AvBL)) Then
                          Unit = Trim(Unit) // '/m3'
                        End If
                      Else
                        If (FieldReq%AvTAvInt == P_Int) Then
                          Unit = Trim(Unit) // ' s'
                        End If
                        If (iHGrid == 0 .and. (iZGrid /= 0 .or. FieldReq%AvBL)) Then
                          Unit = Trim(Unit) // ' / m'
                        Else If (iHGrid /= 0 .and. (iZGrid == 0 .and. .not. FieldReq%AvBL)) Then
                          Unit = Trim(Unit) // ' / m^2'
                        Else If (iHGrid /= 0 .and. (iZGrid /= 0 .or. FieldReq%AvBL)) Then
                          Unit = Trim(Unit) // ' / m^3'
                        End If
                      End If
                    End If
                    Line = Trim(Line) // FormatChar(Trim(Unit), FCW, .true., 'R')
                  Else If (FieldReq%iQuantity == Q_MixingRatio ) Then
                    Unit = MaterialUnits%MaterialUnits(FieldReq%iMaterialUnit)%Name
                    If (FieldReq%NameII) Then
                      If (FieldReq%AvTAvInt == P_Int) Then
                        Unit = Trim(Unit) // 's'
                      End If
                      If ( ( iHGrid /= 0 ) .and. ( iZGrid /= 0 ) ) Then
                        Unit = Trim(Unit)
                      Else
                        Call Message('ERROR: Both horizontal and vertical grid have to be ' // & 
                                     'specified in mixing ratio calculation.',                 &
                                     3                                                         &
                             )
                      End If
                    Else
                      If (FieldReq%AvTAvInt == P_Int) Then
                        Unit = Trim(Unit) // ' s'
                      End If
                      If (iHGrid == 0 .and. (iZGrid /= 0 .or. FieldReq%AvBL)) Then
                        Unit = Trim(Unit) // ' * m^2'
                      Else If (iHGrid /= 0 .and. (iZGrid == 0 .and. .not. FieldReq%AvBL)) Then
                        Unit = Trim(Unit) // ' * m'
                      Else If (iHGrid /= 0 .and. (iZGrid /= 0 .or. FieldReq%AvBL)) Then
                        Unit = Trim(Unit)
                      End If
                    End If
                    Line = Trim(Line) // FormatChar(Trim(Unit), FCW, .true., 'R')
                  Else If (FieldReq%iQuantity == Q_EMixingRatio ) Then
                    Unit = MaterialUnits%MaterialUnits(FieldReq%iMaterialUnit)%Name
                    If (FieldReq%NameII) Then
                      If (FieldReq%AvTAvInt == P_Int) Then
                        Unit = Trim(Unit) // 's'
                      End If
                      If ( ( iHGrid /= 0 ) .and. ( iZGrid /= 0 ) ) Then
                        Unit = Trim(Unit)
                      Else
                        Call Message('ERROR: Both horizontal and vertical grid have to be ' // & 
                                     'specified in E mixing ratio calculation.',               &
                                     3                                                         &
                             )
                      End If
                    Else
                      If (FieldReq%AvTAvInt == P_Int) Then
                        Unit = Trim(Unit) // ' s'
                      End If
                      If (iHGrid == 0 .and. (iZGrid /= 0 .or. FieldReq%AvBL)) Then
                        Unit = Trim(Unit) // ' * m^2'
                      Else If (iHGrid /= 0 .and. (iZGrid == 0 .and. .not. FieldReq%AvBL)) Then
                        Unit = Trim(Unit) // ' * m'
                      Else If (iHGrid /= 0 .and. (iZGrid /= 0 .or. FieldReq%AvBL)) Then
                        Unit = Trim(Unit)
                      End If
                    End If
                    Line = Trim(Line) // FormatChar(Trim(Unit), FCW, .true., 'R')
                  Else If (                             &
                    FieldReq%iQuantity == Q_DryDep .or. &
                    FieldReq%iQuantity == Q_WetDep .or. &
                    FieldReq%iQuantity == Q_Dep         &
                  ) Then
                    Unit = MaterialUnits%MaterialUnits(FieldReq%iMaterialUnit)%Name
                    If (FieldReq%NameII) Then
                      If (FieldReq%AvTAvInt == P_Int) Then
                        If (iHGrid /= 0) Then
                          Unit = Trim(Unit) // '/m2'
                        End If
                      Else
                        If (iHGrid == 0) Then
                          Unit = Trim(Unit) // '/s'
                        Else
                          Unit = Trim(Unit) // '/(m2 s)'
                        End If
                      End If
                    Else
                      If (FieldReq%AvTAvInt == P_Int) Then
                        If (iHGrid /= 0) Then
                          Unit = Trim(Unit) // ' / m^2'
                        End If
                      Else
                        If (iHGrid == 0) Then
                          Unit = Trim(Unit) // ' / s'
                        Else
                          Unit = Trim(Unit) // ' / (m^2 s)'
                        End If
                      End If
                    End If
                    Line = Trim(Line) // FormatChar(Trim(Unit), FCW, .true., 'R')
                  Else If (                        &
                    FieldReq%iQuantity == Q_X .or. &
                    FieldReq%iQuantity == Q_Y      &
                  ) Then
                    Line = Trim(Line) // FormatChar(Trim(FieldReq%HCoordName), FCW, .true., 'R')
                  Else If (                             &
                    FieldReq%iQuantity == Q_MeanZ  .or. &
                    FieldReq%iQuantity == Q_SigmaZ      &
                  ) Then
                    Line = Trim(Line) // FormatChar(Trim(FieldReq%ZCoordName), FCW, .true., 'R')
                  Else If (FieldReq%iQuantity == Q_XStats) Then
                    Line = Trim(Line) // FormatChar(' ', FCW, .true., 'R') ! $$
                  Else If (FieldReq%iQuantity == Q_AduEffCloudGammaDose) Then
                    Unit = 'Sv s^-1'
                    If (FieldReq%AvTAvInt == P_Int) Then
                      Unit = 'Sv'
                    End If
                    Line = Trim(Line) // FormatChar(Trim(Unit), FCW, .true., 'R')
                  Else If (                                           &
                    FieldReq%iQuantity == Q_AduLunCloudGammaDose .or. &
                    FieldReq%iQuantity == Q_AduThyCloudGammaDose .or. &
                    FieldReq%iQuantity == Q_AduBoSCloudGammaDose      &
                  ) Then
                    Unit = 'Gy s^-1'
                    If (FieldReq%AvTAvInt == P_Int) Then
                      Unit = 'Gy'
                    End If
                    Line = Trim(Line) // FormatChar(Trim(Unit), FCW, .true., 'R')
                Else If (                                             &
                  FieldReq%iQuantity == Q_OriginalSourceStrength .or. &
                  FieldReq%iQuantity == Q_RevisedSourceStrength       &
                ) Then
                  Unit = MaterialUnits%MaterialUnits(FieldReq%iMaterialUnit)%Name
                  If (FieldReq%AvTAvInt /= P_Int) Then
                    Unit = Trim(Unit) // '/s'
                  End If
                  Line = Trim(Line) // FormatChar(Trim(Unit), FCW, .true., 'R')
                Else
                  Line = Trim(Line)                          // &
                         FormatChar(                            &
                           Trim(QUnits(FieldReq%iQuantity)),    &
                           FCW, .true., 'R'                     &
                         )
                End If

                ! Source/source group.
                Case (5)

                  If (Scan(QInfo(FieldReq%iQuantity), 'r') == 0) Then
                    Line = Trim(Line) // FormatChar(' ', FCW, .true., 'R')
                  Else If (FieldReq%iSource == 0 .and. FieldReq%iSourceGroup == 0) Then
                    Line = Trim(Line) // FormatChar('All sources', FCW, .true., 'R')
                  Else If (FieldReq%iSource /= 0) Then
                    Line = Trim(Line) // FormatChar(Sources%Sources(FieldReq%iSource)%Name, FCW, .true., 'R')
                  Else If (FieldReq%iSourceGroup /= 0) Then
                    Line = Trim(Line)                                                                     // &
                           FormatChar(Sources%SourceGroups(FieldReq%iSourceGroup)%Name, FCW, .true., 'R')
                  End If

                ! Particle size range.
                Case (23)

                  If ((Scan(QInfo(FieldReq%iQuantity), 'd') == 0) .or. (FieldReq%iSizeDist == 0)) Then
                    Line = Trim(Line) // FormatChar(' ', FCW, .true., 'R')
                  Else
                    Line = Trim(Line)                                                                   // &
                           FormatChar(                                                                     &
                             Trim(Std2Char(                                                                &
                               SizeDists%SizeDists(FieldReq%iSizeDist)%DiameterRangeBoundary(iRange),      &
                               9, .false., 'L', 'F9.2'                                                     &
                             ))                                                                         // &
                             ' - '                                                                      // &
                             Trim(Std2Char(                                                                &
                               SizeDists%SizeDists(FieldReq%iSizeDist)%DiameterRangeBoundary(iRange+1),    &
                               9, .false., 'L', 'F9.2'                                                     &
                             ))                                                                         // &
                             ' um',                                                                        &
                             FCW, .true., 'R'                                                              &
                           )
                  End If

                ! Ensemble averaging information.
                Case (20)

                  If (FieldReq%ComplexProc) Then
                    Header = FormatChar('Complex processing',    FCW, .true., 'R')
                  Else If (.not. FieldReq%AvEnsemble) Then
                    Header = FormatChar('No ensemble averaging', FCW, .true., 'R')
                  Else
                    Header = FormatChar('Ensemble averaging',    FCW, .true., 'R')
                  End If

                  Line = Trim(Line) // Header

                ! Time averaging/integrating information.
                Case (6)

                  Line = Trim(Line) // Header

                ! Horizontal averaging/integrating information.
                Case (7)

                  Line = Trim(Line) // Header

                ! Vertical averaging/integrating information.
                Case (8)

                  Line = Trim(Line) // Header

                ! Probabilities and percentiles.
                Case (16)

                  If (FieldReq%PCode == P_Prob) Then
                    Header = FormatChar('Prob {value > D}', FCW, .true., 'R')
                  Else If (FieldReq%PCode == P_Percent .or. FieldReq%PCode == P_MaxMin) Then
                    Header = FormatChar('D-th percentile', FCW, .true., 'R')
                  Else
                    Header = FormatChar(' ', FCW, .true., 'R')
                  End If

                  Line = Trim(Line) // Header

                ! Probabilities and percentiles - over ensemble.
                Case (19)

                  If (FieldReq%PCode /= 0 .and. FieldReq%PEnsemble) Then
                    Header = FormatChar('ensemble of cases', FCW, .true., 'R')
                  Else
                    Header = FormatChar(' ', FCW, .true., 'R')
                  End If

                  Line = Trim(Line) // Header

                ! Probabilities and percentiles - over time.
                Case (17)

                  If (FieldReq%PCode /= 0 .and. FieldReq%PCode /= P_AvInt) Then
                    Header = FormatChar(                        &
                               Trim(                            &
                                 Time2Char(                     &
                                   FieldReq%PT,                 &
                                   OutputOpts%Seconds,          &
                                   OutputOpts%DecimalPlaces,    &
                                   .true.                       &
                                 )                              &
                               )                             // &
                               'interval',                      &
                               FCW, .true., 'R'                 &
                             )
                  Else
                    Header = FormatChar(' ', FCW, .true., 'R')
                  End If

                  Line = Trim(Line) // Header

                ! T when this is 'across'. If otherwise blank, this repeats time averaging/integrating
                ! information.
                Case (9)

                  If (iTGrid == 0) Then
                    Line = Trim(Line) // Header
                  Else
                    Line = Trim(Line) // FormatChar(T, FCW, .true., 'R')
                  End If

                ! X-Y location name when this is 'across'. If otherwise blank, this repeats horizontal
                ! averaging/integrating information.
                Case (22)

                  If (iHGrid == 0) Then
                    Line = Trim(Line) // Header
                  Else If (HGrid%Named .and. FieldReq%XAcross) Then
                    Line = Trim(Line) // FormatChar(HGrid%Names(kX), FCW, .true., 'R')
                  Else
                    Line = Trim(Line) // Header
                  End If

                ! X when this is 'across'. If otherwise blank, this repeats horizontal averaging/integrating
                ! information.
                Case (10)

                  If (iHGrid == 0) Then
                    Line = Trim(Line) // Header
                  Else
                    Line = Trim(Line) // FormatChar(                      &
                                           'X = '                      // &
                                           Trim(Std2Char(HGrid%X(kX))) // &
                                           ' '                         // &
                                           Trim(HGrid%HCoordName),        &
                                           FCW, .true., 'R'               &
                                         )
                  End If

                ! Y when this is 'across'. If otherwise blank, this repeats horizontal averaging/integrating
                ! information.
                Case (11)

                  If (iHGrid == 0) Then
                    Line = Trim(Line) // Header
                  Else If (HGrid%Unstructured) Then
                    Line = Trim(Line) // FormatChar(                      &
                                           'Y = '                      // &
                                           Trim(Std2Char(HGrid%Y(kX))) // &
                                           ' '                         // &
                                           Trim(HGrid%HCoordName),        &
                                           FCW, .true., 'R'               &
                                         )

                  Else
                    Line = Trim(Line) // FormatChar(                      &
                                           'Y = '                      // &
                                           Trim(Std2Char(HGrid%Y(kY))) // &
                                           ' '                         // &
                                           Trim(HGrid%HCoordName),        &
                                           FCW, .true., 'R'               &
                                         )
                  End If

                ! Z when this is 'across'. If otherwise blank, this repeats vertical averaging/integrating
                ! information.
                Case (12)

                  If (iZGrid == 0) Then
                    Line = Trim(Line) // Header
                  Else
                    Line = Trim(Line) // FormatChar(                      &
                                           'Z = '                      // &
                                           Trim(Std2Char(ZGrid%Z(kZ))) // &
                                           ' '                         // &
                                           Trim(ZGrid%ZCoordName),        &
                                           FCW, .true., 'R'               &
                                         )
                  End If

                ! D when this is 'across'.
                Case (18)

                  If (iDGrid == 0) Then
                    Header = FormatChar(' ', FCW, .true., 'R')
                  Else If (FieldReq%DGridFloating) Then
                    Header = FormatChar('D = variable', FCW, .true., 'R')
                  Else
                    Header = FormatChar('D = ' // Trim(Std2Char(DInDGrid(DGrid, kD, 0))), FCW, .true., 'R')
                  End If

                  Line = Trim(Line) // Header

                ! Commas.
                Case (13)

                  Line = Trim(Line) // FormatChar(' ', FCW, .true., 'R')

                ! Z plus Z averaging information when Z is 'across'.
                Case (14)

                  If (FieldReq%NameII) Then
                    If (Scan(QInfo(FieldReq%iQuantity), 'Z') == 0) Then
                      Line = Trim(Line) // FormatChar('Boundary layer', FCW, .true., 'R')
                    Else If (FieldReq%AvBL) Then
                      Line = Trim(Line) // FormatChar('Boundary layer', FCW, .true., 'R')
                    Else If (FieldReq%iZGrid == 0) Then
                      Line = Trim(Line) // FormatChar('Vertical integral', FCW, .true., 'R')
                    Else If (FieldReq%IntrinsicZAv) Then
                      If (ZGrid%ZCoordName .CIEq. 'FL') Then
                        Line = Trim(Line)                                                         // &
                               FormatChar(                                                           &
                                 'From FL'                                                        // &
                                 Trim(Int2Char(                                                      &
                                   NInt(ZGrid%Z(kZ) - 0.5*ZGrid%AvZ(kZ)), 3, .false., 'R', 'I3.3'    &
                                 ))                                                               // &
                                 ' - FL'                                                          // &
                                 Trim(Int2Char(                                                      &
                                   NInt(ZGrid%Z(kZ) + 0.5*ZGrid%AvZ(kZ)), 3, .false., 'R', 'I3.3'    &
                                 )),                                                                 &
                                 FCW, .true., 'R'                                                    &
                               )
                      Else
                        Line = Trim(Line)                                                       // &
                               FormatChar(                                                         &
                                 'From '                                                        // &
                                 Trim(Int2Char(                                                    &
                                   NInt(ZGrid%Z(kZ) - 0.5*ZGrid%AvZ(kZ)), 5, .false., 'R', 'I5'    &
                                 ))                                                             // &
                                 ' - '                                                          // &
                                 Trim(Int2Char(                                                    &
                                   NInt(ZGrid%Z(kZ) + 0.5*ZGrid%AvZ(kZ)), 5, .false., 'R', 'I5'    &
                                 ))                                                             // &
                                 Trim(ZGrid%ZCoordName),                                           &
                                 FCW, .true., 'R'                                                  &
                               )
                      End If
                    Else
                      If (ZGrid%ZCoordName .CIEq. 'FL') Then
                        Line = Trim(Line)                                                         // &
                               FormatChar(                                                           &
                                 'From FL'                                                        // &
                                 Trim(Int2Char(                                                      &
                                   NInt(ZGrid%Z(kZ) - 0.5*ZGrid%AvZ(kZ)), 3, .false., 'R', 'I3.3'    &
                                 ))                                                               // &
                                 ' - FL'                                                          // &
                                 Trim(Int2Char(                                                      &
                                   NInt(ZGrid%Z(kZ) + 0.5*ZGrid%AvZ(kZ)), 3, .false., 'R', 'I3.3'    &
                                 )),                                                                 &
                                 FCW, .true., 'R'                                                    &
                               )
                      Else
                        Line = Trim(Line)                                                       // &
                               FormatChar(                                                         &
                                 'From '                                                        // &
                                 Trim(Int2Char(                                                    &
                                   NInt(ZGrid%Z(kZ) - 0.5*ZGrid%AvZ(kZ)), 5, .false., 'R', 'I5'    &
                                 ))                                                             // &
                                 ' - '                                                          // &
                                 Trim(Int2Char(                                                    &
                                   NInt(ZGrid%Z(kZ) + 0.5*ZGrid%AvZ(kZ)), 5, .false., 'R', 'I5'    &
                                 ))                                                             // &
                                 Trim(ZGrid%ZCoordName),                                           &
                                 FCW, .true., 'R'                                                  &
                               )
                      End If
                    End If
                  Else
                    If (Scan(QInfo(FieldReq%iQuantity), 'Z') == 0) Then
                      Line = Trim(Line) // FormatChar(' ', FCW, .true., 'R')
                    Else If (FieldReq%AvBL) Then
                      Line = Trim(Line) // FormatChar('Boundary layer average', FCW, .true., 'R')
                    Else If (FieldReq%iZGrid == 0) Then
                      Line = Trim(Line) // FormatChar('Vertical integral', FCW, .true., 'R')
                    Else If (FieldReq%IntrinsicZAv) Then
                      Line = Trim(Line)                                        // &
                             FormatChar(                                          &
                               'From '                                         // &
                               Trim(Std2Char(ZGrid%Z(kZ) - 0.5*ZGrid%AvZ(kZ))) // &
                               ' - '                                           // &
                               Trim(Std2Char(ZGrid%Z(kZ) + 0.5*ZGrid%AvZ(kZ))) // &
                               Trim(ZGrid%ZCoordName),                            &
                               FCW, .true., 'R'                                   &
                             )
                    Else
                      Line = Trim(Line) // FormatChar(                      &
                                             'Z = '                      // &
                                             Trim(Std2Char(ZGrid%Z(kZ))) // &
                                             ' '                         // &
                                             Trim(ZGrid%ZCoordName),        &
                                             FCW, .true., 'R'               &
                                           )
                    End If
                  End If

                ! Blank line.
                Case (15)

              End Select

            End Do kDLoopHeaders
            End Do kZLoopHeaders
            End Do kYLoopHeaders
            End Do kXLoopHeaders
            End Do kTLoopHeaders
            End Do  jLoopHeaders

          End Do iLoopHeaders

          ! Check for no output, check line length, open files/windows and write headers.
          If (k == 1) Then

            ! Disk output.
            If (iDS == 1) Then

              ! If no output, cycle iDS Loop.
              ! $$ Use of the "Cycle" statement fails at three points in this subroutine when using
              !    the Cray compiler - this has been reported as a compiler bug. As an interim fix,
              !    the "Cycle" statements have been replaced with "Go To" statements. This change
              !    could be reverted once the compiler supports the construct.
              
              !If (nFieldCols == 0) Cycle iDSLoopHeaders
              If (nFieldCols == 0) Go To 100

              ! Check line length.
              ! $$ note doesn't work for PCW & FCW = -1 and, for PCW FCW /= -1, doesn't
              ! check PCW & FCW adequate.
              ! Could check PCW & FCW by looking for *** - check only on headers as can
              ! ensure PCW, FCW > output length of real or integer (but think about
              ! real*8).
              ! Could check line length by increasing Line by MaxCharLength and
              ! checking trimmed length doesn't exceed maxoutputlinelength.
              ! Alternatively could try to estimate required col width and line length
              ! at Init/SetUp stage
              ! Same applies to screen below.
              If (nPrelimCols*(PCW + 1) + nFieldCols*(FCW + 1) > MaxOutputLineLength) Then
                Call Message(                                              &
                       'FATAL ERROR: too long an output line required ' // &
                       'for output group "'                             // &
                       Trim(FieldReq%OutputGroup)                       // &
                       '"',                                                &
                       3                                                   &
                     )
              End If

              ! File name.
              If (FieldReq%NameII) Then
                File = Trim(OutputOpts%Folder) // Trim(FieldReq%OutputGroup)
              Else
                If (FieldReq%Ensemble) Then
                  File = Trim(OutputOpts%Folder) // Trim(FieldReq%OutputGroup) // '_CC'
                Else
                  File = Trim(OutputOpts%Folder) // Trim(FieldReq%OutputGroup) // '_C' // Int2Char(iCase)
                End If
              End If
              If (FieldReq%TSeparateFile .and. iTGrid /= 0) Then
                If (FieldReq%NameII) Then
                  T = FileNameTime(ChangeTimeZone(ShortTime2Time(TInTGrid(TGrid, iT)), TGrid%T0))
                  File = Trim(File) // '_' // T
                Else
                  T = FileNameTime(ChangeTimeZone(ShortTime2Time(TInTGrid(TGrid, iT)), TGrid%T0))
                  File = Trim(File) // '_T' // Trim(Int2Char(iT)) // '_' // Trim(T)
                End If
              End If
              If (FieldReq%XSeparateFile .and. iHGrid /= 0) Then
                If (HGrid%Unstructured) Then
                  If (HGrid%Named) Then
                    File = Trim(File) // '_' // Trim(HGrid%Names(iX))
                  Else
                    File = Trim(File) // '_L' // Int2Char(iX)
                  End If
                Else
                  File = Trim(File) // '_X' // Int2Char(iX)
                End If
              End If
              If (FieldReq%YSeparateFile .and. iHGrid /= 0) Then
                If (.not. HGrid%Unstructured) Then
                  File = Trim(File) // '_Y' // Int2Char(iY)
                End If
              End If
              If (FieldReq%ZSeparateFile .and. iZGrid /= 0) Then
                File = Trim(File) // '_Z' // Int2Char(iZ)
              End If
              If (FieldReq%DSeparateFile .and. iDGrid /= 0) Then
                File = Trim(File) // '_D' // Int2Char(iD)
              End If
              File = Trim(File) // '.txt'
              File = ConvertFileName(File)

              ! If nLines = 0 then either (i) the file hasn't yet been opened or written to, or (ii) the run
              ! has been restarted and any previously written lines are to be ignored.
              If (nLines == 0) Then

                ! Open file.
                DiskUnit = OpenFile(                                &
                             File            = File,                &
                             Units           = Units,               &
                             Status          = 'Replace',           &
                             RecL            = MaxOutputLineLength, &
                             FileDescription = 'output file'        &
                           )

                ! Write headers.
                Call FieldNumHeaders(                                   &
                       nPrelimCols, nFieldCols, DiskUnit,               &
                       RadioactiveDecay, Chemistry, StartTime, EndTime, &
                       OutputOpts,                                      &
                       Coords, Grids, Flows, Sources,                   &
                       FieldReq,                                        &
                       Results,                                         &
                       nLines                                           &
                     )

              ! If nLines /= 0 then the run has been restarted and the first nLines previously written lines
              ! are to be preserved and the rest are to be ignored.
              Else

                ! Open file.
                DiskUnit = OpenFile(                                &
                             File            = File,                &
                             Units           = Units,               &
                             Status          = 'Old',               &
                             RecL            = MaxOutputLineLength, &
                             FileDescription = 'output file'        &
                           )

                ! Read first nLines lines.
                Do i = 1, nLines
                  Read (DiskUnit, *, IOStat = IOStat)
                  If (IOStat > 0) Then
                    Call Message(                                                 &
                           'FATAL ERROR: An error occurred reading the file "' // &
                           Trim(File)                                          // &
                           '" following restart',                                 &
                           3                                                      &
                         )
                  Else If (IOStat < 0) Then
                    Call Message(                                             &
                           'FATAL ERROR: The file "'                       // &
                           Trim(File)                                      // &
                           '" is shorter than expected following restart',    &
                           3                                                  &
                         )
                  End If
                End Do

                ! Avoid writing column headers.
                !Cycle iDSLoopHeaders
                Go To 100

              End If

            ! Screen output.
            Else If (iDS == 2) Then

              ! If no output, cycle iDS Loop.
              !If (nFieldCols == 0) Cycle iDSLoopHeaders
              If (nFieldCols == 0) Go To 100

              ! Check line length.
              If (nPrelimCols*(PCW + 1) + nFieldCols*(FCW + 1) > MaxOutputLineLength) Then
                Call Message(                                              &
                       'FATAL ERROR: too long an output line required ' // &
                       'for output group "'                             // &
                       Trim(FieldReq%OutputGroup)                       // &
                       '"',                                                &
                       3                                                   &
                     )
              End If

              ! Window title.
              If (FieldReq%Ensemble) Then
                Title = Trim(FieldReq%OutputGroup) // ': Cases'
              Else
                Title = Trim(FieldReq%OutputGroup) // ': Case ' // Int2Char(iCase)
              End If
              If (FieldReq%TSeparateFile .and. iTGrid /= 0) Then
                Title = Trim(Title) // ': T ' // Int2Char(iT)
              End If
              If (FieldReq%XSeparateFile .and. iHGrid /= 0) Then
                If (HGrid%Unstructured) Then
                  Title = Trim(Title) // ': L ' // Int2Char(iX)
                Else
                  Title = Trim(Title) // ': X ' // Int2Char(iX)
                End If
              End If
              If (FieldReq%YSeparateFile .and. iHGrid /= 0) Then
                If (.not. HGrid%Unstructured) Then
                  Title = Trim(Title) // ': Y ' // Int2Char(iY)
                End If
              End If
              If (FieldReq%ZSeparateFile .and. iZGrid /= 0) Then
                Title = Trim(Title) // ': Z ' // Int2Char(iZ)
              End If
              If (FieldReq%DSeparateFile .and. iDGrid /= 0) Then
                Title = Trim(Title) // ': D ' // Int2Char(iD)
              End If

              ! Open window. In connection with nCols, note that readability is improved by making the window
              ! one column wider than the longest record.
#             ifdef CompaqPCCompiler
                Call GetNewUnit(ScreenUnit, Units)
                Call OpenTextWindow(                   &
                       Unit  = ScreenUnit,             &
                       Title = Trim(Title),            &
                       nRows = 200,                    &
                       nCols = MaxOutputLineLength + 1 &
                     )
#             endif
#             ifdef IntelLinCompiler
                ScreenUnit = 6
#             endif
#             ifdef CrayCLECompiler
                ScreenUnit = 6
#             endif
#             ifdef sun
                ScreenUnit = 6
#             endif

            End If

          End If

          ! Write column headers.
          If (FieldReq%NameII) Then
            If (iDS == 1) Then
              Write (DiskUnit, '(A)') Line(1:Max(Len_Trim(Line),1))
              nLines = nLines + 1
            Else
              Write (ScreenUnit, '(A)') Line(1:Max(Len_Trim(Line),1))
            End If
          Else
            If (iDS == 1) Then
              Write (DiskUnit, '(A)') Trim(Line)
              nLines = nLines + 1
            Else
              Write (ScreenUnit, '(A)') Trim(Line)
            End If
          End If

        End Do kLoop

        100 Continue

      End Do iDSLoopHeaders

      ! Determine range of T-grid indices for looping over output rows.
      If (FieldReq%TSeparateFile) Then
        jTL = iT
        jTU = iT
      Else If (.not. FieldReq%TAcross) Then
        jTL = Field%LastNumOutput + 1
        jTU = Field%LastProc
      Else
        jTL = 1
        jTU = 1
      End If

      ! Determine range of X-grid indices for looping over output rows.
      If (FieldReq%XSeparateFile) Then
        jXL = iX
        jXU = iX
      Else If (.not. FieldReq%XAcross) Then
        jXL = 1
        jXU = nX
      Else
        jXL = 1
        jXU = 1
      End If

      ! Determine range of Y-grid indices for looping over output rows.
      If (FieldReq%YSeparateFile) Then
        jYL = iY
        jYU = iY
      Else If (.not. FieldReq%YAcross) Then
        jYL = 1
        jYU = nY
      Else
        jYL = 1
        jYU = 1
      End If

      ! Determine range of Z-grid indices for looping over output rows.
      If (FieldReq%ZSeparateFile) Then
        jZL = iZ
        jZU = iZ
      Else If (.not. FieldReq%ZAcross) Then
        jZL = 1
        jZU = nZ
      Else
        jZL = 1
        jZU = 1
      End If

      ! Determine range of D-grid indices for looping over output rows.
      If (FieldReq%DSeparateFile) Then
        jDL = iD
        jDU = iD
      Else If (.not. FieldReq%DAcross) Then
        jDL = 1
        jDU = nD
      Else
        jDL = 1
        jDU = 1
      End If

      ! Loop over output rows.
      jTLoop: Do jT = jTL, jTU
      If (iTGrid /= 0) Then
        If (FieldReq%NameII .and. FieldReq%TAcross) Then
          T = Time2Char(                             &
                ShortTime2Time(TInTGrid(TGrid, jT)), &
                OutputOpts%Seconds,                  &
                OutputOpts%DecimalPlaces,            &
                .true.,                              &
                .true.                               &
              )
        Else If (FieldReq%NameII .and. .not. FieldReq%TAcross) Then
          T = Time2Char(                             &
                ShortTime2Time(TInTGrid(TGrid, jT)), &
                OutputOpts%Seconds,                  &
                OutputOpts%DecimalPlaces,            &
                .true.,                              &
                .false.,                             &
                .true.                               &
              )
        Else
          T = Time2Char(                             &
                ShortTime2Time(TInTGrid(TGrid, jT)), &
                OutputOpts%Seconds,                  &
                OutputOpts%DecimalPlaces,            &
                .true.                               &
              )
        End If
      End If
      jXLoop: Do jX = jXL, jXU
      jYLoop: Do jY = jYL, jYU
      jZLoop: Do jZ = jZL, jZU
      jDLoop: Do jD = jDL, jDU

        ! Determine range of T/X/Y/Z/D-grid indices for row except if these vary between fields in the output
        ! group.
        If (FieldReq%TSeparateFile .or. .not. FieldReq%TAcross) Then
          kTL = jT
          kTU = jT
        End If
        If (FieldReq%XSeparateFile .or. .not. FieldReq%XAcross) Then
          kXL = jX
          kXU = jX
        End If
        If (FieldReq%YSeparateFile .or. .not. FieldReq%YAcross) Then
          kYL = jY
          kYU = jY
        End If
        If (FieldReq%ZSeparateFile .or. .not. FieldReq%ZAcross) Then
          kZL = jZ
          kZU = jZ
        End If
        If (FieldReq%DSeparateFile .or. .not. FieldReq%DAcross) Then
          kDL = jD
          kDU = jD
        End If

        ! Loop over disk output and screen output.
        iDSLoop: Do iDS = 1, 2

          If (nNum(iDS) == 0) Cycle

          ! If no output, cycle iDS Loop.
          If (iDS == 1 .and. DiskUnit   == 0) Cycle
          If (iDS == 2 .and. ScreenUnit == 0) Cycle

          Line = REPEAT(' ', MaxOutputLineLength)
          CharPos = 1

          ! Start of line.

          ! Grid indices.
          If (FieldReq%GridIndices) Then
            If (.not. FieldReq%TAcross .and. iTGrid /= 0) Then
              Line(CharPos:CharPos+PCW) = Int2Char(jT, PCW, .true., 'R')
              CharPos = CharPos + PCW + 1
            End If
            If (.not. FieldReq%XAcross .and. iHGrid /= 0) Then
              If (FieldReq%NameII) Then
                Line(CharPos:CharPos+PCW) = Std2Char(Real(jX, Std) + 0.5, PCW, .true., 'R', 'F13.6')
                CharPos = CharPos + PCW + 1
              Else
                Line(CharPos:CharPos+PCW) = Int2Char(jX, PCW, .true., 'R')
                CharPos = CharPos + PCW + 1
              End If
            End If
            If (.not. FieldReq%YAcross .and. iHGrid /= 0) Then
              If (.not. HGrid%Unstructured) Then
                If (FieldReq%NameII) Then
                  Line(CharPos:CharPos+PCW) = Std2Char(Real(jY, Std) + 0.5, PCW, .true., 'R', 'F13.6')
                  CharPos = CharPos + PCW + 1
                Else
                  Line(CharPos:CharPos+PCW) = Int2Char(jY, PCW, .true., 'R')
                  CharPos = CharPos + PCW + 1
                End If
              End If
            End If
            If (.not. FieldReq%ZAcross .and. iZGrid /= 0) Then
              Line(CharPos:CharPos+PCW) = Int2Char(jZ, PCW, .true., 'R')
              CharPos = CharPos + PCW + 1
            End If
            If (.not. FieldReq%DAcross .and. iDGrid /= 0) Then
              Line(CharPos:CharPos+PCW) = Int2Char(jD, PCW, .true., 'R')
              CharPos = CharPos + PCW + 1
            End If
          End If

          ! Grid values.
          If (.not. FieldReq%TAcross .and. iTGrid /= 0) Then
            Line(CharPos:CharPos+PCW) = FormatChar(T, PCW, .true., 'R')
            CharPos = CharPos + PCW + 1
          End If
          If (.not. FieldReq%XAcross .and. iHGrid /= 0) Then
            If (HGrid%Named) Then
              Line(CharPos:CharPos+PCW) = FormatChar(HGrid%Names(jX), PCW, .true., 'R')
              CharPos = CharPos + PCW + 1
            Else
              If (FieldReq%NameII) Then
                Line(CharPos:CharPos+PCW) = Std2Char(HGrid%X(jX), PCW, .true., 'R', 'F13.6')
                CharPos = CharPos + PCW + 1
              Else
                Line(CharPos:CharPos+PCW) = Std2Char(HGrid%X(jX), PCW, .true., 'R')
                CharPos = CharPos + PCW + 1
              End If
            End If
          End If
          If (.not. FieldReq%YAcross .and. iHGrid /= 0) Then
            If (.not. HGrid%Named) Then
              If (HGrid%Unstructured) Then
                Line(CharPos:CharPos+PCW) = Std2Char(HGrid%Y(jX), PCW, .true., 'R')
                CharPos = CharPos + PCW + 1
              Else
                If (FieldReq%NameII) Then
                  Line(CharPos:CharPos+PCW) = Std2Char(HGrid%Y(jY), PCW, .true., 'R', 'F13.6')
                  CharPos = CharPos + PCW + 1
                Else
                  Line(CharPos:CharPos+PCW) = Std2Char(HGrid%Y(jY), PCW, .true., 'R')
                  CharPos = CharPos + PCW + 1
                End If
              End If
            End If
          End If
          If (.not. FieldReq%ZAcross .and. iZGrid /= 0) Then
            Line(CharPos:CharPos+PCW) = Std2Char(ZGrid%Z(jZ), PCW, .true., 'R')
            CharPos = CharPos + PCW + 1
          End If
          If (.not. FieldReq%DAcross .and. iDGrid /= 0) Then
            If (FieldReq%DGridFloating) Then
              Line(CharPos:CharPos+PCW) = FormatChar('variable', PCW, .true., 'R')
              CharPos = CharPos + PCW + 1
            Else
              Line(CharPos:CharPos+PCW) = Std2Char(DInDGrid(DGrid, jD, 0), PCW, .true., 'R')
              CharPos = CharPos + PCW + 1
            End If
          End If

          If (.not. FieldReq%NameII .and. .not. OutputOpts%Pre65Format) Then
            Line(CharPos:CharPos+PCW) = FormatChar(' ', PCW, .true., 'R')
            CharPos = CharPos + PCW + 1
          End If

          NonZeroLine = .false.

          ! Rest of line.

          ! Loop over fields in the output group.
          iLoop: Do i = 1, nNum(iDS)

            ! Abbreviations for field requirement, field, grid indices, grids and grid sizes. Grid sizes are
            ! set to 1 if the grid is missing.
            FieldReq => Reqs%FieldReqs(NumList(i, iDS))
            Field    => Results%Fields(NumList(i, iDS))
            iTGrid = FieldReq%iTGrid
            iHGrid = FieldReq%iHGrid
            iZGrid = FieldReq%iZGrid
            iDGrid = FieldReq%iDGrid
            If (iTGrid == 0) Then
              nT = 1
            Else
              TGrid => Grids%TGrids(iTGrid)
              nT = TGrid%nT
            End If
            If (iHGrid == 0) Then
              nX = 1
              nY = 1
            Else
              HGrid => Grids%HGrids(iHGrid)
              nX = HGrid%nX
              nY = HGrid%nY
            End If
            If (iZGrid == 0) Then
              nZ = 1
            Else
              ZGrid => Grids%ZGrids(iZGrid)
              nZ = ZGrid%nZ
            End If
            If (iDGrid == 0) Then
              nD = 1
            Else
              DGrid => Grids%DGrids(iDGrid)
              nD = nDInDGrid(DGrid)
            End If

            ! Determine range of T/X/Y/Z/D-grid indices for row if these vary between fields in the output
            ! group.
            If (.not. FieldReq%TSeparateFile .and. FieldReq%TAcross) Then
              kTL = Field%LastNumOutput + 1
              kTU = Field%LastProc
            End If
            If (.not. FieldReq%XSeparateFile .and. FieldReq%XAcross) Then
              kXL = 1
              kXU = nX
            End If
            If (.not. FieldReq%YSeparateFile .and. FieldReq%YAcross) Then
              kYL = 1
              kYU = nY
            End If
            If (.not. FieldReq%ZSeparateFile .and. FieldReq%ZAcross) Then
              kZL = 1
              kZU = nZ
            End If
            If (.not. FieldReq%DSeparateFile .and. FieldReq%DAcross) Then
              kDL = 1
              kDU = nD
            End If

            ! Abbreviation for field.
            If (FieldReq%iExt /= 0) Then
              Field => Results%Fields(FieldReq%iExt)
            Else
              Field => Results%Fields(NumList(i, iDS))
            End If

            ! Loop over output columns.
             jLoop: Do j  = 1,   FieldReq%nComponents * FieldReq%nRanges
            kTLoop: Do kT = kTL, kTU
            iTA = CalciTA(Field, kT)
            kXLoop: Do kX = kXL, kXU
            kYLoop: Do kY = kYL, kYU
            kZLoop: Do kZ = kZL, kZU
            kDLoop: Do kD = kDL, kDU
              If (FieldReq%DGridFloating) Then
                Line(CharPos:CharPos+FCW) = Std2Char(DInDGrid(DGrid, kD, Field%S(kX, kY, kZ, 1, iTA, j)), &
                      FCW, .true., 'R')
                CharPos = CharPos + FCW + 1
              End If
              If (FieldReq%TypeType == ' ') Then
                If (FieldReq%NameII) Then
                  Line(CharPos:CharPos+FCW) = Std2Char(Field%Std(kD, kX, kY, kZ, 1, iTA, j), FCW, &
                        .true., 'R', 'E15.8')
                  CharPos = CharPos + FCW + 1
                Else
                  Line(CharPos:CharPos+FCW) = Std2Char(Field%Std(kD, kX, kY, kZ, 1, iTA, j), FCW, .true., 'R')
                  CharPos = CharPos + FCW + 1
                End If
                If (Field%Std(kD, kX, kY, kZ, 1, iTA, j) /= 0.0) NonZeroLine = .true.
              Else If (FieldReq%TypeType == '8') Then
                Line(CharPos:CharPos+FCW) = P642Char(Field%P64(kD, kX, kY, kZ, 1, iTA, j), FCW, .true., 'R')
                CharPos = CharPos + FCW + 1
                If (Field%P64(kD, kX, kY, kZ, 1, iTA, j) /= 0.0) NonZeroLine = .true.
              Else If (FieldReq%TypeType == ':') Then
                Line(CharPos:CharPos+FCW) = FormatChar(                             &
                                              Time2Char(                            &
                                                Field%T(kD, kX, kY, kZ, 1, iTA, j), &
                                                OutputOpts%Seconds,                 &
                                                OutputOpts%DecimalPlaces,           &
                                                .true.                              &
                                              ),                                    &
                                              FCW,                                  &
                                              .true.,                               &
                                              'R'                                   &
                                            )
                CharPos = CharPos + FCW + 1
                If (IsTimeInterval(Field%T(kD, kX, kY, kZ, 1, iTA, j))) Then
                  If (Field%T(kD, kX, kY, kZ, 1, iTA, j) /= ZeroTime()) NonZeroLine = .true.
                Else
                  NonZeroLine = .true.
                End If
              End If
            End Do kDLoop
            End Do kZLoop
            End Do kYLoop
            End Do kXLoop
            End Do kTLoop
            End Do  jLoop

            LastNumOutput(NumList(i, iDS)) = kTU

          End Do iLoop

          If (NonZeroLine .or. FieldReq%ZeroLines) Then
            If (iDS == 1) Then
              Write (DiskUnit, '(A)') Trim(Line)
              nLines = nLines + 1
            Else
              Write (ScreenUnit, '(A)') Trim(Line)
            End If
          End If

        End Do iDSLoop

      End Do jDLoop
      End Do jZLoop
      End Do jYLoop
      End Do jXLoop
      End Do jTLoop

      ! For those files and windows which will not be written to again: (i) close files and (ii) remove
      ! information on files and windows from Results. Also note that (i) the test for DiskUnit /= 0 is needed
      ! because of the possibility that there was no output written (see 'If no output, cycle iDS Loop'), and
      ! (ii) the ScreenUnit isn't closed in order to leave the window open.
      If (                                                                               &
        (      FieldReq%TSeparateFile                                             ) .or. &
        (.not. FieldReq%TSeparateFile .and.       FieldReq%TAcross                ) .or. &
        (.not. FieldReq%TSeparateFile .and. .not. FieldReq%TAcross .and. jTU == nT) .or. &
        EndOfCase                                                                        &
      ) Then

        If (nNum(1) /= 0 .and. DiskUnit /= 0) Then
          Call CloseUnit(DiskUnit, Units)
          nLines   = 0
          DiskUnit = 0
        End If
        If (nNum(2) /= 0) Then
          ScreenUnit = 0
        End If

      ! For other files, flush buffer if FieldReq%Flush set.
      Else

        If (FieldReq%Flush) Then
          If (nNum(1) /= 0 .and. DiskUnit /= 0) Then
            Flush(DiskUnit)
          End If
        End If

      End If

    End Do iDLoop
    End Do iZLoop
    End Do iYLoop
    End Do iXLoop
    End Do iTLoop


  End Do OutputGroupLoop
  !$OMP END PARALLEL DO

  ! For those files and windows which were in the 'may need to be written to again' category but will now not
  ! be written to again: (i) close files and (ii) remove information on files and windows from Results.
  If (EndOfCase) Then
    Do iOutputGroup = 1, Reqs%nFieldGroups
      Do iFile = 1, MaxFieldOutputFiles
        If (Results%FieldDiskUnits(iOutputGroup, iFile) /= 0) Then
          Call CloseUnit(Results%FieldDiskUnits(iOutputGroup, iFile), Units)
        End If
        Results%FieldnLines     (iOutputGroup, iFile) = 0
        Results%FieldDiskUnits  (iOutputGroup, iFile) = 0
        Results%FieldScreenUnits(iOutputGroup, iFile) = 0
      End Do
    End Do
  End If

  !-------------------!
  ! Graphical output. !
  !-------------------!

  ! Loop over fields.
  Do iField = 1, Reqs%nFieldReqs

    ! Abbreviations for field requirement and field.
    FieldReq => Reqs%FieldReqs(iField)
    Field    => Results%Fields(iField)

    ! Determine which fields are to be output.
    If (FieldReq%Ensemble .and. .not.EndOfCases) Cycle
    If (.not. FieldReq%Graph) Cycle

    ! Abbreviations for grid indices, for grids and for grid sizes. Grid sizes are set to 1 if the grid is
    ! missing.
    iTGrid = FieldReq%iTGrid
    iHGrid = FieldReq%iHGrid
    iZGrid = FieldReq%iZGrid
    iDGrid = FieldReq%iDGrid
    If (iTGrid == 0) Then
      nT = 1
    Else
      TGrid => Grids%TGrids(iTGrid)
      nT = TGrid%nT
    End If
    If (iHGrid == 0) Then
      nX = 1
      nY = 1
    Else
      HGrid => Grids%HGrids(iHGrid)
      nX = HGrid%nX
      nY = HGrid%nY
    End If
    If (iZGrid == 0) Then
      nZ = 1
    Else
      ZGrid => Grids%ZGrids(iZGrid)
      nZ = ZGrid%nZ
    End If
    If (iDGrid == 0) Then
      nD = 1
    Else
      DGrid => Grids%DGrids(iDGrid)
      nD = nDInDGrid(DGrid)
    End If

    ! Determine range of T-grid indices for looping over output times.
    iTL = Field%LastGraphOutput + 1
    iTU = Field%LastProc

    ! Abbreviation for graph unit.
    GraphUnit => Field%GraphUnit

    ! Open window and write headers if not already open (tested by GraphUnit = 0).
    If (GraphUnit == 0 .and. iTU >= iTL) Then

      ! Window title.
      If (FieldReq%Ensemble) Then
        Title = Trim(FieldReq%Quantity) // ': Cases'
      Else
        Title = Trim(FieldReq%Quantity) // ': Case ' // Int2Char(iCase)
      End If

      ! Open window.
      Call GetNewUnit(GraphUnit, Units)
      Call OpenGraphWindow(                       &
             Unit  = GraphUnit,                   &
             Title = Trim(Title),                 &
             X1    = -0.02_Std * FieldReq%XScale, &  ! $$ numbers ?
             Y1    =   1.1_Std * FieldReq%YScale, &
             X2    =             FieldReq%XScale, &
             Y2    =  -0.6_Std * FieldReq%YScale  &
           )

      ! Write headers.
      Call FieldGraphHeaders( &
             GraphUnit,       &
             Grids, Sources,  &
             FieldReq,        &
             Results          &
           )

    End If

    ! Abbreviation for field.
    If (FieldReq%iExt /= 0) Then
      Field => Results%Fields(FieldReq%iExt)
    Else
      Field => Results%Fields(iField)
    End If

    ! Loop over output times.
    Do iT = iTL, iTU

      iTA = CalciTA(Field, iT)

      ! Note graphical output only supports the following cases
      !   (i)  nZ > 1, nD = nX = nY = nS = 1, FieldReq%nComponents = 1 (nT can be 1 or > 1)
      !   (ii) nS > 1, nD = nX = nY = nZ = 1, FieldReq%nComponents = 1 (nT can be 1 or > 1).
      If (nZ > 1 .and. nD == 1 .and. nX == 1 .and. nY == 1) Then ! $$  .and. nS == 1
        Call Plot(                                  &
               GraphUnit,                           &
               Field%Std(1, 1, 1, 1:nZ, 1, iTA, 1), &
               ZGrid%Z(1:nZ)                        &
             )
      ! Else If (nS > 1) Then
      !   Call Plot(                                  &
      !          GraphUnit,                           &
      !          Field%Std(1, 1, 1, 1, 1:nS, iTA, 1), &
      !          SGrid%S(1:nS)                        &
      !        ) ! Add support for this when adding S grids to this routine. also plots v t (e.g. sig z)?$$
      End If

      LastGraphOutput(iField) = iT

    End Do

    ! For those windows which will not be written to again, remove information from Results.
    If (iTU == nT .or. EndOfCase) Then
      GraphUnit = 0
    End If

  End Do

  !-------------------------------------------------------------------------!
  ! Update Results%Fields%LastNumOutput and Results%Fields%LastGraphOutput. !
  !-------------------------------------------------------------------------!

  Do iField = 1, Reqs%nFieldReqs
    If (LastNumOutput(iField) /= 0) Then
      Results%Fields(iField)%LastNumOutput = LastNumOutput(iField)
    End If
    If (LastGraphOutput(iField) /= 0) Then
      Results%Fields(iField)%LastGraphOutput = LastGraphOutput(iField)
    End If
  End Do

End Subroutine OutputFields

!-------------------------------------------------------------------------------------------------------------

Subroutine OutputPdfs(Time, Coords, Grids, Reqs, OutputOpts, Results, Units, iCase, Sources)
! Outputs pdfs to screen and disk (by groups).

  Implicit None
  ! Argument list:
  Type(Time_),       Intent(In)    :: Time    ! Current synchronisation time.
  Type(Coords_),     Intent(In)    :: Coords  !
  Type(Grids_),      Intent(In)    :: Grids   !
  Type(Reqs_),       Intent(In)    :: Reqs    !
  Type(OutputOpts_), Intent(In)    :: OutputOpts
  Type(Results_),    Intent(InOut) :: Results !
  Type(Units_),      Intent(InOut) :: Units   !
  Integer,           Intent(In)    :: iCase   !
  Type(Sources_),    Intent(In)    :: Sources !
  ! Locals:
  Integer        :: nScreen                !} Number of pdfs from output group that
  Integer        :: nDisk                  !} are written to the screen and disk.
  Integer        :: ScreenList(MaxPdfReqs) !] Indices of the pdfs from an output group that
  Integer        :: DiskList(MaxPdfReqs)   !] are written to screen and disk.
  Integer        :: iG                     ! Index for output group.
  Integer        :: i                      ! Index for pdf.
  Integer        :: iPdf                   ! Index of a pdf in the output group.
  Integer        :: iHGrid                 ! Index for horizontal grid.
  Integer        :: iZGrid                 ! Index for vertical grid.
  Integer        :: iTGrid                 ! Index for time grid.
  Integer        :: nX                     !} Size of horizontal grid.
  Integer        :: nY                     !}
  Integer        :: nZ                     ! Size of vertical grid.
  Integer        :: nT                     ! Size of time grid.
  Integer        :: iStart                 ! Start index of output times for output group.
  Integer        :: iT                     ! Index for T.
  Integer        :: iX                     ! Index for X.
  Integer        :: iY                     ! Index for Y.
  Integer        :: iZ                     ! Index for Z.
  Integer        :: iP                     ! Index for pdf.
  Integer        :: j                      ! Loop counter.
  Type(Time_)    :: T                      ! Output time.
  Character(MaxOutputLineLength) :: Line       !} Character strings for writing output.
  Character(MaxOutputLineLength) :: ScreenLine !}
  Character(MaxOutputLineLength) :: DiskLine   !}

  ! Output pdfs (by groups).
  Do iG = 1, Reqs%nPdfGroups

    ! Determine how pdfs in each output group are to be output.
    nScreen  = 0
    nDisk    = 0
    Do i = 1, Reqs%nPdfs
      If ((Reqs%PdfReqs(i)%iOutputGroup /= iG) .or. Reqs%PdfReqs(i)%AvEnsemble) Cycle
      If (Reqs%PdfReqs(i)%Screen) Then
        nScreen             = nScreen + 1
        ScreenList(nScreen) = i
        iPdf                = i
      End If
      If (Reqs%PdfReqs(i)%Disk) Then
        nDisk           = nDisk + 1
        DiskList(nDisk) = i
        iPdf            = i
      End If
    End Do
    If (nScreen == 0 .and. nDisk == 0) Cycle   ! No output required for this output group.

    ! Index and size of grids for this output group.
    iHGrid = Reqs%PdfReqs(iPdf)%iHGrid
    iZGrid = Reqs%PdfReqs(iPdf)%iZGrid
    iTGrid = Reqs%PdfReqs(iPdf)%iTGrid
    If (iHGrid == 0) Then
      nX = 1
      nY = 1
    Else
      nX = Grids%HGrids(iHGrid)%nX
      nY = Grids%HGrids(iHGrid)%nY
    End If
    If (iZGrid == 0) Then
      nZ = 1
    Else
      nZ = Grids%ZGrids(iZGrid)%nZ
    End If
    If (iTGrid == 0) Then     ! $$ in fact, TGrid must exist at present.
      nT = 1
    Else
      nT = Grids%TGrids(iTGrid)%nT
    End If

    ! Loop over times in temporal grid starting at the next required output time
    ! for this output group.
    iStart = Results%PdfsNextOutput(iG)
    Do iT = iStart, nT
      If (Time2ShortTime(Time) < TInTGrid(Grids%TGrids(iTGrid), iT)) Exit

      ! Numerical screen and disk output.
      If (nDisk /= 0) Then
        If (Reqs%PdfReqs(DiskList(1))%SeparateFile == 'T' .or. iT == 1) Then
          T = ShortTime2Time(TInTGrid(Grids%TGrids(iTGrid), iT))
          Call PdfHeaders(Reqs, OutputOpts, iCase, iG, Units, Results, T, Grids, Sources)
        End If
      End If

      Do iX = 1, nX
      Do iY = 1, nY
      Do iZ = 1, nZ
        Line = ' '
        T = ShortTime2Time(TInTGrid(Grids%TGrids(iTGrid), iT))

        If (iTGrid /= 0) Line =                              &
                Trim(Line)                                // &
                FormatChar(Time2Char(T, OutputOpts%Seconds, OutputOpts%DecimalPlaces, .true.), &
                      OutputPrelimColumnWidth, .true., 'R')
        If (iHGrid /= 0) Line = &
                Trim(Line)                                                                 // &
                Std2Char(Grids%HGrids(iHGrid)%X(iX), OutputPrelimColumnWidth, .true., 'R')
        If (iHGrid /= 0) Line =                                          &
                Trim(Line)                                            // &
                Std2Char(Grids%HGrids(iHGrid)%Y(iY), OutputPrelimColumnWidth, .true., 'R')
        If (iZGrid /= 0) Line =                                          &
                Trim(Line)                                            // &
                Std2Char(Grids%ZGrids(iZGrid)%Z(iZ), OutputPrelimColumnWidth, .true., 'R')
        If (nScreen /= 0) Then
          ScreenLine = Line
          Do j = 1, nScreen
            Do iP = 1, Results%Pdfs(ScreenList(j))%PdfSize
              ScreenLine = Trim(ScreenLine)                                    // &
                Std2Char(Results%Pdfs(ScreenList(j))%Data(iX, iY, iZ, iT, iP),    &
                15, .true., 'R')
            End Do
            If (.not. Results%Pdfs(ScreenList(j))%FixedThresholds)             &
              ScreenLine = Trim(ScreenLine)                                 // &
                Int2Char(Results%Pdfs(ScreenList(j))%Scale(iX, iY, iZ, iT),    & ! $$
                15, .true., 'R')
          End Do
          Write (Results%PdfScreenUnits(iG), *) ScreenLine(1:Len_Trim(ScreenLine) - 1)
        End If
        If (nDisk /= 0) Then
          DiskLine = Line
          Do j = 1, nDisk
            Do iP = 1, Results%Pdfs(DiskList(j))%PdfSize
              DiskLine = Trim(DiskLine)                                      // &
                Std2Char(Results%Pdfs(DiskList(j))%Data(iX, iY, iZ, iT, iP),    &
                15, .true., 'R')
            End Do
            If (.not. Results%Pdfs(DiskList(j))%FixedThresholds)             &
              DiskLine = Trim(DiskLine)                                   // &
                Int2Char(Results%Pdfs(DiskList(j))%Scale(iX, iY, iZ, iT),    & ! $$
                15, .true., 'R')
          End Do
          Write (Results%PdfDiskUnits(iG), *) DiskLine(1:Len_Trim(DiskLine) - 1)
        End If
      End Do
      End Do
      End Do

      If (nDisk /= 0) Then
        If (Reqs%PdfReqs(DiskList(1))%SeparateFile == 'T') Then
          Call CloseUnit(Results%PdfDiskUnits(iG), Units)
        End If
      End If

      Results%PdfsNextOutput(iG) = Results%PdfsNextOutput(iG) + 1

    End Do

  End Do

End Subroutine OutputPdfs

!-------------------------------------------------------------------------------------------------------------

Subroutine OutputPPInfos(                                    &
             Line,                                           &
             X, XOld, T, TOld,                               &
             iCase, iPPInfo, iT, iUP, Puff,                  &
             OutputOpts,                                     &
             Coords, Grids, Flows, Specieses, Sources, Reqs, &
             Units, Results                                  &
           )
! Outputs particle/puff information.

! Note that DiskUnit and ScreenUnit determine whether a file/window needs to be opened (a zero value indicates
! it needs to be opened) while nLines determines whether, after opening a file, headers need to be written or
! whether the file needs to be positioned after the first nLines lines. Normally opening is followed by
! writing the headers - however this may not be the case after a restart. On some operating systems all screen
! output goes to a single terminal window which is associated with unit 6 (in which case the question of
! opening the window doesn't arise).

  Implicit None
  ! Argument list:
  Character(MaxOutputLineLength), Intent(In)            :: Line
  Real(Std),                      Intent(In)            :: X(3)
  Real(Std),                      Intent(In)            :: XOld(3)
  Type(ShortTime_),               Intent(In)            :: T
  Type(ShortTime_),               Intent(In)            :: TOld
  Integer,                        Intent(In)            :: iCase
  Integer,                        Intent(In)            :: iPPInfo
  Integer,                        Intent(In)            :: iT
  Integer,                        Intent(In)            :: iUP
  Logical,                        Intent(In)            :: Puff
  Type(OutputOpts_),              Intent(In)            :: OutputOpts
  Type(Coords_),                  Intent(In)            :: Coords
  Type(Grids_),                   Intent(In),    Target :: Grids
  Type(Flows_),                   Intent(In)            :: Flows
  Type(Specieses_),               Intent(In)            :: Specieses
  Type(Sources_),                 Intent(In)            :: Sources
  Type(Reqs_),                    Intent(In),    Target :: Reqs
  Type(Units_),                   Intent(InOut)         :: Units
  Type(Results_),                 Intent(InOut), Target :: Results
  ! Line       :: Character string with information to be output.
  ! iCase      :: Number of case.
  ! iPPInfo    :: Index of set of particle/puff information.
  ! iT         :: Index of time in T-grid.
  ! iUP        :: Unique index of particle/puff (i.e. an index unique to the particle/puff which is not reused
  !            :: despite particle recycling).
  ! Puff       :: Indicates information on a puff is to be output.
  ! OutputOpts :: Output options.
  ! Coords     :: Collection of coord systems.
  ! Grids      :: Collection of grids.
  ! Flows      :: Collection of flow module instance states.
  ! Specieses  :: Collection of specieses.
  ! Sources    :: Collection of sources.
  ! Reqs       :: Collection of requirements.
  ! Units      :: Collection of information on input/output unit numbers.
  ! Results    :: Collection of results.
  ! Locals:
  Integer                              :: iFile
  Integer,                     Pointer :: nLines
  Integer,                     Pointer :: DiskUnit
  Integer,                     Pointer :: ScreenUnit
  Integer,                     Pointer :: GraphUnit
  Character(MaxFileNameLength)         :: File
  Character(MaxCharLength3)            :: Title
  Type(TGrid_),                Pointer :: TGrid
  Type(PPInfoReq_),            Pointer :: PPInfoReq
  Type(PPInfo_),               Pointer :: PPInfo
  Integer                              :: iTA
  Integer                              :: iUPA
  Integer                              :: i
  Integer                              :: IOStat
  Integer                              :: DummynLines
  ! iFile       :: Index of file, corresponding to the index of the arrays PPInfo%nLines, PPInfo%DiskUnits,
  !                PPInfo%ScreenUnits, PPInfo%ParticleFiles and PPInfo%PuffFiles.
  ! nLines      :} Abbreviations for number of lines written to file, disk unit and screen unit for numerical
  ! DiskUnit    :} output, and unit for graphical output.
  ! ScreenUnit  :}
  ! GraphUnit   :}
  ! File        :: File name for numerical output.
  ! Title       :: Window title.
  ! TGrid       :} Abbreviations for grids, requirements and results.
  ! PPInfoReq   :}
  ! PPInfo      :}
  ! iTA         :: Adjusted index of time in T-grid.
  ! iUPA        :: Adjusted unique index of particle/puff.
  ! i           :: Loop index.
  ! IOStat      :: Error code for read statements.
  ! DummynLines :: Dummy variable needed for an argument list.

  ! Abbreviations.
  PPInfoReq => Reqs%PPInfoReqs(iPPInfo)
  PPInfo    => Results%PPInfos(iPPInfo)

  !-----------------------------------!
  ! Numerical disk and screen output. !
  !-----------------------------------!

  If (PPInfoReq%Disk .or. PPInfoReq%Screen) Then

    ! Calculate iFile.
    If (PPInfoReq%PSeparateFile) Then

      If (Puff) Then
        If (PPInfoReq%Puffs == 'S') Then
          iUPA = iUP - PPInfoReq%FirstPuff + 1
        Else
          iUPA = iUP
        End If
        If (iUPA > MaxPPInfoOutputFiles) Then
          Call Message(                                                                     &
                 'FATAL ERROR: too many output files/windows required at the same time ' // &
                 'for the set of particle/puff information "'                            // &
                 Trim(PPInfoReq%Name)                                                    // &
                 '"',                                                                       &
                 3                                                                          &
               )
        End If
        If (PPInfo%PuffFiles(iUPA) == 0) Then
          If (PPInfo%LastFile == MaxPPInfoOutputFiles) Then
            Call Message(                                                                     &
                   'FATAL ERROR: too many output files/windows required at the same time ' // &
                   'for the set of particle/puff information "'                            // &
                   Trim(PPInfoReq%Name)                                                    // &
                   '"',                                                                       &
                   3                                                                          &
                 )
          End If
          PPInfo%LastFile        = PPInfo%LastFile + 1
          PPInfo%PuffFiles(iUPA) = PPInfo%LastFile
        End If
        iFile = PPInfo%PuffFiles(iUPA)
      Else
        If (PPInfoReq%Particles == 'S') Then
          iUPA = iUP - PPInfoReq%FirstParticle + 1
        Else
          iUPA = iUP
        End If
        If (iUPA > MaxPPInfoOutputFiles) Then
          Call Message(                                                                     &
                 'FATAL ERROR: too many output files/windows required at the same time ' // &
                 'for the set of particle/puff information "'                            // &
                 Trim(PPInfoReq%Name)                                                    // &
                 '"',                                                                       &
                 3                                                                          &
               )
        End If
        If (PPInfo%ParticleFiles(iUPA) == 0) Then
          If (PPInfo%LastFile == MaxPPInfoOutputFiles) Then
            Call Message(                                                                     &
                   'FATAL ERROR: too many output files/windows required at the same time ' // &
                   'for the set of particle/puff information "'                            // &
                   Trim(PPInfoReq%Name)                                                    // &
                   '"',                                                                       &
                   3                                                                          &
                 )
          End If
          PPInfo%LastFile            = PPInfo%LastFile + 1
          PPInfo%ParticleFiles(iUPA) = PPInfo%LastFile
        End If
        iFile = PPInfo%ParticleFiles(iUPA)
      End If

    Else

      If (PPInfoReq%TSeparateFile) Then

        If (iT < PPInfo%iTL .or. iT > PPInfo%iTL + PPInfo%nT - 1) Then
          Call Message('UNEXPECTED FATAL ERROR in OutputPPInfos', 4)
        End If

        iTA = Mod(iT, PPInfo%nT)
        If (iTA <= 0) iTA = iTA + PPInfo%nT

      Else

        iTA = 1

      End If

      iFile = iTA

      PPInfo%LastFile = Max(PPInfo%LastFile, iFile)

    End If

    ! Abbreviations.
    nLines     => PPInfo%nLines     (iFile)
    DiskUnit   => PPInfo%DiskUnits  (iFile)
    ScreenUnit => PPInfo%ScreenUnits(iFile)

    ! Open window and write headers for disk output.
    If (PPInfoReq%Disk .and. DiskUnit == 0) Then

      ! File name.
      File = Trim(OutputOpts%Folder) // Trim(PPInfoReq%Name) // '_C' // Trim(Int2Char(iCase))
      If (PPInfoReq%PSeparateFile) Then
        If (Puff) Then
          File = Trim(File) // '_Puff' // Trim(Int2Char(iUP))
        Else
          File = Trim(File) // '_Particle' // Trim(Int2Char(iUP))
        End If
      Else If (PPInfoReq%TSeparateFile) Then
        TGrid => Grids%TGrids(PPInfoReq%iTGrid)
        File = Trim(File) // '_T' // Trim(Int2Char(iT)) // '_' //  &
               Trim(                                               &
                 FileNameTime(                                     &
                   ChangeTimeZone(                                 &
                     ShortTime2Time(TInTGrid(TGrid, iT)), TGrid%T0 &
                   )                                               &
                 )                                                 &
               )
      End If
      File = Trim(File) // '.txt'
      File = ConvertFileName(File)

      ! If nLines = 0 then either (i) the file hasn't yet been opened or written to, or (ii) the run has been
      ! restarted and any previously written lines are to be ignored.
      If (nLines == 0) Then

        ! Open file.
        DiskUnit = OpenFile(                                &
                     File            = File,                &
                     Units           = Units,               &
                     Status          = 'Replace',           &
                     RecL            = MaxOutputLineLength, &
                     FileDescription = 'output file'        &
                   )

        ! Write headers.
        Call PPInfoNumHeaders(        &
               DiskUnit,              &
               OutputOpts,            &
               Grids, Flows, Sources, &
               PPInfoReq,             &
               Results,               &
               nLines                 &
             )

        ! Write column headers.
        Call PPInfoColumnHeaders( &
               DiskUnit,          &
               OutputOpts,        &
               Coords, Specieses, &
               PPInfoReq,         &
               nLines             &
             )

      ! If nLines /= 0 then the run has been restarted and the first nLines previously written lines are to be
      ! preserved and the rest are to be ignored.
      Else

        ! Open file.
        DiskUnit = OpenFile(                                &
                     File            = File,                &
                     Units           = Units,               &
                     Status          = 'Old',               &
                     RecL            = MaxOutputLineLength, &
                     FileDescription = 'output file'        &
                   )

        ! Read first nLines lines.
        Do i = 1, nLines
          Read (DiskUnit, *, IOStat = IOStat)
          If (IOStat > 0) Then
            Call Message(                                                 &
                   'FATAL ERROR: An error occurred reading the file "' // &
                   Trim(File)                                          // &
                   '" following restart',                                 &
                   3                                                      &
                 )
          Else If (IOStat < 0) Then
            Call Message(                                             &
                   'FATAL ERROR: The file "'                       // &
                   Trim(File)                                      // &
                   '" is shorter than expected following restart',    &
                   3                                                  &
                 )
          End If
        End Do

      End If

    End If

    ! Open window and write headers for screen output.
    If (PPInfoReq%Screen .and. ScreenUnit == 0) Then

      ! Window title.
      Title = Trim(PPInfoReq%Name) // ': Case ' // Int2Char(iCase)
      If (PPInfoReq%PSeparateFile) Then
        If (Puff) Then
          Title = Trim(Title) // ': Puff ' // Trim(Int2Char(iUP))
        Else
          Title = Trim(Title) // ': Particle ' // Trim(Int2Char(iUP))
        End If
      Else If (PPInfoReq%TSeparateFile) Then
        Title = Trim(Title) // ': T ' // Trim(Int2Char(iT))
      End If

      ! Open window. In connection with nCols, note that readability is improved by making the window one
      ! column wider than the longest record.
#     ifdef CompaqPCCompiler
        Call GetNewUnit(ScreenUnit, Units)
        Call OpenTextWindow(                   &
               Unit  = ScreenUnit,             &
               Title = Trim(Title),            &
               nRows = 200,                    & ! $$ Control over nrows (param initially)?
               nCols = MaxOutputLineLength + 1 & ! $$ different values for disk and screen?
             )                                   ! $$ checking lines long enough (disk too)
#     endif
#     ifdef IntelLinCompiler
        ScreenUnit = 6
#     endif
#     ifdef CrayCLECompiler
        ScreenUnit = 6
#     endif
#     ifdef sun
        ScreenUnit = 6
#     endif

      ! Write column headers.
      Call PPInfoColumnHeaders( &
             ScreenUnit,        &
             OutputOpts,        &
             Coords, Specieses, &
             PPInfoReq,         &
             DummynLines        &
           )

    End If

    ! Write data - disk output.
    If (PPInfoReq%Disk) Then
      Write (DiskUnit, *) Line(1:Len_Trim(Line) - 1)
      nLines = nLines + 1
    End If

    ! Write data - screen output.
    If (PPInfoReq%Screen) Then
      Write (ScreenUnit, *) Line(1:Len_Trim(Line) - 1)
    End If

  End If

  !-------------------!
  ! Graphical output. !
  !-------------------!

  If (PPInfoReq%Graph) Then

    ! Abbreviation for graph unit.
    GraphUnit => PPInfo%GraphUnit

    ! Open window and write headers if not already open (tested by GraphUnit = 0).
    If (GraphUnit == 0) Then

      ! Window title.
      Title = 'Particle/puff trajectories: Case ' // Int2Char(iCase)

      ! Open window.
      Call GetNewUnit(GraphUnit, Units)
      Call OpenGraphWindow(             &
             Unit  = GraphUnit,         &
             Title = Trim(Title),       &
             X1    = -0.02_Std * 200.0, &  ! $$ numbers ? Input scales as for fields
             Y1    =   1.1_Std * 200.0, &
             X2    =             200.0, &
             Y2    =  -0.6_Std * 200.0  &
           )

      ! Write headers.
 !     Call FieldGraphHeaders( & $$ replace with analogous PPInfoGraphHeaders
 !            GraphUnit,       &
 !            Grids, Sources,  &
 !            FieldReq,        &
 !            Results          &
 !          )

    End If

      ! Note graphical output only supports the following cases
      !   (i) Z v S.
      ! $$ add other cases? Make checks that output options chosen are consistent with limitations?
      ! $$ add scaling options
      Call Plot(                                                      &
               GraphUnit,                                             &
               (/ ShortTime2RealTime(TOld), ShortTime2RealTime(T) /), &
               (/ XOld(3), X(3) /)                                    &
             )

  End If

End Subroutine OutputPPInfos

!-------------------------------------------------------------------------------------------------------------

Subroutine ClosePPInfos(    &
             EndOfCase,     &
             Time, IncTime, &
             Grids, Reqs,   &
             Units, Results &
           )
! Closes and/or flushes files used for output of particle/puff information.

  Implicit None
  ! Argument list:
  Logical,          Intent(In)            :: EndOfCase
  Type(ShortTime_), Intent(In)            :: Time
  Logical,          Intent(In)            :: IncTime
  Type(Grids_),     Intent(In),    Target :: Grids
  Type(Reqs_),      Intent(In),    Target :: Reqs
  Type(Units_),     Intent(InOut)         :: Units
  Type(Results_),   Intent(InOut), Target :: Results
  ! EndOfCase :: Indicates the end of the case.
  ! Time      :: Time up to which output can be completed.
  ! IncTime   :: Indicates period for which output can be completed includes Time.
  ! Grids     :: Collection of grids.
  ! Reqs      :: Collection of requirements.
  ! Units     :: Collection of information on input/output unit numbers.
  ! Results   :: Collection of results.
  ! Locals:
  Integer                   :: iFile
  Integer,          Pointer :: nLines
  Integer,          Pointer :: DiskUnit
  Integer,          Pointer :: ScreenUnit
  Integer,          Pointer :: GraphUnit
  Type(PPInfoReq_), Pointer :: PPInfoReq
  Type(PPInfo_),    Pointer :: PPInfo
  Integer                   :: iPPInfo
  Integer                   :: nT
  Integer                   :: iTFromTime
  Integer                   :: iTFromFile
  Type(ShortTime_)          :: DummyT
  ! iFile      :: Index of file, corresponding to the index of the arrays PPInfo%nLines, PPInfo%DiskUnits,
  !               PPInfo%ScreenUnits, PPInfo%ParticleFiles and PPInfo%PuffFiles.
  ! nLines     :} Abbreviations for number of lines written to file, disk unit and screen unit for numerical
  ! DiskUnit   :} output, and unit for graphical output.
  ! ScreenUnit :}
  ! GraphUnit  :}
  ! PPInfoReq  :] Abbreviations for requirements and results.
  ! PPInfo     :]
  ! iPPInfo    :: Index of set of particle/puff information.
  ! nT         :: Number of points in time grid.
  ! iTFromTime :: Index of last point in time grid which is (possibly strictly) before Time.
  ! iTFromFile :: Index of point in time grid corresponding to a particular file.
  ! DummyT     :: Dummy variable needed for an argument list.

  Do iPPInfo = 1, Reqs%nPPInfoReqs

    ! Abbreviations.
    PPInfoReq => Reqs%PPInfoReqs(iPPInfo)
    PPInfo    => Results%PPInfos(iPPInfo)

    !-----------------------------------!
    ! Numerical disk and screen output. !
    !-----------------------------------!

    If (PPInfoReq%Disk .or. PPInfoReq%Screen) Then

      ! nT and iTFromTime.
      If (PPInfoReq%iTGrid /= 0) Then
        Call LastTBeforeT(Grids%TGrids(PPInfoReq%iTGrid), Time, .not.IncTime, DummyT, iTFromTime)
        nT = Grids%TGrids(PPInfoReq%iTGrid)%nT
      End If

      Do iFile = 1, PPInfo%LastFile

        ! Abbreviations.
        nLines     => PPInfo%nLines     (iFile)
        DiskUnit   => PPInfo%DiskUnits  (iFile)
        ScreenUnit => PPInfo%ScreenUnits(iFile)

        ! iTFromFile.
        If (PPInfoReq%TSeparateFile) Then
          iTFromFile = iFile + PPInfo%iTL - Mod(PPInfo%iTL, PPInfo%nT)
          Do While (iTFromFile < PPInfo%iTL)
            iTFromFile = iTFromFile + PPInfo%nT
          End Do
          Do While (iTFromFile > PPInfo%iTL + PPInfo%nT - 1)
            iTFromFile = iTFromFile - PPInfo%nT
          End Do
        End If

        ! For those files and windows which will not be written to again: (i) close files and (ii) remove
        ! information on files and windows from Results. Also note that (i) the test for DiskUnit /= 0 is
        ! needed because of the possibility that there was no output written, and (ii) the ScreenUnit isn't
        ! closed in order to leave the window open.
        If (                                                                                              &
          (PPInfoReq%iTGrid /= 0 .and.       PPInfoReq%TSeparateFile .and. iTFromFile <= iTFromTime) .or. &
          (PPInfoReq%iTGrid /= 0 .and. .not. PPInfoReq%TSeparateFile .and. nT         <= iTFromTime) .or. &
          EndOfCase                                                                                       &
        ) Then

          If (DiskUnit /= 0) Then
            Call CloseUnit(DiskUnit, Units)
            nLines   = 0
            DiskUnit = 0
          End If
          ScreenUnit = 0

        ! For other files, flush buffer if PPInfoReq%Flush set.
        Else

          If (PPInfoReq%Flush) Then
            If (DiskUnit /= 0) Then
              Flush(DiskUnit)
            End If
          End If

        End If

      End Do

    End If

    !-------------------!
    ! Graphical output. !
    !-------------------!

    If (PPInfoReq%Graph) Then

      ! Abbreviation for graph unit.
      GraphUnit => PPInfo%GraphUnit

      ! For those windows which will not be written to again, remove information from Results.
      If (EndOfCase) Then
        GraphUnit = 0
      End If

    End If

  End Do

End Subroutine ClosePPInfos

!-------------------------------------------------------------------------------------------------------------

Subroutine FieldNumHeaders(                                   &
             nPrelimCols, nFieldCols, Unit,                   &
             RadioactiveDecay, Chemistry, StartTime, EndTime, &
             OutputOpts,                                      &
             Coords, Grids, Flows, Sources,                   &
             FieldReq,                                        &
             Results,                                         &
             nLines                                           &
           )
! Writes headers for numerical output of fields.

  Implicit None
  ! Argument List:
  Integer,           Intent(In)           :: nPrelimCols
  Integer,           Intent(In)           :: nFieldCols
  Integer,           Intent(In)           :: Unit
  Logical,           Intent(In)           :: RadioactiveDecay
  Logical,           Intent(In)           :: Chemistry
  Type(ShortTime_),  Intent(In)           :: StartTime
  Type(ShortTime_),  Intent(In)           :: EndTime
  Type(OutputOpts_), Intent(In)           :: OutputOpts
  Type(Coords_),     Intent(In),   Target :: Coords
  Type(Grids_),      Intent(In),   Target :: Grids
  Type(Flows_),      Intent(In)           :: Flows
  Type(Sources_),    Intent(In)           :: Sources
  Type(FieldReq_),   Intent(In)           :: FieldReq
  Type(Results_),    Intent(In)           :: Results
  Integer,           Intent(InOut)        :: nLines
  ! nPrelimCols      :} Number of preliminary and field columns (columns of comma
  ! nFieldCols       :} separated data not columns of characters) in numerical output.
  ! Unit             :: Output unit.
  ! RadioactiveDecay :: Indicates that radioactive decay is modelled.
  ! Chemistry        :: Indicates that chemistry is modelled.
  ! StartTime        :} Start time and expected end time of the run.
  ! EndTime          :}
  ! OutputOpts       :: Output options.
  ! Coords           :: Collection of coord systems.
  ! Grids            :: Collection of grids.
  ! Flows            :: Collection of flow module instance states.
  ! Sources          :: Collection of sources.
  ! FieldReq         :: Field requirement.
  ! Results          :: Collection of results.
  ! nLines           :: Number of lines written to file.
  ! Locals:
  Type(ZCoord_),                 Pointer :: ZCoord !
  Type(HGrid_),                  Pointer :: HGrid  !
  Integer                                :: i      ! Loop index.
  Character(1)                           :: EW     !
  Character(1)                           :: NS     !
  Character(MaxOutputLineLength)         :: Line   ! Character string for writing output.
  Integer                  :: nFlows
  Character(MaxCharLength) :: FlowNames(MaxFlows)
  Type(Time_) :: SourceStartTime
  Type(Time_) :: SourceStopTime

  ! Note PV-Wave require some characters (e.g. blanks) in the second column of headers.

  If (FieldReq%NameII) Then

    Write (Unit, '(A)') ModelVersion

    Write (Unit, *) 'Title:               ' // Trim(Results%RunName)

    Line = 'Run time:            ' // Trim(Time2Char(Results%StartClockTime, .true., 3, .true., .true.))
    Write (Unit, *) Trim(Line)

    Call FlowList(Flows, nFlows, FlowNames)
    Line = ' '
    Do i = 1, Min(nFlows, 2)
      If (i == 1) Then
        Line = FlowNames(i)
      Else
        Line = Trim(Line) // '; ' // FlowNames(i)
      End If
    End Do
    Write (Unit, *) 'Met data:            ', Trim(Line)

    ! Source Information

    SourceStartTime = InfFutureTime()
    SourceStopTime  = InfPastTime  ()
    Do i = 1, Sources%nSources
      SourceStartTime = TMin(SourceStartTime, Sources%Sources(i)%StartTime)
      SourceStopTime  = TMax(SourceStopTime,  Sources%Sources(i)%StopTime)
    End Do
    Line = 'Start of release:     '       // & ! Note this extra space needed by
           Trim(                            &  ! pv-wave for name II format output
             Time2Char(                     &  ! with backwards runs.
               SourceStartTime,             &
               OutputOpts%Seconds,          &
               OutputOpts%DecimalPlaces,    &
               .true.,                      &
               .true.                       &
             )                              &
           )
    Write (Unit, *) Trim(Line)
    Line = 'End of release:       '       // &
           Trim(                            &
             Time2Char(                     &
               SourceStopTime,              &
               OutputOpts%Seconds,          &
               OutputOpts%DecimalPlaces,    &
               .true.,                      &
               .true.                       &
             )                              &
           )
    Write (Unit, *) Trim(Line)

    ! No sources
    If (Sources%nSources == 0) Then

      Write (Unit, *) 'Release rate:        ', 'No Sources'
      Write (Unit, *) 'Release location:    ', 'No Sources'
      Write (Unit, *) 'Release height:      ', 'No Sources'

    ! Multiple sources
    Else If (Sources%nSources > 1) Then

      Write (Unit, *) 'Release rate:        ', 'Multiple Sources'
      Write (Unit, *) 'Release location:    ', 'Multiple Sources'
      Write (Unit, *) 'Release height:      ', 'Multiple Sources'

    ! Single source
    Else

      If (Sources%Sources(1)%nSpecieses /= 1) Then
        Write (Unit, *) 'Release rate:              ', &
                        'Multiple Species'
      Else
        Line = 'Release rate:        '                               // &
               Trim(Std2Char(Sources%Sources(1)%SourceStrengths(1))) // &
               Trim(Sources%Sources(1)%MaterialUnitName(1)) 
        If (Sources%Sources(1)%SourceRate(1)) Line = Trim(Line) // '/s'
        Write (Unit, *) Trim(Line)
      End If

      ! Source lat-long position
      If (Sources%Sources(1)%X(1) >= 0.0) Then
        EW = 'E'
      Else
        EW = 'W'
      EndIf
      If (Sources%Sources(1)%X(2) >= 0.0) Then
        NS = 'N'
      Else
        NS = 'S'
      EndIf

      Line = 'Release location:    '                                             // &
             Trim(Std2Char(Abs(Sources%Sources(1)%X(1)), FormatString = 'F8.4')) // &
             EW                                                                  // &
             '   '                                                               // &
             Trim(Std2Char(Abs(Sources%Sources(1)%X(2)), FormatString = 'F8.4')) // &
             NS
      Write (Unit, *) Trim(Line)

      ! Source height
      ZCoord => Coords%ZCoords(Sources%Sources(1)%iZCoord)
      If (Sources%Sources(1)%dZMetres) Then
        Line =                                                                  &
          'Release height:      '                                            // &
          Trim(Std2Char(Sources%Sources(1)%X(3), FormatString = 'f10.3'))    // &
          Trim(ZCoord%Name)                                                  // &
          ' +/- '                                                            // &
          Trim(Std2Char(Sources%Sources(1)%DX(3)/2, FormatString = 'f10.3')) // &
          'm'
      Else
        Line =                                                                                            &
          'Release height:      '                                                                      // &
          Trim(Std2Char(Sources%Sources(1)%X(3) - Sources%Sources(1)%DX(3)/2, FormatString = 'f10.3')) // &
          ' to '                                                                                       // &
          Trim(Std2Char(Sources%Sources(1)%X(3) + Sources%Sources(1)%DX(3)/2, FormatString = 'f10.3')) // &
          Trim(ZCoord%Name)
      End If
      Write (Unit, *) Trim(Line)

    EndIf

    If (IsInfFuture(EndTime - StartTime)) Then
      Line = 'infinity'
    Else
      Line = Trim(Int2Char(NInt(ShortTime2RealTime(EndTime - StartTime) / 3600.0))) // ' hours'
    End If
    Write (Unit, *) 'Forecast duration:   ' // Trim(Line)

    ! Output grid information ! note iHgrid may vary if across $$
    If (FieldReq%iHGrid /= 0) Then
      HGrid => Grids%HGrids(FieldReq%iHGrid)
      If (FieldReq%NameII .and. .not. FieldReq%TAcross) Then
        ! NameII time series
        Write (Unit, *) 'X grid resolution:   ', HGrid%dx
        Write (Unit, *) 'Y grid resolution:   ', HGrid%dy
      Else
        If (.not. HGrid%Unstructured) Then
          Write (Unit, *) 'X grid origin:       ', HGrid%x0 - 0.5 * HGrid%dx
          Write (Unit, *) 'Y grid origin:       ', HGrid%y0 - 0.5 * HGrid%dy
        Else
          Write (Unit, *) 'X grid origin:       ', ' '
          Write (Unit, *) 'Y grid origin:       ', ' '
        End If
        Write (Unit, *) 'X grid size:         ', HGrid%nx
        Write (Unit, *) 'Y grid size:         ', HGrid%ny
        Write (Unit, *) 'X grid resolution:   ', HGrid%dx
        Write (Unit, *) 'Y grid resolution:   ', HGrid%dy
      End If
    Else
      If (FieldReq%NameII .and. .not. FieldReq%TAcross) Then
        Write (Unit, *) 'X grid resolution:          ', 'No H Grid'
        Write (Unit, *) 'Y grid resolution:          ', 'No H Grid'
      Else
        Write (Unit, *) 'X grid origin:              ', 'No H Grid'
        Write (Unit, *) 'Y grid origin:              ', 'No H Grid'
        Write (Unit, *) 'X grid size:                ', 'No H Grid'
        Write (Unit, *) 'Y grid size:                ', 'No H Grid'
        Write (Unit, *) 'X grid resolution:          ', 'No H Grid'
        Write (Unit, *) 'Y grid resolution:          ', 'No H Grid'
      End If
    End If

    ! Number of fields in output file
    If (.not. FieldReq%TAcross) Then
      Write (Unit, *) 'Number of series:    ', nFieldCols
    Else
      Write (Unit, *) 'Number of fields:    ', nFieldCols
    End If

    ! Blank line
    Write (Unit, *)

    If (.not. FieldReq%TAcross) Then
      nLines = nLines + 14
    Else
      nLines = nLines + 18
    End If

  Else

    Write (Unit, *) ModelVersion

    Write (Unit, *) 'Run name:                   ' // Trim(Results%RunName)

    Line = 'Run time:                   ' // Trim(Time2Char(Results%StartClockTime, .true., 3, .true.))
    Write (Unit, *) Trim(Line)

    Call FlowList(Flows, nFlows, FlowNames)
    Line = ' '
    Do i = 1, nFlows
      If (i == 1) Then
        Line = FlowNames(i)
      Else
        Line = Trim(Line) // '; ' // FlowNames(i)
      End If
    End Do
    Write (Unit, *) 'Met data:                   ', Trim(Line)

    ! Source Information

    SourceStartTime = InfFutureTime()
    SourceStopTime  = InfPastTime  ()
    Do i = 1, Sources%nSources
      SourceStartTime = TMin(SourceStartTime, Sources%Sources(i)%StartTime)
      SourceStopTime  = TMax(SourceStopTime,  Sources%Sources(i)%StopTime)
    End Do
    Line = 'Start of release:           ' // &
           Trim(                             &
             Time2Char(                      &
               SourceStartTime,              &
               OutputOpts%Seconds,           &
               OutputOpts%DecimalPlaces,     &
               .true.                        &
             )                               &
           )
    Write (Unit, *) Trim(Line)
    Line = 'End of release:             ' // &
           Trim(                             &
             Time2Char(                      &
               SourceStopTime,               &
               OutputOpts%Seconds,           &
               OutputOpts%DecimalPlaces,     &
               .true.                        &
             )                               &
           )
    Write (Unit, *) Trim(Line)

    ! No sources
    If (Sources%nSources == 0) Then

      Write (Unit, *) 'Source strength:            ', 'No Sources'
      Write (Unit, *) 'Release location:           ', 'No Sources'
      Write (Unit, *) 'Release height:             ', 'No Sources'

    ! Multiple sources
    Else If (Sources%nSources > 1) Then

      Write (Unit, *) 'Source strength:            ', 'Multiple Sources'
      Write (Unit, *) 'Release location:           ', 'Multiple Sources'
      Write (Unit, *) 'Release height:             ', 'Multiple Sources'

    ! Single source
    Else

      If (Sources%Sources(1)%nSpecieses /= 1) Then
        Write (Unit, *) 'Source strength:            ', &
                        'Multiple Species'
      Else
        Line = 'Source strength:            '                        // &
               Trim(Std2Char(Sources%Sources(1)%SourceStrengths(1))) // &
               ' '                                                   // &
               Trim(Sources%Sources(1)%MaterialUnitName(1)) 
        If (Sources%Sources(1)%SourceRate(1)) Line = Trim(Line) // ' / s'
        Write (Unit, *) Trim(Line)
      End If

      ! Source lat-long position
      If (Sources%Sources(1)%X(1) >= 0.0) Then
        EW = 'E'
      Else
        EW = 'W'
      EndIf
      If (Sources%Sources(1)%X(2) >= 0.0) Then
        NS = 'N'
      Else
        NS = 'S'
      EndIf

      Line = 'Release location:           '                                      // &
             Trim(Std2Char(Abs(Sources%Sources(1)%X(1)), FormatString = 'F8.4')) // &
             EW                                                                  // &
             '   '                                                               // &
             Trim(Std2Char(Abs(Sources%Sources(1)%X(2)), FormatString = 'F8.4')) // &
             NS
      Write (Unit, *) Trim(Line)

      ! Source height
      ZCoord => Coords%ZCoords(Sources%Sources(1)%iZCoord)
      If (Sources%Sources(1)%dZMetres) Then
        Line =                                                                  &
          'Release height:             '                                     // &
          Trim(Std2Char(Sources%Sources(1)%X(3), FormatString = 'f10.3'))    // &
          Trim(ZCoord%Name)                                                  // &
          ' +/- '                                                            // &
          Trim(Std2Char(Sources%Sources(1)%DX(3)/2, FormatString = 'f10.3')) // &
          'm'
      Else
        Line =                                                                                            &
          'Release height:             '                                                               // &
          Trim(Std2Char(Sources%Sources(1)%X(3) - Sources%Sources(1)%DX(3)/2, FormatString = 'f10.3')) // &
          ' to '                                                                                       // &
          Trim(Std2Char(Sources%Sources(1)%X(3) + Sources%Sources(1)%DX(3)/2, FormatString = 'f10.3')) // &
          Trim(ZCoord%Name)
      End If
      Write (Unit, *) Trim(Line)

    EndIf

    Line = Time2Char(                             &
             ShortTime2Time(EndTime - StartTime), &
             OutputOpts%Seconds,                  &
             OutputOpts%DecimalPlaces,            &
             .true.                               &
           )
    Write (Unit, *) 'Run duration:               ' // Trim(Line)

    ! Output grid information ! note iHgrid may vary if across $$
    If (FieldReq%iHGrid /= 0) Then
      HGrid => Grids%HGrids(FieldReq%iHGrid)
      If (.not. HGrid%Unstructured) Then
        Write (Unit, *) 'X grid origin:              ', HGrid%x0
        Write (Unit, *) 'Y grid origin:              ', HGrid%y0
      Else
        Write (Unit, *) 'X grid origin:              ', ' '
        Write (Unit, *) 'Y grid origin:              ', ' '
      End If
      Write (Unit, *) 'X grid size:                ', HGrid%nx
      Write (Unit, *) 'Y grid size:                ', HGrid%ny
      Write (Unit, *) 'X grid resolution:          ', HGrid%dx
      Write (Unit, *) 'Y grid resolution:          ', HGrid%dy
    Else
      Write (Unit, *) 'X grid origin:              ', 'No H Grid'
      Write (Unit, *) 'Y grid origin:              ', 'No H Grid'
      Write (Unit, *) 'X grid size:                ', 'No H Grid'
      Write (Unit, *) 'Y grid size:                ', 'No H Grid'
      Write (Unit, *) 'X grid resolution:          ', 'No H Grid'
      Write (Unit, *) 'Y grid resolution:          ', 'No H Grid'
    EndIf

    ! Number of fields in output file
    Write (Unit, *) 'Number of preliminary cols: ', nPrelimCols
    Write (Unit, *) 'Number of field cols:       ', nFieldCols

    ! Blank line
    Write (Unit, *)

    Write (Unit, *) 'Fields:'

    nLines = nLines + 20

  End If

End Subroutine FieldNumHeaders

!-------------------------------------------------------------------------------------------------------------

Subroutine FieldGraphHeaders(Unit, Grids, Sources, FieldReq, Results)
! Writes headers (i.e. axes etc) for graphical output of fields.

  Implicit None
  ! Argument List:
  Integer,         Intent(In) :: Unit     ! Output unit.
  Type(Grids_),    Intent(In) :: Grids    ! Collection of grids.
  Type(Sources_),  Intent(In) :: Sources  ! Collection of sources.
  Type(FieldReq_), Intent(In) :: FieldReq ! Field requirement.
  Type(Results_),  Intent(In) :: Results  ! Collection of results.
  ! Locals:
  Integer   :: iHGrid !
  Integer   :: iZGrid !
  Integer   :: iTGrid !
  Real(Std) :: YMax   !} Coord range for graph.
  Real(Std) :: XMax   !}

  iZGrid = FieldReq%iZGrid

  ! Set scales for graph.
  YMax = FieldReq%YScale
  XMax = FieldReq%XScale

  ! Draw axes.
  Call Plot(Unit, (/ 0.0, 0.5*XMax /), (/ 0.0, 0.0           /))
  Call Plot(Unit, (/ 0.0, 0.0      /), (/ 0.0, 1.05_Std*YMax /))

  ! BL top and well mixed profile.
  If (Grids%ZGrids(iZGrid)%nZ /= 109) Then
    Call Plot(Unit, (/          0.0, 0.1_Std*XMax /), (/ 1000.0_Std, 1000.0_Std /))
    Call Plot(Unit, (/ 0.1_Std*XMax, 0.1_Std*XMax /), (/        0.0, 1000.0_Std /))
  End If

End Subroutine FieldGraphHeaders

!-------------------------------------------------------------------------------------------------------------

Subroutine PdfHeaders(Reqs, OutputOpts, iCase, iG, Units, Results, T, Grids, Sources)
! Sets up a pdf file and writes out the pdf header.

  Implicit None
  ! Argument List:
  Type(Sources_),    Intent(In)    :: Sources !
  Type(Reqs_),       Intent(In)    :: Reqs    !
  Type(OutputOpts_), Intent(In)    :: OutputOpts
  Type(Results_),    Intent(InOut) :: Results !
  Type(Units_),      Intent(InOut) :: Units   !
  Type(Grids_),      Intent(In)    :: Grids   !
  Integer,           Intent(In)    :: iCase   ! Index of case.
  Integer,           Intent(In)    :: iG      ! Index of pdf output group.
  Type(Time_),       Intent(In)    :: T       ! Current time.
  ! Locals:
  Integer        :: i            ! Loop index.
  Logical        :: ScreenFlag   ! Indicates numerical screen output required.
  Logical        :: DiskFlag     ! Indicates numerical disk output required.
  Character(200) :: Title        ! Window title.
  Character(200) :: ScreenHeader ! Header for numerical screen output.
  Character(200) :: DiskHeader   ! Header for numerical disk output.
                                 ! $$ use MaxLineLength or a new parameter, not 200.
  Integer        :: iHGrid       ! Index of horizontal grid.
  Integer        :: iZGrid       ! Index of vertical grid.
  Integer        :: iTGrid       ! Index of temporal grid.
  Integer        :: nP           ! Size of pdf.
  Integer        :: iP           ! Index for pdf.
  Real(Std)      :: BaseVal      ! Base threshold values in auto-mode calculation.
  Character(MaxOutputLineLength) :: Line ! Character string for writing output.

  ScreenFlag   = .false.
  DiskFlag     = .false.
  ScreenHeader = ' '
  DiskHeader   = ' '

  Do i = 1, Reqs%nPdfs

    If (Reqs%PdfReqs(i)%iOutputGroup /= iG) Cycle

    iHGrid = Reqs%PdfReqs(i)%iHGrid
    iZGrid = Reqs%PdfReqs(i)%iZGrid
    iTGrid = Reqs%PdfReqs(i)%iTGrid
    nP     = Results%Pdfs(i)%PdfSize

    ! Numerical screen output.
    If (Reqs%PdfReqs(i)%Screen) Then
      If (.not.ScreenFlag) Then
        ScreenFlag = .true.
        Call GetNewUnit(Results%PdfScreenUnits(iG), Units)
        Title = 'Pdfs '                           // &
                Trim(Int2Char(iG, Justify = 'L')) // & !$$ Use group name here ?
                ': Case '                         // &
                (Int2Char(iCase, Justify = 'L'))
        If (iTGrid /= 0) ScreenHeader =             &
                Trim(ScreenHeader)               // &
                FormatChar('T', 19, .true., 'R')
        If (iHGrid /= 0) ScreenHeader =                   &
                Trim(ScreenHeader)                     // &
                Trim(FormatChar('X', 15, .true., 'R')) // &
                FormatChar('Y', 15, .true., 'R')
        If (iZGrid /= 0) ScreenHeader =             &
                Trim(ScreenHeader)               // &
                FormatChar('Z', 15, .true., 'R')
      End If
      If (Results%Pdfs(i)%FixedThresholds) Then
        Do iP = 1, nP
          ScreenHeader = Trim(ScreenHeader) //                                     &
                         Std2Char(Results%Pdfs(i)%Thresholds(iP), 15, .true., 'R')
        End Do
      Else
        Do iP = 1, nP
          BaseVal  = 10.0 ** (Float(iP - nP) / Float(AutoPdfResolutionPerDecade))
          ScreenHeader = Trim(ScreenHeader)                 // &
                         Std2Char(BaseVal, 15, .true., 'R')
        End Do
        ScreenHeader = Trim(ScreenHeader) //                         &
                       FormatChar('(Scale Factor)', 15, .true., 'R')
      End If
    End If

    ! Numerical disk output.
    If (Reqs%PdfReqs(i)%Disk) Then
      If (.not.DiskFlag) Then
        DiskFlag = .true.

        If (.not. (Reqs%PdfReqs(i)%SeparateFile == 'T')) Then
          Results%PdfFiles(iG) = Trim(OutputOpts%Folder)              // &
                                 Trim(Reqs%PdfReqs(i)%OutputGroup)    // &
                                 'PdfC'                               // &
                                 Trim(Int2Char(iCase, Justify = 'L')) // &
                                 '.txt'
        Else
          Results%PdfFiles(iG) = Trim(OutputOpts%Folder)              // &
                                 Trim(Reqs%PdfReqs(i)%OutputGroup)    // &
                                 'PdfC'                               // &
                                 Trim(Int2Char(iCase, Justify = 'L')) // &
                                 '_'                                  // &
                                 Trim(FileNameTime(T))                // & ! $$ adjust time zone
                                 '.txt'
        End If

        If (iTGrid /= 0) DiskHeader =               &
                Trim(DiskHeader)                 // &
                FormatChar('Date and Time', 19, .true., 'R')
        If (iHGrid /= 0) DiskHeader =                     &
                Trim(DiskHeader)                       // &
                Trim(FormatChar('X', 15, .true., 'R')) // &
                FormatChar('Y', 15, .true., 'R')
        If (iZGrid /= 0) DiskHeader =               &
                Trim(DiskHeader)                 // &
                FormatChar('Z', 15, .true., 'R')
      End If
      If (Results%Pdfs(i)%FixedThresholds) Then
        Do iP = 1, nP
          DiskHeader = Trim(DiskHeader)                                          // &
                       Std2Char(Results%Pdfs(i)%Thresholds(iP), 15, .true., 'R')
        End Do
      Else
        Do iP = 1, nP
          BaseVal    = 10.0 ** (Float(iP - nP) / Float(AutoPdfResolutionPerDecade))
          DiskHeader = Trim(DiskHeader)                   // &
                       Std2Char(BaseVal, 15, .true., 'R')
        End Do
        DiskHeader = Trim(DiskHeader) //                           &
                     FormatChar('(Scale Factor)', 15, .true., 'R')
      End If
    End If

  End Do

  ! Numerical screen output.
  If (ScreenFlag) Then
    Call OpenTextWindow(Results%PdfScreenUnits(iG), Trim(Title), nRows = 200)
    Write (Results%PdfScreenUnits(iG), *)              &
            ScreenHeader(1:Len_Trim(ScreenHeader) - 1)
  End If

  ! Numerical disk output.
  If (DiskFlag) Then
    Results%PdfDiskUnits(iG) = OpenFile(                      &
                                 Results%PdfFiles(iG),        &
                                 Units,                       &
                                 Status = 'Replace',          &
                                 RecL   = MaxOutputLineLength &
                               )

    ! Header block.
    Write (Results%PdfDiskUnits(iG), *) 'NAME III'
    Write (Results%PdfDiskUnits(iG), *) ModelVersion
    Line = 'Run time:            ' // Trim(Time2Char(Results%StartClockTime, .true., 3, .true., .true.))
    Write (Results%PdfDiskUnits(iG), *) Trim(Line)

    ! List of pdfs in output file.
    Write (Results%PdfDiskUnits(iG), *) 'Number of pdfs: ', 1 ! $$

    ! Blank line.
    Write (Results%PdfDiskUnits(iG), *) ' '

    ! Pdf header line.
    Write (Results%PdfDiskUnits(iG), *) DiskHeader(1:Len_Trim(DiskHeader) - 1)
  End If

End Subroutine PdfHeaders

!-------------------------------------------------------------------------------------------------------------

Subroutine PPInfoNumHeaders(        &
             Unit,                  &
             OutputOpts,            &
             Grids, Flows, Sources, &
             PPInfoReq,             &
             Results,               &
             nLines                 &
           )
! Writes headers for numerical output of sets of particle/puff information.

  Implicit None
  ! Argument List:
  Integer,           Intent(In)    :: Unit       ! Output unit.
  Type(OutputOpts_), Intent(In)    :: OutputOpts ! Output options.
  Type(Grids_),      Intent(In)    :: Grids      ! Collection of grids.
  Type(Flows_),      Intent(In)    :: Flows      ! Collection of flow module instance states.
  Type(Sources_),    Intent(In)    :: Sources    ! Collection of sources.
  Type(PPInfoReq_),  Intent(In)    :: PPInfoReq  ! Requirement for a set of
                                                 ! particle/puff information.
  Type(Results_),    Intent(InOut) :: Results    ! Collection of results.
  Integer,           Intent(InOut) :: nLines     ! Number of lines written to file.
  ! Locals:
  Character(MaxOutputLineLength) :: Line                ! Character string for writing output.
  Integer                        :: i                   ! Loop index variable.
  Integer                        :: nFlows
  Character(MaxCharLength)       :: FlowNames(MaxFlows)

  Write (Unit, *) ModelVersion
  Write (Unit, *) 'Run name:            ' // Trim(Results%RunName)
  Line = 'Run time:            ' // Trim(Time2Char(Results%StartClockTime, .true., 3, .true.))
  Write (Unit, *) Trim(Line)
  Call FlowList(Flows, nFlows, FlowNames)
  Line = ' '
  Do i = 1, nFlows
    If (i == 1) Then
      Line = FlowNames(i)
    Else
      Line = Trim(Line) // '; ' // FlowNames(i)
    End If
  End Do
  Write (Unit, *) 'Met data:            ', Trim(Line)
  If (IsBackwards()) Then
    Line = 'Backward trajectories'
  Else
    Line = 'Forward trajectories'
  End If
  Write (Unit, *) Trim(Line)
  Write (Unit, *)

  nLines = nLines + 6

  If (PPInfoReq%iSource /= 0) Then
    Line = 'Particle/puff information restricted to source ' // Trim(PPInfoReq%SourceName)
  Else
    Line = 'Particle/puff information for all sources'
  End If
  Write (Unit, *) Trim(Line)
  nLines = nLines + 1

  If (PPInfoReq%Puffs == 'A' .or. PPInfoReq%Puffs == 'S') Then
    If (PPInfoReq%Puffs == 'S') Then
      Line = 'Puff details for puff indices in the range ' // &
             Trim(Int2Char(PPInfoReq%FirstPuff))           // &
             ' to '                                        // &
             Trim(Int2Char(PPInfoReq%LastPuff))
    Else
      Line = 'Puff details for all puffs'
    End If
    Write (Unit, *) Trim(Line)
    nLines = nLines + 1
  End If
  If (PPInfoReq%Particles == 'A' .or. PPInfoReq%Particles == 'S') Then
    If (PPInfoReq%Particles == 'S') Then
      Line = 'Particle details for particle indices in the range ' // &
             Trim(Int2Char(PPInfoReq%FirstParticle))               // &
             ' to '                                                // &
             Trim(Int2Char(PPInfoReq%LastParticle))
    Else
      Line = 'Particle details for all particles'
    End If
    Write (Unit, *) Trim(Line)
    nLines = nLines + 1
  End If

  Write (Unit, *)
  nLines = nLines + 1

  Write (Unit, *) 'Set of particle/puff information: ' // Trim(PPInfoReq%Name)
  nLines = nLines + 1

End Subroutine PPInfoNumHeaders

!-------------------------------------------------------------------------------------------------------------

Subroutine PPInfoColumnHeaders( &
             Unit,              &
             OutputOpts,        &
             Coords, Specieses, &
             PPInfoReq,         &
             nLines             &
           )
! Writes column headers for numerical output of sets of particle/puff information.

  Implicit None
  ! Argument List:
  Integer,           Intent(In)    :: Unit       ! Output unit.
  Type(OutputOpts_), Intent(In)    :: OutputOpts ! Output options.
  Type(Coords_),     Intent(In)    :: Coords     ! Collection of coord systems.
  Type(Specieses_),  Intent(In)    :: Specieses  ! Collection of specieses.
  Type(PPInfoReq_),  Intent(In)    :: PPInfoReq  ! Requirement for a set of
                                                 ! particle/puff information.
  Integer,           Intent(InOut) :: nLines     ! Number of lines written to file.
  ! Locals:
  Character(MaxOutputLineLength) :: Line ! Character string for writing output.
  Integer                        :: i    ! Loop index.

  ! 1. Basic information.
  Line = Trim(FormatChar('Puff?',        15, .true., 'R')) // &
         Trim(FormatChar('PP Index',     15, .true., 'R')) // &
         Trim(FormatChar('Source',       15, .true., 'R')) // &
         Trim(FormatChar('Release Time', 25, .true., 'R')) // &
         Trim(FormatChar('Time',         25, .true., 'R')) // &
         Trim(FormatChar('Travel Time',  15, .true., 'R'))
  Line = Trim(Line) // FormatChar(                                   &
                         'X (' // Trim(PPInfoReq%HCoordName) // ')', &
                         15, .true., 'R'                             &
                       )
  Line = Trim(Line) // FormatChar(                                   &
                         'Y (' // Trim(PPInfoReq%HCoordName) // ')', &
                         15, .true., 'R'                             &
                       )
  Line = Trim(Line) // FormatChar(                                   &
                         'Z (' // Trim(PPInfoReq%ZCoordName) // ')', &
                         15, .true., 'R'                             &
                       )
  Line = Trim(Line)                                          // &
         Trim(FormatChar('U Turb',         15, .true., 'R')) // &
         Trim(FormatChar('V Turb',         15, .true., 'R')) // &
         Trim(FormatChar('W Turb',         15, .true., 'R')) // &
         Trim(FormatChar('Puff Spread XX', 15, .true., 'R')) // &
         Trim(FormatChar('Puff Spread YY', 15, .true., 'R')) // &
         Trim(FormatChar('Puff Spread ZZ', 15, .true., 'R'))

  ! 2. Met information.
  If (PPInfoReq%Met) Then
    Line = Trim(Line)                                          // &
           Trim(FormatChar('U Ambient',      15, .true., 'R')) // &
           Trim(FormatChar('V Ambient',      15, .true., 'R')) // &
           Trim(FormatChar('W Ambient',      15, .true., 'R')) // &
           Trim(FormatChar('Sigma UU',       15, .true., 'R')) // &
           Trim(FormatChar('Sigma VV',       15, .true., 'R')) // &
           Trim(FormatChar('Sigma WW',       15, .true., 'R')) // &
           Trim(FormatChar('Temperature',    15, .true., 'R')) // &
           Trim(FormatChar('Pressure',       15, .true., 'R')) // &
           Trim(FormatChar('Potential Temp', 15, .true., 'R')) // &
           Trim(FormatChar('BL Depth',       15, .true., 'R')) // &
           Trim(FormatChar('Cloud (oktas)',  15, .true., 'R')) // &
           Trim(FormatChar('Rel Humidity',   15, .true., 'R')) // &
           Trim(FormatChar('Wind Speed',     15, .true., 'R')) // &
           Trim(FormatChar('Wind Direction', 15, .true., 'R'))
  End If

  ! 3. Mass information
  If (PPInfoReq%Mass) Then
    Do i = 1, Specieses%nParticleSpecieses
      Line = Trim(Line)                                                                     // &
             FormatChar(                                                                       &
               Trim(Specieses%Specieses( Specieses%iParticle2Species(i) )%Name)             // &
               ' ('                                                                         // &
               Trim(Specieses%Specieses( Specieses%iParticle2Species(i) )%MaterialUnitName) // &
               ')',                                                                            &
               15, .true., 'R'                                                                 &
             )
    End Do
  End If

  ! 4. Plume rise information.
  If (PPInfoReq%PlumeRise) Then
    Line = Trim(Line)                                       // &
           Trim(FormatChar('FMass',       15, .true., 'R')) // &
           Trim(FormatChar('FMass_0',     15, .true., 'R')) // &
           Trim(FormatChar('FMomentum_X', 15, .true., 'R')) // &
           Trim(FormatChar('FMomentum_Y', 15, .true., 'R')) // &
           Trim(FormatChar('FMomentum_Z', 15, .true., 'R')) // &
           Trim(FormatChar('FHeat',       15, .true., 'R'))
  End If

  ! 5. Dispersion scheme information.
  If (PPInfoReq%DispersionScheme) Then
    Line = Trim(Line)                                          // &
           Trim(FormatChar('Skew Turb',      15, .true., 'R')) // &
           Trim(FormatChar('X-U RandomWalk', 15, .true., 'R')) // &
           Trim(FormatChar('Inhomg Turb',    15, .true., 'R'))
  End If

  ! 6. Puff family information.
  If (PPInfoReq%PuffFamily) Then
    Line = Trim(Line)                                         // &
           Trim(FormatChar('Original Puff', 15, .true., 'R')) // &
           Trim(FormatChar('N Splits',      15, .true., 'R')) // &
           Trim(FormatChar('Parent',        15, .true., 'R')) // &
           Trim(FormatChar('Sibling',       15, .true., 'R'))
  End If

  Write (Unit, *) Line(1:Len_Trim(Line) - 1)
  nLines = nLines + 1

End Subroutine PPInfoColumnHeaders

!-------------------------------------------------------------------------------------------------------------

End Module OutputModule

! Module:  Global Parameters Module

Module GlobalParametersModule

! This module provides code which defines various parameters.

! Module overview
! ---------------

! The module defines a number of parameters. Many of these are associated with particular modules and, if code
! modularity and reuse were the highest priority, would be better defined within those modules (possibly with
! the public attribute or with a routine provided to return their value in cases where modules which use the
! module in question may need to know the value). However it is more convenient to group them all together
! here to make it easy to change their values without knowledge of the rest of the code.

! Module use
! ----------

! Nothing needed other than a use statement.

! Module call tree
! ----------------

! None

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Public

!-------------------------------------------------------------------------------------------------------------

! Precision parameters (parameters for two types of real variable are defined to allow the possibility of
! altering the precision of horizontal coordinates independently from other real variables).
Integer, Parameter :: I16 = 2     ! Kind parameter for 16-bit integers.
Integer, Parameter :: I32 = 4     ! Kind parameter for 32-bit integers.
Integer, Parameter :: I64 = 8     ! Kind parameter for 64-bit integers.
Integer, Parameter :: P32 = 4     ! Kind parameter for 32-bit precision.
Integer, Parameter :: P64 = 8     ! Kind parameter for 64-bit precision.
#ifdef StdIs64Bit
  Integer, Parameter :: Std = P64 ! Kind parameter for standard reals.
#else
  Integer, Parameter :: Std = P32 ! Kind parameter for standard reals.
#endif
#ifdef PosIs64Bit
  Integer, Parameter :: Pos = P64 ! Kind parameter for position coords.
#else
  Integer, Parameter :: Pos = Std ! Kind parameter for position coords.
#endif

!-------------------------------------------------------------------------------------------------------------

! Lengths for various character strings.
Integer, Parameter :: MaxCharLength     =   64 ! Maximum length for various character strings.
Integer, Parameter :: MaxCharLength2    =  128 !} Maximum length for various character strings for which
Integer, Parameter :: MaxCharLength3    = 1024 !} MaxCharLength is too short.
Integer, Parameter :: MaxFileNameLength =  768 ! Maximum length for names of files and folders.

!-------------------------------------------------------------------------------------------------------------

! Model version.
!Character(MaxCharLength), Parameter :: ModelVersion = 'NAME III (development version 6.5+)'  ! Development.
!Character(MaxCharLength), Parameter :: ModelVersion = 'NAME III (version 6.5 - beta release)' ! Beta release.
Character(MaxCharLength), Parameter :: ModelVersion = 'NAME III (version 6.5)'               ! Main release.

!-------------------------------------------------------------------------------------------------------------

! Parameters for error and message module.
Integer, Parameter :: LogFileUnit       = 7  ! Unit number for the log file (must not be 0, 5, or 6).
Integer, Parameter :: ErrorFileUnit     = 8  ! Unit number for the error file (must not be 0, 5, or 6).
Integer, Parameter :: MessageLineLength = 79 ! Maximum length of message to be output on a given line.
!!!AJM!!!
Integer, Parameter :: ParticleEndLocationFileUnit       = 9  ! Unit number for the log file (must not be 0, 5, or 6).
!!!AJM!!!

!-------------------------------------------------------------------------------------------------------------

! Parameters for unit module.
Integer, Parameter :: FirstUnit = 20   ! Lowest unit number to be used for I/O (must be > 0).
Integer, Parameter :: LastUnit  = 1050 ! Highest unit number to be used for I/O.

!-------------------------------------------------------------------------------------------------------------

! Parameters for headed file module.
Integer, Parameter :: MaxFilesPerSet      = 36   ! Maximum number of headed files per set of headed files.
Integer, Parameter :: MaxBlockForms       = 26   ! Maximum number of block form descriptions that can be used
                                                 ! with a set of headed files in one pass.
Integer, Parameter :: MaxColumnKeys       = 95   ! Maximum number of column keywords that can be defined in
                                                 ! connection with a block in a headed file.
Integer, Parameter :: MaxColumnKeysTwoD   = 11   ! Maximum number of column keywords that can be defined in
                                                 ! connection with a block in a headed file which is read as a
                                                 ! single '2-D' block of data (rather than a line at a time).
                                                 ! Must be <= MaxColumnKeys.
Integer, Parameter :: MaxSetsOfColumnKeys = 3    ! Maximum number of sets of column keywords that can be
                                                 ! defined in connection with a block in a headed file.
Integer, Parameter :: MaxColumns          = 60   ! Maximum number of columns in a block in a headed file.
Integer, Parameter :: MaxLinesTwoD        = 2048 ! Maximum number of lines of data items in a block in a
                                                 ! headed file which is read as a single '2-D' block of data
                                                 ! (rather than a line at a time).
Integer, Parameter :: MaxLineLength       = 2048 ! Maximum line length in headed files. Must be (strictly) >
                                                 ! MaxTokenLength.
Integer, Parameter :: MaxTokenLength      = 300  ! Maximum length for block keywords, block names, column
                                                 ! keywords and data items in headed files.
Integer, Parameter :: MaxArrays           = 12   !} Maximum number and length of input arrays.
Integer, Parameter :: MaxArrayLength      = 85   !}

!-------------------------------------------------------------------------------------------------------------

! Parameters for error and message II module.
Integer, Parameter :: MaxMessageControls = 60 ! Maximum number message controls.

!-------------------------------------------------------------------------------------------------------------

! Parameters for coordinate system module.
Integer, Parameter :: MaxEtaDefns    = 6  ! Maximum number of eta definitions.
Integer, Parameter :: MaxHCoords     = 10 ! Maximum number of horizontal coord systems.
Integer, Parameter :: MaxZCoords     = 10 ! Maximum number of vertical coord systems.

!-------------------------------------------------------------------------------------------------------------

! Parameters for grid and domain module.
Integer, Parameter :: MaxLocationses = 5   ! Maximum number of sets of locations.
Integer, Parameter :: MaxHGrids      = 5000  ! Maximum number of horizontal grids.
Integer, Parameter :: MaxZGrids      = 5000  ! Maximum number of vertical grids.
Integer, Parameter :: MaxTGrids      = 5000 ! Maximum number of temporal grids.
Integer, Parameter :: MaxDGrids      = 5   ! Maximum number of data grids.
Integer, Parameter :: MaxDomains     = 32  ! Maximum number of domains.

!-------------------------------------------------------------------------------------------------------------

! Parameters for met and flow modules.
Integer, Parameter :: MaxTurbLayers          = 2
Integer, Parameter :: MaxFieldQualifiers     = 4
Integer, Parameter :: MaxNWPMetFields        = 90
Integer, Parameter :: MaxRadarMetFields      = 10
Integer, Parameter :: MaxAncillaryMetFields  = 10
Integer, Parameter :: MaxAncillaryMetsPerNWPFlow = 6
Integer, Parameter :: MaxNWPMetDefns         = 32
Integer, Parameter :: MaxNWPMetDefn2s        = 32
Integer, Parameter :: MaxRadarMetDefns       = 2
Integer, Parameter :: MaxRadarMetDefn2s      = 2
Integer, Parameter :: MaxAncillaryMetDefns   = 8
Integer, Parameter :: MaxAncillaryMetDefn2s  = 8
Integer, Parameter :: MaxMetMods             = 5
Integer, Parameter :: MaxPrototypeMets       = 3
Integer, Parameter :: MaxSingleSiteMets      = 1
Integer, Parameter :: MaxNWPMets             = 32
Integer, Parameter :: MaxRadarMets           = 2
Integer, Parameter :: MaxAncillaryMets       = 32
Integer, Parameter :: MaxMetsPerMod          = Max(                 &
                                                 MaxPrototypeMets,  &
                                                 MaxSingleSiteMets, &
                                                 MaxNWPMets,        &
                                                 MaxRadarMets,      &
                                                 MaxAncillaryMets   &
                                               )
Integer, Parameter :: MaxMets                = MaxPrototypeMets  + &
                                               MaxSingleSiteMets + &
                                               MaxNWPMets        + &
                                               MaxRadarMets      + &
                                               MaxAncillaryMets
Integer, Parameter :: MaxMetAttribs          = 6
Integer, Parameter :: MaxFlowMods            = 4
Integer, Parameter :: MaxPrototypeFlows      = 3
Integer, Parameter :: MaxSingleSiteFlows     = 2
Integer, Parameter :: MaxNWPFlows            = 20
Integer, Parameter :: MaxBuildingFlows       = 2
Integer, Parameter :: MaxRadarFlows          = 2
Integer, Parameter :: MaxLincomFlows         = 2
Integer, Parameter :: MaxFlowsPerMod         = Max(                  &
                                                 MaxPrototypeFlows,  &
                                                 MaxSingleSiteFlows, &
                                                 MaxNWPFlows,        &
                                                 MaxBuildingFlows,   &
                                                 MaxLincomFlows,     &
                                                 MaxRadarFlows       &
                                               )
Integer, Parameter :: MaxFlows               = MaxPrototypeFlows  + &
                                               MaxSingleSiteFlows + &
                                               MaxNWPFlows        + &
                                               MaxBuildingFlows   + &
                                               MaxLincomFlows     + &
                                               MaxRadarFlows
Integer, Parameter :: MaxFlowOrders          = 5
Integer, Parameter :: MaxFlowSubsets         = 10
Integer, Parameter :: MaxFlowAttribs         = 8
! MaxTurbLayers              :: Maximum number of turbulent layers allowed when representing the atmosphere in
!                               homogeneous turbulent layers.
! MaxNWPMetFields            :: Maximum number of NWP met fields in any met file structure.
! MaxAncillaryMetFields      :: Maximum number of Ancillary met fields in any met file structure.
! MaxAncillaryMetsPerNWPFlow ::
! MaxNWPMetDefns             :: Maximum number of NWP met definitions (part 1 - basic information).
! MaxNWPMetDefn2s            :: Maximum number of NWP met definitions (part 2 - met file structure).
! MaxFieldQualifiers         :: Maximum number of qualifiers per field in NWP met definitions (part 2).
! MaxAncillaryMetDefns       :: Maximum number of Ancillary met definitions (part 1 - basic information).
! MaxAncillaryMetDefn2s      :: Maximum number of Ancillary met definitions (part 2 - met file structure).
! MaxMetMods                 :: Maximum number of met modules.
! MaxPrototypeMets           :: Maximum number of Prototype Met module instances.
! MaxSingleSiteMets          :: Maximum number of Single Site Met module instances.
! MaxNWPMets                 :: Maximum number of NWP Met module instances.
! MaxRadarMets               :: Maximum number of Radar Met module instances.
! MaxAncillaryMets           :: Maximum number of Ancillary Met module instances.
! MaxMetsPerMod              :: Maximum number of instances of any given met module. Must be the largest of
!                               the maximum number of instances for each of the various met modules.
! MaxMets                    :: Maximum number of met module instances. Must be the sum of the maximum number
!                               of instances for the various met modules.
! MaxMetAttribs              :: Maximum number of met module attributes.
! MaxFlowMods                :: Maximum number of flow modules.
! MaxPrototypeFlows          :: Maximum number of Prototype Flow module instances.
! MaxSingleSiteFlows         :: Maximum number of Single Site Flow module instances.
! MaxNWPFlows                :: Maximum number of NWP Flow module instances.
! MaxBuildingFlows           :: Maximum number of Building Flow module instances.
! MaxRadarFlows              :: Maximum number of Radar Flow module instances.
! MaxLincomFlows             :: Maximum number of Lincom Flow module instances.
! MaxFlowsPerMod             :: Maximum number of instances of any given flow module. Must be the largest of
!                               the maximum number of instances for each of the various flow modules.
! MaxFlows                   :: Maximum number of flow module instances. Must be the sum of the maximum number
!                               of instances for the various flow modules.
! MaxFlowOrders              :: Maximum number of ordered subsets of flow module instances.
! MaxFlowSubsets             :: Maximum number of subsets of flow module instances.
! MaxFlowAttribs             :: Maximum number of flow module attributes.

!-------------------------------------------------------------------------------------------------------------

! Parameters for particle size distribution module.
Integer, Parameter :: MaxSizeDists      =  5 ! Maximum number of particle size distributions.
Integer, Parameter :: MaxDiameterRanges = 10 ! Maximum number of diameter ranges per particle size 
                                             ! distribution.
                                                      
!-------------------------------------------------------------------------------------------------------------
                                                      
! Parameters for species module.
Integer, Parameter :: MaxSpecieses        = 52 ! Maximum number of specieses.
Integer, Parameter :: MaxDecayChainLength = 15 ! Maximum length of radiological decay chains.
Integer, Parameter :: MaxEnergies         = 20 ! Maximum number of photon energies for every species.
Integer, Parameter :: MaxCGSpecies        = 16 ! Maximum number of species for consideration in a cloud gamma
                                               ! calculation - also the max no. of cloudgammaparamses.

!-------------------------------------------------------------------------------------------------------------

! Parameters for source module.
Integer, Parameter :: MaxSourceGroups          = 100  ! Maximum number of source groups.
Integer, Parameter :: MaxSourceGroupsPerSource = 4    ! Maximum number of source groups per source.
Integer, Parameter :: MaxTimeDeps              = 5    ! Maximum number of source time dependencies.
Integer, Parameter :: MaxTimeDepFactors        = 2048 ! Maximum number of factors per source time dependency.
Integer, Parameter :: MaxTimeDepPairs          = 5    ! Maximum number of wild-card time pairs per source time
                                                      ! dependency factor.

!-------------------------------------------------------------------------------------------------------------

! Parameters for output module.
Integer, Parameter :: MaxProcesses            = 100
Integer, Parameter :: MaxProcessesPerFieldReq = 3
Integer, Parameter :: MaxPdfReqs              = 3
Integer, Parameter :: MaxPPInfoReqs           = 10
Integer, Parameter :: MaxFieldOutputFiles     = 400
Integer, Parameter :: MaxPPInfoOutputFiles    = 1000
! Note compaq doesn't like long lines written to screen - even if trimmed to reasonable
! length causes unexpected out of memory failure. $$
#ifdef CompaqPCCompiler
  Integer, Parameter :: MaxOutputLineLength   = 1000
#else
  Integer, Parameter :: MaxOutputLineLength   = 100000
#endif
Integer, Parameter :: OutputPrelimColumnWidth = 24
Integer, Parameter :: OutputFieldColumnWidth  = 26
!Integer, Parameter :: MaxFieldReqsPerReqInfo  = 600
Integer, Parameter :: MaxFieldReqsPerReqInfo  = 5000
Integer, Parameter :: MaxPPInfoReqsPerReqInfo = 10
Integer, Parameter :: MaxReqInfoTimes         = 45
Integer, Parameter :: MaxReqInfoTravelTimes   = 45
! MaxProcesses            :: Maximum number of processing steps.
! MaxProcessesPerFieldReq :: Maximum number of processing steps per field requirements.
! MaxPdfReqs              :: Maximum number of pdfs requirements.
! MaxPPInfoReqs           :: Maximum number of requirements for a set of particle/puff information.
! MaxFieldOutputFiles     :: Maximum number of open output files/windows for fields output which may need to
!                            be written to again (per output group and for each of screen and disk output
!                            considered separately).
! MaxPPInfoOutputFiles    :: Maximum number of open output files/windows for particle/puff information output
!                            which may need to be written to again (for each of screen and disk output
!                            considered separately).
! MaxOutputLineLength     :: Maximum length for output lines.
! OutputPrelimColumnWidth :: Output preliminary column width (excluding separating commas).
! OutputFieldColumnWidth  :: Output field column width (excluding separating commas).
! MaxFieldReqsPerReqInfo  :: Maximum number of field requirements in ReqInfo type.
! MaxPPInfoReqsPerReqInfo :: Maximum number of requirements for a set of particle/puff information in ReqInfo
!                            type.
! MaxReqInfoTimes         :: Maximum number of times in ReqInfo array.
! MaxReqInfoTravelTimes   :: Maximum number of travel-times in ReqInfo array.

!-------------------------------------------------------------------------------------------------------------

! Parameters for fluctuation module. $$ should be under output module?
Integer, Parameter :: MaxPdfSize = 40 ! Maximum number of points in any single probability distribution (i.e.
                                      ! percentiles or exceedence probabilities).

!-------------------------------------------------------------------------------------------------------------

! Parameters for chemistry module.
Integer,   Parameter :: MaxSTOCHEMnX = 72 !] Maximum grid size of STOCHEM fields.
Integer,   Parameter :: MaxSTOCHEMnY = 36 !]
Integer,   Parameter :: MaxSTOCHEMnZ = 9  ! Maximum number of STOCHEM vertical levels.

!-------------------------------------------------------------------------------------------------------------

! Parameters for case module.
Integer, Parameter :: MaxDispOptses = 6 ! Maximum number of sets of dispersion model options.

!-------------------------------------------------------------------------------------------------------------

! Default values used in fluctuations scheme. $$ revisit when stat output stuff written.
! ## Dave - not certain about setting up the default percentiles in this way
! ## using an array as a parameter - better to initialise values in a function?
! ## This would also allow us to check consistency between the array length
! ## and the specified value of AutoPdfNumPercentiles. Any views ??
Integer, Parameter   :: AutoPdfNumProbabilities    = 20
Integer, Parameter   :: AutoPdfResolutionPerDecade = 4
Integer, Parameter   :: AutoPdfNumPercentiles      = 13
Real(Std), Parameter :: AutoPdfPercentiles(AutoPdfNumPercentiles) = (/ &
                          0.01_Std,  &
                          0.10_Std,  &
                          1.00_Std,  &
                          5.00_Std,  &
                         10.00_Std,  &
                         25.00_Std,  &
                         50.00_Std,  &
                         75.00_Std,  &
                         90.00_Std,  &
                         95.00_Std,  &
                         99.00_Std,  &
                         99.90_Std,  &
                         99.99_Std  /)
Real(Std), Parameter :: MaxAvTime = 1800.0_Std
! AutoPdfNumProbabilities    :: Total number of exceedence probabilities calculated in auto-pdf mode.
! AutoPdfResolutionPerDecade :: Number of exceedence probabilities per decade calculated in auto-pdf mode.
! AutoPdfNumPercentiles      :: Number of percentiles calculated in auto-pdf mode.
! AutoPdfPercentiles         :: Default percentile values calculated in auto-pdf mode.
! MaxAvTime                  :: Maximum averaging time permitted in the fluctuations scheme (relevant for case
!                               of continuous sources only).

!-------------------------------------------------------------------------------------------------------------

! Parameters for input module.
Integer, Parameter :: MaxArgLength = 120 ! Max length of command line arguments.

!-------------------------------------------------------------------------------------------------------------

End Module GlobalParametersModule

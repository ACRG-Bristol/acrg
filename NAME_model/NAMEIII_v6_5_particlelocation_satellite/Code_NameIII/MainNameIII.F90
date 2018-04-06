! Module:  Main Name III Module

!-------------------------------------------------------------------------------------------------------------

Program NameIII
! Main routine for the NAME III model.

  Use TimerModule
  Use ServiceModule
  Use FlowsModule, Only: FlowOrder_, FlowSubset_, FlowAttrib_, Flows_,                              &
                         SetUpCoordsEtc_Flows, SetUpFlows_MetsCoordsEtc, SetUpFlowOrders_iFlows,    &
                         SetUpFlowSubsets_iFlows, SetUpFlows_iFlowAttribs, SetUpFlows_iFlowSubsets, &
                         SetUpFlows_iMets, SetUpFlows_iCoordsEtc, CheckFlows_FlowUpdateSubsets,     &
                         ResetMetsFlows,                                                            &
                         PrototypeFlow_, SingleSiteFlow_, NWPFlow_,                                 &
                         BuildingFlow_, RadarFlow_, LINCOMFlow_,                                    &
                         Mets_,                                                                     &
                         SetUpCoordsEtc_Mets, SetUpMets_CoordsEtc, SetUpMets_iCoordsEtc,            &
                         PrototypeMet_, SingleSiteMet_, NWPMet_, RadarMet_, AncillaryMet_,          &
                         CreateFlowMetLocks, DestroyFlowMetLocks
  !$ Use FlowsModule, Only: FlowLocks, MetLocks 
  Use SizeDistModule
  Use SpeciesModule
  Use SourceModule
  Use OutputModule, Only: OutputOpts_, Reqs_, Results_,                                                     &
                          CheckOutputOpts, SetUpReqs, SetUpiReqs, InitResults, RestartAdjustmentForResults, &
                          OutputModuleTimerSummary, OutputModuleTimerInitialise
  Use ParticleModule
  Use PuffModule
  Use ChemistryModule
  Use CaseModule
  Use RestartModule
  Use InputModule
  Use IOThreadModule,  Only: IOThreadModuleTimerSummary, IOThreadModuleTimerInitialise
  Use NWPMetModule,    Only: NWPMetModuleTimerSummary, NWPMetModuleTimerInitialise, EulerianModelGlobal

  Use OpenMPModule, Only : FinaliseOpenMPConfiguration   &
                          ,CheckOpenMPOpts               &
                          ,OpenMPOpts_
  Use CompilerInfoModule
  Use EulerianInterfaceModule, Only: SetUpEulerianGrid

  Implicit None
  ! Locals:
  Real(Std)                    :: BlendingHeights(MaxFlows,MaxFlows) ! $$
  Type(MainOpts_)              :: MainOpts                           !
  Type(RestartOpts_)           :: RestartOpts                        !
  Type(MultiCaseOpts_)         :: MultiCaseOpts                      !
  Type(ChemOpts_)              :: ChemOpts                           !
  Type(OutputOpts_)            :: OutputOpts                         !
  Type(EtaDefns_)              :: EtaDefns                           !
  Type(Coords_)                :: Coords                             !
  Type(Locationses_)           :: Locationses                        !
  Type(Grids_)                 :: Grids                              !
  Type(Domains_)               :: Domains                            !
  Type(Mets_)                  :: Mets                               !
  Type(Flows_)                 :: Flows                              !
  Type(Specieses_)             :: Specieses                          !
  Type(CloudGammaParamses_)    :: CloudGammaParamses                 ! A collection of sets of cloud
                                                                     ! gamma parameters.
  Type(MassLimits_)            :: MassLimits                         !
  Type(Sources_)               :: Sources                            !
  Type(SizeDists_)             :: SizeDists                          !
  Type(TimeDeps_)              :: TimeDeps                           !
  Type(Reqs_)                  :: Reqs                               !
  Type(DispOptses_)            :: DispOptses                         !
  Type(Units_)                 :: Units                              !
  Type(Results_)               :: Results                            !
  Type(DispState_)             :: DispState                          !
  Type(ChemistryDefn_)         :: ChemistryDefn                      !
  Type(ChemistryState_)        :: ChemistryState                     !
  Type(MaterialUnits_)         :: MaterialUnits                      !
  Type(OpenMPOpts_)            :: OpenMPOpts                         !
  Type(TimerOpts_)             :: TimerOpts                          !
  Logical                      :: EndOfCase                          !
  Logical                      :: EndOfCases                         !
  Integer                      :: iCase                              !
  Character(MaxFileNameLength) :: MainInputFile                      !
  Character(MaxFileNameLength) :: LogFile                            !
  Character(MaxFileNameLength) :: ErrorFile                          !
  Integer                      :: ClosePromptLevel                   !
  Character(MaxFileNameLength) :: RestartFileStem                    !
  Logical                      :: Restart                            ! Run will restart from a previous
                                                                     ! restart file.
  Character(MaxFileNameLength) :: RestartFileSuffix
  Logical                      :: NewSourcesDetected                 ! Indicates that one or more new sources
                                                                     ! have been detected in the input file
                                                                     ! during a restart run.
  Integer                      :: iFirstNewSource                    ! Index of first new source in a modified
                                                                     ! input file for a restart run.
  Character(MaxCharLength)     :: CharRunForTime    !
  Character(MaxCharLength)     :: CharRunToTime     ! Blank = no run-to limit
  Type(Time_)                  :: RunForTime
  Type(Time_)                  :: RunToTime
  Integer                      :: RunForCase
  Integer                      :: RunToCase
  Type(Time_)                  :: ClockTime                          !
  Logical                      :: Suspend  !
  Integer                      :: i
  Logical :: Exist1, Exist2
  Type(Timer_)                 :: TotalTimeTimer

!-------------------------------------------------------------------------------------------------------------

  ! Initialise error and message module.
  Call InitErrorAndMessageModule

  ! Set up unit 6 window.
  Call SetTitleTextLinesUnit6('Message window', 302)

  ! Process command line arguments.
  Call ProcessCommandLine(           &
         MainInputFile,              &
         LogFile, ErrorFile,         &
         ClosePromptLevel,           &
         RestartFileStem,            &
         Restart, RestartFileSuffix, &
         CharRunForTime,             &
         CharRunToTime               &
       )

  ! Pre-initialise the restart mechanism.
  Call PreInitRestart(RestartFileStem, Restart)

  ! Set up error and message module. After this call any closing message box will
  ! include the model version, messages are written to the log file (as well as unit
  ! 6), and ClosePromptLevel takes effect.
  Call SetUpErrorAndMessageModule(ModelVersion, LogFile, ErrorFile, ClosePromptLevel, Restart)

  ! Start time of model run.
  ClockTime = CurrentClockTime()

  ! Model version message.
  Call Message(' ')
  Call Message(ModelVersion)

  ! Compiler version message.
  Call Message(' ')
  Call Message('Compile time: '//Trim(CompileTime))

  ! Compiler information message.
  Call Message(' ')
  Call Message('Compiler version: '//Trim(CompilerVersion))

  ! Start time message.
  If (Restart) Then
    Call Message(' ')
    Call Message('Run restarted at ' // Trim(Time2Char(ClockTime, .true., 3, .true.)))
  Else
    Call Message(' ')
    Call Message('Run started at ' // Trim(Time2Char(ClockTime, .true., 3, .true.)))
  End If

  ! Log file message.
  Call Message(' ')
  Call Message('The log file for this run is ' // Trim(LogFile) // '.')

  ! Restart file message.
  Call Message(' ')
  Call Message(                                                   &
         'The restart files (if any) for this run have names ' // &
         'that begin with '                                    // &
         Trim(RestartFileStem)                                 // &
         'Restart.'                                               &
       )

  ! Initialise Units.
  Units = InitUnits((/ 11, 71, 100 /)) ! $$ reserve unit codes which are hard wired.
                                       ! 11 in building
                                       ! 71 in Radar Met
                                       ! 100 in Restart (for restart files)

  ! Create Locks in Flow Module
  Call CreateFlowMetLocks()

  ! Initialise Material Units
  Call MaterialUnitsInitialise(MaterialUnits)

  !-------------------------------!
  ! Read and process input files. !
  !-------------------------------!

  Call Message(' ')
  Call Message('Reading and checking input file(s)')

  ! Read the main set of headed input files.
  Call ReadInputFiles(                                               &
         MainInputFile,                                              &
         MainOpts, RestartOpts, MultiCaseOpts, ChemOpts, OutputOpts, &
         OpenMPOpts, TimerOpts,                                      &
         EtaDefns, Coords, Locationses, Grids, Domains,              &
         Mets, Flows,                                                &
         Specieses, CloudGammaParamses, MassLimits,                  &
         Sources, SizeDists, TimeDeps,                               &
         Reqs,                                                       &
         DispOptses,                                                 &
         BlendingHeights,                                            &
         MaterialUnits,                                              &
         Units                                                       &
       )

  ! Message controls.
  Call AddValidNameToMessageControls(                     &
         ValidName        = 'Inconsistent validity time', &
         CheckUniqueNames = .true.,                       &
         MessageControls  = GlobalMessageControls         &
       )
  Call AddValidNameToMessageControls(                  &
         ValidName        = 'Inconsistent NWP header', &
         CheckUniqueNames = .true.,                    &
         MessageControls  = GlobalMessageControls      &
       )
  Call AddValidNameToMessageControls(             &
         ValidName        = 'Missing NWP field',  &
         CheckUniqueNames = .true.,               &
         MessageControls  = GlobalMessageControls &
       )
  Call AddValidNameToMessageControls(                           &
         ValidName        = 'Inconsistent radar validity time', &
         CheckUniqueNames = .true.,                             &
         MessageControls  = GlobalMessageControls               &
       )
  Call AddValidNameToMessageControls(                        &
         ValidName        = 'Inconsistent radar met header', &
         CheckUniqueNames = .true.,                          &
         MessageControls  = GlobalMessageControls            &
       )
  Call AddValidNameToMessageControls(                  &
         ValidName        = 'Missing radar met field', &
         CheckUniqueNames = .true.,                    &
         MessageControls  = GlobalMessageControls      &
       )
  Call AddValidNameToMessageControls(             &
         ValidName        = 'No flow',            &
         CheckUniqueNames = .true.,               &
         MessageControls  = GlobalMessageControls &
       )
  Call AddValidNameToMessageControls(                    &
         ValidName        = 'No flow for particle/puff', &
         CheckUniqueNames = .true.,                      &
         MessageControls  = GlobalMessageControls        &
       )
  Call AddValidNameToMessageControls(                               &
         ValidName        = 'No time zone in single site met file', &
         CheckUniqueNames = .true.,                                 &
         MessageControls  = GlobalMessageControls                   &
       )
  Call AddValidNameToMessageControls(               &
         ValidName        = 'Zero source strength', &
         CheckUniqueNames = .true.,                 &
         MessageControls  = GlobalMessageControls   &
       )
  Call SetupMessageControls(CheckNamesValid = .true., MessageControls = GlobalMessageControls)


  ! Process some command line variables:
  If (CharRunForTime == ' ') Then
    RunForTime = InfFutureTime(Interval = .true.)
    RunForCase = Huge(RunForCase) / 2
  Else If (Verify(CharRunForTime, '1234567890') /= 0) Then
    RunForTime = Char2Time(CharRunForTime, Interval = .true.)
    RunForCase = 1
  Else
    RunForTime = InfFutureTime(Interval = .true.)
    RunForCase = Char2Int(CharRunForTime)
  End If
  If (CharRunToTime == ' ') Then
    RunToTime = InfFutureTime()
    RunToCase = Huge(RunForCase) / 2
  Else If (Verify(CharRunToTime, '1234567890') /= 0) Then
    RunToTime  = Char2Time(CharRunToTime)
    RunToCase  = Huge(RunForCase) / 2
    RunForCase = 1
  Else
    RunToTime = InfFutureTime()
    RunToCase = Char2Int(CharRunToTime)
  End If

  ! Check the main, dispersion model, output and OpenMP options.
  Call CheckMainOpts     (MainOpts)
  Call CheckRestartOpts  (Restart, RestartOpts)
  Call CheckMultiCaseOpts(MultiCaseOpts, DispOptses)
  Call CheckOutputOpts   (OutputOpts)
  Call CheckOpenMPOpts   (OpenMPOpts)
  Call CheckTimerOpts    (TimerOpts)

  Call TimerCreate(TotalTimeTimer,"TotalTime",TimerOpts)
  Call TimerOn(TotalTimeTimer)

! Initialise timers in the individual modules
  Call CaseModuleTimerInitialise(TimerOpts)
  Call OutputModuleTimerInitialise(TimerOpts)
  Call NWPMetModuleTimerInitialise(TimerOpts, OpenMPOpts)
  Call IOThreadModuleTimerInitialise(TimerOpts)
  Call ChemistryModuleTimerInitialise(TimerOpts)

  ! Set up RestartOpts.
  Call SetUpRestartOpts(RestartFileStem, RestartOpts)

  ! Set up collection of requirements, including adding any coord systems or grids required.
  Call SetUpReqs(Sources, Coords, Grids, Specieses, Reqs, MaterialUnits)

  ! Set up Coords and Grids by adding any extra coords and grids which Mets wants to define.
  Call SetUpCoordsEtc_Mets(Mets, Coords, Grids)

  ! $$ this is called twice as existing grids need to be finalised before SetUpCoordsEtc_Flows,
  ! while the latter may add more grids.
  Call SetUpGrids_CoordsEtc(Locationses, Coords, Grids)

  ! Set up Coords and Grids by adding any extra coords and grids which Flows wants to define.
  Call SetUpCoordsEtc_Flows(Flows, Coords, Grids)

  Call SetUpGrids_CoordsEtc(Locationses, Coords, Grids)

  Call SetUpDomains_CoordsEtc(Locationses, Coords, Domains)

  ! Set up Mets using information from EtaDefns, Coords and Grids.
  Call SetUpMets_CoordsEtc(EtaDefns, Coords, Grids, MultiCaseOpts%MetEnsembleSize, Mets)

  ! Set up Flows using information from EtaDefns, Coords, Grids and Mets.
  Call SetUpFlows_MetsCoordsEtc(EtaDefns, Coords, Grids, Mets, Flows)

  ! Set up Specieses.
  Call SetUpSpecieses(SizeDists, Specieses)

  ! Set up Specieses using information from MassLimits.
  Call SetUpSpecieses_MassLimits(MassLimits, Specieses)

  ! Set up Sources.
  Call SetUpSources(Grids, SizeDists, Specieses, MaterialUnits, Sources)

  ! Add Lat Long coord system to Coords
  Call AddHCoord(HCoord_LatLong(), Coords)

  ! Set up indices in various data types for referring to items in other data types. SetUpFlowOrders_iFlows
  ! and SetUpFlowSubsets_iFlows also set up Flows%Orders%Included and Flows%Subsets%Included respectively
  ! while SetUpFlows_iFlowAttribs also checks that the flow module instances have the required attributes.
  Call SetUpGrids_iCoords              (Coords, Grids)
  Call SetUpDomains_iCoords            (Coords, Domains)
  Call SetUpFlowOrders_iFlows          (Flows)
  Call SetUpFlowSubsets_iFlows         (Flows)
  Call SetUpFlows_iFlowAttribs         (Flows)
  Call SetUpFlows_iFlowSubsets         (Flows)
  Call SetUpMets_iCoordsEtc            (Coords, Grids, Mets)
  Call SetUpFlows_iMets                (Mets, Flows)
  Call SetUpFlows_iCoordsEtc           (Coords, Grids, Domains, Flows)
  Call SetUpiSpecieses                 (SizeDists, Specieses)
  Call SetUpiReqs                      (Coords, Grids, Specieses, Sources, SizeDists, Reqs, MaterialUnits)
  Call SetUpSpecieses_iCloudGammaParams(Specieses, CloudGammaParamses)
  Call SetUpSpecieses_iMaterialUnits   (Specieses, MaterialUnits)
  Call SetUpiSources                   (Coords, Grids, Specieses, MaterialUnits, SizeDists, TimeDeps, Sources)

  ! Sets up full decay chains in Specieses.
  Call SetUpSpecieses_DecayChains(Specieses)    

  ! Check that the update subsets used within the flow module instances only include
  ! flow module instances which are higher up the update priority order.
  Call CheckFlows_FlowUpdateSubsets(Flows)

  Call Message('Input file(s) successfully read and checked')

  !----------------------------------------------------------------!
  ! Initialise various variables which will evolve during the run. !
  !----------------------------------------------------------------!

  ! Initialise a collection of results.
  Results = InitResults(MainOpts%RunName, ClockTime, Reqs)

  ! Initialise a set of information on the state of the dispersion calculation.
  Call InitDispState(DispState)

  ! Initialise state of the chemistry.
  Call PreInitChemistryState(ChemistryState)

  ! Initialise Eulerian model
  If (Specieses%nFields > 0 .or. Any(DispOptses%DispOptses(1:DispOptses%nDispOptses)%Chemistry)) Then
    Call SetUpEulerianGrid(DispState%EulerianField, & 
                           Coords,                  &  
                           Grids,                   &
                           Specieses)
  End If
  ! flag for converting NWP w to eta dot.
  EulerianModelGlobal = Specieses%AdvectedFields

  ! Initialise random number generator.
  Call SetRandomState(MainOpts%RandomMode, GlobalRandomState)

  ! Initialise iCase, EndOfCases and EndOfCase (EndOfCase is set to true to trigger
  ! next case).
  iCase      = 0
  EndOfCases = .false.
  EndOfCase  = .true.

  !------------------------!
  ! Read the restart file. !
  !------------------------!

  If (RestartOpts%RestartType /= 0) Then

    ! Note InitRestart can switch Restart off.
    Call InitRestart(RestartOpts, Restart, Units)

    If (Restart) Then

      ! Note ReadRestartFile can switch Restart off.
      Call ReadRestartFile(           &
             RestartFileSuffix,       &
             RestartOpts, OutputOpts, &
             Sources%nSources,        &
             NewSourcesDetected,      &
             iFirstNewSource,         &
             iCase,                   &
             Results,                 &
             DispState,               &
             ChemistryState,          &
             EndOfCase, EndOfCases,   &
             Restart,                 &
             Units                    &
           )

      Call RestartAdjustmentForInputChanges(      &
             Grids, Domains,                      &
             Sources, SizeDists,                  &
             Reqs,                                &
             NewSourcesDetected, iFirstNewSource, &
             DispState                            &
           )

      Call RestartAdjustmentForResults(Reqs, Results)

    End If

  End If
  
!!! AJM
  !------------------!
  ! Open File to write out particle location when particle leaves domain !
  !------------------!
  WRITE(6,*)'Open ParticleEndLocation File unit:',ParticleEndLocationFileUnit
  OPEN(UNIT=ParticleEndLocationFileUnit,FILE='particle_location.txt',Status = 'Unknown')
  WRITE(ParticleEndLocationFileUnit,'(2A7,A6,A8,A3)')' Long ',' Lat ',' Ht ',' Age(hr)',' Id'
!!! AJM

  !------------------!
  ! Loop over cases. !
  !------------------!

  Suspend = .false.
  If (EndOfCase) Then
    RunToCase = Min(RunToCase, iCase + RunForCase)
  Else
    RunToCase = Min(RunToCase, iCase - 1 + RunForCase)
  End If

  CaseLoop: Do

    If (EndOfCase) Then

      iCase = iCase + 1

      ! Set flag to false.
      EndOfCase = .false.

      ! Prepare for next case.
      Call PrepareForNextCase(                                            &
             iCase,                                                       &
             Coords, Grids, Domains, Specieses, Sources, SizeDists, Reqs, &
             MainOpts, MultiCaseOpts, DispOptses,                         &
             Results, DispState,                                          &
             EndOfCase, EndOfCases                                        &
           )
      If (EndOfCases) Exit CaseLoop

      Call Message(' ')
      Call Message('Case ' // Trim(Int2Char(iCase)) // ' started')
      Call Message(' ') ! $$ should preceed some error meassages in preparefornextcase
                        ! but should not appear if no next case.

      ! Reset met and flow module instances.
      Call ResetMetsFlows(Mets, Flows)

      ! Initially mark chemistry definition and chemistry state as being uninitialised.
      ChemistryDefn%Initialised  = .false.
      ChemistryState%Initialised = .false.

      Call StartCase(                &
             iCase,                  &
             Coords, Grids, Domains, &
             Specieses,              &
             CloudGammaParamses,     &
             Sources,                &
             Reqs,                   &
             MaterialUnits,          &
             OutputOpts,             &
             SizeDists,               &
             OpenMPOpts,             &
             Units,                  &
             Mets,Flows,             &
             DispState,              &
             Results                 &
           )

    End If

    Call InitialiseRandomSeeds(Max(DispState%MaxPuffs, DispState%MaxParticles), MainOpts%RandomMode)

    RunToTime = TMin(RunToTime, DispState%Time + RunForTime)

    ! Set up chemistry definition (for chemistry runs only).
    ! $$ Ultimately we might define chemistry defns via an input file - the set up
    ! $$ routine could then be moved outside the case loop. However, it currently
    ! $$ needs to access the chemistry flag in DispOpts for each case.
    If (DispOptses%DispOptses(DispState%iDispOpts)%Chemistry) Then

      If (.not.ChemOpts%Initialised) Then
        Call Message(                                                                          &
               'FATAL ERROR: A chemistry run is requested but the chemistry options block ' // &
               'is missing from the input file(s)',                                            &
               3                                                                               &
             )
      End If

      Call SetUpChemistryDefns(Coords, Grids, Specieses, ChemOpts, ChemistryDefn)

    End If

    ! $$ restartable runs - add extra reqs and call setup results to set these up

      Call RunToSuspendOrEndOfCase(                                            &
             RunToTime,                                                        &
             iCase,                                                            &
             Coords, Grids, Domains,                                           &
             Specieses, CloudGammaParamses, SizeDists,                         &
             Sources, TimeDeps,                                                &
             Reqs, MaterialUnits,                                              &
             MainOpts, OutputOpts, DispOptses%DispOptses(DispState%iDispOpts), &
             Units,                                                            &
             Mets, Flows,                                                      &
             Results,                                                          &
             DispState,                                                        &
             ChemOpts, ChemistryDefn, ChemistryState,                          &
             OpenMPOpts,                                                       &
             Suspend, EndOfCase, EndOfCases,                                   &
             RestartOpts,RunToCase                                             &
           )

    If (EndOfCases .or. Suspend) Exit CaseLoop

  End Do CaseLoop

! Output timer summaries
  Call CaseModuleTimerSummary()
  Call OutputModuleTimerSummary()
  Call NWPMetModuleTimerSummary(OpenMPOpts)
  Call IOThreadModuleTimerSummary()
  Call ChemistryModuleTimerSummary()

  Call TimerOff(TotalTimeTimer)
  Call TimerWriteSummary(TotalTimeTimer)

  Call FinaliseTimers(TimerOpts,Units)

  ! Destroy locks in Flow Module
  Call DestroyFlowMetLocks()

  Call FinaliseOpenMPConfiguration(OpenMPOpts)

!!! AJM
  !------------------!
  WRITE(6,*)'Close Particle End Location File'
  !------------------!
  CLOSE(ParticleEndLocationFileUnit)
  WRITE(6,*)'Particle End Location File Closed'
!!! AJM

  !-------------!
  ! End of run. !
  !-------------!

  ! End time of model run.
  ClockTime = CurrentClockTime()

  ! End time message. $$ The 'see log file' bit needs an if test as in routine message itself.
  If (.not.Suspend) Then
    If (GetErrorCode() == 0) Then
      Call Message(' ')
      Call Message(                                          &
             'Run completed successfully at '             // &
             Trim(Time2Char(ClockTime, .true., 3, .true.)),  &
             0,                                              &
             .true.                                          &
           )
    Else If (GetErrorCode() == 1) Then
      Call Message(' ')
      Call Message(                                                   &
             'Run completed with some warnings (see log file) at ' // &
             Trim(Time2Char(ClockTime, .true., 3, .true.)),           &
             0,                                                       &
             .true.                                                   &
           )
    Else If (GetErrorCode() == 2) Then
      Call Message(' ')
      Call Message(                                                 &
             'Run completed with some errors (see log file) at ' // &
             Trim(Time2Char(ClockTime, .true., 3, .true.)),         &
             0,                                                     &
             .true.                                                 &
           )
    End If
  Else
    If (GetErrorCode() == 0) Then
      Call Message(' ')
      Call Message(                                          &
             'Run suspended successfully at '             // &
             Trim(Time2Char(ClockTime, .true., 3, .true.)),  &
             0,                                              &
             .false.                                         &
           )
    Else If (GetErrorCode() == 1) Then
      Call Message(' ')
      Call Message(                                                   &
             'Run suspended with some warnings (see log file) at ' // &
             Trim(Time2Char(ClockTime, .true., 3, .true.)),           &
             0,                                                       &
             .false.                                                  &
           )
    Else If (GetErrorCode() == 2) Then
      Call Message(' ')
      Call Message(                                                 &
             'Run suspended with some errors (see log file) at ' // &
             Trim(Time2Char(ClockTime, .true., 3, .true.)),         &
             0,                                                     &
             .false.                                                &
           )
    End If
  End If

End Program NameIII

!-------------------------------------------------------------------------------------------------------------

Subroutine RunToSuspendOrEndOfCase(                                            &
             RunToTime,                                                        &
             iCase,                                                            &
             Coords, Grids, Domains,                                           &
             Specieses, CloudGammaParamses, SizeDists,                         &
             Sources, TimeDeps,                                                &
             Reqs, MaterialUnits,                                              &             
             MainOpts, OutputOpts, DispOpts,                                   &
             Units,                                                            &
             Mets, Flows,                                                      &
             Results,                                                          &
             DispState,                                                        &
             ChemOpts, ChemistryDefn, ChemistryState,                          &
             OpenMPOpts,                                                       &
             Suspend, EndOfCase, EndOfCases,                                   &
             RestartOpts,RunToCase                                             &
           )

  !$ Use omp_lib

  Use ServiceModule
  Use FlowsModule, Only: Mets_, Flows_
  Use SizeDistModule
  Use SpeciesModule
  Use SourceModule
  Use OutputModule, Only: OutputOpts_, Reqs_, Results_
  Use ParticleModule
  Use PuffModule
  Use ChemistryModule
  Use CaseModule
  Use RestartModule
  Use InputModule

  Use OpenMPModule,   only : GetNumberOfParallelIOThreads          &
                            ,WorkerThreadID                        &
                            ,IOThreadID                            &
                            ,InitIOSynchronisation                 &
                            ,LookaheadFileReadCompleteWait         &
                            ,LookaheadFileReadRequest              &
                            ,IOLookaheadComplete                   &
                            ,IOMaxLookahead                        &
                            ,OpenMPOpts_

  Use IOThreadModule, only : StartIOThread

  Implicit None
   ! Argument list:
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
  Type(RestartOpts_),        Intent(In)    :: RestartOpts
  Integer,                   Intent(In)    :: RunToCase

! Local vars
  Type(Time_)                  :: TargetTime
  Type(Time_)                  :: WriteRestartAtTime
  Integer                      :: myid
  Integer                      :: NumberOfParallelIOThreads

  If (OpenMPOpts%ParallelMetRead) Then
     IOMaxLookahead = DispState%sLastDispTime + DispState%MaxTSpread
  End if

  NumberOfParallelIOThreads = GetNumberOfParallelIOThreads()

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(MYID) NUM_THREADS(NumberOfParallelIOThreads)

  myid = WorkerThreadID
  !$ myid=omp_get_thread_num()
  Call InitIOSynchronisation(OpenMPOpts)

!$OMP BARRIER

  If (myid.eq.IOThreadID) then ! I am the IO Thread

    If (GetNumberOfParallelIOThreads() < 2) Then
      Call Message('ERROR: Trying to use parallel MetRead with one thread', 4)
    End If
    Call StartIOThread(Coords,Grids,Mets,Units,OpenMPOpts)

  Else If (myid.eq.WorkerThreadID) Then! I am the computation thread

    Do While (.not.(EndOfCase .or. EndOfCases .or. Suspend))

      WriteRestartAtTime = FindNextRestartWriteTime(DispState%Time, RestartOpts)

      TargetTime = TMin(RunToTime, WriteRestartAtTime)

      Call RunToRestartDumpOrSuspendOrEndOfCase(                               &
             TargetTime, RunToTime,                                            &
             iCase,                                                            &
             Coords, Grids, Domains,                                           &
             Specieses, CloudGammaParamses, SizeDists,                         &
             Sources, TimeDeps,                                                &
             Reqs, MaterialUnits,                                              &
             MainOpts, OutputOpts, DispOpts,                                   &
             Units,                                                            &
             Mets, Flows,                                                      &
             Results,                                                          &
             DispState,                                                        &
             ChemOpts, ChemistryDefn, ChemistryState,                          &
             OpenMPOpts,                                                       &
             Suspend, EndOfCase, EndOfCases                                    &
           )

      ! Note DispState%Time can be infinite if a case doesn't need to start. $$ move to case module?
      If (IsInfFuture(DispState%Time)) EndOfCase = .true.

      If (.not.EndOfCases) Then
        If (                                                                   &
          (.not.IsInfFuture(DispState%Time) .and. DispState%Time >= RunToTime) &
          .or.                                                                 &
          (EndOfCase .and. iCase >= RunToCase)                                 &
        ) Then
          Suspend = .true.
        End If
      End If

    ! Write restart file.
      If (                                                                            &
        (.not.IsInfFuture(DispState%Time) .and. DispState%Time >= WriteRestartAtTime) &
        .or.                                                                          &
        (EndOfCase .and. Mod(iCase, RestartOpts%dCase) == 0)                          &
        .or.                                                                          &
        (RestartOpts%WriteOnSuspend .and. Suspend)                                    &
      ) Then

        Call WriteRestartFile(        &
               RestartOpts,           &
               iCase,                 &
               Results,               &
               DispState,             &
               DispOpts,              &
               ChemistryState,        &
               EndOfCase, EndOfCases, &
               Units,                 &
               Reqs                   &
             )

      End If

    End Do

! Tell the IOThread to stop by setting a termination flag
    If (OpenMPOpts%ParallelMetRead) Then
      Call LookaheadFileReadCompleteWait()
      IOLookaheadComplete=.true.
      Call LookaheadFileReadRequest()
    End If

  Else
    Call Message('ERROR: somehow we have requested more than two threads', 4)
  End If

!$OMP END PARALLEL

End Subroutine RunToSuspendOrEndOfCase

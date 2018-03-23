! Module:  Input Module

Module InputModule

! This module provides code for coordinating input from the command line arguments and
! from the main set of headed input files.

!-----------------------------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowsModule, Only: FlowOrder_, FlowSubset_, FlowAttrib_, Flows_, InitFlows,                           &
                       AddPrototypeFlow, AddSingleSiteFlow, AddNWPFlow,                                   &
                       AddBuildingFlow, AddRadarFlow, AddLINCOMFlow,                                      &
                       InitFlowOrder, AddFlowOrder, InitFlowSubset, AddFlowSubset, InitAttrib, AddAttrib, &
                       PrototypeFlow_, SingleSiteFlow_, NWPFlow_,                                         &
                       BuildingFlow_, RadarFlow_, LINCOMFlow_,                                            &
                       InitPrototypeFlow, InitSingleSiteFlow, InitNWPFlow,                                &
                       InitBuildingFlow, InitRadarFlow, InitLINCOMFlow,                                   &
                       Mets_, InitMets,                                                                   &
                       AddPrototypeMet, AddSingleSiteMet, AddNWPMet, AddRadarMet, AddAncillaryMet,        &
                       PrototypeMet_, SingleSiteMet_, NWPMet_, RadarMet_, AncillaryMet_,                  &
                       InitPrototypeMet, InitSingleSiteMet, InitNWPMet, InitRadarMet, InitAncillaryMet
Use NWPMetModule       ! Needed for met defns - move to Mets.F90
Use RadarMetModule     ! Needed for met defns - move to Mets.F90
Use AncillaryMetModule ! Needed for met defns - move to Mets.F90
Use SizeDistModule
Use SpeciesModule
Use SourceModule
Use OutputModule, Only: Quantities, QInfo, P_AvInt, P_Prob, P_Percent, P_MaxMin, P_Av, P_Int,               &
                        OutputOpts_, FieldReq_, PdfReq_, PPInfoReq_, Reqs_,                                 &
                        PreInitOutputOpts, InitOutputOpts, InitReqs, InitProcess, AddProcess, InitFieldReq, &
                        AddFieldReq, InitPdfReq, AddPdfReq, InitPPInfoReq, AddPPInfoReq
Use ParticleModule
Use PuffModule
Use ChemistryModule
Use CaseModule
Use RestartModule
Use OpenMPModule, Only : PreInitOpenMPOpts, InitOpenMPOpts, OpenMPOpts_
Use TimerModule, Only : TimerOpts_, PreInitTimerOpts, InitTimerOpts

!-----------------------------------------------------------------------------------------------------------------------------------

Implicit None

!-----------------------------------------------------------------------------------------------------------------------------------

Private
Public :: ProcessCommandLine !
Public :: ReadInputFiles     !

!-----------------------------------------------------------------------------------------------------------------------------------

Type :: InputFiles_ ! A collection of input files. ! $$ move to headed file?
  Integer                      :: nInputFiles
  Character(MaxFileNameLength) :: InputFiles(MaxFilesPerSet)
End Type InputFiles_

!-----------------------------------------------------------------------------------------------------------------------------------

Contains

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine ProcessCommandLine(           &
             MainInputFile,              &
             LogFile, ErrorFile,         &
             ClosePromptLevel,           &
             RestartFileStem,            &
             Restart, RestartFileSuffix, &
             CharRunForTime,             &
             CharRunToTime               &
           )
!

  Implicit None
  ! Argument list:
  Character(MaxFileNameLength), Intent(Out) :: MainInputFile     !
  Character(MaxFileNameLength), Intent(Out) :: LogFile           !
  Character(MaxFileNameLength), Intent(Out) :: ErrorFile         !
  Integer,                      Intent(Out) :: ClosePromptLevel  !
  Character(MaxFileNameLength), Intent(Out) :: RestartFileStem   !
  Logical,                      Intent(Out) :: Restart           !
  Character(MaxFileNameLength), Intent(Out) :: RestartFileSuffix ! Blank = latest restart file
  Character(MaxCharLength),     Intent(Out) :: CharRunForTime    ! Blank = no run-for limit
  Character(MaxCharLength),     Intent(Out) :: CharRunToTime     ! Blank = no run-to limit
  ! Locals:
  Integer                   :: ErrorCode
  Integer                   :: i
  Integer                   :: iDot
  Integer                   :: s
  Character(MaxArgLength) :: Keywords(8)
  Character(MaxArgLength) :: Tokens(8)
  Logical                   :: Set(8)

  Keywords = (/                    &
               'InputFile       ', & ! 1
               'LogFolder       ', & ! 2
               'ClosePromptLevel', & ! 3
               'RestartFolder   ', & ! 4
               'Restart         ', & ! 5
               'NewOutputFilesxx', & ! 6 ! $$ remove
               'RunFor          ', & ! 7
               'RunTo           '  & ! 8
             /)

  ! Defaults.
  Tokens = (/            &
             '        ', & ! 1
             '        ', & ! 2
             '       0', & ! 3
             '        ', & ! 4
             '        ', & ! 5
             '        ', & ! 6
             '        ', & ! 7 ! Time defaults set later after know if run is backwards
             '        '  & ! 8 !
           /)

  Call GetProcessedCommandLineArguments(1, Keywords, Tokens, Set)

  If (Len_Trim(Tokens(1)) > MaxFileNameLength) Then
    Call Message(                                               &
           'Error in ProcessCommandLine: '                   // &
           'the name of the main input file is longer than ' // &
           Trim(Int2Char(MaxFileNameLength))                 // &
           ' characters.',                                      &
           3                                                    &
         )
  End If
  MainInputFile = Tokens(1)

  If (Len_Trim(Tokens(2)) > MaxFileNameLength) Then
    Call Message(                                          &
           'Error in ProcessCommandLine: '              // &
           'the name of the log folder is longer than ' // &
           Trim(Int2Char(MaxFileNameLength))            // &
           ' characters.',                                 &
           3                                               &
         )
  End If
  LogFile = Tokens(2)

  ClosePromptLevel = Char2Int(Tokens(3), ErrorCode)
  If (ErrorCode /= 0) Then
    Call Message('FATAL ERROR reading command line arguments', 3)
  End If
  ! Checks - use error code $$

  If (Len_Trim(Tokens(4)) > MaxFileNameLength) Then
    Call Message(                                  &
           'Error in ProcessCommandLine: '      // &
           'the restart folder is longer than ' // &
           Trim(Int2Char(MaxFileNameLength))    // &
           ' characters.',                         &
           3                                       &
         )
  End If
  RestartFileStem = Tokens(4)

  If (Len_Trim(Tokens(5)) > MaxFileNameLength) Then
    Call Message(                                       &
           'Error in ProcessCommandLine: '           // &
           'the restart file suffix is longer than ' // &
           Trim(Int2Char(MaxFileNameLength))         // &
           ' characters.',                              &
           3                                            &
         )
  End If
  Restart           = Set(5)
  RestartFileSuffix = Tokens(5)

  CharRunForTime = Tokens(7)

  CharRunToTime = Tokens(8)

  ! Determine the main input file. ! $$ Should somewhere write command line options plus
  !                                     response to 'what file do you want to use?'
  If (MainInputFile == ' ') Then
    Call SetActiveFocus(6)
    Write (6, *)
    Write (6, *) 'No input file found on command line - what file do you want to use?'
    Read  (5, *) MainInputFile
  End If

  ErrorFile = LogFile

  iDot = Scan(MainInputFile, '.', Back = .true.)
  i    = Scan(MainInputFile, '\/', Back = .true.)

  ! Determine the log file.
  If (LogFile == ' ') Then
    If (iDot == 0 .or. iDot < i) Then
      LogFile = Trim(MainInputFile) // 'Log.txt'
    Else
      LogFile = MainInputFile(1:iDot - 1) // 'Log.txt'
    End If
  Else
    If (Scan(LogFile, '\/', Back = .true.) /= Len_Trim(LogFile)) Then
      LogFile = Trim(LogFile) // '\'
    End If
    If (iDot == 0 .or. iDot < i) Then
      LogFile = Trim(LogFile) // Trim(MainInputFile(i + 1: )) // 'Log.txt'
    Else
      LogFile = Trim(LogFile) // MainInputFile(i + 1:iDot - 1) // 'Log.txt'
    End If
  End If
  LogFile = ConvertFileName(LogFile)

  ! Determine the error file.
  If (ErrorFile == ' ') Then
    If (iDot == 0 .or. iDot < i) Then
      ErrorFile = Trim(MainInputFile) // 'Error.txt'
    Else
      ErrorFile = MainInputFile(1:iDot - 1) // 'Error.txt'
    End If
  Else
    If (Scan(ErrorFile, '\/', Back = .true.) /= Len_Trim(ErrorFile)) Then
      ErrorFile = Trim(ErrorFile) // '\'
    End If
    If (iDot == 0 .or. iDot < i) Then
      ErrorFile = Trim(ErrorFile) // Trim(MainInputFile(i + 1: )) // 'Error.txt'
    Else
      ErrorFile = Trim(ErrorFile) // MainInputFile(i + 1:iDot - 1) // 'Error.txt'
    End If
  End If
  ErrorFile = ConvertFileName(ErrorFile)

  ! Determine the restart file stem.
  If (RestartFileStem == ' ') Then
    If (iDot == 0 .or. iDot < i) Then
      RestartFileStem = Trim(MainInputFile)
    Else
      RestartFileStem = MainInputFile(1:iDot - 1)
    End If
  Else
    If (Scan(RestartFileStem, '\/', Back = .true.) /= Len_Trim(RestartFileStem)) Then
      RestartFileStem = Trim(RestartFileStem) // '\'
    End If
    If (iDot == 0 .or. iDot < i) Then
      RestartFileStem = Trim(RestartFileStem) // Trim(MainInputFile(i + 1: ))
    Else
      RestartFileStem = Trim(RestartFileStem) // MainInputFile(i + 1:iDot - 1)
    End If
  End If

  ! $$ check above logic water tight - e.g. if maininputfile is ./input (OK I think)

  ! Check main input and log files are not the same.
  If (MainInputFile == LogFile) Then
    Call Message(                                                         &
           'ERROR in ProcessCommandLine: '                             // &
           'the main input file and the log file have the same name.',    &
           3                                                              &
         )
  End If

  ! Check restart files won't clash (probably elsewhere)

  ! Replace _ by space in times.
  Do
    s = Scan(CharRunForTime, '_')
    If (s == 0) Exit
    CharRunForTime(s:s) = ' '
  End Do

  Do
    s = Scan(CharRunToTime, '_')
    If (s == 0) Exit
    CharRunToTime(s:s) = ' '
  End Do

End Subroutine ProcessCommandLine

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine GetProcessedCommandLineArguments(MaxNonKeywordArguments, Keywords, Tokens, Set)
!

  Implicit None
  ! Argument list:
  Integer,      Intent(In)    :: MaxNonKeywordArguments
  Character(*), Intent(In)    :: Keywords(:)
  Character(*), Intent(InOut) :: Tokens(:) ! Contains defaults on input
  Logical,      Intent(Out)   :: Set(:)
  ! Locals:
  Integer                  :: nArguments
  Character(Len(Keywords)) :: Arguments(Size(Keywords))
  Logical                  :: KeywordArgUsed
  Integer                  :: s
  Integer                  :: i
  Integer                  :: j

  ! Check keywords not CIEq or blank. check value of MaxNonKeywordArguments $$
  ! Arrays the same size?
  ! Align with fortran 2003 convention?
  ! Note -keyword=blank, -keyword and no argument all return default value. For first two however Set is true.
  ! Set is like Present.
  ! Move routine to System.F90

  ! Get command line arguments.
  Call GetCommandLineArguments(nArguments, Arguments)

  Set(:)         = .false.
  KeywordArgUsed = .false.

  Do i = 1, nArguments

    ! Keyword argument.
    If (Arguments(i)(1:1) == '-') Then

      s = Scan(Arguments(i), '=')
      If (s == 0) s = Len_Trim(Arguments(i)) + 1
      Do j = 1, Size(Keywords)
        If (Arguments(i)(2:s - 1) .CIEq. Keywords(j)) Then
          If (Set(j)) Then
            Call Message(                                                  &
              'FATAL ERROR in processing argument list - the argument ' // &
              Trim(Keywords(j))                                         // &
              ' appears more than once',                                   &
              3                                                            &
            )
          End If
          If (s + 1 <= Len(Arguments(i))) Then
            If (Arguments(i)(s + 1:) /= ' ') Then
              Tokens(j) = Arguments(i)(s + 1:)
            End If
          End If
          Set(j)         = .true.
          KeywordArgUsed = .true.
          Exit
        End If
        If (j == Size(Keywords)) Then
          Call Message(                                                  &
            'FATAL ERROR in processing argument list - the argument ' // &
            Arguments(i)(2:s - 1)                                     // &
            ' has not been recognised',                                  &
            3                                                            &
          )
        End If
      End Do

    ! Non-keyword argument.
    Else If (i <= MaxNonKeywordArguments .and. .not.KeywordArgUsed) Then

      Tokens(i) = Arguments(i)
      Set(i)    = .true.

    ! Non-keyword argument where keyword needed.
    Else

      Call Message(                                                  &
        'FATAL ERROR in processing argument list - the argument ' // &
        Trim(Arguments(i))                                        // &
        ' does not start with "-"',                                  &
        3                                                            &
      )

    End If

  End Do

End Subroutine GetProcessedCommandLineArguments

!-----------------------------------------------------------------------------------------------------------------------------------

Function InitInputFiles()
! Initialises a collection of input files.

  Implicit None
  ! Function result:
  Type(InputFiles_) :: InitInputFiles !
  ! Locals:
  Type(InputFiles_) :: InputFiles !

  InputFiles%nInputFiles = 0

  InitInputFiles = InputFiles

End Function InitInputFiles

!-----------------------------------------------------------------------------------------------------------------------------------

Function InitInputFile(InputFileName)
! Initialises an input file name.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: InputFileName ! Input file name.
  ! Function result:
  Character(MaxFileNameLength) :: InitInputFile ! Input file name.
  ! Locals:
  Type(InputFiles_) :: InputFiles !

  If (Len_Trim(InputFileName) > MaxFileNameLength .or. &
      Len_Trim(InputFileName) == 0) Then
    Call Message('Error in InitInputFile', 3)
  End If

  InitInputFile = InputFileName

End Function InitInputFile

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine AddInputFile(InputFile, InputFiles)
! Adds an input file to the collection of input files.

  Implicit None
  ! Argument list:
  Character(*),      Intent(In)    :: InputFile  ! Input file to be added.
  Type(InputFiles_), Intent(InOut) :: InputFiles !
  ! Locals:
  Integer :: iArray   ! Index of the array.
  Integer :: iElement ! Index of the array element.
  Integer :: i        !

  If (InputFiles%nInputFiles >= MaxFilesPerSet) Then ! $$ check duplicate names
    Call Message('Error in AddInputFile', 3)
  End If
  InputFiles%nInputFiles                        = InputFiles%nInputFiles + 1
  InputFiles%InputFiles(InputFiles%nInputFiles) = InputFile

End Subroutine AddInputFile

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine ReadInputFiles(                                               &
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
! .

  Implicit None
  ! Argument list:
  Character(MaxFileNameLength), Intent(In)    :: MainInputFile
  Type(MainOpts_),              Intent(Out)   :: MainOpts                           !
  Type(RestartOpts_),           Intent(Out)   :: RestartOpts                        !
  Type(MultiCaseOpts_),         Intent(Out)   :: MultiCaseOpts                      !
  Type(ChemOpts_),              Intent(Out)   :: ChemOpts                           !
  Type(OutputOpts_),            Intent(Out)   :: OutputOpts                         !
  Type(OpenMPOpts_),            Intent(Out)   :: OpenMPOpts                         !
  Type(TimerOpts_),             Intent(Out)   :: TimerOpts                          !
  Type(EtaDefns_),              Intent(Out)   :: EtaDefns                           !
  Type(Coords_),                Intent(Out)   :: Coords                             !
  Type(Locationses_),           Intent(Out)   :: Locationses                        !
  Type(Grids_),                 Intent(Out)   :: Grids                              !
  Type(Domains_),               Intent(Out)   :: Domains                            !
  Type(Mets_),                  Intent(Out)   :: Mets                               !  
  Type(Flows_),                 Intent(Out)   :: Flows                              !
  Type(Specieses_),             Intent(Out)   :: Specieses                          !
  Type(CloudGammaParamses_),    Intent(Out)   :: CloudGammaParamses                 !
  Type(MassLimits_),            Intent(Out)   :: MassLimits                         !
  Type(Sources_),               Intent(Out)   :: Sources                            !
  Type(SizeDists_),             Intent(Out)   :: SizeDists                          !
  Type(TimeDeps_),              Intent(Out)   :: TimeDeps                           !
  Type(Reqs_),                  Intent(Out)   :: Reqs                               !
  Type(DispOptses_),            Intent(Out)   :: DispOptses                         !
  Real(Std),                    Intent(Out)   :: BlendingHeights(MaxFlows,MaxFlows) ! $$
  Type(MaterialUnits_),         Intent(InOut) :: MaterialUnits 
  Type(Units_),                 Intent(InOut) :: Units
  ! Locals:
  Type(HFBlockForms_)      :: HFBlockForms
  Type(HFState_)           :: HFState
  Type(Tokens_)            :: Tokens
  Type(Arrays_)            :: Arrays
  Type(MetDefns_)          :: MetDefns
  Type(RadarMetDefns_)     :: RadarMetDefns
  Type(AncillaryMetDefns_) :: AncillaryMetDefns
  Character(MaxLineLength) :: Line
  Type(InputFiles_)        :: InputFiles
  Integer                  :: i
  Logical                  :: EndOfHFs

  Call Message('Main input file used is ' // Trim(MainInputFile))

  ! $$ list other files below.

  ! Initialise collection of names of input files.

  InputFiles%nInputFiles   = 1
  InputFiles%InputFiles(1) = MainInputFile

  ! Preliminary pass through main input file to read names of input files. Note the
  ! main input file is always included in the list of input files and should not be
  ! listed in the main input file.

  Call InitHFState(InputFiles%InputFiles(1:InputFiles%nInputFiles), HFState)

  ! Initialise HFBlockForms.
  Call InitHFBlockForms   (HFBlockForms)
  Call InputFileInputNames(HFBlockForms)

  ! Read file and sort out tokens.
  Do
    Call ReadHFs(HFBlockForms, Tokens, EndOfHFs, Units, HFState)
    If (EndOfHFs) Exit
    If (Tokens%BlockKey .CIEq. 'Input Files:') Then
      Call Tokens2InputFiles(Tokens, InputFiles)
    End If
  End Do

  ! Dependencies: MainOpts needed before any time operations that need to know if
  !               forwards or backwards, what time frame.
  ! $$ others

  ! First pass through main input file and input files named therein.
  ! Read: main options.

  MainOpts = PreInitMainOpts()

  ! Initialise HFState and HFBlockForms.
  Call InitHFState(InputFiles%InputFiles(1:InputFiles%nInputFiles), HFState)
  Call InitHFBlockForms  (HFBlockForms)
  Call MainOptsInputNames(HFBlockForms)

  ! Read file and sort out tokens.
  Do

    Call ReadHFs(HFBlockForms, Tokens, EndOfHFs, Units, HFState)

    If (EndOfHFs) Exit

    If (Tokens%BlockKey .CIEq. 'Main Options:') Then

      Call Tokens2MainOpts(Tokens, MainOpts)

    End If

  End Do

  ! Initialise time module.
  Call InitTimeModule(MainOpts%TimeFrameName, MainOpts%Backwards)

  ! Second pass through main input file and input files named therein.
  ! Read: restart file options,
  !       multiple case options,
  !       output options,
  !       arrays,
  !       eta definitions,
  !       locations.

  RestartOpts           = PreInitRestartOpts()
  MultiCaseOpts         = PreInitMultiCaseOpts()
  ChemOpts              = PreInitChemOpts()
  OutputOpts            = PreInitOutputOpts()
  OpenMPOpts            = PreInitOpenMPOpts()
  TimerOpts             = PreInitTimerOpts()
  GlobalMessageControls = InitMessageControls()
  Arrays                = InitArrays()
  EtaDefns              = InitEtaDefns()
  Locationses           = InitLocationses()

  ! Initialise HFState and HFBlockForms.
  Call InitHFState(InputFiles%InputFiles(1:InputFiles%nInputFiles), HFState)
  Call InitHFBlockForms        (HFBlockForms)
  Call RestartOptsInputNames   (HFBlockForms)
  Call MultiCaseOptsInputNames (HFBlockForms)
  Call MessageControlInputNames(HFBlockForms)
  Call ChemOptsInputNames      (HFBlockForms)
  Call OutputOptsInputNames    (HFBlockForms)
  Call OpenMPOptsInputNames    (HFBlockForms)
  Call TimerOptsInputNames     (HFBlockForms)
  Call ArrayInputNames         (HFBlockForms)
  Call EtaDefnInputNames       (HFBlockForms)
  Call LocationsInputNames     (HFBlockForms)

  ! Read file and sort out tokens.
  Do

    Call ReadHFs(HFBlockForms, Tokens, EndOfHFs, Units, HFState)

    If (EndOfHFs) Exit

    If      (Tokens%BlockKey .CIEq. 'Restart File Options:') Then

      Call Tokens2RestartOpts   (Tokens, MainOpts, RestartOpts)

    Else If (Tokens%BlockKey .CIEq. 'Multiple Case Options:') Then

      Call Tokens2MultiCaseOpts (Tokens, MultiCaseOpts)

    Else If (Tokens%BlockKey .CIEq. 'Message Controls:') Then

      Call Tokens2MessageControl(Tokens, GlobalMessageControls)

    Else If (Tokens%BlockKey .CIEq. 'Chemistry Options:') Then

      Call Tokens2ChemOpts      (Tokens, ChemOpts)

    Else If (Tokens%BlockKey .CIEq. 'Output Options:') Then

      Call Tokens2OutputOpts    (Tokens, OutputOpts)

    Else If (Tokens%BlockKey .CIEq. 'OpenMP Options:') Then

      Call Tokens2OpenMPOpts    (Tokens, OpenMPOpts)

    Else If (Tokens%BlockKey .CIEq. 'Timer Options:') Then

      Call Tokens2TimerOpts    (Tokens, TimerOpts, Units, MainInputFile)

    Else If (Tokens%BlockKey .CIEq. 'Array:') Then

      Call Tokens2Array         (Tokens, Arrays)

    Else If (Tokens%BlockKey .CIEq. 'Eta Definition:') Then

      Call Tokens2EtaDefn       (Tokens, EtaDefns)

    Else If (Tokens%BlockKey .CIEq. 'Locations:') Then

      Call Tokens2Locations     (Tokens, Locationses)

    End If

  End Do

  ! Third pass through main input file and input files named therein.
  ! Read: domains,
  !       met definitions,
  !       coord systems.

  Domains           = InitDomains()
  MetDefns          = InitMetDefns()
  RadarMetDefns     = InitRadarMetDefns()
  AncillaryMetDefns = InitAncillaryMetDefns()
  Coords            = InitCoords()

  ! Initialise HFState and HFBlockForms.
  Call InitHFState(InputFiles%InputFiles(1:InputFiles%nInputFiles), HFState)
  Call InitHFBlockForms           (HFBlockForms)
  Call DomainInputNames           (HFBlockForms)
  Call NWPMetDefnInputNames       (HFBlockForms)
  Call NWPMetDefn2InputNames      (HFBlockForms)
  Call RadarMetDefnInputNames     (HFBlockForms)
  Call RadarMetDefn2InputNames    (HFBlockForms)
  Call AncillaryMetDefnInputNames (HFBlockForms)
  Call AncillaryMetDefn2InputNames(HFBlockForms)
  Call HCoordInputNames           (HFBlockForms)
  Call ZCoordInputNames           (HFBlockForms)

  ! Read file and sort out tokens.
  Do

    Call ReadHFs(HFBlockForms, Tokens, EndOfHFs, Units, HFState)

    If (EndOfHFs) Exit

    If      (Tokens%BlockKey .CIEq. 'Domains:') Then

      Call Tokens2Domain           (Tokens, Domains)

    Else If (Tokens%BlockKey .CIEq. 'NWP Met Definitions:') Then

      Call Tokens2NWPMetDefn       (Tokens, MetDefns)

    Else If (Tokens%BlockKey .CIEq. 'NWP Met File Structure Definition:') Then

      Call Tokens2NWPMetDefn2      (Tokens, MetDefns)

    Else If (Tokens%BlockKey .CIEq. 'Radar Met Definitions:') Then

      Call Tokens2RadarMetDefn       (Tokens, RadarMetDefns)

    Else If (Tokens%BlockKey .CIEq. 'Radar Met File Structure Definition:') Then

      Call Tokens2RadarMetDefn2      (Tokens, RadarMetDefns)

    Else If (Tokens%BlockKey .CIEq. 'Ancillary Met Definitions:') Then

      Call Tokens2AncillaryMetDefn (Tokens, AncillaryMetDefns)

    Else If (Tokens%BlockKey .CIEq. 'Ancillary Met File Structure Definition:') Then

      Call Tokens2AncillaryMetDefn2(Tokens, AncillaryMetDefns)

    Else If (Tokens%BlockKey .CIEq. 'Horizontal Coordinate Systems:') Then

      Call Tokens2HCoord           (Tokens, Coords)

    Else If (Tokens%BlockKey .CIEq. 'Vertical Coordinate Systems:') Then

      Call Tokens2ZCoord           (Tokens, EtaDefns, Coords)

    End If

  End Do

  ! Forth pass through main input file and input files named therein.
  ! Read: grids.

  Grids = InitGrids()

  ! Initialise HFState and HFBlockForms.
  Call InitHFState(InputFiles%InputFiles(1:InputFiles%nInputFiles), HFState)
  Call InitHFBlockForms(HFBlockForms)
  Call HGridInputNames (HFBlockForms)
  Call ZGridInputNames (HFBlockForms)
  Call TGridInputNames (HFBlockForms)

  ! Read file and sort out tokens.
  Do

    Call ReadHFs(HFBlockForms, Tokens, EndOfHFs, Units, HFState)

    If (EndOfHFs) Exit

    If      (Tokens%BlockKey .CIEq. 'Horizontal Grids:') Then

      Call Tokens2HGrid(Tokens, Arrays, Locationses, Grids)

    Else If (Tokens%BlockKey .CIEq. 'Vertical Grids:') Then

      Call Tokens2ZGrid(Tokens, Arrays, Coords, Grids)

    Else If (Tokens%BlockKey .CIEq. 'Temporal Grids:') Then

      Call Tokens2TGrid(Tokens, Arrays, Grids)

    End If

  End Do

  ! Fifth pass through main input file and input files named therein.
  ! Read: met and flow modules instances,
  !       species,
  !       CloudGammaParams,
  !       particle mass limits,
  !       sources,
  !       particle size distributions.
  !       source time dependencies,
  !       output requirements,
  !       dispersion options,

  Mets               = InitMets()
  Flows              = InitFlows(MainOpts%SameResultsWithUpdateOnDemand)
  SizeDists          = InitSizeDists()
  Specieses          = InitSpecieses()
  CloudGammaParamses = InitCloudGammaParamses()
  MassLimits         = InitMassLimits()
  Sources            = InitSources(MainOpts%MaxSources)
  TimeDeps           = InitTimeDeps()
  Reqs               = InitReqs(MainOpts%MaxFieldReqs, MainOpts%MaxFieldOutputGroups)
  DispOptses         = InitDispOptses()

  ! Initialise HFState and HFBlockForms.
  Call InitHFState(InputFiles%InputFiles(1:InputFiles%nInputFiles), HFState)
  Call InitHFBlockForms          (HFBlockForms)
  Call PrototypeMetInputNames    (HFBlockForms)
  Call SingleSiteMetInputNames   (HFBlockForms)
  Call NWPMetInputNames          (HFBlockForms)
  Call RadarMetInputNames        (HFBlockForms)
  Call AncillaryMetInputNames    (HFBlockForms)
  Call PrototypeFlowInputNames   (HFBlockForms)
  Call SingleSiteFlowInputNames  (HFBlockForms)
  Call NWPFlowInputNames         (HFBlockForms)
  Call BuildingFlowInputNames    (HFBlockForms)
  Call RadarFlowInputNames       (HFBlockForms)
  Call LINCOMFlowInputNames      (HFBlockForms)
  Call FlowOrderInputNames       (HFBlockForms)
  Call FlowSubsetInputNames      (HFBlockForms)
  Call FlowAttribInputNames      (HFBlockForms)
  Call SizeDistsInputNames       (HFBlockForms)
  Call SpeciesInputNames         (HFBlockForms)
  Call SpeciesUsesInputNames     (HFBlockForms)
  Call CloudGammaParamsInputNames(HFBlockForms)
  Call MassLimitInputNames       (HFBlockForms)
  Call SourcesInputNames         (HFBlockForms) ! Singular plural not consistent. $$
  Call TimeDepsInputNames        (HFBlockForms)
  Call FieldReqsInputNames       (HFBlockForms)
  Call PdfReqsInputNames         (HFBlockForms)
  Call PPInfoReqsInputNames      (HFBlockForms)
  Call DispOptsesInputNames      (HFBlockForms)

  ! Read file and sort out tokens.
  Do

    Call ReadHFs(HFBlockForms, Tokens, EndOfHFs, Units, HFState)

    If (EndOfHFs) Exit

    If      (Tokens%BlockKey .CIEq. 'Prototype Met Module Instances:') Then

      Call Tokens2PrototypeMet  (Tokens, MainOpts, Mets)

    Else If (Tokens%BlockKey .CIEq. 'Single Site Met Module Instances:') Then

      Call Tokens2SingleSiteMet (Tokens, MainOpts, Mets)

    Else If (Tokens%BlockKey .CIEq. 'NWP Met Module Instances:') Then

      Call Tokens2NWPMet        (Tokens, MainOpts, OpenMPOpts, Arrays, MetDefns, Mets)

    Else If (Tokens%BlockKey .CIEq. 'Radar Met Module Instances:') Then

      Call Tokens2RadarMet      (Tokens, MainOpts, RadarMetDefns, Mets)

    Else If (Tokens%BlockKey .CIEq. 'Ancillary Met Module Instances:') Then

      Call Tokens2AncillaryMet  (Tokens, MainOpts, Arrays, AncillaryMetDefns, Mets)

    Else If (Tokens%BlockKey .CIEq. 'Prototype Flow Module Instances:') Then

      Call Tokens2PrototypeFlow (Tokens, MainOpts, Flows)

    Else If (Tokens%BlockKey .CIEq. 'Single Site Flow Module Instances:') Then

      Call Tokens2SingleSiteFlow(Tokens, MainOpts, Flows)

    Else If (Tokens%BlockKey .CIEq. 'NWP Flow Module Instances:') Then

      Call Tokens2NWPFlow       (Tokens, MainOpts, Flows)

    Else If (Tokens%BlockKey .CIEq. 'Building Flow Module Instances:') Then

      Call Tokens2BuildingFlow  (Tokens, MainOpts, Flows)

    Else If (Tokens%BlockKey .CIEq. 'Radar Flow Module Instances:') Then

      Call Tokens2RadarFlow     (Tokens, MainOpts, Flows)

    Else If (Tokens%BlockKey .CIEq. 'LINCOM Flow Module Instances:') Then

      Call Tokens2LINCOMFlow    (Tokens, MainOpts, Flows)

    Else If (Tokens%BlockKey .CIEq. 'Flow Order:') Then

      Call Tokens2FlowOrder     (Tokens, Flows)

    Else If (Tokens%BlockKey .CIEq. 'Flow Subset:') Then

      Call Tokens2FlowSubset    (Tokens, Flows)

    Else If (Tokens%BlockKey .CIEq. 'Flow Attributes:') Then

      Call Tokens2FlowAttrib    (Tokens, Flows)

    Else If (Tokens%BlockKey .CIEq. 'Particle Size Distribution:') Then

      Call Tokens2SizeDists     (Tokens, SizeDists)

    Else If (Tokens%BlockKey .CIEq. 'Species:') Then

      Call Tokens2Species       (Tokens, Specieses, MaterialUnits)

    Else If (Tokens%BlockKey .CIEq. 'Species Uses:') Then

      Call Tokens2SpeciesUses   (Tokens, Specieses)

    Else If (Tokens%BlockKey .CIEq. 'Cloud Gamma Parameters:') Then

      Call Tokens2CloudGammaParams(Tokens, CloudGammaParamses)

    Else If (Tokens%BlockKey .CIEq. 'Particle Mass Limits:') Then

      Call Tokens2MassLimit     (Tokens, MassLimits)

    Else If (Tokens%BlockKey .CIEq. 'Sources:') Then

      Call Tokens2Sources       (Tokens, MainOpts%FixedMet, Locationses, Sources)

    Else If (Tokens%BlockKey .CIEq. 'Source Time Dependency:') Then

      Call Tokens2TimeDeps      (Tokens, TimeDeps)

    Else If (Tokens%BlockKey .CIEq. 'Output Requirements - Fields:') Then

      Call Tokens2FieldReqs     (Tokens, Grids, Reqs, MaterialUnits)

    Else If (Tokens%BlockKey .CIEq. 'Output Requirements - Pdfs:') Then

      Call Tokens2PdfReqs       (Tokens, Grids, Reqs)

    Else If (Tokens%BlockKey .CIEq. 'Output Requirements - Sets of Particle/Puff Information:') Then

      Call Tokens2PPInfoReqs    (Tokens, Reqs, OpenMPOpts)

    Else If (Tokens%BlockKey .CIEq. 'Sets of Dispersion Options:') Then

      Call Tokens2DispOptses    (Tokens, MainOpts%FixedMet, DispOptses)

    End If

  End Do

  BlendingHeights = 10000.0 ! $$

End Subroutine ReadInputFiles

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine MainOptsInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                                      &
         BlockKey             = 'Main Options:',                                   &
         NamedBlock           = .false.,                                           &
         MatchMultipleLines   = .false.,                                           &
         nHeaderLines         = 1,                                                 &
         ColumnKeys           = Reshape(                                           &
                                  [                                                &
                                    'Absolute Or Relative Time?                 ', & ! 1 ! rename 'Time Frame' $$
                                    'Fixed Met?                                 ', & ! 2
                                    'Time Of Fixed Metxxx                       ', & ! 3 ! remove $$
                                    'Flat Earth?                                ', & ! 4
                                    'Run Name                                   ', & ! 5
                                    'Run-To File                                ', & ! 6
                                    'xxxxxxxxxxxxxxxxxxxxxx                     ', & ! 7 ! remove $$
                                    'Random Seed                                ', & ! 8 ! move to dispopts $$
                                    'Backwards?                                 ', & ! 9
                                    'Max # Sources                              ', & ! 10
                                    'Same Results With/Without Update On Demand?', & ! 11
                                    'Max # Field Reqs                           ', & ! 12
                                    'Max # Field Output Groups                  '  & ! 13
                                  ],                                               &
                                  [ 13, 1 ]                                        &
                                ),                                                 &
         Defaults             = [                                                  &
                                  '   ', & ! 1
                                  '   ', & ! 2
                                  '   ', & ! 3 ! remove $$
                                  '   ', & ! 4
                                  '   ', & ! 5
                                  '   ', & ! 6
                                  '   ', & ! 7 ! remove $$
                                  '   ', & ! 8 ! move to dispopts $$
                                  '   ', & ! 9
                                  '100',                                           & ! 10
                                  'No ',                                           & ! 11
                                  '650',                                           & ! 12
                                  '100'                                            & ! 13
                                ],                                                 &
         TwoD                 = .false.,                                           &
         UnrecognisedMessages = .true.,                                            &
         DisjointColSpecs     = .true.,                                            &
         HFBlockForms         = HFBlockForms                                       &
       )

End Subroutine MainOptsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2MainOpts(Tokens, MainOpts)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)  :: Tokens
  Type(MainOpts_), Intent(Out) :: MainOpts
  ! Local parameters:
  Integer, Parameter :: i_AbsoluteOrRelativeTime               = 1 
  Integer, Parameter :: i_FixedMet                             = 2
  Integer, Parameter :: i_TimeOfFixedMetxxx                    = 3 
  Integer, Parameter :: i_FlatEarth                            = 4
  Integer, Parameter :: i_RunName                              = 5
  Integer, Parameter :: i_RunToFile                            = 6
  Integer, Parameter :: i_xxxxxxxxxxxxxxxxxxxxxx               = 7 
  Integer, Parameter :: i_RandomSeed                           = 8 
  Integer, Parameter :: i_Backwards                            = 9
  Integer, Parameter :: i_MaxNSources                          = 10
  Integer, Parameter :: i_SameResultsWithWithoutUpdateOnDemand = 11
  Integer, Parameter :: i_MaxNFieldReqs                        = 12
  Integer, Parameter :: i_MaxNFieldOutputGroups                = 13
  ! Locals:
  Character(MaxTokenLength) :: TimeFrameName
  Logical                   :: FixedMet
  Logical                   :: FlatEarth
  Logical                   :: Backwards
  Character(1)              :: RandomMode
  Integer                   :: MaxSources
  Logical                   :: SameResultsWithUpdateOnDemand
  Integer                   :: MaxFieldReqs
  Integer                   :: MaxFieldOutputGroups

  If (Tokens%Tokens(i_Backwards) == ' ') Then
    Backwards = .false.
  Else
    Backwards = Token2Log(                                &
                  C         = Tokens%Tokens(i_Backwards), &
                  BlockKey  = 'Main Options',             & ! $$ Could use Tokens%BlockKey ?
                  Item      = ' ',                        &
                  ColumnKey = 'Backwards?'                &
                )
  End If

  TimeFrameName = Tokens%Tokens(i_AbsoluteOrRelativeTime)
  If (TimeFrameName .CIEq. 'Absolute') TimeFrameName = 'Gregorian' ! allow 'absolute' for
                                                                   ! backwards compatability

  If      (TimeFrameName .CIEq. 'Gregorian'   ) Then
  Else If (TimeFrameName .CIEq. '360-day year') Then
  Else If (TimeFrameName .CIEq. 'Relative'    ) Then
  Else
    Call Message('Error in Tokens2MainOpts', 3) ! $$ not quite clear where this test is best placed.
                                                ! Could go in InitMainOpts with time module making list of
                                                ! allowed time frames available.
  End If

  FixedMet = Token2Log(                               &
               C         = Tokens%Tokens(i_FixedMet), &
               BlockKey  = 'Main Options',            &
               Item      = ' ',                       &
               ColumnKey = 'Fixed Met?'               &
             )

  FlatEarth = Token2Log(                                &
                C         = Tokens%Tokens(i_FlatEarth), &
                BlockKey  = 'Main Options',             &
                Item      = ' ',                        &
                ColumnKey = 'Flat Earth?'               &
              )

  If (Tokens%Tokens(i_RandomSeed) == ' ') Then
    RandomMode = 'F'
  Else If (Tokens%Tokens(i_RandomSeed) .CIEq. 'Fixed') Then
    RandomMode = 'F'
  Else If (Tokens%Tokens(i_RandomSeed) .CIEq. 'Random') Then
    RandomMode = 'R'
  Else If (Tokens%Tokens(i_RandomSeed) .CIEq. 'Fixed (Parallel)') Then
    RandomMode = 'P'
  Else If (Tokens%Tokens(i_RandomSeed) .CIEq. 'Input (Parallel)') Then
    RandomMode = 'I'
  Else If (Tokens%Tokens(i_RandomSeed) .CIEq. 'Random (Parallel)') Then
    RandomMode = 'A'
  End If

  MaxSources = Token2Int(                                  &
                 C         = Tokens%Tokens(i_MaxNSources), &
                 BlockKey  = 'Main Options',               &
                 Item      = ' ',                          &
                 ColumnKey = 'Max # Sources'               &
               )

  SameResultsWithUpdateOnDemand = Token2Log(                                                           &
                                    C         = Tokens%Tokens(i_SameResultsWithWithoutUpdateOnDemand), &
                                    BlockKey  = 'Main Options',                                        &
                                    Item      = ' ',                                                   &
                                    ColumnKey = 'Same Results With Update On Demand?'                  &
                                  )

  MaxFieldReqs = Token2Int(                                    &
                   C         = Tokens%Tokens(i_MaxNFieldReqs), &
                   BlockKey  = 'Main Options',                 &
                   Item      = ' ',                            &
                   ColumnKey = 'Max # Field Reqs'              &
                 )
                
  MaxFieldOutputGroups = Token2Int(                                             &
                            C         = Tokens%Tokens(i_MaxNFieldOutputGroups), &
                            BlockKey  = 'Main Options',                         &
                            Item      = ' ',                                    &
                            ColumnKey = 'Max # Field Output Groups'             &
                         )             

  Call InitMainOpts(                                                    &
         RunName                       = Tokens%Tokens(i_RunName),      &
         TimeFrameName                 = TimeFrameName,                 &
         FlatEarth                     = FlatEarth,                     &
         Backwards                     = Backwards,                     &
         FixedMet                      = FixedMet,                      &
         RandomMode                    = RandomMode,                    &
         RunToFile                     = Tokens%Tokens(i_RunToFile),    &
         MaxSources                    = MaxSources,                    &
         SameResultsWithUpdateOnDemand = SameResultsWithUpdateOnDemand, &
         MaxFieldReqs                  = MaxFieldReqs,                  &
         MaxFieldOutputGroups          = MaxFieldOutputGroups,          &
         MainOpts                      = MainOpts                       &
       )

End Subroutine Tokens2MainOpts

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine RestartOptsInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                 &
         BlockKey             = 'Restart File Options:',      &
         NamedBlock           = .false.,                      &
         MatchMultipleLines   = .false.,                      &
         nHeaderLines         = 1,                            &
         ColumnKeys           = Reshape(                      &
                                  (/                          &
                                    '# Cases Between Writes', & ! 1
                                    'Time Between Writes   ', & ! 2
                                    'Delete Old Files?     ', & ! 3
                                    'Write On Suspend?     '  & ! 4
                                  /),                         &
                                  (/ 4, 1 /)                  &
                                ),                            &
         Defaults             = (/                            &
                                  '   ', & ! 1
                                  '   ', & ! 2
                                  '   ', & ! 3
                                  '   '  & ! 4
                                /),                           &
         TwoD                 = .false.,                      &
         UnrecognisedMessages = .true.,                       &
         DisjointColSpecs     = .true.,                       &
         HFBlockForms         = HFBlockForms                  &
       )

End Subroutine RestartOptsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2RestartOpts(Tokens, MainOpts, RestartOpts)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_),      Intent(In)  :: Tokens
  Type(MainOpts_),    Intent(In)  :: MainOpts
  Type(RestartOpts_), Intent(Out) :: RestartOpts
  ! Local parameters:
  Integer, Parameter :: i_NCasesBetweenWrites = 1
  Integer, Parameter :: i_TimeBetweenWrites   = 2
  Integer, Parameter :: i_DeleteOldFiles      = 3
  Integer, Parameter :: i_WriteOnSuspend      = 4
  ! Locals:
  Integer      :: dCase
  Type(Time_)  :: T0
  Type(Time_)  :: dT
  Integer      :: RestartType
  Logical      :: DeleteOldFiles
  Logical      :: WriteOnSuspend

  If (Tokens%Tokens(i_NCasesBetweenWrites) /= ' ' .and. Tokens%Tokens(i_TimeBetweenWrites) /= ' ') Then
    Call Message('Error in Tokens2RestartOpts', 3)
  End If

  If (Tokens%Tokens(i_NCasesBetweenWrites) == ' ' .and. Tokens%Tokens(i_TimeBetweenWrites) == ' ') Then
    RestartType = 0
    DeleteOldFiles = .false.
  Else If (Tokens%Tokens(i_NCasesBetweenWrites) /= ' ') Then
    dCase = Token2Int(                                          &
              C         = Tokens%Tokens(i_NCasesBetweenWrites), &
              BlockKey  = 'Restart File Options',               &
              Item      = ' ',                                  &
              ColumnKey = '# Cases Between Writes'              &
            )
    RestartType = 1
    DeleteOldFiles = Token2Log(                                     &
                       C         = Tokens%Tokens(i_DeleteOldFiles), &
                       BlockKey  = 'Restart File Options',          &
                       Item      = ' ',                             &
                       ColumnKey = 'Delete Old Files?'              &
                     )
  Else If (Tokens%Tokens(i_TimeBetweenWrites) /= ' ') Then
    dT = Token2Time(                                       &
           C         = Tokens%Tokens(i_TimeBetweenWrites), &
           BlockKey  = 'Restart File Options',             &
           Item      = ' ',                                &
           ColumnKey = 'Time Between Writes',              &
           Interval  = .true.                              &
         )
    RestartType = 2
    DeleteOldFiles = Token2Log(                                     &
                       C         = Tokens%Tokens(i_DeleteOldFiles), &
                       BlockKey  = 'Restart File Options',          &
                       Item      = ' ',                             &
                       ColumnKey = 'Delete Old Files?'              &
                     )
  End If

  If (Tokens%Tokens(i_WriteOnSuspend) /= ' ') Then
    WriteOnSuspend = Token2Log(                                     &
                       C         = Tokens%Tokens(i_WriteOnSuspend), &
                       BlockKey  = 'Restart File Options',          &
                       Item      = ' ',                             &
                       ColumnKey = 'Write On Suspend?'              &
                     )
  Else
    WriteOnSuspend = .false.
  End If
  If (WriteOnSuspend .and. RestartType == 0) RestartType = 3

  T0 = ReferenceTime()

  Call InitRestartOpts(                   &
         RestartType    = RestartType,    &
         dCase          = dCase,          &
         T0             = T0,             &
         dT             = dT,             &
         DeleteOldFiles = DeleteOldFiles, &
         WriteOnSuspend = WriteOnSuspend, &
         RestartOpts    = RestartOpts     &
       )

End Subroutine Tokens2RestartOpts

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine MultiCaseOptsInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                           &
         BlockKey             = 'Multiple Case Options:',               &
         NamedBlock           = .false.,                                &
         MatchMultipleLines   = .false.,                                &
         nHeaderLines         = 1,                                      &
         ColumnKeys           = Reshape(                                &
                                  (/                                    &
                                    'Dispersion Options Ensemble Size', & ! 1
                                    'Met Ensemble Size               '  & ! 2
                                  /),                                   &
                                  (/ 2, 1 /)                            &
                                ),                                      &
         Defaults             = (/                                      &
                                  '1',                                  & ! 1
                                  '1'                                   & ! 2
                                /),                                     &
         TwoD                 = .false.,                                &
         UnrecognisedMessages = .true.,                                 &
         DisjointColSpecs     = .true.,                                 &
         HFBlockForms         = HFBlockForms                            &
       )

End Subroutine MultiCaseOptsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2MultiCaseOpts(Tokens, MultiCaseOpts)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_),        Intent(In)  :: Tokens
  Type(MultiCaseOpts_), Intent(Out) :: MultiCaseOpts
  ! Local parameters:
  Integer, Parameter :: i_DispersionOptionsEnsembleSize = 1
  Integer, Parameter :: i_MetEnsembleSize               = 2
  ! Locals:
  Integer :: DispOptsesEnsembleSize
  Integer :: MetEnsembleSize

  DispOptsesEnsembleSize = Token2Int(                                                    &
                             C         = Tokens%Tokens(i_DispersionOptionsEnsembleSize), &
                             BlockKey  = 'Multiple Case Options',                        &
                             Item      = ' ',                                            &
                             ColumnKey = 'Dispersion Options Ensemble Size'              &
                           )
  MetEnsembleSize        = Token2Int(                                      &
                             C         = Tokens%Tokens(i_MetEnsembleSize), &
                             BlockKey  = 'Multiple Case Options',          &
                             Item      = ' ',                              &
                             ColumnKey = 'Met Ensemble Size'               &
                           )

  Call InitMultiCaseOpts(                                 &
         DispOptsesEnsembleSize = DispOptsesEnsembleSize, &
         MetEnsembleSize        = MetEnsembleSize,        &
         MultiCaseOpts          = MultiCaseOpts           &
       )

End Subroutine Tokens2MultiCaseOpts

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine ChemOptsInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                 &
         BlockKey             = 'Chemistry Options:',         &
         NamedBlock           = .false.,                      &
         MatchMultipleLines   = .false.,                      &
         nHeaderLines         = 1,                            &
         ColumnKeys           = Reshape(                      &
                                  (/                          &
                                    'Chemistry Folder      '  & ! 1
                                  /),                         &
                                  (/ 1, 1 /)                  &
                                ),                            &
         Defaults             = (/                            &
                                  '../Resources/Stochem/'     & ! 1
                                /),                           &
         TwoD                 = .false.,                      &
         UnrecognisedMessages = .true.,                       &
         DisjointColSpecs     = .true.,                       &
         HFBlockForms         = HFBlockForms                  &
       )

End Subroutine ChemOptsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2ChemOpts(Tokens, ChemOpts)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)  :: Tokens
  Type(ChemOpts_), Intent(Out) :: ChemOpts
  ! Local parameters:
  Integer, Parameter :: i_ChemistryFolder = 1

  Call InitChemOpts(                                     &
         ChemFolder  = Tokens%Tokens(i_ChemistryFolder), &
         ChemOpts    = ChemOpts                          &
       )

End Subroutine Tokens2ChemOpts

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine OutputOptsInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                              &
         BlockKey             = 'Output Options:',         &
         NamedBlock           = .false.,                   &
         MatchMultipleLines   = .false.,                   &
         nHeaderLines         = 1,                         &
         ColumnKeys           = Reshape(                   &
                                  (/                       &
                                    'Folder             ', & ! 1
                                    'Seconds?           ', & ! 2
                                    'Time Decimal Places', & ! 3
                                    'Pre 6.5 Format?    '  & ! 4 $$ remove this when not needed
                                  /),                      &
                                  (/ 4, 1 /)               &
                                ),                         &
         Defaults             = (/                         &
                                  '           ', & ! 1
                                  '           ', & ! 2
                                  '           ', & ! 3
                                  'Yes        '  & ! 4
                                /),                        &
         TwoD                 = .false.,                   &
         UnrecognisedMessages = .true.,                    &
         DisjointColSpecs     = .true.,                    &
         HFBlockForms         = HFBlockForms               &
       )

End Subroutine OutputOptsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2OutputOpts(Tokens, OutputOpts)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),     Intent(In)    :: Tokens
  Type(OutputOpts_), Intent(InOut) :: OutputOpts
  ! Local parameters:
  Integer, Parameter :: i_Folder            = 1
  Integer, Parameter :: i_Seconds           = 2
  Integer, Parameter :: i_TimeDecimalPlaces = 3
  Integer, Parameter :: i_Pre65Format       = 4 ! $$ remove this when not needed
  ! Locals:
  Logical :: Seconds
  Integer :: DecimalPlaces
  Logical :: Pre65Format

  If (Tokens%Tokens(i_Seconds) /= ' ') Then
    Seconds = Token2Log(                              &
                C         = Tokens%Tokens(i_Seconds), &
                BlockKey  = 'Output Options',         &
                Item      = ' ',                      &
                ColumnKey = 'Seconds?'                &
              )
  Else
    Seconds = .false.
  End If

  If (Seconds) Then
    If (Tokens%Tokens(i_TimeDecimalPlaces) /= ' ') Then
      DecimalPlaces = Token2Int(                                        &
                        C         = Tokens%Tokens(i_TimeDecimalPlaces), &
                        BlockKey  = 'Output Options',                   &
                        Item      = ' ',                                &
                        ColumnKey = 'Time Decimal Places'               &
                      )
    Else
      DecimalPlaces = 0
    End If
  End If

  Pre65Format = Token2Log(                                  &
                  C         = Tokens%Tokens(i_Pre65Format), &
                  BlockKey  = 'Output Options',             &
                  Item      = ' ',                          &
                  ColumnKey = 'Pre 6.5 Format?'             &
                )
    
  Call InitOutputOpts(Tokens%Tokens(i_Folder), Seconds, DecimalPlaces, Pre65Format, OutputOpts)

End Subroutine Tokens2OutputOpts

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine OpenMPOptsInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                     &
         BlockKey             = 'OpenMP Options:',                &
         NamedBlock           = .false.,                          &
         MatchMultipleLines   = .false.,                          &
         nHeaderLines         = 1,                                &
         ColumnKeys           = Reshape(                          &
                                  (/                              &
                                    'Use OpenMP?               ', & ! 1
                                    'Threads                   ', & ! 2
                                    'Particle Threads          ', & ! 3
                                    'Particle Update Threads   ', & ! 4
                                    'Chemistry Threads         ', & ! 5
                                    'Output Group Threads      ', & ! 6
                                    'Output Process Threads    ', & ! 7
                                    'Parallel MetRead          ', & ! 8
                                    'Parallel MetProcess       '  & ! 9
                                  /),                             &
                                  (/ 9, 1 /)                      &
                                ),                                &
         Defaults             = (/                                &
                                  'No                          ', & ! 1
                                  '1                           ', & ! 2
                                  '                            ', & ! 3
                                  '                            ', & ! 4
                                  '                            ', & ! 5
                                  '                            ', & ! 6
                                  '                            ', & ! 7
                                  'No                          ', & ! 8
                                  '                            '  & ! 9
                                /),                               &
         TwoD                 = .false.,                          &
         UnrecognisedMessages = .true.,                           &
         DisjointColSpecs     = .true.,                           &
         HFBlockForms         = HFBlockForms                      &
       )

End Subroutine OpenMPOptsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2OpenMPOpts(Tokens, OpenMPOpts)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),     Intent(In)    :: Tokens
  Type(OpenMPOpts_), Intent(InOut) :: OpenMPOpts
  ! Local parameters:
  Integer, Parameter :: i_UseOpenMP             = 1
  Integer, Parameter :: i_Threads               = 2
  Integer, Parameter :: i_ParticleThreads       = 3
  Integer, Parameter :: i_ParticleUpdateThreads = 4
  Integer, Parameter :: i_ChemistryThreads      = 5
  Integer, Parameter :: i_OutputGroupThreads    = 6
  Integer, Parameter :: i_OutputProcessThreads  = 7
  Integer, Parameter :: i_ParallelMetRead       = 8
  Integer, Parameter :: i_ParallelMetProcess    = 9
  ! Locals:
  Integer :: nThreads

  If (Tokens%Tokens(i_UseOpenMP) /= ' ') Then
    OpenMPOpts%UseOpenMP = Token2Log(                                &
                             C         = Tokens%Tokens(i_UseOpenMP), &
                             BlockKey  = 'OpenMP Options',           &
                             Item      = ' ',                        &
                             ColumnKey = 'Use OpenMP?'               &
                           )
  Else
    OpenMPOpts%UseOpenMP = .false.
  End If

  If (Tokens%Tokens(i_Threads) /= ' ') Then
    nThreads = Token2Int(                              &
                 C         = Tokens%Tokens(i_Threads), &
                 BlockKey  = 'OpenMP Options',         &
                 Item      = ' ',                      &
                 ColumnKey = 'Threads'                 &
               )
  Else
    nThreads = 1
  End If

  OpenMPOpts%nParticleThreads=nThreads
  OpenMPOpts%nParticleUpdateThreads=nThreads
  OpenMPOpts%nChemistryThreads=nThreads
  OpenMPOpts%nOutputGroupThreads=nThreads
  OpenMPOpts%nOutputProcessThreads=nThreads

  If (Tokens%Tokens(i_ParticleThreads) /= ' ') Then
    OpenMPOpts%nParticleThreads = Token2Int(                                      &
                                    C         = Tokens%Tokens(i_ParticleThreads), &
                                    BlockKey  = 'OpenMP Options',                 &
                                    Item      = ' ',                              &
                                    ColumnKey = 'Particle Threads'                &
                                  )
  End If

  If (Tokens%Tokens(i_ParticleUpdateThreads) /= ' ') Then
    OpenMPOpts%nParticleUpdateThreads = Token2Int(                                            &
                                          C         = Tokens%Tokens(i_ParticleUpdateThreads), &
                                          BlockKey  = 'OpenMP Options',                       &
                                          Item      = ' ',                                    &
                                          ColumnKey = 'Particle Update Threads'               &
                                  )
  End If

  If (Tokens%Tokens(i_ChemistryThreads) /= ' ') Then
    OpenMPOpts%nChemistryThreads = Token2Int(                                       &
                                     C         = Tokens%Tokens(i_ChemistryThreads), &
                                     BlockKey  = 'OpenMP Options',                  &
                                     Item      = ' ',                               &
                                     ColumnKey = 'Chemistry Threads'                &
                                   )
  End If

  If (Tokens%Tokens(i_OutputGroupThreads) /= ' ') Then
    OpenMPOpts%nOutputGroupThreads = Token2Int(                                         &
                                       C         = Tokens%Tokens(i_OutputGroupThreads), &
                                       BlockKey  = 'OpenMP Options',                    &
                                       Item      = ' ',                                 &
                                       ColumnKey = 'Output Group Threads'               &
                                     )
  End If
  
  If (Tokens%Tokens(i_OutputProcessThreads) /= ' ') Then
    OpenMPOpts%nOutputProcessThreads = Token2Int(                                           &
                                         C         = Tokens%Tokens(i_OutputProcessThreads), &
                                         BlockKey  = 'OpenMP Options',                      &
                                         Item      = ' ',                                   &
                                         ColumnKey = 'Output Process Threads'               &
                                       )
  End If

  If (Tokens%Tokens(i_ParallelMetRead) /= ' ') Then
    OpenMPOpts%ParallelMetRead = Token2Log(                                      &
                                   C         = Tokens%Tokens(i_ParallelMetRead), &
                                   BlockKey  = 'OpenMP Options',                 &
                                   Item      = ' ',                              &
                                   ColumnKey = 'Parallel MetRead'                &
                                 )
  Else
    OpenMPOpts%ParallelMetRead = .false.
  End If

  OpenMPOpts%ParallelMetProcess = OpenMPOpts%ParallelMetRead

  If (Tokens%Tokens(i_ParallelMetProcess) /= ' ') Then
    OpenMPOpts%ParallelMetProcess = Token2Log(                                         &
                                      C         = Tokens%Tokens(i_ParallelMetProcess), &
                                      BlockKey  = 'OpenMP Options',                    &
                                      Item      = ' ',                                 &
                                      ColumnKey = 'Parallel MetProcess'                &
                                    )
  End If

  Call InitOpenMPOpts(OpenMPOpts)

End Subroutine Tokens2OpenMPOpts

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine TimerOptsInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                     &
         BlockKey             = 'Timer Options:',                 &
         NamedBlock           = .false.,                          &
         MatchMultipleLines   = .false.,                          &
         nHeaderLines         = 1,                                &
         ColumnKeys           = Reshape(                          &
                                  (/                              &
                                    'Use Timers?               ', & ! 1
                                    'Summary Only?             '  & ! 2
                                  /),                             &
                                  (/ 2, 1 /)                      &
                                ),                                &
         Defaults             = (/                                &
                                  'No                          ', & ! 1
                                  'Yes                         '  & ! 2
                                /),                               &
         TwoD                 = .false.,                          &
         UnrecognisedMessages = .true.,                           &
         DisjointColSpecs     = .true.,                           &
         HFBlockForms         = HFBlockForms                      &
       )

End Subroutine TimerOptsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2TimerOpts(Tokens, TimerOpts, Units, MainInputFile)
!.

  Implicit None

  ! Argument list:
  Type(Tokens_),      Intent(In)              :: Tokens
  Type(TimerOpts_),   Intent(InOut)           :: TimerOpts
  Type(Units_),       Intent(InOut)           :: Units
  Character(MaxFileNameLength), Intent(In)    :: MainInputFile
  ! Local parameters:
  Integer, Parameter :: i_UseTimers   = 1
  Integer, Parameter :: i_SummaryOnly = 2
  ! Locals:
  Integer                      :: nThreads
  Character(MaxFileNameLength) :: TimingFile
  Integer                      :: i, iDot

  If (Tokens%Tokens(i_UseTimers) /= ' ') Then
    TimerOpts%UseTimers = Token2Log(                                &
                            C         = Tokens%Tokens(i_UseTimers), &
                            BlockKey  = 'Timer Options',            &
                            Item      = ' ',                        &
                            ColumnKey = 'Use Timers?'               &
                          )
  Else
    TimerOpts%UseTimers = .false.
  End If

  If (Tokens%Tokens(i_SummaryOnly) /= ' ') Then
    TimerOpts%SummaryOnly = Token2Log(                                  &
                              C         = Tokens%Tokens(i_SummaryOnly), &
                              BlockKey  = 'Timer Options',              &
                              Item      = ' ',                          &
                              ColumnKey = 'Summary Only?'               &
                            )
  Else
    TimerOpts%SummaryOnly = .true.
  End If

  iDot = Scan(MainInputFile, '.', Back = .true.)
  i    = Scan(MainInputFile, '\/', Back = .true.)

  ! Determine the timing file.
  If (iDot == 0 .or. iDot < i) Then
    TimingFile = Trim(MainInputFile) // 'Timing.txt'
  Else
    TimingFile = MainInputFile(1:iDot - 1) // 'Timing.txt'
  End If

  TimingFile = ConvertFileName(TimingFile)

  Call InitTimerOpts(TimerOpts, Units, Trim(TimingFile))

End Subroutine Tokens2TimerOpts

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine InputFileInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                    &
         BlockKey             = 'Input Files:',  &
         NamedBlock           = .false.,         &
         MatchMultipleLines   = .false.,         &
         nHeaderLines         = 1,               &
         ColumnKeys           = Reshape(         &
                                  (/             &
                                    'File Names' & ! 1
                                  /),            &
                                  (/ 1, 1 /)     &
                                ),               &
         Defaults             = (/               &
                                  ' '            & ! 1
                                /),              &
         TwoD                 = .false.,         &
         UnrecognisedMessages = .true.,          &
         DisjointColSpecs     = .true.,          &
         HFBlockForms         = HFBlockForms     &
       )

End Subroutine InputFileInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2InputFiles(Tokens, InputFiles)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_),     Intent(In)    :: Tokens
  Type(InputFiles_), Intent(InOut) :: InputFiles
  ! Local parameters:
  Integer, Parameter :: i_FileNames = 1

  Call AddInputFile(                                             &
         InputFile  = InitInputFile(Tokens%Tokens(i_FileNames)), &
         InputFiles = InputFiles                                 &
       )

End Subroutine Tokens2InputFiles

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine ArrayInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                      &
         BlockKey             = 'Array:',          &
         NamedBlock           = .true.,            &
         MatchMultipleLines   = .false.,           &
         nHeaderLines         = 1,                 &
         ColumnKeys           = Reshape(           &
                                  (/               &
                                    'Array Values' & ! 1
                                  /),              &
                                  (/ 1, 1 /)       &
                                ),                 &
         Defaults             = (/                 &
                                  ' '              & ! 1
                                /),                &
         TwoD                 = .true.,            &
         UnrecognisedMessages = .true.,            &
         DisjointColSpecs     = .true.,            &
         HFBlockForms         = HFBlockForms       &
       )

End Subroutine ArrayInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2Array(Tokens, Arrays)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_), Intent(In)    :: Tokens
  Type(Arrays_), Intent(InOut) :: Arrays
  ! Local parameters:
  Integer, Parameter :: i_ArrayValues = 1
  ! Locals:
  Character(MaxTokenLength) :: ArrayValues(MaxArrayLength) ! Array entries.
  Type(Array_)              :: Array                       ! Local copy of array.

  If (Tokens%nLines > MaxArrayLength) Then
    Call Message('ERROR in Tokens2Array: too many entries in the array ' // &
                 Trim(Tokens%BlockName),                                    &
                 3)
  End If

  ArrayValues(1:Tokens%nLines) = Tokens%Tokens2d(i_ArrayValues, 1:Tokens%nLines)

  Array = InitArray(                                      &
            Name          = Tokens%BlockName,             &
            ArrayElements = ArrayValues(1:Tokens%nLines), &
            BlockKey      = 'Array'                       &
          )

  Call AddArray(Array, Arrays)

End Subroutine Tokens2Array

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine MessageControlInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                           &
         BlockKey             = 'Message Controls:',                    &
         NamedBlock           = .false.,                                &
         MatchMultipleLines   = .false.,                                &
         nHeaderLines         = 1,                                      &
         ColumnKeys           = Reshape(                                &
                                  (/                                    &
                                    'Name                            ', & ! 1
                                    'Action                          ', & ! 2
                                    '# of Messages Before Failure    ', & ! 3
                                    '# of Messages Before Suppression'  & ! 4
                                  /),                                   &
                                  (/ 4, 1 /)                            &
                                ),                                      &
         Defaults             = (/                                      &
                                  ' ',                                  & ! 1
                                  ' ',                                  & ! 2
                                  ' ',                                  & ! 3
                                  ' '                                   & ! 4
                                /),                                     &
         TwoD                 = .false.,                                &
         UnrecognisedMessages = .true.,                                 &
         DisjointColSpecs     = .true.,                                 &
         HFBlockForms         = HFBlockForms                            &
       )

End Subroutine MessageControlInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2MessageControl(Tokens, MessageControls)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_),          Intent(In)    :: Tokens
  Type(MessageControls_), Intent(InOut) :: MessageControls
  ! Local parameters:
  Integer, Parameter :: i_Name                       = 1
  Integer, Parameter :: i_Action                     = 2
  Integer, Parameter :: i_NMessagesBeforeFailure     = 3
  Integer, Parameter :: i_NMessagesBeforeSuppression = 4
  ! Locals:
  Integer :: Fail
  Integer :: Suppress

  If (Tokens%Tokens(i_NMessagesBeforeFailure) /= ' ') Then
    Fail = Token2Int(                                             &
             C         = Tokens%Tokens(i_NMessagesBeforeFailure), &
             BlockKey  = 'Message Controls',                      &
             Item      = Tokens%Tokens(i_Name),                   &
             ColumnKey = '# of Messages Before Failure'           &
           )
  End If

  If (Tokens%Tokens(i_NMessagesBeforeSuppression) /= ' ') Then
    Suppress = Token2Int(                                                 &
                 C         = Tokens%Tokens(i_NMessagesBeforeSuppression), &
                 BlockKey  = 'Message Controls',                          &
                 Item      = Tokens%Tokens(i_Name),                       &
                 ColumnKey = '# of Messages Before Suppression'           &
               )
  End If

  Call AddMessageControl(                                                   &
         InitMessageControl(                                                &
           Name       = Tokens%Tokens(i_Name),                              &
           Action     = Tokens%Tokens(i_Action),                            &
           NoFail     = Tokens%Tokens(i_NMessagesBeforeFailure) == ' ',     &
           Fail       = Fail,                                               &
           NoSuppress = Tokens%Tokens(i_NMessagesBeforeSuppression) == ' ', &
           Suppress   = Suppress                                            &
         ),                                                                 &
         MessageControls                                                    &
       )

End Subroutine Tokens2MessageControl

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine EtaDefnInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                      &
         BlockKey             = 'Eta Definition:', &
         NamedBlock           = .true.,            &
         MatchMultipleLines   = .false.,           &
         nHeaderLines         = 1,                 &
         ColumnKeys           = Reshape(           &
                                  (/               &
                                    'Eta',         & ! 1
                                    'A  ',         & ! 2
                                    'B  '          & ! 3
                                  /),              &
                                  (/ 3, 1 /)       &
                                ),                 &
         Defaults             = (/                 &
                                  ' ',             & ! 1
                                  ' ',             & ! 2
                                  ' '              & ! 3
                                /),                &
         TwoD                 = .true.,            &
         UnrecognisedMessages = .true.,            &
         DisjointColSpecs     = .true.,            &
         HFBlockForms         = HFBlockForms       &
       )

End Subroutine EtaDefnInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2EtaDefn(Tokens, EtaDefns)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)    :: Tokens
  Type(EtaDefns_), Intent(InOut) :: EtaDefns
  ! Local parameters:
  Integer, Parameter :: i_Eta = 1
  Integer, Parameter :: i_A   = 2
  Integer, Parameter :: i_B   = 3
  ! Locals:
  Integer        :: i                 ! Loop index.
  Real(Std)      :: Eta(MaxLinesTwoD) ! Eta values.
  Real(Std)      :: A(MaxLinesTwoD)   ! A values.
  Real(Std)      :: B(MaxLinesTwoD)   ! B values.
  Type(EtaDefn_) :: EtaDefn           ! Local copy of eta definition.
  Character(1)   :: EtaDefnType

  ! Eta and A given.
  If (                                                                                                         &
    Tokens%Tokens2d(i_Eta, 1) /= ' ' .and. Tokens%Tokens2d(i_A, 1) /= ' ' .and. Tokens%Tokens2d(i_B, 1) == ' ' &
  ) Then
    EtaDefnType = 'B'
  ! Eta and B given.
  Else If (                                                                                                    &
    Tokens%Tokens2d(i_Eta, 1) /= ' ' .and. Tokens%Tokens2d(i_A, 1) == ' ' .and. Tokens%Tokens2d(i_B, 1) /= ' ' &
  ) Then
    EtaDefnType = 'A'
  ! A and B given.
  Else If (                                                                                                    &
    Tokens%Tokens2d(i_Eta, 1) == ' ' .and. Tokens%Tokens2d(i_A, 1) /= ' ' .and. Tokens%Tokens2d(i_B, 1) /= ' ' &
  ) Then
    EtaDefnType = 'E'
  Else
    Call Message(                                                            &
           'FATAL ERROR in Tokens2EtaDefn: precisely two of Eta, A and B' // &
           ' must be given in an Eta Definition',                            &
           3                                                                 &
         )
  End If

  If (EtaDefnType == 'E') Then

    Do i = 1, Tokens%nLines
      ! Check same two variables out of Eta, A and B present on each line $$
      A(i) = Token2Std(                             &
               C         = Tokens%Tokens2d(i_A, i), &
               BlockKey  = 'Eta Definition',        &
               Item      = Tokens%BlockName,        &
               ColumnKey = 'A',                     &
               i         = i                        &
             )
      B(i) = Token2Std(                             &
               C         = Tokens%Tokens2d(i_B, i), &
               BlockKey  = 'Eta Definition',        &
               Item      = Tokens%BlockName,        &
               ColumnKey = 'B',                     &
               i         = i                        &
             )
    End Do

    EtaDefn = InitEtaDefn(                       &
                Name       = Tokens%BlockName,   &
                nEtaLevels = Tokens%nLines,      &
                A          = A(1:Tokens%nLines), &
                B          = B(1:Tokens%nLines)  &
              )

  Else If (EtaDefnType == 'A') Then

    Do i = 1, Tokens%nLines
      ! Check same two variables out of Eta, A and B present on each line $$
      Eta(i) = Token2Std(                               &
                 C         = Tokens%Tokens2d(i_Eta, i), &
                 BlockKey  = 'Eta Definition',          &
                 Item      = Tokens%BlockName,          &
                 ColumnKey = 'Eta',                     &
                 i         = i                          &
               )
      B(i)   = Token2Std(                               &
                 C         = Tokens%Tokens2d(i_B,   i), &
                 BlockKey  = 'Eta Definition',          &
                 Item      = Tokens%BlockName,          &
                 ColumnKey = 'B',                       &
                 i         = i                          &
               )
    End Do

    EtaDefn = InitEtaDefn(                         &
                Name       = Tokens%BlockName,     &
                nEtaLevels = Tokens%nLines,        &
                Eta        = Eta(1:Tokens%nLines), &
                B          = B  (1:Tokens%nLines)  &
              )

  Else If (EtaDefnType == 'B') Then

    Do i = 1, Tokens%nLines
      ! Check same two variables out of Eta, A and B present on each line $$
      Eta(i) = Token2Std(                               &
                 C         = Tokens%Tokens2d(i_Eta, i), &
                 BlockKey  = 'Eta Definition',          &
                 Item      = Tokens%BlockName,          &
                 ColumnKey = 'Eta',                     &
                 i         = i                          &
               )
      A(i)   = Token2Std(                               &
                 C         = Tokens%Tokens2d(i_A,   i), &
                 BlockKey  = 'Eta Definition',          &
                 Item      = Tokens%BlockName,          &
                 ColumnKey = 'A',                       &
                 i         = i                          &
               )
    End Do

    EtaDefn = InitEtaDefn(                         &
                Name       = Tokens%BlockName,     &
                nEtaLevels = Tokens%nLines,        &
                Eta        = Eta(1:Tokens%nLines), &
                A          = A(1:Tokens%nLines)    &
              )

  End If

  Call AddEtaDefn(EtaDefn, EtaDefns)

End Subroutine Tokens2EtaDefn

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine LocationsInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                  &
         BlockKey             = 'Locations:',  &
         NamedBlock           = .true.,        &
         MatchMultipleLines   = .false.,       &
         nHeaderLines         = 1,             &
         ColumnKeys           = Reshape(       &
                                  (/           &
                                    'Name   ', & ! 1
                                    'H-Coord', & ! 2
                                    'X      ', & ! 3
                                    'Y      '  & ! 4
                                  /),          &
                                  (/ 4, 1 /)   &
                                ),             &
         Defaults             = (/             &
                                  '   ',  & ! 1
                                  '   ',  & ! 2
                                  '   ',  & ! 3
                                  '   '   & ! 4
                                /),            &
         TwoD                 = .true.,        &
         UnrecognisedMessages = .true.,        &
         DisjointColSpecs     = .true.,        &
         HFBlockForms         = HFBlockForms   &
       )

End Subroutine LocationsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2Locations(Tokens, Locationses)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_),      Intent(In)    :: Tokens
  Type(Locationses_), Intent(InOut) :: Locationses
  ! Local parameters:
  Integer, Parameter :: i_Name   = 1
  Integer, Parameter :: i_HCoord = 2
  Integer, Parameter :: i_X      = 3
  Integer, Parameter :: i_Y      = 4
  ! Locals:
  Integer    :: i               ! Loop index.
  Real(Std)  :: X(MaxLinesTwoD) !
  Real(Std)  :: Y(MaxLinesTwoD) !

  Do i = 1, Tokens%nLines
    X(i) = Token2Std(                             &
             C         = Tokens%Tokens2d(i_X, i), &
             BlockKey  = 'Locations',             &
             Item      = Tokens%BlockName,        &
             ColumnKey = 'X',                     &
             i         = i                        &
           )
    Y(i) = Token2Std(                             &
             C         = Tokens%Tokens2d(i_Y, i), &
             BlockKey  = 'Locations',             &
             Item      = Tokens%BlockName,        &
             ColumnKey = 'Y',                     &
             i         = i                        &
           )
  End Do

  Call AddLocations(                                                 &
         InitLocations(                                              &
           Name        = Tokens%BlockName,                           &
           Names       = Tokens%Tokens2d(i_Name,   1:Tokens%nLines), &
           HCoordNames = Tokens%Tokens2d(i_HCoord, 1:Tokens%nLines), &
           X           = X(1:Tokens%nLines),                         &
           Y           = Y(1:Tokens%nLines)                          &
         ),                                                          &
         Locationses                                                 &
       )

End Subroutine Tokens2Locations

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine HCoordInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                     &
         BlockKey             = 'Horizontal Coordinate Systems:', &
         NamedBlock           = .false.,                          &
         MatchMultipleLines   = .false.,                          &
         nHeaderLines         = 1,                                &
         ColumnKeys           = Reshape(                          &
                                  (/                              &
                                    'Name        ',               & ! 1
                                    'Type        ',               & ! 2
                                    'Pole Long   ',               & ! 3
                                    'Pole Lat    ',               & ! 4
                                    'Angle       ',               & ! 5
                                    'X-Origin    ',               & ! 6
                                    'Y-Origin    ',               & ! 7
                                    'X-Unit      ',               & ! 8
                                    'Y-Unit      ',               & ! 9
                                    'Scale Factor'                & ! 10
                                  /),                             &
                                  (/ 10, 1 /)                     &
                                ),                                &
         Defaults             = (/                                &
                                  '   ',                 & ! 1
                                  '   ',                 & ! 2
                                  '   ',                 & ! 3
                                  '   ',                 & ! 4
                                  '   ',                 & ! 5
                                  '   ',                 & ! 6
                                  '   ',                 & ! 7
                                  '   ',                 & ! 8
                                  '   ',                 & ! 9
                                  '   '                  & ! 10
                                /),                               &
         TwoD                 = .false.,                          &
         UnrecognisedMessages = .true.,                           &
         DisjointColSpecs     = .true.,                           &
         HFBlockForms         = HFBlockForms                      &
       )

End Subroutine HCoordInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2HCoord(Tokens, Coords)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_), Intent(In)    :: Tokens
  Type(Coords_), Intent(InOut) :: Coords
  ! Local parameters:
  Integer, Parameter :: i_Name        = 1
  Integer, Parameter :: i_Type        = 2
  Integer, Parameter :: i_PoleLong    = 3
  Integer, Parameter :: i_PoleLat     = 4
  Integer, Parameter :: i_Angle       = 5
  Integer, Parameter :: i_XOrigin     = 6
  Integer, Parameter :: i_YOrigin     = 7
  Integer, Parameter :: i_XUnit       = 8
  Integer, Parameter :: i_YUnit       = 9
  Integer, Parameter :: i_ScaleFactor = 10
  ! Locals:
  Integer       :: iTokens     !
  Integer       :: CoordType   !} Temporary variables, used to convert Tokens into the
  Real(Std)     :: Pole(2)     !} right form.
  Real(Std)     :: Angle       !}
  Real(Std)     :: Origin(2)   !}
  Real(Std)     :: Unit(2)     !}
  Real(Std)     :: ScaleFactor !}
  Type(HCoord_) :: HCoord      !

  ! Standard coord systems.
  If (Tokens%Tokens(i_Name) .CIEq. 'Lat-Long') Then
    Do iTokens = 2, 10
      If (Tokens%Tokens(iTokens) /= ' ') Then
        Call Message(                                                        &
               'FATAL ERROR in processing input: the coordinate system "' // &
               Trim(Tokens%Tokens(i_Name))                                // &
               '" is a standard coordinate system and must be '           // &
               'defined by name only',                                       &
               3                                                             &
             )
      End If
    End Do
    HCoord = HCoord_LatLong()
  Else If (Tokens%Tokens(i_Name) .CIEq. 'EMEP 50km Grid') Then
    Do iTokens = 2, 10
      If (Tokens%Tokens(iTokens) /= ' ') Then
        Call Message(                                                        &
               'FATAL ERROR in processing input: the coordinate system "' // &
               Trim(Tokens%Tokens(i_Name))                                // &
               '" is a standard coordinate system and must be '           // &
               'defined by name only',                                       &
               3                                                             &
             )
      End If
    End Do
    HCoord = HCoord_EMEP50kmGrid()
  Else If (Tokens%Tokens(i_Name) .CIEq. 'EMEP 150km Grid') Then
    Do iTokens = 2, 10
      If (Tokens%Tokens(iTokens) /= ' ') Then
        Call Message(                                                        &
               'FATAL ERROR in processing input: the coordinate system "' // &
               Trim(Tokens%Tokens(i_Name))                                // &
               '" is a standard coordinate system and must be '           // &
               'defined by name only',                                       &
               3                                                             &
             )
      End If
    End Do
    HCoord = HCoord_EMEP150kmGrid()
  Else If (Tokens%Tokens(i_Name) .CIEq. 'UK National Grid (m)') Then
    Do iTokens = 2, 10
      If (Tokens%Tokens(iTokens) /= ' ') Then
        Call Message(                                                        &
               'FATAL ERROR in processing input: the coordinate system "' // &
               Trim(Tokens%Tokens(i_Name))                                // &
               '" is a standard coordinate system and must be '           // &
               'defined by name only',                                       &
               3                                                             &
             )
      End If
    End Do
    HCoord = HCoord_UKNationalGridM()
  Else If (Tokens%Tokens(i_Name) .CIEq. 'UK National Grid (100m)') Then
    Do iTokens = 2, 10
      If (Tokens%Tokens(iTokens) /= ' ') Then
        Call Message(                                                        &
               'FATAL ERROR in processing input: the coordinate system "' // &
               Trim(Tokens%Tokens(i_Name))                                // &
               '" is a standard coordinate system and must be '           // &
               'defined by name only',                                       &
               3                                                             &
             )
      End If
    End Do
    HCoord = HCoord_UKNationalGrid100M()

  ! Non standard coord systems.
  Else

    CoordType = Token2Int(                                     &
                  C         = Tokens%Tokens(i_Type),           &
                  BlockKey  = 'Horizontal Coordinate Systems', &
                  Item      = Tokens%Tokens(i_Name),           &
                  ColumnKey = 'Type'                           &
                )
    Pole(1)   = Token2Std(                                     &
                  C         = Tokens%Tokens(i_PoleLong),       &
                  BlockKey  = 'Horizontal Coordinate Systems', &
                  Item      = Tokens%Tokens(i_Name),           &
                  ColumnKey = 'Pole Long'                      &
                )
    Pole(2)   = Token2Std(                                     &
                  C         = Tokens%Tokens(i_PoleLat),        &
                  BlockKey  = 'Horizontal Coordinate Systems', &
                  Item      = Tokens%Tokens(i_Name),           &
                  ColumnKey = 'Pole Lat'                       &
                )
    Angle     = Token2Std(                                     &
                  C         = Tokens%Tokens(i_Angle),          &
                  BlockKey  = 'Horizontal Coordinate Systems', &
                  Item      = Tokens%Tokens(i_Name),           &
                  ColumnKey = 'Angle'                          &
                )
    Origin(1) = Token2Std(                                     &
                  C         = Tokens%Tokens(i_XOrigin),        &
                  BlockKey  = 'Horizontal Coordinate Systems', &
                  Item      = Tokens%Tokens(i_Name),           &
                  ColumnKey = 'X-Origin'                       &
                )
    Origin(2) = Token2Std(                                     &
                  C         = Tokens%Tokens(i_YOrigin),        &
                  BlockKey  = 'Horizontal Coordinate Systems', &
                  Item      = Tokens%Tokens(i_Name),           &
                  ColumnKey = 'Y-Origin'                       &
                )
    Unit(1)   = Token2Std(                                     &
                  C         = Tokens%Tokens(i_XUnit),          &
                  BlockKey  = 'Horizontal Coordinate Systems', &
                  Item      = Tokens%Tokens(i_Name),           &
                  ColumnKey = 'X-Unit'                         &
                )
    Unit(2)   = Token2Std(                                     &
                  C         = Tokens%Tokens(i_YUnit),          &
                  BlockKey  = 'Horizontal Coordinate Systems', &
                  Item      = Tokens%Tokens(i_Name),           &
                  ColumnKey = 'Y-Unit'                         &
                )
    If (Tokens%Tokens(i_ScaleFactor) == ' ') Then
      ScaleFactor = 1.0
    Else
      ScaleFactor = Token2Std(                                     &
                      C         = Tokens%Tokens(i_ScaleFactor),    &
                      BlockKey  = 'Horizontal Coordinate Systems', &
                      Item      = Tokens%Tokens(i_Name),           &
                      ColumnKey = 'Scale Factor'                   &
                    )
    End If
    HCoord = InitHCoord(                            &
               Name        = Tokens%Tokens(i_Name), &
               CoordType   = CoordType,             &
               Pole        = Pole,                  &
               Angle       = Angle,                 &
               Origin      = Origin,                &
               Unit        = Unit,                  &
               ThetaOrigin = 0.0,                   & ! $$
               ScaleFactor = ScaleFactor            &
             )

  End If

  Call AddHCoord(HCoord, Coords)

End Subroutine Tokens2HCoord

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine ZCoordInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                   &
         BlockKey             = 'Vertical Coordinate Systems:', &
         NamedBlock           = .false.,                        &
         MatchMultipleLines   = .false.,                        &
         nHeaderLines         = 1,                              &
         ColumnKeys           = Reshape(                        &
                                  (/                            &
                                    'Name            ',         & ! 1
                                    'Type            ',         & ! 2
                                    'Unit            ',         & ! 3
                                    'Eta Definition  ',         & ! 4
                                    'Model Top Height',         & ! 5
                                    'Interface Height'          & ! 6
                                  /),                           &
                                  (/ 6, 1 /)                    &
                                ),                              &
         Defaults             = (/                              &
                                  '  ',           & ! 1
                                  '  ',           & ! 2
                                  '  ',           & ! 3
                                  '  ',           & ! 4
                                  '  ',           & ! 5
                                  '  '            & ! 6
                                /),                             &
         TwoD                 = .false.,                        &
         UnrecognisedMessages = .true.,                         &
         DisjointColSpecs     = .true.,                         &
         HFBlockForms         = HFBlockForms                    &
       )

End Subroutine ZCoordInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2ZCoord(Tokens, EtaDefns, Coords)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)    :: Tokens
  Type(EtaDefns_), Intent(In)    :: EtaDefns
  Type(Coords_),   Intent(InOut) :: Coords
  ! Local parameters:
  Integer, Parameter :: i_Name            = 1 
  Integer, Parameter :: i_Type            = 2
  Integer, Parameter :: i_Unit            = 3
  Integer, Parameter :: i_EtaDefinition   = 4
  Integer, Parameter :: i_ModelTopHeight  = 5
  Integer, Parameter :: i_InterfaceHeight = 6
  ! Locals:
  Integer       :: CoordType       !} Temporary variables, used to convert Tokens into
  Real(Std)     :: Unit            !} the right form.
  Real(Std)     :: ModelTopHeight  !}
  Real(Std)     :: InterfaceHeight !}
  Type(ZCoord_) :: ZCoord          !
  Integer       :: i               !
  Integer       :: iTokens

  ! $$ Check for blanks

  ! Standard coord systems.
  If (Tokens%Tokens(i_Name) .CIEq. 'm agl') Then
    Do iTokens = 2, 6
      If (Tokens%Tokens(iTokens) /= ' ') Then
        Call Message(                                                        &
               'FATAL ERROR in processing input: the coordinate system "' // &
               Trim(Tokens%Tokens(i_Name))                                // &
               '" is a standard coordinate system and must be '           // &
               'defined by name only',                                       &
               3                                                             &
             )
      End If
    End Do
    ZCoord = ZCoord_m_agl()
  Else If (Tokens%Tokens(i_Name) .CIEq. 'm asl') Then
    Do iTokens = 2, 6
      If (Tokens%Tokens(iTokens) /= ' ') Then
        Call Message(                                                        &
               'FATAL ERROR in processing input: the coordinate system "' // &
               Trim(Tokens%Tokens(i_Name))                                // &
               '" is a standard coordinate system and must be '           // &
               'defined by name only',                                       &
               3                                                             &
             )
      End If
    End Do
    ZCoord = ZCoord_m_asl()
  Else If (Tokens%Tokens(i_Name) .CIEq. 'Pa') Then
    Do iTokens = 2, 6
      If (Tokens%Tokens(iTokens) /= ' ') Then
        Call Message(                                                        &
               'FATAL ERROR in processing input: the coordinate system "' // &
               Trim(Tokens%Tokens(i_Name))                                // &
               '" is a standard coordinate system and must be '           // &
               'defined by name only',                                       &
               3                                                             &
             )
      End If
    End Do
    ZCoord = ZCoord_Pa()
  Else If (Tokens%Tokens(i_Name) .CIEq. 'FL') Then
    Do iTokens = 2, 6
      If (Tokens%Tokens(iTokens) /= ' ') Then
        Call Message(                                                        &
               'FATAL ERROR in processing input: the coordinate system "' // &
               Trim(Tokens%Tokens(i_Name))                                // &
               '" is a standard coordinate system and must be '           // &
               'defined by name only',                                       &
               3                                                             &
             )
      End If
    End Do
    ZCoord = ZCoord_FL()

  ! Non-standard coord systems.
  Else

    CoordType = Token2Int(                                   &
                  C         = Tokens%Tokens(i_Type),         &
                  BlockKey  = 'Vertical Coordinate Systems', &
                  Item      = Tokens%Tokens(i_Name),         &
                  ColumnKey = 'Type'                         &
                )
    Unit      = Token2Std(                                   &
                  C         = Tokens%Tokens(i_Unit),         &
                  BlockKey  = 'Vertical Coordinate Systems', &
                  Item      = Tokens%Tokens(i_Name),         &
                  ColumnKey = 'Unit'                         &
                )

    ! $$ Further checks might be useful here - e.g. warning message if the user tries to
    ! $$ specify a pressure-eta defn (Token 4) and height-eta defn (Tokens 5+6) together.
    If (                                             &
      Tokens%Tokens(i_EtaDefinition)   == ' '  .and. &
      Tokens%Tokens(i_ModelTopHeight)  == ' '  .and. &
      Tokens%Tokens(i_InterfaceHeight) == ' '        &
    ) Then
      ZCoord = InitZCoord(Tokens%Tokens(i_Name), CoordType, Unit)
    Else If (Tokens%Tokens(i_EtaDefinition) /= ' ') Then
      i = FindEtaDefnIndex(Tokens%Tokens(i_EtaDefinition), EtaDefns)
      ZCoord = InitZCoord(                          &
                 Name      = Tokens%Tokens(i_Name), &
                 CoordType = CoordType,             &
                 Unit      = Unit,                  &
                 EtaDefn   = EtaDefns%EtaDefns(i)   &
               )
    Else
      ModelTopHeight  = Token2Std(                                      &
                          C         = Tokens%Tokens(i_ModelTopHeight),  &
                          BlockKey  = 'Vertical Coordinate Systems',    &
                          Item      = Tokens%Tokens(i_Name),            &
                          ColumnKey = 'Model Top Height'                &
                        )
      InterfaceHeight = Token2Std(                                      &
                          C         = Tokens%Tokens(i_InterfaceHeight), &
                          BlockKey  = 'Vertical Coordinate Systems',    &
                          Item      = Tokens%Tokens(i_Name),            &
                          ColumnKey = 'Interface Height'                &
                        )
      ZCoord = InitZCoord(                                &
                 Name            = Tokens%Tokens(i_Name), &
                 CoordType       = CoordType,             &
                 Unit            = Unit,                  &
                 ModelTopHeight  = ModelTopHeight,        &
                 InterfaceHeight = InterfaceHeight        &
               )
    End If

  End If

  Call AddZCoord(ZCoord, Coords)

End Subroutine Tokens2ZCoord

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine HGridInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                             &
         BlockKey             = 'Horizontal Grids:',      &
         NamedBlock           = .false.,                  &
         MatchMultipleLines   = .false.,                  &
         nHeaderLines         = 1,                        &
         ColumnKeys           = Reshape(                  &
                                  (/                      &
                                    'Name              ', & ! 1
                                    'H-Coord           ', & ! 2
                                    'nX                ', & ! 3
                                    'nY                ', & ! 4
                                    'dX                ', & ! 5
                                    'dY                ', & ! 6
                                    'X Min             ', & ! 7
                                    'Y Min             ', & ! 8
                                    'X-Array           ', & ! 9
                                    'Y-Array           ', & ! 10
                                    'Wrap?             ', & ! 11
                                    'Set Of Locations  ', & ! 12
                                    'Location Of Centre', & ! 13
                                    'X Centre          ', & ! 14
                                    'Y Centre          ', & ! 15
                                    'X Max             ', & ! 16
                                    'Y Max             ', & ! 17
                                    'X Range           ', & ! 18
                                    'Y Range           '  & ! 19
                                  /),                     &
                                  (/ 19, 1 /)             &
                                ),                        &
         Defaults             = (/                        &
                                  '          ', & ! 1
                                  '          ', & ! 2
                                  '          ', & ! 3
                                  '          ', & ! 4
                                  '          ', & ! 5
                                  '          ', & ! 6
                                  '          ', & ! 7
                                  '          ', & ! 8
                                  '          ', & ! 9
                                  '          ', & ! 10
                                  '          ', & ! 11
                                  '          ', & ! 12
                                  '          ', & ! 13
                                  '          ', & ! 14
                                  '          ', & ! 15
                                  '          ', & ! 16
                                  '          ', & ! 17
                                  '          ', & ! 18
                                  '          '  & ! 19
                                /),                       &
         TwoD                 = .false.,                  &
         UnrecognisedMessages = .true.,                   &
         DisjointColSpecs     = .true.,                   &
         HFBlockForms         = HFBlockForms              &
       )

End Subroutine HGridInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2HGrid(Tokens, Arrays, Locationses, Grids)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_),      Intent(In)    :: Tokens
  Type(Arrays_),      Intent(In)    :: Arrays
  Type(Locationses_), Intent(In)    :: Locationses ! $$ remove from here in favour of setup routine
  Type(Grids_),       Intent(InOut) :: Grids
  ! Local parameters:
  Integer, Parameter :: i_Name             = 1
  Integer, Parameter :: i_HCoord           = 2
  Integer, Parameter :: i_nX               = 3
  Integer, Parameter :: i_nY               = 4
  Integer, Parameter :: i_dX               = 5
  Integer, Parameter :: i_dY               = 6
  Integer, Parameter :: i_XMin             = 7
  Integer, Parameter :: i_YMin             = 8
  Integer, Parameter :: i_XArray           = 9
  Integer, Parameter :: i_YArray           = 10
  Integer, Parameter :: i_Wrap             = 11
  Integer, Parameter :: i_SetOfLocations   = 12
  Integer, Parameter :: i_LocationOfCentre = 13
  Integer, Parameter :: i_XCentre          = 14
  Integer, Parameter :: i_YCentre          = 15
  Integer, Parameter :: i_XMax             = 16
  Integer, Parameter :: i_YMax             = 17
  Integer, Parameter :: i_XRange           = 18
  Integer, Parameter :: i_YRange           = 19
  ! Locals:
  Integer   :: nX    !} Temporary variables, used to convert Tokens into the
  Integer   :: nY    !} right form.
  Real(Std) :: dX    !}
  Real(Std) :: dY    !}
  Real(Std) :: X0    !}
  Real(Std) :: Y0    !}
  Real(Std) :: XCen, XMax, XRange, YCen, YMax, YRange
  Integer   :: i
  Integer   :: j
  Real(Std) :: X(MaxArrayLength)
  Real(Std) :: Y(MaxArrayLength)
  Logical   :: Wrap
  Logical   :: Error

  ! Unstructured grids.
  If (                                             &
    Tokens%Tokens(i_HCoord          ) == ' ' .and. &
    Tokens%Tokens(i_nX              ) == ' ' .and. &
    Tokens%Tokens(i_nY              ) == ' ' .and. &
    Tokens%Tokens(i_dX              ) /= ' ' .and. &
    Tokens%Tokens(i_dY              ) /= ' ' .and. &
    Tokens%Tokens(i_XMin            ) == ' ' .and. &
    Tokens%Tokens(i_YMin            ) == ' ' .and. &
    Tokens%Tokens(i_XArray          ) == ' ' .and. &
    Tokens%Tokens(i_YArray          ) == ' ' .and. &
    Tokens%Tokens(i_Wrap            ) == ' ' .and. &
    Tokens%Tokens(i_SetOfLocations  ) /= ' ' .and. &
    Tokens%Tokens(i_LocationOfCentre) == ' '       &
  ) Then

    dX = Token2Std(                           &
           C         = Tokens%Tokens(i_dX),   &
           BlockKey  = 'Horizontal Grids',    &
           Item      = Tokens%Tokens(i_Name), &
           ColumnKey = 'dX'                   &
         )
    dY = Token2Std(                           &
           C         = Tokens%Tokens(i_dY),   &
           BlockKey  = 'Horizontal Grids',    &
           Item      = Tokens%Tokens(i_Name), &
           ColumnKey = 'dY'                   &
         )
    i  = FindLocationsIndex(Tokens%Tokens(i_SetOfLocations), Locationses)

    Call AddHGrid(                                  &
           InitHGrid(                               &
             Name      = Tokens%Tokens(i_Name),     &
             dX        = dX,                        &
             dY        = dY,                        &
             Locations = Locationses%Locationses(i) &
           ),                                       &
           Grids                                    &
         )

  ! Structured variable grids.
  Else If (                                        &
    Tokens%Tokens(i_HCoord          ) /= ' ' .and. &
    Tokens%Tokens(i_nX              ) == ' ' .and. &
    Tokens%Tokens(i_nY              ) == ' ' .and. &
    Tokens%Tokens(i_dX              ) == ' ' .and. &
    Tokens%Tokens(i_dY              ) == ' ' .and. &
    Tokens%Tokens(i_XMin            ) == ' ' .and. &
    Tokens%Tokens(i_YMin            ) == ' ' .and. &
    Tokens%Tokens(i_XArray          ) /= ' ' .and. &
    Tokens%Tokens(i_YArray          ) /= ' ' .and. &
    Tokens%Tokens(i_SetOfLocations  ) == ' ' .and. &
    Tokens%Tokens(i_LocationOfCentre) == ' '       &
  ) Then

    If (Tokens%Tokens(i_Wrap) == ' ') Then
      Wrap = .false.
    Else
      Wrap = Token2Log(                           &
               C         = Tokens%Tokens(i_Wrap), &
               BlockKey  = 'Horizontal Grids',    &
               Item      = Tokens%Tokens(i_Name), &
               ColumnKey = 'Wrap?'                &
             )
    End If
    i  = FindArrayIndex(Tokens%Tokens(i_XArray), Arrays)
    nX = Arrays%Arrays(i)%n
    Do j = 1, nX
      X(j) = Token2Std(                               &
               C         = Arrays%Arrays(i)%Array(j), &
               BlockKey  = 'Array',                   &
               Item      = Arrays%Arrays(i)%Name,     &
               ColumnKey = 'Array Values',            &
               i         = j                          &
             )
    End Do
    i  = FindArrayIndex(Tokens%Tokens(i_YArray), Arrays)
    nY = Arrays%Arrays(i)%n
    Do j = 1, nY
      Y(j) = Token2Std(                               &
               C         = Arrays%Arrays(i)%Array(j), &
               BlockKey  = 'Array',                   &
               Item      = Arrays%Arrays(i)%Name,     &
               ColumnKey = 'Array Values',            &
               i         = j                          &
             )
    End Do

    Call AddHGrid(                                 &
           InitHGrid(                              &
             Name       = Tokens%Tokens(i_Name),   &
             HCoordName = Tokens%Tokens(i_HCoord), &
             Wrap       = Wrap,                    &
             X          = X(1:nX),                 &
             Y          = Y(1:nY)                  &
           ),                                      &
           Grids                                   &
         )

  ! Structured non-variable grids, specified by X0, Y0.
  Else If (                                        &
    Tokens%Tokens(i_HCoord          ) /= ' ' .and. &
    Tokens%Tokens(i_nX              ) /= ' ' .and. &
    Tokens%Tokens(i_nY              ) /= ' ' .and. &
    Tokens%Tokens(i_XArray          ) == ' ' .and. &
    Tokens%Tokens(i_YArray          ) == ' ' .and. &
    Tokens%Tokens(i_SetOfLocations  ) == ' ' .and. &
    Tokens%Tokens(i_LocationOfCentre) == ' '       &
  ) Then

    If (Tokens%Tokens(i_Wrap) == ' ') Then
      Wrap = .false.
    Else
      Wrap = Token2Log(                           &
               C         = Tokens%Tokens(i_Wrap), &
               BlockKey  = 'Horizontal Grids',    &
               Item      = Tokens%Tokens(i_Name), &
               ColumnKey = 'Wrap?'                &
             )
    End If
    nX = Token2Int(                           &
           C         = Tokens%Tokens(i_nX),   &
           BlockKey  = 'Horizontal Grids',    &
           Item      = Tokens%Tokens(i_Name), &
           ColumnKey = 'nX'                   &
         )
    nY = Token2Int(                           &
           C         = Tokens%Tokens(i_nY),   &
           BlockKey  = 'Horizontal Grids',    &
           Item      = Tokens%Tokens(i_Name), &
           ColumnKey = 'nY'                   &
         )
    If (Tokens%Tokens(i_XMin   ) /= ' ') X0     = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_XMin   ), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'X Min'                   &
                                                  )
    If (Tokens%Tokens(i_XMax   ) /= ' ') XMax   = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_XMax   ), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'X Max'                   &
                                                  )
    If (Tokens%Tokens(i_XCentre) /= ' ') XCen   = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_XCentre), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'X Centre'                &
                                                  )
    If (Tokens%Tokens(i_XRange ) /= ' ') XRange = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_XRange ), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'X Range'                 &
                                                  )
    If (Tokens%Tokens(i_dX     ) /= ' ') dX     = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_dX     ), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'dX'                      &
                                                  )
    If (Tokens%Tokens(i_YMin   ) /= ' ') Y0     = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_YMin   ), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'Y Min'                   &
                                                  )
    If (Tokens%Tokens(i_YMax   ) /= ' ') YMax   = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_YMax   ), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'Y Max'                   &
                                                  )
    If (Tokens%Tokens(i_YCentre) /= ' ') YCen   = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_YCentre), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'Y Centre'                &
                                                  )
    If (Tokens%Tokens(i_YRange ) /= ' ') YRange = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_YRange ), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'Y Range'                 &
                                                  )
    If (Tokens%Tokens(i_dY     ) /= ' ') dY     = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_dY     ), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'dY'                      &
                                                  )

    Call MinMaxCentre(                                                                                                           &
           nX,                                                                                                                   &
           Tokens%Tokens(i_XMin), Tokens%Tokens(i_XMax), Tokens%Tokens(i_XCentre), Tokens%Tokens(i_XRange), Tokens%Tokens(i_dX), &
           X0, XMax, XCen, XRange, dX,                                                                                           &
           Error                                                                                                                 &
         )

    If (Error) Then
      Call Message(                                                      &
             'FATAL ERROR in reading the horizontal grid "'           // &
             Trim(Tokens%Tokens(i_Name))                              // &
             '" from the input file(s): '                             // &
             'exactly two of dX, X Range, X Min, X Centre and X Max ' // &
             'must be specified, but not dX and X Range',                &
             3                                                           &
           )
    End If

    Call MinMaxCentre(                                                                                                           &
           nY,                                                                                                                   &
           Tokens%Tokens(i_YMin), Tokens%Tokens(i_YMax), Tokens%Tokens(i_YCentre), Tokens%Tokens(i_YRange), Tokens%Tokens(i_dY), &
           Y0, YMax, YCen, YRange, dY,                                                                                           &
           Error                                                                                                                 &
         )

    If (Error) Then
      Call Message(                                                      &
             'FATAL ERROR in reading the horizontal grid "'           // &
             Trim(Tokens%Tokens(i_Name))                              // &
             '" from the input file(s): '                             // &
             'exactly two of dY, Y Range, Y Min, Y Centre and Y Max ' // &
             'must be specified, but not dY and Y Range',                &
             3                                                           &
           )
    End If

    Call AddHGrid(                                 &
           InitHGrid(                              &
             Name       = Tokens%Tokens(i_Name),   &
             HCoordName = Tokens%Tokens(i_HCoord), &
             Wrap       = Wrap,                    &
             nX         = nX,                      &
             nY         = nY,                      &
             dX         = dX,                      &
             dY         = dY,                      &
             X0         = X0,                      &
             Y0         = Y0                       &
           ),                                      &
           Grids                                   &
         )

  ! Structured non-variable grids, specified by LocationsName, CentreName.
  Else If (                                        &
    Tokens%Tokens(i_HCoord          ) /= ' ' .and. &
    Tokens%Tokens(i_nX              ) /= ' ' .and. &
    Tokens%Tokens(i_nY              ) /= ' ' .and. &
    Tokens%Tokens(i_XMin            ) == ' ' .and. &
    Tokens%Tokens(i_YMin            ) == ' ' .and. &
    Tokens%Tokens(i_XArray          ) == ' ' .and. &
    Tokens%Tokens(i_YArray          ) == ' ' .and. &
    Tokens%Tokens(i_SetOfLocations  ) /= ' ' .and. &
    Tokens%Tokens(i_LocationOfCentre) /= ' ' .and. &
    Tokens%Tokens(i_XCentre         ) == ' ' .and. &
    Tokens%Tokens(i_YCentre         ) == ' ' .and. &
    Tokens%Tokens(i_XMax            ) == ' ' .and. &
    Tokens%Tokens(i_YMax            ) == ' '       &
  ) Then

    If (Tokens%Tokens(i_Wrap) == ' ') Then
      Wrap = .false.
    Else
      Wrap = Token2Log(                           &
               C         = Tokens%Tokens(i_Wrap), &
               BlockKey  = 'Horizontal Grids',    &
               Item      = Tokens%Tokens(i_Name), &
               ColumnKey = 'Wrap?'                &
             )
    End If
    nX = Token2Int(                           &
           C         = Tokens%Tokens(i_nX),   &
           BlockKey  = 'Horizontal Grids',    &
           Item      = Tokens%Tokens(i_Name), &
           ColumnKey = 'nX'                   &
         )
    nY = Token2Int(                           &
           C         = Tokens%Tokens(i_nY),   &
           BlockKey  = 'Horizontal Grids',    &
           Item      = Tokens%Tokens(i_Name), &
           ColumnKey = 'nY'                   &
         )
    If (Tokens%Tokens(i_XRange ) /= ' ') XRange = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_XRange ), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'X Range'                 &
                                                  )
    If (Tokens%Tokens(i_dX     ) /= ' ') dX     = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_dX     ), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'dX'                      &
                                                  )
    XCen = 0.0
    If (Tokens%Tokens(i_YRange ) /= ' ') YRange = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_YRange ), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'Y Range'                 &
                                                  )
    If (Tokens%Tokens(i_dY     ) /= ' ') dY     = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_dY     ), &
                                                    BlockKey  = 'Horizontal Grids',       &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'dY'                      &
                                                  )
    YCen = 0.0

    Call MinMaxCentre(                                                  &
           nX,                                                          &
           ' ', ' ', '0', Tokens%Tokens(i_XRange), Tokens%Tokens(i_dX), &
           X0, XMax, XCen, XRange, dX,                                  &
           Error                                                        &
         )

    If (Error) Then
      Call Message(                                                     &
             'FATAL ERROR in reading the horizontal grid "'          // &
             Trim(Tokens%Tokens(i_Name))                             // &
             '" from the input file(s): '                            // &
             'one and only one of dX and X Range must be specified',    &
             3                                                          &
           )
    End If

    Call MinMaxCentre(                                                  &
           nY,                                                          &
           ' ', ' ', '0', Tokens%Tokens(i_YRange), Tokens%Tokens(i_dY), &
           Y0, YMax, YCen, YRange, dY,                                  &
           Error                                                        &
         )

    If (Error) Then
      Call Message(                                                     &
             'FATAL ERROR in reading the horizontal grid "'          // &
             Trim(Tokens%Tokens(i_Name))                             // &
             '" from the input file(s): '                            // &
             'one and only one of dY and Y Range must be specified',    &
             3                                                          &
           )
    End If

    Call AddHGrid(                                             &
           InitHGrid(                                          &
             Name          = Tokens%Tokens(i_Name),            &
             HCoordName    = Tokens%Tokens(i_HCoord),          &
             Wrap          = Wrap,                             &
             nX            = nX,                               &
             nY            = nY,                               &
             dX            = dX,                               &
             dY            = dY,                               &
             LocationsName = Tokens%Tokens(i_SetOfLocations),  &
             CentreName    = Tokens%Tokens(i_LocationOfCentre) &
           ),                                                  &
           Grids                                               &
         )

  Else

    Call Message(                                                          &
           'FATAL ERROR in reading the horizontal grid "'               // &
           Trim(Tokens%Tokens(i_Name))                                  // &
           '" from the input file(s): '                                 // &
           'an incorrect combination of grid variables has been given',    &
           3                                                               &
         )

  End If

End Subroutine Tokens2HGrid

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine MinMaxCentre(                      &
             n,                               &
             CMin, CMax, CCentre, CRange, Cd, &
              Min,  Max,  Centre,  Range,  d, &
             Error                            &
           )
! .

! The reals corresponding to the non-blank character strings should be computed before calling this routine.
! They are not computed in this routine because the calling routine will be able to give a better error 
! message when the character string is inappropriate.

  Implicit None
  ! Argument list:
  Integer,      Intent(In)    :: n
  Character(*), Intent(In)    :: CMin
  Character(*), Intent(In)    :: CMax
  Character(*), Intent(In)    :: CCentre
  Character(*), Intent(In)    :: CRange
  Character(*), Intent(In)    :: Cd
  Real(Std),    Intent(InOut) :: Min
  Real(Std),    Intent(InOut) :: Max
  Real(Std),    Intent(InOut) :: Centre ! $$ these aren't always computed - should they be?
  Real(Std),    Intent(InOut) :: Range  !
  Real(Std),    Intent(InOut) :: d
  Logical,      Intent(Out)   :: Error
  ! Locals:
  Integer :: i

  i = 0
  If (CMin    /= ' ') i = i + 1
  If (CMax    /= ' ') i = i + 1
  If (CCentre /= ' ') i = i + 1
  If (CRange  /= ' ') i = i + 1
  If (Cd      /= ' ') i = i + 1
  Error = i /= 2 .or. (CRange /= ' ' .and. Cd /= ' ')

  ! If (Cd == ' ' .and. n == 1) Then ! $$ this combination will result in d = ? / 0

  If (Cd /= ' ' .and. CMin /= ' ') Then

    ! d   used as is.
    ! Min used as is.
    Max = Min + Real(n - 1, Std) * d

  Else If (Cd /= ' ' .and. CCentre /= ' ') Then

    ! d      used as is.
    ! Centre used as is.
    Min    = Centre - Real(n - 1, Std) * d / 2.0
    Max    = Min + Real(n - 1, Std) * d

  Else If (Cd /= ' ' .and. CMax /= ' ') Then

    ! d   used as is.
    ! Max used as is.
    Min = Max - Real(n - 1, Std) * d

  Else If (CMin /= ' ' .and. CCentre /= ' ') Then

    ! Min    used as is.
    ! Centre used as is.
    d      = (Centre - Min) * 2.0 / Real(n - 1, Std)
    Max    = Min + Real(n - 1, Std) * d

  Else If (CMin /= ' ' .and. CMax /= ' ') Then

    ! Min used as is.
    ! Max used as is.
    d   = (Max - Min) / Real(n - 1, Std)

  Else If (CMin /= ' ' .and. CRange /= ' ') Then

    ! Min   used as is.
    ! Range used as is.
    d     = Range / Real(n - 1, Std)
    Max   = Min + Real(n - 1, Std) * d

  Else If (CCentre /= ' ' .and. CMax /= ' ') Then

    ! Centre used as is.
    ! Max    used as is.
    d      = (Max - Centre) * 2.0 / Real(n - 1, Std)
    Min    = Max - Real(n - 1, Std) * d

  Else If (CCentre /= ' ' .and. CRange /= ' ') Then

    ! Centre used as is.
    ! Range  used as is.
    d      = Range / Real(n - 1, Std)
    Min    = Centre - Real(n - 1, Std) * d / 2.0
    Max    = Min + Real(n - 1, Std) * d

  Else If (CMax /= ' ' .and. CRange /= ' ') Then

    ! Max   used as is.
    ! Range used as is.
    d     = Range / Real(n - 1, Std)
    Min   = Max - Real(n - 1, Std) * d

  End If

End Subroutine MinMaxCentre

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine ZGridInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                      &
         BlockKey             = 'Vertical Grids:', &
         NamedBlock           = .false.,           &
         MatchMultipleLines   = .false.,           &
         nHeaderLines         = 1,                 &
         ColumnKeys           = Reshape(           &
                                  (/               &
                                    'Name       ', & ! 1
                                    'Z-Coord    ', & ! 2
                                    'nZ         ', & ! 3
                                    'dZ         ', & ! 4
                                    'Z0         ', & ! 5
                                    'Z-Array    ', & ! 6
                                    'Av Z-Array ', & ! 7
                                    'Index-Array'  & ! 8
                                  /),              &
                                  (/ 8, 1 /)       &
                                ),                 &
         Defaults             = (/                 &
                                  '   ',       & ! 1
                                  '   ',       & ! 2
                                  '   ',       & ! 3
                                  '   ',       & ! 4
                                  '   ',       & ! 5
                                  '   ',       & ! 6
                                  '   ',       & ! 7
                                  '   '        & ! 8
                                /),                &
         TwoD                 = .false.,           &
         UnrecognisedMessages = .true.,            &
         DisjointColSpecs     = .true.,            &
         HFBlockForms         = HFBlockForms       &
       )

End Subroutine ZGridInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2ZGrid(Tokens, Arrays, Coords, Grids)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_), Intent(In)    :: Tokens
  Type(Arrays_), Intent(In)    :: Arrays
  Type(Coords_), Intent(In)    :: Coords
  Type(Grids_),  Intent(InOut) :: Grids
  ! Local parameters:
  Integer, Parameter :: i_Name       = 1
  Integer, Parameter :: i_ZCoord     = 2
  Integer, Parameter :: i_nZ         = 3
  Integer, Parameter :: i_dZ         = 4
  Integer, Parameter :: i_Z0         = 5
  Integer, Parameter :: i_ZArray     = 6
  Integer, Parameter :: i_AvZArray   = 7
  Integer, Parameter :: i_IndexArray = 8
  ! Locals:
  Integer      :: nZ    !} Temporary variables, used to convert Tokens into the
  Real(Std)    :: dZ    !} right form.
  Real(Std)    :: Z0    !}
  Type(ZGrid_) :: ZGrid !
  Integer      :: i
  Integer      :: j
  Real(Std)    :: Z(MaxArrayLength)
  Real(Std)    :: AvZ(MaxArrayLength)
  Integer      :: iZ(MaxArrayLength)

  If (                                     &
    Tokens%Tokens(i_nZ      ) /= ' ' .and. &
    Tokens%Tokens(i_dZ      ) /= ' ' .and. &
    Tokens%Tokens(i_Z0      ) /= ' ' .and. &
    Tokens%Tokens(i_ZArray  ) == ' ' .and. &
    Tokens%Tokens(i_AvZArray) == ' '       &
  ) Then

    nZ = Token2Int(                           &
           C         = Tokens%Tokens(i_nZ),   &
           BlockKey  = 'Vertical Grids',      &
           Item      = Tokens%Tokens(i_Name), &
           ColumnKey = 'nZ'                   &
         )
    dZ = Token2Std(                           &
           C         = Tokens%Tokens(i_dZ),   &
           BlockKey  = 'Vertical Grids',      &
           Item      = Tokens%Tokens(i_Name), &
           ColumnKey = 'dZ'                   &
         )
    Z0 = Token2Std(                           &
           C         = Tokens%Tokens(i_Z0),   &
           BlockKey  = 'Vertical Grids',      &
           Item      = Tokens%Tokens(i_Name), &
           ColumnKey = 'Z0'                   &
         )

    If (Tokens%Tokens(i_IndexArray) /= ' ') Then
      i = FindArrayIndex(Tokens%Tokens(i_IndexArray), Arrays)
      If (Arrays%Arrays(i)%n /= nZ) Then
        Call Message(                                                       &
               'Error in Tokens2ZGrid: Index-Array should have nZ entries', &
               3                                                            &
             )
      End If
      Do j = 1, nZ
        iZ(j) = Token2Int(                                 &
                  C         = Arrays%Arrays(i)%Array(j),   &
                  BlockKey  = 'Array',                     &
                  Item      = Tokens%Tokens(i_IndexArray), & ! $$ use array name or name from tokens consistently?
                  ColumnKey = 'Array Values',              &
                  i         = j                            &
                )
      End Do
      ZGrid = InitZGrid(                              &
                Name       = Tokens%Tokens(i_Name),   &
                ZCoordName = Tokens%Tokens(i_ZCoord), &
                nZ         = nZ,                      &
                dZ         = dZ,                      &
                Z0         = Z0,                      &
                iZ         = iZ(1:nZ)                 &
              )
    Else
      ZGrid = InitZGrid(                              &
                Name       = Tokens%Tokens(i_Name),   &
                ZCoordName = Tokens%Tokens(i_ZCoord), &
                nZ         = nZ,                      &
                dZ         = dZ,                      &
                Z0         = Z0                       &
              )
    End If

  Else If (                              &
    Tokens%Tokens(i_nZ    ) == ' ' .and. &
    Tokens%Tokens(i_dZ    ) == ' ' .and. &
    Tokens%Tokens(i_Z0    ) == ' ' .and. &
    Tokens%Tokens(i_ZArray) /= ' '       &
  ) Then

    If (Tokens%Tokens(i_ZArray) .CIEq. 'Use Eta Levels') Then
      i = FindZCoordIndex(Tokens%Tokens(i_ZCoord), Coords)
      If (Coords%ZCoords(i)%CoordType == 5) Then
        nZ      = Coords%ZCoords(i)%nEtaLevels
        Z(1:nZ) = Coords%ZCoords(i)%Eta(1:nZ)
      Else
        Call Message(                                                                &
               'FATAL ERROR in Tokens2ZGrid: Use Eta Levels specified but coord ' // &
               'system is not an eta coord system',                                  &
               3                                                                     &
             )
      End If
    Else
      i  = FindArrayIndex(Tokens%Tokens(i_ZArray), Arrays)
      nZ = Arrays%Arrays(i)%n
      Do j = 1, nZ
        Z(j) = Token2Std(                               &
                 C         = Arrays%Arrays(i)%Array(j), &
                 BlockKey  = 'Array',                   &
                 Item      = Arrays%Arrays(i)%Name,     &
                 ColumnKey = 'Array Values',            &
                 i         = j                          &
               )
      End Do
    End If

    AvZ(:) = 0.0

    If (Tokens%Tokens(i_AvZArray) .CIEq. 'Z on Boundaries') Then

      If (Tokens%Tokens(i_IndexArray) /= ' ') Then
        i = FindArrayIndex(Tokens%Tokens(i_IndexArray), Arrays)
        If (Arrays%Arrays(i)%n /= nZ) Then
          Call Message(                                                                                  &
                 'Error in Tokens2ZGrid: Index-Array should have the same number of entries as Z-Array', &
                 3                                                                                       &
               )
        End If
        Do j = 1, nZ
          iZ(j) = Token2Int(                                 &
                    C         = Arrays%Arrays(i)%Array(j),   &
                    BlockKey  = 'Array',                     &
                    Item      = Tokens%Tokens(i_IndexArray), &
                    ColumnKey = 'Array Values',              &
                    i         = j                            &
                )
        End Do
        ZGrid = InitZGrid(                                 &
                  Name          = Tokens%Tokens(i_Name),   &
                  ZCoordName    = Tokens%Tokens(i_ZCoord), &
                  nZ            = nZ,                      &
                  Z             = Z(1:nZ),                 &
                  AvZ           = AvZ(1:nZ),               & ! $$ needed?
                  ZOnBoundaries = .true.,                  &
                  iZ            = iZ(1:nZ)                 &
                )
      Else
        ZGrid = InitZGrid(                                 &
                  Name          = Tokens%Tokens(i_Name),   &
                  ZCoordName    = Tokens%Tokens(i_ZCoord), &
                  nZ            = nZ,                      &
                  Z             = Z(1:nZ),                 &
                  AvZ           = AvZ(1:nZ),               & ! $$ needed?
                  ZOnBoundaries = .true.                   &
                )
      End If

    Else ! $$ note if av z-array blank, AvZ is still passed to InitZGrid but is zero. Is this OK?

      If (Tokens%Tokens(i_AvZArray) /= ' ') Then
        i = FindArrayIndex(Tokens%Tokens(i_AvZArray), Arrays)
        If (Arrays%Arrays(i)%n /= nZ) Then
          Call Message(                                                                                 &
                 'Error in Tokens2ZGrid: Av Z-Array should have the same number of entries as Z-Array', &
                 3                                                                                      &
               )
        End If
        Do j = 1, nZ
          AvZ(j) = Token2Std(                               &
                     C         = Arrays%Arrays(i)%Array(j), &
                     BlockKey  = 'Array',                   &
                     Item      = Arrays%Arrays(i)%Name,     &
                     ColumnKey = 'Array Values',            &
                     i         = j                          &
                   )
        End Do
      End If

      If (Tokens%Tokens(i_IndexArray) /= ' ') Then
        i = FindArrayIndex(Tokens%Tokens(i_IndexArray), Arrays)
        If (Arrays%Arrays(i)%n /= nZ) Then
          Call Message(                                                                             &
                 'Error in Tokens2ZGrid: Index-Array should have the same number of entries as ' // &
                 'there are grid points',                                                           &
                 3                                                                                  &
               )
        End If
        Do j = 1, nZ
          iZ(j) = Token2Int(                                 &
                    C         = Arrays%Arrays(i)%Array(j),   &
                    BlockKey  = 'Array',                     &
                    Item      = Tokens%Tokens(i_IndexArray), &
                    ColumnKey = 'Array Values',              &
                    i         = j                            &
                )
        End Do
        ZGrid = InitZGrid(                                 &
                  Name          = Tokens%Tokens(i_Name),   &
                  ZCoordName    = Tokens%Tokens(i_ZCoord), &
                  nZ            = nZ,                      &
                  Z             = Z(1:nZ),                 &
                  AvZ           = AvZ(1:nZ),               &
                  ZOnBoundaries = .false.,                 &
                  iZ            = iZ(1:nZ)                 &
                )
      Else
        ZGrid = InitZGrid(                                 &
                  Name          = Tokens%Tokens(i_Name),   &
                  ZCoordName    = Tokens%Tokens(i_ZCoord), &
                  nZ            = nZ,                      &
                  Z             = Z(1:nZ),                 &
                  AvZ           = AvZ(1:nZ),               &
                  ZOnBoundaries = .false.                  &
                )
      End If

    End If

  Else

    Call Message(                                                                             &
           'Error in Tokens2ZGrid: either nZ, dZ and Z0 or '                               // &
           'Names of Z and Av Z Arrays (but not both) must be given in defining a Z Grid',    &
           3                                                                                  &
         )

  End If

  Call AddZGrid(ZGrid, Grids)

End Subroutine Tokens2ZGrid

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine TGridInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                      &
         BlockKey             = 'Temporal Grids:', &
         NamedBlock           = .false.,           &
         MatchMultipleLines   = .false.,           &
         nHeaderLines         = 1,                 &
         ColumnKeys           = Reshape(           &
                                  (/               &
                                    'Name   ',     & ! 1
                                    'nT     ',     & ! 2
                                    'dT     ',     & ! 3
                                    'T0     ',     & ! 4
                                    'T-Array'      & ! 5
                                  /),              &
                                  (/ 5, 1 /)       &
                                ),                 &
         Defaults             = (/                 &
                                  '   ',       & ! 1
                                  '   ',       & ! 2
                                  '   ',       & ! 3
                                  '   ',       & ! 4
                                  '   '        & ! 5
                                /),                &
         TwoD                 = .false.,           &
         UnrecognisedMessages = .true.,            &
         DisjointColSpecs     = .true.,            &
         HFBlockForms         = HFBlockForms       &
       )

End Subroutine TGridInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2TGrid(Tokens, Arrays, Grids)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_), Intent(In)    :: Tokens
  Type(Arrays_), Intent(In)    :: Arrays
  Type(Grids_),  Intent(InOut) :: Grids
  ! Local parameters:
  Integer, Parameter :: i_Name   = 1
  Integer, Parameter :: i_nT     = 2
  Integer, Parameter :: i_dT     = 3
  Integer, Parameter :: i_T0     = 4
  Integer, Parameter :: i_TArray = 5
  ! Locals:
  Integer      :: nT       !} Temporary variables, used to convert Tokens into the
  Type(Time_)  :: dT       !} right form.
  Type(Time_)  :: T0       !}
  Type(TGrid_) :: TGrid    !
  Integer      :: i
  Integer      :: j
  Type(Time_)  :: T(MaxArrayLength)

  If (                                   &
    Tokens%Tokens(i_nT    ) /= ' ' .and. &
    Tokens%Tokens(i_dT    ) /= ' ' .and. &
    Tokens%Tokens(i_T0    ) /= ' ' .and. &
    Tokens%Tokens(i_TArray) == ' '       &
  ) Then

    nT = Token2Int(                           &
           C         = Tokens%Tokens(i_nT),   &
           BlockKey  = 'Temporal Grids',      &
           Item      = Tokens%Tokens(i_Name), &
           ColumnKey = 'nT'                   &
         )
    dT = Token2Time(                          &
           C         = Tokens%Tokens(i_dT),   &
           BlockKey  = 'Temporal Grids',      &
           Item      = Tokens%Tokens(i_Name), &
           ColumnKey = 'dT',                  &
           Interval  = .true.                 &
         )
    T0 = Token2Time(                          &
           C         = Tokens%Tokens(i_T0),   &
           BlockKey  = 'Temporal Grids',      &
           Item      = Tokens%Tokens(i_Name), &
           ColumnKey = 'T0'                   &
         )

    TGrid = InitTGrid(                      &
              Name = Tokens%Tokens(i_Name), &
              nT   = nT,                    &
              dT   = dT,                    &
              T0   = T0                     &
            )

  Else If (                              &
    Tokens%Tokens(i_nT    ) == ' ' .and. &
    Tokens%Tokens(i_dT    ) == ' ' .and. &
    Tokens%Tokens(i_T0    ) == ' ' .and. &
    Tokens%Tokens(i_TArray) /= ' '       &
  ) Then

    i  = FindArrayIndex(Tokens%Tokens(i_TArray), Arrays)
    nT = Arrays%Arrays(i)%n
    Do j = 1, nT
      T(j) = Token2Time(                              &
               C         = Arrays%Arrays(i)%Array(j), &
               BlockKey  = 'Array',                   &
               Item      = Tokens%Tokens(i_TArray),   &
               ColumnKey = 'Array Values',            &
               i         = j                          &
             )
    End Do

    TGrid = InitTGrid(                      &
              Name = Tokens%Tokens(i_Name), &
              nT   = nT,                    &
              T    = T                      &
            )

  Else

    Call Message(                                                                  &
           'Error in Tokens2TGrid: either nT, dT and T0 or '                    // &
           'Name of T Array (but not both) must be given in defining a T Grid',    &
           3                                                                       &
         )

  End If

  Call AddTGrid(TGrid, Grids)

End Subroutine Tokens2TGrid

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine DomainInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                             &
         BlockKey             = 'Domains:',               &
         NamedBlock           = .false.,                  &
         MatchMultipleLines   = .false.,                  &
         nHeaderLines         = 1,                        &
         ColumnKeys           = Reshape(                  &
                                  (/                      &
                                    'Name              ', & ! 1
                                    'H-Coord           ', & ! 2
                                    'X Min             ', & ! 3
                                    'X Max             ', & ! 4
                                    'Y Min             ', & ! 5
                                    'Y Max             ', & ! 6
                                    'H Unbounded?      ', & ! 7
                                    'Z-Coord           ', & ! 8
                                    'Z Max             ', & ! 9 $$ rename Z Top
                                    'Z Unbounded?      ', & ! 10
                                    'Start Time        ', & ! 11
                                    'End Time          ', & ! 12
                                    'Max Travel Time   ', & ! 13
                                    'T Unbounded?      ', & ! 14
                                    'X Centre          ', & ! 15
                                    'Y Centre          ', & ! 16
                                    'X Range           ', & ! 17
                                    'Y Range           ', & ! 18
                                    'Duration          ', & ! 19
                                    'Set of locations  ', & ! 20
                                    'Location of centre', & ! 21
                                    'X Unbounded?      ', & ! 22
                                    'Y Unbounded?      '  & ! 23
                                  /),                     &
                                  (/ 23, 1 /)             &
                                ),                        &
         Defaults             = (/                        &
                                  '   ', & ! 1
                                  '   ', & ! 2
                                  '   ', & ! 3
                                  '   ', & ! 4
                                  '   ', & ! 5
                                  '   ', & ! 6
                                  '   ', & ! 7
                                  '   ', & ! 8
                                  '   ', & ! 9 $$ rename Z Top
                                  '   ', & ! 10
                                  '   ', & ! 11
                                  '   ', & ! 12
                                  '   ', & ! 13
                                  '   ', & ! 14
                                  '   ', & ! 15
                                  '   ', & ! 16
                                  '   ', & ! 17
                                  '   ', & ! 18
                                  '   ', & ! 19
                                  '   ', & ! 20
                                  '   ', & ! 21
                                  '   ', & ! 22
                                  '   '  & ! 23
                                /),                       &
         TwoD                 = .false.,                  &
         UnrecognisedMessages = .true.,                   &
         DisjointColSpecs     = .true.,                   &
         HFBlockForms         = HFBlockForms              &
       )

End Subroutine DomainInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2Domain(Tokens, Domains)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_),  Intent(In)    :: Tokens
  Type(Domains_), Intent(InOut) :: Domains
  ! Local parameters:
  Integer, Parameter :: i_Name             =  1
  Integer, Parameter :: i_HCoord           =  2
  Integer, Parameter :: i_XMin             =  3
  Integer, Parameter :: i_XMax             =  4
  Integer, Parameter :: i_YMin             =  5
  Integer, Parameter :: i_YMax             =  6
  Integer, Parameter :: i_HUnbounded       =  7
  Integer, Parameter :: i_ZCoord           =  8
  Integer, Parameter :: i_ZMax             =  9
  Integer, Parameter :: i_ZUnbounded       =  10
  Integer, Parameter :: i_StartTime        =  11
  Integer, Parameter :: i_EndTime          =  12
  Integer, Parameter :: i_MaxTravelTime    =  13
  Integer, Parameter :: i_TUnbounded       =  14
  Integer, Parameter :: i_XCentre          =  15
  Integer, Parameter :: i_YCentre          =  16
  Integer, Parameter :: i_XRange           =  17
  Integer, Parameter :: i_YRange           =  18
  Integer, Parameter :: i_Duration         =  19
  Integer, Parameter :: i_SetOfLocations   =  20
  Integer, Parameter :: i_LocationOfCentre =  21
  Integer, Parameter :: i_XUnbounded       =  22
  Integer, Parameter :: i_YUnbounded       =  23
  ! Locals:
  Real(Std)     :: XMin
  Real(Std)     :: XMax
  Real(Std)     :: YMin
  Real(Std)     :: YMax
  Logical       :: HUnbounded
  Logical       :: XUnbounded
  Logical       :: YUnbounded
  Real(Std)     :: ZMax
  Logical       :: ZUnbounded
  Type(Time_)   :: StartTime
  Type(Time_)   :: EndTime
  Type(Time_)   :: Duration
  Logical       :: TUnbounded
  Type(Time_)   :: MaxTravelTime
  Type(Domain_) :: Domain
  Real(Std)     :: XCen, YCen, XRange, YRange, dX, dY
  Logical       :: Error

  ! Either H Unbounded? should be given or X Unbounded? and Y Unbounded?
  If (                                       &
    Tokens%Tokens(i_HUnbounded) /= ' ' .and. &
    Tokens%Tokens(i_XUnbounded) == ' ' .and. &
    Tokens%Tokens(i_YUnbounded) == ' '       &
  ) Then

    HUnbounded = Token2Log(                                 &
                   C         = Tokens%Tokens(i_HUnbounded), &
                   BlockKey  = 'Domains',                   &
                   Item      = Tokens%Tokens(i_Name),       &
                   ColumnKey = 'H Unbounded?'               &
                 )

    XUnbounded = HUnbounded
    YUnbounded = HUnbounded

  Else If (                                  &
    Tokens%Tokens(i_HUnbounded) == ' ' .and. &
    Tokens%Tokens(i_XUnbounded) /= ' ' .and. &
    Tokens%Tokens(i_YUnbounded) /= ' '       &
  ) Then

    XUnbounded = Token2Log(                                 &
                   C         = Tokens%Tokens(i_XUnbounded), &
                   BlockKey  = 'Domains',                   &
                   Item      = Tokens%Tokens(i_Name),       &
                   ColumnKey = 'X Unbounded?'               &
                 )

    YUnbounded = Token2Log(                                 &
                   C         = Tokens%Tokens(i_YUnbounded), &
                   BlockKey  = 'Domains',                   &
                   Item      = Tokens%Tokens(i_Name),       &
                   ColumnKey = 'Y Unbounded?'               &
                 )

  Else

    Call Message(                                                                 &
           'Error in Tokens2Domain: either X Unbounded? and Y Unbounded?, or ' // &
           'H Unbounded? (but not both) must be given in defining a Domain',      &
           3                                                                      &
         )

  End If


  If (XUnbounded) Then

    XMin = 0.0
    XMax = 0.0
    If (                                            &
      Tokens%Tokens(i_XMin)             /= ' ' .or. &
      Tokens%Tokens(i_XMax)             /= ' ' .or. &
      Tokens%Tokens(i_XCentre)          /= ' ' .or. &
      Tokens%Tokens(i_XRange)           /= ' '      &
    ) Then
      Call Message(                                              &
             'FATAL ERROR in processing input for domain "'   // &
             Trim(Tokens%Tokens(i_Name))                      // &
             '": X Min, X Max, etc. should not be given for ' // &
             'domains unbounded in the X direction', 3           &
           )
    End If

  Else If (Tokens%Tokens(i_SetOfLocations) == ' ' .and. Tokens%Tokens(i_LocationOfCentre) == ' ') Then

    If (Tokens%Tokens(i_XMin   ) /= ' ') XMin   = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_XMin   ), &
                                                    BlockKey  = 'Domains',                &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'X Min'                   &
                                                  )
    If (Tokens%Tokens(i_XMax   ) /= ' ') XMax   = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_XMax   ), &
                                                    BlockKey  = 'Domains',                &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'X Max'                   &
                                                  )
    If (Tokens%Tokens(i_XCentre) /= ' ') XCen   = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_XCentre), &
                                                    BlockKey  = 'Domains',                &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'X Centre'                &
                                                  )
    If (Tokens%Tokens(i_XRange ) /= ' ') XRange = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_XRange ), &
                                                    BlockKey  = 'Domains',                &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'X Range'                 &
                                                  )
    
    Call MinMaxCentre(                                                                                    &
           2,                                                                                             &
           Tokens%Tokens(i_XMin), Tokens%Tokens(i_XMax), Tokens%Tokens(i_XCentre),                        &
           Tokens%Tokens(i_XRange), ' ',                                                                  &
           XMin, XMax, XCen, XRange, dX,                                                                  &
           Error                                                                                          &
         )

    If (Error) Then
      Call Message(                                                                    &
             'FATAL ERROR in reading the domain "'                                  // &
             Trim(Tokens%Tokens(i_Name))                                            // &
             '" from the input file(s): '                                           // &
             'exactly two of X Range, X Min, X Centre and X Max must be specified',    &
             3                                                                         &
           )
    End If


  Else If (                                        &
    Tokens%Tokens(i_XMin)             == ' ' .and. &
    Tokens%Tokens(i_XMax)             == ' ' .and. &
    Tokens%Tokens(i_XCentre)          == ' ' .and. &
    Tokens%Tokens(i_XRange)           /= ' ' .and. &
    Tokens%Tokens(i_SetOfLocations)   /= ' ' .and. &
    Tokens%Tokens(i_LocationOfCentre) /= ' '       &
  ) Then

    XMax = Token2Std(                             &
             C         = Tokens%Tokens(i_XRange), &
             BlockKey  = 'Domains',               &
             Item      = Tokens%Tokens(i_Name),   &
             ColumnKey = 'X Range'                &
           )
    XMax = XMax/2.0
    XMin = -XMax ! $$ Note XMin/Max set to +/- XRange / 2.0 etc for the moment.

  Else

    Call Message(                                                       &
           'FATAL ERROR in processing input for domain "'            // &
           Trim(Tokens%Tokens(i_Name))                               // &
           '": an inappropriate combination of variables is given '  // &
           'for the X direction of the domain.',                        &
           3                                                            &
         )

  End If

  If (YUnbounded) Then

    YMin = 0.0
    YMax = 0.0
    If (                                            &
      Tokens%Tokens(i_YMin)             /= ' ' .or. &
      Tokens%Tokens(i_YMax)             /= ' ' .or. &
      Tokens%Tokens(i_YCentre)          /= ' ' .or. &
      Tokens%Tokens(i_YRange)           /= ' '      &
    ) Then
      Call Message(                                              &
             'FATAL ERROR in processing input for domain "'   // &
             Trim(Tokens%Tokens(i_Name))                      // &
             '": Y Min, Y Max, etc. should not be given for ' // &
             'domains unbounded in the Y direction', 3           &
           )
    End If

  Else If (Tokens%Tokens(i_SetOfLocations) == ' ' .and. Tokens%Tokens(i_LocationOfCentre) == ' ') Then

    If (Tokens%Tokens(i_YMin   ) /= ' ') YMin   = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_YMin   ), &
                                                    BlockKey  = 'Domains',                &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'Y Min'                   &
                                                  )
    If (Tokens%Tokens(i_YMax   ) /= ' ') YMax   = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_YMax   ), &
                                                    BlockKey  = 'Domains',                &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'Y Max'                   &
                                                  )
    If (Tokens%Tokens(i_YCentre) /= ' ') YCen   = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_YCentre), &
                                                    BlockKey  = 'Domains',                &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'Y Centre'                &
                                                  )
    If (Tokens%Tokens(i_YRange ) /= ' ') YRange = Token2Std(                              &
                                                    C         = Tokens%Tokens(i_YRange ), &
                                                    BlockKey  = 'Domains',                &
                                                    Item      = Tokens%Tokens(i_Name),    &
                                                    ColumnKey = 'Y Range'                 &
                                                  )
    
    Call MinMaxCentre(                                                                                     &
           2,                                                                                              &
           Tokens%Tokens(i_YMin), Tokens%Tokens(i_YMax), Tokens%Tokens(i_YCentre),                         &
           Tokens%Tokens(i_YRange), ' ',                                                                   &
           YMin, YMax, YCen, YRange, dY,                                                                   &
           Error                                                                                           &
         )

    If (Error) Then
      Call Message(                                                                    &
             'FATAL ERROR in reading the domain "'                                  // &
             Trim(Tokens%Tokens(i_Name))                                            // &
             '" from the input file(s): '                                           // &
             'exactly two of Y Range, Y Min, Y Centre and Y Max must be specified',    &
             3                                                                         &
           )
    End If


  Else If (                                        &
    Tokens%Tokens(i_YMin)             == ' ' .and. &
    Tokens%Tokens(i_YMax)             == ' ' .and. &
    Tokens%Tokens(i_YCentre)          == ' ' .and. &
    Tokens%Tokens(i_YRange)           /= ' ' .and. &
    Tokens%Tokens(i_SetOfLocations)   /= ' ' .and. &
    Tokens%Tokens(i_LocationOfCentre) /= ' '       &
  ) Then

    YMax = Token2Std(                             &
             C         = Tokens%Tokens(i_YRange), &
             BlockKey  = 'Domains',               &
             Item      = Tokens%Tokens(i_Name),   &
             ColumnKey = 'Y Range'                &
           )
    YMax = YMax/2.0
    YMin = -YMax ! $$ Note YMin/Max set to +/- YRange / 2.0 etc for the moment.

  Else

    Call Message(                                                       &
           'FATAL ERROR in processing input for domain "'            // &
           Trim(Tokens%Tokens(i_Name))                               // &
           '": an inappropriate combination of variables is given '  // &
           'for the Y direction of the domain.',                        &
           3                                                            &
         )

  End If

  If (XUnbounded .and. YUnbounded) Then

    If (                                            &
      Tokens%Tokens(i_HCoord)           /= ' ' .or. &
      Tokens%Tokens(i_SetOfLocations)   /= ' ' .or. &
      Tokens%Tokens(i_LocationOfCentre) /= ' '      &
    ) Then
      Call Message(                                                     &
             'FATAL ERROR in processing input for domain "'          // &
             Trim(Tokens%Tokens(i_Name))                             // &
             '": H-Coord, Set of Locations, and Location of Centre ' // &
             'should not be given for unbounded domains', 3             &
           )
    End If

  End If

  ZUnbounded = Token2Log(                                 &
                 C         = Tokens%Tokens(i_ZUnbounded), &
                 BlockKey  = 'Domains',                   &
                 Item      = Tokens%Tokens(i_Name),       &
                 ColumnKey = 'Z Unbounded?'               &
               )
  If (ZUnbounded) Then
    ZMax = 0.0
    If (Tokens%Tokens(i_ZCoord) /= ' ' .or. Tokens%Tokens(i_ZMax) /= ' ') Then
      Call Message(                                              &
             'FATAL ERROR in processing input for domain "'   // &
             Trim(Tokens%Tokens(i_Name))                           // &
             '": Z Max, Z-Coord etc should not be given for ' // &
             'unbounded domains', 3                              &
           )
    End If
  Else
    ZMax = Token2Std(                           &
             C         = Tokens%Tokens(i_ZMax), &
             BlockKey  = 'Domains',             &
             Item      = Tokens%Tokens(i_Name), &
             ColumnKey = 'Z Max'                &
           )
  End If

  TUnbounded = Token2Log(                                 &
                 C         = Tokens%Tokens(i_TUnbounded), &
                 BlockKey  = 'Domains',                   &
                 Item      = Tokens%Tokens(i_Name),       &
                 ColumnKey = 'T Unbounded?'               &
               )
  If (TUnbounded) Then
    StartTime = InfPastTime  ()
    EndTime   = InfFutureTime()
    If (Tokens%Tokens(i_StartTime) /= ' ' .or. Tokens%Tokens(i_EndTime) /= ' ' &
                                          .or. Tokens%Tokens(i_Duration) /= ' ') Then
      Call Message(                                                                            &
             'FATAL ERROR: Domain start and end time and duration should not be given for ' // &
             'unbounded domains', 3                                                            &
           )
    End If
  Else
    If (Tokens%Tokens(i_StartTime) /= ' ' .and. Tokens%Tokens(i_EndTime) /= ' ' &
                                          .and. Tokens%Tokens(i_Duration) /= ' ') Then
      Call Message(                                                                   &
             'FATAL ERROR in reading the domain "'                                 // &
             Trim(Tokens%Tokens(i_Name))                                           // &
             '" from the input file(s): '                                          // &
             'exactly two of Start Time, End Time and Duration must be specified',    &
             3                                                                        &
           )
    Else If (Tokens%Tokens(i_StartTime) /= ' ' .and. Tokens%Tokens(i_EndTime) /= ' ') Then
      StartTime = Token2Time(                               &
                    C         = Tokens%Tokens(i_StartTime), &
                    BlockKey  = 'Domains',                  &
                    Item      = Tokens%Tokens(i_Name),      &
                    ColumnKey = 'Start Time'                &
                  )
      EndTime   = Token2Time(                             &
                    C         = Tokens%Tokens(i_EndTime), &
                    BlockKey  = 'Domains',                &
                    Item      = Tokens%Tokens(i_Name)  ,  &
                    ColumnKey = 'End Time'                &
                  )
    Else If (Tokens%Tokens(i_StartTime) /= ' ' .and. Tokens%Tokens(i_Duration) /= ' ') Then
      StartTime = Token2Time(                               &
                    C         = Tokens%Tokens(i_StartTime), &
                    BlockKey  = 'Domains',                  &
                    Item      = Tokens%Tokens(i_Name),      &
                    ColumnKey = 'Start Time'                &
                  )
      Duration  = Token2Time(                              &
                    C         = Tokens%Tokens(i_Duration), &
                    BlockKey  = 'Domains',                 &
                    Item      = Tokens%Tokens(i_Name),     &
                    ColumnKey = 'Duration',                &
                    Interval  = .true.                     &
                  )
      EndTime   = StartTime + Duration
    Else If (Tokens%Tokens(i_EndTime) /= ' ' .and. Tokens%Tokens(i_Duration) /= ' ') Then
      EndTime   = Token2Time(                             &
                    C         = Tokens%Tokens(i_EndTime), &
                    BlockKey  = 'Domains',                &
                    Item      = Tokens%Tokens(i_Name),    &
                    ColumnKey = 'End Time'                &
                  )
      Duration  = Token2Time(                              &
                    C         = Tokens%Tokens(i_Duration), &
                    BlockKey  = 'Domains',                 &
                    Item      = Tokens%Tokens(i_Name),     &
                    ColumnKey = 'Duration',                &
                    Interval  = .true.                     &
                  )
      StartTime = EndTime - Duration
    Else
      Call Message(                                                                   &
             'FATAL ERROR in reading the domain "'                                 // &
             Trim(Tokens%Tokens(i_Name))                                           // &
             '" from the input file(s): '                                          // &
             'exactly two of Start Time, End Time and Duration must be specified',    &
             3                                                                        &
           )
    End If
  End If

  MaxTravelTime = Token2Time(                                   &
                    C         = Tokens%Tokens(i_MaxTravelTime), &
                    BlockKey  = 'Domains',                      &
                    Item      = Tokens%Tokens(i_Name),          &
                    ColumnKey = 'Max Travel Time',              &
                    Interval  = .true.                          &
                  )

  Domain = InitDomain(                                          &
             Name          = Tokens%Tokens(i_Name),             &
             HCoordName    = Tokens%Tokens(i_HCoord),           &
             XMin          = XMin,                              &
             XMax          = XMax,                              &
             YMin          = YMin,                              &
             YMax          = YMax,                              &
             XUnbounded    = XUnbounded,                        &
             YUnbounded    = YUnbounded,                        &
             ZCoordName    = Tokens%Tokens(i_ZCoord),           &
             ZMax          = ZMax,                              &
             ZUnbounded    = ZUnbounded,                        &
             StartTime     = StartTime,                         &
             EndTime       = EndTime,                           &
             MaxTravelTime = MaxTravelTime,                     &
             LocationsName = Tokens%Tokens(i_SetOfLocations),   &
             CentreName    = Tokens%Tokens(i_LocationOfCentre)  &
           )

  Call AddDomain(Domain, Domains)

End Subroutine Tokens2Domain

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine PrototypeMetInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                      &
         BlockKey             = 'Prototype Met Module Instances:', &
         NamedBlock           = .false.,                           &
         MatchMultipleLines   = .false.,                           &
         nHeaderLines         = 1,                                 &
         ColumnKeys           = Reshape(                           &
                                  (/                               &
                                    'Name             ',           & ! 1
                                    'H-Coord          ',           & ! 2
                                    'Z-Coord          ',           & ! 3
                                    'Met File         ',           & ! 4
                                    'Update On Demand?'            & ! 5
                                  /),                              &
                                  (/ 5, 1 /)                       &
                                ),                                 &
         Defaults             = (/                                 &
                                  '  ',                      & ! 1
                                  '  ',                      & ! 2
                                  '  ',                      & ! 3
                                  '  ',                      & ! 4
                                  'No'                             & ! 5
                                /),                                &
         TwoD                 = .false.,                           &
         UnrecognisedMessages = .true.,                            &
         DisjointColSpecs     = .true.,                            &
         HFBlockForms         = HFBlockForms                       &
       )

End Subroutine PrototypeMetInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2PrototypeMet(Tokens, MainOpts, Mets)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)    :: Tokens
  Type(MainOpts_), Intent(In)    :: MainOpts
  Type(Mets_),     Intent(InOut) :: Mets
  ! Local parameters:
  Integer, Parameter :: i_Name           = 1
  Integer, Parameter :: i_HCoord         = 2
  Integer, Parameter :: i_ZCoord         = 3
  Integer, Parameter :: i_MetFile        = 4
  Integer, Parameter :: i_UpdateOnDemand = 5
  ! Locals:
  Logical :: UpdateOnDemand ! Indicates the met module instance is to be updated using update-on-demand.

  UpdateOnDemand = Token2Log(                                      &
                     C         = Tokens%Tokens(i_UpdateOnDemand),  &
                     BlockKey  = 'Prototype Met Module Instances', &
                     Item      = Tokens%Tokens(i_Name),            &
                     ColumnKey = 'Update On Demand?'               &
                   )

  Call AddPrototypeMet(                              &
         InitPrototypeMet(                           &
           MetName        = Tokens%Tokens(i_Name),   &
           HCoordName     = Tokens%Tokens(i_HCoord), &
           ZCoordName     = Tokens%Tokens(i_ZCoord), &
           FixedMet       = MainOpts%FixedMet,       &
           UpdateOnDemand = UpdateOnDemand,          &
           MetFile        = Tokens%Tokens(i_MetFile) &
         ),                                          &
         Mets                                        &
       )

End Subroutine Tokens2PrototypeMet

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine SingleSiteMetInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                        &
         BlockKey             = 'Single Site Met Module Instances:', &
         NamedBlock           = .false.,                             &
         MatchMultipleLines   = .false.,                             &
         nHeaderLines         = 1,                                   &
         ColumnKeys           = Reshape(                             &
                                  (/                                 &
                                    'Name                  ',        & ! 1
                                    'H-Coord               ',        & ! 2
                                    'Long                  ',        & ! 3
                                    'Lat                   ',        & ! 4
                                    'Height                ',        & ! 5
                                    'z0                    ',        & ! 6
                                    'z0D                   ',        & ! 7
                                    'Representative?       ',        & ! 8
                                    'Met File              ',        & ! 9
                                    'Ignore Fixed Met Time?',        & ! 10
                                    'Mesoscale SigU        ',        & ! 11
                                    'Mesoscale TauU        ',        & ! 12
                                    'Free Trop SigU        ',        & ! 13
                                    'Free Trop SigW        ',        & ! 14
                                    'Free Trop TauU        ',        & ! 15
                                    'Free Trop TauW        ',        & ! 16
                                    'Update On Demand?     '         & ! 17
                                  /),                                &
                                  (/ 17, 1 /)                        &
                                ),                                   &
         Defaults             = (/                                   &
                                  '       ',          & ! 1
                                  '       ',          & ! 2
                                  '       ',          & ! 3
                                  '       ',          & ! 4
                                  '       ',          & ! 5
                                  '       ',          & ! 6
                                  '       ',          & ! 7
                                  '       ',          & ! 8
                                  '       ',          & ! 9
                                  '       ',          & ! 10
                                  '0.5    ',                         & ! 11
                                  '7200.0 ',                         & ! 12
                                  '0.25   ',                         & ! 13
                                  '0.1    ',                         & ! 14
                                  '300.0  ',                         & ! 15
                                  '100.0  ',                         & ! 16
                                  'No     '                          & ! 17
                                /),                                  &
         TwoD                 = .false.,                             &
         UnrecognisedMessages = .true.,                              &
         DisjointColSpecs     = .true.,                              &
         HFBlockForms         = HFBlockForms                         &
       )

End Subroutine SingleSiteMetInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2SingleSiteMet(Tokens, MainOpts, Mets)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)    :: Tokens
  Type(MainOpts_), Intent(In)    :: MainOpts
  Type(Mets_),     Intent(InOut) :: Mets
  ! Local parameters:
  Integer, Parameter :: i_Name               = 1
  Integer, Parameter :: i_HCoord             = 2
  Integer, Parameter :: i_Long               = 3
  Integer, Parameter :: i_Lat                = 4
  Integer, Parameter :: i_Height             = 5
  Integer, Parameter :: i_z0                 = 6
  Integer, Parameter :: i_z0D                = 7
  Integer, Parameter :: i_Representative     = 8
  Integer, Parameter :: i_MetFile            = 9
  Integer, Parameter :: i_IgnoreFixedMetTime = 10
  Integer, Parameter :: i_MeanderSigU        = 11
  Integer, Parameter :: i_MeanderTauU        = 12
  Integer, Parameter :: i_FreeTropSigU       = 13
  Integer, Parameter :: i_FreeTropSigW       = 14
  Integer, Parameter :: i_FreeTropTauU       = 15
  Integer, Parameter :: i_FreeTropTauW       = 16
  Integer, Parameter :: i_UpdateOnDemand     = 17
  ! Locals:
  Real(Std) :: Long
  Real(Std) :: Lat
  Real(Std) :: Z
  Real(Std) :: Z0
  Real(Std) :: Z0D
  Real(Std) :: SigUUM
  Real(Std) :: TauUUM
  Real(Std) :: SigU2HPlus
  Real(Std) :: SigW2HPlus
  Real(Std) :: TauUHPlus
  Real(Std) :: TauWHPlus
  Logical   :: Representative
  Logical   :: IgnoreFixedMetTime
  Logical   :: UpdateOnDemand     ! Indicates the met module instance is to be updated using update-on-demand.

  Long       = Token2Std(                                        &
                 C         = Tokens%Tokens(i_Long),              &
                 BlockKey  = 'Single Site Met Module Instances', &
                 Item      = Tokens%Tokens(i_Name),              &
                 ColumnKey = 'Long'                              &
               )
  Lat        = Token2Std(                                        &
                 C         = Tokens%Tokens(i_Lat),               &
                 BlockKey  = 'Single Site Met Module Instances', &
                 Item      = Tokens%Tokens(i_Name),              &
                 ColumnKey = 'Lat'                               &
               )
  Z          = Token2Std(                                        &
                 C         = Tokens%Tokens(i_Height),            &
                 BlockKey  = 'Single Site Met Module Instances', &
                 Item      = Tokens%Tokens(i_Name),              &
                 ColumnKey = 'Height'                            &
               )
  Z0         = Token2Std(                                        &
                 C         = Tokens%Tokens(i_z0),                &
                 BlockKey  = 'Single Site Met Module Instances', &
                 Item      = Tokens%Tokens(i_Name),              &
                 ColumnKey = 'z0'                                &
               )
  Z0D        = Token2Std(                                        &
                 C         = Tokens%Tokens(i_z0D),               &
                 BlockKey  = 'Single Site Met Module Instances', &
                 Item      = Tokens%Tokens(i_Name),              &
                 ColumnKey = 'z0D'                               &
               )
  SigUUM     = (                                                   &
                 Token2Std(                                        &
                   C         = Tokens%Tokens(i_MeanderSigU),       &
                   BlockKey  = 'Single Site Met Module Instances', &
                   Item      = Tokens%Tokens(i_Name),              &
                   ColumnKey = 'Meander SigU'                      &
                 )                                                 &
               )**2
  TauUUM     = Token2Std(                                        &
                 C         = Tokens%Tokens(i_MeanderTauU),       &
                 BlockKey  = 'Single Site Met Module Instances', &
                 Item      = Tokens%Tokens(i_Name),              &
                 ColumnKey = 'Meander TauU'                      &
               )
  SigU2HPlus = (                                                   &
                 Token2Std(                                        &
                   C         = Tokens%Tokens(i_FreeTropSigU),      &
                   BlockKey  = 'Single Site Met Module Instances', &
                   Item      = Tokens%Tokens(i_Name),              &
                   ColumnKey = 'Free Trop SigU'                    &
                 )                                                 &
               )**2
  SigW2HPlus = (                                                   &
                 Token2Std(                                        &
                   C         = Tokens%Tokens(i_FreeTropSigW),      &
                   BlockKey  = 'Single Site Met Module Instances', &
                   Item      = Tokens%Tokens(i_Name),              &
                   ColumnKey = 'Free Trop SigW'                    &
                 )                                                 &
               )**2
  TauUHPlus  = Token2Std(                                        &
                 C         = Tokens%Tokens(i_FreeTropTauU),      &
                 BlockKey  = 'Single Site Met Module Instances', &
                 Item      = Tokens%Tokens(i_Name),              &
                 ColumnKey = 'Free Trop TauU'                    &
               )
  TauWHPlus  = Token2Std(                                        &
                 C         = Tokens%Tokens(i_FreeTropTauW),      &
                 BlockKey  = 'Single Site Met Module Instances', &
                 Item      = Tokens%Tokens(i_Name),              &
                 ColumnKey = 'Free Trop TauW'                    &
               )
  Representative = Token2Log(                                        &
                     C         = Tokens%Tokens(i_Representative),    &
                     BlockKey  = 'Single Site Met Module Instances', &
                     Item      = Tokens%Tokens(i_Name),              &
                     ColumnKey = 'Representative?'                   &
                   )
  IgnoreFixedMetTime = Token2Log(                                         &
                         C         = Tokens%Tokens(i_IgnoreFixedMetTime), &
                         BlockKey  = 'Single Site Met Module Instances',  &
                         Item      = Tokens%Tokens(i_Name),               &
                         ColumnKey = 'Ignore Fixed Met Time?'             &
                       )

  UpdateOnDemand = Token2Log(                                        &
                     C         = Tokens%Tokens(i_UpdateOnDemand),    &
                     BlockKey  = 'Single Site Met Module Instances', &
                     Item      = Tokens%Tokens(i_Name),              &
                     ColumnKey = 'Update On Demand?'                 &
                   )

  Call AddSingleSiteMet(                                               &
         InitSingleSiteMet(                                            &
           MetName            = Tokens%Tokens(i_Name),                 &
           HCoordName         = Tokens%Tokens(i_HCoord),               &
           FixedMet           = MainOpts%FixedMet,                     &
           UpdateOnDemand     = UpdateOnDemand,                        &
           Sequential         = .true.,                                & !$$
           IgnoreFixedMetTime = IgnoreFixedMetTime,                    &
           TSample            = Char2Time('01:00', Interval = .true.), & !$$
           Long               = Long,                                  &
           Lat                = Lat,                                   &
           Z                  = Z,                                     &
           Z0                 = Z0,                                    &
           Z0D                = Z0D,                                   &
           Alpha              = 1.0,                                   & !$$
           AlphaD             = 1.0,                                   & !$$
           R                  = 0.23,                                  & !$$
           RD                 = 0.23,                                  & !$$
           RecipLMOMax        = 1.0,                                   & !$$
           RecipLMOMaxD       = 1.0,                                   & !$$
           PrecipFactor       = 1.0,                                   & !$$
           Representative     = Representative,                        &
           SigUUM             = SigUUM,                                &
           TauUUM             = TauUUM,                                &
           SigU2HPlus         = SigU2HPlus,                            &
           SigW2HPlus         = SigW2HPlus,                            &
           TauUHPlus          = TauUHPlus,                             &
           TauWHPlus          = TauWHPlus,                             &
           MetFile            = Tokens%Tokens(i_MetFile),              &
           MetLimitsFile      = 'dummy'                                &
         ),                                                            &
         Mets                                                          &
       )

End Subroutine Tokens2SingleSiteMet

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine NWPMetDefnInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                        &
         BlockKey             = 'NWP Met Definitions:',              &
         NamedBlock           = .false.,                             &
         MatchMultipleLines   = .false.,                             &
         nHeaderLines         = 1,                                   &
         ColumnKeys           = Reshape(                             &
                                  (/                                 &
                                    'Name                         ', & ! 1
                                    'Binary Format                ', & ! 2
                                    'dT                           ', & ! 3
                                    'T0                           ', & ! 4
                                    'Day Per File                 ', & ! 5
                                    'Suffix                       ', & ! 6
                                    'Topography File              ', & ! 7
                                    'Next Heat Flux               ', & ! 8
                                    'Prefix                       ', & ! 9
                                    'File Type                    ', & ! 10
                                    'Met File Structure Definition', & ! 11
                                    'Z-Coord - W                  ', & ! 12
                                    'Z-Coord - Cloud Height       ', & ! 13
                                    'Z-Grid                       ', & ! 14
                                    'Z-Grid - UV                  ', & ! 15
                                    'Z-Grid - W                   ', & ! 16
                                    'Z-Grid - P                   ', & ! 17
                                    'H-Grid                       ', & ! 18
                                    'H-Grid - U                   ', & ! 19
                                    'H-Grid - V                   ', & ! 20
                                    'Next Precipitation           ', & ! 21
                                    'Next Cloud                   ', & ! 22
                                    'Mesoscale SigU               ', & ! 23
                                    'Mesoscale TauU               '  & ! 24
                                  /),                                &
                                  (/ 24, 1 /)                        &
                                ),                                   &
         Defaults             = (/                                   &
                                  '              ', & ! 1
                                  '              ', & ! 2
                                  '              ', & ! 3
                                  '1/1/2000 00:00',                  & ! 4
                                  '              ', & ! 5
                                  '              ', & ! 6
                                  '              ', & ! 7
                                  '              ', & ! 8
                                  '              ', & ! 9
                                  '              ', & ! 10
                                  '              ', & ! 11
                                  '              ', & ! 12
                                  '              ', & ! 13
                                  '              ', & ! 14
                                  '              ', & ! 15
                                  '              ', & ! 16
                                  '              ', & ! 17
                                  '              ', & ! 18
                                  '              ', & ! 19
                                  '              ', & ! 20
                                  'No            ',                  & ! 21
                                  'No            ',                  & ! 22
                                  '              ',                  & ! 23
                                  '              '                   & ! 24
                                /),                                  &
         TwoD                 = .false.,                             &
         UnrecognisedMessages = .true.,                              &
         DisjointColSpecs     = .true.,                              &
         HFBlockForms         = HFBlockForms                         &
       )

End Subroutine NWPMetDefnInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2NWPMetDefn(Tokens, MetDefns)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)    :: Tokens
  Type(MetDefns_), Intent(InOut) :: MetDefns
  ! Local parameters:
  Integer, Parameter :: i_Name                       = 1
  Integer, Parameter :: i_BinaryFormat               = 2
  Integer, Parameter :: i_dT                         = 3
  Integer, Parameter :: i_T0                         = 4
  Integer, Parameter :: i_DayPerFile                 = 5
  Integer, Parameter :: i_Suffix                     = 6
  Integer, Parameter :: i_TopographyFile             = 7
  Integer, Parameter :: i_NextHeatFlux               = 8
  Integer, Parameter :: i_Prefix                     = 9
  Integer, Parameter :: i_FileType                   = 10
  Integer, Parameter :: i_MetFileStructureDefinition = 11
  Integer, Parameter :: i_ZCoordW                    = 12
  Integer, Parameter :: i_ZCoordCloudHeight          = 13
  Integer, Parameter :: i_ZGrid                      = 14
  Integer, Parameter :: i_ZGridUV                    = 15
  Integer, Parameter :: i_ZGridW                     = 16
  Integer, Parameter :: i_ZGridP                     = 17
  Integer, Parameter :: i_HGrid                      = 18
  Integer, Parameter :: i_HGridU                     = 19
  Integer, Parameter :: i_HGridV                     = 20
  Integer, Parameter :: i_NextPrecipitation          = 21
  Integer, Parameter :: i_NextCloud                  = 22
  Integer, Parameter :: i_MesoscaleSigU              = 23
  Integer, Parameter :: i_MesoscaleTauU              = 24
  ! Locals:
  Type(Time_)       :: T0           ! Reference time for the first met field
  Type(Time_)       :: Dt           ! Time interval between met fields
  Logical           :: DayPerFile
  Logical           :: NextHeatFlux ! Use next time step for heat flux and ustar
  Logical           :: NextPrecip   ! Use next time step for precip
  Logical           :: NextCloud    ! Use next time step for cloud
  Type(NWPMetDefn_) :: NWPMetDefn
  Integer           :: nSuffixs
  Integer           :: nMetDefn2Names
  Real(Std)         :: SigUUM
  Real(Std)         :: TauUUM
  Character(MaxCharLength) :: Suffixs(MaxNWPMetDefn2s)
  Character(MaxCharLength) :: MetDefn2Names(MaxNWPMetDefn2s)

  T0           = Token2Time(                          &
                   C         = Tokens%Tokens(i_T0),   &
                   BlockKey  = 'NWP Met Definitions', &
                   Item      = Tokens%Tokens(i_Name), &
                   ColumnKey = 'T0'                   &
                 )

  Dt           = Token2Time(                          &
                   C         = Tokens%Tokens(i_dT),   &
                   BlockKey  = 'NWP Met Definitions', &
                   Item      = Tokens%Tokens(i_Name), &
                   ColumnKey = 'dT',                  &
                   Interval  = .true.                 &
                 )

  DayPerFile   = Token2Log(                                 &
                   C         = Tokens%Tokens(i_DayPerFile), &
                   BlockKey  = 'NWP Met Definitions',       &
                   Item      = Tokens%Tokens(i_Name),       &
                   ColumnKey = 'Day Per File'               &
                 )

  NextHeatFlux = Token2Log(                                   &
                   C         = Tokens%Tokens(i_NextHeatFlux), &
                   BlockKey  = 'NWP Met Definitions',         &
                   Item      = Tokens%Tokens(i_Name),         &
                   ColumnKey = 'Next Heat Flux'               &
                 )

  NextPrecip   = Token2Log(                                        &
                   C         = Tokens%Tokens(i_NextPrecipitation), &
                   BlockKey  = 'NWP Met Definitions',              &
                   Item      = Tokens%Tokens(i_Name),              &
                   ColumnKey = 'Next Precipitation'                &
                 )

  NextCloud    = Token2Log(                                &
                   C         = Tokens%Tokens(i_NextCloud), &
                   BlockKey  = 'NWP Met Definitions',      &
                   Item      = Tokens%Tokens(i_Name),      &
                   ColumnKey = 'Next Cloud'                &
                 )

  Call ParseListChar(                         &
         C         = Tokens%Tokens(i_Suffix), &
         Delim     = ';',                     &
         BlockKey  = 'NWP Met Definitions',   &
         Item      = Tokens%Tokens(i_Name),   &
         ColumnKey = 'Suffix',                &
         nValues   = nSuffixs,                &
         Values    = Suffixs                  &
       )

  Call ParseListChar(                                             &
         C         = Tokens%Tokens(i_MetFileStructureDefinition), &
         Delim     = ';',                                         &
         BlockKey  = 'NWP Met Definitions',                       &
         Item      = Tokens%Tokens(i_Name),                       &
         ColumnKey = 'Met File Structure Definition',             &
         nValues   = nMetDefn2Names,                              &
         Values    = MetDefn2Names                                &
       )

  SigUUM = (Token2Std(                                    &
              C         = Tokens%Tokens(i_MesoscaleSigU), &
              BlockKey  = 'NWP Met Definitions',          &
              Item      = Tokens%Tokens(i_Name),          &
              ColumnKey = 'Mesoscale SigU'                &
            ))**2
                        
  TauUUM = Token2Std(                                    &
             C         = Tokens%Tokens(i_MesoscaleTauU), &
             BlockKey  = 'NWP Met Definitions',          &
             Item      = Tokens%Tokens(i_Name),          &
             ColumnKey = 'Mesoscale TauU'                &
           )

  NWPMetDefn = InitNWPMetDefn(                                        &
                 Name           = Tokens%Tokens(i_Name),              &
                 BinaryFormat   = Tokens%Tokens(i_BinaryFormat),      &
                 FileType       = Tokens%Tokens(i_FileType),          &
                 Prefix         = Tokens%Tokens(i_Prefix),            &
                 Suffixs        = Suffixs(1:nSuffixs),                &
                 TopogFile      = Tokens%Tokens(i_TopographyFile),    &
                 DayPerFile     = DayPerFile,                         &
                 T0             = T0,                                 &
                 Dt             = Dt,                                 &
                 NextHeatFlux   = NextHeatFlux,                       &
                 NextPrecip     = NextPrecip,                         &
                 NextCloud      = NextCloud,                          &
                 MetDefn2Names  = MetDefn2Names(1:nMetDefn2Names),    &
                 ZCoordNameW    = Tokens%Tokens(i_ZCoordW),           &
                 ZCoordNameCl   = Tokens%Tokens(i_ZCoordCloudHeight), &
                 HGridName      = Tokens%Tokens(i_HGrid),             &
                 HGridNameU     = Tokens%Tokens(i_HGridU),            &
                 HGridNameV     = Tokens%Tokens(i_HGridV),            &
                 ZGridName      = Tokens%Tokens(i_ZGrid),             &
                 ZGridNameUV    = Tokens%Tokens(i_ZGridUV),           &
                 ZGridNameW     = Tokens%Tokens(i_ZGridW),            &
                 ZGridNameP     = Tokens%Tokens(i_ZGridP),            &
                 SigUUM         = SigUUM,                             &
                 TauUUM         = TauUUM                              &
               )

  Call AddNWPMetDefn(NWPMetDefn, MetDefns)

End Subroutine Tokens2NWPMetDefn

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine NWPMetDefn2InputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                         &
         BlockKey             = 'NWP Met File Structure Definition:', &
         NamedBlock           = .true.,                               &
         MatchMultipleLines   = .false.,                              &
         nHeaderLines         = 1,                                    &
         ColumnKeys           = Reshape(                              &
                                  (/                                  &
                                    'Field Name      ',               & ! 1
                                    'Lowest Level    ',               & ! 2
                                    'Highest Level   ',               & ! 3
                                    'Field Code      ',               & ! 4
                                    '3-d?            ',               & ! 5
                                    'Field Qualifiers',               & ! 6
                                    'NC Field Name'                   & ! 7
                                  /),                                 &
                                  (/ 7, 1 /)                          &
                                ),                                    &
         Defaults             = (/                                    &
                                  '                   ', & ! 1
                                  '                   ', & ! 2
                                  '                   ', & ! 3
                                  '                   ', & ! 4
                                  '                   ', & ! 5
                                  '                   ', & ! 6
                                  '                   '  & ! 7
                                /),                                   &
         TwoD                 = .true.,                               &
         UnrecognisedMessages = .true.,                               &
         DisjointColSpecs     = .true.,                               &
         HFBlockForms         = HFBlockForms                          &
       )

End Subroutine NWPMetDefn2InputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2NWPMetDefn2(Tokens, MetDefns)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)    :: Tokens
  Type(MetDefns_), Intent(InOut) :: MetDefns
  ! Local parameters:
  Integer, Parameter :: i_FieldName       = 1
  Integer, Parameter :: i_LowestLevel     = 2
  Integer, Parameter :: i_HighestLevel    = 3
  Integer, Parameter :: i_FieldCode       = 4
  Integer, Parameter :: i_3d              = 5
  Integer, Parameter :: i_FieldQualifiers = 6
  Integer, Parameter :: i_NCFieldName     = 7
  ! Locals:
  Integer                   :: i
  Integer                   :: iFieldQualifier
  Character(MaxTokenLength) :: FieldNames(MaxNWPMetFields)
  Integer                   :: LowestLevels(MaxNWPMetFields)
  Integer                   :: HighestLevels(MaxNWPMetFields)
  Logical                   :: Top(MaxNWPMetFields)
  Integer                   :: FieldCodes(MaxNWPMetFields)
  Logical                   :: ThreeD(MaxNWPMetFields)
  Logical                   :: Total(MaxNWPMetFields)
  Logical                   :: Dynamic(MaxNWPMetFields)
  Type(NWPMetDefn2_)        :: NWPMetDefn2
  Integer                   :: nFieldQualifiers
  Character(MaxCharLength)  :: FieldQualifiers(MaxFieldQualifiers)
  Character(MaxTokenLength) :: NCFieldNames(MaxNWPMetFields)
  ! i                :: Loop index.
  ! iFieldQualifier  :: Loop index
  ! FieldNames       :: Names of met fields.
  ! FieldCodes       :: Codes of met fields. 0 indicates no code given.
  ! LowestLevels     :: Lowest model levels of met fields.
  ! HighestLevels    :: Highest model levels of met fields.
  ! Top              :: Indicates that highest level of a met field is top of the grid.
  ! Total            :: Indicates that the total or dynamic cloud stashes are total
  ! Dynamic          :: Indicates that the total or dynamic cloud stashes are dynamic
  ! NWPMetDefn2      :: NWP met definition (part 2 - met file structure definition).
  ! nFieldQualifiers :: number of Field Qualifiers for field
  ! FieldQualifiers  :: Field Qualifiers for field
  ! NCFieldNames     :: Names of met fields in the NetCDF files.

  If (Tokens%nLines > MaxNWPMetFields) Then ! $$ replace by maxRows... Then no need for if-test
    Call Message('ERROR in Tokens2NWPMetDefn2: too many fields ' // &
                 'declared in the met definition block '         // &
                 Trim(Tokens%BlockName),                            &
                 3)
  End If

  Do i = 1, Tokens%nLines

    FieldNames(i)   = Tokens%Tokens2d(i_FieldName, i)

    LowestLevels(i) = Token2Int(                                         &
                        C         = Tokens%Tokens2d(i_LowestLevel, i),   &
                        BlockKey  = 'NWP Met File Structure Definition', &
                        Item      = Tokens%BlockName,                    &
                        ColumnKey = 'Lowest Level'                       &
                      )

    If (Tokens%Tokens2d(i_HighestLevel, i) .CIEq. 'Top') Then
      HighestLevels(i) = 0
      Top(i)           = .true.
    Else
      HighestLevels(i) = Token2Int(                                         &
                           C         = Tokens%Tokens2d(i_HighestLevel, i),  &
                           BlockKey  = 'NWP Met File Structure Definition', &
                           Item      = Tokens%BlockName,                    &
                           ColumnKey = 'Highest Level'                      &
                         )
      Top(i)           = .false.
    End If

    FieldCodes(i) = Token2Int(                                         &
                      C         = Tokens%Tokens2d(i_FieldCode, i),     &
                      BlockKey  = 'NWP Met File Structure Definition', &
                      Item      = Tokens%BlockName,                    &
                      ColumnKey = 'Field Code'                         &
                    )

    ThreeD(i) = Token2Log(                                         &
                  C         = Tokens%Tokens2d(i_3d, i),            &
                  BlockKey  = 'NWP Met File Structure Definition', &
                  Item      = Tokens%BlockName,                    &
                  ColumnKey = '3-d?'                               &
                )

    Call ParseListChar(                                       &
           C         = Tokens%Tokens2d(i_FieldQualifiers, i), &
           Delim     = ';',                                   &
           BlockKey  = 'NWP Met File Structure Definition',   &
           Item      = Tokens%BlockName,                      &
           ColumnKey = 'Field Qualifiers',                    &
           nValues   = nFieldQualifiers,                      &
           Values    = FieldQualifiers                        &
         )
         
    Total(i)   = .false.
    Dynamic(i) = .false.
    Do iFieldQualifier = 1, nFieldQualifiers
      If (FieldQualifiers(iFieldQualifier) .CIEq. 'Total') Then
        Total(i) = .true.
      Else If (FieldQualifiers(iFieldQualifier) .CIEq. 'Dynamic') Then
        Dynamic(i) = .true.
      Else
        Call Message('FATAL ERROR in Tokens2NWPMetDefn2: Field qualifier "'// &
                     Trim(FieldQualifiers(iFieldQualifier))                // &
                     '" not recognised in the met definition block '       // &
                     Trim(Tokens%BlockName),                                  &
                     3)
      End If
    End Do
    
    If (Total(i) .and. Dynamic(i)) Then
      Call Message('FATAL ERROR in Tokens2NWPMetDefn2: Fields '     // &
                   'declared in the met definition block '          // &
                   Trim(Tokens%BlockName)                           // &
                   ' cannot be both "Total" and "Dynamic"',            &
                   3)
    End If

    NCFieldNames(i) = Tokens%Tokens2d(i_NCFieldName, i)

  End Do

  NWPMetDefn2 = InitNWPMetDefn2(                                  &
                  Name          = Tokens%BlockName,               &
                  FieldNames    = FieldNames   (1:Tokens%nLines), &
                  LowestLevels  = LowestLevels (1:Tokens%nLines), &
                  HighestLevels = HighestLevels(1:Tokens%nLines), &
                  Top           = Top          (1:Tokens%nLines), &
                  FieldCodes    = FieldCodes   (1:Tokens%nLines), &
                  ThreeD        = ThreeD       (1:Tokens%nLines), &
                  Total         = Total        (1:Tokens%nLines), &
                  NCFieldNames  = NCFieldNames (1:Tokens%nLines)  &
                )

  Call AddNWPMetDefn2(NWPMetDefn2, MetDefns)

End Subroutine Tokens2NWPMetDefn2

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine NWPMetInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                  &
         BlockKey             = 'NWP Met Module Instances:',   &
         NamedBlock           = .false.,                       &
         MatchMultipleLines   = .false.,                       &
         nHeaderLines         = 1,                             &
         ColumnKeys           = Reshape(                       &
                                  (/                           &
                                    'Name                   ', & ! 1
                                    'Min B L Depth          ', & ! 2
                                    'Use NWP BL Depth?      ', & ! 3
                                    'Restore Met Script     ', & ! 4
                                    'Delete Met?            ', & ! 5
                                    'Met Folder             ', & ! 6
                                    'Met Folder Stem        ', & ! 7
                                    'Ensemble Met Folder    ', & ! 8
                                    'Met Folders            ', & ! 9
                                    'Met Definition Name    ', & ! 10
                                    'Max B L Depth          ', & ! 11
                                    'Topography Folder      ', & ! 12
                                    'Mesoscale SigU         ', & ! 13
                                    'Mesoscale TauU         ', & ! 14
                                    'Free Trop SigU         ', & ! 15
                                    'Free Trop SigW         ', & ! 16
                                    'Free Trop TauU         ', & ! 17
                                    'Free Trop TauW         ', & ! 18
                                    'Update On Demand?      '  & ! 19
                                  /),                          &
                                  (/ 19, 1 /)                  &
                                ),                             &
         Defaults             = (/                             &
                                  '       ',     & ! 1
                                  '       ',     & ! 2
                                  '       ',     & ! 3
                                  '       ',     & ! 4
                                  '       ',     & ! 5
                                  '       ',     & ! 6
                                  '       ',     & ! 7
                                  '       ',     & ! 8
                                  '       ',     & ! 9
                                  '       ',     & ! 10
                                  '       ',     & ! 11
                                  '       ',     & ! 12
                                  '       ',                   & ! 13
                                  '       ',                   & ! 14
                                  '0.25   ',                   & ! 15
                                  '0.1    ',                   & ! 16
                                  '300.0  ',                   & ! 17
                                  '100.0  ',                   & ! 18
                                  'No     '                    & ! 19
                                /),                            &
         TwoD                 = .false.,                       &
         UnrecognisedMessages = .true.,                        &
         DisjointColSpecs     = .true.,                        &
         HFBlockForms         = HFBlockForms                   &
       )

End Subroutine NWPMetInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2NWPMet(Tokens, MainOpts, OpenMPOpts, Arrays, MetDefns, Mets)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),        Intent(In)           :: Tokens
  Type(MainOpts_),      Intent(In)           :: MainOpts
  Type(OpenMPOpts_),    Intent(In)           :: OpenMPOpts
  Type(Arrays_),        Intent(In),   Target :: Arrays
  Type(MetDefns_),      Intent(In)           :: MetDefns
  Type(Mets_),          Intent(InOut)        :: Mets
  ! Local parameters:
  Integer, Parameter :: i_Name              = 1
  Integer, Parameter :: i_MinBLDepth        = 2
  Integer, Parameter :: i_UseNWPBLDepth     = 3
  Integer, Parameter :: i_RestoreMetScript  = 4
  Integer, Parameter :: i_DeleteMet         = 5
  Integer, Parameter :: i_MetFolder         = 6
  Integer, Parameter :: i_MetFolderStem     = 7
  Integer, Parameter :: i_EnsembleMetFolder = 8
  Integer, Parameter :: i_MetFolders        = 9
  Integer, Parameter :: i_MetDefinitionName = 10
  Integer, Parameter :: i_MaxBLDepth        = 11
  Integer, Parameter :: i_TopographyFolder  = 12
  Integer, Parameter :: i_MesoscaleSigU     = 13
  Integer, Parameter :: i_MesoscaleTauU     = 14
  Integer, Parameter :: i_FreeTropSigU      = 15
  Integer, Parameter :: i_FreeTropSigW      = 16
  Integer, Parameter :: i_FreeTropTauU      = 17
  Integer, Parameter :: i_FreeTropTauW      = 18
  Integer, Parameter :: i_UpdateOnDemand    = 19
  ! Locals:
  Real(Std)                 :: HMin
  Real(Std)                 :: HMax
  Logical                   :: UseNWPH
  Real(Std)                 :: SigUUM
  Real(Std)                 :: TauUUM
  Real(Std)                 :: SigU2HPlus
  Real(Std)                 :: SigW2HPlus
  Real(Std)                 :: TauUHPlus
  Real(Std)                 :: TauWHPlus
  Logical                   :: DeleteMet
  Integer                   :: iMetDefn
  Integer                   :: iMetDefn2s(MaxNWPMetDefn2s)
  Type(NWPMetDefn2_)        :: SelNWPMetDefn2s(MaxNWPMetDefn2s)
  Integer                   :: iFoldersArray
  Type(Array_), Pointer     :: FoldersArray
  Integer                   :: ArraySize
  Character(MaxTokenLength) :: MetFolders(MaxArrayLength)
  Integer                   :: i
  Integer                   :: j
  Logical                   :: UpdateOnDemand             ! Indicates the met module instance is to be updated
                                                          ! using update-on-demand.
                                                           
  HMin       = Token2Std(                                 &
                 C         = Tokens%Tokens(i_MinBLDepth), &
                 BlockKey  = 'NWP Met Module Instances',  &
                 Item      = Tokens%Tokens(i_Name),       &
                 ColumnKey = 'Min B L Depth'              &
               )

  HMax       = Token2Std(                                 &
                 C         = Tokens%Tokens(i_MaxBLDepth), &
                 BlockKey  = 'NWP Met Module Instances',  &
                 Item      = Tokens%Tokens(i_Name),       &
                 ColumnKey = 'Max B L Depth'              &
               )

  UseNWPH    = Token2Log(                                    &
                 C         = Tokens%Tokens(i_UseNWPBLDepth), &
                 BlockKey  = 'NWP Met Module Instances',     &
                 Item      = Tokens%Tokens(i_Name),          &
                 ColumnKey = 'Use NWP BL Depth?'             &
               )

  SigU2HPlus = (                                              &
                 Token2Std(                                   &
                   C         = Tokens%Tokens(i_FreeTropSigU), &
                   BlockKey  = 'NWP Met Module Instances',    &
                   Item      = Tokens%Tokens(i_Name),         &
                   ColumnKey = 'Free Trop SigU'               &
                 )                                            &
               )**2

  SigW2HPlus = (                                              &
                 Token2Std(                                   &
                   C         = Tokens%Tokens(i_FreeTropSigW), &
                   BlockKey  = 'NWP Met Module Instances',    &
                   Item      = Tokens%Tokens(i_Name),         &
                   ColumnKey = 'Free Trop SigW'               &
                 )                                            &
               )**2

  TauUHPlus  = Token2Std(                                   &
                 C         = Tokens%Tokens(i_FreeTropTauU), &
                 BlockKey  = 'NWP Met Module Instances',    &
                 Item      = Tokens%Tokens(i_Name),         &
                 ColumnKey = 'Free Trop TauU'               &
               )

  TauWHPlus  = Token2Std(                                   &
                 C         = Tokens%Tokens(i_FreeTropTauW), &
                 BlockKey  = 'NWP Met Module Instances',    &
                 Item      = Tokens%Tokens(i_Name),         &
                 ColumnKey = 'Free Trop TauW'               &
               )

  DeleteMet  = Token2Log(                                &
                 C         = Tokens%Tokens(i_DeleteMet), &
                 BlockKey  = 'NWP Met Module Instances', &
                 Item      = Tokens%Tokens(i_Name),      &
                 ColumnKey = 'Delete Met?'               &
               )

  ArraySize     = 0
  MetFolders(:) = ' '

  If (Tokens%Tokens(i_MetFolders) /= ' ') Then

    iFoldersArray =  FindArrayIndex(Tokens%Tokens(i_MetFolders), Arrays)
    FoldersArray  => Arrays%Arrays(iFoldersArray)
    ArraySize     =  FoldersArray%n
    If (FoldersArray%n < 1) Then ! $$ needed? can arrays have 0 elements?
      Call Message(                                                                       &
             'FATAL ERROR in reading the NWP met module instance "'                    // &
             Trim(Tokens%Tokens(i_Name))                                               // &
             '" from the input file(s): the specified array of met folder names has '  // &
             'no entries',                                                                &
             3                                                                            &
           )
    End If
    Do j = 1, FoldersArray%n
      MetFolders(j) = FoldersArray%Array(j)
    End Do

  End If

  iMetDefn  = FindNWPMetDefnIndex(Tokens%Tokens(i_MetDefinitionName), MetDefns)

  Do i = 1, MetDefns%NWPMetDefns(iMetDefn)%nMetDefn2Names
    iMetDefn2s(i)      = FindNWPMetDefn2Index(                             &
                          MetDefns%NWPMetDefns(iMetDefn)%MetDefn2Names(i), &
                          MetDefns                                         &
                        )
    SelNWPMetDefn2s(i) = MetDefns%NWPMetDefn2s(iMetDefn2s(i))
  End Do              

  If ((Tokens%Tokens(i_MesoscaleSigU) == ' ') .neqv. (Tokens%Tokens(i_MesoscaleTauU) == ' ')) Then
    Call Message(                                                                        &
           'FATAL ERROR in reading the NWP met module instance "'                     // &
           Trim(Tokens%Tokens(i_Name))                                                // &
           '" from the input file(s): only one of Mesoscale SigU or Mesoscale TauU '  // &
           'is given',                                                                   &
           3                                                                             &
         )
  End If

  If (Tokens%Tokens(i_MesoscaleSigU) /= ' ') Then
    SigUUM = (Token2Std(                                    &
                C         = Tokens%Tokens(i_MesoscaleSigU), &
                BlockKey  = 'NWP Met Module Instances',     &
                Item      = Tokens%Tokens(i_Name),          &
                ColumnKey = 'Mesoscale SigU'                &
              ))**2
  Else
    SigUUM = MetDefns%NWPMetDefns(iMetDefn)%SigUUM
  End If
  
  If (Tokens%Tokens(i_MesoscaleTauU) /= ' ') Then
    TauUUM = Token2Std(                                    &
               C         = Tokens%Tokens(i_MesoscaleTauU), &
               BlockKey  = 'NWP Met Module Instances',     &
               Item      = Tokens%Tokens(i_Name),          &
               ColumnKey = 'Mesoscale TauU'                &
             )
  Else
    TauUUM = MetDefns%NWPMetDefns(iMetDefn)%TauUUM
  End If
        
  UpdateOnDemand = Token2Log(                                     &
                     C         = Tokens%Tokens(i_UpdateOnDemand), &
                     BlockKey  = 'NWP Met Module Instances',      &
                     Item      = Tokens%Tokens(i_Name),           &
                     ColumnKey = 'Update On Demand?'              &
                   )

  Call AddNWPMet(                                                                                         &
         NWPMet = InitNWPMet(                                                                             &
                    MetName           = Tokens%Tokens(i_Name),                                            &
                    FixedMet          = MainOpts%FixedMet,                                                &
                    UpdateOnDemand    = UpdateOnDemand,                                                   &
                    NWPMetDefn        = MetDefns%NWPMetDefns(iMetDefn),                                   &
                    NWPMetDefn2s      = SelNWPMetDefn2s(1:MetDefns%NWPMetDefns(iMetDefn)%nMetDefn2Names), &
                    HMin              = HMin,                                                             &
                    HMax              = HMax,                                                             &
                    UseNWPH           = UseNWPH,                                                          &
                    SigUUM            = SigUUM,                                                           &
                    TauUUM            = TauUUM,                                                           &
                    SigU2HPlus        = SigU2HPlus,                                                       &
                    SigW2HPlus        = SigW2HPlus,                                                       &
                    TauUHPlus         = TauUHPlus,                                                        &
                    TauWHPlus         = TauWHPlus,                                                        &
                    MetFolder         = Tokens%Tokens(i_MetFolder),                                       &
                    MetFolderStem     = Tokens%Tokens(i_MetFolderStem),                                   &
                    MetFolders        = MetFolders(1:ArraySize),                                          &
                    EnsembleMetFolder = Tokens%Tokens(i_EnsembleMetFolder),                               &
                    TopogFolder       = Tokens%Tokens(i_TopographyFolder),                                &
                    RestoreMetScript  = Tokens%Tokens(i_RestoreMetScript),                                &
                    DeleteMet         = DeleteMet,                                                        &
                    OpenMPOpts        = OpenMPOpts                                                        &
                  ),                                                                                      &
         Mets   = Mets                                                                                    &
       )

End Subroutine Tokens2NWPMet

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine RadarMetDefnInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                        &
         BlockKey             = 'Radar Met Definitions:',            &
         NamedBlock           = .false.,                             &
         MatchMultipleLines   = .false.,                             &
         nHeaderLines         = 1,                                   &
         ColumnKeys           = Reshape(                             &
                                  (/                                 &
                                    'Name                         ', & ! 1
                                    'Binary Format                ', & ! 2
                                    'dT                           ', & ! 3
                                    'T0                           ', & ! 4
                                    'Day Per File                 ', & ! 5
                                    'Suffix                       ', & ! 6
                                    'Prefix                       ', & ! 7
                                    'File Type                    ', & ! 8
                                    'Met File Structure Definition', & ! 9
                                    'H-Grid                       ', & ! 10
                                    'Next Precipitation           '  & ! 11
                                  /),                                &
                                  (/ 11, 1 /)                        &
                                ),                                   &
         Defaults             = (/                                   &
                                  '              ', & ! 1
                                  '              ', & ! 2
                                  '              ', & ! 3
                                  '1/1/2000 00:00',                  & ! 4
                                  '              ', & ! 5
                                  '              ', & ! 6
                                  '              ', & ! 7
                                  '              ', & ! 8
                                  '              ', & ! 9
                                  '              ', & ! 10
                                  'No            '                   & ! 11
                                /),                                  &
         TwoD                 = .false.,                             &
         UnrecognisedMessages = .true.,                              &
         DisjointColSpecs     = .true.,                              &
         HFBlockForms         = HFBlockForms                         &
       )

End Subroutine RadarMetDefnInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2RadarMetDefn(Tokens, MetDefns)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),        Intent(In)    :: Tokens
  Type(RadarMetDefns_), Intent(InOut) :: MetDefns
  ! Local parameters:
  Integer, Parameter :: i_Name                       = 1
  Integer, Parameter :: i_BinaryFormat               = 2
  Integer, Parameter :: i_dT                         = 3
  Integer, Parameter :: i_T0                         = 4
  Integer, Parameter :: i_DayPerFile                 = 5
  Integer, Parameter :: i_Suffix                     = 6
  Integer, Parameter :: i_Prefix                     = 7
  Integer, Parameter :: i_FileType                   = 8
  Integer, Parameter :: i_MetFileStructureDefinition = 9
  Integer, Parameter :: i_HGrid                      = 10
  Integer, Parameter :: i_NextPrecipitation          = 11
  ! Locals:
  Type(Time_)              :: T0           ! Reference time for the first met field
  Type(Time_)              :: Dt           ! Time interval between met fields
  Logical                  :: DayPerFile
  Logical                  :: NextPrecip   ! Use next time step for precip
  Integer                  :: nSuffixs
  Integer                  :: nMetDefn2Names
  Character(MaxCharLength) :: Suffixs(MaxRadarMetDefn2s)
  Character(MaxCharLength) :: MetDefn2Names(MaxRadarMetDefn2s)

  Type(RadarMetDefn_)      :: RadarMetDefn

  T0           = Token2Time(                            &
                   C         = Tokens%Tokens(i_T0),     &
                   BlockKey  = 'Radar Met Definitions', &
                   Item      = Tokens%Tokens(i_Name),   &
                   ColumnKey = 'T0'                     &
                 )

  Dt           = Token2Time(                            &
                   C         = Tokens%Tokens(i_dT),     &
                   BlockKey  = 'Radar Met Definitions', &
                   Item      = Tokens%Tokens(i_Name),   &
                   ColumnKey = 'dT',                    &
                   Interval  = .true.                   &
                 )

  DayPerFile   = Token2Log(                                 &
                   C         = Tokens%Tokens(i_DayPerFile), &
                   BlockKey  = 'Radar Met Definitions',     &
                   Item      = Tokens%Tokens(i_Name),       &
                   ColumnKey = 'Day Per File'               &
                 )

  NextPrecip   = Token2Log(                                        &
                   C         = Tokens%Tokens(i_NextPrecipitation), &
                   BlockKey  = 'Radar Met Definitions',            &
                   Item      = Tokens%Tokens(i_Name),              &
                   ColumnKey = 'Next Precipitation'                &
                 )

  Call ParseListChar(                         &
         C         = Tokens%Tokens(i_Suffix), &
         Delim     = ';',                     &
         BlockKey  = 'Radar Met Definitions', &
         Item      = Tokens%Tokens(i_Name),   &
         ColumnKey = 'Suffix',                &
         nValues   = nSuffixs,                &
         Values    = Suffixs                  &
       )

  Call ParseListChar(                                             &
         C         = Tokens%Tokens(i_MetFileStructureDefinition), &
         Delim     = ';',                                         &
         BlockKey  = 'Radar Met Definitions',                     &
         Item      = Tokens%Tokens(i_Name),                       &
         ColumnKey = 'Met File Structure Definition',             &
         nValues   = nMetDefn2Names,                              &
         Values    = MetDefn2Names                                &
       )

  RadarMetDefn = InitRadarMetDefn(                                   &
                   Name           = Tokens%Tokens(i_Name),           &
                   BinaryFormat   = Tokens%Tokens(i_BinaryFormat),   &
                   FileType       = Tokens%Tokens(i_FileType),       &
                   Prefix         = Tokens%Tokens(i_Prefix),         &
                   Suffixs        = Suffixs(1:nSuffixs),             &
                   DayPerFile     = DayPerFile,                      &
                   T0             = T0,                              &
                   Dt             = Dt,                              &
                   NextPrecip     = NextPrecip,                      &
                   MetDefn2Names  = MetDefn2Names(1:nMetDefn2Names), &
                   HGridName      = Tokens%Tokens(i_HGrid)           &
                 )

  Call AddRadarMetDefn(RadarMetDefn, MetDefns)

End Subroutine Tokens2RadarMetDefn

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine RadarMetDefn2InputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                           &
         BlockKey             = 'Radar Met File Structure Definition:', &
         NamedBlock           = .true.,                                 &
         MatchMultipleLines   = .false.,                                &
         nHeaderLines         = 1,                                      &
         ColumnKeys           = Reshape(                                &
                                  (/                                    &
                                    'Field Name      ',                 & ! 1
                                    'Lowest Level    ',                 & ! 2
                                    'Highest Level   ',                 & ! 3
                                    'Field Code      ',                 & ! 4
                                    '3-d?            ',                 & ! 5
                                    'Field Qualifiers'                  & ! 6
                                  /),                                   &
                                  (/ 6, 1 /)                            &
                                ),                                      &
         Defaults             = (/                                      &
                                  '                   ', & ! 1
                                  '                   ', & ! 2
                                  '                   ', & ! 3
                                  '                   ', & ! 4
                                  '                   ', & ! 5
                                  '                   '  & ! 6
                                /),                                     &
         TwoD                 = .true.,                                 &
         UnrecognisedMessages = .true.,                                 &
         DisjointColSpecs     = .true.,                                 &
         HFBlockForms         = HFBlockForms                            &
       )

End Subroutine RadarMetDefn2InputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2RadarMetDefn2(Tokens, MetDefns)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),        Intent(In)    :: Tokens
  Type(RadarMetDefns_), Intent(InOut) :: MetDefns
  ! Local parameters:
  Integer, Parameter :: i_FieldName       = 1
  Integer, Parameter :: i_LowestLevel     = 2
  Integer, Parameter :: i_HighestLevel    = 3
  Integer, Parameter :: i_FieldCode       = 4
  Integer, Parameter :: i_3d              = 5
  Integer, Parameter :: i_FieldQualifiers = 6
  ! Locals:
  Integer                   :: i
  Integer                   :: iFieldQualifier
  Character(MaxTokenLength) :: FieldNames(MaxRadarMetFields)
  Integer                   :: LowestLevels(MaxRadarMetFields)
  Integer                   :: HighestLevels(MaxRadarMetFields)
  Logical                   :: Top(MaxRadarMetFields)
  Integer                   :: FieldCodes(MaxRadarMetFields)
  Logical                   :: ThreeD(MaxRadarMetFields)
  Logical                   :: Total(MaxRadarMetFields)
  Logical                   :: Dynamic(MaxRadarMetFields)
  Integer                   :: nFieldQualifiers
  Character(MaxCharLength)  :: FieldQualifiers(MaxFieldQualifiers)
  Type(RadarMetDefn2_)      :: RadarMetDefn2
  ! i                :: Loop index.
  ! iFieldQualifier  :: Loop index
  ! FieldNames       :: Names of met fields.
  ! LowestLevels     :: Lowest model levels of met fields.
  ! HighestLevels    :: Highest model levels of met fields.
  ! Top              :: Indicates that highest level of a met field is top of the grid.
  ! FieldCodes       :: Codes of met fields. 0 indicates no code given.
  ! Total            :: Indicates that the total or dynamic cloud stashes are total
  ! Dynamic          :: Indicates that the total or dynamic cloud stashes are dynamic
  ! nFieldQualifiers :: Number of Field Qualifiers for field
  ! FieldQualifiers  :: Collection of Field Qualifiers for field
  ! RadarMetDefn2    :: A radar met definition (part 2 - met file structure definition).
  
  ! $$ The distinction between Total / Dynamic is not used for radar met fields but the code has been
  !    retained to assist in any future merging of these radar met routines with the NWP met routines.

  If (Tokens%nLines > MaxRadarMetFields) Then ! $$ replace by maxRows... Then no need for if-test
    Call Message('ERROR in Tokens2RadarMetDefn2: too many fields ' // &
                 'declared in the met definition block '           // &
                 Trim(Tokens%BlockName),                              &
                 3)
  End If

  Do i = 1, Tokens%nLines

    FieldNames(i)   = Tokens%Tokens2d(i_FieldName, i)

    LowestLevels(i) = Token2Int(                                           &
                        C         = Tokens%Tokens2d(i_LowestLevel, i),     &
                        BlockKey  = 'Radar Met File Structure Definition', &
                        Item      = Tokens%BlockName,                      &
                        ColumnKey = 'Lowest Level'                         &
                      )

    If (Tokens%Tokens2d(i_HighestLevel, i) .CIEq. 'Top') Then
      HighestLevels(i) = 0
      Top(i)           = .true.
    Else
      HighestLevels(i) = Token2Int(                                           &
                           C         = Tokens%Tokens2d(i_HighestLevel, i),    &
                           BlockKey  = 'Radar Met File Structure Definition', &
                           Item      = Tokens%BlockName,                      &
                           ColumnKey = 'Highest Level'                        &
                         )
      Top(i)           = .false.
    End If

    FieldCodes(i) = Token2Int(                                           &
                      C         = Tokens%Tokens2d(i_FieldCode, i),       &
                      BlockKey  = 'Radar Met File Structure Definition', &
                      Item      = Tokens%BlockName,                      &
                      ColumnKey = 'Field Code'                           &
                    )

    ThreeD(i) = Token2Log(                                           &
                  C         = Tokens%Tokens2d(i_3d, i),              &
                  BlockKey  = 'Radar Met File Structure Definition', &
                  Item      = Tokens%BlockName,                      &
                  ColumnKey = '3-d?'                                 &
                )

    Call ParseListChar(                                       &
           C         = Tokens%Tokens2d(i_FieldQualifiers, i), &
           Delim     = ';',                                   &
           BlockKey  = 'Radar Met File Structure Definition', &
           Item      = Tokens%BlockName,                      &
           ColumnKey = 'Field Qualifiers',                    &
           nValues   = nFieldQualifiers,                      &
           Values    = FieldQualifiers                        &
         )
         
    Total(i)   = .false.
    Dynamic(i) = .false.
    Do iFieldQualifier = 1, nFieldQualifiers
      If (FieldQualifiers(iFieldQualifier) .CIEq. 'Total') Then
        Total(i) = .true.
      Else If (FieldQualifiers(iFieldQualifier) .CIEq. 'Dynamic') Then
        Dynamic(i) = .true.
      Else
        Call Message('FATAL ERROR in Tokens2RadarMetDefn2: Field qualifier '// &
                      Trim(FieldQualifiers(iFieldQualifier))                // & 
                     ' not recognised in the met definition block '         // &
                     Trim(Tokens%BlockName),                                   &
                     3)
      End If      
    End Do
    
    If (Total(i) .and. Dynamic(i)) Then
      Call Message('FATAL ERROR in Tokens2RadarMetDefn2: Fields '   // &
                   'declared in the met definition block '          // &
                   Trim(Tokens%BlockName)                           // &
                   ' cannot be both Total and Dynamic',                &
                   3)
    End If
    
  End Do

  RadarMetDefn2 = InitRadarMetDefn2(                                &
                    Name          = Tokens%BlockName,               &
                    FieldNames    = FieldNames   (1:Tokens%nLines), &
                    LowestLevels  = LowestLevels (1:Tokens%nLines), &
                    HighestLevels = HighestLevels(1:Tokens%nLines), &
                    Top           = Top          (1:Tokens%nLines), &
                    FieldCodes    = FieldCodes   (1:Tokens%nLines), &
                    ThreeD        = ThreeD       (1:Tokens%nLines), &
                    Total         = Total        (1:Tokens%nLines)  &
                  )

  Call AddRadarMetDefn2(RadarMetDefn2, MetDefns)

End Subroutine Tokens2RadarMetDefn2

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine RadarMetInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                  &
         BlockKey             = 'Radar Met Module Instances:', &
         NamedBlock           = .false.,                       &
         MatchMultipleLines   = .false.,                       &
         nHeaderLines         = 1,                             &
         ColumnKeys           = Reshape(                       &
                                  (/                           &
                                    'Name                   ', & ! 1
                                    'Restore Met Script     ', & ! 2
                                    'Delete Met?            ', & ! 3
                                    'Met Folder             ', & ! 4
                                    'Met Definition Name    ', & ! 5
                                    'Update On Demand?      '  & ! 6
                                  /),                          &
                                  (/ 6, 1 /)                   &
                                ),                             &
         Defaults             = (/                             &
                                  '      ', & ! 1
                                  '      ', & ! 2
                                  '      ', & ! 3
                                  '      ', & ! 4
                                  '      ', & ! 5
                                  'No    '                     & ! 6
                                /),                            &
         TwoD                 = .false.,                       &
         UnrecognisedMessages = .true.,                        &
         DisjointColSpecs     = .true.,                        &
         HFBlockForms         = HFBlockForms                   &
       )

End Subroutine RadarMetInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2RadarMet(Tokens, MainOpts, MetDefns, Mets)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),        Intent(In)    :: Tokens
  Type(MainOpts_),      Intent(In)    :: MainOpts
  Type(RadarMetDefns_), Intent(In)    :: MetDefns
  Type(Mets_),          Intent(InOut) :: Mets
  ! Local parameters:
  Integer, Parameter :: i_Name              = 1
  Integer, Parameter :: i_RestoreMetScript  = 2
  Integer, Parameter :: i_DeleteMet         = 3
  Integer, Parameter :: i_MetFolder         = 4
  Integer, Parameter :: i_MetDefinitionName = 5
  Integer, Parameter :: i_UpdateOnDemand    = 6
  ! Locals:
  Integer              :: i
  Integer              :: iMetDefn
  Integer              :: nMetDefn2s
  Integer              :: iMetDefn2s(MaxRadarMetDefn2s)
  Type(RadarMetDefn2_) :: SelRadarMetDefn2s(MaxRadarMetDefn2s)
  Logical              :: DeleteMet
  Logical              :: UpdateOnDemand
  ! i                                    :: Loop index
  ! iMetDefn                             :: Index of radar met definition in collection of all met definitions
  ! nMetDefn2s                           :: Number of radar met definition part 2s occuring in the definition
  ! iMetDefn2s(MaxRadarMetDefn2s)        :: Indexes of radar met definition part 2s occuring in the definition
  !                                         with respect to the collection of all met definition part 2s
  ! SelRadarMetDefn2s(MaxRadarMetDefn2s) :: Radar met definition part 2s occuring in the met definition
  ! DeleteMet                            :: Indicates met files are to be deleted from the file system when
  !                                         they are no longer required by this met module instance
  ! UpdateOnDemand                       :: Indicates the met module instance is to be updated
  !                                         using update-on-demand

  DeleteMet      = Token2Log(                                  &
                     C         = Tokens%Tokens(i_DeleteMet),   &
                     BlockKey  = 'Radar Met Module Instances', &
                     Item      = Tokens%Tokens(i_Name),        &
                     ColumnKey = 'Delete Met?'                 &
                   )

  UpdateOnDemand = Token2Log(                                     &
                     C         = Tokens%Tokens(i_UpdateOnDemand), &
                     BlockKey  = 'Radar Met Module Instances',    &
                     Item      = Tokens%Tokens(i_Name),           &
                     ColumnKey = 'Update On Demand?'              &
                   )

  iMetDefn   = FindRadarMetDefnIndex(Tokens%Tokens(i_MetDefinitionName), MetDefns)

  nMetDefn2s = MetDefns%RadarMetDefns(iMetDefn)%nMetDefn2Names

  Do i=1, nMetDefn2s
    iMetDefn2s(i)        = FindRadarMetDefn2Index(                              &
                             MetDefns%RadarMetDefns(iMetDefn)%MetDefn2Names(i), &
                             MetDefns                                           &
                           )
    SelRadarMetDefn2s(i) = MetDefns%RadarMetDefn2s(iMetDefn2s(i))
  End Do   

  Call AddRadarMet(                                               &
         InitRadarMet(                                            &
           MetName           = Tokens%Tokens(i_Name),             &
           FixedMet          = MainOpts%FixedMet,                 &
           UpdateOnDemand    = UpdateOnDemand,                    &
           RadarMetDefn      = MetDefns%RadarMetDefns(iMetDefn),  &
           RadarMetDefn2s    = SelRadarMetDefn2s(1:nMetDefn2s),   &
           MetFolder         = Tokens%Tokens(i_MetFolder),        &
           RestoreMetScript  = Tokens%Tokens(i_RestoreMetScript), &
           DeleteMet         = DeleteMet                          &
         ),                                                       &
         Mets                                                     &
       )

End Subroutine Tokens2RadarMet

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine AncillaryMetDefnInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                        &
         BlockKey             = 'Ancillary Met Definitions:',        &
         NamedBlock           = .false.,                             &
         MatchMultipleLines   = .false.,                             &
         nHeaderLines         = 1,                                   &
         ColumnKeys           = Reshape(                             &
                                  (/                                 &
                                    'Name                         ', & ! 1
                                    'Binary Format                ', & ! 2
                                    'dT                           ', & ! 3
                                    'T0                           ', & ! 4
                                    'Suffix                       ', & ! 5
                                    'Next Heat Flux               ', & ! 6
                                    'Prefix                       ', & ! 7
                                    'File Type                    ', & ! 8
                                    'Met File Structure Definition', & ! 9
                                    'H-Grid                       ', & ! 10
                                    'Next Precipitation           '  & ! 11
                                  /),                                &
                                  (/ 11, 1 /)                        &
                                ),                                   &
         Defaults             = (/                                   &
                                  '              ',                  & ! 1
                                  '              ',                  & ! 2
                                  '              ',                  & ! 3
                                  '1/1/2000 00:00',                  & ! 4
                                  '              ',                  & ! 5
                                  'No            ',                  & ! 6  $$ review default when 'next'
                                  '              ',                  & ! 7     sorted for ancil met
                                  '              ',                  & ! 8
                                  '              ',                  & ! 9
                                  '              ',                  & ! 10
                                  'No            '                   & ! 11
                                /),                                  &
         TwoD                 = .false.,                             &
         UnrecognisedMessages = .true.,                              &
         DisjointColSpecs     = .true.,                              &
         HFBlockForms         = HFBlockForms                         &
       )

End Subroutine AncillaryMetDefnInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2AncillaryMetDefn(Tokens, AncillaryMetDefns)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),            Intent(In)    :: Tokens
  Type(AncillaryMetDefns_), Intent(InOut) :: AncillaryMetDefns
  ! Local parameters:
  Integer, Parameter :: i_Name                       = 1
  Integer, Parameter :: i_BinaryFormat               = 2
  Integer, Parameter :: i_dT                         = 3
  Integer, Parameter :: i_T0                         = 4
  Integer, Parameter :: i_Suffix                     = 5
  Integer, Parameter :: i_NextHeatFlux               = 6
  Integer, Parameter :: i_Prefix                     = 7
  Integer, Parameter :: i_FileType                   = 8
  Integer, Parameter :: i_MetFileStructureDefinition = 9
  Integer, Parameter :: i_HGrid                      = 10
  Integer, Parameter :: i_NextPrecipitation          = 11
  ! Locals:
  Type(Time_)             :: T0           ! Reference time for the first met field
  Type(Time_)             :: Dt           ! Time interval between met fields
  Logical                 :: NextHeatFlux ! Use next time step for heat flux and ustar
  Logical                 :: NextPrecip   ! Use next time step for precip
  Logical                 :: Monthly      !
  Type(AncillaryMetDefn_) :: AncillaryMetDefn

  T0           = Token2Time(                                &
                   C         = Tokens%Tokens(i_T0),         &
                   BlockKey  = 'Ancillary Met Definitions', &
                   Item      = Tokens%Tokens(i_Name),       &
                   ColumnKey = 'T0'                         &
                 )

  If (Tokens%Tokens(i_dT) .CIEq. 'Monthly') Then
    Monthly = .true.
  Else
    Monthly = .false.
    Dt           = Token2Time(                                &
                     C         = Tokens%Tokens(i_dT),         &
                     BlockKey  = 'Ancillary Met Definitions', &
                     Item      = Tokens%Tokens(i_Name),       &
                     ColumnKey = 'dT',                        &
                     Interval  = .true.                       &
                   )
  End If

  NextHeatFlux = Token2Log(                                   &
                   C         = Tokens%Tokens(i_NextHeatFlux), &
                   BlockKey  = 'Ancillary Met Definitions',   &
                   Item      = Tokens%Tokens(i_Name),         &
                   ColumnKey = 'Next Heat Flux'               &
                 )

  NextPrecip   = Token2Log(                                        &
                   C         = Tokens%Tokens(i_NextPrecipitation), &
                   BlockKey  = 'Ancillary Met Definitions',        &
                   Item      = Tokens%Tokens(i_Name),              &
                   ColumnKey = 'Next Precipitation'                &
                 )

  AncillaryMetDefn = InitAncillaryMetDefn(                                         &
                       Name         = Tokens%Tokens(i_Name),                       &
                       BinaryFormat = Tokens%Tokens(i_BinaryFormat),               &
                       FileType     = Tokens%Tokens(i_FileType),                   &
                       Prefix       = Tokens%Tokens(i_Prefix),                     &
                       Suffix       = Tokens%Tokens(i_Suffix),                     &
                       T0           = T0,                                          &
                       Dt           = Dt,                                          &
                       Monthly      = Monthly,                                     &
                       NextHeatFlux = NextHeatFlux,                                &
                       NextPrecip   = NextPrecip,                                  &
                       MetDefn2Name = Tokens%Tokens(i_MetFileStructureDefinition), &
                       HGridName    = Tokens%Tokens(i_HGrid)                       &
                     )

  Call AddAncillaryMetDefn(AncillaryMetDefn, AncillaryMetDefns)

End Subroutine Tokens2AncillaryMetDefn

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine AncillaryMetDefn2InputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                               &
         BlockKey             = 'Ancillary Met File Structure Definition:', &
         NamedBlock           = .true.,                                     &
         MatchMultipleLines   = .false.,                                    &
         nHeaderLines         = 1,                                          &
         ColumnKeys           = Reshape(                                    &
                                  (/                                        &
                                    'Field Name   ',                        & ! 1
                                    'Lowest Level ',                        & ! 2
                                    'Highest Level',                        & ! 3
                                    'Field Code   ',                        & ! 4
                                    '3-d?         '                         & ! 5
                                  /),                                       &
                                  (/ 5, 1 /)                                &
                                ),                                          &
         Defaults             = (/                                          &
                                  '                   ', & ! 1
                                  '                   ', & ! 2
                                  '                   ', & ! 3
                                  '                   ', & ! 4
                                  '                   '  & ! 5
                                /),                                         &
         TwoD                 = .true.,                                     &
         UnrecognisedMessages = .true.,                                     &
         DisjointColSpecs     = .true.,                                     &
         HFBlockForms         = HFBlockForms                                &
       )

End Subroutine AncillaryMetDefn2InputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2AncillaryMetDefn2(Tokens, AncillaryMetDefns)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),            Intent(In)    :: Tokens
  Type(AncillaryMetDefns_), Intent(InOut) :: AncillaryMetDefns
  ! Local parameters:
  Integer, Parameter :: i_FieldName    = 1
  Integer, Parameter :: i_LowestLevel  = 2
  Integer, Parameter :: i_HighestLevel = 3
  Integer, Parameter :: i_FieldCode    = 4
  Integer, Parameter :: i_3d           = 5
  ! Locals:
  Integer                   :: i
  Character(MaxTokenLength) :: FieldNames(MaxAncillaryMetFields)
  Integer                   :: LowestLevels(MaxAncillaryMetFields)
  Integer                   :: HighestLevels(MaxAncillaryMetFields)
  Logical                   :: Top(MaxAncillaryMetFields)
  Integer                   :: FieldCodes(MaxAncillaryMetFields)
  Logical                   :: ThreeD(MaxAncillaryMetFields)
  Type(AncillaryMetDefn2_)  :: AncillaryMetDefn2
  ! i                 :: Loop index.
  ! FieldNames        :: Names of met fields.
  ! FieldCodes        :: Codes of met fields. 0 indicates no code given.
  ! LowestLevels      :: Lowest model levels of met fields.
  ! HighestLevels     :: Highest model levels of met fields.
  ! Top               :: Indicates that highest level of a met field is top of the grid.
  ! AncillaryMetDefn2 :: Ancillary met definition (part 2 - met file structure definition).

  If (Tokens%nLines > MaxAncillaryMetFields) Then ! $$ replace by maxRows... Then no need for if-test
    Call Message('ERROR in Tokens2AncillaryMetDefn2: too many fields ' // &
                 'declared in the met definition block '               // &
                 Trim(Tokens%BlockName),                                  &
                 3)
  End If

  Do i = 1, Tokens%nLines

    FieldNames(i)   = Tokens%Tokens2d(i_FieldName, i)

    LowestLevels(i) = Token2Int(                                               &
                        C         = Tokens%Tokens2d(i_LowestLevel, i),         &
                        BlockKey  = 'Ancillary Met File Structure Definition', &
                        Item      = Tokens%BlockName,                          &
                        ColumnKey = 'Lowest Level'                             &
                      )

    If (Tokens%Tokens2d(i_HighestLevel, i) .CIEq. 'Top') Then
      HighestLevels(i) = 0
      Top(i)           = .true.
    Else
      HighestLevels(i) = Token2Int(                                               &
                           C         = Tokens%Tokens2d(i_HighestLevel, i),        &
                           BlockKey  = 'Ancillary Met File Structure Definition', &
                           Item      = Tokens%BlockName,                          &
                           ColumnKey = 'Highest Level'                            &
                         )
      Top(i)           = .false.
    End If

    FieldCodes(i) = Token2Int(                                               &
                      C         = Tokens%Tokens2d(i_FieldCode, i),           &
                      BlockKey  = 'Ancillary Met File Structure Definition', &
                      Item      = Tokens%BlockName,                          &
                      ColumnKey = 'Field Code'                               &
                    )

    ThreeD(i) = Token2Log(                                               &
                  C         = Tokens%Tokens2d(i_3d, i),                  &
                  BlockKey  = 'Ancillary Met File Structure Definition', &
                  Item      = Tokens%BlockName,                          &
                  ColumnKey = '3-d?'                                     &
                )

  End Do

  AncillaryMetDefn2 = InitAncillaryMetDefn2(                            &
                        Name          = Tokens%BlockName,               &
                        FieldNames    = FieldNames   (1:Tokens%nLines), &
                        LowestLevels  = LowestLevels (1:Tokens%nLines), &
                        HighestLevels = HighestLevels(1:Tokens%nLines), &
                        Top           = Top          (1:Tokens%nLines), &
                        FieldCodes    = FieldCodes   (1:Tokens%nLines), &
                        ThreeD        = ThreeD       (1:Tokens%nLines)  &
                      )

  Call AddAncillaryMetDefn2(AncillaryMetDefn2, AncillaryMetDefns)

End Subroutine Tokens2AncillaryMetDefn2

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine AncillaryMetInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                      &
         BlockKey             = 'Ancillary Met Module Instances:', &
         NamedBlock           = .false.,                           &
         MatchMultipleLines   = .false.,                           &
         nHeaderLines         = 1,                                 &
         ColumnKeys           = Reshape(                           &
                                  (/                               &
                                    'Name                   ',     & ! 1
                                    'Met Folder             ',     & ! 2
                                    'Met Folder Stem        ',     & ! 3
                                    'Ensemble Met Folder    ',     & ! 4
                                    'Met Folders            ',     & ! 5
                                    'Met Definition Name    ',     & ! 6
                                    'Update On Demand?      '      & ! 7
                                  /),                              &
                                  (/ 7, 1 /)                       &
                                ),                                 &
         Defaults             = (/                                 &
                                  '       ',     & ! 1
                                  '       ',     & ! 2
                                  '       ',     & ! 3
                                  '       ',     & ! 4
                                  '       ',     & ! 5
                                  '       ',     & ! 6
                                  'No     '                        & ! 7
                                /),                                &
         TwoD                 = .false.,                           &
         UnrecognisedMessages = .true.,                            &
         DisjointColSpecs     = .true.,                            &
         HFBlockForms         = HFBlockForms                       &
       )

End Subroutine AncillaryMetInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2AncillaryMet(Tokens, MainOpts, Arrays, AncillaryMetDefns, Mets)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),            Intent(In)           :: Tokens
  Type(MainOpts_),          Intent(In)           :: MainOpts
  Type(Arrays_),            Intent(In),   Target :: Arrays
  Type(AncillaryMetDefns_), Intent(In)           :: AncillaryMetDefns
  Type(Mets_),              Intent(InOut)        :: Mets
  ! Local parameters:
  Integer, Parameter :: i_Name              = 1
  Integer, Parameter :: i_MetFolder         = 2
  Integer, Parameter :: i_MetFolderStem     = 3
  Integer, Parameter :: i_EnsembleMetFolder = 4
  Integer, Parameter :: i_MetFolders        = 5
  Integer, Parameter :: i_MetDefinitionName = 6
  Integer, Parameter :: i_UpdateOnDemand    = 7
  ! Locals:
  Integer                   :: iMetDefn
  Integer                   :: iMetDefn2
  Integer                   :: iFoldersArray
  Type(Array_), Pointer     :: FoldersArray
  Integer                   :: ArraySize
  Character(MaxTokenLength) :: MetFolders(MaxArrayLength)
  Integer                   :: j
  Logical                   :: UpdateOnDemand             ! Indicates the met module instance is to be updated
                                                          ! using update-on-demand.

  ArraySize     = 0
  MetFolders(:) = ' '

  If (Tokens%Tokens(i_MetFolders) /= ' ') Then

    iFoldersArray =  FindArrayIndex(Tokens%Tokens(i_MetFolders), Arrays)
    FoldersArray  => Arrays%Arrays(iFoldersArray)
    ArraySize     =  FoldersArray%n
    If (FoldersArray%n < 1) Then ! $$ needed? can arrays have 0 elements?
      Call Message(                                                                       &
             'FATAL ERROR in reading the Ancillary met module instance "'              // &
             Trim(Tokens%Tokens(i_Name))                                               // &
             '" from the input file(s): the specified array of met folder names has '  // &
             'no entries',                                                                &
             3                                                                            &
           )
    End If
    Do j = 1, FoldersArray%n
      MetFolders(j) = FoldersArray%Array(j)
    End Do

  End If

  iMetDefn  = FindAncillaryMetDefnIndex(Tokens%Tokens(i_MetDefinitionName), AncillaryMetDefns)
  iMetDefn2 = FindAncillaryMetDefn2Index(                             &
                AncillaryMetDefns%NWPMetDefns(iMetDefn)%MetDefn2Name, &
                AncillaryMetDefns                                     &
              )

  UpdateOnDemand = Token2Log(                                      &
                     C         = Tokens%Tokens(i_UpdateOnDemand),  &
                     BlockKey  = 'Ancillary Met Module Instances', &
                     Item      = Tokens%Tokens(i_Name),            &
                     ColumnKey = 'Update On Demand?'               &
                   )

  Call AddAncillaryMet(                                                                  &
         AncillaryMet = InitAncillaryMet(                                                &
                          MetName           = Tokens%Tokens(i_Name),                     &
                          FixedMet          = MainOpts%FixedMet,                         &
                          UpdateOnDemand    = UpdateOnDemand,                            &
                          NWPMetDefn        = AncillaryMetDefns%NWPMetDefns(iMetDefn),   &
                          NWPMetDefn2       = AncillaryMetDefns%NWPMetDefn2s(iMetDefn2), &
                          MetFolder         = Tokens%Tokens(i_MetFolder),                &
                          MetFolderStem     = Tokens%Tokens(i_MetFolderStem),            &
                          MetFolders        = MetFolders(1:ArraySize),                   &
                          EnsembleMetFolder = Tokens%Tokens(i_EnsembleMetFolder)         &
                        ),                                                               &
         Mets   = Mets                                                                   &
       )

End Subroutine Tokens2AncillaryMet

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine PrototypeFlowInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                       &
         BlockKey             = 'Prototype Flow Module Instances:', &
         NamedBlock           = .false.,                            &
         MatchMultipleLines   = .false.,                            &
         nHeaderLines         = 1,                                  &
         ColumnKeys           = Reshape(                            &
                                  (/                                &
                                    'Name             ',            & ! 1
                                    'Met Module       ',            & ! 2
                                    'Met              ',            & ! 3
                                    'H-Coord          ',            & ! 4
                                    'Z-Coord          ',            & ! 5
                                    'Domain           ',            & ! 6
                                    'Update On Demand?'             & ! 7
                                  /),                               &
                                  (/ 7, 1 /)                        &
                                ),                                  &
         Defaults             = (/                                  &
                                  '   ',                     & ! 1
                                  '   ',                     & ! 2
                                  '   ',                     & ! 3
                                  '   ',                     & ! 4
                                  '   ',                     & ! 5
                                  '   ',                     & ! 6
                                  'No '                             & ! 7
                                /),                                 &
         TwoD                 = .false.,                            &
         UnrecognisedMessages = .true.,                             &
         DisjointColSpecs     = .true.,                             &
         HFBlockForms         = HFBlockForms                        &
       )

End Subroutine PrototypeFlowInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2PrototypeFlow(Tokens, MainOpts, Flows)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)    :: Tokens
  Type(MainOpts_), Intent(In)    :: MainOpts
  Type(Flows_),    Intent(InOut) :: Flows
  ! Local parameters:
  Integer, Parameter :: i_Name           = 1
  Integer, Parameter :: i_MetModule      = 2
  Integer, Parameter :: i_Met            = 3
  Integer, Parameter :: i_HCoord         = 4
  Integer, Parameter :: i_ZCoord         = 5
  Integer, Parameter :: i_Domain         = 6
  Integer, Parameter :: i_UpdateOnDemand = 7
  ! Locals:
  Logical :: UpdateOnDemand ! Indicates the flow module instance is to be updated using update-on-demand.

  UpdateOnDemand = Token2Log(                                       &
                     C         = Tokens%Tokens(i_UpdateOnDemand),   &
                     BlockKey  = 'Prototype Flow Module Instances', &
                     Item      = Tokens%Tokens(i_Name),             &
                     ColumnKey = 'Update On Demand?'                &
                   )

  Call AddPrototypeFlow(                                                &
         PrototypeFlow = InitPrototypeFlow(                             &
                           FlowName       = Tokens%Tokens(i_Name),      &
                           MetModName     = Tokens%Tokens(i_MetModule), &
                           MetName        = Tokens%Tokens(i_Met),       &
                           HCoordName     = Tokens%Tokens(i_HCoord),    &
                           ZCoordName     = Tokens%Tokens(i_ZCoord),    &
                           DomainName     = Tokens%Tokens(i_Domain),    &
                           FixedMet       = MainOpts%FixedMet,          &
                           UpdateOnDemand = UpdateOnDemand              &
                         ),                                             &
         Flows         = Flows                                          &
       )

End Subroutine Tokens2PrototypeFlow

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine SingleSiteFlowInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                         &
         BlockKey             = 'Single Site Flow Module Instances:', &
         NamedBlock           = .false.,                              &
         MatchMultipleLines   = .false.,                              &
         nHeaderLines         = 1,                                    &
         ColumnKeys           = Reshape(                              &
                                  (/                                  &
                                    'Name             ',              & ! 1
                                    'Met Module       ',              & ! 2
                                    'Met              ',              & ! 3
                                    'Domain           ',              & ! 4
                                    'Update On Demand?'               & ! 5
                                  /),                                 &
                                  (/ 5, 1 /)                          &
                                ),                                    &
         Defaults             = (/                                    &
                                  '    ',                       & ! 1
                                  '    ',                       & ! 2
                                  '    ',                       & ! 3
                                  '    ',                       & ! 4
                                  'No  '                              & ! 5
                                /),                                   &
         TwoD                 = .false.,                              &
         UnrecognisedMessages = .true.,                               &
         DisjointColSpecs     = .true.,                               &
         HFBlockForms         = HFBlockForms                          &
       )

End Subroutine SingleSiteFlowInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2SingleSiteFlow(Tokens, MainOpts, Flows)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)    :: Tokens
  Type(MainOpts_), Intent(In)    :: MainOpts
  Type(Flows_),    Intent(InOut) :: Flows
  ! Local parameters:
  Integer, Parameter :: i_Name           = 1
  Integer, Parameter :: i_MetModule      = 2
  Integer, Parameter :: i_Met            = 3
  Integer, Parameter :: i_Domain         = 4
  Integer, Parameter :: i_UpdateOnDemand = 5
  ! Locals:
  Logical :: UpdateOnDemand ! Indicates the flow module instance is to be updated using update-on-demand.

  UpdateOnDemand = Token2Log(                                         &
                     C         = Tokens%Tokens(i_UpdateOnDemand),     &
                     BlockKey  = 'Single Site Flow Module Instances', &
                     Item      = Tokens%Tokens(i_Name),               &
                     ColumnKey = 'Update On Demand?'                  &
                   )

  Call AddSingleSiteFlow(                                                &
         SingleSiteFlow = InitSingleSiteFlow(                            &
                            FlowName       = Tokens%Tokens(i_Name),      &
                            MetModName     = Tokens%Tokens(i_MetModule), &
                            MetName        = Tokens%Tokens(i_Met),       &
                            DomainName     = Tokens%Tokens(i_Domain),    &
                            FixedMet       = MainOpts%FixedMet,          &
                            UpdateOnDemand = UpdateOnDemand              &
                          ),                                             &
         Flows          = Flows                                          &
       )

End Subroutine Tokens2SingleSiteFlow

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine NWPFlowInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                         &
         BlockKey             = 'NWP Flow Module Instances:',         &
         NamedBlock           = .false.,                              &
         MatchMultipleLines   = .false.,                              &
         nHeaderLines         = 1,                                    &
         ColumnKeys           = Reshape(                              &
                                  (/                                  &
                                    'Name                          ', & ! 1
                                    'Met Module                    ', & ! 2
                                    'Met                           ', & ! 3
                                    'Domain                        ', & ! 4
                                    'Ancillary Met Module Instances', & ! 5
                                    'Update On Demand?             ', & ! 6
                                    'Urban Canopy?                 '  & ! 7
                                  /),                                 &
                                  (/ 7, 1 /)                          &
                                ),                                    &
         Defaults             = (/                                    &
                                  '    ',                      & ! 1
                                  '    ',                      & ! 2
                                  '    ',                      & ! 3
                                  '    ',                      & ! 4
                                  '    ',                      & ! 5
                                  'No  ',                             & ! 6
                                  'No  '                              & ! 7
                                /),                                   &
         TwoD                 = .false.,                              &
         UnrecognisedMessages = .true.,                               &
         DisjointColSpecs     = .true.,                               &
         HFBlockForms         = HFBlockForms                          &
       )

End Subroutine NWPFlowInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2NWPFlow(Tokens, MainOpts, Flows)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)    :: Tokens
  Type(MainOpts_), Intent(In)    :: MainOpts
  Type(Flows_),    Intent(InOut) :: Flows
  ! Local parameters:
  Integer, Parameter :: i_Name                        = 1
  Integer, Parameter :: i_MetModule                   = 2
  Integer, Parameter :: i_Met                         = 3
  Integer, Parameter :: i_Domain                      = 4
  Integer, Parameter :: i_AncillaryMetModuleInstances = 5
  Integer, Parameter :: i_UpdateOnDemand              = 6
  Integer, Parameter :: i_UrbanCanopy                 = 7
  ! Locals:
  Integer                  :: nValues
  Character(MaxCharLength) :: Values(MaxAncillaryMetsPerNWPFlow)
  Logical                  :: UpdateOnDemand                     ! Indicates the flow module instance is to be
                                                                 ! updated using update-on-demand.
  Logical                  :: UrbanCanopy

  Call ParseListChar(                                              &
         C         = Tokens%Tokens(i_AncillaryMetModuleInstances), &
         Delim     = ';',                                          &
         BlockKey  = 'NWP Flow Module Instances',                  &
         Item      = Tokens%Tokens(i_Name),                        &
         ColumnKey = 'Ancillary Met Module Instances',             &
         nValues   = nValues,                                      &
         Values    = Values                                        &
       )

  UpdateOnDemand = Token2Log(                                     &
                     C         = Tokens%Tokens(i_UpdateOnDemand), &
                     BlockKey  = 'NWP Flow Module Instances',     &
                     Item      = Tokens%Tokens(i_Name),           &
                     ColumnKey = 'Update On Demand?'              &
                   )

  UrbanCanopy    = Token2Log(                                  &
                     C         = Tokens%Tokens(i_UrbanCanopy), &
                     BlockKey  = 'NWP Flow Module Instances',  &
                     Item      = Tokens%Tokens(i_Name),        &
                     ColumnKey = 'Urban Canopy?'               &
                   )

  Call AddNWPFlow(                                                   &
         NWPFlow = InitNWPFlow(                                      &
                     FlowName          = Tokens%Tokens(i_Name),      &
                     MetModName        = Tokens%Tokens(i_MetModule), &
                     MetName           = Tokens%Tokens(i_Met),       &
                     AncillaryMetNames = Values(1:nValues),          &
                     DomainName        = Tokens%Tokens(i_Domain),    &
                     FixedMet          = MainOpts%FixedMet,          &
                     UpdateOnDemand    = UpdateOnDemand,             &
                     UrbanCanopy       = UrbanCanopy                 &
                   ),                                                &
         Flows   = Flows                                             &
       )

End Subroutine Tokens2NWPFlow

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine BuildingFlowInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                      &
         BlockKey             = 'Building Flow Module Instances:', &
         NamedBlock           = .false.,                           &
         MatchMultipleLines   = .false.,                           &
         nHeaderLines         = 1,                                 &
         ColumnKeys           = Reshape(                           &
                                  (/                               &
                                    'Name              ',          & ! 1
                                    'H-Coord           ',          & ! 2
                                    'Z-Coord           ',          & ! 3
                                    'Domain            ',          & ! 4
                                    'Length            ',          & ! 5
                                    'Width             ',          & ! 6
                                    'Height            ',          & ! 7
                                    'X                 ',          & ! 8
                                    'Y                 ',          & ! 9
                                    'Angle             ',          & ! 10
                                    'Update Flow Subset',          & ! 11
                                    'Update On Demand? '           & ! 12
                                  /),                              &
                                  (/ 12, 1 /)                      &
                                ),                                 &
         Defaults             = (/                                 &
                                  '           ',            & ! 1
                                  '           ',            & ! 2
                                  '           ',            & ! 3
                                  '           ',            & ! 4
                                  '           ',            & ! 5
                                  '           ',            & ! 6
                                  '           ',            & ! 7
                                  '           ',            & ! 8
                                  '           ',            & ! 9
                                  '           ',            & ! 10
                                  '           ',            & ! 11
                                  'No         '                    & ! 12
                                /),                                &
         TwoD                 = .false.,                           &
         UnrecognisedMessages = .true.,                            &
         DisjointColSpecs     = .true.,                            &
         HFBlockForms         = HFBlockForms                       &
       )

End Subroutine BuildingFlowInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2BuildingFlow(Tokens, MainOpts, Flows)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)    :: Tokens
  Type(MainOpts_), Intent(In)    :: MainOpts
  Type(Flows_),    Intent(InOut) :: Flows
  ! Local parameters:
  Integer, Parameter :: i_Name             = 1
  Integer, Parameter :: i_HCoord           = 2
  Integer, Parameter :: i_ZCoord           = 3
  Integer, Parameter :: i_Domain           = 4
  Integer, Parameter :: i_Length           = 5
  Integer, Parameter :: i_Width            = 6
  Integer, Parameter :: i_Height           = 7
  Integer, Parameter :: i_X                = 8
  Integer, Parameter :: i_Y                = 9
  Integer, Parameter :: i_Angle            = 10
  Integer, Parameter :: i_UpdateFlowSubset = 11
  Integer, Parameter :: i_UpdateOnDemand   = 12
  ! Locals:
  Real(Std) :: Lb
  Real(Std) :: Wb
  Real(Std) :: Hb
  Real(Std) :: Xb
  Real(Std) :: Yb
  Real(Std) :: BuildingOrientation
  Logical   :: UpdateOnDemand      ! Indicates the flow module instance is to be updated using
                                   ! update-on-demand.

  Lb                  = Token2Std(                                      &
                          C         = Tokens%Tokens(i_Length),          &
                          BlockKey  = 'Building Flow Module Instances', &
                          Item      = Tokens%Tokens(i_Name),            &
                          ColumnKey = 'Length'                          &
                        )
  Wb                  = Token2Std(                                      &
                          C         = Tokens%Tokens(i_Width),           &
                          BlockKey  = 'Building Flow Module Instances', &
                          Item      = Tokens%Tokens(i_Name),            &
                          ColumnKey = 'Width'                           &
                        )
  Hb                  = Token2Std(                                      &
                          C         = Tokens%Tokens(i_Height),          &
                          BlockKey  = 'Building Flow Module Instances', &
                          Item      = Tokens%Tokens(i_Name),            &
                          ColumnKey = 'Height'                          &
                        )
  Xb                  = Token2Std(                                      &
                          C         = Tokens%Tokens(i_X),               &
                          BlockKey  = 'Building Flow Module Instances', &
                          Item      = Tokens%Tokens(i_Name),            &
                          ColumnKey = 'X'                               &
                        )
  Yb                  = Token2Std(                                      &
                          C         = Tokens%Tokens(i_Y),               &
                          BlockKey  = 'Building Flow Module Instances', &
                          Item      = Tokens%Tokens(i_Name),            &
                          ColumnKey = 'Y'                               &
                        )
  BuildingOrientation = Token2Std(                                      &
                          C         = Tokens%Tokens(i_Angle),           &
                          BlockKey  = 'Building Flow Module Instances', &
                          Item      = Tokens%Tokens(i_Name),            &
                          ColumnKey = 'Angle'                           &
                        )

  UpdateOnDemand = Token2Log(                                      &
                     C         = Tokens%Tokens(i_UpdateOnDemand),  &
                     BlockKey  = 'Building Flow Module Instances', &
                     Item      = Tokens%Tokens(i_Name),            &
                     ColumnKey = 'Update On Demand?'               &
                   )

  Call AddBuildingFlow(                                                           &
         BuildingFlow = InitBuildingFlow(                                         &
                          FlowName            = Tokens%Tokens(i_Name),            &
                          HCoordName          = Tokens%Tokens(i_HCoord),          &
                          ZCoordName          = Tokens%Tokens(i_ZCoord),          &
                          DomainName          = Tokens%Tokens(i_Domain),          &
                          FixedMet            = MainOpts%FixedMet,                &
                          UpdateOnDemand      = UpdateOnDemand,                   &
                          Lb                  = Lb,                               &
                          Wb                  = Wb,                               &
                          Hb                  = Hb,                               &
                          Xb                  = Xb,                               &
                          Yb                  = Yb,                               &
                          BuildingOrientation = BuildingOrientation,              &
                          UpdateSubsetName    = Tokens%Tokens(i_UpdateFlowSubset) &
                        ),                                                        &
         Flows        = Flows                                                     &
       )

End Subroutine Tokens2BuildingFlow

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine RadarFlowInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                   &
         BlockKey             = 'Radar Flow Module Instances:', &
         NamedBlock           = .false.,                        &
         MatchMultipleLines   = .false.,                        &
         nHeaderLines         = 1,                              &
         ColumnKeys           = Reshape(                        &
                                  (/                            &
                                    'Name             ',        & ! 1
                                    'Met Module       ',        & ! 2
                                    'Met              ',        & ! 3
                                    'Domain           ',        & ! 4
                                    'Update On Demand?'         & ! 5
                                  /),                           &
                                  (/ 5, 1 /)                    &
                                ),                              &
         Defaults             = (/                              &
                                  '   ',  & ! 1
                                  '   ',  & ! 2
                                  '   ',  & ! 3
                                  '   ',  & ! 4
                                  'No '                         & ! 5
                                /),                             &
         TwoD                 = .false.,                        &
         UnrecognisedMessages = .true.,                         &
         DisjointColSpecs     = .true.,                         &
         HFBlockForms         = HFBlockForms                    &
       )

End Subroutine RadarFlowInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2RadarFlow(Tokens, MainOpts, Flows)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)    :: Tokens
  Type(MainOpts_), Intent(In)    :: MainOpts
  Type(Flows_),    Intent(InOut) :: Flows
  ! Local parameters:
  Integer, Parameter :: i_Name           = 1
  Integer, Parameter :: i_MetModule      = 2
  Integer, Parameter :: i_Met            = 3
  Integer, Parameter :: i_Domain         = 4
  Integer, Parameter :: i_UpdateOnDemand = 5
  ! Locals:
  Logical :: UpdateOnDemand ! Indicates the flow module instance is to be updated using update-on-demand.

  UpdateOnDemand = Token2Log(                                     &
                     C         = Tokens%Tokens(i_UpdateOnDemand), &
                     BlockKey  = 'Radar Flow Module Instances',   &
                     Item      = Tokens%Tokens(i_Name),           &
                     ColumnKey = 'Update On Demand?'              &
                   )

  Call AddRadarFlow(                                                &
         RadarFlow = InitRadarFlow(                                 &
                       FlowName       = Tokens%Tokens(i_Name),      &
                       MetModName     = Tokens%Tokens(i_MetModule), &
                       MetName        = Tokens%Tokens(i_Met),       &
                       DomainName     = Tokens%Tokens(i_Domain),    &
                       FixedMet       = MainOpts%FixedMet,          &
                       UpdateOnDemand = UpdateOnDemand              &
                     ),                                             &
         Flows     = Flows                                          &
       )

End Subroutine Tokens2RadarFlow

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine LINCOMFlowInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                    &
         BlockKey             = 'LINCOM Flow Module Instances:', &
         NamedBlock           = .false.,                         &
         MatchMultipleLines   = .false.,                         &
         nHeaderLines         = 1,                               &
         ColumnKeys           = Reshape(                         &
                                  (/                             &
                                    'Name                   ',   & ! 1
                                    'Terrain File           ',   & ! 2
                                    'Roughness File         ',   & ! 3
                                    'Terrain Pertubation?   ',   & ! 4
                                    'Roughness Perturbation?',   & ! 5
                                    'xxxxxxxxxxxxxxxxxxxxxx ',   & ! 6 ! remove $$
                                    'H-Grid                 ',   & ! 7
                                    'Z-Grid                 ',   & ! 8
                                    'X Comp Grid            ',   & ! 9
                                    'Y Comp Grid            ',   & ! 10
                                    'Domain                 ',   & ! 11
                                    'Update Flow Subset     ',   & ! 12
                                    'Lincom Executable      ',   & ! 13
                                    'Update On Demand?      '    & ! 14
                                  /),                            &
                                  (/ 14, 1 /)                    &
                                ),                               &
         Defaults             = (/                               &
                                  '   ',     & ! 1
                                  '   ',     & ! 2
                                  '   ',     & ! 3
                                  '   ',     & ! 4
                                  '   ',     & ! 5
                                  '   ',     & ! 6 ! remove $$
                                  '   ',     & ! 7
                                  '   ',     & ! 8
                                  '   ',     & ! 9
                                  '   ',     & ! 10
                                  '   ',     & ! 11
                                  '   ',     & ! 12
                                  '   ',     & ! 13
                                  'No '                          & ! 14
                                /),                              &
         TwoD                 = .false.,                         &
         UnrecognisedMessages = .true.,                          &
         DisjointColSpecs     = .true.,                          &
         HFBlockForms         = HFBlockForms                     &
       )

End Subroutine LINCOMFlowInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2LINCOMFlow(Tokens, MainOpts, Flows)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)    :: Tokens
  Type(MainOpts_), Intent(In)    :: MainOpts
  Type(Flows_),    Intent(InOut) :: Flows
  ! Local parameters:
  Integer, Parameter :: i_Name                  = 1
  Integer, Parameter :: i_TerrainFile           = 2
  Integer, Parameter :: i_RoughnessFile         = 3
  Integer, Parameter :: i_TerrainPertubation    = 4
  Integer, Parameter :: i_RoughnessPerturbation = 5
  Integer, Parameter :: i_xxxxxxxxxxxxxxxxxxxxx = 6 
  Integer, Parameter :: i_HGrid                 = 7
  Integer, Parameter :: i_ZGrid                 = 8
  Integer, Parameter :: i_XCompGrid             = 9
  Integer, Parameter :: i_YCompGrid             = 10
  Integer, Parameter :: i_Domain                = 11
  Integer, Parameter :: i_UpdateFlowSubset      = 12
  Integer, Parameter :: i_LincomExecutable      = 13
  Integer, Parameter :: i_UpdateOnDemand        = 14
  ! Locals:
  Logical :: TerrainPert
  Logical :: RoughnessPert
  Integer :: nXCompGrid
  Integer :: nYCompGrid
  Logical :: UpdateOnDemand ! Indicates the flow module instance is to be updated using update-on-demand.

  TerrainPert   = Token2Log(                                         &
                    C         = Tokens%Tokens(i_TerrainPertubation), &
                    BlockKey  = 'LINCOM Flow Module Instances',      &
                    Item      = Tokens%Tokens(i_Name),               &
                    ColumnKey = 'Terrain Pertubation?'               &
                  )
  RoughnessPert = Token2Log(                                            &
                    C         = Tokens%Tokens(i_RoughnessPerturbation), &
                    BlockKey  = 'LINCOM Flow Module Instances',         &
                    Item      = Tokens%Tokens(i_Name),                  &
                    ColumnKey = 'Roughness Perturbation?'               &
                  )
  nXCompGrid    = Token2Int(                                    &
                    C         = Tokens%Tokens(i_XCompGrid),     &
                    BlockKey  = 'LINCOM Flow Module Instances', &
                    Item      = Tokens%Tokens(i_Name),          &
                    ColumnKey = 'X Comp Grid'                   &
                  )
  nYCompGrid    = Token2Int(                                    &
                    C         = Tokens%Tokens(i_YCompGrid),     &
                    BlockKey  = 'LINCOM Flow Module Instances', &
                    Item      = Tokens%Tokens(i_Name),          &
                    ColumnKey = 'Y Comp Grid'                   &
                  )

  UpdateOnDemand = Token2Log(                                     &
                     C         = Tokens%Tokens(i_UpdateOnDemand), &
                     BlockKey  = 'LINCOM Flow Module Instances',  &
                     Item      = Tokens%Tokens(i_Name),           &
                     ColumnKey = 'Update On Demand?'              &
                   )

  Call AddLINCOMFlow(                                                         &
         LINCOMFlow = InitLINCOMFlow(                                         &
                        FlowName         = Tokens%Tokens(i_Name),             &
                        TerrainFile      = Tokens%Tokens(i_TerrainFile),      &
                        RoughnessFile    = Tokens%Tokens(i_RoughnessFile),    &
                        HGrid            = Tokens%Tokens(i_HGrid),            &
                        ZGrid            = Tokens%Tokens(i_ZGrid),            &
                        TerrainPert      = TerrainPert,                       &
                        RoughnessPert    = RoughnessPert,                     &
                        nXCompGrid       = nXCompGrid,                        &
                        nYCompGrid       = nYCompGrid,                        &
                        LincomExe        = Tokens%Tokens(i_LincomExecutable), &
                        DomainName       = Tokens%Tokens(i_Domain),           &
                        UpdateSubsetName = Tokens%Tokens(i_UpdateFlowSubset), &
                        FixedMet         = MainOpts%FixedMet,                 &
                        UpdateOnDemand   = UpdateOnDemand                     &
                      ),                                                      &
         Flows      = Flows                                                   &
       )

End Subroutine Tokens2LINCOMFlow

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine FlowOrderInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                      &
         BlockKey             = 'Flow Order:',     &
         NamedBlock           = .true.,            &
         MatchMultipleLines   = .false.,           &
         nHeaderLines         = 1,                 &
         ColumnKeys           = Reshape(           &
                                  (/               &
                                    'Flow Module', & ! 1
                                    'Flow       '  & ! 2
                                  /),              &
                                  (/ 2, 1 /)       &
                                ),                 &
         Defaults             = (/                 &
                                  ' ', & ! 1
                                  ' '  & ! 2
                                /),                &
         TwoD                 = .true.,            &
         UnrecognisedMessages = .true.,            &
         DisjointColSpecs     = .true.,            &
         HFBlockForms         = HFBlockForms       &
       )

End Subroutine FlowOrderInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2FlowOrder(Tokens, Flows)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),  Intent(In)    :: Tokens
  Type(Flows_),   Intent(InOut) :: Flows
  ! Local parameters:
  Integer, Parameter :: i_FlowModule = 1
  Integer, Parameter :: i_Flow       = 2
  ! Locals:
  Type(FlowOrder_) :: FlowOrder

  FlowOrder = InitFlowOrder(                                                  &
                Name        = Tokens%BlockName,                               &
                FlowModName = Tokens%Tokens2d(i_FlowModule, 1:Tokens%nLines), &
                FlowName    = Tokens%Tokens2d(i_Flow,       1:Tokens%nLines)  &
              )
  Call AddFlowOrder(FlowOrder, Flows)

End Subroutine Tokens2FlowOrder

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine FlowSubsetInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                      &
         BlockKey             = 'Flow Subset:',    &
         NamedBlock           = .true.,            &
         MatchMultipleLines   = .false.,           &
         nHeaderLines         = 1,                 &
         ColumnKeys           = Reshape(           &
                                  (/               &
                                    'Flow Module', & ! 1
                                    'Flow       '  & ! 2
                                  /),              &
                                  (/ 2, 1 /)       &
                                ),                 &
         Defaults             = (/                 &
                                  ' ', & ! 1
                                  ' '  & ! 2
                                /),                &
         TwoD                 = .true.,            &
         UnrecognisedMessages = .true.,            &
         DisjointColSpecs     = .true.,            &
         HFBlockForms         = HFBlockForms       &
       )

End Subroutine FlowSubsetInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2FlowSubset(Tokens, Flows)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),  Intent(In)    :: Tokens
  Type(Flows_),   Intent(InOut) :: Flows
  ! Local parameters:
  Integer, Parameter :: i_FlowModule = 1
  Integer, Parameter :: i_Flow       = 2
  ! Locals:
  Type(FlowSubset_) :: FlowSubset

  FlowSubset = InitFlowSubset(                                                 &
                 Name        = Tokens%BlockName,                               &
                 FlowModName = Tokens%Tokens2d(i_FlowModule, 1:Tokens%nLines), &
                 FlowName    = Tokens%Tokens2d(i_Flow      , 1:Tokens%nLines)  &
               )
  Call AddFlowSubset(FlowSubset, Flows)

End Subroutine Tokens2FlowSubset

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine FlowAttribInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                       &
         BlockKey             = 'Flow Attributes:', &
         NamedBlock           = .false.,            &
         MatchMultipleLines   = .false.,            &
         nHeaderLines         = 1,                  &
         ColumnKeys           = Reshape(            &
                                  (/                &
                                    'Name      ',   & ! 1
                                    'Flow Order'    & ! 2
                                  /),               &
                                  (/ 2, 1 /)        &
                                ),                  &
         Defaults             = (/                  &
                                  '  ',     & ! 1
                                  '  '      & ! 2
                                /),                 &
         TwoD                 = .false.,            &
         UnrecognisedMessages = .true.,             &
         DisjointColSpecs     = .true.,             &
         HFBlockForms         = HFBlockForms        &
       )

End Subroutine FlowAttribInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2FlowAttrib(Tokens, Flows)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),  Intent(In)    :: Tokens
  Type(Flows_),   Intent(InOut) :: Flows
  ! Local parameters:
  Integer, Parameter :: i_Name      = 1
  Integer, Parameter :: i_FlowOrder = 2

  Call AddAttrib(                                          &
         Attrib = InitAttrib(                              &
                    Name      = Tokens%Tokens(i_Name),     &
                    OrderName = Tokens%Tokens(i_FlowOrder) &
                  ),                                       &
         Flows  = Flows                                    &
       )

End Subroutine Tokens2FlowAttrib

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine SizeDistsInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                  &
         BlockKey             = 'Particle Size Distribution:', &
         NamedBlock           = .true.,                        &
         MatchMultipleLines   = .false.,                       &
         nHeaderLines         = 1,                             &
         ColumnKeys           = Reshape(                       &
                                  (/                           &
                                    'Diameter Range Boundary', & ! 1
                                    'Cumulative Fraction    ', & ! 2
                                    'Particle Density       ', & ! 3
                                    'Representative Diameter', & ! 4
                                    'Particle Shape         '  & ! 5 
                                  /),                          &
                                  (/ 5, 1 /)                   &
                                ),                             &
         Defaults             = (/                             &
                                  ' ',                         & ! 1
                                  ' ',                         & ! 2
                                  ' ',                         & ! 3
                                  ' ',                         & ! 4
                                  ' '                          & ! 5
                                /),                            &
         TwoD                 = .true.,                        &
         UnrecognisedMessages = .true.,                        &
         DisjointColSpecs     = .true.,                        &
         HFBlockForms         = HFBlockForms                   &
       )

End Subroutine SizeDistsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2SizeDists(Tokens, SizeDists)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_),    Intent(In)    :: Tokens
  Type(SizeDists_), Intent(InOut) :: SizeDists
  ! Local parameters:
  Integer, Parameter :: i_DiameterRangeBoundary  = 1
  Integer, Parameter :: i_CumulativeFrac         = 2
  Integer, Parameter :: i_Density                = 3
  Integer, Parameter :: i_RepresentativeDiameter = 4
  Integer, Parameter :: i_ParticleShape          = 5
  ! Locals:
  Real(Std)       :: DiameterRangeBoundary(MaxLinesTwoD)
  Real(Std)       :: CumulativeFrac(MaxLinesTwoD)
  Real(Std)       :: Density(MaxLinesTwoD)
  Real(Std)       :: RepresentativeDiameter(MaxLinesTwoD)
  Real(Std)       :: ParticleShape(MaxLinesTwoD)
  Logical         :: CumulativeFracPresent
  Logical         :: DensityPresent
  Logical         :: RepresentativeDiameterPresent
  Logical         :: ParticleShapePresent
  Integer         :: i

  If (Tokens%nLines < 2) Then
    Call Message(                                                     &
           'FATAL ERROR in inputting particle size distribution "' // &
           Trim(Tokens%BlockName)                                  // &
           '". There must be at least two lines of data.',            &
           3                                                          &
         )
  End If
  
  Do i = 1, Tokens%nLines
    DiameterRangeBoundary(i) = Token2Std(                                                 &
                                 C         = Tokens%Tokens2d(i_DiameterRangeBoundary, i), &
                                 BlockKey  = 'Particle Size Distribution',                &
                                 Item      = Tokens%BlockName,                            &
                                 ColumnKey = 'Diameter Range Boundary',                   &
                                 i         = i                                            &
                               )
  End Do

  CumulativeFracPresent         = Tokens%Tokens2d(i_CumulativeFrac,         1) /= ' '
  DensityPresent                = Tokens%Tokens2d(i_Density,                1) /= ' '
  RepresentativeDiameterPresent = Tokens%Tokens2d(i_RepresentativeDiameter, 1) /= ' '
  ParticleShapePresent          = Tokens%Tokens2d(i_ParticleShape,          1) /= ' '

  If (CumulativeFracPresent) Then
    Do i = 1, Tokens%nLines
      CumulativeFrac(i) = Token2Std(                                          &
                            C         = Tokens%Tokens2d(i_CumulativeFrac, i), &
                            BlockKey  = 'Particle Size Distribution',         &
                            Item      = Tokens%BlockName,                     &
                            ColumnKey = 'Cumulative Fraction',                &
                            i         = i                                     &
                          )
    End Do
  Else
    Do i = 2, Tokens%nLines
      If (Tokens%Tokens2d(i_CumulativeFrac, i) /= ' ') Then
        Call Message(                                                                   &
               'FATAL ERROR in inputting particle size distribution "'               // &
               Trim(Tokens%BlockName)                                                // &
               '". The value of "Cumulative Fraction" should be blank on all lines ' // &
               'if it is blank on the first line.',                                     &
               3                                                                        &
             )
      End If
    End Do
  End If

  If (DensityPresent) Then
    Do i = 1, Tokens%nLines - 1
      Density(i) = Token2Std(                                   &
                     C         = Tokens%Tokens2d(i_Density, i), &
                     BlockKey  = 'Particle Size Distribution',  &
                     Item      = Tokens%BlockName,              &
                     ColumnKey = 'Particle Density',            &
                     i         = i                              &
                   )
    End Do
    If (Tokens%Tokens2d(i_Density, Tokens%nLines) /= ' ') Then
      Call Message(                                                            &
             'FATAL ERROR in inputting particle size distribution "'        // &
             Trim(Tokens%BlockName)                                         // &
             '". The value of "Density" on the last line should be blank.',    &
             3                                                                 &
           )
    End If
  Else
    Do i = 2, Tokens%nLines
      If (Tokens%Tokens2d(i_Density, i) /= ' ') Then
        Call Message(                                                       &
               'FATAL ERROR in inputting particle size distribution "'   // &
               Trim(Tokens%BlockName)                                    // &
               '". The value of "Density" should be blank on all lines ' // &
               'if it is blank on the first line.',                         &
               3                                                            &
             )
      End If
    End Do
  End If

  If (RepresentativeDiameterPresent) Then
    Do i = 1, Tokens%nLines - 1
      RepresentativeDiameter(i) = Token2Std(                                                  &
                                    C         = Tokens%Tokens2d(i_RepresentativeDiameter, i), &
                                    BlockKey  = 'Particle Size Distribution',                 &
                                    Item      = Tokens%BlockName,                             &
                                    ColumnKey = 'Representative Diameter',                    &
                                    i         = i                                             &
                                  )
    End Do
    If (Tokens%Tokens2d(i_RepresentativeDiameter, Tokens%nLines) /= ' ') Then
      Call Message(                                                                            &
             'FATAL ERROR in inputting particle size distribution "'                        // &
             Trim(Tokens%BlockName)                                                         // &
             '". The value of "Representative Diameter" on the last line should be blank.',    &
             3                                                                                 &
           )
    End If
  Else
    Do i = 2, Tokens%nLines
      If (Tokens%Tokens2d(i_RepresentativeDiameter, i) /= ' ') Then
        Call Message(                                                                       &
               'FATAL ERROR in inputting particle size distribution "'                   // &
               Trim(Tokens%BlockName)                                                    // &
               '". The value of "Representative Diameter" should be blank on all lines ' // &
               'if it is blank on the first line.',                                         &
               3                                                                            &
             )
      End If
    End Do
  End If

  If (ParticleShapePresent) Then            
    Do i = 1, Tokens%nLines - 1
      ParticleShape(i) = Token2Std(                                         &
                           C         = Tokens%Tokens2d(i_ParticleShape, i), &
                           BlockKey  = 'Particle Size Distribution',        &
                           Item      = Tokens%BlockName,                    &
                           ColumnKey = 'Particle Shape',                    &
                           i         = i                                    &
                         )
    End Do
    If (Tokens%Tokens2d(i_ParticleShape, Tokens%nLines) /= ' ') Then
      Call Message(                                                                   &
             'FATAL ERROR in inputting particle size distribution "'               // &
             Trim(Tokens%BlockName)                                                // &
             '". The value of "Particle Shape" on the last line should be blank.',    &
             3                                                                        &
           )
    End If
  Else
    Do i = 2, Tokens%nLines
      If (Tokens%Tokens2d(i_ParticleShape, i) /= ' ') Then
        Call Message(                                                              &
               'FATAL ERROR in inputting particle size distribution "'          // &
               Trim(Tokens%BlockName)                                           // &
               '". The value of "Particle Shape" should be blank on all lines ' // &
               'if it is blank on the first line.',                                &
               3                                                                   &
             )
      End If
    End Do
  End If

  Call AddSizeDist(                                                                                 &
         SizeDist  = InitSizeDist(                                                                  &
                       Name                          = Tokens%BlockName,                            &
                       DiameterRangeBoundary         = DiameterRangeBoundary (1:Tokens%nLines    ), &
                       CumulativeFrac                = CumulativeFrac        (1:Tokens%nLines    ), &
                       Density                       = Density               (1:Tokens%nLines - 1), &
                       ParticleShape                 = ParticleShape         (1:Tokens%nLines - 1), &
                       RepresentativeDiameter        = RepresentativeDiameter(1:Tokens%nLines - 1), &
                       CumulativeFracPresent         = CumulativeFracPresent,                       &
                       DensityPresent                = DensityPresent,                              &
                       ParticleShapePresent          = ParticleShapePresent,                        &
                       RepresentativeDiameterPresent = RepresentativeDiameterPresent                &
                     ),                                                                             &
         SizeDists = SizeDists                                                                      &
       )

End Subroutine Tokens2SizeDists

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine SpeciesInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                        &
         BlockKey             = 'Species:',                          &
         NamedBlock           = .false.,                             &
         MatchMultipleLines   = .false.,                             &
         nHeaderLines         = 1,                                   &
         ColumnKeys           = Reshape(                             &
                                  (/                                 &
                                    'Name                         ', & ! 1
                                    'Typexx                       ', & ! 2 ! $$ remove
                                    'Category                     ', & ! 3
                                    'Molecular Weight             ', & ! 4
                                    'Particle Massxxx             ', & ! 5 ! $$ remove
                                    'Densityxx                    ', & ! 6 ! $$ remove
                                    'Half Life                    ', & ! 7
                                    'UV Loss Rate                 ', & ! 8
                                    'Surface Resistance           ', & ! 9
                                    'Deposition Velocity          ', & ! 10
                                    'Material Unit                ', & ! 11
                                    'Daughter                     ', & ! 12
                                    'Branching Ratio              ', & ! 13
                                    'Set Of Cloud Gamma Parameters', & ! 14
                                    'Land Use Dependent Dry Dep   ', & ! 15
                                    'Mean Aerosol Diameter        ', & ! 16
                                    'A Rain - BC                  ', & ! 17 ! below cloud, rain } wet 
                                    'B Rain - BC                  ', & ! 18 ! below cloud, rain } deposition
                                    'A Snow - BC                  ', & ! 19 ! below cloud, snow } scavenging
                                    'B Snow - BC                  ', & ! 20 ! below cloud, snow } coefficients
                                    'A Rain - IC                  ', & ! 21 ! in cloud, rain    }
                                    'B Rain - IC                  ', & ! 22 ! in cloud, rain    }
                                    'A Snow - IC                  ', & ! 23 ! in cloud, snow    }
                                    'B Snow - IC                  ', & ! 24 ! in cloud, snow    }
                                    'Power Law Decay Exponent     ', & ! 25
                                    'Power Law Decay Delay        '  & ! 26
                                  /),                                &
                                  (/ 26, 1 /)                        &
                                ),                                   &
         Defaults             = (/                                   &
                                  '      ', & ! 1
                                  '      ', & ! 2 ! $$ remove
                                  '      ', & ! 3
                                  '      ', & ! 4
                                  '      ', & ! 5 ! $$ remove
                                  '      ', & ! 6 ! $$ remove
                                  'Stable',                & ! 7
                                  'Stable',                & ! 8
                                  '      ', & ! 9
                                  '      ', & ! 10
                                  '      ', & ! 11
                                  '      ', & ! 12
                                  '      ', & ! 13
                                  '      ', & ! 14
                                  'No    ',                & ! 15
                                  '      ', & ! 16
                                  '      ',                & ! 17
                                  '      ',                & ! 18
                                  '      ',                & ! 19
                                  '      ',                & ! 20
                                  '      ',                & ! 21
                                  '      ',                & ! 22
                                  '      ',                & ! 23
                                  '      ',                & ! 24
                                  '      ',                & ! 25
                                  '      '                 & ! 26
                                /),                        &
         TwoD                 = .false.,                   &
         UnrecognisedMessages = .true.,                    &
         DisjointColSpecs     = .true.,                    &
         HFBlockForms         = HFBlockForms               &
       )

End Subroutine SpeciesInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2Species(Tokens, Specieses, MaterialUnits)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),        Intent(In)    :: Tokens        !
  Type(Specieses_),     Intent(InOut) :: Specieses     ! A collection of specieses.
  Type(MaterialUnits_), Intent(InOut) :: MaterialUnits ! Collection of material units.
  ! Local parameters:
  Integer, Parameter :: i_Name                      = 1
  Integer, Parameter :: i_Typexx                    = 2 
  Integer, Parameter :: i_Category                  = 3
  Integer, Parameter :: i_MolecularWeight           = 4
  Integer, Parameter :: i_ParticleMassxxx           = 5 
  Integer, Parameter :: i_Densityxx                 = 6 
  Integer, Parameter :: i_HalfLife                  = 7
  Integer, Parameter :: i_UVLossRate                = 8
  Integer, Parameter :: i_SurfaceResistance         = 9
  Integer, Parameter :: i_DepositionVelocity        = 10
  Integer, Parameter :: i_MaterialUnit              = 11
  Integer, Parameter :: i_Daughter                  = 12
  Integer, Parameter :: i_BranchingRatio            = 13
  Integer, Parameter :: i_SetOfCloudGammaParameters = 14
  Integer, Parameter :: i_LandUseDependentDryDep    = 15
  Integer, Parameter :: i_MeanAerosolDiameter       = 16
  Integer, Parameter :: i_ARainBC                   = 17
  Integer, Parameter :: i_BRainBC                   = 18
  Integer, Parameter :: i_ASnowBC                   = 19
  Integer, Parameter :: i_BSnowBC                   = 20
  Integer, Parameter :: i_ARainIC                   = 21
  Integer, Parameter :: i_BRainIC                   = 22
  Integer, Parameter :: i_ASnowIC                   = 23
  Integer, Parameter :: i_BSnowIC                   = 24
  Integer, Parameter :: i_PowerLawDecayExponent     = 25
  Integer, Parameter :: i_PowerLawDecayDelay        = 26
  ! Locals:
  Real(Std)      :: MolecularWeight       ! Molecular weight.
  Real(Std)      :: InvHalfLife           ! Reciprocal of the half life.
  Logical        :: HasDaughter           ! Determines whether or not the radionuclide
                                          ! concerned decays to a daughter product.
  Real(Std)      :: BranchingRatio        ! Fraction of the activity of the parent which is
                                          ! available for decay to the daughter product.
  Real(Std)      :: UVLossRate            ! Loss rate (per hour) from UV decay.
  Logical        :: UsePowerLawDecay      !
  Real(Std)      :: PowerLawDecayExponent !
  Real(Std)      :: PowerLawDecayDelay    !
  Logical        :: UseVd                 ! Indicates Vd rather than Rs is to be used.
  Real(Std)      :: Rs                    ! Surface resistance.
  Real(Std)      :: Vd                    ! Deposition velocity.
  Real(Std)      :: AerosolDiameter       ! Aerosol mean diameter (um) for dry deposition parametrisation
  Real(Std)      :: ABelowRain            !] Wet deposition coefficients. Scavenging coefficient is given
  Real(Std)      :: BBelowRain            !] by lambda = A * r^B where r is the rainfall rate in mm/hr.
  Real(Std)      :: ABelowSnow            !] Different coefficients used in-cloud / below-cloud and for
  Real(Std)      :: BBelowSnow            !] snow / rain.
  Real(Std)      :: AInRain               !] 
  Real(Std)      :: BInRain               !] 
  Real(Std)      :: AInSnow               !] 
  Real(Std)      :: BInSnow               !] 
  Logical        :: LandUseDryDep         ! Indicates land use dependent dry deposition scheme to be used
  Type(Species_) :: Species               ! A species.

  MolecularWeight = Token2Std(                                      &
                      C         = Tokens%Tokens(i_MolecularWeight), &
                      BlockKey  = 'Species',                        &
                      Item      = Tokens%Tokens(i_Name),            &
                      ColumnKey = 'Molecular Weight'                &
                    )

  If (Tokens%Tokens(i_HalfLife) .CIEq. 'Stable') Then
    InvHalfLife = 0.0
  Else
    InvHalfLife = Token2Std(                               &
                    C         = Tokens%Tokens(i_HalfLife), &
                    BlockKey  = 'Species',                 &
                    Item      = Tokens%Tokens(i_Name),     &
                    ColumnKey = 'Half Life'                &
                  )
    If (InvHalfLife == 0.0) Call Message('FATAL ERROR: Species half life cannot be given as zero', 3)
    InvHalfLife = 1.0/InvHalfLife
  End If

  If (Tokens%Tokens(i_UVLossRate) .CIEq. 'Stable') Then
    UVLossRate = 0.0
  Else
    UVLossRate = Token2Std(                                 &
                   C         = Tokens%Tokens(i_UVLossRate), &
                   BlockKey  = 'Species',                   &
                   Item      = Tokens%Tokens(i_Name),       &
                   ColumnKey = 'UV Loss Rate'               &
                 )
  End If

  LandUseDryDep = Token2Log(                                             &
                    C         = Tokens%Tokens(i_LandUseDependentDryDep), &
                    BlockKey  = 'Species',                               &
                    Item      = Tokens%Tokens(i_Name),                   &
                    ColumnKey = 'Land use dependent dry dep'             &
                  )
  If (                                                     &
   LandUseDryDep .and. .not.                               &
   (                                                       &
     (Tokens%Tokens(i_Name) .CIEq. 'O3'             ) .or. &
     (Tokens%Tokens(i_Name) .CIEq. 'NO'             ) .or. &
     (Tokens%Tokens(i_Name) .CIEq. 'SULPHUR-DIOXIDE') .or. &
     (Tokens%Tokens(i_Name) .CIEq. 'HYDROGEN'       ) .or. &
     (Tokens%Tokens(i_Name) .CIEq. 'AMMONIA'        ) .or. &
     (Tokens%Tokens(i_Name) .CIEq. 'CH4'            ) .or. &
     (Tokens%Tokens(i_Name) .CIEq. 'HCHO'           ) .or. &
     (Tokens%Tokens(i_Name) .CIEq. 'PAN'            ) .or. &
     (Tokens%Tokens(i_Name) .CIEq. 'NO2'            ) .or. &
     (Tokens%Tokens(i_Name) .CIEq. 'CO'             ) .or. &
     (Tokens%Tokens(i_Name) .CIEq. 'HNO3'           )      &
   )                                                       &
  ) Then

    Call Message(                                                                                      &
           'FATAL ERROR: Land Use dependent dry deposition scheme not available for this species: ' // &
           Trim(Tokens%Tokens(i_Name)),                                                                &
           3                                                                                           &
         )

  End If

  If (LandUseDryDep .and. Tokens%Tokens(i_SurfaceResistance) /= ' ') Then
    Call Message(                                                                                     &
           'FATAL ERROR: Surface resistance should not be given when the land use dry deposition ' // &
           'scheme is invoked for a species',                                                         &
           3                                                                                          &
         )
  End If

  If (LandUseDryDep .and. Tokens%Tokens(i_DepositionVelocity) /= ' ') Then
    Call Message(                                                                                      &
           'FATAL ERROR: Deposition velocity should not be given when the land use dry deposition ' // &
           'scheme is invoked for a species',                                                          &
           3                                                                                           &
         )
  End If

  If (LandUseDryDep .and. Tokens%Tokens(i_MeanAerosolDiameter) /= ' ') Then
    Call Message(                                                                              &
           'FATAL ERROR: Incompatible combination: A mean aerosol diameter has been given ' // &
           'and the land use dry deposition scheme has been invoked for a species',            &
           3                                                                                   &
         )
  End If


  AerosolDiameter = 0.0

  If (LandUseDryDep) Then
    UseVd = .false.
  Else
    If (                                                                                              &
      (                                                                                               &
        Tokens%Tokens(i_SurfaceResistance) /= ' ' .and.                                               &
        (Tokens%Tokens(i_DepositionVelocity) /= ' ' .or. Tokens%Tokens(i_MeanAerosolDiameter) /= ' ') &
      )                                                                                               &
      .or.                                                                                            &
      (Tokens%Tokens(i_DepositionVelocity) /= ' ' .and. Tokens%Tokens(i_MeanAerosolDiameter) /= ' ')  &
    ) Then
      Call Message(                                                                                  &
             'FATAL ERROR: Only one of surface resistance, deposition velocity and mean aerosol ' // &
             'diameter should be given for a species',                                               &
             3                                                                                       &
           )
    Else If (                                           &
      Tokens%Tokens(i_SurfaceResistance)   == ' ' .and. &
      Tokens%Tokens(i_DepositionVelocity)  == ' ' .and. &
      Tokens%Tokens(i_MeanAerosolDiameter) == ' '       &
    ) Then
      UseVd = .true.
      Vd    = 0.0
    Else If (Tokens%Tokens(i_SurfaceResistance) /= ' ') Then
      UseVd = .false.
      Rs = Token2Std(                                        &
             C         = Tokens%Tokens(i_SurfaceResistance), &
             BlockKey  = 'Species',                          &
             Item      = Tokens%Tokens(i_Name),              &
             ColumnKey = 'Surface Resistance'                &
           )
    Else If (Tokens%Tokens(i_DepositionVelocity) /= ' ') Then
      UseVd = .true.
      Vd = Token2Std(                                         &
             C         = Tokens%Tokens(i_DepositionVelocity), &
             BlockKey  = 'Species',                           &
             Item      = Tokens%Tokens(i_Name),               &
             ColumnKey = 'Deposition velocity'                &
           )
    Else If (Tokens%Tokens(i_MeanAerosolDiameter) /= ' ') Then
      UseVd = .false.
      AerosolDiameter = Token2Std(                                          &
                          C         = Tokens%Tokens(i_MeanAerosolDiameter), &
                          BlockKey  = 'Species',                            &
                          Item      = Tokens%Tokens(i_Name),                &
                          ColumnKey = 'Mean aerosol diameter'               &
                        )
    End If
  End If

  If (Tokens%Tokens(i_Daughter) /= ' ') Then
    HasDaughter = .true.
    If (Tokens%Tokens(i_BranchingRatio) /= ' ') Then
      BranchingRatio = Token2Std(                                     &
                         C         = Tokens%Tokens(i_BranchingRatio), &
                         BlockKey  = 'Species',                       &
                         Item      = Tokens%Tokens(i_Name),           &
                         ColumnKey = 'Branching Ratio'                &
                       )
    Else
      Call Message(                                                               &
             'FATAL ERROR: Branching ratio must be given for a daughter product', &
             3                                                                    &
           )
    End If
  Else
    HasDaughter = .false.
    If (Tokens%Tokens(i_BranchingRatio) /= ' ') Then
      Call Message(                                                                                 &
             'FATAL ERROR: Branching ratio must not be given when a daughter product is not given', &
             3                                                                                      &
           )
    End If
    BranchingRatio = 0.0
  End If


  If (HasDaughter .and. InvHalfLife == 0.0) Then
    Call Message(                                                                  &
           'FATAL ERROR: A daughter product cannot be given for a stable species', &
           3                                                                       &
         )
  End If

  If (                                      &
    (Tokens%Tokens(i_ARainBC) == ' ') .and. &
    (Tokens%Tokens(i_BRainBC) == ' ') .and. &
    (Tokens%Tokens(i_ASnowBC) == ' ') .and. &
    (Tokens%Tokens(i_BSnowBC) == ' ') .and. &
    (Tokens%Tokens(i_ARainIC) == ' ') .and. &
    (Tokens%Tokens(i_BRainIC) == ' ') .and. &
    (Tokens%Tokens(i_ASnowIC) == ' ') .and. &
    (Tokens%Tokens(i_BSnowIC) == ' ')       &
  ) Then
    ABelowRain = -1.0
    BBelowRain = -1.0
    ABelowSnow = -1.0
    BBelowSnow = -1.0
    AInRain    = -1.0
    BInRain    = -1.0
    AInSnow    = -1.0
    BInSnow    = -1.0
  Else If (                                 &
    (Tokens%Tokens(i_ARainBC) /= ' ') .and. &
    (Tokens%Tokens(i_BRainBC) /= ' ') .and. &
    (Tokens%Tokens(i_ASnowBC) /= ' ') .and. &
    (Tokens%Tokens(i_BSnowBC) /= ' ') .and. &
    (Tokens%Tokens(i_ARainIC) /= ' ') .and. &
    (Tokens%Tokens(i_BRainIC) /= ' ') .and. &
    (Tokens%Tokens(i_ASnowIC) /= ' ') .and. &
    (Tokens%Tokens(i_BSnowIC) /= ' ')       &
  ) Then
    ABelowRain = Token2Std(                              &
                   C         = Tokens%Tokens(i_ARainBC), &
                   BlockKey  = 'Species',                &
                   Item      = Tokens%Tokens(i_Name),    &
                   ColumnKey = 'A rain - BC'             &
                 )
    BBelowRain = Token2Std(                              &
                   C         = Tokens%Tokens(i_BRainBC), &
                   BlockKey  = 'Species',                &
                   Item      = Tokens%Tokens(i_Name),    &
                   ColumnKey = 'B rain - BC'             &
                 )
    ABelowSnow = Token2Std(                              &
                   C         = Tokens%Tokens(i_ASnowBC), &
                   BlockKey  = 'Species',                &
                   Item      = Tokens%Tokens(i_Name),    &
                   ColumnKey = 'A snow - BC'             &
                 )
    BBelowSnow = Token2Std(                              &
                   C         = Tokens%Tokens(i_BSnowBC), &
                   BlockKey  = 'Species',                &
                   Item      = Tokens%Tokens(i_Name),    &
                   ColumnKey = 'B snow - BC'             &
                 )
    AInRain    = Token2Std(                              &
                   C         = Tokens%Tokens(i_ARainIC), &
                   BlockKey  = 'Species',                &
                   Item      = Tokens%Tokens(i_Name),    &
                   ColumnKey = 'A rain - IC'             &
                 )
    BInRain    = Token2Std(                              &
                   C         = Tokens%Tokens(i_BRainIC), &
                   BlockKey  = 'Species',                &
                   Item      = Tokens%Tokens(i_Name),    &
                   ColumnKey = 'B rain - IC'             &
                 )
    AInSnow    = Token2Std(                              &
                   C         = Tokens%Tokens(i_ASnowIC), &
                   BlockKey  = 'Species',                &
                   Item      = Tokens%Tokens(i_Name),    &
                   ColumnKey = 'A snow - IC'             &
                 )
    BInSnow    = Token2Std(                              &
                   C         = Tokens%Tokens(i_BSnowIC), &
                   BlockKey  = 'Species',                &
                   Item      = Tokens%Tokens(i_Name),    &
                   ColumnKey = 'B snow - IC'             &
                 )
    If (Any((/ ABelowRain, BBelowRain, ABelowSnow, BBelowSnow, AInRain, BInRain, AInSnow, BInSnow /) < 0.0)) Then
      Call Message(                                                           &
             'FATAL ERROR: Negative wet deposition coefficients (A or B) ' // &
             ' specified for specie: '                                     // &
             Trim(Tokens%Tokens(i_Name)),                                     &
             3                                                                &
           )
    End If
  Else
    Call Message(                                                             &
           'FATAL ERROR: Not all wet deposition coefficients (As and Bs) ' // &
           ' have been specified for specie: '                             // &
           Trim(Tokens%Tokens(i_Name)),                                       &
           3                                                                  &
         )
  End If

  UsePowerLawDecay = Tokens%Tokens(i_PowerLawDecayExponent) /= ' ' .or. Tokens%Tokens(i_PowerLawDecayDelay) /= ' '
  If (UsePowerLawDecay) Then
    PowerLawDecayExponent = Token2Std(                                            &
                              C         = Tokens%Tokens(i_PowerLawDecayExponent), &
                              BlockKey  = 'Species',                              &
                              Item      = Tokens%Tokens(i_Name),                  &
                              ColumnKey = 'Power Law Decay Exponent'              &
                            )
    PowerLawDecayDelay    = Token2Std(                                            &
                              C         = Tokens%Tokens(i_PowerLawDecayDelay),    &
                              BlockKey  = 'Species',                              &
                              Item      = Tokens%Tokens(i_Name),                  &
                              ColumnKey = 'Power Law Decay Delay'                 &
                            )
  End If
 
  Call AddSpecies(                                                             &
         InitSpecies(                                                          &
           Name                  = Tokens%Tokens(i_Name),                      &
           Category              = Tokens%Tokens(i_Category),                  &
           MaterialUnits         = MaterialUnits,                              &
           MaterialUnitName      = Tokens%Tokens(i_MaterialUnit),              &
           MolecularWeight       = MolecularWeight,                            &
           InvHalfLife           = InvHalfLife,                                &
           HasDaughter           = HasDaughter,                                &
           Daughter              = Tokens%Tokens(i_Daughter),                  &
           BranchingRatio        = BranchingRatio,                             &
           CloudGammaParams      = Tokens%Tokens(i_SetOfCloudGammaParameters), &
           UVLossRate            = UVLossRate,                                 &
           UsePowerLawDecay      = UsePowerLawDecay,                           &
           PowerLawDecayExponent = PowerLawDecayExponent,                      &
           PowerLawDecayDelay    = PowerLawDecayDelay,                         &
           UseVd                 = UseVd,                                      &
           Rs                    = Rs,                                         &
           Vd                    = Vd,                                         &
           ABelowRain            = ABelowRain,                                 &
           BBelowRain            = BBelowRain,                                 &
           ABelowSnow            = ABelowSnow,                                 &
           BBelowSnow            = BBelowSnow,                                 &
           AInRain               = AInRain,                                    &
           BInRain               = BInRain,                                    &
           AInSnow               = AInSnow,                                    &
           BInSnow               = BInSnow,                                    &
           LandUseDryDep         = LandUseDryDep,                              &
           AerosolDiameter       = AerosolDiameter                             &
         ),                                                                    &
         Specieses                                                             &
       )

End Subroutine Tokens2Species

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine SpeciesUsesInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                                &
         BlockKey             = 'Species Uses:',                             &
         NamedBlock           = .false.,                                     &
         MatchMultipleLines   = .false.,                                     &
         nHeaderLines         = 1,                                           &
         ColumnKeys           = Reshape(                                     &
                                  (/                                         &
                                    'Species                              ', & ! 1
                                    'On Particles?                        ', & ! 2
                                    'On Fields?                           ', & ! 3
                                    'Advect Field?                        ', & ! 4
                                    'Particle Size Distribution For Fields'  & ! 5
                                  /),                                        &
                                  (/ 5, 1 /)                                 &
                                ),                                           &
         Defaults             = (/                                           &
                                  '                                       ', & ! 1
                                  'No                                     ', & ! 2 
                                  'No                                     ', & ! 3 
                                  'No                                     ', & ! 4 
                                  '                                       '  & ! 5
                                /),                                          &
         TwoD                 = .false.,                                     &
         UnrecognisedMessages = .true.,                                      &
         DisjointColSpecs     = .true.,                                      &
         HFBlockForms         = HFBlockForms                                 &
       )

End Subroutine SpeciesUsesInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2SpeciesUses(Tokens, Specieses)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),        Intent(In)    :: Tokens    !
  Type(Specieses_),     Intent(InOut) :: Specieses ! A collection of specieses.
  ! Local parameters:
  Integer, Parameter :: i_Species                           = 1
  Integer, Parameter :: i_OnParticles                       = 2
  Integer, Parameter :: i_OnFields                          = 3
  Integer, Parameter :: i_AdvectField                       = 4
  Integer, Parameter :: i_ParticleSizeDistributionForFields = 5
  ! Locals:
  Logical :: OnParticles !
  Logical :: OnFields    !
  Logical :: AdvectField !

  OnParticles = Token2Log(                                  &
                  C         = Tokens%Tokens(i_OnParticles), &
                  BlockKey  = 'Species Uses',               &
                  Item      = Tokens%Tokens(i_Species),     &
                  ColumnKey = 'On Particles?'               &
                )
  OnFields    = Token2Log(                                  &
                  C         = Tokens%Tokens(i_OnFields),    &
                  BlockKey  = 'Species Uses',               &
                  Item      = Tokens%Tokens(i_Species),     &
                  ColumnKey = 'On Fields?'                  &
                )
  AdvectField = Token2Log(                                  &
                  C         = Tokens%Tokens(i_AdvectField), &
                  BlockKey  = 'Species Uses',               &
                  Item      = Tokens%Tokens(i_Species),     &
                  ColumnKey = 'Advect Field?'               &
                )
 
  Call AddSpeciesUses(                                                       &
         InitSpeciesUses(                                                    &
           SpeciesName  = Tokens%Tokens(i_Species),                          &
           OnParticles  = OnParticles,                                       &
           OnFields     = OnFields,                                          &
           AdvectField  = AdvectField,                                       &
           SizeDistName = Tokens%Tokens(i_ParticleSizeDistributionForFields) &
         ),                                                                  &
         Specieses                                                           &
       )

End Subroutine Tokens2SpeciesUses

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine CloudGammaParamsInputNames(HFBlockForms)
! This block contains photon energy dependent parameters (which are all also
! dependent on the Species) for use in the cloud gamma dose calculation

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                                &
         BlockKey             = 'Cloud Gamma Parameters:',                   &
         NamedBlock           = .true.,                                      &
         MatchMultipleLines   = .false.,                                     &
         nHeaderLines         = 1,                                           &
         ColumnKeys           = Reshape(                                     &
                                  (/                                         &
                                    'Photon Energy                        ', & ! 1
                                    'Photon Intensity                     ', & ! 2
                                    'Linear Attenuation Coefficient       ', & ! 3
                                    'B Build-Up Factor A                  ', & ! 4
                                    'B Build-Up Factor B                  ', & ! 5
                                    'Air Kerma Pu Fluence                 ', & ! 6
                                    'Adult Effective Dose Pu Air Kerma    ', & ! 7
                                    'Adult Thyroid Dose Pu Air Kerma      ', & ! 8
                                    'Adult Lung Dose Pu Air Kerma         ', & ! 9
                                    'Adult Bone Surface Dose Pu Air Kerma '  & ! 10
                                  /),                                        &
                                  (/ 10, 1 /)                                &
                                ),                                           &
         Defaults             = (/                                           &
                                  '      ',                                  & ! 1
                                  '      ',                                  & ! 2
                                  '      ',                                  & ! 3
                                  '      ',                                  & ! 4
                                  '      ',                                  & ! 5
                                  '      ',                                  & ! 6
                                  '      ',                                  & ! 7
                                  '      ',                                  & ! 8
                                  '      ',                                  & ! 9
                                  '      '                                   & ! 10
                                /),                                          &
         TwoD                 = .true.,                                      &
         UnrecognisedMessages = .true.,                                      &
         DisjointColSpecs     = .true.,                                      &
         HFBlockForms         = HFBlockForms                                 &
       )

End Subroutine CloudGammaParamsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2CloudGammaParams(Tokens, CloudGammaParamses)
! .

  Implicit None
  ! Argument list:
  Type(Tokens_),             Intent(In)    :: Tokens              !
  Type(CloudGammaParamses_), Intent(InOut) :: CloudGammaParamses  ! A collection of sets cloud gamma
                                                                  ! parameters.
  ! Local parameters:
  Integer, Parameter :: i_PhotonEnergy                   = 1
  Integer, Parameter :: i_PhotonIntensity                = 2
  Integer, Parameter :: i_LinearAttenuationCoefficient   = 3
  Integer, Parameter :: i_BBuildUpFactorA                = 4
  Integer, Parameter :: i_BBuildUpFactorB                = 5
  Integer, Parameter :: i_AirKermaPuFluence              = 6
  Integer, Parameter :: i_AdultEffectiveDosePuAirKerma   = 7
  Integer, Parameter :: i_AdultThyroidDosePuAirKerma     = 8
  Integer, Parameter :: i_AdultLungDosePuAirKerma        = 9
  Integer, Parameter :: i_AdultBoneSurfaceDosePuAirKerma = 10
  ! Locals:
  Integer                 :: i                                    !  Loop index.
  Real(Std)               :: PhotonEnergy(MaxLinesTwoD)           !  The energy of a photon.
  Real(Std)               :: PhotonIntensity(MaxLinesTwoD)        !  The frequency with which the
                                                                  ! photons are emitted.
  Real(Std)               :: LinearAttenCoeff(MaxLinesTwoD)       !  Parameter defining how the
                                                                  ! photon interacts with
                                                                  !  the matter through which it is passing.
  Real(Std)               :: BBuildUpFactora(MaxLinesTwoD)        !} Berger build-up factor. Contribution from
  Real(Std)               :: BBuildUpFactorb(MaxLinesTwoD)        !} scattered photons.
  Real(Std)               :: AirKermapuFluence(MaxLinesTwoD)      !  Air kerma per unit fluence
                                                                  ! is effectively the
                                                                  !  dose in air incident on a unit area.
  Real(Std)               :: AdEffDosepuAirKerma(MaxLinesTwoD)    !} Adult effective or organ dose per unit
  Real(Std)               :: AdThyDosepuAirKerma(MaxLinesTwoD)    !} dose in air.
  Real(Std)               :: AdLunDosepuAirKerma(MaxLinesTwoD)    !}
  Real(Std)               :: AdBoSDosepuAirKerma(MaxLinesTwoD)    !}
  Type(CloudGammaParams_) :: CloudGammaParams                     !  A collection of a single set of
                                                                  ! cloud gamma parameters.

  Do i = 1, Tokens%nlines

    PhotonEnergy(i)         = Token2Std(                                                         &
                                C         = Tokens%Tokens2d(i_PhotonEnergy, i),                  &
                                BlockKey  = 'Cloud Gamma Parameters',                            &
                                Item      = Tokens%BlockName,                                    &
                                ColumnKey = 'Photon Energy'                                      &
                              )

    PhotonIntensity(i)      = Token2Std(                                                         &
                                C         = Tokens%Tokens2d(i_PhotonIntensity, i),               &
                                BlockKey  = 'Cloud Gamma Parameters',                            &
                                Item      = Tokens%BlockName,                                    &
                                ColumnKey = 'Photon Intensity'                                   &
                              )

    LinearAttenCoeff(i)     = Token2Std(                                                         &
                                C         = Tokens%Tokens2d(i_LinearAttenuationCoefficient, i),  &
                                BlockKey  = 'Cloud Gamma Parameters',                            &
                                Item      = Tokens%BlockName,                                    &
                                ColumnKey = 'Linear Attenuation Coefficient'                     &
                              )

    BBuildUpFactora(i)     = Token2Std(                                                          &
                               C         = Tokens%Tokens2d(i_BBuildUpFactorA, i),                &
                               BlockKey  = 'Cloud Gamma Parameters',                             &
                               Item      = Tokens%BlockName,                                     &
                               ColumnKey = 'GP Build-up Factor b'                                &
                             )

    BBuildUpFactorb(i)     = Token2Std(                                                          &
                               C         = Tokens%Tokens2d(i_BBuildUpFactorB, i),                &
                               BlockKey  = 'Cloud Gamma Parameters',                             &
                               Item      = Tokens%BlockName,                                     &
                               ColumnKey = 'GP Build-up Factor c'                                &
                             )

    AirKermapuFluence(i)   = Token2Std(                                                          &
                               C         = Tokens%Tokens2d(i_AirKermaPuFluence, i),              &
                               BlockKey  = 'Cloud Gamma Parameters',                             &
                               Item      = Tokens%BlockName,                                     &
                               ColumnKey = 'Air kerma pu fluence'                                &
                             )

    AdEffDosepuAirKerma(i) = Token2Std(                                                          &
                               C         = Tokens%Tokens2d(i_AdultEffectiveDosePuAirKerma, i),   &
                               BlockKey  = 'Cloud Gamma Parameters',                             &
                               Item      = Tokens%BlockName,                                     &
                               ColumnKey = 'Adult effective dose pu air kerma'                   &
                             )

    AdThyDosepuAirKerma(i) = Token2Std(                                                          &
                               C         = Tokens%Tokens2d(i_AdultThyroidDosePuAirKerma, i),     &
                               BlockKey  = 'Cloud Gamma Parameters',                             &
                               Item      = Tokens%BlockName,                                     &
                               ColumnKey = 'Adult thyroid dose pu air kerma'                     &
                             )

    AdLunDosepuAirKerma(i) = Token2Std(                                                          &
                               C         = Tokens%Tokens2d(i_AdultLungDosePuAirKerma, i),        &
                               BlockKey  = 'Cloud Gamma Parameters',                             &
                               Item      = Tokens%BlockName,                                     &
                               ColumnKey = 'Adult lung dose pu air kerma'                        &
                             )

    AdBoSDosepuAirKerma(i) = Token2Std(                                                          &
                               C         = Tokens%Tokens2d(i_AdultBoneSurfaceDosePuAirKerma, i), &
                               BlockKey  = 'Cloud Gamma Parameters',                             &
                               Item      = Tokens%BlockName,                                     &
                               ColumnKey = 'Adult bone surface dose pu air kerma'                &
                             )

  End Do

    CloudGammaParams = InitCloudGammaParams(                                         &
                         Name                = Tokens%BlockName,                     &
                         nEnergies           = Tokens%nLines,                        &
                         PhotonEnergy        = PhotonEnergy(1:Tokens%nLines),        &
                         PhotonIntensity     = PhotonIntensity(1:Tokens%nLines),     &
                         LinearAttenCoeff    = LinearAttenCoeff(1:Tokens%nLines),    &
                         BBuildUpFactora     = BBuildUpFactora(1:Tokens%nLines),     &
                         BBuildUpFactorb     = BBuildUpFactorb(1:Tokens%nLines),     &
                         AirKermapuFluence   = AirKermapuFluence(1:Tokens%nLines),   &
                         AdEffDosepuAirKerma = AdEffDosepuAirKerma(1:Tokens%nLines), &
                         AdThyDosepuAirKerma = AdThyDosepuAirKerma(1:Tokens%nLines), &
                         AdLunDosepuAirKerma = AdLunDosepuAirKerma(1:Tokens%nLines), &
                         AdBoSDosepuAirKerma = AdBoSDosepuAirKerma(1:Tokens%nLines)  &
                       )

  Call AddCloudGammaParams(CloudGammaParams, CloudGammaParamses)

End Subroutine Tokens2CloudGammaParams

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine MassLimitInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                              &
         BlockKey             = 'Particle Mass Limits:',   &
         NamedBlock           = .false.,                   &
         MatchMultipleLines   = .false.,                   &
         nHeaderLines         = 1,                         &
         ColumnKeys           = Reshape(                   &
                                  (/                       &
                                    'Species Name       ', & !  1
                                    'Particle Mass Limit'  & !  2
                                  /),                      &
                                  (/ 2, 1 /)               &
                                ),                         &
         Defaults             = (/                         &
                                  ' ', & !  1
                                  ' '  & !  2
                                /),                        &
         TwoD                 = .false.,                   &
         UnrecognisedMessages = .true.,                    &
         DisjointColSpecs     = .true.,                    &
         HFBlockForms         = HFBlockForms               &
       )

End Subroutine MassLimitInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2MassLimit(Tokens, MassLimits)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),     Intent(In)    :: Tokens
  Type(MassLimits_), Intent(InOut) :: MassLimits
  ! Local parameters:
  Integer, Parameter :: i_SpeciesName       = 1
  Integer, Parameter :: i_ParticleMassLimit = 2
  ! Locals:
  Real(Std) :: Limit

  Limit = Token2Std(                                        &
            C         = Tokens%Tokens(i_ParticleMassLimit), &
            BlockKey  = 'Particle Mass Limits',             &
            Item      = Tokens%Tokens(i_SpeciesName),       &
            ColumnKey = 'Particle Mass Limit'               &
          )

  Call AddMassLimit(                                   &
         InitMassLimit(                                &
           SpeciesName = Tokens%Tokens(i_SpeciesName), &
           Limit       = Limit                         &
         ),                                            &
         MassLimits                                    &
       )

End Subroutine Tokens2MassLimit

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine SourcesInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms
  ! Locals:
  Integer :: i

  Call InitAndAddHFBlockForm(                                        &
         BlockKey             = 'Sources:',                          &
         NamedBlock           = .false.,                             &
         MatchMultipleLines   = .false.,                             &
         nHeaderLines         = 1,                                   &
         ColumnKeys           = Reshape(                             &
                                  (/                                 &
                                      'Name                       ', & ! 1
                                      'Shape                      ', & ! 2
                                      'H-Coord                    ', & ! 3
                                      'Z-Coord                    ', & ! 4
                                      'X                          ', & ! 5
                                      'Y                          ', & ! 6
                                      'Z                          ', & ! 7
                                      'dX                         ', & ! 8
                                      'dY                         ', & ! 9
                                      'dZ                         ', & ! 10
                                      'Angle                      ', & ! 11 ! label 'degrees' $$
                                      'dH-Metres?                 ', & ! 12
                                      'dZ-Metres?                 ', & ! 13
                                      'Time Dependency            ', & ! 14
                                      'Plume Rise?                ', & ! 15
                                      'Temperature                ', & ! 16
                                      'Volume Flow Rate           ', & ! 17
                                      'Flow Velocity              ', & ! 18
                                      '# Particles                ', & ! 19
                                      'Max Age                    ', & ! 20
                                      'Top Hat                    ', & ! 21
                                      'Start Time                 ', & ! 22
                                      'Stop Time                  ', & ! 23
                                      'Set Of Locations           ', & ! 24
                                      'Location                   ', & ! 25
                                      'Particle Diameter          ', & ! 26
                                      'Particle Size Distribution ', & ! 27
                                      'Particle Density           ', & ! 28
                                      'Uniform Area?              ', & ! 29
                                      'No Reflect?                ', & ! 30
                                      'Met-Dependent Source Type  ', & ! 31
                                      'H-Grid                     ', & ! 32
                                      'Z-Grid                     ', & ! 33
                                      'Lagrangian-Eulerian Time   ', & ! 34
                                      'Particle Shape             ', & ! 35
                                      'Sedimentation Scheme       ', & ! 36
                                      'Source Groups              ', & ! 37
                                      'Z-Grid Iterative Plume Rise', & ! 38
                                      'Dry Gas Mass Fraction      ', & ! 39
                                      'Water Vapour Mass Fraction ', & ! 40
                                      'Plume Rise Height          ', & ! 41
                                      'Summit Elevation           ', & ! 42
                                      'Distal Fine Ash Fraction   ', & ! 43
                                    (                                &
                                      'Source Strength            ', & ! 43 + 1 to 43 + MaxSpecieses
                                      i = 1, MaxSpecieses            &
                                    )                                &
                                  /),                                &
                                  (/ 43 + MaxSpecieses, 1 /)         &
                                ),                                   &
         Defaults             = (/                                   &
                                    '        ', & ! 1
                                    'Cuboid  ',                      & ! 2
                                    '        ', & ! 3
                                    '        ', & ! 4
                                    '        ', & ! 5
                                    '        ', & ! 6
                                    '        ', & ! 7
                                    '        ', & ! 8
                                    '        ', & ! 9
                                    '        ', & ! 10
                                    '        ', & ! 11
                                    '        ', & ! 12
                                    '        ', & ! 13
                                    '        ', & ! 14
                                    'No      ',                      & ! 15
                                    '        ', & ! 16
                                    '        ', & ! 17
                                    '        ', & ! 18
                                    '        ', & ! 19
                                    'infinity',                      & ! 20
                                    'Yes     ',                      & ! 21
                                    '        ', & ! 22
                                    '        ', & ! 23
                                    '        ', & ! 24
                                    '        ', & ! 25
                                    '        ', & ! 26
                                    '        ', & ! 27
                                    '        ', & ! 28
                                    'No      ',                      & ! 29
                                    'No      ',                      & ! 30
                                    '        ', & ! 31
                                    '        ', & ! 32
                                    '        ', & ! 33
                                    'infinity',                      & ! 34
                                    '        ', & ! 35
                                    '        ', & ! 36
                                    '        ', & ! 37
                                    '        ', & ! 38
                                    '        ', & ! 39
                                    '        ', & ! 40
                                    '        ', & ! 41
                                    '        ', & ! 42
                                    '        ', & ! 43
                                  (                                  &
                                    '        ',                      & ! 43 + 1 to 43 + MaxSpecieses [footnote *)]
                                    i = 1, MaxSpecieses              & 
                                  )                                  &
                                /),                                  &
         TwoD                 = .false.,                             &
         UnrecognisedMessages = .true.,                              &
         DisjointColSpecs     = .true.,                              &
         HFBlockForms         = HFBlockForms                         &
       )

!$$ [footnote *)] could have default for species 1?

End Subroutine SourcesInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2Sources(Tokens, FixedMet, Locationses, Sources)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),      Intent(In)    :: Tokens
  Logical,            Intent(In)    :: FixedMet
  Type(Locationses_), Intent(In)    :: Locationses
  Type(Sources_),     Intent(InOut) :: Sources
  ! Local parameters:
  Integer, Parameter :: i_Name                     = 1
  Integer, Parameter :: i_Shape                    = 2
  Integer, Parameter :: i_HCoord                   = 3
  Integer, Parameter :: i_ZCoord                   = 4
  Integer, Parameter :: i_X                        = 5
  Integer, Parameter :: i_Y                        = 6
  Integer, Parameter :: i_Z                        = 7
  Integer, Parameter :: i_dX                       = 8
  Integer, Parameter :: i_dY                       = 9
  Integer, Parameter :: i_dZ                       = 10
  Integer, Parameter :: i_Angle                    = 11
  Integer, Parameter :: i_dHMetres                 = 12
  Integer, Parameter :: i_dZMetres                 = 13
  Integer, Parameter :: i_TimeDependency           = 14
  Integer, Parameter :: i_PlumeRise                = 15
  Integer, Parameter :: i_Temperature              = 16
  Integer, Parameter :: i_VolumeFlowRate           = 17
  Integer, Parameter :: i_FlowVelocity             = 18
  Integer, Parameter :: i_NParticles               = 19
  Integer, Parameter :: i_MaxAge                   = 20
  Integer, Parameter :: i_TopHat                   = 21
  Integer, Parameter :: i_StartTime                = 22
  Integer, Parameter :: i_StopTime                 = 23
  Integer, Parameter :: i_SetOfLocations           = 24
  Integer, Parameter :: i_Location                 = 25
  Integer, Parameter :: i_ParticleDiameter         = 26
  Integer, Parameter :: i_ParticleSizeDistribution = 27
  Integer, Parameter :: i_ParticleDensity          = 28
  Integer, Parameter :: i_UniformArea              = 29
  Integer, Parameter :: i_NoReflect                = 30
  Integer, Parameter :: i_MetDependentSourceType   = 31
  Integer, Parameter :: i_HGrid                    = 32
  Integer, Parameter :: i_ZGrid                    = 33
  Integer, Parameter :: i_LagrangianEulerianTime   = 34
  Integer, Parameter :: i_ParticleShape            = 35
  Integer, Parameter :: i_ShapeScheme              = 36
  Integer, Parameter :: i_SourceGroups             = 37
  Integer, Parameter :: i_ZGridIterativePlumeRise  = 38
  Integer, Parameter :: i_DryGasFraction0          = 39
  Integer, Parameter :: i_WaterVapourFraction0     = 40
  Integer, Parameter :: i_PlumeRiseHeight          = 41
  Integer, Parameter :: i_SummitElevation          = 42
  Integer, Parameter :: i_DistalFineAshFraction    = 43
  Integer, Parameter :: i_ZerothSourceStrength     = 43
  Integer, Parameter :: i_FirstSourceStrength      = i_ZerothSourceStrength + 1
  Integer, Parameter :: i_LastSourceStrength       = i_ZerothSourceStrength + MaxSpecieses

  ! Locals:
  Integer                           :: SourceType     ! Type of source
  Real(Std)                         :: X(3)
  Real(Std)                         :: dX(3)
  Real(Std)                         :: Angle
  Logical                           :: PlumeRise
  Real(Std)                         :: Temperature
  Real(Std)                         :: VolumeFlowRate
  Real(Std)                         :: FlowVelocity
  Type(Time_)                       :: MaxTravelTime
  Logical                           :: TopHat
  Logical                           :: dHMetres
  Logical                           :: dZMetres
  Logical                           :: UniformArea
  Logical                           :: NoReflect
  Type(Time_)                       :: StartTime
  Type(Time_)                       :: StopTime
  Real(Std)                         :: Density
  Character(MaxCharLength)          :: ShapeScheme
  Real(Std)                         :: ParticleShape
  Logical                           :: DensityPresent
  Logical                           :: ParticleShapePresent
  Real(Std)                         :: Diameter
  Logical                           :: Particulate
  Integer                           :: nSpecieses     !
  Integer                           :: i              ! Loop index.
  Integer                           :: j              ! Loop index.
  Character(MaxTokenLength)         :: HCoordName
  Integer                           :: nSourceGroups
  Character(MaxCharLength)          :: SourceGroupNames(MaxSourceGroupsPerSource)
  Type(ShortTime_)                  :: EulerianTime
  Real(Std)                         :: DryGasFraction0
  Real(Std)                         :: WaterVapourFraction0
  Real(Std)                         :: PlumeRiseHeight
  Real(Std)                         :: SummitElevation
  Real(Std)                         :: DistalFineAshFraction

  DensityPresent       = Tokens%Tokens(i_ParticleDensity) /= ' '
  ParticleShapePresent = Tokens%Tokens(i_ParticleShape)   /= ' '
  
  Call ParseListChar(                               &
         C         = Tokens%Tokens(i_SourceGroups), &
         Delim     = ';',                           &
         BlockKey  = 'Sources',                     &
         Item      = Tokens%Tokens(i_Name),         &
         ColumnKey = 'Source Groups',               &
         nValues   = nSourceGroups,                 &
         Values    = SourceGroupNames               &
       )

  ! Determine the type of source.
  ! Note that met-dependent sources (dust and sea-salt) are considered separately from generic sources
  ! because a more limited set of parameters are input for these special source types and the releases
  ! are generally defined over an area through a H-Grid rather than at a specific position.
  If (Tokens%Tokens(i_MetDependentSourceType) .CIEq. ' ') Then

    SourceType = S_Generic

  Else If (Tokens%Tokens(i_MetDependentSourceType) .CIEq. SourceTypeNames(S_Dust)) Then

    SourceType = S_Dust

  Else If (Tokens%Tokens(i_MetDependentSourceType) .CIEq. SourceTypeNames(S_SeaSalt)) Then

    SourceType = S_SeaSalt

  Else If (Tokens%Tokens(i_MetDependentSourceType) .CIEq. SourceTypeNames(S_IterativePlumeModel)) Then

    SourceType = S_IterativePlumeModel

  Else

    Call Message(                                                                            &
           'Error in a source-block in the input file(s): the type of source is unknown.' // &
           'The "Met-dependent Source Type" should be left blank for generic sources or ' // &
           'specified as either "Dust" or "Sea Salt" for met-dependent sources',             &
           3                                                                                 &
         )

  End If

  If (SourceType == S_Generic) Then

    !
    ! ** Generic source **
    !

    ! Check that H-Grid is not given for a generic source.
    If (Tokens%Tokens(i_HGrid) /= ' ') Then
      Call Message('Error in a source-block in the input file(s): a H-Grid must not be specified ' // &
                   'for a generic source', 3)
    End If
    If (Tokens%Tokens(i_ZGrid) /= ' ') Then
      Call Message('Error in a source-block in the input file(s): a Z-Grid must not be specified ' // &
                   'for a generic source', 3)
    End If

    ! Source location given by coords.
    If (                                           &
      Tokens%Tokens(i_X)              /= ' ' .and. &
      Tokens%Tokens(i_Y)              /= ' ' .and. &
      Tokens%Tokens(i_SetOfLocations) == ' ' .and. &
      Tokens%Tokens(i_Location)       == ' '       &
    ) Then

      X(1) = Token2Std(                           &
               C         = Tokens%Tokens(i_X),    &
               BlockKey  = 'Sources',             &
               Item      = Tokens%Tokens(i_Name), &
               ColumnKey = 'X'                    &
             )

      X(2) = Token2Std(                           &
               C         = Tokens%Tokens(i_Y),    &
               BlockKey  = 'Sources',             &
               Item      = Tokens%Tokens(i_Name), &
               ColumnKey = 'Y'                    &
             )

      HCoordName = Tokens%Tokens(i_HCoord)

    ! Source location given by name.
    Else If (                                      &
      Tokens%Tokens(i_X)              == ' ' .and. &
      Tokens%Tokens(i_Y)              == ' ' .and. &
      Tokens%Tokens(i_SetOfLocations) /= ' ' .and. &
      Tokens%Tokens(i_Location)       /= ' '       &
    ) Then

      ! $$ Need to prevent input of HCoordName via Tokens(3) here, so user doesn't think its being used
      ! for dX, dY. If want option to not use the location coord system or metres aligned to that system,
      ! best to add a dHCoordName input.
      i = FindLocationsIndex(Tokens%Tokens(i_SetOfLocations), Locationses)
      j = FindLocationIndex(Tokens%Tokens(i_Location), Locationses%Locationses(i))
      X(1)       = Locationses%Locationses(i)%X(j)
      X(2)       = Locationses%Locationses(i)%Y(j)
      HCoordName = Locationses%Locationses(i)%HCoordNames(j)

    Else

      Call Message(                                                                   &
             'Error in a source-block in the input file(s): '                      // &
             'the source location must be specified by giving '                    // &
             'either "X", "Y" and "H-Coord" or "Set of locations" and "Location"',    &
             3                                                                        &
           )

    End If

    X(3) = Token2Std(                           &
             C         = Tokens%Tokens(i_Z),    &
             BlockKey  = 'Sources',             &
             Item      = Tokens%Tokens(i_Name), &
             ColumnKey = 'Z'                    &
           )

    dX(1) = Token2Std(                           &
              C         = Tokens%Tokens(i_dX),   &
              BlockKey  = 'Sources',             &
              Item      = Tokens%Tokens(i_Name), &
              ColumnKey = 'dX'                   &
            )

    dX(2) = Token2Std(                           &
              C         = Tokens%Tokens(i_dY),   &
              BlockKey  = 'Sources',             &
              Item      = Tokens%Tokens(i_Name), &
              ColumnKey = 'dY'                   &
            )

    dX(3) = Token2Std(                           &
              C         = Tokens%Tokens(i_dZ),   &
              BlockKey  = 'Sources',             &
              Item      = Tokens%Tokens(i_Name), &
              ColumnKey = 'dZ'                   &
            )

    Angle = Token2Std(                            &
              C         = Tokens%Tokens(i_Angle), &
              BlockKey  = 'Sources',              &
              Item      = Tokens%Tokens(i_Name),  &
              ColumnKey = 'Angle'                 &
            )

    dHMetres = Token2Log(                               &
                 C         = Tokens%Tokens(i_dHMetres), &
                 BlockKey  = 'Sources',                 &
                 Item      = Tokens%Tokens(i_Name),     &
                 ColumnKey = 'dH-Metres?'               &
               )

    dZMetres = Token2Log(                               &
                 C         = Tokens%Tokens(i_dZMetres), &
                 BlockKey  = 'Sources',                 &
                 Item      = Tokens%Tokens(i_Name),     &
                 ColumnKey = 'dZ-Metres?'               &
               )

    PlumeRise = Token2Log(                                &
                  C         = Tokens%Tokens(i_PlumeRise), &
                  BlockKey  = 'Sources',                  &
                  Item      = Tokens%Tokens(i_Name),      &
                  ColumnKey = 'Plume Rise?'               &
                )

    If (PlumeRise) Then

      If (Tokens%Tokens(i_Temperature) /= ' ') Then
        Temperature = Token2Std(                                  &
                        C         = Tokens%Tokens(i_Temperature), &
                        BlockKey  = 'Sources',                    &
                        Item      = Tokens%Tokens(i_Name),        &
                        ColumnKey = 'Temperature'                 &
                      )
      Else
        Call Message('Error: must specify temperature of source for plume rise', 3)
      End If
      If (                                            &
        (Tokens%Tokens(i_Shape) .CIEq. 'Cuboid') .or. &
        dX(1) /= dX(2) .or.                           &
        dX(3) /= 0.0                                  &
      ) Then
        Call Message(                                                      &
               'Error: plume rise only possible for Ellipsoid sources ' // &
               'with dX = dY and dZ = 0 at present',                       &
               3                                                           &
             )
      End If

      If (Tokens%Tokens(i_VolumeFlowRate) /= ' ') Then
        If (Tokens%Tokens(i_FlowVelocity) /= ' ') Then
          Call Message('Error: cannot specify VolumeFlowRate and FlowVelocity simultaneously', 3)
        Else
          VolumeFlowRate = Token2Std(                                     &
                             C         = Tokens%Tokens(i_VolumeFlowRate), &
                             BlockKey  = 'Sources',                       &
                             Item      = Tokens%Tokens(i_Name),           &
                             ColumnKey = 'Volume Flow Rate'               &
                           )
          FlowVelocity = VolumeFlowRate/(0.25*Pi*dX(1)*dX(2)) ! DX might not be metres $$
        End If                                                ! Better to convert later?
      Else If (Tokens%Tokens(i_FlowVelocity) /= ' ') Then
        FlowVelocity = Token2Std(                                   &
                         C         = Tokens%Tokens(i_FlowVelocity), &
                         BlockKey  = 'Sources',                     &
                         Item      = Tokens%Tokens(i_Name),         &
                         ColumnKey = 'Flow Velocity'                &
                       )
        VolumeFlowRate = FlowVelocity*(0.25*Pi*dX(1)*dX(2))
      Else
        Call Message('Error: flow rate or velocity needed for plume rise', 3)
      End If

    Else
      Temperature    = -999.0
      FlowVelocity   = -999.0
      VolumeFlowRate = -999.0
    End If

    MaxTravelTime = Token2Time(                            &
                      C         = Tokens%Tokens(i_MaxAge), &
                      BlockKey  = 'Sources',               &
                      Item      = Tokens%Tokens(i_Name),   &
                      ColumnKey = 'Max Age',               &
                      Interval  = .true.                   &
                    )

    TopHat = Token2Log(                             &
               C         = Tokens%Tokens(i_TopHat), &
               BlockKey  = 'Sources',               &
               Item      = Tokens%Tokens(i_Name),   &
               ColumnKey = 'Top Hat'                &
             )

    StartTime = Token2Time(                               &
                  C         = Tokens%Tokens(i_StartTime), &
                  BlockKey  = 'Sources',                  &
                  Item      = Tokens%Tokens(i_Name),      &
                  ColumnKey = 'Start Time'                &
                )

    StopTime  = Token2Time(                              &
                  C         = Tokens%Tokens(i_StopTime), &
                  BlockKey  = 'Sources',                 &
                  Item      = Tokens%Tokens(i_Name),     &
                  ColumnKey = 'Stop Time'                &
                )

    nSpecieses = 0
    Do i = 1, MaxSpecieses
      If (Tokens%Tokens(i_ZerothSourceStrength + i) == ' ') Exit
      nSpecieses = nSpecieses + 1
    End Do

    If (                                                    &
      Tokens%Tokens(i_ParticleSizeDistribution) /= ' ' .or. &
      Tokens%Tokens(i_ParticleDiameter)         /= ' ' .or. &
      Tokens%Tokens(i_ParticleDensity)          /= ' ' .or. &
      Tokens%Tokens(i_ParticleShape)            /= ' ' .or. &
      Tokens%Tokens(i_ShapeScheme)              /= ' '      &
    ) Then
      Particulate = .true.
    Else
      Particulate = .false.
    End If

    If (Particulate) Then                      
                  
        ! If no PSD then diameter must be defined.
        If (Tokens%Tokens(i_ParticleSizeDistribution) == ' ') Then
          Diameter = Token2Std(                                      &
                       C        = Tokens%Tokens(i_ParticleDiameter), &
                       BlockKey = 'Sources',                         &
                       Item     = Tokens%Tokens(i_Name),             &
                       ColumnKey = 'Particle Diameter'               &
                     )
        End If
        
        If (Tokens%Tokens(i_ParticleDensity) /= ' ') Then            
          Density = Token2Std(                                     &
                      C        = Tokens%Tokens(i_ParticleDensity), &
                      BlockKey = 'Sources',                        &
                      Item     = Tokens%Tokens(i_Name),            &
                      ColumnKey = 'Particle Density'               &
                    )
          If (Density <= 0.0) Then
            Call Message(                                   &
                  'FATAL ERROR:  Particle density is <= 0', &
                   3                                        &
                 )
          End If
        End If
        
        If (Tokens%Tokens(i_ParticleShape) /= ' ') Then
          ParticleShape = Token2Std(                                         &
                            C        = Tokens%Tokens(i_ParticleShape),       &
                            BlockKey = 'Sources',                            &
                            Item     = Tokens%Tokens(i_Name),                &
                            ColumnKey = 'Particle Shape'                     &
                          )
          If (ParticleShape > 1.0 .or. ParticleShape <= 0.0) Then
            Call Message(                                        &
                   'FATAL ERROR: Particle shape is > 1 or <= 0', &
                   3                                             &
                 )
          End If
        End If
        
        If (Tokens%Tokens(i_ShapeScheme) .CIEq. 'Ganser') Then
          ShapeScheme = 'Ganser'
        Else If (Tokens%Tokens(i_ShapeScheme) .CIEq. 'WH') Then
          ShapeScheme = 'WH'
        Else If (Tokens%Tokens(i_ShapeScheme) .CIEq. 'White') Then
          ShapeScheme = 'White'
          Call Message(                                                         &
                 'WARNING: When using the "White" sedimentation scheme '     // &
                 'the particle shape is necessarily equal to 1"',               &
                 1                                                              &
               )
        Else If (Tokens%Tokens(i_ShapeScheme) .CIEq. ' ') Then
          ShapeScheme = ' '
        Else
          Call Message(                                                           &
                 'FATAL ERROR: Error in a source-block in the input file(s): ' // &
                 'the type of sedimentation scheme is unknown. '               // &
                 'The "Sedimentation Scheme" should be left blank or '         // &
                 'specified as either "White" or "Ganser" or "WH"',               &
                 3                                                                &
               )
        End If
        
    End If
    
    UniformArea = Token2Log(                                  &
                    C         = Tokens%Tokens(i_UniformArea), &
                    BlockKey  = 'Sources',                    &
                    Item      = Tokens%Tokens(i_Name),        &
                    ColumnKey = 'Uniform Area?'               &
                  )

    NoReflect   = Token2Log(                                &
                    C         = Tokens%Tokens(i_NoReflect), &
                    BlockKey  = 'Sources',                  &
                    Item      = Tokens%Tokens(i_Name),      &
                    ColumnKey = 'No Reflect?'               &
                  )
 
  Else If (SourceType == S_Dust .Or. SourceType == S_SeaSalt) Then

    !
    ! ** Dust/sea-salt source **
    ! Note: the following subset of source parameters is given as user input here
    !   Name
    !   Met-dependent Source Type
    !   H-Grid
    !   Z-Coord (= 'm agl')
    !   dZ-Metres?
    !   Z (= 0)
    !   dZ (= 0)
    !   Shape (= Cuboid, the default)
    !   TopHat (= Yes, the default)
    !   Uniform Area (default = No)
    !   No Reflect (default = No)
    !   Start Time
    !   Stop Time
    !   Max Age
    !   # Particles
    !   Source Strength
    !   Plume Rise (= No, the default)
    !   Particle Density
    !   Particle Size Distribution
    ! Other tokens should be left blank. The corresponding source variables are initialised
    ! but these set values should not be used for met-dependent sources.

    ! No explicit specification of a location. Coord system given via H-Grid.
    If (                                           &
      Tokens%Tokens(i_HCoord        ) == ' ' .and. &
      Tokens%Tokens(i_X             ) == ' ' .and. &
      Tokens%Tokens(i_Y             ) == ' ' .and. &
      Tokens%Tokens(i_dX            ) == ' ' .and. &
      Tokens%Tokens(i_dY            ) == ' ' .and. &
      Tokens%Tokens(i_Angle         ) == ' ' .and. &
      Tokens%Tokens(i_dHMetres      ) == ' ' .and. &
      Tokens%Tokens(i_SetOfLocations) == ' ' .and. &
      Tokens%Tokens(i_Location      ) == ' ' .and. &
      Tokens%Tokens(i_HGrid         ) /= ' '       &
    ) Then

    Else

      Call Message(                                                                                   &
             'Error in a source-block in the input file(s): '                                      // &
             'the location of a dust or sea-salt source can only be specified using a "H-Grid", '  // &
             'with all other location variables ("X", "Y", etc.) kept blank',                         &
             3                                                                                        &
           )

    End If

    HCoordName = Tokens%Tokens(i_HCoord)

    ! Release height.
    X(3)     = Token2Std(                           &
                 C         = Tokens%Tokens(i_Z),    &
                 BlockKey  = 'Sources',             &
                 Item      = Tokens%Tokens(i_Name), &
                 ColumnKey = 'Z'                    &
               )

    dX(3)    = Token2Std(                           &
                 C         = Tokens%Tokens(i_dZ),   &
                 BlockKey  = 'Sources',             &
                 Item      = Tokens%Tokens(i_Name), &
                 ColumnKey = 'dZ'                   &
               )

    dZMetres = Token2Log(                               &
                 C         = Tokens%Tokens(i_dZMetres), &
                 BlockKey  = 'Sources',                 &
                 Item      = Tokens%Tokens(i_Name),     &
                 ColumnKey = 'dZ-Metres?'               &
               )

    If (Tokens%Tokens(i_ZGrid) /= ' ') Then
      Call Message('Error in a source-block in the input file(s): a Z-Grid must not be specified ' // &
                   'for a dust or sea-salt source', 3)
    End If

    ! Time dependency not applicable for dust and sea-salt sources
    If (Tokens%Tokens(i_TimeDependency) /= ' ') Then
      Call Message('Error: "Time Dependency" should be blank for a dust or sea-salt source', 3)
    End If

    ! Plume rise not applicable for dust and sea-salt sources
    PlumeRise = Token2Log(                                &
                  C         = Tokens%Tokens(i_PlumeRise), &
                  BlockKey  = 'Sources',                  &
                  Item      = Tokens%Tokens(i_Name),      &
                  ColumnKey = 'Plume Rise?'               &
                )

    If (PlumeRise) Then
      Call Message('Error: plume rise cannot be used with a dust or sea-salt source', 3)
    End If

    Temperature    = -999.0
    FlowVelocity   = -999.0
    VolumeFlowRate = -999.0

    ! Set defaults (not used, but needed by InitSource routine)
!    X(:)  = 0.0
!    dX(:) = 0.0
!    Angle = 0.0
    dHMetres = .false.
!    dZMetres = .true.

    MaxTravelTime = Token2Time(                            &
                      C         = Tokens%Tokens(i_MaxAge), &
                      BlockKey  = 'Sources',               &
                      Item      = Tokens%Tokens(i_Name),   &
                      ColumnKey = 'Max Age',               &
                      Interval  = .true.                   &
                    )

    TopHat = Token2Log(                             &
               C         = Tokens%Tokens(i_TopHat), &
               BlockKey  = 'Sources',               &
               Item      = Tokens%Tokens(i_Name),   &
               ColumnKey = 'Top Hat'                &
             )

    StartTime = Token2Time(                               &
                  C         = Tokens%Tokens(i_StartTime), &
                  BlockKey  = 'Sources',                  &
                  Item      = Tokens%Tokens(i_Name),      &
                  ColumnKey = 'Start Time'                &
                )

    StopTime  = Token2Time(                              &
                  C         = Tokens%Tokens(i_StopTime), &
                  BlockKey  = 'Sources',                 &
                  Item      = Tokens%Tokens(i_Name),     &
                  ColumnKey = 'Stop Time'                &
                )

    nSpecieses = 0
    Do i = 1, MaxSpecieses
      If (Tokens%Tokens(i_ZerothSourceStrength + i) == ' ') Exit
      nSpecieses = nSpecieses + 1
    End Do

    ! Particulate details
    If (                                                    &
      Tokens%Tokens(i_ParticleSizeDistribution) /= ' ' .or. &
      Tokens%Tokens(i_ParticleDiameter)         /= ' ' .or. &
      Tokens%Tokens(i_ParticleDensity)          /= ' '      &
    ) Then
      Particulate = .true.
    Else
      Particulate = .false.
      Call Message('Error: dust or sea-salt sources must be particulates', 3)
    End If

    If (Tokens%Tokens(i_ParticleSizeDistribution) == ' ') Then
      If (Tokens%Tokens(i_ParticleDiameter) /= ' ') Then
        Diameter = Token2Std(                                      &
                     C        = Tokens%Tokens(i_ParticleDiameter), &
                     BlockKey = 'Sources',                         &
                     Item     = Tokens%Tokens(i_Name),             &
                     ColumnKey = 'Particle Diameter'               &
                   )
      Else
        Call Message('Error: a diameter must be specified if a particle size distribution is not used', 3)
      End If
    Else
      If (Tokens%Tokens(i_ParticleDiameter) /= ' ') Then
        Call Message('Error: a diameter must not be specified if a particle size distribution is used', 3)
      End If
      Diameter = 0.0
    End If

    If (Tokens%Tokens(i_ParticleDensity) /= ' ') Then
      Density = Token2Std(                                     &
                  C        = Tokens%Tokens(i_ParticleDensity), &
                  BlockKey = 'Sources',                        &
                  Item     = Tokens%Tokens(i_Name),            &
                  ColumnKey = 'Particle Density'               &
                )
! else block no longer needed as density could be input via PSD. But just commented out for now. 
!    Else
!      Call Message('Error: must specify particle density for a dust or sea-salt source', 3)
    End If

    If (Tokens%Tokens(i_ParticleShape) /= ' ') Then
      ParticleShape = Token2Std(                                   &
                        C        = Tokens%Tokens(i_ParticleShape), &
                        BlockKey = 'Sources',                      &
                        Item     = Tokens%Tokens(i_Name),          &
                        ColumnKey = 'Particle Shape'               &
                      )
    End If

    UniformArea = Token2Log(                                  &
                    C         = Tokens%Tokens(i_UniformArea), &
                    BlockKey  = 'Sources',                    &
                    Item      = Tokens%Tokens(i_Name),        &
                    ColumnKey = 'Uniform Area?'               &
                  )

    NoReflect   = Token2Log(                                &
                    C         = Tokens%Tokens(i_NoReflect), &
                    BlockKey  = 'Sources',                  &
                    Item      = Tokens%Tokens(i_Name),      &
                    ColumnKey = 'No Reflect?'               &
                  )

  Else    ! Iterative plume rise model

    ! Check that H-Grid is not given for a generic source.
    If (Tokens%Tokens(i_HGrid) /= ' ') Then
      Call Message('Error in a source-block in the input file(s): a H-Grid must not be specified ' // &
                   'for a generic source', 3)
    End If
    If (Tokens%Tokens(i_ZGrid) /= ' ') Then
      Call Message('Error in a source-block in the input file(s): a Z-Grid must not be specified ' // &
                   'for a generic source', 3)
    End If

    ! Source location given by coords.
    If (                                           &
      Tokens%Tokens(i_X)              /= ' ' .and. &
      Tokens%Tokens(i_Y)              /= ' ' .and. &
      Tokens%Tokens(i_SetOfLocations) == ' ' .and. &
      Tokens%Tokens(i_Location      ) == ' '       &
    ) Then

      X(1) = Token2Std(                           &
               C         = Tokens%Tokens(i_X),    &
               BlockKey  = 'Sources',             &
               Item      = Tokens%Tokens(i_Name), &
               ColumnKey = 'X'                    &
             )

      X(2) = Token2Std(                           &
               C         = Tokens%Tokens(i_Y),    &
               BlockKey  = 'Sources',             &
               Item      = Tokens%Tokens(i_Name), &
               ColumnKey = 'Y'                    &
             )

      HCoordName = Tokens%Tokens(i_HCoord)

    ! Source location given by name.
    Else If (                        &
      Tokens%Tokens(i_X)              == ' ' .and. &
      Tokens%Tokens(i_Y)              == ' ' .and. &
      Tokens%Tokens(i_SetOfLocations) /= ' ' .and. &
      Tokens%Tokens(i_Location)       /= ' '       &
    ) Then

      ! $$ Need to prevent input of HCoordName via Tokens(i_HCoord) here, so user doesn't think its being used
      ! for dX, dY. If want option to not use the location coord system or metres aligned to that system,
      ! best to add a dHCoordName input.
      i = FindLocationsIndex(Tokens%Tokens(i_SetOfLocations), Locationses)
      j = FindLocationIndex(Tokens%Tokens(i_Location), Locationses%Locationses(i))
      X(1)       = Locationses%Locationses(i)%X(j)
      X(2)       = Locationses%Locationses(i)%Y(j)
      HCoordName = Locationses%Locationses(i)%HCoordNames(j)

    Else

      Call Message(                                                                   &
             'Error in a source-block in the input file(s): '                      // &
             'the source location must be specified by giving '                    // &
             'either "X", "Y" and "H-Coord" or "Set of locations" and "Location"',    &
             3                                                                        &
           )

    End If

    X(3) = Token2Std(                           &
             C         = Tokens%Tokens(i_Z),    &
             BlockKey  = 'Sources',             &
             Item      = Tokens%Tokens(i_Name), &
             ColumnKey = 'Z'                    &
           )

    dX(1) = Token2Std(                           &
              C         = Tokens%Tokens(i_dX),   &
              BlockKey  = 'Sources',             &
              Item      = Tokens%Tokens(i_Name), &
              ColumnKey = 'dX'                   &
            )

    dX(2) = Token2Std(                           &
              C         = Tokens%Tokens(i_dY),   &
              BlockKey  = 'Sources',             &
              Item      = Tokens%Tokens(i_Name), &
              ColumnKey = 'dY'                   &
            )

    dX(3) = Token2Std(                           &
              C         = Tokens%Tokens(i_dZ),   &
              BlockKey  = 'Sources',             &
              Item      = Tokens%Tokens(i_Name), &
              ColumnKey = 'dZ'                   &
            )

    Angle = Token2Std(                            &
              C         = Tokens%Tokens(i_Angle), &
              BlockKey  = 'Sources',              &
              Item      = Tokens%Tokens(i_Name),  &
              ColumnKey = 'Angle'                 &
            )

    dHMetres = Token2Log(                               &
                 C         = Tokens%Tokens(i_dHMetres), &
                 BlockKey  = 'Sources',                 &
                 Item      = Tokens%Tokens(i_Name),     &
                 ColumnKey = 'dH-Metres?'               &
               )

    dZMetres = Token2Log(                               &
                 C         = Tokens%Tokens(i_dZMetres), &
                 BlockKey  = 'Sources',                 &
                 Item      = Tokens%Tokens(i_Name),     &
                 ColumnKey = 'dZ-Metres?'               &
               )

    PlumeRise = Token2Log(                                &
                  C         = Tokens%Tokens(i_PlumeRise), &
                  BlockKey  = 'Sources',                  &
                  Item      = Tokens%Tokens(i_Name),      &
                  ColumnKey = 'Plume Rise?'               &
                )

    If (Tokens%Tokens(i_Temperature) /= ' ') Then
       Temperature = Token2Std(                                  &
                       C         = Tokens%Tokens(i_Temperature), &
                       BlockKey  = 'Sources',                    &
                       Item      = Tokens%Tokens(i_Name),        &
                       ColumnKey = 'Temperature'                 &
                     )
    Else
        Temperature = 1273.0  ! Default temperature for volcanic plumes
    End If

    If (Tokens%Tokens(i_FlowVelocity) /= ' ') Then
       FlowVelocity = Token2Std(                                   &
                        C         = Tokens%Tokens(i_FlowVelocity), &
                        BlockKey  = 'Sources',                     &
                        Item      = Tokens%Tokens(i_Name),         &
                        ColumnKey = 'Flow Velocity'                &
                      )
    Else
      FlowVelocity = 100.0    ! Default exit velocity for volcanic plumes
    End If

    If (Tokens%Tokens(i_DryGasFraction0) /= ' ') Then
       DryGasFraction0 = Token2Std(                                      &
                           C         = Tokens%Tokens(i_DryGasFraction0), &
                           BlockKey  = 'Sources',                        &
                           Item      = Tokens%Tokens(i_Name),            &
                           ColumnKey = 'Dry Gas Mass Fraction'           &
                          )
    Else
      DryGasFraction0 = 0.03  ! Default initial dry gas fraction 
    End If
    
    If (Tokens%Tokens(i_WaterVapourFraction0) /= ' ') Then
       WaterVapourFraction0 = Token2Std(                                           &
                                C         = Tokens%Tokens(i_WaterVapourFraction0), &
                                BlockKey  = 'Sources',                             &
                                Item      = Tokens%Tokens(i_Name),                 &
                                ColumnKey = 'Water Vapour Mass Fraction'           &
                              )
    Else
      WaterVapourFraction0 = 0.0  ! Default initial water vapour fraction
    End If

    If (Tokens%Tokens(i_PlumeRiseHeight) /= ' ') Then
       PlumeRiseHeight = Token2Std(                                      &
                           C         = Tokens%Tokens(i_PlumeRiseHeight), &
                           BlockKey  = 'Sources',                        &
                           Item      = Tokens%Tokens(i_Name),            &
                           ColumnKey = 'Plume Rise Height'               &
                         )
    Else
      Call Message(                                                                               &
             'FATAL ERROR: Plume rise height must be specificed for iterative plume rise scheme', &
             3                                                                                    &
           )
    End If

    If (Tokens%Tokens(i_SummitElevation) /= ' ') Then
       SummitElevation = Token2Std(                                      &
                           C         = Tokens%Tokens(i_SummitElevation), &
                           BlockKey  = 'Sources',                        &
                           Item      = Tokens%Tokens(i_Name),            &
                           ColumnKey = 'Summit Elevation'                &
                         )
    Else
      Call Message(                                                                              &
             'FATAL ERROR: Summit elevation must be specificed for iterative plume rise scheme', &
             3                                                                                   &
           )
    End If

    If (Tokens%Tokens(i_DistalFineAshFraction) /= ' ') Then
       DistalFineAshFraction = Token2Std(                                      &
                           C         = Tokens%Tokens(i_DistalFineAshFraction), &
                           BlockKey  = 'Sources',                              &
                           Item      = Tokens%Tokens(i_Name),                  &
                           ColumnKey = 'Distal Fine Ash Fraction'              &
                         )
    Else
      DistalFineAshFraction = 0.05 
    End If

    If (PlumeRise) Then

      ! $$ Default temperature for volcanic plume rise model may not be appropriate here
      If (                                            &
        (Tokens%Tokens(i_Shape) .CIEq. 'Cuboid') .or. &
        dX(1) /= dX(2) .or.                           &
        dX(3) /= 0.0                                  &
      ) Then
        Call Message(                                                      &
               'Error: plume rise only possible for Ellipsoid sources ' // &
               'with dX = dY and dZ = 0 at present',                       &
               3                                                           &
             )
      End If

      If (Tokens%Tokens(i_VolumeFlowRate) /= ' ') Then
        If (Tokens%Tokens(i_FlowVelocity) /= ' ') Then
          Call Message('Error: cannot specify VolumeFlowRate and FlowVelocity simultaneously', 3)
        Else
          VolumeFlowRate = Token2Std(                                     &
                             C         = Tokens%Tokens(i_VolumeFlowRate), &
                             BlockKey  = 'Sources',                       &
                             Item      = Tokens%Tokens(i_Name),           &
                             ColumnKey = 'Volume Flow Rate'               &
                           )
          FlowVelocity = VolumeFlowRate/(0.25*Pi*dX(1)*dX(2)) ! DX might not be metres $$
        End If                                                ! Better to convert later?
      Else 
        VolumeFlowRate = FlowVelocity*(0.25*Pi*dX(1)*dX(2))
      End If

    Else
      VolumeFlowRate = -999.0
    End If

    MaxTravelTime = Token2Time(                            &
                      C         = Tokens%Tokens(i_MaxAge), &
                      BlockKey  = 'Sources',               &
                      Item      = Tokens%Tokens(i_Name),   &
                      ColumnKey = 'Max Age',               &
                      Interval  = .true.                   &
                    )

    TopHat = Token2Log(                             &
               C         = Tokens%Tokens(i_TopHat), &
               BlockKey  = 'Sources',               &
               Item      = Tokens%Tokens(i_Name),   &
               ColumnKey = 'Top Hat'                &
             )

    StartTime = Token2Time(                               & 
                  C         = Tokens%Tokens(i_StartTime), &
                  BlockKey  = 'Sources',                  &
                  Item      = Tokens%Tokens(i_Name),      &
                  ColumnKey = 'Start Time'                &
                )
    StopTime  = Token2Time(                              &
                  C         = Tokens%Tokens(i_StopTime), &
                  BlockKey  = 'Sources',                 &
                  Item      = Tokens%Tokens(i_Name),     &
                  ColumnKey = 'Stop Time'                &
                )

    nSpecieses = 0
    Do i = 1, MaxSpecieses
      If (Tokens%Tokens(i_ZerothSourceStrength + i) == ' ') Exit
      nSpecieses = nSpecieses + 1
    End Do

    If (                                                    &
      Tokens%Tokens(i_ParticleSizeDistribution) /= ' ' .or. &
      Tokens%Tokens(i_ParticleDiameter)         /= ' ' .or. &
      Tokens%Tokens(i_ParticleDensity)          /= ' ' .or. &
      Tokens%Tokens(i_ParticleShape)            /= ' ' .or. &
      Tokens%Tokens(i_ShapeScheme)              /= ' '      &
    ) Then
      Particulate = .true.
    Else
      Particulate = .false.
    End If

    If (Particulate) Then                      
                  
        ! If no PSD then diameter must be defined.
        If (Tokens%Tokens(i_ParticleSizeDistribution) == ' ') Then
          Diameter = Token2Std(                                      &
                       C        = Tokens%Tokens(i_ParticleDiameter), &
                       BlockKey = 'Sources',                         &
                       Item     = Tokens%Tokens(i_Name),             &
                       ColumnKey = 'Particle Diameter'               &
                     )
        End If
        
        If (Tokens%Tokens(i_ParticleDensity) /= ' ') Then            
          Density = Token2Std(                                     &
                      C        = Tokens%Tokens(i_ParticleDensity), &
                      BlockKey = 'Sources',                        &
                      Item     = Tokens%Tokens(i_Name),            &
                      ColumnKey = 'Particle Density'               &
                    )
          If (Density <= 0.0) Then
            Call Message(                                   &
                  'FATAL ERROR:  Particle density is <= 0', &
                   3                                        &
                 )
          End If
        End If
        
        If (Tokens%Tokens(i_ParticleShape) /= ' ') Then
          ParticleShape = Token2Std(                                         &
                            C        = Tokens%Tokens(i_ParticleShape),       &
                            BlockKey = 'Sources',                            &
                            Item     = Tokens%Tokens(i_Name),                &
                            ColumnKey = 'Particle Shape'                     &
                          )
        End If
        
        If (Tokens%Tokens(i_ShapeScheme) .CIEq. 'Ganser') Then
          ShapeScheme = 'Ganser'
        Else If (Tokens%Tokens(i_ShapeScheme) .CIEq. 'WH') Then
          ShapeScheme = 'WH'
        Else If (Tokens%Tokens(i_ShapeScheme) .CIEq. 'White') Then
          ShapeScheme = 'White'
        Else If (Tokens%Tokens(i_ShapeScheme) .CIEq. ' ') Then
          ShapeScheme = ' '        
        Else
          Call Message(                                                           &
                 'FATAL ERROR: Error in a source-block in the input file(s): ' // & 
                 'the type of shape scheme is unknown. '                       // &
                 'The "Shape Scheme" should be left blank or '                 // &
                 'specified as either "White" or "Ganser" or "WH"',               &
                 3                                                                &
               )            
        End If
        
                    
        If (Tokens%Tokens(i_ShapeScheme) .CIEq. 'White') Then
          Call Message(                                                         &                              
                 'WARNING: '                                                 // &
                 'When using the "White" shape scheme the particle shape '   // &
                 'is necessarily equal to 1"',                                  &
                 1                                                              &
               )             
        End If
                         
        If (Tokens%Tokens(i_ParticleShape) /= ' ') Then
          If (ParticleShape > 1.0 .or. ParticleShape <= 0.0) Then
            Call Message(                                        &
                   'FATAL ERROR: Particle shape is > 1 or <= 0', &
                   3                                             &
                 )
          End If
        End If 

    End If

    UniformArea = Token2Log(                                  &
                    C         = Tokens%Tokens(i_UniformArea), &
                    BlockKey  = 'Sources',                    &
                    Item      = Tokens%Tokens(i_Name),        &
                    ColumnKey = 'Uniform Area?'               &
                  )
    NoReflect   = Token2Log(                                &
                    C         = Tokens%Tokens(i_NoReflect), &
                    BlockKey  = 'Sources',                  &
                    Item      = Tokens%Tokens(i_Name),      &
                    ColumnKey = 'No Reflect?'               &
                  )

  End If

  EulerianTime = Time2ShortTime(                                          &
                   Token2Time(                                            &
                     C         = Tokens%Tokens(i_LagrangianEulerianTime), &
                     BlockKey  = 'Sources',                               &
                     Item      = Tokens%Tokens(i_Name),                   &
                     ColumnKey = 'Lagrangian-Eulerian Time',              &
                     Interval  = .true.                                   &
                   )                                                      &
                 )

  Call AddSource(                                                                                           &
         InitSource(                                                                                        &
           Name                    = Tokens%Tokens(i_Name),                                                    &
           SourceGroupNames        = SourceGroupNames(1:nSourceGroups),                                        &
           SourceType              = SourceType,                                                               &
           HGridName               = Tokens%Tokens(i_HGrid),                                                   &
           ZGridName               = Tokens%Tokens(i_ZGrid),                                                   &
           HCoordName              = HCoordName,                                                               &
           ZCoordName              = Tokens%Tokens(i_ZCoord),                                                  &
           dHCoordName             = HCoordName,                                                               &
           dZCoordName             = Tokens%Tokens(i_ZCoord),                                                  &
           dHMetres                = dHMetres,                                                                 &
           dZMetres                = dZMetres,                                                                 &
           X                       = X(:),                                                                     &
           Shape                   = Tokens%Tokens(i_Shape),                                                   &
           TopHat                  = TopHat,                                                                   &
           UniformArea             = UniformArea,                                                              &
           NoReflect               = NoReflect,                                                                &
           dX                      = dX(:),                                                                    &
           Angle                   = Angle,                                                                    &
           StartTime               = StartTime,                                                                &
           StopTime                = StopTime,                                                                 &
           MaxTravelTime           = MaxTravelTime,                                                            &
           CharNParticles          = Tokens%Tokens(i_NParticles),                                              &
           CharSourceStrengths     = Tokens%Tokens(i_FirstSourceStrength:i_ZerothSourceStrength + nSpecieses), &
           PlumeRise               = PlumeRise,                                                                &
           VolumeFlowRate          = VolumeFlowRate,                                                           &
           FlowVelocity            = FlowVelocity,                                                             &
           Temperature             = Temperature,                                                              &
           Particulate             = Particulate,                                                              &
           Diameter                = Diameter,                                                                 &
           Density                 = Density,                                                                  &
           ParticleShape           = ParticleShape,                                                            &
           DensityPresent          = DensityPresent,                                                           &
           ParticleShapePresent    = ParticleShapePresent,                                                     &
           ShapeScheme             = ShapeScheme,                                                              &
           SizeDistName            = Tokens%Tokens(i_ParticleSizeDistribution),                                &
           TimeDepName             = Tokens%Tokens(i_TimeDependency),                                          &
           FixedMet                = FixedMet,                                                                 &
           EulerianTime            = EulerianTime,                                                             &                 
           ZGridIterativePlumeRise = Tokens%Tokens(i_ZGridIterativePlumeRise),                                 &
           DryGasFraction0         = DryGasFraction0,                                                          &
           WaterVapourFraction0    = WaterVapourFraction0,                                                     &
           PlumeRiseHeight         = PlumeRiseHeight,                                                          &
           SummitElevation         = SummitElevation,                                                          &
           DistalFineAshFraction   = DistalFineAshFraction                                                     &      
         ),                                                                                                    &
         Sources                                                                                               &
       )

End Subroutine Tokens2Sources

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine TimeDepsInputNames(HFBlockForms)
!.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                              &
         BlockKey             = 'Source Time Dependency:', &
         NamedBlock           = .true.,                    &
         MatchMultipleLines   = .false.,                   &
         nHeaderLines         = 1,                         &
         ColumnKeys           = Reshape(                   &
                                  (/                       &
                                    'From 1',              & ! 1
                                    'From 2',              & ! 2
                                    'From 3',              & ! 3
                                    'From 4',              & ! 4
                                    'From 5',              & ! 5
                                    'To 1  ',              & ! 6
                                    'To 2  ',              & ! 7
                                    'To 3  ',              & ! 8
                                    'To 4  ',              & ! 9
                                    'To 5  ',              & ! 10
                                    'Factor'               & ! 11
                                  /),                      &
                                  (/ 11, 1 /)              &
                                ),                         &
         Defaults             = (/                         &
                                  ' ',                & ! 1
                                  ' ',                & ! 2
                                  ' ',                & ! 3
                                  ' ',                & ! 4
                                  ' ',                & ! 5
                                  ' ',                & ! 6
                                  ' ',                & ! 7
                                  ' ',                & ! 8
                                  ' ',                & ! 9
                                  ' ',                & ! 10
                                  ' '                 & ! 11
                                /),                        &
         TwoD                 = .true.,                    &
         UnrecognisedMessages = .true.,                    &
         DisjointColSpecs     = .true.,                    &
         HFBlockForms         = HFBlockForms               &
       )

  If (MaxTimeDepPairs /= 5) Then
    Call Message(                                          &
           'UNEXPECTED FATAL ERROR in TimeDepsInputNames', &
           4                                               &
         )
  End If

End Subroutine TimeDepsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2TimeDeps(Tokens, TimeDeps)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),   Intent(In)    :: Tokens
  Type(TimeDeps_), Intent(InOut) :: TimeDeps
  ! Locals:
  Real(Std)           :: Factors(MaxTimeDepFactors)
  Integer             :: i, j, nPairs(MaxTimeDepFactors)
  Type(WildTimePair_) :: WildTimePairs(MaxTimeDepPairs, MaxTimeDepFactors)
  Type(WildTime_)     :: FromTime, ToTime

  Do i = 1, Tokens%nLines
    Factors(i) = Token2Std(                                               &
                   C         = Tokens%Tokens2d(2*MaxTimeDepPairs + 1, i), &
                   BlockKey  = 'Source Time Dependency',                  &
                   Item      = Tokens%BlockName,                          &
                   ColumnKey = 'Factor',                                  &
                   i         = i                                          &
                 )
    nPairs(i) = 0
    Do j = 1, MaxTimeDepPairs
      If (                                                   &
        Tokens%Tokens2d(j,                   i) == ' ' .and. &
        Tokens%Tokens2d(j + MaxTimeDepPairs, i) == ' '       &
      ) Exit
      nPairs(i) = j
      FromTime = Char2WildTime(Tokens%Tokens2d(j,                   i)) ! $$  use Errorcodes? Create Token2WildCardTime routine?
      ToTime   = Char2WildTime(Tokens%Tokens2d(j + MaxTimeDepPairs, i)) ! $$  use Errorcodes? Create Token2WildCardTime routine?
      WildTimePairs(j, i) = InitWildTimePair(FromTime, ToTime) ! $$ use Errorcodes? Create Token2WildCardTimePair routine?
    End Do
  End Do

  Call AddTimeDep(                                                            &
         InitTimeDep(                                                         &
           Name           = Tokens%BlockName,                                 &
           Factors        = Factors(1:Tokens%nLines),                         &
           nWildTimePairs = nPairs(1:Tokens%nLines),                          &
           WildTimePairs  = WildTimePairs(1:MaxTimeDepPairs, 1:Tokens%nLines) &
         ),                                                                   &
         TimeDeps                                                             &
       )

End Subroutine Tokens2TimeDeps

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine FieldReqsInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  ! $$ Output makes some assumptions about future additions to column keys.
  ! Check messages in output.F90 agree.

  Call InitAndAddHFBlockForm(                                     &
         BlockKey             = 'Output Requirements - Fields:',  &
         NamedBlock           = .false.,                          &
         MatchMultipleLines   = .false.,                          &
         nHeaderLines         = 1,                                &
         ColumnKeys           = Reshape(                          &
                                  (/                              &
                                    'Quantity                  ', & ! 1
                                    'Species                   ', & ! 2
                                    'Source                    ', & ! 3
                                    'Source Group              ', & ! 4
                                    'H-Grid                    ', & ! 5
                                    'Z-Grid                    ', & ! 6
                                    'T-Grid                    ', & ! 7
                                    'S-Grid                    ', & ! 8
                                    'H-Coord                   ', & ! 9
                                    'Z-Coord                   ', & ! 10
                                    'BL Average                ', & ! 11
                                    'T Av Or Int               ', & ! 12
                                    'Sync?                     ', & ! 13
                                    'Graph?xxxxx               ', & ! 14 $$ use output route
                                    'Screen?xxxxx              ', & ! 15 $$
                                    'Disk?xxxxx                ', & ! 16 $$
                                    'xxxxxxx                   ', & ! 17 $$ remove
                                    'Decay Deposition?         ', & ! 18
                                    'Across                    ', & ! 19
                                    'Separate File             ', & ! 20
                                    'X Scale                   ', & ! 21
                                    'Output Group              ', & ! 22
                                    'Av Time                   ', & ! 23 $$ alt name Av T
                                    '# Av Times                ', & ! 24
                                    'Output Format             ', & ! 25
                                    'Y Scale                   ', & ! 26
                                    'Percentiles               ', & ! 27
                                    'P Time                    ', & ! 28 $$ alt name P T
                                    'P dT                      ', & ! 29
                                    'Probabilities             ', & ! 30
                                    'Ensemble Av?              ', & ! 31
                                    'Ensemble P?               ', & ! 32
                                    'Name                      ', & ! 33
                                    'Output Route              ', & ! 34
                                    'Fluctuations?             ', & ! 35
                                    'Particle Size Distribution', & ! 36
                                    'Semi-Infinite Approx?     ', & ! 37 ! $$ right name?
                                    'Material Unit             '  & ! 38
                                  /),                             &
                                  (/ 38, 1 /)                     &
                                ),                                &
         Defaults             = (/                                &
                                  '    ',           & ! 1
                                  '    ',           & ! 2
                                  '    ',           & ! 3
                                  '    ',           & ! 4
                                  '    ',           & ! 5
                                  '    ',           & ! 6
                                  '    ',           & ! 7
                                  '    ',           & ! 8
                                  '    ',           & ! 9
                                  '    ',           & ! 10
                                  'No  ',                         & ! 11
                                  '    ',           & ! 12
                                  '    ',           & ! 13
                                  '    ',           & ! 14
                                  '    ',           & ! 15
                                  '    ',           & ! 16
                                  '    ',           & ! 17
                                  '    ',           & ! 18
                                  '    ',           & ! 19
                                  '    ',           & ! 20
                                  '1.0 ',                         & ! 21
                                  '    ',           & ! 22
                                  '    ',           & ! 23
                                  '    ',           & ! 24
                                  '    ',           & ! 25
                                  '1.0 ',                         & ! 26
                                  '    ',           & ! 27
                                  '    ',           & ! 28
                                  '    ',           & ! 29
                                  '    ',           & ! 30
                                  'No  ',                         & ! 31
                                  'No  ',                         & ! 32
                                  '    ',           & ! 33
                                  '    ',                         & ! 34
                                  'No  ',                         & ! 35
                                  '    ',                         & ! 36
                                  'No  ',                         & ! 37 ! $$
                                  '    '                          & ! 38
                                /),                               &
         TwoD                 = .false.,                          &
         UnrecognisedMessages = .true.,                           &
         DisjointColSpecs     = .true.,                           &
         HFBlockForms         = HFBlockForms                      &
       )

End Subroutine FieldReqsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2FieldReqs(Tokens, Grids, Reqs, MaterialUnits)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),        Intent(In)    :: Tokens
  Type(Grids_),         Intent(InOut) :: Grids
  Type(Reqs_),          Intent(InOut) :: Reqs
  Type(MaterialUnits_), Intent(InOut) :: MaterialUnits
  ! Local parameters:
  Integer, Parameter :: i_Quantity                 = 1
  Integer, Parameter :: i_Species                  = 2
  Integer, Parameter :: i_Source                   = 3
  Integer, Parameter :: i_SourceGroup              = 4
  Integer, Parameter :: i_HGrid                    = 5
  Integer, Parameter :: i_ZGrid                    = 6
  Integer, Parameter :: i_TGrid                    = 7
  Integer, Parameter :: i_SGrid                    = 8
  Integer, Parameter :: i_HCoord                   = 9
  Integer, Parameter :: i_ZCoord                   = 10
  Integer, Parameter :: i_BLAverage                = 11
  Integer, Parameter :: i_TAvOrInt                 = 12
  Integer, Parameter :: i_Sync                     = 13
  Integer, Parameter :: i_Graphxxxxx               = 14
  Integer, Parameter :: i_Screenxxxxx              = 15
  Integer, Parameter :: i_Diskxxxxx                = 16
  Integer, Parameter :: i_xxxxxxx                  = 17
  Integer, Parameter :: i_DecayDeposition          = 18
  Integer, Parameter :: i_Across                   = 19
  Integer, Parameter :: i_SeparateFile             = 20
  Integer, Parameter :: i_XScale                   = 21
  Integer, Parameter :: i_OutputGroup              = 22
  Integer, Parameter :: i_AvTime                   = 23
  Integer, Parameter :: i_nAvTimes                 = 24
  Integer, Parameter :: i_OutputFormat             = 25
  Integer, Parameter :: i_YScale                   = 26
  Integer, Parameter :: i_Percentiles              = 27
  Integer, Parameter :: i_PTime                    = 28
  Integer, Parameter :: i_PdT                      = 29
  Integer, Parameter :: i_Probabilities            = 30
  Integer, Parameter :: i_EnsembleAv               = 31
  Integer, Parameter :: i_EnsembleP                = 32
  Integer, Parameter :: i_Name                     = 33
  Integer, Parameter :: i_OutputRoute              = 34
  Integer, Parameter :: i_Fluctuations             = 35
  Integer, Parameter :: i_ParticleSizeDistribution = 36
  Integer, Parameter :: i_SemiInfiniteApprox       = 37
  Integer, Parameter :: i_MaterialUnit             = 38
  ! Locals:
  Integer         :: iQuantity
  Logical         :: DecayDep
  Logical         :: AvBL
  Logical         :: AvEnsemble
  Integer         :: TAvInt
  Type(Time_)     :: AvTime
  Integer         :: nAvTimes
  Integer         :: HAvInt
  Integer         :: ZAvInt
  Integer         :: PCode
  Logical         :: EnsembleP
  Type(Time_)     :: PTime
  Type(Time_)     :: PdT
  Integer         :: nPTimes
  Integer         :: nDValues
  Integer         :: nSpeciesNames
  Real(Std)       :: DValues(50)
  Character(MaxCharLength) :: SpeciesNames(50)
  Real(Std)       :: X,Y,Z,dX,dY,dZ
  Logical         :: Sync
  Logical         :: Graph
  Logical         :: Screen
  Logical         :: Disk
  Real(Std)       :: XScale
  Real(Std)       :: YScale
  Type(FieldReq_) :: FieldReq
  Character(MaxTokenLength) :: Name
  Character(MaxCharLength)  :: NamePart2                            ! $$ temp fix up to allow
                                                                    ! multiple species to be given
                                                                    ! in input file
  Character(MaxCharLength) :: ProcessNames(MaxProcessesPerFieldReq) ! $$ need to trap attempt to
                                                                    ! create too many processes?
  Integer                  :: nProcesses
  Character(MaxCharLength) :: DGridName
  Integer                  :: i
  Logical                  :: IntrinsicHAv
  Logical                  :: IntrinsicZAv
  Logical                  :: Fluctuations
  Logical                  :: SemiInfApprox

!  If (Tokens%Tokens(i_Name) == ' ') Then ! $$ remove when names compulsory
!    Name = 'Unnamed requirement ' // Int2Char(Reqs%nFieldReqs + 1)
!  Else
    Name = Tokens%Tokens(i_Name)
!  End If

  If (                                                            &
    (Tokens%Tokens(i_Quantity) .CIEq. 'Dry Deposition Rate') .or. &
    (Tokens%Tokens(i_Quantity) .CIEq. 'Wet Deposition Rate') .or. &
    (Tokens%Tokens(i_Quantity) .CIEq. 'Deposition Rate'    )      &
  ) Then
    If (Tokens%Tokens(i_DecayDeposition) /= ' ') Then
      DecayDep = Token2Log(                                      &
                   C         = Tokens%Tokens(i_DecayDeposition), &
                   BlockKey  = 'Output Requirements - Fields',   &
                   Item      = Tokens%Tokens(i_Name),            &
                   ColumnKey = 'Decay deposition?'               &
                 )
    Else
      DecayDep = .true.
    End If
  Else
    If (Tokens%Tokens(i_DecayDeposition) /= ' ') Then
      Call Message("Error: Can't specify 'Decay Deposition?' for a non-deposition field", 3)
    Else
      DecayDep = .false.
    End If
  End If

  Sync = Token2Log(                                    &
           C         = Tokens%Tokens(i_Sync),          &
           BlockKey  = 'Output Requirements - Fields', &
           Item      = Tokens%Tokens(i_Name),          &
           ColumnKey = 'Sync?'                         &
         )

  XScale = Token2Std(                                    &
             C         = Tokens%Tokens(i_XScale),        &
             BlockKey  = 'Output Requirements - Fields', &
             Item      = Tokens%Tokens(i_Name),          &
             ColumnKey = 'X Scale'                       &
           )
  YScale = Token2Std(                                    &
             C         = Tokens%Tokens(i_YScale),        &
             BlockKey  = 'Output Requirements - Fields', &
             Item      = Tokens%Tokens(i_Name),          &
             ColumnKey = 'Y Scale'                       &
           )

  ! Averaging.
  ! $$ input options to support: no T av
  !                              2 of nT, Av T and Av dT
  !                              nT - av based on tgrid
  !                              intrinsic TAv - intrinsic av based on tgrid
  !                              nT + intrinsic TAv

  AvEnsemble = Token2Log(                                    &
                 C         = Tokens%Tokens(i_EnsembleAv),    &
                 BlockKey  = 'Output Requirements - Fields', &
                 Item      = Tokens%Tokens(i_Name),          &
                 ColumnKey = 'Ensemble Av?'                  &
               )
  If (Tokens%Tokens(i_TAvOrInt) .CIEq. 'Av') Then
    TAvInt  = P_Av
  Else If (Tokens%Tokens(i_TAvOrInt) .CIEq. 'Int') Then
    TAvInt  = P_Int
  Else If ((Tokens%Tokens(i_TAvOrInt) .CIEq. ' ') .or. (Tokens%Tokens(i_TAvOrInt) .CIEq. 'No')) Then
    TAvInt  = 0
  Else
    Call Message('FATAL ERROR: "T Av Or Int" must be "Av", "Int", "No" or blank', 3)
  End If

  If ((Tokens%Tokens(i_TAvOrInt) .CIEq. ' ') .or. (Tokens%Tokens(i_TAvOrInt) .CIEq. 'No')) Then
    AvTime   = ZeroTime()
    nAvTimes = 1
  Else
    AvTime   = Token2Time(                                   &
                 C         = Tokens%Tokens(i_AvTime),        &
                 BlockKey  = 'Output Requirements - Fields', &
                 Item      = Tokens%Tokens(i_Name),          &
                 ColumnKey = 'Av Time',                      &
                 Interval  = .true.                          &
               )
    ! $$ Input Av dT not nAvTimes. Allow blank for no time processing.
    nAvTimes = Token2Int(                                    &
                 C         = Tokens%Tokens(i_nAvTimes),      &
                 BlockKey  = 'Output Requirements - Fields', &
                 Item      = Tokens%Tokens(i_Name),          &
                 ColumnKey = '# Av Times'                    &
               )
  End If

  iQuantity = 0
  Do i = 1, Size(Quantities)
    If (Tokens%Tokens(i_Quantity) .CIEq. Quantities(i)) Then
      iQuantity = i
      Exit
    End If
  End Do
  If (iQuantity == 0) Then
    Call Message(                                                     &
           'FATAL ERROR in reading Output Requirements - Fields: ' // &
           'Quantity "'                                            // &
           Trim(Tokens%Tokens(i_Quantity))                         // &
           '" not recognised',                                        &
           3                                                          &
         )
  End If

  AvBL = Token2Log(                                    &
           C         = Tokens%Tokens(i_BLAverage),     &
           BlockKey  = 'Output Requirements - Fields', &
           Item      = Tokens%Tokens(i_Name),          &
           ColumnKey = 'BL Average'                    &
         )
  ZAvInt = 0
  Z  = 0.0
  dZ = 0.0
  If (Scan(QInfo(iQuantity), 'Z') /= 0 .and. Tokens%Tokens(i_ZGrid) == ' ') Then
    If (AvBL) Then
      ZAvInt = P_Av
    Else
      ZAvInt = P_Int
      Z  = -1.0
      dZ = -1.0
    End If
  End If
  HAvInt = 0
  X  = 0.0
  Y  = 0.0
  dX = 0.0
  dY = 0.0
  If (Scan(QInfo(iQuantity), 'H') /= 0 .and. Tokens%Tokens(i_HGrid) == ' ') Then
    HAvInt = P_Int
    X  = -1.0
    Y  = -1.0
    dX = -1.0
    dY = -1.0
  End If

  ! Probabilities. $$ add support for moments,
  !                   named arrays, named d-grids, named processes.
  ! If processes input, this must include intermediate prob steps in calculating percentiles, or have
  ! option in inputting percent processes.

  If (Tokens%Tokens(i_Percentiles) == ' ' .and. Tokens%Tokens(i_Probabilities) /= ' ') Then
    PCode   = P_Prob
    Call ParseListStd(                                 &
           C         = Tokens%Tokens(i_Probabilities), &
           Delim     = ';',                            &
           BlockKey  = 'Output Requirements - Fields', &
           Item      = Tokens%Tokens(i_Name),          &
           ColumnKey = 'Probabilities',                &
           nValues   = nDValues,                       &
           Values    = DValues                         &
         )
  Else If (Tokens%Tokens(i_Percentiles) /= ' ' .and. Tokens%Tokens(i_Probabilities) == ' ') Then
    PCode   = P_Percent
    Call ParseListStd(                                 &
           C         = Tokens%Tokens(i_Percentiles),   &
           Delim     = ';',                            &
           BlockKey  = 'Output Requirements - Fields', &
           Item      = Tokens%Tokens(i_Name),          &
           ColumnKey = 'Percentiles',                  &
           nValues   = nDValues,                       &
           Values    = DValues                         &
         )
  Else If (Tokens%Tokens(i_Percentiles) == ' ' .and. Tokens%Tokens(i_Probabilities) == ' ') Then
    PCode = 0
  Else
    Call Message ('FATAL ERROR: At most one of "Percentile" and "Probabilities" can be used at once', 3)
  End If

  If (PCode == P_Prob .or. PCode == P_Percent) Then
    ! $$ check PTime, Pdt, EnsembleP present/absent consistently
    ! but allow PTime, Pdt both absent for no time
    ! processing (AvInt stuff too)
    EnsembleP = Token2Log(                                    &
                  C         = Tokens%Tokens(i_EnsembleP),     &
                  BlockKey  = 'Output Requirements - Fields', &
                  Item      = Tokens%Tokens(i_Name),          &
                  ColumnKey = 'Ensemble P?'                   &
                )
    If (Tokens%Tokens(i_PTime) /= ' ' .and. Tokens%Tokens(i_PdT) /= ' ') Then
      PTime = Token2Time(                                   &
                C         = Tokens%Tokens(i_PTime),         &
                BlockKey  = 'Output Requirements - Fields', &
                Item      = Tokens%Tokens(i_Name),          &
                ColumnKey = 'P Time',                       &
                Interval  = .true.                          &
              )
      PdT   = Token2Time(                                   &
                C         = Tokens%Tokens(i_PdT),           &
                BlockKey  = 'Output Requirements - Fields', &
                Item      = Tokens%Tokens(i_Name),          &
                ColumnKey = 'P dT',                         &
                Interval  = .true.                          &
              )
      If (PTime == ZeroTime() .and. PdT == ZeroTime()) Then ! $$ not all possibilities covered yet.
        nPTimes = 1
      Else If (PTime /= ZeroTime() .and. PdT /= ZeroTime()) Then
        nPTimes = PTime/PdT
      End If
    Else
      PTime   = ZeroTime()
      PdT     = ZeroTime()
      nPTimes = 1
    End If
    If (PCode == P_Percent .and. All(DValues(1:nDValues) == 100.0 .or. DValues(1:nDValues) == 0.0)) Then
      PCode = P_MaxMin
    End If
  End If

  Fluctuations = Token2Log(                                    &
                   C         = Tokens%Tokens(i_Fluctuations),  &
                   BlockKey  = 'Output Requirements - Fields', &
                   Item      = Tokens%Tokens(i_Name),          &
                   ColumnKey = 'Fluctuations?'                 &
                 )

  SemiInfApprox = Token2Log(                                         &
                    C         = Tokens%Tokens(i_SemiInfiniteApprox), &
                    BlockKey  = 'Output Requirements - Fields',      &
                    Item      = Tokens%Tokens(i_Name),               &
                    ColumnKey = 'Semi-infinite approx?'              &
                  )

  IntrinsicHAv = HAvInt /= 0
  IntrinsicZAv = ZAvInt /= 0

  ! Add required D-grids and processing steps.

  nProcesses = 0

  If (TAvInt /= 0 .or. HAvInt /= 0 .or. ZAvInt /= 0) Then
    If (TAvInt == 0) TAvInt = P_Av
    If (HAvInt == 0) HAvInt = P_Av
    If (ZAvInt == 0) ZAvInt = P_Av
    nProcesses = nProcesses + 1
    Call AddProcess(                                     &
           InitProcess(                                  &
             Name      = ' ',                            &
             Code      = P_AvInt,                        &
             Ensemble  = AvEnsemble,                     &
             TAvInt    = TAvInt,                         &
             nT        = nAvTimes,                       &
             T         = AvTime,                         &
             dT        = AvTime/nAvTimes,                &
             UseTGrid  = .false.,                        &
             T0        = ReferenceTime(),                &
             HAvInt    = HAvInt,                         &
             nX        = 1,                              &
             nY        = 1,                              &
             X         = X,                              &
             Y         = Y,                              &
             dX        = dX,                             &
             dY        = dY,                             &
             UseHGrid  = .false.,                        &
             ZAvInt    = ZAvInt,                         &
             nZ        = 1,                              &
             Z         = Z,                              &
             dZ        = dZ,                             &
             UseZGrid  = .false.,                        &
             BL        = AvBL,                           &
             DGridName = ' ',                            &
             BlockKey  = 'Output Requirements - Fields', &
             Item      = ' '                             &
           ),                                            &
           Reqs,                                         &
           Name = ProcessNames(nProcesses)               &
         )
  End If

  If (PCode == P_Prob) Then

    nProcesses = nProcesses + 1
    DGridName = Trim(Int2Char(PCode)) // Trim(FileNameTime(PTime)) // Trim(FileNameTime(PdT)) &
             // Trim(Std2Char(DValues(1)))
    Call AddDGrid(                      &
           InitDGrid(                   &
             Name = DGridName,          &
             nD   = nDValues,           &
             D    = DValues(1:nDValues) &
           ),                           &
           Grids                        &
         )
    Call AddProcess(                                     &
           InitProcess(                                  &
             Name      = ' ',                            &
             Code      = PCode,                          &
             Ensemble  = EnsembleP,                      &
             TAvInt    = 0,                              &
             nT        = nPTimes,                        &
             T         = PTime,                          &
             dT        = PdT,                            &
             UseTGrid  = .false.,                        &
             T0        = ReferenceTime(),                &
             HAvInt    = 0,                              &
             nX        = 1,                              &
             nY        = 1,                              &
             X         = 0.0,                            &
             Y         = 0.0,                            &
             dX        = 0.0,                            &
             dY        = 0.0,                            &
             UseHGrid  = .false.,                        &
             ZAvInt    = 0,                              &
             nZ        = 1,                              &
             Z         = 0.0,                            &
             dZ        = 0.0,                            &
             UseZGrid  = .false.,                        &
             BL        = .false.,                        &
             DGridName = DGridName,                      &
             BlockKey  = 'Output Requirements - Fields', &
             Item      = ' '                             &
           ),                                            &
           Reqs,                                         &
           Name = ProcessNames(nProcesses)               &
         )

  End If

  If (PCode == P_Percent) Then

    nProcesses = nProcesses + 1
    ! Note 'P' needed to keep name unique - there might be a prob request that generates the same name.
    ! Eventually rationalise with AddDGrid returning unique name. $$
    DGridName = 'P' // Trim(Int2Char(P_Prob)) // Trim(FileNameTime(PTime)) // Trim(FileNameTime(PdT)) &
             // Trim(Std2Char(DValues(1)))
    Call AddDGrid(                   &
           InitDGrid(                &
             Name      = DGridName,  &
             nD        = 20,         &
             dD        = 10.0**0.25, &
             D0        = 1.0,        &
             Floating  = .true.,     &
             Geometric = .true.,     &
             Zero      = .true.      &
           ),                        &
           Grids                     &
         )
    Call AddProcess(                                     &
           InitProcess(                                  &
             Name      = ' ',                            &
             Code      = P_Prob,                         &
             Ensemble  = EnsembleP,                      &
             TAvInt    = 0,                              &
             nT        = nPTimes,                        &
             T         = PTime,                          &
             dT        = PdT,                            &
             UseTGrid  = .false.,                        &
             T0        = ReferenceTime(),                &
             HAvInt    = 0,                              &
             nX        = 1,                              &
             nY        = 1,                              &
             X         = 0.0,                            &
             Y         = 0.0,                            &
             dX        = 0.0,                            &
             dY        = 0.0,                            &
             UseHGrid  = .false.,                        &
             ZAvInt    = 0,                              &
             nZ        = 1,                              &
             Z         = 0.0,                            &
             dZ        = 0.0,                            &
             UseZGrid  = .false.,                        &
             BL        = .false.,                        &
             DGridName = DGridName,                      &
             BlockKey  = 'Output Requirements - Fields', &
             Item      = ' '                             &
           ),                                            &
           Reqs,                                         &
           Name = ProcessNames(nProcesses)               &
         )

    nProcesses = nProcesses + 1
    DGridName = Trim(Int2Char(PCode)) // Trim(FileNameTime(PTime)) // Trim(FileNameTime(PdT)) &
              // Trim(Std2Char(DValues(1)))
    Call AddDGrid(                      &
           InitDGrid(                   &
             Name = DGridName,          &
             nD   = nDValues,           &
             D    = DValues(1:nDValues) &
           ),                           &
           Grids                        &
         )
    Call AddProcess(                                     &
           InitProcess(                                  &
             Name      = ' ',                            &
             Code      = PCode,                          &
             Ensemble  = .false.,                        &
             TAvInt    = 0,                              &
             nT        = 1,                              &
             T         = ZeroTime(),                     &
             dT        = ZeroTime(),                     &
             UseTGrid  = .false.,                        &
             T0        = ReferenceTime(),                &
             HAvInt    = 0,                              &
             nX        = 1,                              &
             nY        = 1,                              &
             X         = 0.0,                            &
             Y         = 0.0,                            &
             dX        = 0.0,                            &
             dY        = 0.0,                            &
             UseHGrid  = .false.,                        &
             ZAvInt    = 0,                              &
             nZ        = 1,                              &
             Z         = 0.0,                            &
             dZ        = 0.0,                            &
             UseZGrid  = .false.,                        &
             BL        = .false.,                        &
             DGridName = DGridName,                      &
             BlockKey  = 'Output Requirements - Fields', &
             Item      = ' '                             &
           ),                                            &
           Reqs,                                         &
           Name = ProcessNames(nProcesses)               &
         )

  End If

  If (PCode == P_MaxMin) Then

    nProcesses = nProcesses + 1
    DGridName = Trim(Int2Char(PCode)) // Trim(FileNameTime(PTime)) // Trim(FileNameTime(PdT)) &
              // Trim(Std2Char(DValues(1)))
    Call AddDGrid(                      &
           InitDGrid(                   &
             Name = DGridName,          &
             nD   = nDValues,           &
             D    = DValues(1:nDValues) &
           ),                           &
           Grids                        &
         )
    Call AddProcess(                                     &
           InitProcess(                                  &
             Name      = ' ',                            &
             Code      = PCode,                          &
             Ensemble  = EnsembleP,                      &
             TAvInt    = 0,                              &
             nT        = nPTimes,                        &
             T         = PTime,                          &
             dT        = PdT,                            &
             UseTGrid  = .false.,                        &
             T0        = ReferenceTime(),                &
             HAvInt    = 0,                              &
             nX        = 1,                              &
             nY        = 1,                              &
             X         = 0.0,                            &
             Y         = 0.0,                            &
             dX        = 0.0,                            &
             dY        = 0.0,                            &
             UseHGrid  = .false.,                        &
             ZAvInt    = 0,                              &
             nZ        = 1,                              &
             Z         = 0.0,                            &
             dZ        = 0.0,                            &
             UseZGrid  = .false.,                        &
             BL        = .false.,                        &
             DGridName = DGridName,                      &
             BlockKey  = 'Output Requirements - Fields', &
             Item      = ' '                             & ! $$ set value
           ),                                            &
           Reqs,                                         &
           Name = ProcessNames(nProcesses)               &
         )

  End If

  Call TokenLengthTest(                              & !$$ this is a check for blank strings, such as might 
         C         = Tokens%Tokens(i_OutputRoute),   & !   occur if user was using 'Disk?' instead of 'Output Route'.
         Length    = MaxCharLength,                  & !   Disk? has been removed as an option.
         Zero      = .true.,                         & !   Generally we don't check for blanks in Input.F90
         BlockKey  = 'Output Requirements - Fields', & !   but perhaps should review this.
         Item      = Tokens%Tokens(i_Name),          &
         ColumnKey = 'Output Route'                  &
       ) 

  ! Add an individual field request for each species in the input list.
  ! $$ this is a temporary fix up to allow multiple species to be specified for a field in the input file.
  ! $$ Ideally these species should be stored and processed within a single field request.
  If (Tokens%Tokens(i_Species) /= ' ') Then
    ! Species dependent field request

    Call ParseListChar(                                &
           C         = Tokens%Tokens(i_Species),       &
           Delim     = ';',                            &
           BlockKey  = 'Output Requirements - Fields', &
           Item      = Tokens%Tokens(i_Name),          &
           ColumnKey = 'Species',                      &
           nValues   = nSpeciesNames,                  &
           Values    = SpeciesNames                    &
         )

    Do i = 1, nSpeciesNames

      ! If multiple species are requested then append the species index to the name of the field request
      ! (this ensures that each individual request has a unique name).
      If (nSpeciesNames > 1) Then
        NamePart2 = '_' // Int2Char(i)
      Else
        NamePart2 = ''
      End If

      FieldReq = InitFieldReq(                                                   &
                   Name             = Trim(Name) // NamePart2,                   &
                   Quantity         = Tokens%Tokens(i_Quantity),                 &
                   SpeciesName      = SpeciesNames(i),                           &
                   SourceName       = Tokens%Tokens(i_Source),                   &
                   SourceGroupName  = Tokens%Tokens(i_SourceGroup),              &
                   SizeDistName     = Tokens%Tokens(i_ParticleSizeDistribution), &
                   DecayDep         = DecayDep,                                  &
                   SemiInfApprox    = SemiInfApprox,                             &
                   TGridName        = Tokens%Tokens(i_TGrid),                    &
                   SGridName        = Tokens%Tokens(i_SGrid),                    &
                   HGridName        = Tokens%Tokens(i_HGrid),                    &
                   ZGridName        = Tokens%Tokens(i_ZGrid),                    &
                   HCoordName       = Tokens%Tokens(i_HCoord),                   &
                   ZCoordName       = Tokens%Tokens(i_ZCoord),                   &
                   ProcessNames     = ProcessNames(1:nProcesses),                &
                   IntrinsicTAv     = .false.,                                   &
                   IntrinsicHAv     = IntrinsicHAv,                              &
                   IntrinsicZAv     = IntrinsicZAv,                              &
                   Fluctuations     = Fluctuations,                              &
                   FlAvT            = ZeroShortTime(), & ! $$ Set UseFlAvT etc based on presence in
                   FlAvX            = 0.0,             & ! $$ input file. FlAvT = 0 etc if Use false
                   FlAvY            = 0.0,             & ! $$
                   FlAvZ            = 0.0,             & ! $$
                   UseFlAvT         = .false.,         & ! $$
                   UseFlAvH         = .false.,         & ! $$
                   UseFlAvZ         = .false.,         & ! $$
                   Sync             = Sync,                                      &
                   BlockKey         = 'Output Requirements - Fields',            &
                   OutputRoute      = Tokens%Tokens(i_OutputRoute),              &
                   OutputGroup      = Tokens%Tokens(i_OutputGroup),              &
                   Across           = Tokens%Tokens(i_Across),                   &
                   SeparateFile     = Tokens%Tokens(i_SeparateFile),             &
                   OutputFormat     = Tokens%Tokens(i_OutputFormat),             &
                   XScale           = XScale,                                    &
                   YScale           = YScale,                                    &
                   MaterialUnits    = MaterialUnits,                             & 
                   MaterialUnitName = Tokens%Tokens(i_MaterialUnit)              &
                 )

      Call AddFieldReq(FieldReq, Reqs)

    End Do

  Else
    ! Other field request (no species dependency)

    FieldReq = InitFieldReq(                                                   &
                 Name             = Name,                                      &
                 Quantity         = Tokens%Tokens(i_Quantity),                 &
                 SpeciesName      = Tokens%Tokens(i_Species),                  &
                 SourceName       = Tokens%Tokens(i_Source),                   &
                 SourceGroupName  = Tokens%Tokens(i_SourceGroup),              &
                 SizeDistName     = Tokens%Tokens(i_ParticleSizeDistribution), &
                 DecayDep         = DecayDep,                                  &
                 SemiInfApprox    = SemiInfApprox,                             &
                 TGridName        = Tokens%Tokens(i_TGrid),                    &
                 SGridName        = Tokens%Tokens(i_SGrid),                    &
                 HGridName        = Tokens%Tokens(i_HGrid),                    &
                 ZGridName        = Tokens%Tokens(i_ZGrid),                    &
                 HCoordName       = Tokens%Tokens(i_HCoord),                   &
                 ZCoordName       = Tokens%Tokens(i_ZCoord),                   &
                 ProcessNames     = ProcessNames(1:nProcesses),                &
                 IntrinsicTAv     = .false.,                                   &
                 IntrinsicHAv     = IntrinsicHAv,                              &
                 IntrinsicZAv     = IntrinsicZAv,                              &
                 Fluctuations     = Fluctuations,                              &
                 FlAvT            = ZeroShortTime(), & ! $$ Set UseFlAvT etc based on presence in
                 FlAvX            = 0.0,             & ! $$ input file. FlAvT = 0 etc if Use false
                 FlAvY            = 0.0,             & ! $$
                 FlAvZ            = 0.0,             & ! $$
                 UseFlAvT         = .false.,         & ! $$
                 UseFlAvH         = .false.,         & ! $$
                 UseFlAvZ         = .false.,         & ! $$
                 Sync             = Sync,                                      &
                 MaterialUnits    = MaterialUnits,                             & 
                 MaterialUnitName = '',                                        &
                 BlockKey         = 'Output Requirements - Fields',            &
                 OutputRoute      = Tokens%Tokens(i_OutputRoute),              &
                 OutputGroup      = Tokens%Tokens(i_OutputGroup),              &
                 Across           = Tokens%Tokens(i_Across),                   &
                 SeparateFile     = Tokens%Tokens(i_SeparateFile),             &
                 OutputFormat     = Tokens%Tokens(i_OutputFormat),             &
                 XScale           = XScale,                                    &
                 YScale           = YScale                                     &
               )

    Call AddFieldReq(FieldReq, Reqs)

  End If

End Subroutine Tokens2FieldReqs

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine PdfReqsInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  ! $$ Currently PdfMode not input - set to 1. Also CalcType not input.

  Call InitAndAddHFBlockForm(                                  &
         BlockKey             = 'Output Requirements - Pdfs:', &
         NamedBlock           = .false.,                       &
         MatchMultipleLines   = .false.,                       &
         nHeaderLines         = 1,                             &
         ColumnKeys           = Reshape(                       &
                                  (/                           &
                                    'Name         ',           & ! 1
                                    'Species      ',           & ! 2
                                    'Source       ',           & ! 3
                                    'Source Group ',           & ! 4
                                    'Pdf Type     ',           & ! 5
                                    'H-Grid       ',           & ! 6
                                    'Z-Grid       ',           & ! 7
                                    'T-Grid       ',           & ! 8
                                    'T Av Or Int  ',           & ! 9
                                    'Average Time ',           & ! 10
                                    'Screen?      ',           & ! 11
                                    'Disk?        ',           & ! 12
                                    'Stat?        '            & ! 13
                                  /),                          &
                                  (/ 13, 1 /)                  &
                                ),                             &
         Defaults             = (/                             &
                                  '   ',             & ! 1
                                  '   ',             & ! 2
                                  '   ',             & ! 3
                                  '   ',             & ! 4
                                  '   ',             & ! 5
                                  '   ',             & ! 6
                                  '   ',             & ! 7
                                  '   ',             & ! 8
                                  '   ',             & ! 9
                                  '   ',             & ! 10
                                  '   ',             & ! 11
                                  '   ',             & ! 12
                                  'No '                        & ! 13
                                /),                            &
         TwoD                 = .false.,                       &
         UnrecognisedMessages = .true.,                        &
         DisjointColSpecs     = .true.,                        &
         HFBlockForms         = HFBlockForms                   &
       )

End Subroutine PdfReqsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2PdfReqs (Tokens, Grids, Reqs)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_), Intent(In)    :: Tokens
  Type(Grids_),  Intent(In)    :: Grids
  Type(Reqs_),   Intent(InOut) :: Reqs
  ! Local parameters:
  Integer, Parameter :: i_Name        = 1
  Integer, Parameter :: i_Species     = 2
  Integer, Parameter :: i_Source      = 3
  Integer, Parameter :: i_SourceGroup = 4
  Integer, Parameter :: i_PdfType     = 5
  Integer, Parameter :: i_HGrid       = 6
  Integer, Parameter :: i_ZGrid       = 7
  Integer, Parameter :: i_TGrid       = 8
  Integer, Parameter :: i_TAvOrInt    = 9
  Integer, Parameter :: i_AverageTime = 10
  Integer, Parameter :: i_Screen      = 11
  Integer, Parameter :: i_Disk        = 12
  Integer, Parameter :: i_Stat        = 13
  ! Locals:
  Integer                  :: PdfType
  Integer                  :: TAvOrInt
  Real(Std)                :: AvTime
  Real(Std)                :: PdfThresholds(MaxPdfSize)
  Logical                  :: Screen
  Logical                  :: Disk
  Logical                  :: Stat
  Type(PdfReq_)            :: PdfReq
  Integer                  :: iHGrid
  Integer                  :: iZGrid
  Character(MaxCharLength) :: HCoordName
  Character(MaxCharLength) :: ZCoordName

  PdfType = Token2Int(                                  &
              C         = Tokens%Tokens(i_PdfType),     &
              BlockKey  = 'Output Requirements - Pdfs', &
              Item      = Tokens%Tokens(i_Name),        &
              ColumnKey = 'Pdf Type'                    &
            )
  If (Tokens%Tokens(i_TAvOrInt) .CIEq. 'Av') Then
    TAvOrInt  = P_Av
  Else If (Tokens%Tokens(i_TAvOrInt) .CIEq. 'Int') Then
    TAvOrInt  = P_Int
  Else If ((Tokens%Tokens(i_TAvOrInt) .CIEq. ' ') .or. (Tokens%Tokens(i_TAvOrInt) .CIEq. 'No')) Then
    TAvOrInt  = 0
  Else
    Call Message('FATAL ERROR: "T Av Or Int" must be "Av", "Int", "No" or blank', 3)
  End If
  AvTime = Token2Std(                                  &
             C         = Tokens%Tokens(i_AverageTime), &
             BlockKey  = 'Output Requirements - Pdfs', &
             Item      = Tokens%Tokens(i_Name),        &
             ColumnKey = 'Average Time'                &
           )
  PdfThresholds(:) = 0.0
  Screen = Token2Log(                                  &
             C         = Tokens%Tokens(i_Screen),      &
             BlockKey  = 'Output Requirements - Pdfs', &
             Item      = Tokens%Tokens(i_Name),        &
             ColumnKey = 'Screen?'                     &
           )
  Disk   = Token2Log(                                  &
             C         = Tokens%Tokens(i_Disk),        &
             BlockKey  = 'Output Requirements - Pdfs', &
             Item      = Tokens%Tokens(i_Name),        &
             ColumnKey = 'Disk?'                       &
           )
  Stat   = Token2Log(                                  &
             C         = Tokens%Tokens(i_Stat),        &
             BlockKey  = 'Output Requirements - Pdfs', &
             Item      = Tokens%Tokens(i_Name),        &
             ColumnKey = 'Stat?'                       &
           )

  iHGrid     = FindHGridIndex(Tokens%Tokens(i_HGrid), Grids)
  HCoordName = Grids%HGrids(iHGrid)%HCoordName
  iZGrid     = FindZGridIndex(Tokens%Tokens(i_ZGrid), Grids)
  ZCoordName = Grids%ZGrids(iZGrid)%ZCoordName

  PdfReq = InitPdfReq(                                       &
             Name            = Tokens%Tokens(i_Name),        &
             SpeciesName     = Tokens%Tokens(i_Species),     &
             SourceName      = Tokens%Tokens(i_Source),      &
             SourceGroupName = Tokens%Tokens(i_SourceGroup), &
             PdfType         = PdfType,                      &
             HGridName       = Tokens%Tokens(i_HGrid),       &
             ZGridName       = Tokens%Tokens(i_ZGrid),       &
             TGridName       = Tokens%Tokens(i_TGrid),       &
             HCoordName      = HCoordName,                   &
             ZCoordName      = ZCoordName,                   &
             CalcType        = 1,                            & ! Hard wired for now $$
             TAvOrInt        = TAvOrInt,                     &
             AvTime          = AvTime,                       &
             PdfMode         = 1,                            & ! Hard wired for now $$
             PdfSize         = -1,                           &
             PdfThresholds   = PdfThresholds,                &
             Screen          = Screen,                       &
             Disk            = Disk,                         &
             AvEnsemble      = Stat                          &
           )

  Call AddPdfReq(PdfReq, Reqs)

End Subroutine Tokens2PdfReqs

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine PPInfoReqsInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                                               &
         BlockKey             = 'Output Requirements - Sets of Particle/Puff Information:', &
         NamedBlock           = .false.,                                                    &
         MatchMultipleLines   = .false.,                                                    &
         nHeaderLines         = 1,                                                          &
         ColumnKeys           = Reshape(                                                    &
                                  (/                                                        &
                                    'Output Name       ',                                   & ! 1 $$ call Name
                                    'Particles?        ',                                   & ! 2
                                    'Puffs?            ',                                   & ! 3
                                    'First Particle    ',                                   & ! 4
                                    'Last Particle     ',                                   & ! 5
                                    'First Puff        ',                                   & ! 6
                                    'Last Puff         ',                                   & ! 7
                                    'Met?              ',                                   & ! 8
                                    'Plume Rise?       ',                                   & ! 9
                                    'Dispersion Scheme?',                                   & ! 10
                                    'Puff Family?      ',                                   & ! 11
                                    'H-Coord           ',                                   & ! 12
                                    'Z-Coord           ',                                   & ! 13
                                    'T-Grid            ',                                   & ! 14
                                    'Sync?             ',                                   & ! 15
                                    'Screen?xxxx       ',                                   & ! 16
                                    'Disk?xxxx         ',                                   & ! 17
                                    'Mass?             ',                                   & ! 18
                                    'Source            ',                                   & ! 19
                                    'Output Route      ',                                   & ! 20
                                    'Output Format     '                                    & ! 21
                                  /),                                                       &
                                  (/ 21, 1 /)                                               &
                                ),                                                          &
         Defaults             = (/                                                          &
                                  '    ',             & ! 1
                                  '    ',             & ! 2
                                  '    ',             & ! 3
                                  '0   ',                                                   & ! 4
                                  '0   ',                                                   & ! 5
                                  '0   ',                                                   & ! 6
                                  '0   ',                                                   & ! 7
                                  '    ',             & ! 8
                                  '    ',             & ! 9
                                  '    ',             & ! 10
                                  '    ',             & ! 11
                                  '    ',             & ! 12
                                  '    ',             & ! 13
                                  '    ',             & ! 14
                                  '    ',             & ! 15
                                  '    ',             & ! 16
                                  '    ',             & ! 17
                                  '    ',             & ! 18
                                  '    ',             & ! 19
                                  '    ',                                                   & ! 20
                                  '    '                                                    & ! 21
                                /),                                                         &
         TwoD                 = .false.,                                                    &
         UnrecognisedMessages = .true.,                                                     &
         DisjointColSpecs     = .true.,                                                     &
         HFBlockForms         = HFBlockForms                                                &
       )

End Subroutine PPInfoReqsInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2PPInfoReqs(Tokens, Reqs, OpenMPOpts)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_), Intent(In)     :: Tokens
  Type(Reqs_),   Intent(InOut)  :: Reqs
  Type(OpenMPOpts_), Intent(In) :: OpenMPOpts
  ! Local parameters:
  Integer, Parameter :: i_OutputName       = 1 
  Integer, Parameter :: i_Particles        = 2
  Integer, Parameter :: i_Puffs            = 3
  Integer, Parameter :: i_FirstParticle    = 4
  Integer, Parameter :: i_LastParticle     = 5
  Integer, Parameter :: i_FirstPuff        = 6
  Integer, Parameter :: i_LastPuff         = 7
  Integer, Parameter :: i_Met              = 8
  Integer, Parameter :: i_PlumeRise        = 9
  Integer, Parameter :: i_DispersionScheme = 10
  Integer, Parameter :: i_PuffFamily       = 11
  Integer, Parameter :: i_HCoord           = 12
  Integer, Parameter :: i_ZCoord           = 13
  Integer, Parameter :: i_TGrid            = 14
  Integer, Parameter :: i_Sync             = 15
  Integer, Parameter :: i_Screenxxxx       = 16
  Integer, Parameter :: i_Diskxxxx         = 17
  Integer, Parameter :: i_Mass             = 18
  Integer, Parameter :: i_Source           = 19
  Integer, Parameter :: i_OutputRoute      = 20
  Integer, Parameter :: i_OutputFormat     = 21
  ! Locals:
  Logical          :: Particles
  Logical          :: Puffs
  Integer          :: FirstParticle
  Integer          :: LastParticle
  Integer          :: FirstPuff
  Integer          :: LastPuff
  Logical          :: Met
  Logical          :: Mass
  Logical          :: PlumeRise
  Logical          :: DispersionScheme
  Logical          :: PuffFamily
  Logical          :: Sync
  Logical          :: Graph
  Logical          :: Screen
  Logical          :: Disk
  Character(1)     :: CParticles
  Character(1)     :: CPuffs
  Type(PPInfoReq_) :: PPInfoReq !

  Particles = Token2Log(                                                               &
                C         = Tokens%Tokens(i_Particles),                                &
                BlockKey  = 'Output Requirements - Sets of Particle/Puff Information', &
                Item      = Tokens%Tokens(i_OutputName),                               &
                ColumnKey = 'Particles?'                                               &
              )
  Puffs     = Token2Log(                                                               &
                C         = Tokens%Tokens(i_Puffs),                                    &
                BlockKey  = 'Output Requirements - Sets of Particle/Puff Information', &
                Item      = Tokens%Tokens(i_OutputName),                               &
                ColumnKey = 'Puffs?'                                                   &
              )

  ! $$ avoid need to give zero for no/all particles
  ! $$ Ideally move to a single variable: all, none or n:m
  FirstParticle = Token2Int(                                                               &
                    C         = Tokens%Tokens(i_FirstParticle),                            &
                    BlockKey  = 'Output Requirements - Sets of Particle/Puff Information', &
                    Item      = Tokens%Tokens(i_OutputName),                               &
                    ColumnKey = 'First Particle'                                           &
                  )
  LastParticle  = Token2Int(                                                               &
                    C         = Tokens%Tokens(i_LastParticle),                             &
                    BlockKey  = 'Output Requirements - Sets of Particle/Puff Information', &
                    Item      = Tokens%Tokens(i_OutputName),                               &
                    ColumnKey = 'Last Particle'                                            &
                  )
  FirstPuff     = Token2Int(                                                               &
                    C         = Tokens%Tokens(i_FirstPuff),                                &
                    BlockKey  = 'Output Requirements - Sets of Particle/Puff Information', &
                    Item      = Tokens%Tokens(i_OutputName),                               &
                    ColumnKey = 'First Puff'                                               &
                  )
  LastPuff      = Token2Int(                                                               &
                    C         = Tokens%Tokens(i_LastPuff),                                 &
                    BlockKey  = 'Output Requirements - Sets of Particle/Puff Information', &
                    Item      = Tokens%Tokens(i_OutputName),                               &
                    ColumnKey = 'Last Puff'                                                &
                  )

  If (Particles .and. FirstParticle == 0 .and. LastParticle == 0) Then
    CParticles = 'A' ! All
  Else If (Particles) Then
    CParticles = 'S' ! Some
  Else
    CParticles = 'N' ! None
  End If
  If (Puffs .and. FirstPuff == 0 .and. LastPuff == 0) Then
    CPuffs = 'A' ! All
  Else If (Puffs) Then
    CPuffs = 'S' ! Some
  Else
    CPuffs = 'N' ! None
  End If

  Met              = Token2Log(                                                               &
                       C         = Tokens%Tokens(i_Met),                                      &
                       BlockKey  = 'Output Requirements - Sets of Particle/Puff Information', &
                       Item      = Tokens%Tokens(i_OutputName),                               &
                       ColumnKey = 'Met?'                                                     &
                     )
  Mass             = Token2Log(                                                               &
                       C         = Tokens%Tokens(i_Mass),                                     &
                       BlockKey  = 'Output Requirements - Sets of Particle/Puff Information', &
                       Item      = Tokens%Tokens(i_OutputName),                               &
                       ColumnKey = 'Mass?'                                                    &
                     )
  PlumeRise        = Token2Log(                                                               &
                       C         = Tokens%Tokens(i_PlumeRise),                                &
                       BlockKey  = 'Output Requirements - Sets of Particle/Puff Information', &
                       Item      = Tokens%Tokens(i_OutputName),                               &
                       ColumnKey = 'Plume Rise?'                                              &
                     )
  DispersionScheme = Token2Log(                                                               &
                       C         = Tokens%Tokens(i_DispersionScheme),                         &
                       BlockKey  = 'Output Requirements - Sets of Particle/Puff Information', &
                       Item      = Tokens%Tokens(i_OutputName),                               &
                       ColumnKey = 'Dispersion Scheme?'                                       &
                     )
  If (Puffs) Then
    PuffFamily = Token2Log(                                                               &
                   C         = Tokens%Tokens(i_PuffFamily),                               &
                   BlockKey  = 'Output Requirements - Sets of Particle/Puff Information', &
                   Item      = Tokens%Tokens(i_OutputName),                               &
                   ColumnKey = 'Puff Family?'                                             &
                 )
  Else
    PuffFamily = .false.
  End If

  Sync = Token2Log(                                                               &
           C         = Tokens%Tokens(i_Sync),                                     &
           BlockKey  = 'Output Requirements - Sets of Particle/Puff Information', &
           Item      = Tokens%Tokens(i_OutputName),                               &
           ColumnKey = 'Sync?'                                                    &
         )

  ! When compiled with Intel and running on more than one thread, 
  ! trajectory output can only be requested at sync times.
  ! Write a warning message to the user if Sync is set to false in any
  ! request for particle-puff information.
 
#ifdef IntelLinCompiler
  If (OpenMPOpts%nParticleThreads > 1 .and. (.not.Sync)) Then
    Call Message('WARNING: When running with more than one particle thread, '      // &
                 'Particle/Puff Information can only be requested at Sync times. ' // &
                 'Otherwise, the code will crash with a segmentation fault.',         &
                 1)
  End If
#endif

  Call TokenLengthTest(                              & !$$ this is a check for blank strings, such as might 
         C         = Tokens%Tokens(i_OutputRoute),   & !   occur if user was using 'Disk?' instead of 'Output Route'.
         Length    = MaxCharLength,                  & !   Disk? has been removed as an option.
         Zero      = .true.,                         & !   Generally we don't check for blanks in Input.F90
         BlockKey  = 'Output Requirements - Fields', & !   but perhaps should review this.
         Item      = Tokens%Tokens(i_OutputName),    &
         ColumnKey = 'Output Route'                  &
       ) 

  PPInfoReq = InitPPInfoReq(                                     &
                Name             = Tokens%Tokens(i_OutputName),  &
                Particles        = CParticles,                   &
                Puffs            = CPuffs,                       &
                FirstParticle    = FirstParticle,                &
                LastParticle     = LastParticle,                 &
                FirstPuff        = FirstPuff,                    &
                LastPuff         = LastPuff,                     &
                SourceName       = Tokens%Tokens(i_Source),      &
                Met              = Met,                          &
                Mass             = Mass,                         &
                PlumeRise        = PlumeRise,                    &
                DispersionScheme = DispersionScheme,             &
                PuffFamily       = PuffFamily,                   &
                HCoordName       = Tokens%Tokens(i_HCoord),      &
                ZCoordName       = Tokens%Tokens(i_ZCoord),      &
                TGridName        = Tokens%Tokens(i_TGrid),       &
                Sync             = Sync,                         &
                OutputRoute      = Tokens%Tokens(i_OutputRoute), &
                OutputFormat     = Tokens%Tokens(i_OutputFormat) &
              )

  Call AddPPInfoReq(PPInfoReq, Reqs)

End Subroutine Tokens2PPInfoReqs

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine DispOptsesInputNames(HFBlockForms)
! .

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(InOut) :: HFBlockForms

  Call InitAndAddHFBlockForm(                                         &
         BlockKey             = 'Sets of Dispersion Options:',        &
         NamedBlock           = .false.,                              &
         MatchMultipleLines   = .false.,                              &
         nHeaderLines         = 1,                                    &
         ColumnKeys           = Reshape(                              &
                                  (/                                  &
                                    'Skew Time                     ', & ! 1
                                    'Velocity Memory Time          ', & ! 2
                                    'Inhomogeneous Time            ', & ! 3
                                    'Puff Time                     ', & ! 4
                                    'Sync Time                     ', & ! 5
                                    'Computational Domain          ', & ! 6
                                    'Puff Interval                 ', & ! 7
                                    'DeltaOpt                      ', & ! 8
                                    'Deep Convection?              ', & ! 9
                                    'Radioactive Decay?            ', & ! 10
                                    'Agent Decay?                  ', & ! 11
                                    'Dry Deposition?               ', & ! 12
                                    'Wet Deposition?               ', & ! 13
                                    'Mesoscale Motions?            ', & ! 14
                                    'Chemistry?                    ', & ! 15
                                    'Turbulence?                   ', & ! 16
                                    'Mesoscale Velocity Memory Time', & ! 17
                                    'Damping?                      ', & ! 18
                                    'Max # Particles               ', & ! 19
                                    'Max # Full Particles          ', & ! 20
                                    'Max # Puffs                   ', & ! 21
                                    'Max # Original Puffs          ', & ! 22
                                    'Time Of Fixed Met             ', & ! 23
                                    'Particle Ceiling              ', & ! 24
                                    'Particle Factor               ', & ! 25
                                    'A1                            ', & ! 26
                                    'A5                            ', & ! 27
                                    'A7                            ', & ! 28
                                    'Vertical Velocity?            ', & ! 29
                                    'Max Deposition Height         '  & ! 30
                                  /),                                 &
                                  (/ 30, 1 /)                         &
                                ),                                    &
         Defaults             = (/                                    &
                                  '       ', & ! 1
                                  '       ', & ! 2
                                  '       ', & ! 3
                                  '       ', & ! 4
                                  '       ', & ! 5
                                  '       ', & ! 6
                                  '       ', & ! 7
                                  '       ', & ! 8
                                  '       ', & ! 9
                                  '       ', & ! 10
                                  '       ', & ! 11
                                  '       ', & ! 12
                                  '       ', & ! 13
                                  '       ', & ! 14
                                  '       ', & ! 15
                                  '       ', & ! 16
                                  '       ', & ! 17
                                  '.true. ',                        & ! 18
                                  '       ', & ! 19
                                  '       ', & ! 20
                                  '       ', & ! 21
                                  '       ', & ! 22
                                  '       ', & ! 23
                                  '       ', & ! 24
                                  '1.0    ',                        & ! 25
                                  '50.0   ',                        & ! 26
                                  '3.0    ',                        & ! 27
                                  '2.0    ',                        & ! 28
                                  '.true. ',                        & ! 29
                                  '1.0E20 '                         & ! 30
                                /),                                 &
         TwoD                 = .false.,                            &
         UnrecognisedMessages = .true.,                             &
         DisjointColSpecs     = .true.,                             &
         HFBlockForms         = HFBlockForms                        &
       )

End Subroutine DispOptsesInputNames

!-----------------------------------------------------------------------------------------------------------------------------------

Subroutine Tokens2DispOptses(Tokens, FixedMet, DispOptses)
!.

  Implicit None
  ! Argument list:
  Type(Tokens_),     Intent(In)    :: Tokens
  Logical,           Intent(In)    :: FixedMet
  Type(DispOptses_), Intent(InOut) :: DispOptses
  ! Local parameters:
  Integer, Parameter :: i_SkewTime                    = 1
  Integer, Parameter :: i_VelocityMemoryTime          = 2
  Integer, Parameter :: i_InhomogeneousTime           = 3
  Integer, Parameter :: i_PuffTime                    = 4
  Integer, Parameter :: i_SyncTime                    = 5
  Integer, Parameter :: i_ComputationalDomain         = 6
  Integer, Parameter :: i_PuffInterval                = 7
  Integer, Parameter :: i_DeltaOpt                    = 8
  Integer, Parameter :: i_DeepConvection              = 9
  Integer, Parameter :: i_RadioactiveDecay            = 10
  Integer, Parameter :: i_AgentDecay                  = 11
  Integer, Parameter :: i_DryDeposition               = 12
  Integer, Parameter :: i_WetDeposition               = 13
  Integer, Parameter :: i_MesoscaleMotions            = 14
  Integer, Parameter :: i_Chemistry                   = 15
  Integer, Parameter :: i_Turbulence                  = 16
  Integer, Parameter :: i_MesoscaleVelocityMemoryTime = 17
  Integer, Parameter :: i_Damping                     = 18
  Integer, Parameter :: i_MaxNParticles               = 19
  Integer, Parameter :: i_MaxNFullParticles           = 20
  Integer, Parameter :: i_MaxNPuffs                   = 21
  Integer, Parameter :: i_MaxNOriginalPuffs           = 22
  Integer, Parameter :: i_TimeOfFixedMet              = 23
  Integer, Parameter :: i_ParticleCeiling             = 24
  Integer, Parameter :: i_ParticleFactor              = 25
  Integer, Parameter :: i_A1                          = 26
  Integer, Parameter :: i_A5                          = 27
  Integer, Parameter :: i_A7                          = 28
  Integer, Parameter :: i_VerticalVelocity            = 29
  Integer, Parameter :: i_MaxDepositionHeight         = 30
  ! Locals:
  Integer                  :: MaxParticles
  Integer                  :: MaxFullParticles
  Integer                  :: MaxPuffs
  Integer                  :: MaxOriginalPuffs
  Integer                  :: ParticleCeiling
  Real(Std)                :: ParticleFactor
  Type(Time_)              :: SkewTime
  Type(Time_)              :: VelMemTime
  Type(Time_)              :: InhomogTime
  Type(Time_)              :: MVelMemTime
  Type(Time_)              :: PuffTime
  Type(Time_)              :: SyncTime
  Type(Time_)              :: PuffInterval
  Integer                  :: DeltaOpt
  Character(MaxCharLength) :: DeepConvection
  Logical                  :: RadioactiveDecay
  Logical                  :: AgentDecay
  Logical                  :: DryDep
  Logical                  :: WetDep
  Logical                  :: MesoscaleMotions
  Logical                  :: Turbulence
  Logical                  :: VerticalVelocity
  Logical                  :: Damping
  Logical                  :: Chemistry
  Type(Time_)              :: FixedMetTime
  Real(Std)                :: A1
  Real(Std)                :: A5
  Real(Std)                :: A7
  Real(Std)                :: Zs
  Type(DispOpts_)          :: DispOpts

  If (FixedMet) Then
    If (Tokens%Tokens(i_TimeOfFixedMet) == ' ') Then
      Call Message('Error in Tokens2DispOptses', 3)
    End If
    FixedMetTime = Token2Time(                                    &
                     C         = Tokens%Tokens(i_TimeOfFixedMet), &
                     BlockKey  = 'Sets of Dispersion Options',    &
                     Item      = ' ',                             &
                     ColumnKey = 'Time of Fixed Met'              &
                   )
  Else
    If (Tokens%Tokens(i_TimeOfFixedMet) /= ' ') Then
      Call Message('Error in Tokens2DispOptses', 3)
    End If
  End If

  DeltaOpt = Token2Int(                                  &
               C         = Tokens%Tokens(i_DeltaOpt),    &
               BlockKey  = 'Sets of Dispersion Options', &
               Item      = ' ',                          &
               ColumnKey = 'DeltaOpt'                    &
             )

  MaxParticles     = Token2Int(                                    &
                       C         = Tokens%Tokens(i_MaxNParticles), &
                       BlockKey  = 'Sets of Dispersion Options',   &
                       Item      = ' ',                            &
                       ColumnKey = 'Max # Particles'               &
                     )

  MaxFullParticles = Token2Int(                                        &
                       C         = Tokens%Tokens(i_MaxNFullParticles), &
                       BlockKey  = 'Sets of Dispersion Options',       &
                       Item      = ' ',                                &
                       ColumnKey = 'Max # Full Particles'              &
                     )

  MaxPuffs         = Token2Int(                                  &
                       C         = Tokens%Tokens(i_MaxNPuffs),   &
                       BlockKey  = 'Sets of Dispersion Options', &
                       Item      = ' ',                          &
                       ColumnKey = 'Max # Puffs'                 &
                     )

  MaxOriginalPuffs = Token2Int(                                        &
                       C         = Tokens%Tokens(i_MaxNOriginalPuffs), &
                       BlockKey  = 'Sets of Dispersion Options',       &
                       Item      = ' ',                                &
                       ColumnKey = 'Max # Original Puffs'              &
                     )

  If (Tokens%Tokens(i_ParticleCeiling) == ' ') Then
    ParticleCeiling = MaxParticles
  Else
    ParticleCeiling = Token2Int(                                      &
                        C         = Tokens%Tokens(i_ParticleCeiling), &
                        BlockKey  = 'Sets of Dispersion Options',     &
                        Item      = ' ',                              &
                        ColumnKey = 'Particle Ceiling'                &
                      )
  End If

  ParticleFactor = Token2Std(                                     &
                     C         = Tokens%Tokens(i_ParticleFactor), &
                     BlockKey  = 'Sets of Dispersion Options',    &
                     Item      = ' ',                             &
                     ColumnKey = 'Particle Factor'                &
                   )

  SkewTime     = Token2Time(                                 &
                   C         = Tokens%Tokens(i_SkewTime),    &
                   BlockKey  = 'Sets of Dispersion Options', &
                   Item      = ' ',                          &
                   ColumnKey = 'Skew Time',                  &
                   Interval  = .true.                        &
                 )

  VelMemTime   = Token2Time(                                        &
                   C         = Tokens%Tokens(i_VelocityMemoryTime), &
                   BlockKey  = 'Sets of Dispersion Options',        &
                   Item      = ' ',                                 &
                   ColumnKey = 'Velocity Memory Time',              &
                   Interval  = .true.                               &
                 )

  InhomogTime  = Token2Time(                                       &
                   C         = Tokens%Tokens(i_InhomogeneousTime), &
                   BlockKey  = 'Sets of Dispersion Options',       &
                   Item      = ' ',                                &
                   ColumnKey = 'Inhomogeneous Time',               &
                   Interval  = .true.                              &
                 )

  PuffTime     = Token2Time(                                 &
                   C         = Tokens%Tokens(i_PuffTime),    &
                   BlockKey  = 'Sets of Dispersion Options', &
                   Item      = ' ',                          &
                   ColumnKey = 'Puff Time',                  &
                   Interval  = .true.                        &
                 )

  SyncTime     = Token2Time(                                 &
                   C         = Tokens%Tokens(i_SyncTime),    &
                   BlockKey  = 'Sets of Dispersion Options', &
                   Item      = ' ',                          &
                   ColumnKey = 'Sync Time',                  &
                   Interval  = .true.                        &
                 )

  MVelMemTime  = Token2Time(                                                 &
                   C         = Tokens%Tokens(i_MesoscaleVelocityMemoryTime), &
                   BlockKey  = 'Sets of Dispersion Options',                 &
                   Item      = ' ',                                          &
                   ColumnKey = 'Mesoscale Velocity Memory Time',             &
                   Interval  = .true.                                        &
                 )

  PuffInterval = Token2Time(                                  &
                   C         = Tokens%Tokens(i_PuffInterval), &
                   BlockKey  = 'Sets of Dispersion Options',  &
                   Item      = ' ',                           &
                   ColumnKey = 'Puff Interval',               &
                   Interval  = .true.                         &
                 )
  ! $$ currently must be defined as used in MaxOutputMemoryTime.

  DeepConvection   = Tokens%Tokens(i_DeepConvection)
                      
  RadioactiveDecay = Token2Log(                                       &
                       C         = Tokens%Tokens(i_RadioactiveDecay), &
                       BlockKey  = 'Sets of Dispersion Options',      &
                       Item      = ' ',                               &
                       ColumnKey = 'Radioactive Decay?'               &
                     )

  AgentDecay       = Token2Log(                                  &
                       C         = Tokens%Tokens(i_AgentDecay),  &
                       BlockKey  = 'Sets of Dispersion Options', &
                       Item      = ' ',                          &
                       ColumnKey = 'Agent Decay?'                &
                     )

  DryDep           = Token2Log(                                    &
                       C         = Tokens%Tokens(i_DryDeposition), &
                       BlockKey  = 'Sets of Dispersion Options',   &
                       Item      = ' ',                            &
                       ColumnKey = 'Dry Deposition?'               &
                     )

  WetDep           = Token2Log(                                    &
                       C         = Tokens%Tokens(i_WetDeposition), &
                       BlockKey  = 'Sets of Dispersion Options',   &
                       Item      = ' ',                            &
                       ColumnKey = 'Wet Deposition?'               &
                     )

  MesoscaleMotions = Token2Log(                                       &
                       C         = Tokens%Tokens(i_MesoscaleMotions), &
                       BlockKey  = 'Sets of Dispersion Options',      &
                       Item      = ' ',                               &
                       ColumnKey = 'Mesoscale Motions?'               &
                     )

  Chemistry        = Token2Log(                                  &
                       C         = Tokens%Tokens(i_Chemistry),   &
                       BlockKey  = 'Sets of Dispersion Options', &
                       Item      = ' ',                          &
                       ColumnKey = 'Chemistry?'                  &
                     )

  Turbulence       = Token2Log(                                  &
                       C         = Tokens%Tokens(i_Turbulence),  &
                       BlockKey  = 'Sets of Dispersion Options', &
                       Item      = ' ',                          &
                       ColumnKey = 'Turbulence?'                 &
                     )

  Damping          = Token2Log(                                  &
                       C         = Tokens%Tokens(i_Damping),     &
                       BlockKey  = 'Sets of Dispersion Options', &
                       Item      = ' ',                          &
                       ColumnKey = 'Damping?'                    &
                     )

  VerticalVelocity = Token2Log(                                       &
                       C         = Tokens%Tokens(i_VerticalVelocity), &
                       BlockKey  = 'Sets of Dispersion Options',      &
                       Item      = ' ',                               &
                       ColumnKey = 'Vertical Velocity?'               &
                     )

  A1 = Token2Std(                                  &
         C         = Tokens%Tokens(i_A1),          &
         BlockKey  = 'Sets of Dispersion Options', &
         Item      = ' ',                          &
         ColumnKey = 'A1'                          &
       )

  A5 = Token2Std(                                  &
         C         = Tokens%Tokens(i_A5),          &
         BlockKey  = 'Sets of Dispersion Options', &
         Item      = ' ',                          &
         ColumnKey = 'A5'                          &
       )

  A7 = Token2Std(                                  &
         C         = Tokens%Tokens(i_A7),          &
         BlockKey  = 'Sets of Dispersion Options', &
         Item      = ' ',                          &
         ColumnKey = 'A7'                          &
       )

  Zs = Token2Std(                                          &
         C         = Tokens%Tokens(i_MaxDepositionHeight), &
         BlockKey  = 'Sets of Dispersion Options',         &
         Item      = ' ',                                  &
         ColumnKey = 'Max Deposition Height'               &
       )

  Call InitDispOpts(                                              &
         FixedMetTime     = FixedMetTime,                         &
         MaxParticles     = MaxParticles,                         &
         MaxFullParticles = MaxFullParticles,                     &
         MaxPuffs         = MaxPuffs,                             &
         MaxOriginalPuffs = MaxOriginalPuffs,                     &
         ParticleCeiling  = ParticleCeiling,                      &
         ParticleFactor   = ParticleFactor,                       &
         SkewTime         = SkewTime,                             &
         VelMemTime       = VelMemTime,                           &
         InhomogTime      = InhomogTime,                          &
         MVelMemTime      = MVelMemTime,                          &
         PuffTime         = PuffTime,                             &
         SyncdT           = SyncTime,                             &
         DomainName       = Tokens%Tokens(i_ComputationalDomain), &
         PuffInterval     = PuffInterval,                         &
         DeltaOpt         = DeltaOpt,                             &
         DeepConvection   = DeepConvection,                       &
         RadioactiveDecay = RadioactiveDecay,                     &
         AgentDecay       = AgentDecay,                           &
         DryDep           = DryDep,                               &
         WetDep           = WetDep,                               &
         Turbulence       = Turbulence,                           &
         MesoscaleMotions = MesoscaleMotions,                     &
         Damping          = Damping,                              &
         VerticalVelocity = VerticalVelocity,                     &
         Chemistry        = Chemistry,                            &
         A1               = A1,                                   &
         A5               = A5,                                   &
         A7               = A7,                                   &
         Zs               = Zs,                                   &
         DispOpts         = DispOpts                              &
       )
  Call AddDispOpts(DispOpts, DispOptses)

End Subroutine Tokens2DispOptses

!-----------------------------------------------------------------------------------------------------------------------------------

End Module InputModule

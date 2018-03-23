! Module:  Restart Module

Module RestartModule

! This module provides code to control the restartability of the model.

! Three restart catalogues are used to keep information on available valid restart
! files and restart catalogue validity files are used to keep information on the
! validity of the restart catalogues.

! The catalogues can be 'valid', 'invalid' or 'valid only as a list of files which
! might need to be deleted' ('valid for delete' for short). If more than one is valid,
! the lowest numbered one counts. If more than one is valid for delete, the lowest
! numbered one counts.

! Being the lowest numbered valid catalogue usually means it is a correct catalogue.
! However the catalogue can contain files which have been deleted by the user.

! A non-existent validity file is equivalent to one indicating an invalid catalogue.
! A non-existent catalogue is equivalent to an empty catalogue (which may or may not
! be a valid and/or correct catalogue).

! Currently any alterations to the command line argument list (other than the -Restart
! flag), input, output or restart files (other than deleting them) between restarts
! can have unpredictable effects and is not recommended.

! Deleting any of the restart files is safe (in that it should not result in incorrect
! answers). However the effects may not be predictable without a detailed
! understanding of the restart code. In particular deleting a restart catalogue or
! restart catalogue validity file can have the effect of making the restart files
! unusable and may result in some of them being deleted when Name III is started or
! restarted. Its best to either delete all restart files (catalogues, catalogue
! validity files and individual restart files) which will force the model to start at
! the beginning, or to delete individual restart files which can be useful to save
! space and will do nothing except preventing a restart from that file.

! Any restart files that are not in a catalogue (e.g. from a previous run) are ignored
! but may be overwritten. Also, notwithstanding any comments about safety above,
! neither the restart mechanism nor Name III as a whole is robust against problems
! caused by two different runs attempting to use the same files.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use OutputModule, Only: OutputOpts_, Results_, WriteRestartFileResults, ReadRestartFileResults, Reqs_
Use ChemistryModule
Use CaseModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: RestartOpts_             ! Restart options.
Public :: PreInitRestartOpts       ! Pre-initialises the restart options.
Public :: InitRestartOpts          ! Initialises the restart options.
Public :: CheckRestartOpts         ! Checks the restart options.
Public :: SetUpRestartOpts         ! Sets up RestartOpts.
Public :: PreInitRestart           ! Pre-initialises the restart mechanism.
Public :: InitRestart              ! Initialises the restart mechanism.
Public :: FindNextRestartWriteTime
Public :: WriteRestartFile
Public :: ReadRestartFile

!-------------------------------------------------------------------------------------------------------------

Type :: RestartOpts_ ! Restart options.
  Logical                      :: Initialised
  Integer                      :: RestartType
  Character(MaxFileNameLength) :: FileStem
  Integer                      :: dCase
  Type(Time_)                  :: T0
  Type(Time_)                  :: dT
  Logical                      :: WriteOnSuspend
  Logical                      :: DeleteOldFiles
  ! Initialised    :: Indicates the restart options have been initialised.
  ! RestartType    :: 0 = none, 1 = every so many cases, 2 = within case, 3 = on suspend, but not 1 or 2
  ! FileStem       :: File stem for restart files.
  ! dCase          :: Number of cases between writing restart files.
  ! dT             :: Time interval between writing restart files.
  ! WriteOnSuspend :: Write restart file when run suspended.
End Type RestartOpts_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function PreInitRestartOpts() Result(RestartOpts)
! Pre-initialises the restart options.

  Implicit None
  ! Function result:
  Type(RestartOpts_) :: RestartOpts ! Restart options.

  RestartOpts%Initialised = .false.

End Function PreInitRestartOpts

!-------------------------------------------------------------------------------------------------------------

Subroutine InitRestartOpts(  &
             RestartType,    &
             dCase, T0, dT,  &
             DeleteOldFiles, &
             WriteOnSuspend, &
             RestartOpts     &
           )
! Initialises the restart options.

  Implicit None
  ! Argument list:
  Integer,            Intent(In)    :: RestartType
  Integer,            Intent(In)    :: dCase
  Type(Time_),        Intent(In)    :: T0
  Type(Time_),        Intent(In)    :: dT
  Logical,            Intent(In)    :: DeleteOldFiles
  Logical,            Intent(In)    :: WriteOnSuspend
  Type(RestartOpts_), Intent(InOut) :: RestartOpts
  ! RestartType ::
  ! dCase       ::
  ! dT          ::
  ! RestartOpts :: Restart options.

  If (RestartOpts%Initialised) Then
    Call Message(                                                             &
           'FATAL ERROR: an attempt has been made to '                     // &
           're-initialise the restart options - check that there is only ' // &
           'one set of restart options specified in the input',               &
           3                                                                  &
         )
  End If

  RestartOpts%RestartType = RestartType
  If (RestartType == 0 .or. RestartType == 3) Then ! $$ do we need restarttype? Could some of this be
                                                   ! $$ in input?
    RestartOpts%dCase = Huge(RestartOpts%dCase) / 2
  Else If (RestartType == 2) Then
    RestartOpts%dCase = 1
  Else
    RestartOpts%dCase = dCase
  End If
  RestartOpts%T0             = T0 ! $$ Infinite for restarttype 1?
  RestartOpts%dT             = dT ! $$
  RestartOpts%WriteOnSuspend = WriteOnSuspend
  RestartOpts%DeleteOldFiles = DeleteOldFiles
 ! If (WriteOnSuspend .and. RestartOpts%RestartType == 0) RestartOpts%RestartType = 3
  If (RestartType == 1) Then
    If (dCase <= 0) Then
      Call Message(                                                &
             'FATAL ERROR: # Cases Between Writes cannot be <= 0', &
             3                                                     &
           )
    End If
  End If
  If (RestartType == 2) Then
    If (dT <= ZeroTime()) Then
      Call Message(                                             &
             'FATAL ERROR: Time Between Writes cannot be <= 0', &
             3                                                  &
           )
    End If
  End If

  RestartOpts%Initialised = .true.

End Subroutine InitRestartOpts

!-------------------------------------------------------------------------------------------------------------

Subroutine CheckRestartOpts(Restart, RestartOpts)
! Checks the restart options.

  Implicit None
  ! Argument list:
  Logical,            Intent(In) :: Restart     ! Indicates the run is being restarted.
  Type(RestartOpts_), Intent(In) :: RestartOpts ! Restart options.

  If (.not.RestartOpts%Initialised) Then
    Call Message('FATAL ERROR: No Restart options block found in input files', 3)
  End If
  If (Restart .and. RestartOpts%RestartType == 0) Then
    Call Message(                                                           &
           'FATAL ERROR: The -Restart flag is used on the command line ' // &
           'but restart options are not set in the input file',             &
           3                                                                &
         )
  End If

End Subroutine CheckRestartOpts

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpRestartOpts(FileStem, RestartOpts)
! Sets up RestartOpts.

  Implicit None
  ! Argument list:
  Character(MaxFileNameLength), Intent(In)    :: FileStem
  Type(RestartOpts_),           Intent(InOut) :: RestartOpts ! Restart options.

 RestartOpts%FileStem = FileStem

End Subroutine SetUpRestartOpts

!-------------------------------------------------------------------------------------------------------------

Subroutine PreInitRestart(FileStem, Restart)
! Pre-initialises the restart mechanism.

! Checks whether Restart is consistent with existence of restart files (more precisely
! restart catalogue validity files). Restart set but restart files not present is
! allowed - in this case the run 'restarts' at the beginning. Restart not set but
! restart files present is not allowed. If a restart was intended this prevents
! overwriting of the restart files and, if called before SetUpErrorAndMessageModule,
! prevents overwriting of the log file. If a restart was not intended, calling before
! SetUpErrorAndMessageModule prevents the possibility of the following sequence:
! overwriting of the log file, a fatal error or user action to stop the run before the
! first call to InitRestart is completed, and a restart which restarts from the
! previous run's restart files but which uses the current input files and looks, from
! the log file, like a restart of the new run.

  Implicit None
  ! Argument list:
  Logical,                      Intent(In) :: Restart
  Character(MaxFileNameLength), Intent(In) :: FileStem
  ! Locals:
  Logical :: Exist1
  Logical :: Exist2
  Logical :: Exist3

  Inquire(                                                     &
    File  = Trim(FileStem) // 'RestartCatalogueValidity1.txt', &
    Exist = Exist1                                             &
  )
  Inquire(                                                     &
    File  = Trim(FileStem) // 'RestartCatalogueValidity2.txt', &
    Exist = Exist2                                             &
  )
  Inquire(                                                     &
    File  = Trim(FileStem) // 'RestartCatalogueValidity3.txt', &
    Exist = Exist3                                             &
  )
  If (.not.Restart .and. (Exist1 .or. Exist2 .or. Exist3)) Then
    Call Message(                                                               &
      'FATAL ERROR: Restart catalogue validity files exist but the restart ' // &
      'command line flag is not set',                                           &
      3                                                                         &
    )
  End If

End Subroutine PreInitRestart

!-------------------------------------------------------------------------------------------------------------

Subroutine InitRestart(RestartOpts, Restart, Units)
! Initialises the restart mechanism. by ensuring a consistent set of restart
! catalogues.

  Implicit None
  ! Argument list:
  Type(RestartOpts_), Intent(In)    :: RestartOpts
  Logical,            Intent(InOut) :: Restart
  Type(Units_),       Intent(InOut) :: Units
  ! Locals:
  Logical :: Valid1
  Logical :: Valid2
  Logical :: Valid3
  Logical :: ValidForDelete1
  Logical :: ValidForDelete2
  Logical :: ValidForDelete3

  If (.not.Restart) Then

    Call EmptyCatalogue    (RestartOpts%FileStem, '1', Units)
    Call MarkCatalogueValid(RestartOpts%FileStem, '1', Units)
    Call EmptyCatalogue    (RestartOpts%FileStem, '2', Units)
    Call MarkCatalogueValid(RestartOpts%FileStem, '2', Units)
    Call EmptyCatalogue    (RestartOpts%FileStem, '3', Units)
    Call MarkCatalogueValid(RestartOpts%FileStem, '3', Units)

  Else

    Valid1          = IsCatalogueValid         (RestartOpts%FileStem, '1', Units)
    Valid2          = IsCatalogueValid         (RestartOpts%FileStem, '2', Units)
    Valid3          = IsCatalogueValid         (RestartOpts%FileStem, '3', Units)
    ValidForDelete1 = IsCatalogueValidForDelete(RestartOpts%FileStem, '1', Units)
    ValidForDelete2 = IsCatalogueValidForDelete(RestartOpts%FileStem, '2', Units)
    ValidForDelete3 = IsCatalogueValidForDelete(RestartOpts%FileStem, '3', Units)

    If (.not.Valid1 .and. .not.Valid2 .and. .not.Valid2) Then

      If (ValidForDelete1) Then
        Call DeleteCataloguedFiles(RestartOpts%FileStem, '0', '1', Units)
      End If
      If (ValidForDelete2) Then
        Call DeleteCataloguedFiles(RestartOpts%FileStem, '0', '2', Units)
      End If
      If (ValidForDelete3) Then
        Call DeleteCataloguedFiles(RestartOpts%FileStem, '0', '3', Units)
      End If
      Call EmptyCatalogue    (RestartOpts%FileStem, '1', Units)
      Call MarkCatalogueValid(RestartOpts%FileStem, '1', Units)
      Call EmptyCatalogue    (RestartOpts%FileStem, '2', Units)
      Call MarkCatalogueValid(RestartOpts%FileStem, '2', Units)
      Call EmptyCatalogue    (RestartOpts%FileStem, '3', Units)
      Call MarkCatalogueValid(RestartOpts%FileStem, '3', Units)
      Call Message(' ')
      Call Message(                                  &
             'No valid restart catalogue found. ' // &
             'Run restarting at the begining.'       &
           )
      Restart = .false.

    Else If (Valid1) Then

      If (.not.Valid2 .and. ValidForDelete2) Then
        Call DeleteCataloguedFiles(RestartOpts%FileStem, '1', '2', Units)
      End If
      If (.not.Valid3 .and. ValidForDelete3) Then
        Call DeleteCataloguedFiles(RestartOpts%FileStem, '1', '3', Units)
      End If
      Call MarkCatalogueInValid(RestartOpts%FileStem,      '2',      Units)
      Call CopyCatalogue       (RestartOpts%FileStem, '1', '2', ' ', Units)
      Call MarkCatalogueValid  (RestartOpts%FileStem,      '2',      Units)
      Call MarkCatalogueInValid(RestartOpts%FileStem,      '3',      Units)
      Call CopyCatalogue       (RestartOpts%FileStem, '1', '3', ' ', Units)
      Call MarkCatalogueValid  (RestartOpts%FileStem,      '3',      Units)

    Else If (Valid2) Then

      If (.not.Valid3 .and. ValidForDelete3) Then
        Call DeleteCataloguedFiles(RestartOpts%FileStem, '2', '3', Units)
      End If
      If (.not.Valid1 .and. ValidForDelete1) Then
        Call DeleteCataloguedFiles(RestartOpts%FileStem, '2', '1', Units)
      End If
      Call MarkCatalogueInValid(RestartOpts%FileStem,      '3',      Units)
      Call CopyCatalogue       (RestartOpts%FileStem, '2', '3', ' ', Units)
      Call MarkCatalogueValid  (RestartOpts%FileStem,      '3',      Units)
      Call MarkCatalogueInValid(RestartOpts%FileStem,      '1',      Units)
      Call CopyCatalogue       (RestartOpts%FileStem, '2', '1', ' ', Units)
      Call MarkCatalogueValid  (RestartOpts%FileStem,      '1',      Units)

    Else If (Valid3) Then

      If (.not.Valid2 .and. ValidForDelete2) Then
        Call DeleteCataloguedFiles(RestartOpts%FileStem, '3', '2', Units)
      End If
      If (.not.Valid1 .and. ValidForDelete1) Then
        Call DeleteCataloguedFiles(RestartOpts%FileStem, '3', '1', Units)
      End If
      Call MarkCatalogueInValid(RestartOpts%FileStem,      '2',      Units)
      Call CopyCatalogue       (RestartOpts%FileStem, '3', '2', ' ', Units)
      Call MarkCatalogueValid  (RestartOpts%FileStem,      '2',      Units)
      Call MarkCatalogueInValid(RestartOpts%FileStem,      '1',      Units)
      Call CopyCatalogue       (RestartOpts%FileStem, '3', '1', ' ', Units)
      Call MarkCatalogueValid  (RestartOpts%FileStem,      '1',      Units)

    End If

  End If

End Subroutine InitRestart

!-------------------------------------------------------------------------------------------------------------

Function FindNextRestartWriteTime(Time, RestartOpts) Result(NextRestartWriteTime)

  Implicit None
  ! Argument list:
  Type(Time_),        Intent(In) :: Time
  Type(RestartOpts_), Intent(In) :: RestartOpts
  ! Function result:
  Type(Time_) :: NextRestartWriteTime

  If (RestartOpts%RestartType == 2) Then
    NextRestartWriteTime = Round(            &
                             Time,           &
                             RestartOpts%T0, &
                             RestartOpts%dT, &
                             Up = .false.    &
                           )
    NextRestartWriteTime = Time + RestartOpts%dT
  Else
    NextRestartWriteTime = InfFutureTime()
  End If

End Function FindNextRestartWriteTime

!-------------------------------------------------------------------------------------------------------------

Subroutine WriteRestartFile(        &
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
!

  Implicit None
  ! Argument list:
  Type(RestartOpts_),    Intent(In)    :: RestartOpts
  Integer,               Intent(In)    :: iCase
  Type(Results_),        Intent(InOut) :: Results
  Type(DispState_),      Intent(InOut) :: DispState
  Type(DispOpts_),       Intent(In)    :: DispOpts
  Type(ChemistryState_), Intent(In)    :: ChemistryState
  Logical,               Intent(In)    :: EndOfCase
  Logical,               Intent(In)    :: EndOfCases
  Type(Units_),          Intent(InOut) :: Units
  Type(Reqs_),           Intent(In)    :: Reqs
  ! Locals:
  Integer                      :: Unit
  Integer                      :: iOutputGroup
  Integer                      :: iPPInfoReq
  Integer                      :: iFile
  Character(MaxFileNameLength) :: RestartFile
  Integer                      :: IOStat
  Integer                      :: ValueOfErrorCode

  ! Restart file name.
  RestartFile = Trim(RestartOpts%FileStem)         // &
                'Restart_C'                        // &
                Trim(Int2Char(iCase))              // &
                '_T'                               // &
                Trim(FileNameTime(DispState%Time)) // & ! Add time index$$
                '.dat'

  ! Add restart file to catalogue 3 as a candidate for deletion if run crashes or is
  ! stopped during writing the restart file.
  Call MarkCatalogueInvalid       (RestartOpts%FileStem, '3', Units)
  Call AddToCatalogue(RestartFile, RestartOpts%FileStem, '3', Units)
  Call MarkCatalogueValidForDelete(RestartOpts%FileStem, '3', Units)

  ! Open restart file.
  ! NOTE: Little-endian to big-endian conversion not supported for derived types in
  ! unformatted I/O. The restart file must therefore be written in the native endian.
  ! This can be achieved by reserving a unit number for which conversion is not done.
  ! On the Linux system, endian conversion is not performed on logical unit N by
  ! setting the environment variable F_UFMTENDIAN = "big;little:N"
  ! Call GetNewUnit(Unit, Units)
  !$$ Should this be made compiler specific??
  ! $$ use openfile in due course once convert sorted
  Unit = 100
  Open(                                          &
    Unit   = Unit,                               &
    File   = Trim(ConvertFileName(RestartFile)), &
    Form   = 'Unformatted',                      &
    IOStat = IOStat                              &
  )
  If (IOStat /= 0) Then
    Call Message(                                       &
           'FATAL ERROR: Cannot open restart file "' // &
           Trim(RestartFile)                         // &
           '"',                                         &
           3                                            &
         )
  End If

  ! Write restart file.
  ValueOfErrorCode = GetErrorCode()
  Write (Unit) ValueOfErrorCode
  Write (Unit) iCase
  Write (Unit) EndOfCase
  Write (Unit) EndOfCases
  Call WriteRandomState(Unit, GlobalRandomState)
  Call WriteRestartFileResults(Unit, Results)
  Call WriteRestartFileDispState(Unit, DispState)

  Write (Unit) DispOpts%Chemistry
  If (DispOpts%Chemistry) Call WriteRestartFileChemistry(Unit, ChemistryState)

  ! Close restart file.
  !Call CloseUnit(Unit, Units)
  Close (Unit)

  ! Flush buffers of files currently active for field output.
  ! $$ Use nFieldGroups and nFiles instead here?
  ! $$ Could move to WriteRestartFileResults?
  Do iOutputGroup = 1, Reqs%MaxFieldOutputGroups
    Do iFile = 1, MaxFieldOutputFiles
      If (Results%FieldDiskUnits(iOutputGroup, iFile) /= 0) Then
        Flush(Results%FieldDiskUnits(iOutputGroup, iFile))
      End If
    End Do
  End Do
  Do iPPInfoReq = 1, MaxPPInfoReqs
    Do iFile = 1, MaxPPInfoOutputFiles
      If (Results%PPInfos(iPPInfoReq)%DiskUnits(iFile) /= 0) Then
        Flush(Results%PPInfos(iPPInfoReq)%DiskUnits(iFile))
      End If
    End Do
  End Do

  If (RestartOpts%DeleteOldFiles) Then

    ! Catalogue 1: Mark as invalid, update, and mark as valid.
    Call MarkCatalogueInvalid       (RestartOpts%FileStem, '1', Units)
    Call EmptyCatalogue             (RestartOpts%FileStem, '1', Units)
    Call AddToCatalogue(RestartFile, RestartOpts%FileStem, '1', Units)
    Call MarkCatalogueValid         (RestartOpts%FileStem, '1', Units)

  Else

    ! Catalogue 1: Mark as invalid, update, and mark as valid.
    Call MarkCatalogueInvalid       (RestartOpts%FileStem, '1', Units)
    Call AddToCatalogue(RestartFile, RestartOpts%FileStem, '1', Units)
    Call MarkCatalogueValid         (RestartOpts%FileStem, '1', Units)

  End If

  ! Delete old restart files.
  Call DeleteCataloguedFiles(RestartOpts%FileStem, '1', '3', Units)

  ! Catalogue 2: Mark as invalid, copy catalogue 1 to 2, and mark as valid.
  Call MarkCatalogueInvalid(RestartOpts%FileStem,      '2',      Units)
  Call CopyCatalogue       (RestartOpts%FileStem, '1', '2', ' ', Units)
  Call MarkCatalogueValid  (RestartOpts%FileStem,      '2',      Units)

  ! Catalogue 3: Mark as invalid, copy catalogue 1 to 3, and mark as valid.
  Call MarkCatalogueInvalid(RestartOpts%FileStem,      '3',      Units)
  Call CopyCatalogue       (RestartOpts%FileStem, '1', '3', ' ', Units)
  Call MarkCatalogueValid  (RestartOpts%FileStem,      '3',      Units)

End Subroutine WriteRestartFile

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadRestartFile(           &
             RestartFileSuffix,       &
             RestartOpts, OutputOpts, &
             nSources,                &
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
!

  Implicit None
  ! Argument list:
  Character(MaxFileNameLength), Intent(In)    :: RestartFileSuffix
  Type(RestartOpts_),           Intent(In)    :: RestartOpts
  Type(OutputOpts_),            Intent(In)    :: OutputOpts ! $$ replace with restart options
  Integer,                      Intent(In)    :: nSources
  Logical,                      Intent(Out)   :: NewSourcesDetected
  Integer,                      Intent(Out)   :: iFirstNewSource
  Integer,                      Intent(InOut) :: iCase
  Type(Results_),               Intent(InOut) :: Results
  Type(DispState_),             Intent(InOut) :: DispState
  Type(ChemistryState_),        Intent(InOut) :: ChemistryState
  Logical,                      Intent(InOut) :: EndOfCase
  Logical,                      Intent(InOut) :: EndOfCases
  Logical,                      Intent(InOut) :: Restart
  Type(Units_),                 Intent(InOut) :: Units
  ! nSources           :: Number of sources defined in input file.
  ! NewSourcesDetected :: Indicates that one or more new sources have been detected in the input file
  !                       during a restart run.
  ! iFirstNewSource    :: Index of first new source in a modified input file for a restart run.
  
  ! Locals:
  Logical            :: Exist
  Integer            :: Unit
  Integer            :: IOStat
  Logical            :: IsChemistry
  Character(MaxFileNameLength) :: FileName
  Integer                      :: ValueOfErrorCode

  ! Initialise flag for detecting new sources in the input file.
  NewSourcesDetected = .false.
  iFirstNewSource      = -1
  
  ! Get name of restart file and find out if its in the catalogue and exists.
  ! Note that if the file isn't the latest file, its necessary to delete the catalogue
  ! beyond FileName
  ! because of output files. Once the run has been restarted the output files will
  ! reflect this and be too short to restart from a later restart file.
  If (RestartFileSuffix == ' ') Then

    Call LastFileInCatalogue(RestartOpts%FileStem, '1', FileName, Units)
    If (FileName == ' ') Then
      Call Message(' ')
      Call Message(                                                  &
             'No valid restart file found in restart catalogue. ' // &
             'Run restarting at the begining.'                       &
           )
      Restart = .false.
      Return
    End If
    Inquire(File = FileName, Exist = Exist)
    If (.not. Exist) Then
      Call Message(' ')
      Call Message(                          &
             'FATAL ERROR: Restart file ' // &
             Trim(FileName)               // &
             ' not found.',                  &
             3                               &
           )
    End If

  Else

    FileName = Trim(RestartOpts%FileStem) // &
               'Restart_'                 // &
               Trim(RestartFileSuffix)    // &
               '.dat'
    Call InCatalogue(RestartOpts%FileStem, '1', FileName, Exist, Units)
    If (.not. Exist) Then
      Call Message(' ')
      Call Message(                          &
             'FATAL ERROR: Restart file ' // &
             Trim(FileName)               // &
             ' not in catalogue.',           &
             3                               &
           )
    End If
    Inquire(File = FileName, Exist = Exist)
    If (.not. Exist) Then
      Call Message(' ')
      Call Message(                          &
             'FATAL ERROR: Restart file ' // &
             Trim(FileName)               // &
             ' not found.',                  &
             3                               &
           )
    End If

    ! Catalogue 2: Mark as invalid, copy catalogue 1 to 2, and mark as valid.
    Call MarkCatalogueInvalid(RestartOpts%FileStem,      '2',      Units)
    Call CopyCatalogue       (RestartOpts%FileStem, '1', '2', ' ', Units)
    Call MarkCatalogueValid  (RestartOpts%FileStem,      '2',      Units)

    ! Catalogue 1: Mark as invalid, copy catalogue 2 to 1 up to FileName, and mark as
    ! valid.
    Call MarkCatalogueInvalid(RestartOpts%FileStem,      '1',           Units)
    Call CopyCatalogue       (RestartOpts%FileStem, '2', '1', FileName, Units)
    Call MarkCatalogueValid  (RestartOpts%FileStem,      '1',           Units)

    Call DeleteCataloguedFiles(RestartOpts%FileStem, '1', '2', Units)

    ! Catalogue 2: Mark as invalid, copy catalogue 1 to 2, and mark as valid.
    Call MarkCatalogueInvalid(RestartOpts%FileStem,      '2',      Units)
    Call CopyCatalogue       (RestartOpts%FileStem, '1', '2', ' ', Units)
    Call MarkCatalogueValid  (RestartOpts%FileStem,      '2',      Units)

    ! Catalogue 3: Mark as invalid, copy catalogue 1 to 3, and mark as valid.
    Call MarkCatalogueInvalid(RestartOpts%FileStem,      '3',      Units)
    Call CopyCatalogue       (RestartOpts%FileStem, '1', '3', ' ', Units)
    Call MarkCatalogueValid  (RestartOpts%FileStem,      '3',      Units)

  End If

  ! Open restart file.
  ! NOTE: Little-endian to big-endian conversion not supported for derived types in
  ! unformatted I/O. The restart file must therefore be written in the native endian.
  ! This can be achieved by reserving a unit number for which conversion is not done.
  ! Call GetNewUnit(Unit, Units)
  !$$ Should this be made compiler specific??
  ! $$ use openfile in due course once convert sorted
  Unit = 100
  Open(                                       &
    Unit   = Unit,                            &
    File   = Trim(ConvertFileName(FileName)), &
    Status = 'Old',                           &
    Form   = 'Unformatted',                   &
    Action = 'Read',                          &
    IOStat = IOStat                           &
  )
  If (IOStat /= 0) Then
    Call Message(                                       &
           'FATAL ERROR: Cannot open restart file "' // &
           Trim(FileName)                            // &
           '"',                                         &
           3                                            &
         )
  End If

  ! Read restart file.
  Read (Unit) ValueOfErrorCode
  Call SetMinErrorCode(ValueOfErrorCode)
  Read (Unit) iCase
  Read (Unit) EndOfCase
  Read (Unit) EndOfCases
  Call ReadRandomState(Unit, GlobalRandomState)
  Call ReadRestartFileResults(Unit, Results)
  Call ReadRestartFileDispState(nSources, iFirstNewSource, NewSourcesDetected, Unit, DispState)

  Read (Unit) IsChemistry
  If (IsChemistry) Call ReadRestartFileChemistry(Unit, ChemistryState)

  ! Close restart file.
  !Call CloseUnit(Unit, Units)
  Close (Unit)

  If (EndOfCase) Then
    Call Message(' ')
    Call Message('Restart file '            // &
           Trim(FileName)                   // &
           ' read succesfully. '            // &
           'Run restarting at end of case ' // &
           Trim(Int2Char(iCase))            // &
           '.'                                 &
         )
  Else
    Call Message(' ')
    Call Message(                           &
           'Restart file '               // &
           Trim(FileName)                // &
           ' read succesfully. '         // &
           'Run restarting with case '   // &
           Trim(Int2Char(iCase))         // &
           ' at '                        // &
           Trim(                            &
             Time2Char(                     &
               DispState%Time,              &
               OutputOpts%Seconds,          &
               OutputOpts%DecimalPlaces,    &
               .true.                       &
             )                              &
           )                             // &
           '.'                              &
         )
    Call Message(' ')
  End If

End Subroutine ReadRestartFile

!-------------------------------------------------------------------------------------------------------------

Function IsCatalogueValid(FileStem, Index, Units) Result(Valid)
!

  Implicit None
  ! Argument list:
  Character(*), Intent(In)    :: FileStem
  Character(1), Intent(In)    :: Index
  Type(Units_), Intent(InOut) :: Units
  ! Function result:
  Logical :: Valid
  ! Locals:
  Character(MaxFileNameLength) :: File
  Integer                      :: Unit
  Character(MaxCharLength)     :: String
  Integer                      :: IOStat

  File = Trim(FileStem) // 'RestartCatalogueValidity' // Index // '.txt'

  Inquire(File = File, Exist = Valid)

  If (Valid) Then

    Unit = OpenFile(         &
             File   = File,  &
             Units  = Units, &
             Status = 'Old', &
             Action = 'Read' &
           )

    Read (Unit, '(A)', IOStat = IOStat) String

    If (IOStat > 0) Call Message('UNEXPECTED FATAL ERROR in IsCatalogueValid', 4)

    If (IOStat < 0 .or. String /= 'Valid') Valid = .false.

    Call CloseUnit(Unit, Units)

  End If

End Function IsCatalogueValid

!-------------------------------------------------------------------------------------------------------------

Function IsCatalogueValidForDelete(FileStem, Index, Units) Result(Valid)
!

  Implicit None
  ! Argument list:
  Character(*), Intent(In)    :: FileStem
  Character(1), Intent(In)    :: Index
  Type(Units_), Intent(InOut) :: Units
  ! Function result:
  Logical :: Valid
  ! Local parameters:
  Character(60) :: ValidForDelete = 'Valid only as a list of files ' // &
                                    'which might need to be deleted'
  ! Locals:
  Character(MaxFileNameLength) :: File
  Integer                      :: Unit
  Character(MaxCharLength)     :: String
  Integer                      :: IOStat

  File = Trim(FileStem) // 'RestartCatalogueValidity' // Index // '.txt'

  Inquire(File = File, Exist = Valid)

  If (Valid) Then

    Unit = OpenFile(         &
             File   = File,  &
             Units  = Units, &
             Status = 'Old', &
             Action = 'Read' &
           )

    Read (Unit, '(A)', IOStat = IOStat) String

    If (IOStat > 0) Call Message('UNEXPECTED FATAL ERROR in IsCatalogueValidForDelete', 4)

    If (IOStat < 0 .or. String /= ValidForDelete) Valid = .false.

    Call CloseUnit(Unit, Units)

  End If

End Function IsCatalogueValidForDelete

!-------------------------------------------------------------------------------------------------------------

Subroutine MarkCatalogueValid(FileStem, Index, Units)
!

  Implicit None
  ! Argument list:
  Character(*), Intent(In)    :: FileStem
  Character(1), Intent(In)    :: Index
  Type(Units_), Intent(InOut) :: Units
  ! Locals:
  Character(MaxFileNameLength) :: File
  Integer                      :: Unit

  File = Trim(FileStem) // 'RestartCatalogueValidity' // Index // '.txt'

  Unit = OpenFile(       &
           File  = File, &
           Units = Units &
         )

  Write (Unit, '(A)') 'Valid'

  Call CloseUnit(Unit = Unit, Units = Units)

End Subroutine MarkCatalogueValid

!-------------------------------------------------------------------------------------------------------------

Subroutine MarkCatalogueValidForDelete(FileStem, Index, Units)
! Marks catalogue as valid as a list of files which might need to be deleted.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)    :: FileStem
  Character(1), Intent(In)    :: Index
  Type(Units_), Intent(InOut) :: Units
  ! Local parameters:
  Character(60) :: ValidForDelete = 'Valid only as a list of files ' // &
                                    'which might need to be deleted'
  ! Locals:
  Character(MaxFileNameLength) :: File
  Integer                      :: Unit

  File = Trim(FileStem) // 'RestartCatalogueValidity' // Index // '.txt'

  Unit = OpenFile(       &
           File  = File, &
           Units = Units &
         )

  Write (Unit, '(A)') ValidForDelete

  Call CloseUnit(Unit = Unit, Units = Units)

End Subroutine MarkCatalogueValidForDelete

!-------------------------------------------------------------------------------------------------------------

Subroutine MarkCatalogueInvalid(FileStem, Index, Units)
!

  Implicit None
  ! Argument list:
  Character(*), Intent(In)    :: FileStem
  Character(1), Intent(In)    :: Index
  Type(Units_), Intent(InOut) :: Units
  ! Locals:
  Character(MaxFileNameLength) :: File
  Integer                      :: Unit

  File = Trim(FileStem) // 'RestartCatalogueValidity' // Index // '.txt'

  Unit = OpenFile(       &
           File  = File, &
           Units = Units &
         )

  Write (Unit, '(A)') 'Invalid'

  Call CloseUnit(Unit = Unit, Units = Units)

End Subroutine MarkCatalogueInvalid

!-------------------------------------------------------------------------------------------------------------

Subroutine EmptyCatalogue(FileStem, Index, Units)
!

  Implicit None
  ! Argument list:
  Character(*), Intent(In)    :: FileStem
  Character(1), Intent(In)    :: Index
  Type(Units_), Intent(InOut) :: Units
  ! Locals:
  Character(MaxFileNameLength) :: File
  Integer                      :: Unit

  File = Trim(FileStem) // 'RestartCatalogue' // Index // '.txt'

  Unit = OpenFile(       &
           File  = File, &
           Units = Units &
         )
  Call CloseUnit(Unit = Unit, Units = Units, DeleteFile = .true.)

  Unit = OpenFile(         &
           File   = File,  &
           Units  = Units, &
           Status = 'New'  &
         )
  Call CloseUnit(Unit = Unit, Units = Units)

End Subroutine EmptyCatalogue

!-------------------------------------------------------------------------------------------------------------

Subroutine AddToCatalogue(String, FileStem, Index, Units)
!

  Implicit None
  ! Argument list:
  Character(*), Intent(In)    :: String
  Character(*), Intent(In)    :: FileStem
  Character(1), Intent(In)    :: Index
  Type(Units_), Intent(InOut) :: Units
  ! Locals:
  Character(MaxFileNameLength) :: File
  Integer                      :: Unit

  File = Trim(FileStem) // 'RestartCatalogue' // Index // '.txt'

  Unit = OpenFile(             &
           File     = File,    &
           Units    = Units,   &
           Position = 'Append' &
         )

  Write (Unit, '(A)') Trim(String)

  Call CloseUnit(Unit = Unit, Units = Units)

End Subroutine AddToCatalogue

!-------------------------------------------------------------------------------------------------------------

Subroutine DeleteCataloguedFiles(FileStem, ValidIndex, ValidForDeleteIndex, Units)
! Deletes files in ValidForDeleteIndex catalogue but not in ValidIndex catalogue

  Implicit None
  ! Argument list:
  Character(*), Intent(In)    :: FileStem
  Character(1), Intent(In)    :: ValidIndex
  Character(1), Intent(In)    :: ValidForDeleteIndex
  Type(Units_), Intent(InOut) :: Units
  ! Locals:
  Character(MaxFileNameLength) :: File
  Integer                      :: Unit
  Character(MaxFileNameLength) :: File1
  Integer                      :: Unit1
  Logical                      :: Exist
  Logical                      :: Found

  File = Trim(FileStem) // 'RestartCatalogue' // ValidForDeleteIndex // '.txt'

  Inquire(File = File, Exist = Exist)

  If (Exist) Then

    Unit = OpenFile(         &
             File   = File,  &
             Units  = Units, &
             Status = 'Old', &
             Action = 'Read' &
           )

    Do
      Read (Unit, '(A)', End = 5) File1
      Call InCatalogue(FileStem, ValidIndex, File1, Found, Units)
      If (.not.Found) Then
        Unit1 = OpenFile(               &
                  File  = File1,        &
                  Units = Units,        &
                  Form  = 'Unformatted' &
                )
        Call CloseUnit(Unit = Unit1, Units = Units, DeleteFile = .true.)
      End If
    End Do
5   Continue

    Call CloseUnit(Unit = Unit, Units = Units)

  End If

End Subroutine DeleteCataloguedFiles

!-------------------------------------------------------------------------------------------------------------

Subroutine CopyCatalogue(FileStem, FromIndex, ToIndex, LastFile, Units)
!

  Implicit None
  ! Argument list:
  Character(*), Intent(In)    :: FileStem
  Character(1), Intent(In)    :: FromIndex
  Character(1), Intent(In)    :: ToIndex
  Character(*), Intent(In)    :: LastFile       ! If non-blank, copies up to this file
  Type(Units_), Intent(InOut) :: Units
  ! Locals:
  Character(MaxFileNameLength) :: FromFile
  Integer                      :: FromUnit
  Character(MaxFileNameLength) :: ToFile
  Integer                      :: ToUnit
  Logical                      :: Exist
  Character(MaxFileNameLength) :: File

  FromFile = Trim(FileStem) // 'RestartCatalogue' // FromIndex // '.txt'
  ToFile   = Trim(FileStem) // 'RestartCatalogue' // ToIndex   // '.txt'

  Inquire(File = FromFile, Exist = Exist)

  If (Exist) Then

    FromUnit = OpenFile(            &
                 File   = FromFile, &
                 Units  = Units,    &
                 Status = 'Old',    &
                 Action = 'Read'    &
               )
    ToUnit   = OpenFile(          & ! $$ just wondering if status = old is correct here. Probably is,
                 File   = ToFile, & !    but needs checking and explaining in a comment.
                 Units  = Units,  &
                 Status = 'Old'   &
               )

    Do
      Read  (FromUnit, '(A)', End = 5) File
      Write (ToUnit,   '(A)'         ) File
      If (LastFile /= ' ' .and. (LastFile .CIEq. File)) Exit
    End Do
5   Continue

    Call CloseUnit(Unit = FromUnit, Units = Units)
    Call CloseUnit(Unit = ToUnit,   Units = Units)

  Else

    Call EmptyCatalogue(FileStem, ToIndex, Units)

  End If

End Subroutine CopyCatalogue

!-------------------------------------------------------------------------------------------------------------

Subroutine InCatalogue(FileStem, Index, FileName, Found, Units)
!

  Implicit None
  ! Argument list:
  Character(*), Intent(In)    :: FileStem
  Character(1), Intent(In)    :: Index
  Character(*), Intent(In)    :: FileName
  Logical,      Intent(Out)   :: Found
  Type(Units_), Intent(InOut) :: Units
  ! Locals:
  Character(MaxFileNameLength) :: File
  Integer                      :: Unit
  Logical                      :: Exist

  File = Trim(FileStem) // 'RestartCatalogue' // Index // '.txt'

  Inquire(File = File, Exist = Exist)

  Found = .false.

  If (Exist) Then

    Unit = OpenFile(         &
             File   = File,  &
             Units  = Units, &
             Status = 'Old', &
             Action = 'Read' &
           )

    Do
      Read (Unit, '(A)', End = 5) File
      If (File .CIEq. FileName) Found = .true.
    End Do
5   Continue

    Call CloseUnit(Unit, Units)

  End If

End Subroutine InCatalogue

!-------------------------------------------------------------------------------------------------------------

Subroutine LastFileInCatalogue(FileStem, Index, LastFile, Units)
!

  Implicit None
  ! Argument list:
  Character(*), Intent(In)    :: FileStem
  Character(1), Intent(In)    :: Index
  Character(*), Intent(Out)   :: LastFile   ! Blank means catalogue empty. $$ (*)
  Type(Units_), Intent(InOut) :: Units
  ! Locals:
  Character(MaxFileNameLength) :: File
  Integer                      :: Unit
  Logical                      :: Exist

  File = Trim(FileStem) // 'RestartCatalogue' // Index // '.txt'

  Inquire(File = File, Exist = Exist)

  LastFile = ' '

  If (Exist) Then

    Unit = OpenFile(         &
             File   = File,  &
             Units  = Units, &
             Status = 'Old', &
             Action = 'Read' &
           )
    Do
      Read (Unit, '(A)', End = 5) LastFile
    End Do
5   Continue
    Call CloseUnit(Unit, Units)

  End If

End Subroutine LastFileInCatalogue

!-------------------------------------------------------------------------------------------------------------

End Module RestartModule

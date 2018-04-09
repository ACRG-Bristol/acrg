! Module: Error and message Module

Module ErrorAndMessageModule

! This module provides code for handling errors and messages.

! Module overview
! ---------------

! $$
!
! Error codes are as follows:
!   0 indicates a message (i.e. not a warning or error),
!   1 indicates a warning,
!   2 indicates a non-fatal error,
!   3 indicates a fatal error,
!   4 indicates an unexpected fatal error.
! Unexpected errors are those which require the code or system to be changed to correct (i.e. not due to
! errors in the user input). Untrapped errors are not considered here but are handled by the compiler.

! Module use
! ----------

! $$
!
! $$ include comment along lines of the following:
! After call to SetUpErrorAndMessageModule any closing message box will
! include the model version, messages are written to the log file (as well as unit
! 6), and ClosePromptLevel takes effect.

! Module call tree
! ----------------

! $$

!-------------------------------------------------------------------------------------------------------------

#ifdef CompaqPCCompiler
  Use DFLib
#endif
#ifdef IntelWin
  Use DFLib
#endif
#ifdef IntelLinCompiler
  Use IFCore, Only: TracebackQQ
#endif
Use GlobalParametersModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: InitErrorAndMessageModule  ! Initialises the error and message module.
Public :: SetUpErrorAndMessageModule ! Sets up the error and message module by setting the caption for the
                                     ! message box used when the run terminates, opening the log and error
                                     ! files, and setting the error code at or above which the user is
                                     ! prompted to close the window when the run terminates.
Public :: GetErrorCode               ! Returns the worst (highest) error code to have occurred.
Public :: SetMinErrorCode            ! Imposes a minimum on the variable recording the worst (highest) error
                                     ! code to have occurred (to be used only following a restart).
Public :: Message                    ! Outputs a message to the screen and, if possible and appropriate, to
                                     ! the log and error files, updates the record of the worst (highest)
                                     ! error code to have occurred, and, if appropriate, ends or suspends the
                                     ! run.

!-------------------------------------------------------------------------------------------------------------

Character(MaxCharLength2), Private, Save :: EndOfRunCaption
Logical,                   Private, Save :: LogFileOpened
Logical,                   Private, Save :: ErrorFileOpened
Integer,                   Private, Save :: PromptErrorCode
Integer,                   Private, Save :: WorstErrorCode
! EndOfRunCaption :: Caption for the message box used when the run terminates.
! LogFileOpened   :: Indicates the log file has been opened.
! ErrorFileOpened :: Indicates the error file has been opened.
! PromptErrorCode :: Error code at or above which the user is prompted to close the window when the run
!                    terminates. Below this value the window closes automatically.
! WorstErrorCode  :: Worst (highest) error code to have occured.

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine InitErrorAndMessageModule
! Initialises the error and message module.

  Implicit None

  EndOfRunCaption = ' '
  LogFileOpened   = .false.
  ErrorFileOpened = .false.
  PromptErrorCode = 0
  WorstErrorCode  = 0

End Subroutine InitErrorAndMessageModule

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpErrorAndMessageModule(EndRunCaption, LogFile, ErrorFile, ClosePromptLevel, Restart)
! Sets up the error and message module by setting the caption for the message box used when the run
! terminates, opening the log and error files, and setting the error code at or above which the user is
! prompted to close the window when the run terminates.
!
! Note the code is arranged so that none of EndRunCaption, LogFile, ErrorFile and ClosePromptLevel take effect
! until after this routine is completed.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: EndRunCaption    ! Caption for the message box used when the run terminates.
  Character(*), Intent(In) :: LogFile          ! Log file (blank indicates no log file is to be used).
  Character(*), Intent(In) :: ErrorFile        ! Error file (blank indicates no error file is to be used).
  Integer,      Intent(In) :: ClosePromptLevel ! Error code at or above which the user is prompted to close
                                               ! the window when the run terminates. Below this value the
                                               ! window closes automatically.
  Logical,      Intent(In) :: Restart          ! Indicates the run is being restarted.
  ! Locals:
  Integer :: IOStat      ! Status code for open statement.
  Logical :: LogOpened   ! Indicates the log file has been opened.
  Logical :: ErrorOpened ! Indicates the error file has been opened.

  If (Len_Trim(EndRunCaption) > Len(EndOfRunCaption)) Then
    Call Message(                                            &
           'UNEXPECTED FATAL ERROR in SetUpErrorAndMessage', &
           4                                                 &
         )
  End If

  LogOpened = .false.
  If (LogFile /= ' ') Then
    If (Restart) Then
      Open (                    &
        Unit     = LogFileUnit, &
        File     = LogFile,     &
        Status   = 'Old',       &
        Position = 'Append',    &
        IOStat   = IOStat       &
      )
      If (IOStat /= 0) Then
        Call Message(                                                               &
               'FATAL ERROR: unable to open log file (possibly because restart ' // &
               'selected but no log file yet exists)',                              &
               3                                                                    &
             )
      End If
    Else
      Open (                  &
        Unit   = LogFileUnit, &
        File   = LogFile,     &
        Status = 'Replace',   &
        IOStat = IOStat       &
      )
      If (IOStat /= 0) Then
        Call Message('FATAL ERROR: unable to open log file', 3)
      End If
    End If
    LogOpened = .true.
  End If

  ErrorOpened = .false.
  If (ErrorFile /= ' ') Then
    If (Restart) Then
      Open (                      &
        Unit     = ErrorFileUnit, &
        File     = ErrorFile,     &
        Status   = 'Old',         &
        Position = 'Append',      &
        IOStat   = IOStat         &
      )
      If (IOStat /= 0) Then
        Call Message(                                                                 &
               'FATAL ERROR: unable to open error file (possibly because restart ' // &
               'selected but no error file yet exists)',                              &
               3                                                                      &
             )
      End If
    Else
      Open (                    &
        Unit   = ErrorFileUnit, &
        File   = ErrorFile,     &
        Status = 'Replace',     &
        IOStat = IOStat         &
      )
      If (IOStat /= 0) Then
        Call Message('FATAL ERROR: unable to open error file', 3)
      End If
    End If
    ErrorOpened = .true.
  End If

  LogFileOpened   = LogOpened

  ErrorFileOpened = ErrorOpened

  EndOfRunCaption = EndRunCaption

  PromptErrorCode = ClosePromptLevel

End Subroutine SetUpErrorAndMessageModule

!-------------------------------------------------------------------------------------------------------------

Function GetErrorCode()
! Returns the worst (highest) error code to have occurred.

  Implicit None
  ! Function result:
  Integer :: GetErrorCode ! The worst (highest) error code to have occurred.

  GetErrorCode = WorstErrorCode

End Function GetErrorCode

!-------------------------------------------------------------------------------------------------------------

Subroutine SetMinErrorCode(MinErrorCode)
! Imposes a minimum on the variable recording the worst (highest) error code to have occurred (to be used only
! following a restart).

  Implicit None
  ! Argument list:
  Integer, Intent(In) :: MinErrorCode ! The minimum value to be imposed on the variable recording the worst
                                      ! (highest) error code to have occurred.

  WorstErrorCode = Max(WorstErrorCode, MinErrorCode)

End Subroutine SetMinErrorCode

!-------------------------------------------------------------------------------------------------------------

Subroutine Message(Char, ErrorCode, EndRun)
! Outputs a message to the screen and, if possible and appropriate, to the log and error files, updates the
! record of the worst (highest) error code to have occurred, and, if appropriate, ends or suspends the run.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: Char
  Integer,      Intent(In), Optional :: ErrorCode
  Logical,      Intent(In), Optional :: EndRun
  ! Char      :: Message to output.
  ! ErrorCode :: Error code. If omitted, zero is assumed (indicating no error or warning).
  ! EndRun    :: Affects whether and how EndOfRun is called.
  !              a) If present and true, EndOfRun is called with the 'run completed' message.
  !              b) If present and false, EndOfRun is called with the 'run suspended' message.
  !              c) If absent and worst error code to date >= 3, EndOfRun is called with the 'run completed'
  !                 message.
  !              d) If absent and worst error code to date < 3, EndOfRun is not called.
  ! Locals:
  Character(MaxCharLength2) :: Msg        ! Message to be displayed via EndOfRun.
  Integer                   :: ErrorCodeL ! Local copy of ErrorCode, set to zero if ErrorCode missing.

  If (Present(ErrorCode)) Then
    ErrorCodeL = ErrorCode
  Else
    ErrorCodeL = 0
  End If

  If (ErrorCodeL < 0 .or. ErrorCodeL > 4) Then
    ErrorCodeL = 4
    Call WriteMessage('UNEXPECTED FATAL ERROR in Message', ErrorCodeL)
  Else
    Call WriteMessage(Char, ErrorCodeL)
  End If

  WorstErrorCode = Max(WorstErrorCode, ErrorCodeL)

# ifdef ExtraChecks
    If (WorstErrorCode >= 3) Then
#     ifdef IntelLinCompiler
        Call WriteMessage(                                                                    &
               'Traceback initiated by NAME (note it will not appear in log or error file).', &
               ErrorCodeL                                                                     &
             )
        Call TracebackQQ(User_Exit_Code = -1)
#     endif
    End If
# endif

  If (Present(EndRun) .or. WorstErrorCode >= 3) Then

    If (Present(EndRun)) Then
      If (EndRun) Then
        Msg = 'Run completed'
      Else
        Msg = 'Run suspended'
      End If
    Else
      Msg = 'Run completed'
    End If

    If (WorstErrorCode == 0) Then
      Msg = Trim(Msg) // ' successfully'
    Else If (WorstErrorCode == 1) Then
      Msg = Trim(Msg) // ' with some warnings'
    Else If (WorstErrorCode == 2) Then
      Msg = Trim(Msg) // ' with some errors'
    Else If (WorstErrorCode == 3) Then
      Msg = Trim(Msg) // ' with a fatal error'
    Else If (WorstErrorCode == 4) Then
      Msg = Trim(Msg) // ' with an unexpected fatal error'
    End If

    If (WorstErrorCode > 0) Then
      If (LogFileOpened .and. ErrorFileOpened) Then
        Msg = Trim(Msg) // ' (see log or error file)'
      Else If (LogFileOpened) Then
        Msg = Trim(Msg) // ' (see log file)'
      Else If (ErrorFileOpened) Then
        Msg = Trim(Msg) // ' (see error file)'
      End If
    End If

    Msg = Trim(Msg) // ': Close Window?'

    Call EndOfRun(WorstErrorCode >= PromptErrorCode, Msg, EndOfRunCaption)

  End If

End Subroutine Message

!-------------------------------------------------------------------------------------------------------------

Subroutine WriteMessage(Char, ErrorCode)
! Outputs a message to the screen and, if possible and appropriate, to the log and error files.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Char      ! Message to output.
  Integer,      Intent(In) :: ErrorCode ! Error code.
  ! Locals:
  Integer :: Length ! Trimmed length of Char.
  Integer :: i      !} Start and end locations of the part of Char to be output on a given line.
  Integer :: j      !}

  Length = Len_Trim(Char)

  If (Length == 0) Then

    Write (6, *)
    If (LogFileOpened) Write (LogFileUnit, *)
    If (ErrorFileOpened .and. ErrorCode >= 1) Write (ErrorFileUnit, *)

  Else

    j = 0

    Do

      i = Verify(Char(j + 1: ), ' ') + j
      If (i + MessageLineLength > Length) Then
        j = Length
      Else
        j = Scan(Char(i:i + MessageLineLength), ' ', .true.)
        If (j == 0) j = MessageLineLength + 1
        j = i + j - 2
      End If

      Write (6, *) Char(i:j)
      If (LogFileOpened) Write (LogFileUnit, *) Char(i:j)
      If (ErrorFileOpened .and. ErrorCode >= 1) Write (ErrorFileUnit, *) Char(i:j)

      If (j == Length) Exit

    End Do

  End If

End Subroutine WriteMessage

!-------------------------------------------------------------------------------------------------------------

Subroutine EndOfRun(Prompt, Msg, Caption)
! Terminates the run, possibly using a message box.
!
! Note that, in order to avoid possible recursion, WriteMessage is called instead of Message when an error
! occurs.

  Implicit None
  ! Argument list:
  Logical,      Intent(In) :: Prompt  ! Indicates the message and caption are to be displayed and the user
                                      ! prompted to close the window. Otherwise the window is closed without
                                      ! first prompting the user.
  Character(*), Intent(In) :: Msg     ! Message. This should include the question 'close window?' if Prompt is
                                      ! true.
  Character(*), Intent(In) :: Caption ! Message box caption.
  ! Locals:
  Integer :: Status ! Results of MsgBox and SetExitQQ.

# ifdef CompaqPCCompiler

    If (Prompt) Then

      Status = MessageBoxQQ(                 &
                 Trim(Msg) // Char(0),       &
                 Trim(Caption) // Char(0),   &
                 MB$YesNo .or. MB$DefButton2 &
               )
      If (Status == 0) Then
        Call WriteMessage(                                                                           &
               'UNEXPECTED FATAL ERROR in EndOfRun: insufficient memory to display the message box', &
               4                                                                                     &
             )
        Stop
      End If

      If (Status == MB$IDYes) Then
        Status = SetExitQQ(QWin$ExitNoPersist)
      Else If (Status == MB$IDNo) Then
        Status = SetExitQQ(QWin$ExitPersist)
      End If
      If (Status /= 0) Then
        Call WriteMessage('UNEXPECTED FATAL ERROR in EndOfRun: cannot set exit behaviour', 4)
        Stop
      End If

    Else

      Status = SetExitQQ(QWin$ExitNoPersist)
      If (Status /= 0) Then
        Call WriteMessage('UNEXPECTED FATAL ERROR in EndOfRun: cannot set exit behaviour', 4)
        Stop
      End If

    End If

# endif

# ifdef IntelWin

    If (Prompt) Then

      Status = MessageBoxQQ(                 &
                 Trim(Msg) // Char(0),       &
                 Trim(Caption) // Char(0),   &
                 MB$YesNo .or. MB$DefButton2 &
               )
      If (Status == 0) Then
        Call WriteMessage(                                                                           &
               'UNEXPECTED FATAL ERROR in EndOfRun: insufficient memory to display the message box', &
               4                                                                                     &
             )
        Stop
      End If

      If (Status == MB$IDYes) Then
        Status = SetExitQQ(QWin$ExitNoPersist)
      Else If (Status == MB$IDNo) Then
        Status = SetExitQQ(QWin$ExitPersist)
      End If
      If (Status /= 0) Then
        Call WriteMessage('UNEXPECTED FATAL ERROR in EndOfRun: cannot set exit behaviour', 4)
        Stop
      End If

    Else

      Status = SetExitQQ(QWin$ExitNoPersist)
      If (Status /= 0) Then
        Call WriteMessage('UNEXPECTED FATAL ERROR in EndOfRun: cannot set exit behaviour', 4)
        Stop
      End If

    End If

# endif

  Stop

End Subroutine EndOfRun

!-------------------------------------------------------------------------------------------------------------

End Module ErrorAndMessageModule

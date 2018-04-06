! Module:  System Module

Module SystemModule

! This module provides code for interactions with the operating system and presents a consistent platform- and
! compiler-independent interface for such interactions to the rest of the code.

!-------------------------------------------------------------------------------------------------------------

#ifdef IntelLinCompiler
  Use IFLPort
#endif
Use ErrorAndMessageModule
Use StringModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: GetCommandLineArguments ! Gets command line arguments.
Public :: SubmitSystemCommand     ! Submits a system command.

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine GetCommandLineArguments(nArguments, Arguments)
! Gets command line arguments.

  Implicit None
  ! Argument list:
  Integer,      Intent(Out) :: nArguments   ! Number of command line arguments.
  Character(*), Intent(Out) :: Arguments(:) ! Command line arguments.
  ! Locals:
  Integer :: iArgument      ! Index for command line arguments.
  Integer :: ArgumentLength ! Length of a command line argument.
  Integer :: Status         ! Status code.
  
  nArguments   = 0
  Arguments(:) = ' '

  nArguments = Command_Argument_Count()
  If (nArguments > Size(Arguments)) Then
    Call Message(                                              &
           'FATAL ERROR: '                                  // &
           'number of command line arguments, '             // &
           Trim(Int2Char(nArguments))                       // &
           ', is greater than the maximum number allowed, ' // &
           Trim(Int2Char(Size(Arguments)))                  // &
           '.',                                                &
           3                                                   &
         )
  End If
  Do iArgument = 1, nArguments
    Call Get_Command_Argument(iArgument, Arguments(iArgument), ArgumentLength, Status) 
    If (Status == -1) Then
      Call Message(                        &
             'FATAL ERROR: '            // &
             'command line argument "'  // &
             Trim(Arguments(iArgument)) // &
             '..." is too long.',          &
             3                             &
           )
    Else If (Status /= 0) Then
      Call Message(                                                   &
             'UNEXPECTED FATAL ERROR in GetCommandLineArguments: ' // &
             'error in retrieving the '                            // &
             Trim(Int2Char(iArgument))                             // &
             'th command line argument',                              &
             4                                                        &
           )
    End If
    ! $$ probably unnecessary as a truncated argument should be picked up in the test above
    If (ArgumentLength > Len(Arguments)) Then
      Call Message(                        &
             'FATAL ERROR: '            // &
             'command line argument "'  // &
             Trim(Arguments(iArgument)) // &
             '..." is too long.',          &
             3                             &
           )
    End If
  End Do

End Subroutine GetCommandLineArguments

!-------------------------------------------------------------------------------------------------------------

Subroutine SubmitSystemCommand(SystemCommand)
! Submits a system command.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: SystemCommand ! Command to be submitted.
  ! Locals:
  Integer :: ErrorCode  ! Error code.

# ifdef IntelLinCompiler
    Integer :: ErrNum   ! Error handling.
# endif

  Call Message('System command submitted: ' // Trim(SystemCommand))

# ifdef IntelLinCompiler

    ! Retained portability function for current version of Intel compiler,
    ! but the intrinsic subroutine 'Execute_Command_Line' is available in v15.
    
    ErrorCode = System(SystemCommand)
    If (ErrorCode == -1) Then
      ErrNum = ierrno( )
      ! $$ Check for other errors
      If (ErrNum == ENOMEM) Then
        Call Message(                                                      &
               'FATAL ERROR in SubmitSystemCommand: '                   // &
               'unable to submit system command - insufficient memory',    &
               3                                                           &
             )
      End If
    End If

# else

  Call Execute_Command_Line(Trim(SystemCommand), ExitStat=ErrorCode)
  
  If (ErrorCode /= 0) Then
    ! $$ Check for specific error codes here
    Call Message(                                                                                         &
           'FATAL ERROR in SubmitSystemCommand: '                                                      // &
           'unable to submit system command or an error has occurred while executing the system command', &
           3                                                                                              &
         )
  End If

# endif

End Subroutine SubmitSystemCommand

!-------------------------------------------------------------------------------------------------------------

End Module SystemModule

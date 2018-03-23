! Module: Error and message II Module

Module ErrorAndMessageIIModule

! This module provides code for handling errors and messages in a slightly more sophisticated way than is
! available through the error and message module.

! Module overview
! ---------------

! $$

! Module use
! ----------

! To use the module one needs a variable of type MessageControls_ for storing a collection of message
! controls. This could be either the variable GlobalMessageControls which is made available by this module or
! a locally declared variable. InitMessageControls must be used to initialise this variable.
!
! Then the required message controls of type MessageControl_ must be initialised using InitMessageControl and
! stored in the collection of message controls using AddMessageControl.
!
! If it is desired to check that the message controls all have appropriate names (e.g. to trap for misspelling
! when the names are part of the user input), AddValidNameToMessageControls must be used to add the valid
! names to the collection of message controls.
!
! If it is desired to check whether a certain message control is present in the collection of message controls
! this can be done using FindMessageControlIndex.
!
! If it is desired to check whether a certain message control's properties are appropriate, this can be done
! using GetMessageControlProperties.
!
! These three types of check can be done in any order. The first type can be done before, during or after
! initialising and adding the message controls (but after initialsing the collection of message controls).
! Note however the check is not actually done until SetupMessageControls is called (see below). The other two
! types can be done after adding the message control in question.
!
! Then SetupMessageControls should be called.
!
! To output a message under the control of a message control, call ControlledMessage. If the required message
! control is not present or SetupMessageControls has not yet been called, then the default message action is
! followed. In calling ControlledMessage one can either pass the name of the message control or the index of
! the message control (which can be found using FindMessageControlIndex if the required message control exists
! and which should be set to zero if the required message control doesn't exist or if it is to be ignored).
! Note messages will appear if ControlledMessage is invoked before SetupMessageControls is called, but the
! default message action is followed and the message will not be counted in deciding when to suppress messages
! or generate a fatal error.

! Module call tree
! ----------------

! InitMessageControls
!
!                      (ErrorAndMessage.Message
! InitMessageControl---(String.CaseInsensitiveEq(.CIEq.)
!                      (HeadedFile.TokenLengthTest
!
!                     (MessageControlEq---String.CaseInsensitiveEq(.CIEq.)
! AddMessageControl---(ErrorAndMessage.Message
!                     (String.CaseInsensitiveEq(.CIEq.)
!
! FindMessageControlIndex---(ErrorAndMessage.Message
!                           (String.CaseInsensitiveEq(.CIEq.)
!
! AddValidNameToMessageControls---(ErrorAndMessage.Message
!                                 (String.CaseInsensitiveEq(.CIEq.)
!
! GetMessageControlProperties---FindMessageControlIndex
!
! SetupMessageControls---(ErrorAndMessage.Message
!                        (String.CaseInsensitiveEq(.CIEq.)
!
! ControlledMessage---(FindMessageControlIndex
!                     (ErrorAndMessage.Message
!
! Calls to routines in other modules are included (indicated by ModuleName.RoutineName), but subsequent calls
! from the called routine are not.

!-------------------------------------------------------------------------------------------------------------

Use GlobalParametersModule
Use ErrorAndMessageModule
Use StringModule
Use HeadedFileModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: MessageControl_               ! A message control.
Public :: MessageControls_              ! A collection of message controls.
Public :: GlobalMessageControls         ! A globally available collection of message controls.
Public :: InitMessageControls           ! Initialises a collection of message controls.
Public :: InitMessageControl            ! Initialises a message control.
Public :: AddMessageControl             ! Adds a message control to a collection of message controls.
Public :: FindMessageControlIndex       ! Finds the index of a message control.
Public :: AddValidNameToMessageControls ! Adds a valid name to a collection of message controls.
Public :: GetMessageControlProperties   ! Gets message control properties.
Public :: SetupMessageControls          ! Sets up MessageControls.
Public :: ControlledMessage             ! Outputs a message to the screen and, if possible and appropriate, to
                                        ! the log and error files, updates the record of the worst (highest)
                                        ! error code to have occurred, and, if appropriate, ends or suspends
                                        ! the run, all under the control of a message control.

!-------------------------------------------------------------------------------------------------------------

Type :: MessageControl_ ! A message control.
  Character(MaxCharLength) :: Name     ! Name of message control.
  Character(MaxCharLength) :: Action   ! Action required for the message.
  Integer                  :: iAction  ! Code for action required.
  Integer                  :: Fail     ! After this number of messages the message will be escalated to a
                                       ! fatal error. A huge value indicates no escalation.
  Integer                  :: Suppress ! After this number of messages the messages will be suppressed. A huge
                                       ! value indicates no suppression. Fatal error messages will not be
                                       ! suppressed however, nor will the first occurrance of non-fatal errors
                                       ! and warnings, nor will the first occurance of a message if this
                                       ! message may be later escalated to a fatal error. .

  ! The possible actions and codes are as follows (where the actions are case insensitive):
  !   Action                  Code  Action taken when message invoked
  !   -------------------------------------------------------------
  !   Default or blank        -1    Same as if no message control was in place.
  !   Message                  0    } Message issued as a message, warning, error, fatal error or unexpected
  !   Warning                  1    } fatal error.
  !   Error                    2    }
  !   Fatal Error              3    }
  !   Unexpected Fatal Error   4    }

End Type

!-------------------------------------------------------------------------------------------------------------

Type :: MessageControls_ ! A collection of message controls.
  Integer                  :: nMessageControls
  Type(MessageControl_)    :: MessageControls(MaxMessageControls)
  Integer                  :: nValidNames
  Character(MaxCharLength) :: ValidNames(MaxMessageControls)
  Logical                  :: SetupComplete
  Integer                  :: Defaults(MaxMessageControls)
  Integer                  :: nMessages(MaxMessageControls)
  ! nMessageControls :: Number of message controls.
  ! MessageControls  :: Collection of message controls.
  ! nValidNames      :: Number of valid names for message controls.
  ! ValidNames       :: Valid names for message controls.
  ! SetupComplete    :: Indicates SetupMessageControls has been called.
  ! Defaults         :: If a message control specifies default action, the default action is stored here on
  !                     the first use of the message control.
  ! nMessages        :: Number of messages generated under the control of message controls.
End Type

!-------------------------------------------------------------------------------------------------------------

Type(MessageControls_), Save :: GlobalMessageControls ! A globally available collection of message controls.

! $$ add interactive user decision option
! $$ write evolving stuff to restart file?
! $$ output of error stats?
! $$ to allow internal generation of message controls, add 'BlockKey' to control the error messages in
!    InitMessageControl etc.
! $$ ability to set default action in ControlledMessage to no message.

!-------------------------------------------------------------------------------------------------------------

Interface Operator(==) ! Equality of message controls.
  Module Procedure MessageControlEq
End Interface

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitMessageControls() Result (MessageControls)
! Initialises a collection of message controls.

  Implicit None
  ! Function result:
  Type(MessageControls_) :: MessageControls ! Initialised collection of message controls.

  MessageControls%nMessageControls = 0
  MessageControls%nValidNames      = 0
  MessageControls%SetupComplete    = .false.
  MessageControls%nMessages(:)     = 0

End Function InitMessageControls

!-------------------------------------------------------------------------------------------------------------

Function InitMessageControl(Name, Action, NoFail, Fail, NoSuppress, Suppress) Result (MessageControl)
! Initialises a message control.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name       ! Name of message control.
  Character(*), Intent(In) :: Action     ! Action required for the message.
  Logical,      Intent(In) :: NoFail     ! Indicates that the message is never escalated to a fatal error.
  Integer,      Intent(In) :: Fail       ! After this number of messages the message will be escalated to a
                                         ! fatal error. Value not significant if NoFail is true.
  Logical,      Intent(In) :: NoSuppress ! Indicates that the message is never suppressed.
  Integer,      Intent(In) :: Suppress   ! After this number of messages the messages will be suppressed.
                                         ! Value not significant if NoSuppress is true.
  ! Function result:
  Type (MessageControl_) :: MessageControl ! Initialised message control.

  Call TokenLengthTest(Name, MaxCharLength, .true., 'Message Controls', ' ', 'Name')
  MessageControl%Name = Name

  Call TokenLengthTest(Name, MaxCharLength, .false., 'Message Controls', Name, 'Action')
  MessageControl%Action = Action

  If (Action == ' ' .or. (Action .CIEq. 'Default')) Then
    MessageControl%iAction = -1
  Else If (Action .CIEq. 'Message') Then
    MessageControl%iAction = 0
  Else If (Action .CIEq. 'Warning') Then
    MessageControl%iAction = 1
  Else If (Action .CIEq. 'Error') Then
    MessageControl%iAction = 2
  Else If (Action .CIEq. 'Fatal Error') Then
    MessageControl%iAction = 3
  Else If (Action .CIEq. 'Unexpected Fatal Error') Then
    MessageControl%iAction = 4
  Else
    Call Message(                                                   &
           'FATAL ERROR: the action given for message control "' // &
           Trim(Name)                                            // &
           '" is not understood.',                                  &
           3                                                        &
         )
  End If

  If (NoFail) Then
    MessageControl%Fail = Huge(Fail)
  Else
    If (Fail < 0) Then
      Call Message(                                                                &
             'FATAL ERROR: the message control "'                               // &
             Trim(Name)                                                         // &
             '" specifies a negative value for "# of Messages Before Failure"',    &
             3                                                                     &
           )
    End If
    MessageControl%Fail = Fail
  End If

  If (NoSuppress) Then
    MessageControl%Suppress = Huge(Suppress)
  Else
    If (Suppress < 0) Then
      Call Message(                                                                    &
             'FATAL ERROR: the message control "'                                   // &
             Trim(Name)                                                             // &
             '" specifies a negative value for "# of Messages Before Suppression"',    &
             3                                                                         &
           )
    End If
    MessageControl%Suppress = Suppress
  End If

End Function InitMessageControl

!-------------------------------------------------------------------------------------------------------------

Subroutine AddMessageControl(MessageControl, MessageControls)
! Adds a message control to a collection of message controls.

  Implicit None
  ! Argument list:
  Type(MessageControl_),  Intent(In)    :: MessageControl  ! Message control to be added.
  Type(MessageControls_), Intent(InOut) :: MessageControls ! Collection of message controls.
  ! Locals:
  Integer :: i ! Loop index.

  If (MessageControls%SetupComplete) Then
    Call Message('UNEXPECTED FATAL ERROR in AddValidNameToMessageControls', 4)
  End If

  Do i = 1, MessageControls%nMessageControls
    If (MessageControl%Name .CIEq. MessageControls%MessageControls(i)%Name) Then
      If (MessageControl == MessageControls%MessageControls(i)) Then
        Return
      Else
        Call Message(                                                                 &
               'FATAL ERROR in adding the message control "'                       // &
               Trim(MessageControl%Name)                                           // &
               '": a different message control with the same name already exists',    &
               3                                                                      &
             )
      End If
    End If
  End Do

  If (MessageControls%nMessageControls >= MaxMessageControls) Then
    Call Message(                                           &
           'FATAL ERROR in adding the message control "' // &
           Trim(MessageControl%Name)                     // &
           '": there are too many message controls',        &
           3                                                &
         )
  End If

  MessageControls%nMessageControls                                  = MessageControls%nMessageControls + 1
  MessageControls%MessageControls(MessageControls%nMessageControls) = MessageControl

End Subroutine AddMessageControl

!-------------------------------------------------------------------------------------------------------------

Function FindMessageControlIndex(Name, MessageControls, Error)
! Finds the index of a message control.

  Implicit None
  ! Argument list:
  Character(*),           Intent(In)            :: Name            ! Name of message control.
  Type(MessageControls_), Intent(In)            :: MessageControls ! Collection of message controls.
  Logical,                Intent(Out), Optional :: Error           ! Error flag for message control not found.
  ! Function result:
  Integer :: FindMessageControlIndex ! Index of the message control.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, MessageControls%nMessageControls
    If (Name .CIEq. MessageControls%MessageControls(i)%Name) Then
      FindMessageControlIndex = i
      If (Present(Error)) Error = .false.
      Return
    End If
  End Do

  If (Present(Error)) Then
    FindMessageControlIndex = 0
    Error = .true.
  Else
    Call Message(                              &
           'FATAL ERROR: message control "' // &
           Trim(Name)                       // &
           '" not found',                      &
           3                                   &
         )
  End If

End Function FindMessageControlIndex

!-------------------------------------------------------------------------------------------------------------

Function MessageControlEq(MessageControl1, MessageControl2)
! Tests for equality of message controls.

  Implicit None
  ! Argument list:
  Type(MessageControl_), Intent(In) :: MessageControl1 !} The two message controls.
  Type(MessageControl_), Intent(In) :: MessageControl2 !}
  ! Function result:
  Logical :: MessageControlEq ! Indicates if message controls are equal.

  MessageControlEq = (MessageControl1%Name     .CIEq. MessageControl2%Name)   .and. &
                     (MessageControl1%Action   .CIEq. MessageControl2%Action) .and. &
                      MessageControl1%Suppress   ==   MessageControl2%Suppress

End Function MessageControlEq

!-------------------------------------------------------------------------------------------------------------

Subroutine AddValidNameToMessageControls(ValidName, CheckUniqueNames, MessageControls)
! Adds a valid name to a collection of message controls.

  Implicit None
  ! Argument list:
  Character(*),           Intent(In)    :: ValidName        ! Valid name.
  Logical,                Intent(In)    :: CheckUniqueNames ! Indicates the added names are checked for
                                                            ! uniqueness.
  Type(MessageControls_), Intent(InOut) :: MessageControls  ! Collection of message controls.
  ! Locals:
  Integer :: i ! Loop index.

  If (MessageControls%SetupComplete) Then
    Call Message('UNEXPECTED FATAL ERROR in AddValidNameToMessageControls', 4)
  End If

  If (ValidName == ' ') Then
    Call Message('UNEXPECTED FATAL ERROR in AddValidNameToMessageControls', 4)
  End If

  If (CheckUniqueNames) Then
    Do i = 1, MessageControls%nValidNames
      If (ValidName .CIEq. MessageControls%ValidNames(i)) Then
        Call Message('UNEXPECTED FATAL ERROR in AddValidNameToMessageControls', 4)
      End If
    End Do
  End If

  If (MessageControls%nValidNames >= MaxMessageControls) Then
    Call Message(                                                                    &
           'UNEXPECTED FATAL ERROR in AddValidNameToMessageControls: '            // &
           'Too many valid names specified for message controls. '                // &
           'This can be fixed by increasing MaxMessageControls and recompiling.',    &
           4                                                                         &
         )
  End If

  MessageControls%nValidNames                             = MessageControls%nValidNames + 1
  MessageControls%ValidNames(MessageControls%nValidNames) = ValidName

End Subroutine AddValidNameToMessageControls

!-------------------------------------------------------------------------------------------------------------

Subroutine GetMessageControlProperties(Name, MessageControls, Action, iAction, Fail, Suppress, Error)
! Gets message control properties.

  Implicit None
  ! Argument list:
  Character(*),             Intent(In)            :: Name            ! Name of message control.
  Type(MessageControls_),   Intent(In)            :: MessageControls ! Collection of message controls.
  Character(MaxCharLength), Intent(Out)           :: Action          ! Action required for the message.
  Integer,                  Intent(Out)           :: iAction         ! Code for action required for the
                                                                     ! message.
  Integer,                  Intent(Out)           :: Fail            ! After this number of messages the
                                                                     ! message will be escalated to a fatal
                                                                     ! error. A huge value indicates no
                                                                     ! escalation.
  Integer,                  Intent(Out)           :: Suppress        ! After this number of messages the
                                                                     ! messages will be suppressed. A huge
                                                                     ! value indicates no suppression.
  Logical,                  Intent(Out), Optional :: Error           ! Error flag for message control not
                                                                     ! found.
  ! Locals:
  Integer :: i ! Index of message control.

  i = FindMessageControlIndex(Name, MessageControls, Error)

  If (Present(Error)) Then
    If (Error) Return
  End If

  Action   = MessageControls%MessageControls(i)%Action
  iAction  = MessageControls%MessageControls(i)%iAction
  Fail     = MessageControls%MessageControls(i)%Fail
  Suppress = MessageControls%MessageControls(i)%Suppress

End Subroutine GetMessageControlProperties

!-------------------------------------------------------------------------------------------------------------

Subroutine SetupMessageControls(CheckNamesValid, MessageControls)
! Sets up MessageControls.

  Implicit None
  ! Argument list:
  Logical,                Intent(In)    :: CheckNamesValid ! Indicates the validity of the message control
                                                           ! names will be checked
  Type(MessageControls_), Intent(InOut) :: MessageControls ! Collection of message controls.
  ! Locals:
  Integer :: i !} Loop indices.
  Integer :: j !}

  If (CheckNamesValid) Then
    Do i = 1, MessageControls%nMessageControls
      Do j = 1, MessageControls%nValidNames
        If (MessageControls%MessageControls(i)%Name .CIEq. MessageControls%ValidNames(j)) Exit
        If (j == MessageControls%nValidNames) Then
          Call Message(                                                                                     &
                 'FATAL ERROR: The message control "'                                                    // &
                 Trim(MessageControls%MessageControls(i)%Name)                                           // &
                 '" (which is likely to have been specified in the input but might have been generated ' // &
                 'internally) has an invalid name.',                                                        &
                 3                                                                                          &
               )
        End If
      End Do
    End Do
  End If

  MessageControls%SetupComplete = .true.

End Subroutine SetupMessageControls

!-------------------------------------------------------------------------------------------------------------

Subroutine ControlledMessage(                                      &
             Char,                                                 &
             MessageControls, MessageControlName, iMessageControl, &
             ErrorCode, EndRun                                     &
           )
! Outputs a message to the screen and, if possible and appropriate, to the log and error files, updates the
! record of the worst (highest) error code to have occurred, and, if appropriate, ends or suspends the run,
! all under the control of a message control.

  Implicit None
  ! Argument list:
  Character(*),           Intent(In)             :: Char
  Type(MessageControls_), Intent(InOut)          :: MessageControls
  Character(*),           Intent(In),   Optional :: MessageControlName
  Integer,                Intent(In),   Optional :: iMessageControl
  Integer,                Intent(In),   Optional :: ErrorCode
  Logical,                Intent(In),   Optional :: EndRun
  ! Char               :: Message to output.
  ! MessageControls    :: Collection of message controls.
  ! MessageControlName :} Name and index of message control. Only one of these should be specified.
  ! iMessageControl    :}
  ! ErrorCode          :: Error code specifying the default action to be taken. The values correspond to the
  !                       error codes used by the error and message module and to the action codes discussed
  !                       above, with the exception that -1 is not allowed. If omitted, zero is assumed
  !                       (indicating no error or warning). ErrorCode should be the same for all calls to this
  !                       routine using the same message control.
  ! EndRun             :: Affects whether and how EndOfRun is called.
  !                       a) If present and true, EndOfRun is called with the 'run completed' message.
  !                       b) If present and false, EndOfRun is called with the 'run suspended' message.
  !                       c) If absent and worst error code to date >= 3, EndOfRun is called with the 'run
  !                          completed' message.
  !                       d) If absent and worst error code to date < 3, EndOfRun is not called.

  ! Locals:
  Integer :: iMessageControlL ! Index of message control to be used.
  Integer :: iAction          ! Code specifying the action to be taken. The values correspond to the error
                              ! codes used by the error and message module and to the action codes discussed
                              ! above, with the exception that -1 is not allowed.
  Logical :: Fail             ! Indicates a message/warning/error escalated to a fatal error,
  Logical :: SuppressNextTime ! Indicates a message/warning/error which will be suppressed next time.
  Logical :: Error            ! Error flag for message control not found.

  ! Check only one of MessageControlName and iMessageControl present.
  If (Present(MessageControlName) .eqv. Present(iMessageControl)) Then
    Call Message('UNEXPECTED FATAL ERROR in ControlledMessage', 4)
  End If

  ! Set iAction = ErrorCode.
  If (Present(ErrorCode)) Then
    iAction = ErrorCode
  Else
    iAction = 0
  End If

  ! Initialise Fail and SuppressNextTime.
  Fail             = .false.
  SuppressNextTime = .false.

  ! Find index of control to use.
  If (MessageControls%SetupComplete) Then
    If (Present(MessageControlName)) Then
      iMessageControlL = FindMessageControlIndex(MessageControlName, MessageControls, Error)
      If (Error) iMessageControlL = 0
    Else
      iMessageControlL = iMessageControl
      If (iMessageControl < 0 .or. iMessageControl > MessageControls%nMessageControls) Then
        Call Message('UNEXPECTED FATAL ERROR in ControlledMessage', 4)
      End If
    End If
  Else
    iMessageControlL = 0
  End If

  ! If message control available, revise values of iAction, Fail and SuppressNextTime.
  If (iMessageControlL /= 0) Then

    MessageControls%nMessages(iMessageControlL) = MessageControls%nMessages(iMessageControlL) + 1

    ! Record default for first message and check default on subsequent messages. The purpose of this is to
    ! trap the possibility of the message control action being 'default' and the message being invoked with
    ! different defaults in different places. This could allow a warning or error to occur without any warning
    ! or error message being issued.
    If (MessageControls%nMessages(iMessageControlL) == 1) Then
      MessageControls%Defaults(iMessageControlL) = iAction
    Else
      If (MessageControls%Defaults(iMessageControlL) /= iAction) Then
        Call Message('UNEXPECTED FATAL ERROR in ControlledMessage', 4)
      End If
    End If

    ! Determine action to be taken.
    If (MessageControls%MessageControls(iMessageControlL)%iAction /= -1) Then
      iAction = MessageControls%MessageControls(iMessageControlL)%iAction
    End If

    ! Escalate message. Note this is done in a different way if nMessages = 1 to avoid writing 'FATAL ERROR
    ! due to an excessive number of occurrances of the following message' when there have been no previous
    ! messages actually written.
    If (                                                                    &
      MessageControls%MessageControls(iMessageControlL)%Fail       /=       &
      Huge(MessageControls%MessageControls(iMessageControlL)%Fail)    .and. &
      MessageControls%nMessages(iMessageControlL)            >              &
      MessageControls%MessageControls(iMessageControlL)%Fail                &
    ) Then
      If (iAction <= 2) Then
        If (MessageControls%nMessages(iMessageControlL) == 1) Then
          iAction = 3
        Else
          Fail = .true.
        End If
      End If
    End If

    ! Deal with possible suppression of message.
    If (.not. Fail) Then

      ! Suppress message. Note fatal error messages are not suppressed, nor is the first occurrance of
      ! non-fatal errors and warnings, nor is the first occurance of a message if this message may be later
      ! escalated to a fatal error.
      If (                                                                        &
        MessageControls%MessageControls(iMessageControlL)%Suppress       /=       &
        Huge(MessageControls%MessageControls(iMessageControlL)%Suppress)    .and. &
        MessageControls%nMessages(iMessageControlL)                >              &
        MessageControls%MessageControls(iMessageControlL)%Suppress                &
      ) Then
        If (                                                                       &
          iAction > 2                                                         .or. &
          (                                                                        &
            iAction > 0                                      .and.                 &
            MessageControls%nMessages(iMessageControlL) == 1                       &
          )                                                                   .or. &
          (                                                                        &
            iAction == 0                                                .and.      &
            MessageControls%nMessages(iMessageControlL) == 1            .and.      &
            MessageControls%MessageControls(iMessageControlL)%Fail /=              &
            Huge(MessageControls%MessageControls(iMessageControlL)%Fail)           &
          )                                                                        &
        ) Then
        Else
          Return
        End If
      End If

      ! Determine if message will be suppressed next time. Note the 'iAction > 0 ...' and 'iAction == 0 ...'
      ! expressions can never be true, but they are retained for consistency with the if test above. Note also
      ! that the message won't be suppressed next time if it will be escalated to a fatal error.
      If (                                                                        &
        MessageControls%MessageControls(iMessageControlL)%Suppress       /=       &
        Huge(MessageControls%MessageControls(iMessageControlL)%Suppress)    .and. &
        MessageControls%nMessages(iMessageControlL) + 1            >              &
        MessageControls%MessageControls(iMessageControlL)%Suppress                &
      ) Then
        If (                                                                       &
          iAction > 2                                                         .or. &
          (                                                                        &
            iAction > 0 .and.                                                      &
            MessageControls%nMessages(iMessageControlL) + 1 == 1                   &
          )                                                                   .or. &
          (                                                                        &
            iAction == 0                                                .and.      &
            MessageControls%nMessages(iMessageControlL) + 1 == 1        .and.      &
            MessageControls%MessageControls(iMessageControlL)%Fail /=              &
            Huge(MessageControls%MessageControls(iMessageControlL)%Fail)           &
          )                                                                        &
        ) Then
        Else If (                                                               &
          MessageControls%MessageControls(iMessageControlL)%Fail       /=       &
          Huge(MessageControls%MessageControls(iMessageControlL)%Fail)    .and. &
          MessageControls%nMessages(iMessageControlL) + 1        >              &
          MessageControls%MessageControls(iMessageControlL)%Fail                &
        ) Then
        Else
          SuppressNextTime = .true.
        End If
      End If

    End If

  End If

  If (Fail) Then

    Select Case (iAction)

      Case (0)

        Call Message(                                                         &
               'FATAL ERROR due to an excessive number of occurrances ('   // &
               Trim(Int2Char(MessageControls%nMessages(iMessageControlL))) // &
               ') of the following message: '                              // &
               Char,                                                          &
               3,                                                             &
               EndRun                                                         &
             )

      Case (1)

        Call Message(                                                         &
               'FATAL ERROR due to an excessive number of occurrances ('   // &
               Trim(Int2Char(MessageControls%nMessages(iMessageControlL))) // &
               ') of the following warning: '                              // &
               Char,                                                          &
               3,                                                             &
               EndRun                                                         &
             )

      Case (2)

        Call Message(                                                         &
               'FATAL ERROR due to an excessive number of occurrances ('   // &
               Trim(Int2Char(MessageControls%nMessages(iMessageControlL))) // &
               ') of the following error: '                                // &
               Char,                                                          &
               3,                                                             &
               EndRun                                                         &
             )

      Case Default

        Call Message('UNEXPECTED FATAL ERROR in ControlledMessage.', 4)

    End Select

  Else If (SuppressNextTime) Then

    Select Case (iAction)

      Case (0)

        Call Message(                                                                            &
                              Trim(Char) // ' Future messages of this type will be suppressed.', &
               0,                                                                                &
               EndRun                                                                            &
             )

      Case (1)

        Call Message(                                                                            &
               'WARNING: ' // Trim(Char) // ' Future messages of this type will be suppressed.', &
               1,                                                                                &
               EndRun                                                                            &
             )

      Case (2)

        Call Message(                                                                            &
               'ERROR: '   // Trim(Char) // ' Future messages of this type will be suppressed.', &
               2,                                                                                &
               EndRun                                                                            &
             )

      Case Default

        Call Message('UNEXPECTED FATAL ERROR in ControlledMessage.', 4)

    End Select

  Else

    Select Case (iAction)

      Case (0)

        Call Message(                              Char, 0, EndRun)

      Case (1)

        Call Message('WARNING: '                // Char, 1, EndRun)

      Case (2)

        Call Message('ERROR: '                  // Char, 2, EndRun)

      Case (3)

        Call Message('FATAL ERROR: '            // Char, 3, EndRun)

      Case (4)

        Call Message('UNEXPECTED FATAL ERROR: ' // Char, 4, EndRun)

      Case Default

        Call Message('UNEXPECTED FATAL ERROR in ControlledMessage.', 4)

    End Select

  End If

End Subroutine ControlledMessage

!-------------------------------------------------------------------------------------------------------------

End Module ErrorAndMessageIIModule

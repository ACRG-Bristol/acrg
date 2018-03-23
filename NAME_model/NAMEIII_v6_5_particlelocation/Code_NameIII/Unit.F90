! Module: Unit Module

Module UnitModule

! This module provides code for handling Fortran unit numbers for input and output.

!-------------------------------------------------------------------------------------------------------------

#ifdef CompaqPCCompiler
  Use DFPort
#endif
#ifdef IntelLinCompiler
  Use IFLPort
#endif
Use GlobalParametersModule
Use ErrorAndMessageModule
Use StringModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: Units_     ! A collection of information on input/output unit numbers.
Public :: InitUnits  ! Initialises a collection of information on input/output unit numbers.
Public :: GetNewUnit ! Returns a new unit number.
Public :: OpenFile   ! Opens a file.
Public :: CloseUnit  ! Closes and frees up a unit, with option to delete the associated file.

!-------------------------------------------------------------------------------------------------------------

Type :: Units_ ! A collection of information on input/output unit numbers.
  Private
  Logical :: Assigned(LastUnit) ! For each unit number, indicates whether the unit number has been assigned by
                                ! GetNewUnit and not yet unassigned by CloseUnit.
  Logical :: Reserved(LastUnit) ! For each unit number, indicates whether the unit number is reserved for use
                                ! elsewhere. This includes units 5, 6, LogFileUnit and ErrorFileUnit. (Unit 0,
                                ! which also has a special role in Fortran, is not considered since the arrays
                                ! Assigned and Reserved start at 1).
End Type Units_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitUnits(ReservedUnits) Result(Units)
! Initialises a collection of information on input/output unit numbers.

  Implicit None
  ! Argument list:
  Integer, Intent(In), Optional :: ReservedUnits(:)
  ! Function result:
  Type(Units_) :: Units ! Collection of unit numbers to be initialised.
  ! Local:
  Integer :: i ! Loop index.

  Units%Assigned(:)             = .false.
  Units%Reserved(:)             = .false.
  Units%Reserved(5)             = .true.
  Units%Reserved(6)             = .true.
  Units%Reserved(LogFileUnit)   = .true.
  Units%Reserved(ErrorFileUnit) = .true.
  Units%Reserved(ParticleEndLocationFileUnit) = .true.
  If (Present(ReservedUnits)) Then
    Do i = 1, Size(ReservedUnits)
      If (                                     &
        ReservedUnits(i) <  1             .or. &
        ReservedUnits(i) >  LastUnit      .or. &
        ReservedUnits(i) == 5             .or. &
        ReservedUnits(i) == 6             .or. &
        ReservedUnits(i) == LogFileUnit   .or. &
        ReservedUnits(i) == ParticleEndLocationFileUnit   .or. &
        ReservedUnits(i) == ErrorFileUnit      &
      ) Then
        Call Message(                                    &
               'UNEXPECTED ERROR in InitUnits: Unit ' // &
               Trim(Int2Char(ReservedUnits(i)))       // &
               ' cannot be reserved',                    &
               4                                         &
             )
      End If
      Units%Reserved(ReservedUnits(i)) = .true.
    End Do
  End If

End Function InitUnits

!-------------------------------------------------------------------------------------------------------------

Subroutine GetNewUnit(Unit, Units) ! $$ could be made more efficient with stack of unallocated units
! Returns a new unit number.

  Implicit None
  ! Argument list:
  Integer,      Intent(Out)   :: Unit  ! New unit number.
  Type(Units_), Intent(InOut) :: Units ! Collection of unit numbers.
  ! Local:
  Integer :: i ! Loop index.

!$OMP CRITICAL
  Do i = FirstUnit, LastUnit
    If (.not.Units%Assigned(i) .and. .not.Units%Reserved(i)) Then
      Units%Assigned(i) = .true.
      Unit = i
      Exit
    End If
  End Do
!$OMP END CRITICAL
  Return

  Call Message('FATAL ERROR: no spare input/output unit numbers left', 3)

End Subroutine GetNewUnit

!-------------------------------------------------------------------------------------------------------------

Function OpenFile(                                       &
           File,                                         &
           Units,                                        &
           Status, Access, Form, Action, Position, RecL, &
           FileDescription                               &
         )                                               &
Result(Unit)
! Opens a file.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)             :: File            ! File.
  Type(Units_), Intent(InOut)          :: Units           ! Collection of information on input/output unit
                                                          ! numbers.
  Character(*), Intent(In),   Optional :: Status          !} Options for open statement.
  Character(*), Intent(In),   Optional :: Access          !}
  Character(*), Intent(In),   Optional :: Form            !}
  Character(*), Intent(In),   Optional :: Action          !}
  Character(*), Intent(In),   Optional :: Position        !}
  Integer,      Intent(In),   Optional :: RecL            !}
  Character(*), Intent(In),   Optional :: FileDescription ! File description for error message.
  ! Function result:
  Integer :: Unit
  ! Locals:
  Character(MaxCharLength)  :: StatusL   !} Local copies of Status, Access, Form, Action and Position.
  Character(MaxCharLength)  :: AccessL   !}
  Character(MaxCharLength)  :: FormL     !}
  Character(MaxCharLength)  :: ActionL   !}
  Character(MaxCharLength)  :: PositionL !}
  Integer                   :: IOStat    ! Error code for open statement.
  Character(MaxCharLength3) :: IOMsg     ! Location for compiler error message.

  If (Present(Status)) Then
    StatusL = Status
  Else
    StatusL = 'Unknown'
  End If

  If (Present(Access)) Then
    AccessL = Access
  Else
    AccessL = 'Sequential'
  End If

  If (Present(Form)) Then
    FormL = Form
  Else
    If (AccessL .CIEq. 'Sequential') Then
      FormL = 'Formatted'
    Else If (AccessL .CIEq. 'Direct') Then
      FormL = 'Unformatted'
    Else
      Call Message('UNEXPECTED FATAL ERROR in OpenFile', 4) ! $$ trap other invalid values?
    End If
  End If

  If (Present(Action)) Then
    ActionL = Action
  Else
    ActionL = 'ReadWrite'
  End If

  If (Present(Position)) Then
    PositionL = Position
  Else
    PositionL = 'AsIs'
  End If

  Call GetNewUnit(Unit, Units)

  If (Present(RecL)) Then
    Open (                                    &
      Unit     = Unit,                        &
      File     = Trim(ConvertFileName(File)), &
      Status   = StatusL,                     &
      Access   = AccessL,                     &
      Form     = FormL,                       &
      Action   = ActionL,                     &
      Position = PositionL,                   &
      RecL     = RecL,                        &
      IOStat   = IOStat,                      &
      IOMsg    = IOMsg                        &
    )
  Else
    Open (                                    &
      Unit     = Unit,                        &
      File     = Trim(ConvertFileName(File)), &
      Status   = StatusL,                     &
      Access   = AccessL,                     &
      Form     = FormL,                       &
      Action   = ActionL,                     &
      Position = PositionL,                   &
      IOStat   = IOStat,                      &
      IOMsg    = IOMsg                        &
    )
  End If

  If (IOStat /= 0) Then ! $$ improve error messages for issues related to open options
                        !    other than status. Option to return errorcode (and free up unit number).

    If (StatusL .CIEq. 'Old') Then

      If (Present(FileDescription)) Then
        Call Message(                                                               &
               'FATAL ERROR in opening '                                         // &
               Trim(FileDescription)                                             // &
               ' "'                                                              // &
               Trim(ConvertFileName(File))                                       // &
               '" (possibly because the file is expected to exist but does not ' // &
               'or is not in the appropriate directory). '                       // &
               'Error message from compiler is: '                                // &
               Trim(IOMsg),                                                         &
               3                                                                    &
             )
      Else
        Call Message(                                                               &
               'FATAL ERROR in opening file "'                                   // &
               Trim(ConvertFileName(File))                                       // &
               '" (possibly because the file is expected to exist but does not ' // &
               'or is not in the appropriate directory). '                       // &
               'Error message from compiler is: '                                // &
               Trim(IOMsg),                                                         &
               3                                                                    &
             )
      End If

    Else If (StatusL .CIEq. 'New') Then

      If (Present(FileDescription)) Then
        Call Message(                                                               &
               'FATAL ERROR in opening '                                         // &
               Trim(FileDescription)                                             // &
               ' "'                                                              // &
               Trim(ConvertFileName(File))                                       // &
               '" (possibly because the file is not expected to exist but does ' // &
               'or because the directory does not exist). '                      // &
               'Error message from compiler is: '                                // &
               Trim(IOMsg),                                                         &
               3                                                                    &
             )
      Else
        Call Message(                                                               &
               'FATAL ERROR in opening file "'                                   // &
               Trim(ConvertFileName(File))                                       // &
               '" (possibly because the file is not expected to exist but does ' // &
               'or because the directory does not exist). '                      // &
               'Error message from compiler is: '                                // &
               Trim(IOMsg),                                                         &
               3                                                                    &
             )
      End If

    Else

      If (Present(FileDescription)) Then
        Call Message(                                                   &
               'FATAL ERROR in opening '                             // &
               Trim(FileDescription)                                 // &
               ' "'                                                  // &
               Trim(ConvertFileName(File))                           // &
               '" (possibly because the directory does not exist). ' // &
               'Error message from compiler is: '                    // &
               Trim(IOMsg),                                             &
               3                                                        &
             )
      Else
        Call Message(                                                   &
               'FATAL ERROR in opening file "'                       // &
               Trim(ConvertFileName(File))                           // &
               '" (possibly because the directory does not exist). ' // &
               'Error message from compiler is: '                    // &
               Trim(IOMsg),                                             &
               3                                                        &
             )
      End If

    End If

  End If

End Function OpenFile

!-------------------------------------------------------------------------------------------------------------

Subroutine CloseUnit(Unit, Units, DeleteFile)
! Closes and frees up a unit, with option to delete the associated file.

  Implicit None
  ! Argument list:
  Integer,      Intent(In)             :: Unit       ! Unit number to be closed.
  Type(Units_), Intent(InOut)          :: Units      ! Collection of unit numbers.
  Logical,      Intent(In),   Optional :: DeleteFile ! Indicates associated file is to be deleted.

!$OMP CRITICAL
  If (Unit <  FirstUnit .or. Unit >  LastUnit) Then
    Call Message(                                                     &
           'UNEXPECTED FATAL ERROR in CloseUnit: the unit number ' // &
           Trim(Int2Char(Unit))                                    // &
           ' lies outside the permitted range '                    // &
           Trim(Int2Char(FirstUnit))                               // &
           ' to '                                                  // &
           Trim(Int2Char(LastUnit)),                                  &
           4                                                          &
         )
  End If

  If (Units%Reserved(Unit)) Then
    Call Message(                                                     &
           'UNEXPECTED FATAL ERROR in CloseUnit: the unit number ' // &
           Trim(Int2Char(Unit))                                    // &
           ' is not permitted',                                       &
           4                                                          &
         )
  End If

  If (Units%Assigned(Unit)) Then
    Units%Assigned(Unit) = .false.
    If (Present(DeleteFile)) Then
      If (DeleteFile) Then
        Close(Unit, Status = 'Delete')
      Else
        Close(Unit)
      End If
    Else
      Close(Unit)
    End If
  Else
    Call Message(                                                     &
           'UNEXPECTED FATAL ERROR in CloseUnit: the unit number ' // &
           Trim(Int2Char(Unit))                                    // &
           ' is not assigned',                                        &
           4                                                          &
         )
  End If
!$OMP END CRITICAL

End Subroutine CloseUnit

!-------------------------------------------------------------------------------------------------------------

End Module UnitModule

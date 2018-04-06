! Module:  Sort Module

Module SortModule

! This module provides code for sorting.

! Notes
! -----

! Note QSort doesn't pass character string lengths correctly to the ordering subroutine (on both Compaq and
! Intel fortran, although the precise behaviour regarding what works or doesn't work or works only without
! bounds checking is different for the two compilers). Hence receiving the character strings as Character(*)
! is unpredictable. This is why the length is passed via the module variable CharLengthForSort, bypassing
! QSort. $$

!-------------------------------------------------------------------------------------------------------------

#ifdef CompaqPCCompiler
  USE DFPort, Only: QSort
#endif
#ifdef IntelLinCompiler
  USE IFLPort, Only: QSort, SizeOf_Size_T
#endif
Use ErrorAndMessageModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: CharLengthForSort ! Length of character strings to be sorted (used to communicate this information
                            ! between the routine calling QSort and the ordering routine called by QSort).
Public :: SortInPlace       ! Sorts in place.

!-------------------------------------------------------------------------------------------------------------

Interface SortInPlace ! Sorts in place.
  Module Procedure SortInPlaceChar
End Interface

!-------------------------------------------------------------------------------------------------------------

Integer :: CharLengthForSort ! Length of character strings to be sorted (used to communicate this information
                             ! between the routine calling QSort and the ordering routine called by QSort).

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine SortInPlaceChar(OrderingType, C)
! Sorts character strings in place.

  ! Argument list:
  Character(*), Intent(In)    :: OrderingType ! Type of ordering:
                                              !   CI: Case insensitive ordering.
  Character(*), Intent(InOut) :: C(:)         ! Character strings to be sorted.
  ! Locals:
  Integer(2), External :: OrderCharCI ! Routine used to compute order of character strings.

  CharLengthForSort = Len(C(1))

  If (OrderingType == 'CI') Then
#   ifdef CompaqPCCompiler
      Call QSort(                 &
             Array  = C,          &
             Len    = Size(C),    &
             iSize  = Len(C(1)),  &
             Compar = OrderCharCI &
           )
#   endif
#   ifdef IntelLinCompiler
      Call QSort(                                      &
             Array  = C,                               &
             Len    = Size(C, Kind = SizeOf_Size_T),   &
             iSize  = Len(C(1), Kind = SizeOf_Size_T), &
             Compar = OrderCharCI                      &
           )
#   endif
  ! $$ Test for no option given for other compilers.
  Else
    Call Message('UNEXPECTED FATAL ERROR in SortInPlaceChar.', 4)
  End If

End Subroutine SortInPlaceChar

!-------------------------------------------------------------------------------------------------------------

End Module SortModule

!-------------------------------------------------------------------------------------------------------------

Function OrderCharCI(C1, C2)
! Computes order of character strings using a case insensitive ordering.

  Use GlobalParametersModule
  Use SortModule

  Implicit None
  ! Argument list:
  Character(CharLengthForSort), Intent(In) :: C1 !} Character strings to be ordered.
  Character(CharLengthForSort), Intent(In) :: C2 !}
  ! Function result:
  Integer(2) :: OrderCharCI ! -1 if C1 before C2, 0 if C1 = C2, and 1 if C1 after C2.
  ! Locals:
  Integer      :: i  ! Loop index.
  Character(1) :: A1 !} Characters from the strings C1 and C2 respectively.
  Character(1) :: A2 !}

  OrderCharCI = 0
  Do i = 1, Len(C1)
    A1 = C1(i:i)
    If ('a' <= A1 .and. A1 <= 'z') A1 = Char(IChar(A1) - IChar('a') + IChar('A'))
    A2 = C2(i:i)
    If ('a' <= A2 .and. A2 <= 'z') A2 = Char(IChar(A2) - IChar('a') + IChar('A'))
    If (A1 > A2) Then
      OrderCharCI = 1
      Exit
    Else If (A1 < A2) Then
      OrderCharCI = - 1
      Exit
    End If
  End Do

End Function OrderCharCI


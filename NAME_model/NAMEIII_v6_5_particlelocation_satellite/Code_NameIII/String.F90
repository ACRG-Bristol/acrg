! Module:  String Module

Module StringModule

! This module provides code for handling character strings and converting numbers to and from character
! strings in particular formats.

!-------------------------------------------------------------------------------------------------------------

Use GlobalParametersModule
Use ErrorAndMessageModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: Int2Char         ! Converts an integer to a character string.
Public :: Std2Char         ! Converts a real of kind Std to a character string.
Public :: P642Char         ! Converts a real of kind P64 to a character string. $$ other kinds?
Public :: FormatChar       ! Formats a character string.
Public :: Char2Int         ! Converts a character string to an integer.
Public :: Char2Std         ! Converts a character string to a real of kind Std.
Public :: Char2Log         ! Converts a character string to a logical variable.
Public :: Int2Ordinal      ! Converts an integer to the ordinal character string 'st', 'nd', 'rd' or 'th' as
                           ! appropriate.
Public :: ConvertFileName  ! Converts a file name so as to use '\' or '/' as appropriate.
Public :: Operator(.CIEq.) ! Case-insensitive equality of trimmed character strings.

!-------------------------------------------------------------------------------------------------------------

Interface Operator(.CIEq.) ! Case-insensitive equality of trimmed character strings.
  Module Procedure CaseInsensitiveEq
End Interface

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function Int2Char(I, Length, Comma, Justify, FormatString)
! Converts an integer to a character string.

  Implicit None
  ! Argument list:
  Integer,      Intent(In)           :: I            ! Number to be formatted.
  Integer,      Intent(In), Optional :: Length       ! Length of desired output (excluding any comma). If
                                                     ! missing, the length is taken to be the trimmed length
                                                     ! of the left justified formatted I.
  Logical,      Intent(In), Optional :: Comma        ! Indicates a comma is to be added to the output. If
                                                     ! missing, Comma is taken to be false.
  Character(1), Intent(In), Optional :: Justify      ! This can be L, l, C, c, R or r indicating left, centre
                                                     ! or right justification. If missing, Justify is taken to
                                                     ! be L.
  Character(*), Intent(In), Optional :: FormatString ! Character string giving the format to be used. If
                                                     ! missing, list-directed formatting is used.
  ! Function result:
  Character(MaxCharLength) :: Int2Char ! The formatted number padded with blanks.

  If (Present(FormatString)) Then
    Write (Int2Char, '(' // FormatString // ')') I
  Else
    Write (Int2Char, *) I
  End If

  Int2Char = FormatChar(Int2Char, Length, Comma, Justify)

End Function Int2Char

!-------------------------------------------------------------------------------------------------------------

Function Std2Char(R, Length, Comma, Justify, FormatString)
! Converts a real of kind Std to a character string.

  Implicit None
  ! Argument list:
  Real(Std),    Intent(In)           :: R            ! Number to be formatted.
  Integer,      Intent(In), Optional :: Length       ! Length of desired output (excluding any comma). If
                                                     ! missing, the length is taken to be the trimmed length
                                                     ! of the left justified formatted R.
  Logical,      Intent(In), Optional :: Comma        ! Indicates a comma is to be added to the output. If
                                                     ! missing, Comma is taken to be false.
  Character(1), Intent(In), Optional :: Justify      ! This can be L, l, C, c, R or r indicating left, centre
                                                     ! or right justification. If missing, Justify is taken to
                                                     ! be L.
  Character(*), Intent(In), Optional :: FormatString ! Character string giving the format to be used. If
                                                     ! missing, list-directed formatting is used.
  ! Function result:
  Character(MaxCharLength) :: Std2Char ! The formatted number padded with blanks.

  If (Present(FormatString)) Then
    Write (Std2Char, '(' // FormatString // ')') R
  Else
    Write (Std2Char, *) R
  End If

  Std2Char = FormatChar(Std2Char, Length, Comma, Justify)

End Function Std2Char

!-------------------------------------------------------------------------------------------------------------

Function P642Char(R, Length, Comma, Justify, FormatString)
! Converts a real of kind Std to a character string.

  Implicit None
  ! Argument list:
  Real(P64),    Intent(In)           :: R            ! Number to be formatted.
  Integer,      Intent(In), Optional :: Length       ! Length of desired output (excluding any comma). If
                                                     ! missing, the length is taken to be the trimmed length
                                                     ! of the left justified formatted R.
  Logical,      Intent(In), Optional :: Comma        ! Indicates a comma is to be added to the output. If
                                                     ! missing, Comma is taken to be false.
  Character(1), Intent(In), Optional :: Justify      ! This can be L, l, C, c, R or r indicating left, centre
                                                     ! or right justification. If missing, Justify is taken to
                                                     ! be L.
  Character(*), Intent(In), Optional :: FormatString ! Character string giving the format to be used. If
                                                     ! missing, list-directed formatting for type Real(Std) is
                                                     ! used.
  ! Function result:
  Character(MaxCharLength) :: P642Char ! The formatted number padded with blanks.
  ! Locals:
  Real(Std) :: RL ! Local copy of R.

  If (Present(FormatString)) Then
    Write (P642Char, '(' // FormatString // ')') R
  Else
    RL = R
    Write (P642Char, *) RL
  End If

  P642Char = FormatChar(P642Char, Length, Comma, Justify)

End Function P642Char

!-------------------------------------------------------------------------------------------------------------

Function FormatChar(C, Length, Comma, Justify, Truncate)
! Formats a character string.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: C        ! Character string to be formatted.
  Integer,      Intent(In), Optional :: Length   ! Length of desired output (excluding any comma). If missing
                                                 ! or negative, the length is taken to be the trimmed length
                                                 ! of the left justified C, suitably limited by MaxCharLength
                                                 ! and Comma.
  Logical,      Intent(In), Optional :: Comma    ! Indicates a comma is to be added to the output. If missing,
                                                 ! Comma is taken to be false.
  Character(1), Intent(In), Optional :: Justify  ! This can be L, l, C, c, R or r indicating left, centre or
                                                 ! right justification. If missing, Justify is taken to be L.
  Logical,      Intent(In), Optional :: Truncate ! Indicates an over-long string is to be truncated rather
                                                 ! than replaced with asterisks.
  ! Function result:
  Character(MaxCharLength) :: FormatChar ! The formatted character string padded with blanks.
  ! Locals:
  Integer :: LengthL   !} Local copies of Length, Comma and Truncate, fixed up for missing values or negative
  Logical :: CommaL    !} Length.
  Logical :: TruncateL !}
  Integer :: LenC      ! Trimmed length of left adjusted C.
  Integer :: Shift     ! Amount the character string needs to be shifted by to justify it correctly.

  LenC = Len_Trim(AdjustL(C))

  If (Present(Comma)) Then
    CommaL = Comma
  Else
    CommaL = .false.
  End If

  If (Present(Truncate)) Then
    TruncateL = Truncate
  Else
    TruncateL = .false.
  End If

  If (Present(Length)) Then
    If (                                          &
      (Length >  MaxCharLength)              .or. &
      (Length == MaxCharLength .and. CommaL)      &
    ) Then
      Call Message('FATAL ERROR in FormatChar: Value of "Length" is inappropriate', 4)
    End If
    If (Length >= 0) Then
      LengthL = Length
    Else
      If (CommaL) Then
        LengthL = Min(LenC, MaxCharLength - 1)
      Else
        LengthL = Min(LenC, MaxCharLength)
      End If
    End If
  Else
    If (CommaL) Then
      LengthL = Min(LenC, MaxCharLength - 1)
    Else
      LengthL = Min(LenC, MaxCharLength)
    End If
  End If

  If (LenC > LengthL) Then
    If (TruncateL) Then
      ! $$ Note FormatChar(1:LengthL) = AdjustL(C) causes strange behaviour on Compaq compiler. The strange
      !    behaviour is actually in a subsequent call when LenC <= LengthL.
      FormatChar = AdjustL(C)
      FormatChar = FormatChar(1:LengthL)
    Else
      FormatChar = Repeat('*', LengthL)
    End If
  Else If (.not.Present(Justify)) Then
    FormatChar = AdjustL(C)
  Else If (Justify == 'L' .or. Justify == 'l') Then
    FormatChar = AdjustL(C)
  Else If (Justify == 'R' .or. Justify == 'r') Then
    Shift = LengthL - LenC
    FormatChar(1:Shift)     = ' '
    FormatChar(Shift + 1: ) = AdjustL(C)
  Else If (Justify == 'C' .or. Justify == 'c') Then
    Shift = (LengthL - LenC)/2
    FormatChar(1:Shift)     = ' '
    FormatChar(Shift + 1: ) = AdjustL(C)
  Else
    Call Message('FATAL ERROR in FormatChar: Value of "Justify" is inappropriate', 4)
  End If

  If (CommaL) FormatChar(LengthL + 1:LengthL + 1) = ','

End Function FormatChar

!-------------------------------------------------------------------------------------------------------------

Function Char2Int(C, Error)
! Converts a character string to an integer.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)            :: C     ! Character string to be converted. Acceptable values for the
                                               ! character string are those accepted by a list directed read
                                               ! of an integer.
  Integer,      Intent(Out), Optional :: Error ! Error code. 0 indicates no error while other values
                                               ! correspond to the IOStat error codes for a Read statement. If
                                               ! not present, the routine stops when there is an error.
  ! Function result:
  Integer :: Char2Int ! The integer result of the conversion.
  ! Locals:
  Integer :: IOStat ! IOStat error code.

  Read (C, *, IOStat = IOStat) Char2Int

  If (Present(Error)) Then
    Error = IOStat
  Else
    If (IOStat /= 0) Then
      Call Message('FATAL ERROR in Char2Int', 4)
    End If
  End If

End Function Char2Int

!-------------------------------------------------------------------------------------------------------------

Function Char2Std(C, Error)
! Converts a character string to a real of kind Std.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)            :: C     ! Character string to be converted. Acceptable values for the
                                               ! character string are those accepted by a list directed read
                                               ! of a real of kind Std.
  Integer,      Intent(Out), Optional :: Error ! Error code. 0 indicates no error while other values
                                               ! correspond to the IOStat error codes for a Read statement. If
                                               ! not present, the routine stops when there is an error.
  ! Function result:
  Real(Std) :: Char2Std ! The real result of the conversion.
  ! Locals:
  Integer :: IOStat ! IOStat error code.

  Read (C, *, IOStat = IOStat) Char2Std

  If (Present(Error)) Then
    Error = IOStat
  Else
    If (IOStat /= 0) Then
      Call Message('FATAL ERROR in Char2Std', 4)
    End If
  End If

End Function Char2Std

!-------------------------------------------------------------------------------------------------------------

Function Char2Log(C, Error)
! Converts a character string to a logical variable.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)            :: C     ! Character string to be converted. Acceptable values for the
                                               ! character string are those accepted by a list directed read
                                               ! of a logical (i.e. those begining T, .T, t, .t, F, .F, f or
                                               ! .f, possibly preceeded by blanks) and those begining Y, y, N
                                               ! or n (again possibly preceeded by blanks).
  Integer,      Intent(Out), Optional :: Error ! Error code. 0 indicates no error while 1 indicates that C has
                                               ! an inappropriate value. If not present, the routine stops
                                               ! when there is an error.
  ! Function result:
  Logical :: Char2Log ! The logical result of the conversion.
  ! Locals:
  Character(2) :: TempChar ! Temporary variable.

  TempChar = AdjustL(C)

  If (                               &
    (TempChar(1:2) .CIEq. '.t') .or. &
    (TempChar(1:1) .CIEq. 't' ) .or. &
    (TempChar(1:1) .CIEq. 'y' )      &
  ) Then
    Char2Log = .true.
    If (Present(Error)) Error = 0
  Else If (                           &
    (TempChar(1:2) .CIEq. '.f') .or.  &
    (TempChar(1:1) .CIEq. 'f' ) .or.  &
    (TempChar(1:1) .CIEq. 'n' )       &
  ) Then
    Char2Log = .false.
    If (Present(Error)) Error = 0
  Else
    If (Present(Error)) Then
      ! Note Char2Log needs to be defined to ensure no extra error on return from this routine (otherwise Result = 
      ! Char2Log(C, Error) uses an uninitialised variable), but the value should not be used to do anything further if Error is 
      ! set).
      Char2Log = .false. 
      Error = 1
    Else
      Call Message('FATAL ERROR in Char2Log', 4)
    End If
  End If

End Function Char2Log

!-------------------------------------------------------------------------------------------------------------

Function Int2Ordinal(I)
! Converts an integer to the ordinal character string 'st', 'nd', 'rd' or 'th' as appropriate.

  Implicit None
  ! Argument list:
  Integer, Intent(In) :: I ! Integer to be converted.
  ! Function result:
  Character(2) :: Int2Ordinal ! The ordinal string.
  ! Locals:
  Integer :: IMod10  ! |I| Mod 10.
  Integer :: IMod100 ! |I| Mod 100.

  IMod10  = Mod(Abs(I),  10)
  IMod100 = Mod(Abs(I), 100)

  If (IMod100 == 11 .or. IMod100 == 12 .or. IMod100 == 13) Then
    Int2Ordinal = 'th'
  Else If (IMod10 == 1) Then
    Int2Ordinal = 'st'
  Else If (IMod10 == 2) Then
    Int2Ordinal = 'nd'
  Else If (IMod10 == 3) Then
    Int2Ordinal = 'rd'
  Else
    Int2Ordinal = 'th'
  End If

End Function Int2Ordinal

!-------------------------------------------------------------------------------------------------------------

Function CaseInsensitiveEq(C1, C2)
! Tests to see if two character strings are case-insensitively equal when trimmed.

  Implicit none
  ! Argument list:
  Character(*), Intent(In) :: C1 !} The two character strings to be compared.
  Character(*), Intent(In) :: C2 !}
  ! Function result:
  Logical :: CaseInsensitiveEq ! Indicates if the two character strings are case-insensitively equal when
                               ! trimmed.
  ! Locals:
  Character(1) :: A1 !} Characters from the strings C1 and C2 respectively.
  Character(1) :: A2 !}
  Integer      :: i  ! Loop index.

  If (Len_Trim(C1) /= Len_Trim(C2)) Then
    CaseInsensitiveEq = .false.
  Else
    CaseInsensitiveEq = .true.
    Do i = 1, Len_Trim(C1)
      A1 = C1(i:i)
      If ('a' <= A1 .and. A1 <= 'z') A1 = Char(IChar(A1) - IChar('a') + IChar('A'))
      A2 = C2(i:i)
      If ('a' <= A2 .and. A2 <= 'z') A2 = Char(IChar(A2) - IChar('a') + IChar('A'))
      If (A1 /= A2) Then
        CaseInsensitiveEq = .false.
        Exit
      End If
    End Do
  End If

End Function CaseInsensitiveEq

!-------------------------------------------------------------------------------------------------------------

Function ConvertFileName(FileName)
! Converts a file name so as to use '\' or '/' as appropriate.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: FileName ! Name of file.
  ! Function result:
  Character(MaxFileNameLength) :: ConvertFileName ! The result of the conversion.
  ! Locals:
  Integer :: i ! Loop index.

  If (Len_Trim(FileName) > MaxFileNameLength) Then
    Call Message(                                           &
           'FATAL ERROR in ConvertFileName: file name "' // &
           Trim(FileName)                                // &
           '" is too long',                                 &
           4                                                &
         )
  End If

  ConvertFileName = FileName

# ifdef CompaqPCCompiler
    Do i = 1, Len_Trim(ConvertFileName)
      If (ConvertFileName(i:i) == '/') ConvertFileName(i:i) = '\'
    End Do
# endif
# ifdef IntelLinCompiler
    Do i = 1, Len_Trim(ConvertFileName)
      If (ConvertFileName(i:i) == '\') ConvertFileName(i:i) = '/'
    End Do
# endif
# ifdef CrayCLECompiler
    Do i = 1, Len_Trim(ConvertFileName)
      If (ConvertFileName(i:i) == '\') ConvertFileName(i:i) = '/'
    End Do
# endif
# ifdef sun
    Do i = 1, Len_Trim(ConvertFileName)
      If (ConvertFileName(i:i) == '\') ConvertFileName(i:i) = '/'
    End Do
# endif

End Function ConvertFileName

!-------------------------------------------------------------------------------------------------------------

End Module StringModule

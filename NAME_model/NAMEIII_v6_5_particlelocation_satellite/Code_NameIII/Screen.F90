! Module: Screen Module

Module ScreenModule

! This module provides code for producing screen output, including graphics.

!-------------------------------------------------------------------------------------------------------------

#ifdef CompaqPCCompiler
  Use DFLib
#endif
Use GlobalParametersModule
Use ErrorAndMessageModule
Use StringModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: SetActiveFocus         ! Makes the window associated with the unit number active and gives it the
                                 ! focus.
Public :: OpenTextWindow         ! Opens a text window for other than units 0, 5 or 6.
Public :: SetTitleTextLinesUnit6 ! Sets title and number of text lines for output to the unit 6 window.
Public :: OpenGraphWindow        ! Opens a graph window and sets the coords of the corners.
Public :: Plot                   ! Draws a line joining a number of points.
Public :: MsgBox                 ! Produces a message box.

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine SetActiveFocus(Unit)
! Makes the window associated with the unit number active and gives it the focus.

  Implicit None
  ! Argument list:
  Integer, Intent(In) :: Unit ! Unit number.
  ! Locals:
  Integer :: Status ! Results of SetActiveQQ and FocusQQ.

# ifdef CompaqPCCompiler
    Status = SetActiveQQ(Unit)
    If (Status /= 1) Then
      Call Message(                                             &
             'Error in SetActiveFocus: unable to make unit ' // &
             Trim(Int2Char(Unit))                            // &
             ' active',                                         &
             4                                                  &
           )
    End If
    Status = FocusQQ(Unit)
    If (Status /= 0) Then
      Call Message(                                             &
             'Error in SetActiveFocus: unable to give unit ' // &
             Trim(Int2Char(Unit))                            // &
             ' the focus',                                      &
             4                                                  &
           )
    End If
# endif

End Subroutine SetActiveFocus

!-------------------------------------------------------------------------------------------------------------

Subroutine OpenTextWindow(Unit, Title, nCols, nRows)
! Opens a text window for other than units 0, 5 or 6.

! Note RecL in the open statement controls the point at which wrapping occurs, while
! WC%NumTextCols controls the maximum size of the window.

  Implicit None
  ! Argument list:
  Integer,      Intent(In)           :: Unit  ! Unit number of window.
  Character(*), Intent(In)           :: Title ! Title of window.
  Integer,      Intent(In), Optional :: nCols ! Number of columns of text allowed.
  Integer,      Intent(In), Optional :: nRows ! Number of rows of text allowed.
  ! Locals:
# ifdef CompaqPCCompiler
    Type(WindowConfig) :: WC      ! Window configuration.
    Integer            :: Status  ! Results of SetActiveQQ.
    Logical            :: Status2 ! Result of SetWindowConfig.
# endif

# ifdef CompaqPCCompiler

    If (             &
      Unit <= 0 .or. &
      Unit == 5 .or. &
      Unit == 6      &
    ) Then
      Call Message(                                         &
             'Error in OpenTextWindow: the unit number ' // &
             Trim(Int2Char(Unit))                        // &
             ' is not permitted',                           &
             4                                              &
           )
    End If

    If (Present(nCols)) Then
      Open (               &
        Unit    = Unit,    &
        File    = 'User',  &
        IOFocus = .false., &
        Title   = Title,   &
        RecL    = nCols    &
      )
    Else
      Open (               &
        Unit    = Unit,    &
        File    = 'User',  &
        IOFocus = .false., &
        Title   = Title    &
      )
    End If

    WC%NumXPixels = -1
    WC%NumYPixels = -1
    If (Present(nCols)) Then
      WC%NumTextCols = nCols + 2 ! The '+ 2' is an empirical correction.
    Else
      WC%NumTextCols = -1
    End If
    If (Present(nRows)) Then
      WC%NumTextRows = nRows
    Else
      WC%NumTextRows = -1
    End If
    WC%NumColors = -1
    WC%Title     = Title // Char(0)
    WC%FontSize  = -1

    Status = SetActiveQQ(Unit)
    If (Status /= 1) Then
      Call Message(                                             &
             'Error in OpenTextWindow: unable to make unit ' // &
             Trim(Int2Char(Unit))                            // &
             ' active',                                         &
             4                                                  &
           )
    End If

    Status2 = SetWindowConfig(WC)
    If (.not.Status2) Then
      Call Message(                                              &
             'Error in OpenTextWindow: unable to set window ' // &
             'characteristics as desired for unit '           // &
             Trim(Int2Char(Unit))                             // &
             ' - nearest permitted values will be used',         &
             2                                                   &
           )
      Status2 = SetWindowConfig(WC)
    End If

# endif

End Subroutine OpenTextWindow

!-------------------------------------------------------------------------------------------------------------

Subroutine SetTitleTextLinesUnit6(Title, nRows)
! Sets title and number of text lines for output to the unit 6 window.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: Title ! Title of window.
  Integer,      Intent(In), Optional :: nRows ! Number of rows of text allowed.
  ! Locals:
# ifdef CompaqPCCompiler
    Type(WindowConfig) :: WC      ! Window configuration.
    Integer            :: Status  ! Results of SetActiveQQ.
    Logical            :: Status2 ! Result of SetWindowConfig.
# endif

# ifdef CompaqPCCompiler

    WC%NumXPixels  = -1
    WC%NumYPixels  = -1
    WC%NumTextCols = -1
    If (Present(nRows)) Then
      WC%NumTextRows = nRows
    Else
      WC%NumTextRows = -1
    End If
    WC%NumColors = -1
    WC%Title     = Title // Char(0)
    WC%FontSize  = -1

    Status = SetActiveQQ(6)
    If (Status /= 1) Then
      Call Message(                                                     &
             'Error in SetTitleTextLinesUnit6: unable to make unit ' // &
             '6'                                                     // &
             ' active',                                                 &
             4                                                          &
           )
    End If

    Status2 = SetWindowConfig(WC)
    If (.not.Status2) Then
      Call Message(                                                      &
             'Error in SetTitleTextLinesUnit6: unable to set window ' // &
             'characteristics as desired for unit '                   // &
             '6'                                                      // &
             ' - nearest permitted values will be used',                 &
             2                                                           &
           )
      Status2 = SetWindowConfig(WC)
    End If

# endif

End Subroutine SetTitleTextLinesUnit6

!-------------------------------------------------------------------------------------------------------------

Subroutine OpenGraphWindow(Unit, Title, X1, Y1, X2, Y2)
! Opens a graph window and sets the coords of the corners.

  Implicit None
  ! Argument list:
  Integer,      Intent(In) :: Unit  ! Unit number
  Character(*), Intent(In) :: Title ! Title of window.
  Real(Std),    Intent(In) :: X1    !} Coords of upper left corner.
  Real(Std),    Intent(In) :: Y1    !}
  Real(Std),    Intent(In) :: X2    !] Coords of lower right corner.
  Real(Std),    Intent(In) :: Y2    !]
  ! Locals:
  Integer(2) :: Status ! Result of SetWindow.
  Real(P64)  :: WX1    !} Coords of upper left corner.
  Real(P64)  :: WY1    !}
  Real(P64)  :: WX2    !] Coords of lower right corner.
  Real(P64)  :: WY2    !]

# ifdef CompaqPCCompiler

    Open (            &
      Unit  = Unit,   &
      File  = 'User', &
      Title = Title   &
    )

    WX1 = X1
    WY1 = Y1
    WX2 = X2
    WY2 = Y2
    Status = SetWindow(.true._2, WX1, WY1, WX2, WY2)
    If (Status == 0) Then
      Call Message(                                               &
             'Error in OpenGraphWindow: unable to set window ' // &
             'characteristics as desired for unit '            // &
             Trim(Int2Char(Unit)),                                &
             4                                                    &
           )
    End If

# endif

End Subroutine OpenGraphWindow

!-------------------------------------------------------------------------------------------------------------

Subroutine Plot(Unit, X, Y)
! Draws a line joining a number of points.

  Implicit None
  ! Argument list:
  Integer,   Intent(In) :: Unit ! Unit number of window in which the line is to be
                                ! drawn.
  Real(Std), Intent(In) :: X(:) !} Coords of points.
  Real(Std), Intent(In) :: Y(:) !}
  ! Locals:
# ifdef CompaqPCCompiler
    Integer        :: Status  ! Result of SetActiveQQ.
    Type(WXYCoord) :: WXY     ! Result of MoveTo_W.
    Integer(2)     :: Status2 ! Result of LineTo_W.
    Real(P64)      :: WX      !} Coords.
    Real(P64)      :: WY      !}
    Integer        :: i       ! Loop index.
# endif


# ifdef CompaqPCCompiler

    If (Size(X) /= Size(Y)) Then
      Call Message(                                                      &
             'UNEXPECTED FATAL ERROR in Plot: Size of arrays not equal', &
             4                                                           &
           )
    End If

    If (Size(X) < 2) Return

    Status = SetActiveQQ(Unit)
    If (Status /= 1) Then
      Call Message(                                   &
             'Error in Plot: unable to make unit ' // &
             Trim(Int2Char(Unit))                  // &
             ' active',                               &
             4                                        &
           )
    End If

    WX = X(1)
    WY = Y(1)
    Call MoveTo_W(WX, WY, WXY)
    Do i = 2, Size(X)
      WX = X(i)
      WY = Y(i)
      Status2 = LineTo_W(WX, WY)
      If (Status2 == 0) Then
        Call Message('Error in Plot', 4)
      End If
    End Do

# endif

End Subroutine Plot

!-------------------------------------------------------------------------------------------------------------

Function MsgBox(Msg, Caption, MsgBoxType)
! Produces a message box.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Msg        ! Message.
  Character(*), Intent(In) :: Caption    ! Message box caption.
  Integer,      Intent(In) :: MsgBoxType ! Indicates buttons and symbols to be
                                         ! displayed, and which button is the default.
                                         ! The possible values are as follows:
                                         !     MB$OK               = #000 } Buttons
                                         !     MB$OKCancel         = #001 } to
                                         !     MB$AbortRetryIgnore = #002 } display.
                                         !     MB$YesNoCancel      = #003 }
                                         !     MB$YesNo            = #004 }
                                         !     MB$RetryCancel      = #005 }
                                         !     MB$IconStop         = #010 ] Icon to
                                         !     MB$IconQuestion     = #020 ] display.
                                         !     MB$IconExclamation  = #030 ]
                                         !     MB$IconInformation  = #040 ]
                                         !     MB$DefButton1       = #000 } Default
                                         !     MB$DefButton2       = #100 } button.
                                         !     MB$DefButton3       = #200 }
                                         ! They can be combined (at least provided at
                                         ! most 1 value of each type is used) by '+'
                                         ! or by '.or.'.
  ! Function result:
  Integer :: MsgBox ! Indicates which button was pressed. The possible values are as
                    ! follows:
                    !     MB$IDOK     = 1
                    !     MB$IDCancel = 2
                    !     MB$IDAbort  = 3
                    !     MB$IDRetry  = 4
                    !     MB$IDIgnore = 5
                    !     MB$IDYes    = 6
                    !     MB$IDNo     = 7.
                    ! 0 is returned on platforms or compilers that don't support
                    ! message boxes.

# ifdef CompaqPCCompiler
    MsgBox = MessageBoxQQ(               &
               Trim(Msg) // Char(0),     &
               Trim(Caption) // Char(0), &
               MsgBoxType                &
             )
    If (MsgBox == 0) Then
      Call Message(                                                             &
             'Error in MsgBox: insufficient memory to display the message box', &
             4                                                                  &
           )
    End If
# else
    MsgBox = 0
# endif

End Function MsgBox

!-------------------------------------------------------------------------------------------------------------

End Module ScreenModule



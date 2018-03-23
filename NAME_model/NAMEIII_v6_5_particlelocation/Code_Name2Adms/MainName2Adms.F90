! Module: Main Name to Adms Module

!-------------------------------------------------------------------------------------------------------------

Program Name2Adms
! .
  
  Use ServiceModule
   
  Implicit None
  ! Locals: 
  Integer             :: nArguments
  Character(80)       :: Arguments(2)
  Type(Units_)        :: Units
  Integer             :: Unit
  Type(HFBlockForms_) :: HFBlockForms
  Type(HFState_)      :: HFState
  Type(Tokens_)       :: Tokens
  Logical             :: EndOfHFs
  Type(Time_)         :: T
  Real                :: Values(9)
  Character(1024)     :: Line
  Integer             :: i

  ! Initialise error and message module.
  Call InitErrorAndMessageModule

  ! Set up unit 6 window.
  Call SetTitleTextLinesUnit6('Message window', 302)  

  Call SetUpErrorAndMessageModule(       &
         EndRunCaption    = 'Name2Adms', &
         LogFile          = ' ',         &
         ErrorFile        = ' ',         &
         ClosePromptLevel = 0,           &
         Restart          = .false.      &
       )

  Call Message('NAME output to ADMS met input conversion') 

  ! Get command line arguments.
  Call GetCommandLineArguments(nArguments, Arguments)
  If (nArguments /= 2) Then
    Call Message('Error: incorrect number of command line arguments', 3)
  End If

  Units = InitUnits() 
  Call InitTimeModule('Gregorian', IsBackwards = .false.)

  ! Initialise HFState.
  Call InitHFState(Arguments(1:1), HFState)

  Call InitHFBlockForms(HFBlockForms)

  Call InitAndAddHFBlockForm(                                     &
         BlockKey             = 'Fields:',                        &
         NamedBlock           = .false.,                          &
         MatchMultipleLines   = .true.,                           &
         nHeaderLines         = 19,                               &
         HeaderLinesToMatch   = (/ 3, 19 /),                      &
         ColumnKeys           = Reshape(                          &
                                  (/                              &
                                    '                          ', & ! 1
                                    'Wind Speed                ', & ! 2
                                    'Wind Direction (degrees)  ', & ! 3
                                    'Cloud Amount (oktas)      ', & ! 4
                                    'Temperature (C)           ', & ! 5
                                    'Sensible heat flux        ', & ! 6
                                    'Boundary layer depth      ', & ! 7
                                    'Precipitation rate (mm/hr)', & ! 8
                                    'Relative Humidity (%)     ', & ! 9
                                    'T                         ', & ! 1
                                    '                          ', & ! 2
                                    '                          ', & ! 3
                                    '                          ', & ! 4
                                    '                          ', & ! 5
                                    '                          ', & ! 6
                                    '                          ', & ! 7
                                    '                          ', & ! 8
                                    '                          '  & ! 9
                                  /),                             &
                                  (/ 9, 2 /)                      &
                                ),                                &
         TwoD                 = .false.,                          &
         UnrecognisedMessages = .false.,                          &
         DisjointColSpecs     = .true.,                           &
         HFBlockForms         = HFBlockForms                      &
       )

  Unit = OpenFile(       &
           Arguments(2), &
           Units,        &
           RecL = 256    &
         )

  ! Write comments.
  Write (Unit, '(A)')
  Write (Unit, '(A)') 'Note the sensible heat fluxes have been limited so that they are >= -99.9.' 
  Write (Unit, '(A)') 'This is because ADMS regards more negative heat fluxes as errors.'

  ! Write keywords.
  Write (Unit, '(A)')
  Write (Unit, '(A)') 'VARIABLES:'
  Write (Unit, '(A)') '13'
  Write (Unit, '(A)') 'YEAR'
  Write (Unit, '(A)') 'MONTH'
  Write (Unit, '(A)') 'DAY OF MONTH'
  Write (Unit, '(A)') 'DAY'
  Write (Unit, '(A)') 'HOUR'
  Write (Unit, '(A)') 'WIND SPEED                  '
  Write (Unit, '(A)') 'WIND DIRECTION (DEGREES)    '
  Write (Unit, '(A)') 'CLOUD AMOUNT (OKTAS)        '
  Write (Unit, '(A)') 'TEMPERATURE (C)             '
  Write (Unit, '(A)') 'SENSIBLE HEAT FLUX          '
  Write (Unit, '(A)') 'BOUNDARY LAYER DEPTH        '
  Write (Unit, '(A)') 'PRECIPITATION RATE (MM/HOUR)'
  Write (Unit, '(A)') 'RELATIVE HUMIDITY (PERCENT) '
  Write (Unit, '(A)') 
  Write (Unit, '(A)') 'DATA:'

  Do    

    Call ReadHFs(HFBlockForms, Tokens, EndOfHFs, Units, HFState)
    If (EndOfHFs) Exit

    T         = Char2Time(Tokens%Tokens(1))
    Values(2) = Char2Std (Tokens%Tokens(2))
    Values(3) = Char2Std (Tokens%Tokens(3))
    Values(4) = Char2Std (Tokens%Tokens(4))
    Values(5) = Char2Std (Tokens%Tokens(5))
    Values(6) = Char2Std (Tokens%Tokens(6))
    Values(7) = Char2Std (Tokens%Tokens(7))
    Values(8) = Char2Std (Tokens%Tokens(8))
    Values(9) = Char2Std (Tokens%Tokens(9))

    ! Limit sensible heat flux so its >= -99.9.
    If (Values(6) < -99.9) Values(6) = -99.9

    ! Convert missing data values.
    Do i = 2, 9
      If (Values(i) < - Huge(Values) / 2.0) Values(i) = -999.0
    End Do

    Line =                                                          &
      Trim(Int2Char(T%Year,             5, .true., 'R', 'I4'  )) // & ! $$ introduce and use functions
      Trim(Int2Char(T%Month,            3, .true., 'R', 'I2'  )) // & !    to return year etc (to keep elements 
      Trim(Int2Char(T%Day,              3, .true., 'R', 'I2'  )) // & !    private). Could have RealDay,
      Trim(Int2Char(JulianDayNumber(T), 4, .true., 'R', 'I3'  )) // & !    RealHour etc too.
      Trim(Std2Char(Real(T%Hour),       6, .true., 'R', 'F5.2')) // &
      Trim(Std2Char(Values(2),          8, .true., 'R', 'F7.2')) // & ! Wind speed                 
      Trim(Std2Char(Values(3),          7, .true., 'R', 'F6.1')) // & ! Wind direction (degrees)   
      Trim(Std2Char(Values(4),          8, .true., 'R', 'F7.2')) // & ! Cloud amount (oktas)       
      Trim(Std2Char(Values(5),          7, .true., 'R', 'F6.1')) // & ! Temperature (C)            
      Trim(Std2Char(Values(6),          7, .true., 'R', 'F6.1')) // & ! Sensible heat flux         
      Trim(Std2Char(Values(7),          7, .true., 'R', 'F6.1')) // & ! Boundary layer depth       
      Trim(Std2Char(Values(8),          8, .true., 'R', 'F7.2')) // & ! Precipitation rate (mm/hr) 
      Trim(Std2Char(Values(9),          7, .true., 'R', 'F6.1'))      ! Relative humidity (percent)

    Write (Unit, '(A)') Trim(Line)

  End Do

  If (GetErrorCode() == 0) Then
    Call Message(' ')
    Call Message(                                            &
           'Name to Adms conversion completed successfully', &
           0,                                                &
           .true.                                            &
         )
  Else If (GetErrorCode() == 1) Then
    Call Message(' ')
    Call Message(                                                  &
           'Name to Adms conversion completed with some warnings', &
           0,                                                      &
           .true.                                                  &
         )
  Else If (GetErrorCode() == 2) Then
    Call Message(' ')
    Call Message(                                                &
           'Name to Adms conversion completed with some errors', &
           0,                                                    &
           .true.                                                &
         )
  End If
    
End Program Name2Adms

!-------------------------------------------------------------------------------------------------------------

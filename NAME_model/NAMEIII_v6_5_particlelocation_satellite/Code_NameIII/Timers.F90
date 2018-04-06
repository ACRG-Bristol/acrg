! Module:  Timer Module

Module TimerModule

!
! Module for timing sections of the code
!
!
!  At the moment this module only works together with openmp as it uses the openmp
!  function omp_get_wtime. If UseTimers is not defined with -DUseTimers, the subroutines
!  in this module do nothing. Otherwise all output information is written to the file
!  <Filestem>Timing.txt, where <Filestem> is the stem of the name Input file, e.g.
!  Inputfile = 'Example1.txt' => timer information written to 'Example1Timing.txt'.
!  If compiled with -DUseTimers the behaviour of the module can be controlled with a
!  structure of Type(TimerOpts_). Output will only be written to the timing log file
!  if the initialisation routine is called with TimerOpts%UseTimers == .true. Detailled timing
!  information is only written if TimerOpts%SummaryOnly is .false.

!-------------------------------------------------------------------------------------------------------------

! OpenMP library
!$ Use omp_lib

! Used to assign a Unit number to the output file
  Use UnitModule
! Used for warning and error messages
  Use ErrorAndMessageModule

!-------------------------------------------------------------------------------------------------------------

  Implicit None

!-------------------------------------------------------------------------------------------------------------

  Private

  Public :: PreInitTimerOpts
  Public :: InitTimerOpts
  Public :: CheckTimerOpts
  Public :: FinaliseTimers
  Public :: TimerCreate
  Public :: TimerOn, TimerOff
  Public :: TimerWriteLast
  Public :: TimerWriteSummary
  Public :: Timer_
  Public :: TimerOpts_

! Maximal length of timer labels
  Integer, Parameter :: TimerLabelLength = 60

! Timer options type
  Type :: TimerOpts_
    Logical :: Initialised                         ! Have the timer module options been initialised?
    Logical :: UseTimers                           ! Do we use timers?
    Logical :: SummaryOnly                         ! Print timer summary only?
  End Type TimerOpts_

! Timer Type
  Type :: Timer_
    Logical :: IsInitialised = .False.             ! Has this Timer been initialised?
    Logical :: On = .False.                        ! On/Off switch
    Character(len=TimerLabelLength) :: Label       ! Label
    Real*8 :: Start                                ! Most recent start time
    Real*8 :: Finish                               ! Most renet finish time
    Integer :: Calls                               ! Total number of calls
    Real*8 :: TotalTime                            ! Total time the timer was active
    Type(TimerOpts_) :: TimerOpts                  ! Timer options
  End Type Timer_

! Unit for Output (set in InitialiseTimers)
  Integer :: TimerUnit

! Reference Time (set in SetReferenceTime)
  Real*8 :: ReferenceTime

! Have we initialised the Library?
  Logical :: TimersIsInitialised = .false.

!-------------------------------------------------------------------------------------------------------------

  Contains

!-------------------------------------------------------------------------------------------------------------

  Function PreInitTimerOpts()
! Return default values for Type(TimerOpts_)

   Implicit None

   Type(TimerOpts_) :: PreInitTimerOpts

   PreInitTimerOpts%Initialised=.true.
   PreInitTimerOpts%UseTimers=.false.
   PreInitTimerOpts%SummaryOnly=.true.

  End Function PreInitTimerOpts

!-------------------------------------------------------------------------------------------------------------

  Subroutine InitTimerOpts(TimerOpts,Units,TimingFileName)
! Initialise library:
! Open an output unit for timer writes, set reference time

    Implicit None

    Type(TimerOpts_), Intent(InOut) :: TimerOpts

    ! Information about available Units
    Type(Units_), Intent(InOut) :: Units
    ! Stem of Restart File. This is used to construct the name of the
    ! timing output file <RestartFileStem>Timing.txt
    Character(len=*), Intent(In) :: TimingFileName
    Character(len=200) :: Filename

    TimerOpts%Initialised=.true.
#ifdef UseTimers
    If (TimerOpts%UseTimers) Then
      ! Set reference time
      Call SetReferenceTime()
      ! Open output unit
      Filename = TRIM(TimingFileName)
      TimerUnit = OpenFile(Filename,                           &
                           Units,                              &
                           Action='WRITE',                     &
                           FileDescription="Timer output file" &
                  )
    End If
#else
    Call Message("WARNING: Trying to use timers with code that has been compiled " // &
                 "without support for timer module",1)
    TimerOpts%UseTimers=.false.
#endif


  End Subroutine InitTimerOpts

!-------------------------------------------------------------------------------------------------------------

  Subroutine CheckTimerOpts(TimerOpts)
! Check correct initialisation of timer module, print info messages

    Implicit None

    Type(TimerOpts_), Intent(InOut) :: TimerOpts

    If (.not. TimerOpts%Initialised) Then
      Call Message('FATAL ERROR: The Timer options have not been initialised', 3)
    End If

#ifdef UseTimers
    If (TimerOpts%UseTimers) Then
      If (TimerOpts%SummaryOnly) Then
        Call Message('Collecting timer information',0)
      Else
        Call Message('Collecting detailled timer information',0)
      End If
    End If
#endif

  End Subroutine CheckTimerOpts


!-------------------------------------------------------------------------------------------------------------

  Subroutine FinaliseTimers(TimerOpts,Units)
! Finalise module, close output unit

    Implicit None

    Type(TimerOpts_), Intent(In) :: TimerOpts

    ! Information about available Units
    Type(Units_), Intent(InOut) :: Units

#ifdef UseTimers
    If (TimerOpts%UseTimers) Then
      Write(TimerUnit,'("Finalised Timer module")')
      Call CloseUnit(TimerUnit,Units)
    End If
#endif

  End Subroutine FinaliseTimers

!-------------------------------------------------------------------------------------------------------------

  Subroutine SetReferenceTime()
! Set Reference Time

    Implicit None

#ifdef UseTimers
    ReferenceTime = omp_get_wtime()
#endif

  End Subroutine SetReferenceTime

!-------------------------------------------------------------------------------------------------------------

  Subroutine TimerCreate(Timer,Label,TimerOpts)
! Create new Timer with name <Label>

    Implicit None

    Type(Timer_),     Intent(InOut) :: Timer
    Character(len=*), Intent(In)    :: Label
    Type(TimerOpts_), Intent(In)    :: TimerOpts

#ifdef UseTimers
    If (TimerOpts%UseTimers) Then
      Timer%IsInitialised = .true.
      Timer%On = .false.
      Timer%Label = Label
      Timer%Start = -1.0
      Timer%Finish = -1.0
      Timer%Calls = 0
      Timer%TotalTime = 0.0
      Timer%TimerOpts = TimerOpts
      Write(TimerUnit,*) "Created Timer ",TRIM(Label)
    End If
#endif

  End Subroutine TimerCreate

!-------------------------------------------------------------------------------------------------------------

  Subroutine TimerOn(Timer)
! Switch Timer on

    Implicit None

    Type(Timer_), Intent(InOut) :: Timer

#ifdef UseTimers
    If (Timer%TimerOpts%UseTimers) Then
      If (Timer%IsInitialised) Then
        If (Timer%On) Then
          Write(TimerUnit,*) "WARNING: Timer ",TRIM(Timer%Label)," is already switched on."
        Else
          Timer%Start = omp_get_wtime()
          Timer%On = .true.
        End If
      Else
        Write(TimerUnit,*) "ERROR: Timer ",TRIM(Timer%Label)," has not been initialised."
      End If
    End If
#endif

  End Subroutine TimerOn

!-------------------------------------------------------------------------------------------------------------

  Subroutine TimerOff(Timer)
! Switch Timer off

    Implicit None

    Type(Timer_), Intent(InOut) :: Timer

#ifdef UseTimers
    If (Timer%TimerOpts%UseTimers) Then
      If (Timer%IsInitialised) Then
        If (.not.Timer%On) Then
          Write(TimerUnit,*) "ERROR: Timer ",TRIM(Timer%Label)," is already switched off."
        Else
          Timer%On = .false.
          Timer%Finish = omp_get_wtime()
          Timer%Calls = Timer%Calls + 1
          Timer%TotalTime = Timer%TotalTime + (Timer%Finish - Timer%Start)
        End If
      Else
        Write(TimerUnit,*) "ERROR: Timer ",TRIM(Timer%Label)," has not been initialised."
      End If
    End IF
#endif

  End Subroutine TimerOff

!-------------------------------------------------------------------------------------------------------------

  Subroutine TimerWriteLast(Timer)
! Print out information about last timer call

    Implicit None

    Type(Timer_), Intent(In) :: Timer

#ifdef UseTimers
    If (Timer%TimerOpts%UseTimers .and. (.not. Timer%TimerOpts%SummaryOnly)) Then
      If (Timer%IsInitialised) Then
          If (Timer%On) Then
          Write(TimerUnit,*) "ERROR: Timer ",TRIM(Timer%Label)," is currently switched on."
        Else
          Write(TimerUnit,'("Timer ",A60," ")',advance='no') TRIM(Timer%Label)
          Write(TimerUnit,'(" call ",I10," time ",F20.3)',advance='no') Timer%Calls, &
            (Timer%Finish - Timer%Start)
          Write(TimerUnit,'(" interval = ",F20.3," ",F20.3)') &
            Timer%Start-ReferenceTime, Timer%Finish-ReferenceTime
        End If
      Else
        Write(TimerUnit,*) "ERROR: Timer ",TRIM(Timer%Label)," has not been initialied."
      End If
    End If
#endif

  End Subroutine TimerWriteLast

!-------------------------------------------------------------------------------------------------------------

  Subroutine TimerWriteSummary(Timer)
! Print full Timer information

    Implicit None

    Type(Timer_), Intent(In) :: Timer

#ifdef UseTimers
    If (Timer%TimerOpts%UseTimers) Then
      If (Timer%IsInitialised) Then
        If (Timer%On) Then
          Write(TimerUnit,*) "ERROR: Timer ",TRIM(Timer%Label)," is currently switched on."
        Else
          Write(TimerUnit,*) "Timer ",TRIM(Timer%Label), " summary:"
          Write(TimerUnit, '("   Total time     ",F20.2)') Timer%TotalTime
          Write(TimerUnit, '("   Calls          ",I10)') Timer%Calls
          If ( Timer%Calls > 0 ) Then
            Write(TimerUnit, '("   Time per call  ",F20.4)') Timer%TotalTime/Timer%Calls
          Else
            Write(TimerUnit, '("   Time per call  ",F20.4)') 0.0
          End If
        End If
      Else
        Write(TimerUnit,*) "ERROR: Timer ",TRIM(Timer%Label)," has not been initialied."
      End If
    End If
#endif

  End Subroutine TimerWriteSummary

!-------------------------------------------------------------------------------------------------------------

End Module TimerModule

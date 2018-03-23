! Module:  Time Module

Module TimeModule

! This module provides code for handling time.

! There are three derived types used to store times. Time_ and ShortTime_ both store similar information - the
! differences are that ShortTime_ stores it in a more efficient but more cryptic form and doesn't store time
! zone information. WildTime_ is like Time_ but allows wild cards for some of the elements and can include day
! of the week information.

! Times can be times of events or time intervals.
! Times can be model times (related to the time within the model) or clock times (related to the time when the
! simulation is performed).
! Time_ and ShortTime_ can be any of these times while WildTime_ can only be a model time of an event.

! The module can be used in various 'time frames', although only one time frame can be used in any one run.
! The time frames correspond to the values of the variables Frame and iFrame which is declared globally in the
! module. The time frames are:
! (1) Gregorian:    times defined using the Gregorian calendar.
! (2) 360-day year: times defined using a year of 12 30-day months; days of the week are defined by taking
!                   1/1/2000 to be a Saturday (as it is in the Gregorian calendar).
! (3) relative:     times defined relative to some fixed but unknown reference time; this time frame makes no
!                   use of months, years or time zones.
! Time frames (1) and (2) are collectively called 'calendar time frames'.

! The module can be used in forwards or backwards mode, although only one mode can be used in any one run.
! The modes correspond to the values of the variable Backwards which is declared globally in the module.
! In forwards mode time within the model progresses in the normal chronological direction.
! In backwards mode time within the model progresses in the opposite direction to the normal chronological
! direction.
!
! In backwards mode we distinguish between the concepts of 'chronological' and 'run direction' ordering. The
! time of an event stored in the derived type Time_, ShortTime_ or WildTime_ still represents the true time
! and a time interval dT = Time1 - Time2 stored in the derived type Time_ or ShortTime_ is still stored as
! positive if Time1 is chronologically later than Time2. The operators +, -, *, and / also work as in forwards
! mode. However the following refer to the run direction:
!     (1) >=, >, <=, <, TMin, TMax, and rounding up and down in Round, RoundToYear, RoundToMonth
!     (2) IsInfFuture, IsInfPast, InfFutureTime, InfPastTime, InfFutureShortTime and InfPastShortTime
!     (3) Real time intervals created with ShortTime2RealTime or used in RealTime2ShortTime.
!     (4) Character time intervals used in Char2Time or created with Time2Char or FileNameTime.
!     (5) InWildTimeInterval.

! Note clock times aren't affected by the time frame or the mode. They always use the Gregorian time frame and
! the forwards mode.

! Notes on +, -, *, / and rounding for times:
! Suppose
!     Time1, Time2 and Time3 are the same type as regards model/clock times or times of events/time intervals,
!     dT1 and dT2 are time intervals and the same type as regards model/clock times,
!     n is an integer.
! Then the following types of calculation are valid:
!     Time2 = Time1 + dT1   (but not +infinity + -infinity or -infinity + +infinity)
!     Time2 = Time1 - dT1   (but not +infinity - +infinity or -infinity - -infinity)
!     dT1   = Time1 - Time2 (but not +infinity - +infinity or -infinity - -infinity)
!     dT2   = dT1 * n       (but not +/-infinity * 0)
!     dT2   = dT1 / n       (but not with infinite times or zero denominators)
!     n     = dT1 / dT2     (but not with infinite times or zero denominators)
!     Time3 = RoundTime(Time1, Time2, dT) (but not with infinite times or zero dT).

!-------------------------------------------------------------------------------------------------------------

Use GlobalParametersModule
Use ErrorAndMessageModule
Use StringModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: Time_               ! Information defining a time.
Public :: ShortTime_          ! Abbreviated form of the derived type Time_.
Public :: WildTime_           ! A wild-card time, i.e. a non-interval model time, possibly with wild cards for
                              ! some of the elements or, for calendar time frames, with a day of the week
                              ! specified.
Public :: WildTimePair_       ! A pair of wild-card times, used to specify a collection of time intervals.
Public :: InitTimeModule      ! Initialises the time module. This must be called before other parts of the
                              ! module are invoked in connection with model times. Initialisation is not
                              ! necessary when dealing with clock times.
Public :: GetFrame            ! Gets the time frame used.
Public :: IsGregorian         ! Indicates whether the run uses the Gregorian time frame.
Public :: Is360DayYear        ! Indicates whether the run uses the 360-day year time frame.
Public :: IsRelative          ! Indicates whether the run uses the relative time frame.
Public :: IsCalendar          ! Indicates whether the run uses a calendar time frame.
Public :: IsBackwards         ! Indicates whether the run is backwards.
Public :: Operator(==)        ! Equality of times.
Public :: Operator(/=)        ! Inequality of times.
Public :: Operator(>=)        ! Comparison of times (>=), based on the run direction.
Public :: Operator(>)         ! Comparison of times (>), based on the run direction.
Public :: Operator(<=)        ! Comparison of times (<=), based on the run direction.
Public :: Operator(<)         ! Comparison of times (<), based on the run direction.
Public :: Operator(+)         ! Addition of times.
Public :: Operator(-)         ! Subtraction of times.
Public :: Operator(*)         ! Multiplication of a time by an integer.
Public :: Operator(/)         ! Division of a time by an integer or of one time by another (rounded towards
                              ! zero).
Public :: Round               ! Rounds time.
Public :: RoundToYear         ! Rounds time to the start or end of the year.
Public :: RoundToMonth        ! Rounds time to the start or end of the month.
Public :: TMin                ! Minimum of two times.
Public :: TMax                ! Maximum of two times.
Public :: IsInfFuture         ! Finds out whether a time is in the infinite future (in the run direction for
                              ! model times).
Public :: IsInfPast           ! Finds out whether a time is in the infinite past (in the run direction for
                              ! model times) .
Public :: IsTimeInterval      ! Finds out whether a time is a time interval.
Public :: IsClockTime         ! Finds out whether a time is a clock time.
Public :: InfFutureTime       ! Returns a time in the infinite future (in run direction for model times).
Public :: InfPastTime         ! Returns a time in the infinite past (in run direction for model times).
Public :: ZeroTime            ! Returns a zero time interval.
Public :: ReferenceTime       ! Returns a reference time (for when one wants a time variable set to a fixed
                              ! time, but the time itself doesn't matter).
Public :: Char2Time           ! Converts a character string to a time.
Public :: Time2Char           ! Converts a time to a character string.
Public :: FileNameTime        ! Converts a time to a character string in a form suitable for use in a file
                              ! name.
Public :: ChangeTimeZone      ! Changes the time zone of a time, while keeping the time the same.
Public :: DayOfWeek           ! Evaluates the day of the week (Monday = 1, Sunday = 7) for finite
                              ! calendar-frame times and returns zero for relative-frame times or infinite
                              ! times.
Public :: JulianDayNumber     ! Evaluates the Julian day number for finite calendar-frame times and returns
                              ! zero for relative-frame times or infinite times.
Public :: CurrentClockTime    ! Returns the current clock time.
Public :: InfFutureShortTime  ! Returns a model time in the infinite future (in run direction for model
                              ! times).
Public :: InfPastShortTime    ! Returns a model time in the infinite past (in run direction for model times).
Public :: ZeroShortTime       ! Returns a zero short time interval.
Public :: ReferenceShortTime  ! Returns a reference short time (for when one wants a short time variable set
                              ! to a fixed time, but the time itself doesn't matter).
Public :: Time2ShortTime      ! Converts a time to a short time.
Public :: ShortTime2Time      ! Converts a short time to a time.
Public :: ShortTime2RealTime  ! Converts a short time interval to a duration in seconds.
Public :: RealTime2ShortTime  ! Converts a duration in seconds to a short time interval.
Public :: Char2WildTime       ! Converts a character string to a wild-card time.
Public :: SubstituteWildCards ! Substitutes for the wild cards in a wild-card time using a given time.
Public :: InitWildTimePair    ! Initialises a wild-card time pair.
Public :: InWildTimeInterval  ! Checks whether a time lies in the intervals defined by a pair of wild-card
                              ! times.

!-------------------------------------------------------------------------------------------------------------

! Parameters associated with times:
Integer, Parameter :: YearOrigin   = 2000     !} Year, month and day origin for times held as short times in
Integer, Parameter :: MonthOrigin  = 1        !} calendar time frames. Also values used as substitutes for
Integer, Parameter :: DayOrigin    = 1        !} wild cards in WildTime_. For the latter role in the Gregorian
                                              !} time frame, YearOrigin must be a leap year and MonthOrigin a
                                              !} 31 day month, with Month in [1, 12] and Day in [1, 28]. For
                                              !} the latter role in the 360-day year time frame, Month must be
                                              !} in [1, 12] and Day in [1, 30]. This prevents the non-wild
                                              !} card elements being corrupted by justification.
Integer, Parameter :: FracsPerSec  = 10000000 ! Number of fractions per second.
Integer, Parameter :: F_Gregorian  = 1        !} Codes for time frames "Gregorian", "360-day year", and
Integer, Parameter :: F_360DayYear = 2        !} "Relative".
Integer, Parameter :: F_Relative   = 3        !}

!-------------------------------------------------------------------------------------------------------------

Type :: Time_ ! Information defining a time.
  Logical(1) :: Interval   ! Indicates whether the time is a time interval or the time of an event.
  Logical(1) :: Clock      ! Indicates whether the time is a clock time or a model time.
  Integer(2) :: Infinite   ! 1 indicates an infinite time, 0 a finite time and -1 a negative infinite time.
                           ! The remaining variables in this type are only defined if the time is finite.
  Integer    :: Year       ! Year.                         } Year and Month are defined for non-interval
  Integer    :: Month      ! Month.                        } calendar-frame times only. For such times
  Integer    :: Day        ! Day.                          } midnight at the start of 2000 AD corresponds to
  Integer    :: Hour       ! Hour.                         } values of 2000, 1, 1, 0, 0, 0, 0. For other
  Integer    :: Minute     ! Minute.                       } times, the reference time corresponds to zero
  Integer    :: Second     ! Second.                       } values (in particular, in contrast to
  Integer    :: Fraction   ! Fraction of a second in units } non-interval calendar-frame times, Day = 0 on the
                           !   defined by FracsPerSec.     } first day).
  Integer    :: ZoneHour   !} The time zone in which year, month, day and hour are expressed, specified as
  Integer    :: ZoneMinute !} hours and minutes ahead of UTC (so that time zones to the east of the Greenwich
                           !} meridian will tend to have a positive value and time zones to the west will tend
                           !} to have a negative value). These variables are defined for non-interval
                           !} calendar-frame times only.

  ! Note the values should be 'justified', i.e. Month, Day, Hour, Minute, Second, Fraction, and ZoneMinute
  ! should lie in [1, 12], [1, number of days in the month], [0, 23], [0, 59], [0, 59], [0, FracsPerSec - 1]
  ! and [-59, 59] respectively and ZoneHour and ZoneMinute should have the same sign. An exception concerns
  ! time intervals and relative-frame times where Day has no upper or lower limit.
End Type Time_

!-------------------------------------------------------------------------------------------------------------

Type :: ShortTime_ ! Abbreviated form of the derived type Time_. From this abbreviated form, it is possible to
                   ! reconstruct the full time, although, for calendar-frame times, the reconstructed time
                   ! will always be expressed in UTC and so any time zone information is lost.
  Logical(1)   :: Interval ! Indicates whether the time is a time interval or the time of an event.
  Logical(1)   :: Clock    ! Indicates whether the time is a clock time or a model time.
  Integer(2)   :: Infinite ! 1 indicates an infinite time, 0 a finite time and -1 a negative infinite time.
  Integer(I64) :: Fraction ! Time expressed in fractions of a second (in units defined by FracsPerSec)
                           ! relative to midnight at the start of the day defined in UTC by YearOrigin,
                           ! MonthOrigin and DayOrigin for non-interval calendar-frame times and relative to
                           ! zero for time intervals or relative-frame times. Defined only if the time is
                           ! finite.
End Type ShortTime_

!-------------------------------------------------------------------------------------------------------------

Type :: WildTime_ ! A wild-card time, i.e. a non-interval model time, possibly with wild cards for some of the
                  ! elements or, for calendar time frames, with a day of the week specified.
  Private
  Integer     :: nNonWildCards ! Defines which elements (out of year, month, day, hour, minute, second and
                               ! fraction for calendar time frames when day of week isn't given, out of hour,
                               ! minute, second and fraction for calendar time frames when day of week is
                               ! given, and out of day, hour, minute, second and fraction for relative time
                               ! frames) are not wild cards by giving the number of such elements. The wild
                               ! card elements are always the first so many (possibly none but not all)
                               ! elements from the above list. If there are no wild cards then nNonWildCards =
                               ! 7, 4 and 5 for calendar time frames when day of week isn't given, for
                               ! calendar time frames when day of week is given, and for relative time frames.
  Integer     :: DayOfWeek     ! Day of week (Monday = 1, Sunday = 7) with 0 indicating value unspecified.
  Type(Time_) :: Time          ! A non-interval model time. For calendar time frames when day of week isn't
                               ! given, values of year, month, day, hour, minute and second corresponding to
                               ! wild cards are set equal to YearOrigin, MonthOrigin, DayOrigin, 0, 0 and 0.
                               ! For calendar time frames when day of week is given, values of year, month and
                               ! day are set equal to YearOrigin, MonthOrigin and DayOrigin (wild cards are
                               ! not allowed). For relative time frames, values of day, hour, minute and
                               ! second corresponding to wild cards are set equal to 0, 0, 0 and 0.

  ! The following combinations are not allowed:
  !    Day of week specified and relative time frame,
  !    Day of week specified and wild cards,
  !    Day of week specified and infinite time,
  !    Wild cards and infinite time.

  ! If day of week is not specified, values from another time can be substituted for wild cards using the
  ! routine SubstituteWildCards. When this is done, impossible times can appear in the Gregorian time frame,
  ! such as 29/2/2001 hh:mm or 31/4/2000 hh:mm. These are interpreted as midnight at the start of the
  ! following month.
End Type WildTime_

!-------------------------------------------------------------------------------------------------------------

Type :: WildTimePair_ ! A pair of wild-card times, used to specify a collection of time intervals.
  Private
  Type(WildTime_) :: FromTime !} The wild-card times. Both times must be the same as regards their time zone,
  Type(WildTime_) :: ToTime   !} the number of wild cards, and whether they specify a day of the week.

  ! If no days of the week are used, FromTime and ToTime define a collection of time intervals as follows.
  ! When values are substituted for the wild cards we get a sequence of times SubFromTime_i and SubToTime_i.
  ! Now we consider two cases, FromTime <= ToTime (i.e. SubFromTime_i <= SubToTime_i for all i) and FromTime >
  ! ToTime. In the 1st case the substituted times satisfy
  !            SubFromTime_i <= SubToTime_i <= SubFromTime_i+1
  ! and the time intervals are
  !                    [SubFromTime_i, SubToTime_i),
  ! while in the 2nd case the substituted times satisfy
  !            SubFromTime_i-1 <= SubToTime_i <= SubFromTime_i
  ! and the time intervals are the complement of the time intervals
  !                    [SubToTime_i, SubFromTime_i).
  ! The second case leads to a cyclic interpretation (even if no wild cards are specified), so that, assuming
  ! the model is running forwards, 1/1/2000 22:00 lies in the collection of time intervals defined by
  ! (* / * / * 18:00, * / * / * 06:00).

  ! If days of the week are used, ... $$

  ! [Note the fortran preprocessor fpp doesn't like * / * / * without spaces - hence the spaces above.]
End Type WildTimePair_

!-------------------------------------------------------------------------------------------------------------

Interface Operator(==) ! Equality of times.
  Module Procedure TimeEq
  Module Procedure ShortTimeEq
  Module Procedure WildTimeEq
  Module Procedure WildTimePairEq
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface Operator(/=) ! Inequality of times.
  Module Procedure TimeNE
  Module Procedure ShortTimeNE
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface Operator(>=) ! Comparison of times (>=), based on the run direction.
  Module Procedure TimeGE
  Module Procedure ShortTimeGE
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface Operator(>) ! Comparison of times (>), based on the run direction.
  Module Procedure TimeGT
  Module Procedure ShortTimeGT
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface Operator(<=) ! Comparison of times (<=), based on the run direction.
  Module Procedure TimeLE
  Module Procedure ShortTimeLE
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface Operator(<) ! Comparison of times (<), based on the run direction.
  Module Procedure TimeLT
  Module Procedure ShortTimeLT
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface Operator(+) ! Addition of times.
  Module Procedure AddTime
  Module Procedure AddShortTime
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface Operator(-) ! Subtraction of times.
  Module Procedure SubtractTime
  Module Procedure SubtractShortTime
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface Operator(*) ! Multiplication of a time by an integer.
  Module Procedure MultiplyTime
  Module Procedure MultiplyShortTime
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface Operator(/) ! Division of a time by an integer or of one time by another (rounded towards zero).
  Module Procedure DivideTime
  Module Procedure RatioOfTimes
  Module Procedure DivideShortTime
  Module Procedure RatioOfShortTimes
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface Round ! Rounds time.
  Module Procedure RoundTime
  Module Procedure RoundShortTime
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface RoundToYear ! Rounds time to the start or end of the year.
  Module Procedure RoundTimeToYear
  Module Procedure RoundShortTimeToYear
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface RoundToMonth ! Rounds time to the start or end of the month.
  Module Procedure RoundTimeToMonth
  Module Procedure RoundShortTimeToMonth
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface TMin ! Minimum of two times.
  Module Procedure TimeMin
  Module Procedure ShortTimeMin
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface TMax ! Maximum of two times.
  Module Procedure TimeMax
  Module Procedure ShortTimeMax
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface IsInfFuture ! Finds out whether a time is in the infinite future (in the run direction for model
                      ! times).
  Module Procedure TimeIsInfFuture
  Module Procedure ShortTimeIsInfFuture
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface IsInfPast ! Finds out whether a time is in the infinite past (in the run direction for model times).
  Module Procedure TimeIsInfPast
  Module Procedure ShortTimeIsInfPast
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface IsTimeInterval ! Finds out whether a time is a time interval.
  Module Procedure TimeIsTimeInterval
  Module Procedure ShortTimeIsTimeInterval
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface IsClockTime ! Finds out whether a time is a clock time.
  Module Procedure TimeIsClockTime
  Module Procedure ShortTimeIsClockTime
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface UsesYearAndMonth ! Finds out whether a time uses year, month and time zone information.
  Module Procedure TimeUsesYearAndMonth
  Module Procedure ShortTimeUsesYearAndMonth
End Interface

!-------------------------------------------------------------------------------------------------------------

! Variable module properties:
Logical,                  Private, Save :: Initialised = .false. ! Indicates the module has been initialised.
Character(MaxCharLength), Private, Save :: Frame                 ! The time frame used.
Integer,                  Private, Save :: iFrame                ! Code for the time frame used.
Logical,                  Private, Save :: Backwards             ! Indicates that model time runs backwards.
! American time formats? $$

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine InitTimeModule(FrameName, IsBackwards)
! Initialises the time module. This must be called before other parts of the module are invoked in connection
! with model times. Initialisation is not necessary when dealing with clock times.

! Module initialisation is checked for in the time frame query routines (IsGregorian etc), in IsBackwards,
! and, when creating model times but not clock times, in the time creation routines InfFutureTime,
! InfPastTime, ZeroTime, ReferenceTime, Char2Time, InfFutureShortTime, InfPastShortTime, ZeroShortTime,
! ReferenceShortTime, RealTime2ShortTime and Char2WildTime. In the time creation routines this checking is
! sometimes indirect via a call to a time frame query routine or to IsBackwards. Module initialisation isn't
! necessary for clock times.

! This is sufficient checking provided the time types are treated as if their elements are private. It isn't
! possible however to declare them private because they need to be written to and read from restart files, and
! Fortran doesn't permit this with private elements.
! $$ query Intel and fortran standards group on this?
! $$ if can set private, add functions to return Year, Month, Day, Hour ... ZoneMinute and RealDay,
! $$ RealHour ...
! RealZoneHour.

! Note initialisation might not be checked for if uninitialised times are used as these won't have been
! created by a creation routine. However in such cases results would be unpredictable even if initialisation
! was checked for.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: FrameName   ! The time frame used.
  Logical,      Intent(In) :: IsBackwards ! Indicates that the run is backwards.

  If (Initialised) Call Message('UNEXPECTED FATAL ERROR in InitTimeModule', 4)

  If (Len_Trim(FrameName) > MaxCharLength) Call Message('UNEXPECTED FATAL ERROR in InitTimeModule', 4)
  Frame = FrameName

  If (FrameName .CIEq. 'Gregorian') Then
    iFrame = F_Gregorian
  Else If (FrameName .CIEq. '360-day year') Then
    iFrame = F_360DayYear
  Else If (FrameName .CIEq. 'Relative') Then
    iFrame = F_Relative
  Else
    Call Message('FATAL ERROR: Time frame "' // Trim(FrameName) // '" not recognised', 3)
  End If

  Backwards = IsBackwards

  Initialised = .true.

End Subroutine InitTimeModule

!-------------------------------------------------------------------------------------------------------------

Function GetFrame()
! Gets the time frame used.

  Implicit None
  ! Function result:
  Character(MaxCharLength) :: GetFrame ! The time frame used.

  If (.not.Initialised) Call Message('UNEXPECTED FATAL ERROR in GetFrame', 4)

  GetFrame = Frame

End Function GetFrame

!-------------------------------------------------------------------------------------------------------------

Function IsGregorian()
! Indicates whether the run uses the Gregorian time frame.

  Implicit None
  ! Function result:
  Logical :: IsGregorian ! Indicates that the run uses the Gregorian time frame.

  If (.not.Initialised) Call Message('UNEXPECTED FATAL ERROR in IsGregorian', 4)

  IsGregorian = iFrame == F_Gregorian

End Function IsGregorian

!-------------------------------------------------------------------------------------------------------------

Function Is360DayYear()
! Indicates whether the run uses the 360-day year time frame.

  Implicit None
  ! Function result:
  Logical :: Is360DayYear ! Indicates that the run uses the 360-day year time frame.

  If (.not.Initialised) Call Message('UNEXPECTED FATAL ERROR in Is360DayYear', 4)

  Is360DayYear = iFrame == F_360DayYear

End Function Is360DayYear

!-------------------------------------------------------------------------------------------------------------

Function IsRelative()
! Indicates whether the run uses the relative time frame.

  Implicit None
  ! Function result:
  Logical :: IsRelative ! Indicates that the run uses the relative time frame.

  If (.not.Initialised) Call Message('UNEXPECTED FATAL ERROR in IsRelative', 4)

  IsRelative = iFrame == F_Relative

End Function IsRelative

!-------------------------------------------------------------------------------------------------------------

Function IsCalendar()
! Indicates whether the run uses a calendar time frame.

  Implicit None
  ! Function result:
  Logical :: IsCalendar ! Indicates that the run uses a calendar time frame.

  If (.not.Initialised) Call Message('UNEXPECTED FATAL ERROR in IsCalendar', 4)

  IsCalendar = iFrame == F_Gregorian .or. iFrame == F_360DayYear

End Function IsCalendar

!-------------------------------------------------------------------------------------------------------------

Function IsBackwards()
! Indicates whether the run is backwards.

  Implicit None
  ! Function result:
  Logical :: IsBackwards ! Indicates that the run is backwards.

  If (.not.Initialised) Call Message('UNEXPECTED FATAL ERROR in IsBackwards', 4)

  IsBackwards = Backwards

End Function IsBackwards

!-------------------------------------------------------------------------------------------------------------

Function TimeEq(Time1, Time2)
! Tests for Time1 = Time2.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time1 !} Times to be compared. They must both be of the same type.
  Type(Time_), Intent(In) :: Time2 !}
  ! Function result:
  Logical :: TimeEq ! Indicates whether Time1 = Time2.

  TimeEq = Time2ShortTime(Time1) == Time2ShortTime(Time2)

End Function TimeEq

!-------------------------------------------------------------------------------------------------------------

Function TimeNE(Time1, Time2)
! Tests for Time1 /= Time2.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time1 !} Times to be compared. They must both be of the same type.
  Type(Time_), Intent(In) :: Time2 !}
  ! Function result:
  Logical :: TimeNE ! Indicates whether Time1 /= Time2.

  TimeNE = .not. (Time2 == Time1)

End Function TimeNE

!-------------------------------------------------------------------------------------------------------------

Function TimeGE(Time1, Time2)
! Tests for Time1 >= Time2.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time1 !} Times to be compared. They must both be of the same type.
  Type(Time_), Intent(In) :: Time2 !}
  ! Function result:
  Logical :: TimeGE ! Indicates whether Time1 >= Time2.

  TimeGE = Time2ShortTime(Time1) >= Time2ShortTime(Time2)

End Function TimeGE

!-------------------------------------------------------------------------------------------------------------

Function TimeGT(Time1, Time2)
! Tests for Time1 > Time2.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time1 !} Times to be compared. They must both be of the same type.
  Type(Time_), Intent(In) :: Time2 !}
  ! Function result:
  Logical :: TimeGT ! Indicates whether Time1 > Time2.

  TimeGT = .not. Time2 >= Time1

End Function TimeGT

!-------------------------------------------------------------------------------------------------------------

Function TimeLE(Time1, Time2)
! Tests for Time1 <= Time2.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time1 !} Times to be compared. They must both be of the same type.
  Type(Time_), Intent(In) :: Time2 !}
  ! Function result:
  Logical :: TimeLE ! Indicates whether Time1 <= Time2.

  TimeLE = Time2 >= Time1

End Function TimeLE

!-------------------------------------------------------------------------------------------------------------

Function TimeLT(Time1, Time2)
! Tests for Time1 < Time2.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time1 !} Times to be compared. They must both be of the same type.
  Type(Time_), Intent(In) :: Time2 !}
  ! Function result:
  Logical :: TimeLT ! Indicates whether Time1 < Time2.

  TimeLT = Time2 > Time1

End Function TimeLT

!-------------------------------------------------------------------------------------------------------------

Function AddTime(Time, dT)
! Adds two times.

! +infinity + -infinity and -infinity + +infinity are not allowed.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time !} Times to be added. dT must be a time interval. Both must be either model
  Type(Time_), Intent(In) :: dT   !} or clock times.
  ! Function result:
  Type(Time_) :: AddTime ! Time + dT.

  AddTime = ShortTime2Time(Time2ShortTime(Time) + Time2ShortTime(dT))

End Function AddTime

!-------------------------------------------------------------------------------------------------------------

Function SubtractTime(Time1, Time2)
! Subtracts one time from another.

! +infinity - +infinity or -infinity - -infinity are not allowed.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time1 !} Times to be subtracted. They must both be or both not be time intervals
  Type(Time_), Intent(In) :: Time2 !} or Time2 must be a time interval. Both must be either model or clock
                                   !} times.
  ! Function result:
  Type(Time_) :: SubtractTime ! Time1 - Time2.

  SubtractTime = ShortTime2Time(Time2ShortTime(Time1) - Time2ShortTime(Time2))

End Function SubtractTime

!-------------------------------------------------------------------------------------------------------------

Function MultiplyTime(dT, n)
! Multiplies a time by an integer.

! +/-infinity * 0 is not allowed.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: dT ! Time to be multiplied by n. Must be a time interval.
  Integer,     Intent(In) :: n  ! Quantity by which dT is to be multiplied.
  ! Function result:
  Type(Time_) :: MultiplyTime ! dT * n.

  MultiplyTime = ShortTime2Time(Time2ShortTime(dT) * n)

End Function MultiplyTime

!-------------------------------------------------------------------------------------------------------------

Function DivideTime(dT, n)
! Divides a time by an integer (rounded towards zero).

! Infinite times or zero denominators are not allowed.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: dT ! Time to be divided by n. Must be a time interval.
  Integer,     Intent(In) :: n  ! Quantity by which dT is to be divided.
  ! Function result:
  Type(Time_) :: DivideTime ! dT / n.

  DivideTime = ShortTime2Time(Time2ShortTime(dT) / n)

End Function DivideTime

!-------------------------------------------------------------------------------------------------------------

Function RatioOfTimes(dT1, dT2)
! Finds the ratio of two times (rounded towards zero).

! Infinite times or zero denominators are not allowed.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: dT1 !} Times for which the ratio is required. Must be time intervals. Both must
  Type(Time_), Intent(In) :: dT2 !} be either model or clock times.
  ! Function result:
  Integer :: RatioOfTimes ! dT1 / dT2.

  RatioOfTimes = Time2ShortTime(dT1) / Time2ShortTime(dT2)

End Function RatioOfTimes

!-------------------------------------------------------------------------------------------------------------

Function RoundTime(Time, T0, dT, Up)
! Rounds a time to a time of the form T0 + n*dT.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In)           :: Time ! Time. } Must be finite. Time and T0 must both be or both not be
  Type(Time_), Intent(In)           :: T0   ! T0.   } time intervals. dT must be non-zero and a time interval.
  Type(Time_), Intent(In)           :: dT   ! dT.   } All must be either model or clock times.
  Logical,     Intent(In), Optional :: Up   ! Indicates if rounding is up or down. If absent, rounds towards
                                            ! T0.
  ! Function result:
  Type(Time_) :: RoundTime ! Rounded time.

  RoundTime = ShortTime2Time(           &
                RoundShortTime(         &
                  Time2ShortTime(Time), &
                  Time2ShortTime(T0),   &
                  Time2ShortTime(dT),   &
                  Up                    &
                )                       &
              )

End Function RoundTime

!-------------------------------------------------------------------------------------------------------------

Function RoundTimeToYear(Time, Up, Strictly)
! Rounds a time to the start or end of the year.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time     ! Time. Must use year and month, must be finite and must not be a time
                                      ! interval.
  Logical,     Intent(In) :: Up       ! Indicates if rounding is up or down.
  Logical,     Intent(In) :: Strictly ! Ensures the rounded time doesn't equal time.
  ! Function result:
  Type(Time_) :: RoundTimeToYear ! Rounded time.

  If (.not. TimeUsesYearAndMonth(Time) .or. Time%Infinite /= 0 .or. IsTimeInterval(Time)) Then
    Call Message('UNEXPECTED FATAL ERROR in RoundTimeToYear', 4)
  End If

  RoundTimeToYear          = Time
  RoundTimeToYear%Month    = 1
  RoundTimeToYear%Day      = 1
  RoundTimeToYear%Hour     = 0
  RoundTimeToYear%Minute   = 0
  RoundTimeToYear%Second   = 0
  RoundTimeToYear%Fraction = 0

  ! If rounding up, ensure time <= required time. Then increase until time found.
  If (Up) Then

    If (IsClockTime(Time)) Then
      RoundTimeToYear%Year = RoundTimeToYear%Year - 1
    Else If (IsBackwards()) Then
      RoundTimeToYear%Year = RoundTimeToYear%Year + 1
    Else
      RoundTimeToYear%Year = RoundTimeToYear%Year - 1
    End If
    Call JustifyTime(RoundTimeToYear)
    Do
      If (RoundTimeToYear > Time .or. (RoundTimeToYear == Time .and. .not. Strictly)) Exit
      If (IsClockTime(Time)) Then
        RoundTimeToYear%Year = RoundTimeToYear%Year + 1
      Else If (IsBackwards()) Then
        RoundTimeToYear%Year = RoundTimeToYear%Year - 1
      Else
        RoundTimeToYear%Year = RoundTimeToYear%Year + 1
      End If
      Call JustifyTime(RoundTimeToYear)
    End Do

  ! If rounding down ensure time >= required time. Then reduce until time found.
  Else

    If (IsClockTime(Time)) Then
      RoundTimeToYear%Year = RoundTimeToYear%Year + 1
    Else If (IsBackwards()) Then
      RoundTimeToYear%Year = RoundTimeToYear%Year - 1
    Else
      RoundTimeToYear%Year = RoundTimeToYear%Year + 1
    End If
    Call JustifyTime(RoundTimeToYear)
    Do
      If (RoundTimeToYear < Time .or. (RoundTimeToYear == Time .and. .not. Strictly)) Exit
      If (IsClockTime(Time)) Then
        RoundTimeToYear%Year = RoundTimeToYear%Year - 1
      Else If (IsBackwards()) Then
        RoundTimeToYear%Year = RoundTimeToYear%Year + 1
      Else
        RoundTimeToYear%Year = RoundTimeToYear%Year - 1
      End If
      Call JustifyTime(RoundTimeToYear)
    End Do

  End If

End Function RoundTimeToYear

!-------------------------------------------------------------------------------------------------------------

Function RoundTimeToMonth(Time, Up, Strictly)
! Rounds a time to the start or end of the month.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time     ! Time. Must use year and month, must be finite and must not be a time
                                      ! interval.
  Logical,     Intent(In) :: Up       ! Indicates if rounding is up or down.
  Logical,     Intent(In) :: Strictly ! Ensures the rounded time doesn't equal time.
  ! Function result:
  Type(Time_) :: RoundTimeToMonth ! Rounded time.

  If (.not. TimeUsesYearAndMonth(Time) .or. Time%Infinite /= 0 .or. IsTimeInterval(Time)) Then
    Call Message('UNEXPECTED FATAL ERROR in RoundTimeToMonth', 4)
  End If

  RoundTimeToMonth          = Time
  RoundTimeToMonth%Day      = 1
  RoundTimeToMonth%Hour     = 0
  RoundTimeToMonth%Minute   = 0
  RoundTimeToMonth%Second   = 0
  RoundTimeToMonth%Fraction = 0

  ! If rounding up, ensure time <= required time. Then increase until time found.
  If (Up) Then

    If (IsClockTime(Time)) Then
      RoundTimeToMonth%Month = RoundTimeToMonth%Month - 1
    Else If (IsBackwards()) Then
      RoundTimeToMonth%Month = RoundTimeToMonth%Month + 1
    Else
      RoundTimeToMonth%Month = RoundTimeToMonth%Month - 1
    End If
    Call JustifyTime(RoundTimeToMonth)
    Do
      If (RoundTimeToMonth > Time .or. (RoundTimeToMonth == Time .and. .not. Strictly)) Exit
      If (IsClockTime(Time)) Then
        RoundTimeToMonth%Month = RoundTimeToMonth%Month + 1
      Else If (IsBackwards()) Then
        RoundTimeToMonth%Month = RoundTimeToMonth%Month - 1
      Else
        RoundTimeToMonth%Month = RoundTimeToMonth%Month + 1
      End If
      Call JustifyTime(RoundTimeToMonth)
    End Do

  ! If rounding down ensure time >= required time. Then reduce until time found.
  Else

    If (IsClockTime(Time)) Then
      RoundTimeToMonth%Month = RoundTimeToMonth%Month + 1
    Else If (IsBackwards()) Then
      RoundTimeToMonth%Month = RoundTimeToMonth%Month - 1
    Else
      RoundTimeToMonth%Month = RoundTimeToMonth%Month + 1
    End If
    Call JustifyTime(RoundTimeToMonth)
    Do
      If (RoundTimeToMonth < Time .or. (RoundTimeToMonth == Time .and. .not. Strictly)) Exit
      If (IsClockTime(Time)) Then
        RoundTimeToMonth%Month = RoundTimeToMonth%Month - 1
      Else If (IsBackwards()) Then
        RoundTimeToMonth%Month = RoundTimeToMonth%Month + 1
      Else
        RoundTimeToMonth%Month = RoundTimeToMonth%Month - 1
      End If
      Call JustifyTime(RoundTimeToMonth)
    End Do

  End If

End Function RoundTimeToMonth

!-------------------------------------------------------------------------------------------------------------

Function TimeMin(Time1, Time2)
! Returns the minimum of two times.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time1 !} The two times. They must both be of the same type.
  Type(Time_), Intent(In) :: Time2 !}
  ! Function result:
  Type(Time_) :: TimeMin ! Minimum of Time1 and Time2.

  If (Time2 >= Time1) Then
    TimeMin = Time1
  Else
    TimeMin = Time2
  End If

End Function TimeMin

!-------------------------------------------------------------------------------------------------------------

Function TimeMax(Time1, Time2)
! Returns the maximum of two times.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time1 !} The two times. They must both be of the same type.
  Type(Time_), Intent(In) :: Time2 !}
  ! Function result:
  Type(Time_) :: TimeMax ! Maximum of Time1 and Time2.

  If (Time2 >= Time1) Then
    TimeMax = Time2
  Else
    TimeMax = Time1
  End If

End Function TimeMax

!-------------------------------------------------------------------------------------------------------------

Function TimeIsInfFuture(Time)
! Finds out whether a time is in the infinite future (in the run direction for model times).

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time ! Time.
  ! Function result:
  Logical :: TimeIsInfFuture ! Indicates the time is in the infinite future.

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (IsClockTime(Time)) Then
    TimeIsInfFuture = Time%Infinite == 1
  Else If (IsBackwards()) Then
    TimeIsInfFuture = Time%Infinite == -1
  Else
    TimeIsInfFuture = Time%Infinite == 1
  End If

End Function TimeIsInfFuture

!-------------------------------------------------------------------------------------------------------------

Function TimeIsInfPast(Time)
! Finds out whether a time is in the infinite past (in the run direction for model times).

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time ! Time.
  ! Function result:
  Logical :: TimeIsInfPast ! Indicates the time is in the infinite past.

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (IsClockTime(Time)) Then
    TimeIsInfPast = Time%Infinite == -1
  Else If (IsBackwards()) Then
    TimeIsInfPast = Time%Infinite == 1
  Else
    TimeIsInfPast = Time%Infinite == -1
  End If

End Function TimeIsInfPast

!-------------------------------------------------------------------------------------------------------------

Function TimeIsTimeInterval(Time)
! Finds out whether a time is a time interval.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time ! Time.
  ! Function result:
  Logical :: TimeIsTimeInterval ! Indicates the time is a time interval.

  TimeIsTimeInterval = Time%Interval

End Function TimeIsTimeInterval

!-------------------------------------------------------------------------------------------------------------

Function TimeIsClockTime(Time)
! Finds out whether a time is a clock time.

! Note this routine can be used on partially constructed times provided the Clock part of the time is correct.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time ! Time.
  ! Function result:
  Logical :: TimeIsClockTime ! Indicates the time is a clock time.

  TimeIsClockTime = Time%Clock

End Function TimeIsClockTime

!-------------------------------------------------------------------------------------------------------------

Function InfFutureTime(Interval, Clock) Result(Time)
! Returns a time in the infinite future (in run direction for model times).

  Implicit None
  ! Argument list:
  Logical, Intent(In), Optional :: Interval ! If present and true, the time is a time interval.
  Logical, Intent(In), Optional :: Clock    ! If present and true, the time is a clock time.
  ! Function result:
  Type(Time_) :: Time ! The time in the infinite future.

  Time%Interval = .false.
  If (Present(Interval)) Then
    If (Interval) Time%Interval = .true.
  End If

  Time%Clock = .false.
  If (Present(Clock)) Then
    If (Clock) Time%Clock = .true.
  End If

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (IsClockTime(Time)) Then
    Time%Infinite = 1
  Else If (IsBackwards()) Then
    Time%Infinite = -1
  Else
    Time%Infinite = 1
  End If

End Function InfFutureTime

!-------------------------------------------------------------------------------------------------------------

Function InfPastTime(Interval, Clock) Result(Time)
! Returns a time in the infinite past (in run direction for model times).

  Implicit None
  ! Argument list:
  Logical, Intent(In), Optional :: Interval ! If present and true, the time is a time interval.
  Logical, Intent(In), Optional :: Clock    ! If present and true, the time is a clock time.
  ! Function result:
  Type(Time_) :: Time ! The time in the infinite past.

  Time%Interval = .false.
  If (Present(Interval)) Then
    If (Interval) Time%Interval = .true.
  End If

  Time%Clock = .false.
  If (Present(Clock)) Then
    If (Clock) Time%Clock = .true.
  End If

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (IsClockTime(Time)) Then
    Time%Infinite = -1
  Else If (IsBackwards()) Then
    Time%Infinite = 1
  Else
    Time%Infinite = -1
  End If

End Function InfPastTime

!-------------------------------------------------------------------------------------------------------------

Function ZeroTime(Clock) Result(Time)
! Returns a zero time interval.

  Implicit None
  ! Argument list:
  Logical, Intent(In), Optional :: Clock ! If present and true, the time is a clock time.
  ! Function result:
  Type(Time_) :: Time ! The zero time interval.

  Time%Interval = .true.

  Time%Clock = .false.
  If (Present(Clock)) Then
    If (Clock) Time%Clock = .true.
  End If

  Time%Infinite = 0
  Time%Day      = 0
  Time%Hour     = 0
  Time%Minute   = 0
  Time%Second   = 0
  Time%Fraction = 0

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (.not. IsClockTime(Time) .and. .not.Initialised) Call Message('UNEXPECTED FATAL ERROR in ZeroTime', 4)

End Function ZeroTime

!-------------------------------------------------------------------------------------------------------------

Function ReferenceTime(Clock) Result(Time)
! Returns a reference time (for when one wants a time variable set to a fixed time, but the time itself
! doesn't matter).

  Implicit None
  ! Argument list:
  Logical, Intent(In), Optional :: Clock ! If present and true, the time is a clock time.
  ! Function result:
  Type(Time_) :: Time ! The reference time.

  Time%Clock = .false.
  If (Present(Clock)) Then
    If (Clock) Time%Clock = .true.
  End If

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (IsClockTime(Time)) Then
    Time = Char2Time(CharTime = '1/1/2000 00:00', Clock = .true.)
  Else If (.not. IsCalendar()) Then
    Time = Char2Time('00:00')
  Else
    Time = Char2Time('1/1/2000 00:00')
  End If

End Function ReferenceTime

!-------------------------------------------------------------------------------------------------------------

Function Char2Time(CharTime, Interval, Clock, AllowWildCards, ErrorCode) Result(Time)
! Converts a character string to a time.

! The format of the character string is as follows:
! (1) finite times involving year and month (i.e. finite non-interval times in a calendar time frame):
!                1/5/2001 15:30{:20{.34}}{UTC{s01:00}}
! (2) finite times not involving year and month (i.e. finite time intervals or finite times in the relative
!     time frame):
!                {-}{5d}15:30{:20{.34}}
!     (non-descriptive format) or
!                {-}{5day}{15hr}{30min}{20{.34}sec}
!     (descriptive format),
! (3) infinite times:
!                {-}Infinity

! Here {...} indicates optional components; if omitted, the numerical parts of these optional components are
! taken as zero. s indicates a sign, + or -. For descriptive format, at least one of the optional components
! corresponding to 'day', 'hr', 'min' and 'sec' must be present.
!
! Extra spaces are permitted except within 'UTC', 'day', 'hr', 'min', 'sec' and within the numerical
! components (a numerical component being a consecutive sequence of the digits 0 to 9, but not including any
! sign, decimal point, colon etc). The only compulsory space is that after the year (for finite times
! involving year and month). Leading zeros are not significant in any of the numerical components except after
! the decimal point, or except where the value is zero when a zero must be present. There is no set number of
! digits for each numerical component so that, for example, a time interval of 1000s could be indicated by
! '0:0:1000' as well as by '0:16:40'. However values can only lie outside their appropriate range where we
! have finite times not involving year and month and where the value is the most significant non-zero element
! of the time (this condition could be relaxed a bit if desired - see comments in code below). 'UTC', 'd',
! 'day', 'hr', 'min' and 'sec' are case insensitive.
!
! The value following 'UTC' indicates the time zone in which the time is expressed, specified as hours and
! minutes ahead of UTC (so that time zones to the east of the Greenwich meridian will tend to have a positive
! value and time zones to the west will tend to have a negative value). The sign associated with a time zone
! applies to the whole time zone offset (hours and minutes). The time '1/5/2001 15:30 UTC+01:00' can be
! thought of as 1/5/2001 15:30 expressed in the time code 'UTC + 1 hour'. If there is nothing after UTC (or if
! UTC itself is absent) then the time must be expressed in UTC.
!
! For finite times not involving year and month, the optional number followed by 'd' indicates the number of
! days. Also any minus sign associated with such a time applies to the whole number, so that, e.g. '-1:30' is
! one and a half hours before the relative time frame reference time or a time interval of minus one and a
! half hours, and '-1d12:00' or '- 1day 12hr 0min' is minus one and a half days.
!
! Positive and negative infinite times are indicated by 'infinity' and '-infinity' (with 'infinity' being case
! insensitive and with optional spaces allowed between '-' and 'infinity').

  Implicit None
  ! Argument list:
  Character(*), Intent(In)            :: CharTime
  Logical,      Intent(In),  Optional :: Interval
  Logical,      Intent(In),  Optional :: Clock
  Logical,      Intent(In),  Optional :: AllowWildCards
  Integer,      Intent(Out), Optional :: ErrorCode
  ! CharTime       :: Character string to be converted. Character formats are as described above.
  ! Interval       :: If present and true, the time is a time interval.
  ! Clock          :: If present and true, the time is a clock time.
  ! AllowWildCards :: If present and true, indicates that wild cards (*) are allowed in CharTime for some
  !                   elements out of year, month, day, hour, minute and second (but not fraction) for
  !                   non-interval calendar-frame times, and out of day, hour, minute, and second (but not
  !                   fraction) for non-interval relative-frame times. Interval times and clock times are not
  !                   allowed if AllowWildCards is true. Any wild card elements must be the first so many
  !                   (possibly none or all) elements from the above list. In the output wild cards are
  !                   replaced by YearOrigin, MonthOrigin, DayOrigin, 0, 0, 0 for calendar-frame times and by
  !                   0, 0, 0, 0 for relative-frame times. However wild card values count as non-zero with
  !                   regard to the statement above that 'values can only lie outside their appropriate range
  !                   where we have finite times not involving year and month and where the value is the most
  !                   significant non-zero element of the time'.
  ! ErrorCode      :: Non-zero values indicate that an error has occurred. If not present, the routine stops
  !                   when there is an error. If an unexpected fatal error occurs the routine stops whether or
  !                   not ErrorCode is present. The error codes are as follows:
  !                     $$ list error codes here
  !                     $$ split error 1 into 2 (too long/blank)
  !                        and 3 into 4 (relative/interval, descriptive/non-descriptive)
  !                        or possibly reduce to just 3 error types or even 1: blank, too long and incorrect format
  ! Function result:
  Type(Time_) :: Time ! The converted character string.
  ! Locals:
  Integer                      :: Year            ! Year.
  Integer                      :: Month           ! Month.
  Integer                      :: Day             ! Day.
  Integer                      :: Hour            ! Hour.
  Integer                      :: Minute          ! Minute.
  Integer                      :: Second          ! Second.
  Integer                      :: Fraction        ! Fraction of a second in units defined by FracsPerSec.
  Integer                      :: ZoneHour        ! Time zone hour.
  Integer                      :: ZoneMinute      ! Time zone minute.
  Character(MaxCharLength)     :: CYear           ! Year as a character string.
  Character(MaxCharLength)     :: CMonth          ! Month as a character string.
  Character(MaxCharLength)     :: CDay            ! Day as a character string.
  Character(MaxCharLength)     :: CHour           ! Hour as a character string.
  Character(MaxCharLength)     :: CMinute         ! Minute as a character string.
  Character(MaxCharLength)     :: CSecond         ! Second as a character string.
  Character(MaxCharLength)     :: CFraction       ! Fraction as a character string.
  Character(MaxCharLength)     :: CZoneHour       ! ZoneHour as a character string.
  Character(MaxCharLength)     :: CZoneMinute     ! ZoneMinute as a character string.
  Logical                      :: SecondNext      !} Indicates that the next bit of CharTemp is the number of
  Logical                      :: FractionNext    !} seconds, the number of fractions, the string 'UTC', and time
  Logical                      :: UTCNext         !} zone information respectively.
  Logical                      :: ZoneNext        !}
  Logical                      :: Minus           ! Indicates a minus sign present.
  Integer                      :: i               ! Character index in CharTemp.
  Integer                      :: MostSig         ! Indicates, for finite times not involving year and month,
                                                  ! which element of the time is most significant (Day = 1,
                                                  ! Fraction = 5).
  Logical                      :: WildCards       ! Indicates, for finite times not involving year and month,
                                                  ! that wild cards are present.
  Character(MaxCharLength + 1) :: CharTemp        ! Temporary copy of CharTime for processing. Note the length is
                                                  ! declared as MaxCharLength + 1 to exclude the possibility that
                                                  ! any of the CharTemp(i + 1:) expressions might be undefined.
  Integer(I64)                 :: FractionI64     ! I64 version of Fraction.
  Integer                      :: ErrorInChar2Int ! Indicates an error occurred in Char2Int.

  ! $$ provide wrapper routine Char2Time and rename and make private this routine to prevent use of
  !    AllowWildCards except from within this module

  If (Present(ErrorCode)) ErrorCode = 0

  If (                                               &
    Len_Trim(AdjustL(CharTime)) > MaxCharLength .or. &
    Len_Trim(CharTime) == 0                          &
  ) Go To 1

  CharTemp = AdjustL(CharTime)

  Time%Interval = .false.
  If (Present(Interval)) Then
    If (Interval) Time%Interval = .true.
  End If

  Time%Clock = .false.
  If (Present(Clock)) Then
    If (Clock) Time%Clock = .true.
  End If

  If (Present(AllowWildCards)) Then
    If ((IsTimeInterval(Time) .or. IsClockTime(Time)) .and. AllowWildCards) Then
      Call Message('UNEXPECTED FATAL ERROR in Char2Time', 4)
    End If
  End If

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (.not. IsClockTime(Time) .and. .not.Initialised) Call Message('UNEXPECTED FATAL ERROR in Char2Time', 4)

  ! Infinite times.
  If (CharTemp .CIEq. 'infinity') Then

    ! Note this is designed so clock times can be processed before InitTimeModule is called.
    If (IsClockTime(Time)) Then
      Time%Infinite = 1
    Else If (IsBackwards() .and. IsTimeInterval(Time)) Then
      Time%Infinite = -1
    Else
      Time%Infinite = 1
    End If
    Return

  Else If (CharTemp(1:1) == '-') Then
    If (AdjustL(CharTemp(2:)) .CIEq. 'infinity') Then

      ! Note this is designed so clock times can be processed before InitTimeModule is called.
      If (IsClockTime(Time)) Then
        Time%Infinite = -1
      Else If (IsBackwards() .and. IsTimeInterval(Time)) Then
        Time%Infinite = 1
      Else
        Time%Infinite = -1
      End If
      Return

    End If
  End If

  Time%Infinite = 0

  ! Check interval/non-interval and time frame consistent with character format for finite times.
  i = Scan(CharTemp, '/')
  If (i /= 0) Then
    If (.not. UsesYearAndMonth(Time)) Go To 2
  Else
    If (UsesYearAndMonth(Time)) Go To 2
  End If

  ! Finite times involving year and month.
  If (UsesYearAndMonth(Time)) Then

    SecondNext   = .false.
    FractionNext = .false.
    UTCNext      = .false.
    ZoneNext     = .false.
    Minus        = .false.

    i = Scan(CharTemp, '/')
    If (i == 0 .or. i == 1) Go To 2
    CDay     = CharTemp(1:i - 1)
    CharTemp = AdjustL(CharTemp(i + 1:))

    i = Scan(CharTemp, '/')
    If (i == 0 .or. i == 1) Go To 2
    CMonth   = CharTemp(1:i - 1)
    CharTemp = AdjustL(CharTemp(i + 1:))

    i = Scan(CharTemp, ' ')
    If (i == 0 .or. i == 1) Go To 2
    CYear    = CharTemp(1:i - 1)
    CharTemp = AdjustL(CharTemp(i + 1:))

    i = Scan(CharTemp, ':')
    If (i == 0 .or. i == 1) Go To 2
    CHour    = CharTemp(1:i - 1)
    CharTemp = AdjustL(CharTemp(i + 1:))

    i = Scan(CharTemp, ':Uu')
    If (i == 1) Go To 2
    If (i > 1) Then
      CMinute  = CharTemp(1:i - 1)
      If (CharTemp(i:i) == ':') Then
        SecondNext = .true.
        CharTemp = AdjustL(CharTemp(i + 1:))
      Else If (CharTemp(i:i) .CIEq. 'U') Then
        UTCNext  = .true.
        CharTemp = AdjustL(CharTemp(i:))
      End If
    Else
      CMinute = CharTemp
    End If

    If (SecondNext) Then
      i = Scan(CharTemp, '.Uu')
      If (i == 1) Go To 2
      If (i > 1) Then
        CSecond  = CharTemp(1:i - 1)
        If (CharTemp(i:i) == '.') Then
          FractionNext = .true.
          CharTemp = AdjustL(CharTemp(i + 1:))
        Else If (CharTemp(i:i) .CIEq. 'U') Then
          UTCNext  = .true.
          CharTemp = AdjustL(CharTemp(i:))
        End If
      Else
        CSecond = CharTemp
      End If
    Else
      CSecond = '0'
    End If

    If (FractionNext) Then
      i = Scan(CharTemp, 'Uu')
      If (i == 1) Go To 2
      If (i > 1) Then
        CFraction = CharTemp(1:i - 1)
        UTCNext  = .true.
        CharTemp = AdjustL(CharTemp(i:))
      Else
        CFraction = CharTemp
      End If
    Else
      CFraction = '0'
    End If

    If (UTCNext) Then
      If (.not. (CharTemp(1:3) .CIEq. 'UTC')) Go To 2
      CharTemp = AdjustL(CharTemp(4:))
      If (CharTemp /= ' ') Then
        ZoneNext = .true.
      End If
    End If

    If (ZoneNext) Then
      If (CharTemp(1:1) == '-') Then
        Minus    = .true.
        CharTemp = AdjustL(CharTemp(2:))
      Else If (CharTemp(1:1) == '+') Then
        Minus    = .false.
        CharTemp = AdjustL(CharTemp(2:))
      Else
        Go To 2
      End If
      i = Scan(CharTemp, ':')
      If (i == 0 .or. i == 1) Go To 2
      CZoneHour   = CharTemp(1:i - 1)
      CharTemp    = AdjustL(CharTemp(i + 1:))
      CZoneMinute = CharTemp
    Else
      CZoneHour   = '0'
      CZoneMinute = '0'
    End If

    If (Present(AllowWildCards)) Then
      If (AllowWildCards) Then
        If (CYear == '*') Then
          CYear = Int2Char(YearOrigin)
          If (CMonth == '*') Then
            CMonth = Int2Char(MonthOrigin)
            If (CDay == '*') Then
              CDay = Int2Char(DayOrigin)
              If (CHour == '*') Then
                CHour = '0'
                If (CMinute == '*') Then
                  CMinute = '0'
                  If (CSecond == '*') Then
                    CSecond = '0'
                  End If
                End If
              End If
            End If
          End If
        End If
      End If
    End If

    ! $$ Trim(Day) etc should have non-zero length, no spaces and only
    ! numbers (no '.' to ensure '00:00.1' gives error). Also they, and in particular
    ! Fraction, shouldn't be too long.
    ! Note internal read error wont catch '01 rubbish' and might not catch 00.1, but we could alter Char2Int
    ! to trap some of these. Could perhaps alter format definition to allow times to be "too long", e.g. 001/012/00002014, but Fraction
    ! still needs special care.
    ! Need to check carefully what is trapped currently and what isn't before altering code.

    Year       = Char2Int(CYear,       ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Month      = Char2Int(CMonth,      ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Day        = Char2Int(CDay,        ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Hour       = Char2Int(CHour,       ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Minute     = Char2Int(CMinute,     ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Second     = Char2Int(CSecond,     ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Fraction   = Char2Int(CFraction,   ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    ZoneHour   = Char2Int(CZoneHour,   ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    ZoneMinute = Char2Int(CZoneMinute, ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10

    FractionI64 = Int(Fraction, I64) * Int(FracsPerSec, I64)
    FractionI64 = FractionI64 / Int(10**Len_Trim(CFraction), I64)

    ! Trap values outside appropriate range. Note this is not essential if no wild cards or if only ZoneMinute
    ! is outside range. Could remove error trap for these cases (or change to warning) to allow a wider range
    ! of time formats, although with greater possibilities for user error. Fraction out of range should be
    ! impossible.
    If      (Month < 1 .or. Month > 12                      ) Then
      Go To 4
    Else If (Day   < 1 .or. Day   > DaysInMonth(Year, Month, IsClockTime(Time))) Then
      Go To 5
    Else If (Hour        >= 24         ) Then
      Go To 6
    Else If (Minute      >= 60         ) Then
      Go To 7
    Else If (Second      >= 60         ) Then
      Go To 8
    Else If (FractionI64 >= FracsPerSec) Then
      Call Message('UNEXPECTED FATAL ERROR in Char2Time', 4)
    Else If (ZoneMinute  >= 60         ) Then
      Go To 9
    End If

    If (Minus) Then
      ZoneHour   = -ZoneHour
      ZoneMinute = -ZoneMinute
    End If

    Time = InitTime(                             &
             Interval   = IsTimeInterval(Time),  &
             Clock      = IsClockTime(Time),     &
             Infinite   = 0,                     &
             Year       = Year,                  &
             Month      = Month,                 &
             Day        = Day,                   &
             Hour       = Hour,                  &
             Minute     = Minute,                &
             Second     = Second,                &
             Fraction   = Int(FractionI64, I32), &
             ZoneHour   = ZoneHour,              &
             ZoneMinute = ZoneMinute             &
           )

  ! Finite times not involving year and month (non-descriptive format).
  Else If (Scan(CharTemp, ':') /= 0) Then

    SecondNext   = .false.
    FractionNext = .false.
    Minus        = .false.

    i = Scan(CharTemp, '-')
    If (i > 1) Go To 2
    If (i == 1) Then
      Minus    = .true.
      CharTemp = AdjustL(CharTemp(i + 1:))
    End If

    i = Scan(CharTemp, 'Dd')
    If (i == 1) Go To 2
    If (i > 1) Then
      CDay     = CharTemp(1:i - 1)
      CharTemp = AdjustL(CharTemp(i + 1:))
    Else
      CDay = '0'
    End If

    i = Scan(CharTemp, ':')
    If (i == 0 .or. i == 1) Go To 2
    CHour    = CharTemp(1:i - 1)
    CharTemp = AdjustL(CharTemp(i + 1:))

    i = Scan(CharTemp, ':')
    If (i == 1) Go To 2
    If (i > 1) Then
      CMinute    = CharTemp(1:i - 1)
      SecondNext = .true.
      CharTemp   = AdjustL(CharTemp(i + 1:))
    Else
      CMinute = CharTemp
    End If

    If (SecondNext) Then
      i = Scan(CharTemp, '.')
      If (i == 1) Go To 2
      If (i > 1) Then
        CSecond      = CharTemp(1:i - 1)
        FractionNext = .true.
        CharTemp     = AdjustL(CharTemp(i + 1:))
      Else
        CSecond = CharTemp
      End If
    Else
      CSecond = '0'
    End If

    If (FractionNext) Then
      CFraction = CharTemp
    Else
      CFraction = '0'
    End If

    WildCards = .false.
    If (Present(AllowWildCards)) Then
      If (AllowWildCards) Then
        If (CDay == '*') Then
          CDay = '0'
          WildCards = .true.
          If (CHour == '*') Then
            CHour = '0'
            WildCards = .true.
            If (CMinute == '*') Then
              CMinute = '0'
              WildCards = .true.
              If (CSecond == '*') Then
                CSecond = '0'
                WildCards = .true.
              End If
            End If
          End If
        End If
      End If
    End If

    Day      = Char2Int(CDay,      ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Hour     = Char2Int(CHour,     ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Minute   = Char2Int(CMinute,   ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Second   = Char2Int(CSecond,   ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Fraction = Char2Int(CFraction, ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10

    FractionI64 = Int(Fraction, I64) * Int(FracsPerSec, I64)
    FractionI64 = FractionI64 / Int(10**Len_Trim(CFraction), I64)

    ! Trap values outside appropriate range. Note this is not essential if no wild cards. Could remove error
    ! trap for these cases (or change to warning) to allow a wider range of time formats, although with
    ! greater possibilities for error. Fraction out of range should be impossible.
    If (WildCards) Then
      MostSig = 1
    Else If (Day         > 0) Then
      MostSig = 1
    Else If (Hour        > 0) Then
      MostSig = 2
    Else If (Minute      > 0) Then
      MostSig = 3
    Else If (Second      > 0) Then
      MostSig = 4
    Else If (FractionI64 > 0) Then
      MostSig = 5
    Else
      MostSig = 0
    End If
    If      (Hour        >= 24          .and. MostSig /= 2) Then
      Go To 6
    Else If (Minute      >= 60          .and. MostSig /= 3) Then
      Go To 7
    Else If (Second      >= 60          .and. MostSig /= 4) Then
      Go To 8
    Else If (FractionI64 >= FracsPerSec .and. MostSig /= 5) Then
      Call Message('UNEXPECTED FATAL ERROR in Char2Time', 4)
    End If

    ! Note this is designed so clock times can be processed before InitTimeModule is called.
    If (.not. IsClockTime(Time)) Then
      If (IsBackwards() .and. IsTimeInterval(Time)) Then
        Minus = .not.Minus
      End If
    End If

    If (Minus) Then
      Day         = -Day
      Hour        = -Hour
      Minute      = -Minute
      Second      = -Second
      FractionI64 = -FractionI64
    End If

    Time = InitTime(                          &
             Interval = IsTimeInterval(Time), &
             Clock    = IsClockTime(Time),    &
             Infinite = 0,                    &
             Day      = Day,                  &
             Hour     = Hour,                 &
             Minute   = Minute,               &
             Second   = Second,               &
             Fraction = Int(FractionI64, I32) &
           )

    ! Set day = 0 if wild cards present (if minus is set it may have been set to -1 by InitTime).
    If (WildCards) Time%Day = 0

  ! Finite times not involving year and month (descriptive format).
  Else

    Minus = .false.

    If (Scan(CharTemp, 'DdHhMmSs') == 0) Go To 2

    i = Scan(CharTemp, '-')
    If (i > 1) Go To 2
    If (i == 1) Then
      Minus    = .true.
      CharTemp = AdjustL(CharTemp(i + 1:))
    End If

    i = Scan(CharTemp, 'Dd')
    If (i == 1) Go To 2
    If (i > 1) Then
      If (.not. (CharTemp(i:i + 2) .CIEq. 'day')) Go To 2
      CDay     = CharTemp(1:i - 1)
      CharTemp = AdjustL(CharTemp(i + 3:))
    Else
      CDay = '0'
    End If

    i = Scan(CharTemp, 'Hh')
    If (i == 1) Go To 2
    If (i > 1) Then
      If (.not. (CharTemp(i:i + 1) .CIEq. 'hr')) Go To 2
      CHour    = CharTemp(1:i - 1)
      CharTemp = AdjustL(CharTemp(i + 2:))
    Else
      CHour = '0'
    End If

    i = Scan(CharTemp, 'Mm')
    If (i == 1) Go To 2
    If (i > 1) Then
      If (.not. (CharTemp(i:i + 2) .CIEq. 'min')) Go To 2
      CMinute  = CharTemp(1:i - 1)
      CharTemp = AdjustL(CharTemp(i + 3:))
    Else
      CMinute = '0'
    End If

    i = Scan(CharTemp, 'Ss')
    If (i == 1) Go To 2
    If (i > 1) Then
      If (.not. (CharTemp(i:i + 2) .CIEq. 'sec')) Go To 2
      CSecond  = CharTemp(1:i - 1)
      CharTemp = AdjustL(CharTemp(i + 3:))
    Else
      CSecond = '0'
    End If

    If (CharTemp /= ' ') Go To 2

    i = Scan(CSecond, '.')
    If (i == 1) Go To 2
    If (i > 1) Then
      CFraction = AdjustL(CSecond(i + 1:))
      CSecond   = CSecond(1:i - 1)
    Else
      CFraction = '0'
    End If

    WildCards = .false.
    If (Present(AllowWildCards)) Then
      If (AllowWildCards) Then
        If (CDay == '*') Then
          CDay = '0'
          WildCards = .true.
          If (CHour == '*') Then
            CHour = '0'
            WildCards = .true.
            If (CMinute == '*') Then
              CMinute = '0'
              WildCards = .true.
              If (CSecond == '*') Then
                CSecond = '0'
                WildCards = .true.
              End If
            End If
          End If
        End If
      End If
    End If

    Day      = Char2Int(CDay,      ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Hour     = Char2Int(CHour,     ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Minute   = Char2Int(CMinute,   ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Second   = Char2Int(CSecond,   ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10
    Fraction = Char2Int(CFraction, ErrorInChar2Int); If (ErrorInChar2Int /= 0) Go To 10

    FractionI64 = Int(Fraction, I64) * Int(FracsPerSec, I64)
    FractionI64 = FractionI64 / Int(10**Len_Trim(CFraction), I64)

    ! Trap values outside appropriate range. Note this is not essential if no wild cards. Could remove error
    ! trap for these cases (or change to warning) to allow a wider range of time formats, although with
    ! greater possibilities for error. Fraction out of range should be impossible.
    If (WildCards) Then
      MostSig = 1
    Else If (Day         > 0) Then
      MostSig = 1
    Else If (Hour        > 0) Then
      MostSig = 2
    Else If (Minute      > 0) Then
      MostSig = 3
    Else If (Second      > 0) Then
      MostSig = 4
    Else If (FractionI64 > 0) Then
      MostSig = 5
    Else
      MostSig = 0
    End If
    If      (Hour        >= 24          .and. MostSig /= 2) Then
      Go To 6
    Else If (Minute      >= 60          .and. MostSig /= 3) Then
      Go To 7
    Else If (Second      >= 60          .and. MostSig /= 4) Then
      Go To 8
    Else If (FractionI64 >= FracsPerSec .and. MostSig /= 5) Then
      Call Message('UNEXPECTED FATAL ERROR in Char2Time', 4)
    End If

    ! Note this is designed so clock times can be processed before InitTimeModule is called.
    If (.not. IsClockTime(Time)) Then
      If (IsBackwards() .and. IsTimeInterval(Time)) Then
        Minus = .not.Minus
      End If
    End If

    If (Minus) Then
      Day         = -Day
      Hour        = -Hour
      Minute      = -Minute
      Second      = -Second
      FractionI64 = -FractionI64
    End If

    Time = InitTime(                          &
             Interval = IsTimeInterval(Time), &
             Clock    = IsClockTime(Time),    &
             Infinite = 0,                    &
             Day      = Day,                  &
             Hour     = Hour,                 &
             Minute   = Minute,               &
             Second   = Second,               &
             Fraction = Int(FractionI64, I32) &
           )

    ! Set day = 0 if wild cards present (if minus is set it may have been set to -1 by InitTime).
    If (WildCards) Time%Day = 0

  End If

  Return

1 Continue
  If (Present(ErrorCode)) Then
    ErrorCode = 1
    Return
  Else
    Call Message(                       &
           'FATAL ERROR: Time "'     // &
           Trim(CharTime)            // &
           '" is too long or blank',    &
           3                            &
         )
  End If

2 Continue
  If (UsesYearAndMonth(Time)) Then
    If (Present(ErrorCode)) Then
      ErrorCode = 2
      Return
    Else
      Call Message(                                                                           &
             'FATAL ERROR: Time "'                                                         // &
             Trim(CharTime)                                                                // &
             '" has an incorrect format for a non-interval time in a calendar time frame',    &
             3                                                                                &
           )
    End If
  Else
    If (Present(ErrorCode)) Then
      ErrorCode = 3
      Return
    Else
      Call Message(                                                                                   &
             'FATAL ERROR: Time "'                                                                 // &
             Trim(CharTime)                                                                        // &
             '" has an incorrect format for a time interval or a time in the relative time frame',    &
             3                                                                                        &
           )
    End If
  End If

4 Continue
  If (Present(ErrorCode)) Then
    ErrorCode = 4
    Return
  Else
    Call Message(                         &
           'FATAL ERROR: Time "'       // &
           Trim(CharTime)              // &
           '" has month out of range',    &
           3                              &
         )
  End If

5 Continue
  If (Present(ErrorCode)) Then
    ErrorCode = 5
    Return
  Else
    Call Message(                       &
           'FATAL ERROR: Time "'     // &
           Trim(CharTime)            // &
           '" has day out of range',    &
           3                            &
         )
  End If

6 Continue
  If (Present(ErrorCode)) Then
    ErrorCode = 6
    Return
  Else
    Call Message(                        &
           'FATAL ERROR: Time "'      // &
           Trim(CharTime)             // &
           '" has hour out of range',    &
           3                             &
         )
  End If

7 Continue
  If (Present(ErrorCode)) Then
    ErrorCode = 7
    Return
  Else
    Call Message(                          &
           'FATAL ERROR: Time "'        // &
           Trim(CharTime)               // &
           '" has minute out of range',    &
           3                               &
         )
  End If

8 Continue
  If (Present(ErrorCode)) Then
    ErrorCode = 8
    Return
  Else
    Call Message(                          &
           'FATAL ERROR: Time "'        // &
           Trim(CharTime)               // &
           '" has second out of range',    &
           3                               &
         )
  End If

9 Continue
  If (Present(ErrorCode)) Then
    ErrorCode = 9
    Return
  Else
    Call Message(                                    &
           'FATAL ERROR: Time "'                  // &
           Trim(CharTime)                         // &
           '" has time zone minute out of range',    &
           3                                         &
         )
  End If

10 Continue
  If (Present(ErrorCode)) Then
    ErrorCode = 10
    Return
  Else
    Call Message(                        &
           'FATAL ERROR: Time "'      // &
           Trim(CharTime)             // &
           '" has an invalid format',    &
           3                             &
         )
  End If

End Function Char2Time

!-------------------------------------------------------------------------------------------------------------

Function Time2Char(Time, Seconds, DecimalPlaces, Describe, NameII, NameIITS)
! Converts a time to a character string.

! The format of the character string is as described in Char2Time.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In)           :: Time
  Logical,     Intent(In)           :: Seconds
  Integer,     Intent(In)           :: DecimalPlaces
  Logical,     Intent(In)           :: Describe
  Logical,     Intent(In), Optional :: NameII
  Logical,     Intent(In), Optional :: NameIITS
  ! Time          :: Time to be converted.
  ! Seconds       :: Indicates seconds are to be included in result.
  ! DecimalPlaces :: Number of decimal places to be included for seconds in result.
  ! Describe      :: For finite times not involving year and month, indicates that the time is to be formatted
  !                  in the descriptive form (e.g. '5hr 30min' or '5day 12hr 30min 30.55sec', not '05:30' or
  !                  '5d12:30:30.55').
  ! NameII        :: Indicates Name II format for finite times involving year and month.
  ! NameIITS      :: Indicates Name II time series format for finite times involving year and month.
  ! Function result:
  Character(MaxCharLength) :: Time2Char ! The converted time.
  ! Locals:
  Character(MaxCharLength) :: FormatString ! Format string for fractions of a second.
  Integer                  :: Day          !} Values from Time for finite times not involving year and month.
  Integer                  :: Hour         !} Used to convert to output format (where sign applies to all
  Integer                  :: Minute       !} elements) from internal format (where sign applies only to Day).
  Integer                  :: Second       !}
  Integer                  :: Fraction     !}
  Character(1)             :: Sign         ! Sign (' ' or '-') for finite times not involving year and month.
  Type(Time_)              :: TempTime     ! Copy of Time with time zone changed to UTC.
  Type(Time_)              :: TimeForZone  ! A time with UTC time zone.

  If (Time%Infinite == 1) Then

    ! Note this is designed so clock times can be processed before InitTimeModule is called.
    If (IsClockTime(Time)) Then
      Time2Char = 'infinity'
    Else If (IsBackwards() .and. IsTimeInterval(Time)) Then
      Time2Char = '-infinity'
    Else
      Time2Char = 'infinity'
    End If

  Else If (Time%Infinite == -1) Then

    ! Note this is designed so clock times can be processed before InitTimeModule is called.
    If (IsClockTime(Time)) Then
      Time2Char = '-infinity'
    Else If (IsBackwards() .and. IsTimeInterval(Time)) Then
      Time2Char = 'infinity'
    Else
      Time2Char = '-infinity'
    End If

  Else If (UsesYearAndMonth(Time)) Then

    If (Present(NameII)) Then
      If (NameII) Then

        TimeForZone%Interval   = Time%Interval
        TimeForZone%Clock      = Time%Clock
        TimeForZone%Infinite   = 0
        TimeForZone%ZoneHour   = 0
        TimeForZone%ZoneMinute = 0
        TempTime = ChangeTimeZone(Time, TimeForZone)

        Time2Char =                                                   &
          Trim(Int2Char(TempTime%Hour,   2, .false., 'L', 'i2.2')) // &
          Trim(Int2Char(TempTime%Minute, 2, .false., 'L', 'i2.2')) // &
          'UTC '                                                   // &
          Trim(Int2Char(TempTime%Day,    2, .false., 'L', 'i2.2')) // &
          '/'                                                      // &
          Trim(Int2Char(TempTime%Month,  2, .false., 'L', 'i2.2')) // &
          '/'                                                      // &
          Trim(Int2Char(TempTime%Year,   4, .false., 'L', 'i4.4'))

        Return
      End If
    End If

    If (Present(NameIITS)) Then
      If (NameIITS) Then

        TimeForZone%Interval   = Time%Interval
        TimeForZone%Clock      = Time%Clock
        TimeForZone%Infinite   = 0
        TimeForZone%ZoneHour   = 0
        TimeForZone%ZoneMinute = 0
        TempTime = ChangeTimeZone(Time, TimeForZone)

        Time2Char =                                                   &
          Trim(Int2Char(TempTime%Day,    2, .false., 'L', 'i2.2')) // &
          '/'                                                      // &
          Trim(Int2Char(TempTime%Month,  2, .false., 'L', 'i2.2')) // &
          '/'                                                      // &
          Trim(Int2Char(TempTime%Year,   4, .false., 'L', 'i4.4')) // &
          ',     '                                                 // &
          Trim(Int2Char(TempTime%Hour,   2, .false., 'L', 'i2.2')) // &
          ':'                                                      // &
          Trim(Int2Char(TempTime%Minute, 2, .false., 'L', 'i2.2')) // &
          ':'                                                      // &
          Trim(Int2Char(TempTime%Second, 2, .false., 'L', 'i2.2'))

        Return
      End If
    End If

    If (Seconds .and. DecimalPlaces > 0) Then

      FormatString = 'i'                           // &
                     Trim(Int2Char(DecimalPlaces)) // &
                     '.'                           // &
                     Trim(Int2Char(DecimalPlaces))

      Time2Char =                                                                                            &
        Trim(Int2Char(Time%Day,    2, .false., 'L', 'i2.2'))                                              // &
        '/'                                                                                               // &
        Trim(Int2Char(Time%Month,  2, .false., 'L', 'i2.2'))                                              // &
        '/'                                                                                               // &
        Trim(Int2Char(Time%Year,   4, .false., 'L', 'i4.4'))                                              // &
        ' '                                                                                               // &
        Trim(Int2Char(Time%Hour,   2, .false., 'L', 'i2.2'))                                              // &
        ':'                                                                                               // &
        Trim(Int2Char(Time%Minute, 2, .false., 'L', 'i2.2'))                                              // &
        ':'                                                                                               // &
        Trim(Int2Char(Time%Second, 2, .false., 'L', 'i2.2'))                                              // &
        '.'                                                                                               // &
        Trim(                                                                                                &
          Int2Char(                                                                                          &
            Int(Real(Time%Fraction, P64) * Real(10.0, P64)**DecimalPlaces / Real(FracsPerSec, P64), I32),    &
            DecimalPlaces, .false., 'L', FormatString                                                        &
          )                                                                                                  &
        )                                                                                                 // &
        ' UTC'

    Else If (Seconds) Then

      Time2Char =                                               &
        Trim(Int2Char(Time%Day,    2, .false., 'L', 'i2.2')) // &
        '/'                                                  // &
        Trim(Int2Char(Time%Month,  2, .false., 'L', 'i2.2')) // &
        '/'                                                  // &
        Trim(Int2Char(Time%Year,   4, .false., 'L', 'i4.4')) // &
        ' '                                                  // &
        Trim(Int2Char(Time%Hour,   2, .false., 'L', 'i2.2')) // &
        ':'                                                  // &
        Trim(Int2Char(Time%Minute, 2, .false., 'L', 'i2.2')) // &
        ':'                                                  // &
        Trim(Int2Char(Time%Second, 2, .false., 'L', 'i2.2')) // &
        ' UTC'

    Else

      Time2Char =                                               &
        Trim(Int2Char(Time%Day,    2, .false., 'L', 'i2.2')) // &
        '/'                                                  // &
        Trim(Int2Char(Time%Month,  2, .false., 'L', 'i2.2')) // &
        '/'                                                  // &
        Trim(Int2Char(Time%Year,   4, .false., 'L', 'i4.4')) // &
        ' '                                                  // &
        Trim(Int2Char(Time%Hour,   2, .false., 'L', 'i2.2')) // &
        ':'                                                  // &
        Trim(Int2Char(Time%Minute, 2, .false., 'L', 'i2.2')) // &
        ' UTC'

    End If

    If (Time%ZoneHour /= 0 .or. Time%ZoneMinute /= 0) Then
      If (Time%ZoneHour < 0 .or. Time%ZoneMinute < 0) Then
        Time2Char =                                                    &
          Trim(Time2Char)                                           // &
          '-'                                                       // &
          Trim(Int2Char(-Time%ZoneHour,   2, .false., 'L', 'i2.2')) // &
          ':'                                                       // &
          Trim(Int2Char(-Time%ZoneMinute, 2, .false., 'L', 'i2.2'))
      Else
        Time2Char =                                                   &
          Trim(Time2Char)                                          // &
          '+'                                                      // &
          Trim(Int2Char(Time%ZoneHour,   2, .false., 'L', 'i2.2')) // &
          ':'                                                      // &
          Trim(Int2Char(Time%ZoneMinute, 2, .false., 'L', 'i2.2'))
      End If
    End If

  Else

    Day      = Time%Day
    Hour     = Time%Hour
    Minute   = Time%Minute
    Second   = Time%Second
    Fraction = Time%Fraction
    Sign     = ' '

    If (Day < 0) Then
      If (Fraction > 0) Then
        Fraction = Fraction - FracsPerSec
        Second   = Second + 1
      End If
      If (Second > 0) Then
        Second = Second - 60
        Minute = Minute + 1
      End If
      If (Minute > 0) Then
        Minute = Minute - 60
        Hour   = Hour + 1
      End If
      If (Hour > 0) Then
        Hour = Hour - 24
        Day  = Day + 1
      End If
      Day      = -Day
      Hour     = -Hour
      Minute   = -Minute
      Second   = -Second
      Fraction = -Fraction
      Sign     = '-'
    End If

    ! Note this is designed so clock times can be processed before InitTimeModule is called.
    If (.not. IsClockTime(Time)) Then
      If (IsBackwards() .and. IsTimeInterval(Time)) Then
        If (Sign == ' ') Then
          Sign = '-'
        Else
          Sign = ' '
        End If
      End If
    End If

    If (Describe) Then

      If (Seconds .and. DecimalPlaces > 0) Then

        FormatString = 'i'                           // &
                       Trim(Int2Char(DecimalPlaces)) // &
                       '.'                           // &
                       Trim(Int2Char(DecimalPlaces))

        Time2Char = Sign
        If (Day /= 0) Then
          Time2Char = Trim(Time2Char) // ' ' // Trim(Int2Char(Day   )) // 'day'
        End If
        If (Hour /= 0 .or. Day /= 0) Then
          Time2Char = Trim(Time2Char) // ' ' // Trim(Int2Char(Hour  )) // 'hr'
        End If
        If (Minute /= 0 .or. Hour /= 0 .or. Day /= 0) Then
          Time2Char = Trim(Time2Char) // ' ' // Trim(Int2Char(Minute)) // 'min'
        End If
        Time2Char =                                                                                       &
          Trim(Time2Char)                                                                              // &
          ' '                                                                                          // &
          Trim(Int2Char(Second))                                                                       // &
          '.'                                                                                          // &
          Trim(                                                                                           &
            Int2Char(                                                                                     &
              Int(Real(Fraction, P64) * Real(10.0, P64)**DecimalPlaces / Real(FracsPerSec, P64), I32),    &
              DecimalPlaces, .false., 'L', FormatString                                                   &
            )                                                                                             &
          )                                                                                            // &
          'sec'
        Time2Char = AdjustL(Time2Char)

      Else If (Seconds) Then

        Time2Char = Sign
        If (Day /= 0) Then
          Time2Char = Trim(Time2Char) // ' ' // Trim(Int2Char(Day   )) // 'day'
        End If
        If (Hour /= 0 .or. Day /= 0) Then
          Time2Char = Trim(Time2Char) // ' ' // Trim(Int2Char(Hour  )) // 'hr'
        End If
        If (Minute /= 0 .or. Hour /= 0 .or. Day /= 0) Then
          Time2Char = Trim(Time2Char) // ' ' // Trim(Int2Char(Minute)) // 'min'
        End If
        Time2Char = Trim(Time2Char) // ' ' // Trim(Int2Char(Second)) // 'sec'
        Time2Char = AdjustL(Time2Char)

      Else

        Time2Char = Sign
        If (Day /= 0) Then
          Time2Char = Trim(Time2Char) // ' ' // Trim(Int2Char(Day   )) // 'day'
        End If
        If (Hour /= 0 .or. Day /= 0) Then
          Time2Char = Trim(Time2Char) // ' ' // Trim(Int2Char(Hour  )) // 'hr'
        End If
        Time2Char = Trim(Time2Char) // ' ' // Trim(Int2Char(Minute)) // 'min'
        Time2Char = AdjustL(Time2Char)

      End If

    Else

      If (Seconds .and. DecimalPlaces > 0) Then

        FormatString = 'i'                           // &
                       Trim(Int2Char(DecimalPlaces)) // &
                       '.'                           // &
                       Trim(Int2Char(DecimalPlaces))

        Time2Char = Sign
        If (Day /= 0) Then
          Time2Char = Trim(Time2Char) // Trim(Int2Char(Day)) // 'd'
        End If
        Time2Char =                                                                                       &
          Trim(Time2Char)                                                                              // &
          Trim(Int2Char(Hour,   2, .false., 'L', 'i2.2'))                                              // &
          ':'                                                                                          // &
          Trim(Int2Char(Minute, 2, .false., 'L', 'i2.2'))                                              // &
          ':'                                                                                          // &
          Trim(Int2Char(Second, 2, .false., 'L', 'i2.2'))                                              // &
          '.'                                                                                          // &
          Trim(                                                                                           &
            Int2Char(                                                                                     &
              Int(Real(Fraction, P64) * Real(10.0, P64)**DecimalPlaces / Real(FracsPerSec, P64), I32),    &
              DecimalPlaces, .false., 'L', FormatString                                                   &
            )                                                                                             &
          )
        Time2Char = AdjustL(Time2Char)

      Else If (Seconds) Then

        Time2Char = Sign
        If (Day /= 0) Then
          Time2Char = Trim(Time2Char) // Trim(Int2Char(Day)) // 'd'
        End If
        Time2Char =                                          &
          Trim(Time2Char)                                 // &
          Trim(Int2Char(Hour,   2, .false., 'L', 'i2.2')) // &
          ':'                                             // &
          Trim(Int2Char(Minute, 2, .false., 'L', 'i2.2')) // &
          ':'                                             // &
          Trim(Int2Char(Second, 2, .false., 'L', 'i2.2'))
        Time2Char = AdjustL(Time2Char)

      Else

        Time2Char = Sign
        If (Day /= 0) Then
          Time2Char = Trim(Time2Char) // Trim(Int2Char(Day)) // 'd'
        End If
        Time2Char =                                          &
          Trim(Time2Char)                                 // &
          Trim(Int2Char(Hour,   2, .false., 'L', 'i2.2')) // &
          ':'                                             // &
          Trim(Int2Char(Minute, 2, .false., 'L', 'i2.2'))
        Time2Char = AdjustL(Time2Char)

      End If

    End If

  End If

End Function Time2Char

!-------------------------------------------------------------------------------------------------------------

Function FileNameTime(Time)
! Converts a time to a character string in a form suitable for use in a file name.

! The format of the character string is as follows:
! (1) finite times involving year and month (i.e. finite non-interval times in a calendar time frame):
!                YYYYMMDDHHMM
! (2) finite times not involving year and month (i.e. finite time intervals or finite times in the relative
!     time frame):
!                {-}DayHHMM
! (3) infinite times:
!                {-}Infinity

! Here {} indicates optional components. In (2) Day has a variable length and the sign applies to Day, HH and
! MM.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time ! Time to be converted.
  ! Function result:
  Character(MaxCharLength) :: FileNameTime ! The converted time.
  ! Locals:
  Integer      :: Day      !} Values from Time for finite times not involving year and month. Used to convert
  Integer      :: Hour     !} to output format (where sign applies to all elements) from internal format
  Integer      :: Minute   !} (where sign applies only to Day).
  Integer      :: Second   !}
  Integer      :: Fraction !}
  Character(1) :: Sign     ! Sign (' ' or '-') for finite times not involving year and month.

  If (Time%Infinite == 1) Then

    ! Note this is designed so clock times can be processed before InitTimeModule is called.
    If (IsClockTime(Time)) Then
      FileNameTime = 'Infinity'
    Else If (IsBackwards() .and. IsTimeInterval(Time)) Then
      FileNameTime = '-Infinity'
    Else
      FileNameTime = 'Infinity'
    End If

  Else If (Time%Infinite == -1) Then

    ! Note this is designed so clock times can be processed before InitTimeModule is called.
    If (IsClockTime(Time)) Then
      FileNameTime = '-Infinity'
    Else If (IsBackwards() .and. IsTimeInterval(Time)) Then
      FileNameTime = 'Infinity'
    Else
      FileNameTime = '-Infinity'
    End If

  Else If (UsesYearAndMonth(Time)) Then

    FileNameTime =                                            &
      Trim(Int2Char(Time%Year,   4, .false., 'L', 'i4.4')) // &
      Trim(Int2Char(Time%Month,  2, .false., 'L', 'i2.2')) // &
      Trim(Int2Char(Time%Day,    2, .false., 'L', 'i2.2')) // &
      Trim(Int2Char(Time%Hour,   2, .false., 'L', 'i2.2')) // &
      Trim(Int2Char(Time%Minute, 2, .false., 'L', 'i2.2'))

  Else

    Day      = Time%Day
    Hour     = Time%Hour
    Minute   = Time%Minute
    Second   = Time%Second
    Fraction = Time%Fraction
    Sign     = ' '

    If (Day < 0) Then
      If (Fraction > 0) Then
        Fraction = Fraction - FracsPerSec
        Second   = Second + 1
      End If
      If (Second > 0) Then
        Second = Second - 60
        Minute = Minute + 1
      End If
      If (Minute > 0) Then
        Minute = Minute - 60
        Hour   = Hour + 1
      End If
      If (Hour > 0) Then
        Hour = Hour - 24
        Day  = Day + 1
      End If
      Day      = -Day
      Hour     = -Hour
      Minute   = -Minute
      Second   = -Second
      Fraction = -Fraction
      Sign     = '-'
    End If

    ! Note this is designed so clock times can be processed before InitTimeModule is called.
    If (.not. IsClockTime(Time)) Then
      If (IsBackwards() .and. IsTimeInterval(Time)) Then
        If (Sign == ' ') Then
          Sign = '-'
        Else
          Sign = ' '
        End If
      End If
    End If

    FileNameTime =                                       &
      Trim(AdjustL(Sign // Int2Char(Day)))            // &
      Trim(Int2Char(Hour,   2, .false., 'L', 'i2.2')) // &
      Trim(Int2Char(Minute, 2, .false., 'L', 'i2.2'))

  End If

End Function FileNameTime

!-------------------------------------------------------------------------------------------------------------

Function ChangeTimeZone(Time, TimeForZone) Result(ChangedTime)
! Changes the time zone of a time, while keeping the time the same.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time        ! Time.
  Type(Time_), Intent(In) :: TimeForZone ! A time giving the required time zone (must be finite).
  ! Function result:
  Type(Time_) :: ChangedTime ! Time with changed time zone.

  If (                                                               &
    ( IsTimeInterval(Time) .neqv. IsTimeInterval(TimeForZone) ) .or. &
    ( IsClockTime   (Time) .neqv. IsClockTime   (TimeForZone) )      &
  ) Then
    Call Message('UNEXPECTED FATAL ERROR in ChangeTimeZone', 4)
  End If

  If (TimeForZone%Infinite /= 0) Call Message('UNEXPECTED FATAL ERROR in ChangeTimeZone', 4)

  If (UsesYearAndMonth(Time)) Then
    If (                                              &
      Time%ZoneHour   == TimeForZone%ZoneHour   .and. &
      Time%ZoneMinute == TimeForZone%ZoneMinute       &
    ) Then
      ChangedTime = Time
    Else
      ChangedTime            = Time
      ChangedTime%Hour       = ChangedTime%Hour -     &
                               ChangedTime%ZoneHour + &
                               TimeForZone%ZoneHour
      ChangedTime%Minute     = ChangedTime%Minute -     &
                               ChangedTime%ZoneMinute + &
                               TimeForZone%ZoneMinute
      ChangedTime%ZoneHour   = TimeForZone%ZoneHour
      ChangedTime%ZoneMinute = TimeForZone%ZoneMinute
      Call JustifyTime(ChangedTime)
    End If
  Else
    ChangedTime = Time
  End If

End Function ChangeTimeZone

!-------------------------------------------------------------------------------------------------------------

Function DayOfWeek(Time)
! Evaluates the day of the week (Monday = 1, Sunday = 7) for finite calendar-frame times and returns zero for
! relative-frame times or infinite times.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time ! Time. Must not be a time interval.
  ! Function result:
  Integer :: DayOfWeek ! Day of the week (Monday = 1, Sunday = 7).
  ! Locals:
  Type(Time_) :: TempTime ! Temporary variable. TempTime%Day is days since end of 1999 so that 1/1/2000 -
                          ! which is a Saturday - has TempTime%Day = 1.

  If (IsTimeInterval(Time)) Call Message('UNEXPECTED FATAL ERROR in DayOfWeek', 4)

  If (UsesYearAndMonth(Time)) Then
    TempTime = Time
    Call ConvertTime(Year = 2000, Month = 1, Time = TempTime)
    DayOfWeek = Mod(TempTime%Day + 5, 7)
    If (DayOfWeek <= 0) DayOfWeek = DayOfWeek + 7
  Else
    DayOfWeek = 0
  End If

End Function DayOfWeek

!-------------------------------------------------------------------------------------------------------------

Function JulianDayNumber(Time)
! Evaluates the Julian day number for finite calendar-frame times and returns zero for relative-frame times or
! infinite times.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time ! Time. Must not be a time interval.
  ! Function result:
  Integer :: JulianDayNumber ! Julian day number.
  ! Locals:
  Type(Time_) :: TempTime ! Temporary variable. TempTime%Day is days since end of previous year so that
                          ! 1/1/Time%Year has TempTime%Day = 1.

  If (IsTimeInterval(Time)) Call Message('UNEXPECTED FATAL ERROR in DayOfWeek', 4)

  If (UsesYearAndMonth(Time)) Then
    TempTime = Time
    Call ConvertTime(Year = Time%Year, Month = 1, Time = TempTime)
    JulianDayNumber = TempTime%Day
  Else
    JulianDayNumber = 0
  End If

End Function JulianDayNumber

!-------------------------------------------------------------------------------------------------------------

Function CurrentClockTime()
! Returns the current clock time.

  Implicit None
  ! Function result:
  Type(Time_) :: CurrentClockTime ! Current clock time.
  ! Locals:
  Integer      :: Values(8) ! Values returned by Date_And_Time.
  Integer(I64) :: Fraction  ! Fractions of a second.

  Call Date_And_Time(Values = Values)

  ! Values(8) is in milli-seconds - convert to fractions as defined by FracsPerSec.
  Fraction = Int(Values(8), I64) * Int(FracsPerSec, I64)
  Fraction = Fraction / 1000_I64

  CurrentClockTime = InitTime(                          &
                       Interval   = .false.,            &
                       Clock      = .true.,             &
                       Infinite   = 0,                  &
                       Year       = Values(1),          &
                       Month      = Values(2),          &
                       Day        = Values(3),          &
                       Hour       = Values(5),          &
                       Minute     = Values(6),          &
                       Second     = Values(7),          &
                       Fraction   = Int(Fraction, I32), &
                       ZoneHour   = 0,                  &
                       ZoneMinute = Values(4)           &
                     )

End Function CurrentClockTime

!-------------------------------------------------------------------------------------------------------------

Function InitTime(                                           &
           Interval, Clock,                                  &
           Infinite,                                         &
           Year, Month, Day, Hour, Minute, Second, Fraction, &
           ZoneHour, ZoneMinute                              &
         )                                                   &
Result(Time)
! Initialises a time. The inputs to the routine need not be justified but the result returned is always a
! justified time.

  Implicit None
  ! Argument list:
  Logical, Intent(In)           :: Interval   ! Indicates whether the time is a time interval or the time of
                                              ! an event.
  Logical, Intent(In)           :: Clock      ! Indicates whether the time is a clock time or a model time.
  Integer, Intent(In)           :: Infinite   ! 1, 0 and -1 indicate +infinity, a finite time and -infinity.
  Integer, Intent(In), Optional :: Year       ! Year.
  Integer, Intent(In), Optional :: Month      ! Month.
  Integer, Intent(In), Optional :: Day        ! Day.
  Integer, Intent(In), Optional :: Hour       ! Hour.
  Integer, Intent(In), Optional :: Minute     ! Minute.
  Integer, Intent(In), Optional :: Second     ! Second.
  Integer, Intent(In), Optional :: Fraction   ! Fraction of a second in units defined by FracsPerSec.
  Integer, Intent(In), Optional :: ZoneHour   ! Time zone hour.
  Integer, Intent(In), Optional :: ZoneMinute ! Time zone minute.
  ! Function result:
  Type(Time_) :: Time ! Initialised time.

  If (Infinite < -1 .or. Infinite > 1) Call Message('UNEXPECTED FATAL ERROR in InitTime', 4)

  Time%Interval = Interval
  Time%Clock    = Clock

  Time%Infinite = Infinite

  If (Infinite == 0) Then
    If (UsesYearAndMonth(Time)) Then
      If (                             &
        .not. Present(Year)       .or. &
        .not. Present(Month)      .or. &
        .not. Present(Day)        .or. &
        .not. Present(Hour)       .or. &
        .not. Present(Minute)     .or. &
        .not. Present(Second)     .or. &
        .not. Present(Fraction)   .or. &
        .not. Present(ZoneHour)   .or. &
        .not. Present(ZoneMinute)      &
      ) Then
        Call Message('UNEXPECTED FATAL ERROR in InitTime', 4)
      End If
    Else
      If (                           &
        .not. Present(Day)      .or. &
        .not. Present(Hour)     .or. &
        .not. Present(Minute)   .or. &
        .not. Present(Second)   .or. &
        .not. Present(Fraction)      &
      ) Then
        Call Message('UNEXPECTED FATAL ERROR in InitTime', 4)
      End If
    End If
  End If

  If (Infinite == 0) Then
    Time%Day      = Day
    Time%Hour     = Hour
    Time%Minute   = Minute
    Time%Second   = Second
    Time%Fraction = Fraction
    If (UsesYearAndMonth(Time)) Then
      Time%Year       = Year
      Time%Month      = Month
      Time%ZoneHour   = ZoneHour
      Time%ZoneMinute = ZoneMinute
    End If
  End If

  Call JustifyTime(Time)

End Function InitTime

!-------------------------------------------------------------------------------------------------------------

Subroutine JustifyTime(Time)
! Justifies a time so that
!     (i)   the meaning is unchanged,
!     (ii)  the meaning of the time zone is unchanged,
!     (iii) Month, Day, Hour, Minute, Second and Fraction are all as small as possible for non-interval
!           calendar-frame times (with Month >= 1, Day >= 1, Hour >= 0, Minute >= 0, Second >= 0, Fraction >=
!           0) and Hour, Minute, Second and Fraction are all as small as possible for time intervals and
!           relative-frame times (with Hour >= 0, Minute >= 0, Second >= 0, Fraction >= 0), and
!     (iv)  ZoneHour and ZoneMinute have the same sign and ZoneMinute lies in [-59, 59] for non-interval
!           calendar-frame times.
! Has no effect on infinite times. Unlike most times, the input to this routine does not of course have to be
! justified.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(InOut) :: Time ! Time to be justified.
  ! Locals:
  Integer :: YearShift   !} Number of years, days, hours, minutes and seconds the time needs to be adjusted
  Integer :: DayShift    !} by.
  Integer :: HourShift   !}
  Integer :: MinuteShift !}
  Integer :: SecondShift !}
  Integer :: MonthDays   ! Number of days in month.
  Integer :: ZoneMinute  ! Time zone in minutes.

  If (Time%Infinite == 0) Then

    SecondShift     = Time%Fraction / FracsPerSec
    Time%Second     = Time%Second + SecondShift
    Time%Fraction   = Time%Fraction - FracsPerSec * SecondShift
    If (Time%Fraction < 0) Then
      Time%Fraction = Time%Fraction + FracsPerSec
      Time%Second   = Time%Second - 1
    End If

    MinuteShift     = Time%Second / 60
    Time%Minute     = Time%Minute + MinuteShift
    Time%Second     = Time%Second - 60 * MinuteShift
    If (Time%Second < 0) Then
      Time%Second   = Time%Second + 60
      Time%Minute   = Time%Minute - 1
    End If

    HourShift       = Time%Minute / 60
    Time%Hour       = Time%Hour + HourShift
    Time%Minute     = Time%Minute - 60 * HourShift
    If (Time%Minute < 0) Then
      Time%Minute   = Time%Minute + 60
      Time%Hour     = Time%Hour - 1
    End If

    DayShift        = Time%Hour / 24
    Time%Day        = Time%Day + DayShift
    Time%Hour       = Time%Hour - 24 * DayShift
    If (Time%Hour < 0) Then
      Time%Hour     = Time%Hour + 24
      Time%Day      = Time%Day - 1
    End If

    If (UsesYearAndMonth(Time)) Then

      YearShift    = (Time%Month - 1) / 12
      Time%Year    = Time%Year + YearShift
      Time%Month   = Time%Month - 12 * YearShift
      If (Time%Month - 1 < 0) Then
        Time%Month = Time%Month + 12
        Time%Year  = Time%Year - 1
      End If

      MonthDays = DaysInMonth(Time%Year, Time%Month, IsClockTime(Time))
      Do While (Time%Day > MonthDays)
        Time%Day   = Time%Day - MonthDays
        Time%Month = Time%Month + 1
        If (Time%Month > 12) Then
          Time%Month = Time%Month - 12
          Time%Year  = Time%Year + 1
        End If
        MonthDays = DaysInMonth(Time%Year, Time%Month, IsClockTime(Time))
      End Do

      Do While (Time%Day < 1)
        MonthDays  = DaysInMonth(Time%Year, Time%Month - 1, IsClockTime(Time))
        Time%Day   = Time%Day + MonthDays
        Time%Month = Time%Month - 1
        If (Time%Month < 1) Then
          Time%Month = Time%Month + 12
          Time%Year  = Time%Year - 1
        End If
      End Do

      ZoneMinute      = 60 * Time%ZoneHour + Time%ZoneMinute
      Time%ZoneHour   = ZoneMinute / 60
      Time%ZoneMinute = ZoneMinute - 60 * Time%ZoneHour

    End If

  End If

End Subroutine JustifyTime

!-------------------------------------------------------------------------------------------------------------

Subroutine ConvertTime(Year, Month, Time)
! Converts a time so that the meaning is the same, but the time's values of Year and Month are equal to
! specified values. Has no effect except on finite non-interval calendar-frame times. Unlike most times, the
! result is of course not necessarily justified as regards the value of Day.

  Implicit None
  ! Argument list:
  Integer,     Intent(In)    :: Year  !} Values of Year and Month required. Month must lie in [1,12].
  Integer,     Intent(In)    :: Month !}
  Type(Time_), Intent(InOut) :: Time  ! Time.
  ! Locals:
  Integer :: DayShift ! Number of days the time needs to be adjusted by.
  Integer :: i        ! Loop index.

  If (Month < 1 .or. Month > 12) Call Message('UNEXPECTED FATAL ERROR in ConvertTime', 4)

  If (UsesYearAndMonth(Time)) Then

    DayShift = 0
    Do i = 2, Time%Month
      DayShift = DayShift + DaysInMonth(Time%Year, i - 1, IsClockTime(Time))
    End Do
    Do i = 2, Month
      DayShift = DayShift - DaysInMonth(Year, i - 1, IsClockTime(Time))
    End Do
    If (Year < Time%Year) Then
      Do i = Year, Time%Year - 1
        DayShift = DayShift + DaysInYear(i, IsClockTime(Time))
      End Do
    End If
    If (Time%Year < Year) Then
      Do i = Time%Year, Year - 1
        DayShift = DayShift - DaysInYear(i, IsClockTime(Time))
      End Do
    End If

    Time%Year  = Year
    Time%Month = Month
    Time%Day   = Time%Day + DayShift

  End If

End Subroutine ConvertTime

!-------------------------------------------------------------------------------------------------------------

Function DaysInMonth(Year, Month, Clock)
! Evaluates number of days in a month.

  Implicit None
  ! Argument list:
  Integer, Intent(In)           :: Year  ! Year.
  Integer, Intent(In)           :: Month ! Month. Must lie in [0,12] (it is convenient to allow 0 for use in
                                         ! JustifyTime).
  Logical, Intent(In), Optional :: Clock ! If present and true, calculates result for clock times.
  ! Function result:
  Integer :: DaysInMonth ! Days in the month.
  ! Local parameters:
  Integer, Parameter :: MonthDays(12) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  ! MonthDays :: Number of days  in each month for non leap years.
  ! Locals:
  Logical :: Gregorian ! Indicates the Gregorian calendar is to be used.

  If (Month < 0 .or. Month > 12) Call Message('UNEXPECTED FATAL ERROR in DaysInMonth', 4)

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  Gregorian = .false.
  If (Present(Clock)) Then
    If (Clock) Then
      Gregorian = .true.
    Else If (IsGregorian()) Then
      Gregorian = .true.
    End If
  Else If (IsGregorian()) Then
    Gregorian = .true.
  End If

  If (Gregorian) Then

    If (Month == 0) Then
      DaysInMonth = MonthDays(12)
    Else If (Month == 2) Then
      If (IsGregorianLeapYear(Year)) Then
        DaysInMonth = 29
      Else
        DaysInMonth = MonthDays(2)
      End If
    Else
      DaysInMonth = MonthDays(Month)
    End If

  Else If (Is360DayYear()) Then

    DaysInMonth = 30

  Else

    Call Message('UNEXPECTED FATAL ERROR in DaysInMonth', 4)

  End If

End Function DaysInMonth

!-------------------------------------------------------------------------------------------------------------

Function DaysInYear(Year, Clock)
! Evaluates number of days in a year.

  Implicit None
  ! Argument list:
  Integer, Intent(In)           :: Year  ! Year.
  Logical, Intent(In), Optional :: Clock ! If present and true, calculates result for clock times.
  ! Function result:
  Integer :: DaysInYear ! Days in the year.
  ! Locals:
  Logical :: Gregorian ! Indicates the Gregorian calendar is to be used.

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  Gregorian = .false.
  If (Present(Clock)) Then
    If (Clock) Then
      Gregorian = .true.
    Else If (IsGregorian()) Then
      Gregorian = .true.
    End If
  Else If (IsGregorian()) Then
    Gregorian = .true.
  End If

  If (Gregorian) Then

    If (IsGregorianLeapYear(Year)) Then
      DaysInYear = 366
    Else
      DaysInYear = 365
    End If

  Else If (Is360DayYear()) Then

    DaysInYear = 360

  Else

    Call Message('UNEXPECTED FATAL ERROR in DaysInYear', 4)

  End If

End Function DaysInYear

!-------------------------------------------------------------------------------------------------------------

Function IsGregorianLeapYear(Year)
! Indicates whether it is a leap year in the Gregorian calendar.

  Implicit None
  ! Argument list:
  Integer, Intent(In) :: Year ! Year.
  ! Function result:
  Logical :: IsGregorianLeapYear ! Indicates whether it is a leap year in the Gregorian calendar.

  IsGregorianLeapYear =                           Year ==   4 * (Year /   4)
  IsGregorianLeapYear = IsGregorianLeapYear .and. Year /= 100 * (Year / 100)
  IsGregorianLeapYear = IsGregorianLeapYear .or.  Year == 400 * (Year / 400)

End Function IsGregorianLeapYear

!-------------------------------------------------------------------------------------------------------------

Function TimeUsesYearAndMonth(Time)
! Finds out whether a time uses year, month and time zone information.

! Note this routine can be used on partially constructed times provided the Interval, Clock and Infinite parts
! of the time are correct.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time ! The time.
  ! Function result:
  Logical :: TimeUsesYearAndMonth ! Indicates Year, Month, ZoneHour and ZoneMinute are used. The result is
                                  ! true only for finite non-interval calendar-frame times.

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (Time%Infinite /= 0) Then
    TimeUsesYearAndMonth = .false.
  Else If (IsClockTime(Time)) Then
    If (IsTimeInterval(Time)) Then
      TimeUsesYearAndMonth = .false.
    Else
      TimeUsesYearAndMonth = .true.
    End If
  Else If (IsCalendar()) Then
    If (IsTimeInterval(Time)) Then
      TimeUsesYearAndMonth = .false.
    Else
      TimeUsesYearAndMonth = .true.
    End If
  Else
    TimeUsesYearAndMonth = .false.
  End If

End Function TimeUsesYearAndMonth

!-------------------------------------------------------------------------------------------------------------

Function ShortTimeEq(ShortTime1, ShortTime2)
! Tests for ShortTime1 = ShortTime2.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime1 !} Times to be compared. They must both be of the same type.
  Type(ShortTime_), Intent(In) :: ShortTime2 !}
  ! Function result:
  Logical :: ShortTimeEq ! Indicates whether ShortTime1 = ShortTime2.

  If (                                                                    &
    ( IsTimeInterval(ShortTime1) .neqv. IsTimeInterval(ShortTime2) ) .or. &
    ( IsClockTime   (ShortTime1) .neqv. IsClockTime   (ShortTime2) )      &
  ) Then
    Call Message('UNEXPECTED FATAL ERROR in ShortTimeEq', 4)
  End If

  ShortTimeEq = ShortTime1%Infinite == ShortTime2%Infinite   .and. &
                (                                                  &
                  ShortTime1%Infinite /= 0                   .or.  &
                  ShortTime1%Fraction == ShortTime2%Fraction       &
                )

End Function ShortTimeEq

!-------------------------------------------------------------------------------------------------------------

Function ShortTimeNE(ShortTime1, ShortTime2)
! Tests for ShortTime1 /= ShortTime2.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime1 !} Times to be compared. They must both be of the same type.
  Type(ShortTime_), Intent(In) :: ShortTime2 !}
  ! Function result:
  Logical :: ShortTimeNE ! Indicates whether ShortTime1 /= ShortTime2.

  ShortTimeNE = .not. (ShortTime2 == ShortTime1)

End Function ShortTimeNE

!-------------------------------------------------------------------------------------------------------------

Function ShortTimeGE(ShortTime1, ShortTime2)
! Tests for ShortTime1 >= ShortTime2.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime1 !} Times to be compared. They must both be of the same type.
  Type(ShortTime_), Intent(In) :: ShortTime2 !}
  ! Function result:
  Logical :: ShortTimeGE ! Indicates whether ShortTime1 >= ShortTime2.

  If (                                                                    &
    ( IsTimeInterval(ShortTime1) .neqv. IsTimeInterval(ShortTime2) ) .or. &
    ( IsClockTime   (ShortTime1) .neqv. IsClockTime   (ShortTime2) )      &
  ) Then
    Call Message('UNEXPECTED FATAL ERROR in ShortTimeGE', 4)
  End If

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (IsClockTime(ShortTime1)) Then
    ShortTimeGE = ShortTime1%Infinite > ShortTime2%Infinite      .or.  &
                  (                                                    &
                    ShortTime1%Infinite == ShortTime2%Infinite   .and. &
                    (                                                  &
                      ShortTime1%Infinite /= 0                   .or.  &
                      ShortTime1%Fraction >= ShortTime2%Fraction       &
                    )                                                  &
                  )
  Else If (IsBackwards()) Then
    ShortTimeGE = ShortTime1%Infinite < ShortTime2%Infinite      .or.  &
                  (                                                    &
                    ShortTime1%Infinite == ShortTime2%Infinite   .and. &
                    (                                                  &
                      ShortTime1%Infinite /= 0                   .or.  &
                      ShortTime1%Fraction <= ShortTime2%Fraction       &
                    )                                                  &
                  )
  Else
    ShortTimeGE = ShortTime1%Infinite > ShortTime2%Infinite      .or.  &
                  (                                                    &
                    ShortTime1%Infinite == ShortTime2%Infinite   .and. &
                    (                                                  &
                      ShortTime1%Infinite /= 0                   .or.  &
                      ShortTime1%Fraction >= ShortTime2%Fraction       &
                    )                                                  &
                  )
  End If

End Function ShortTimeGE

!-------------------------------------------------------------------------------------------------------------

Function ShortTimeGT(ShortTime1, ShortTime2)
! Tests for ShortTime1 > ShortTime2.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime1 !} Times to be compared. They must both be of the same type.
  Type(ShortTime_), Intent(In) :: ShortTime2 !}
  ! Function result:
  Logical :: ShortTimeGT ! Indicates whether ShortTime1 > ShortTime2.

  ShortTimeGT = .not. ShortTime2 >= ShortTime1

End Function ShortTimeGT

!-------------------------------------------------------------------------------------------------------------

Function ShortTimeLE(ShortTime1, ShortTime2)
! Tests for ShortTime1 <= ShortTime2.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime1 !} Times to be compared. They must both be of the same type.
  Type(ShortTime_), Intent(In) :: ShortTime2 !}
  ! Function result:
  Logical :: ShortTimeLE ! Indicates whether ShortTime1 <= ShortTime2.

  ShortTimeLE = ShortTime2 >= ShortTime1

End Function ShortTimeLE

!-------------------------------------------------------------------------------------------------------------

Function ShortTimeLT(ShortTime1, ShortTime2)
! Tests for ShortTime1 < ShortTime2.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime1 !} Times to be compared. They must both be of the same type.
  Type(ShortTime_), Intent(In) :: ShortTime2 !}
  ! Function result:
  Logical :: ShortTimeLT ! Indicates whether ShortTime1 < ShortTime2.

  ShortTimeLT = ShortTime2 > ShortTime1

End Function ShortTimeLT

!-------------------------------------------------------------------------------------------------------------

Function AddShortTime(ShortTime, dT) Result(Sum)
! Adds two short times.

! +infinity + -infinity and -infinity + +infinity are not allowed.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime !} Times to be added. dT must be a time interval. Both must be
  Type(ShortTime_), Intent(In) :: dT        !} either model or clock times.
  ! Function result:
  Type(ShortTime_) :: Sum ! ShortTime + dT.

  If (                                                               &
    ( .not. IsTimeInterval(dT)                                ) .or. &
    ( IsClockTime(ShortTime)           .neqv. IsClockTime(dT) ) .or. &
    ( ShortTime%Infinite * dT%Infinite   ==   -1              )      &
  ) Then
    Call Message('UNEXPECTED FATAL ERROR in AddShortTime', 4)
  End If

  Sum%Interval = ShortTime%Interval

  Sum%Clock = ShortTime%Clock

  If (ShortTime%Infinite /= 0) Then
    Sum%Infinite = ShortTime%Infinite
  Else If (dT%Infinite /= 0) Then
    Sum%Infinite = dT%Infinite
  Else
    Sum%Infinite = 0
    Sum%Fraction = ShortTime%Fraction + dT%Fraction
  End If

End Function AddShortTime

!-------------------------------------------------------------------------------------------------------------

Function SubtractShortTime(ShortTime1, ShortTime2) Result(Difference)
! Subtracts one short time from another.

! +infinity - +infinity or -infinity - -infinity are not allowed.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime1 !} Times to be subtracted. They must both be or both not be time
  Type(ShortTime_), Intent(In) :: ShortTime2 !} intervals or ShortTime2 must be a time interval. Both must be
                                             !} either model or clock times.
  ! Function result:
  Type(ShortTime_) :: Difference ! ShortTime1 - ShortTime2.

  If (                                                                                         &
    ( IsTimeInterval(ShortTime1)                .and.  .not. IsTimeInterval(ShortTime2) ) .or. &
    ( IsClockTime(ShortTime1)                   .neqv. IsClockTime(ShortTime2)          ) .or. &
    ( ShortTime1%Infinite * ShortTime2%Infinite   == 1                                  )      &
  ) Then
    Call Message('UNEXPECTED FATAL ERROR in SubtractShortTime', 4)
  End If

  If (IsTimeInterval(ShortTime1) .eqv. IsTimeInterval(ShortTime2)) Then
    Difference%Interval = .true.
  Else If (IsTimeInterval(ShortTime2)) Then
    Difference%Interval = ShortTime1%Interval
  End If

  Difference%Clock = ShortTime1%Clock

  If (ShortTime1%Infinite /= 0) Then
    Difference%Infinite = ShortTime1%Infinite
  Else If (ShortTime2%Infinite /= 0) Then
    Difference%Infinite = - ShortTime2%Infinite
  Else
    Difference%Infinite = 0
    Difference%Fraction = ShortTime1%Fraction - ShortTime2%Fraction
  End If

End Function SubtractShortTime

!-------------------------------------------------------------------------------------------------------------

Function MultiplyShortTime(dT, n) Result(Product)
! Multiplies a short time by an integer.

! +/-infinity * 0 is not allowed.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: dT ! Time to be multiplied by n. Must be a time interval.
  Integer,          Intent(In) :: n  ! Quantity by which dT is to be multiplied.
  ! Function result:
  Type(ShortTime_) :: Product ! dT * n.

  If (                                             &
    ( .not. IsTimeInterval(dT)              ) .or. &
    ( dT%Infinite /= 0         .and. n == 0 )      &
  ) Then
    Call Message('UNEXPECTED FATAL ERROR in MultiplyShortTime', 4)
  End If

  Product%Interval = dT%Interval

  Product%Clock = dT%Clock

  Product%Infinite = dT%Infinite * Sign(1, n)
  If (Product%Infinite == 0) Then
    Product%Fraction = dT%Fraction * Int(n, I64)
  End If

End Function MultiplyShortTime

!-------------------------------------------------------------------------------------------------------------

Function DivideShortTime(dT, n) Result(Quotient)
! Divides a short time by an integer (rounded towards zero).

! Infinite times or zero denominators are not allowed.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: dT ! Time to be divided by n. Must be a time interval.
  Integer,          Intent(In) :: n  ! Quantity by which dT is to be divided.
  ! Function result:
  Type(ShortTime_) :: Quotient ! dT / n.

  If (                                            &
    ( .not. IsTimeInterval(dT)             ) .or. &
    ( dT%Infinite /= 0         .or. n == 0 )      &
  ) Then
    Call Message('UNEXPECTED FATAL ERROR in DivideShortTime', 4)
  End If

  Quotient%Interval = dT%Interval

  Quotient%Clock = dT%Clock

  Quotient%Infinite = 0
  Quotient%Fraction = dT%Fraction / Int(n, I64)

End Function DivideShortTime

!-------------------------------------------------------------------------------------------------------------

Function RatioOfShortTimes(dT1, dT2) Result(Quotient)
! Finds the ratio of two short times (rounded towards zero).

! Infinite times or zero denominators are not allowed.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: dT1 !} Times for which the ratio is required. Must be time intervals. Both
  Type(ShortTime_), Intent(In) :: dT2 !} must be either model or clock times.
  ! Function result:
  Integer :: Quotient ! dT1 / dT2.

  If (                                                                  &
    ( .not. IsTimeInterval(dT1)  .or.  .not. IsTimeInterval(dT2) ) .or. &
    ( IsClockTime(dT1)          .neqv. IsClockTime(dT2)          ) .or. &
    ( dT1%Infinite /= 0                                          ) .or. &
    ( dT2%Infinite /= 0                                          ) .or. &
    ( dT2%Infinite == 0         .and.  dT2%Fraction == 0         )      &
  ) Then
    Call Message('UNEXPECTED FATAL ERROR in RatioOfShortTimes', 4)
  End If

  Quotient = Int(dT1%Fraction / dT2%Fraction, I32)

End Function RatioOfShortTimes

!-------------------------------------------------------------------------------------------------------------

Function RoundShortTime(Time, T0, dT, Up)
! Rounds a short time to a time of the form T0 + n*dT.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In)           :: Time ! Time. } Must be finite. Time and T0 must both be or both
  Type(ShortTime_), Intent(In)           :: T0   ! T0.   } not be time intervals. dT must be non-zero and a
  Type(ShortTime_), Intent(In)           :: dT   ! dT.   } time interval. All must be either model or clock
                                                 !       } times.
  Logical,          Intent(In), Optional :: Up   ! Indicates if rounding is up or down. If absent, rounds
                                                 ! towards T0.
  ! Function result:
  Type(ShortTime_) :: RoundShortTime ! Rounded short time.

  RoundShortTime = T0 + dT * ((Time - T0) / dT)

  If (Present(Up)) Then
    If (Up) Then
      If (dT > ZeroShortTime(IsClockTime(Time))) Then
        If (RoundShortTime < Time) RoundShortTime = RoundShortTime + dT
      Else
        If (RoundShortTime < Time) RoundShortTime = RoundShortTime - dT
      End If
    Else
      If (dT > ZeroShortTime(IsClockTime(Time))) Then
        If (RoundShortTime > Time) RoundShortTime = RoundShortTime - dT
      Else
        If (RoundShortTime > Time) RoundShortTime = RoundShortTime + dT
      End If
    End If
  End If

End Function RoundShortTime

!-------------------------------------------------------------------------------------------------------------

Function RoundShortTimeToYear(Time, Up, Strictly)
! Rounds a short time to the start or end of the year.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: Time     ! Time. Must be finite and not a time interval.
  Logical,          Intent(In) :: Up       ! Indicates if rounding is up or down.
  Logical,          Intent(In) :: Strictly ! Ensures the rounded time doesn't equal time.
  ! Function result:
  Type(ShortTime_) :: RoundShortTimeToYear ! Rounded time.

  RoundShortTimeToYear = Time2ShortTime(RoundTimeToYear(ShortTime2Time(Time), Up, Strictly))

End Function RoundShortTimeToYear

!-------------------------------------------------------------------------------------------------------------

Function RoundShortTimeToMonth(Time, Up, Strictly)
! Rounds a short time to the start or end of the month.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: Time     ! Time. Must be finite and not a time interval.
  Logical,          Intent(In) :: Up       ! Indicates if rounding is up or down.
  Logical,          Intent(In) :: Strictly ! Ensures the rounded time doesn't equal time.
  ! Function result:
  Type(ShortTime_) :: RoundShortTimeToMonth ! Rounded time.

  RoundShortTimeToMonth = Time2ShortTime(RoundTimeToMonth(ShortTime2Time(Time), Up, Strictly))

End Function RoundShortTimeToMonth

!-------------------------------------------------------------------------------------------------------------

Function ShortTimeMin(ShortTime1, ShortTime2)
! Returns the minimum of two short times.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime1 !} The two short times. They must both be of the same type.
  Type(ShortTime_), Intent(In) :: ShortTime2 !}
  ! Function result:
  Type(ShortTime_) :: ShortTimeMin ! Minimum of ShortTime1 and ShortTime2.

  If (ShortTime2 >= ShortTime1) Then
    ShortTimeMin = ShortTime1
  Else
    ShortTimeMin = ShortTime2
  End If

End Function ShortTimeMin

!-------------------------------------------------------------------------------------------------------------

Function ShortTimeMax(ShortTime1, ShortTime2)
! Returns the maximum of two short times.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime1 !} The two short times. They must both be of the same type.
  Type(ShortTime_), Intent(In) :: ShortTime2 !}
  ! Function result:
  Type(ShortTime_) :: ShortTimeMax ! Maximum of ShortTime1 and ShortTime2.

  If (ShortTime2 >= ShortTime1) Then
    ShortTimeMax = ShortTime2
  Else
    ShortTimeMax = ShortTime1
  End If

End Function ShortTimeMax

!-------------------------------------------------------------------------------------------------------------

Function ShortTimeIsInfFuture(ShortTime)
! Finds out whether a short time is in the infinite future (in the run direction for model times).

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime ! The short time.
  ! Function result:
  Logical :: ShortTimeIsInfFuture ! Indicates the short time is in the infinite future.

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (IsClockTime(ShortTime)) Then
    ShortTimeIsInfFuture = ShortTime%Infinite == 1
  Else If (IsBackwards()) Then
    ShortTimeIsInfFuture = ShortTime%Infinite == -1
  Else
    ShortTimeIsInfFuture = ShortTime%Infinite == 1
  End If

End Function ShortTimeIsInfFuture

!-------------------------------------------------------------------------------------------------------------

Function ShortTimeIsInfPast(ShortTime)
! Finds out whether a short time is in the infinite past (in the run direction for model times).

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime ! The short time.
  ! Function result:
  Logical :: ShortTimeIsInfPast ! Indicates the short time is in the infinite past.

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (IsClockTime(ShortTime)) Then
    ShortTimeIsInfPast = ShortTime%Infinite == -1
  Else If (IsBackwards()) Then
    ShortTimeIsInfPast = ShortTime%Infinite == 1
  Else
    ShortTimeIsInfPast = ShortTime%Infinite == -1
  End If

End Function ShortTimeIsInfPast

!-------------------------------------------------------------------------------------------------------------

Function ShortTimeIsTimeInterval(ShortTime)
! Finds out whether a time is a time interval.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime ! Time.
  ! Function result:
  Logical :: ShortTimeIsTimeInterval ! Indicates the time is a time interval.

  ShortTimeIsTimeInterval = ShortTime%Interval

End Function ShortTimeIsTimeInterval

!-------------------------------------------------------------------------------------------------------------

Function ShortTimeIsClockTime(ShortTime)
! Finds out whether a time is a clock time.

! Note this routine can be used on partially constructed times provided the Clock part of the time is correct.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime ! Time.
  ! Function result:
  Logical :: ShortTimeIsClockTime ! Indicates the time is a clock time.

  ShortTimeIsClockTime = ShortTime%Clock

End Function ShortTimeIsClockTime

!-------------------------------------------------------------------------------------------------------------

Function InfFutureShortTime(Interval, Clock) Result(ShortTime)
! Returns a short time in the infinite future (in run direction for model times).

  Implicit None
  ! Argument list:
  Logical, Intent(In), Optional :: Interval ! If present and true, the time is a time interval.
  Logical, Intent(In), Optional :: Clock    ! If present and true, the time is a clock time.
  ! Function result:
  Type(ShortTime_) :: ShortTime ! The short time in the infinite future.

  ShortTime%Interval = .false.
  If (Present(Interval)) Then
    If (Interval) ShortTime%Interval = .true.
  End If

  ShortTime%Clock = .false.
  If (Present(Clock)) Then
    If (Clock) ShortTime%Clock = .true.
  End If

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (IsClockTime(ShortTime)) Then
    ShortTime%Infinite = 1
  Else If (IsBackwards()) Then
    ShortTime%Infinite = -1
  Else
    ShortTime%Infinite = 1
  End If

End Function InfFutureShortTime

!-------------------------------------------------------------------------------------------------------------

Function InfPastShortTime(Interval, Clock) Result(ShortTime)
! Returns a short time in the infinite past (in run direction for model times).

  Implicit None
  ! Argument list:
  Logical, Intent(In), Optional :: Interval ! If present and true, the time is a time interval.
  Logical, Intent(In), Optional :: Clock    ! If present and true, the time is a clock time.
  ! Function result:
  Type(ShortTime_) :: ShortTime ! The short time in the infinite past.

  ShortTime%Interval = .false.
  If (Present(Interval)) Then
    If (Interval) ShortTime%Interval = .true.
  End If

  ShortTime%Clock = .false.
  If (Present(Clock)) Then
    If (Clock) ShortTime%Clock = .true.
  End If

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (IsClockTime(ShortTime)) Then
    ShortTime%Infinite = -1
  Else If (IsBackwards()) Then
    ShortTime%Infinite = 1
  Else
    ShortTime%Infinite = -1
  End If

End Function InfPastShortTime

!-------------------------------------------------------------------------------------------------------------

Function ZeroShortTime(Clock) Result(ShortTime)
! Returns a zero short time interval.

  Implicit None
  ! Argument list:
  Logical, Intent(In), Optional :: Clock ! If present and true, the time is a clock time.
  ! Function result:
  Type(ShortTime_) :: ShortTime ! The zero short time interval.

  ShortTime%Interval = .true.

  ShortTime%Clock = .false.
  If (Present(Clock)) Then
    If (Clock) ShortTime%Clock = .true.
  End If

  ShortTime%Infinite = 0
  ShortTime%Fraction = 0

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (.not. IsClockTime(ShortTime) .and. .not.Initialised) Then
    Call Message('UNEXPECTED FATAL ERROR in ZeroShortTime', 4)
  End If

End Function ZeroShortTime

!-------------------------------------------------------------------------------------------------------------

Function ReferenceShortTime(Clock)
! Returns a reference short time (for when one wants a short time variable set to a fixed time, but the time
! itself doesn't matter).

  Implicit None
  ! Argument list:
  Logical, Intent(In), Optional :: Clock ! If present and true, the time is a clock time.
  ! Function result:
  Type(ShortTime_) :: ReferenceShortTime ! The reference short time.

  ReferenceShortTime = Time2ShortTime(ReferenceTime(Clock))

End Function ReferenceShortTime

!-------------------------------------------------------------------------------------------------------------

Function Time2ShortTime(Time) Result(ShortTime)
! Converts a time to a short time.

  Implicit None
  ! Argument list:
  Type(Time_), Intent(In) :: Time ! Time to be converted.
  ! Function result:
  Type(ShortTime_) :: ShortTime ! The converted time.
  ! Locals:
  Type(Time_) :: TempTime ! Temporary variable.

  ShortTime%Interval = Time%Interval
  ShortTime%Clock    = Time%Clock

  ShortTime%Infinite = Time%Infinite

  If (Time%Infinite == 0) Then

    If (UsesYearAndMonth(Time)) Then
      TempTime = Time
      Call ConvertTime(YearOrigin, MonthOrigin, TempTime)
      TempTime%Hour     = TempTime%Hour   - TempTime%ZoneHour
      TempTime%Minute   = TempTime%Minute - TempTime%ZoneMinute
      TempTime%Day      = TempTime%Day    - DayOrigin
    Else
      TempTime = Time
    End If

    ShortTime%Fraction = (TempTime%Day * 24 + TempTime%Hour) * 60 + TempTime%Minute
    ShortTime%Fraction = ShortTime%Fraction * 60_I64 + Int(TempTime%Second, I64)
    ShortTime%Fraction = ShortTime%Fraction * Int(FracsPerSec, I64) + Int(TempTime%Fraction, I64)

  End If

End Function Time2ShortTime

!-------------------------------------------------------------------------------------------------------------

Function ShortTime2Time(ShortTime) Result(Time)
! Converts a short time to a time.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime ! Short time to be converted.
  ! Function result:
  Type(Time_) :: Time ! The converted short time.
  ! Local parameters:
  Integer(I64), Parameter :: FracsPerMin = Int(FracsPerSec, I64) * 60_I64 ! Fractions per minute.
  ! Locals:
  Integer(I64) :: Minute   !} Time expressed as minutes and fractions.
  Integer(I64) :: Fraction !}

  ! Infinite times.
  If (ShortTime%Infinite /= 0) Then

    Time%Interval = ShortTime%Interval
    Time%Clock    = ShortTime%Clock

    Time%Infinite = ShortTime%Infinite

  ! Finite times.
  Else

    Minute   = ShortTime%Fraction/FracsPerMin
    Fraction = ShortTime%Fraction - Minute * FracsPerMin

    If (UsesYearAndMonth(ShortTime)) Then
      Time = InitTime(                                 &
               Interval   = IsTimeInterval(ShortTime), &
               Clock      = IsClockTime   (ShortTime), &
               Infinite   = 0,                         &
               Year       = YearOrigin,                &
               Month      = MonthOrigin,               &
               Day        = DayOrigin,                 &
               Hour       = 0,                         &
               Minute     = Int(Minute, I32),          &
               Second     = 0,                         &
               Fraction   = Int(Fraction, I32),        &
               ZoneHour   = 0,                         &
               ZoneMinute = 0                          &
             )
    Else
      Time = InitTime(                               &
               Interval = IsTimeInterval(ShortTime), &
               Clock    = IsClockTime   (ShortTime), &
               Infinite = 0,                         &
               Day      = 0,                         &
               Hour     = 0,                         &
               Minute   = Int(Minute, I32),          &
               Second   = 0,                         &
               Fraction = Int(Fraction, I32)         &
             )
    End If

    Call JustifyTime(Time)

  End If

End Function ShortTime2Time

!-------------------------------------------------------------------------------------------------------------

Function ShortTime2RealTime(ShortTime) Result(RealTime)
! Converts a short time interval to a duration in seconds.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime ! Short time to be converted. Must be a finite time interval.
  ! Function result:
  Real(Std) :: RealTime ! The converted short time.

  If (.not. IsTimeInterval(ShortTime) .or. ShortTime%Infinite /= 0) Then
    Call Message('UNEXPECTED FATAL ERROR in ShortTime2RealTime', 4)
  End If

  RealTime = Real(ShortTime%Fraction, P32) / Real(FracsPerSec, P32)

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (.not. IsClockTime(ShortTime)) Then
    If (IsBackwards()) Then
      RealTime = - RealTime
    End If
  End If

End Function ShortTime2RealTime

!-------------------------------------------------------------------------------------------------------------

Function RealTime2ShortTime(RealTime, Clock) Result(ShortTime)
! Converts a duration in seconds to a short time interval.

  Implicit None
  ! Argument list:
  Real(Std), Intent(In)           :: RealTime ! Time in seconds to be converted.
  Logical,   Intent(In), Optional :: Clock    ! If present and true, the real time is a clock time.
  ! Function result:
  Type(ShortTime_) :: ShortTime ! The converted time in seconds.

  ShortTime%Interval = .true.

  ShortTime%Clock = .false.
  If (Present(Clock)) Then
    If (Clock) Then
      ShortTime%Clock = .true.
    End If
  End If

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (.not. IsClockTime(ShortTime) .and. .not.Initialised) Then
    Call Message('UNEXPECTED FATAL ERROR in RealTime2ShortTime', 4)
  End If

  ShortTime%Infinite = 0
  ShortTime%Fraction = Int(Real(RealTime, P64) * Real(FracsPerSec, P64), I64)

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (.not. IsClockTime(ShortTime)) Then
    If (IsBackwards()) Then
      ShortTime%Fraction = - ShortTime%Fraction
    End If
  End If

End Function RealTime2ShortTime

!-------------------------------------------------------------------------------------------------------------

Function ShortTimeUsesYearAndMonth(ShortTime)
! Finds out whether a short time uses year, month and time zone information.

! Note this routine can be used on partially constructed times provided the Interval, Clock and Infinite parts
! of the time are correct.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: ShortTime ! The short time.
  ! Function result:
  Logical :: ShortTimeUsesYearAndMonth ! Indicates Year, Month, ZoneHour and ZoneMinute are used. The result
                                       ! is true only for finite non-interval calendar-frame times.

  ! Note this is designed so clock times can be processed before InitTimeModule is called.
  If (ShortTime%Infinite /= 0) Then
    ShortTimeUsesYearAndMonth = .false.
  Else If (IsClockTime(ShortTime)) Then
    If (IsTimeInterval(ShortTime)) Then
      ShortTimeUsesYearAndMonth = .false.
    Else
      ShortTimeUsesYearAndMonth = .true.
    End If
  Else If (IsCalendar()) Then
    If (IsTimeInterval(ShortTime)) Then
      ShortTimeUsesYearAndMonth = .false.
    Else
      ShortTimeUsesYearAndMonth = .true.
    End If
  Else
    ShortTimeUsesYearAndMonth = .false.
  End If

End Function ShortTimeUsesYearAndMonth

!-------------------------------------------------------------------------------------------------------------

Function Char2WildTime(CharTime, ErrorCode) Result (WildTime)
! Converts a character string to a wild-card time.

! Character formats for wild-card times are similar to those used for times. There are three possibilities:
! 1) Same format as used for non-interval times.
! 2) Same format as used for non-infinite non-interval times but with a number of elements replaced by * to
!    indicate a wild card. The elements which can be wild cards are year, month, day, hour, minute and second
!    for calendar time frames, and day, hour, minute, and second for relative time frames. The wild card
!    elements are always the first so many (possibly none or all) elements from the above list.
! 3) Same format as used for non-infinite non-interval calendar-frame times but with the year, month and day
!    replaced by a day of the week. The day of the week is case insensitive and can be abbreviated provided
!    this is unambiguous.
! Wild card values count as non-zero with regard to the statement in Char2Time (in defining the format of
! times) that 'values can only lie outside their appropriate range where we have finite times not involving
! year and month and where the value is the most significant non-zero element of the time'.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)            :: CharTime
  Integer,      Intent(Out), Optional :: ErrorCode
  ! CharTime  :: Character string to be converted.
  ! ErrorCode :: Non-zero values indicate that an error has occurred. If not present, the routine stops when
  !              there is an error. If an unexpected fatal error occurs the routine stops whether or not
  !              ErrorCode is present. The error codes are as follows:
  !                 $$ list error codes here
  ! Function result:
  Type(WildTime_) :: WildTime ! Wild-card time.
  ! Local parameters:
  Character(MaxCharLength) :: DaysOfWeek(7) = (/             & ! Names of days of the week.
                                                'Monday   ', & !
                                                'Tuesday  ', & !
                                                'Wednesday', & !
                                                'Thursday ', & !
                                                'Friday   ', & !
                                                'Saturday ', & !
                                                'Sunday   '  & !
                                              /)               !
  ! Locals:
  Integer :: i          ! Index of first non-blank character in CharTime.
  Integer :: j          ! Index of first character after day of week in CharTime.
  Integer :: k          ! Loop index.
  Integer :: ErrorCodeL ! Indicates an error has occurred in a subroutine or function called from here.

  If (.not.Initialised) Call Message('UNEXPECTED FATAL ERROR in Char2WildTime', 4)

  If (Present(ErrorCode)) ErrorCode = 0

  i = Verify(CharTime, ' ')

  If (i == 0) Go To 1

  ! Times with day of week specified.
  If (Scan(CharTime(i:i), 'MTWFSmtwfs') /= 0) Then

    If (.not. IsCalendar()) Go To 1

    If (Scan(CharTime, '*') /= 0) Go To 1

    j = Scan(CharTime(i + 1:), ' ') + i
    If (j == i) Go To 1
    WildTime%nNonWildCards = 4
    WildTime%Time = Char2Time(                                  &
                      CharTime       = '*/*/*' // CharTime(j:), &
                      AllowWildCards = .true.,                  &
                      ErrorCode      = ErrorCodeL               &
                    )
    If (ErrorCodeL > 0) Go To 1
    WildTime%DayOfWeek = 0
    Do k = 1, 7
      ! Identify day of week, noting that need to check at least two characters for
      ! Tues, Thurs, Sat and Sun.
      If (                                                                     &
        (CharTime(i:j - 1) .CIEq. DaysOfWeek(k)(1: Min(j - i, MaxCharLength))) &
        .and.                                                                  &
        (k == 1 .or. k == 3 .or. k == 5 .or. j - i >= 2)                       &
      ) Then
        WildTime%DayOfWeek = k
      End If
    End Do
    If (WildTime%DayOfWeek == 0) Go To 1

  ! Times without day of week specified.
  Else

    If (IsCalendar()) Then
      WildTime%nNonWildCards = 7
    Else
      WildTime%nNonWildCards = 5
    End If
    Do k = 1, Len_Trim(CharTime)
      If (CharTime(k:k) == '*') Then
        WildTime%nNonWildCards = WildTime%nNonWildCards - 1
      End If
    End Do
    WildTime%Time = Char2Time(                    &
                      CharTime       = CharTime,  &
                      AllowWildCards = .true.,    &
                      ErrorCode      = ErrorCodeL &
                    )
    If (ErrorCodeL > 0) Go To 1
    WildTime%DayOfWeek = 0

  End If

  Return

1 Continue
  If (Present(ErrorCode)) Then
    ErrorCode = 1
    Return
  Else
    Call Message(                                              & ! better messages? $$
           'FATAL ERROR: Time (possibly with wild cards) "' // & ! use errorcodeL
           Trim(CharTime)                                   // & ! to give useful messages
           '" has an incorrect format',                        &
           3                                                   &
         )
  End If

End Function Char2WildTime

!-------------------------------------------------------------------------------------------------------------

Function WildTimeEq(WildTime1, WildTime2)
! Tests for equality of wild-card times.

! Note that, unlike in TimeEq, equality here requires equality of the time zones.

  Implicit None
  ! Argument list:
  Type(WildTime_), Intent(In) :: WildTime1 !} The wild-card times.
  Type(WildTime_), Intent(In) :: WildTime2 !}
  ! Function result:
  Logical :: WildTimeEq ! Indicates if wild-card times are equal.

  WildTimeEq = WildTime1%nNonWildCards == WildTime2%nNonWildCards .and. &
               WildTime1%DayOfWeek     == WildTime2%DayOfWeek     .and. &
               WildTime1%Time          == WildTime2%Time
  If (IsCalendar() .and. WildTime1%Time%Infinite == 0) Then
    WildTimeEq = WildTimeEq                                             .and. &
                 WildTime1%Time%ZoneHour   == WildTime2%Time%ZoneHour   .and. &
                 WildTime1%Time%ZoneMinute == WildTime2%Time%ZoneMinute
  End If

End Function WildTimeEq

!-------------------------------------------------------------------------------------------------------------

Function SubstituteWildCards(Time, WildTime) Result(SubTime)
! Substitutes for the wild cards in a wild-card time using a given time.

  Implicit None
  ! Argument list:
  Type(Time_),     Intent(In) :: Time     ! The time.
  Type(WildTime_), Intent(In) :: WildTime ! The wild-card time. This must not specify a day of the week.
  ! Function result:
  Type(Time_) :: SubTime ! The time part of WildTime with the wild cards substituted.
  ! Locals:
  Integer     :: nNonWildCards ! Number of non-wild cards in WildTime.
  Type(Time_) :: TimeL         ! Local copy of Time.

  SubTime = WildTime%Time

  nNonWildCards = WildTime%nNonWildCards

  ! Convert time zone of Time.
  TimeL = ChangeTimeZone(Time, SubTime)

  ! Check no days of week.
  If (WildTime%DayOfWeek /= 0) Call Message('UNEXPECTED FATAL ERROR in SubstituteWildCards', 4)

  ! Calendar time frames.
  If (IsCalendar()) Then

    If (nNonWildCards <= 6) Then
      SubTime%Year = TimeL%Year
      If (nNonWildCards <= 5) Then
        SubTime%Month = TimeL%Month
        If (nNonWildCards <= 4) Then
          SubTime%Day = TimeL%Day
          If (nNonWildCards <= 3) Then
            SubTime%Hour = TimeL%Hour
            If (nNonWildCards <= 2) Then
              SubTime%Minute = TimeL%Minute
              If (nNonWildCards == 1) Then
                SubTime%Second = TimeL%Second
              End If
            End If
          End If
        End If
      End If
    End If

    ! Ensure impossible times in the Gregorian time frame, such as 29/2/2001 hh:mm or 31/4/2000 hh:mm are
    ! interpreted as midnight at the start of the following month. Note December has 31 days so there's no
    ! need to check whether the year needs incrementing too.
    If (SubTime%Day > DaysInMonth(SubTime%Year, SubTime%Month, IsClockTime(SubTime))) Then
      SubTime%Month    = SubTime%Month + 1
      SubTime%Day      = 1
      SubTime%Hour     = 0
      SubTime%Minute   = 0
      SubTime%Second   = 0
      SubTime%Fraction = 0
    End If

  ! Relative time frame.
  Else

    If (nNonWildCards <= 4) Then
      SubTime%Day = TimeL%Day
      If (nNonWildCards <= 3) Then
        SubTime%Hour = TimeL%Hour
        If (nNonWildCards <= 2) Then
          SubTime%Minute = TimeL%Minute
          If (nNonWildCards == 1) Then
            SubTime%Second = TimeL%Second
          End If
        End If
      End If
    End If

  End If

End Function SubstituteWildCards

!-------------------------------------------------------------------------------------------------------------

Function InitWildTimePair(FromTime, ToTime, ErrorCode) Result (WildTimePair)
! Initialises a wild-card time pair.

  Implicit None
  ! Argument list:
  Type(WildTime_), Intent(In)            :: FromTime  !} The pair of wild-card times to use in the wild-card
  Type(WildTime_), Intent(In)            :: ToTime    !} time pair.
  Integer,         Intent(Out), Optional :: ErrorCode ! Non-zero values indicate that an error has occurred.
                                                      ! If not present, the routine stops when there is an
                                                      ! error. The error codes are as follows:
                                                      !  $$ list codes
  ! Function result:
  Type(WildTimePair_) :: WildTimePair ! Initialised wild-card time pair.

  If (Present(ErrorCode)) ErrorCode = 0

  ! Check same time zone.
  If (IsCalendar()) Then
    If (FromTime%Time%Infinite == 0 .and. ToTime%Time%Infinite == 0) Then
      If (                                                       &
        FromTime%Time%ZoneHour   /= ToTime%Time%ZoneHour   .and. &
        FromTime%Time%ZoneMinute /= ToTime%Time%ZoneMinute       &
      ) Go To 1
    End If
  End If

  ! Check same number of wild cards.
  If (FromTime%nNonWildCards /= ToTime%nNonWildCards) Go To 2

  ! Check same as regards whether they specify a day of the week.
  If (                                                         &
    (FromTime%DayOfWeek == 0 .and. ToTime%DayOfWeek /= 0) .or. &
    (FromTime%DayOfWeek /= 0 .and. ToTime%DayOfWeek == 0)      &
  ) Go To 3

  WildTimePair%FromTime = FromTime
  WildTimePair%ToTime   = ToTime

  Return

1 Continue
  If (Present(ErrorCode)) Then
    ErrorCode = 1
    Return
  Else
    Call Message(                                                                &
           'FATAL ERROR: Two wild-card times to be used together to specify ' // &
           'a collection of time intervals have different time zones',           &
           3                                                                     &
         )
  End If

2 Continue
  If (Present(ErrorCode)) Then
    ErrorCode = 2
    Return
  Else
    Call Message(                                                                &
           'FATAL ERROR: Two wild-card times to be used together to specify ' // &
           'a collection of time intervals have a different number of wild '  // &
           'cards',                                                              &
           3                                                                     &
         )
  End If

3 Continue
  If (Present(ErrorCode)) Then
    ErrorCode = 3
    Return
  Else
    Call Message(                                                                &
           'FATAL ERROR: Two wild-card times to be used together to specify ' // &
           'a collection of time intervals are not the same with regard to '  // &
           'whether they specify a day of the week',                             &
           3                                                                     &
         )
  End If

End Function InitWildTimePair

!-------------------------------------------------------------------------------------------------------------

Function WildTimePairEq(WildTimePair1, WildTimePair2)
! Tests for equality of wild-card time pairs.

  Implicit None
  ! Argument list:
  Type(WildTimePair_), Intent(In) :: WildTimePair1 !} The wild-card time pairs.
  Type(WildTimePair_), Intent(In) :: WildTimePair2 !}
  ! Function result:
  Logical :: WildTimePairEq ! Indicates if wild-card time pairs are equal.

  WildTimePairEq = WildTimePair1%FromTime == WildTimePair2%FromTime .and. &
                   WildTimePair1%ToTime   == WildTimePair2%ToTime

End Function WildTimePairEq

!-------------------------------------------------------------------------------------------------------------

Subroutine InWildTimeInterval(Time, WildTimePair, In, NextTime)
! Checks whether a time lies in the time intervals defined by a pair of wild-card times.

  Implicit None
  ! Argument list:
  Type(Time_),         Intent(In)            :: Time
  Type(WildTimePair_), Intent(In)            :: WildTimePair
  Logical,             Intent(Out)           :: In
  Type(Time_),         Intent(Out), Optional :: NextTime
  ! Time         :: Time.
  ! WildTimePair :: The wild-card time pair.
  ! In           :: Indicates if the time lies in the time intervals.
  ! NextTime     :: Next time that a change might occur in whether the time lies in the time intervals defined
  !                 by the pair of wild-card times.
  ! Locals:
  Type(Time_) :: TimeL         !} Local copies of Time and the time-component of the two wild-card times in
  Type(Time_) :: ToTimeL       !} the wild-card time pair.
  Type(Time_) :: FromTimeL     !}
  Type(Time_) :: AlteredTime   ! Altered version of TimeL.
  Integer     :: nNonWildCards ! Number of non-wild cards in the wild-card times.
  Integer     :: Sign          ! 1 for forwards run and -1 for backwards run.

  ! Calculate Sign.
  If (IsBackwards()) Then
    Sign = -1
  Else
    Sign = 1
  End If

  ! Convert time zone of Time.
  TimeL = ChangeTimeZone(Time, WildTimePair%FromTime%Time)

  ! Wild-card times not involving the day of the week.
  If (WildTimePair%FromTime%DayOfWeek == 0) Then

    nNonWildCards = WildTimePair%FromTime%nNonWildCards

    AlteredTime = TimeL

    ! If running backwards, time is 00:00 on the 1st of a month other than January, and year and month are
    ! wild cards, need to compute FromTimeL and ToTime in the previous month to avoid problems with impossible
    ! times in the Gregorian time frame (e.g. 29/1/2001 hh:mm or 31/4/2000 hh:mm).
    If (IsCalendar() .and. IsBackwards()) Then
      If (                        &
        nNonWildCards  == 5 .and. &
        TimeL%Month    /= 1 .and. &
        TimeL%Day      == 1 .and. &
        TimeL%Hour     == 0 .and. &
        TimeL%Minute   == 0 .and. &
        TimeL%Second   == 0 .and. &
        TimeL%Fraction == 0       &
      ) Then
        AlteredTime%Month = AlteredTime%Month - 1
      End If
    End If

    FromTimeL = SubstituteWildCards(AlteredTime, WildTimePair%FromTime)
    ToTimeL   = SubstituteWildCards(AlteredTime, WildTimePair%ToTime)

    ! Note this if test is not 'if ToTimeL >= FromTimeL' because we need to consider the possibility that
    ! SubstituteWildCards has converted impossible times in the Gregorian time frame (e.g. 29/1/2001 hh:mm or
    ! 31/4/2000 hh:mm) to the start of the next month.
    If (WildTimePair%ToTime%Time >= WildTimePair%FromTime%Time) Then
      In = ToTimeL > TimeL .and. TimeL >= FromTimeL
    Else
      In = ToTimeL > TimeL .or.  TimeL >= FromTimeL
    End If

    If (Present(NextTime)) Then

      ! One-time loop.
      Do

        ! No wild cards.
        If (nNonWildCards == 7) Then
          NextTime = InfFutureTime()
          Exit
        Else If (nNonWildCards == 5 .and. .not. IsCalendar()) Then
          NextTime = InfFutureTime()
          Exit
        End If

        ! NextTime is FromTime or ToTime.
        If (ToTimeL >= FromTimeL) Then
          If (FromTimeL > TimeL) Then
            NextTime = FromTimeL
            Exit
          Else If (ToTimeL > TimeL) Then
            NextTime = ToTimeL
            Exit
          End If
        Else
          If (ToTimeL > TimeL) Then
            NextTime = ToTimeL
            Exit
          Else If (FromTimeL > TimeL) Then
            NextTime = FromTimeL
            Exit
          End If
        End If

        ! Increment time, compute the next FromTime and ToTime using the incremented time, and then see if
        ! NextTime is FromTime or ToTime. Note setting of Day to avoid justification altering meaning of
        ! year-month combination.
        If (nNonWildCards == 6) Then
          AlteredTime%Year = AlteredTime%Year + Sign
          AlteredTime%Day  = 1
        Else If (nNonWildCards == 5) Then
          AlteredTime%Month = AlteredTime%Month + Sign
          AlteredTime%Day   = 1
        Else If (nNonWildCards == 4) Then
          AlteredTime%Day = AlteredTime%Day + Sign
        Else If (nNonWildCards == 3) Then
          AlteredTime%Hour = AlteredTime%Hour + Sign
        Else If (nNonWildCards == 2) Then
          AlteredTime%Minute = AlteredTime%Minute + Sign
        Else If (nNonWildCards == 1) Then
          AlteredTime%Second = AlteredTime%Second + Sign
        End If
        Call JustifyTime(AlteredTime)
        FromTimeL = SubstituteWildCards(AlteredTime, WildTimePair%FromTime)
        ToTimeL   = SubstituteWildCards(AlteredTime, WildTimePair%ToTime)
        If (ToTimeL >= FromTimeL) Then
          If (FromTimeL > TimeL) Then
            NextTime = FromTimeL
            Exit
          Else If (ToTimeL > TimeL) Then
            NextTime = ToTimeL
            Exit
          End If
        Else
          If (ToTimeL > TimeL) Then
            NextTime = ToTimeL
            Exit
          Else If (FromTimeL > TimeL) Then
            NextTime = FromTimeL
            Exit
          End If
        End If

        Call Message('UNEXPECTED FATAL ERROR in InWildTimeInterval', 4)

      End Do

    End If

  ! Wild-card times involving the day of the week. Here we (mis-)use the Time_ type by using the day component
  ! to represent the day of the week.
  Else

    AlteredTime     = TimeL
    AlteredTime%Day = DayOfWeek(AlteredTime)

    FromTimeL = WildTimePair%FromTime%Time
    ToTimeL   = WildTimePair%ToTime%Time

    FromTimeL%Year = TimeL%Year
    ToTimeL%Year   = TimeL%Year

    FromTimeL%Month = TimeL%Month
    ToTimeL%Month   = TimeL%Month

    FromTimeL%Day = WildTimePair%FromTime%DayOfWeek
    ToTimeL%Day   = WildTimePair%ToTime%DayOfWeek

    If (ToTimeL >= FromTimeL) Then
      In = ToTimeL > AlteredTime .and. AlteredTime >= FromTimeL
    Else
      In = ToTimeL > AlteredTime .or.  AlteredTime >= FromTimeL
    End If

    If (Present(NextTime)) Then

      If (ToTimeL >= FromTimeL) Then
        If (FromTimeL > AlteredTime) Then
          NextTime     = FromTimeL
          NextTime%Day = TimeL%Day + FromTimeL%Day - AlteredTime%Day
        Else If (ToTimeL > AlteredTime) Then
          NextTime     = ToTimeL
          NextTime%Day = TimeL%Day + ToTimeL%Day - AlteredTime%Day
        Else
          NextTime     = FromTimeL
          NextTime%Day = TimeL%Day + FromTimeL%Day - AlteredTime%Day + 7 * Sign
        End If
      Else
        If (ToTimeL > AlteredTime) Then
          NextTime = ToTimeL
          NextTime%Day = TimeL%Day + ToTimeL%Day - AlteredTime%Day
        Else If (FromTimeL > AlteredTime) Then
          NextTime = FromTimeL
          NextTime%Day = TimeL%Day + FromTimeL%Day - AlteredTime%Day
        Else
          NextTime = FromTimeL
          NextTime%Day = TimeL%Day + ToTimeL%Day - AlteredTime%Day + 7 * Sign
        End If
      End If

      Call JustifyTime(NextTime)

    End If

  End If

End Subroutine InWildTimeInterval

!-------------------------------------------------------------------------------------------------------------

End Module TimeModule

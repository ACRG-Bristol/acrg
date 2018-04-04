! Module: Single Site Module

! Notes: drawing on code from the ADMS 3.0 met input module by Dave Thomson - the code which
!        was used is that in the files MetInput.For, Met.Inc and MetDemo.Met, before the final changes made
!        by Dave Thomson for CERC as given in the files MetInput_For_CERC.For and MetDemo_For_CERC.Met.

! Only the following limited changes have been made to the ADMS code:
!
! Structures/records replaced by types.
! . replaced by %.
! Types and routines put in a single module.
! Logical*1 replaced by Logical.
! Local variable LMetSite removed from the MetInput routine and .true./.false. used
!     directly in the call to SiteInit.
! Saved variables in the MetInput routine returned through the argument list, except
!     for TSampleL, PCorrL, LSeqL and LMetSiteRepresentativeL which are available
!     anyway in the calling routine and so were deleted.
! MetInput routine split into InitMetInput and UpdateMetInput. The ability to print
!     error and warning statistics via CALL Message(50, 0.0, RMessage, 0, 0, ' ') has
!     not been retained in either routine - instead this call can be invoked directly
!     as required. The argument lists of the two routines have been reduced to the
!     variables needed. In particular LFirstTime and LMessage are no longer needed
!     (and have been replaced by .true. and .false. in the argument list to the
!     routine Check).
! Met_ renamed SSMet_.
! Char in argument list of message routine changed from length 60 to length *.
! NDetailedFr in MetFr_ changed to an integer.
! Calculation of MetOut%LocalMeanTime from MetOut%Hour, MetOut%TimeZone and
!     Site%Long/15.0 corrected. Both the time zone and longitude corrections had the
!     wrong sign.
! Private and public attributes added.
! Adaptions made so as to work in calm conditions. Much of the code already worked in
!     calms, but four final changes needed to be made. These are (i) removing message
!     for calms in ProcessMet and the following Return, (ii) altering WStarCalc and
!     its argument list to work in calms, (iii) altering SigmaThetaCalc to work in
!     calms by putting U10 = 0.5 in calms, and (iv) not applying RecipLMOMax in calms
!     (to maintain the convention RecipLMO = +/-200000). These changes don't alter
!     results for non-calms. (iii) could be improved and made continuous, at expense
!     of altering existing results. However the calculated SigmaTheta is not used in
!     Name III.
!     $$ Its just possible that problems could occur for very small non-zero winds.
! Sevice module now used.
! Message routine substantially revised and renamed SSMessage. Consequent changes on
!     rest of code: (i) test for variable name 'STATION DCNN' added before call to
!     message routine for message 14 (previously this test was in the message
!     routine). Calls to messages 24 and 25 removed (these messages have been
!     removed in the message routine).
! Length of character string MetFile changed from 92 to the global parameter
!     MaxFileNameLength.
! Minor changes to open statements.
! Wind direction calculation altered so that the wind direction at the input height
!     equals the input wind direction. Previously the surface wind direction was set
!     equal to the input wind direction and so the wind direction at the input height
!     was different if there was turning of wind with height.

Module SingleSiteModule

! This module provides code for reading single site met files and doing some initial
! processing.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: SSMet_         ! Met data.
Public  :: MetFr_         ! Met frequency data.
Public  :: Site_          ! Site details.
Public  :: MetFileState_  ! Information on the met file and the parameters associated
                          ! with reading from it.
Public  :: InitMetInput   !
Public  :: UpdateMetInput !

!-------------------------------------------------------------------------------------------------------------

!----------------------------- Parameters -----------------------------!

! MaxVariables: maximum number of variables which can be read from the
!    met file.
! NVariableNames: number of recognized met-file variable names, excluding
!    names of the detailed frequency variables.
! NReadAtATime: Number of lines of data read from the met file at a time.

      INTEGER*4 MaxVariables, NVariableNames, NReadAtATime
      PARAMETER (MaxVariables = 60, NVariableNames = 35, &
                 NReadAtATime = 55)

!-------------------------------- Met ---------------------------------!

! Structure containing met data.

! For all the structure elements, -999.0 indicates missing data.
! In the dispersion area (as opposed to at the met site) the 'wind
! measurement' refers to a hypothetical measurement of the geostrophic
! wind. It is the geostropic wind which is assumed to be the same
! at the met site and in the dispersion area, and so it is this wind
! which is the input into the dispersion area met processing (if
! dispersion area differs from the met site).

! 1) Met variables:
! U: Wind speed at measurement height (friction velocity if measurement
!    height = 0; geostophic wind if measurement  height = 1000.0).
! UGStar: Geostrophic wind speed/friction velocity.
! UStar: Friction velocity.
! UG: Geostrophic wind speed.
! Phi: Wind direction (angle wind is coming from in degrees clockwise from
!    north) at measurement height (surface wind direction if measurement
!    height = 0; geostophic wind direction if measurement height = 1000.0).
! DeltaPhi: Geostrophic wind direction minus surface wind direction (both
!    directions measured in degrees clockwise from north).
! Phi0: Surface wind direction (angle wind is coming from in degrees
!    clockwise from north).
! PhiG: Geostrophic wind direction (angle wind is coming from in degrees
!    clockwise from north).
! FTheta0: Surface sensible heat flux.
! RecipLMO: 1/Monin-Obukhov length (set to +/-200000 in calms).
! ThetaStar: Temperature scale (positive for FTheta0 < 0 and set to +/-200000
!    in calms).
! H: Boundary layer depth.
! WStar: Convective velocity scale if FTheta0 > 0, zero if FTheta0 <or= 0.
! S: Sine of the solar elevation.
! Cl: Cloud amount (oktas).
! K: Incoming solar radiation.
! T0C: Near surface temperature (degrees C).
! T0K: Near surface temperature (degrees K).
! NU: Buoyancy frequency above the boundary layer.
! DeltaTheta: Temperature jump across the boundary layer top.
! P: Precipitation rate (mm/hour).
! TSea: Sea surface temperature (degrees C).
! DeltaT: Near surface temperature over land minus sea surface temperature.
! SigmaThetaDeg: Standard deviation of changes in mean wind direction
!    (degrees).
! SigmaThetaRad: Standard deviation of changes in mean wind direction
!    (radians).
! Q0: Near surface specific humidity.
! RH0: Near surface relative humidity (percent).
! RHU: Relative humidity just above the boundary layer (percent).
! DRHDZU: d(Relative humidity)/dz above the boundary layer (percent/m).
! LambdaE: Surface latent heat flux.

! 2) Site variables:
! Z: Height of wind measurement at met site (1,000.0 is used to indicate
!    geostrophic wind, 0.0 to indicate friction velocity).
! ZD: Height of wind measurement in dispersion area.
! Z0: Roughness length at met site.
! Z0D: Roughness length in dispersion area.
! R: Surface albedo at met site.
! RD: Surface albedo in dispersion area.
! Alpha: Modified Priestley-Taylor parameter at met site (as defined in
!    Holtslag and van Ulden, 1983, J. Clim. Appl. Met., vol 22, 517-529).
! AlphaD: Modified Priestley-Taylor parameter in dispersion area (as defined
!    in Holtslag and van Ulden, 1983, J. Clim. Appl. Met., vol 22, 517-529).
! RecipLMOMax: Maximum value of 1/LMO at met site, with values > 1.0E6
!    indicating no maximum.
! RecipLMOMaxD: Maximum value of 1/LMO in dispersion area, with values > 1.0E6
!    indicating no maximum.

! 3) Time variables:
! Hour: Hour of day (e.g. 5.30am = 5.5).
! Day: Day of year (e.g. 1st Jan = 1.0).
! Year: Year (full 4 digit number, e.g. 1999.0).
! TimeZone: The time zone in which the above time data are expressed,
!    specified as hours ahead of GMT (so that time zones to the east of
!    the Greenwich meridian will tend to have a positive value and time
!    zones to the west will tend to have a negative value). If this
!    variable is missing local mean time is assumed, i.e. GMT + longitude
!    in degrees (east positive)/15. The value specified is the time zone
!    applicable to all of the variables HOUR, DAY and YEAR (to see how the
!    time zone can be important for identifying the day and year note
!    that, for example, in one time zone the time may be 11.30pm on 31st
!    December 1999, while in another it may be 12.30am on 1st January
!    2000).
! LocalMeanTime: Hour of day calculated as GMT + longitude in degrees (east
!    positive)/15).

      Type :: SSMet_
        REAL*4 U, UGStar, UStar, UG, Phi, DeltaPhi, Phi0, PhiG, FTheta0, &
               RecipLMO, ThetaStar, H, WStar, S, Cl, K, T0C, T0K, NU,    &
               DeltaTheta, P, TSea, DeltaT, SigmaThetaDeg,               &
               SigmaThetaRad, Q0, RH0, RHU, DRHDZU, LambdaE
        REAL*4 Z, ZD, Z0, Z0D, R, RD, Alpha, AlphaD, RecipLMOMax, &
               RecipLMOMaxD
        REAL*4 Hour, Day, Year, TimeZone, LocalMeanTime
      End Type

!------------------------------- MetFr --------------------------------!

! Structure containing met frequency data.

! Fr: Frequency with which a given set of met conditions occurs
!    (arbitrary units, e.g. percentage of occasions or number of hours
!    per year). -999.0 indicates missing data.
! FrPresent: Indicates that frequency information is present in the met
!    file.
! NDetailedFr: Number of variables in the met file giving detailed
!    frequency information.
! DetailedFr: Detailed frequency information giving the frequency with
!    which a given set of met conditions occurs for certain months of the
!    year or times of day (same units as for Fr).
! HourLowerLimits, HourUpperLimits: Hour limits for detailed frequency
!    data, for example (1.50, 4.50) indicates 1.30am to 4.30am.
! MonthLowerLimits, MonthUpperLimits: Month limits for detailed frequency
!    data in inclusive whole months, for example (1.0, 3.0) indicates
!    January to March inclusive.
! TimeZone: The time zone in which the hour and month limits are
!    expressed, specified as hours ahead of GMT (so that time zones to the
!    east of the Greenwich meridian will tend to have a positive value and
!    time zones to the west will tend to have a negative value). -999.0
!    indicates local mean time, i.e. GMT + longitude in degrees (east
!    positive)/15.

      Type :: MetFr_
        REAL*4 Fr
        LOGICAL FrPresent
        INTEGER*4 NDetailedFr
        REAL*4 DetailedFr(MaxVariables),       &
               HourLowerLimits(MaxVariables),  &
               HourUpperLimits(MaxVariables),  &
               MonthLowerLimits(MaxVariables), &
               MonthUpperLimits(MaxVariables), &
               TimeZone
      End Type

!-------------------------------- Site --------------------------------!

! Structure containing site details.

! 1) Site location:
! Long: Longitude (degrees, east positive).
! Lat: Latitude (degrees, north positive).
! AbsF, SignF: Absolute value and sign of Coriolis parameter.
! MetSite: Indicates whether the structure refers to the met site or the
!    dispersion area.

! 2) Site characteristics:
! Alpha, RecipLMOMax, R, Z, Z0: Values of Alpha, RecipLMOMax, R, Z, and
!    Z0. Missing data is indicated by -999.0. (If values are missing,
!    then values from the met file will be used instead. If there are no
!    values in the met file either, then, for Alpha and R, default values
!    will be used.)

! 3) Site limits for UStar calculation:
! UG2B: minimum allowed value of U**2/(-surface buoyancy flux), for U a
!    geostrophic wind.
! U2BStar: minimum allowed value of U**2/(buoyancy scale), for U a surface
!    layer wind at height ZLast with roughness length Z0Last.
! U3B: minimum allowed value of U**3/(-surface buoyancy flux), for U a
!    surface layer wind at height ZLast with roughness length Z0Last.
! ZLast, Z0Last: Z and Z0 from the last set of met data (used to assess
!    whether UG2B, U2BStar and UB3 need recalculating).

! 4) Site history:
! UStar, FTheta0, S, T0K, NU: Values of UStar, FTheta0, S, T0K and NU
!    for previous hours.

      Type :: Site_
        REAL*4 Long, Lat, AbsF, SignF
        LOGICAL MetSite
        REAL*4 Alpha, RecipLMOMax, R, Z, Z0
        REAL*4 UG2B, U2BStar, U3B, ZLast, Z0Last
        REAL*4 UStar(24), FTheta0(24), S(24), T0K(24), NU(24)
      End Type

!---------------------------- MetFileState ----------------------------!

! Structure containing information on the met file and the parameters
! associated with reading from it.

! Codes: Code numbers to indicate what the data in the met file
!    represents - the code numbers correspond to the items of data listed
!    in VariableNames (see subroutine ReadSetup), with -1 denoting
!    detailed frequency data.
! FrCodes: Code numbers to indicate what the data in the met file
!    represents - the code numbers correspond to the detailed frequency
!    data in the met file.
! NAvailable: Number of hours' data in the arrays DetailedFr and
!    SortedData which have not yet been returned by the subroutine ReadMet.
! NRead: Number of hours' data read from the met file and stored in the
!    arrays DetailedFr and SortedData on the last occasion the file was
!    read.
! NRec: Number of records of the met file which have been read.
! NVar: Number of variables given in the met file for each
!    each hour of data (or MaxVariables if this is smaller).
! DetailedFr: Array containing up to NReadAtATime hours' detailed
!    frequency data from the met file. -999.0 indicates missing values.
! SortedData: Array containing up to NReadAtATime hours' met data and
!    frequency data (but not detailed frequency data) from the met file.
!    -999.0 indicates missing values.
! MetFile: Name of met file.

      Type :: MetFileState_
        INTEGER*4 Codes(MaxVariables), FrCodes(MaxVariables), &
                  NAvailable, NRead, NRec, NVar
        REAL*4 DetailedFr(MaxVariables,NReadAtATime), &
               SortedData(NVariableNames,NReadAtATime)
        CHARACTER(MaxFileNameLength) MetFile
      End Type

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------


!***********************************************************************
! This module reads data from the met file and uses these data to      *
! estimate values of the various meteorological quantitites required   *
! for running the dispersion model. The module should be called once   *
! for each hour's data.                                                *
!                                                                      *
!      Documentation: ADMS papers P05/01, P05/03 and P05/04            *
!        Input Files: Met file - contains the met data to be read in.  *
!                     Met-limits file - contains the limits within     *
!                        which data in the met file and data           *
!                        transferred into the met input module from    *
!                        other parts of ADMS must lie.                 *
!       Output Files: Log file - contains output messages              *
!     Included Files: Met.inc - contains definitions of structures and *
!                        parameters.                                   *
!***********************************************************************
!                                                                      *
! Module structure:                                                    *
!                                                                      *
!            (Check---ReadLimits                                       *
!            (                                                         *
!            (ReadMetSetUp                                             *
! MetInput---(ReadMet---ReadData                                       *
!            (                                                         *
!            (QualityControlMetSetUp                                   *
!            (QualityControlMet                                        *
!            (                                                         *
!            (SiteInit                                                 *
!            (SiteIncrement                                            *
!            (ProcessMet---(see below)                                 *
!                                                                      *
!              (CORIOLISCALC                                           *
!              (              (UCALC2A                                 *
!              (SETUPLIMITS---(UCALC3A                                 *
!              (                                                       *
!              (FTHETA0LIMIT                                           *
!              (                                                       *
!              (SCALC                                                  *
!              (                                                       *
!              (KCalc                                                  *
!              (                                                       *
!              (         (SURFACELAYERSTABLE & UNSTABLE                *
!              (UCALC1---(ROSSBYCOEFFSTABLE & UNSTABLE                 *
!              (                                                       *
!              (         (SURFACELAYERSTABLE & UNSTABLE                *
!              (UCALC2---(ROSSBYCOEFFSTABLE & UNSTABLE                 *
!              (                                                       *
!              (FTheta0Calc---(see below)                              *
!              (                                                       *
!              (UGCALC---ROSSBYCOEFFSTABLE & UNSTABLE                  *
!              (                                                       *
!              (ANGLES                                                 *
! ProcessMet---(                                                       *
!              (HStableCalc---(see below)                              *
!              (                                                       *
!              (                (HStableCalc---(see below)             *
!              (HUnstableCalc---(HStableBasicCalc                      *
!              (                (BLGROWTH                              *
!              (HLimits                                                *
!              (                                                       *
!              (DELTATHETACALC                                         *
!              (                                                       *
!              (WSTARCALC                                              *
!              (                                                       *
!              (LambdaECalc                                            *
!              (                                                       *
!              (         (SatVapourPressure                            *
!              (Q0Calc---(MixingRatio                                  *
!              (                                                       *
!              (          (FTheta0Calc---(see below)                   *
!              (FillSeq---(INTERPOLATION                               *
!              (                                                       *
!              (FillNonSeq---FTheta0Calc---(see below)                 *
!                                                                      *
!               (THETASTARLIMIT                                        *
!               (                                                      *
!               (         (SURFACELAYERSTABLE & UNSTABLE               *
!               (UCALC2---(ROSSBYCOEFFSTABLE & UNSTABLE                *
! FTheta0Calc---(                                                      *
!               (         (SURFACELAYERSTABLE                          *
!               (UCALC3---(ROSSBYCOEFFSTABLE                           *
!               (                                                      *
!               (FTheta0DayTimeCalc                                    *
!                                                                      *
!                         (SURFACELAYERSTABLE & UNSTABLE               *
!               (UCALC2---(ROSSBYCOEFFSTABLE & UNSTABLE                *
! HStableCalc---(                                                      *
!               (HStableBasicCalc                                      *
!                                                                      *
! There is also a message subroutine which can be called from any part *
! of the module.                                                       *
!                                                                      *
! Transfer of data between subroutines is done by argument list.       *
!                                                                      *
! All units are SI units except where stated.                          *
!***********************************************************************

      SUBROUTINE InitMetInput(                                        &
         RLong, RLat, TSample,                                        &
         Z, Z0, Z0D, Alpha, AlphaD, R, RD, RecipLMOMax, RecipLMOMaxD, &
         PCorr,                                                       &
         MetFile, MetLimitsFile,                                      &
         LSeq, LMetSiteRepresentative,                                &
         MetFr,                                                       &
         IHour, NInad, NCalm, NTotal, FInad, FCalm, FTotal,           &
         LFatalError,                                                 &
         MinMet, MaxMet, Site, SiteD, MetFileState) ! ex local saved variables
!-----------------------------------------------------------------------
! Input variables:
! RLong: Longitude (degrees, east positive).
! RLat: Latitude (degrees, north positive).
! TSample: Sampling time.
! Z: Height of the wind measurement at met site (1,000.0 is used to
!    indicate geostrophic wind, 0.0 to indicate friction velocity).
! Z0: Roughness length at met site.
! Z0D: Roughness length in dispersion area.
! Alpha: Modified Priestley-Taylor parameter at met site (as defined in
!    Holtslag and van Ulden, 1983, J. Clim. Appl. Met., vol 22, 517-529).
! AlphaD: Modified Priestley-Taylor parameter in dispersion area (as defined
!    in Holtslag and van Ulden, 1983, J. Clim. Appl. Met., vol 22, 517-529).
! R: Surface albedo at met site.
! RD: Surface albedo in dispersion area.
! RecipLMOMax: Maximum value of 1/LMO at met site, with values > 1.0E6
!    indicating no maximum.
! RecipLMOMaxD: Maximum value of 1/LMO in dispersion area, with values > 1.0E6
!    indicating no maximum.
! PCorr: Precipitation correction factor (ratio of average rainfall in the
!    dispersion area to that at the measurement site).
! MetFile: Name of met file.
! MetLimitsFile: Name of met-limits file.
! LSeq: Indicates that the data provided in the met file are hourly
!    sequential.
! LMetSiteRepresentative: Indicates that no distinction is to be made
!    between met site and dispersion area properties.
! LFirstTime: Indicates that the met input module is being called for the
!    first time.
! LMessage: Indicates that the met input module is to issue a message
!    giving statistics of errors and warnings, but is not to process any
!    hours of data.
! Apart from LFirstTime and LMessage, these input variables only matter
! the first time the module is called.
! For Z, Z0, Z0D, Alpha, AlphaD, R, RD, RecipLMOMax and RecipLMOMaxD,
! missing values are indicated by -999.0. (If values are missing,
! then values from the met file will be used instead. If there are no
! values in the met file either, then, for Alpha and R, default values
! will be used.)
!
! Output variables:
! MetAsRead: Met data as read from the met file.
! ProcMet: Processed met data.
! MetFr: Frequency data associated with the met conditions.
! IHour: Indicates that the data being considered are the IHour-th hour's
!    data from the met file.
! NInad: Cumulative number of hours for which LInadequateData = .TRUE..
! NCalm: Cumulative number of hours for which LCalm = .TRUE..
! NTotal: Cumulative number of hours for which neither LInadequateData nor
!    LCalm is set.
! FInad: Cumulative frequency of occasions for which LInadequateData = .TRUE..
! FCalm: Cumulative frequency of occasions for which LCalm = .TRUE..
! FTotal: Cumulative frequency of occasions for which neither LInadequateData
!    nor LCalm is set.
! LFatalError: Indicates that a fatal error has occurred and that the met
!    input module should not be called again.
! LNoMoreData: Indicates that there are no more data in the met file.
! LInadequateData: Indicates that a particular hours' data from the met
!    file are inadequate or not sensible.
! LCalm: Indicates that U = 0 for a particular hours' data.
!
! Local variables:
! MinMet, MaxMet: Limits within which the met variables must lie to be
!    regarded as sensible.
! MetAsReadD: Hypothetical set of met data as might be measured in the
!    dispersion area.
! Site, SiteD: Site characteristics at the met site and in the dispersion
!    area (including some information on the recent met).
! MetFileState: The state of the met file.
! RMessage: Array required in the argument list of the subroutine SSMessage.
!    It has no significance in this subroutine.
! LMetSite: Indicates that the met site, as opposed to the dispersion
!    area, is being considered.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      REAL*4 RLong, RLat, TSample, Z, Z0, Z0D, Alpha, AlphaD, &
             R, RD, RecipLMOMax, RecipLMOMaxD, PCorr
      CHARACTER(MaxFileNameLength) MetFile
      CHARACTER*32 MetLimitsFile
      LOGICAL LSeq, LMetSiteRepresentative
!..........Output variables. (Some are actually InOut)
      Type(MetFr_) MetFr
      INTEGER*4 IHour, NInad, NCalm, NTotal
      REAL*4 FInad, FCalm, FTotal
      LOGICAL LFatalError
!..........Variables which need to be remembered between calls (originally local saved
!          variables; now InOut variables).
      Type(SSMet_) MinMet, MaxMet
      Type(Site_) Site, SiteD
      Type(MetFileState_) MetFileState
!      SAVE MinMet, MaxMet, Site, SiteD, MetFileState, TSampleL, PCorrL, &
!           LSeqL, LMetSiteRepresentativeL
!..........Local variables.
      REAL*4 RMessage(1)

!==========Initialise messages.

      CALL SSMessage(0, 0.0, RMessage, 0, 0, ' ')

!==========Check values obtained from other parts of ADMS are sensible
!          and not missing when needed.

      CALL Check(                               &
       RLat, TSample, Z, Z0, Z0D, PCorr,        &
       MetLimitsFile,                           &
       LMetSiteRepresentative, .true., .false., &
       LFatalError)
      IF (LFatalError) RETURN

!==========Initialise everything.

      IHour  = 0
      NInad  = 0
      FInad  = 0.0
      NCalm  = 0
      FCalm  = 0.0
      NTotal = 0
      FTotal = 0.0
      CALL ReadMetSetUp(MetFile, MetFr, MetFileState, LFatalError)
      IF (LFatalError) RETURN
      CALL QualityControlMetSetUp(MetFr, MetLimitsFile, LSeq, MinMet, &
                                  MaxMet, LFatalError)
      IF (LFatalError) RETURN
      CALL SiteInit(RLong, RLat, Z, Z0, R, Alpha, &
                    RecipLMOMax, .true., Site)
      CALL SiteInit(RLong, RLat, 1000.0, Z0D, RD, AlphaD, &
                    RecipLMOMaxD, .false., SiteD)

      RETURN
      END Subroutine InitMetInput

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

      SUBROUTINE UpdateMetInput(                                      &
         TSample,                                                     &
         PCorr,                                                       &
         LSeq, LMetSiteRepresentative,                                &
         MetAsRead, ProcMet, MetFr,                                   &
         IHour, NInad, NCalm, NTotal, FInad, FCalm, FTotal,           &
         LFatalError, LNoMoreData, LInadequateData, LCalm,            &
         MinMet, MaxMet, Site, SiteD, MetFileState) ! ex local saved variables
!-----------------------------------------------------------------------
! Input variables:
! RLong: Longitude (degrees, east positive).
! RLat: Latitude (degrees, north positive).
! TSample: Sampling time.
! Z: Height of the wind measurement at met site (1,000.0 is used to
!    indicate geostrophic wind, 0.0 to indicate friction velocity).
! Z0: Roughness length at met site.
! Z0D: Roughness length in dispersion area.
! Alpha: Modified Priestley-Taylor parameter at met site (as defined in
!    Holtslag and van Ulden, 1983, J. Clim. Appl. Met., vol 22, 517-529).
! AlphaD: Modified Priestley-Taylor parameter in dispersion area (as defined
!    in Holtslag and van Ulden, 1983, J. Clim. Appl. Met., vol 22, 517-529).
! R: Surface albedo at met site.
! RD: Surface albedo in dispersion area.
! RecipLMOMax: Maximum value of 1/LMO at met site, with values > 1.0E6
!    indicating no maximum.
! RecipLMOMaxD: Maximum value of 1/LMO in dispersion area, with values > 1.0E6
!    indicating no maximum.
! PCorr: Precipitation correction factor (ratio of average rainfall in the
!    dispersion area to that at the measurement site).
! MetFile: Name of met file.
! MetLimitsFile: Name of met-limits file.
! LSeq: Indicates that the data provided in the met file are hourly
!    sequential.
! LMetSiteRepresentative: Indicates that no distinction is to be made
!    between met site and dispersion area properties.
! LFirstTime: Indicates that the met input module is being called for the
!    first time.
! LMessage: Indicates that the met input module is to issue a message
!    giving statistics of errors and warnings, but is not to process any
!    hours of data.
! Apart from LFirstTime and LMessage, these input variables only matter
! the first time the module is called.
! For Z, Z0, Z0D, Alpha, AlphaD, R, RD, RecipLMOMax and RecipLMOMaxD,
! missing values are indicated by -999.0. (If values are missing,
! then values from the met file will be used instead. If there are no
! values in the met file either, then, for Alpha and R, default values
! will be used.)
!
! Output variables:
! MetAsRead: Met data as read from the met file.
! ProcMet: Processed met data.
! MetFr: Frequency data associated with the met conditions.
! IHour: Indicates that the data being considered are the IHour-th hour's
!    data from the met file.
! NInad: Cumulative number of hours for which LInadequateData = .TRUE..
! NCalm: Cumulative number of hours for which LCalm = .TRUE..
! NTotal: Cumulative number of hours for which neither LInadequateData nor
!    LCalm is set.
! FInad: Cumulative frequency of occasions for which LInadequateData = .TRUE..
! FCalm: Cumulative frequency of occasions for which LCalm = .TRUE..
! FTotal: Cumulative frequency of occasions for which neither LInadequateData
!    nor LCalm is set.
! LFatalError: Indicates that a fatal error has occurred and that the met
!    input module should not be called again.
! LNoMoreData: Indicates that there are no more data in the met file.
! LInadequateData: Indicates that a particular hours' data from the met
!    file are inadequate or not sensible.
! LCalm: Indicates that U = 0 for a particular hours' data.
!
! Local variables:
! MinMet, MaxMet: Limits within which the met variables must lie to be
!    regarded as sensible.
! MetAsReadD: Hypothetical set of met data as might be measured in the
!    dispersion area.
! Site, SiteD: Site characteristics at the met site and in the dispersion
!    area (including some information on the recent met).
! MetFileState: The state of the met file.
! RMessage: Array required in the argument list of the subroutine SSMessage.
!    It has no significance in this subroutine.
! LMetSite: Indicates that the met site, as opposed to the dispersion
!    area, is being considered.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      REAL*4 TSample, PCorr
      LOGICAL LSeq, LMetSiteRepresentative
!..........Output variables. (Some are actually InOut)
      Type(SSMet_) MetAsRead, ProcMet
      Type(MetFr_) MetFr
      INTEGER*4 IHour, NInad, NCalm, NTotal
      REAL*4 FInad, FCalm, FTotal
      LOGICAL LFatalError, LNoMoreData, LInadequateData, LCalm
!..........Variables which need to be remembered between calls (originally local saved
!          variables; now InOut variables).
      Type(SSMet_) MinMet, MaxMet
      Type(Site_) Site, SiteD
      Type(MetFileState_) MetFileState
!      SAVE MinMet, MaxMet, Site, SiteD, MetFileState, TSampleL, PCorrL, &
!           LSeqL, LMetSiteRepresentativeL
!..........Local variables.
      Type(SSMet_) MetAsReadD
      REAL*4 RMessage(1)

!==========Increment IHour and store in message routine.

      IHour = IHour + 1
      CALL SSMessage(-1, 0.0, RMessage, IHour, 0, ' ')

!==========Read Met.

      CALL ReadMet(MetFr, MetFileState, MetAsRead, LNoMoreData, &
                   LFatalError)
      IF (LFatalError) RETURN
      IF (LNoMoreData) THEN
        CALL SSMessage(26, 0.0, RMessage, 0,0, ' ')
        RETURN
      ENDIF

!==========Increment Site and SiteD to reflect the fact that time has
!          advanced by one hour.

      CALL SiteIncrement(LSeq,Site)
      CALL SiteIncrement(LSeq,SiteD)

!==========Quality control data (and set Frequency = 1 if frequency is
!          missing).

      CALL QualityControlMet(MinMet, MaxMet, MetAsRead, MetFr, &
                             LInadequateData)
      IF (MetFr%Fr.EQ.-999.0) MetFr%Fr = 1.0
      IF (LInadequateData) THEN
        NInad = NInad + 1
        FInad = FInad + MetFr%Fr
        RETURN
      ENDIF

!==========Call BLP submodule.

      CALL ProcessMet(MetAsRead, TSample, LSeq, Site, ProcMet, &
                      LInadequateData, LCalm)
      IF (LInadequateData) THEN
        NInad = NInad + 1
        FInad = FInad + MetFr%Fr
        RETURN
      ENDIF
      IF (LCalm) THEN
        NCalm = NCalm + 1
        FCalm = FCalm + MetFr%Fr
        RETURN
      ENDIF
      NTotal = NTotal + 1
      FTotal = FTotal + MetFr%Fr

!==========Correct for differences between met site and dispersion area.

      IF (.NOT.LMetSiteRepresentative) THEN
        MetAsReadD     = MetAsRead
        MetAsReadD%U   = ProcMet%UG
        MetAsReadD%Phi = ProcMet%PhiG
        MetAsReadD%ZD  = 1000.0
        IF (MetAsRead%P.NE.-999.0 .AND. PCorr.NE.-999.0) MetAsReadD%P = &
           MetAsRead%P*PCorr
        IF (MetAsRead%FTheta0.NE.-999.0 .OR.   &
            MetAsRead%RecipLMO.NE.-999.0) THEN
          IF (ProcMet%FTheta0.GT.0.0) THEN
            MetAsReadD%RecipLMO  = -999.0
            MetAsReadD%FTheta0   = ProcMet%FTheta0
            MetAsReadD%ThetaStar = -999.0
          ELSE
            MetAsReadD%RecipLMO  = -999.0
            MetAsReadD%FTheta0   = -999.0
            MetAsReadD%ThetaStar = ProcMet%ThetaStar
          ENDIF
        ENDIF
        CALL ProcessMet(MetAsReadD, TSample, LSeq, SiteD, ProcMet, &
                        LInadequateData, LCalm)
      ENDIF

      RETURN
      END Subroutine UpdateMetInput

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

Subroutine SSMessage(N, R, RMessage, I1, I2, Char)
! This subroutine generates messages. When called with N = 0 it initialises the arrays MessageCount and
! MessageLimits. When called with N = 50 it issues a message giving statistics on errors and warnings.

  Implicit None
  ! Argument list:
  Integer,      Intent(In) :: N           ! Number of message to be generated.
  Real(Std),    Intent(In) :: R           ! A real variable to be included in the message.
  Real(Std),    Intent(In) :: RMessage(:) ! An array of real variables to be included in the message.
  Integer,      Intent(In) :: I1          ! An integer variable to be included in the message.
  Integer,      Intent(In) :: I2          ! An integer variable to be included in the message or to indicate
                                          ! the number of characters in Char to be printed.
  Character(*), Intent(In) :: Char        ! A character variable to be included in the message.
  ! Local parameters:
  Integer, Parameter :: NumberOfMessages = 31 ! The number of message types.
  ! Locals:
  Integer                   :: MessageCount(NumberOfMessages)
  Integer                   :: MessageLimits(NumberOfMessages)
  Integer                   :: IHour
  Integer                   :: Loop
  Character(MaxCharLength3) :: Line
  Logical                   :: LMessages
  Save :: MessageCount
  Save :: MessageLimits
  Save :: IHour
  ! MessageCount  :: Array giving the number of messages of each type issued so far.
  ! MessageLimits :: Array giving the maximum number of messages of each type to be issued.
  ! IHour         :: Hour of met data being considered.
  ! Loop          :: Loop index.
  ! Line          :: Text of message.
  ! LMessages     :: Indicates whether any error or warning messages have been issued.

  Data MessageLimits/                            &
          1,  1,  1,  1,  1,  1,  1,  1,  1,  5, & !  1-10
          5,  5,  5, 50,  1,  0,  5,  5,  0,  5, & ! 11-20
          1,  1,  5,  0,  0,  0,  1,  1,  1,  1, & ! 21-30
         50                                      & ! 31
       /

!***********************************************************************
!          When this subroutine is called with N = 0, initialise       *
!          messages.                                                   *
!***********************************************************************

  If (N == 0) Then

    Do Loop = 1, NumberOfMessages
      MessageCount(Loop) = 0
    End Do

  Else If (N == -1) Then

    IHour = I1

!***********************************************************************
!          When this subroutine is called with N = 50, issue message   *
!          giving statistics on errors and warnings.                   *
!***********************************************************************

  Else If (N == 50) Then

    Write (6, *) 'Met input module. Summary of errors and warnings:' ! $$ use call message or remove
    Write (6, *) 'Error/warning number         Number of occasions'
    LMessages = .FALSE.
    Do Loop = 1,23
      If (MessageCount(Loop) > 0) Then
        If (Loop <= 9) Then
          Line = '    MET'                                               // &
                 Trim(Int2Char(Loop, FormatString = 'I1'))               // &
                 '                        '                              // &
                 Trim(Int2Char(MessageCount(Loop), FormatString = 'I8'))
          Write (6, *) Trim(Line)
        Else
          Line = '    MET'                                               // &
                 Trim(Int2Char(Loop, FormatString = 'I2'))               // &
                 '                       '                               // &
                 Trim(Int2Char(MessageCount(Loop), FormatString = 'I8'))
          Write (6, *) Trim(Line)
        End If
        LMessages = .TRUE.
      End If
    End Do
    If (.not.LMessages) Then
      Write (6, *) 'No errors or warnings issued.'
    End If

!***********************************************************************
!          Update MessageCount.                                        *
!***********************************************************************

  Else

    MessageCount(N) = MessageCount(N) + 1

!***********************************************************************
!          Write message.                                              *
!***********************************************************************

    If (MessageCount(N) <= MessageLimits(N)) Then

      If (N == 1) Then

        Line = 'ERROR: Error ' // Char(1:I2)

      Else If (N == 2) Then

        Line = 'ERROR: The value of '   // &
               Char(1:I2)               // &
               ' is inappropriate, '    // &
               Char(1:I2)               // &
               ' = '                    // &
               Trim(Std2Char(R))        // &
               ' - program terminated.'

      Else If (N == 3) Then

        Line = 'ERROR: The value of '                    // &
               Char(1:I2)                                // &
               ' is not available - program terminated.'

      Else If (N == 4) Then

        Line = 'ERROR: The time zones used in defining the detailed '          // &
               'frequency met data are not all the same - program terminated.'

      Else If (N == 5) Then

        Line = 'ERROR: The sequential data flag is set but there are ' // &
               'frequency data in the met file - program terminated.'

      Else If (N == 6) Then

        Line = 'ERROR: The month or hour limits for the detailed ' // &
               'frequency data lie outside the range 0 to 12 or '  // &
               '0 to 24 respectively - program terminated.'

      Else If (N == 7) Then

        Line = 'ERROR: The values of PCORR and the met-site-in-dispersion-area ' // &
               'flag are inconsistent - program terminated.'

      Else If (N == 8) Then

        Line = 'ERROR: The message and first-time flags are both set ' // &
               '- program terminated.'

      Else If (N == 9) Then

        Line = 'ERROR: The met file contains both of the variables ' // &
               'TIME ZONE and THOUR - program terminated.'

      Else If (N == 10) Then

        Line = 'ERROR: Hour number '                                           // &
               Trim(Int2Char(IHour))                                           // &
               '. The detailed frequency data obtained from the met file are ' // &
               'inappropriate. The '                                           // &
               Trim(Int2Char(I1)) // Int2Ordinal(I1)                           // &
               ' detailed frequency variables has the value '                  // &
               Trim(Std2Char(R))                                               // &
               ' - it is impossible to use this data.'

      Else If (N == 11) Then

        Line = 'ERROR: Hour number '                            // &
               Trim(Int2Char(IHour))                            // &
               '. The value of '                                // &
               Char(1:I2)                                       // &
               ' obtained from the met file is inappropriate: ' // &
               Char(1:I2)                                       // &
               ' = '                                            // &
               Trim(Std2Char(R))                                // &
               ' - it is impossible to use this data.'

      Else If (N == 12) Then

        Line = 'ERROR: Hour number '                  // &
               Trim(Int2Char(IHour))                  // &
               '. '                                   // &
               Char(1:I2)                             // &
               ' is not available from the met file ' // &
               '- it is impossible to use this data.'

      Else If (N == 13) Then

        Line = 'ERROR: Hour number '                                            // &
               Trim(Int2Char(IHour))                                            // &
               '. There is insufficient information to estimate the heat flux ' // &
               '- it is impossible to use this data.'

      Else If (N == 14) Then

        Line = 'WARNING: The met file contains an unrecognised variable name: ' // &
               Char(1:I2)                                                       // &
               '. The values given for this variable will be ignored.'

      Else If (N == 15) Then

        Line = 'WARNING: The number of variables given in the met file for each hour is ' // &
               Trim(Int2Char(I1))                                                         // &
               ' and is too large. Except for the first '                                 // &
               Trim(Int2Char(I2))                                                         // &
               ' variables, the data are ignored.'

      Else If (N == 16) Then

        Line = 'WARNING: The Coriolis parameter has been calculated as '                   // &
               Trim(Std2Char(RMessage(1)))                                                 // &
               '. This is less than (and has been replaced by) the minimum allowed value ' // &
               Trim(Std2Char(RMessage(2)))                                                 // &
               '.'

      Else If (N == 17) Then

        Line = 'WARNING: Hour number '                                                    // &
               Trim(Int2Char(IHour))                                                      // &
               '. The USTAR calculation has not converged. The calculation will proceed ' // &
               'using the value of USTAR from the final iteration. Please report the '    // &
               'following values: U = '                                                   // &
               Trim(Std2Char(RMessage(1)))                                                // &
               ', Z = '                                                                   // &
               Trim(Std2Char(RMessage(2)))                                                // &
               ', '                                                                       // &
               Char(1:I2)                                                                 // &
               ' = '                                                                      // &
               Trim(Std2Char(RMessage(3)))                                                // &
               ', Z0 = '                                                                  // &
               Trim(Std2Char(RMessage(4)))                                                // &
               ', ABSF = '                                                                // &
               Trim(Std2Char(RMessage(5)))                                                // &
               '.'

      Else If (N == 18) Then

        Line = 'WARNING: Hour number '                                  // &
               Trim(Int2Char(IHour))                                    // &
               '. The surface sensible heat flux from the met file is ' // &
               Trim(Std2Char(RMessage(1)))                              // &
               '. This is less than (and has been replaced by) '        // &
               'the minimum allowed value for the conditions, namely '  // &
               Trim(Std2Char(RMessage(2)))                              // &
               '.'

      Else If (N == 19) Then

        Line = 'WARNING: Hour number '            // &
               Trim(Int2Char(IHour))              // &
               '. The friction Rossby number is ' // &
               Trim(Std2Char(RMessage(1)))        // &
               ' and is less than '               // &
               Trim(Std2Char(RMessage(2)))        // &
               '.'

      Else If (N == 20) Then

        Line = 'WARNING: Hour number '                                                          // &
               Trim(Int2Char(IHour))                                                            // &
               '. There is insufficient information to estimate the depth to which the '        // &
               'growing daytime boundary layer has grown. As a result, an estimate of '         // &
               'boundary layer depth based on the assumption of neutral conditions is adopted.'

      Else If (N == 21) Then

        Line = 'ERROR: The met file contains detailed frequency data ' // &
               'but not basic frequency data - program terminated.'

      Else If (N == 22) Then

        Line = 'ERROR: The time zone used in defining the detailed frequency met data (UTC + ' // &
               Trim(Std2Char(R, FormatString = 'F8.2'))                                        // &
               ') is inappropriate - program terminated.'


      Else If (N == 23) Then

        Line = 'WARNING: Hour number '                                               // &
               Trim(Int2Char(IHour))                                                 // &
               '. The value of U obtained from the met file is 0.0 indicating calm ' // &
               'conditions - it is impossible to use this data.'

      Else If (N == 24) Then ! There's no message 24

      Else If (N == 25) Then ! There's no message 25

      Else If (N == 26) Then

        Line = 'Number of hours/lines of met data: ' // &
               Trim(Int2Char(IHour - 1))             // &
               '.'

      Else If (N == 27) Then

        Line = 'ERROR: Error opening met file. Please check the met ' // &
               'file name is specified correctly in the input.'

      Else If (N == 28) Then

        Line = 'ERROR: Error reading met file. Please check the met file is in the correct format.'

      Else If (N == 29) Then

        Line = 'WARNING: The met file contains the variable R. If you are doing a calculation '    // &
               'where "met site not representative of source site" is selected, you should note '  // &
               'that separate values of R can be specified in the met file for the met site and '  // &
               'for the dispersion area. The variables R and ALBEDO(M) will only affect the met '  // &
               'site value of albedo. The dispersion area value can be set through the variable '  // &
               'ALBEDO(D). If no value is given the value from the menu is used and, if no value ' // &
               'is given there, a default is assumed. Renaming R as ALBEDO(M) will suppress this ' // &
               'warning message.'

      Else If (N == 30) Then

        Line = 'WARNING: The met file contains the variable ALPHA. If you are doing a calculation '   // &
               'where "met site not representative of source site" is selected, you should note '     // &
               'that separate values of ALPHA can be specified in the met file for the met site and ' // &
               'for the dispersion area. The variables ALPHA and ALPHA(M) will only affect the met '  // &
               'site value of alpha. The dispersion area value can be set through the variable '      // &
               'ALPHA(D). If no value is given the value from the menu is used and, if no value '     // &
               'is given there, a default is assumed. Renaming ALPHA as ALPHA(M) will suppress this ' // &
               'warning message.'

      Else If (N == 31) Then

        Line = 'WARNING: The '                                                // &
               Trim(Int2Char(I1)) // Int2Ordinal(I1)                          // &
               ' and '                                                        // &
               Trim(Int2Char(I2)) // Int2Ordinal(I2)                          // &
               ' variables in the met file represent the same quantity'       // &
               Trim(Char)                                                     // &
               '. The last value of the quantity which occurs, will be used.'

      End If

      ! Add warning when further messages will be suppressed, except for messages
      ! 1-9, 15, 21, 22, 26-30 (these are only relevant once, so there's no need to
      ! warn of further suppression of messages).
      If (                                        &
        MessageCount(N) == MessageLimits(N) .and. &
        N > 9                               .and. &
        N /= 15                             .and. &
        N /= 21                             .and. &
        N /= 22                             .and. &
        N /= 26                             .and. &
        N /= 27                             .and. &
        N /= 28                             .and. &
        N /= 29                             .and. &
        N /= 30                                   &
      ) Then
        Line = Trim(Line) // ' Further messages of this type will be suppressed.'
      End If

      Call Message(Line)

    End If

  End If

End Subroutine

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

      SUBROUTINE Check(                                &
         RLat, TSample, Z, Z0, Z0D, PCorr,             &
         MetLimitsFile,                                &
         LMetSiteRepresentative, LFirstTime, LMessage, &
         LFatalError)
!***********************************************************************
! This subroutine checks that the data transferred into the met input  *
! module from other parts of ADMS are sensible and not missing when    *
! needed.                                                              *
!                                                                      *
!                                                                      *
! Input variables:                                                     *
!                                                                      *
! Local variables:                                                     *
! RLOWERLIMITS ) limits within which various variables from other      *
! UPPERLIMITS  )    parts of ADMS must lie in order to be counted as   *
!                   sensible.                                          *
! RMESSAGE       array required in the argument list of the subroutine *
!                   SSMessage. It has no significance in this          *
!                   subroutine.                                        *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Items in argument list.
      REAL*4        Z,PCORR,Z0,Z0D,RLat,TSample
      CHARACTER*32  MetLimitsFile
      LOGICAL     LMetSiteRepresentative
      LOGICAL       LFIRSTTIME
      LOGICAL       LMESSAGE
      LOGICAL LFatalError
!..........Local variables.
      REAL*4 RLOWERLIMITS(7),UPPERLIMITS(7),RMESSAGE(1)

!***********************************************************************
!          If LFIRSTTIME = .TRUE., read limits within which various    *
!          variables from other parts of ADMS must lie in order to be  *
!          counted as sensible.                                        *
!***********************************************************************

      LFatalError = .FALSE.
      IF (LFIRSTTIME) THEN
        CALL READLIMITS(RLOWERLIMITS,UPPERLIMITS,LFatalError, &
           MetLimitsFile)
        IF (LFatalError) RETURN
      ELSE
        LFatalError = .FALSE.
      ENDIF

!***********************************************************************
!          Case 1: LFIRSTTIME = .TRUE..                                *
!***********************************************************************

      IF (LFIRSTTIME) THEN
        IF (Z.NE.-999.0.AND.Z.NE.0.0.AND.Z.NE.1000.0.AND.       &
            (Z.GT.UPPERLIMITS(1).OR.Z.LT.RLOWERLIMITS(1))) THEN
          LFatalError = .TRUE.
          CALL SSMessage(2,Z,RMESSAGE,0,5,'ZWIND')
          RETURN
        ENDIF
        IF (PCORR.NE.-999.0.AND.                                        &
            (PCORR.GT.UPPERLIMITS(2).OR.PCORR.LT.RLOWERLIMITS(2))) THEN
          LFatalError = .TRUE.
          CALL SSMessage(2,PCORR,RMESSAGE,0,5,'PCORR')
          RETURN
        ENDIF
        IF (PCORR.NE.-999.0.AND.PCORR.NE.1.0.AND. &
            LMetSiteRepresentative) THEN
          LFatalError = .TRUE.
          CALL SSMessage(7,0.0,RMESSAGE,0,0,' ')
          RETURN
        ENDIF
        IF (Z0.NE.-999.0.AND.                                     &
            (Z0.GT.UPPERLIMITS(3).OR.Z0.LT.RLOWERLIMITS(3))) THEN
          LFatalError = .TRUE.
          CALL SSMessage(2,Z0,RMESSAGE,0,2,'Z0')
          RETURN
        ENDIF
        IF (Z0D.NE.-999.0.AND.                                      &
            (Z0D.GT.UPPERLIMITS(3).OR.Z0D.LT.RLOWERLIMITS(3))) THEN
          LFatalError = .TRUE.
          CALL SSMessage(2,Z0D,RMESSAGE,0,3,'Z0D')
          RETURN
        ENDIF
        IF (RLat.NE.-999.0.AND.                                       &
            (RLat.GT.UPPERLIMITS(5).OR.RLat.LT.RLOWERLIMITS(5))) THEN
          LFatalError = .TRUE.
          CALL SSMessage(2,RLat,RMESSAGE,0,4,'RLat')
          RETURN
        ENDIF
        IF (RLat.EQ.-999.0) THEN
          LFatalError = .TRUE.
          CALL SSMessage(3,0.0,RMESSAGE,0,4,'RLat')
          RETURN
        ENDIF
        IF (LMESSAGE) THEN
          LFatalError = .TRUE.
          CALL SSMessage(8,0.0,RMESSAGE,0,0,' ')
          RETURN
        ENDIF
      ENDIF

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE READLIMITS(RLOWERLIMITS,UPPERLIMITS,LFatalError, &
                            MetLimitsFile)
!***********************************************************************
! This subroutine reads the met-limits file for use in the             *
! subroutine CHECK.                                                    *
!                                                                      *
! Input Common Blocks: METIN.                                          *
! Output Common Blocks: METOUT.                                        *
!                                                                      *
! Output variables:                                                    *
! RLOWERLIMITS ) limits within which various variables from other      *
! UPPERLIMITS  )    parts of ADMS must lie in order to be counted as   *
!                   sensible.                                          *
!                                                                      *
! Local variables:                                                     *
! RMESSAGE array required in the argument list of the subroutine       *
!             SSMessage. It has no significance in this subroutine.    *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Items in argument list.
      REAL*4 RLOWERLIMITS(7),UPPERLIMITS(7)
      CHARACTER*32  MetLimitsFile
      LOGICAL LFatalError
!..........Local variables.
!      REAL*4 RMESSAGE(1)

!==========Avoid compiler warning.
      IF (MetLimitsFile .EQ. ' ') RLowerLimits(1) = 1.0

      RLOWERLIMITS(1) =    1.0
      RLOWERLIMITS(2) =    0.2
      RLOWERLIMITS(3) =    1.0E-10
      RLOWERLIMITS(4) =    0.0
      RLOWERLIMITS(5) =  -90.0
      RLOWERLIMITS(6) =    0.0
      RLOWERLIMITS(7) = -360.0
      UPPERLIMITS(1)  = 1000.0
      UPPERLIMITS(2)  =    5.0
      UPPERLIMITS(3)  =  100.0
      UPPERLIMITS(4)  =  360.0
      UPPERLIMITS(5)  =   90.0
      UPPERLIMITS(6)  =  100.0
      UPPERLIMITS(7)  =  720.0

      LFatalError = .FALSE.

!..........Read RLOWERLIMITS and UPPERLIMITS (initialise before
!          list-directed read).
!      DO I = 1,7
!        RLowerLimits = 1.0
!        RUpperLimits = 0.0
!      ENDDO
!      OPEN (UNIT=11,FILE=MetLimitsFile,STATUS='OLD', Action = 'Read', ERR=100) ! $$ Use OpenFile
!      READ (11,*,ERR=200,END=200) RLOWERLIMITS
!      READ (11,*,ERR=200,END=200) UPPERLIMITS
!      CLOSE (UNIT=11)

!..........Normal return.
      RETURN

!..........Jump points for read errors.
!  100 CONTINUE
!      LFatalError = .TRUE.
!      CALL SSMessage(1,0.0,RMESSAGE,0,23,'opening met-limits file')
!      CLOSE (UNIT=11)
!      RETURN

!  200 CONTINUE
!      LFatalError = .TRUE.
!      CALL SSMessage(1,0.0,RMESSAGE,0,23,'reading met-limits file')
!      CLOSE (UNIT=11)
!      RETURN

      END Subroutine

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

      SUBROUTINE ReadMetSetUp(MetFile, MetFr, State, LFatalError)
!-----------------------------------------------------------------------
! This subroutine reads the met file headers and initialises MetFr and
! State.
!
! Input variables:
! MetFile: Name of met file.
!
! Output variables:
! MetFr: Met frequency data.
! State: The state of the met file.
! LFatalError: Indicates that a fatal error has occurred.
!
! Local variables:
! I, J: Loop counters.
! IOCheck1, IOCheck2, IOCheck3, IOCheck4, IOCheck5: Variables reporting
!    status of read statements.
! RMessage: Array required in the argument list of the subroutine SSMessage.
!    It has no significance in this subroutine.
! FrTimeZone: Time zone associated with detailed frequency data.
! LTHour, LTimeZone: Indicate the presence of THour and Time Zone in the
!    met file.
! Blank: A blank character string.
! DataLine: A charater string read from the met file.
! KeyWordData, KeyWordVARIABLES: Keywords used in the met file format.
! VariableNames: Short names. These agree with names used for messages
!    produced within the module.
! VariableNamesAlt1: Alternative long names.
! VariableNamesAlt2: Alternative names needed for historical reasons.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      CHARACTER(MaxFileNameLength) MetFile
!..........Output variables.
      Type(MetFr_) MetFr
      Type(MetFileState_) State
      LOGICAL LFatalError
!..........Local variables.
      INTEGER*4 I, J, IOCheck1, IOCheck2, IOCheck3, IOCheck4, IOCheck5
      REAL*4 RMessage(1), FrTimeZone
      LOGICAL LTHour, LTimeZone
      CHARACTER*80 Blank, DataLine, KeyWordData, KeyWordVARIABLES, &
                   VariableNames(NVariableNames + 1),              &
                   VariableNamesAlt1(NVariableNames + 1),          &
                   VariableNamesAlt2(NVariableNames + 1)

!==========Initialise Blank, KeyWordDATA, KeyWordVARIABLES and VariableNames.
!          Note the order of these names should not be changed without
!          understanding all the ReadMet routines and the message
!          routine. In particular note the hard wired numbers 25, 27, 31
!          and 34 for R, Alpha, THour and Time Zone in this routine.

      DATA Blank,KeyWordDATA,KeyWordVARIABLES/                           &
         '                                                            ', &
         'DATA:                                                       ', &
         'VARIABLES:                                                  '/

      DATA VariableNames/                                                &
         'WIND SPEED                                                  ', &
         'UG/USTAR                                                    ', &
         'WIND DIRN                                                   ', &
         'DIRN CHANGE                                                 ', &
         'HEAT FLUX                                                   ', &
         '1/LMO                                                       ', &
         'BL DEPTH                                                    ', &
         'CLOUD                                                       ', &
         'SOLAR RAD                                                   ', &
         'TEMPERATURE                                                 ', &
         'N ABOVE BL                                                  ', &
         'DELTA THETA                                                 ', &
         'PRECIP                                                      ', &
         'SEA TEMP                                                    ', &
         'DELTA T                                                     ', &
         'SIGMA THETA                                                 ', &
         'S HUMIDITY                                                  ', &
         'R HUMIDITY                                                  ', &
         'RH ABOVE BL                                                 ', &
         'DRH/DZ                                                      ', &
         'LAT HT FLUX                                                 ', &
         'WIND HEIGHT                                                 ', &
         'Z0 (M)                                                      ', &
         'Z0 (D)                                                      ', &
         'ALBEDO (M)                                                  ', &
         'ALBEDO (D)                                                  ', &
         'ALPHA (M)                                                   ', &
         'ALPHA (D)                                                   ', &
         '1/LMO < (M)                                                 ', &
         '1/LMO < (D)                                                 ', &
         'HOUR                                                        ', &
         'DAY                                                         ', &
         'YEAR                                                        ', &
         'TIME ZONE                                                   ', &
         'FREQUENCY                                                   ', &
         'FREQUENCY FOR MONTHS xx TO xx, HOURS xxxxx TO xxxxx (GMT + xxxxx)'/
      DATA VariableNamesAlt1/                                            &
         'WIND SPEED                                                  ', &
         'GEOSTROPHIC WIND SPEED/FRICTION VELOCITY                    ', &
         'WIND DIRECTION (DEGREES)                                    ', &
         'GEOSTROPHIC MINUS SURFACE WIND DIRECTION (DEGREES)          ', &
         'SENSIBLE HEAT FLUX                                          ', &
         '1/MONIN-OBUKHOV LENGTH                                      ', &
         'BOUNDARY LAYER DEPTH                                        ', &
         'CLOUD AMOUNT (OKTAS)                                        ', &
         'INCOMING SOLAR RADIATION                                    ', &
         'TEMPERATURE (C)                                             ', &
         'BUOYANCY FREQUENCY ABOVE BOUNDARY LAYER                     ', &
         'TEMPERATURE JUMP ACROSS BOUNDARY LAYER TOP                  ', &
         'PRECIPITATION RATE (MM/HOUR)                                ', &
         'SEA SURFACE TEMPERATURE (C)                                 ', &
         'TEMPERATURE OVER LAND MINUS SEA SURFACE TEMPERATURE         ', &
         'SIGMA THETA (DEGREES)                                       ', &
         'SPECIFIC HUMIDITY                                           ', &
         'RELATIVE HUMIDITY (PERCENT)                                 ', &
         'RELATIVE HUMIDITY ABOVE BOUNDARY LAYER (PERCENT)            ', &
         'D(RELATIVE HUMIDITY)/DZ ABOVE BOUNDARY LAYER (PERCENT/M)    ', &
         'LATENT HEAT FLUX                                            ', &
         'WIND MEASUREMENT HEIGHT                                     ', &
         'ROUGHNESS LENGTH (MET SITE)                                 ', &
         'ROUGHNESS LENGTH (DISPERSION AREA)                          ', &
         'ALBEDO (MET SITE)                                           ', &
         'ALBEDO (DISPERSION AREA)                                    ', &
         'MODIFIED PRIESTLEY-TAYLOR PARAMETER (MET SITE)              ', &
         'MODIFIED PRIESTLEY-TAYLOR PARAMETER (DISPERSION AREA)       ', &
         'MAX 1/MONIN-OBUKHOV LENGTH (MET SITE)                       ', &
         'MAX 1/MONIN-OBUKHOV LENGTH (DISPERSION AREA)                ', &
         'HOUR                                                        ', &
         'DAY                                                         ', &
         'YEAR                                                        ', &
         'TIME ZONE                                                   ', &
         'FREQUENCY                                                   ', &
         'FREQUENCY FOR MONTHS xx TO xx, HOURS xxxxx TO xxxxx (GMT + xxxxx)'/
      DATA VariableNamesAlt2/                                            &
         'U                                                           ', &
         'UGSTAR                                                      ', &
         'PHI                                                         ', &
         'DELTAPHI                                                    ', &
         'FTHETA0                                                     ', &
         'RECIPLMO                                                    ', &
         'H                                                           ', &
         'CL                                                          ', &
         'SOLAR RAD                                                   ', &
         'T0C                                                         ', &
         'NU                                                          ', &
         'DELTATHETA                                                  ', &
         'P                                                           ', &
         'TSEA                                                        ', &
         'DELTAT                                                      ', &
         'SIGMATHETA                                                  ', &
         'S HUMIDITY                                                  ', &
         'RHUM                                                        ', &
         'RH ABOVE BL                                                 ', &
         'DRH/DZ                                                      ', &
         'LAT HT FLUX                                                 ', &
         'WIND HEIGHT                                                 ', &
         'Z0 (M)                                                      ', &
         'Z0 (D)                                                      ', &
         'R                                                           ', &
         'ALBEDO (D)                                                  ', &
         'ALPHA                                                       ', &
         'ALPHA (D)                                                   ', &
         '1/LMO < (M)                                                 ', &
         '1/LMO < (D)                                                 ', &
         'THOUR                                                       ', &
         'TDAY                                                        ', &
         'YEAR                                                        ', &
         'TIME ZONE                                                   ', &
         'FR                                                          ', &
         'MONTHS xx TO xx, HOURS xx TO xx                             '/

!==========Initialise LFatalError.

      LFatalError = .FALSE.

!==========Initialise State.

      DO I = 1, MaxVariables
        State%Codes(I)   = 0
        State%FrCodes(I) = 0
      ENDDO
      State%NAvailable = 0
      State%NRec       = 0
      State%NVar       = 0
      State%MetFile    = MetFile

!==========Initialise MetFr.

      MetFr%FrPresent   = .FALSE.
      MetFr%NDetailedFr = 0

!==========Initialise LTimeZone and LTHour.

      LTimeZone = .FALSE.
      LTHour    = .FALSE.

!==========Open met file.

      OPEN (UNIT = 10, FILE = State%MetFile, STATUS = 'OLD', Action = 'Read', ERR = 100) ! $$ Use OpenFile

!==========Read comments.

   10 CONTINUE
      READ (10, '(A60)', ERR = 200, END = 200) DataLine
      State%NRec = State%NRec + 1
      DO J = 1, 60
        IF (LLE('a',DataLine(J:J)) .AND. LLE(DataLine(J:J),'z'))  &
           DataLine(J:J) = CHAR(ICHAR(DataLine(J:J)) + ICHAR('A') &
           - ICHAR('a'))
      ENDDO
      IF (DataLine.NE.KeyWordVARIABLES) GOTO 10

!==========Read number of variables.

      READ (10, *, ERR = 200, END = 200) State%NVar
      State%NRec = State%NRec + 1

!==========Read list of variables and calculate codes and hour and month
!          frequency limits.

      DO I = 1, State%NVar
        READ (10, '(A60)', ERR = 200, END = 200) DataLine
        State%NRec = State%NRec + 1
        IF (I.LE.MaxVariables) THEN
          DO J = 1, 60
            IF (LLE('a',DataLine(J:J)) .AND. LLE(DataLine(J:J),'z'))  &
               DataLine(J:J) = CHAR(ICHAR(DataLine(J:J)) + ICHAR('A') &
               - ICHAR('a'))
          ENDDO
          DO J = 1, NVariableNames
            IF (DataLine.EQ.VariableNames(J) .OR.      &
                DataLine.EQ.VariableNamesAlt1(J) .OR.  &
                DataLine.EQ.VariableNamesAlt2(J)) THEN
              State%Codes(I) = J
              IF (J.EQ.NVariableNames) MetFr%FrPresent = .TRUE.
            ENDIF
          ENDDO
          IF (DataLine( 1:21).EQ.                             &
              VariableNames(NVariableNames + 1)( 1:21) .AND.  &
              DataLine(24:27).EQ.                             &
              VariableNames(NVariableNames + 1)(24:27) .AND.  &
              DataLine(30:37).EQ.                             &
              VariableNames(NVariableNames + 1)(30:37) .AND.  &
              DataLine(43:46).EQ.                             &
              VariableNames(NVariableNames + 1)(43:46) .AND.  &
              ((                                              &
              DataLine(52:59).EQ.                             &
              VariableNames(NVariableNames + 1)(52:59) .AND.  &
              DataLine(65:80).EQ.                             &
              VariableNames(NVariableNames + 1)(65:80)        &
              ) .OR. DataLine(52:80).EQ.Blank(52:80) ) ) THEN
            READ (DataLine(22:23), *, IOSTAT = IOCheck1)     &
               MetFr%MonthLowerLimits(MetFr%NDetailedFr + 1)
            READ (DataLine(28:29), *, IOSTAT = IOCheck2)     &
               MetFr%MonthUpperLimits(MetFr%NDetailedFr + 1)
            READ (DataLine(38:42), *, IOSTAT = IOCheck3)    &
               MetFr%HourLowerLimits(MetFr%NDetailedFr + 1)
            READ (DataLine(47:51), *, IOSTAT = IOCheck4)    &
               MetFr%HourUpperLimits(MetFr%NDetailedFr + 1)
            IF (DataLine(52:80).EQ.Blank(52:80)) THEN
              FrTimeZone = -999.0
              IOCheck5 = 0
            ELSE
              READ (DataLine(60:64), *, IOSTAT = IOCheck5) FrTimeZone
            ENDIF
            IF (IOCheck1.EQ.0 .AND. IOCheck2.EQ.0 .AND. &
                IOCheck3.EQ.0 .AND. IOCheck4.EQ.0 .AND. &
                IOCheck5.EQ.0) THEN
              MetFr%NDetailedFr = MetFr%NDetailedFr + 1
              State%Codes(I) = -1
              State%FrCodes(I) = MetFr%NDetailedFr
              IF (MetFr%NDetailedFr .EQ. 1) THEN
                MetFr%TimeZone = FrTimeZone
              ELSEIF (MetFr%TimeZone .NE. FrTimeZone) THEN
                LFatalError = .TRUE.
                CALL SSMessage(4, 0.0, RMessage, 0, 0, ' ')
                RETURN
              ENDIF
            ENDIF
          ENDIF
          IF (DataLine( 1: 7).EQ.                                &
              VariableNamesAlt2(NVariableNames + 1)( 1: 7) .AND. &
              DataLine(10:13).EQ.                                &
              VariableNamesAlt2(NVariableNames + 1)(10:13) .AND. &
              DataLine(16:23).EQ.                                &
              VariableNamesAlt2(NVariableNames + 1)(16:23) .AND. &
              DataLine(26:29).EQ.                                &
              VariableNamesAlt2(NVariableNames + 1)(26:29) .AND. &
              DataLine(32:60).EQ.                                &
              VariableNamesAlt2(NVariableNames + 1)(32:60)) THEN
            READ (DataLine( 8: 9), *, IOSTAT = IOCheck1)     &
               MetFr%MonthLowerLimits(MetFr%NDetailedFr + 1)
            READ (DataLine(14:15), *, IOSTAT = IOCheck2)     &
               MetFr%MonthUpperLimits(MetFr%NDetailedFr + 1)
            READ (DataLine(24:25), *, IOSTAT = IOCheck3)    &
               MetFr%HourLowerLimits(MetFr%NDetailedFr + 1)
            READ (DataLine(30:31), *, IOSTAT = IOCheck4)    &
               MetFr%HourUpperLimits(MetFr%NDetailedFr + 1)
            MetFr%HourLowerLimits(MetFr%NDetailedFr + 1) =        &
               MetFr%HourLowerLimits(MetFr%NDetailedFr + 1) - 1.0
            FrTimeZone = 0.0
            IF (IOCheck1.EQ.0 .AND. IOCheck2.EQ.0 .AND. &
                IOCheck3.EQ.0 .AND. IOCheck4.EQ.0) THEN
              MetFr%NDetailedFr = MetFr%NDetailedFr + 1
              State%Codes(I) = -1
              State%FrCodes(I) = MetFr%NDetailedFr
              IF (MetFr%NDetailedFr .EQ. 1) THEN
                MetFr%TimeZone = FrTimeZone
              ELSEIF (MetFr%TimeZone .NE. FrTimeZone) THEN
                LFatalError = .TRUE.
                CALL SSMessage(4, 0.0, RMessage, 0, 0, ' ')
                RETURN
              ENDIF
            ENDIF
          ENDIF
          If (State%Codes(I) == 0 .and. .not.(DataLine .CIEq. 'STATION DCNN')) Then
            CALL SSMessage(14, 0.0, RMessage, 0, 60, DataLine)
          End If
          IF (DataLine.EQ.VariableNamesAlt2(25)) CALL SSMessage(29, &
             0.0, RMessage, 0, 0, ' ')
          IF (DataLine.EQ.VariableNamesAlt2(27)) CALL SSMessage(30, &
             0.0, RMessage, 0, 0, ' ')
          IF (DataLine.EQ.VariableNamesAlt2(31)) LTHour = .TRUE.
          IF (DataLine.EQ.VariableNames(34) .OR.                   &
              DataLine.EQ.VariableNamesAlt1(34) .OR.               &
              DataLine.EQ.VariableNamesAlt2(34)) LTimeZone = .TRUE.
        ENDIF
      ENDDO

!==========Error message for THour and Time Zone both present in the met file.

      IF (LTHour .AND. LTimeZone) THEN
        LFatalError = .TRUE.
        CALL SSMessage(9, 0.0, RMessage, 0, 0, ' ')
        RETURN
      ENDIF

!==========Message to warn that NVar > MaxVariables.

      IF (State%NVar.GT.MaxVariables) THEN
        CALL SSMessage(15, 0.0, RMessage, State%NVar, MaxVariables, ' ')
        State%NVar = MaxVariables
      ENDIF

!==========Message to warn of duplicate variable names.

      DO I = 1, State%NVar - 1
        IF (State%Codes(I).GT.0) THEN
          DO J = I + 1, State%NVar
            IF (State%Codes(I).EQ.State%Codes(J)) CALL SSMessage(31, 0.0, &
               RMessage, I, J, VariableNames(State%Codes(I)))
          ENDDO
        ENDIF
      ENDDO

!==========Read comments.

   30 CONTINUE
      READ (10, '(A60)', ERR = 200, END = 200) DataLine
      State%NRec = State%NRec + 1
      DO J = 1,60
        IF (LLE('a',DataLine(J:J)) .AND. LLE(DataLine(J:J),'z'))  &
           DataLine(J:J) = CHAR(ICHAR(DataLine(J:J)) + ICHAR('A') &
           - ICHAR('a'))
      ENDDO
      IF (DataLine.NE.KeyWordDATA) GOTO 30

!..........Normal return.
      CLOSE (UNIT = 10)
      RETURN

!..........Jump points for errors.
  100 CONTINUE
      LFatalError = .TRUE.
      CALL SSMessage(27, 0.0, RMessage, 0, 16, ' ')
      CLOSE (UNIT = 10)
      RETURN

  200 CONTINUE
      LFatalError = .TRUE.
      CALL SSMessage(28, 0.0, RMessage, 0, 16, ' ')
      CLOSE (UNIT = 10)
      RETURN

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE ReadMet(MetFr, State, Met, LNoMoreData, LFatalError)
!-----------------------------------------------------------------------
! This subroutine controls the reading of the met data from the met file.
!
! Input/output variables:
! MetFr: Met frequency data.
! State: The state of the met file.
!
! Output variables:
! Met: Met data.
! LNoMoreData: Indicates that there are no more data in the met file.
! LFatalError: Indicates that a fatal error has occurred.
!
! Local variables:
! I: Loop counter.
! MetData: Array containing the met data.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input/output variables.
      Type(MetFr_) MetFr
      Type(MetFileState_) State
!..........Output variables.
      Type(SSMet_) Met
      LOGICAL LNoMoreData,LFatalError
!..........Local variables.
      INTEGER*4 I
      REAL*4 MetData(NVariableNames)

!==========Initialise LFatalError.

      LFatalError = .FALSE.

!==========If no data available, try to read some more.

      IF (State%NAvailable.EQ.0) THEN
        CALL ReadData(State, LFatalError)
        IF (LFatalError) RETURN
      ENDIF

!==========If data available, return next hour's data. Otherwise set
!          LNoMoreData = .TRUE..

      IF (State%NAvailable.EQ.0) THEN
        LNoMoreData = .TRUE.
      ELSE
        LNoMoreData = .FALSE.
        State%NAvailable = State%NAvailable - 1
! Vectorisation of following loop gives wrong results with certain complier options (see
! LinuxIntelRelease_EMARC - at least its called that in v4.1). Unclear why. $$
# ifdef IntelLinCompiler
!DEC$ NOVECTOR
        DO I = 1, NVariableNames
          MetData(I) = State%SortedData(I,State%NRead - State%NAvailable)
        ENDDO
# else
        DO I = 1, NVariableNames
          MetData(I) = State%SortedData(I,State%NRead - State%NAvailable)
        ENDDO
# endif
        DO I = 1, MaxVariables
          MetFr%DetailedFr(I) = State%DetailedFr(I, State%NRead - State%NAvailable)
        ENDDO
        Met%U             = MetData(1)
        Met%UGStar        = MetData(2)
        Met%Phi           = MetData(3)
        Met%DeltaPhi      = MetData(4)
        Met%FTheta0       = MetData(5)
        Met%RecipLMO      = MetData(6)
        Met%H             = MetData(7)
        Met%Cl            = MetData(8)
        Met%K             = MetData(9)
        Met%T0C           = MetData(10)
        Met%NU            = MetData(11)
        Met%DeltaTheta    = MetData(12)
        Met%P             = MetData(13)
        Met%TSea          = MetData(14)
        Met%DeltaT        = MetData(15)
        Met%SigmaThetaDeg = MetData(16)
        Met%Q0            = MetData(17)
        Met%RH0           = MetData(18)
        Met%RHU           = MetData(19)
        Met%DRHDZU        = MetData(20)
        Met%LambdaE       = MetData(21)
        Met%Z             = MetData(22)
        Met%Z0            = MetData(23)
        Met%Z0D           = MetData(24)
        Met%R             = MetData(25)
        Met%RD            = MetData(26)
        Met%Alpha         = MetData(27)
        Met%AlphaD        = MetData(28)
        Met%RecipLMOMax   = MetData(29)
        Met%RecipLMOMaxD  = MetData(30)
        Met%Hour          = MetData(31)
        Met%Day           = MetData(32)
        Met%Year          = MetData(33)
        Met%TimeZone      = MetData(34)
        Met%UStar         = -999.0
        Met%UG            = -999.0
        Met%Phi0          = -999.0
        Met%PhiG          = -999.0
        Met%ThetaStar     = -999.0
        Met%WStar         = -999.0
        Met%S             = -999.0
        Met%T0K           = -999.0
        Met%ZD            = -999.0
        Met%SigmaThetaRad = -999.0
        MetFr%Fr          = MetData(NVariableNames)

      ENDIF

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE ReadData(State, LFatalError)
!-----------------------------------------------------------------------
! This subroutine reads the met data from the met file.
!
! Input/output variables:
! State: The state of the met file.
!
! Output variables:
! LFatalError: Indicates that a fatal error has occurred.
!
! Local variables:
! I, J: Loop counters.
! RMessage: Array required in the argument list of the subroutine Message.
!    It has no significance in this subroutine.
! UnsortedData: Array containing the unsorted met data.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input/output variables.
      Type(MetFileState_) State
!..........Output variables.
      LOGICAL LFatalError
!..........Local variables.
      INTEGER*4 I, J
      REAL*4 RMessage(1), UnsortedData(MaxVariables)

!==========Set read-error flag, open met file and read to current position.

      LFatalError = .False.
      OPEN (UNIT = 10, FILE = State%MetFile, STATUS = 'OLD', Action = 'Read', ERR = 100) ! $$ Use OpenFile
      DO I = 1, State%NRec
        READ (10, *, ERR = 200, END = 200)
      ENDDO

!==========Set DetailedFr and SortedData equal to -999.0.

      DO I = 1, MaxVariables
        DO J = 1, NReadAtATime
          State%DetailedFr(I,J) = -999.0
        ENDDO
      ENDDO
      DO I = 1, NVariableNames
        DO J = 1, NReadAtATime
          State%SortedData(I,J) = -999.0
        ENDDO
      ENDDO

!==========Read and sort met data.

      State%NRead = 0
      DO J = 1, NReadAtATime
        DO I = 1, State%NVar
          UnsortedData(I) = -999.0
        ENDDO
        READ (10, *, END = 1000, ERR = 200) (UnsortedData(I), I = 1, State%NVar)
        State%NRec = State%NRec + 1
        State%NRead = State%NRead + 1
        DO I = 1, State%NVar
          IF (State%Codes(I).GT.0)                                &
             State%SortedData(State%Codes(I),J) = UnsortedData(I)
          IF (State%Codes(I).EQ.-1)                                 &
             State%DetailedFr(State%FrCodes(I),J) = UnsortedData(I)
        ENDDO
      ENDDO
 1000 CONTINUE

      State%NAvailable = State%NRead

!..........Normal return.
      CLOSE (UNIT = 10)
      RETURN

!..........Jump points for read errors.
  100 CONTINUE
      LFatalError = .TRUE.
      CALL SSMessage(27, 0.0, RMessage, 0, 16, ' ')
      CLOSE (UNIT = 10)
      RETURN

  200 CONTINUE
      LFatalError = .TRUE.
      CALL SSMessage(28, 0.0, RMessage, 0, 16, ' ')
      CLOSE (UNIT = 10)
      RETURN

      END Subroutine

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

      SUBROUTINE QualityControlMetSetUp(MetFr, MetLimitsFile, LSeq,  &
                                        MinMet, MaxMet, LFatalError)
!-----------------------------------------------------------------------
! This subroutine checks the value of MetFr and reads the met-limits file
! for use in the subroutine QualityControlMet.
!
! Input variables:
! MetFr: Met frequency data to be checked.
! MetLimitsFile: Name of met-limits file.
! LSeq: Indicates that the data provided in the met file are hourly
!    sequential.
!
! Output variables:
! MinMet, MaxMet: Limits within which the met variables must lie to be
!    regarded as sensible.
! LFatalError: Indicates that a fatal error has occurred.
!
! Local variables:
! I: Loop counter.
! RMessage: Array required in the argument list of the subroutine Message.
!    It has no significance in this subroutine.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      Type(MetFr_) MetFr
      CHARACTER*32 MetLimitsFile
      LOGICAL LSeq
!..........Output variables.
      Type(SSMet_) MinMet, MaxMet
      LOGICAL LFatalError
!..........Local variables.
      INTEGER*4 I
      REAL*4 RMessage(1)

!..........To avoid compiler warning:
      IF (MetLimitsFile .EQ. ' ') MinMet%U = 0.0

!==========Initialise LFatalError.

      LFatalError = .FALSE.

!==========Check frequency variables in met file are OK.

      IF (LSeq) THEN
        IF (MetFr%FrPresent .OR. MetFr%NdetailedFr.NE.0) THEN
          LFatalError = .TRUE.
          CALL SSMessage(5, 0.0, RMessage, 0, 0, ' ')
          RETURN
        ENDIF
        IF (.NOT.MetFr%FrPresent .AND. MetFr%NdetailedFr.NE.0) THEN
          LFatalError = .TRUE.
          CALL SSMessage(21, 0.0, RMessage, 0, 0, ' ')
          RETURN
        ENDIF
      ELSE
        DO I = 1, MetFr%NDetailedFr
          IF (MetFr%HourLowerLimits(I).LT. 0.0 .OR.  &
              MetFr%HourLowerLimits(I).GT.24.0 .OR.  &
              MetFr%HourUpperLimits(I).LT. 0.0 .OR.  &
              MetFr%HourUpperLimits(I).GT.24.0) THEN
            LFatalError = .TRUE.
            CALL SSMessage(6, 0.0, RMessage, 0, 0, ' ')
            RETURN
          ENDIF
          IF (MetFr%MonthLowerLimits(I).LT. 0.0 .OR.  &
              MetFr%MonthLowerLimits(I).GT.12.0 .OR.  &
              MetFr%MonthUpperLimits(I).LT. 0.0 .OR.  &
              MetFr%MonthUpperLimits(I).GT.12.0) THEN
            LFatalError = .TRUE.
            CALL SSMessage(6, 0.0, RMessage, 0, 0, ' ')
            RETURN
          ENDIF
          IF (MetFr%TimeZone.LT.-15.0 .OR. MetFr%TimeZone.GT. 15.0) THEN
            LFatalError = .TRUE.
            CALL SSMessage(22, MetFr%TimeZone, RMessage, 0, 0, ' ')
            RETURN
          ENDIF
        ENDDO
      ENDIF

!==========Initialise MinMet and MaxMet.

      MinMet%U             =       0.0
      MinMet%UGStar        =       5.0
      MinMet%Phi           =       0.0
      MinMet%DeltaPhi      =     -60.0
      MinMet%FTheta0       =    -100.0
      MinMet%RecipLMO      =     -10.0
      MinMet%H             =       0.0
      MinMet%Cl            =       0.0
      MinMet%K             =       0.0
      MinMet%T0C           =    -100.0
      MinMet%NU            =       0.0
      MinMet%DeltaTheta    =       0.0
      MinMet%P             =       0.0
      MinMet%TSea          =     -10.0
      MinMet%DeltaT        =     -40.0
      MinMet%SigmaThetaDeg =       0.0
      MinMet%Q0            =       0.0
      MinMet%RH0           =       0.0
      MinMet%RHU           =       0.0
      MinMet%DRHDZU        =     -10.0
      MinMet%LambdaE       =    -100.0
      MinMet%Z             =       0.0
      MinMet%Z0            =   1.0E-10
      MinMet%Z0D           =   1.0E-10
      MinMet%R             =       0.0
      MinMet%RD            =       0.0
      MinMet%Alpha         =       0.0
      MinMet%AlphaD        =       0.0
      MinMet%RecipLMOMax   =       0.0
      MinMet%RecipLMOMaxD  =       0.0
      MinMet%Hour          =       0.0
      MinMet%Day           =       1.0
      MinMet%Year          =       0.0
      MinMet%TimeZone      =     -15.0

      MaxMet%U             =     100.0
      MaxMet%UGStar        =    1000.0
      MaxMet%Phi           =     360.0
      MaxMet%DeltaPhi      =      60.0
      MaxMet%FTheta0       =    1000.0
      MaxMet%RecipLMO      =      10.0
      MaxMet%H             =    4000.0
      MaxMet%Cl            =       8.0
      MaxMet%K             =    1500.0
      MaxMet%T0C           =      60.0
      MaxMet%NU            =       0.1
      MaxMet%DeltaTheta    =      25.0
      MaxMet%P             =     100.0
      MaxMet%TSea          =      40.0
      MaxMet%DeltaT        =      40.0
      MaxMet%SigmaThetaDeg =      90.0
      MaxMet%Q0            =       0.1
      MaxMet%RH0           =     100.0
      MaxMet%RHU           =     100.0
      MaxMet%DRHDZU        =      10.0
      MaxMet%LambdaE       =    1000.0
      MaxMet%Z             =    1000.0
      MaxMet%Z0            =     100.0
      MaxMet%Z0D           =     100.0
      MaxMet%R             =       1.0
      MaxMet%RD            =       1.0
      MaxMet%Alpha         =       3.0
      MaxMet%AlphaD        =       3.0
      MaxMet%RecipLMOMax   =    1.0E20
      MaxMet%RecipLMOMaxD  =    1.0E20
      MaxMet%Hour          =      24.0
      MaxMet%Day           =     366.0
      MaxMet%Year          =    2500.0
      MaxMet%TimeZone      =      15.0

! If restoring the following code, initialise MinMet and MaxMet before
! reading by replacing above values with zeros.
!      OPEN (UNIT = 11, FILE = MetLimitsFile, STATUS = 'OLD', Action = 'Read', ERR = 100) ! $$ Use OpenFile
!      READ (11, *, ERR = 200, END = 200)
!      READ (11, *, ERR = 200, END = 200)
!      READ (11, *, ERR = 200, END = 200)
!         MinMet.U,             MinMet.UGStar,      MinMet.Phi,
!         MinMet.DeltaPhi,      MinMet.FTheta0,     MinMet.RecipLMO,
!         MinMet.H,             MinMet.Cl,          MinMet.K,
!         MinMet.T0C,           MinMet.NU,          MinMet.DeltaTheta,
!         MinMet.P,             MinMet.TSea,        MinMet.DeltaT,
!         MinMet.SigmaThetaDeg, MinMet.Q0,          MinMet.RH0,
!         MinMet.RHU,           MinMet.DRHDZU,      MinMet.LambdaE,
!         MinMet.Z,             MinMet.Z0,          MinMet.Z0D,
!         MinMet.R,             MinMet.RD,          MinMet.Alpha,
!         MinMet.AlphaD,        MinMet.RecipLMOMax, MinMet.RecipLMOMaxD,
!         MinMet.Hour,          MinMet.Day,         MinMet.Year,
!         MinMet.TimeZone
!      READ (11, *, ERR = 200, END = 200)
!         MaxMet.U,             MaxMet.UGStar,      MaxMet.Phi,
!         MaxMet.DeltaPhi,      MaxMet.FTheta0,     MaxMet.RecipLMO,
!         MaxMet.H,             MaxMet.Cl,          MaxMet.K,
!         MaxMet.T0C,           MaxMet.NU,          MaxMet.DeltaTheta,
!         MaxMet.P,             MaxMet.TSea,        MaxMet.DeltaT,
!         MaxMet.SigmaThetaDeg, MaxMet.Q0,          MaxMet.RH0,
!         MaxMet.RHU,           MaxMet.DRHDZU,      MaxMet.LambdaE,
!         MaxMet.Z,             MaxMet.Z0,          MaxMet.Z0D,
!         MaxMet.R,             MaxMet.RD,          MaxMet.Alpha,
!         MaxMet.AlphaD,        MaxMet.RecipLMOMax, MaxMet.RecipLMOMaxD,
!         MaxMet.Hour,          MaxMet.Day,         MaxMet.Year,
!         MaxMet.TimeZone
!      CLOSE (UNIT = 11)

!..........Normal return.
      RETURN

!..........Jump points for read errors.
!  100 CONTINUE
!      LFatalError = .TRUE.
!      CALL SSMessage(1, 0.0, RMessage, 0, 23, 'opening met-limits file')
!      CLOSE (UNIT = 11)
!      RETURN

!  200 CONTINUE
!      LFatalError = .TRUE.
!      CALL SSMessage(1, 0.0, RMessage, 0, 23, 'reading met-limits file')
!      CLOSE (UNIT = 11)
!      RETURN

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE QualityControlMet(MinMet, MaxMet, Met, MetFr, LInadequateData)
!-----------------------------------------------------------------------
! This subroutine checks that the data obtained from the met file are
! sensible.
!
! Input variables:
! MinMet, MaxMet: Limits within which the met variables must lie to be
!    regarded as sensible.
! Met: Met data to be checked.
! MetFr: Met frequency data to be checked.
!
! Output variables:
! LInadequateData: Indicates that the data are not sensible.
!
! Local variables:
! RMessage: Array required in the argument list of the subroutine SSMessage.
!    It has no significance in this subroutine.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      Type(SSMet_) MinMet, MaxMet, Met
      Type(MetFr_) MetFr
!..........Output variables.
      LOGICAL LInadequateData
!..........Local variables.
      INTEGER*4 I
      REAL*4 RMessage(1)

!==========Initialise LInadequateData.

      LInadequateData = .FALSE.

!==========Check Met values are sensible.

      IF (Met%U.NE.-999.0 .AND.    &
         (Met%U.GT.MaxMet%U .OR.   &
          Met%U.LT.MinMet%U)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%U, RMessage, 0, 10, 'WIND SPEED')
        RETURN
      ENDIF

      IF (Met%UGStar.NE.-999.0 .AND.         &
         (Met%UGStar.GT.MaxMet%UGStar .OR.   &
          Met%UGStar.LT.MinMet%UGStar)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%UGStar, RMessage, 0, 8, 'UG/USTAR')
        RETURN
      ENDIF

      IF (Met%Phi.NE.-999.0 .AND.      &
         (Met%Phi.GT.MaxMet%Phi .OR.   &
          Met%Phi.LT.MinMet%Phi)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%Phi, RMessage, 0, 9, 'WIND DIRN')
        RETURN
      ENDIF

      IF (Met%DeltaPhi.NE.-999.0 .AND.           &
         (Met%DeltaPhi.GT.MaxMet%DeltaPhi .OR.   &
          Met%DeltaPhi.LT.MinMet%DeltaPhi)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%DeltaPhi, RMessage, 0, 11, 'DIRN CHANGE')
        RETURN
      ENDIF

      IF (Met%FTheta0.NE.-999.0 .AND.          &
         (Met%FTheta0.GT.MaxMet%FTheta0 .OR.   &
          Met%FTheta0.LT.MinMet%FTheta0)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%FTheta0, RMessage, 0, 9, 'HEAT FLUX')
        RETURN
      ENDIF

      IF (Met%RecipLMO.NE.-999.0 .AND.           &
         (Met%RecipLMO.GT.MaxMet%RecipLMO .OR.   &
          Met%RecipLMO.LT.MinMet%RecipLMO)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%RecipLMO, RMessage, 0, 5, '1/LMO')
        RETURN
      ENDIF

      IF (Met%H.NE.-999.0 .AND.    &
         (Met%H.GT.MaxMet%H .OR.   &
          Met%H.LT.MinMet%H)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%H, RMessage, 0, 8, 'BL DEPTH')
        RETURN
      ENDIF

      IF (Met%Cl.NE.-999.0 .AND.     &
         (Met%Cl.GT.MaxMet%Cl .OR.   &
          Met%Cl.LT.MinMet%Cl)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%Cl, RMessage, 0, 5, 'CLOUD')
        RETURN
      ENDIF

      IF (Met%K.NE.-999.0 .AND.    &
         (Met%K.GT.MaxMet%K .OR.   &
          Met%K.LT.MinMet%K)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%K, RMessage, 0, 9, 'SOLAR RAD')
        RETURN
      ENDIF

      IF (Met%T0C.NE.-999.0 .AND.      &
         (Met%T0C.GT.MaxMet%T0C .OR.   &
          Met%T0C.LT.MinMet%T0C)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%T0C, RMessage, 0, 11, 'TEMPERATURE')
        RETURN
      ENDIF

      IF (Met%NU.NE.-999.0 .AND.     &
         (Met%NU.GT.MaxMet%NU .OR.   &
          Met%NU.LT.MinMet%NU)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%NU, RMessage, 0, 10, 'N ABOVE BL')
        RETURN
      ENDIF

      IF (Met%DeltaTheta.NE.-999.0 .AND.             &
         (Met%DeltaTheta.GT.MaxMet%DeltaTheta .OR.   &
          Met%DeltaTheta.LT.MinMet%DeltaTheta)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%DeltaTheta, RMessage, 0, 11, 'DELTA THETA')
        RETURN
      ENDIF

      IF (Met%P.NE.-999.0 .AND.    &
         (Met%P.GT.MaxMet%P .OR.   &
          Met%P.LT.MinMet%P)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%P, RMessage, 0, 6, 'PRECIP')
        RETURN
      ENDIF

      IF (Met%TSea.NE.-999.0 .AND.       &
         (Met%TSea.GT.MaxMet%TSea .OR.   &
          Met%TSea.LT.MinMet%TSea)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%TSea, RMessage, 0, 8, 'SEA TEMP')
        RETURN
      ENDIF

      IF (Met%DeltaT.NE.-999.0 .AND.         &
         (Met%DeltaT.GT.MaxMet%DeltaT .OR.   &
          Met%DeltaT.LT.MinMet%DeltaT)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%DeltaT, RMessage, 0, 7, 'DELTA T')
        RETURN
      ENDIF

      IF (Met%SigmaThetaDeg.NE.-999.0 .AND.                &
         (Met%SigmaThetaDeg.GT.MaxMet%SigmaThetaDeg .OR.   &
          Met%SigmaThetaDeg.LT.MinMet%SigmaThetaDeg)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%SigmaThetaDeg, RMessage, 0, 11, 'SIGMA THETA')
        RETURN
      ENDIF

      IF (Met%Q0.NE.-999.0 .AND.     &
         (Met%Q0.GT.MaxMet%Q0 .OR.   &
          Met%Q0.LT.MinMet%Q0)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%Q0, RMessage, 0, 10, 'S HUMIDITY')
        RETURN
      ENDIF

      IF (Met%RH0.NE.-999.0 .AND.      &
         (Met%RH0.GT.MaxMet%RH0 .OR.   &
          Met%RH0.LT.MinMet%RH0)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%RH0, RMessage, 0, 10, 'R HUMIDITY')
        RETURN
      ENDIF

      IF (Met%RHU.NE.-999.0 .AND.      &
         (Met%RHU.GT.MaxMet%RHU .OR.   &
          Met%RHU.LT.MinMet%RHU)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%RHU, RMessage, 0, 11, 'RH ABOVE BL')
        RETURN
      ENDIF

      IF (Met%DRHDZU.NE.-999.0 .AND.         &
         (Met%DRHDZU.GT.MaxMet%DRHDZU .OR.   &
          Met%DRHDZU.LT.MinMet%DRHDZU)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%DRHDZU, RMessage, 0, 6, 'DRH/DZ')
        RETURN
      ENDIF

      IF (Met%LambdaE.NE.-999.0 .AND.          &
         (Met%LambdaE.GT.MaxMet%LambdaE .OR.   &
          Met%LambdaE.LT.MinMet%LambdaE)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%LambdaE, RMessage, 0, 11, 'LAT HT FLUX')
        RETURN
      ENDIF

      IF (Met%Z.NE.-999.0 .AND.    &
         (Met%Z.GT.MaxMet%Z .OR.   &
          Met%Z.LT.MinMet%Z)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%Z, RMessage, 0, 11, 'WIND HEIGHT')
        RETURN
      ENDIF

      IF (Met%Z0.NE.-999.0 .AND.     &
         (Met%Z0.GT.MaxMet%Z0 .OR.   &
          Met%Z0.LT.MinMet%Z0)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%Z0, RMessage, 0, 6, 'Z0 (M)')
        RETURN
      ENDIF

      IF (Met%Z0D.NE.-999.0 .AND.      &
         (Met%Z0D.GT.MaxMet%Z0D .OR.   &
          Met%Z0D.LT.MinMet%Z0D)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%Z0D, RMessage, 0, 6, 'Z0 (D)')
        RETURN
      ENDIF

      IF (Met%R.NE.-999.0 .AND.    &
         (Met%R.GT.MaxMet%R .OR.   &
          Met%R.LT.MinMet%R)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%R, RMessage, 0, 10, 'ALBEDO (M)')
        RETURN
      ENDIF

      IF (Met%RD.NE.-999.0 .AND.     &
         (Met%RD.GT.MaxMet%RD .OR.   &
          Met%RD.LT.MinMet%RD)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%RD, RMessage, 0, 10, 'ALBEDO (D)')
        RETURN
      ENDIF

      IF (Met%Alpha.NE.-999.0 .AND.        &
         (Met%Alpha.GT.MaxMet%Alpha .OR.   &
          Met%Alpha.LT.MinMet%Alpha)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%Alpha, RMessage, 0, 9, 'ALPHA (M)')
        RETURN
      ENDIF

      IF (Met%AlphaD.NE.-999.0 .AND.         &
         (Met%AlphaD.GT.MaxMet%AlphaD .OR.   &
          Met%AlphaD.LT.MinMet%AlphaD)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%AlphaD, RMessage, 0, 9, 'ALPHA (D)')
        RETURN
      ENDIF

      IF (Met%RecipLMOMax.NE.-999.0 .AND.              &
         (Met%RecipLMOMax.GT.MaxMet%RecipLMOMax .OR.   &
          Met%RecipLMOMax.LT.MinMet%RecipLMOMax)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%RecipLMOMax, RMessage, 0, 11, '1/LMO < (M)')
        RETURN
      ENDIF

      IF (Met%RecipLMOMaxD.NE.-999.0 .AND.               &
         (Met%RecipLMOMaxD.GT.MaxMet%RecipLMOMaxD .OR.   &
          Met%RecipLMOMaxD.LT.MinMet%RecipLMOMaxD)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%RecipLMOMaxD, RMessage, 0, 11, '1/LMO < (D)')
        RETURN
      ENDIF

      IF (Met%Hour.NE.-999.0 .AND.       &
         (Met%Hour.GT.MaxMet%Hour .OR.   &
          Met%Hour.LT.MinMet%Hour)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%Hour, RMessage, 0, 4, 'HOUR')
        RETURN
      ENDIF

      IF (Met%Day.NE.-999.0 .AND.      &
         (Met%Day.GT.MaxMet%Day .OR.   &
          Met%Day.LT.MinMet%Day)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%Day, RMessage, 0, 3, 'DAY')
        RETURN
      ENDIF

      IF (Met%Year.NE.-999.0 .AND.       &
         (Met%Year.GT.MaxMet%Year .OR.   &
          Met%Year.LT.MinMet%Year)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%Year, RMessage, 0, 4, 'YEAR')
        RETURN
      ENDIF

      IF (Met%TimeZone.NE.-999.0 .AND.           &
         (Met%TimeZone.GT.MaxMet%TimeZone .OR.   &
          Met%TimeZone.LT.MinMet%TimeZone)) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, Met%TimeZone, RMessage, 0, 9, 'TIME ZONE')
        RETURN
      ENDIF

!==========Check MetFr values are sensible.

      IF (MetFr%FrPresent .AND. MetFr%Fr.LT.0.0) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(11, MetFr%Fr, RMessage, 0, 9, 'FREQUENCY')
        RETURN
      ENDIF
      DO I = 1, MetFr%NDetailedFr
        IF (MetFr%DetailedFr(I).LT.0.0) THEN
          LInadequateData = .TRUE.
          CALL SSMessage(10, MetFr%DetailedFr(I), RMessage, I, 0, ' ')
          RETURN
        ENDIF
      ENDDO

      RETURN
      END Subroutine

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

      SUBROUTINE SiteInit(RLong, RLat, Z, Z0, R, Alpha, RecipLMOMax, &
                          LMetSite, Site)
!-----------------------------------------------------------------------
! This subroutine initialises Site.
!
! Input variables:
! RLong: Longitude (degrees, east positive).
! RLat: Latitude (degrees, north positive).
! Z, Z0, R, Alpha, RecipLMOMax: Values of Z, Z0, R, Alpha and RecipLMOMax.
!    Missing data is indicated by -999.0. (If values are missing, then
!    values from the met file will be used instead. If there are no values
!    in the met file either, then, for Alpha and R, default values will be
!    used.)
! LMetSite: Indicates that the met site, as opposed to the dispersion
!    area, is being considered.
!
! Output variables:
! Site: Site characteristics (including some information on the recent met).
!
! Local variables:
! I: Loop counter.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      REAL*4 RLong, RLat, Z, Z0, R, Alpha, RecipLMOMax
      LOGICAL LMetSite
!..........Output variables.
      Type(Site_) Site
!..........Local variables.
      INTEGER*4 I

!==========Initialise Long, Lat and MetSite and calculate AbsF and SignF.

      Site%Long    = RLong
      Site%Lat     = RLat
      Site%MetSite = LMetSite
      CALL CoriolisCalc(Site%Lat, Site%AbsF, Site%SignF)

!==========Initialise Defaults.

      Site%Z           = Z
      Site%Z0          = Z0
      Site%R           = R
      Site%Alpha       = Alpha
      Site%RecipLMOMax = RecipLMOMax

!==========Initialise UG2B, U2BStar, U3B, ZLast and Z0Last.

      Site%UG2B    = -999.0
      Site%U2BStar = -999.0
      Site%U3B     = -999.0
      Site%ZLast   = -999.0
      Site%Z0Last  = -999.0

!==========Initialise UStar, FTheta0, RNU, T0 and S.

      DO I = 1,24
        Site%UStar(I)   = -999.0
        Site%FTheta0(I) = -999.0
        Site%S(I)       = -999.0
        Site%T0K(I)     = -999.0
        Site%NU(I)      = -999.0
      ENDDO

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE SiteIncrement(LSeq, Site)
!-----------------------------------------------------------------------
! This subroutine increments Site to reflect the fact that time has
! advanced by one hour.
!
! Input variables:
! LSeq: Indicates that the data provided in the met file are hourly
!    sequential.
!
! Input/output variables:
! Site: Site characteristics (including some information on the recent met).
!
! Local variables:
! I: Loop counter.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      LOGICAL LSeq
!..........Input/output variables.
      Type(Site_) Site
!..........Local variables.
      INTEGER*4 I

!==========If sequential data, shift the history arrays to reflect the
!          fact that time has advanced by one hour. If not sequential
!          data, set arrays to -999.0.

      IF (LSeq) THEN
        DO I = 1, 23
          Site%UStar(I)   = Site%UStar(I + 1)
          Site%FTheta0(I) = Site%FTheta0(I + 1)
          Site%S(I)       = Site%S(I + 1)
          Site%T0K(I)     = Site%T0K(I + 1)
          Site%NU(I)      = Site%NU(I + 1)
        ENDDO
        Site%UStar(24)   = -999.0
        Site%FTheta0(24) = -999.0
        Site%S(24)       = -999.0
        Site%T0K(24)     = -999.0
        Site%NU(24)      = -999.0
      ELSE
        DO I = 1, 24
          Site%UStar(I)   = -999.0
          Site%FTheta0(I) = -999.0
          Site%S(I)       = -999.0
          Site%T0K(I)     = -999.0
          Site%NU(I)      = -999.0
        ENDDO
      ENDIF

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE ProcessMet(MetIn, TSample, LSeq, Site, MetOut, &
                            LInadequateData, LCalm)
!-----------------------------------------------------------------------
! This subroutine processes the met data.
!
! Parameters:
! Pi: Pi.
! CP: Specific heat capacity of dry air at constant pressure.
! RhoA: Density of air.
! TCToK: Absolute temperature of the zero of the Celsius scale.
! RDef: Default value for R.
! AlphaDef: Default value for Alpha.
! ClDef: Default value for Cl.
! T0CDef: Default value for T0C.
! RNUDef: Default value for NU.
! RHUDef: Default value for RHU.
! DRHDZUDef: Default value for DRHDZU.
!
! Input variables:
! MetIn: Met data to be processed.
! TSample: Sampling time.
! LSeq: Indicates that the data provided in the met file are hourly
!    sequential.
!
! Input/output variables:
! Site: Site characteristics (including some information on the recent met).
!
! Output variables:
! MetOut: Processed met data.
! LInadequateData: Indicates that a particular hours' data from the met
!    file are inadequate or not sensible.
! LCalm: Indicates that U = 0 for a particular hours' data.
!
! Local variables:
! IterationCount: Number of iterations used in subroutine UCalc1, UCalc2 or
!    UCalc3 (used for testing purposes only).
! Z, Z0, R, Alpha, RecipLMOMax: Values of Z, Z0, R, Alpha and RecipLMOMax
!    to be used in this routine.
! Cl: Cloud for use in estimating FTheta0 (for given RK).
! DeltaThetaConv: Temperature jump across the boundary layer top estimated
!    from the model of boundary layer growth if the value of MetOut.H is
!    the value obtained from this model. -999.0 otherwise.
! H: Temporary value of boundary layer depth.
! RK: Incoming solar radiation for use in estimating FTheta0.
! RMessage: Array required in the argument list of the subroutine SSMessage.
!    It has no significance in this subroutine.
! UStarDum, RecipLMODum: Quantities which are returned from subroutines
!    and which we have no interest in in this subroutine.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Parameters.
      REAL*4 Pi, CP, RhoA, TCToK
      PARAMETER (CP = 1012.0, RhoA = 1.225, TCToK = 273.15)
      REAL*4 RDef, AlphaDef, ClDef, T0CDef, RNUDef, RHUDef, DRHDZUDef
      PARAMETER (RDef = 0.23, AlphaDef = 1.0, ClDef = 5.0,     &
                 T0CDef = 15.0, RNUDef = 0.013, RHUDef = 50.0, &
                 DRHDZUDef = 0.0)
!..........Input variables.
      Type(SSMet_) MetIn
      REAL*4 TSample
      LOGICAL LSeq
!..........Input/output variables.
      Type(Site_) Site
!..........Output variables.
      Type(SSMet_) MetOut
      LOGICAL LInadequateData, LCalm
!..........Local variables.
      INTEGER*4 IterationCount
      REAL*4 Z, Z0, R, Alpha, RecipLMOMax, Cl, DeltaThetaConv, H, RK, &
             RMessage(1), UStarDum, RecipLMODum
      LOGICAL LRecipLMOMaxUsed
!..........Pi.
      Pi = 4.0*ATAN(1.0)

!==========Initialise LInadequateData.

      LInadequateData = .FALSE.

!==========Initialise MetOut.

      MetOut = MetIn

!==========Process site variables in met data.

      IF (Site%MetSite) THEN
        IF (Site%Z          .NE.-999.0) MetOut%Z     = Site%Z
        IF (Site%Z0         .NE.-999.0) MetOut%Z0    = Site%Z0
        IF (Site%R          .NE.-999.0) MetOut%R     = Site%R
        IF (Site%Alpha      .NE.-999.0) MetOut%Alpha = Site%Alpha
        IF (Site%RecipLMOMax.NE.-999.0) MetOut%RecipLMOMax = Site%RecipLMOMax
        IF (MetOut%R    .EQ.-999.0) MetOut%R     = RDef
        IF (MetOut%Alpha.EQ.-999.0) MetOut%Alpha = AlphaDef
      ELSE
        IF (Site%Z          .NE.-999.0) MetOut%ZD     = Site%Z
        IF (Site%Z0         .NE.-999.0) MetOut%Z0D    = Site%Z0
        IF (Site%R          .NE.-999.0) MetOut%RD     = Site%R
        IF (Site%Alpha      .NE.-999.0) MetOut%AlphaD = Site%Alpha
        IF (Site%RecipLMOMax.NE.-999.0) MetOut%RecipLMOMaxD = Site%RecipLMOMax
        IF (MetOut%RD    .EQ.-999.0) MetOut%RD     = RDef
        IF (MetOut%AlphaD.EQ.-999.0) MetOut%AlphaD = AlphaDef
      ENDIF

!==========Calculate site variables for use in this routine.

      IF (Site%MetSite) THEN
        Z           = MetOut%Z
        Z0          = MetOut%Z0
        R           = MetOut%R
        Alpha       = MetOut%Alpha
        RecipLMOMax = MetOut%RecipLMOMax
      ELSE
        Z           = MetOut%ZD
        Z0          = MetOut%Z0D
        R           = MetOut%RD
        Alpha       = MetOut%AlphaD
        RecipLMOMax = MetOut%RecipLMOMaxD
      ENDIF

!==========Set elements of MetOut which shouldn't be input through MetIn
!          or which we want to prevent being output to -999.0.

      MetOut%UStar         = -999.0
      MetOut%UG            = -999.0
      MetOut%Phi0          = -999.0
      MetOut%PhiG          = -999.0
      MetOut%WStar         = -999.0
      MetOut%S             = -999.0
      MetOut%T0K           = -999.0
      MetOut%SigmaThetaRad = -999.0
      MetOut%LocalMeanTime = -999.0

      MetOut%U             = -999.0
      MetOut%Phi           = -999.0
      MetOut%TSea          = -999.0
      MetOut%RH0           = -999.0
      MetOut%Z             = -999.0
      MetOut%ZD            = -999.0

!==========Calculate LocalMeanTime.

      IF (MetOut%Hour.NE.-999.0) THEN
        IF (MetOut%TimeZone.NE.-999.0) THEN
          MetOut%LocalMeanTime = MetOut%Hour - MetOut%TimeZone + Site%Long/15.0
        ELSE
          MetOut%LocalMeanTime = MetOut%Hour
        ENDIF
      ENDIF

!==========Calculate DeltaT.

      IF (MetOut%DeltaT.EQ.-999.0) THEN
        IF (MetIn%T0C.NE.-999.0 .AND. MetIn%TSea.NE.-999.0) &
           MetOut%DeltaT = MetIn%T0C - MetIn%TSea
      ENDIF

!==========Estimate NU.

      IF (MetOut%NU.EQ.-999.0) MetOut%NU = RNUDef

!==========Estimate T0.

      IF (MetOut%T0C.EQ.-999.0) MetOut%T0C = T0CDef
      MetOut%T0K = MetOut%T0C + TCToK

!==========Estimate S.

      IF (MetOut%S.EQ.-999.0) THEN
        IF (MetOut%Day.NE.-999.0 .AND. MetOut%LocalMeanTime.NE.-999.0) THEN
          CALL SCalc(Site%Lat, MetOut%Day, MetOut%LocalMeanTime, MetOut%S)
        ENDIF
      ENDIF

!==========Estimate K.

      IF (MetOut%K.EQ.-999.0 .AND. MetOut%Cl.NE.-999.0 .AND. &
          MetOut%S.NE.-999.0) THEN
        CALL KCalc(MetOut%Cl, MetOut%S, RK)
        MetOut%K = MAX(RK, 0.0)
      ELSE
        RK = MetOut%K
      ENDIF

!==========Estimate Cl.

      IF (MetOut%Cl.NE.-999.0) THEN
        Cl = MetOut%Cl
      ELSE
        Cl = ClDef
      ENDIF

!==========Check Z, Z0 and U available.

      IF (Z.EQ.-999.0) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(12, 0.0, RMessage, 0, 11, 'WIND HEIGHT')
      ENDIF
      IF (Z0.EQ.-999.0) THEN
        LInadequateData = .TRUE.
        IF (Site%MetSite) THEN
          CALL SSMessage(12, 0.0, RMessage, 0, 6, 'Z0 (M)')
        ELSE
          CALL SSMessage(12, 0.0, RMessage, 0, 6, 'Z0 (D)')
        ENDIF
      ENDIF
      IF (MetIn%U.EQ.-999.0) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(12, 0.0, RMessage, 0, 10, 'WIND SPEED')
      ENDIF
      IF (LInadequateData) THEN
        IF (LSeq) THEN
          Site%NU(24)  = MetOut%NU
          Site%T0K(24) = MetOut%T0K
          Site%S(24)   = MetOut%S
          IF (MetOut%FTheta0.NE.-999.0) THEN
            Site%FTheta0(24) = MetOut%FTheta0
          ELSEIF (MetOut%S.NE.-999.0 .AND. RK.NE.-999.0) THEN
!           In the following call we use an input friction velocity of
!           1.0. Because the result is used only when the surface
!           sensible heat flux is non-negative, the value adopted has
!           no effect on the results.
            CALL FTheta0Calc(1.0, 0.0, MetOut%T0K, MetOut%S, Cl,         &
                             RK, R, Alpha, 0.0, 0.0, 0.0, 0.0, UStarDum, &
                             Site%FTheta0(24), RecipLMODum)
            IF (Site%FTheta0(24).LT.0.0) Site%FTheta0(24) = -999.0
          ENDIF
        ENDIF
        RETURN
      ENDIF

!==========Calculate UG2B, U2BStar and U3B and update ZLast and
!          Z0Last.

!3.0 CAM 19/11/98  SAVE TIME IN MET PROCESSOR
!      IF (Site.UG2B   .EQ.-999.0 .OR.
!     +    Site.U2BSTAR.EQ.-999.0 .OR.
!     +    Site.U3B    .EQ.-999.0 .OR.
!     +    Z .NE.Site.ZLast .OR.
!     +    Z0.NE.Site.Z0Last) THEN
!        CALL SetUpLimits(Z, Z0, Site.AbsF, Site.UG2B, Site.U2BStar,
!     +                   Site.U3B)
!        Site.ZLast  = Z
!        Site.Z0Last = Z0
!      ENDIF
      IF ((Site%UG2B   .EQ.-999.0 .AND. Z.EQ.1000.0) .OR.                &
          (Site%U2BSTAR.EQ.-999.0 .AND. (Z.LT.1000.0.AND.Z.GT.0.0)) .OR. &
          (Site%U3B    .EQ.-999.0 .AND. (Z.LT.1000.0.AND.Z.GT.0.0)) .OR. &
          Z .NE.Site%ZLast .OR.                                          &
          Z0.NE.Site%Z0Last) THEN
        CALL SetUpLimits(Z, Z0, Site%AbsF, Site%UG2B, Site%U2BStar, Site%U3B)
        Site%ZLast  = Z
        Site%Z0Last = Z0
      ENDIF

!==========Estimate UStar, FTheta0 and RecipLMO.

!..........Case 1: RecipLMO available.

      IF (MetOut%RecipLMO.NE.-999.0) THEN
        CALL UCalc1(MetIn%U, Z, MetOut%RecipLMO, MetOut%T0K, Z0, &
                    MetIn%UGStar, Site%AbsF, MetOut%UStar,       &
                    MetOut%FTheta0, IterationCount)

!..........Case 2: FTheta0 available, but not RecipLMO.

      ELSEIF (MetOut%FTheta0.NE.-999.0) THEN
        CALL FTheta0Limit(MetIn%U, Z, MetOut%T0K, MetIn%UGStar, &
                          Site%UG2B, Site%U3B, MetOut%FTheta0)
        CALL UCalc2(MetIn%U, Z, MetOut%FTheta0, MetOut%T0K, Z0, &
                    MetIn%UGStar, Site%AbsF, MetOut%UStar,      &
                    MetOut%RecipLMO, IterationCount)

!..........Case 3: RecipLMO and FTheta0 not available, but ThetaStar
!                  available and >or= 0.

      ELSEIF (MetOut%ThetaStar.GE.0.0) THEN
        CALL UCalc3(MetIn%U, Z, MetOut%ThetaStar, MetOut%T0K, Z0,    &
                    MetIn%UGStar, Site%AbsF, MetOut%UStar,           &
                    MetOut%FTheta0, MetOut%RecipLMO, IterationCount)

!..........Case 4: RecipLMO and FTheta0 not available, but S and
!          RK are.

      ELSEIF (MetOut%S.NE.-999.0 .AND. RK.NE.-999.0) THEN
        CALL FTheta0Calc(MetIn%U, Z, MetOut%T0K, MetOut%S, Cl,       &
                         RK, R, Alpha, Z0, MetIn%UGStar, Site%AbsF,  &
                         Site%U2BStar, MetOut%UStar, MetOut%FTheta0, &
                         MetOut%RecipLMO)

!..........Case 5: Neither RecipLMO nor FTheta0 nor S and
!          RK available.

      ELSE
        IF (LSeq) THEN
          Site%NU(24)  = MetOut%NU
          Site%T0K(24) = MetOut%T0K
          Site%S(24)   = MetOut%S
          IF (Z.EQ.0.0) Site%UStar(24) = MetIn%U
        ENDIF
        LInadequateData = .TRUE.
        CALL SSMessage(13, 0.0, RMessage, 0, 0, ' ')
        RETURN
      ENDIF

!==========Revise UStar, FTheta0 and RecipLMO if RecipLMO > RecipLMOMax.

      If (MetOut%UStar /= 0.0) Then
        IF (RecipLMOMax.NE.-999.0 .AND. RecipLMOMax.LE.1.0E6 .AND. &
           MetOut%RecipLMO.GT.RecipLMOMax) THEN
          LRecipLMOMaxUsed = .TRUE.
          MetOut%RecipLMO  = RecipLMOMax
          CALL UCalc1(MetIn%U, Z, MetOut%RecipLMO, MetOut%T0K, Z0, &
                      MetIn%UGStar, Site%AbsF, MetOut%UStar,       &
                      MetOut%FTheta0, IterationCount)
        ELSE
          LRecipLMOMaxUsed = .FALSE.
        ENDIF
      End If

!==========Calculate ThetaStar.

      IF (MetOut%UStar.NE.0.0) THEN
        MetOut%ThetaStar = -MetOut%FTheta0/(RHOA*CP*MetOut%UStar)
      ELSEIF (MetOut%FTheta0.GT.0.0) THEN
        MetOut%ThetaStar = -200000.0
      ELSE
        MetOut%ThetaStar = 200000.0
      ENDIF

!==========Add latest values to arrays of UStar, FTheta0, S, NU and T0K
!          values and fill in gaps as far as possible.

      Site%UStar(24)   = MetOut%UStar
      Site%FTheta0(24) = MetOut%FTheta0
      Site%S(24)       = MetOut%S
      Site%T0K(24)     = MetOut%T0K
      Site%NU(24)      = MetOut%NU
      IF (LSeq) THEN
        CALL FillSeq(MetOut%Cl, R, Alpha, MetOut%Day,                &
           MetOut%LocalMeanTime, Site%Lat, Site%UStar, Site%FTheta0, &
           Site%S, Site%NU, Site%T0K)
      ELSE
        CALL FillNonSeq(MetOut%Cl, R, Alpha, MetOut%Day,             &
           MetOut%LocalMeanTime, Site%Lat, Site%UStar, Site%FTheta0, &
           Site%S, Site%NU, Site%T0K)
      ENDIF

!==========Estimate UG, UGStar and DeltaPhi.

      CALL UGCalc(MetOut%UStar, MetIn%U, Z, MetOut%RecipLMO, Z0,       &
                  MetIn%UGStar, MetIn%DeltaPhi, Site%AbsF, Site%SignF, &
                  MetOut%UG, MetOut%UGStar, MetOut%DeltaPhi)

!==========Check for calm conditions.
!      write(6,*) 'checking for calm conditions'
      IF (MetIn%U.EQ.0.0) THEN
        LCalm = .TRUE.
 !       CALL SSMessage(23, 0.0, RMessage, 0, 0, ' ')
 !       RETURN
      ELSE
        LCalm = .FALSE.
      ENDIF

!==========Set DeltaThetaConv.

      DeltaThetaConv = -999.0

!==========Estimate H and amend DeltaThetaConv if appropriate.

      IF (MetOut%H.EQ.-999.0) THEN
!..........Stable conditions.
        IF (MetOut%FTheta0.LE.0.0) THEN
          CALL HStableCalcDawnDusk(Site%AbsF, Site%UStar, Site%FTheta0, &
                                   Site%S, Site%NU, Site%T0K, MetOut%H)
        ELSE
!..........Unstable conditions.
          CALL HUnstableCalc(Site%AbsF, Site%UStar, Site%FTheta0, &
                             Site%S, Site%NU, Site%T0K, MetOut%H, &
                             DeltaThetaConv)
        ENDIF
!..........Limit H.
        H = MetOut%H
        CALL HLimits(MetOut%H)
        IF (MetOut%H.NE.H) DeltaThetaConv = -999.0
!..........Correct for RecipLMOMax.
      ELSEIF (LRecipLMOMaxUsed) THEN
        CALL HStableCalcDawnDusk(Site%AbsF, Site%UStar, Site%FTheta0, &
                                 Site%S, Site%NU, Site%T0K, H)
        CALL HLimits(H)
        MetOut%H = MAX(MetOut%H, H)
      ENDIF

!==========Estimate DeltaTheta.

      IF (MetOut%DeltaTheta.EQ.-999.0)                            &
         CALL DeltaThetaCalc(MetOut%FTheta0, MetOut%H, MetOut%NU, &
                             MetOut%T0K, DeltaThetaConv,          &
                             MetOut%DeltaTheta)

!==========Check Phi available.

      IF (MetIn%Phi.EQ.-999.0) THEN
        LInadequateData = .TRUE.
        CALL SSMessage(12, 0.0, RMessage, 0, 9, 'WIND DIRN')
        RETURN
      ENDIF

!==========Calculate Phi0 and PhiG.

      CALL Angles(MetIn%Phi, Z, MetOut%DeltaPhi, MetOut%H, MetOut%Phi0, MetOut%PhiG)

!==========Calculate WStar.

      CALL WStarCalc(MetOut%UStar, MetOut%FTheta0, MetOut%RecipLMO, &
                     MetOut%T0K, MetOut%H,                          &
                     MetOut%WStar)

!==========Calculate SigmaTheta.

      IF (MetOut%SigmaThetaDeg.EQ.-999.0) THEN
        CALL SigmaThetaCalc(TSample, MetOut%UStar, MetOut%RecipLMO, Z0, &
                            MetOut%SigmaThetaRad)
        MetOut%SigmaThetaDeg = (180.0/PI)*MetOut%SigmaThetaRad
      ELSE
        MetOut%SigmaThetaRad = (PI/180.0)*MetOut%SigmaThetaDeg
      ENDIF

!==========Estimate humidity variables.

      IF (MetOut%Q0.EQ.-999.0) THEN
        IF (MetIn%RH0.NE.-999.0 .AND. MetIn%T0C.NE.-999.0) &
           CALL Q0Calc(MetOut%T0K, MetIn%RH0, MetOut%Q0)
      ENDIF
      IF (MetOut%RHU    .EQ.-999.0) MetOut%RHU    = RHUDef
      IF (MetOut%DRHDZU .EQ.-999.0) MetOut%DRHDZU = DRHDZUDef
      IF (MetOut%LambdaE.EQ.-999.0) CALL LambdaECalc(MetOut%FTheta0, &
         Alpha, MetOut%T0K, MetOut%LambdaE)

      RETURN
      END Subroutine

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

      SUBROUTINE CORIOLISCALC(RLAT,ABSF,SIGNF)
!***********************************************************************
! This subroutine calculates ABSF and SIGNF.                           *
!                                                                      *
! Parameters:                                                          *
! ABSFMIN minimum value imposed on the Coriolis parameter.             *
! PI      pi.                                                          *
!                                                                      *
! Input Variables:                                                     *
! RLAT latitude of source (degrees, being positive in the northern     *
!         hemisphere).                                                 *
!                                                                      *
! Output Variables:                                                    *
! ABSF  absolute value of the Coriolis parameter.                      *
! SIGNF sign of the Coriolis parameter.                                *
!                                                                      *
! Local Variables:                                                     *
! F        Coriolis parameter.                                         *
! RMESSAGE array used to pass information to the message subroutine.   *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      REAL*4 ABSFMIN,PI
      PARAMETER (ABSFMIN = 0.00005)
!..........Items in argument list.
      REAL*4 RLAT,ABSF,SIGNF
!..........Local variables.
      REAL*4 F,RMESSAGE(2)
!..........Pi.
      PI = 4.0*ATAN(1.0)

!==========Calculate ABSF and SIGNF.

      F = (4.0*PI/86400.0)*SIN(2.0*PI*RLAT/360.0)
      ABSF = ABS(F)
      SIGNF = SIGN(1.0,F)
      IF (ABSF.LT.ABSFMIN) THEN
        RMESSAGE(1) = ABSF
        ABSF = ABSFMIN
        RMESSAGE(2) = ABSF
        CALL SSMessage(16,0.0,RMESSAGE,0,0,' ')
      ENDIF

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE SETUPLIMITS(Z,Z0,ABSF,UG2B,U2BSTAR,U3B)
!***********************************************************************
! This subroutine calculates the minimum allowed value of UG2B,        *
! U2BSTAR and U3B.                                                     *
!                                                                      *
! Parameters:                                                          *
! AA )           parameters defining the assumed surface layer profile *
! BB )           in stable conditions.                                 *
! CC )                                                                 *
! CP             specific heat capacity of dry air at constant         *
!                   pressure.                                          *
! G              gravitational acceleration.                           *
! ITERATIONLIMIT maximum number of iterations.                         *
! RHOA           density of air.                                       *
! T0TEST         temperature adopted to test convergence of USTAR      *
!                   calculation.                                       *
! UTEST          wind speed adopted to test convergence of USTAR       *
!                   calculation.                                       *
! VK             von Karman's constant.                                *
!                                                                      *
! Input Variables:                                                     *
! Z    height of the wind measurement (1000.0 is used to indicate      *
!         geostrophic wind, 0.0 to indicate friction velocity).        *
! Z0   roughness length.                                               *
! ABSF absolute value of the Coriolis parameter.                       *
!                                                                      *
! Output Variables:                                                    *
! UG2B    minimum allowed value of U**2/(-surface buoyancy flux), for  *
!            U a geostrophic wind.                                     *
! U2BSTAR minimum allowed value of U**2/(buoyancy scale), for U a      *
!            surface layer wind at height Z with roughness length Z0.  *
! U3B     minimum allowed value of U**3/(-surface buoyancy flux), for  *
!            U a surface layer wind at height Z with roughness length  *
!            Z0.                                                       *
!                                                                      *
! Local Variables:                                                     *
! B              surface buoyancy flux.                                *
! BSTAR          buoyancy scale.                                       *
! FTHETA0        surface sensible heat flux.                           *
! ITERATIONCOUNT count of number of iterations required.               *
! LCONVERGE      indicates that the calculation of USTAR in UCALC2 or  *
!                   UCALC3 has converged.                              *
! TEMP           intermediate quantity in calculation, stored          *
!                         temporarily.                                 *
! THETASTAR      temperature scale.                                    *
! U2BSTARLOWER ) lower and upper bounds on U2BSTAR in bisection        *
! U2BSTARUPPER )    calculation.                                       *
! U3BUPPER )     lower and upper bounds on U3B in bisection            *
! U3BLOWER )        calculation.                                       *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      INTEGER*4 ITERATIONLIMIT
      REAL*4 AA,BB,CC,CP,G,RHOA,T0TEST,UTEST,VK
      PARAMETER (AA = 0.7, BB = 0.75, CC = 5.0, CP = 1012.0, G = 9.807, &
                 ITERATIONLIMIT = 30, RHOA = 1.225, T0TEST = 288.15,    &
                 UTEST = 40.0, VK = 0.4)
!..........Items in argument list.
      REAL*4 Z,Z0,ABSF,UG2B,U2BSTAR,U3B
!..........Local variables.
      INTEGER*4 ITERATIONCOUNT
      REAL*4 B,BSTAR,FTHETA0,TEMP,THETASTAR,U2BSTARUPPER,U2BSTARLOWER, &
             U3BUPPER,U3BLOWER
      LOGICAL LCONVERGE

!==========Calculate UG2B, U2BSTAR and U3B.

      IF (Z.EQ.1000.0) THEN

!..........Calculate UG2B.
        UG2B = 1.0/(0.8*0.145*ABSF)

      ELSEIF (Z.NE.0.0) THEN

!..........Calculate U2BSTAR.
        TEMP = 4.0*Z*ALOG((Z + Z0)/Z0)/VK
        U2BSTARLOWER = TEMP*(AA - BB*EXP(-CC-2.0))
        U2BSTARUPPER = TEMP*(AA + BB + BB*CC)
        DO ITERATIONCOUNT=1,ITERATIONLIMIT
          U2BSTAR = (U2BSTARLOWER + U2BSTARUPPER)*0.5
          BSTAR = UTEST**2/U2BSTAR
          THETASTAR = BSTAR*T0TEST/G
          CALL UCALC3A(UTEST,Z,THETASTAR,T0TEST,Z0,LCONVERGE)
          IF (LCONVERGE) THEN
            U2BSTARUPPER = U2BSTAR
          ELSE
            U2BSTARLOWER = U2BSTAR
          ENDIF
        ENDDO
        U2BSTAR = U2BSTARUPPER

!..........Calculate U3B.
        TEMP = 27.0*Z*(ALOG((Z + Z0)/Z0))**2/(4*VK**2)
        U3BLOWER = TEMP*(AA - BB*EXP(-CC-2.0))
        U3BUPPER = TEMP*(AA + BB + BB*CC)
        DO ITERATIONCOUNT=1,ITERATIONLIMIT
          U3B = (U3BLOWER + U3BUPPER)*0.5
          B = -UTEST**3/U3B
          FTHETA0 = B*RHOA*CP*T0TEST/G
          CALL UCALC2A(UTEST,Z,FTHETA0,T0TEST,Z0,LCONVERGE)
          IF (LCONVERGE) THEN
            U3BUPPER = U3B
          ELSE
            U3BLOWER = U3B
          ENDIF
        ENDDO
        U3B = U3BUPPER
      ENDIF

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE UCALC2A(U,Z,FTHETA0,T0,Z0,LCONVERGE)
!***********************************************************************
! This subroutine estimates USTAR from a surface layer wind as in      *
! UCALC2 and identifies whether or not the calculation converges.      *
!                                                                      *
! Parameters:                                                          *
! AA )             parameters defining the assumed surface layer       *
! BB )                profile in stable conditions.                    *
! CC )                                                                 *
! ACCURACYREQUIRED accuracy required in solving for USTAR.             *
! CP               specific heat capacity of dry air at constant       *
!                     pressure.                                        *
! G                gravitational acceleration.                         *
! ITERATIONLIMIT   maximum number of iterations.                       *
! RHOA             density of air.                                     *
! VK               von Karman's constant.                              *
!                                                                      *
! Input Variables:                                                     *
! U       wind speed.                                                  *
! Z       height of wind speed U (1000.0 is used to indicate           *
!            geostrophic wind, 0.0 to indicate friction velocity).     *
! FTHETA0 surface sensible heat flux.                                  *
! T0      near surface temperature (degrees K).                        *
! Z0      roughness length.                                            *
!                                                                      *
! Output Variables:                                                    *
! LCONVERGE indicates that the calculation of USTAR has converged.     *
!                                                                      *
! Local Variables:                                                     *
! B              surface buoyancy flux.                                *
! ITERATIONCOUNT number of iterations required.                        *
! RECIPLMO       reciprocal of the Monin-Obukhov length.               *
! TEMP           intermediate quantity in calculation, stored          *
!                   temporarily.                                       *
! USTAR          friction velocity.                                    *
! USTAROLD       previous value of USTAR in iteration.                 *
! UUSTAR         U normalised by USTAR.                                *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      INTEGER*4 ITERATIONLIMIT
      REAL*4 AA,BB,CC,ACCURACYREQUIRED,CP,G,RHOA,VK
      PARAMETER (AA = 0.7, BB = 0.75, CC = 5.0,                     &
                 ACCURACYREQUIRED = 1.0E-4, CP = 1012.0, G = 9.807, &
                 ITERATIONLIMIT = 50, RHOA = 1.225, VK = 0.4)
!..........Items in argument list.
      REAL*4 U,Z,FTHETA0,T0,Z0
      LOGICAL LCONVERGE
!..........Local variables.
      INTEGER*4 ITERATIONCOUNT
      REAL*4 B,RECIPLMO,TEMP,USTAR,USTAROLD,UUSTAR

!==========Estimate USTAR.

!..........Initialise LCONVERGE.
      LCONVERGE = .FALSE.
!..........Calculate B.
      B = FTHETA0*G/(RHOA*CP*T0)
!..........Only stable cases treated. We use an iteration technique.
      USTAR = U*VK/ALOG((Z + Z0)/Z0)
      DO ITERATIONCOUNT=1,ITERATIONLIMIT
        USTAROLD = USTAR
        RECIPLMO = -VK*B/(USTAR**3)
        TEMP = Z*RECIPLMO*(AA - BB*EXP(-CC-2.0))
        IF ((TEMP + ALOG((Z + Z0)/Z0)).GT.VK*U/USTAR.AND. &
            TEMP.GT.0.5*ALOG((Z + Z0)/Z0)) GOTO 20
        CALL SURFACELAYERSTABLE(RECIPLMO,Z,Z0,UUSTAR)
        USTAR = U/UUSTAR
        IF (ABS(USTAR - USTAROLD).LT.ACCURACYREQUIRED) THEN
          LCONVERGE = .TRUE.
          GOTO 20
        ENDIF
      ENDDO
   20 CONTINUE

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE UCALC3A(U,Z,THETASTAR,T0,Z0,LCONVERGE)
!***********************************************************************
! This subroutine estimates USTAR from a surface layer wind as in      *
! UCALC3 and identifies whether or not the calculation converges.      *
!                                                                      *
! Parameters:                                                          *
! AA )             parameters defining the assumed surface layer       *
! BB )                profile in stable conditions.                    *
! CC )                                                                 *
! ACCURACYREQUIRED accuracy required in solving for USTAR.             *
! CP               specific heat capacity of dry air at constant       *
!                     pressure.                                        *
! G                gravitational acceleration.                         *
! ITERATIONLIMIT   maximum number of iterations.                       *
! RHOA             density of air.                                     *
! VK               von Karman's constant.                              *
!                                                                      *
! Input Variables:                                                     *
! U         wind speed.                                                *
! Z         height of wind speed U (1000.0 is used to indicate         *
!              geostrophic wind, 0.0 to indicate friction velocity).   *
! THETASTAR temperature scale.                                         *
! T0        near surface temperature (degrees K).                      *
! Z0        roughness length.                                          *
!                                                                      *
! Output Variables:                                                    *
! LCONVERGE indicates that the calculation of USTAR has converged.     *
!                                                                      *
! Local Variables:                                                     *
! BSTAR          buoyancy scale.                                       *
! ITERATIONCOUNT number of iterations required.                        *
! RECIPLMO       reciprocal of the Monin-Obukhov length.               *
! TEMP           intermediate quantity in calculation, stored          *
!                   temporarily.                                       *
! USTAR          friction velocity.                                    *
! USTAROLD       previous value of USTAR in iteration.                 *
! UUSTAR         U normalised by USTAR.                                *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      INTEGER*4 ITERATIONLIMIT
      REAL*4 AA,BB,CC,ACCURACYREQUIRED,CP,G,RHOA,VK
      PARAMETER (AA = 0.7, BB = 0.75, CC = 5.0,                     &
                 ACCURACYREQUIRED = 1.0E-4, CP = 1012.0, G = 9.807, &
                 ITERATIONLIMIT = 50, RHOA = 1.225, VK = 0.4)
!..........Items in argument list.
      REAL*4 U,Z,THETASTAR,T0,Z0
      LOGICAL LCONVERGE
!..........Local variables.
      INTEGER*4 ITERATIONCOUNT
      REAL*4 BSTAR,RECIPLMO,TEMP,USTAR,USTAROLD,UUSTAR

!==========Estimate USTAR.

!..........Initialise LCONVERGE.
      LCONVERGE = .FALSE.
!..........Calculate BSTAR.
      BSTAR = THETASTAR*G/T0
!..........Only stable cases treated. We use an iteration technique.
      USTAR = U*VK/ALOG((Z + Z0)/Z0)
      DO ITERATIONCOUNT=1,ITERATIONLIMIT
        USTAROLD = USTAR
        RECIPLMO = VK*BSTAR/USTAR**2
        TEMP = Z*RECIPLMO*(AA - BB*EXP(-CC-2.0))
        IF ((TEMP + ALOG((Z + Z0)/Z0)).GT.VK*U/USTAR.AND. &
            TEMP.GT.ALOG((Z + Z0)/Z0)) GOTO 20
        CALL SURFACELAYERSTABLE(RECIPLMO,Z,Z0,UUSTAR)
        USTAR = U/UUSTAR
        IF (ABS(USTAR - USTAROLD).LT.ACCURACYREQUIRED) THEN
          LCONVERGE = .TRUE.
          GOTO 20
        ENDIF
      ENDDO
   20 CONTINUE

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE FTHETA0LIMIT(U,Z,T0,UGSTARIN,UG2B,U3B,FTHETA0)
!***********************************************************************
! This subroutine calculates the minimum allowed value of FTHETA0 and  *
! amends FTHETA0 if necessary.                                         *
!                                                                      *
! Parameters:                                                          *
! CP   specific heat capacity of dry air at constant pressure.         *
! G    gravitational acceleration.                                     *
! RHOA density of air.                                                 *
!                                                                      *
! Input Variables:                                                     *
! U        wind speed.                                                 *
! Z        height of wind speed U (1000.0 is used to indicate          *
!             geostrophic wind, 0.0 to indicate friction velocity).    *
! T0       near surface temperature (degrees K).                       *
! UGSTARIN geostrophic wind speed normalised by the friction velocity  *
!             as given in the met file. -999.0 indicates no value      *
!             given.                                                   *
! UG2B     minimum allowed value of U**2/(-surface buoyancy flux),     *
!             for U a geostrophic wind.                                *
! U3B      minimum allowed value of U**3/(-surface buoyancy flux),     *
!             for U a surface layer wind at height Z with appropriate  *
!             roughness length.                                        *
!                                                                      *
! Input/Output Variables:                                              *
! FTHETA0 surface sensible heat flux.                                  *
!                                                                      *
! Local Variables:                                                     *
! BMIN       minimum allowed value of surface buoyancy flux.           *
! FTHETA0MIN minimum allowed value of surface sensible heat flux.      *
! RMESSAGE   array used to pass information to the message subroutine. *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      REAL*4 CP,G,RHOA
      PARAMETER (CP = 1012.0, G = 9.807, RHOA = 1.225)

!..........Items in argument list.
      REAL*4 U,Z,T0,UGSTARIN,UG2B,U3B,FTHETA0
!..........Local variables.
      REAL*4 BMIN,FTHETA0MIN,RMESSAGE(2)

!==========Calculate minimum value of the surface sensible heat flux
!          and amend FTHETA0 if necessary.

      IF (Z.EQ.1000.0) THEN
        IF (UGSTARIN.EQ.-999.0) THEN
          BMIN = -U*U/UG2B
          FTHETA0MIN = BMIN*RHOA*CP*T0/G
          IF (FTHETA0.LT.FTHETA0MIN) THEN
            RMESSAGE(1) = FTHETA0
            FTHETA0 = FTHETA0MIN
            RMESSAGE(2) = FTHETA0
            CALL SSMessage(18,0.0,RMESSAGE,0,0,' ')
          ENDIF
        ENDIF
      ELSEIF (Z.NE.0.0) THEN
        BMIN = -U*U*U/U3B
        FTHETA0MIN = BMIN*RHOA*CP*T0/G
        IF (FTHETA0.LT.FTHETA0MIN) FTHETA0 = FTHETA0MIN
      ENDIF

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE THETASTARLIMIT(U,Z,T0,U2BSTAR,THETASTAR)
!***********************************************************************
! This subroutine calculates the maximum allowed value of THETASTAR    *
! and amends THETASTAR if necessary.                                   *
!                                                                      *
! Parameters:                                                          *
! G gravitational acceleration.                                        *
!                                                                      *
! Input Variables:                                                     *
! U       wind speed.                                                  *
! Z       height of wind speed U (1000.0 is used to indicate           *
!            geostrophic wind, 0.0 to indicate friction velocity).     *
! T0      near surface temperature (degrees K).                        *
! U2BSTAR minimum allowed value of U**2/(buoyancy scale), for U a      *
!            surface layer wind at height Z with appropriate roughness *
!            length.                                                   *
!                                                                      *
! Input/Output Variables:                                              *
! THETASTAR temperature scale.                                         *
!                                                                      *
! Local Variables:                                                     *
! BSTARMAX     maximum allowed value of buoyancy scale.                *
! THETASTARMAX maximum allowed value of temperature scale.             *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      REAL*4 G
      PARAMETER (G = 9.807)
!..........Items in argument list.
      REAL*4 U,Z,T0,U2BSTAR,THETASTAR
!..........Local variables.
      REAL*4 BSTARMAX,THETASTARMAX

!==========Calculate maximum value of the temperature scale
!          and amend THETASTAR if necessary.

      IF (Z.NE.1000.0.AND.Z.NE.0.0) THEN
        BSTARMAX = U*U/U2BSTAR
        THETASTARMAX = BSTARMAX*T0/G
        IF (THETASTAR.GT.THETASTARMAX) THETASTAR = THETASTARMAX
      ENDIF

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE UCALC1(U,Z,RECIPLMO,T0,Z0,UGSTARIN,ABSF,USTAR,FTHETA0, &
                        ITERATIONCOUNT)
!***********************************************************************
! This subroutine estimates USTAR and FTHETA0 from U, Z, RECIPLMO, T0, *
! Z0, UGSTARIN and ABSF.                                               *
!                                                                      *
! Parameters:                                                          *
! ACCURACYREQUIRED accuracy required in solving for USTAR.             *
! CP               specific heat capacity of dry air at constant       *
!                     pressure.                                        *
! G                gravitational acceleration.                         *
! ITERATIONLIMIT   maximum number of iterations.                       *
! RHOA             density of air.                                     *
! SIGNF            sign of the Coriolis parameter used in calls to     *
!                     ROSSBYCOEFFSTABLE and ROSSBYCOEFUNSTABLE (for    *
!                     the purposes of this routine the value is        *
!                     irrelevant provided it is either 1 or -1).       *
! VK               von Karman's constant.                              *
!                                                                      *
! Input Variables:                                                     *
! U        wind speed.                                                 *
! Z        height of wind speed U (1000.0 is used to indicate          *
!             geostrophic wind, 0.0 to indicate friction velocity).    *
! RECIPLMO reciprocal of the Monin-Obukhov length.                     *
! T0       near surface temperature (degrees K).                       *
! Z0       roughness length.                                           *
! UGSTARIN geostrophic wind speed normalised by the friction velocity  *
!             as given in the met file. -999.0 indicates no value      *
!             given.                                                   *
! ABSF     absolute value of the Coriolis parameter.                   *
!                                                                      *
! Output Variables:                                                    *
! USTAR          friction velocity.                                    *
! FTHETA0        surface sensible heat flux.                           *
! ITERATIONCOUNT number of iterations required.                        *
!                                                                      *
! Local Variables:                                                     *
! DELTAPHI   geostrophic wind direction minus surface wind direction   *
!               (degrees). This quantity is returned by the            *
!               subroutines ROSSBYCOEFFSTABLE and ROSSBYCOEFFUNSTABLE  *
!               but is not used here.                                  *
! RMESSAGE   array used to pass information to the message subroutine. *
! UGSTAR     geostrophic wind speed normalised by the friction         *
!               velocity.                                              *
! USTARLOWER lower bound on USTAR in bisection calculation.            *
! USTAROLD   previous value of USTAR in iteration.                     *
! USTARUPPER upper bound on USTAR in bisection calculation.            *
! UTEST      U calculated from an estimate of USTAR.                   *
! UUSTAR     U normalised by USTAR.                                    *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      INTEGER*4 ITERATIONLIMIT
      REAL*4 ACCURACYREQUIRED,CP,G,RHOA,SIGNF,VK
      PARAMETER (ACCURACYREQUIRED = 1.0E-4, CP = 1012.0, G = 9.807, &
                 ITERATIONLIMIT = 70, RHOA = 1.225, SIGNF = 1.0,    &
                 VK = 0.4)
!..........Items in argument list.
      INTEGER*4 ITERATIONCOUNT
      REAL*4 U,Z,RECIPLMO,T0,Z0,UGSTARIN,ABSF,USTAR,FTHETA0
!..........Local variables.
      REAL*4 DELTAPHI,RMESSAGE(5),UGSTAR,USTARLOWER,USTAROLD,USTARUPPER, &
             UTEST,UUSTAR

!***********************************************************************
!          Case 1: U is a friction velocity or U = 0.0.                *
!***********************************************************************

      IF (Z.EQ.0.0.OR.U.EQ.0.0) THEN

!==========Set USTAR = U.

        USTAR = U

!***********************************************************************
!          Case 2: U is a geostrophic wind.                            *
!***********************************************************************

      ELSEIF (Z.EQ.1000.0) THEN

!==========Estimate USTAR.

        IF (UGSTARIN.NE.-999.0) THEN
          USTAR = U/UGSTARIN
        ELSE
!..........Stable case. We use an iteration technique.
          IF (RECIPLMO.GE.0.0) THEN
            USTAR = U*VK/ALOG((250.0 + Z0)/Z0)
            DO ITERATIONCOUNT=1,ITERATIONLIMIT
              USTAROLD = USTAR
              CALL ROSSBYCOEFFSTABLE(USTAR,RECIPLMO,Z0,ABSF,SIGNF, UGSTAR,DELTAPHI)
              USTAR = U/UGSTAR
              IF (ABS(USTAR - USTAROLD).LT.ACCURACYREQUIRED) GOTO 10
            ENDDO
!..........Unstable case. We use a bisection technique.
          ELSE
            USTARUPPER = U*VK/SQRT(1.38*1.38 + (ALOG(100.0) - 3.69)**2)
            CALL ROSSBYCOEFFUNSTABLE(USTARUPPER,0.0,Z0,ABSF,SIGNF, UGSTAR,DELTAPHI)
            USTARLOWER = U/UGSTAR
            DO ITERATIONCOUNT=1,ITERATIONLIMIT
              USTAR = 2.0*(USTARUPPER*USTARLOWER)/(USTARUPPER + USTARLOWER)
              CALL ROSSBYCOEFFUNSTABLE(USTAR,RECIPLMO,Z0,ABSF,SIGNF, UGSTAR,DELTAPHI)
              UTEST = USTAR*UGSTAR
              IF (UTEST.LT.U) THEN
                USTARLOWER = USTAR
              ELSE
                USTARUPPER = USTAR
              ENDIF
              IF (ABS(USTARUPPER - USTARLOWER).LT.ACCURACYREQUIRED) GOTO 10
            ENDDO
          ENDIF
!..........Message to indicate lack of convergence of iteration.
          RMESSAGE(1) = U
          RMESSAGE(2) = Z
          RMESSAGE(3) = RECIPLMO
          RMESSAGE(4) = Z0
          RMESSAGE(5) = ABSF
          CALL SSMessage(17,0.0,RMESSAGE,0,5,'1/LMO')
   10     CONTINUE
        ENDIF

!***********************************************************************
!          Case 3: U is a surface layer wind.                          *
!***********************************************************************

      ELSE

!==========Estimate USTAR.

!..........Stable case.
        IF (RECIPLMO.GE.0.0) THEN
          CALL SURFACELAYERSTABLE(RECIPLMO,Z,Z0,UUSTAR)
          USTAR = U/UUSTAR
!..........Unstable case.
        ELSE
          CALL SURFACELAYERUNSTABLE(RECIPLMO,Z,Z0,UUSTAR)
          USTAR = U/UUSTAR
        ENDIF
      ENDIF

!***********************************************************************
!          Calculate FTHETA0.                                          *
!***********************************************************************

      FTHETA0 = -RHOA*CP*T0*RECIPLMO*USTAR**3/(VK*G)

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE UCALC2(U,Z,FTHETA0,T0,Z0,UGSTARIN,ABSF,USTAR,RECIPLMO, &
                        ITERATIONCOUNT)
!***********************************************************************
! This subroutine estimates USTAR and RECIPLMO from U, Z, FTHETA0, T0, *
! Z0, UGSTARIN and ABSF.                                               *
!                                                                      *
! Parameters:                                                          *
! ACCURACYREQUIRED accuracy required in solving for USTAR.             *
! CP               specific heat capacity of dry air at constant       *
!                     pressure.                                        *
! G                gravitational acceleration.                         *
! ITERATIONLIMIT   maximum number of iterations.                       *
! RHOA             density of air.                                     *
! SIGNF            sign of the Coriolis parameter used in calls to     *
!                     ROSSBYCOEFFSTABLE and ROSSBYCOEFUNSTABLE (for    *
!                     the purposes of this routine the value is        *
!                     irrelevant provided it is either 1 or -1).       *
! VK               von Karman's constant.                              *
!                                                                      *
! Input Variables:                                                     *
! U        wind speed.                                                 *
! Z        height of wind speed U (1000.0 is used to indicate          *
!             geostrophic wind, 0.0 to indicate friction velocity).    *
! FTHETA0  surface sensible heat flux.                                 *
! T0       near surface temperature (degrees K).                       *
! Z0       roughness length.                                           *
! UGSTARIN geostrophic wind speed normalised by the friction velocity  *
!             as given in the met file. -999.0 indicates no value      *
!             given.                                                   *
! ABSF     absolute value of the Coriolis parameter.                   *
!                                                                      *
! Output Variables:                                                    *
! USTAR          friction velocity.                                    *
! RECIPLMO       reciprocal of the Monin-Obukhov length.               *
! ITERATIONCOUNT number of iterations required.                        *
!                                                                      *
! Local Variables:                                                     *
! B          surface buoyancy flux.                                    *
! DELTAPHI   geostrophic wind direction minus surface wind direction   *
!               (degrees). This quantity is returned by the            *
!               subroutines ROSSBYCOEFFSTABLE and ROSSBYCOEFFUNSTABLE  *
!               but is not used here.                                  *
! RMESSAGE   array used to pass information to the message subroutine. *
! UGSTAR     geostrophic wind speed normalised by the friction         *
!               velocity.                                              *
! USTARLOWER lower bound on USTAR in bisection calculation.            *
! USTAROLD   previous value of USTAR in iteration.                     *
! USTARUPPER upper bound on USTAR in bisection calculation.            *
! UTEST      U calculated from an estimate of USTAR.                   *
! UUSTAR     U normalised by USTAR.                                    *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      INTEGER*4 ITERATIONLIMIT
      REAL*4 ACCURACYREQUIRED,CP,G,RHOA,SIGNF,VK
      PARAMETER (ACCURACYREQUIRED = 1.0E-4, CP = 1012.0, G = 9.807, &
                 ITERATIONLIMIT = 70, RHOA = 1.225, SIGNF = 1.0,    &
                 VK = 0.4)
!..........Items in argument list.
      INTEGER*4 ITERATIONCOUNT
      REAL*4 U,Z,FTHETA0,T0,Z0,UGSTARIN,ABSF,USTAR,RECIPLMO
!..........Local variables.
      REAL*4 B,DELTAPHI,RMESSAGE(5),UGSTAR,USTARLOWER,USTAROLD, &
             USTARUPPER,UTEST,UUSTAR

!***********************************************************************
!          Calculate B.                                                *
!***********************************************************************

      B = FTHETA0*G/(RHOA*CP*T0)

!***********************************************************************
!          Case 1: U is a friction velocity or U = 0.0.                *
!***********************************************************************

      IF (Z.EQ.0.0.OR.U.EQ.0.0) THEN

!==========Set USTAR = U.

        USTAR = U

!***********************************************************************
!          Case 2: U is a geostrophic wind.                            *
!***********************************************************************

      ELSEIF (Z.EQ.1000.0) THEN

!==========Estimate USTAR.

        IF (UGSTARIN.NE.-999.0) THEN
          USTAR = U/UGSTARIN
        ELSE
!..........Stable case. We use an iteration technique.
          IF (B.LE.0.0) THEN
            USTAR = U*VK/ALOG((250.0 + Z0)/Z0)
            DO ITERATIONCOUNT=1,ITERATIONLIMIT
              USTAROLD = USTAR
              RECIPLMO = -VK*B/USTAR**3
              CALL ROSSBYCOEFFSTABLE(USTAR,RECIPLMO,Z0,ABSF,SIGNF, UGSTAR,DELTAPHI)
              USTAR = U/UGSTAR
              IF (ABS(USTAR - USTAROLD).LT.ACCURACYREQUIRED) GOTO 10
            ENDDO
!..........Unstable case. We use a bisection technique.
          ELSE
            USTARUPPER = U*VK/SQRT(1.38*1.38 + (ALOG(100.0) - 3.69)**2)
            CALL ROSSBYCOEFFUNSTABLE(USTARUPPER,0.0,Z0,ABSF,SIGNF, UGSTAR,DELTAPHI)
            USTARLOWER = U/UGSTAR
            DO ITERATIONCOUNT=1,ITERATIONLIMIT
              USTAR = 2.0*(USTARUPPER*USTARLOWER)/(USTARUPPER + USTARLOWER)
              RECIPLMO = -VK*B/USTAR**3
              CALL ROSSBYCOEFFUNSTABLE(USTAR,RECIPLMO,Z0,ABSF,SIGNF, UGSTAR,DELTAPHI)
              UTEST = USTAR*UGSTAR
              IF (UTEST.LT.U) THEN
                USTARLOWER = USTAR
              ELSE
                USTARUPPER = USTAR
              ENDIF
              IF (ABS(USTARUPPER - USTARLOWER).LT.ACCURACYREQUIRED) GOTO 10
            ENDDO
          ENDIF
!..........Message to indicate lack of convergence of iteration.
          RMESSAGE(1) = U
          RMESSAGE(2) = Z
          RMESSAGE(3) = B
          RMESSAGE(4) = Z0
          RMESSAGE(5) = ABSF
          CALL SSMessage(17,0.0,RMESSAGE,0,1,'B')
   10     CONTINUE
        ENDIF

!***********************************************************************
!          Case 3: U is a surface layer wind.                          *
!***********************************************************************

      ELSE

!==========Estimate USTAR.

!..........Stable case. We use an iteration technique.
        IF (B.LE.0.0) THEN
          USTAR = U*VK/ALOG((Z + Z0)/Z0)
          DO ITERATIONCOUNT=1,ITERATIONLIMIT
            USTAROLD = USTAR
            RECIPLMO = -VK*B/USTAR**3
            CALL SURFACELAYERSTABLE(RECIPLMO,Z,Z0,UUSTAR)
            USTAR = U/UUSTAR
            IF (ABS(USTAR - USTAROLD).LT.ACCURACYREQUIRED) GOTO 20
          ENDDO
!..........Unstable case. We use a bisection technique.
        ELSE
          USTARLOWER = U*VK/ALOG((Z + Z0)/Z0)
          RECIPLMO = -VK*B/USTARLOWER**3
          CALL SURFACELAYERUNSTABLE(RECIPLMO,Z,Z0,UUSTAR)
          USTARUPPER = U/UUSTAR
          DO ITERATIONCOUNT=1,ITERATIONLIMIT
            USTAR = 2.0*(USTARUPPER*USTARLOWER)/(USTARUPPER + USTARLOWER)
            RECIPLMO = -VK*B/USTAR**3
            CALL SURFACELAYERUNSTABLE(RECIPLMO,Z,Z0,UUSTAR)
            UTEST = USTAR*UUSTAR
            IF (UTEST.LT.U) THEN
              USTARLOWER = USTAR
            ELSE
              USTARUPPER = USTAR
            ENDIF
            IF (ABS(USTARUPPER - USTARLOWER).LT.ACCURACYREQUIRED) GOTO 20
          ENDDO
        ENDIF
!..........Message to indicate lack of convergence of iteration.
        RMESSAGE(1) = U
        RMESSAGE(2) = Z
        RMESSAGE(3) = B
        RMESSAGE(4) = Z0
        RMESSAGE(5) = ABSF
        CALL SSMessage(17,0.0,RMESSAGE,0,1,'B')
   20   CONTINUE
      ENDIF

!***********************************************************************
!          Calculate RECIPLMO.                                         *
!***********************************************************************

      IF (USTAR.NE.0.0) THEN
        RECIPLMO = -VK*B/USTAR**3
      ELSEIF (B.GT.0.0) THEN
        RECIPLMO = -200000.0
      ELSE
        RECIPLMO = 200000.0
      ENDIF

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE UCALC3(U,Z,THETASTAR,T0,Z0,UGSTARIN,ABSF,USTAR,FTHETA0, &
                        RECIPLMO,ITERATIONCOUNT)
!***********************************************************************
! This subroutine estimates USTAR, FTHETA0 and RECIPLMO from U, Z,     *
! THETASTAR, T0, Z0, UGSTARIN and ABSF provided THETASTAR >or= 0.0.    *
!                                                                      *
! Parameters:                                                          *
! ACCURACYREQUIRED accuracy required in solving for USTAR.             *
! CP               specific heat capacity of dry air at constant       *
!                     pressure.                                        *
! G                gravitational acceleration.                         *
! ITERATIONLIMIT   maximum number of iterations.                       *
! RHOA             density of air.                                     *
! SIGNF            sign of the Coriolis parameter used in calls to     *
!                     ROSSBYCOEFFSTABLE and ROSSBYCOEFUNSTABLE (for    *
!                     the purposes of this routine the value is        *
!                     irrelevant provided it is either 1 or -1).       *
! VK               von Karman's constant.                              *
!                                                                      *
! Input Variables:                                                     *
! U         wind speed.                                                *
! Z         height of wind speed U (1000.0 is used to indicate         *
!              geostrophic wind, 0.0 to indicate friction velocity).   *
! THETASTAR temperature scale.                                         *
! T0        near surface temperature (degrees K).                      *
! Z0        roughness length.                                          *
! UGSTARIN geostrophic wind speed normalised by the friction velocity  *
!             as given in the met file. -999.0 indicates no value      *
!             given.                                                   *
! ABSF      absolute value of the Coriolis parameter.                  *
!                                                                      *
! Output Variables:                                                    *
! USTAR          friction velocity.                                    *
! FTHETA0        surface sensible heat flux.                           *
! RECIPLMO       reciprocal of the Monin-Obukhov length.               *
! ITERATIONCOUNT number of iterations required.                        *
!                                                                      *
! Local Variables:                                                     *
! BSTAR    buoyancy scale.                                             *
! DELTAPHI geostrophic wind direction minus surface wind direction     *
!             (degrees). This quantity is returned by the subroutine   *
!             ROSSBYCOEFFSTABLE but is not used here.                  *
! RMESSAGE array used to pass information to the message subroutine.   *
! UGSTAR   geostrophic wind speed normalised by the friction velocity. *
! USTAROLD previous value of USTAR in iteration.                       *
! UUSTAR   U normalised by USTAR.                                      *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      INTEGER*4 ITERATIONLIMIT
      REAL*4 ACCURACYREQUIRED,CP,G,RHOA,SIGNF,VK
      PARAMETER (ACCURACYREQUIRED = 1.0E-4, CP = 1012.0, G = 9.807, &
                 ITERATIONLIMIT = 70, RHOA = 1.225, SIGNF = 1.0,    &
                 VK = 0.4)
!..........Items in argument list.
      INTEGER*4 ITERATIONCOUNT
      REAL*4 U,Z,THETASTAR,T0,Z0,UGSTARIN,ABSF,USTAR,FTHETA0,RECIPLMO
!..........Local variables.
      REAL*4 BSTAR,DELTAPHI,RMESSAGE(5),UGSTAR,USTAROLD,UUSTAR

!***********************************************************************
!          Calculate BSTAR.                                            *
!***********************************************************************

      BSTAR = THETASTAR*G/T0

!***********************************************************************
!          Case 1: U is a friction velocity or U = 0.0.                *
!***********************************************************************

      IF (Z.EQ.0.0.OR.U.EQ.0.0) THEN

!==========Set USTAR = U.

        USTAR = U

!***********************************************************************
!          Case 2: U is a geostrophic wind.                            *
!***********************************************************************

      ELSEIF (Z.EQ.1000.0) THEN

!==========Estimate USTAR.

        IF (UGSTARIN.NE.-999.0) THEN
          USTAR = U/UGSTARIN
        ELSE
!..........Only stable cases treated. We use an iteration technique.
          USTAR = U*VK/ALOG((250.0 + Z0)/Z0)
          DO ITERATIONCOUNT=1,ITERATIONLIMIT
            USTAROLD = USTAR
            RECIPLMO = VK*BSTAR/USTAR**2
            CALL ROSSBYCOEFFSTABLE(USTAR,RECIPLMO,Z0,ABSF,SIGNF, UGSTAR,DELTAPHI)
            USTAR = U/UGSTAR
            IF (ABS(USTAR - USTAROLD).LT.ACCURACYREQUIRED) GOTO 10
          ENDDO
!..........Message to indicate lack of convergence of iteration.
          RMESSAGE(1) = U
          RMESSAGE(2) = Z
          RMESSAGE(3) = BSTAR
          RMESSAGE(4) = Z0
          RMESSAGE(5) = ABSF
          CALL SSMessage(17,0.0,RMESSAGE,0,5,'BSTAR')
   10     CONTINUE
        ENDIF

!***********************************************************************
!          Case 3: U is a surface layer wind.                          *
!***********************************************************************

      ELSE

!==========Estimate USTAR.

!..........Only stable cases treated. We use an iteration technique.
        USTAR = U*VK/ALOG((Z + Z0)/Z0)
        DO ITERATIONCOUNT=1,ITERATIONLIMIT
          USTAROLD = USTAR
          RECIPLMO = VK*BSTAR/USTAR**2
          CALL SURFACELAYERSTABLE(RECIPLMO,Z,Z0,UUSTAR)
          USTAR = U/UUSTAR
          IF (ABS(USTAR - USTAROLD).LT.ACCURACYREQUIRED) GOTO 20
        ENDDO
!..........Message to indicate lack of convergence of iteration.
        RMESSAGE(1) = U
        RMESSAGE(2) = Z
        RMESSAGE(3) = BSTAR
        RMESSAGE(4) = Z0
        RMESSAGE(5) = ABSF
        CALL SSMessage(17,0.0,RMESSAGE,0,5,'BSTAR')
   20   CONTINUE
      ENDIF

!***********************************************************************
!          Calculate FTHETA0 and RECIPLMO.                             *
!***********************************************************************

      IF (USTAR.NE.0.0) THEN
        RECIPLMO = VK*BSTAR/USTAR**2
      ELSE
        RECIPLMO = 200000.0
      ENDIF
      FTHETA0 = -RHOA*CP*USTAR*THETASTAR

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE FTheta0Calc(U, Z, T0K, S, Cl, RK, R, Alpha, Z0,      &
                             UGStarIn, AbsF, U2BStar, UStar, FTheta0, &
                             RecipLMO)
!-----------------------------------------------------------------------
! This subroutine estimates UStar, FTheta0 and RecipLMO from U, Z, T0K,
! Cl, R, Alpha, S, RK, Z0, UGStarIn, ABSF, U2BStar and RLat. The
! values of Z0, UGStarIn, ABSF and U2BStar given on calling this
! routine are irrelevant if Z = 0.0.
!
! Parameters:
! Vk: von Karman's constant.
!
! Input variables:
! U: Wind speed.
! Z: Height of wind speed U (1000.0 is used to indicate geostrophic wind,
!    0.0 to indicate friction velocity).
! T0K: Near surface temperature (degrees K).
! Cl: Cloud amount (oktas).
! R: Surface albedo.
! Alpha: Alpha.
! S: Sine of the solar elevation.
! RK: Incoming solar radiation.
! Z0: Roughness length.
! UGStarIn: Geostrophic wind speed normalised by the friction velocity
!    as given in the met file. -999.0 indicates no value given.
! AbsF: Absolute value of the Coriolis parameter.
! U2BStar: Minimum allowed value of U**2/(buoyancy scale), for U a surface
!    layer wind at height Z with roughness length Z0.
!
! Output variables:
! UStar: Friction velocity.
! FTheta0: Surface sensible heat flux.
! RecipLMO: Reciprocal of the Monin-Obukhov length.
!
! Local variables:
! IterationCount: Number of iterations used in subroutines UCalc2 and
!    UCalc3 (used for testing purposes only).
! FTheta01: Surface sensible heat flux estimated from ThetaStar.
! RecipLMO1: Reciprocal of the Monin-Obukhov length estimated from ThetaStar.
! ThetaStar: Temperature scale.
! UStar1: Friction velocity estimated from ThetaStar.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Parameters.
      REAL*4 VK
      PARAMETER (VK = 0.4)
!..........Input variables.
      REAL*4 U, Z, T0K, Cl, R, Alpha, S, RK, Z0, UGStarIn, AbsF, U2BStar
!..........Output variables.
      REAL*4 UStar, FTheta0, RecipLMO
!..........Local variables.
      INTEGER*4 IterationCount
      REAL*4 FTheta01, RecipLMO1, ThetaStar, UStar1

!==========Estimate UStar, FTheta0 and RecipLMO.

!..........Case 1: Daytime.
      IF (S.GT.0.0) THEN
        CALL FTheta0DayTimeCalc(T0K, Cl, R, Alpha, RK, FTheta0)
!..........Case 1.1: FTheta0 >or= 0.0.
        IF (FTheta0.GE.0.0) THEN
          CALL UCalc2(U, Z, FTheta0, T0K, Z0, UGStarIn, AbsF, UStar, &
                      RecipLMO, IterationCount)
!..........Case 1.2: FTheta0 < 0.0.
        ELSE
          ThetaStar = 0.09*(1.0 - 0.5*(Cl/8.0)**2)
          CALL ThetaStarLimit(U, Z, T0K, U2BStar, ThetaStar)
          CALL UCalc3(U, Z, ThetaStar, T0K, Z0, UGStarIn, AbsF, UStar1, &
                      FTheta01, RecipLMO1, IterationCount)
          IF (FTheta0.LE.FTheta01) THEN
            UStar    = UStar1
            FTheta0  = FTheta01
            RecipLMO = RecipLMO1
          ELSE
            CALL UCalc2(U, Z, FTheta0, T0K, Z0, UGStarIn, AbsF, UStar, &
                        RecipLMO, IterationCount)
          ENDIF
        ENDIF
!..........Case 2: Nighttime.
      ELSE
        ThetaStar = 0.09*(1.0 - 0.5*(Cl/8.0)**2)
        CALL ThetaStarLimit(U, Z, T0K, U2BStar, ThetaStar)
        CALL UCalc3(U, Z, ThetaStar, T0K, Z0, UGStarIn, AbsF, UStar, &
                    FTheta0, RecipLMO, IterationCount)
      ENDIF

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE UGCALC(USTAR,U,Z,RECIPLMO,Z0,UGSTARIN,DELTAPHIIN, &
                        ABSF,SIGNF,UG,UGSTAR,DELTAPHI)
!***********************************************************************
! This subroutine estimates UG, UGSTAR and DELTAPHI from USTAR, U, Z,  *
! RECIPLMO, Z0, UGSTARIN, DELTAPHIIN, ABSF and SIGNF.                  *
!                                                                      *
! Parameters:                                                          *
! ROMIN minimum value of RO for which no warning message is generated. *
!                                                                      *
! Input Variables:                                                     *
! USTAR      friction velocity.                                        *
! U          wind speed.                                               *
! Z          height of wind speed U (1000.0 is used to indicate        *
!               geostrophic wind, 0.0 to indicate friction velocity).  *
! RECIPLMO   reciprocal of the Monin-Obukhov length.                   *
! Z0         roughness length.                                         *
! UGSTARIN   geostrophic wind speed normalised by the friction         *
!               velocity as given in the met file. -999.0 indicates no *
!               value given.                                           *
! DELTAPHIIN geostrophic wind direction minus surface wind direction   *
!               (degrees) as given in the met file. -999.0 indicates   *
!               no value given.                                        *
! ABSF       absolute value of the Coriolis parameter.                 *
! SIGNF      sign of the Coriolis parameter.                           *
!                                                                      *
! Output Variables:                                                    *
! UG       geostrophic wind speed.                                     *
! UGSTAR   geostrophic wind speed normalised by the friction velocity. *
! DELTAPHI geostrophic wind direction minus surface wind direction     *
!             (degrees).                                               *
!                                                                      *
! Local Variables:                                                     *
! RMESSAGE array used to pass information to the message subroutine.   *
! RO       friction Rossby number.                                     *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      REAL*4 ROMIN
      PARAMETER (ROMIN = 500.0)
!..........Items in argument list.
      REAL*4 USTAR,U,Z,RECIPLMO,Z0,UGSTARIN,DELTAPHIIN,ABSF,SIGNF,UG, &
             UGSTAR,DELTAPHI
!..........Local variables.
      REAL*4 RMESSAGE(2),RO

!***********************************************************************
!          Case 1: U is a geostrophic wind.                            *
!***********************************************************************

      IF (Z.EQ.1000.0) THEN

!==========Set UG = U.

        UG = U

!==========Estimate UGSTAR and DELTAPHI.

        IF (RECIPLMO.GE.0.0) THEN
          CALL ROSSBYCOEFFSTABLE(USTAR,RECIPLMO,Z0,ABSF,SIGNF,UGSTAR, DELTAPHI)
        ELSE
          CALL ROSSBYCOEFFUNSTABLE(USTAR,RECIPLMO,Z0,ABSF,SIGNF,UGSTAR, DELTAPHI)
        ENDIF
        IF (UGSTARIN.NE.-999.0) UGSTAR = UGSTARIN
        IF (UGSTARIN.EQ.-999.0.AND.USTAR.NE.0.0) UGSTAR = UG/USTAR
        IF (DELTAPHIIN.NE.-999.0) DELTAPHI = DELTAPHIIN

!***********************************************************************
!          Case 2: U is a friction velocity or a surface layer wind.   *
!***********************************************************************

      ELSE

!==========Estimate UGSTAR and DELTAPHI.

        IF (RECIPLMO.GE.0.0) THEN
          CALL ROSSBYCOEFFSTABLE(USTAR,RECIPLMO,Z0,ABSF,SIGNF,UGSTAR, DELTAPHI)
        ELSE
          CALL ROSSBYCOEFFUNSTABLE(USTAR,RECIPLMO,Z0,ABSF,SIGNF,UGSTAR, DELTAPHI)
        ENDIF
        IF (UGSTARIN.NE.-999.0) UGSTAR = UGSTARIN
        IF (DELTAPHIIN.NE.-999.0) DELTAPHI = DELTAPHIIN

!==========Set UG = USTAR*UGSTAR.

        UG = USTAR*UGSTAR
      ENDIF

!***********************************************************************
!          Message to warn that the friction Rossby number is less     *
!          than ROMIN.                                                 *
!***********************************************************************

      IF (UGSTARIN.EQ.-999.0.OR.DELTAPHIIN.EQ.-999.0) THEN
        RO = USTAR/(ABSF*Z0)
        IF (RO.LT.ROMIN) THEN
          RMESSAGE(1) = RO
          RMESSAGE(2) = ROMIN
          CALL SSMessage(19,0.0,RMESSAGE,0,0,' ')
        ENDIF
      ENDIF

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE ANGLES(PHI,Z,DELTAPHI,H,PHI0,PHIG)
!***********************************************************************
! This subroutine calculates PHI0 and PHIG from PHI, Z and DELTAPHI.   *
!                                                                      *
! Input Variables:                                                     *
! PHI      wind direction (angle wind is coming from in degrees        *
!             clockwise from north).                                   *
! Z        height of wind direction PHI (1000.0 is used to indicate    *
!             geostrophic wind, 0.0 to indicate friction velocity).    *
! DELTAPHI geostrophic wind direction minus surface wind direction     *
!             (degrees).                                               *
! H        boundary layer depth.                                       *
!                                                                      *
! Output Variables:                                                    *
! PHI0 surface wind direction (angle wind is coming from in degrees    *
!         clockwise from north).                                       *
! PHIG geostrophic wind direction (angle wind is coming from in        *
!         degrees clockwise from north).                               *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Items in argument list.
      REAL*4 PHI,Z,DELTAPHI,H,PHI0,PHIG

!==========Calculate PHI0, PHIG.

      IF (Z.EQ.1000.0) THEN
        PHIG = PHI
        PHI0 = PHIG - DELTAPHI
        IF (PHI0.LT.0.0) THEN
          PHI0 = PHI0 + 360.0
        ENDIF
        IF (PHI0.GE.360.0) THEN
          PHI0 = PHI0 - 360.0
        ENDIF
      ELSE IF (Z.EQ.0.0) THEN
        PHI0 = PHI
        PHIG = PHI0 + DELTAPHI
        IF (PHIG.LT.0.0) THEN
          PHIG = PHIG + 360.0
        ENDIF
        IF (PHIG.GE.360.0) THEN
          PHIG = PHIG - 360.0
        ENDIF
      ELSE
        PHI0 = PHI - DELTAPHI * (Z / H)
        PHIG = PHI0 + DELTAPHI
        IF (PHIG.LT.0.0) THEN
          PHIG = PHIG + 360.0
        ENDIF
        IF (PHIG.GE.360.0) THEN
          PHIG = PHIG - 360.0
        ENDIF
      ENDIF

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE HStableCalcDawnDusk(AbsF, UStarArray, FTheta0Array, &
                                     SArray, RNUArray, T0KArray, H)
!-----------------------------------------------------------------------
! This subroutine estimates H from AbsF, UStarArray, FTheta0Array,
! SArray, RNUArray and T0KArray in stable cases, with near dawn and near
! dusk corrections.
!
! Input variables:
! AbsF: Absolute value of the Coriolis parameter.
! UStarArray, FTheta0Array, SArray, RNUArray, T0KArray: Arrays containing
!    the values of UStar, FTheta0, S, RNU and T0K over 24 hours. -999.0
!    indicates missing values. We assume UStarArray(24), FTheta0Array(24)
!    and T0KArray(24) are present and that FTheta0Array(24) <or= 0.
!    We also assume that, if FTheta0Array(23) is present and > 0, then the
!    arrays of UStar, FTheta0, NU and T0K are complete as far back as the
!    next time that FTheta0 is <or= 0 or missing.
!
! Output variables:
! H: Boundary layer depth.
!
! Local variables:
! I: Loop counter.
! DeltaThetaConv: DeltaTheta as returned by HUnstableCalc.
! FThetaArray1: Copy of FTheta0Array with last FTheta0 value set to 0.1.
! HConv: H as returned by HUnstableCalc.
!
! Note it would be logical to combine this routine with HStableCalcDawn,
! but there is then a risk of recursion and its not clear how the
! compiler will handle it.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      REAL*4 AbsF, UStarArray(24), FTheta0Array(24), SArray(24), &
             RNUArray(24), T0KArray(24)
!..........Output variables.
      REAL*4 H
!..........Local variables.
      INTEGER*4 I
      REAL*4 HConv, DeltaThetaConv, FTheta0Array1(24)

!==========Value with near dawn correction.

      CALL HStableCalcDawn(AbsF, 24, UStarArray, FTheta0Array, &
                           SArray, T0KArray, H)

!==========Near dusk correction.

      IF (FTheta0Array(23).NE.-999.0 .AND. FTheta0Array(23).GT.0.0) THEN
        DO I = 1, 23
          FTheta0Array1(I) = FTheta0Array(I)
        ENDDO
        FTheta0Array1(24) = 0.1
        CALL HUnstableCalc(AbsF, UStarArray, FTheta0Array1, SArray,   &
                           RNUArray, T0KArray, HConv, DeltaThetaConv)
        H = MIN(H, HCONV)
      ENDIF

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE HStableCalcDawn(AbsF, Index, UStarArray, FTheta0Array, &
                                 SArray, T0KArray, H)
!-----------------------------------------------------------------------
! This subroutine estimates H from AbsF, Index, UStarArray, FTheta0Array,
! SArray and T0KArray in stable cases, with near dawn correction.
!
! Input variables:
! AbsF: Absolute value of the Coriolis parameter.
! Index: Index of array elements corresponding to the time for which we
!    wish to calculate H.
! UStarArray, FTheta0Array, SArray, T0KArray: Arrays containing the values
!    of UStar, FTheta0, S and T0K over 24 hours. -999.0 indicates missing
!    values. We assume UStarArray(Index), FTheta0Array(Index) and
!    T0KArray(Index) are present and that FTheta0Array(Index) <or= 0.
!
! Output variables:
! H: Boundary layer depth.
!
! Local variables:
! I: Loop counter.
! IterationCount: Quantity returned by UCalc2 but not used here.
! RecipLMOForH: Value of RecipLMO used in calculating H.
! UStar: Quantity returned by UCalc2 but not used here.
! NearDawn: Indicates that the near dawn correction is to be made.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      INTEGER*4 Index
      REAL*4 AbsF, UStarArray(24), FTheta0Array(24), SArray(24), T0KArray(24)
!..........Output variables.
      REAL*4 H
!..........Local variables.
      INTEGER*4 I, IterationCount
      REAL*4 RecipLMOForH, UStar
      LOGICAL NearDawn

!==========Find out whether near dawn correction needed.

!..........Set NearDawn = .TRUE.
      NearDawn = .TRUE.
!..........Set NearDawn = .FALSE. if S missing.
      IF (SArray(Index).EQ.-999.0) NearDawn = .FALSE.
!..........Set NearDawn = .FALSE. if sun below horizon.
      IF (SArray(Index).LE.0.0) NearDawn = .FALSE.
!..........Set NearDawn = .FALSE. if Index = 1.
      IF (Index.EQ.1) NearDawn = .FALSE.
!..........Set NearDawn = .FALSE. if can't calculate RecipLMO at
!          previous hour.
      IF (NearDawn) THEN
        IF (UStarArray(Index - 1)  .EQ.-999.0 .OR. &
            FTheta0Array(Index - 1).EQ.-999.0 .OR. &
            T0KArray(Index - 1)    .EQ.-999.0)     &
            NearDawn = .FALSE.
      ENDIF
!..........Set NearDawn = .FALSE. if, going backwards in time, FTheta0
!          doesn't stay negative until the sun is below the horizon.
      IF (NearDawn) THEN
        DO I = Index - 1, 1, -1
          IF (FTheta0Array(I).EQ.-999.0 .OR. &
              SArray(I)      .EQ.-999.0) THEN
            NearDawn = .FALSE.
            GOTO 10
          ENDIF
          IF (FTheta0Array(I).GT.0.0) THEN
            NearDawn = .FALSE.
            GOTO 10
          ENDIF
          IF (SArray(I).LE.0.0) GOTO 10
          IF (I.EQ.1) NearDawn = .FALSE.
        ENDDO
   10   CONTINUE
      ENDIF

!==========Calculate value of RecipLMO to use in calculating H. Note
!          Z0 and UGStarIn are set to dummy values in the calls to UCalc2
!          because they are not known and because they are irrelevent
!          for this application of UCalc2.

      IF (NearDawn) THEN
        CALL UCalc2(UStarArray(Index - 1), 0.0, FTheta0Array(Index - 1), &
                    T0KArray(Index - 1), 0.0, 0.0, AbsF, UStar,          &
                    RecipLMOForH, IterationCount)
      ELSE
        CALL UCalc2(UStarArray(Index), 0.0, FTheta0Array(Index), &
                    T0KArray(Index), 0.0, 0.0, AbsF, UStar,      &
                    RecipLMOForH, IterationCount)
      ENDIF

!==========Calculate H.

      CALL HStableCalc(UStarArray(Index), RecipLMOForH, AbsF, H)

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE HStableCalc(UStar, RecipLMO, AbsF, H)
!-----------------------------------------------------------------------
! This subroutine estimates H from UStar, RecipLMO and AbsF when
! RecipLMO >or= 0.0, but without the near-dawn correction.
!
! Input variables:
! UStar: Friction velocity.
! RecipLMO: Reciprocal of the Monin-Obukhov length.
! AbsF: Absolute value of the Coriolis parameter.
!
! Output variables:
! H: Boundary layer depth.
!
! Local variables:
! Temp: Intermediate quantity in calculation, stored temporarily.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      REAL*4 UStar, RecipLMO, AbsF
!..........Output variables.
      REAL*4 H
!..........Local variables.
      REAL*4 Temp

!==========Estimate H.

      Temp = 2.28*UStar*RecipLMO/AbsF
      H = 0.6*UStar/(AbsF*(SQRT(1.0 + Temp) + 1.0))

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE HUnstableCalc(AbsF, UStarArray, FTheta0Array, SArray, &
                               RNUArray, T0KArray, H, DeltaThetaConv)
!-----------------------------------------------------------------------
! This subroutine estimates H and DeltaThetaConv from AbsF, UStarArray,
! FTheta0Array, RNUArray and T0KArray in unstable cases.
!
! Input variables:
! AbsF: Absolute value of the Coriolis parameter.
! UStarArray, FTheta0Array, SArray, RNUArray, T0KArray: Arrays containing
!    the values of UStar, FTheta0, S, NU and T0K over the 24 hours up to
!    the current hour. -999.0 indicates missing values. We assume the
!    arrays of UStar, FTheta0, NU and T0K are complete as far back as the
!    first time that FTheta0 is <or= 0 or missing.
!
! Output variables:
! H: Boundary layer depth.
! DeltaThetaConv: Temperature jump across the boundary layer top estimated
!    from the model of boundary layer growth if the value of H is the value
!    obtained from this model, -999.0 otherwise.
!
! Local variables:
! DeltaTheta0: Temperature jump across the boundary layer top used as
!    initial conditions in the subroutine BLGrowth.
! DeltaTheta1: Temperature jump across the boundary layer top returned
!    by the subroutine BLGrowth.
! HS: Boundary layer depth in last hour with FTheta0 <or= 0.0.
! H0: Boundary layer depth used as initial conditions in the subroutine
!    BLGrowth.
! H1: Boundary layer depth returned by the subroutine BLGrowth.
! LastLoop: Subscript corresponding to the most recent hour with FTheta0
!    <or= 0.0 or for which FTheta0 is missing.
! I: Loop counter.
! RMessage: Array required in the argument list of the subroutine SSMessage.
!    It has no significance in this subroutine.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      REAL*4 AbsF, UStarArray(24), FTheta0Array(24), SArray(24), &
             RNUArray(24), T0KArray(24)
!..........Output variables.
      REAL*4 H, DeltaThetaConv
!..........Local variables.
      INTEGER*4 I, LastLoop
      REAL*4 DeltaTheta0, DeltaTheta1, HS, H0, H1, RMessage(1)

!==========Find most recent hour with FTheta0 <or= 0.0 or for which FTheta0
!          is missing.

      DO I = 23, 1, -1
        LastLoop = I
!..........If FTheta0Array(I) = -999.0 or FTheta0Array(I) is
!          non-positive, exit loop.
        IF (FTheta0Array(I).LE.0.0) GOTO 10
      ENDDO
   10 CONTINUE

!==========Estimate H.

!..........If it has not been possible to fill the arrays (this is the
!          case iff FTheta0Array(LastLoop) = -999.0) or all values in
!          FTheta0Array are positive (this is the case iff
!          FTheta0Array(LastLoop) > 0.0), use neutral value for H.
      IF (FTheta0Array(LastLoop).EQ.-999.0 .OR. &
          FTheta0Array(LastLoop).GT.0.0) THEN
        CALL HStableCalc(UStarArray(24), 0.0, AbsF, H)
        DeltaThetaConv = -999.0
        CALL SSMessage(20, 0.0, RMessage, 0, 0, ' ')
      ELSE
!..........Estimate boundary layer depth at last hour with FTheta0 <or=
!          0.0 using HStableCalcDawn.
        CALL HStableCalcDawn(AbsF, LastLoop, UStarArray, FTheta0Array, &
                             SArray, T0KArray, HS)
!..........Estimate values of boundary layer depth and temperature jump
!          across the boundary layer top for the current hour using the
!          model of boundary layer growth.
        H0 = 0.0
        DeltaTheta0 = 0.0
        DO I = LastLoop + 1, 23
        CALL BLGrowth(UStarArray(I), FTheta0Array(I), RNUArray(I), &
                      T0KArray(I), H0, DeltaTheta0, 3600.0, H1,    &
                      DeltaTheta1)
        H0 = H1
        DeltaTheta0 = DeltaTheta1
        ENDDO
        CALL BLGrowth(UStarArray(24), FTheta0Array(24), RNUArray(24), &
                      T0KArray(24), H0, DeltaTheta0, 1800.0, H1,      &
                      DeltaTheta1)
!..........Estimate value of boundary layer depth for the current hour.
        IF (H1.GE.HS) THEN
          H = H1
          DeltaThetaConv = DeltaTheta1
        ELSE
          H = HS
          DeltaThetaConv = -999.0
        ENDIF
      ENDIF

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE BLGROWTH(USTAR,FTHETA0,RNU,T0,H0,DELTATHETA0,TSTEP,H1, &
                          DELTATHETA1)
!***********************************************************************
! This subroutine estimates the change over a time TSTEP in the        *
! boundary layer depth and the temperature jump across the boundary    *
! layer top using the model of boundary layer growth.                  *
!                                                                      *
! Parameters:                                                          *
! A  )               constants used in the boundary layer growth       *
! CF )                  model.                                         *
! CP                 specific heat capacity of dry air at constant     *
!                       pressure.                                      *
! G                  gravitational acceleration.                       *
! HMAX               maximum allowed increase in boundary layer depth. *
! RHOA               density of air.                                   *
! NUMBEROFBISECTIONS number of bisections used in finding H1.          *
!                                                                      *
! Input Variables:                                                     *
! USTAR       friction velocity.                                       *
! FTHETA0     surface sensible heat flux.                              *
! RNU         buoyancy frequency above the boundary layer.             *
! T0          near surface temperature (degrees K).                    *
! H0          initial boundary layer depth.                            *
! DELTATHETA0 initial temperature jump across the boundary layer top.  *
! TSTEP       time over which the evolution of boundary layer depth    *
!                and temperature jump across the boundary layer top    *
!                are estimated.                                        *
!                                                                      *
! Output Variables:                                                    *
! H1          final boundary layer depth.                              *
! DELTATHETA1 final temperature jump across the boundary layer top.    *
!                                                                      *
! Local Variables:                                                     *
! AA    )        intermediate quantities in calculation, stored        *
! BB    )           temporarily.                                       *
! CFP1  )                                                              *
! CF2P1 )                                                              *
! S1    )                                                              *
! S2    )                                                              *
! S2HAT )                                                              *
! DELTATHETA0D ) double precision versions of DELTATHETA0, H0 and H1.  *
! H0D          )                                                       *
! H1D          )                                                       *
! GAMMA          potential temperature lapse rate.                     *
! HLOWER )       upper and lower limits within which it is known that  *
! HUPPER )          H1 must lie.                                       *
! LOOP           loop counter.                                         *
! T              time required for boundary layer depth to reach a     *
!                   given value.                                       *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      INTEGER*4 NUMBEROFBISECTIONS
      REAL*4 A,CP,G,HMAX,RHOA
      REAL*8 CF
      PARAMETER (A = 5.0, CF = 0.2D0, CP = 1012.0, G = 9.807,          &
                 HMAX = 4000.0, NUMBEROFBISECTIONS = 12, RHOA = 1.225)
!..........Items in argument list.
      REAL*4 USTAR,FTHETA0,RNU,T0,H0,DELTATHETA0,TSTEP,H1,DELTATHETA1
!..........Local variables.
      INTEGER*4 LOOP
      REAL*4 S2,HLOWER,HUPPER
      REAL*8 AA,BB,CFP1,CF2P1,S1,S2HAT,DELTATHETA0D,H0D,H1D,GAMMA,T

!==========Set up variables needed for iteration.

      GAMMA = RNU*RNU*T0/G
      S1 = MAX(FTHETA0,1.0)/(RHOA*CP)
      S2 = T0*(USTAR*USTAR*USTAR)/G
      S2HAT = A*S2/S1
      H0D = H0
      DELTATHETA0D = DELTATHETA0
      CFP1 = 1.0D0 + CF
      CF2P1 = CFP1 + CF
      AA = H0D*DELTATHETA0D - 0.5D0*GAMMA*H0D*H0D
      BB = AA/S1 + GAMMA*H0D*H0D/(2.0D0*S1*CF2P1)      &
           - GAMMA*S2HAT*(H0D - S2HAT)/(S1*CF2P1*CFP1)

!==========Iteration to estimate H1.

      HUPPER = H0 + HMAX
      HLOWER = H0
      DO LOOP = 1,NUMBEROFBISECTIONS
        H1 = (HUPPER + HLOWER)*0.5
        H1D = H1
        T = AA/S1 + GAMMA*H1D*H1D/(2.0D0*S1*CF2P1)               &
            - GAMMA*S2HAT*(H1D - S2HAT)/(S1*CF2P1*CFP1)          &
            - BB*((CF*H0D + S2HAT)/(CF*H1D + S2HAT))**(1.0D0/CF)
        IF(T.LT.TSTEP)THEN
          HLOWER = H1
        ELSE
          HUPPER = H1
        ENDIF
      ENDDO

!==========Estimate DELTATHETA1.

      DELTATHETA1 = (0.5*GAMMA*H1*H1 - S1*T + AA)/H1

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE HLimits(H)
!-----------------------------------------------------------------------
! This subroutine ensures H lies between HMin and HMax.
!
! Parameters:
! HMax, HMin: Maximum and minimum limits imposed on estimates of H.
!
! Input/output variables:
! H: Boundary layer depth.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Parameters.
      REAL*4 HMax, HMin
      PARAMETER (HMax = 4000.0, HMin = 50.0)
!..........Input/output variables.
      REAL*4 H

!==========Ensure H lies between HMin and HMax.

      IF (H.LT.HMin) THEN
        H = HMin
      ENDIF
      IF (H.GT.HMax) THEN
        H = HMax
      ENDIF

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE DELTATHETACALC(FTHETA0,H,RNU,T0,DELTATHETACONV, DELTATHETA)
!***********************************************************************
! This subroutine estimates DELTATHETA from FTHETA0, H, RNU, T0 and    *
! DELTATHETACONV.                                                      *
!                                                                      *
! Parameters:                                                          *
! CF constant used in the boundary layer growth model.                 *
! G  gravitational acceleration.                                       *
!                                                                      *
! Input Variables:                                                     *
! FTHETA0        surface sensible heat flux.                           *
! H              boundary layer depth.                                 *
! RNU            buoyancy frequency above the boundary layer.          *
! T0             near surface temperature (degrees K).                 *
! DELTATHETACONV temperature jump across the boundary layer top        *
!                   estimated from the model of boundary layer growth  *
!                   if the value of H is the value obtained from this  *
!                   model, -999.0 otherwise.                           *
!                                                                      *
! Output Variables:                                                    *
! DELTATHETA temperature jump across the boundary layer top.           *
!                                                                      *
! Local Variables:                                                     *
! GAMMA potential temperature lapse rate.                              *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      REAL*4 CF,G
      PARAMETER (CF = 0.2, G = 9.807)
!..........Items in argument list.
      REAL*4 FTHETA0,H,RNU,T0,DELTATHETACONV,DELTATHETA
!..........Local variables.
      REAL*4 GAMMA

!==========Estimate DELTATHETA.

      IF (DELTATHETACONV.NE.-999.0) THEN
        DELTATHETA = DELTATHETACONV
      ELSEIF (FTHETA0.GT.0.0) THEN
        GAMMA = RNU*RNU*T0/G
        DELTATHETA = GAMMA*H*CF/(1.0 + 2.0*CF)
      ELSE
        DELTATHETA = 0.0
      ENDIF

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE WSTARCALC(USTAR,FTHETA0,RECIPLMO,T0,H,WSTAR)
!***********************************************************************
! This subroutine calculates WSTAR from USTAR, RECIPLMO and H.         *
!                                                                      *
! Parameters:                                                          *
! CP   specific heat capacity of dry air at constant pressure.         *
! G    gravitational acceleration.                                     *
! RHOA density of air.                                                 *
! VK   von Karman's constant.                                          *
!                                                                      *
! Input Variables:                                                     *
! USTAR    friction velocity.                                          *
! FTHETA0  surface sensible heat flux.                                 *
! RECIPLMO reciprocal of the Monin-Obukhov length.                     *
! T0       near surface temperature (degrees K).                       *
! H        boundary layer depth.                                       *
!                                                                      *
! Output Variables:                                                    *
! WSTAR convective velocity scale if FTHETA0 > 0.0, zero if FTHETA0    *
!          <or= 0.0.                                                   *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      REAL*4 CP,G,RHOA,VK
      PARAMETER (CP = 1012.0, G = 9.807, RHOA = 1.225, VK = 0.4)
!..........Items in argument list.
      REAL*4 USTAR,FTHETA0,RECIPLMO,T0,H,WSTAR

!==========Calculate WSTAR.

      IF (RECIPLMO.LT.0.0) THEN
        If (UStar == 0.0) Then
          WSTAR = (H*FTHETA0*G/(RHOA*CP*T0))**(1.0/3.0)
        Else
          WSTAR = USTAR*(-H*RECIPLMO/VK)**(1.0/3.0)
        End If
      ELSE
        WSTAR = 0.0
      ENDIF

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE FillSeq(Cl, R, Alpha, Day, Hour, RLat, UStarArray, &
                         FTheta0Array, SArray, RNUArray, T0KArray)
!-----------------------------------------------------------------------
! This subroutine attempts to fill in the arrays UStarArray,
! FTheta0Array, RNUArray and T0KArray in the case of sequential data.
!
! Parameters:
! InterpolationLimit: Limit to the number of consecutive hours with
!    FTheta0 missing which can be filled in by interpolation in the case
!    of sequential data.
! RNUDef: Default value for RNU.
! T0CDef: Default value for T0C.
! TCToK: Absolute temperature of the zero of the Celsius scale.
!
! Input variables:
! Cl: Cloud amount (oktas).
! R: Surface albedo.
! Alpha: Alpha.
! Day:  Day of year.
! Hour: Hour of day (local time, hours).
! RLat: Latitude of source (degrees, being positive in the northern
!    hemisphere).
!
! Input/output variables:
! UStarArray, FTheta0Array, SArray, RNUArray, T0KArray: Arrays containing
!    the values of UStar, FTheta0, NU and T0K over the 24 hours up to the
!    current hour. -999.0 indicates missing values.
!
! Local variables:
! LInterpolate: Indicates that interpolation of FTheta0 is permitted.
! I: Loop counter.
! RecipLMO: Reciprocal of the Monin-Obukhov length.
! RK: Incoming solar radiation.
! UStar: Friction velocity.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Parameters.
      INTEGER*4 InterpolationLimit
      REAL*4 RNUDef, TCToK, T0CDef
      PARAMETER (InterpolationLimit = 2, RNUDef = 0.013, &
                 TCToK = 273.15, T0CDef = 15.0)
!..........Input variables:.
      REAL*4 Cl, R, Alpha, Day, Hour, RLat
!..........Input/output variables:.
      REAL*4 UStarArray(24), FTheta0Array(24), SArray(24), RNUArray(24), T0KArray(24)
!..........Local variables.
      INTEGER*4 I
      REAL*4 RecipLMO, RK, UStar
      LOGICAL LInterpolate

!==========Initialise LInterpolate.

      LInterpolate = .TRUE.

!==========Loop to fill in arrays.

      DO I = 23, 1, -1

!==========UStarArray.

!..........If value not available, try interpolation.
        IF (UStarArray(I).EQ.-999.0) THEN
          CALL Interpolation(I, 24, UStarArray)
!..........If interpolation unsuccesful, set equal to first later value.
          IF (UStarArray(I).EQ.-999.0) UStarArray(I) = UStarArray(I + 1)
        ENDIF

!==========RNUArray.

        IF (RNUArray(I).EQ.-999.0) RNUArray(I) = RNUDef

!==========T0KArray.

        IF (T0KArray(I).EQ.-999.0) T0KArray(I) = T0CDef + TCTOK

!==========SArray.

!..........Try estimating from Day and Hour.
        IF (SArray(I).EQ.-999.0 .AND. Day.NE.-999.0 .AND. &
            Hour.NE.-999.0) CALL SCalc(RLat, Day,         &
            Hour - 24.0 + FLOAT(I), SArray(I))

!==========FTheta0Array.

!..........If value available, enable subsequent interpolations using
!          this value.
        IF (FTheta0Array(I).NE.-999.0) THEN
          LInterpolate = .TRUE.
!..........If value not available, try interpolation; if interpolation
!          fails, prevent subsequent interpolations until a new
!          value of FTheta0 is found - this is in order to prevent
!          interpolations using values of FTheta0 which are estimated
!          below.
        ELSE
          IF (LInterpolate) THEN
            CALL Interpolation(I, InterpolationLimit, FTheta0Array)
            IF (FTheta0Array(I).EQ.-999.0) LInterpolate = .FALSE.
          ENDIF
!..........If FTheta0 still unknown, try estimating from Cl and S.
          IF (Cl.NE.-999.0 .AND. SArray(I).NE.-999.0 .AND. &
              FTheta0Array(I).EQ.-999.0) THEN
            CALL KCalc(Cl, SArray(I), RK)
            CALL FTheta0Calc(UStarArray(I), 0.0, T0KArray(I), SArray(I), &
                             Cl, RK, R, Alpha, 0.0, 0.0, 0.0, 0.0,       &
                             UStar, FTheta0Array(I), RecipLMO)
          ENDIF
        ENDIF

!==========End of loop.

!..........If FTheta0Array(I) = -999.0, or FTheta0Array(I)
!          and SArray(I) non-positive, or FTheta0Array(I)
!          non-positive and SArray(I) = -999.0, exit loop.
        IF (FTheta0Array(I).EQ.-999.0 .OR.                          &
           (FTheta0Array(I).LE.0.0 .AND. SArray(I).LE.0.0)) GOTO 10
      ENDDO
   10 CONTINUE

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE INTERPOLATION(LOOP,INTERPOLATIONLIMIT,ARRAY)
!***********************************************************************
! This subroutine fills in ARRAY(LOOP) by interpolation if possible.   *
!                                                                      *
! Input Variables:                                                     *
! LOOP               subscript of array element to be interpolated.    *
! INTERPOLATIONLIMIT limit on the number of missing hours over which   *
!                       to interpolate.                                *
!                                                                      *
! Input/Output Variables:                                              *
! ARRAY array of values to be interpolated.                            *
!                                                                      *
! Local Variables:                                                     *
! LOOPLOWER ) subscripts corresponding to the array elements to be     *
! LOOPUPPER )    used in the interpolation.                            *
! LOOP1       loop counter.                                            *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Items in argument list.
      INTEGER*4 LOOP,INTERPOLATIONLIMIT
      REAL*4 ARRAY(24)
!..........Local variables.
      INTEGER*4 LOOPLOWER,LOOPUPPER,LOOP1

!==========Interpolate.

      LOOPUPPER = -999
      LOOPLOWER = -999
      DO LOOP1 = LOOP + 1,24
        IF (ARRAY(LOOP1).NE.-999.0) THEN
          LOOPUPPER = LOOP1
          GOTO 10
        ENDIF
      ENDDO
   10 CONTINUE
      DO LOOP1 = LOOP - 1,1,-1
        IF (ARRAY(LOOP1).NE.-999.0) THEN
          LOOPLOWER = LOOP1
          GOTO 20
        ENDIF
      ENDDO
   20 CONTINUE
      IF (LOOPUPPER.NE.-999.AND.LOOPLOWER.NE.-999                   &
          .AND.(LOOPUPPER - LOOPLOWER).LE.INTERPOLATIONLIMIT+1)     &
          ARRAY(LOOP) = (FLOAT(LOOP - LOOPLOWER)*ARRAY(LOOPUPPER) + &
          FLOAT(LOOPUPPER - LOOP)*ARRAY(LOOPLOWER))/                &
          FLOAT(LOOPUPPER - LOOPLOWER)

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE FillNonSeq(Cl, R, Alpha, Day, Hour, RLat, UStarArray, &
                            FTheta0Array, SArray, RNUArray, T0KArray)
!-----------------------------------------------------------------------
! This subroutine attempts to fill in the arrays UStarArray, FTheta0Array,
! RNUArray and T0KArray in the case of non-sequential data.
!
! Input variables:
! Cl: Cloud amount (oktas).
! R: Surface albedo.
! Alpha: Alpha.
! Day: Day of year.
! Hour: Hour of day (local time, hours).
! RLat: Latitude of source (degrees, being positive in the northern
!    hemisphere).
!
! Input/output variables:
! UStarArray, FTheta0Array, SArray, RNUArray, T0KArray: Arrays containing
!    the values of UStar, FTheta0, NU and T0K over the 24 hours up to the
!    current hour. -999.0 indicates missing values.
!
! Local variables:
! I: Loop counter.
! RecipLMO: Reciprocal of the Monin-Obukhov length.
! RK: Incoming solar radiation.
! UStar: Friction velocity.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      REAL*4 Cl, R, Alpha, Day, Hour, RLat
!..........Input/output variables.
      REAL*4 UStarArray(24), FTheta0Array(24), SArray(24), RNUArray(24), T0KArray(24)
!..........Local variables.
      INTEGER*4 I
      REAL*4 RecipLMO, RK, UStar

!==========Loop to fill in arrays.

      DO I = 23, 1, -1

!==========UStarArray, RNUArray and T0KArray.

!..........Set equal to current values.
        UStarArray(I) = UStarArray(24)
        RNUArray(I)   = RNUArray(24)
        T0KArray(I)   = T0KArray(24)

!==========SArray.

!..........Try estimating from Day and Hour.
        IF (Day.NE.-999.0 .AND. Hour.NE.-999.0) THEN
          CALL SCalc(RLat, Day, Hour - 24.0 + FLOAT(I), SArray(I))
        ELSE
          SArray(I) = -999.0
        ENDIF

!==========FTheta0Array.

!..........Try estimating from Cl and S.
        IF (Cl.NE.-999.0 .AND. SArray(I).NE.-999.0) THEN
          CALL KCalc(Cl, SArray(I), RK)
          CALL FTheta0Calc(UStarArray(I), 0.0, T0KArray(I), SArray(I), &
                           Cl, RK, R, Alpha, 0.0, 0.0, 0.0, 0.0,       &
                           UStar, FTheta0Array(I), RecipLMO)
        ELSE
          FTheta0Array(I) = -999.0
        ENDIF

!==========End of loop.

!..........If FTheta0Array(I) = -999.0, or FTheta0Array(I)
!          and SArray(I) non-positive, or FTheta0Array(I)
!          non-positive and SArray(I) = -999.0, exit loop.
        IF (FTheta0Array(I).EQ.-999.0 .OR.                          &
           (FTheta0Array(I).LE.0.0 .AND. SArray(I).LE.0.0)) GOTO 10
      ENDDO
   10 CONTINUE

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE SURFACELAYERSTABLE(RECIPLMO,Z,Z0,UUSTAR)
!***********************************************************************
! This subroutine estimates UUSTAR in stable conditions.               *
!                                                                      *
! Parameters:                                                          *
! AA ) parameters defining the assumed surface layer profile in stable *
! BB )    conditions.                                                  *
! CC )                                                                 *
! DD )                                                                 *
! VK   von Karman's constant.                                          *
!                                                                      *
! Input Variables:                                                     *
! RECIPLMO reciprocal of the Monin-Obukhov length.                     *
! Z        height of surface layer wind.                               *
! Z0       roughness length.                                           *
!                                                                      *
! Output Variables:                                                    *
! UUSTAR wind at height Z normalised by the friction velocity.         *
!                                                                      *
! Local Variables:                                                     *
! TEMP ) intermediate quantities in calculation, stored temporarily.   *
! X    )                                                               *
! X0   )                                                               *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      REAL*4 AA,BB,CC,DD,VK
      PARAMETER (AA = 0.7, BB = 0.75, CC = 5.0, DD = 0.35, VK = 0.4)
!..........Items in argument list.
      REAL*4 RECIPLMO,Z,Z0,UUSTAR
!..........Local variables.
      REAL*4 TEMP,X,X0

!==========Estimate UUSTAR.

      TEMP = (Z + Z0)*RECIPLMO
      X = AA*TEMP + BB*(TEMP - CC/DD)*EXP(-DD*TEMP)
      TEMP = Z0*RECIPLMO
      X0 = AA*TEMP + BB*(TEMP - CC/DD)*EXP(-DD*TEMP)
      UUSTAR = (X - X0 + ALOG((Z + Z0)/Z0))/VK

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE SURFACELAYERUNSTABLE(RECIPLMO,Z,Z0,UUSTAR)
!***********************************************************************
! This subroutine estimates UUSTAR in unstable conditions.             *
!                                                                      *
! Parameters:                                                          *
! VK von Karman's constant.                                            *
!                                                                      *
! Input Variables:                                                     *
! RECIPLMO reciprocal of the Monin-Obukhov length.                     *
! Z        height of surface layer wind.                               *
! Z0       roughness length.                                           *
!                                                                      *
! Output Variables:                                                    *
! UUSTAR wind at height Z normalised by the friction velocity.         *
!                                                                      *
! Local Variables:                                                     *
! TEMP ) intermediate quantities in calculation, stored temporarily.   *
! X    )                                                               *
! X0   )                                                               *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      REAL*4 VK
      PARAMETER (VK = 0.4)
!..........Items in argument list.
      REAL*4 RECIPLMO,Z,Z0,UUSTAR
!..........Local variables.
      REAL*4 TEMP,X,X0

!==========Estimate UUSTAR.

      X = (1.0 - 16.0*(Z + Z0)*RECIPLMO)**0.25
      X0 = (1.0 - 16.0*Z0*RECIPLMO)**0.25
      TEMP = 2.0*(ATAN(X) - ATAN(X0)) -      &
             ALOG((1.0 + X**2)*(1.0 + X)**2/ &
             ((1.0 + X0**2)*(1.0 + X0)**2))
      UUSTAR = (TEMP + ALOG((Z + Z0)/Z0))/VK

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE ROSSBYCOEFFSTABLE(USTAR,RECIPLMO,Z0,ABSF,SIGNF,UGSTAR, DELTAPHI)
!***********************************************************************
! This subroutine estimates UGSTAR and DELTAPHI in stable conditions.  *
!                                                                      *
! Parameters:                                                          *
! PI pi.                                                               *
! VK von Karman's constant.                                            *
!                                                                      *
! Input Variables:                                                     *
! USTAR    friction velocity.                                          *
! RECIPLMO reciprocal of the Monin-Obukhov length.                     *
! Z0       roughness length.                                           *
! ABSF     absolute value of the Coriolis parameter.                   *
! SIGNF    sign of the Coriolis parameter.                             *
!                                                                      *
! Output Variables:                                                    *
! UGSTAR   geostrophic wind speed normalised by the friction velocity. *
! DELTAPHI geostrophic wind direction minus surface wind direction     *
!             (degrees).                                               *
!                                                                      *
! Local Variables:                                                     *
! H     value of boundary layer depth used in calculation.             *
! TEMP  intermediate quantity in calculation, stored temporarily.      *
! UGX ) UGX and -UGY are the components of the geostrophic wind        *
! UGY )    normalised by the friction velocity in coordinates with x   *
!          aligned with the surface stress.                            *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      REAL*4 PI,VK
      PARAMETER (VK = 0.4)
!..........Items in argument list.
      REAL*4 USTAR,RECIPLMO,Z0,ABSF,SIGNF,UGSTAR,DELTAPHI
!..........Local variables.
      REAL*4 H,TEMP,UGX,UGY
!..........Pi.
      PI = 4.0*ATAN(1.0)

!==========Estimate UGSTAR and DELTAPHI.

      TEMP = 2.28*USTAR*RECIPLMO/ABSF
      H = 0.6*USTAR/(ABSF*(SQRT(1.0 + TEMP) + 1.0))
      UGX = (2.2*H*RECIPLMO + ALOG(H/Z0 + 30.0) + 0.19)/VK
      UGY = MAX((3.55*H*RECIPLMO + 1.87),5.14)*SIGNF/VK

      UGSTAR = SQRT(UGX**2 + UGY**2)
      DELTAPHI = ATAN(UGY/UGX)*360.0/(2.0*PI)

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE ROSSBYCOEFFUNSTABLE(USTAR,RECIPLMO,Z0,ABSF,SIGNF, &
                                     UGSTAR,DELTAPHI)
!***********************************************************************
! This subroutine estimates UGSTAR and DELTAPHI in unstable            *
! conditions.                                                          *
!                                                                      *
! Parameters:                                                          *
! PI pi.                                                               *
! VK von Karman's constant.                                            *
!                                                                      *
! Input Variables:                                                     *
! USTAR    friction velocity.                                          *
! RECIPLMO reciprocal of the Monin-Obukhov length.                     *
! Z0       roughness length.                                           *
! ABSF     absolute value of the Coriolis parameter.                   *
! SIGNF    sign of the Coriolis parameter.                             *
!                                                                      *
! Output Variables:                                                    *
! UGSTAR   geostrophic wind speed normalised by the friction velocity. *
! DELTAPHI geostrophic wind direction minus surface wind direction     *
!             (degrees).                                               *
!                                                                      *
! Local Variables:                                                     *
! AA )  Rossby number similarity coefficients.                         *
! BB )                                                                 *
! RMU   VK*USTAR*RECIPLMO/ABSF                                         *
! UGX ) UGX and -UGY are the components of the geostrophic wind        *
! UGY )    normalised by the friction velocity in coordinates with x   *
!          aligned with the surface stress.                            *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      REAL*4 PI,VK
      PARAMETER (VK = 0.4)
!..........Items in argument list.
      REAL*4 USTAR,RECIPLMO,Z0,ABSF,SIGNF,UGSTAR,DELTAPHI
!..........Local variables.
      REAL*4 AA,BB,RMU,UGX,UGY
!..........Pi.
      PI = 4.0*ATAN(1.0)

!==========Estimate UGSTAR and DELTAPHI.

      RMU = VK*USTAR*RECIPLMO/ABSF
      IF (RMU.GT.-50.0) THEN
        AA = 1.01 + RMU*(-0.105 + RMU*(-9.9E-4 + RMU*8.1E-7))
        BB = 5.14 + RMU*(0.142 + RMU*(1.17E-3 - RMU*3.3E-6))
      ELSE
        AA = 3.69
        BB = 1.38
      ENDIF
      UGX = (ALOG(USTAR/(Z0*ABSF) + 100.0) - AA)/VK
      UGY = BB*SIGNF/VK

      UGSTAR = SQRT(UGX**2 + UGY**2)
      DELTAPHI = ATAN(UGY/UGX)*360.0/(2.0*PI)

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE SCALC(RLAT,TDAY,THOUR,S)
!***********************************************************************
! This subroutine estimates S.                                         *
!                                                                      *
! Parameters:                                                          *
! PI pi.                                                               *
!                                                                      *
! Input Variables:                                                     *
! RLAT  latitude of source (degrees, being positive in the northern    *
!          hemisphere).                                                *
! TDAY  Julian day number.                                             *
! THOUR local time (hours).                                            *
!                                                                      *
! Output Variables:                                                    *
! S sine of the solar elevation.                                       *
!                                                                      *
! Local Variables:                                                     *
! A solar declination (radians).                                       *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Parameters.
      REAL*4 PI
!..........Items in argument list.
      REAL*4 RLAT,TDAY,THOUR,S
!..........Local variables.
      REAL*4 A
!..........Pi.
      PI = 4.0*ATAN(1.0)

!==========Estimate S.

      A = 23.45*(2.0*PI/360.0)*SIN((TDAY + 284.0)*2.0*PI/365.0)
      S = SIN(2.0*PI*RLAT/360.0)*SIN(A) +                               &
          COS(2.0*PI*RLAT/360.0)*COS(A)*COS((THOUR - 12.0)*2.0*PI/24.0)

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE FTheta0DayTimeCalc(T0K, Cl, R, Alpha, RK, FTheta0)
!-----------------------------------------------------------------------
! This subroutine estimates FTheta0 in during the daytime.
!
! Input variables:
! T0K: Near surface temperature (degrees K).
! Cl: Cloud amount (oktas).
! R: Surface albedo.
! Alpha: Alpha.
! RK: Incoming solar radiation.
!
! Output variables:
! FTheta0: Surface sensible heat flux.
!
! Local variables:
! RN: Net radiation.
! S: (lamda/CP)d(qs)/dT evaluated at T = T0K, where lamda is the specific
!    latent heat of vaporization of water and qs is the saturated specific
!    humidity at temperature T.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      REAL*4 T0K, Cl, R, Alpha, RK
!..........Output variables.
      REAL*4 FTheta0
!..........Local variables.
      REAL*4 RN, S

!==========Estimate FTheta0.

      RN = ((1.0 - R)*RK + 5.31E-13*T0K**6 - 5.67E-8*T0K**4 + &
           60.0*Cl/8.0)/1.12
      S = EXP(0.055*(T0K - 279.0))
      FTheta0 = 0.9*RN*((1.0 - Alpha)*S + 1.0)/(S + 1.0) - 20.0*Alpha

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE KCalc(Cl,S,RK)
!-----------------------------------------------------------------------
! This subroutine estimates RK in during the daytime.
!
! Input variables:
! Cl: Cloud amount (oktas).
! S: Sine of the solar elevation.
!
! Output variables:
! RK: Incoming solar radiation.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      REAL*4 Cl,S
!..........Output variables.
      REAL*4 RK

!==========Estimate RK.

      RK = (990.0*S - 30.0)*(1.0 - 0.75*(Cl/8.0)**3.4)

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE SIGMATHETACALC(TS,USTAR,RECIPLMO,Z0,SIGMATHETA)
!***********************************************************************
! This subroutine estimates SIGMATHETA.                                *
!                                                                      *
! Input Variables:                                                     *
! TS       sampling time (hours).                                      *
! USTAR    friction velocity.                                          *
! RECIPLMO reciprocal of the Monin-Obukhov length.                     *
! Z0       roughness length.                                           *
!                                                                      *
! Output Variables:                                                    *
! SIGMATHETA standard deviation of the mean wind direction.            *
!                                                                      *
! Local Variables:                                                     *
! UUSTAR wind at height Z normalised by the friction velocity.         *
! U10 wind speed at 10m.                                               *
!***********************************************************************

!==========Declarations.

      IMPLICIT NONE
!..........Items in argument list.
      REAL*4 TS,USTAR,RECIPLMO,Z0,SIGMATHETA
!..........Local variables.
      REAL*4 UUSTAR,U10

!==========Estimate SIGMATHETA.

      If (UStar == 0.0) Then
        U10 = 0.5
      Else
        IF (RECIPLMO.GE.0.0) THEN
          CALL SURFACELAYERSTABLE(RECIPLMO,10.0,Z0,UUSTAR)
        ELSE
          CALL SURFACELAYERUNSTABLE(RECIPLMO,10.0,Z0,UUSTAR)
        ENDIF
        U10 = USTAR*UUSTAR
      End If
      SIGMATHETA = 0.065*SQRT(7.0*TS/U10)

      RETURN
      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE Q0Calc(T0K, RH0, Q0)
!-----------------------------------------------------------------------
! This subroutine estimates Q0.
!
! Parameters:
! P: Value of pressure assumed (millibars).
!
! Input variables:
! T0K: T0K.
! RH0: RH0.
!
! Output variables:
! Q0: Specific humidity.
!
! Local variables:
! P_kPa: Pressure (kPa).
! SatVP: Saturated vapour pressure (kPa).
! SatMR: Saturated mixing ratio.
! MR: Mixing ratio.
!
! Functions:
! SatVapourPressure: Saturated vapour pressure (kPa).
! MixingRatio: Mixing ratio.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Parameters.
      REAL*4 P
      PARAMETER (P = 1013.25)
!..........Input variables.
      REAL*4 RH0, T0K
!..........Ouput variables.
      REAL*4 Q0
!..........Local variables.
      REAL*4 P_kPa, SatVP, SatMR, MR
!..........Functions.
  !    REAL*4 SatVapourPressure, MixingRatio   ! Not needed now in module
  !    EXTERNAL SatVapourPressure, MixingRatio !

!==========Calculate specific humidity.

      P_kPa = P/10.0
      SatVP = SatVapourPressure(T0K)
      SatMR = MixingRatio(P_kPa, SatVP)
      MR = SatMR*RH0/100.0
      Q0 = MR/(1.0 + MR)

      END Subroutine

!----------------------------------------------------------------------!

      SUBROUTINE LambdaECalc(FTheta0, Alpha, T0K, LambdaE)
!-----------------------------------------------------------------------
! This subroutine estimates LambdaE.
!
! Input variables:
! FTheta0: FTheta0.
! Alpha: Alpha.
! T0K: T0K.
!
! Output variables:
! LambdaE: LambdaE.
!
! Local variables:
! S: (lamda/CP)d(qs)/dT evaluated at T = T0K, where lamda is the specific
!    latent heat of vaporization of water and qs is the saturated specific
!    humidity at temperature T.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!..........Input variables.
      REAL*4 FTheta0, Alpha, T0K
!..........Ouput variables.
      REAL*4 LambdaE
!..........Local variables.
      REAL*4 S

!==========Estimate latent heat flux.

      IF (FTheta0.GT.0.0) THEN
        S = EXP(0.055*(T0K - 279.0))
        LambdaE = (FTheta0 + 20.0*Alpha)*Alpha*S/((1.0 - Alpha)*S + &
                  1.0) + 20.0*Alpha
      ELSE
        LambdaE = 0.0
      ENDIF

      END Subroutine

!----------------------------------------------------------------------!

      REAL FUNCTION SatVapourPressure(T)
!________________________________________________________________________
!
!  Routine:   SatVapourPressure
!  Purpose:   Calculates the saturation vapour pressure in kPa at temperature T
!             using Wexler's formula (1976) (from Bolton 1979 Mon Wea Rev)
!  Passed:    T temperature (K)
!  Returns:
!  Remarks:
!  Author:    Steph Dyster, 04/06/98
!  Revisions: --/--/-- ***:
!  (C) CERC 1998
!________________________________________________________________________
      IMPLICIT NONE
!.....Arguments
      REAL*4        T
!.....Local variables
      REAL*4    G(8),SUMGT
      INTEGER*4        I
!.....Functions
!.....Constants
      DATA G/-2.9912729E3,-6.0170128E3,18.87643854,-2.8354721E-2, &
      1.7838301E-5,-8.4150417E-10,4.4412543E-13,2.858487 /
!.....Debug/Log

!D     CALL Log_Screen('Entered SatVapourPressure routine')
!D     CALL Set_Routine('SatVapourPressure')
      SUMGT=0.
      DO I=1,7
          SUMGT=SUMGT+( G(I)*(T**(I-3)) )
      END DO
      SUMGT=SUMGT+G(8)*ALOG(T)
      IF (SUMGT.GT.85.) SUMGT=85.

      SatVapourPressure=EXP(SUMGT)/1000.


!D     CALL Clear_Routine()
      RETURN
      END Function



      REAL FUNCTION MixingRatio(P,VapP)
!________________________________________________________________________
!
!  Routine:   MixingRatio
!  Purpose:   Calculates the mixing ratio from pressure and vapour pressure
!  Passed:    P = pressure; VapP = vapour pressure
!  Returns:
!  Remarks:
!  Author:    Steph Dyster, 04/06/98
!  Revisions: --/--/-- ***:
!  (C) CERC 1998
!________________________________________________________________________
      IMPLICIT NONE
!.....Arguments
      REAL*4        P, VapP
!.....Local variables
      REAL*4        eps
!DJT      INCLUDE 'ADMScnst.inc'
!.....Functions
!.....Constants
!.....Debug/Log

!D     CALL Log_Screen('Entered MixingRatio routine')
!D     CALL Set_Routine('MixingRatio')

       eps = 0.62197
!DJT       eps = MW_water/MW_air
       MixingRatio = eps*VapP/(P - VapP)

!D     CALL Clear_Routine()
      RETURN
      END Function

End Module SingleSiteModule

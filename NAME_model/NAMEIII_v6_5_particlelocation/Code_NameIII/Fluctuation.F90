! Module: Fluctuation Module

! Notes:  draws on code from the ADMS 3.0 fluctuations module by Dave Thomson
!         - the code which was used is that in the file Fluct.For, as sent to
!         CERC for incorporation in ADMS 3.0 but before the final changes made by CERC
!         in integrating the code into ADMS 3.0

Module FluctuationModule

! This module provides code for estimating statistics of concentration fluctuations.

! Module structure:
!
! CalculatePlumeSigmas
!
! CalculateSigmaC---(C2BarCalc---MuCalc
!                   (CCBarCalc---MuCalc
!
! CalculatePdf---(                         (CRatioCalcNeg
!                (GammaCalc---CRatioCalc---(CRatioCalcPos
!                (                         (CRatioCalcInt---ErrorFunc
!                (
!                (            (SigmaCalcNeg
!                (SigmaCalc---(SigmaCalcPos
!                (            (SigmaCalcInt---ErrorFunc
!                (
!                (DoseCalc---DoseContribCalc
!                (CHatCalc
!                (ProbCalc---ErrorFunc
!                (ConcCalc---Std_Normal_Inv

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private

! Items from this module which need to be made available outside the module:
Public  :: CalculatePlumeSigmas  ! Calculate estimates of sigma_y and sigma_z.
Public  :: CalculateSigmaC       ! Calculate estimate of sigma_c.
Public  :: CalculatePdf          ! Calculate estimates of various pdf statistics.

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine CalculatePlumeSigmas(DualTimes, InterpCoeff,  &
                                SXX1, SXY1, SYY1, SZZ1,  &
                                SXX2, SXY2, SYY2, SZZ2,  &
                                U1, V1, U2, V2,          &
                                SigmaY, SigmaZ)
! Calculates estimates of the standard deviations SigmaY and SigmaZ of particle
! displacements at a given travel time. The values are calculated by processing
! X-Stats information into plume-oriented coordinate frames.

  Implicit None
  ! Argument list:
  Logical,   Intent(In)  :: DualTimes   ! Indicates that values of SigmaY and SigmaZ
                                        ! are calculated for two consecutive travel
                                        ! times in the travel time grid (with linear
                                        ! interpolation used between these points).
                                        ! Otherwise, the value of SigmaY and SigmaZ
                                        ! is calculated for a single travel time in
                                        ! the travel time grid.
  Real(Std), Intent(In)  :: InterpCoeff ! Interpolation coefficient used when values
                                        ! are interpolated between two travel times.
  Real(Std), Intent(In)  :: SXX1        !} Entries of the particle displacement
  Real(Std), Intent(In)  :: SXY1        !} covariance matrix (in default Cartesian
  Real(Std), Intent(In)  :: SYY1        !} coord system) for earlier grid time s1.
  Real(Std), Intent(In)  :: SZZ1        !}
  Real(Std), Intent(In)  :: SXX2        !] Entries of the particle displacement
  Real(Std), Intent(In)  :: SXY2        !] covariance matrix (in default Cartesian
  Real(Std), Intent(In)  :: SYY2        !] coord system) for later grid time s2
  Real(Std), Intent(In)  :: SZZ2        !] (processed only if DualTimes is true).
  Real(Std), Intent(In)  :: U1          !} Horizontal components of mean wind vector
  Real(Std), Intent(In)  :: V1          !} at plume centre for travel time s1.
  Real(Std), Intent(In)  :: U2          !] Horizontal components of mean wind vector
  Real(Std), Intent(In)  :: V2          !] at plume centre for travel time s2.
  Real(Std), Intent(Out) :: SigmaY      ! Calculated SigmaY value.
  Real(Std), Intent(Out) :: SigmaZ      ! Calculated SigmaZ value.

  ! Locals:
  Real(Std) :: HorzWind      ! Modulus of mean horizontal wind vector.
  Real(Std) :: SinPhi        ! Sine of rotation angle in SigmaY calculation.
  Real(Std) :: CosPhi        ! Cosine of rotation angle in SigmaY calculation.
  Real(Std) :: SigmaYSquared ! Square of calculated SigmaY values.
  Real(Std) :: SigmaY1       !} Calculated SigmaY and SigmaZ for travel time s1.
  Real(Std) :: SigmaZ1       !}
  Real(Std) :: SigmaY2       !] Calculated SigmaY and SigmaZ for travel time s2.
  Real(Std) :: SigmaZ2       !]

  ! Process X-Stats information for travel time s1.
  HorzWind = Sqrt(U1**2 + V1**2)
  ! Establish that the wind direction is well-defined at the plume centroid, and
  ! write a warning if mean wind speed is less than 1.0 m/s.
  If ((U1 == 0.0).and.(V1 == 0.0)) Then
    Call Message('Error in CalculatePlumeSigmas: undefined wind direction', 3)
  Else If (HorzWind < 1.0) Then
    Call Message('Warning: calculation of SigmaY may be unreliable in light winds', 1)
  End If
  ! Calculate entries in rotation matrix.
  SinPhi = V1 / HorzWind
  CosPhi = U1 / HorzWind
  ! Apply rotation of coordinate frame to the input covariance matrix
  ! (note that we only need to calculate the SigmaY entry here!).
  SigmaYSquared =   CosPhi * (CosPhi*SYY1 - SinPhi*SXY1)  &
                  - SinPhi * (CosPhi*SXY1 - SinPhi*SXX1)
  SigmaY1       = Sqrt(SigmaYSquared)
  ! Calculate SigmaZ value.
  SigmaZ1 = Sqrt(SZZ1)

  ! Process X-Stats information for travel time s2 (if applicable).
  If (DualTimes) Then
    HorzWind = Sqrt(U2**2 + V2**2)
    ! Establish that the wind direction is well-defined at the plume centroid, and
    ! write a warning if mean wind speed is less than 1.0 m/s.
    If ((U2 == 0.0).and.(V2 == 0.0)) Then
      Call Message('Error in CalculatePlumeSigmas: undefined wind direction', 3)
    Else If (HorzWind < 1.0) Then
      Call Message('Warning: calculation of SigmaY may be unreliable in light winds', 1)
    End If
    ! Calculate entries in rotation matrix.
    SinPhi = V2 / HorzWind
    CosPhi = U2 / HorzWind
    ! Apply rotation of coordinate frame to the input covariance matrix
    ! (note that we only need to calculate the SigmaY entry here!).
    SigmaYSquared =   CosPhi * (CosPhi*SYY2 - SinPhi*SXY2)  &
                    - SinPhi * (CosPhi*SXY2 - SinPhi*SXX2)
    SigmaY2       = Sqrt(SigmaYSquared)
    ! Calculate SigmaZ value.
    SigmaZ2 = Sqrt(SZZ2)
    ! Interpolate between s1 and s2 values.
    SigmaY = (1.0 - InterpCoeff)*SigmaY1 + InterpCoeff*SigmaY2
    SigmaZ = (1.0 - InterpCoeff)*SigmaZ1 + InterpCoeff*SigmaZ2
  Else
    ! Use s1 values.
    SigmaY = SigmaY1
    SigmaZ = SigmaZ1
  End If

End Subroutine CalculatePlumeSigmas

!-------------------------------------------------------------------------------------------------------------

Subroutine CalculateSigmaC(CalcType,                       &
                           TravelTime, YOffset, ZOffset,   &
                           MeanConc, SigmaY, SigmaZ,       &
                           U, SigmaVel2, Epsilon,          &
                           SourceDiameter, ReleaseRate,    &
                           SigmaPlumeRise2,                &
                           TAverage, AveragingTime,        &
                           ReleaseTime,                    &
                           FluctsAdvectVel,                &
                           SigmaC)
! Calculates an estimate of the standard deviation of the concentration fluctuations.
! Currently restricted to a single continuous point source only.

  Implicit None
  ! Argument list:
  Integer,   Intent(In)  :: CalcType
  Real(Std), Intent(In)  :: TravelTime
  Real(Std), Intent(In)  :: YOffset
  Real(Std), Intent(In)  :: ZOffset
  Real(Std), Intent(In)  :: MeanConc
  Real(Std), Intent(In)  :: SigmaY
  Real(Std), Intent(In)  :: SigmaZ
  Real(Std), Intent(In)  :: U
  Real(Std), Intent(In)  :: SigmaVel2
  Real(Std), Intent(In)  :: Epsilon
  Real(Std), Intent(In)  :: SourceDiameter
  Real(Std), Intent(In)  :: ReleaseRate
  Real(Std), Intent(In)  :: SigmaPlumeRise2
  Logical,   Intent(In)  :: TAverage
  Real(Std), Intent(In)  :: AveragingTime
  Real(Std), Intent(In)  :: ReleaseTime
  Real(Std), Intent(In)  :: FluctsAdvectVel
  Real(Std), Intent(Out) :: SigmaC
  ! CalcType        :: Indicates type of fluctuations calculation to be carried out
  !                    (currently must be set to 1 - continuous releases only).
  ! TravelTime      :: Travel time from source to receptor grid point (actually an
  !                    average over all particles in the receptor grid-box).
  ! YOffset         :: Cross-wind offset of grid point relative to the centre-line
  !                    of the plume (in metres).
  ! ZOffset         :: Vertical offset of grid point relative to the mean plume height
  !                    (in metres). $$ not needed at present?
  ! MeanConc        :: Mean concentration at receptor grid point.
  ! SigmaY          :} Lateral and vertical spread of the plume at the given travel
  ! SigmaZ          :} time (this includes contribution from source size).
  ! U               :: Mean wind speed evaluated at the plume's centre of mass P.
  ! SigmaVel2       :: Average of the three velocity variance components at P.
  ! Epsilon         :: Turbulent kinetic energy dissipation rate at P.
  ! SourceDiameter  :: Effective diameter of point source.
  ! ReleaseRate     :: Source release rate.
  ! SigmaPlumeRise2 :: Variance of the plume spread due to plume rise
  !                    (currently not supported - set to zero).
  ! TAverage        :: Indicates time-averaged concentrations are considered
  !                    (defined for CalcType = 1 case only).
  ! AveragingTime   :: Averaging time in seconds used for time-averaged concentrations
  !                    (defined for CalcType = 1 case if TAverage is true).
  ! ReleaseTime     :: Duration of release in seconds (defined for CalcType = 2 only).
  ! FluctsAdvectVel :: Effective advection velocity for the fluctuations.
  ! SigmaC          :: Calculated standard deviation of the fluctuations.

  ! Locals:
  Real(Std) :: SmoothingTime
  Real(Std) :: XAv
  Real(Std) :: Sigma0X2
  Real(Std) :: Sigma0Y2
  Real(Std) :: Sigma0Z2
  Real(Std) :: Sigma1X2
  Real(Std) :: Sigma1Y2
  Real(Std) :: Sigma1Z2
  Real(Std) :: SigmaYD2
  Real(Std) :: MaxConcYZ
  Real(Std) :: MaxConcZ
  Real(Std) :: HX
  Real(Std) :: HY
  Real(Std) :: HZ
  Real(Std) :: C2Bar
  Real(Std) :: SigmaDeltaX2
  Real(Std) :: SigmaDeltaY2
  Real(Std) :: SigmaDeltaZ2
  Real(Std) :: SigmaC2
  ! SmoothingTime :: Smoothing time in calculation of time-averaged concentrations:
  !                    SmoothingTime is {AveragingTime for CalcType = 1,
  !                                     {ReleaseTime   for CalcType = 2.
  ! XAv           :: Averaging length appropriate for time-averaged concentrations.
  ! Sigma0X2      :} Variance components of the source's spatial distribution.
  ! Sigma0Y2      :}
  ! Sigma0Z2      :}
  ! Sigma1X2      :] Variance components of the displacement of a single particle.
  ! Sigma1Y2      :]
  ! Sigma1Z2      :]
  ! SigmaYD2      :: Variance of the displacement of a single particle in the
  !                  cross-wind direction used in the SigmaDelta calculation.
  ! MaxConcYZ     :: Estimate of the maximum plume concentration for the given
  !                  travel time (spatial peak over Y and Z).
  ! MaxConcZ      :: Estimate of the maximum vertical concentration at the receptor
  !                  grid point (spatial peak over Z).
  ! HX            :} Dilution factors of mean concentration field.
  ! HY            :}
  ! HZ            :}
  ! C2Bar         :: Calculated mean-square concentration.
  ! SigmaDeltaX2  :} Calculated variance components of the change in separation
  ! SigmaDeltaY2  :} of a particle pair which are initially close together.
  ! SigmaDeltaZ2  :}
  ! SigmaC2       :: Estimated concentration variance.

  ! Consistency checks on the type of calculation.
  ! $$ currently only CalcType = 1 case is supported
  If (CalcType /= 1) Then
    Call Message('Error in CalculateSigmaC: calculation type not currently supported', 3)
  End If

  ! Check for zero mean concentration (in this case SigmaC is also zero).
  If (MeanConc == 0.0) Then
    SigmaC = 0.0
    Return
  End If

  ! $$ Need to consider the Sigma0 definition and the calculation of Sigma1
  ! w.r.t. source size contribution. Also take care with the subtraction here
  ! since it could lead to negative values - check needed ?

  ! Sigma0 terms: variance components of the source's spatial distribution.
  Sigma0Y2 = SourceDiameter**2
  Sigma0Z2 = SourceDiameter**2

  ! Sigma1 terms: variance components of the displacement of a single particle.
  Sigma1Y2 = SigmaY**2 - Sigma0Y2
  Sigma1Z2 = SigmaZ**2 - Sigma0Z2

  ! Calculate estimate of maximum plume concentration (spatial peak over Y and Z).
  MaxConcYZ = ReleaseRate / (2*Pi*U*SigmaY*SigmaZ)

  ! Calculate estimate of maximum vertical concentration at the receptor grid point
  ! (assuming the plume is Gaussian in the lateral direction).
  MaxConcZ = MaxConcYZ * Exp(- YOffset**2 / (2*SigmaY**2))

  ! Adjust these estimates of the maximum concentrations to be consistent with the
  ! receptor mean concentration.
  MaxConcZ  = Max(MaxConcZ, MeanConc)
  MaxConcYZ = Max(MaxConcYZ, MaxConcZ)

  ! Calculate the dilution factors of the mean concentration field.
  HY = MaxConcZ / MaxConcYZ
  HZ = MeanConc / MaxConcZ

  ! Set smoothing time to be used in time-averaging calculations.
  ! $$ note that CalcType = 2 case *is* considered here for use in the future.
  If (CalcType == 1) Then
    ! Set smoothing time to the averaging time (for time-averaged concentrations),
    ! or zero otherwise.
    If (TAverage) Then
      SmoothingTime = AveragingTime
    Else
      SmoothingTime = 0.0
    End If
  Else If (CalcType == 2) Then
    ! Set smoothing time equal to the release duration.
    SmoothingTime = ReleaseTime
  End If

  ! Set averaging length scale appropriate for time-averaged concentrations.
  XAv = SmoothingTime * FluctsAdvectVel

  ! Set the value to be used for Sigma1Y2 in the SigmaDelta calculation.
  ! $$ currently we use the same value of Sigma1Y2 as elsewhere (i.e. there is no
  ! reduction for unresolved mesoscale motions as in ADMS scheme) - consider this issue in the future ??
  SigmaYD2 = Sigma1Y2

  ! Call C2BarCalc to calculate an estimate of the mean-square concentration.
  ! Note that the subroutine also computes the three components of SigmaDelta2
  ! (which are not used again in any single-source calculation but would be
  ! needed for any multiple-source calculations).
  Call C2BarCalc(CalcType, SigmaVel2, Epsilon, TravelTime, XAv,           &
                 Sigma0X2, Sigma0Y2, Sigma0Z2, MaxConcYZ, HX, HY, HZ,     &
                 Sigma1X2, Sigma1Y2, Sigma1Z2, SigmaYD2, SigmaPlumeRise2, &
                 C2Bar, SigmaDeltaX2, SigmaDeltaY2, SigmaDeltaZ2)

  ! Calculate the standard deviation of the fluctuations from the estimated C2Bar.
  SigmaC2 = Max(0.0, C2Bar - MeanConc**2)
  SigmaC  = Sqrt(SigmaC2)

End Subroutine CalculateSigmaC

!-------------------------------------------------------------------------------------------------------------

Subroutine CalculatePdf(PdfType, FixedThresholds, PdfSize, PdfThresholds, &
                        MeanC, SigmaC, PdfScale, PdfValues)
! Calculates estimates of various statistics of the concentration pdf.

  Implicit None
  ! Argument list:
  Integer,   Intent(In)            :: PdfType
  Logical,   Intent(In)            :: FixedThresholds
  Integer,   Intent(In)            :: PdfSize
  Real(Std), Intent(In), Optional  :: PdfThresholds(MaxPdfSize)
  Real(Std), Intent(In)            :: MeanC
  Real(Std), Intent(In)            :: SigmaC
  Integer,   Intent(Out), Optional :: PdfScale
  Real(Std), Intent(Out)           :: PdfValues(:)
  ! PdfType         :: Type of pdf calculation: 1 = exceedence probs, 2 = percentiles
  ! FixedThresholds :: Calculation uses fixed concentration or probability thresholds
  ! PdfSize         :: Number of points in the computed pdf
  ! PdfThresholds   :: Fixed concentration or probability thresholds (optional)
  ! MeanC           :: Mean of the distribution
  ! SigmaC          :: Standard deviation of the distribution
  ! PdfScale        :: Scale factor for automated calculation of exceedence probs
  ! PdfValues       :: Calculated values of the exceedence probs or percentiles

  ! Locals:
  Real(Std) :: CRatio !  Mean square concentration divided by mean concentration squared
  Real(Std) :: Gamma  !} Parameters of the clipped normal distribution
  Real(Std) :: Sigma  !}
  Real(Std) :: Thresholds(MaxPdfSize)

  ! Check that mean and standard deviation of the concentration are positive.
# ifdef ExtraChecks
  If (MeanC <= 0.0) Then
    Call Message('Error in CalculatePdf: mean concentration is not positive', 3)
  End If
  If (SigmaC <= 0.0) Then
    Call Message('Error in CalculatePdf: standard deviation is not positive', 3)
  End If
# endif

  ! Calculate the parameters Gamma and Sigma of the clipped normal distribution.
  CRatio = 1.0 + (SigmaC/MeanC)**2
  Call GammaCalc(CRatio, Gamma)
  Call SigmaCalc(MeanC, Gamma, Sigma)

  If (PdfType == 1) Then
    ! Calculate exceedence probabilities
    If (FixedThresholds) Then
      ! Use fixed concentration thresholds (prescribed by the user)
      Thresholds = PdfThresholds
    Else
      ! Auto-mode calculation of concentration thresholds
      Call CHatCalc(Gamma, Sigma, PdfSize, Thresholds, PdfScale)
    End If
    Call ProbCalc(PdfSize, Gamma, Sigma, Thresholds, PdfValues)
  Else If (PdfType == 2) Then
    ! Calculate percentiles
    Call ConcCalc(PdfSize, Gamma, Sigma, Thresholds, PdfValues)
  Else
    Call Message('Error in CalculatePdf: invalid PdfType', 3)
  End If

End Subroutine CalculatePdf

!-------------------------------------------------------------------------------------------------------------

!#### ADMS CODE STARTS
!## Note that some changes to the format of the code have been necessary to run it
!## under Fortran 90 (comment lines, continuation lines, End statements, etc.).
!## Locally declared parameters (e.g. Pi) have been removed and the values of
!## global parameters used in their place.
!## Any material changes to the ADMS code are noted at the start of each routine.

!
! ######################################################################
!
      SUBROUTINE C2BarCalc(ICalcType,SigmaVel2,Epsilon,TMS,XAv,          &
                           Sigma0X2,Sigma0Y2,Sigma0Z2,                   &
                           DBarMax,HX,HY,HZ,                             &
                           Sigma1X2,Sigma1Y2,Sigma1Z2,SigmaYD2,SigmaPR2, &
                           D2Bar,                                        &
                           SigmaDeltaX2,SigmaDeltaY2,SigmaDeltaZ2)
!***********************************************************************
! This subroutine calculates the mean-square concentration.            *
!                                                                      *
! Changes to ADMS code:                                                *
! - the arguments DBarMax, HX, HY, HZ, D2Bar changed from type REAL*8  *
!   to type REAL*4 in order to run consistently at standard precision. *
!                                                                      *
! Parameters:                                                          *
! CKol Kolmogorov constant in the Eulerian velocity structure          *
!         function.                                                    *
! A0 ) Tunable constants.                                              *
! A4 )                                                                 *
! A5 )                                                                 *
!                                                                      *
! Input variables:                                                     *
! ICalcType    Indicates the calculation type (1 denotes a calculation *
!                 for a continuous release, 2 denotes the calculation  *
!                 of time integrated quantities for a finite duration  *
!                 release, and 3 denotes a calculation of quantities   *
!                 as a function of time for a finite duration          *
!                 release).                                            *
! SigmaVel2    Average of the three velocity component variances.      *
! Epsilon      Dissipation rate.                                       *
! TMS          Travel time, t - s.                                     *
! XAv          Averaging length.                                       *
! Sigma0X/Y/Z2 Components of the variance of the source's spatial      *
!                 distribution.                                        *
! DBarMax      ICalcType = 1: maximum value of the mean concentration  *
!                 at downwind distance where output is required (the   *
!                 maximum is over Y and Z); ICalcType = 2: maximum     *
!                 value of the time-integral of the mean concentration *
!                 at downwind distance where output is required (the   *
!                 maximum is over Y and Z); ICalcType = 3: maximum     *
!                 value of the mean concentration at downwind distance *
!                 for which output is required (the maximum is over Y, *
!                 Z and T).                                            *
! HX/Y/Z       Dilution factors.                                       *
! Sigma1X/Y/Z2 Components of the variance of the displacement of a     *
!                 single particle.                                     *
! SigmaYD2     Variance of the displacement of a single particle in    *
!                 Y direction as used to calculate SigmaDelta.         *
! SigmaPR2     Variance of the plume spread due to plume rise.         *
!                                                                      *
! Output variables:                                                    *
! D2Bar            Mean-square instantaneous concentration.            *
! SigmaDeltaX/Y/Z2 These are defined in the documentation.             *
!                                                                      *
! Local variables:                                                     *
! R2                    Epsilon*TMS**3/3.0 plus time averaging terms.  *
! RMu                   Mu.                                            *
! SigmaSigmaX/Y/Z2    ) Apart from the _02, these are defined in the   *
! Sigma1X/Y/Z2_02     )    documentation; the _02 indicates the        *
! SigmaDeltaX/Y/Z2_02 )    addition of Sigma0X/Y/Z2.                   *
! SigmaSigmaX/Y/Z2_02 )                                                *
! ChiX/Y/Z,GX/Y/Z     )                                                *
!***********************************************************************
      IMPLICIT NONE
!..........Parameters.
      REAL*4 A0,A4,A5,CKol
      PARAMETER (A0 = 0.2, A4 = 0.01, A5 = 0.02, CKol = 2.0)       ! Tune
!..........Variables in argument list.
      INTEGER*4 ICalcType
      REAL*4 SigmaVel2,Epsilon,TMS,XAv,                     &
             Sigma0X2,Sigma0Y2,Sigma0Z2,                    &
             Sigma1X2,Sigma1Y2,Sigma1Z2,SigmaYD2,SigmaPR2,  &
             SigmaDeltaX2,SigmaDeltaY2,SigmaDeltaZ2
      REAL*4 DBarMax,HX,HY,HZ,D2Bar
!..........Local variables.
      REAL*4 R2,RMu,                                          &
             SigmaSigmaX2,SigmaSigmaY2,SigmaSigmaZ2,          &
             Sigma1X2_02,Sigma1Y2_02,Sigma1Z2_02,             &
             SigmaDeltaX2_02,SigmaDeltaY2_02,SigmaDeltaZ2_02, &
             SigmaSigmaX2_02,SigmaSigmaY2_02,SigmaSigmaZ2_02, &
             ChiX,ChiY,ChiZ,GX,GY,GZ

!..........Sigmas.
      R2 = Epsilon*TMS**3/3.0 + SigmaPR2 +                    &
           (4.0*CKol/27.0)*(Epsilon*XAv)**(2.0/3.0)*TMS**2 +  &
           (A4*SigmaVel2*TMS**2 + A5*Epsilon*TMS**3)*         &
           (XAv*Epsilon)**2/SigmaVel2**3
      SigmaDeltaY2 = R2*SigmaYD2/ &
                     (R2 + A0*SQRT(R2*SigmaYD2) + SigmaYD2)
      SigmaSigmaY2 = 2.0*Sigma1Y2 - SigmaDeltaY2
      Sigma1Y2_02 = Sigma1Y2 + Sigma0Y2
      SigmaDeltaY2_02 = SigmaDeltaY2 + Sigma0Y2
      SigmaSigmaY2_02 = SigmaSigmaY2 + Sigma0Y2
      SigmaDeltaZ2 = R2*Sigma1Z2/ &
                     (R2 + A0*SQRT(R2*Sigma1Z2) + Sigma1Z2)
      SigmaSigmaZ2 = 2.0*Sigma1Z2 - SigmaDeltaZ2
      Sigma1Z2_02 = Sigma1Z2 + Sigma0Z2
      SigmaDeltaZ2_02 = SigmaDeltaZ2 + Sigma0Z2
      SigmaSigmaZ2_02 = SigmaSigmaZ2 + Sigma0Z2
      IF (ICalcType.EQ.3) THEN
        SigmaDeltaX2 = R2*Sigma1X2/ &
                       (R2 + A0*SQRT(R2*Sigma1X2) + Sigma1X2)
        SigmaSigmaX2 = 2.0*Sigma1X2 - SigmaDeltaX2
        Sigma1X2_02 = Sigma1X2 + Sigma0X2
        SigmaDeltaX2_02 = SigmaDeltaX2 + Sigma0X2
        SigmaSigmaX2_02 = SigmaSigmaX2 + Sigma0X2
      ENDIF

!..........Chi.
      ChiY = 2.0*Sigma1Y2_02/SigmaSigmaY2_02
      ChiZ = 2.0*Sigma1Z2_02/SigmaSigmaZ2_02
      IF (ICalcType.EQ.3) ChiX = 2.0*Sigma1X2_02/SigmaSigmaX2_02

!..........G.
      GY = Sigma1Y2_02/SQRT(SigmaSigmaY2_02*SigmaDeltaY2_02)
      GZ = Sigma1Z2_02/SQRT(SigmaSigmaZ2_02*SigmaDeltaZ2_02)
      IF (ICalcType.EQ.3) GX = Sigma1X2_02/ &
                          SQRT(SigmaSigmaX2_02*SigmaDeltaX2_02)

!..........RMu.
      CALL MuCalc(ICalcType,SigmaVel2,Epsilon,TMS,XAv,    &
                  Sigma0X2,Sigma0Y2,Sigma0Z2,             &
                  SigmaDeltaX2,SigmaDeltaY2,SigmaDeltaZ2, &
                  RMu)

!..........D2Bar.
      D2Bar = RMu*(DBarMax**2)*GY*GZ
      IF (ICalcType.EQ.3) D2Bar = D2Bar*GX
      D2Bar = D2Bar*(HY**ChiY)*(HZ**ChiZ)
      IF (ICalcType.EQ.3) D2Bar = D2Bar*(HX**ChiX)

      END SUBROUTINE C2BarCalc
!
! ######################################################################
!
      SUBROUTINE CCBarCalc(SigmaVel2_1,Epsilon_1,TMS_1,XAv_1, &
                           Sigma0Y2_1,Sigma0Z2_1,             &
                           DBarMax_1,HY_1,HZ_1,               &
                           Sigma1Y2_1,Sigma1Z2_1,             &
                           SigmaDeltaY2_1,SigmaDeltaZ2_1,     &
                           SigmaVel2_2,Epsilon_2,TMS_2,XAv_2, &
                           Sigma0Y2_2,Sigma0Z2_2,             &
                           DBarMax_2,HY_2,HZ_2,               &
                           Sigma1Y2_2,Sigma1Z2_2,             &
                           SigmaDeltaY2_2,SigmaDeltaZ2_2,     &
                           YOutside,ZOutside,                 &
                           YSeparation,ZSeparation,           &
                           DDBar)
!***********************************************************************
! This subroutine calculates the mean of the product of the            *
! concentrations from two sources, for continuous releases only.       *
!                                                                      *
! Changes to ADMS code: nil                                            *
!                                                                      *
! Parameters:                                                          *
! A1 A tunable constant.                                               *
!                                                                      *
! Input variables (except for Y/ZOutside, _1 and _2 refer to different *
! sources):                                                            *
! SigmaVel2      Average of the three velocity component variances.    *
! Epsilon        Dissipation rate.                                     *
! TMS            Travel time, t - s.                                   *
! XAv            Averaging length.                                     *
! Sigma0Y/Z2     Components of the variance of the source's spatial    *
!                   distribution.                                      *
! DBarMax        Maximum value of the mean concentration at downwind   *
!                   distance where output is required (the maximum is  *
!                   over Y and Z).                                     *
! HY/Z           Dilution factors.                                     *
! Sigma1Y/Z2     Components of the variance of the displacement of a   *
!                   single particle.                                   *
! SigmaDeltaY/Z2 These are defined in the documentation.               *
! Y/ZOutside     1 if receptor lies outside the sources in the Y/Z     *
!                   direction and -1 otherwise.                        *
! Y/ZSeparation  Distance between the two sources in the Y/Z           *
!                   direction.                                         *
!                                                                      *
! Output variables:                                                    *
! DDBar Mean of the product of the concentrations from two sources.    *
!                                                                      *
! Local variables:                                                     *
! A,B,Beta,C,Chi ) These are defined in the documentation.             *
! FY/Z,GY/Z      )                                                     *
! Sigma0Y/Z2_eff   Effective Sigma0Y/Z2 for calculating two source mu. *
! RDum             Dummy variable of no interest.                      *
! RMu_12           Two source mu.                                      *
! RMu_1,RMu_2      Temporary variables used in calculating mu.         *
!***********************************************************************
      IMPLICIT NONE
!..........Parameters.
      REAL*4 A1
      PARAMETER (A1 = 1.0)                            ! Tune
!..........Variables in argument list.
      REAL*4 SigmaVel2_1,Epsilon_1,TMS_1,XAv_1, &
             Sigma0Y2_1,Sigma0Z2_1,             &
             Sigma1Y2_1,Sigma1Z2_1,             &
             SigmaDeltaY2_1,SigmaDeltaZ2_1,     &
             SigmaVel2_2,Epsilon_2,TMS_2,XAv_2, &
             Sigma0Y2_2,Sigma0Z2_2,             &
             Sigma1Y2_2,Sigma1Z2_2,             &
             SigmaDeltaY2_2,SigmaDeltaZ2_2,     &
             YOutside,ZOutside,                 &
             YSeparation,ZSeparation
      REAL*8 DBarMax_1,HY_1,HZ_1,DBarMax_2,HY_2,HZ_2,DDBar
!..........Local variables.
      REAL*4 A,B,Beta,C,Chi,GY,GZ,      &
             Sigma0Y2_eff,Sigma0Z2_eff, &
             RDum,RMu_12,RMu_1,RMu_2
      REAL*8 FY,FZ

!..........A, B, C, Chi, Beta, G and F in y direction.
      A = Sigma1Y2_1 + Sigma0Y2_1
      C = Sigma1Y2_2 + Sigma0Y2_2
      B = (Sigma1Y2_1 - SigmaDeltaY2_1)*(Sigma1Y2_2 - SigmaDeltaY2_2)
      B = SQRT(B)
      Chi = A*C/(A*C - B**2)
      Beta = 2.0*B*SQRT(A*C)*YOutside/(A*C - B**2)
      GY = SQRT(Chi)
      FY = DEXP(Chi*DLOG(HY_1*HY_2) + Beta*SQRT(DLOG(HY_1)*DLOG(HY_2)))

!..........A, B, C, Chi, Beta, G and F in z direction.
      A = Sigma1Z2_1 + Sigma0Z2_1
      C = Sigma1Z2_2 + Sigma0Z2_2
      B = (Sigma1Z2_1 - SigmaDeltaZ2_1)*(Sigma1Z2_2 - SigmaDeltaZ2_2)
      B = SQRT(B)
      Chi = A*C/(A*C - B**2)
      Beta = 2.0*B*SQRT(A*C)*ZOutside/(A*C - B**2)
      GZ = SQRT(Chi)
      FZ = DEXP(Chi*DLOG(HZ_1*HZ_2) + Beta*SQRT(DLOG(HZ_1)*DLOG(HZ_2)))

!..........Sigma0Y/Z2_eff.
      Sigma0Y2_eff = 0.5*(SQRT(Sigma1Y2_1) - SQRT(Sigma1Y2_2))**2 + &
                     0.5*(Sigma0Y2_1 + Sigma0Y2_2) +                &
                     0.5*A1*YSeparation**2
      Sigma0Z2_eff = 0.5*(SQRT(Sigma1Z2_1) - SQRT(Sigma1Z2_2))**2 + &
                     0.5*(Sigma0Z2_1 + Sigma0Z2_2) +                &
                     0.5*A1*ZSeparation**2

!..........RMu_12.
      CALL MuCalc(1,SigmaVel2_1,Epsilon_1,TMS_1,XAv_1, &
                  RDum,Sigma0Y2_eff,Sigma0Z2_eff,      &
                  RDum,SigmaDeltaY2_1,SigmaDeltaZ2_1,  &
                  RMu_1)
      CALL MuCalc(1,SigmaVel2_2,Epsilon_2,TMS_2,XAv_2, &
                  RDum,Sigma0Y2_eff,Sigma0Z2_eff,      &
                  RDum,SigmaDeltaY2_2,SigmaDeltaZ2_2,  &
                  RMu_2)
      RMu_12 = 1.0 + SQRT((RMu_1 - 1.0)*(RMu_2 - 1.0))

!..........DDBar.
      DDBar = RMu_12*DBarMax_1*DBarMax_2*GY*GZ*FY*FZ

      END SUBROUTINE CCBarCalc
!
! ######################################################################
!
      SUBROUTINE MuCalc(ICalcType,SigmaVel2,Epsilon,TMS,XAv,    &
                        Sigma0X2,Sigma0Y2,Sigma0Z2,             &
                        SigmaDeltaX2,SigmaDeltaY2,SigmaDeltaZ2, &
                        RMu)
!***********************************************************************
! This subroutine calculates mu.                                       *
!                                                                      *
! Changes to ADMS code: nil                                            *
!                                                                      *
! Parameters:                                                          *
! A2 ) Tunable constants.                                              *
! A3 )                                                                 *
!                                                                      *
! Input variables:                                                     *
! ICalcType        Indicates the calculation type (1 denotes a         *
!                     calculation for a continuous release, 2 denotes  *
!                     the calculation of time integrated quantities    *
!                     for a finite duration release, and 3 denotes a   *
!                     calculation of quantities as a function of time  *
!                     for a finite duration release).                  *
! SigmaVel2        Average of the three velocity component variances.  *
! Epsilon          Dissipation rate.                                   *
! TMS              Travel time, t - s.                                 *
! XAv              Averaging length.                                   *
! Sigma0X/Y/Z2     Components of the variance of the source's spatial  *
!                     distribution.                                    *
! SigmaDeltaX/Y/Z2 These are defined in the documentation.             *
!                                                                      *
! Output variables:                                                    *
! RMu Mu.                                                              *
!                                                                      *
! Local variables:                                                     *
! SigRatioA/L/P (Sigma0/SigmaDelta)**2 appropriate for area,line and   *
!                  point sources (in Taylor transformed frame).        *
! Temp          Temporary variable.                                    *
! Tstar         TMS*Epsilon/SigmaVel2.                                 *
! RMuA/L/P      Values of mu appropriate for area,line and point       *
!                  sources (in Taylor transformed frame).              *
! RLC           Length scale for damping mu.                           *
!***********************************************************************
      IMPLICIT NONE
!..........Parameters.
      REAL*4 A2,A3
      PARAMETER (A2 = 1.0, A3 = 1.0)                       ! Tune
!..........Variables in argument list.
      INTEGER*4 ICalcType
      REAL*4 SigmaVel2,Epsilon,TMS,XAv,              &
             Sigma0X2,Sigma0Y2,Sigma0Z2,             &
             SigmaDeltaX2,SigmaDeltaY2,SigmaDeltaZ2, &
             RMu
!..........Local variables.
      REAL*4 SigRatioA,SigRatioL,SigRatioP,Temp,Tstar,RMuA,RMuL,RMuP,RLC

!..........Tstar.
      Tstar = TMS*Epsilon/SigmaVel2

!..........SigRatioA/L/P (these should be increasing in size).
      SigRatioA = SQRT(Sigma0Y2/SigmaDeltaY2)
      SigRatioL = SQRT(Sigma0Z2/SigmaDeltaZ2)
      IF (ICalcType.EQ.3) SigRatioP = SQRT(Sigma0X2/SigmaDeltaX2)
      IF (SigRatioL.LT.SigRatioA) THEN
        Temp = SigRatioL
        SigRatioL = SigRatioA
        SigRatioA = Temp
      ENDIF
      IF (ICalcType.EQ.3) THEN
        IF (SigRatioP.LT.SigRatioA) THEN
          Temp = SigRatioP
          SigRatioP = SigRatioL
          SigRatioL = SigRatioA
          SigRatioA = Temp
        ELSEIF (SigRatioP.LT.SigRatioL) THEN
          Temp = SigRatioP
          SigRatioP = SigRatioL
          SigRatioL = Temp
        ENDIF
      ENDIF

!..........RMuA.
      IF (SigRatioA.EQ.0.0) THEN
        RMuA = 1.4
      ELSE
        SigRatioA = 1.0/SigRatioA
        IF (SigRatioA.GE.9.0) THEN
          RMuA = 1.4
        ELSEIF (SigRatioA.LE.0.7) THEN
          RMuA = 1.0
        ELSE
          RMuA = 1.0 + 0.4*ALOG(SigRatioA/0.7)/ALOG(9.0/0.7)
        ENDIF
      ENDIF
      RMuA = MIN(RMuA,MAX(1.0,1.4 - (0.4/3)*Tstar))

!..........RMuL.
      IF (SigRatioL.EQ.0.0) THEN
        RMuL = 2.8
      ELSE
        SigRatioL = 1.0/SigRatioL
        IF (SigRatioL.GE.17.0) THEN
          RMuL = 2.8
        ELSEIF (SigRatioL.LE.0.9) THEN
          RMuL = 1.0
        ELSE
          RMuL = 1.0 + 1.8*ALOG(SigRatioL/0.9)/ALOG(17.0/0.9)
        ENDIF
      ENDIF
      RMuL = MIN(RMuL,MAX(1.0,2.8 - 0.6*Tstar))

      RMu = MAX(RMuA,RMuL)

!..........RMuP.
      IF (ICalcType.EQ.3) THEN
        IF (SigRatioP.EQ.0.0) THEN
          RMuP = 12.0
        ELSE
          SigRatioP = 1.0/SigRatioP
          IF (SigRatioP.GE.100.0) THEN
            RMuP = 12.0
          ELSEIF (SigRatioP.LE.1.0) THEN
            RMuP = 1.0
          ELSE
            RMuP = 1.0 + 11.0*ALOG(SigRatioP/1.0)/ALOG(100.0/1.0)
          ENDIF
        ENDIF
        RMuP = MIN(RMuP,MAX(1.0,12.0 - (11.0/3.0)*Tstar))
        RMu = MAX(RMu,RMuP)
      ENDIF

!..........Time-averaging damping.
      RLC = 1.0/(A2*SQRT(Epsilon*TMS**3)) + Epsilon/(A3*SigmaVel2**1.5)
      RLC = 1.0/RLC
      RMu = 1.0 + (RMu - 1.0)/(1.0 + XAv/RLC)

      END SUBROUTINE MuCalc
!
! ######################################################################
!
      SUBROUTINE GAMMACALC(CRATIO,GAMMA)
!***********************************************************************
! This subroutine calculates GAMMA as a function of CRATIO.            *
!                                                                      *
! Changes to ADMS code: removed local parameter PI                     *
!                                                                      *
! Input variables:                                                     *
! CRATIO mean square concentration divided by the mean concentration   *
!            squared.                                                  *
!                                                                      *
! Output variables:                                                    *
! GAMMA a parameter of the clipped normal distribution.                *
!                                                                      *
! Local variables:                                                     *
! EST      )   values of CRATIO corresponding to GAMMA, GAMMALOWER     *
! ESTLOWER )      and GAMMAUPPER.                                      *
! ESTUPPER )                                                           *
! GAMMALOWER ) temporary values of GAMMA.                              *
! GAMMAUPPER )                                                         *
! LOOP          loop counter.                                          *
! ROOT2PI     SQRT(2*PI).                                              *
!***********************************************************************
      IMPLICIT NONE
!..........Variables in argument list.
      REAL*4 CRATIO,GAMMA
!..........Local variables.
      INTEGER*4 LOOP
      REAL*4 EST,ESTLOWER,ESTUPPER,GAMMALOWER,GAMMAUPPER,ROOT2PI

!==========Case 1: CRATIO close to unity.

      IF (CRATIO.LT.1.0816) THEN

        GAMMA = SQRT(1.0/(CRATIO - 1.0))

!==========Case 2: CRATIO large.

      ELSEIF (CRATIO.GT.8.35E3) THEN

        ROOT2PI = SQRT(2.0*PI)
        GAMMA = -3.5
        DO LOOP = 1,4
          GAMMA = - CRATIO/ &
                    (2.0*ROOT2PI*GAMMA*(1.0 + 6.0/(GAMMA**4)))
          GAMMA = - SQRT(2.0*ALOG(GAMMA))
        ENDDO

!==========Case 3: CRATIO not close to unity.

      ELSE

!..........Determine limits within which GAMMA lies.
        GAMMALOWER = -3.6
        GAMMAUPPER = 3.6

!..........Determine CRATIO corresponding to these limits.
        CALL CRATIOCALC(GAMMALOWER,ESTLOWER)
        CALL CRATIOCALC(GAMMAUPPER,ESTUPPER)

!..........Determine GAMMA so that the corresponding value of CRATIO
!          is accurate to 1 in 10**3.
   30   GAMMA = (GAMMALOWER + GAMMAUPPER)*0.5
        CALL CRATIOCALC(GAMMA,EST)
        IF (EST.LT.CRATIO) THEN
          GAMMAUPPER = GAMMA
          ESTUPPER = EST
        ELSEIF (EST.GE.CRATIO) THEN
          GAMMALOWER = GAMMA
          ESTLOWER = EST
        ENDIF
        IF (ABS(ESTLOWER - ESTUPPER).GT.1.0E-3*CRATIO) THEN
          GOTO 30
        ELSE
          GAMMA = (GAMMALOWER + GAMMAUPPER)*0.5
        ENDIF
      ENDIF

      END SUBROUTINE GAMMACALC
!
! ######################################################################
!
      SUBROUTINE CRATIOCALC(GAMMA,CRATIO)
!***********************************************************************
! This subroutine calculates CRATIO as a function of GAMMA and, where  *
! appropriate, makes use of the asymptotic expressions which are valid *
! for extreme values of GAMMA.                                         *
!                                                                      *
! Changes to ADMS code: nil                                            *
!                                                                      *
! Input variables:                                                     *
! GAMMA a parameter of the clipped normal distribution.                *
!                                                                      *
! Output variables:                                                    *
! CRATIO mean square concentration divided by the mean concentration   *
!            squared.                                                  *
!***********************************************************************
      IMPLICIT NONE
!..........Variables in argument list.
      REAL*4 GAMMA,CRATIO

!..........Calculate value for large -ve GAMMA.
      IF (GAMMA.LT.-3.5) THEN
        CALL CRATIOCALCNEG(GAMMA,CRATIO)

!..........Calculate value for large +ve GAMMA.
      ELSEIF (GAMMA.GT.3.5) THEN
        CALL CRATIOCALCPOS(GAMMA,CRATIO)

!..........Calculate value for intermediate GAMMA.
      ELSE
        CALL CRATIOCALCINT(GAMMA,CRATIO)
      ENDIF

      END SUBROUTINE CRATIOCALC
!
! ######################################################################
!
      SUBROUTINE CRATIOCALCNEG(GAMMA,CRATIO)
!***********************************************************************
! This subroutine calculates CRATIO as a function of GAMMA using the   *
! asymptotic expression which is valid for large negative GAMMA.       *
!                                                                      *
! Changes to ADMS code: removed local parameter PI                     *
!                                                                      *
! Input variables:                                                     *
! GAMMA a parameter of the clipped normal distribution.                *
!                                                                      *
! Output variables:                                                    *
! CRATIO mean square concentration divided by the mean concentration   *
!            squared.                                                  *
!                                                                      *
! Local variables:                                                     *
! ROOT2PI SQRT(2*PI).                                                  *
!***********************************************************************
      IMPLICIT NONE
!..........Variables in argument list.
      REAL*4 GAMMA,CRATIO
!..........Local variables.
      REAL*4 ROOT2PI

!..........Calculate CRATIO.
      ROOT2PI = SQRT(2.0*PI)
      CRATIO = -2.0*ROOT2PI*GAMMA*EXP(GAMMA*GAMMA*0.5)* &
                (1.0 + 6.0/(GAMMA**4))

      END SUBROUTINE CRATIOCALCNEG
!
! ######################################################################
!
      SUBROUTINE CRATIOCALCPOS(GAMMA,CRATIO)
!***********************************************************************
! This subroutine calculates CRATIO as a function of GAMMA using the   *
! asymptotic expression which is valid for large positive GAMMA.       *
!                                                                      *
! Changes to ADMS code: nil                                            *
!                                                                      *
! Input variables:                                                     *
! GAMMA a parameter of the clipped normal distribution.                *
!                                                                      *
! Output variables:                                                    *
! CRATIO mean square concentration divided by the mean concentration   *
!            squared.                                                  *
!***********************************************************************
      IMPLICIT NONE
!..........Variables in argument list.
      REAL*4 GAMMA,CRATIO

!..........Calculate CRATIO.
      CRATIO = 1.0 + 1.0/(GAMMA*GAMMA)

      END SUBROUTINE CRATIOCALCPOS
!
! ######################################################################
!
      SUBROUTINE CRATIOCALCINT(GAMMA,CRATIO)
!***********************************************************************
! This subroutine calculates CRATIO as a function of GAMMA without     *
! making use of the asymptotic expressions which are valid for extreme *
! values of GAMMA. As a result, the routine is inaccurate for large    *
! negative GAMMA and inefficient for large positive GAMMA.             *
!                                                                      *
! Changes to ADMS code: removed local parameter PI                     *
!                                                                      *
! Input variables:                                                     *
! GAMMA a parameter of the clipped normal distribution.                *
!                                                                      *
! Output variables:                                                    *
! CRATIO mean square concentration divided by the mean concentration   *
!            squared.                                                  *
!                                                                      *
! Local variables:                                                     *
! ROOT2PI SQRT(2*PI).                                                  *
! EFUNC    value of ErrorFunc.                                         *
! CBAR  )  mean and mean square concentrations corresponding to the    *
! C2BAR )     given value of GAMMA and SIGMA (the other parameter of   *
!              the clipped normal distribution) equal to unity.        *
!                                                                      *
! Function subprograms:                                                *
! ErrorFunc error function.                                            *
!***********************************************************************
      IMPLICIT NONE
!..........Variables in argument list.
      REAL*4 GAMMA,CRATIO
!..........Local variables.
      REAL*4 ROOT2PI,EFUNC,CBAR,C2BAR

!..........Calculate CRATIO.
      ROOT2PI = SQRT(2.0*PI)
      EFUNC = GAMMA/SQRT(2.0)
      EFUNC = ErrorFunc(EFUNC)
      CBAR = EXP(-GAMMA*GAMMA*0.5)/ROOT2PI + &
              0.5*GAMMA*(1.0 + EFUNC)
      C2BAR = GAMMA*EXP(-GAMMA*GAMMA*0.5)/ROOT2PI + &
               0.5*(1.0 + GAMMA*GAMMA)*(1.0 + EFUNC)
      CRATIO = C2BAR/(CBAR*CBAR)

      END SUBROUTINE CRATIOCALCINT
!
! ######################################################################
!
      SUBROUTINE SIGMACALC(CBAR,GAMMA,SIGMA)
!***********************************************************************
! This subroutine calculates SIGMA as a function of CBAR and GAMMA     *
! and, where appropriate, makes use of the asymptotic expressions      *
! which are valid for extreme values of GAMMA.                         *
!                                                                      *
! Changes to ADMS code: nil                                            *
!                                                                      *
! Input variables:                                                     *
! CBAR mean concentration.                                             *
! GAMMA a parameter of the clipped normal distribution.                *
!                                                                      *
! Output variables:                                                    *
! SIGMA a parameter of the clipped normal distribution.                *
!***********************************************************************
      IMPLICIT NONE
!..........Variables in argument list.
      REAL*4 GAMMA,CBAR,SIGMA

!..........Calculate value for large -ve GAMMA.
      IF (GAMMA.LT.-3.5) THEN
        CALL SIGMACALCNEG(CBAR,GAMMA,SIGMA)

!..........Calculate value for large +ve GAMMA.
      ELSEIF (GAMMA.GT.3.5) THEN
        CALL SIGMACALCPOS(CBAR,GAMMA,SIGMA)

!..........Calculate value for intermediate GAMMA.
      ELSE
        CALL SIGMACALCINT(CBAR,GAMMA,SIGMA)
      ENDIF

      END SUBROUTINE SIGMACALC
!
! ######################################################################
!
      SUBROUTINE SIGMACALCNEG(CBAR,GAMMA,SIGMA)
!***********************************************************************
! This subroutine calculates SIGMA as a function of CBAR and GAMMA     *
! using the asymptotic expression which is valid for large negative    *
! GAMMA.                                                               *
!                                                                      *
! Changes to ADMS code: removed local parameter PI                     *
!                                                                      *
! Input variables:                                                     *
! CBAR mean concentration.                                             *
! GAMMA a parameter of the clipped normal distribution.                *
!                                                                      *
! Output variables:                                                    *
! SIGMA a parameter of the clipped normal distribution.                *
!                                                                      *
! Local variables:                                                     *
! RATIO  CBAR/SIGMA.                                                   *
!***********************************************************************
      IMPLICIT NONE
!..........Variables in argument list.
      REAL*4 GAMMA,CBAR,SIGMA
!..........Local variables.
      REAL*4 RATIO

!..........Calculate SIGMA.
      RATIO = (1.0/GAMMA**2 - 3.0/GAMMA**4 + 15.0/GAMMA**6)* &
               EXP(-0.5*GAMMA*GAMMA)/SQRT(2.0*PI)
      SIGMA = CBAR/RATIO

      END SUBROUTINE SIGMACALCNEG
!
! ######################################################################
!
      SUBROUTINE SIGMACALCPOS(CBAR,GAMMA,SIGMA)
!***********************************************************************
! This subroutine calculates SIGMA as a function of CBAR and GAMMA     *
! using the asymptotic expression which is valid for large positive    *
! GAMMA.                                                               *
!                                                                      *
! Changes to ADMS code: nil                                            *
!                                                                      *
! Input variables:                                                     *
! CBAR mean concentration.                                             *
! GAMMA a parameter of the clipped normal distribution.                *
!                                                                      *
! Output variables:                                                    *
! SIGMA a parameter of the clipped normal distribution.                *
!***********************************************************************
      IMPLICIT NONE
!..........Variables in argument list.
      REAL*4 GAMMA,CBAR,SIGMA

!..........Calculate SIGMA.
      SIGMA = CBAR/GAMMA

      END SUBROUTINE SIGMACALCPOS
!
! ######################################################################
!
      SUBROUTINE SIGMACALCINT(CBAR,GAMMA,SIGMA)
!***********************************************************************
! This subroutine calculates SIGMA as a function of CBAR and GAMMA     *
! without making use of the asymptotic expressions which are valid for *
! extreme values of GAMMA. As a result, the routine is inaccurate for  *
! large negative GAMMA and inefficient for large positive GAMMA.       *
!                                                                      *
! Changes to ADMS code: removed local parameter PI                     *
!                                                                      *
! Input variables:                                                     *
! CBAR mean concentration.                                             *
! GAMMA a parameter of the clipped normal distribution.                *
!                                                                      *
! Output variables:                                                    *
! SIGMA a parameter of the clipped normal distribution.                *
!                                                                      *
! Local variables:                                                     *
! EFUNC value of ErrorFunc.                                            *
! RATIO  CBAR/SIGMA.                                                   *
!                                                                      *
! Function subprograms:                                                *
! ErrorFunc error function.                                            *
!***********************************************************************
      IMPLICIT NONE
!..........Variables in argument list.
      REAL*4 GAMMA,CBAR,SIGMA
!..........Local variables.
      REAL*4 EFUNC,RATIO

!..........Calculate SIGMA.
      EFUNC = GAMMA/SQRT(2.0)
      EFUNC = ErrorFunc(EFUNC)
      RATIO = 0.5*GAMMA*(1.0 + EFUNC) + EXP(-0.5*GAMMA*GAMMA)/ &
              SQRT(2.0*PI)
      SIGMA = CBAR/RATIO

      END SUBROUTINE SIGMACALCINT
!
! ######################################################################
!
      SUBROUTINE DOSECALC(PDOSE,GAMMA,SIGMA,DOSE)
!***********************************************************************
! This subroutine calculates DOSE by numerical integration using       *
! Simpson's rule.                                                      *
!                                                                      *
! Changes to ADMS code: nil                                            *
!                                                                      *
! Parameters:                                                          *
! NUMBEROFPOINTS  number of concentration values used in the           *
!                      numerical integration.                          *
! HOWFARINTOTAIL a measure of how far into the tail of the             *
!                      distribution the numerical integration is to    *
!                      proceed.                                        *
!                                                                      *
! Input variables:                                                     *
! PDOSE power used in calculating DOSE.                                *
! GAMMA  a parameter of the clipped normal distribution.               *
! SIGMA  a parameter of the clipped normal distribution.               *
!                                                                      *
! Output variables:                                                    *
! DOSE the ensemble average of the PDOSEth power of the                *
!         concentration.                                               *
!                                                                      *
! Local variables:                                                     *
! A )          end points of integration.                              *
! B )                                                                  *
! C            concentration.                                          *
! H            concentration interval used in integration.             *
! DOSECONTRIB contribution to the dose from a particular               *
!                 concentration value.                                 *
! LOOP         loop counter.                                           *
!***********************************************************************
      IMPLICIT NONE
!..........Parameters.
      INTEGER*4 NUMBEROFPOINTS
      REAL*4 HOWFARINTOTAIL
      PARAMETER (NUMBEROFPOINTS = 30, HOWFARINTOTAIL = 6.0)
!..........Variables in argument list.
      REAL*4 PDOSE,GAMMA,SIGMA
      REAL*8 DOSE
!..........Local variables.
      INTEGER*4 LOOP
      REAL*4 A,B,C,H
      REAL*8 DOSECONTRIB

!..........Set up A, B, H.
      IF (GAMMA.GE.0.0) THEN
        A = MAX(0.0,(GAMMA - HOWFARINTOTAIL)*SIGMA)
        B = (GAMMA + HOWFARINTOTAIL)*SIGMA
      ELSE
        A = 0.0
        B = HOWFARINTOTAIL*MIN(SIGMA,-3.0*SIGMA/GAMMA)
      ENDIF
      H = (B - A)/FLOAT(NUMBEROFPOINTS)

!..........First point.
      C = A
      CALL DOSECONTRIBCALC(PDOSE,GAMMA,SIGMA,C,DOSECONTRIB)
      DOSE = DOSECONTRIB

!..........Last point.
      C = B
      CALL DOSECONTRIBCALC(PDOSE,GAMMA,SIGMA,C,DOSECONTRIB)
      DOSE = DOSE + DOSECONTRIB

!..........Intermediate points.
      DO LOOP = 1,NUMBEROFPOINTS/2 - 1
        C = A + H*FLOAT(2*LOOP - 1)
        CALL DOSECONTRIBCALC(PDOSE,GAMMA,SIGMA,C,DOSECONTRIB)
        DOSE = DOSE + 4.0D0*DOSECONTRIB
        C = A + H*FLOAT(2*LOOP)
        CALL DOSECONTRIBCALC(PDOSE,GAMMA,SIGMA,C,DOSECONTRIB)
        DOSE = DOSE + 2.0D0*DOSECONTRIB
      ENDDO

!..........Penultimate point.
      LOOP = NUMBEROFPOINTS/2
      C = A + H*FLOAT(2*LOOP - 1)
      CALL DOSECONTRIBCALC(PDOSE,GAMMA,SIGMA,C,DOSECONTRIB)
      DOSE = DOSE + 4.0D0*DOSECONTRIB

!..........Calculate result.
      DOSE = DOSE*H/3.0

      END SUBROUTINE DOSECALC
!
! ######################################################################
!
      SUBROUTINE DOSECONTRIBCALC(PDOSE,GAMMA,SIGMA,C,DOSECONTRIB)
!***********************************************************************
! This subroutine calculates the p.d.f. of the concentration at the    *
! value C multiplied by C to the power PDOSE.                          *
!                                                                      *
! Changes to ADMS code: removed local parameter PI                     *
!                                                                      *
! Input variables:                                                     *
! PDOSE power used in calculating DOSE.                                *
! C      concentration.                                                *
! GAMMA  a parameter of the clipped normal distribution.               *
! SIGMA  a parameter of the clipped normal distribution.               *
!                                                                      *
! Output variables:                                                    *
! DOSECONTRIB contribution to the dose from a particular               *
!                 concentration value.                                 *
!                                                                      *
! Local variables:                                                     *
! CPDOSE C to the power PDOSE.                                         *
!***********************************************************************
      IMPLICIT NONE
!..........Variables in argument list.
      REAL*4 PDOSE,GAMMA,SIGMA,C
      REAL*8 DOSECONTRIB
!..........Local variables.
      REAL*8 CPDOSE

!..........Calculate DOSECONTRIB.
      IF (PDOSE.NE.0.0) THEN
        CPDOSE = DBLE(C)**PDOSE
      ELSE
        CPDOSE = 1.0D0
      ENDIF
      DOSECONTRIB = CPDOSE*EXP(DBLE(-0.5*(C/SIGMA - GAMMA)**2))/ &
                     DBLE(SIGMA*SQRT(2.0*PI))

      END SUBROUTINE DOSECONTRIBCALC
!
! ######################################################################
!
      SUBROUTINE CHatCalc(Gamma,Sigma,NPdf,CHat,Scale)
!***********************************************************************
! This subroutine calculates the values of concentration for which     *
! the exceedence probability is calculated.                            *
!                                                                      *
! Changes to ADMS code:                                                *
! - local parameter MaxNPdf replaced by the global parameter MaxPdfSize*
! - local variable CHatMax renamed CMax (for clarity)                  *
! - local variable Res added (to store pdf resolution as a real)       *
! - calculation of threshold values is modified to give standardised   *
!   concentrations on log10 scale                                      *
! - the Scale factor is returned via argument list                     *
!                                                                      *
! Input variables:                                                     *
! Gamma A parameter of the clipped normal distribution.                *
! Sigma A parameter of the clipped normal distribution.                *
! NPdf  Number of concentration values                                 *
!                                                                      *
! Output variables:                                                    *
! CHat  Concentration values for which the exceedence probability is   *
!         calculated.                                                  *
! Scale Index of the largest concentration threshold (essentially acts *
!         as a scale factor for the concentrations).                   *
!                                                                      *
! Local variables:                                                     *
! Loop Loop counter.                                                   *
! CMax Maximum concentration value considered in the calculation.      *
! Res  Resolution per decade in pdf (as a real value).                 *
!***********************************************************************
      IMPLICIT NONE
!..........Variables in argument list.
      INTEGER*4 NPdf,Scale
      REAL*4 Gamma,Sigma
      REAL*4 CHat(MaxPdfSize)
!..........Local variables.
      INTEGER*4 Loop
      REAL*4 CMax,Res

!..........Calculate CHat.
      IF (Gamma .GE. 0.0) THEN
        CMax = (6.0 + Gamma)*Sigma
      ELSEIF (Gamma .GE. -3.0) THEN
        CMax = 6.0*Sigma
      ELSE
        CMax = -18.0*Sigma/Gamma
      ENDIF
      Res   = FLOAT(AutoPdfResolutionPerDecade)
      Scale = NINT(Res * LOG10(CMax))
      DO Loop = 1,NPdf
        CHat(Loop) = 10.0**(FLOAT(Scale + Loop - NPdf)/Res)
      ENDDO

      END SUBROUTINE CHatCalc
!
! ######################################################################
!
      SUBROUTINE ProbCalc(NPdf,Gamma,Sigma,CHat,Prob)
!***********************************************************************
! This subroutine calculates the probability distribution.             *
!                                                                      *
! Changes to ADMS code:                                                *
! - local parameter MaxNPdf replaced by the global parameter MaxPdfSize*
! - argument Prob changed to assumed shape array                       *
!                                                                      *
! Input variables:                                                     *
! NPdf    Number of concentration values in CHat.                      *
! Gamma ) Parameters of the clipped normal distribution.               *
! Sigma )                                                              *
! CHat    Concentration values for which the exceedence probability is *
!            calculated.                                               *
!                                                                      *
! Output variables:                                                    *
! Prob Probability of concentration exceeding the value given in the   *
!         corresponding array element in CHat.                         *
!                                                                      *
! Local variables:                                                     *
! Loop Loop counter.                                                   *
! P    Argument of error function.                                     *
! DP   Double precision version of P.                                  *
!                                                                      *
! Function subprograms:                                                *
! ErrorFunc Error function.                                            *
!***********************************************************************
      IMPLICIT NONE
!..........Variables in argument list.
      INTEGER*4 NPdf
      REAL*4 Gamma,Sigma
      REAL*4 Prob(:),CHat(MaxPdfSize)
!..........Local variables.
      INTEGER*4 Loop
      REAL*4 P
      REAL*8 DP

!..........Calculate Prob (setting 100th percentile equal to 99.999th
!          percentile to avoid infinite concentrations).
      DO Loop = 1,NPdf
        DP = (CHat(Loop)/Sigma - Gamma)/SQRT(2.0)
        IF (DP.GT.1.0E20) DP = 1.0E20
        IF (DP.LT.-1.0E20) DP = -1.0E20
        P = DP
        Prob(Loop) = 0.5*(1.0 - ErrorFunc(P))
        IF (Prob(Loop).LT.0.00001) Prob(Loop) = 0.0
      ENDDO

      END SUBROUTINE ProbCalc
!
! ######################################################################
!
      SUBROUTINE ConcCalc(NPerc,Gamma,Sigma,Percentile,Conc)
!***********************************************************************
! This subroutine calculates some percentiles of the concentration     *
! probability distribution.                                            *
!                                                                      *
! Changes to ADMS code:                                                *
! - local parameter MaxNPdf replaced by the global parameter MaxPdfSize*
! - argument Conc changed to assumed shape array                       *
!                                                                      *
! Input variables:                                                     *
! NPerc      Number of percentiles required.                           *
! Gamma )    Parameters of the clipped normal distribution.            *
! Sigma )                                                              *
! Percentile Percentiles required.                                     *
!                                                                      *
! Output variables:                                                    *
! Conc Concentrations corresponding to the desired percentiles.        *
!                                                                      *
! Local variables:                                                     *
! Loop Loop counter.                                                   *
! P    Argument of function Std_Normal_Inv.                            *
!                                                                      *
! Function subprograms:                                                *
! Std_Normal_Inv Inverse of the upper-tail integral of the standard    *
!                   normal pdf.                                        *
!***********************************************************************
      IMPLICIT NONE
!..........Variables in argument list.
      INTEGER*4 NPerc
      REAL*4 Gamma, Sigma
      REAL*4 Percentile(MaxPdfSize), Conc(:)
!..........Local variables.
      INTEGER*4 Loop
      REAL*4 P

!..........Calculate Conc (setting 100th percentile equal to 99.999th
!          percentile to avoid infinite concentrations)
      DO Loop = 1,NPerc
        IF (Percentile(Loop).EQ.100.0) THEN
          P = 1.0 - 99.999/100.0
        ELSE
          P = 1.0 - Percentile(Loop)/100.0
        ENDIF
        IF (P .LE. 0.0) THEN
          Conc(Loop) = 1.0E24
        ELSEIF (P .GE. 1.0) THEN
          Conc(Loop) = 0.0
        ELSE
          Conc(Loop) = MAX((Std_Normal_Inv(P) + Gamma)*Sigma, 0.0)
        ENDIF
      ENDDO

      END SUBROUTINE ConcCalc
!
! ######################################################################
!
      FUNCTION ErrorFunc(X)
!***********************************************************************
! This function calculates the error function using an approximation   *
! given by Abramowitz and Stegun, Handbook of mathematical functions,  *
! Dover Publications (1965).                                           *
!                                                                      *
! Changes to ADMS code: nil                                            *
!                                                                      *
! Parameters:                                                          *
! P  ) numbers used in the numerical approximation of the error        *
! A1 )    function.                                                    *
! A2 )                                                                 *
! A3 )                                                                 *
! A4 )                                                                 *
! A5 )                                                                 *
!                                                                      *
! Input variables:                                                     *
! X argument of error function.                                        *
!                                                                      *
! Local variables:                                                     *
! Y absolute value of X                                                *
! T 1.0/(1.0 + P*Y)                                                    *
!                                                                      *
! Function name:                                                       *
! ErrorFunc error function.                                            *
!***********************************************************************
      IMPLICIT NONE
!..........Parameters.
      REAL*4 P,A1,A2,A3,A4,A5
      PARAMETER(P = 0.3275911, A1 = 0.254829592, A2 = -0.284496736, &
                A3 = 1.421413741, A4 = -1.453152027, A5 = 1.061405429)
!..........Variables in argument list.
      REAL*4 X
!..........Local variables.
      REAL*4 Y,T
!..........Function name.
      REAL*4 ErrorFunc

      Y = ABS(X)
      T = 1.0/(1.0 + P*Y)
      ErrorFunc = 1.0 - EXP(-Y*Y)* &
                   (A1*T + A2*T**2 + A3*T**3 + A4*T**4 + A5*T**5)
      IF (X.LT.0.0) ErrorFunc = -ErrorFunc

      END FUNCTION ErrorFunc
!
! ######################################################################
!
      FUNCTION Std_Normal_Inv(P)
!***********************************************************************
! This function calculates the inverse of the upper-tail integral of the
! standard normal distribution using an approximation given by
! Abramowitz and Stegun, Handbook of mathematical functions, Dover
! Publications (1965).
!
! Changes to ADMS code: nil
!
! Parameters:
! C0 ) numbers used in the numerical approximation of the function.
! C1 )
! C2 )
! D1 )
! D2 )
! D3 )
!
! Input variables:
! P argument of function.
!
! Local variables:
! Y 1 - P if P > 0.5; otherwise P
! T SQRT(ALOG(1.0/(Y*Y)))
!
! Function name:
! Std_Normal_Inv required function.
!***********************************************************************
      IMPLICIT NONE
!..........Parameters.
      REAL*4 C0, C1, C2, D1, D2, D3
      PARAMETER(C0 = 2.515517, C1 = 0.802853, C2 = 0.010328, &
                D1 = 1.432788, D2 = 0.189269, D3 = 0.001308)
!..........Variables in argument list.
      REAL*4 P
!..........Local variables.
      REAL*4 T
      REAL*8 Y
      REAL*4 Z
!..........Function name.
      REAL*4 Std_Normal_Inv

      IF (P .GT. 0.5) THEN
        Y = 1.0 - P
      ELSE
        Y = P
      ENDIF
      T = SQRT(LOG(1.0D0/(Y*Y)))
      Z = T - (C0 + (C1 + C2*T)*T)/(1.0 + (D1 + (D2 + D3*T)*T)*T)
      IF (P .GT. 0.5) THEN
        Std_Normal_Inv = -Z
      ELSE
        Std_Normal_Inv = Z
      ENDIF

      END FUNCTION Std_Normal_Inv



!#### ADMS CODE ENDS





!-------------------------------------------------------------------------------------------------------------

End Module FluctuationModule

! Module:  Plume Rise Module

Module PlumeRiseModule

! This module provides code for calculating plume rise.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: PlumeRise

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine PlumeRise(U, V, W, T, TH, PAMBIENT, UOLD, VOLD, WOLD, TOLD,       &
                     THOLD, PAMBIENTOLD, TAUW, SIGW, EPS,                    &
                     IPARTAGE, IDELT, DEBUG,                                 &
                     UP, VP, WP, RANDPLUMSIG,                                &
                     FMASS, FMX, FMY, FMZ, FH, FMASS0, B0OLD)

! DZBYDETA, R0, TSTACK, and VS removed from NAME arg list.
! RANDPLUMX/Y/Z replaced by RANDPLUMSIG, the standard deviation of RANDPLUMX/Y/Z

! Purpose
! Determine plume rise
! Drag included
!
!
! Method
!   Using equations of conservation of mass, momentum and energy, the
!    new position of each particle including the effects of plume rise
!    is calculated.
!
!   Original code  H N Webster


  ! ambient variables at particle position
  REAL,    Intent(In)    :: U         ! AMBIENT U VELOCITY AT PARTICLE POSITION
  REAL,    Intent(In)    :: V         ! AMBIENT V VELOCITY AT PARTICLE POSITION
  REAL,    Intent(In)    :: W         ! AMBIENT W VELOCITY AT PARTICLE POSITION
  REAL,    Intent(In)    :: PAMBIENT  ! ATMOSPHERIC PRESSURE AT PARTICLE POSITION
  REAL,    Intent(In)    :: T         ! ATMOSPHERIC TEMPERATURE AT PARTICLE POSITION
  REAL,    Intent(In)    :: TH        ! AMBIENT POTENTIAL TEMPERATURE AT PARTICLE
                                      ! POSITION
  REAL,    Intent(In)    :: UOLD      ! AMBIENT U VELOCITY AT PARTICLE POSITION
                                      ! AT PREVIOUS TIME STEP
  REAL,    Intent(In)    :: VOLD      ! AMBIENT V VELOCITY AT PARTICLE POSITION
                                      ! AT PREVIOUS TIME STEP
  REAL,    Intent(In)    :: WOLD      ! AMBIENT W VELOCITY AT PARTICLE POSITION
                                      ! AT PREVIOUS TIME STEP
  REAL,    Intent(In)    :: PAMBIENTOLD  ! ATMOSPHERIC PRESSURE AT PARTICLE POSITION
                                         ! AT PREVIOUS TIME STEP
  REAL,    Intent(In)    :: TOLD      ! ATMOSPHERIC TEMPERATURE AT PARTICLE POSITION
                                      ! AT PREVIOUS TIME STEP
  REAL,    Intent(In)    :: THOLD     ! AMBIENT POTENTIAL TEMPERATURE AT PARTICLE
                                      ! POSITION AT PREVIOUS TIME STEP

!  REAL,    Intent(In)    :: DTHETADX  ! DTHETA/DX
!  REAL,    Intent(In)    :: DTHETADY  ! DTHETA/DY
!  REAL,    Intent(In)    :: DTHETADZ  ! DTHETA/DZ
  !REAL,    Intent(In)    :: DZBYDETA  ! DZ/DETA ! Not needed if w not eta-dot.
!  REAL,    Intent(In)    :: DWDX      ! DW/DX
!  REAL,    Intent(In)    :: DWDY      ! DW/DY
!  REAL,    Intent(In)    :: DWDZ      ! DW/DZ
!  REAL,    Intent(In)    :: DUDX      ! DU/DX
!  REAL,    Intent(In)    :: DUDY      ! DU/DY
!  REAL,    Intent(In)    :: DUDZ      ! DU/DZ
!  REAL,    Intent(In)    :: DVDX      ! DV/DX
!  REAL,    Intent(In)    :: DVDY      ! DV/DY
!  REAL,    Intent(In)    :: DVDZ      ! DV/DZ
  REAL,    Intent(In)    :: TAUW      ! VERTICAL LAGRANGIAN TIMESCALE
  REAL,    Intent(In)    :: SIGW      ! VELOCITY VARIANCE
  REAL,    Intent(In)    :: EPS       ! TURBULENT DECAY

  ! stack variables - not needed as particle initialised elsewhere.
  !REAL,    Intent(In)    :: R0        ! STACK RADIUS
  !REAL,    Intent(In)    :: VS        ! EMISSION VELOCITY
  !REAL,    Intent(In)    :: TSTACK    ! EMISSION TEMPERATURE

  ! other variables
  REAL, Intent(In)    :: IPARTAGE  ! AGE OF PARTICLE
  REAL, Intent(In)    :: IDELT     ! TIMESTEP
  LOGICAL, Intent(In)    :: DEBUG     ! DEBUG FLAG

  ! plume variables
  REAL,    Intent(Out)   :: UP        ! U VELOCITY OF THE PARTICLE/PLUME
  REAL,    Intent(Out)   :: VP        ! V VELOCITY OF THE PARTICLE/PLUME
  REAL,    Intent(Out)   :: WP        ! W VELOCITY OF THE PARTICLE/PLUME
  REAL,    Intent(Out)   :: RANDPLUMSIG ! Standard deviation of random displacement in
                                        ! each coord direction due to plume induced
                                        ! turbulence.
  ! fluxes
  REAL,    Intent(InOut) :: FH        ! HEAT FLUX AT TIME IPARTAGE
  REAL,    Intent(InOut) :: FMASS     ! MASS FLUX AT TIME IPARTAGE
  REAL,    Intent(InOut) :: FMX       ! X COMPONENT OF MOMENTUM FLUX AT TIME IPARTAGE
  REAL,    Intent(InOut) :: FMY       ! Y COMPONENT OF MOMENTUM FLUX AT TIME IPARTAGE
  REAL,    Intent(InOut) :: FMZ       ! VERTICAL COMPONENT OF MOMENTUM FLUX AT TIME
                                      ! IPARTAGE
  REAL,    Intent(InOut) :: FMASS0    ! MASS FLUX DUE TO ENTRAINMENT DUE TO PLUME RISE
                                      ! ONLY

  ! plume variables
  REAL,    Intent(InOut) :: B0OLD     ! PLUME RADIUS FROM ENTRAINMENT DUE TO
                                      ! PLUME RISE ONLY AT TIME IPARTAGE - IDELT

  ! Note:
  ! FMass = 0 switches off plume rise as a whole. This is set if
  ! Abs(w_plume - w_ambient) < 0.01 or travel time > 3600.
  ! B0Old = -999.0 switches off plume rise induced turbulence. This is set if B0 <
  ! B0Old.

  ! Locals:

  REAL          RHOA      ! AMBIENT DENSITY AT PARTICLE POSITION
!  REAL          DTHETADT  ! DTHETA/DT
!  REAL          DWDT      ! DW/DT
!  REAL          DUDT      ! DU/DT
!  REAL          DVDT      ! DV/DT

  ! plume variables
  REAL          UZETA     ! MAGNITUDE OF THE PARTICLE/PLUME VELOCITY
  REAL          THP       ! POTENTIAL TEMPERATURE OF PARTICLE/PLUME
  REAL          TP        ! TEMPERATURE OF PARTICLE/PLUME
  REAL          RHOP      ! DENSITY OF PARTICLE/PLUME
  REAL          B         ! RADIUS OF PLUME
  REAL          UE        ! ENTRAINMENT VELOCITY
  REAL          UETURB    ! ENTRAINMENT DUE TO AMBIENT TURBULENCE
  REAL          UERISE    ! ENTRAINMENT DUE TO PLUME RISE
  REAL          B0        ! PLUME RADIUS FROM ENTRAINMENT DUE TO
                          ! PLUME RISE ONLY AT TIME IPARTAGE

  ! fluxes at new times
  REAL          FHNEW     ! HEAT FLUX AT TIME (IPARTAGE + IDELT)
  REAL          FMASSNEW  ! MASS FLUX AT TIME (IPARTAGE + IDELT)
  REAL          FMXNEW    ! X COMPONENT OF MOMENTUM FLUX AT TIME
                          ! (IPARTAGE + IDELT)
  REAL          FMYNEW    ! Y COMPONENT OF MOMENTUM FLUX AT TIME
                          ! (IPARTAGE + IDELT)
  REAL          FMZNEW    ! VERTICAL COMPONENT OF MOMENTUM FLUX
                          ! AT TIME (IPARTAGE + IDELT)
  REAL          FMASS0NEW ! MASS FLUX DUE TO ENTRAINMENT FROM PLUME
                          ! RISE ONLY AT TIME (IPARTAGE + IDELT)

  ! other variables
  REAL          DELTAUZETA ! COMPONENT OF THE MAGNITUDE OF THE
                           ! DIFFERENCE BETWEEN THE PARTICLE/
                           ! PLUME AND AMBIENT VELOCITY IN THE
                           ! DIRECTION OF THE PLUME AXIS
  REAL          DELTAUNX   ! X COMPONENT OF THE DIFFERENCE
                           ! BETWEEN THE PARTICLE/ PLUME AND
                           ! AMBIENT VELOCITY
  REAL          DELTAUNY   ! Y COMPONENT OF THE DIFFERENCE
                           ! BETWEEN THE PARTICLE/ PLUME AND
                           ! AMBIENT VELOCITY
  REAL          DELTAUNZ   ! Z COMPONENT OF THE DIFFERENCE
                           ! BETWEEN THE PARTICLE/ PLUME AND
                           ! AMBIENT VELOCITY
  REAL          DELTAUN    ! COMPONENT OF THE MAGNITUDE OF THE
                           ! DIFFERENCE BETWEEN THE PARTICLE/
                           ! PLUME AND AMBIENT VELOCITY
                           ! PERPENDICULAR TO THE PLUME AXIS
  REAL          DRAG       ! TEMPORARY VARIABLE TO CALCULATE DRAG FORCE

  ! external functions - replaced by Gauss in MathsModule
  ! REAL         GRANUM           !RANDOM NUMBERS

  ! parameters
  REAL          PI            ! PI
  PARAMETER     (PI = 3.1415926)
  REAL          CP            ! SPECIFIC HEAT AT CONSTANT PRESSURE
  PARAMETER     (CP = 1004.0)
  REAL          GASCONST      ! SPECIFIC GAS CONSTANT
  PARAMETER     (GASCONST = 287.05)
  REAL          RCP           ! R/CP
  PARAMETER     (RCP = 0.2856)
  REAL          G             ! GRAVITATIONAL CONSTANT
  PARAMETER     (G = 9.81)
  REAL          ALPHA1        ! ENTRAINMENT PARAMETER
  PARAMETER     (ALPHA1 = 0.11)
  REAL          ALPHA2        ! ENTRAINMENT PARAMETER
  PARAMETER     (ALPHA2 = 0.5)
  REAL          ALPHA3        ! ENTRAINMENT PARAMETER
  PARAMETER     (ALPHA3 = 0.655)
  REAL          CD            ! DRAG COEFFICIENT
  PARAMETER     (CD = 0.21)

!
!Set inital plume variables
!

! If IPartAge = 0 block executed, routine exited, routine re-entered with IPartAge /= 0
! block executed, effect is unchanged. Hence IPartAge = 0 block not needed if
! particle initialised elsewhere.
   !   IF (IPARTAGE.EQ.0.0) THEN
   !     TP = TSTACK
   !     THP = TSTACK*(1.0E5/(PAMBIENT))**RCP
   !     B = R0
   !     UP = 0.0
   !     VP = 0.0
   !     WP = VS
   !     UZETA = SQRT(UP**2 + VP**2 + WP**2)
   !     RHOP = (PAMBIENT)/(GASCONST*TP)
   !     RHOA = (PAMBIENT)/(GASCONST*T)
   !     FMASS = PI*B**2*RHOP*UZETA
   !     FMX = (UP - U)*FMASS
   !     FMY = (VP - V)*FMASS
   !     FMZ = (WP - W)*FMASS
   !     FH = CP*(THP - TH)*FMASS
   !     IF (B0OLD.GE.0.0) THEN
   !       B0 = R0
   !       B0OLD = R0
   !       FMASS0 = PI*B0**2*RHOP*UZETA
   !     ENDIF
   !   ELSE

! calculate variables at t - dt
        UP = UOLD + FMX/FMASS
        VP = VOLD + FMY/FMASS
        WP = WOLD + FMZ/FMASS
        UZETA = SQRT(UP**2 + VP**2 + WP**2)
        THP = FH/(CP*FMASS) + THOLD
        RHOA = (PAMBIENT)/(GASCONST*T)
        TP = THP*(PAMBIENTOLD/1.0E5)**RCP
        RHOP = (PAMBIENTOLD)/(GASCONST*TP)
!        If ((FMASS/(PI*RHOP*UZETA)) > 0.0) Then
!        Else
!          write(6,*)'pambientold=',PAMBIENTOLD ! $$
!          write(6,*)'tp=',TP
!          write(6,*)'fmass=',FMASS
!          write(6,*)'rhop=',RHOP
!          write(6,*)'uzeta=',UZETA
!        End If
        B = SQRT(FMASS/(PI*RHOP*UZETA))

   !   ENDIF
!
! Calculate velocity of particle at time IPARTAGE
!

!
! Calculate the magnitude of the components of u(p) - u parallel and
! perpendicular to the plume's axis at time IPARTAGE
!
      DELTAUZETA = ((UP - UOLD)*UP + (VP - VOLD)*VP + (WP - WOLD)*WP)/UZETA
      DELTAUNX = (1.0 - DELTAUZETA/UZETA)*UP - UOLD
      DELTAUNY = (1.0 - DELTAUZETA/UZETA)*VP - VOLD
      DELTAUNZ = (1.0 - DELTAUZETA/UZETA)*WP - WOLD
      DELTAUN = SQRT((DELTAUNX)**2 + (DELTAUNY)**2 + (DELTAUNZ)**2)

!
! Calculate fluxes and positions at new timestep (IPARTAGE + IDELT)
!
!      DTHETADT = UP*DTHETADX + VP*DTHETADY + WP*DTHETADZ
      FHNEW = FH - FMASS*CP*(TH - THOLD)
!      DWDT = UP*DWDX + VP*DWDY + WP*DWDZ
      DRAG = RHOA*PI*B*DELTAUN*CD
      FMZNEW = FMZ + IDELT*(FH*G/(CP*TH) - UZETA*DRAG*DELTAUNZ) - FMASS*(W - WOLD)
      UETURB = ALPHA3 * MIN((EPS*B)**0.333333,SIGW/SQRT(1.0 + IPARTAGE/(2.0*TAUW)))
      UERISE = ALPHA1*Abs(DELTAUZETA) + ALPHA2*DELTAUN
      UE = UERISE + UETURB
      FMASSNEW = FMASS + 2.0*IDELT*PI*B*RHOA*UE*UZETA
!      DUDT = UP*DUDX + VP*DUDY + WP*DUDZ
      FMXNEW = FMX - FMASS*(U - UOLD) - IDELT*UZETA*DRAG*DELTAUNX
!      DVDT = UP*DVDX + VP*DVDY + WP*DVDZ
      FMYNEW = FMY - FMASS*(V - VOLD) - IDELT*UZETA*DRAG*DELTAUNY
      IF (B0OLD.GE.0.0) THEN
        FMASS0NEW = FMASS0 + 2.0*PI*B0OLD*RHOA*IDELT*UERISE*UZETA
      ELSE
        FMASS0NEW = -999.0
      ENDIF

!      If ((FMASS0NEW < 0.0) .And. (FMASS0NEW /= -999.0)) Then
!        write(6,*)'fmass0new=',FMASS0NEW ! $$
!        write(6,*)'fmass0=',FMASS0
!        write(6,*)'rhoa=',RHOa
!        write(6,*)'uzeta=',uzeta
!        write(6,*)'b0=',B0
!        write(6,*)'dt=',IDELT
!      End If

! update fluxes to return
      FMX = FMXNEW
      FMY = FMYNEW
      FMZ = FMZNEW
      FMASS = FMASSNEW
      FH = FHNEW
      FMASS0 = FMASS0NEW


!   If (((FMZ/FMASS) > WP) .And. ((FMZ/FMASS) > 3.0)) Then
!     write(6,*)'wp=',WP
!   write(6,*)'fmz/fmass=',FMZ/FMASS
!   End If

!   If ((FMZ/FMASS) > WP) Then
!     write(6,*)'wp=',WP
!   write(6,*)'fmz/fmass=',FMZ/FMASS
!   End If


! calculate variables at t
        UP = U + FMX/FMASS
        VP = V + FMY/FMASS
        WP = W + FMZ/FMASS
        UZETA = SQRT(UP**2 + VP**2 + WP**2)
        THP = FH/(CP*FMASS) + TH
        TP = THP*(PAMBIENT/1.0E5)**RCP
        RHOP = (PAMBIENT)/(GASCONST*TP)
        B = SQRT(FMASS/(PI*RHOP*UZETA))

        IF (B0OLD.GE.0.0) THEN
          B0 = SQRT(FMASS0/(PI*RHOP*UZETA))
        ELSE
          B0 = -999.0
        ENDIF


! calculate standard deviation for turbulence calculation for plume
! induced turbulence - -999.0 is a flag to terminate plume induced
! turbulence for this particle
      IF (B0OLD.GE.0.0) THEN
        IF (B0.GE.B0OLD) THEN
          RANDPLUMSig = SQRT(B0**2/4.0 - B0OLD**2/4.0)
          B0OLD = B0
        ELSE
          RANDPLUMSig = 0.0
          B0OLD = -999.0
        ENDIF
      Else
        RANDPLUMSig = 0.0
        B0OLD = -999.0
      ENDIF

!      IF (DEBUG) THEN
!        WRITE(83, '( F7.2, 4X, 5(F6.2, 2X), E12.4, 2X, 2(F8.2, 2X), 5(F6.2, 2X), &
!        2(F8.6, 2X), 3(E12.6, 2X), 21(E12.4, 2X))')                               &
!        IPARTAGE, UP, VP, WP,                                                     &
!        U, V, W, B, B0, UZETA, TP, T, THP, TH, RHOA, RHOP,                        &
!        RANDPLUMSig,                                                              &
!        FMASS, FMX, FMY, FMZ, FH, FMASS0,                                         &
!        DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DWDX, DWDY, DWDZ, DWDX, DWDY, DWDZ,   &
!        DTHETADX, DTHETADY, DTHETADZ
!      ENDIF


! Flag to stop plumerise
      IF ((FMZ/FMASS).LT.0.01.OR.IPARTAGE.GT.3600.0) THEN
!      IF ((FMZ/FMASS).LT.0.5.OR.IPARTAGE.GT.3600.0) THEN
        FMASS = 0.0
!        IF (DEBUG) THEN
!          WRITE(84, '(A30, I4)')'PLUME RISE TERMINATED AT ', IPARTAGE
!        ENDIF
      ENDIF

End Subroutine PlumeRise

!-------------------------------------------------------------------------------------------------------------

End Module PlumeRiseModule

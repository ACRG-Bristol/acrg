! Module:  Puff Module

Module PuffModule

! This module provides code to treat puffs.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowsModule, Only: Flows_, Flow_, Cloud_, Rain_
Use SpeciesModule
Use SourceModule
Use OutputModule, Only: OutputOpts_, Reqs_, Results_, OutputPPInfos
Use PlumeRiseModule
Use ParticleModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: A6
Public :: PuffOpts_
Public :: Puff_
Public :: AssignPuff
Public :: GetNewPuffIndex
Public :: GetNewOriginalPuffIndex
Public :: InactivatePuff
Public :: NextPuff
Public :: InitPuff
Public :: InitPuff2
Public :: Step2
Public :: Step3
Public :: Step4
Public :: Step5
Public :: Step6
Public :: Step7
Public :: PuffFraction
Public :: PuffConc
Public :: PuffCentres
Public :: PuffMass
Public :: PuffMeanZ
Public :: PuffSigZ2
Public :: PuffInfo

!-------------------------------------------------------------------------------------------------------------

Real(Std), Parameter :: A2 = 1.2_Std ! Large values make puff-recombination more
                                     ! difficult.
Real(Std), Parameter :: A3 = 1.0     ! Large values reduce error in Gaussian
                                     ! approximation.
Real(Std), Parameter :: A4 = 0.8_Std ! In (x,u) Markov model we aim to keep Sigma_p
                                     ! between A4*Delta and Delta.
Real(Std), Parameter :: A6 = 0.5     ! Multiple of sigma_h used in limiting Delta.

!-------------------------------------------------------------------------------------------------------------

Type :: PuffOpts_ ! Puff options.
  Real(Std) :: A1 ! Large values reduce statistical noise.
  Real(Std) :: A5 ! How far into tail of distribution in calculating concentration.
  Real(Std) :: A7 ! Critical distance moved/sigma_x^2 for deciding whether to use error function or Gaussian
                  ! approx.
End Type PuffOpts_

!-------------------------------------------------------------------------------------------------------------

Type :: Puff_ ! Information defining a puff.
  ! If changing this type, remember restart read and write routines.
  Type(Particle_)  :: P            ! $$ note XOld not yet used.
  Integer          :: iUOP         ! } Unique index of original puff (i.e. an index unique to the
                                   ! } original puff which is
                                   ! not reused despite puff recycling). $$ Integer(8)
  Real(Std)        :: XXh(3)       !} Homogeneous spread.
  Real(Std)        :: XUh(3)       !}
  Real(Std)        :: UUh(3)       !}
  Real(Std)        :: XXp(3)       !] Puff spread.
  Real(Std)        :: XUp(3)       !]
  Real(Std)        :: UUp(3)       !]
  Real(Std)        :: XXpOld(3)    ! Puff spread at end of previous time step.
  Real(Std)        :: SigT         !
  Integer          :: N            ! Number of times puff has been split.
  Integer          :: Sibling      ! Index of sibling puff; 0 = no sibling puff.
  Integer          :: Parent       ! Index of parent puff; 0 = no parent puff.
  Integer          :: Child1       ! Index of first child puff; 0 = no child puff.
  Integer          :: Child2       ! Index of second child puff; 0 = no child puff.
  Integer          :: OriginalPuff ! Index of original puff
  Logical          :: G            ! Indicates puff is Gaussian.
  Logical          :: R            ! Indicates puff position may have a random component.
  Real(Std)        :: dXdT(3)      ! rate of change of puff centroid with time as time
                                   ! varies across the time extent of the puff.
End Type Puff_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine AssignPuff(Puff1, Puff2)
! Set Puff2 = Puff1 except for relational values.

  Implicit None
  ! Argument list:
  Type(Puff_), Intent(In)    :: Puff1 !
  Type(Puff_), Intent(InOut) :: Puff2 !
  ! Locals:
  Integer :: Sibling ! Number of sibling puff; 0 = no sibling puff.
  Integer :: Parent  ! Number of parent puff; 0 = no parent puff.
  Integer :: Child1  ! Number of first child puff; 0 = no child puff.
  Integer :: Child2  ! Number of second child puff; 0 = no child puff.
  Integer :: iUP

  iUP           = Puff2%P%iUP
  Sibling       = Puff2%Sibling
  Parent        = Puff2%Parent
  Child1        = Puff2%Child1
  Child2        = Puff2%Child2
  Puff2         = Puff1
  Puff2%P%iUP   = iUP
  Puff2%Sibling = Sibling
  Puff2%Parent  = Parent
  Puff2%Child1  = Child1
  Puff2%Child2  = Child2

End Subroutine AssignPuff

!-------------------------------------------------------------------------------------------------------------

Subroutine GetNewPuffIndex(iUP, iP, nPuffs, LastPuff, FreePuffStack, iULastPuff)
!

  Implicit None
  ! Argument list:
  Integer, Intent(Out)   :: iUP              ! Unique index of puff (i.e. an index unique to the puff which is
                                             ! not reused despite puff recycling). $$ Integer(8)
  Integer, Intent(Out)   :: iP               ! Puff index.
  Integer, Intent(InOut) :: nPuffs           !
  Integer, Intent(InOut) :: LastPuff         !
  Integer, Intent(InOut) :: FreePuffStack(:) !
  Integer, Intent(InOut) :: iULastPuff       !

  iULastPuff = iULastPuff + 1
  iUP        = iULastPuff

  If (nPuffs + 1 > Size(FreePuffStack)) Then
    Call Message('FATAL ERROR: too many puffs', 3)
  End If

  nPuffs = nPuffs + 1
  iP = FreePuffStack(nPuffs)
  If (iP > LastPuff) LastPuff = iP

End Subroutine GetNewPuffIndex

!-------------------------------------------------------------------------------------------------------------

Subroutine GetNewOriginalPuffIndex(                                                       &
             iUOP, iUP, iOP, iP,                                                          &
             nOriginalPuffs, LastOriginalPuff, FreeOriginalPuffStack, iULastOriginalPuff, &
             nPuffs, LastPuff, FreePuffStack, iULastPuff,                                 &
             iPuffs                                                                       &
           )
!

  Implicit None
  ! Argument list:
  Integer, Intent(Out)   :: iUOP                     ! Unique index of original puff (i.e. an index unique
                                                     ! to the
                                                     ! original puff which is not reused despite puff
                                                     ! recycling). $$ Integer(8)
  Integer, Intent(Out)   :: iUP                      ! Unique index of puff (i.e. an index unique to the puff
                                                     ! which is not reused despite puff recycling).
                                                     ! $$ Integer(8)
  Integer, Intent(Out)   :: iOP                      !
  Integer, Intent(Out)   :: iP                       !
  Integer, Intent(InOut) :: nOriginalPuffs           !
  Integer, Intent(InOut) :: LastOriginalPuff         !
  Integer, Intent(InOut) :: FreeOriginalPuffStack(:) !
  Integer, Intent(InOut) :: iULastOriginalPuff       !
  Integer, Intent(InOut) :: nPuffs                   !
  Integer, Intent(InOut) :: LastPuff                 !
  Integer, Intent(InOut) :: FreePuffStack(:)         !
  Integer, Intent(InOut) :: iULastPuff               !
  Integer, Intent(InOut) :: iPuffs(:)                !

  iULastOriginalPuff = iULastOriginalPuff + 1
  iUOP               = iULastOriginalPuff

  If (nOriginalPuffs + 1 > Size(FreeOriginalPuffStack)) Then
    Call Message('FATAL ERROR: too many original puffs', 3)
  End If

  nOriginalPuffs = nOriginalPuffs + 1
  iOP = FreeOriginalPuffStack(nOriginalPuffs)
  If (iOP > LastOriginalPuff) LastOriginalPuff = iOP

  iULastPuff = iULastPuff + 1
  iUP        = iULastPuff

  If (nPuffs + 1 > Size(FreePuffStack)) Then
    Call Message('FATAL ERROR: too many puffs', 3)
  End If

  nPuffs = nPuffs + 1
  iP = FreePuffStack(nPuffs)
  If (iP > LastPuff) LastPuff = iP

  iPuffs(iOP) = iP

End Subroutine GetNewOriginalPuffIndex

!-------------------------------------------------------------------------------------------------------------

Subroutine InactivatePuff(                                            &
             iP,                                                      &
             Puffs,                                                   &
             nPuffs, LastPuff, FreePuffStack,                         &
             nOriginalPuffs, LastOriginalPuff, FreeOriginalPuffStack, &
             iPuffs                                                   &
           )
! Inactivates a puff, and propagates inactivation to parent puffs until (i) a puff which is active is found or
! (ii) a puff with an active child puff is found. (i) is used when combining puffs. By making the parent puff
! active before calling this routine, the propagation is prevented.

  Implicit None
  ! Argument list:
  Integer,     Intent(In)    :: iP                       !
  Type(Puff_), Intent(InOut) :: Puffs(:)                 !
  Integer,     Intent(InOut) :: nPuffs                   !
  Integer,     Intent(InOut) :: LastPuff                 !
  Integer,     Intent(InOut) :: FreePuffStack(:)         !
  Integer,     Intent(InOut) :: nOriginalPuffs           !
  Integer,     Intent(InOut) :: LastOriginalPuff         !
  Integer,     Intent(InOut) :: FreeOriginalPuffStack(:) !
  Integer,     Intent(InOut) :: iPuffs(:)                !
  ! Locals:
  Integer :: iPuff
  Integer :: iParent

  ! Inactivate puff.
  Puffs(iP)%P%Active = .false.

  ! Set indices for puff and parent puff.
  iPuff   = iP
  iParent = Puffs(iPuff)%Parent

  Do

    ! Return puff to the stack of free puffs and reset number of active puffs.
    FreePuffStack(nPuffs) = iPuff
    nPuffs = nPuffs - 1

    ! If no parent puff, return original puff to the stack of free original puffs and
    ! reset number of active original puffs.
    If (iParent == 0) Then
      FreeOriginalPuffStack(nOriginalPuffs) = Puffs(iPuff)%OriginalPuff
      nOriginalPuffs = nOriginalPuffs - 1
      iPuffs(Puffs(iPuff)%OriginalPuff) = 0
      Return
    End If

    ! Set child index in parent puff to 0.
    If (Puffs(iParent)%Child1 == iPuff) Then
      Puffs(iParent)%Child1 = 0
    Else If (Puffs(iParent)%Child2 == iPuff) Then
      Puffs(iParent)%Child2 = 0
    Else
      Call Message('UNEXPECTED FATAL ERROR in InactivatePuff', 4)
    End If

    ! If parent puff active or has an active child puff, return.
    If (                                    &
      ParticleActive(Puffs(iParent)%P) .or. &
      Puffs(iParent)%Child1 /= 0       .or. &
      Puffs(iParent)%Child2 /= 0            &
    ) Return

    ! Set indices for puff and parent puff.
    iPuff   = iParent
    iParent = Puffs(iPuff)%Parent

  End Do

End Subroutine InactivatePuff

!-------------------------------------------------------------------------------------------------------------

Subroutine NextPuff(LastOriginalPuff, Puffs, iPuffs, iP, iOP, NoMorePuffs)

  Implicit None
  ! Argument list:
  Integer,     Intent(In)    :: LastOriginalPuff
  Type(Puff_), Intent(In)    :: Puffs(:)    !
  Integer,     Intent(In)    :: iPuffs(:)
  Integer,     Intent(InOut) :: iP          !} Set to zero on first call
  Integer,     Intent(InOut) :: iOP         !}
  Logical,     Intent(Out)   :: NoMorePuffs !
  ! Locals:
  Integer :: iPLast

  NoMorePuffs = .false.

  ! Find next puff down the tree (if possible).
  If (iP == 0) Then

    Do
      iOP = iOP + 1
      If (iOP > LastOriginalPuff) Then
        NoMorePuffs = .true.
        Return
      End If
      iP = iPuffs(iOP)
      If (iP /= 0) Return
    End Do

  Else If (Puffs(iP)%Child1 /= 0) Then

    iP = Puffs(iP)%Child1
    Return

  Else If (Puffs(iP)%Child2 /= 0) Then

    iP = Puffs(iP)%Child2
    Return

  End If

  ! If no more puffs down the tree, back up and look for new branch.
  Do

    ! Note last value of iP.
    iPLast = iP

    ! Move one step up the tree.
    iP = Puffs(iPLast)%Parent

    ! Find next puff down the tree (if possible) taking care to exclude puffs already visited.
    If (iP == 0) Then

      Do
        iOP = iOP + 1
        If (iOP > LastOriginalPuff) Then
          NoMorePuffs = .true.
          Return
        End If
        iP = iPuffs(iOP)
        If (iP /= 0) Return
      End Do

    Else If (Puffs(iP)%Child2 /= 0 .and. Puffs(iP)%Child1 == iPLast) Then

      iP = Puffs(iP)%Child2
      Return

    End If

  End Do

End Subroutine NextPuff

!-------------------------------------------------------------------------------------------------------------

Subroutine InitPuff(                                         &
             iUP,                                            &
             TimeM, Time, TimeP, Instantaneous,              &
             SkewTime, VelMemTime, InhomogTime, MVelMemTime, &
             Coords, Specieses,                              &
             iSource, Source, H3,                            &
             iUOP, iPuff,                                    &
             Puff, Extra                                     &
           )
! Initialise a puff.

  Implicit None
  ! Argument list:
  Integer,          Intent(In)  :: iUP           !
  Type(ShortTime_), Intent(In)  :: TimeM         !
  Type(ShortTime_), Intent(In)  :: Time          !
  Type(ShortTime_), Intent(In)  :: TimeP         !
  Logical,          Intent(In)  :: Instantaneous !
  Type(ShortTime_), Intent(In)  :: SkewTime      !
  Type(ShortTime_), Intent(In)  :: VelMemTime    !
  Type(ShortTime_), Intent(In)  :: InhomogTime   !
  Type(ShortTime_), Intent(In)  :: MVelMemTime   !
  Type(Coords_),    Intent(In)  :: Coords        !
  Type(Specieses_), Intent(In)  :: Specieses     !
  Integer,          Intent(In)  :: iSource       !
  Type(Source_),    Intent(In)  :: Source        !
  Real(Std),        Intent(In)  :: H3            !
  Integer,          Intent(In)  :: iUOP          !
  Integer,          Intent(In)  :: iPuff         !
  Type(Puff_),      Intent(Out) :: Puff          !
  Type(Extra_),     Intent(Out) :: Extra         !
  ! Locals:
  Real(Std) :: HMax            ! Max metric coefficient.
  Real(Std) :: H(3)   !
  Real(Std) :: Factor !
  Integer   :: SourceDim       ! Number of elements of dX(:) which are non-zero.
  Real(Std) :: GeometricFactor !
  Integer   :: i      !

  Puff%P%iUP = iUP
  Puff%iUOP  = iUOP

  Puff%P%iSource = iSource

  Extra%Skew    = SkewTime    > ZeroShortTime()
  Extra%VelMem  = VelMemTime  > ZeroShortTime()
  Extra%Inhomog = InhomogTime > ZeroShortTime()
  Extra%MVelMem = MVelMemTime > ZeroShortTime()

  Puff%P%iHCoord = Source%iHCoord
  Puff%P%iZCoord = Source%iZCoord

  Call MetricCoeffs(                     &
         Coords%HCoords(Puff%P%iHCoord), &
         Source%X(1:2),                  &
         HMax, H(1), H(2)                &
       )
  H(3) = H3

  Puff%P%X(1)  = Source%X(1)
  Puff%P%X(2)  = Source%X(2)
  Puff%P%X(3)  = Source%X(3)
  Extra%U(1:3) = 0.0

  Puff%P%XOld(:)    = Puff%P%X(:)
  Puff%P%iHCoordOld = Puff%P%iHCoord
  Puff%P%iZCoordOld = Puff%P%iZCoord

  If ( Source%Shape .CIEq. 'Ellipsoid' ) Then
    SourceDim = 0
    Do i = 1, 3
      If (Source%dX(i) > 0.0) SourceDim = SourceDim + 1
    End Do
    GeometricFactor = Sqrt( 1.0 / ( 4.0 * (Real(SourceDim) + 2.0) ) )
    Do i = 1, 3
      Puff%XXh(i) = (GeometricFactor * Source%dX(i)) ** 2 ! Off diag bits? $$
      Puff%XUh(i) = 0.0
      Puff%UUh(i) = 0.0   !$$
      Puff%XXp(i) = (GeometricFactor * Source%dX(i)) ** 2
      Puff%XUp(i) = 0.0
      Puff%UUp(i) = 0.0   !$$
    End Do

    Puff%SigT   = 0.0
  Else
    Do i = 1, 3
      Puff%XXh(i) = (1.0/12.0) * Source%dX(i) ** 2 ! Off diag bits? $$
      Puff%XUh(i) = 0.0
      Puff%UUh(i) = 0.0                                     !$$
      Puff%XXp(i) = (1.0/12.0) * Source%dX(i) ** 2
      Puff%XUp(i) = 0.0
      Puff%UUp(i) = 0.0                                     !$$
    End Do

    Puff%SigT   = Sqrt(1.0/12.0) * Source%dX(3) ! $$
  End If
  If (.not.Source%dHMetres) Then
    Do i = 1, 2
      Puff%XXh(i) = Puff%XXh(i) * H(i) ** 2
      Puff%XXp(i) = Puff%XXp(i) * H(i) ** 2
    End Do
  End If
  If (.not.Source%dZMetres) Then
    Puff%XXh(3) = Puff%XXh(3) * H(3) ** 2
    Puff%XXp(3) = Puff%XXp(3) * H(3) ** 2
  End If

  Puff%dXdT(:) = 0.0

  Puff%XXpOld(:) = Puff%XXp(:)

  Puff%P%T     = ZeroShortTime()
  Puff%P%TOld  = Puff%P%T
  Puff%P%T0    = Time
  Extra%TPlus  = TimeP - Puff%P%T0
  Extra%TMinus = TimeM - Puff%P%T0
  If (Instantaneous) Then
    Extra%TInst = InfPastShortTime(Interval = .true.)
  Else
    Extra%TInst = InfFutureShortTime(Interval = .true.)
  End If

  Puff%N            = 0
  Puff%Sibling      = 0
  Puff%Parent       = 0
  Puff%Child1       = 0
  Puff%Child2       = 0
  Puff%G            = .not.Source%TopHat
  Puff%R            = .false.
  Puff%OriginalPuff = iPuff

  Puff%P%Active       = .true.
  Puff%P%Marked       = 0
  Puff%P%NeedToSetVel = .false.

  If (Source%PlumeRise) Then
    Extra%FMass  = Source%VolumeFlowRate
    Extra%FM(1)  = 0.0
    Extra%FM(2)  = 0.0
    Extra%FM(3)  = Source%FlowVelocity
    Extra%FH     = Source%Temperature
    Extra%FMass0 = Source%VolumeFlowRate
    Extra%B0Old  = Source%dX(1)/2.0
  Else
    Extra%FMass  = 0.0
  End If

  Puff%P%WSed     = 0.0 ! $$
  Puff%P%Diameter = 0.0 ! $$
  Puff%P%Density  = 0.0 ! $$

End Subroutine InitPuff

!-------------------------------------------------------------------------------------------------------------

Subroutine InitPuff2(DeltaOpt, Flow, Puff, Extra, iP)
! .

  Implicit None
  ! Argument list:
  Integer,      Intent(In)    :: DeltaOpt !
  Type(Flow_),  Intent(In)    :: Flow     !
  Type(Puff_),  Intent(InOut) :: Puff     !
  Type(Extra_), Intent(InOut) :: Extra    !
  Integer,      Intent(In)    :: iP       ! Puff index
  ! Locals:
  Real(Std) :: R       !
  Real(Std) :: Rho_p   !
  Real(Std) :: Theta_p !
  Real(Std) :: Beta
  Real(Std) :: Delta

  If (Extra%VelMem) Then

    If (DeltaOpt == 0) Then
      Delta = Huge(Delta)
    Else
      Delta = Flow%DeltaI
    End If
    If (Puff%XXp(3) >= Delta**2) Then
      Beta = 1.0
    Else If (DeltaOpt == 2) Then
      Beta = 1.0 - A6**2
    Else
      Beta = 0.0
    End If
    Extra%U(3) = InitW(Flow, Beta, Extra%Skew, iP)

    Puff%UUp(1) = Flow%SigUU(1)
    Puff%UUp(2) = Flow%SigUU(2)
    Puff%UUp(3) = Flow%SigUU(3) * (1.0 - Beta)
    Puff%UUh(1) = Flow%SigUU(1)
    Puff%UUh(2) = Flow%SigUU(2)
    Puff%UUh(3) = Flow%SigUU(3)

    Extra%SigUU(1)  = Flow%SigUU(1)
    Extra%SigUU(2)  = Flow%SigUU(2)
    Extra%SigUU(3)  = Flow%SigUU(3)
    Extra%UAOLD(1)  = Flow%U(1)
    Extra%UAOLD(2)  = Flow%U(2)
    Extra%UAOLD(3)  = Flow%U(3)
    Extra%TAOLD     = Flow%T
    Extra%ThetaAOLD = Flow%Theta
    Extra%PAOLD     = Flow%P

  End If

  If (Extra%FMass /= 0.0) Then ! test on real dangerous - hence else block
    ! Plume density and potential temperature.
    Rho_p   = Flow%P/(GasConstant * Extra%FH)
    Theta_p = Extra%FH * (PRef/Flow%P)**(GasConstant/Cp)
    ! Plume fluxes.
    Extra%FMass  = Rho_p * Extra%FMass
    Extra%FM(1)  = (Extra%FM(1) - Flow%U(1)) * Extra%FMass
    Extra%FM(2)  = (Extra%FM(2) - Flow%U(2)) * Extra%FMass
    Extra%FM(3)  = (Extra%FM(3) - Flow%U(3)) * Extra%FMass
    Extra%FH     = Cp * (Theta_p - Flow%Theta) * Extra%FMass
    Extra%FMass0 = Rho_p * Extra%FMass0
  Else
    Extra%FMass  = 0.0
    Extra%FM(1)  = 0.0
    Extra%FM(2)  = 0.0
    Extra%FM(3)  = 0.0
    Extra%FH     = 0.0
    Extra%FMass0 = 0.0
  End If

  Puff%P%H = Flow%H

End Subroutine InitPuff2

!-------------------------------------------------------------------------------------------------------------

Subroutine Step2(LastPuff, nOriginalPuffs, Puffs, Masses, SigA2, T)
! Compute sigma_a.

  Implicit None
  ! Argument list:
  Integer,          Intent(In)  :: LastPuff
  Integer,          Intent(In)  :: nOriginalPuffs !
  Type(Puff_),      Intent(In)  :: Puffs(:)
  Real(Std),        Intent(In)  :: Masses(:, :)
  Real(Std),        Intent(Out) :: SigA2(:)       ! size MaxOriginalPuffs
  Type(ShortTime_), Intent(Out) :: T(:)
  ! Locals:
  Integer   :: i               !
  Integer   :: iOP             !
  Real(Std) :: S0(Size(SigA2)) !
  Real(Std) :: S1(Size(SigA2)) !
  Real(Std) :: S2(Size(SigA2)) !
  Real(Std) :: Temp            !

  S0(1:nOriginalPuffs) = 0.0
  S1(1:nOriginalPuffs) = 0.0
  S2(1:nOriginalPuffs) = 0.0

  Do i = 1, LastPuff

    iOP = Puffs(i)%OriginalPuff

    If (.not.Puffs(i)%P%Active) Cycle
    If (Puffs(i)%XXp(3) == 0.0) Then
      Temp = Puffs(i)%P%X(3)
    Else
      Temp    = Puffs(i)%P%X(3) / Sqrt(2.0*Puffs(i)%XXp(3))
      Temp    = Puffs(i)%P%X(3)*Erf(Temp) + Sqrt(2.0*Puffs(i)%XXp(3)/Pi)*Exp(-Temp**2)
    End If
    S0(iOP) = S0(iOP) + Masses(1, i) / 2**Puffs(i)%N
    S1(iOP) = S1(iOP) + Temp * Masses(1, i) / 2**Puffs(i)%N
    S2(iOP) = S2(iOP) + (Puffs(i)%XXp(3) + Puffs(i)%P%X(3)**2) * &
                        Masses(1, i) / 2**Puffs(i)%N ! $$ note use of species 1 here
    T(iOP) = ParticleTime(Puffs(i)%P)

  End Do

  Do iOP = 1, nOriginalPuffs

    If (S0(iOP) /= 0.0) Then
      S1(iOP)    = S1(iOP)/S0(iOP)
      S2(iOP)    = S2(iOP)/S0(iOP)
      SigA2(iOP) = S2(iOP) - S1(iOP)**2
    End If

  End Do

End Subroutine Step2

!-------------------------------------------------------------------------------------------------------------

Subroutine Step3(                                                     &
             PuffOpts,                                                &
             SigA2, Puffs, Extras, Masses,                            &
             nPuffs, LastPuff, FreePuffStack,                         &
             nOriginalPuffs, LastOriginalPuff, FreeOriginalPuffStack, &
             iPuffs                                                   &
           )
! Re-combine puffs.

  Implicit None
  ! Argument list:
  Type(PuffOpts_), Intent(In)    :: PuffOpts
  Real(Std),       Intent(In)    :: SigA2(:)                 ! size MaxOriginalPuffs
  Type(Puff_),     Intent(InOut) :: Puffs(:)                 !
  Type(Extra_),    Intent(InOut) :: Extras(:)                !
  Real(Std),       Intent(InOut) :: Masses(:, :)
  Integer,         Intent(InOut) :: nPuffs                   !
  Integer,         Intent(InOut) :: LastPuff                 !
  Integer,         Intent(InOut) :: FreePuffStack(:)         !
  Integer,         Intent(InOut) :: nOriginalPuffs           !
  Integer,         Intent(InOut) :: LastOriginalPuff         !
  Integer,         Intent(InOut) :: FreeOriginalPuffStack(:) !
  Integer,         Intent(InOut) :: iPuffs(:)                !
  ! Locals:
  Integer   :: i    !
  Integer   :: j    !
  Integer   :: k    !
  Real(Std) :: R    !
  Real(Std) :: Temp !

  Do i = 1, Size(Puffs)       !$$ What if one is inactive but other active? May still
    If (i > LastPuff) Exit ! wish to combine.
    If (.not.Puffs(i)%P%Active) Cycle
    If (Puffs(i)%Sibling == 0) Cycle
    If (.not.Puffs(Puffs(i)%Sibling)%P%Active) Cycle

    Temp = (A2 * PuffOpts%A1 / 2.0**(Puffs(i)%N - 1))**2 * &
           Min(                                            &
             Puffs(i)%XXh(3),                              &
             SigA2(Puffs(i)%OriginalPuff)                  &
           )
    If (Puffs(i)%G .and. Puffs(i)%R .and. Puffs(i)%XXp(3) > Temp) Then
      j = Puffs(i)%Sibling
      k = Puffs(i)%Parent
      If (j > i) Then
        Temp = (A2 * PuffOpts%A1 / 2.0**(Puffs(j)%N - 1))**2 * &
               Min(                                            &
                 Puffs(j)%XXh(3),                              &
                 SigA2(Puffs(i)%OriginalPuff)                  &
               )
        If (Puffs(j)%G .and. Puffs(j)%R .and. Puffs(j)%XXp(3) > Temp) Then
          Call GetRandomNumber(R, i)
          If (R < 0.5) Then
            Call AssignPuff(Puffs(i), Puffs(k))
            Masses(:, k) = Masses(:, i)
            Extras(   k) = Extras(   i)
          Else
            Call AssignPuff(Puffs(j), Puffs(k))
            Masses(:, k) = Masses(:, j)
            Extras(   k) = Extras(   j)
          End If
          Puffs(k)%N      = Puffs(k)%N - 1
          Puffs(k)%P%Active = .true.

          Call InactivatePuff(                                            &
                 i,                                                       &
                 Puffs,                                                   &
                 nPuffs, LastPuff, FreePuffStack,                         &
                 nOriginalPuffs, LastOriginalPuff, FreeOriginalPuffStack, &
                 iPuffs                                                   &
               )
          Call InactivatePuff(                                            &
                 j,                                                       &
                 Puffs,                                                   &
                 nPuffs, LastPuff, FreePuffStack,                         &
                 nOriginalPuffs, LastOriginalPuff, FreeOriginalPuffStack, &
                 iPuffs                                                   &
               )

        End If
      End If
    End If

  End Do

End Subroutine Step3

!-------------------------------------------------------------------------------------------------------------

Subroutine Step4(Puff)
! Approximate by Gaussian.

  Implicit None
  ! Argument list:
  Type(Puff_), Intent(InOut) :: Puff !

!  If (.not.Puff%G .and. Puff%XXp(3) > (A3 * Puff%SigT)**2) Then ! $$ currently puff scheme doesn't
!  treat top hat source
    Puff%G    = .true.
    Puff%SigT = 0.0
!  End If

End Subroutine Step4

!-------------------------------------------------------------------------------------------------------------

Subroutine Step5(PuffOpts, DeltaOpt, Delta, Flow, SigA2, SigPGtDelta, PuffMoved, Puff, Extra, iP)
! Reduce sigma_p.

  Implicit None
  ! Argument list:
  Type(PuffOpts_), Intent(In)    :: PuffOpts
  Integer,         Intent(In)    :: DeltaOpt    !
  Real(Std),       Intent(In)    :: Delta       !
  Type(Flow_),     Intent(In)    :: Flow
  Real(Std),       Intent(In)    :: SigA2       !
  Logical,         Intent(Out)   :: SigPGtDelta !
  Logical,         Intent(Out)   :: PuffMoved   !
  Type(Puff_),     Intent(InOut) :: Puff        !
  Type(Extra_),    Intent(InOut) :: Extra       !
  Integer,         Intent(In)    :: iP          ! Puff index
  ! Locals:
  Real(Std) :: Lambda2 !
  Real(Std) :: Temp    !
  Real(Std) :: Corr    !
  Real(Std) :: R1      !
  Real(Std) :: R2      !
  Real(Std) :: DeltaL  !
  Real(Std) :: Beta

 ! $$ note commenting this out degrades in3.txt comparison with particles. But its more logical to
 ! $$ comment out.
 ! Could probably adjust in3.txt to agree better, or introduce a 'split deterministically at source (into n)
 ! even if not otherwise required' option (n being user determined) - could help with tophat issue too.
 ! If (DeltaOpt == 2 .and. TravelTime(Puff%P) == ZeroShortTime()) Then
 !   DeltaL = Min(Delta, A6 * Sqrt(Puff%XXh(3)))
 ! Else
    DeltaL = Delta
 ! End If

  SigPGtDelta = Puff%XXp(3) > DeltaL**2
  Temp = (PuffOpts%A1 / 2.0**Puff%N)**2 * Min(Puff%XXh(3), SigA2)
  If (Puff%R .and. SigPGtDelta .and. Puff%XXp(3) > Temp) Then
    If (Temp <= DeltaL**2) Then
      Lambda2     = DeltaL**2/Puff%XXp(3)
      SigPGtDelta = .false.
    Else
      Lambda2     = Temp/Puff%XXp(3)
      SigPGtDelta = .true.
    End If
    R1          = Gauss(iP)
    Puff%P%X(3) = Puff%P%X(3) + Sqrt((1.0 - Lambda2)*Puff%XXp(3))*R1
    If (Extra%VelMem) Then
      R2   = Gauss(iP)

      Corr = Sqrt(Puff%XXp(3)*Puff%UUp(3))
      If (Corr == 0.0) Then
        Corr = 1.0
      Else
        Corr = Puff%XUp(3)/Sqrt(Puff%XXp(3)*Puff%UUp(3))
      End If
      Corr = Min(Corr, 1.0)

!      Corr = Puff%XUp(3)/Sqrt(Puff%XXp(3)*Puff%UUp(3))

      Extra%U(3) = Extra%U(3) +                      &
                   Sqrt(                             &
                     (1.0 - Lambda2)*Puff%UUp(3)     &
                   )                                 &
                   *(Corr*R1 + Sqrt(Max(1.0 - Corr**2, 0.0))*R2)

      Puff%XUp(3) = Puff%XUp(3)*Lambda2
      Puff%UUp(3) = Puff%UUp(3)*Lambda2

     ! If (TravelTime(Puff%P) == ZeroShortTime()) Then
     !  ! $$ note splitting large sources on release is a situation where recalling flow would be worthwhile.
     !  ! $$ currently reinitialisation only ensures random numbers different.
     !  ! $$ should consider reinitialisation when deterministic splitting used too.
     !   If (Puff%XXp(3) >= Delta**2) Then ! $$ Only first option should be relevant here.
     !     Beta = 1.0 ! - A6**2
     !   Else If (DeltaOpt == 2) Then
     !     Beta = 1.0 - A6**2
     !   Else
     !     Beta = 0.0
     !   End If
     !   Extra%U(3) = InitW(Flow, Beta, Extra%Skew, iP)
     !   Puff%UUp(3) = Flow%SigUU(3) * (1.0 - Beta)
    !  Else
    !    Puff%XXp(3) = Puff%XXp(3)
     ! End If

    End If

    Puff%XXp(3) = Puff%XXp(3)*Lambda2 ! $$ Should ensure hitting delta**2 (no precision error)

    PuffMoved = .true.
  Else
    PuffMoved = .false.
  End If

  If (Extra%VelMem) Then
    If (TravelTime(Puff%P) == ZeroShortTime()) Then
     ! $$ note splitting large sources on release is a situation where recalling flow would be worthwhile.
     ! $$ currently reinitialisation only ensures random numbers different.
      If (.not. Puff%R) Then ! $$ shouldn't be needed except through rounding errors - ensure isn't needed
        Beta = 0.0
      Else If (Puff%XXp(3) >= Delta**2) Then
        Beta = 1.0 ! - A6**2
      Else If (DeltaOpt == 2) Then
        Beta = 1.0 - A6**2
      Else
        Beta = 0.0
      End If
      Extra%U(3) = InitW(Flow, Beta, Extra%Skew, iP)
      Puff%UUp(3) = Flow%SigUU(3) * (1.0 - Beta)
    End If
  End If

End Subroutine Step5

!-------------------------------------------------------------------------------------------------------------

Subroutine Step6(                                        &
             PuffOpts,                                   &
             i, DeltaOpt, Delta, SigA2, SigPGtDelta,     &
             Puffs, Extras, Masses,                      &
             nPuffs, LastPuff, FreePuffStack, iULastPuff &
           )
! Split puffs.

  Implicit None
  ! Argument list:
  Type(PuffOpts_), Intent(In)    :: PuffOpts
  Integer,         Intent(In)    :: i                !
  Integer,         Intent(In)    :: DeltaOpt         !
  Real(Std),       Intent(In)    :: Delta            ! $$ not used
  Real(Std),       Intent(In)    :: SigA2            !
  Logical,         Intent(In)    :: SigPGtDelta      !
  Type(Puff_),     Intent(InOut) :: Puffs(:)         !
  Type(Extra_),    Intent(InOut) :: Extras(:)        !
  Real(Std),       Intent(InOut) :: Masses(:, :)     !
  Integer,         Intent(InOut) :: nPuffs           !
  Integer,         Intent(InOut) :: LastPuff         !
  Integer,         Intent(InOut) :: FreePuffStack(:) !
  Integer,         Intent(InOut) :: iULastPuff       !
  ! Locals:
  Integer   :: Child1 !
  Integer   :: Child2 !
  Integer   :: iUP    !
  Real(Std) :: Temp   !

  Temp = (PuffOpts%A1 / 2.0**Puffs(i)%N)**2 * Min(Puffs(i)%XXh(3), SigA2)
  If (                                                        &
    SigPGtDelta                                          .or. &
    (Puffs(i)%XXp(3) < Temp      .and.       Puffs(i)%R) .or. &
 !   (Puffs(i)%XXp(3) >= Delta**2 .and. .not. Puffs(i)%R) .or. &
    (DeltaOpt == 2               .and. .not. Puffs(i)%R)      &
  ) Then

    If (Puffs(i)%Child1 == 0) Then
      Call GetNewPuffIndex(iUP, Child1, nPuffs, LastPuff, FreePuffStack, iULastPuff)
      Puffs(i     )%Child1 = Child1
      Puffs(Child1)%Parent = i
      Puffs(Child1)%Child1 = 0
      Puffs(Child1)%Child2 = 0
      Puffs(Child1)%P%iUP  = iUP
    Else
      Child1 = Puffs(i)%Child1
    End If
    If (Puffs(i)%Child2 == 0) Then
      Call GetNewPuffIndex(iUP, Child2, nPuffs, LastPuff, FreePuffStack, iULastPuff)
      Puffs(i     )%Child2 = Child2
      Puffs(Child2)%Parent = i
      Puffs(Child2)%Child1 = 0
      Puffs(Child2)%Child2 = 0
      Puffs(Child2)%P%iUP  = iUP
    Else
      Child2 = Puffs(i)%Child2
    End If
    Puffs(Child1)%Sibling = Child2
    Puffs(Child2)%Sibling = Child1

    If (Puffs(i)%G) Then
      Call AssignPuff(Puffs(i), Puffs(Child1))
      Masses(:, Child1) = Masses(:, i)
      Extras(   Child1) = Extras(   i)
      Puffs(Child1)%N = Puffs(Child1)%N + 1
      Puffs(Child1)%R = .true.
      Call AssignPuff(Puffs(Child1), Puffs(Child2))
      Masses(:, Child2) = Masses(:, Child1)
      Extras(   Child2) = Extras(   Child1)
    Else
      Call AssignPuff(Puffs(i), Puffs(Child1))
      Masses(:, Child1) = Masses(:, i)
      Extras(   Child1) = Extras(   i)
      Extras(Child1)%U(3) = 0.0 ! $$ ? This looks odd
      Puffs(Child1)%XXp(3)  = Puffs(Child1)%XXp(3) - 0.75*Puffs(Child1)%SigT**2
      Puffs(Child1)%N    = Puffs(Child1)%N + 1
      Puffs(Child1)%R    = .false.
      Puffs(Child1)%G    = .false.
      Puffs(Child1)%SigT = 0.5*Puffs(Child1)%SigT
      Call AssignPuff(Puffs(Child1), Puffs(Child2))
      Masses(:, Child2) = Masses(:, Child1)
      Extras(   Child2) = Extras(   Child1)
      Puffs(Child1)%P%X(3) = Puffs(Child1)%P%X(3) + Sqrt(12.0)*0.25*Puffs(Child1)%SigT
      Puffs(Child2)%P%X(3) = Puffs(Child2)%P%X(3) - Sqrt(12.0)*0.25*Puffs(Child2)%SigT
    End If
    Puffs(i     )%P%Active = .false.
    Puffs(Child1)%P%Active = .true.
    Puffs(Child2)%P%Active = .true.

  End If

End Subroutine Step6

!-------------------------------------------------------------------------------------------------------------

Subroutine Step7(DeltaOpt, Delta, Flow, Dt, H1, H2, Zs, WMeanAndTurb, Turbulence, MesoscaleMotions, Puff, Extra, iP)
! Time-step puffs.

  Implicit None
  ! Argument list:
  Integer,      Intent(In)    :: DeltaOpt         !
  Real(Std),    Intent(In)    :: Delta            !
  Type(Flow_),  Intent(In)    :: Flow             !
  Real(Std),    Intent(In)    :: Dt               !
  Real(Std),    Intent(In)    :: H1               ! Use for XY
  Real(Std),    Intent(In)    :: H2               ! use for xy
  Real(Std),    Intent(In)    :: Zs               ! Height over which dry deposition occurs
  Real(Std),    Intent(Out)   :: WMeanAndTurb     ! mean and turbulent vertical velocity
  Logical,      Intent(In)    :: Turbulence       !
  Logical,      Intent(In)    :: MesoscaleMotions !
  Type(Puff_),  Intent(InOut) :: Puff             !
  Type(Extra_), Intent(InOut) :: Extra            !
  Integer,      Intent(In)    :: iP               ! Puff index
  ! Locals:
  Real(Std) :: AA(3)            !
  Real(Std) :: B(3)             !
  Real(Std) :: U(3)             ! Velocity due to mean flow and mean and random plume rise.
  Real(Std) :: Beta             !
  Real(Std) :: Gamma            !
  Real(Std) :: R                !
  Integer   :: i                !
  Real(Std) :: UP               ! U VELOCITY OF THE PARTICLE/PLUME
  Real(Std) :: VP               ! V VELOCITY OF THE PARTICLE/PLUME
  Real(Std) :: WP               ! W VELOCITY OF THE PARTICLE/PLUME
  Real(Std) :: RANDPLUMX        ! RANDOM PLUME INDUCED TURBULENT COMPONENT
  Real(Std) :: RANDPLUMY        ! RANDOM PLUME INDUCED TURBULENT COMPONENT
  Real(Std) :: RANDPLUMZ        ! RANDOM PLUME INDUCED TURBULENT COMPONENT
  Real(Std) :: RANDPLUMSIG      ! Standard deviation of random displacement in
                                ! each coord direction due to plume induced
                                ! turbulence.
  Real(Std) :: IPartAge
  Real(Std) :: KFactor(3)
  Real(Std) :: KFactorM
  Real(Std) :: PuffTravelTime
  Real(Std) :: ZGradientDamping !
  Real(Std) :: ZTemp
  Real(Std) :: DtLeft

  Puff%P%XOld       = Puff%P%X
  Puff%P%iHCoordOld = Puff%P%iHCoord
  Puff%P%iZCoordOld = Puff%P%iZCoord
  Puff%XXpOld(:)    = Puff%XXp(:)

  If (Extra%FMass > 0.0) Then ! move out to Case.f90 and remove use stmt?
    IPartAge = ShortTime2RealTime(Puff%P%T)
    If (IPartAge /= 0.0) Then
      If (Extra%Inhomog) Then
        Call PlumeRise(                                                     &
               Flow%U(1), Flow%U(2), Flow%U(3), Flow%T, Flow%Theta, Flow%P, &
               Extra%UAOLD(1), Extra%UAOLD(2), Extra%UAOLD(3),              &
               Extra%TAOLD, Extra%ThetaAOLD, Extra%PAOLD,                   &
               Flow%TauUU(3), Sqrt(Flow%SigUU(3)), Flow%Eps,                &
               IPartAge, dT, .false.,                                       &
               UP, VP, WP, RANDPLUMSig,                                     &
               Extra%FMass, Extra%FM(1), Extra%FM(2), Extra%FM(3),          &
               Extra%FH, Extra%FMASS0, Extra%B0OLD                          &
             )
      Else
        Call PlumeRise(                                                     &
               Flow%U(1), Flow%U(2), Flow%U(3), Flow%T, Flow%Theta, Flow%P, &
               Extra%UAOLD(1), Extra%UAOLD(2), Extra%UAOLD(3),              &
               Extra%TAOLD, Extra%ThetaAOLD, Extra%PAOLD,                   &
               Flow%HTauUU(3), Sqrt(Flow%HSigUU(3)), Flow%HEps,             &
               IPartAge, dT, .false.,                                       &
               UP, VP, WP, RANDPLUMSig,                                     &
               Extra%FMass, Extra%FM(1), Extra%FM(2), Extra%FM(3),          &
               Extra%FH, Extra%FMASS0, Extra%B0OLD                          &
             )
    End If
    Else
      UP=0.0
      VP=0.0
      WP = Flow%U(3) + Extra%FM(3) / Extra%FMass
      RandPlumSig = 0.0
    End If
    U(1) = UP ! $$ + RandPlumSig*Gauss(iP)/Dt
    U(2) = VP ! $$ + RandPlumSig*Gauss(iP)/Dt
    U(3) = WP ! $$ + RandPlumSig*Gauss(iP)/Dt
  Else
    U(1) = Flow%U(1)
    U(2) = Flow%U(2)
    U(3) = Flow%U(3)
    RandPlumSig = 0.0
  End If

  ! Turbulence.
  If (Turbulence) Then

    If (Puff%P%X(3) * Puff%P%X(3) < 9.0 * Puff%XXp(3)) Then ! Note 9 hardwired $$
      ZGradientDamping = Erf(Puff%P%X(3) / Sqrt(2.0 * Puff%XXp(3)))
    Else
      ZGradientDamping = 1.0
    End If

    ! Velocity memory.
    If (Extra%VelMem) Then

      If (.not. Puff%R) Then ! $$ shouldn't be needed except through rounding errors - ensure isn't needed
        Beta = 0.0
      Else If (Puff%XXp(3) >= Delta**2) Then
        Beta = 1.0
      Else If (DeltaOpt == 2) Then
        Beta = 1.0 - A6**2
      Else
        Beta = 0.0
      End If
      Gamma = 1.0 - Beta

      Call XUSdeTerms(Puff%P, Extra, Flow, Dt, U(3), AA, B, Zs, ZGradientDamping, Beta)
 !     Call XUSdeTerms(Puff%P, Extra, Flow, Dt, AA, B, 1.0, Beta)

      R = Gauss(iP)
      Extra%U(1) = Extra%U(1) + AA(1)
      R = Gauss(iP)
      Extra%U(2) = Extra%U(2) + AA(2)
      Extra%U(3) = Extra%U(3) + AA(3) + B(3) * Gauss(iP)
   !   Extra%U(3) = Extra%U(3) + Extra%U(3) * dt * (Flow%dSigUUdX(3) / (2.0 * Flow%SigUU(3))) * Sqrt(2.0/Pi) &
   !                * Exp(- Puff%P%X(3) * Puff%P%X(3) / (2.0 * Puff%XXp(3))) * Puff%XUp(3) / Sqrt(Puff%XXp(3))

      ! Drift term  dt [ ZWp (d tau_w / dz) / tau**2 + (gamma / 2 + WWp / 2 sigma_w^2) d sigma_w^2 / dz ].
      Extra%U(3) = Extra%U(3)                                                                           &
                   + dt * ZGradientDamping * (Puff%XUp(3) / Flow%K(3)) &
                   * (Flow%dKdX(3) / Flow%TauUU(3) - Flow%dSigUUdX(3)) &
                   + dt * ZGradientDamping * (Gamma / 2.0 + Puff%UUp(3) / (2.0 * Flow%SigUU(3))) &
                   * Flow%dSigUUdX(3)

      Puff%P%X(1) = Puff%P%X(1) + (U(1) +  Extra%U(1)) * dt / H1
      Puff%P%X(2) = Puff%P%X(2) + (U(2) +  Extra%U(2)) * dt / H2

      If (RandPlumSig > 0.0) Then
        WMeanAndTurb = U(3) +  Extra%U(3) + (Sqrt(Beta) * RandPlumSig * Gauss(iP)) / dt
      Else
        WMeanAndTurb = U(3) +  Extra%U(3)
      End If

      If (Puff%P%WSed > 0.0 ) Then
        If (Puff%P%X(3) < Zs) Then
          If (Puff%P%WSed * dt / Zs < 0.01) Then
            ZTemp = Puff%P%X(3) + (WMeanAndTurb - Puff%P%WSed * Puff%P%X(3) / Zs) * dt
          Else
            ZTemp = WMeanAndTurb * Zs / Puff%P%WSed +            &
              (Puff%P%X(3) - WMeanAndTurb * Zs / Puff%P%WSed) * Exp(-Puff%P%WSed * dt / Zs)
          End If
          If (ZTemp > Zs) Then
            If (Puff%P%WSed *(Puff%P%X(3) - Zs) / (Puff%P%X(3) * Puff%P%WSed - WMeanAndTurb * Zs) < 0.01) Then
              DtLeft = dt + Zs * (Zs - Puff%P%X(3)) / (Puff%P%X(3) * Puff%P%WSed - WMeanAndTurb * Zs)
            Else
              DtLeft = dt + Zs / Puff%P%WSed * Log(Zs * (Puff%P%WSed - WMeanAndTurb) /    &
                (Puff%P%X(3) * Puff%P%WSed - WMeanAndTurb * Zs))
            End If
            ZTemp = Zs + (WMeanAndTurb - Puff%P%WSed) * DtLeft
          End If
        Else
          ZTemp = Puff%P%X(3) + (WMeanAndTurb - Puff%P%WSed) * dt
          If (ZTemp < Zs) Then
            DtLeft = dt - (Zs - Puff%P%X(3)) / (WMeanAndTurb - Puff%P%WSed)
            If (Puff%P%WSed * DtLeft / Zs > 0.01) Then
              ZTemp = WMeanAndTurb * Zs / Puff%P%WSed +                                      &
                (Zs  - WMeanAndTurb * Zs / Puff%P%WSed) * Exp(-Puff%P%WSed * DtLeft / Zs)
            End If
          End If
        End If
      Else
        ZTemp = Puff%P%X(3) + WMeanAndTurb * dt
      End If

      Puff%P%X(3) = ZTemp

      Do i = 1, 3
        Puff%UUh(i) = Flow%SigUU(i)
        Puff%XUh(i) = Puff%XUh(i) + (Puff%UUh(i) - Puff%XUh(i) / Flow%TauUU(i)) * dt
        Puff%XXh(i) = Puff%XXh(i) + 2.0 * Puff%XUh(i) * dt + RandPlumSig**2
      End Do

      If (Puff%XUh(3) > Sqrt(Puff%UUh(3)*Puff%XXh(3))) Then
        Puff%XUh(3) = Sqrt(Puff%UUh(3)*Puff%XXh(3))
      Else If (Puff%XUh(3) < - Sqrt(Puff%UUh(3)*Puff%XXh(3))) Then
        Puff%XUh(3) = - Sqrt(Puff%UUh(3)*Puff%XXh(3))
      End If

      Do i = 1, 2
        Puff%UUp(i) = Puff%UUh(i)
        Puff%XUp(i) = Puff%XUh(i)
        Puff%XXp(i) = Puff%XXh(i)
      End Do
      Puff%UUp(3) = Puff%UUp(3) / Extra%SigUU(3)
      Puff%UUp(3) = Puff%UUp(3) + (2.0 * Gamma - 2.0 * Puff%UUp(3)) * dt / Flow%TauUU(3)
      Puff%UUp(3) = Puff%UUp(3) * Flow%SigUU(3)
      Puff%XUp(3) = Puff%XUp(3) / Sqrt(Extra%SigUU(3))
      Puff%XUp(3) = Puff%XUp(3) + (Puff%UUp(3) / Sqrt(Flow%SigUU(3)) - Puff%XUp(3) / Flow%TauUU(3)) * dt
      Puff%XUp(3) = Puff%XUp(3) * Sqrt(Flow%SigUU(3))
      Puff%XXp(3) = Puff%XXp(3) + 2.0 * Puff%XUp(3) * dt + (1.0 - Beta) * RandPlumSig**2

      ! Note: Once puff stops growing the xu correlation will tend towards 1. With
      ! numerical errors it can gradually creep up to values significantly bigger than 1.
      ! Hence the following. The similar bit above for XUh may not be necessary, but is
      ! currently included for safety.
      If (Puff%XUp(3) > Sqrt(Puff%UUp(3)*Puff%XXp(3))) Then
        Puff%XUp(3) = Sqrt(Puff%UUp(3)*Puff%XXp(3))
      Else If (Puff%XUp(3) < - Sqrt(Puff%UUp(3)*Puff%XXp(3))) Then
        Puff%XUp(3) = - Sqrt(Puff%UUp(3)*Puff%XXp(3))
      End If

      Extra%SigUU(1) = Flow%SigUU(1)
      Extra%SigUU(2) = Flow%SigUU(2)
      Extra%SigUU(3) = Flow%SigUU(3)

    ! No velocity memory.
    Else

      If (Extra%Inhomog) Then
        Call XInhomogSdeTerms(Puff%P, Flow, Dt, AA, B, ZGradientDamping)
      Else
        Call XHomogSdeTerms(Puff%P, Flow, Dt, AA, B)
      End If

      ! Near source correction. Not applied to 3rd component because can upset
      ! profiles (probably because it changes dK/dZ at constant travel time if tau
      ! varies with height). $$
      ! Note the algebraic expression is significantly faster than the exp. Change
      ! unresolved mesoscale motions too? $$
      PuffTravelTime = ShortTime2RealTime(Puff%P%T) + Dt/2.0
      If (Extra%Inhomog) Then
        KFactor(1:2) = PuffTravelTime / Sqrt(PuffTravelTime**2 + Flow%TauUU(1:2)**2)
 !       KFactor(1:2) = 1.0 - Exp(- PuffTravelTime / Flow%TauUU(1:2))
      Else
        KFactor(1:2) = PuffTravelTime / Sqrt(PuffTravelTime**2 + Flow%HTauUU(1:2)**2)
 !       KFactor(1:2) = 1.0 - Exp(- PuffTravelTime / Flow%HTauUU(1:2))
      End If
      AA(1:2) = AA(1:2) * KFactor(1:2)
      B(1:2) = B(1:2) * Sqrt(KFactor(1:2))

      If (.not. Puff%R) Then
        Beta = 0.0
      Else If (Puff%XXp(3) >= Delta**2) Then
        Beta = 1.0
      Else If (DeltaOpt == 2) Then
        Beta = 1.0 - A6**2
      Else
        Beta = 0.0
      End If
      If (Puff%R) Then ! $$ Could do similar for vel memory case? Or remove this bit?
        Beta = Max(Beta, 1.0 - (Delta**2 - Puff%XXp(3)) / B(3)**2)
        Beta = Max(0.0, Min(1.0, Beta))
      End If

      R = Gauss(iP)
      Puff%P%X(1) = Puff%P%X(1) + (U(1) * Dt + AA(1)) / H1
      R = Gauss(iP)
      Puff%P%X(2) = Puff%P%X(2) + (U(2) * Dt + AA(2)) / H2

      If (RandPlumSig > 0.0) Then ! $$ is two gauss calls in one line OK?
        WMeanAndTurb = U(3) + (AA(3) + Sqrt(Beta) * B(3) * Gauss(iP) + Sqrt(Beta) * RandPlumSig &
                     * Gauss(iP)) / dt
      Else
        WMeanAndTurb = U(3) + (AA(3) + Sqrt(Beta) * B(3) * Gauss(iP)) / dt
      End If

      If (Puff%P%WSed > 0.0 ) Then
        If (Puff%P%X(3) < Zs) Then
          If (Puff%P%WSed * dt / Zs < 0.01) Then
            ZTemp = Puff%P%X(3) + (WMeanAndTurb - Puff%P%WSed * Puff%P%X(3) / Zs) * dt
          Else
            ZTemp = WMeanAndTurb * Zs / Puff%P%WSed +            &
              (Puff%P%X(3) - WMeanAndTurb * Zs / Puff%P%WSed) * Exp(-Puff%P%WSed * dt / Zs)
          End If
          If (ZTemp > Zs) Then
            If (Puff%P%WSed *(Puff%P%X(3) - Zs) / (Puff%P%X(3) * Puff%P%WSed - WMeanAndTurb * Zs) < 0.01) Then
              DtLeft = dt + Zs * (Zs - Puff%P%X(3)) / (Puff%P%X(3) * Puff%P%WSed - WMeanAndTurb * Zs)
            Else
              DtLeft = dt + Zs / Puff%P%WSed * Log(Zs * (Puff%P%WSed - WMeanAndTurb) /    &
                (Puff%P%X(3) * Puff%P%WSed - WMeanAndTurb * Zs))
            End If
            ZTemp = Zs + (WMeanAndTurb - Puff%P%WSed) * DtLeft
          End If
        Else
          ZTemp = Puff%P%X(3) + (WMeanAndTurb - Puff%P%WSed) * dt
          If (ZTemp < Zs) Then
            DtLeft = dt - (Zs - Puff%P%X(3)) / (WMeanAndTurb - Puff%P%WSed)
            If (Puff%P%WSed * DtLeft / Zs > 0.01) Then
              ZTemp = WMeanAndTurb * Zs / Puff%P%WSed +                                      &
                (Zs  - WMeanAndTurb * Zs / Puff%P%WSed) * Exp(-Puff%P%WSed * DtLeft / Zs)
            End If
          End If
        End If
      Else
        ZTemp = Puff%P%X(3) + WMeanAndTurb * dt
      End If

      Puff%P%X(3) = ZTemp

      Do i = 1, 3
        Puff%XXh(i) = Puff%XXh(i) + B(i)**2 + RandPlumSig**2
      End Do

      Do i = 1, 2
        Puff%XXp(i) = Puff%XXh(i)
      End Do
      Puff%XXp(3) = Puff%XXp(3) + (1.0 - Beta) * B(3)**2 + (1.0 - Beta) * RandPlumSig**2

    End If

  ! No turbulence.
  Else

    ! Output from puffs for this case needs some thought (to avoid infinities) $$
    Puff%P%X(1) = Puff%P%X(1) + U(1) * Dt / H1
    Puff%P%X(2) = Puff%P%X(2) + U(2) * Dt / H2

    WMeanAndTurb = U(3)

    If (Puff%P%WSed > 0.0 ) Then
      If (Puff%P%X(3) < Zs) Then
        If (Puff%P%WSed * dt / Zs < 0.01) Then
          ZTemp = Puff%P%X(3) + (WMeanAndTurb - Puff%P%WSed * Puff%P%X(3) / Zs) * dt
        Else
          ZTemp = WMeanAndTurb * Zs / Puff%P%WSed +            &
            (Puff%P%X(3) - WMeanAndTurb * Zs / Puff%P%WSed) * Exp(-Puff%P%WSed * dt / Zs)
        End If
        If (ZTemp > Zs) Then
          If (Puff%P%WSed *(Puff%P%X(3) - Zs) / (Puff%P%X(3) * Puff%P%WSed - WMeanAndTurb * Zs) < 0.01) Then
            DtLeft = dt + Zs * (Zs - Puff%P%X(3)) / (Puff%P%X(3) * Puff%P%WSed - WMeanAndTurb * Zs)
          Else
            DtLeft = dt + Zs / Puff%P%WSed * Log(Zs * (Puff%P%WSed - WMeanAndTurb) /    &
              (Puff%P%X(3) * Puff%P%WSed - WMeanAndTurb * Zs))
          End If
          ZTemp = Zs + (WMeanAndTurb - Puff%P%WSed) * DtLeft
        End If
      Else
        ZTemp = Puff%P%X(3) + (WMeanAndTurb - Puff%P%WSed) * dt
        If (ZTemp < Zs) Then
          DtLeft = dt - (Zs - Puff%P%X(3)) / (WMeanAndTurb - Puff%P%WSed)
          If (Puff%P%WSed * DtLeft / Zs > 0.01) Then
            ZTemp = WMeanAndTurb * Zs / Puff%P%WSed +                                      &
              (Zs  - WMeanAndTurb * Zs / Puff%P%WSed) * Exp(-Puff%P%WSed * DtLeft / Zs)
          End If
        End If
      End If
    Else
      ZTemp = Puff%P%X(3) + WMeanAndTurb * dt
    End If

    Puff%P%X(3) = ZTemp

  End If

  If (MesoscaleMotions) Then
    ! Note the Dt/2 term gives the correct result for the first time step, and correct
    ! leading order behaviour for all travel times if Dt << Flow%TauUUM.
    KFactorM =                                                             &
      2.0 * Flow%SigUUM * Flow%TauUUM * Dt *                               &
      (1.0 - Exp(- (ShortTime2RealTime(Puff%P%T) + Dt/2.0) / Flow%TauUUM))
    Puff%XXh(1) = Puff%XXh(1) + KFactorM
    Puff%XXh(2) = Puff%XXh(2) + KFactorM
    Puff%XXp(1) = Puff%XXp(1) + KFactorM
    Puff%XXp(2) = Puff%XXp(2) + KFactorM
  End If

  Extra%UAOLD(1)  = Flow%U(1)
  Extra%UAOLD(2)  = Flow%U(2)
  Extra%UAOLD(3)  = Flow%U(3)
  Extra%TAOLD     = Flow%T
  Extra%ThetaAOLD = Flow%Theta
  Extra%PAOLD     = Flow%P

  Puff%P%H = Flow%H

  Puff%dXdT(:) = Puff%dXdT(:) + Flow%dUdT(:) * Dt

! Eventually use the following to turn off time spread. $$

!  ! Consider setting TInst to a finite value
!  If (Extra%TInst == InfFutureShortTime(Interval = .true.)) Then
!
!    If (                                                           &
!      Puff%XXp(1) >=                                               &
!      0.25 * (ShortTime2RealTime(Extra%TPlus - Extra%TMinus))**2 * &
!      (Flow%U(1)**2 + Flow%U(2)**2)                                &
!    ) Then
!      Extra%TInst = Puff%P%T + Extra%TPlus
!    End If
!
!  End If

End Subroutine Step7

!-------------------------------------------------------------------------------------------------------------

Subroutine PuffFraction(Puff, Extra, Time, RdT, Frac, XOffset)
! Calculates the fraction of the puff mass between Time and Time + RdT and the
! centroid of this part of the puff.

  Implicit None
  ! Argument list:
  Type(Puff_),      Intent(In)  :: Puff       !
  Type(Extra_),     Intent(In)  :: Extra      !
  Type(ShortTime_), Intent(In)  :: Time       !
  Real(Std),        Intent(In)  :: RdT        !
  Real(Std),        Intent(Out) :: Frac       !
  Real(Std),        Intent(Out) :: XOffset(3) !
  ! Locals:
  Real(Std) :: T1  ! Time relative to puff centre, limited to lie within [PTM, PTP].
  Real(Std) :: T2  ! Time + RdT relative to puff centre, limited to lie within [PTM,
                   ! PTP].
  Real(Std) :: PTM ! Puff trailing time-edge relative to puff centre.
  Real(Std) :: PTP ! Puff leading time-edge relative to puff centre.

  ! Puffs without time extent at this output time.
  If (Time >= Puff%P%T0 + Extra%TInst) Then

    Frac       = 1.0
    XOffset(:) = 0.0

  ! Puffs with time extent at this output time.
  Else

    PTM = ShortTime2RealTime(Extra%TMinus)
    PTP = ShortTime2RealTime(Extra%TPlus)
    T1  = ShortTime2RealTime(Time - (Puff%P%T0 + Puff%P%T))
    T2  = T1 + RdT
    T2  = Max(Min(T2, PTP), PTM)
    T1  = Max(Min(T1, PTP), PTM)
    If (T1 >= 0.0) Then
      If (PTP == 0.0) Then ! $$ need check on small values?
        Frac = 0.0
      Else
        Frac = (T2 - T1) * (PTP - (T1 + T2)/2.0)/PTP
      End If
    Else If (T2 <= 0.0) Then
      If (PTM == 0.0) Then
        Frac = 0.0
      Else
        Frac = (T2 - T1) * (PTM - (T1 + T2)/2.0)/PTM
      End If
    Else
      Frac = T2 * (PTP - T2/2.0)/PTP - T1 * (PTM - T1/2.0)/PTM
    End If
    Frac = Frac * 2.0 / (PTP - PTM)
    XOffset(:) = 0.5 * (T1 + T2) * Puff%dXdT

  End If

End Subroutine PuffFraction

!-------------------------------------------------------------------------------------------------------------

Subroutine PuffConc(                               &
             PuffOpts,                             &
             Puff, Puff2, Mass, Frac, XOffset, H3, &
             Coords, Grids, Reqs,                  &
             iField, iT,                           &
             Results                               &
           )
! Calculate concentration due to puffs.

  Implicit None
  ! Argument list:
  Type(PuffOpts_),   Intent(In)           :: PuffOpts
  Type(Puff_),       Intent(In)           :: Puff       !
  Type(Puff_),       Intent(In)           :: Puff2      !
  Real(Std),         Intent(In)           :: Mass(:)    ! mass or
                                                        ! dry, wet or total dep rate
  Real(Std),         Intent(In)           :: Frac       !
  Real(Std),         Intent(In)           :: XOffset(3) !
  Real(Std),         Intent(In)           :: H3
  Type(Coords_),     Intent(In)           :: Coords     !
  Type(Grids_),      Intent(In),   Target :: Grids      !
  Type(Reqs_),       Intent(In)           :: Reqs       !
  Integer,           Intent(In)           :: iField     !
  Integer,           Intent(In)           :: iT         ! $$ rename iTA (other routines below too)
  Type(Results_),    Intent(InOut)        :: Results    !
  ! Locals:
  Integer               :: iParticleSpecies !
  Integer               :: iHGrid           !
  Integer               :: iZGrid           !
  Type(HGrid_), Pointer :: HGrid
  Type(ZGrid_), Pointer :: ZGrid
  Integer               :: nX(3)            !
  Real(Std)             :: X0(3)            !
  Real(Std)             :: dX(3)            !
  Real(Std)             :: Q
  Real(Std)             :: Xp(3)
  Real(Std)             :: Xq(3)
  Real(Std)             :: XXp(3)
  Real(Std)             :: XXq(3)
  Real(Std)             :: HMax             ! Max metric coefficient.
  Real(Std)             :: H(3)             !
  Real(Std)             :: Upper(4)         ! 1,2,3,4 => x, y, direct z, reflected z
  Real(Std)             :: Lower(4)         !
  Integer               :: iUpper(4)        !
  Integer               :: iLower(4)        !
        ! Note iUpper and iLower obtained by rounding down and up respectively, and then
        ! adding 1 because the origin is at the first grid point.
  Real(Std)             :: Temp       !
  Real(Std)             :: TempX      !
  Real(Std)             :: TempY      !
  Real(Std)             :: TempZ      !
  Integer               :: i
  Integer               :: j          !
  Integer               :: k          !
  Real(Std)             :: QCopy      ! Copies used for substeps
  Real(Std)             :: XpCopy(3)  !
  Real(Std)             :: XqCopy(3)  !
  Real(Std)             :: XXpCopy(3) !
  Real(Std)             :: XXqCopy(3) !
  Integer               :: nSubStep
  Integer               :: iSubStep
  Logical               :: Gaussian
  Real(Std)             :: XL(2)
  Real(Std)             :: XT(2)
  Real(Std)             :: XLT(2)
  Real(Std)             :: XC(2)
  Real(Std)             :: A
  Real(Std)             :: BL
  Real(Std)             :: BT
  Real(Std)             :: QQ
  Real(Std)             :: SubStepFactor
  Real(Std)             :: SubStepFactorAcc ! Accumulated factor
  Real(Std), Parameter :: MaxXXRatio = 1.4

  iParticleSpecies = Reqs%FieldReqs(iField)%iParticleSpecies
  iHGrid           = Reqs%FieldReqs(iField)%iHGrid
  iZGrid           = Reqs%FieldReqs(iField)%iZGrid

  If (iHGrid /= 0) Then
    HGrid => Grids%HGrids(iHGrid)
    dX(1) = HGrid%dX
    X0(1) = HGrid%X0
    nX(1) = HGrid%nX
    dX(2) = HGrid%dY
    X0(2) = HGrid%Y0
    nX(2) = HGrid%nY
  End If
  If (iZGrid /= 0) Then
    ZGrid => Grids%ZGrids(iZGrid)
    If (ZGrid%Variable) Then
      Call Message('FATAL ERROR: variable vertical grids are not supported for puffs', 3)
    End If
    dX(3) = ZGrid%dZ
    X0(3) = ZGrid%Z0
    nX(3) = ZGrid%nZ
  End If

  If (dX(1) == 0.0) dX(1) = 0.00001 ! $$ temp fix up for fact that Floor(a/0) for a>0 seems to be -2^31 + 1
  If (dX(2) == 0.0) dX(2) = 0.00001 ! $$ which upsets calculation of iUpper and iLower.
  If (dX(3) == 0.0) dX(3) = 0.00001 !

  Call MetricCoeffs(Coords%HCoords(Puff2%P%iHCoord), Puff2%P%X(:), HMax, H(1), H(2))
  H(3) = H3

  Xp(:)  = Puff2%P%X(:)
  Xq (:) = Puff2%P%XOld(:)

  ! $$ Note this helps with cyclic grids, but wont help if a puff contributes to points on both
  ! sides of a grid join. Also question marks over unstructured grids (e.g. some use of XCycle
  ! needed but puffs may still contrib to points on both sides of join and different points in grid
  ! may not all be reduced mod XCycle or reduced mod XCycle relative to XCentre if XCentre defined).
  ! OK for most current uses of the puff scheme.
  If (iHGrid /= 0 .and. .not. HGrid%Unstructured) Then
    If (HGrid%XCycle > 0) Then
      Do While (Xp(1) < HGrid%XCentre - 0.5 * HGrid%XCycle)
        Xp(1) = Xp(1) + HGrid%XCycle
      End Do
      Do While (Xp(1) > HGrid%XCentre + 0.5 * HGrid%XCycle)
        Xp(1) = Xp(1) - HGrid%XCycle
      End Do

      Do While (Xq(1) < Xp(1) - 0.5 * HGrid%XCycle)
        Xq(1) = Xq(1) + HGrid%XCycle
      End Do
      Do While (Xq(1) > Xp(1) + 0.5 * HGrid%XCycle)
        Xq(1) = Xq(1) - HGrid%XCycle
      End Do

    End If
    If (HGrid%YCycle > 0) Then
      Do While (Xp(2) < HGrid%YCentre - 0.5 * HGrid%YCycle)
        Xp(2) = Xp(2) + HGrid%YCycle
      End Do
      Do While (Xp(2) > HGrid%YCentre + 0.5 * HGrid%YCycle)
        Xp(2) = Xp(2) - HGrid%YCycle
      End Do

      Do While (Xq(2) < Xp(2) - 0.5 * HGrid%YCycle)
        Xq(2) = Xq(2) + HGrid%YCycle
      End Do
      Do While (Xq(2) > Xp(2) + 0.5 * HGrid%YCycle)
        Xq(2) = Xq(2) - HGrid%YCycle
      End Do

    End If
  End If

  Q      = Mass(iParticleSpecies) * Frac / 2.0**Puff%N
  Xp (:) = Xp(:) + XOffset(:) / H(:)
  Xq (:) = Xq(:) + XOffset(:) / H(:)
  XXp(:) = Puff%XXp(:)
  XXq(:) = Puff%XXpOld(:)

  If (iHGrid /= 0) Then
    Gaussian = ((Xp(1) - Xq(1)) * H(1) )**2 + ((Xp(2) - Xq(2)) * H(2) )**2 < PuffOpts%A7 * XXp(1)
       ! $$ XXp(1) hard wired
    If (Gaussian) Then
      nSubStep = 0
    Else If (XXp(1) == 0.0 .or. XXq(1) == 0.0) Then
      nSubStep = 10
    Else
      nSubStep = ALog(Max(XXp(1) / XXq(1), 1.0)) / 2.0 * ALog(MaxXXRatio)
      nSubStep = Min(nSubStep, 10)
    End If
  Else
    nSubStep = 0
  End If
  ! nSubStep = 0

  SubStepLoop: Do iSubStep = 0, nSubStep

  If (nSubStep > 0) Then

    If (iSubStep == 0) Then
      QCopy      = Q
      XpCopy (:) = Xp (:)
      XqCopy (:) = Xq (:)
      XXpCopy(:) = XXp(:)
      XXqCopy(:) = XXq(:)

      Xq (:) = Xp (:)
      XXq(:) = XXp(:)
      SubStepFactorAcc = 1.0
    End If

    If (SubStepFactorAcc == 0.0) Exit SubStepLoop

    Xp (:) = Xq (:)
    XXp(:) = XXq(:)
    SubStepFactor = (XXp(1) / MaxXXRatio - XXqCopy(1)) / (XXp(1) - XXqCopy(1)) ! Factor SubStepFactorAcc
                                                                               ! will reduce by this substep
    If (SubStepFactor <= 0.0 .or. iSubStep == nSubStep) SubStepFactor = 0.0
    Xq (:) = XqCopy (:) + (Xp(:)  - XqCopy (:)) * SubStepFactor
    XXq(:) = XXqCopy(:) + (XXp(:) - XXqCopy(:)) * SubStepFactor
    Q      = QCopy * SubStepFactorAcc * (1.0 - SubStepFactor)
    SubStepFactorAcc = SubStepFactorAcc * SubStepFactor ! Fraction of step left at end of substep

  End If

  ! No grids.
  If (iHGrid == 0 .and. iZGrid == 0) Then

    If (Reqs%FieldReqs(iField)%AvBL) Then
      BL = (Xp(3) - Puff%P%H) / Sqrt(2.0 * XXp(3))
      BT = (Xp(3) + Puff%P%H) / Sqrt(2.0 * XXp(3))
      TempZ = 0.5*(Erf(BL) - Erf(BT)) / Puff%P%H
      TempZ = Abs(TempZ)
    Else
      TempZ = 1.0
    End If

    Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, 1) =           &
      Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, 1) + Q*TempZ

  ! H grid.
  Else If (iHGrid /= 0 .and. iZGrid == 0) Then

    QQ = Q / Sqrt(4.0 * Pi**2 * XXp(1) * XXp(2))

    Gaussian = ((Xp(1) - Xq(1)) * H(1) )**2 + ((Xp(2) - Xq(2)) * H(2) )**2 < &
               PuffOpts%A7 * XXp(1)              ! XXp(1) hard wired

    If (Reqs%FieldReqs(iField)%AvBL) Then
      BL = (Xp(3) - Puff%P%H) / Sqrt(2.0 * XXp(3))
      BT = (Xp(3) + Puff%P%H) / Sqrt(2.0 * XXp(3))
      TempZ = 0.5*(Erf(BL) - Erf(BT)) / Puff%P%H
      TempZ = Abs(TempZ)
    Else
      TempZ = 1.0
    End If

    ! H grid, H unstructured.
    If (HGrid%Unstructured) Then

      ! H grid, H unstructured, Gaussian.
      If (Gaussian) Then

        Do i = 1, 2
          Upper(i)  = Xp(i) + PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          Lower(i)  = Xp(i) - PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          iUpper(i) = Min(Floor((Upper(i) - X0(i))/dX(i)) + 1, nX(i))
          iLower(i) = Max(Floor((Lower(i) - X0(i))/dX(i)) + 2,     1)
        End Do

        Do i = 1, HGrid%nX
          TempX = ((HGrid%X(i) - Xp(1))**2) * H(1)**2 / (2.0 * XXp(1))
          TempY = ((HGrid%Y(i) - Xp(2))**2) * H(2)**2 / (2.0 * XXp(2))
          If (TempX + TempY > 0.5 * PuffOpts%A5 ** 2) Cycle
          Temp = Exp(- (TempX + TempY)) * TempZ
          Results%Fields(iField)%Std(1, i, 1, 1, 1, iT, 1) =           &
            Results%Fields(iField)%Std(1, i, 1, 1, 1, iT, 1) + QQ*Temp
        End Do

      ! H grid, H unstructured, non-Gaussian.
      Else

        Do i = 1, 2
          Upper(i)  = Max(Xp(i), Xq(i)) + PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          Lower(i)  = Min(Xp(i), Xq(i)) - PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          iUpper(i) = Min(Floor((Upper(i) - X0(i))/dX(i)) + 1, nX(i))
          iLower(i) = Max(Floor((Lower(i) - X0(i))/dX(i)) + 2,     1)
        End Do

        Do i = 1, HGrid%nX
          If (.not.(                                                           &
                     Upper(1) >= HGrid%X(i) .and. Lower(1) <= HGrid%X(i) .and. &
                     Upper(2) >= HGrid%Y(i) .and. Lower(2) <= HGrid%Y(i)       &
                   )                                                           &
          ) Cycle
          XL(1)  = (HGrid%X(i) - Xp(1)) * H(1)
          XT(1)  = (HGrid%X(i) - Xq(1)) * H(1)
          XLT(1) = XL(1) - XT(1)
          XL(2)  = (HGrid%Y(i) - Xp(2)) * H(2)
          XT(2)  = (HGrid%Y(i) - Xq(2)) * H(2)
          XLT(2) = XL(2) - XT(2)

          XC(:)  = (XL(1)*XLT(1) + XL(2)*XLT(2)) * XT(:) - &
                   (XT(1)*XLT(1) + XT(2)*XLT(2)) * XL(:)
          XC(:)  = XC(:) / (XLT(1)**2 + XLT(2)**2)

          A     = Sqrt(XC(1)**2 + XC(2)**2)
          TempX = A**2 / (2.0 * XXp(1))
          If (TempX > 0.5 * PuffOpts%A5 ** 2) Cycle

          BL = Sqrt((XL(1) - XC(1))**2 + (XL(2) - XC(2))**2)
          BT = Sqrt((XT(1) - XC(1))**2 + (XT(2) - XC(2))**2)
          If (Max(BL, BT) < Sqrt(XLT(1)**2 + XLT(2)**2)) BT = - BT
          BL = BL / Sqrt(2.0 * XXp(1))
          BT = BT / Sqrt(2.0 * XXp(1))

          TempY = 0.5*(Erf(BL) - Erf(BT)) / Sqrt(XLT(1)**2 + XLT(2)**2)
          TempY = Abs(TempY) * Sqrt(2.0 * Pi * XXp(1))

          Temp = Exp(- TempX) * TempY * TempZ
          Results%Fields(iField)%Std(1, i, 1, k, 1, iT, 1) =           &
            Results%Fields(iField)%Std(1, i, 1, k, 1, iT, 1) + QQ*Temp
        End Do

      End If

    ! H grid, H structured.
    Else

      ! H grid, H structured, Gaussian.
      If (Gaussian) Then

        Do i = 1, 2
          Upper(i)  = Xp(i) + PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          Lower(i)  = Xp(i) - PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          iUpper(i) = Min(Floor((Upper(i) - X0(i))/dX(i)) + 1, nX(i))
          iLower(i) = Max(Floor((Lower(i) - X0(i))/dX(i)) + 2,     1)
        End Do

        Do i = iLower(1), iUpper(1)
          TempX = ((HGrid%X(i) - Xp(1))**2) * H(1)**2 / (2.0 * XXp(1))
          Do j = iLower(2), iUpper(2)
            TempY = ((HGrid%Y(j) - Xp(2))**2) * H(2)**2 / (2.0 * XXp(2))
            If (TempX + TempY > 0.5 * PuffOpts%A5 ** 2) Cycle
            Temp = Exp(- (TempX + TempY)) * TempZ
            Results%Fields(iField)%Std(1, i, j, 1, 1, iT, 1) =           &
              Results%Fields(iField)%Std(1, i, j, 1, 1, iT, 1) + QQ*Temp
          End Do
        End Do

      ! H grid, H structured, non-Gaussian.
      Else

        Do i = 1, 2
          Upper(i)  = Max(Xp(i), Xq(i)) + PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          Lower(i)  = Min(Xp(i), Xq(i)) - PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          iUpper(i) = Min(Floor((Upper(i) - X0(i))/dX(i)) + 1, nX(i))
          iLower(i) = Max(Floor((Lower(i) - X0(i))/dX(i)) + 2,     1)
        End Do

        Do i = iLower(1), iUpper(1)
          XL(1)  = (HGrid%X(i) - Xp(1)) * H(1)
          XT(1)  = (HGrid%X(i) - Xq(1)) * H(1)
          XLT(1) = XL(1) - XT(1)
          Do j = iLower(2), iUpper(2)
            XL(2)  = (HGrid%Y(j) - Xp(2)) * H(2)
            XT(2)  = (HGrid%Y(j) - Xq(2)) * H(2)
            XLT(2) = XL(2) - XT(2)

            XC(:)  = (XL(1)*XLT(1) + XL(2)*XLT(2)) * XT(:) - &
                     (XT(1)*XLT(1) + XT(2)*XLT(2)) * XL(:)
            XC(:)  = XC(:) / (XLT(1)**2 + XLT(2)**2)

            A     = Sqrt(XC(1)**2 + XC(2)**2)
            TempX = A**2 / (2.0 * XXp(1))
            If (TempX > 0.5 * PuffOpts%A5 ** 2) Cycle

            BL = Sqrt((XL(1) - XC(1))**2 + (XL(2) - XC(2))**2)
            BT = Sqrt((XT(1) - XC(1))**2 + (XT(2) - XC(2))**2)
            If (Max(BL, BT) < Sqrt(XLT(1)**2 + XLT(2)**2)) BT = - BT
            BL = BL / Sqrt(2.0 * XXp(1))
            BT = BT / Sqrt(2.0 * XXp(1))

            TempY = 0.5*(Erf(BL) - Erf(BT)) / Sqrt(XLT(1)**2 + XLT(2)**2)
            TempY = Abs(TempY) * Sqrt(2.0 * Pi * XXp(1))

            Temp = Exp(- TempX) * TempY * TempZ
            Results%Fields(iField)%Std(1, i, j, 1, 1, iT, 1) =           &
              Results%Fields(iField)%Std(1, i, j, 1, 1, iT, 1) + QQ*Temp
          End Do
        End Do

      End If

    End If

  ! Z grid.
  Else If (iHGrid == 0 .and. iZGrid /= 0) Then

    QQ = Q / Sqrt(2.0 * Pi * XXp(3))

    Upper(3)  = Xp(3) + PuffOpts%A5 * Sqrt(XXp(3))
    Lower(3)  = Xp(3) - PuffOpts%A5 * Sqrt(XXp(3))
    iUpper(3) = Min(Int((Upper(3) - X0(3))/dX(3)) + 2, nX(3)) ! $$ Make consistent with
    iLower(3) = Max(Int((Lower(3) - X0(3))/dX(3)) + 1,     1) ! 3-d grid calc
    Upper(4)  = - Xp(3) + PuffOpts%A5 * Sqrt(XXp(3))
    Lower(4)  = - Xp(3) - PuffOpts%A5 * Sqrt(XXp(3))
    iUpper(4) = Min(Int((Upper(4) - X0(3))/dX(3)) + 2, nX(3)) ! $$ Make consistent with
    iLower(4) = Max(Int((Lower(4) - X0(3))/dX(3)) + 1,     1) ! 3-d grid calc

    Do k = iLower(3), iUpper(3)
      TempZ = ((ZGrid%Z(k) - Xp(3))**2) / (2.0*XXp(3))
      Temp  = Exp(- TempZ)
      Results%Fields(iField)%Std(1, 1, 1, k, 1, iT, 1) =           &
        Results%Fields(iField)%Std(1, 1, 1, k, 1, iT, 1) + QQ*Temp
    End Do
    Do k = iLower(4), iUpper(4)
      TempZ = ((ZGrid%Z(k) + Xp(3))**2) / (2.0*XXp(3))
      Temp  = Exp(- TempZ)
      Results%Fields(iField)%Std(1, 1, 1, k, 1, iT, 1) =           &
        Results%Fields(iField)%Std(1, 1, 1, k, 1, iT, 1) + QQ*Temp
    End Do

  ! H and Z grids.
  Else If (iHGrid /= 0 .and. iZGrid /= 0) Then

    QQ = Q / Sqrt(8.0 * Pi**3 * XXp(1) * XXp(2) * XXp(3))

    Gaussian = ((Xp(1) - Xq(1)) * H(1) )**2 + ((Xp(2) - Xq(2)) * H(2) )**2 < &
               PuffOpts%A7 * XXp(1)              ! XXp(1) hard wired

    ! H and Z grids, H unstructured.
    If (HGrid%Unstructured) Then

      ! H and Z grids, H unstructured, Gaussian.
      If (Gaussian) Then

        Do i = 1, 3
          Upper(i)  = Xp(i) + PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          Lower(i)  = Xp(i) - PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          iUpper(i) = Min(Floor((Upper(i) - X0(i))/dX(i)) + 1, nX(i))
          iLower(i) = Max(Floor((Lower(i) - X0(i))/dX(i)) + 2,     1)
        End Do
        Upper(4)  = - Xp(3) + PuffOpts%A5 * Sqrt(XXp(3))
        Lower(4)  = - Xp(3) - PuffOpts%A5 * Sqrt(XXp(3))
        iUpper(4) = Min(Floor((Upper(4) - X0(3))/dX(3)) + 1, nX(3))
        iLower(4) = Max(Floor((Lower(4) - X0(3))/dX(3)) + 2,     1)

        Do i = 1, HGrid%nX
          If (.not.(                                                           &
                     Upper(1) >= HGrid%X(i) .and. Lower(1) <= HGrid%X(i) .and. &
                     Upper(2) >= HGrid%Y(i) .and. Lower(2) <= HGrid%Y(i)       &
                   )                                                           &
          ) Cycle
          TempX = ((HGrid%X(i) - Xp(1))**2) * H(1)**2 / (2.0 * XXp(1))
          TempY = ((HGrid%Y(i) - Xp(2))**2) * H(2)**2 / (2.0 * XXp(2))
          If (TempX + TempY > 0.5 * PuffOpts%A5 ** 2) Cycle
          Do k = iLower(3), iUpper(3)
            TempZ = ((ZGrid%Z(k) - Xp(3))**2) / (2.0 * XXp(3))
            If (TempX + TempY + TempZ > 0.5 * PuffOpts%A5 ** 2) Cycle
            Temp = Exp(- (TempX + TempY + TempZ))
            Results%Fields(iField)%Std(1, i, 1, k, 1, iT, 1) =           &
              Results%Fields(iField)%Std(1, i, 1, k, 1, iT, 1) + QQ*Temp
          End Do
          Do k = iLower(4), iUpper(4)
            TempZ = ((ZGrid%Z(k) + Xp(3))**2) / (2.0 * XXp(3))
            If (TempX + TempY + TempZ > 0.5 * PuffOpts%A5 ** 2) Cycle
            Temp = Exp(- (TempX + TempY + TempZ))
            Results%Fields(iField)%Std(1, i, 1, k, 1, iT, 1) =           &
              Results%Fields(iField)%Std(1, i, 1, k, 1, iT, 1) + QQ*Temp
          End Do
        End Do

      ! H and Z grids, H unstructured, non-Gaussian.
      Else

        Do i = 1, 3
          Upper(i)  = Max(Xp(i), Xq(i)) + PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          Lower(i)  = Min(Xp(i), Xq(i)) - PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          iUpper(i) = Min(Floor((Upper(i) - X0(i))/dX(i)) + 1, nX(i))
          iLower(i) = Max(Floor((Lower(i) - X0(i))/dX(i)) + 2,     1)
        End Do
        Upper(4)  = - Xp(3) + PuffOpts%A5 * Sqrt(XXp(3))
        Lower(4)  = - Xp(3) - PuffOpts%A5 * Sqrt(XXp(3))
        iUpper(4) = Min(Floor((Upper(4) - X0(3))/dX(3)) + 1, nX(3))
        iLower(4) = Max(Floor((Lower(4) - X0(3))/dX(3)) + 2,     1)

        Do i = 1, HGrid%nX
          If (.not.(                                                           &
                     Upper(1) >= HGrid%X(i) .and. Lower(1) <= HGrid%X(i) .and. &
                     Upper(2) >= HGrid%Y(i) .and. Lower(2) <= HGrid%Y(i)       &
                   )                                                           &
          ) Cycle
          XL(1)  = (HGrid%X(i) - Xp(1)) * H(1)
          XT(1)  = (HGrid%X(i) - Xq(1)) * H(1)
          XLT(1) = XL(1) - XT(1)
          XL(2)  = (HGrid%Y(i) - Xp(2)) * H(2)
          XT(2)  = (HGrid%Y(i) - Xq(2)) * H(2)
          XLT(2) = XL(2) - XT(2)

          XC(:)  = (XL(1)*XLT(1) + XL(2)*XLT(2)) * XT(:) - &
                   (XT(1)*XLT(1) + XT(2)*XLT(2)) * XL(:)
          XC(:)  = XC(:) / (XLT(1)**2 + XLT(2)**2)

          A     = Sqrt(XC(1)**2 + XC(2)**2)
          TempX = A**2 / (2.0 * XXp(1))
          If (TempX > 0.5 * PuffOpts%A5 ** 2) Cycle

          BL = Sqrt((XL(1) - XC(1))**2 + (XL(2) - XC(2))**2)
          BT = Sqrt((XT(1) - XC(1))**2 + (XT(2) - XC(2))**2)
          If (Max(BL, BT) < Sqrt(XLT(1)**2 + XLT(2)**2)) BT = - BT
          BL = BL / Sqrt(2.0 * XXp(1))
          BT = BT / Sqrt(2.0 * XXp(1))

          TempY = 0.5*(Erf(BL) - Erf(BT)) / Sqrt(XLT(1)**2 + XLT(2)**2)
          TempY = Abs(TempY) * Sqrt(2.0 * Pi * XXp(1))

          Do k = iLower(3), iUpper(3)
            TempZ = ((ZGrid%Z(k) - Xp(3))**2) / (2.0 * XXp(3))
            If (TempX + TempZ > 0.5 * PuffOpts%A5 ** 2) Cycle
            Temp = Exp(- (TempX + TempZ)) * TempY
            Results%Fields(iField)%Std(1, i, 1, k, 1, iT, 1) =           &
              Results%Fields(iField)%Std(1, i, 1, k, 1, iT, 1) + QQ*Temp
          End Do
          Do k = iLower(4), iUpper(4)
            TempZ = ((ZGrid%Z(k) + Xp(3))**2) / (2.0 * XXp(3))
            If (TempX + TempZ > 0.5 * PuffOpts%A5 ** 2) Cycle
            Temp = Exp(- (TempX + TempZ)) * TempY
            Results%Fields(iField)%Std(1, i, 1, k, 1, iT, 1) =           &
              Results%Fields(iField)%Std(1, i, 1, k, 1, iT, 1) + QQ*Temp
          End Do
        End Do

      End If

    ! H and Z grids, H structured.
    Else

      ! H and Z grids, H structured, Gaussian.
      If (Gaussian) Then

        Do i = 1, 3
          Upper(i)  = Xp(i) + PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          Lower(i)  = Xp(i) - PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          iUpper(i) = Min(Floor((Upper(i) - X0(i))/dX(i)) + 1, nX(i))
          iLower(i) = Max(Floor((Lower(i) - X0(i))/dX(i)) + 2,     1)
        End Do
        Upper(4)  = - Xp(3) + PuffOpts%A5 * Sqrt(XXp(3))
        Lower(4)  = - Xp(3) - PuffOpts%A5 * Sqrt(XXp(3))
        iUpper(4) = Min(Floor((Upper(4) - X0(3))/dX(3)) + 1, nX(3))
        iLower(4) = Max(Floor((Lower(4) - X0(3))/dX(3)) + 2,     1)

        Do i = iLower(1), iUpper(1)
          TempX = ((HGrid%X(i) - Xp(1))**2) * H(1)**2 / (2.0 * XXp(1))
          Do j = iLower(2), iUpper(2)
            TempY = ((HGrid%Y(j) - Xp(2))**2) * H(2)**2 / (2.0 * XXp(2))
            If (TempX + TempY > 0.5 * PuffOpts%A5 ** 2) Cycle
            Do k = iLower(3), iUpper(3)
              TempZ = ((ZGrid%Z(k) - Xp(3))**2) / (2.0 * XXp(3))
              If (TempX + TempY + TempZ > 0.5 * PuffOpts%A5 ** 2) Cycle
              Temp = Exp(- (TempX + TempY + TempZ))
              Results%Fields(iField)%Std(1, i, j, k, 1, iT, 1) =           &
                Results%Fields(iField)%Std(1, i, j, k, 1, iT, 1) + QQ*Temp
            End Do
            Do k = iLower(4), iUpper(4)
              TempZ = ((ZGrid%Z(k) + Xp(3))**2) / (2.0 * XXp(3))
              If (TempX + TempY + TempZ > 0.5 * PuffOpts%A5 ** 2) Cycle
              Temp = Exp(- (TempX + TempY + TempZ))
              Results%Fields(iField)%Std(1, i, j, k, 1, iT, 1) =           &
                Results%Fields(iField)%Std(1, i, j, k, 1, iT, 1) + QQ*Temp
            End Do
          End Do
        End Do

      ! H and Z grids, H structured, non-Gaussian.
      Else

        Do i = 1, 3
          Upper(i)  = Max(Xp(i), Xq(i)) + PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          Lower(i)  = Min(Xp(i), Xq(i)) - PuffOpts%A5 * Sqrt(XXp(i)) / H(i)
          iUpper(i) = Min(Floor((Upper(i) - X0(i))/dX(i)) + 1, nX(i))
          iLower(i) = Max(Floor((Lower(i) - X0(i))/dX(i)) + 2,     1)
        End Do
        Upper(4)  = - Xp(3) + PuffOpts%A5 * Sqrt(XXp(3))
        Lower(4)  = - Xp(3) - PuffOpts%A5 * Sqrt(XXp(3))
        iUpper(4) = Min(Floor((Upper(4) - X0(3))/dX(3)) + 1, nX(3))
        iLower(4) = Max(Floor((Lower(4) - X0(3))/dX(3)) + 2,     1)

        Do i = iLower(1), iUpper(1)
          XL(1)  = (HGrid%X(i) - Xp(1)) * H(1)
          XT(1)  = (HGrid%X(i) - Xq(1)) * H(1)
          XLT(1) = XL(1) - XT(1)
          Do j = iLower(2), iUpper(2)
            XL(2)  = (HGrid%Y(j) - Xp(2)) * H(2)
            XT(2)  = (HGrid%Y(j) - Xq(2)) * H(2)
            XLT(2) = XL(2) - XT(2)

            XC(:)  = (XL(1)*XLT(1) + XL(2)*XLT(2)) * XT(:) - &
                     (XT(1)*XLT(1) + XT(2)*XLT(2)) * XL(:)
            XC(:)  = XC(:) / (XLT(1)**2 + XLT(2)**2)

            A     = Sqrt(XC(1)**2 + XC(2)**2)
            TempX = A**2 / (2.0 * XXp(1))
            If (TempX > 0.5 * PuffOpts%A5 ** 2) Cycle

            BL = Sqrt((XL(1) - XC(1))**2 + (XL(2) - XC(2))**2)
            BT = Sqrt((XT(1) - XC(1))**2 + (XT(2) - XC(2))**2)
            If (Max(BL, BT) < Sqrt(XLT(1)**2 + XLT(2)**2)) BT = - BT
            BL = BL / Sqrt(2.0 * XXp(1))
            BT = BT / Sqrt(2.0 * XXp(1))

            TempY = 0.5*(Erf(BL) - Erf(BT)) / Sqrt(XLT(1)**2 + XLT(2)**2)
            TempY = Abs(TempY) * Sqrt(2.0 * Pi * XXp(1))

            Do k = iLower(3), iUpper(3)
              TempZ = ((ZGrid%Z(k) - Xp(3))**2) / (2.0 * XXp(3))
              If (TempX + TempZ > 0.5 * PuffOpts%A5 ** 2) Cycle
              Temp = Exp(- (TempX + TempZ)) * TempY
              Results%Fields(iField)%Std(1, i, j, k, 1, iT, 1) =           &
                Results%Fields(iField)%Std(1, i, j, k, 1, iT, 1) + QQ*Temp
            End Do
            Do k = iLower(4), iUpper(4)
              TempZ = ((ZGrid%Z(k) + Xp(3))**2) / (2.0 * XXp(3))
              If (TempX + TempZ > 0.5 * PuffOpts%A5 ** 2) Cycle
              Temp = Exp(- (TempX + TempZ)) * TempY
              Results%Fields(iField)%Std(1, i, j, k, 1, iT, 1) =           &
                Results%Fields(iField)%Std(1, i, j, k, 1, iT, 1) + QQ*Temp
            End Do
          End Do
        End Do

      End If

    End If

  End If

  EndDo SubStepLoop

End Subroutine PuffConc

!-------------------------------------------------------------------------------------------------------------

Subroutine PuffCentres(Puff, Puff2, Mass, Coords, Grids, Reqs, iField, iT, Results)
! Calculate concentration of puff centres.

  Implicit None
  ! Argument list:
  Type(Puff_),    Intent(In)    :: Puff    !
  Type(Puff_),    Intent(In)    :: Puff2   !
  Real(Std),      Intent(In)    :: Mass(:) !
  Type(Coords_),  Intent(In)    :: Coords  !
  Type(Grids_),   Intent(In)    :: Grids   !
  Type(Reqs_),    Intent(In)    :: Reqs    !
  Integer,        Intent(In)    :: iField  !
  Integer,        Intent(In)    :: iT      !
  Type(Results_), Intent(InOut) :: Results !
  ! Locals:
  Real(Std) :: ZOrigin          !
  Real(Std) :: DeltaZ           !
  Integer   :: iZ               !
  Integer   :: NZ               !
  Integer   :: iParticleSpecies !
  Integer   :: iHGrid           !
  Integer   :: iZGrid           !
  Real(Std) :: Volume           !

  iParticleSpecies = Reqs%FieldReqs(iField)%iParticleSpecies
  iHGrid           = Reqs%FieldReqs(iField)%iHGrid
  iZGrid           = Reqs%FieldReqs(iField)%iZGrid

  DeltaZ  = Grids%ZGrids(iZGrid)%dZ
  ZOrigin = Grids%ZGrids(iZGrid)%Z0
  NZ      = Grids%ZGrids(iZGrid)%nZ

  iZ = NInt((Puff2%P%X(3) - ZOrigin)/DeltaZ + 1.0) ! $$ need horizontal dependence
  If (iZ > Grids%ZGrids(iZGrid)%nZ .or. iZ < 1) Return
  Volume = Min(                                                         &
             Grids%ZGrids(iZGrid)%dZ,                                   &
             Grids%ZGrids(iZGrid)%Z(iZ) + 0.5 * Grids%ZGrids(iZGrid)%dZ &
           ) ! $$ assumes ground at z = 0.0

  Results%Fields(iField)%Std(1, 1, 1, iZ, 1, iT, 1) =   &
    Results%Fields(iField)%Std(1, 1, 1, iZ, 1, iT, 1) + &
    Mass(iParticleSpecies) / (Volume * 2.0**Puff%N)

End Subroutine PuffCentres

!-------------------------------------------------------------------------------------------------------------

Subroutine PuffMass(Puff, Puff2, Mass, Frac, Reqs, iField, iT, Results)
! Compute sigma_z^2.

  Implicit None
  ! Argument list:
  Type(Puff_),    Intent(In)    :: Puff    !
  Type(Puff_),    Intent(In)    :: Puff2   !
  Real(Std),      Intent(In)    :: Mass(:) !
  Real(Std),      Intent(In)    :: Frac    !
  Type(Reqs_),    Intent(In)    :: Reqs    !
  Integer,        Intent(In)    :: iField  !
  Integer,        Intent(In)    :: iT      !
  Type(Results_), Intent(InOut) :: Results !
  ! Locals:
  Integer   :: iParticleSpecies !

  iParticleSpecies = Reqs%FieldReqs(iField)%iParticleSpecies

  Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, 1) =    &
    Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, 1) +  &
    Frac * Mass(iParticleSpecies) / 2.0**Puff%N

End Subroutine PuffMass

!-------------------------------------------------------------------------------------------------------------

Subroutine PuffMeanZ(Puff, Puff2, Mass, Frac, Reqs, iField, iT, Results)
! Compute sigma_z^2.

  Implicit None
  ! Argument list:
  Type(Puff_),    Intent(In)    :: Puff    !
  Type(Puff_),    Intent(In)    :: Puff2   !
  Real(Std),      Intent(In)    :: Mass(:) !
  Real(Std),      Intent(In)    :: Frac    !
  Type(Reqs_),    Intent(In)    :: Reqs    !
  Integer,        Intent(In)    :: iField  !
  Integer,        Intent(In)    :: iT      !
  Type(Results_), Intent(InOut) :: Results !
  ! Locals:
  Integer   :: iParticleSpecies !
  Real(Std) :: Temp             !

  iParticleSpecies = Reqs%FieldReqs(iField)%iParticleSpecies

  Temp = Puff2%P%X(3) / Sqrt(2.0*Puff%XXp(3))
  Temp = Puff2%P%X(3)*erf(Temp) + Sqrt(2.0*Puff%XXp(3)/Pi)*Exp(-Temp**2)
  Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, 1) =   &
    Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, 1) + &
    Temp * Frac * Mass(iParticleSpecies) / 2**Puff%N

End Subroutine PuffMeanZ

!-------------------------------------------------------------------------------------------------------------

Subroutine PuffSigZ2(Puff, Puff2, Mass, Frac, Reqs, iField, iT, Results)
! Compute sigma_z^2.

  Implicit None
  ! Argument list:
  Type(Puff_),    Intent(In)    :: Puff    !
  Type(Puff_),    Intent(In)    :: Puff2   !
  Real(Std),      Intent(In)    :: Mass(:) !
  Real(Std),      Intent(In)    :: Frac    !
  Type(Reqs_),    Intent(In)    :: Reqs    !
  Integer,        Intent(In)    :: iField  !
  Integer,        Intent(In)    :: iT      !
  Type(Results_), Intent(InOut) :: Results !
  ! Locals:
  Integer   :: iParticleSpecies !

  iParticleSpecies = Reqs%FieldReqs(iField)%iParticleSpecies

  Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, 1) =                            &
    Results%Fields(iField)%Std(1, 1, 1, 1, 1, iT, 1) +                          &
    (Puff%XXp(3) + Puff2%P%X(3)**2) * Frac * Mass(iParticleSpecies) / 2**Puff%N

End Subroutine PuffSigZ2

!-------------------------------------------------------------------------------------------------------------

Subroutine PuffInfo(                                         &
             iCase,                                          &
             Puff, Puff2, Extra, Mass,                       &
             Flow, Cloud, Rain,                              &
             OutputOpts,                                     &
             Coords, Grids, Flows, Specieses, Sources, Reqs, &
             iPPInfo, iT,                                    &
             Units, Results                                  &
           )
! Outputs puff info and infomation on flow properties etc which have been used to
! time-step the particle.

  Implicit None
  ! Argument list:
  Integer,           Intent(In)    :: iCase      ! Number of case.
  Type(Puff_),       Intent(In)    :: Puff       !
  Type(Puff_),       Intent(In)    :: Puff2      !
  Type(Extra_),      Intent(In)    :: Extra      !
  Real(Std),         Intent(In)    :: Mass(:)    !
  Type(Flow_),       Intent(In)    :: Flow       !
  Type(Cloud_),      Intent(In)    :: Cloud      !
  Type(Rain_),       Intent(In)    :: Rain       !
  Type(OutputOpts_), Intent(In)    :: OutputOpts ! Output options.
  Type(Coords_),     Intent(In)    :: Coords
  Type(Grids_),      Intent(In)    :: Grids      !
  Type(Flows_),      Intent(In)    :: Flows
  Type(Specieses_),  Intent(In)    :: Specieses  !
  Type(Sources_),    Intent(In)    :: Sources    !
  Type(Reqs_),       Intent(In)    :: Reqs       !
  Integer,           Intent(In)    :: iPPInfo    !
  Integer,           Intent(In)    :: iT         !
  Type(Units_),      Intent(InOut) :: Units      ! Collection of information on input/output unit numbers.
  Type(Results_),    Intent(InOut) :: Results    !
  ! Locals:
  Integer                        :: iTGrid
  Type(Time_)                    :: T
  Type(Time_)                    :: T0
  Character(MaxOutputLineLength) :: Line
  Character(8)                   :: Skew
  Character(8)                   :: VelMem
  Character(8)                   :: Inhomog
  Integer                        :: i
  Real(Std) :: RH             ! Relative humidity (%).
  Real(Std) :: CloudOktas     ! Cloud Cover in Oktas.
  Real(Std) :: WindSpeed      ! Wind Speed.
  Real(Std) :: WindDirection  ! Wind Direction (degrees).

  ! Index of TGrid.
  iTGrid = Reqs%PPInfoReqs(iPPInfo)%iTGrid

  ! Numerical screen and disk output.
  If (Reqs%PPInfoReqs(iPPInfo)%Screen .or. &
      Reqs%PPInfoReqs(iPPInfo)%Disk) Then

    If (iTGrid == 0) Then
      T = ShortTime2Time(Puff%P%T0 + Puff%P%T)
    Else
      T = ShortTime2Time(TInTGrid(Grids%TGrids(iTGrid), iT))
    End If
    T0 = ShortTime2Time(Puff%P%T0)

    ! Construct output line.

    ! 1. Basic information.
    Line =                                                                         &
        Trim(FormatChar('Y',                                  15, .true., 'R')) // &
        Trim(  Int2Char(Puff%P%iUP,                           15, .true., 'R')) // &
        Trim(FormatChar(Sources%Sources(Puff%P%iSource)%Name, 15, .true., 'R')) // &
        Trim(FormatChar(Time2Char(T0, .true., 0, .false.),    25, .true., 'R')) // &
        Trim(FormatChar(Time2Char(T,  .true., 0, .false.),    25, .true., 'R')) // &
        Trim(  Std2Char(ShortTime2RealTime(Puff%P%T),         15, .true., 'R')) // &
        Trim(  Std2Char(Puff2%P%X(1),                         15, .true., 'R')) // &
        Trim(  Std2Char(Puff2%P%X(2),                         15, .true., 'R')) // &
        Trim(  Std2Char(Puff2%P%X(3),                         15, .true., 'R')) // &
        Trim(  Std2Char(Extra%U(1),                           15, .true., 'R')) // &
        Trim(  Std2Char(Extra%U(2),                           15, .true., 'R')) // &
        Trim(  Std2Char(Extra%U(3),                           15, .true., 'R')) // &
        Trim(  Std2Char(Puff%XXp(1),                          15, .true., 'R')) // &
        Trim(  Std2Char(Puff%XXp(2),                          15, .true., 'R')) // &
        Trim(  Std2Char(Puff%XXp(3),                          15, .true., 'R'))

    ! 2. Met information. $$ Currently from puff (i.e. previous time step) - use Flow instead?
    If (Reqs%PPInfoReqs(iPPInfo)%Met) Then

      ! Limit relative humidity.
      RH = CalcRH(Flow%Q, Flow%T, Flow%P)
      If (RH <   0.0) RH =   0.0
      If (RH > 100.0) RH = 100.0

      CloudOktas = Cloud%Cloud * 8.0

      ! Calculate wind speed and direction at particle location
      WindSpeed = Sqrt(Extra%UAOLD(1)**2 + Extra%UAOLD(2)**2)

      WindDirection = ATan2ZeroTest(-Extra%UAOLD(1), -Extra%UAOLD(2)) * 180.0 / Pi
      If (WindDirection < 0.0) WindDirection = WindDirection + 360.0

      Line = Trim(Line)                                       // & ! $$ check old values
             Trim(Std2Char(Extra%UAOLD(1),  15, .true., 'R')) // & ! always defined.
             Trim(Std2Char(Extra%UAOLD(2),  15, .true., 'R')) // &
             Trim(Std2Char(Extra%UAOLD(3),  15, .true., 'R')) // &
             Trim(Std2Char(Extra%SigUU(1),  15, .true., 'R')) // &
             Trim(Std2Char(Extra%SigUU(2),  15, .true., 'R')) // &
             Trim(Std2Char(Extra%SigUU(3),  15, .true., 'R')) // &
             Trim(Std2Char(Extra%TAOLD,     15, .true., 'R')) // &
             Trim(Std2Char(Extra%PAOLD,     15, .true., 'R')) // &
             Trim(Std2Char(Extra%ThetaAOLD, 15, .true., 'R')) // &
             Trim(Std2Char(Puff%P%H,        15, .true., 'R')) // &
             Trim(Std2Char(CloudOktas,      15, .true., 'R')) // &
             Trim(Std2Char(RH,              15, .true., 'R')) // &
             Trim(Std2Char(WindSpeed,       15, .true., 'R')) // &
             Trim(Std2Char(WindDirection,   15, .true., 'R'))

    End If

    ! 3. Mass information
    If (Reqs%PPInfoReqs(iPPInfo)%Mass) Then
      Do i = 1, Specieses%nParticleSpecieses
        Line = Trim(Line) // Std2Char(Mass(i), 15, .true., 'R')
      End Do
    End If

    ! 4. Plume rise information. $$ Need to check what these values are actually?
    If (Reqs%PPInfoReqs(iPPInfo)%PlumeRise) Then
      Line = Trim(Line)                                      // &
             Trim(Std2Char(Extra%FMass,    15, .true., 'R')) // &
             Trim(Std2Char(Extra%FMass0,   15, .true., 'R')) // &
             Trim(Std2Char(Extra%FM(1),    15, .true., 'R')) // &
             Trim(Std2Char(Extra%FM(2),    15, .true., 'R')) // &
             Trim(Std2Char(Extra%FM(3),    15, .true., 'R')) // &
             Trim(Std2Char(Extra%FH,       15, .true., 'R'))
    End If

    ! 5. Dispersion scheme information.
    If (Reqs%PPInfoReqs(iPPInfo)%DispersionScheme) Then
      Write (Skew,    *) Extra%Skew
      Write (VelMem,  *) Extra%VelMem
      Write (Inhomog, *) Extra%Inhomog
      Line = Trim(Line)                                 // &
             Trim(FormatChar(Skew,    15, .true., 'R')) // &
             Trim(FormatChar(VelMem,  15, .true., 'R')) // &
             Trim(FormatChar(Inhomog, 15, .true., 'R'))
    End If

    ! 6. Puff family information.
    If (Reqs%PPInfoReqs(iPPInfo)%PuffFamily) Then
      Line = Trim(Line)                                          // &
             Trim(Int2Char(Puff%OriginalPuff,  15, .true., 'R')) // &
             Trim(Int2Char(Puff%N,             15, .true., 'R')) // &
             Trim(Int2Char(Puff%Parent,        15, .true., 'R')) // &
             Trim(Int2Char(Puff%Sibling,       15, .true., 'R'))
    End If

    Call OutputPPInfos(                                    &
           Line,                                           &
           Puff2%P%X, Puff%P%XOld, Puff%P%T, Puff%P%TOld,  & ! $$ coords for XOld may be suspect
           iCase, iPPInfo, iT, Puff%P%iUP, .true.,         &
           OutputOpts,                                     &
           Coords, Grids, Flows, Specieses, Sources, Reqs, &
           Units, Results                                  &
         )

  End If

End Subroutine PuffInfo

!-------------------------------------------------------------------------------------------------------------

End Module PuffModule



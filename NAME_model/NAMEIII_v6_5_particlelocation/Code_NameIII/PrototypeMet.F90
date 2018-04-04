! Module: Prototype Met Module

Module PrototypeMetModule

! This is a simple prototype met module which reads some basic met data from a file.
! The met is constant in time and in the horizontal and doesn't vary between
! cases. The topographic height is assumed to be zero.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use CommonMetModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: PrototypeMet_
Public  :: InitPrototypeMet
Public  :: PrepareForUpdatePrototypeMet ! Prepares for updating an instance of the prototype met module.
Public  :: UpdatePrototypeMet
Public  :: Levels_ ! (exported for prototypeflow)
Public  :: Met_    ! (exported for prototypeflow)

!-------------------------------------------------------------------------------------------------------------

Type :: Levels_ ! Values of a met parameter at three particular heights.
  Real(Std) Z     ! Value at z = 0.
  Real(Std) H     ! Value at the boundary layer top.
  Real(Std) HPlus ! Value at the top of the transition between boundary layer and free
                  ! atmosphere.
End Type Levels_

!-------------------------------------------------------------------------------------------------------------

Type :: Met_ ! Met data.
  Real(Std)     H      ! Boundary layer depth.
  Real(Std)     HPlus  ! Top of transition between boundary layer and free atmosphere.
  Type(Levels_) SigW2  ! SigW2.
  Type(Levels_) TauW   ! TauW.
  Type(Levels_) KZ     ! KZ.
  Type(Levels_) SigU2  ! SigU2.
  Type(Levels_) TauU   ! TauU.
  Type(Levels_) KX     ! KX.
  Type(Levels_) SigV2  ! SigV2.
  Type(Levels_) TauV   ! TauV.
  Type(Levels_) KY     ! KY.
  Type(Levels_) U      ! U.
  Type(Levels_) V      ! V.
  Type(Levels_) DeltaI ! Limit on Delta due to inhomogeneity.
  Real(Std)     Z00    ! Height below which turbulence properties are constant.
End Type Met_

!-------------------------------------------------------------------------------------------------------------

Type :: PrototypeMet_ ! Information describing the state of an instance of the
                      ! prototype met module.
  Type(CommonMet_)             C       ! The part of the met state common to all met modules.
  Character(MaxFileNameLength) MetFile ! Name of file containing met data.
  Type(Met_)                   Met     ! Met data.
End Type PrototypeMet_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitPrototypeMet(         &
           MetName,                &
           HCoordName, ZCoordName, &
           FixedMet,               &
           UpdateOnDemand,         &
           MetFile                 &
         )
! Initialises an instance of the prototype met module.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: MetName        ! Name of met module instance.
  Character(*), Intent(In) :: HCoordName
  Character(*), Intent(In) :: ZCoordName
  Logical,      Intent(In) :: FixedMet
  Logical,      Intent(In) :: UpdateOnDemand ! Indicates the met module instance is to be updated using
                                             ! update-on-demand.
  Character(*), Intent(In) :: MetFile
  ! Function result:
  Type(PrototypeMet_) InitPrototypeMet !
  ! Locals:
  Type(PrototypeMet_) PrototypeMet !
  Character(MaxCharLength) HCoordNames(MaxHCoords)
  Character(MaxCharLength) ZCoordNames(MaxZCoords)

  If (Len_Trim(MetFile) > MaxFileNameLength) Then
    Call Message('Error in InitPrototypeMet', 3)
  End If

  If (Len_Trim(HCoordName) > MaxCharLength .or. Len_Trim(HCoordName) == 0) Then
    Call Message('Error', 3)
  End If
  If (Len_Trim(ZCoordName) > MaxCharLength .or. Len_Trim(ZCoordName) == 0) Then
    Call Message('Error', 3)
  End If
  HCoordNames(1) = HCoordName
  ZCoordNames(1) = ZCoordName

  PrototypeMet%C = InitCommonMet('Prototype Met', MetName, FixedMet, UpdateOnDemand)
  Call AddCoordsToCommonMet(             &
         1, HCoordNames, 1, ZCoordNames, &
         PrototypeMet%C                  &
       )
  PrototypeMet%MetFile = MetFile

  InitPrototypeMet = PrototypeMet

End Function InitPrototypeMet

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForUpdatePrototypeMet( &
             Coords, Grids,              &
             iCase, iMetCase,            &
             Time,                       &
             TValid, UpdateNow,          &
             PrototypeMet,               &
             Units                       &
           )
! Prepares for updating an instance of the prototype met module.

! This routine must set TValid and UpdateNow but must not alter the validity of the prototype met module
! instance.
! $$ this comment should be with generic definition.

  Implicit None
  ! Argument List:
  Type(Coords_),       Intent(In)           :: Coords
  Type(Grids_),        Intent(In),   Target :: Grids
  Integer,             Intent(In)           :: iCase
  Integer,             Intent(In)           :: iMetCase
  Type(Time_),         Intent(In)           :: Time
  Type(Time_),         Intent(Out)          :: TValid
  Logical,             Intent(Out)          :: UpdateNow
  Type(PrototypeMet_), Intent(InOut)        :: PrototypeMet
  Type(Units_),        Intent(InOut)        :: Units
  ! Coords       :: Collection of coord systems.
  ! Grids        :: Collection of grids.
  ! iCase        :: Number of case.
  ! iMetCase     :: Number of the met realisation in the met ensemble.
  ! Time         :: Time for which the prototype met module instance might be updated.
  ! TValid       :: Earliest time that the validity (overall or for any single attribute) of the met module
  !                 instance might change, assuming the met module instance is updated now. The value is that
  !                 determined at the end of this routine (the actual time may be later).
  ! UpdateNow    :: Indicates the met module instance must be updated now (even if update-on-demand is
  !                 specified). If set, TValid need not be set to any particular time.
  ! PrototypeMet :: State of a prototype met module instance.
  ! Units        :: Collection of information on input/output unit numbers.

  TValid = InfFutureTime()
  UpdateNow = .false.

End Subroutine PrepareForUpdatePrototypeMet

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdatePrototypeMet( &
             Coords, Grids,    &
             iCase, iMetCase,  &
             Time,             &
             PrototypeMet,     &
             Units             &
           )
! Updates an instance of the prototype met module.

  Implicit None
  ! Argument list:
  Type(Coords_),       Intent(In)    :: Coords       ! Collection of coord systems.
  Type(Grids_),        Intent(In)    :: Grids        ! Collection of grids.
  Integer,             Intent(In)    :: iCase        ! Number of case.
  Integer,             Intent(In)    :: iMetCase     ! Number of the met realisation in the met ensemble.
  Type(Time_),         Intent(In)    :: Time         ! Time for which the prototype met module instance is to
                                                     ! be updated.
  Type(PrototypeMet_), Intent(InOut) :: PrototypeMet ! State of a prototype met module instance.
  Type(Units_),        Intent(InOut) :: Units        ! Collection of information on input/output unit numbers.
  ! Locals:
  Integer Unit   !
  Integer IOStat !

  Unit = OpenFile(                 &
    File   = PrototypeMet%MetFile, &
    Units  = Units,                &
    Status = 'Old',                &
    Action = 'Read'                &
  )

  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%H, PrototypeMet%Met%HPlus
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)
  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%SigW2
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)
  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%TauW
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)
  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%KZ
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)
  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%SigU2
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)
  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%TauU
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)
  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%KX
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)
  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%SigV2
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)
  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%TauV
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)
  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%KY
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)
  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%U
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)
  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%V
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)
  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%DeltaI
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)
  Read (Unit, *, IOStat = IOStat) PrototypeMet%Met%Z00
  If (IOStat /= 0) Call ReadError(PrototypeMet%MetFile, 'UpdatePrototypeMet', IOStat)

  Call CloseUnit(Unit, Units)

  PrototypeMet%C%Valid  = .true.
  PrototypeMet%C%TValid = InfFutureTime()

End Subroutine UpdatePrototypeMet

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadError(File, Routine, IOStat)
! Read error handler.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: File    !
  Character(*), Intent(In) :: Routine !
  Integer,      Intent(In) :: IOStat  !
  ! Local:
  Character(20) AccessL !
  Character(20) StatusL !
  Character(20) FormL   !
  Character(20) ActionL !

  !$$ need to provide interpretation of IOStat

  If (IOStat /= 0) Then
    Call Message(                     &
           'Error in '             // &
           Trim(Routine)           // &
           ': error reading file ' // &
           Trim(File),                &
           3                          &
         )
  End If

End Subroutine ReadError

!-------------------------------------------------------------------------------------------------------------

End Module PrototypeMetModule

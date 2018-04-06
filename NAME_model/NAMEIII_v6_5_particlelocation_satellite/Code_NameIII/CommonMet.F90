! Module: Common Met Module

Module CommonMetModule

! This module provides code common to all met modules.

! Met modules must not alter any elements of the CommonMet_ type after initialising it
! and setting it up using the public routines below, except for alterations to the
! elements Valid, ValidAttribs and TValid made in the module's Update????Met routine.

! Met module attributes indicate what the module can supply to the flow modules. Unlike for flow modules, the
! attributes are specific to each met module.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: CommonMet_           ! The part of the met state common to all met modules.
Public  :: InitCommonMet        ! Initialises an instance of the part of the met state
                                ! common to all met modules.
Public  :: AddCoordsToCommonMet ! Adds names of coord systems in Coords (external to
                                ! this module) which the met module instance needs to
                                ! refer to to an instance of the part of the met state
                                ! common to all met modules.
Public  :: AddGridsToCommonMet  ! Adds names of grids in Grids (external to this
                                ! module) which the met module instance needs to refer
                                ! to to an instance of the part of the met state
                                ! common to all met modules.

!-------------------------------------------------------------------------------------------------------------

Type :: CommonMet_ ! The part of the met state common to all met modules.
  Character(MaxCharLength) :: MetModName
  Character(MaxCharLength) :: MetName
  Integer                  :: iMetMod
  Integer                  :: iMet
  Integer                  :: nHCoords
  Character(MaxCharLength) :: HCoordNames(MaxHCoords)
  Integer                  :: iHCoords(MaxHCoords)
  Integer                  :: nZCoords
  Character(MaxCharLength) :: ZCoordNames(MaxZCoords)
  Integer                  :: iZCoords(MaxZCoords)
  Integer                  :: nHGrids
  Character(MaxCharLength) :: HGridNames(MaxHGrids)
  Integer                  :: iHGrids(MaxHGrids)
  Integer                  :: nZGrids
  Character(MaxCharLength) :: ZGridNames(MaxZGrids)
  Integer                  :: iZGrids(MaxZGrids)
  Logical                  :: FixedMet
  Logical                  :: UpdateOnDemand
  Logical                  :: Valid
  Logical                  :: ValidAttribs(MaxMetAttribs)
  Type(Time_)              :: TValid
  Logical                  :: Prefetch = .false.
  Integer                  :: iCase
  Integer                  :: iMetCase
  Type(Time_)              :: MetTime
  Logical                  :: Error
  Integer                  :: PrefetchBufferIndex
  Integer                  :: NewBufferIndex
  Logical                  :: DueForUpdate
  ! MetModName     :: Name of met module.
  ! MetName        :: Name of met module instance.
  ! iMetMod        :: Index of met module.
  ! iMet           :: Index of met module instance.
  ! nHCoords       :} Names and indices of the horizontal/vertical coord systems in
  ! HCoordNames    :} Coords (external to this module) which the met module instance
  ! iHCoords       :} needs to refer to and the number of such coord systems.
  ! nZCoords       :}
  ! ZCoordNames    :}
  ! iZCoords       :}
  ! nHGrids        :] Names and indices of the horizontal/vertical grids in Grids
  ! HGridNames     :] (external to this module) which the met module instance needs to
  ! iHGrids        :] refer to and the number of such grids.
  ! nZGrids        :]
  ! ZGridNames     :]
  ! iZGrids        :]
  ! FixedMet       :: Indicates that the met is fixed.
  ! UpdateOnDemand :: Indicates the met module instance is to be updated using update-on-demand.
  ! Valid          :: Indicates the met module instance was valid immediately after the
  !                   last call to Update????Met.
  ! ValidAttribs   :: Indicates, for each of the module's attributes, the met module instance was valid
  !                   immediately after the last call to Update????Met.
  ! TValid         :: Earliest time that the validity (overall or for any single attribute) of the met module
  !                   instance might change. The value is that determined at the time Update????Met was last
  !                   called (the actual time may be later).
  ! DueForUpdate   :: Indicates PrepareForUpdateMet has been called but the update has not yet been done. This
  !                   is useful only in connection with update-on-demand.

  ! MetName must identify the met module instance uniquely, without the need to specify MetModName.

  ! The times for which a met module instance is valid for an attribute are determined as follows. We
  ! consider non-fixed met first.
  !
  ! After a successful call to Update????Met at time t (this requires sufficient data to be available for the
  ! module to operate successfully for the attribute over [t, t + epsilon)), the module is valid until the end
  ! of the period it can operate successfully over. Note the module must provide a consistent set of data to
  ! the dispersion calculation for the period of its validity and needs to retain validity at the start of the
  ! period even after it is used later in the period.
  !
  ! After an unsuccessful call to Update????Met at time t, the module is invalid until the earliest time when
  ! sufficient data might be available for the module to operate successfully for the attribute over [t, t +
  ! epsilon).
  !
  ! At the end of any such period of validity or invalidity Update????Met is called again.
  !
  ! For fixed met, the module only needs to be able to operate successfully for the attribute at the fixed met
  ! time for the update to be successful, and the module validity subsequently never changes and is never
  ! updated (within a single case).

  ! $$ currently we don't have attributes for met (or we only have 1). These would be particular to met module
  ! types, not universal like for flow modules.

  ! In the above Update????Met refers to any particular met module's version of this routine.

End Type CommonMet_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitCommonMet(         &
           MetModName, MetName, &
           FixedMet,            &
           UpdateOnDemand       &
           )                    &
Result(CommonMet)
! Initialises an instance of the part of the met state common to all met modules.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: MetModName
  Character(*), Intent(In) :: MetName
  Logical,      Intent(In) :: FixedMet
  Logical,      Intent(In) :: UpdateOnDemand
  ! MetModName     :: Name of met module.
  ! MetName        :: Name of met module instance.
  ! FixedMet       :: Indicates that the met is fixed.
  ! UpdateOnDemand :: Indicates the met module instance is to be updated using update-on-demand.
  ! Function result:
  Type(CommonMet_) :: CommonMet ! Initialised instance of the part of the met state
                                ! common to all met modules.

  ! MetModName and MetName.
  If (Len_Trim(MetModName) > MaxCharLength) Then
    Call Message(                                               &
           'FATAL ERROR in InitCommonMet: met module name (' // &
           Trim(MetModName)                                  // &
           ') is too long',                                     &
           3                                                    &
         )
  End If
  If (Len_Trim(MetModName) <= 0) Then
    Call Message('FATAL ERROR in InitCommonMet: met module name is blank', 3)
  End If
  If (Len_Trim(MetName) > MaxCharLength) Then
    Call Message(                                                        &
           'FATAL ERROR in InitCommonMet: met module instance name (' // &
           Trim(MetName)                                              // &
           ') is too long',                                              &
           3                                                             &
         )
  End If
  If (Len_Trim(MetName) <= 0) Then
    Call Message('FATAL ERROR in InitCommonMet: met module instance name is blank', 3)
  End If
  CommonMet%MetModName = MetModName
  CommonMet%MetName    = MetName

  ! nHCoords and nZCoords.
  CommonMet%nHCoords = 0
  CommonMet%nZCoords = 0

  ! nHGrids and nZGrids.
  CommonMet%nHGrids = 0
  CommonMet%nZGrids = 0

  ! FixedMet.
  CommonMet%FixedMet = FixedMet

  ! UpdateOnDemand.
  CommonMet%UpdateOnDemand = UpdateOnDemand

  ! Validity.
  CommonMet%Valid           = .false.
  CommonMet%ValidAttribs(:) = .false.
  CommonMet%TValid          = InfPastTime()
  CommonMet%DueForUpdate    = .false.

End Function InitCommonMet

!-------------------------------------------------------------------------------------------------------------

Subroutine AddCoordsToCommonMet(                           &
             nHCoords, HCoordNames, nZCoords, ZCoordNames, &
             CommonMet                                     &
           )
! Adds names of coord systems in Coords (external to this module) which the met module
! instance needs to refer to to an instance of the part of the met state common to all
! met modules.

  Implicit None
  ! Argument list:
  Integer,          Intent(In)    :: nHCoords
  Character(*),     Intent(In)    :: HCoordNames(:)
  Integer,          Intent(In)    :: nZCoords
  Character(*),     Intent(In)    :: ZCoordNames(:)
  Type(CommonMet_), Intent(InOut) :: CommonMet
  ! nHCoords    :} Names of horizontal/vertical coord systems to be added and the
  ! HCoordNames :} number of such coord systems.
  ! nZCoords    :}
  ! ZCoordNames :}
  ! CommonMet   :: Instance of the part of the met state common to all met modules.
  ! Locals:
  Integer :: i ! Loop index.

  ! nHCoords and HCoordNames.
  If (nHCoords < 0) Then
    Call Message(                                                                  &
           'Error in AddCoordToCommonMet: number of horizontal coord systems (' // &
           Trim(Int2Char(nHCoords))                                             // &
           ') is given as less than zero',                                         &
           3                                                                       &
         )
  End If
  If (CommonMet%nHCoords + nHCoords > MaxHCoords) Then
    Call Message(                                                                  &
           'Error in AddCoordToCommonMet: number of horizontal coord systems (' // &
           Trim(Int2Char(CommonMet%nHCoords + nHCoords))                        // &
           ') is greater than the maximum allowed ('                            // &
           Trim(Int2Char(MaxHCoords))                                           // &
           ')',                                                                    &
           3                                                                       &
         )
  End If
  If (nHCoords > Size(HCoordNames)) Then
    Call Message(                                                                  &
           'Error in AddCoordToCommonMet: number of horizontal coord systems (' // &
           Trim(Int2Char(nHCoords))                                             // &
           ') is given as greater than the size of the array of names ('        // &
           Trim(Int2Char(Size(HCoordNames)))                                    // &
           ')',                                                                    &
           3                                                                       &
         )
  End If
  Do i = 1, nHCoords
    If (Len_Trim(HCoordNames(i)) > MaxCharLength) Then
      Call Message(                                                            &
             'Error in AddCoordToCommonMet: horizontal coord system name (' // &
             Trim(HCoordNames(i))                                           // &
             ') is too long',                                                  &
             3                                                                 &
           )
    End If
    If (Len_Trim(HCoordNames(i)) <= 0) Then
      Call Message('Error in AddCoordToCommonMet: a horizontal coord system name is blank', 3)
    End If
  End Do
  CommonMet%HCoordNames(CommonMet%nHCoords + 1:CommonMet%nHCoords + nHCoords) = &
    HCoordNames(1:nHCoords)
  CommonMet%nHCoords = CommonMet%nHCoords + nHCoords

  ! nZCoords and ZCoordNames.
  If (nZCoords < 0) Then
    Call Message(                                                                &
           'Error in AddCoordToCommonMet: number of vertical coord systems (' // &
           Trim(Int2Char(nZCoords))                                           // &
           ') is given as less than zero',                                       &
           3                                                                     &
         )
  End If
  If (CommonMet%nZCoords + nZCoords > MaxZCoords) Then
    Call Message(                                                                &
           'Error in AddCoordToCommonMet: number of vertical coord systems (' // &
           Trim(Int2Char(CommonMet%nZCoords + nZCoords))                      // &
           ') is greater than the maximum allowed ('                          // &
           Trim(Int2Char(MaxZCoords))                                         // &
           ')',                                                                  &
           3                                                                     &
         )
  End If
  If (nZCoords > Size(ZCoordNames)) Then
    Call Message(                                                                &
           'Error in AddCoordToCommonMet: number of vertical coord systems (' // &
           Trim(Int2Char(nZCoords))                                           // &
           ') is given as greater than the size of the array of names ('      // &
           Trim(Int2Char(Size(ZCoordNames)))                                  // &
           ')',                                                                  &
           3                                                                     &
         )
  End If
  Do i = 1, nZCoords
    If (Len_Trim(ZCoordNames(i)) > MaxCharLength) Then
      Call Message(                                                          &
             'Error in AddCoordToCommonMet: vertical coord system name (' // &
             Trim(ZCoordNames(i))                                         // &
             ') is too long',                                                &
             3                                                               &
           )
    End If
    If (Len_Trim(ZCoordNames(i)) <= 0) Then
      Call Message('Error in AddCoordToCommonMet: a vertical coord system name is blank', 3)
    End If
  End Do
  CommonMet%ZCoordNames(CommonMet%nZCoords + 1:CommonMet%nZCoords + nZCoords) = &
    ZCoordNames(1:nZCoords)
  CommonMet%nZCoords = CommonMet%nZCoords + nZCoords

End Subroutine AddCoordsToCommonMet

!-------------------------------------------------------------------------------------------------------------

Subroutine AddGridsToCommonMet(                        &
             nHGrids, HGridNames, nZGrids, ZGridNames, &
             CommonMet                                 &
           )
! Adds names of grids in Grids (external to this module) which the met module instance
! needs to refer to to an instance of the part of the met state common to all met
! modules.

  Implicit None
  ! Argument list:
  Integer,          Intent(In)    :: nHGrids
  Character(*),     Intent(In)    :: HGridNames(:)
  Integer,          Intent(In)    :: nZGrids
  Character(*),     Intent(In)    :: ZGridNames(:)
  Type(CommonMet_), Intent(InOut) :: CommonMet
  ! nHGrids    :} Names of horizontal/vertical grids to be added and the number of
  ! HGridNames :} such grids.
  ! nZGrids    :}
  ! ZGridNames :}
  ! CommonMet  :: Instance of the part of the met state common to all met modules.
  ! Locals:
  Integer :: i ! Loop index.

  ! nHGrids and HGridNames.
  If (nHGrids < 0) Then
    Call Message(                                                         &
           'Error in AddGridToCommonMet: number of horizontal grids (' // &
           Trim(Int2Char(nHGrids))                                     // &
           ') is given as less than zero',                                &
           3                                                              &
         )
  End If
  If (CommonMet%nHGrids + nHGrids > MaxHGrids) Then
    Call Message(                                                         &
           'Error in AddGridToCommonMet: number of horizontal grids (' // &
           Trim(Int2Char(CommonMet%nHGrids + nHGrids))                 // &
           ') is greater than the maximum allowed ('                   // &
           Trim(Int2Char(MaxHGrids))                                   // &
           ')',                                                           &
           3                                                              &
         )
  End If
  If (nHGrids > Size(HGridNames)) Then
    Call Message(                                                           &
           'Error in AddGridToCommonMet: number of horizontal grids ('   // &
           Trim(Int2Char(nHGrids))                                       // &
           ') is given as greater than the size of the array of names (' // &
           Trim(Int2Char(Size(HGridNames)))                              // &
           ')',                                                             &
           3                                                                &
         )
  End If
  Do i = 1, nHGrids
    If (Len_Trim(HGridNames(i)) > MaxCharLength) Then
      Call Message(                                                   &
             'Error in AddGridToCommonMet: horizontal grid name (' // &
             Trim(HGridNames(i))                                   // &
             ') is too long',                                         &
             3                                                        &
           )
    End If
    If (Len_Trim(HGridNames(i)) <= 0) Then
      Call Message('Error in AddGridToCommonMet: a horizontal grid name is blank', 3)
    End If
  End Do
  CommonMet%HGridNames(CommonMet%nHGrids + 1:CommonMet%nHGrids + nHGrids) = &
    HGridNames(1:nHGrids)
  CommonMet%nHGrids = CommonMet%nHGrids + nHGrids

  ! nZGrids and ZGridNames.
  If (nZGrids < 0) Then
    Call Message(                                                       &
           'Error in AddGridToCommonMet: number of vertical grids (' // &
           Trim(Int2Char(nZGrids))                                   // &
           ') is given as less than zero',                              &
           3                                                            &
         )
  End If
  If (CommonMet%nZGrids + nZGrids > MaxZGrids) Then
    Call Message(                                                       &
           'Error in AddGridToCommonMet: number of vertical grids (' // &
           Trim(Int2Char(CommonMet%nZGrids + nZGrids))               // &
           ') is greater than the maximum allowed ('                 // &
           Trim(Int2Char(MaxZGrids))                                 // &
           ')',                                                         &
           3                                                            &
         )
  End If
  If (nZGrids > Size(ZGridNames)) Then
    Call Message(                                                           &
           'Error in AddGridToCommonMet: number of vertical grids ('     // &
           Trim(Int2Char(nZGrids))                                       // &
           ') is given as greater than the size of the array of names (' // &
           Trim(Int2Char(Size(ZGridNames)))                              // &
           ')',                                                             &
           3                                                                &
         )
  End If
  Do i = 1, nZGrids
    If (Len_Trim(ZGridNames(i)) > MaxCharLength) Then
      Call Message(                                                 &
             'Error in AddGridToCommonMet: vertical grid name (' // &
             Trim(ZGridNames(i))                                 // &
             ') is too long',                                       &
             3                                                      &
           )
    End If
    If (Len_Trim(ZGridNames(i)) <= 0) Then
      Call Message('Error in AddGridToCommonMet: a vertical grid name is blank', 3)
    End If
  End Do
  CommonMet%ZGridNames(CommonMet%nZGrids + 1:CommonMet%nZGrids + nZGrids) = &
    ZGridNames(1:nZGrids)
  CommonMet%nZGrids = CommonMet%nZGrids + nZGrids

End Subroutine AddGridsToCommonMet

!-------------------------------------------------------------------------------------------------------------

End Module CommonMetModule

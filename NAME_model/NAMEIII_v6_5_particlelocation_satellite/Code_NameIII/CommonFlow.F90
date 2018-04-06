! Module: Common Flow Module

Module CommonFlowModule

! This module provides code common to all flow modules.

! Flow modules must not alter any elements of the CommonFlow_ type after initialising
! it and setting it up using the public routines below, except for alterations to the
! elements Valid, ValidAttribParams, ValidityUnimprovable and TValid made in the module's Update????Flow
! routine.

! Flow module attributes indicate what the module can do, including being updated, doing vertical coord
! conversions, and supplying information of various sorts to the flow modules. Unlike for met modules, these
! must be drawn from the attributes listed in AttribParamNames.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: AttribParamNames             ! Names of attributes (in the parameter list).
Public :: A_Update                     !} Codes for attributes.
Public :: A_Convert                    !}
Public :: A_Flow                       !}
Public :: A_Cloud                      !}
Public :: A_Rain                       !}
Public :: A_Surface                    !}
Public :: A_Soil                       !}
Public :: A_Plant                      !}
Public :: CommonFlow_                  ! The part of the flow state common to all flow modules.
Public :: InitCommonFlow               ! Initialises an instance of the part of the flow state common to all
                                       ! flow modules.
Public :: AddCoordsToCommonFlow        ! Adds names of coord systems in Coords (external to this module) which
                                       ! the flow module instance needs to refer to to an instance of the part
                                       ! of the flow state common to all flow modules.
Public :: AddGridsToCommonFlow         ! Adds names of grids in Grids (external to this module) which the flow
                                       ! module instance needs to refer to to an instance of the part of the
                                       ! flow state common to all flow modules.
Public :: AddUpdateSubsetsToCommonFlow ! Adds names of update subsets to an instance of the part of the flow
                                       ! state common to all flow modules.

!-------------------------------------------------------------------------------------------------------------

Character(MaxCharLength), Parameter :: AttribParamNames(8) = &
                                       (/                    &
                                         'Update ',          & ! 1
                                         'Convert',          & ! 2
                                         'Flow   ',          & ! 3
                                         'Cloud  ',          & ! 4
                                         'Rain   ',          & ! 5
                                         'Surface',          & ! 6
                                         'Soil   ',          & ! 7
                                         'Plant  '           & ! 8
                                       /)
! AttribParamNames :: Names of attributes (in the parameter list).

! Codes for attributes.
Integer, Parameter :: A_Update  = 1 !} Codes for attributes.
Integer, Parameter :: A_Convert = 2 !}
Integer, Parameter :: A_Flow    = 3 !}
Integer, Parameter :: A_Cloud   = 4 !}
Integer, Parameter :: A_Rain    = 5 !}
Integer, Parameter :: A_Surface = 6 !}
Integer, Parameter :: A_Soil    = 7 !}
Integer, Parameter :: A_Plant   = 8 !}

!-------------------------------------------------------------------------------------------------------------

Type :: CommonFlow_ ! The part of the flow state common to all flow modules.
  Character(MaxCharLength) :: FlowModName
  Character(MaxCharLength) :: FlowName
  Integer                  :: iFlowMod
  Integer                  :: iFlow
  Integer                  :: nMets
  Character(MaxCharLength) :: MetModNames(MaxMets)
  Character(MaxCharLength) :: MetNames(MaxMets)
  Integer                  :: iMetMods(MaxMets)
  Integer                  :: iMets(MaxMets)
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
  Character(MaxCharLength) :: DomainName
  Integer                  :: iDomain
  Integer                  :: nAttribs
  Character(MaxCharLength) :: AttribNames(MaxFlowAttribs)
  Integer                  :: iAttribs(MaxFlowAttribs)
  Integer                  :: nUpdateSubsets
  Character(MaxCharLength) :: UpdateSubsetNames(MaxFlowSubsets)
  Integer                  :: iUpdateSubsets(MaxFlowSubsets)
  Logical                  :: FixedMet
  Logical                  :: UpdateOnDemand
  Logical                  :: Valid
  Logical                  :: ValidAttribParams(Size(AttribParamNames))
  Logical                  :: ValidityUnimprovable
  Type(Time_)              :: TValid
  Logical                  :: DueForUpdate
  ! FlowModName          :: Name of flow module.
  ! FlowName             :: Name of flow module instance.
  ! iFlowMod             :: Index of flow module.
  ! iFlow                :: Index of flow module instance.
  ! nMets                :: Number of met modules used as data sources by the flow
  !                         module.
  ! MetModNames          :: Names of used met modules.
  ! MetNames             :: Names of used met module instances.
  ! iMetMods             :: Indices of used met modules.
  ! iMets                :: Indices of used met module instances.
  ! nHCoords             :} Names and indices of the horizontal/vertical coord systems in
  ! HCoordNames          :} Coords (external to this module) which the flow module
  ! iHCoords             :} instance needs to refer to and the number of such coord
  ! nZCoords             :} systems.
  ! ZCoordNames          :}
  ! iZCoords             :}
  ! nHGrids              :] Names and indices of the horizontal/vertical grids in Grids
  ! HGridNames           :] (external to this module) which the flow module instance
  ! iHGrids              :] needs to refer to and the number of such grids.
  ! nZGrids              :]
  ! ZGridNames           :]
  ! iZGrids              :]
  ! DomainName           :} Name and index of the domain of the flow module instance.
  ! iDomain              :} This domain must be unbounded in time if the met is fixed.
  ! nAttribs             :: Number of attributes of the flow module.
  ! AttribNames          :: Names and indices (in the parameters in this module) of attributes of the flow
  ! iAttribs             :: module.
  ! nUpdateSubsets       :} Names and indices of the subsets of flow module instances
  ! UpdateSubsetNames    :} used in updating the flow module instance, and the number of
  ! iUpdateSubsets       :} such subsets. These are used to specify which (other) flow
  !                         module instances are searched for information.
  ! FixedMet             :: Indicates that the met is fixed.
  ! UpdateOnDemand       :: Indicates the flow module instance is to be updated using update-on-demand.
  ! Valid                :: Indicates the flow module instance was valid for some attributes immediately after
  !                         the last call to Update????Flow.
  ! ValidAttribParams    :: Indicates, for each of the attributes in AttribParamNames other than update, the
  !                         flow module instance was valid immediately after the last call to Update????Flow.
  ! ValidityUnimprovable :: Indicates the validity of the flow module instance cannot be improved before time
  !                         TValid.
  ! TValid               :: Earliest time that the validity (overall or for any single attribute) of the flow
  !                         module instance might change, except perhaps for a change from invalid to valid
  !                         (overall or for any single attribute) caused by a change in the validity of the
  !                         met and flow module instances acting as data sources. If ValidityUnimprovable is
  !                         true, this indicates that this exception will not occur. The value is that
  !                         determined at the time Update????Flow was last called (the actual time may be
  !                         later).
  ! DueForUpdate         :: Indicates PrepareForUpdateFlow has been called but the update has not yet been
  !                         done. This is useful only in connection with update-on-demand.

  ! FlowName must identify the flow module instance uniquely, without the need to specify FlowModName.

  ! The times for which a flow module instance is valid for an attribute are determined as follows. We
  ! consider non-fixed met first.
  !
  ! After a successful call to Update????Flow at time t (this requires [t, t + epsilon) to lie within the
  ! temporal extent of the domain and sufficient of the met/flow module instances used as data sources to be
  ! valid for the module to operate successfully for the attribute over the spatial extent of its domain and
  ! over [t, t + epsilon)), the module is valid until the earliest of the end of the temporal extent of the
  ! domain, the first update of any met/flow module instances acting as data sources, and the time if any when
  ! the module needs to update itself again for reasons unconnected with changes in the data sources. Note the
  ! module must provide a consistent set of data to the dispersion calculation for the period of its validity
  ! and needs to retain validity at the start of the period even after it is used later in the period.
  !
  ! After an unsuccessful call to Update????Flow at time t, the module is invalid until the earliest of the
  ! start of the temporal extent of the domain and the first update of any met/flow module instances acting
  ! as data sources (unless this occurs at or after the end of the temporal extent of the domain).
  !
  ! At the end of any such period of validity or invalidity Update????Flow is called again.
  !
  ! For fixed met, the temporal extent of the domain is ignored (or equivalently regarded as unbounded), the
  ! module only needs to be able to operate successfully for the attribute at the fixed met time for the
  ! update to be successful, and the module validity subsequently never changes and is never updated (within a
  ! single case). The temporal extent of the domain can however be used to control the times for which the
  ! dispersion calculation can obtain information from the module.
  !
  ! $$ Currently the domain must be unbounded in time for fixed met but this should be relaxed - see
  ! description of Domain above and comments in WhichFlow. Note some flow modules may need changing to follow
  ! above ideas if we do relax this.

  ! In the above Update????Flow refers to any particular flow module's version of this routine.

End Type CommonFlow_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitCommonFlow(                 &
           FlowModName, FlowName,        &
           nMets, MetModNames, MetNames, & ! $$ remove nMets and nAttribs and use Size(array).
           DomainName,                   &
           nAttribs, AttribNames,        &
           FixedMet,                     &
           UpdateOnDemand                &
         )                               &
Result(CommonFlow)
! Initialises an instance of the part of the flow state common to all flow modules.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: FlowModName
  Character(*), Intent(In) :: FlowName
  Integer,      Intent(In) :: nMets
  Character(*), Intent(In) :: MetModNames(:)
  Character(*), Intent(In) :: MetNames(:)
  Character(*), Intent(In) :: DomainName
  Integer,      Intent(In) :: nAttribs
  Character(*), Intent(In) :: AttribNames(:)
  Logical,      Intent(In) :: FixedMet
  Logical,      Intent(In) :: UpdateOnDemand
  ! FlowModName    :: Name of flow module.
  ! FlowName       :: Name of flow module instance.
  ! nMets          :: Number of met modules used as data sources by the flow module.
  ! MetModNames    :: Names of used met modules.
  ! MetNames       :: Names of used met module instances.
  ! DomainName     :: Name of the domain of the flow module instance.
  ! nAttribs       :: Number of attributes of the flow module.
  ! AttribNames    :: Names of attributes of the flow module.
  ! FixedMet       :: Indicates that the met is fixed.
  ! UpdateOnDemand :: Indicates the flow module instance is to be updated using update-on-demand.
  ! Function result:
  Type(CommonFlow_) :: CommonFlow ! Initialised instance of the part of the flow state common to all flow
                                  ! modules.
  ! Locals:
  Integer :: i !} Loop indices.
  Integer :: j !}

  ! FlowModName and FlowName.
  If (Len_Trim(FlowModName) > MaxCharLength) Then
    Call Message(                                           &
           'Error in InitCommonFlow: flow module name (' // &
           Trim(FlowModName)                             // &
           ') is too long',                                 &
           3                                                &
         )
  End If
  If (Len_Trim(FlowModName) <= 0) Then
    Call Message('Error in InitCommonFlow: flow module name is blank', 3)
  End If
  If (Len_Trim(FlowName) > MaxCharLength) Then
    Call Message(                                                    &
           'Error in InitCommonFlow: flow module instance name (' // &
           Trim(FlowName)                                         // &
           ') is too long',                                          &
           3                                                         &
         )
  End If
  If (Len_Trim(FlowName) <= 0) Then
    Call Message('Error in InitCommonFlow: flow module instance name is blank', 3)
  End If
  CommonFlow%FlowModName = FlowModName
  CommonFlow%FlowName    = FlowName

  ! nMets, MetModNames and MetNames.
  If (nMets < 0) Then
    Call Message(                                                         &
           'Error in InitCommonFlow: number of met module instances (' // &
           Trim(Int2Char(nMets))                                       // &
           ') is given as less than zero',                                &
           3                                                              &
         )
  End If
  If (nMets > MaxMets) Then
    Call Message(                                                         &
           'Error in InitCommonFlow: number of met module instances (' // &
           Trim(Int2Char(nMets))                                       // &
           ') is greater than the maximum allowed ('                   // &
           Trim(Int2Char(MaxMets))                                     // &
           ')',                                                           &
           3                                                              &
         )
  End If
  If (nMets > Size(MetModNames)) Then
    Call Message(                                                           &
           'Error in InitCommonFlow: number of met modules ('            // &
           Trim(Int2Char(nMets))                                         // &
           ') is given as greater than the size of the array of names (' // &
           Trim(Int2Char(Size(MetModNames)))                             // &
           ')',                                                             &
           3                                                                &
         )
  End If
  If (nMets > Size(MetNames)) Then
    Call Message(                                                           &
           'Error in InitCommonFlow: number of met module instances ('   // &
           Trim(Int2Char(nMets))                                         // &
           ') is given as greater than the size of the array of names (' // &
           Trim(Int2Char(Size(MetNames)))                                // &
           ')',                                                             &
           3                                                                &
         )
  End If
  Do i = 1, nMets
    If (Len_Trim(MetModNames(i)) > MaxCharLength) Then
      Call Message(                                          &
             'Error in InitCommonFlow: met module name (' // &
             Trim(MetModNames(i))                         // &
             ') is too long',                                &
             3                                               &
           )
    End If
    If (Len_Trim(MetModNames(i)) <= 0) Then
      Call Message('Error in InitCommonFlow: a met module name is blank', 3)
    End If
    If (Len_Trim(MetNames(i)) > MaxCharLength) Then
      Call Message(                                                   &
             'Error in InitCommonFlow: met module instance name (' // &
             Trim(MetNames(i))                                     // &
             ') is too long',                                         &
             3                                                        &
           )
    End If
    If (Len_Trim(MetNames(i)) <= 0) Then
      Call Message('Error in InitCommonFlow: a met module instance name is blank', 3)
    End If
  End Do
  CommonFlow%nMets                = nMets
  CommonFlow%MetModNames(1:nMets) = MetModNames(1:nMets)
  CommonFlow%MetNames(1:nMets)    = MetNames(1:nMets)

  ! nHCoords and nZCoords.
  CommonFlow%nHCoords = 0
  CommonFlow%nZCoords = 0

  ! nHGrids and nZGrids.
  CommonFlow%nHGrids = 0
  CommonFlow%nZGrids = 0

  ! DomainName.
  If (Len_Trim(DomainName) > MaxCharLength) Then
    Call Message(                                      &
           'FATAL ERROR: the domain name "'         // &
           Trim(DomainName)                         // &
           '" given for the flow module instance "' // &
           Trim(FlowName)                           // &
           '" is too long.',                           &
           3                                           &
         )
  End If
  If (Len_Trim(DomainName) <= 0) Then
    Call Message(                                                                 &
           'FATAL ERROR: the domain name given for the flow module instance "' // &
           Trim(FlowName)                                                      // &
           '" is blank.',                                                         &
           3                                                                      &
         )
  End If
  CommonFlow%DomainName = DomainName

  ! nAttribs and AttribNames.
  If (nAttribs < 0) Then
    Call Message(                                               &
           'Error in InitCommonFlow: number of attributes (' // &
           Trim(Int2Char(nAttribs))                          // &
           ') is given as less than zero',                      &
           3                                                    &
         )
  End If
  If (nAttribs > MaxFlowAttribs) Then
    Call Message(                                               &
           'Error in InitCommonFlow: number of attributes (' // &
           Trim(Int2Char(nAttribs))                          // &
           ') is greater than the maximum allowed ('         // &
           Trim(Int2Char(MaxFlowAttribs))                    // &
           ')',                                                 &
           3                                                    &
         )
  End If
  If (nAttribs > Size(AttribNames)) Then
    Call Message(                                                           &
           'Error in InitCommonFlow: number of attributes ('             // &
           Trim(Int2Char(nAttribs))                                      // &
           ') is given as greater than the size of the array of names (' // &
           Trim(Int2Char(Size(AttribNames)))                             // &
           ')',                                                             &
           3                                                                &
         )
  End If
  Do i = 1, nAttribs
    If (Len_Trim(AttribNames(i)) > MaxCharLength) Then
      Call Message(                                         &
             'Error in InitCommonFlow: attribute name (' // &
             Trim(AttribNames(i))                        // &
             ') is too long',                               &
             3                                              &
           )
    End If
    If (Len_Trim(AttribNames(i)) <= 0) Then
      Call Message('Error in InitCommonFlow: an attribute name is blank', 3)
    End If
  End Do
  CommonFlow%nAttribs                = nAttribs
  CommonFlow%AttribNames(1:nAttribs) = AttribNames(1:nAttribs)

  ! iAttribs.
  Do j = 1, nAttribs
    CommonFlow%iAttribs(j) = 0
    Do i = 1, Size(AttribParamNames)
      If (AttribNames(j) .CIEq. AttribParamNames(i)) Then
        CommonFlow%iAttribs(j) = i
        Exit
      End If
    End Do
    If (CommonFlow%iAttribs(j) == 0) Then
      Call Message('UNEXPECTED FATAL ERROR in InitCommonFlow', 4)
    End If
  End Do

  ! nUpdateSubsets.
  CommonFlow%nUpdateSubsets = 0

  ! FixedMet.
  CommonFlow%FixedMet = FixedMet

  ! UpdateOnDemand.
  CommonFlow%UpdateOnDemand = UpdateOnDemand

  ! Validity.
  CommonFlow%Valid                = .false.
  CommonFlow%ValidAttribParams(:) = .false.
  CommonFlow%ValidityUnimprovable = .false.
  CommonFlow%TValid               = InfPastTime()
  CommonFlow%DueForUpdate         = .false.

End Function InitCommonFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine AddCoordsToCommonFlow(                          &
             nHCoords, HCoordNames, nZCoords, ZCoordNames, &
             CommonFlow                                    &
           )
! Adds names of coord systems in Coords (external to this module) which the flow
! module instance needs to refer to to an instance of the part of the flow state
! common to all flow modules.

  Implicit None
  ! Argument list:
  Integer,           Intent(In)    :: nHCoords
  Character(*),      Intent(In)    :: HCoordNames(:)
  Integer,           Intent(In)    :: nZCoords
  Character(*),      Intent(In)    :: ZCoordNames(:)
  Type(CommonFlow_), Intent(InOut) :: CommonFlow
  ! nHCoords    :} Names of horizontal/vertical coord systems to be added and the
  ! HCoordNames :} number of such coord systems.
  ! nZCoords    :}
  ! ZCoordNames :}
  ! CommonFlow  :: Instance of the part of the flow state common to all flow modules.
  ! Locals:
  Integer :: i ! Loop index.

  ! nHCoords and HCoordNames.
  If (nHCoords < 0) Then
    Call Message(                                                                   &
           'Error in AddCoordToCommonFlow: number of horizontal coord systems (' // &
           Trim(Int2Char(nHCoords))                                              // &
           ') is given as less than zero',                                          &
           3                                                                        &
         )
  End If
  If (CommonFlow%nHCoords + nHCoords > MaxHCoords) Then
    Call Message(                                                                   &
           'Error in AddCoordToCommonFlow: number of horizontal coord systems (' // &
           Trim(Int2Char(CommonFlow%nHCoords + nHCoords))                        // &
           ') is greater than the maximum allowed ('                             // &
           Trim(Int2Char(MaxHCoords))                                            // &
           ')',                                                                     &
           3                                                                        &
         )
  End If
  If (nHCoords > Size(HCoordNames)) Then
    Call Message(                                                                   &
           'Error in AddCoordToCommonFlow: number of horizontal coord systems (' // &
           Trim(Int2Char(nHCoords))                                              // &
           ') is given as greater than the size of the array of names ('         // &
           Trim(Int2Char(Size(HCoordNames)))                                     // &
           ')',                                                                     &
           3                                                                        &
         )
  End If
  Do i = 1, nHCoords
    If (Len_Trim(HCoordNames(i)) > MaxCharLength) Then
      Call Message(                                                             &
             'Error in AddCoordToCommonFlow: horizontal coord system name (' // &
             Trim(HCoordNames(i))                                            // &
             ') is too long',                                                   &
             3                                                                  &
           )
    End If
    If (Len_Trim(HCoordNames(i)) <= 0) Then
      Call Message('Error in AddCoordToCommonFlow: a horizontal coord system name is blank', 3)
    End If
  End Do
  CommonFlow%HCoordNames(CommonFlow%nHCoords + 1:CommonFlow%nHCoords + nHCoords) = &
    HCoordNames(1:nHCoords)
  CommonFlow%nHCoords = CommonFlow%nHCoords + nHCoords

  ! nZCoords and ZCoordNames.
  If (nZCoords < 0) Then
    Call Message(                                                                 &
           'Error in AddCoordToCommonFlow: number of vertical coord systems (' // &
           Trim(Int2Char(nZCoords))                                            // &
           ') is given as less than zero',                                        &
           3                                                                      &
         )
  End If
  If (CommonFlow%nZCoords + nZCoords > MaxZCoords) Then
    Call Message(                                                                 &
           'Error in AddCoordToCommonFlow: number of vertical coord systems (' // &
           Trim(Int2Char(CommonFlow%nZCoords + nZCoords))                      // &
           ') is greater than the maximum allowed ('                           // &
           Trim(Int2Char(MaxZCoords))                                          // &
           ')',                                                                   &
           3                                                                      &
         )
  End If
  If (nZCoords > Size(ZCoordNames)) Then
    Call Message(                                                                 &
           'Error in AddCoordToCommonFlow: number of vertical coord systems (' // &
           Trim(Int2Char(nZCoords))                                            // &
           ') is given as greater than the size of the array of names ('       // &
           Trim(Int2Char(Size(ZCoordNames)))                                   // &
           ')',                                                                   &
           3                                                                      &
         )
  End If
  Do i = 1, nZCoords
    If (Len_Trim(ZCoordNames(i)) > MaxCharLength) Then
      Call Message(                                                           &
             'Error in AddCoordToCommonFlow: vertical coord system name (' // &
             Trim(ZCoordNames(i))                                          // &
             ') is too long',                                                 &
             3                                                                &
           )
    End If
    If (Len_Trim(ZCoordNames(i)) <= 0) Then
      Call Message('Error in AddCoordToCommonFlow: a vertical coord system name is blank', 3)
    End If
  End Do
  CommonFlow%ZCoordNames(CommonFlow%nZCoords + 1:CommonFlow%nZCoords + nZCoords) = &
    ZCoordNames(1:nZCoords)
  CommonFlow%nZCoords = CommonFlow%nZCoords + nZCoords

End Subroutine AddCoordsToCommonFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine AddGridsToCommonFlow(                       &
             nHGrids, HGridNames, nZGrids, ZGridNames, &
             CommonFlow                                &
           )
! Adds names of grids in Grids (external to this module) which the flow module
! instance needs to refer to to an instance of the part of the flow state common to
! all flow modules.

  Implicit None
  ! Argument list:
  Integer,           Intent(In)    :: nHGrids
  Character(*),      Intent(In)    :: HGridNames(:)
  Integer,           Intent(In)    :: nZGrids
  Character(*),      Intent(In)    :: ZGridNames(:)
  Type(CommonFlow_), Intent(InOut) :: CommonFlow
  ! nHGrids    :} Names of horizontal/vertical grids to be added and the number of
  ! HGridNames :} such grids.
  ! nZGrids    :}
  ! ZGridNames :}
  ! CommonFlow :: Instance of the part of the flow state common to all flow modules.
  ! Locals:
  Integer :: i ! Loop index.

  ! nHGrids and HGridNames.
  If (nHGrids < 0) Then
    Call Message(                                               &
           'Error in AddGridToCommonFlow: number of grids (' // &
           Trim(Int2Char(nHGrids))                           // &
           ') is given as less than zero',                      &
           3                                                    &
         )
  End If
  If (CommonFlow%nHGrids + nHGrids > MaxHGrids) Then
    Call Message(                                               &
           'Error in AddGridToCommonFlow: number of grids (' // &
           Trim(Int2Char(CommonFlow%nHGrids + nHGrids))      // &
           ') is greater than the maximum allowed ('         // &
           Trim(Int2Char(MaxHGrids))                         // &
           ')',                                                 &
           3                                                    &
         )
  End If
  If (nHGrids > Size(HGridNames)) Then
    Call Message(                                                           &
           'Error in AddGridToCommonFlow: number of grids ('             // &
           Trim(Int2Char(nHGrids))                                       // &
           ') is given as greater than the size of the array of names (' // &
           Trim(Int2Char(Size(HGridNames)))                              // &
           ')',                                                             &
           3                                                                &
         )
  End If
  Do i = 1, nHGrids
    If (Len_Trim(HGridNames(i)) > MaxCharLength) Then
      Call Message(                                                    &
             'Error in AddGridToCommonFlow: horizontal grid name (' // &
             Trim(HGridNames(i))                                    // &
             ') is too long',                                          &
             3                                                         &
           )
    End If
    If (Len_Trim(HGridNames(i)) <= 0) Then
      Call Message('Error in AddGridToCommonFlow: a horizontal grid name is blank', 3)
    End If
  End Do
  CommonFlow%HGridNames(CommonFlow%nHGrids + 1:CommonFlow%nHGrids + nHGrids) = &
    HGridNames(1:nHGrids)
  CommonFlow%nHGrids = CommonFlow%nHGrids + nHGrids

  ! nZGrids and ZGridNames.
  If (nZGrids < 0) Then
    Call Message(                                                        &
           'Error in AddGridToCommonFlow: number of vertical grids (' // &
           Trim(Int2Char(nZGrids))                                    // &
           ') is given as less than zero',                               &
           3                                                             &
         )
  End If
  If (CommonFlow%nZGrids + nZGrids > MaxZGrids) Then
    Call Message(                                                        &
           'Error in AddGridToCommonFlow: number of vertical grids (' // &
           Trim(Int2Char(CommonFlow%nZGrids + nZGrids))               // &
           ') is greater than the maximum allowed ('                  // &
           Trim(Int2Char(MaxZGrids))                                  // &
           ')',                                                          &
           3                                                             &
         )
  End If
  If (nZGrids > Size(ZGridNames)) Then
    Call Message(                                                           &
           'Error in AddGridToCommonFlow: number of vertical grids ('    // &
           Trim(Int2Char(nZGrids))                                       // &
           ') is given as greater than the size of the array of names (' // &
           Trim(Int2Char(Size(ZGridNames)))                              // &
           ')',                                                             &
           3                                                                &
         )
  End If
  Do i = 1, nZGrids
    If (Len_Trim(ZGridNames(i)) > MaxCharLength) Then
      Call Message(                                                  &
             'Error in AddGridToCommonFlow: vertical grid name (' // &
             Trim(ZGridNames(i))                                  // &
             ') is too long',                                        &
             3                                                       &
           )
    End If
    If (Len_Trim(ZGridNames(i)) <= 0) Then
      Call Message('Error in AddGridToCommonFlow: a vertical grid name is blank', 3)
    End If
  End Do
  CommonFlow%ZGridNames(CommonFlow%nZGrids + 1:CommonFlow%nZGrids + nZGrids) = &
    ZGridNames(1:nZGrids)
  CommonFlow%nZGrids = CommonFlow%nZGrids + nZGrids

End Subroutine AddGridsToCommonFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine AddUpdateSubsetsToCommonFlow(        &
             nUpdateSubsets, UpdateSubsetNames, &
             CommonFlow                         &
           )
! Adds names of update subsets to an instance of the part of the flow state common to
! all flow modules.

  Implicit None
  ! Argument list:
  Integer,           Intent(In)    :: nUpdateSubsets
  Character(*),      Intent(In)    :: UpdateSubsetNames(:)
  Type(CommonFlow_), Intent(InOut) :: CommonFlow
  ! nUpdateSubsets    :} Names of the update subsets to be added and the number of
  ! UpdateSubsetNames :} such subsets.
  ! CommonFlow        :: Instance of the part of the flow state common to all flow
  !                      modules.
  ! Locals:
  Integer :: i ! Loop index.

  If (nUpdateSubsets < 0) Then
    Call Message(                                                                 &
           'Error in AddUpdateSubsetsToCommonFlow: number of update subsets (' // &
           Trim(Int2Char(nUpdateSubsets))                                      // &
           ') is given as less than zero',                                        &
           3                                                                      &
         )
  End If
  If (CommonFlow%nUpdateSubsets + nUpdateSubsets > MaxFlowSubsets) Then
    Call Message(                                                                 &
           'Error in AddUpdateSubsetsToCommonFlow: number of update subsets (' // &
           Trim(Int2Char(CommonFlow%nUpdateSubsets + nUpdateSubsets))          // &
           ') is greater than the maximum allowed ('                           // &
           Trim(Int2Char(MaxFlowSubsets))                                      // &
           ')',                                                                   &
           3                                                                      &
         )
  End If
  If (nUpdateSubsets > Size(UpdateSubsetNames)) Then
    Call Message(                                                                 &
           'Error in AddUpdateSubsetsToCommonFlow: number of update subsets (' // &
           Trim(Int2Char(nUpdateSubsets))                                      // &
           ') is given as greater than the size of the array of names ('       // &
           Trim(Int2Char(Size(UpdateSubsetNames)))                             // &
           ')',                                                                   &
           3                                                                      &
         )
  End If
  Do i = 1, nUpdateSubsets
    If (Len_Trim(UpdateSubsetNames(i)) > MaxCharLength) Then
      Call Message(                                                           &
             'Error in AddUpdateSubsetsToCommonFlow: update subset name (' // &
             Trim(UpdateSubsetNames(i))                                    // &
             ') is too long',                                                 &
             3                                                                &
           )
    End If
    If (Len_Trim(UpdateSubsetNames(i)) <= 0) Then
      Call Message('Error in AddUpdateSubsetsToCommonFlow: a update subset name is blank', 3)
    End If
  End Do
  CommonFlow%UpdateSubsetNames(                                              &
    CommonFlow%nUpdateSubsets + 1:CommonFlow%nUpdateSubsets + nUpdateSubsets &
  ) = UpdateSubsetNames(1:nUpdateSubsets)
  CommonFlow%nUpdateSubsets = CommonFlow%nUpdateSubsets + nUpdateSubsets

End Subroutine AddUpdateSubsetsToCommonFlow

!-------------------------------------------------------------------------------------------------------------

End Module CommonFlowModule

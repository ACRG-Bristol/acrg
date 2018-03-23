! Module: Flows Module

Module FlowsModule

! This module provides code to make all flow modules look identical.

! Module overview
! ---------------

! $$

! Certain pieces of code need to be repeated for some or all flow modules. This code repetition is achieved by
! prepreocessing Flows.P90 with Preprocessor.exe to produce Flows.F90. Whenever a new flow module is added,
! appropriate $UseBlock directives should be added for use by Preprocessor.exe.

! Note there are two ways of updating met and flow modules, update-at-once (i.e. as soon as UpdateMetsFlows is
! called) and update-on-demand (i.e. only when required).

! Module use
! ----------

! $$

! Module call tree
! ----------------

! Here xxxx ranges over the various types of flow module and yyyy over the various flow attributes.
!
! InitFlows
!
! AddxxxxFlow
!
! InitFlowOrder
!
! AddFlowOrder
!
! InitFlowSubset
!
! AddFlowSubset
!
! InitAttrib
!
! AddAttrib
!
! SetUpCoordsEtc_Flows
!
! SetUpFlows_MetsCoordsEtc
!
! SetUpFlowOrders_iFlows
!
! SetUpFlowSubsets_iFlows
!
! SetUpFlows_iFlowAttribs
!
! SetUpFlows_iFlowSubsets
!
! SetUpFlows_iMets
!
! SetUpFlows_iCoordsEtc
!
! CheckFlows_FlowUpdateSubsets
!
! FlowList
!
! MetsFlowsOverallTValid
!
! $$ above here not done yet.
!
! UpdateMetsFlows
!   |
!   |---MetsModule.UpdateMets[---(update-at-once)---MetsModule.UpdateMet]
!   |
!   |---UpdateFlows
!         |
!         |---PrepareForUpdateFlow
!         |     |
!         |     |---xxxxFlowModule.PrepareForUpdatexxxxFlow
!         |
!         |---(update-at-once)---UpdateFlow
!                                  |
!                                  |---(update-on-demand)---MetsModule.UpdateMet
!                                  |
!                                  |---xxxxFlowModule.xxxxFlowReqs
!                                  |
!                                  |---GetAttribField
!                                  |     |
!                                  |     |---WhichFlow
!                                  |     |     |
!                                  |     |     |---(update-on-demand)---UpdateFlow - -
!                                  |     |     |
!                                  |     |     |---ConvertToZKnownFlow
!                                  |     |     |     |
!                                  |     |     |     |---xxxxFlowModule.xxxxConvertCoordIndices
!                                  |     |     |     |
!                                  |     |     |     |---xxxxFlowModule.xxxxConvertToZ
!                                  |     |     |
!                                  |     |     |---ResetFlowMemorySpecificFlow
!                                  |     |           |
!                                  |     |           |---xxxxFlowModule.ResetxxxxFlowMemory
!                                  |     |           |
!                                  |     |           |---FlowPresent
!                                  |     |
!                                  |     |---GetAttribKnownFlow
!                                  |     |     |
!                                  |     |     |---xxxxFlowModule.xxxxyyyyCoordIndices
!                                  |     |     |
!                                  |     |     |---ConvertToZKnownFlow
!                                  |     |     |     |
!                                  |     |     |     |---xxxxFlowModule.xxxxConvertCoordIndices
!                                  |     |     |     |
!                                  |     |     |     |---xxxxFlowModule.xxxxConvertToZ
!                                  |     |     |
!                                  |     |     |---xxxxFlowModule.Getxxxxyyyy
!                                  |     |
!                                  |     |---ResetFlowMemory - -
!                                  |
!                                  |---xxxxFlowModule.UpdatexxxxFlow
!
! ResetFlowMemory
!   |
!   |---xxxxFlowModule.ResetxxxxFlowMemory
!
! ConvertToZ
!   |
!   |---WhichFlow
!   |     |
!   |     |---(update-on-demand)---UpdateFlow - -
!   |     |
!   |     |---ConvertToZKnownFlow
!   |     |     |
!   |     |     |---xxxxFlowModule.xxxxConvertCoordIndices
!   |     |     |
!   |     |     |---xxxxFlowModule.xxxxConvertToZ
!   |     |
!   |     |---ResetFlowMemorySpecificFlow
!   |           |
!   |           |---xxxxFlowModule.ResetxxxxFlowMemory
!   |           |
!   |           |---FlowPresent
!   |
!   |---ConvertToZKnownFlow
!         |
!         |---xxxxFlowModule.xxxxConvertCoordIndices
!         |
!         |---xxxxFlowModule.xxxxConvertToZ
!
! GetAttrib
!   |
!   |---WhichFlow
!   |     |
!   |     |---(update-on-demand)---UpdateFlow - -
!   |     |
!   |     |---ConvertToZKnownFlow
!   |     |     |
!   |     |     |---xxxxFlowModule.xxxxConvertCoordIndices
!   |     |     |
!   |     |     |---xxxxFlowModule.xxxxConvertToZ
!   |     |
!   |     |---ResetFlowMemorySpecificFlow
!   |           |
!   |           |---xxxxFlowModule.ResetxxxxFlowMemory
!   |           |
!   |           |---FlowPresent
!   |
!   |---GetAttribKnownFlow
!         |
!         |---xxxxFlowModule.xxxxyyyyCoordIndices
!         |
!         |---ConvertToZKnownFlow
!         |     |
!         |     |---xxxxFlowModule.xxxxConvertCoordIndices
!         |     |
!         |     |---xxxxFlowModule.xxxxConvertToZ
!         |
!         |---xxxxFlowModule.Getxxxxyyyy
!
! Reflect
!   |
!   |---WhichFlow
!   |     |
!   |     |---(update-on-demand)---UpdateFlow - -
!   |     |
!   |     |---ConvertToZKnownFlow
!   |     |     |
!   |     |     |---xxxxFlowModule.xxxxConvertCoordIndices
!   |     |     |
!   |     |     |---xxxxFlowModule.xxxxConvertToZ
!   |     |
!   |     |---ResetFlowMemorySpecificFlow
!   |           |
!   |           |---xxxxFlowModule.ResetxxxxFlowMemory
!   |           |
!   |           |---FlowPresent
!   |
!   |---ReflectKnownFlow
!         |
!         |---xxxxFlowModule.xxxxReflectCoordIndices
!         |
!         |---ConvertToZKnownFlow
!         |     |
!         |     |---xxxxFlowModule.xxxxConvertCoordIndices
!         |     |
!         |     |---xxxxFlowModule.xxxxConvertToZ
!         |
!         |---xxxxFlowModule.xxxxReflect
!
! CalcdZdZ
!   |
!   |---CalcdZdZ_ZBased
!   |     |
!   |     |---GetAttrib - -
!   |     |
!   |     |---ConvertToZ - -
!   |
!   |---CalcdZdZ_PBased
!   |     |
!   |     |---GetAttrib - -
!   |     |
!   |     |---ConvertToZ - -
!   |
!   |---GetAttrib - -
!
! XInTheDomain
!   |
!   |---ConvertToZ - -
!
! ResetMetsFlows
!   |
!   |---MetsModule.ResetMets
!   |
!   |---ResetFlows
!
! Calls to routines in the service modules are not included.
! Calls to routines in other modules are included (indicated by ModuleName.RoutineName), but subsequent calls
! from the called routine are not in general.
! Calls to public routines in this module are followed by "- -" indicating that one should refer to the main
! reference for this routine.
! All except one of the calls to UpdateFlow are followed by "- -", indicating that one should refer to the one
! call which isn't followed by "- -".
! Note UpdateFlow, GetAttribField and WhichFlow are potentially recursive.

!-------------------------------------------------------------------------------------------------------------

!$ Use omp_lib

Use ServiceModule
Use FlowAndFlowProfileModule
Use MetsModule
Use CommonFlowModule

Use PrototypeFlowModule
Use NWPFlowModule
Use SingleSiteFlowModule
Use BuildingFlowModule
Use RadarFlowModule
Use LINCOMFlowModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private

! Items from this module which need to be made available:
Public  :: FlowOrder_                   ! An ordered subset of flow module instances.
Public  :: FlowSubset_                  ! A subset of flow module instances.
Public  :: FlowAttrib_                  ! A flow attribute.
Public  :: Flows_                       ! A collection of flow module instance states.
Public  :: FlowMemory_                  ! Information related to a particular
                                        ! space-time point of interest, including
                                        ! information on the domains the point lies
                                        ! within, the choices of flow module instances
                                        ! approprite, and anything the individual flow
                                        ! modules instances wish to record. The
                                        ! purpose of storing this information together
                                        ! is to provide a memory to avoid the need to
                                        ! recalculate the information.
Public  :: InitFlows                    ! Initialises a collection of flow module
                                        ! instance states.

Public  :: AddPrototypeFlow                  ! Adds the state of a Prototype flow module
                                        ! instance to a collection of flow module
                                        ! instance states.
Public  :: AddSingleSiteFlow                  ! Adds the state of a SingleSite flow module
                                        ! instance to a collection of flow module
                                        ! instance states.
Public  :: AddNWPFlow                  ! Adds the state of a NWP flow module
                                        ! instance to a collection of flow module
                                        ! instance states.
Public  :: AddBuildingFlow                  ! Adds the state of a Building flow module
                                        ! instance to a collection of flow module
                                        ! instance states.
Public  :: AddRadarFlow                  ! Adds the state of a Radar flow module
                                        ! instance to a collection of flow module
                                        ! instance states.
Public  :: AddLINCOMFlow                  ! Adds the state of a LINCOM flow module
                                        ! instance to a collection of flow module
                                        ! instance states.

Public  :: InitFlowOrder                ! Initialises an ordered subset of flow module
                                        ! instances.
Public  :: AddFlowOrder                 ! Adds an ordered subset of flow module
                                        ! instances to a collection of flow module
                                        ! instance states.
Public  :: InitFlowSubset               ! Initialises a subset of flow module
                                        ! instances.
Public  :: AddFlowSubset                ! Adds a subset of flow module instances to a
                                        ! collection of flow module instance states.
Public  :: InitAttrib                   ! Initialises an attribute.
Public  :: AddAttrib                    ! Adds an attribute to the collection of flow
                                        ! module instance states.
Public  :: SetUpCoordsEtc_Flows         ! Sets up Coords and Grids by adding any extra
                                        ! coords and grids which Flows wants to
                                        ! define.
Public  :: SetUpFlows_MetsCoordsEtc     ! Sets up Flows using information from
                                        ! EtaDefns, Coords, Grids and Mets.
Public  :: SetUpFlowOrders_iFlows       ! Sets up indices in Flows%Orders for
                                        ! referring to flow module instances. Also
                                        ! sets up Flows%Orders%Included.
Public  :: SetUpFlowSubsets_iFlows      ! Sets up indices in Flows%Subsets for
                                        ! referring to flow module instances. Also
                                        ! sets up Flows%Subsets%Included.
Public  :: SetUpFlows_iFlowAttribs      ! Sets up indices in Flows for referring to
                                        ! flow attributes. Also checks that the flow
                                        ! module instances have the required
                                        ! attributes.
Public  :: SetUpFlows_iFlowSubsets      ! Sets up indices in the flow module instances
                                        ! for referring to flow subsets.
Public  :: SetUpFlows_iMets             ! Sets up indices in the flow module instances
                                        ! for referring to met module instances.
Public  :: SetUpFlows_iCoordsEtc        ! Sets up indices in the flow module instances
                                        ! for referring coord systems and grids.
Public  :: CheckFlows_FlowUpdateSubsets ! Checks that the update subsets used within
                                        ! the flow module instances only include flow
                                        ! module instances which are higher up the
                                        ! update priority order.
Public  :: FlowList                     ! Returns a list of the flow module instances.
Public  :: MetsFlowsOverallTValid       ! Returns the earliest time that the validity of any of the met or
                                        ! flow module instances might change, assuming all the instances which
                                        ! have been prepared for update-on-demand are updated now. The value
                                        ! is that determined at the time UpdateMetsFlows was last called (the
                                        ! actual time may be later).
Public  :: UpdateMetsFlows              ! Updates the state of the met and flow module
                                        ! instances.
Public  :: ResetFlowMemory              ! Resets the flow memory (for use when the
                                        ! space-time point which the memory applies to
                                        ! changes).
Public  :: ConvertToZ                   ! Converts the vertical coordinate of a point
                                        ! to a particular coord system.
Public  :: GetAttrib                    ! Gets attribute information.
Public  :: Reflect                      ! Reflects position of particle or puff
                                        ! centroid.
Public  :: CalcdZdZ                     ! Calculates dZCoord1/dZCoord2 for two
                                        ! vertical coord systems ZCoord1 and dZCoord2.
Public  :: XInTheDomain                 ! Finds out whether a point is in a domain
                                        ! (other than a flow module instance's
                                        ! domain).
Public  :: ResetMetsFlows               ! Resets the state of the met and flow module
                                        ! instances for a new realisation.
!$ Public  :: FlowLocks                 ! OpenMP locks for protecting FlowModule updates
!$ Public  :: MetLocks                  ! OpenMP locks for protecting MetModule updates
Public  :: CreateFlowMetLocks           ! Create Flow locks
Public  :: DestroyFlowMetLocks          ! Destroy Flow locks

! Items from the various flow modules which need to be made available by the flows module:
Public  :: PrototypeFlow_    ! Information describing the state of an instance of the Prototype
                        ! flow module.
Public  :: SingleSiteFlow_    ! Information describing the state of an instance of the SingleSite
                        ! flow module.
Public  :: NWPFlow_    ! Information describing the state of an instance of the NWP
                        ! flow module.
Public  :: BuildingFlow_    ! Information describing the state of an instance of the Building
                        ! flow module.
Public  :: RadarFlow_    ! Information describing the state of an instance of the Radar
                        ! flow module.
Public  :: LINCOMFlow_    ! Information describing the state of an instance of the LINCOM
                        ! flow module.

Public  :: InitPrototypeFlow ! Initialises an instance of PrototypeFlow_.
Public  :: InitSingleSiteFlow ! Initialises an instance of SingleSiteFlow_.
Public  :: InitNWPFlow ! Initialises an instance of NWPFlow_.
Public  :: InitBuildingFlow ! Initialises an instance of BuildingFlow_.
Public  :: InitRadarFlow ! Initialises an instance of RadarFlow_.
Public  :: InitLINCOMFlow ! Initialises an instance of LINCOMFlow_.

! Items from the common flow module which need to be made available by the flows module:
Public  :: A_Update  !} Codes for attributes.
Public  :: A_Convert !}
Public  :: A_Flow    !}
Public  :: A_Cloud   !}
Public  :: A_Rain    !}
Public  :: A_Surface !}
Public  :: A_Soil    !}
Public  :: A_Plant   !}

! Items from the mets module which need to be made available by the flows module:
Public  :: Mets_                ! A collection of met module instance states.
Public  :: InitMets             ! Initialises a collection of met module instance
                                ! states.

Public  :: AddPrototypeMet           ! Adds the state of a Prototype met module instance to a
                                ! collection of met module instance states.
Public  :: AddSingleSiteMet           ! Adds the state of a SingleSite met module instance to a
                                ! collection of met module instance states.
Public  :: AddNWPMet           ! Adds the state of a NWP met module instance to a
                                ! collection of met module instance states.
Public  :: AddRadarMet           ! Adds the state of a Radar met module instance to a
                                ! collection of met module instance states.
Public  :: AddAncillaryMet           ! Adds the state of a Ancillary met module instance to a
                                ! collection of met module instance states.

Public  :: SetUpCoordsEtc_Mets  ! Sets up Coords and Grids by adding any extra coords
                                ! and grids which Mets wants to define.
Public  :: SetUpMets_CoordsEtc  ! Sets up Mets using information from EtaDefns, Coords
                                ! and Grids.
Public  :: SetUpMets_iCoordsEtc ! Sets up indices in the met module instances for
                                ! referring coord systems and grids.

Public  :: PrototypeMet_             ! Information describing the state of an instance of
                                ! the Prototype met module. This is defined in the
                                ! PrototypeMetModule.
Public  :: SingleSiteMet_             ! Information describing the state of an instance of
                                ! the SingleSite met module. This is defined in the
                                ! SingleSiteMetModule.
Public  :: NWPMet_             ! Information describing the state of an instance of
                                ! the NWP met module. This is defined in the
                                ! NWPMetModule.
Public  :: RadarMet_             ! Information describing the state of an instance of
                                ! the Radar met module. This is defined in the
                                ! RadarMetModule.
Public  :: AncillaryMet_             ! Information describing the state of an instance of
                                ! the Ancillary met module. This is defined in the
                                ! AncillaryMetModule.

Public  :: InitPrototypeMet          ! Initialises an instance of PrototypeMet_. This is defined
                                ! in the PrototypeMetModule.
Public  :: InitSingleSiteMet          ! Initialises an instance of SingleSiteMet_. This is defined
                                ! in the SingleSiteMetModule.
Public  :: InitNWPMet          ! Initialises an instance of NWPMet_. This is defined
                                ! in the NWPMetModule.
Public  :: InitRadarMet          ! Initialises an instance of RadarMet_. This is defined
                                ! in the RadarMetModule.
Public  :: InitAncillaryMet          ! Initialises an instance of AncillaryMet_. This is defined
                                ! in the AncillaryMetModule.

! Items from the flow and flow profile module which need to be made available by the flows module:
Public  :: Flow_        ! Flow properties at a particular location, plus some bulk
                        ! boundary layer properties. Used to return flow information
                        ! from the flow modules or, in conjunction with FlowField_, to
                        ! transfer information between the flow modules.
Public  :: Cloud_       ! Cloud properties at a particular location. Used to return
                        ! cloud information from the flow modules or, in conjunction
                        ! with CloudField_, to transfer information between the flow
                        ! modules.
Public  :: Rain_        ! Rain properties at a particular location. Used to return
                        ! rain information from the flow modules or, in conjunction
                        ! with RainField_, to transfer information between the flow
                        ! modules.
Public  :: Surface_     ! Surface properties at a particular location. Used to return surface information from
                        ! the flow modules or, in conjunction with SurfaceField_, to transfer information
                        ! between the flow modules.
Public  :: Soil_        ! Soil properties at a particular location. Used to return soil information from
                        ! the flow modules or, in conjunction with SoilField_, to transfer information
                        ! between the flow modules.
Public  :: Plant_       ! Plant properties at a particular location. Used to return plant information from
                        ! the flow modules or, in conjunction with PlantField_, to transfer information
                        ! between the flow modules.
Public  :: ProfileData_ ! Data needed to contruct idealised analytic mean flow and
                        ! turbulence profiles. Used internally within the flow
                        ! modules, or to return flow information from the flow modules
                        ! or, in conjunction with FlowField_, to transfer information
                        ! between the flow modules.
Public  :: ConvertFlow  ! Converts flow information between different coord systems.
Public  :: ReverseFlow  ! Time reverses the flow.

!-------------------------------------------------------------------------------------------------------------

Type :: FlowOrder_ ! An ordered subset of flow module instances.
  Private
  Character(MaxCharLength) :: Name
  Integer                  :: nFlows
  Character(MaxCharLength) :: FlowModNames(MaxFlows)
  Character(MaxCharLength) :: FlowNames(MaxFlows)
  Integer                  :: iFlowMods(MaxFlows)
  Integer                  :: iFlows(MaxFlows)
  Logical                  :: Included(MaxFlowMods, MaxFlowsPerMod)
  ! Name         :: Name of order.
  ! nFlows       :: Number of flow module instances.
  ! FlowModNames :} Names and indices of the flow modules and flow module instances
  ! FlowNames    :} which are included in the subset, ordered in order of decreasing
  ! iFlowMods    :} priority.
  ! iFlows       :}
  ! Included     :: Indicates which flow module instances are included.
End Type FlowOrder_

!-------------------------------------------------------------------------------------------------------------

Type :: FlowSubset_ ! A subset of flow module instances.
  Private
  Character(MaxCharLength) :: Name
  Integer                  :: nFlows
  Character(MaxCharLength) :: FlowModNames(MaxFlows)
  Character(MaxCharLength) :: FlowNames(MaxFlows)
  Integer                  :: iFlowMods(MaxFlows)
  Integer                  :: iFlows(MaxFlows)
  Logical                  :: Included(MaxFlowMods, MaxFlowsPerMod)
  ! Name         :: Name of subset.
  ! nFlows       :: Number of flow module instances.
  ! FlowModNames :} Names and indices of the flow modules and flow module instances
  ! FlowNames    :} which are included in the subset.
  ! iFlowMods    :}
  ! iFlows       :}
  ! Included     :: Indicates which flow module instances are included.
End Type FlowSubset_

!-------------------------------------------------------------------------------------------------------------

Type :: FlowAttrib_ ! A flow attribute.
  Private
  Character(MaxCharLength) :: Name         ! Name of attribute.
  Integer                  :: iAttribParam ! Index of attribute (in the parameter list in the common flow
                                           ! module).
  Character(MaxCharLength) :: OrderName    !} Name and index of ordered subset of flow module instances giving
  Integer                  :: iOrder       !} the priority order for the attribute.
End Type FlowAttrib_

!-------------------------------------------------------------------------------------------------------------

Type :: CommonFlowP_ ! A pointer to an instance of the type CommonFlow_.
  Private
  Type(CommonFlow_), Pointer :: P ! Pointer to an instance of the type CommonFlow_.
End Type CommonFlowP_

!-------------------------------------------------------------------------------------------------------------

Type :: Flows_ ! A collection of flow module instance states.
  Private
  Integer            :: nFlowMods
  Integer            :: nFlows(MaxFlowMods)
  Type(CommonFlowP_) :: C(MaxFlowMods, MaxFlowsPerMod)
  Integer            :: nOrders
  Type(FlowOrder_)   :: Orders(MaxFlowOrders)
  Integer            :: nSubsets
  Type(FlowSubset_)  :: Subsets(MaxFlowSubsets)
  Integer            :: nAttribs
  Type(FlowAttrib_)  :: Attribs(MaxFlowAttribs)
  Integer            :: iAttribs(Size(AttribParamNames))
  Integer            :: iZCoordMagl
  Integer            :: iZCoordPa
  Logical            :: AnyUpdateOnDemand
  Logical            :: SameResultsWithUpdateOnDemand
  Integer            :: iCase
  Integer            :: iMetCase
  Type(ShortTime_)   :: MetTime
  Type(Time_)        :: OverallTValid

  Integer            :: nPrototypeFlows
  Type(PrototypeFlow_)    :: PrototypeFlows(MaxPrototypeFlows)
  Integer            :: nSingleSiteFlows
  Type(SingleSiteFlow_)    :: SingleSiteFlows(MaxSingleSiteFlows)
  Integer            :: nNWPFlows
  Type(NWPFlow_)    :: NWPFlows(MaxNWPFlows)
  Integer            :: nBuildingFlows
  Type(BuildingFlow_)    :: BuildingFlows(MaxBuildingFlows)
  Integer            :: nRadarFlows
  Type(RadarFlow_)    :: RadarFlows(MaxRadarFlows)
  Integer            :: nLINCOMFlows
  Type(LINCOMFlow_)    :: LINCOMFlows(MaxLINCOMFlows)

  ! nFlowMods                     :: Number of flow modules.
  ! nFlows                        :: Number of instances of each flow module.
  ! C                             :: Pointers to the common parts of the flow module instance states.
  ! nOrders                       :: Number of ordered subsets of flow module instances.
  ! Orders                        :: Ordered subsets of flow module instances.
  ! nSubsets                      :: Number of subsets of flow module instances.
  ! Subsets                       :: Subsets of flow module instances.
  ! nAttribs                      :: Number of attributes.
  ! Attribs                       :: Attributes.
  ! iAttribs                      :: Attributes are indexed in two ways - via the parameters in the common
  !                                  flow module and via the Attribs array. This array gives the index in the
  !                                  Attribs array of a given attribute in the parameters in the common flow
  !                                  module. It is set to zero if the attribute is not in the Attribs array.
  !                                  Note the inverse function is available as part of Attribs.
  ! iZCoordMagl                   :} Indices of vertical coord systems m agl and Pa.
  ! iZCoordPa                     :}
  ! AnyUpdateOnDemand             :: Indicates if any met and flow modules are to use update-on-demand.
  ! SameResultsWithUpdateOnDemand :: Indicates results should be made the same whether the met and flow module
  !                                  instances use update-on-demand or update-at-once (so it should really be
  !                                  called SameResultsWithOrWithoutUpdateOnDemand).
  ! iCase                         :} Values of iCase, iMetCase and MetTime on last call to UpdateMetsFlows.
  ! iMetCase                      :}
  ! MetTime                       :}
  ! OverallTValid                 :: Earliest time that the validity of any of the met or flow module
  !                                  instances might change, assuming all the instances which have been
  !                                  prepared for update-on-demand are updated now. The value is that
  !                                  determined at the time UpdateMetsFlows was last called (the actual time
  !                                  may be later).
  ! n????Flows                    :: Number of instances of the ???? flow module.
  ! ????Flows                     :: States of all instances of the ???? flow module.
End Type Flows_

!-------------------------------------------------------------------------------------------------------------

Type :: FlowMemory_ ! Information related to a particular (space, time)-point or
                    ! (space, time, travel time)-point of interest, including
                    ! information on the flow module instance domains the point lies
                    ! within, the choices of flow module instances appropriate, and
                    ! anything the individual flow module instances wish to record.
                    ! The purpose of storing this information together is to provide a
                    ! memory to avoid the need to recalculate the information.
  Private
  Logical               :: In(MaxFlowMods, MaxFlowsPerMod)
  Logical               :: InValid(MaxFlowMods, MaxFlowsPerMod)
  Integer               :: iFlowMod(MaxFlowAttribs)
  Integer               :: iFlow(MaxFlowAttribs)
  Logical               :: iValid(MaxFlowAttribs)

  Type(NWPFlowMemory_) :: NWPFlowMemories(MaxNWPFlows)

  ! In               :} Indicates whether the point in question lies in the domain of
  ! InValid          :} each of the flow module instances and whether this information
  !                     is valid.
  ! iFlowMod         :] Indices of the flow modules and flow module instances
  ! iFlow            :] appropriate for use at the point in question for each of the
  ! iValid           :] attributes and whether this information is valid.
  ! ????FlowMemories :: Information which the ???? flow module instances wish to
  !                     record.
End Type FlowMemory_

!-------------------------------------------------------------------------------------------------------------

! Maximum number of threads
Integer, Parameter :: MaxThreads = 8

! Locks

!$ Integer(kind=OMP_lock_kind) :: FlowLocks(MaxFlowMods,MaxFlowsPerMod)
!$ Integer(kind=OMP_lock_kind) :: MetLocks(MaxMetMods,MaxMetsPerMod)

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitFlows(SameResultsWithUpdateOnDemand) Result(Flows)
! Initialises a collection of flow module instance states.

  Implicit None
  Logical :: SameResultsWithUpdateOnDemand ! Indicates results should be made the same whether the met and
                                           ! flow module instances use update-on-demand or update-at-once (so
                                           ! it should really be called
                                           ! SameResultsWithOrWithoutUpdateOnDemand).
  ! Function result:
  Type(Flows_) :: Flows ! The initialised collection of flow module instance states.

  Flows%nFlowMods = 0
  Flows%nOrders   = 0
  Flows%nSubsets  = 0
  Flows%nAttribs  = 0

  Flows%SameResultsWithUpdateOnDemand = SameResultsWithUpdateOnDemand

  Flows%nPrototypeFlows = 0
  Flows%nSingleSiteFlows = 0
  Flows%nNWPFlows = 0
  Flows%nBuildingFlows = 0
  Flows%nRadarFlows = 0
  Flows%nLINCOMFlows = 0

End Function InitFlows

!-------------------------------------------------------------------------------------------------------------


Subroutine AddPrototypeFlow(PrototypeFlow, Flows)
! Adds the state of a Prototype flow module instance to a collection of flow module
! instance states.

  Implicit None
  ! Argument list:
  Type(PrototypeFlow_), Intent(In)            :: PrototypeFlow ! The state of the Prototype flow
                                                     ! module instance.
  Type(Flows_),    Intent(InOut), Target :: Flows    ! The collection of flow module
                                                     ! instance states.
  ! Locals:
  Integer :: i        ! Loop index.
  Integer :: j        ! Loop index.
  Integer :: iFlowMod ! Index of flow module to be added.
  Integer :: iFlow    ! Index of flow module instance to be added.

  ! Calculate iFlowMod and, if necessary, check for too many flow modules, update
  ! Flows%nFlowMods and initialise Flows%nFlows(iFlowMod).
  If (Flows%nPrototypeFlows == 0) Then
    If (Flows%nFlowMods == MaxFlowMods) Then
      Call Message('FATAL ERROR in AddPrototypeFlow: Too many flow modules', 3)
    End If
    iFlowMod               = Flows%nFlowMods + 1
    Flows%nFlowMods        = iFlowMod
    Flows%nFlows(iFlowMod) = 0
  Else
    iFlowMod = Flows%PrototypeFlows(1)%C%iFlowMod
  End If

  ! Check for too many flow module instances of this type and for duplicate names.
  If (Flows%nFlows(iFlowMod) == MaxPrototypeFlows) Then
    Call Message('FATAL ERROR in AddPrototypeFlow: Too many instances of the Prototype flow module', 3)
  End If
  Do i = 1, Flows%nFlowMods
  Do j = 1, Flows%nFlows(i)
    If (Flows%C(i, j)%P%FlowName .CIEq. PrototypeFlow%C%FlowName) Then
      Call Message(                                                            &
             'FATAL ERROR in AddPrototypeFlow: a flow module instance with the ' // &
             'same name already exists',                                       &
             3                                                                 &
           )
    End If
  End Do
  End Do

  ! Calculate iFlow.
  iFlow = Flows%nFlows(iFlowMod) + 1

  ! Add flow module instance.
  Flows%nFlows(iFlowMod)            = iFlow
  Flows%nPrototypeFlows                  = iFlow
  Flows%PrototypeFlows(iFlow)            = PrototypeFlow
  Flows%PrototypeFlows(iFlow)%C%iFlowMod = iFlowMod
  Flows%PrototypeFlows(iFlow)%C%iFlow    = iFlow
  Flows%C(iFlowMod, iFlow)%P        => Flows%PrototypeFlows(iFlow)%C

  ! Check TValid is the right type of time.
  If (IsTimeInterval(Flows%C(iFlowMod, iFlow)%P%TValid)) Then
    Call Message('UNEXPECTED FATAL ERROR in AddPrototypeFlow', 4)
  End If

  ! Check consistency of FixedMet between flow module instances.
  If (Flows%C(1, 1)%P%FixedMet .neqv. Flows%C(iFlowMod, iFlow)%P%FixedMet) Then
    Call Message('UNEXPECTED FATAL ERROR in AddPrototypeFlow', 4)
  End If

End Subroutine AddPrototypeFlow


Subroutine AddSingleSiteFlow(SingleSiteFlow, Flows)
! Adds the state of a SingleSite flow module instance to a collection of flow module
! instance states.

  Implicit None
  ! Argument list:
  Type(SingleSiteFlow_), Intent(In)            :: SingleSiteFlow ! The state of the SingleSite flow
                                                     ! module instance.
  Type(Flows_),    Intent(InOut), Target :: Flows    ! The collection of flow module
                                                     ! instance states.
  ! Locals:
  Integer :: i        ! Loop index.
  Integer :: j        ! Loop index.
  Integer :: iFlowMod ! Index of flow module to be added.
  Integer :: iFlow    ! Index of flow module instance to be added.

  ! Calculate iFlowMod and, if necessary, check for too many flow modules, update
  ! Flows%nFlowMods and initialise Flows%nFlows(iFlowMod).
  If (Flows%nSingleSiteFlows == 0) Then
    If (Flows%nFlowMods == MaxFlowMods) Then
      Call Message('FATAL ERROR in AddSingleSiteFlow: Too many flow modules', 3)
    End If
    iFlowMod               = Flows%nFlowMods + 1
    Flows%nFlowMods        = iFlowMod
    Flows%nFlows(iFlowMod) = 0
  Else
    iFlowMod = Flows%SingleSiteFlows(1)%C%iFlowMod
  End If

  ! Check for too many flow module instances of this type and for duplicate names.
  If (Flows%nFlows(iFlowMod) == MaxSingleSiteFlows) Then
    Call Message('FATAL ERROR in AddSingleSiteFlow: Too many instances of the SingleSite flow module', 3)
  End If
  Do i = 1, Flows%nFlowMods
  Do j = 1, Flows%nFlows(i)
    If (Flows%C(i, j)%P%FlowName .CIEq. SingleSiteFlow%C%FlowName) Then
      Call Message(                                                            &
             'FATAL ERROR in AddSingleSiteFlow: a flow module instance with the ' // &
             'same name already exists',                                       &
             3                                                                 &
           )
    End If
  End Do
  End Do

  ! Calculate iFlow.
  iFlow = Flows%nFlows(iFlowMod) + 1

  ! Add flow module instance.
  Flows%nFlows(iFlowMod)            = iFlow
  Flows%nSingleSiteFlows                  = iFlow
  Flows%SingleSiteFlows(iFlow)            = SingleSiteFlow
  Flows%SingleSiteFlows(iFlow)%C%iFlowMod = iFlowMod
  Flows%SingleSiteFlows(iFlow)%C%iFlow    = iFlow
  Flows%C(iFlowMod, iFlow)%P        => Flows%SingleSiteFlows(iFlow)%C

  ! Check TValid is the right type of time.
  If (IsTimeInterval(Flows%C(iFlowMod, iFlow)%P%TValid)) Then
    Call Message('UNEXPECTED FATAL ERROR in AddSingleSiteFlow', 4)
  End If

  ! Check consistency of FixedMet between flow module instances.
  If (Flows%C(1, 1)%P%FixedMet .neqv. Flows%C(iFlowMod, iFlow)%P%FixedMet) Then
    Call Message('UNEXPECTED FATAL ERROR in AddSingleSiteFlow', 4)
  End If

End Subroutine AddSingleSiteFlow


Subroutine AddNWPFlow(NWPFlow, Flows)
! Adds the state of a NWP flow module instance to a collection of flow module
! instance states.

  Implicit None
  ! Argument list:
  Type(NWPFlow_), Intent(In)            :: NWPFlow ! The state of the NWP flow
                                                     ! module instance.
  Type(Flows_),    Intent(InOut), Target :: Flows    ! The collection of flow module
                                                     ! instance states.
  ! Locals:
  Integer :: i        ! Loop index.
  Integer :: j        ! Loop index.
  Integer :: iFlowMod ! Index of flow module to be added.
  Integer :: iFlow    ! Index of flow module instance to be added.

  ! Calculate iFlowMod and, if necessary, check for too many flow modules, update
  ! Flows%nFlowMods and initialise Flows%nFlows(iFlowMod).
  If (Flows%nNWPFlows == 0) Then
    If (Flows%nFlowMods == MaxFlowMods) Then
      Call Message('FATAL ERROR in AddNWPFlow: Too many flow modules', 3)
    End If
    iFlowMod               = Flows%nFlowMods + 1
    Flows%nFlowMods        = iFlowMod
    Flows%nFlows(iFlowMod) = 0
  Else
    iFlowMod = Flows%NWPFlows(1)%C%iFlowMod
  End If

  ! Check for too many flow module instances of this type and for duplicate names.
  If (Flows%nFlows(iFlowMod) == MaxNWPFlows) Then
    Call Message('FATAL ERROR in AddNWPFlow: Too many instances of the NWP flow module', 3)
  End If
  Do i = 1, Flows%nFlowMods
  Do j = 1, Flows%nFlows(i)
    If (Flows%C(i, j)%P%FlowName .CIEq. NWPFlow%C%FlowName) Then
      Call Message(                                                            &
             'FATAL ERROR in AddNWPFlow: a flow module instance with the ' // &
             'same name already exists',                                       &
             3                                                                 &
           )
    End If
  End Do
  End Do

  ! Calculate iFlow.
  iFlow = Flows%nFlows(iFlowMod) + 1

  ! Add flow module instance.
  Flows%nFlows(iFlowMod)            = iFlow
  Flows%nNWPFlows                  = iFlow
  Flows%NWPFlows(iFlow)            = NWPFlow
  Flows%NWPFlows(iFlow)%C%iFlowMod = iFlowMod
  Flows%NWPFlows(iFlow)%C%iFlow    = iFlow
  Flows%C(iFlowMod, iFlow)%P        => Flows%NWPFlows(iFlow)%C

  ! Check TValid is the right type of time.
  If (IsTimeInterval(Flows%C(iFlowMod, iFlow)%P%TValid)) Then
    Call Message('UNEXPECTED FATAL ERROR in AddNWPFlow', 4)
  End If

  ! Check consistency of FixedMet between flow module instances.
  If (Flows%C(1, 1)%P%FixedMet .neqv. Flows%C(iFlowMod, iFlow)%P%FixedMet) Then
    Call Message('UNEXPECTED FATAL ERROR in AddNWPFlow', 4)
  End If

End Subroutine AddNWPFlow


Subroutine AddBuildingFlow(BuildingFlow, Flows)
! Adds the state of a Building flow module instance to a collection of flow module
! instance states.

  Implicit None
  ! Argument list:
  Type(BuildingFlow_), Intent(In)            :: BuildingFlow ! The state of the Building flow
                                                     ! module instance.
  Type(Flows_),    Intent(InOut), Target :: Flows    ! The collection of flow module
                                                     ! instance states.
  ! Locals:
  Integer :: i        ! Loop index.
  Integer :: j        ! Loop index.
  Integer :: iFlowMod ! Index of flow module to be added.
  Integer :: iFlow    ! Index of flow module instance to be added.

  ! Calculate iFlowMod and, if necessary, check for too many flow modules, update
  ! Flows%nFlowMods and initialise Flows%nFlows(iFlowMod).
  If (Flows%nBuildingFlows == 0) Then
    If (Flows%nFlowMods == MaxFlowMods) Then
      Call Message('FATAL ERROR in AddBuildingFlow: Too many flow modules', 3)
    End If
    iFlowMod               = Flows%nFlowMods + 1
    Flows%nFlowMods        = iFlowMod
    Flows%nFlows(iFlowMod) = 0
  Else
    iFlowMod = Flows%BuildingFlows(1)%C%iFlowMod
  End If

  ! Check for too many flow module instances of this type and for duplicate names.
  If (Flows%nFlows(iFlowMod) == MaxBuildingFlows) Then
    Call Message('FATAL ERROR in AddBuildingFlow: Too many instances of the Building flow module', 3)
  End If
  Do i = 1, Flows%nFlowMods
  Do j = 1, Flows%nFlows(i)
    If (Flows%C(i, j)%P%FlowName .CIEq. BuildingFlow%C%FlowName) Then
      Call Message(                                                            &
             'FATAL ERROR in AddBuildingFlow: a flow module instance with the ' // &
             'same name already exists',                                       &
             3                                                                 &
           )
    End If
  End Do
  End Do

  ! Calculate iFlow.
  iFlow = Flows%nFlows(iFlowMod) + 1

  ! Add flow module instance.
  Flows%nFlows(iFlowMod)            = iFlow
  Flows%nBuildingFlows                  = iFlow
  Flows%BuildingFlows(iFlow)            = BuildingFlow
  Flows%BuildingFlows(iFlow)%C%iFlowMod = iFlowMod
  Flows%BuildingFlows(iFlow)%C%iFlow    = iFlow
  Flows%C(iFlowMod, iFlow)%P        => Flows%BuildingFlows(iFlow)%C

  ! Check TValid is the right type of time.
  If (IsTimeInterval(Flows%C(iFlowMod, iFlow)%P%TValid)) Then
    Call Message('UNEXPECTED FATAL ERROR in AddBuildingFlow', 4)
  End If

  ! Check consistency of FixedMet between flow module instances.
  If (Flows%C(1, 1)%P%FixedMet .neqv. Flows%C(iFlowMod, iFlow)%P%FixedMet) Then
    Call Message('UNEXPECTED FATAL ERROR in AddBuildingFlow', 4)
  End If

End Subroutine AddBuildingFlow


Subroutine AddRadarFlow(RadarFlow, Flows)
! Adds the state of a Radar flow module instance to a collection of flow module
! instance states.

  Implicit None
  ! Argument list:
  Type(RadarFlow_), Intent(In)            :: RadarFlow ! The state of the Radar flow
                                                     ! module instance.
  Type(Flows_),    Intent(InOut), Target :: Flows    ! The collection of flow module
                                                     ! instance states.
  ! Locals:
  Integer :: i        ! Loop index.
  Integer :: j        ! Loop index.
  Integer :: iFlowMod ! Index of flow module to be added.
  Integer :: iFlow    ! Index of flow module instance to be added.

  ! Calculate iFlowMod and, if necessary, check for too many flow modules, update
  ! Flows%nFlowMods and initialise Flows%nFlows(iFlowMod).
  If (Flows%nRadarFlows == 0) Then
    If (Flows%nFlowMods == MaxFlowMods) Then
      Call Message('FATAL ERROR in AddRadarFlow: Too many flow modules', 3)
    End If
    iFlowMod               = Flows%nFlowMods + 1
    Flows%nFlowMods        = iFlowMod
    Flows%nFlows(iFlowMod) = 0
  Else
    iFlowMod = Flows%RadarFlows(1)%C%iFlowMod
  End If

  ! Check for too many flow module instances of this type and for duplicate names.
  If (Flows%nFlows(iFlowMod) == MaxRadarFlows) Then
    Call Message('FATAL ERROR in AddRadarFlow: Too many instances of the Radar flow module', 3)
  End If
  Do i = 1, Flows%nFlowMods
  Do j = 1, Flows%nFlows(i)
    If (Flows%C(i, j)%P%FlowName .CIEq. RadarFlow%C%FlowName) Then
      Call Message(                                                            &
             'FATAL ERROR in AddRadarFlow: a flow module instance with the ' // &
             'same name already exists',                                       &
             3                                                                 &
           )
    End If
  End Do
  End Do

  ! Calculate iFlow.
  iFlow = Flows%nFlows(iFlowMod) + 1

  ! Add flow module instance.
  Flows%nFlows(iFlowMod)            = iFlow
  Flows%nRadarFlows                  = iFlow
  Flows%RadarFlows(iFlow)            = RadarFlow
  Flows%RadarFlows(iFlow)%C%iFlowMod = iFlowMod
  Flows%RadarFlows(iFlow)%C%iFlow    = iFlow
  Flows%C(iFlowMod, iFlow)%P        => Flows%RadarFlows(iFlow)%C

  ! Check TValid is the right type of time.
  If (IsTimeInterval(Flows%C(iFlowMod, iFlow)%P%TValid)) Then
    Call Message('UNEXPECTED FATAL ERROR in AddRadarFlow', 4)
  End If

  ! Check consistency of FixedMet between flow module instances.
  If (Flows%C(1, 1)%P%FixedMet .neqv. Flows%C(iFlowMod, iFlow)%P%FixedMet) Then
    Call Message('UNEXPECTED FATAL ERROR in AddRadarFlow', 4)
  End If

End Subroutine AddRadarFlow


Subroutine AddLINCOMFlow(LINCOMFlow, Flows)
! Adds the state of a LINCOM flow module instance to a collection of flow module
! instance states.

  Implicit None
  ! Argument list:
  Type(LINCOMFlow_), Intent(In)            :: LINCOMFlow ! The state of the LINCOM flow
                                                     ! module instance.
  Type(Flows_),    Intent(InOut), Target :: Flows    ! The collection of flow module
                                                     ! instance states.
  ! Locals:
  Integer :: i        ! Loop index.
  Integer :: j        ! Loop index.
  Integer :: iFlowMod ! Index of flow module to be added.
  Integer :: iFlow    ! Index of flow module instance to be added.

  ! Calculate iFlowMod and, if necessary, check for too many flow modules, update
  ! Flows%nFlowMods and initialise Flows%nFlows(iFlowMod).
  If (Flows%nLINCOMFlows == 0) Then
    If (Flows%nFlowMods == MaxFlowMods) Then
      Call Message('FATAL ERROR in AddLINCOMFlow: Too many flow modules', 3)
    End If
    iFlowMod               = Flows%nFlowMods + 1
    Flows%nFlowMods        = iFlowMod
    Flows%nFlows(iFlowMod) = 0
  Else
    iFlowMod = Flows%LINCOMFlows(1)%C%iFlowMod
  End If

  ! Check for too many flow module instances of this type and for duplicate names.
  If (Flows%nFlows(iFlowMod) == MaxLINCOMFlows) Then
    Call Message('FATAL ERROR in AddLINCOMFlow: Too many instances of the LINCOM flow module', 3)
  End If
  Do i = 1, Flows%nFlowMods
  Do j = 1, Flows%nFlows(i)
    If (Flows%C(i, j)%P%FlowName .CIEq. LINCOMFlow%C%FlowName) Then
      Call Message(                                                            &
             'FATAL ERROR in AddLINCOMFlow: a flow module instance with the ' // &
             'same name already exists',                                       &
             3                                                                 &
           )
    End If
  End Do
  End Do

  ! Calculate iFlow.
  iFlow = Flows%nFlows(iFlowMod) + 1

  ! Add flow module instance.
  Flows%nFlows(iFlowMod)            = iFlow
  Flows%nLINCOMFlows                  = iFlow
  Flows%LINCOMFlows(iFlow)            = LINCOMFlow
  Flows%LINCOMFlows(iFlow)%C%iFlowMod = iFlowMod
  Flows%LINCOMFlows(iFlow)%C%iFlow    = iFlow
  Flows%C(iFlowMod, iFlow)%P        => Flows%LINCOMFlows(iFlow)%C

  ! Check TValid is the right type of time.
  If (IsTimeInterval(Flows%C(iFlowMod, iFlow)%P%TValid)) Then
    Call Message('UNEXPECTED FATAL ERROR in AddLINCOMFlow', 4)
  End If

  ! Check consistency of FixedMet between flow module instances.
  If (Flows%C(1, 1)%P%FixedMet .neqv. Flows%C(iFlowMod, iFlow)%P%FixedMet) Then
    Call Message('UNEXPECTED FATAL ERROR in AddLINCOMFlow', 4)
  End If

End Subroutine AddLINCOMFlow


!-------------------------------------------------------------------------------------------------------------

Subroutine FindFlowIndex(FlowModName, FlowName, Flows, iFlowMod, iFlow)
! Finds the indices of a flow module and a flow module instance.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)  :: FlowModName ! Flow module name.
  Character(*), Intent(In)  :: FlowName    ! Flow module instance name.
  Type(Flows_), Intent(In)  :: Flows       ! Collection of flow module instance
                                           ! states.
  Integer,      Intent(Out) :: iFlowMod    ! Flow module index.
  Integer,      Intent(Out) :: iFlow       ! Flow module instance index.
  ! Locals:
  Integer :: i ! Loop index.

  iFlowMod = 0
  Do i = 1, Flows%nFlowMods
    If (Flows%nFlows(i) > 0) Then
      If (FlowModName .CIEq. Flows%C(i, 1)%P%FlowModName) Then
        iFlowMod = i
        Exit
      End If
    End If
  End Do
  If (iFlowMod == 0) Then
    Call Message('FATAL ERROR in FindFlowIndex: flow module not found', 3)
  End If

  iFlow = 0
  Do i = 1, Flows%nFlows(iFlowMod)
    If (FlowName .CIEq. Flows%C(iFlowMod, i)%P%FlowName) Then
      iFlow = i
      Exit
    End If
  End Do
  If (iFlow == 0) Then
    Call Message('FATAL ERROR in FindFlowIndex: flow module instance not found', 3)
  End If

End Subroutine FindFlowIndex

!-------------------------------------------------------------------------------------------------------------

Function InitFlowOrder(Name, FlowModName, FlowName) Result(FlowOrder)
! Initialises an ordered subset of flow module instances.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name           ! Name of ordered subset of flow module
                                             ! instances.
  Character(*), Intent(In) :: FlowModName(:) ! Flow module names.
  Character(*), Intent(In) :: FlowName(:)    ! Flow module instance names.
  ! Function result:
  Type(FlowOrder_) :: FlowOrder ! Ordered subset of flow module instances.
  ! Locals:
  Integer :: i ! Loop index.

  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message('FATAL ERROR in InitFlowOrder', 3)
  End If
  If (Len_Trim(Name) == 0) Then
    Call Message('FATAL ERROR in InitFlowOrder', 3)
  End If

  If (Size(FlowModName) /= Size(FlowName)) Then
    Call Message('UNEXPECTED FATAL ERROR in InitFlowOrder', 4)
  End If
  If (Size(FlowName) > MaxFlows) Then
    Call Message('FATAL ERROR in InitFlowOrder', 3)
  End If

  Do i = 1, Size(FlowName)
    If (Len_Trim(FlowModName(i)) > MaxCharLength) Then
      Call Message('FATAL ERROR in InitFlowOrder', 3)
    End If
    If (Len_Trim(FlowModName(i)) == 0) Then
      Call Message('FATAL ERROR in InitFlowOrder', 3)
    End If
    If (Len_Trim(FlowName(i)) > MaxCharLength) Then
      Call Message('FATAL ERROR in InitFlowOrder', 3)
    End If
    If (Len_Trim(FlowName(i)) == 0) Then
      Call Message('FATAL ERROR in InitFlowOrder', 3)
    End If
  End Do

  FlowOrder%Name                           = Name
  FlowOrder%nFlows                         = Size(FlowName)
  FlowOrder%FlowModNames(1:Size(FlowName)) = FlowModName(1:Size(FlowName))
  FlowOrder%FlowNames   (1:Size(FlowName)) = FlowName   (1:Size(FlowName))

End Function InitFlowOrder

!-------------------------------------------------------------------------------------------------------------

Subroutine AddFlowOrder(FlowOrder, Flows)
! Adds an ordered subset of flow module instances to a collection of flow module
! instance states.

  Implicit None
  ! Argument list:
  Type(FlowOrder_), Intent(In)    :: FlowOrder ! Ordered subset of flow module
                                               ! instances.
  Type(Flows_),     Intent(InOut) :: Flows     ! Collection of flow module instance
                                               ! states.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Flows%nOrders
    If (FlowOrder%Name .CIEq. Flows%Orders(i)%Name) Then
      If (FlowOrderEq(FlowOrder, Flows%Orders(i))) Then
        Return
      Else
        Call Message(                                                             &
               'FATAL ERROR in adding Flow Order "'                            // &
               Trim(FlowOrder%Name)                                            // &
               '": an different flow order with the same name already exists',    &
               3                                                                  &
             )
      End If
    End If
  End Do

  If (Flows%nOrders >= MaxFlowOrders) Then
    Call Message(                                    &
           'FATAL ERROR in adding a Flow Order: ' // &
           'there are too many Flow Orders ',        &
           3                                         &
         )
  End If

  Flows%nOrders               = Flows%nOrders + 1
  Flows%Orders(Flows%nOrders) = FlowOrder

End Subroutine AddFlowOrder

!-------------------------------------------------------------------------------------------------------------

Function FlowOrderEq(FlowOrder1, FlowOrder2)
! Tests for equality of ordered subsets of flow module instances.

  Implicit None
  ! Argument list:
  Type(FlowOrder_), Intent(In) :: FlowOrder1 !} The two ordered subsets of flow module
  Type(FlowOrder_), Intent(In) :: FlowOrder2 !} instances.
  ! Function result:
  Logical :: FlowOrderEq ! Indicates if ordered subsets of flow module instances are
                         ! equal.
  ! Locals:
  Integer :: i ! Loop index.

  FlowOrderEq = (FlowOrder1%Name   .CIEq. FlowOrder2%Name)  .and. &
                 FlowOrder1%nFlows   ==   FlowOrder2%nFlows
  Do i = 1, Min(FlowOrder1%nFlows, FlowOrder2%nFlows)
    FlowOrderEq =                                                          &
      FlowOrderEq                                                    .and. &
      (FlowOrder1%FlowModNames(i) .CIEq. FlowOrder2%FlowModNames(i)) .and. &
      (FlowOrder1%FlowNames(i)    .CIEq. FlowOrder2%FlowNames(i))
  End Do

End Function FlowOrderEq

!-------------------------------------------------------------------------------------------------------------

Function FindOrderIndex(Name, Flows)
! Finds the index of an ordered subset of flow module instances.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name  ! Name of ordered subset of flow module instances.
  Type(Flows_), Intent(In) :: Flows ! Collection of flow module instance states.
  ! Function result:
  Integer :: FindOrderIndex ! Index of the ordered subset of flow module instances.
  ! Locals:
  Integer :: i ! Loop index.

  FindOrderIndex = 0
  Do i = 1, Flows%nOrders
    If (Name .CIEq. Flows%Orders(i)%Name) Then
      FindOrderIndex = i
      Exit
    End If
  End Do
  If (FindOrderIndex == 0) Then
    Call Message('FATAL ERROR in FindOrderIndex: order not found', 3)
  End If

End Function FindOrderIndex

!-------------------------------------------------------------------------------------------------------------

Function InitFlowSubset(Name, FlowModName, FlowName) Result(FlowSubset)
! Initialises a subset of flow module instances.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name           ! Name of subset of flow module
                                             ! instances.
  Character(*), Intent(In) :: FlowModName(:) ! Flow module names.
  Character(*), Intent(In) :: FlowName(:)    ! Flow module instance names.
  ! Function result:
  Type(FlowSubset_) :: FlowSubset ! Subset of flow module instances.
  ! Locals:
  Integer :: i ! Loop index.

  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message('FATAL ERROR in InitFlowSubset', 3)
  End If
  If (Len_Trim(Name) == 0) Then
    Call Message('FATAL ERROR in InitFlowSubset', 3)
  End If

  If (Size(FlowModName) /= Size(FlowName)) Then
    Call Message('UNEXPECTED FATAL ERROR in InitFlowSubset', 4)
  End If
  If (Size(FlowName) > MaxFlows) Then
    Call Message('FATAL ERROR in InitFlowSubset', 3)
  End If

  Do i = 1, Size(FlowName)
    If (Len_Trim(FlowModName(i)) > MaxCharLength) Then
      Call Message('FATAL ERROR in InitFlowSubset', 3)
    End If
    If (Len_Trim(FlowModName(i)) == 0) Then
      Call Message('FATAL ERROR in InitFlowSubset', 3)
    End If
    If (Len_Trim(FlowName(i)) > MaxCharLength) Then
      Call Message('FATAL ERROR in InitFlowSubset', 3)
    End If
    If (Len_Trim(FlowName(i)) == 0) Then
      Call Message('FATAL ERROR in InitFlowSubset', 3)
    End If
  End Do

  FlowSubset%Name                           = Name
  FlowSubset%nFlows                         = Size(FlowName)
  FlowSubset%FlowModNames(1:Size(FlowName)) = FlowModName(1:Size(FlowName))
  FlowSubset%FlowNames   (1:Size(FlowName)) = FlowName   (1:Size(FlowName))

End Function InitFlowSubset

!-------------------------------------------------------------------------------------------------------------

Subroutine AddFlowSubset(FlowSubset, Flows)
! Adds a subset of flow module instances to a collection of flow module instance
! states.

  Implicit None
  ! Argument list:
  Type(FlowSubset_), Intent(In)    :: FlowSubset ! Subset of flow module instances.
  Type(Flows_),      Intent(InOut) :: Flows      ! Collection of flow module instance
                                                 ! states.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Flows%nSubsets
    If (FlowSubset%Name .CIEq. Flows%Subsets(i)%Name) Then
      If (FlowSubsetEq(FlowSubset, Flows%Subsets(i))) Then
        Return
      Else
        Call Message(                                                              &
               'FATAL ERROR in adding Flow Subset "'                            // &
               Trim(FlowSubset%Name)                                            // &
               '": an different flow subset with the same name already exists',    &
               3                                                                   &
             )
      End If
    End If
  End Do

  If (Flows%nSubsets >= MaxFlowSubsets) Then
    Call Message(                                     &
           'FATAL ERROR in adding a Flow Subset: ' // &
           'there are too many Flow Subsets ',        &
           3                                          &
         )
  End If

  Flows%nSubsets                = Flows%nSubsets + 1
  Flows%Subsets(Flows%nSubsets) = FlowSubset

End Subroutine AddFlowSubset

!-------------------------------------------------------------------------------------------------------------

Function FlowSubsetEq(FlowSubset1, FlowSubset2)
! Tests for equality of subsets of flow module instances.

  Implicit None
  ! Argument list:
  Type(FlowSubset_), Intent(In) :: FlowSubset1 !} The two subsets of flow module
  Type(FlowSubset_), Intent(In) :: FlowSubset2 !} instances.
  ! Function result:
  Logical :: FlowSubsetEq ! Indicates if subsets of flow module instances are equal.
  ! Locals:
  Integer :: i ! Loop index.

  FlowSubsetEq = (FlowSubset1%Name   .CIEq. FlowSubset2%Name)  .and. &
                  FlowSubset1%nFlows   ==   FlowSubset2%nFlows
  Do i = 1, Min(FlowSubset1%nFlows, FlowSubset2%nFlows)
    FlowSubsetEq =                                                           &
      FlowSubsetEq                                                     .and. &
      (FlowSubset1%FlowModNames(i) .CIEq. FlowSubset2%FlowModNames(i)) .and. &
      (FlowSubset1%FlowNames(i)    .CIEq. FlowSubset2%FlowNames(i))
  End Do

End Function FlowSubsetEq

!-------------------------------------------------------------------------------------------------------------

Function FindSubsetIndex(Name, Flows)
! Finds the index of a subset of flow module instances.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name  ! Name of subset of flow module instances.
  Type(Flows_), Intent(In) :: Flows ! Collection of flow module instance states.
  ! Function result:
  Integer :: FindSubsetIndex ! Index of the subset of flow module instances.
  ! Locals:
  Integer :: i ! Loop index.

  FindSubsetIndex = 0
  Do i = 1, Flows%nSubsets
    If (Name .CIEq. Flows%Subsets(i)%Name) Then
      FindSubsetIndex = i
      Exit
    End If
  End Do
  If (FindSubsetIndex == 0) Then
    Call Message('FATAL ERROR in FindSubsetIndex: subset not found', 3)
  End If

End Function FindSubsetIndex

!-------------------------------------------------------------------------------------------------------------

Function InitAttrib(Name, OrderName) Result(Attrib)
! Initialises an attribute.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name      ! Name of attribute.
  Character(*), Intent(In) :: OrderName ! Name of ordered subset of flow module instances giving the priority
                                        ! order for the attribute.
  ! Function result:
  Type(FlowAttrib_) :: Attrib ! Initialised attribute.
  ! Locals:
  Integer :: i ! Loop index.

  ! Name.
  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message('FATAL ERROR in InitAttrib', 3)
  End If
  If (Len_Trim(Name) == 0) Then
    Call Message('FATAL ERROR in InitAttrib', 3)
  End If
  Attrib%Name      = Name

  ! iAttribParam.
  Attrib%iAttribParam = 0
  Do i = 1, Size(AttribParamNames)
    If (Name .CIEq. AttribParamNames(i)) Then
      Attrib%iAttribParam = i
      Exit
    End If
  End Do
  If (Attrib%iAttribParam == 0) Then
    Call Message(                                              &
           'FATAL ERROR in reading Flow Attributes: Name "' // &
           Trim(Name)                                       // &
           '" is not a valid attribute name',                  &
           3                                                   &
         )
  End If

  ! OrderName.
  If (Len_Trim(OrderName) > MaxCharLength) Then
    Call Message('FATAL ERROR in InitAttrib', 3)
  End If
  If (Len_Trim(OrderName) == 0) Then
    Call Message('FATAL ERROR in InitAttrib', 3)
  End If
  Attrib%OrderName = OrderName

End Function InitAttrib

!-------------------------------------------------------------------------------------------------------------

Subroutine AddAttrib(Attrib, Flows)
! Adds an attribute to the collection of flow module instance states.

  Implicit None
  ! Argument list:
  Type(FlowAttrib_), Intent(In)    :: Attrib ! Attribute.
  Type(Flows_),      Intent(InOut) :: Flows  ! Collection of flow module instance states.

  If (Flows%nAttribs == MaxFlowAttribs) Then
    Call Message('FATAL ERROR in AddAttrib', 3) ! $$ check for identical attrib names?
  End If
  Flows%nAttribs                = Flows%nAttribs + 1
  Flows%Attribs(Flows%nAttribs) = Attrib

End Subroutine AddAttrib

!-------------------------------------------------------------------------------------------------------------

Function FindAttribIndex(Name, Flows, Error)
! Finds the index of an attribute.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)            :: Name  ! Name of attribute.
  Type(Flows_), Intent(In)            :: Flows ! Collection of flow module instance states.
  Logical,      Intent(Out), Optional :: Error ! Indicates whether attribute found. If
                                               ! not present, the routine stops when
                                               ! the attribute is not found.
  ! Function result:
  Integer :: FindAttribIndex ! Index of attribute.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Flows%nAttribs
    If (Name .CIEq. Flows%Attribs(i)%Name) Then
      FindAttribIndex = i
      If (Present(Error)) Error = .false.
      Return
    End If
  End Do

  If (Present(Error)) Then
    FindAttribIndex = 0
    Error = .true.
  Else
    Call Message(                        &
           'FATAL ERROR: attribute "' // &
           Trim(Name)                 // &
           '" not found',                &
           3                             &
         )
  End If

End Function FindAttribIndex

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpCoordsEtc_Flows(Flows, Coords, Grids)
! Sets up Coords and Grids by adding any extra coords and grids which Flows wants to
! define.

  Implicit None
  ! Argument list:
  Type(Flows_),  Intent(In)    :: Flows  ! Collection of flow module instance states.
  Type(Coords_), Intent(InOut) :: Coords ! Collection of coord systems.
  Type(Grids_),  Intent(InOut) :: Grids  ! Collection of grids.
  ! Locals:
  Integer :: iFlowMod ! Flow module index.
  Integer :: iFlow    ! Flow module instance index.

  ! Add vertical coord systems m agl and Pa to Coords.
  Call AddZCoord(ZCoord_m_agl(), Coords)
  Call AddZCoord(ZCoord_Pa(), Coords)

  Do iFlowMod = 1, Flows%nFlowMods
  Do iFlow    = 1, Flows%nFlows(iFlowMod)

    If (.false.) Then

    Else If (                                   &
      Flows%nNWPFlows /= 0 .and.               &
      iFlowMod == Flows%NWPFlows(1)%C%iFlowMod &
    ) Then

      Call SetUpCoordsEtc_NWPFlow(Flows%NWPFlows(iFlow), Coords, Grids)

    Else If (                                   &
      Flows%nBuildingFlows /= 0 .and.               &
      iFlowMod == Flows%BuildingFlows(1)%C%iFlowMod &
    ) Then

      Call SetUpCoordsEtc_BuildingFlow(Flows%BuildingFlows(iFlow), Coords, Grids)

    Else If (                                   &
      Flows%nLINCOMFlows /= 0 .and.               &
      iFlowMod == Flows%LINCOMFlows(1)%C%iFlowMod &
    ) Then

      Call SetUpCoordsEtc_LINCOMFlow(Flows%LINCOMFlows(iFlow), Coords, Grids)


    ! Check flow module present.
    Else If (.not. FlowPresent(iFlowMod, iFlow, Flows)) Then
      Call Message('UNEXPECTED FATAL ERROR in SetUpCoordsEtc_Flows: flow module instance not found', 4)
    End If

  End Do
  End Do

End Subroutine SetUpCoordsEtc_Flows

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpFlows_MetsCoordsEtc(EtaDefns, Coords, Grids, Mets, Flows)
! Sets up Flows using information from EtaDefns, Coords, Grids and Mets.

  Implicit None
  ! Argument list:
  Type(EtaDefns_), Intent(In)    :: EtaDefns ! Collection of eta definitions.
  Type(Coords_),   Intent(In)    :: Coords   ! Collection of coord systems.
  Type(Grids_),    Intent(In)    :: Grids    ! Collection of grids.
  Type(Mets_),     Intent(In)    :: Mets     ! Collection of met module instance
                                             ! states.
  Type(Flows_),    Intent(InOut) :: Flows    ! Collection of flow module instance
                                             ! states.
  ! Locals:
  Integer :: iFlowMod ! Flow module index.
  Integer :: iFlow    ! Flow module instance index.

  Do iFlowMod = 1, Flows%nFlowMods
  Do iFlow    = 1, Flows%nFlows(iFlowMod)

    If (.false.) Then

    Else If (                                   &
      Flows%nSingleSiteFlows /= 0 .and.               &
      iFlowMod == Flows%SingleSiteFlows(1)%C%iFlowMod &
    ) Then

      Call SetUpSingleSiteFlow_MetsCoordsEtc(     &
             EtaDefns, Coords, Grids, Mets, &
             Flows%SingleSiteFlows(iFlow)         &
           )

    Else If (                                   &
      Flows%nNWPFlows /= 0 .and.               &
      iFlowMod == Flows%NWPFlows(1)%C%iFlowMod &
    ) Then

      Call SetUpNWPFlow_MetsCoordsEtc(     &
             EtaDefns, Coords, Grids, Mets, &
             Flows%NWPFlows(iFlow)         &
           )

    Else If (                                   &
      Flows%nRadarFlows /= 0 .and.               &
      iFlowMod == Flows%RadarFlows(1)%C%iFlowMod &
    ) Then

      Call SetUpRadarFlow_MetsCoordsEtc(     &
             EtaDefns, Coords, Grids, Mets, &
             Flows%RadarFlows(iFlow)         &
           )

    Else If (                                   &
      Flows%nBuildingFlows /= 0 .and.               &
      iFlowMod == Flows%BuildingFlows(1)%C%iFlowMod &
    ) Then

      Call SetUpBuildingFlow_MetsCoordsEtc(     &
             EtaDefns, Coords, Grids, Mets, &
             Flows%BuildingFlows(iFlow)         &
           )

    Else If (                                   &
      Flows%nLINCOMFlows /= 0 .and.               &
      iFlowMod == Flows%LINCOMFlows(1)%C%iFlowMod &
    ) Then

      Call SetUpLINCOMFlow_MetsCoordsEtc(     &
             EtaDefns, Coords, Grids, Mets, &
             Flows%LINCOMFlows(iFlow)         &
           )


    ! Check flow module present.
    Else If (.not. FlowPresent(iFlowMod, iFlow, Flows)) Then
      Call Message('UNEXPECTED FATAL ERROR in SetUpFlows_MetsCoordsEtc: flow module instance not found', 4)
    End If

  End Do
  End Do

  ! AnyUpdateOnDemand.
  Flows%AnyUpdateOnDemand = .false.

  Do iFlowMod = 1, Flows%nFlowMods
  Do iFlow    = 1, Flows%nFlows(iFlowMod)

    If (Flows%C(iFlowMod, iFlow)%P%UpdateOnDemand) Flows%AnyUpdateOnDemand = .true.

  End Do
  End Do

End Subroutine SetUpFlows_MetsCoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpFlowOrders_iFlows(Flows)
! Sets up indices in Flows%Orders for referring to flow module instances. Also sets up
! Flows%Orders%Included.

  Implicit None
  ! Argument list:
  Type(Flows_), Intent(InOut) :: Flows ! Collection of flow module instance states.
  ! Locals:
  Integer :: iOrder ! Index of ordered subset of flow module instances.
  Integer :: iFlow  ! Index of the flow within the ordered subset of flow module instances.

  Do iOrder = 1, Flows%nOrders
    Flows%Orders(iOrder)%Included = .false.
    Do iFlow = 1, Flows%Orders(iOrder)%nFlows
      ! Get flow index.
      Call FindFlowIndex(                              &
             Flows%Orders(iOrder)%FlowModNames(iFlow), &
             Flows%Orders(iOrder)%FlowNames(iFlow),    &
             Flows,                                    &
             Flows%Orders(iOrder)%iFlowMods(iFlow),    &
             Flows%Orders(iOrder)%iFlows(iFlow)        &
           )
      ! Record the flow as included.
      Flows%Orders(iOrder)%Included(           &
        Flows%Orders(iOrder)%iFlowMods(iFlow), &
        Flows%Orders(iOrder)%iFlows(iFlow)     &
      ) = .true.
    End Do
  End Do

End Subroutine SetUpFlowOrders_iFlows

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpFlowSubsets_iFlows(Flows)
! Sets up indices in Flows%Subsets for referring to flow module instances. Also sets
! up Flows%Subsets%Included.

  Implicit None
  ! Argument list:
  Type(Flows_), Intent(InOut) :: Flows ! Collection of flow module instance states.
  ! Locals:
  Integer :: iSubset ! Index of subset of flow module instances.
  Integer :: iFlow   ! Index of the flow within the subset of flow module instances.

  Do iSubset = 1, Flows%nSubsets
    Flows%Subsets(iSubset)%Included = .false.
    Do iFlow = 1, Flows%Subsets(iSubset)%nFlows
      ! Get flow index.
      Call FindFlowIndex(                                &
             Flows%Subsets(iSubset)%FlowModNames(iFlow), &
             Flows%Subsets(iSubset)%FlowNames(iFlow),    &
             Flows,                                      &
             Flows%Subsets(iSubset)%iFlowMods(iFlow),    &
             Flows%Subsets(iSubset)%iFlows(iFlow)        &
           )
      ! Record the flow as included.
      Flows%Subsets(iSubset)%Included(           &
        Flows%Subsets(iSubset)%iFlowMods(iFlow), &
        Flows%Subsets(iSubset)%iFlows(iFlow)     &
      ) = .true.
    End Do
  End Do

End Subroutine SetUpFlowSubsets_iFlows

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpFlows_iFlowAttribs(Flows)
! Sets up indices in Flows for referring to flow attributes. Also checks that the flow
! module instances have the required attributes.

  Implicit None
  ! Argument list:
  Type(Flows_), Intent(InOut) :: Flows ! Collection of flow module instance states.
  ! Locals:
  Integer :: iAttribParam !} Indices of attribute in the parameters in the common flow module and in
  Integer :: iAttrib      !} Flows%Attrib.
  Integer :: iOrder       ! Index of ordered subset of flow module instances.
  Integer :: iFlowInOrder ! Index of the flow within the ordered subset of flow module instances.
  Integer :: iFlowMod     ! Flow module index.
  Integer :: iFlow        ! Flow module instance index.
  Logical :: Flag         ! Indicates whether the flow module instance has the attribute.
  Integer :: i            ! Loop Index.
  Logical :: Error        ! Indicates attribute not found by FindAttribIndex.

  Do iAttrib = 1, Flows%nAttribs

    iOrder = FindOrderIndex(Flows%Attribs(iAttrib)%OrderName,  Flows)
    Flows%Attribs(iAttrib)%iOrder = iOrder

    Do iFlowInOrder = 1, Flows%Orders(iOrder)%nFlows

      ! Get flow index.
      Call FindFlowIndex(                                     &
             Flows%Orders(iOrder)%FlowModNames(iFlowInOrder), &
             Flows%Orders(iOrder)%FlowNames(iFlowInOrder),    &
             Flows,                                           &
             iFlowMod,                                        &
             iFlow                                            &
           )

      ! Check flow has attribute.
      Flag = .false.
      Do i = 1, Flows%C(iFlowMod, iFlow)%P%nAttribs
        If (Flows%C(iFlowMod, iFlow)%P%AttribNames(i) .CIEq. Flows%Attribs(iAttrib)%Name) Then
          Flag = .true.
          Exit
        End If
      End Do
      If (.not.Flag) Then
        Call Message('FATAL ERROR in SetUpFlows_iFlowAttribs', 3)
      End If

    End Do

  End Do

  ! Set attribute indices.

  Do iAttribParam = 1, Size(AttribParamNames)
    Flows%iAttribs(iAttribParam) = FindAttribIndex(AttribParamNames(iAttribParam), Flows, Error)
    If (Error) Flows%iAttribs(iAttribParam) = 0
  End Do

  If (Flows%iAttribs(A_Update) == 0) Then
    Call Message('FATAL ERROR in SetUpFlows_iFlowAttribs: Update attribute not defined', 3)
  End If

  If (Flows%iAttribs(A_Convert) == 0) Then
    Call Message('FATAL ERROR in SetUpFlows_iFlowAttribs: Convert attribute not defined', 3)
  End If

  If (Flows%iAttribs(A_Flow) == 0) Then
    Call Message('FATAL ERROR in SetUpFlows_iFlowAttribs: Flow attribute not defined', 3)
  End If

End Subroutine SetUpFlows_iFlowAttribs

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpFlows_iFlowSubsets(Flows)
! Sets up indices in the flow module instances for referring to flow subsets.

  Implicit None
  ! Argument list:
  Type(Flows_), Intent(InOut) :: Flows ! Collection of flow module instance states.
  ! Locals:
  Integer :: iFlowMod ! Flow module index.
  Integer :: iFlow    ! Flow module instance index.
  Integer :: i        ! Loop index.

  Do iFlowMod = 1, Flows%nFlowMods
  Do iFlow    = 1, Flows%nFlows(iFlowMod)
    Do i = 1, Flows%C(iFlowMod, iFlow)%P%nUpdateSubsets
      Flows%C(iFlowMod, iFlow)%P%iUpdateSubsets(i) = FindSubsetIndex(                                   &
                                                       Flows%C(iFlowMod, iFlow)%P%UpdateSubsetNames(i), &
                                                       Flows                                            &
                                                     )
    End Do
  End Do
  End Do

End Subroutine SetUpFlows_iFlowSubsets

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpFlows_iMets(Mets, Flows)
! Sets up indices in the flow module instances for referring to met module instances.

  Implicit None
  ! Argument list:
  Type(Mets_),  Intent(In)    :: Mets  ! Collection of met module instance states.
  Type(Flows_), Intent(InOut) :: Flows ! Collection of flow module instance states.
  ! Locals:
  Integer :: iFlowMod ! Flow module index.
  Integer :: iFlow    ! Flow module instance index.
  Integer :: i        ! Loop Index.

  Do iFlowMod = 1, Flows%nFlowMods
  Do iFlow    = 1, Flows%nFlows(iFlowMod)
    Do i = 1, Flows%C(iFlowMod, iFlow)%P%nMets
      Call FindMetIndex(                           &
        Flows%C(iFlowMod, iFlow)%P%MetModNames(i), &
        Flows%C(iFlowMod, iFlow)%P%MetNames(i),    &
        Mets,                                      &
        Flows%C(iFlowMod, iFlow)%P%iMetMods(i),    &
        Flows%C(iFlowMod, iFlow)%P%iMets(i)        &
      )
    End Do
  End Do
  End Do

End Subroutine SetUpFlows_iMets

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpFlows_iCoordsEtc(Coords, Grids, Domains, Flows)
! Sets up indices in the flow module instances for referring coord systems and grids.

  Implicit None
  ! Argument list:
  Type(Coords_),  Intent(In)            :: Coords  ! Collection of coord systems.
  Type(Grids_),   Intent(In)            :: Grids   ! Collection of grids.
  Type(Domains_), Intent(In)            :: Domains ! Collection of domains.
  Type(Flows_),   Intent(InOut), Target :: Flows   ! Collection of flow module instance states.
  ! Locals:
  Type(CommonFlow_), Pointer :: CommonFlow ! Abreviation for the part of the flow state common to all flow
                                           ! modules.
  Integer                    :: iFlowMod   ! Flow module index.
  Integer                    :: iFlow      ! Flow module instance index.
  Integer                    :: i          ! Loop index.

  Flows%iZCoordMagl = FindZCoordIndex('m agl', Coords)
  Flows%iZCoordPa   = FindZCoordIndex('Pa',    Coords)

  Do iFlowMod = 1, Flows%nFlowMods
  Do iFlow    = 1, Flows%nFlows(iFlowMod)

    CommonFlow => Flows%C(iFlowMod, iFlow)%P

    Do i = 1, CommonFlow%nHCoords
      CommonFlow%iHCoords(i) = FindHCoordIndex(CommonFlow%HCoordNames(i), Coords)
    End Do

    Do i = 1, CommonFlow%nZCoords
      CommonFlow%iZCoords(i) = FindZCoordIndex(CommonFlow%ZCoordNames(i), Coords)
    End Do

    Do i = 1, CommonFlow%nHGrids
      CommonFlow%iHGrids(i) = FindHGridIndex(CommonFlow%HGridNames(i), Grids)
    End Do

    Do i = 1, CommonFlow%nZGrids
      CommonFlow%iZGrids(i) = FindZGridIndex(CommonFlow%ZGridNames(i), Grids)
    End Do

    CommonFlow%iDomain = FindDomainIndex(CommonFlow%DomainName, Domains)

    If (CommonFlow%FixedMet) Then
      If (                                                                          &
        .not.IsInfPast(StartTimeOfDomain(Domains%Domains(CommonFlow%iDomain))) .or. &
        .not.IsInfFuture(EndTimeOfDomain(Domains%Domains(CommonFlow%iDomain)))      &
      ) Then
        Call Message(                                                                                    &
               'FATAL ERROR: The domain "'                                                            // &
               Trim(CommonFlow%DomainName)                                                            // &
               '" of flow module instance "'                                                          // &
               Trim(CommonFlow%FlowName)                                                              // &
               '" is not unbounded in time but the met is fixed. This combination is not permitted.',    &
               3                                                                                         &
             )
      End If
    End If

  End Do
  End Do

End Subroutine SetUpFlows_iCoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine CheckFlows_FlowUpdateSubsets(Flows)
! Checks that the update subsets used within the flow module instances only include
! flow module instances which are higher up the update priority order.

  Implicit None
  ! Argument list:
  Type(Flows_), Intent(InOut) :: Flows ! Collection of flow module instance states.
  ! Locals:
  Integer :: iFlowMod  !} Indices of the flow module instance whose update subsets are
  Integer :: iFlow     !} being considered.
  Integer :: iOrder    ! Index of priority order for the update attribute.
  Integer :: iSubset   ! Index of update subset being considered.
  Integer :: jFlowMod  !} Indices of the flow module instance whose membership of the
  Integer :: jFlow     !} update subset is being considered.
  Logical :: GotToFlow ! Indicates that, in looping through the flow module instances,
                       ! the flow module instance whose update subsets are being
                       ! considered has been reached.
  Integer :: i         ! Loop index.
  Integer :: j         ! Loop index.

  iOrder = Flows%Attribs(Flows%iAttribs(A_Update))%iOrder

  Do iFlowMod = 1, Flows%nFlowMods
  Do iFlow    = 1, Flows%nFlows(iFlowMod)

    Do i = 1, Flows%C(iFlowMod, iFlow)%P%nUpdateSubsets
      iSubset = Flows%C(iFlowMod, iFlow)%P%iUpdateSubsets(i)
      GotToFlow = .false.

      Do j = 1, Flows%Orders(iOrder)%nFlows
        jFlowMod = Flows%Orders(iOrder)%iFlowMods(j)
        jFlow    = Flows%Orders(iOrder)%iFlows(j)
        If (jFlowMod == iFlowMod .and. jFlow == iFlow) Then
          GotToFlow = .true.
        End If
        If (GotToFlow) Then
          If (Flows%Subsets(iSubset)%Included(jFlowMod, jFlow)) Then
            Call Message('FATAL ERROR in CheckFlows_FlowUpdateSubsets', 3)
          End If
        End If
      End Do

    End Do

  End Do
  End Do

End Subroutine CheckFlows_FlowUpdateSubsets

!-------------------------------------------------------------------------------------------------------------

Subroutine FlowList(Flows, nFlows, FlowNames)
! Returns a list of the flow module instances.

  Implicit None
  ! Argument list:
  Type(Flows_),             Intent(In)  :: Flows
  Integer,                  Intent(Out) :: nFlows
  Character(MaxCharLength), Intent(Out) :: FlowNames(MaxFlows)
  ! Flows     :: Collection of flow module instance states.
  ! nFlows    :: Number of the flow module instance states.
  ! FlowNames :: Names of the flow module instance states.
  ! Locals:
  Integer :: iFlowMod ! Flow module index.
  Integer :: iFlow    ! Flow module instance index.

  nFlows = 0
  Do iFlowMod = 1, Flows%nFlowMods
  Do iFlow    = 1, Flows%nFlows(iFlowMod)
    nFlows = nFlows + 1
    FlowNames(nFlows) = Trim(Flows%C(iFlowMod, iFlow)%P%FlowModName) // &
                        '.'                                          // &
                        Trim(Flows%C(iFlowMod, iFlow)%P%FlowName)
  End Do
  End Do


End Subroutine FlowList

!-------------------------------------------------------------------------------------------------------------

Function FlowPresent(iFlowMod, iFlow, Flows)
! Checks whether a flow module instance with given indices is present.

  Implicit None
  ! Argument list:
  Integer,      Intent(In) :: iFlowMod ! Flow module index.
  Integer,      Intent(In) :: iFlow    ! Flow module instance index.
  Type(Flows_), Intent(In) :: Flows    ! Collection of flow module instance states.
  ! Function result:
  Logical :: FlowPresent ! Indicates whether the flow module instance is present.

  FlowPresent = .false.

  If (.false.) Then

  Else If (                                   &
    Flows%nPrototypeFlows /= 0 .and.               &
    iFlowMod == Flows%PrototypeFlows(1)%C%iFlowMod &
  ) Then

    FlowPresent = .true.

  Else If (                                   &
    Flows%nSingleSiteFlows /= 0 .and.               &
    iFlowMod == Flows%SingleSiteFlows(1)%C%iFlowMod &
  ) Then

    FlowPresent = .true.

  Else If (                                   &
    Flows%nNWPFlows /= 0 .and.               &
    iFlowMod == Flows%NWPFlows(1)%C%iFlowMod &
  ) Then

    FlowPresent = .true.

  Else If (                                   &
    Flows%nBuildingFlows /= 0 .and.               &
    iFlowMod == Flows%BuildingFlows(1)%C%iFlowMod &
  ) Then

    FlowPresent = .true.

  Else If (                                   &
    Flows%nRadarFlows /= 0 .and.               &
    iFlowMod == Flows%RadarFlows(1)%C%iFlowMod &
  ) Then

    FlowPresent = .true.

  Else If (                                   &
    Flows%nLINCOMFlows /= 0 .and.               &
    iFlowMod == Flows%LINCOMFlows(1)%C%iFlowMod &
  ) Then

    FlowPresent = .true.


  End If

End Function FlowPresent

!-------------------------------------------------------------------------------------------------------------

Function MetsFlowsOverallTValid(Flows) Result(OverallTValid)
! Returns the earliest time that the validity of any of the met or flow module instances might change,
! assuming all the instances which have been prepared for update-on-demand are updated now. The value is that
! determined at the time UpdateMetsFlows was last called (the actual time may be later).

  Implicit None
  ! Argument list:
  Type(Flows_), Intent(In) :: Flows ! Collection of flow module instance states.
  ! Function result:
  Type(Time_) :: OverallTValid ! Earliest time that the validity of any of the met or flow module instances
                               ! might change, assuming all the instances which have been prepared for
                               ! update-on-demand are updated now. The value is that determined at the time
                               ! UpdateMetsFlows was last called (the actual time may be later).

  OverallTValid = Flows%OverallTValid

End Function MetsFlowsOverallTValid

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateMetsFlows(          &
             Coords, Grids, Domains, &
             iCase, iMetCase,        &
             MetTime,                &
             Mets, Flows,            &
             Units                   &
           )
! Updates the state of the met and flow module instances.

  Implicit None
  ! Argument list:
  Type(Coords_),  Intent(In)    :: Coords        ! Collection of coord systems.
  Type(Grids_),   Intent(In)    :: Grids         ! Collection of grids.
  Type(Domains_), Intent(In)    :: Domains       ! Collection of domains.
  Integer,        Intent(In)    :: iCase         ! Number of case.
  Integer,        Intent(In)    :: iMetCase      ! Number of the met realisation in the met ensemble.
  Type(Time_),    Intent(In)    :: MetTime       ! Time for which the met and flow module instances are to be
                                                 ! updated (this must be the current time unless we have a
                                                 ! fixed met case, in which case it must be the time of the
                                                 ! fixed met).
  Type(Mets_),    Intent(InOut) :: Mets          ! Collection of met module instance states.
  Type(Flows_),   Intent(InOut) :: Flows         ! Collection of flow module instance states.
  Type(Units_),   Intent(InOut) :: Units         ! Collection of information on input/output unit numbers.

  ! Record iCase, iMetCase and MetTime.
  Flows%iCase    = iCase
  Flows%iMetCase = iMetCase
  Flows%MetTime  = Time2ShortTime(MetTime)
  ! $$ for fixed met should check MetTime the same each time within one case
  ! $$ should check iCase and iMetCase the same each time within one case?

  ! Message to announce updating of met and flow modules.
  If (Flows%AnyUpdateOnDemand) Then
    If (Flows%C(1, 1)%P%FixedMet) Then
      Call Message(' ')
      Call Message(                                                             &
             'Preparing to update met and flow modules at (fixed met) time ' // &
             Trim(Time2Char(MetTime, .false., 0, .true.))                       &
           )
    Else
      Call Message(' ')
      Call Message(                                                 &
             'Preparing to update met and flow modules at time ' // &
             Trim(Time2Char(MetTime, .false., 0, .true.))           &
           )
    End If
  Else
    If (Flows%C(1, 1)%P%FixedMet) Then
      Call Message(' ')
      Call Message(                                                  &
             'Updating met and flow modules at (fixed met) time ' // &
             Trim(Time2Char(MetTime, .false., 0, .true.))            &
           )
    Else
      Call Message(' ')
      Call Message(                                          &
             'Updating met and flow modules at time '     // &
             Trim(Time2Char(MetTime, .false., 0, .true.))    &
           )
    End If
  End If

  ! Initialise OverallTValid.
  Flows%OverallTValid = InfFutureTime()

  ! Update met module instances.
  Call UpdateMets(                            &
         Coords, Grids,                       &
         Flows%SameResultsWithUpdateOnDemand, &
         iCase, iMetCase,                     &
         MetTime,                             &
         Flows%OverallTValid,                 &
         Mets,                                &
         Units                                &
       )

  ! Update flow module instances.
  Call UpdateFlows(                           &
         Coords, Grids, Domains,              &
         Flows%SameResultsWithUpdateOnDemand, &
         iCase, iMetCase,                     &
         MetTime,                             &
         Flows%OverallTValid,                 &
         Mets, Flows,                         &
         Units                                &
       )

End Subroutine UpdateMetsFlows

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateFlows(                     &
             Coords, Grids, Domains,        &
             SameResultsWithUpdateOnDemand, &
             iCase, iMetCase,               &
             MetTime,                       &
             OverallTValid,                 &
             Mets, Flows,                   &
             Units                          &
           )
! Updates the state of the flow module instances.

! Note this routine might also update met module instances if update-on-demand is being used.

  Implicit None
  ! Argument list:
  Type(Coords_),  Intent(In)    :: Coords
  Type(Grids_),   Intent(In)    :: Grids
  Type(Domains_), Intent(In)    :: Domains
  Logical,        Intent(In)    :: SameResultsWithUpdateOnDemand
  Integer,        Intent(In)    :: iCase
  Integer,        Intent(In)    :: iMetCase
  Type(Time_),    Intent(In)    :: MetTime
  Type(Time_),    Intent(InOut) :: OverallTValid
  Type(Mets_),    Intent(InOut) :: Mets
  Type(Flows_),   Intent(InOut) :: Flows
  Type(Units_),   Intent(InOut) :: Units
  ! Coords                        :: Collection of coord systems.
  ! Grids                         :: Collection of grids.
  ! Domains                       :: Collection of domains.
  ! SameResultsWithUpdateOnDemand :: Indicates results should be made the same whether the met and flow module
  !                                  instances use update-on-demand or update-at-once (so it should really be
  !                                  called SameResultsWithOrWithoutUpdateOnDemand).
  ! iCase                         :: Number of case.
  ! iMetCase                      :: Number of the met realisation in the met ensemble.
  ! MetTime                       :: Time for which the flow module instances are to be updated (this must be
  !                                  the current time unless we have a fixed met case, in which case it must
  !                                  be the time of the fixed met).
  ! OverallTValid                 :: Earliest time that the validity of any of the met or flow module
  !                                  instances might change, assuming all the instances which have been
  !                                  prepared for update-on-demand are updated now. The value is that
  !                                  determined at the end of this routine (the actual time may be later).
  ! Mets                          :: Collection of met module instance states.
  ! Flows                         :: Collection of flow module instance states.
  ! Units                         :: Collection of information on input/output unit numbers.
  ! Locals:
  Type(CommonFlow_), Pointer :: CommonFlow
  Integer                    :: iFlowMod
  Integer                    :: iFlow
  Integer                    :: iOrder
  Logical                    :: ValidFlow
  Type(Time_)                :: TValid
  Logical                    :: UpdateNow
  Integer                    :: i
  ! CommonFlow :: Abbreviation for the common part of the flow module instance.
  ! iFlowMod   :: Flow module index.
  ! iFlow      :: Flow module instance index.
  ! iOrder     :: Index of ordered subset of flow module instances.
  ! ValidFlow  :: Indicates there is a valid flow module with the flow attribute.
  ! TValid     :: Earliest time that the validity (overall or for any single attribute) of the flow module
  !               instance might change, except perhaps for a change caused by a change in the validity of the
  !               met and flow module instances acting as data sources, assuming the flow module instance is
  !               updated now. The value is that determined at the end of PrepareForUpdateFlow (the actual
  !               time may be later).
  ! UpdateNow  :: Indicates the flow module instance must be updated now (even if update-on-demand is
  !               specified). If set, the value of TValid is not reliable.
  ! i          :: Loop index.

  ! Initialise ValidFlow.
  ValidFlow = .false.

  ! Index of ordered subset of flow module instances.
  iOrder = Flows%Attribs(Flows%iAttribs(A_Update))%iOrder

  ! Loop over flow module instances in the right order.
  Do i = 1, Flows%Orders(iOrder)%nFlows

    ! Flow module instance indices.
    iFlowMod = Flows%Orders(iOrder)%iFlowMods(i)
    iFlow    = Flows%Orders(iOrder)%iFlows(i)

    CommonFlow => Flows%C(iFlowMod, iFlow)%P

    ! Check MetTime is the right type of time and is finite.
    If (i == 1) Then
      If (IsTimeInterval(MetTime)) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlows', 4)
      End If
      If (IsInfFuture(MetTime) .or. IsInfPast(MetTime)) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlows', 4)
      End If
    End If

    ! Is the validity of the flow module instance improvable?
    If (.not. CommonFlow%ValidityUnimprovable .or. MetTime >= CommonFlow%TValid) Then

        ! $$ Note that, for the same results whether the met and flow module instances use update-on-demand or
        ! update-at-once, two things need to be done. Firstly we need to ensure that (i) the sync times, and
        ! (ii) the times at which flow data is obtained to update other flow modules, are not altered by the
        ! need to break to call UpdateMetsFlows or by the absence of such a need. This requires OverallTValid
        ! to use the local variable TValid for update-at-once as it would for update-on-demand (of course
        ! often local
        ! TValid = CommonFlow%TValid). This is done here. Note that for this to work the UpdateNow flag
        ! returned by PrepareForUpdate????Flow must not depend on whether update-on-demand or update-at-once
        ! is used. Secondly we need to ensure that met and flow values after (the final value at the end of
        ! UpdateMets and UpdateFlows of) time OverallTValid are not altered by whether or not a met or flow
        ! module instance is updated before time OverallTValid. This is considered elsewhere.
        !
        ! $$ UpdateNow flag mustn't depend on MetOnDemand
        !
        ! $$ The second issue is trivial for many met/flow modules but hard where its not trivial
        ! and hasn't yet been
        ! done yet. Suggested approach:
        ! Ensure modules which depend on their history always use UpdateNow.
        ! Ensure modules which depend on the precise update time (i.e. value at 0100 depends on
        ! whether updated at 1100
        !   and valid till 0200 or updated at 1200 and valid till 0200) either
        !   update (or prepare to update) on every call to UpdateMetsFlows by limiting TValid
        !   or go back to latest update of any source data and update forwards from there
        !   (or from latest of any
        !   subsequent members of a sequence of quantised times)
        !   or always use UpdateNow.

      ! Prepare to update flow module instance.
      Call PrepareForUpdateFlow(     &
             Coords, Grids, Domains, &
             Mets,                   &
             iCase,                  &
             iFlowMod, iFlow,        &
             MetTime,                &
             TValid, UpdateNow,      &
             Flows,                  &
             Units                   &
           )

      If (UpdateNow) Then

        ! Try to update flow module instance.
        Call UpdateFlow(               &
               Coords, Grids, Domains, &
               iCase, iMetCase,        &
               iFlowMod, iFlow,        &
               MetTime,                &
               Mets, Flows,            &
               Units                   &
             )
        ! Update OverallTValid.
        OverallTValid = TMin(OverallTValid, CommonFlow%TValid)

      Else If (.not. CommonFlow%UpdateOnDemand) Then

        ! Try to update flow module instance.
        Call UpdateFlow(               &
               Coords, Grids, Domains, &
               iCase, iMetCase,        &
               iFlowMod, iFlow,        &
               MetTime,                &
               Mets, Flows,            &
               Units                   &
             )
        ! Update OverallTValid. For same results whether update-on-demand or update-at-once is used, update
        ! OverallTValid using the local variable TValid for update-at-once, as would be used for
        ! update-on-demand.
        If (.not. SameResultsWithUpdateOnDemand) Then
          OverallTValid = TMin(OverallTValid, CommonFlow%TValid)
        Else
          OverallTValid = TMin(OverallTValid, TValid)
        End If

      Else

        ! Update OverallTValid.
        OverallTValid = TMin(OverallTValid, TValid)

      End If

    Else

      ! Update OverallTValid.
      OverallTValid = TMin(OverallTValid, CommonFlow%TValid)

    End If

    ! Update ValidFlow.
    ValidFlow = ValidFlow .or. CommonFlow%ValidAttribParams(A_Flow)

  End Do

  ! Check for presence of flow module which is valid for the flow attribute.
  ! $$ extend to other attribs?
  If (.not. Flows%AnyUpdateOnDemand) Then
    If (.not. ValidFlow) Then
      Call ControlledMessage(                                                       &
             'There are no flow modules which are valid for the "flow" attribute.', &
             MessageControls    = GlobalMessageControls,                            &
             MessageControlName = 'No Flow',                                        &
             ErrorCode          = 3                                                 &
           )
    End If
  End If

End Subroutine UpdateFlows

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForUpdateFlow(     &
             Coords, Grids, Domains, &
             Mets,                   &
             iCase,                  &
             iFlowMod, iFlow,        &
             MetTime,                &
             TValid, UpdateNow,      &
             Flows,                  &
             Units                   &
           )
! Prepares for updating the state of a flow module instance.

  Implicit None
  ! Argument list:
  Type(Coords_),  Intent(In)    :: Coords    ! Collection of coord systems.
  Type(Grids_),   Intent(In)    :: Grids     ! Collection of grids.
  Type(Domains_), Intent(In)    :: Domains   ! Collection of domains.
  Type(Mets_),    Intent(In)    :: Mets      ! Collection of met module instance states.
  Integer,        Intent(In)    :: iCase     ! Number of case.
  Integer,        Intent(In)    :: iFlowMod  ! Flow module index.
  Integer,        Intent(In)    :: iFlow     ! Flow module instance index.
  Type(Time_),    Intent(In)    :: MetTime   ! Time for which the flow module instance is to be updated (this
                                             ! must be the current time unless we have a fixed met case, in
                                             ! which case it must be the time of the fixed met).
  Type(Time_),    Intent(Out)   :: TValid    ! Earliest time that the validity (overall or for any single
                                             ! attribute) of the flow module instance might change, except
                                             ! perhaps for a change caused by a change in the validity of the
                                             ! met and flow module instances acting as data sources, assuming
                                             ! the flow module instance is updated now. The value is that
                                             ! determined at the end of this routine (the actual time may be
                                             ! later).
  Logical,        Intent(Out)   :: UpdateNow ! Indicates the flow module instance must be updated now (even if
                                             ! update-on-demand is specified). If set, TValid need not be set
                                             ! to any particular time.
  Type(Flows_),   Intent(InOut) :: Flows     ! Collection of flow module instance states.
  Type(Units_),   Intent(InOut) :: Units     ! Collection of information on input/output unit numbers.
  ! Locals:
  Type(CommonFlow_), Pointer :: CommonFlow
  ! CommonFlow :: Abbreviation for the common part of the flow module instance.

  CommonFlow => Flows%C(iFlowMod, iFlow)%P

  If (.false.) Then

  Else If (Flows%nPrototypeFlows /= 0 .and. iFlowMod == Flows%PrototypeFlows(1)%C%iFlowMod) Then

    Call PrepareForUpdatePrototypeFlow( &
           Coords, Grids, Domains, &
           Mets,                   &
           iCase,                  &
           MetTime,                &
           TValid, UpdateNow,      &
           Flows%PrototypeFlows(iFlow), &
           Units                   &
         )

  Else If (Flows%nSingleSiteFlows /= 0 .and. iFlowMod == Flows%SingleSiteFlows(1)%C%iFlowMod) Then

    Call PrepareForUpdateSingleSiteFlow( &
           Coords, Grids, Domains, &
           Mets,                   &
           iCase,                  &
           MetTime,                &
           TValid, UpdateNow,      &
           Flows%SingleSiteFlows(iFlow), &
           Units                   &
         )

  Else If (Flows%nNWPFlows /= 0 .and. iFlowMod == Flows%NWPFlows(1)%C%iFlowMod) Then

    Call PrepareForUpdateNWPFlow( &
           Coords, Grids, Domains, &
           Mets,                   &
           iCase,                  &
           MetTime,                &
           TValid, UpdateNow,      &
           Flows%NWPFlows(iFlow), &
           Units                   &
         )

  Else If (Flows%nBuildingFlows /= 0 .and. iFlowMod == Flows%BuildingFlows(1)%C%iFlowMod) Then

    Call PrepareForUpdateBuildingFlow( &
           Coords, Grids, Domains, &
           Mets,                   &
           iCase,                  &
           MetTime,                &
           TValid, UpdateNow,      &
           Flows%BuildingFlows(iFlow), &
           Units                   &
         )

  Else If (Flows%nRadarFlows /= 0 .and. iFlowMod == Flows%RadarFlows(1)%C%iFlowMod) Then

    Call PrepareForUpdateRadarFlow( &
           Coords, Grids, Domains, &
           Mets,                   &
           iCase,                  &
           MetTime,                &
           TValid, UpdateNow,      &
           Flows%RadarFlows(iFlow), &
           Units                   &
         )

  Else If (Flows%nLINCOMFlows /= 0 .and. iFlowMod == Flows%LINCOMFlows(1)%C%iFlowMod) Then

    Call PrepareForUpdateLINCOMFlow( &
           Coords, Grids, Domains, &
           Mets,                   &
           iCase,                  &
           MetTime,                &
           TValid, UpdateNow,      &
           Flows%LINCOMFlows(iFlow), &
           Units                   &
         )


  Else
    Call Message('UNEXPECTED FATAL ERROR in PrepareForUpdateFlow: flow module instance not found', 4)
  End If

  ! Check TValid.
  If (.not. UpdateNow) Then

    ! Check TValid is the right type of time.
    If (IsTimeInterval(TValid)) Then
      Call Message('UNEXPECTED FATAL ERROR in PrepareForUpdateFlow', 4)
    End If

    ! Check that TValid is in the future, and, for fixed met, is infinite.
    If (TValid <= MetTime) Then
      Call Message('UNEXPECTED FATAL ERROR in PrepareForUpdateFlow', 4)
    End If
    If (CommonFlow%FixedMet) Then
      If (.not. IsInfFuture(TValid)) Then
        Call Message('UNEXPECTED FATAL ERROR in PrepareForUpdateFlow', 4)
      End If
    End If

  End If

  ! Update DueForUpdate.
  CommonFlow%DueForUpdate = .true.

End Subroutine PrepareForUpdateFlow

!-------------------------------------------------------------------------------------------------------------

Recursive Subroutine UpdateFlow(               &
                       Coords, Grids, Domains, &
                       iCase, iMetCase,        &
                       iFlowMod, iFlow,        &
                       MetTime,                &
                       Mets, Flows,            &
                       Units                   &
                     )
! Updates the state of a flow module instance.

! Note this routine might also update met module instances if update-on-demand is being used.

  Implicit None
  ! Argument list:
  Type(Coords_),  Intent(In)    :: Coords   ! Collection of coord systems.
  Type(Grids_),   Intent(In)    :: Grids    ! Collection of grids.
  Type(Domains_), Intent(In)    :: Domains  ! Collection of domains.
  Integer,        Intent(In)    :: iCase    ! Number of case.
  Integer,        Intent(In)    :: iMetCase ! Number of the met realisation in the met ensemble.
  Integer,        Intent(In)    :: iFlowMod ! Flow module index.
  Integer,        Intent(In)    :: iFlow    ! Flow module instance index.
  Type(Time_),    Intent(In)    :: MetTime  ! Time for which the flow module instance is to be updated (this
                                            ! must be the current time unless we have a fixed met case, in
                                            ! which case it must be the time of the fixed met).
  Type(Mets_),    Intent(InOut) :: Mets     ! Collection of met module instance states.
  Type(Flows_),   Intent(InOut) :: Flows    ! Collection of flow module instance states.
  Type(Units_),   Intent(InOut) :: Units    ! Collection of information on input/output unit numbers.
  ! Locals:
  Type(CommonFlow_), Pointer :: CommonFlow
  Logical                    :: WantFlowField
  Logical                    :: WantCloudField
  Logical                    :: WantRainField
  Integer                    :: iFlowUpdateSubset
  Integer                    :: iCloudUpdateSubset
  Integer                    :: iRainUpdateSubset
  Type(FlowField_)           :: FlowField
  Type(CloudField_)          :: CloudField
  Type(RainField_)           :: RainField
  Logical                    :: Error
  Integer                    :: i
  ! CommonFlow         :: Abbreviation for the common part of the flow module
  !                       instance.
  ! WantFlowField      :} Indicates the flow module instance wants information for
  ! WantCloudField     :} various attributes from other flow module instances.
  ! WantRainField      :}
  ! iFlowUpdateSubset  :] Index of update subsets for various attributes in the array
  ! iCloudUpdateSubset :] of update subset information in the common part of the flow
  ! iRainUpdateSubset  :] module instance.
  ! FlowField          :} Field of information for various attributes used (i) to
  ! CloudField         :} specify what information the flow module instance wants from
  ! RainField          :} the other flow module instances, and (ii) to transfer this
  !                       information from the other flow module instances to the flow
  !                       module instance.
  ! Error              :: Indicates no suitable flow module instance was found for
  !                       providing the information wanted from the other flow module
  !                       instances in connection with a particular attribute. If this
  !                       is the case, the updating of the flow module instance is
  !                       still attempted, with the flow module instance deciding
  !                       whether updating is possible.
  ! i                  :: Loop index.

  CommonFlow => Flows%C(iFlowMod, iFlow)%P

  ! Update met module instances used by flow module if these are due for update.
  Do i = 1, CommonFlow%nMets
    If (Mets%C(CommonFlow%iMetMods(i), CommonFlow%iMets(i))%P%DueForUpdate) Then
      Call SetMetLock(CommonFlow%iMetMods(i), CommonFlow%iMets(i))
      Call UpdateMet(                                     &
             Coords, Grids,                               &
             iCase, iMetCase,                             &
             CommonFlow%iMetMods(i), CommonFlow%iMets(i), &
             MetTime,                                     &
             Mets,                                        &
             Units                                        &
           )
      Call UnsetMetLock(CommonFlow%iMetMods(i), CommonFlow%iMets(i))
    End If
  End Do

  If (.false.) Then

  Else If (Flows%nPrototypeFlows /= 0 .and. iFlowMod == Flows%PrototypeFlows(1)%C%iFlowMod) Then

    Call PrototypeFlowReqs(                                               &
           Flows%PrototypeFlows(iFlow),                                   &
           MetTime,                                                  &
           WantFlowField,     WantCloudField,     WantRainField,     &
           iFlowUpdateSubset, iCloudUpdateSubset, iRainUpdateSubset, &
           FlowField,         CloudField,         RainField          &
         )

    If (WantFlowField) Then
      If (iFlowUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iFlowUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                 &
             A_Flow,                                       &
             Coords, Grids, Domains,                       &
             CommonFlow%iUpdateSubsets(iFlowUpdateSubset), &
             MetTime,                                      &
             Units, Mets, Flows,                           &
             FlowField, CloudField, RainField,             &
             Error                                         &
           )
      FlowField%C%Valid = .not.Error
    End If

    If (WantCloudField) Then
      If (iCloudUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iCloudUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                  &
             A_Cloud,                                       &
             Coords, Grids, Domains,                        &
             CommonFlow%iUpdateSubsets(iCloudUpdateSubset), &
             MetTime,                                       &
             Units, Mets, Flows,                            &
             FlowField, CloudField, RainField,              &
             Error                                          &
           )
      CloudField%C%Valid = .not.Error
    End If

    If (WantRainField) Then
      If (iRainUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iRainUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                 &
             A_Rain,                                       &
             Coords, Grids, Domains,                       &
             CommonFlow%iUpdateSubsets(iRainUpdateSubset), &
             MetTime,                                      &
             Units, Mets, Flows,                           &
             FlowField, CloudField, RainField,             &
             Error                                         &
           )
      RainField%C%Valid = .not.Error
    End If

    Call UpdatePrototypeFlow(                     &
           Coords, Grids, Domains,           &
           Mets,                             &
           iCase,                            &
           MetTime,                          &
           FlowField, CloudField, RainField, &
           Flows%PrototypeFlows(iFlow),           &
           Units                             &
         )

  Else If (Flows%nSingleSiteFlows /= 0 .and. iFlowMod == Flows%SingleSiteFlows(1)%C%iFlowMod) Then

    Call SingleSiteFlowReqs(                                               &
           Flows%SingleSiteFlows(iFlow),                                   &
           MetTime,                                                  &
           WantFlowField,     WantCloudField,     WantRainField,     &
           iFlowUpdateSubset, iCloudUpdateSubset, iRainUpdateSubset, &
           FlowField,         CloudField,         RainField          &
         )

    If (WantFlowField) Then
      If (iFlowUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iFlowUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                 &
             A_Flow,                                       &
             Coords, Grids, Domains,                       &
             CommonFlow%iUpdateSubsets(iFlowUpdateSubset), &
             MetTime,                                      &
             Units, Mets, Flows,                           &
             FlowField, CloudField, RainField,             &
             Error                                         &
           )
      FlowField%C%Valid = .not.Error
    End If

    If (WantCloudField) Then
      If (iCloudUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iCloudUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                  &
             A_Cloud,                                       &
             Coords, Grids, Domains,                        &
             CommonFlow%iUpdateSubsets(iCloudUpdateSubset), &
             MetTime,                                       &
             Units, Mets, Flows,                            &
             FlowField, CloudField, RainField,              &
             Error                                          &
           )
      CloudField%C%Valid = .not.Error
    End If

    If (WantRainField) Then
      If (iRainUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iRainUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                 &
             A_Rain,                                       &
             Coords, Grids, Domains,                       &
             CommonFlow%iUpdateSubsets(iRainUpdateSubset), &
             MetTime,                                      &
             Units, Mets, Flows,                           &
             FlowField, CloudField, RainField,             &
             Error                                         &
           )
      RainField%C%Valid = .not.Error
    End If

    Call UpdateSingleSiteFlow(                     &
           Coords, Grids, Domains,           &
           Mets,                             &
           iCase,                            &
           MetTime,                          &
           FlowField, CloudField, RainField, &
           Flows%SingleSiteFlows(iFlow),           &
           Units                             &
         )

  Else If (Flows%nNWPFlows /= 0 .and. iFlowMod == Flows%NWPFlows(1)%C%iFlowMod) Then

    Call NWPFlowReqs(                                               &
           Flows%NWPFlows(iFlow),                                   &
           MetTime,                                                  &
           WantFlowField,     WantCloudField,     WantRainField,     &
           iFlowUpdateSubset, iCloudUpdateSubset, iRainUpdateSubset, &
           FlowField,         CloudField,         RainField          &
         )

    If (WantFlowField) Then
      If (iFlowUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iFlowUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                 &
             A_Flow,                                       &
             Coords, Grids, Domains,                       &
             CommonFlow%iUpdateSubsets(iFlowUpdateSubset), &
             MetTime,                                      &
             Units, Mets, Flows,                           &
             FlowField, CloudField, RainField,             &
             Error                                         &
           )
      FlowField%C%Valid = .not.Error
    End If

    If (WantCloudField) Then
      If (iCloudUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iCloudUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                  &
             A_Cloud,                                       &
             Coords, Grids, Domains,                        &
             CommonFlow%iUpdateSubsets(iCloudUpdateSubset), &
             MetTime,                                       &
             Units, Mets, Flows,                            &
             FlowField, CloudField, RainField,              &
             Error                                          &
           )
      CloudField%C%Valid = .not.Error
    End If

    If (WantRainField) Then
      If (iRainUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iRainUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                 &
             A_Rain,                                       &
             Coords, Grids, Domains,                       &
             CommonFlow%iUpdateSubsets(iRainUpdateSubset), &
             MetTime,                                      &
             Units, Mets, Flows,                           &
             FlowField, CloudField, RainField,             &
             Error                                         &
           )
      RainField%C%Valid = .not.Error
    End If

    Call UpdateNWPFlow(                     &
           Coords, Grids, Domains,           &
           Mets,                             &
           iCase,                            &
           MetTime,                          &
           FlowField, CloudField, RainField, &
           Flows%NWPFlows(iFlow),           &
           Units                             &
         )

  Else If (Flows%nBuildingFlows /= 0 .and. iFlowMod == Flows%BuildingFlows(1)%C%iFlowMod) Then

    Call BuildingFlowReqs(                                               &
           Flows%BuildingFlows(iFlow),                                   &
           MetTime,                                                  &
           WantFlowField,     WantCloudField,     WantRainField,     &
           iFlowUpdateSubset, iCloudUpdateSubset, iRainUpdateSubset, &
           FlowField,         CloudField,         RainField          &
         )

    If (WantFlowField) Then
      If (iFlowUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iFlowUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                 &
             A_Flow,                                       &
             Coords, Grids, Domains,                       &
             CommonFlow%iUpdateSubsets(iFlowUpdateSubset), &
             MetTime,                                      &
             Units, Mets, Flows,                           &
             FlowField, CloudField, RainField,             &
             Error                                         &
           )
      FlowField%C%Valid = .not.Error
    End If

    If (WantCloudField) Then
      If (iCloudUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iCloudUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                  &
             A_Cloud,                                       &
             Coords, Grids, Domains,                        &
             CommonFlow%iUpdateSubsets(iCloudUpdateSubset), &
             MetTime,                                       &
             Units, Mets, Flows,                            &
             FlowField, CloudField, RainField,              &
             Error                                          &
           )
      CloudField%C%Valid = .not.Error
    End If

    If (WantRainField) Then
      If (iRainUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iRainUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                 &
             A_Rain,                                       &
             Coords, Grids, Domains,                       &
             CommonFlow%iUpdateSubsets(iRainUpdateSubset), &
             MetTime,                                      &
             Units, Mets, Flows,                           &
             FlowField, CloudField, RainField,             &
             Error                                         &
           )
      RainField%C%Valid = .not.Error
    End If

    Call UpdateBuildingFlow(                     &
           Coords, Grids, Domains,           &
           Mets,                             &
           iCase,                            &
           MetTime,                          &
           FlowField, CloudField, RainField, &
           Flows%BuildingFlows(iFlow),           &
           Units                             &
         )

  Else If (Flows%nRadarFlows /= 0 .and. iFlowMod == Flows%RadarFlows(1)%C%iFlowMod) Then

    Call RadarFlowReqs(                                               &
           Flows%RadarFlows(iFlow),                                   &
           MetTime,                                                  &
           WantFlowField,     WantCloudField,     WantRainField,     &
           iFlowUpdateSubset, iCloudUpdateSubset, iRainUpdateSubset, &
           FlowField,         CloudField,         RainField          &
         )

    If (WantFlowField) Then
      If (iFlowUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iFlowUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                 &
             A_Flow,                                       &
             Coords, Grids, Domains,                       &
             CommonFlow%iUpdateSubsets(iFlowUpdateSubset), &
             MetTime,                                      &
             Units, Mets, Flows,                           &
             FlowField, CloudField, RainField,             &
             Error                                         &
           )
      FlowField%C%Valid = .not.Error
    End If

    If (WantCloudField) Then
      If (iCloudUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iCloudUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                  &
             A_Cloud,                                       &
             Coords, Grids, Domains,                        &
             CommonFlow%iUpdateSubsets(iCloudUpdateSubset), &
             MetTime,                                       &
             Units, Mets, Flows,                            &
             FlowField, CloudField, RainField,              &
             Error                                          &
           )
      CloudField%C%Valid = .not.Error
    End If

    If (WantRainField) Then
      If (iRainUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iRainUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                 &
             A_Rain,                                       &
             Coords, Grids, Domains,                       &
             CommonFlow%iUpdateSubsets(iRainUpdateSubset), &
             MetTime,                                      &
             Units, Mets, Flows,                           &
             FlowField, CloudField, RainField,             &
             Error                                         &
           )
      RainField%C%Valid = .not.Error
    End If

    Call UpdateRadarFlow(                     &
           Coords, Grids, Domains,           &
           Mets,                             &
           iCase,                            &
           MetTime,                          &
           FlowField, CloudField, RainField, &
           Flows%RadarFlows(iFlow),           &
           Units                             &
         )

  Else If (Flows%nLINCOMFlows /= 0 .and. iFlowMod == Flows%LINCOMFlows(1)%C%iFlowMod) Then

    Call LINCOMFlowReqs(                                               &
           Flows%LINCOMFlows(iFlow),                                   &
           MetTime,                                                  &
           WantFlowField,     WantCloudField,     WantRainField,     &
           iFlowUpdateSubset, iCloudUpdateSubset, iRainUpdateSubset, &
           FlowField,         CloudField,         RainField          &
         )

    If (WantFlowField) Then
      If (iFlowUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iFlowUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                 &
             A_Flow,                                       &
             Coords, Grids, Domains,                       &
             CommonFlow%iUpdateSubsets(iFlowUpdateSubset), &
             MetTime,                                      &
             Units, Mets, Flows,                           &
             FlowField, CloudField, RainField,             &
             Error                                         &
           )
      FlowField%C%Valid = .not.Error
    End If

    If (WantCloudField) Then
      If (iCloudUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iCloudUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                  &
             A_Cloud,                                       &
             Coords, Grids, Domains,                        &
             CommonFlow%iUpdateSubsets(iCloudUpdateSubset), &
             MetTime,                                       &
             Units, Mets, Flows,                            &
             FlowField, CloudField, RainField,              &
             Error                                          &
           )
      CloudField%C%Valid = .not.Error
    End If

    If (WantRainField) Then
      If (iRainUpdateSubset < 1 .or. CommonFlow%nUpdateSubsets < iRainUpdateSubset) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
      End If
      Call GetAttribField(                                 &
             A_Rain,                                       &
             Coords, Grids, Domains,                       &
             CommonFlow%iUpdateSubsets(iRainUpdateSubset), &
             MetTime,                                      &
             Units, Mets, Flows,                           &
             FlowField, CloudField, RainField,             &
             Error                                         &
           )
      RainField%C%Valid = .not.Error
    End If

    Call UpdateLINCOMFlow(                     &
           Coords, Grids, Domains,           &
           Mets,                             &
           iCase,                            &
           MetTime,                          &
           FlowField, CloudField, RainField, &
           Flows%LINCOMFlows(iFlow),           &
           Units                             &
         )


  Else
    Call Message('UNEXPECTED FATAL ERROR in UpdateFlow: flow module instance not found', 4)
  End If

  ! Check TValid is the right type of time.
  If (IsTimeInterval(CommonFlow%TValid)) Then
    Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
  End If

  ! Check that TValid is in the future, and, for fixed met, that TValid is infinite and ValidityUmimprovable
  ! is true.
  If (CommonFlow%TValid <= MetTime) Then
    Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
  End If
  If (CommonFlow%FixedMet) Then
    If (.not. IsInfFuture(CommonFlow%TValid) .or. .not. CommonFlow%ValidityUnimprovable) Then
      Call Message('UNEXPECTED FATAL ERROR in UpdateFlow', 4)
    End If
  End If

  ! Messages for whether flow update was sucessful.
  If (CommonFlow%Valid) Then
    If (CommonFlow%FixedMet) Then
      Call Message(                          &
             'Flow module "'              // &
             Trim(CommonFlow%FlowModName) // &
             '.'                          // &
             Trim(CommonFlow%FlowName)    // &
             '" now valid'                   & ! $$ list valid attributes
           )
    Else
      Call Message(                                                    &
             'Flow module "'                                        // &
             Trim(CommonFlow%FlowModName)                           // &
             '.'                                                    // &
             Trim(CommonFlow%FlowName)                              // &
             '" now valid until '                                   // &
             Trim(Time2Char(CommonFlow%TValid, .false., 0, .true.))    & ! $$ list valid attributes
           )
    End If
  Else
    Call Message(                              &
           'Unable to update flow module "' // &
           Trim(CommonFlow%FlowModName)     // &
           '.'                              // &
           Trim(CommonFlow%FlowName)        // &
           '"'                                 &
         )
  EndIf

  ! Update DueForUpdate.
  CommonFlow%DueForUpdate = .false.

End Subroutine UpdateFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine ResetFlowMemory(Flows, FlowMemory)
! Resets the flow memory (for use when the space-time point which the memory applies
! to changes).

  Implicit None
  ! Argument list:
  Type(Flows_),      Intent(In)    :: Flows      ! Collection of flow module instance
                                                 ! states.
  Type(FlowMemory_), Intent(InOut) :: FlowMemory ! Flow memory.
  ! Locals:
  Integer :: iFlowMod ! Flow module index.
  Integer :: iFlow    ! Flow module instance index.
  Integer :: iAttrib  ! Attribute index.

  Do iFlowMod = 1, Flows%nFlowMods
  Do iFlow    = 1, Flows%nFlows(iFlowMod)
    FlowMemory%InValid(iFlowMod, iFlow) = .false.
  End Do
  End Do

  Do iFlow = 1, Flows%nNWPFlows
    Call ResetNWPFlowMemory(FlowMemory%NWPFlowMemories(iFlow))
  End Do

  Do iAttrib = 1, Flows%nAttribs
    FlowMemory%iValid(iAttrib) = .false.
  End Do

End Subroutine ResetFlowMemory

!-------------------------------------------------------------------------------------------------------------

Subroutine ResetFlowMemorySpecificFlow(Flows, iFlowMod, iFlow, FlowMemory)
! Resets the flow memory of a specific flow module instance (for use in WhichFlow in
! situations where the flow memory may have been set while performing a vertical
! coordinate conversion which turns out to have been performed using the wrong flow
! module instance).

  Implicit None
  ! Argument list:
  Type(Flows_),      Intent(In)    :: Flows      ! Collection of flow module instance
                                                 ! states.
  Integer,           Intent(In)    :: iFlowMod   !} Indices of the flow module
  Integer,           Intent(In)    :: iFlow      !} instance whose memory is to be
                                                 !} reset.
  Type(FlowMemory_), Intent(InOut) :: FlowMemory ! Flow memory.

  If (.false.) Then

  Else If (                                   &
    Flows%nNWPFlows /= 0 .and.               &
    iFlowMod == Flows%NWPFlows(1)%C%iFlowMod &
  ) Then

    Call ResetNWPFlowMemory(FlowMemory%NWPFlowMemories(iFlow))


  ! Check flow module present.
  Else If (.not. FlowPresent(iFlowMod, iFlow, Flows)) Then
    Call Message('UNEXPECTED FATAL ERROR in ResetFlowMemorySpecificFlow: flow module instance not found', 4)
  End If

End Subroutine ResetFlowMemorySpecificFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine ConvertToZ(                                  &
             Coords, Grids, Domains,                    &
             iZCoord,                                   &
             Time, AnyTravelTime, TravelTime, Position, &
             Units, Mets, Flows,                        &
             FlowMemory,                                &
             ErrorCode                                  &
           )
! Converts the vertical coordinate of a point to a particular coord system.

  Implicit None
  ! Argument list:
  Type(Coords_),     Intent(In)    :: Coords
  Type(Grids_),      Intent(In)    :: Grids
  Type(Domains_),    Intent(In)    :: Domains
  Integer,           Intent(In)    :: iZCoord
  Type(ShortTime_),  Intent(In)    :: Time
  Logical,           Intent(In)    :: AnyTravelTime
  Type(ShortTime_),  Intent(In)    :: TravelTime
  Type(Position_),   Intent(InOut) :: Position
  Type(Units_),      Intent(InOut) :: Units
  Type(Mets_),       Intent(InOut) :: Mets
  Type(Flows_),      Intent(InOut) :: Flows
  Type(FlowMemory_), Intent(InOut) :: FlowMemory
  Integer,           Intent(Out)   :: ErrorCode
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! Domains       :: Collection of domains.
  ! iZCoord       :: Index in Coords of coord system to convert to.
  ! Time          :: Current time or time of the met to be used for performing the
  !                  vertical coord conversion (these must be the same unless we have
  !                  a fixed met case; for fixed met cases the value is ignored).
  ! AnyTravelTime :: Indicates that choices of flow module instances are not to be
  !                  restricted by travel time.
  ! TravelTime    :: Travel time.
  ! Position      :: Position.
  ! Units         :: Collection of information on input/output unit numbers.
  ! Mets          :: Collection of met module instance states.
  ! Flows         :: Collection of flow module instance states.
  ! FlowMemory    :: Flow memory.
  ! ErrorCode     :: Error code. 0 = no error, 1 = no flow module for convert attribute.
  ! Locals:
  Type(ShortTime_) :: MetTime  ! Time of the met to be used for performing the
                               ! vertical coord conversion.
  Integer          :: iFlowMod !} Indices of the flow module instance selected to
  Integer          :: iFlow    !} perform the vertical coord conversion.
  Logical          :: Error    ! Indicates no suitable flow module instance was found.

# ifdef ExtraChecks
    If (                        &
      iZCoord < 1 .or.          &
      Coords%nZCoords < iZCoord &
    ) Then
      Call Message('UNEXPECTED FATAL ERROR in ConvertToZ', 4)
    End If
# endif

  ErrorCode = 0

  If (.not.Position%ZValid(iZCoord)) Then

    If (Flows%C(1, 1)%P%FixedMet) Then
      MetTime = Flows%MetTime
    Else
      MetTime = Time
    End If

    Call WhichFlow(                               &
           Coords, Grids, Domains,                &
           Flows%iAttribs(A_Convert), 0, MetTime, &
           AnyTravelTime, TravelTime,             &
           iFlowMod, iFlow, Error,                &
           Position,                              &
           Units, Mets, Flows,                    &
           FlowMemory                             &
         )
    If (Error) Then
      ErrorCode = 1
      Return
    End If

    Call ConvertToZKnownFlow( &
           Flows,             &
           iFlowMod, iFlow,   &
           Coords, Grids,     &
           iZCoord,           &
           MetTime, Position, &
           FlowMemory         &
         )

  End If

End Subroutine ConvertToZ

!-------------------------------------------------------------------------------------------------------------

Subroutine ConvertToZKnownFlow( &
             Flows,             &
             iFlowMod, iFlow,   &
             Coords, Grids,     &
             iZCoord,           &
             MetTime, Position, &
             FlowMemory         &
           )
! Converts the vertical coordinate of a point to a particular coord system in
! situations where the flow module instance to be used for the conversion is known.

! Note Position%ZValid(iZCoord) is set to .true. and FlowMemory may be updated, but
! this may not be appropriate if the point lies outside the domain of the flow module
! instance selected to perform the conversion. The calling routine needs to correct
! this where appropriate.

  Implicit None
  ! Argument list:
  Type(Flows_),      Intent(In)    :: Flows      ! Collection of flow module instance
                                                 ! states.
  Integer,           Intent(In)    :: iFlowMod   !} Indices of the flow module
  Integer,           Intent(In)    :: iFlow      !} instance selected to perform the
                                                 !} conversion.
  Type(Coords_),     Intent(In)    :: Coords     ! A collection of coord systems.
  Type(Grids_),      Intent(In)    :: Grids      ! A collection of grids.
  Integer,           Intent(In)    :: iZCoord    ! Index in Coords of coord system to
                                                 ! convert to.
  Type(ShortTime_),  Intent(In)    :: MetTime    ! Time of the met to be used for
                                                 ! performing the vertical coord
                                                 ! conversion.
  Type(Position_),   Intent(InOut) :: Position   ! Position.
  Type(FlowMemory_), Intent(InOut) :: FlowMemory ! Flow memory.
  ! Locals:
  Integer :: nHCoords             !} Number and indices of the horizontal coords which
  Integer :: iHCoords(MaxHCoords) !} are needed for performing the vertical coord
                                  !} conversion.
  Integer :: i                    ! Loop index.

# ifdef ExtraChecks
    If (                        &
      iZCoord < 1 .or.          &
      Coords%nZCoords < iZCoord &
    ) Then
      Call Message('UNEXPECTED FATAL ERROR in ConvertToZKnownFlow', 4)
    End If
# endif

  If (.not.Position%ZValid(iZCoord)) Then

    If (.false.) Then

    Else If (                                   &
      Flows%nPrototypeFlows /= 0 .and.               &
      iFlowMod == Flows%PrototypeFlows(1)%C%iFlowMod &
    ) Then

      Call PrototypeConvertCoordIndices(Flows%PrototypeFlows(iFlow), nHCoords, iHCoords)
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Call PrototypeConvertToZ(                          &
             Flows%PrototypeFlows(iFlow), Coords, Grids, &
             iZCoord, MetTime, Position             &
           )

    Else If (                                   &
      Flows%nSingleSiteFlows /= 0 .and.               &
      iFlowMod == Flows%SingleSiteFlows(1)%C%iFlowMod &
    ) Then

      Call SingleSiteConvertCoordIndices(Flows%SingleSiteFlows(iFlow), nHCoords, iHCoords)
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Call SingleSiteConvertToZ(                          &
             Flows%SingleSiteFlows(iFlow), Coords, Grids, &
             iZCoord, MetTime, Position             &
           )


    Else If (                                   &
      Flows%nNWPFlows /= 0 .and.               &
      iFlowMod == Flows%NWPFlows(1)%C%iFlowMod &
    ) Then

      Call NWPConvertCoordIndices(Flows%NWPFlows(iFlow), nHCoords, iHCoords)
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Call NWPConvertToZ(                          &
             Flows%NWPFlows(iFlow), Coords, Grids, &
             iZCoord, MetTime, Position,            &
             FlowMemory%NWPFlowMemories(iFlow)     &
           )


    Else
      Call Message('UNEXPECTED FATAL ERROR in ConvertToZKnownFlow: flow module instance not found', 4)
    End If

    Position%ZValid(iZCoord) = .true.

  End If

End Subroutine ConvertToZKnownFlow

!-------------------------------------------------------------------------------------------------------------

Recursive Subroutine GetAttribField(                     &
                       iAttribParam,                     &
                       Coords, Grids, Domains,           &
                       iSubset,                          &
                       MetTime,                          &
                       Units, Mets, Flows,               &
                       FlowField, CloudField, RainField, &
                       Error                             &
                     )
! Gets attribute field information (used for transferring flow information between flow module instances).

! Information on unrequired attributes is left unchanged.

  Implicit None
  ! Argument list:
  Integer,           Intent(In)            :: iAttribParam
  Type(Coords_),     Intent(In)            :: Coords
  Type(Grids_),      Intent(In)            :: Grids
  Type(Domains_),    Intent(In)            :: Domains
  Integer,           Intent(In)            :: iSubset
  Type(Time_),       Intent(In)            :: MetTime
  Type(Units_),      Intent(InOut)         :: Units
  Type(Mets_),       Intent(InOut)         :: Mets
  Type(Flows_),      Intent(InOut)         :: Flows
  Type(FlowField_),  Intent(InOut), Target :: FlowField
  Type(CloudField_), Intent(InOut), Target :: CloudField
  Type(RainField_),  Intent(InOut), Target :: RainField
  Logical,           Intent(Out)           :: Error
  ! iAttribParam :: Index of required attribute (in the parameter list in the common flow module).
  ! Flows        :: Collection of flow module instance states.
  ! Coords       :: Collection of coord systems.
  ! Grids        :: Collection of grids.
  ! Domains      :: Collection of domains.
  ! iSubset      :: Index of subset giving a restricted choice of flow module instances.
  ! MetTime      :: Time of the met to be used for supplying the attribute information.
  ! Units        :: Collection of information on input/output unit numbers.
  ! Mets         :: Collection of met module instance states.
  ! Flows        :: Collection of flow module instance states.
  ! FlowField    :: Field of flow information.
  ! CloudField   :: Field of cloud information.
  ! RainField    :: Field of rain information.
  ! Error        :: Indicates no suitable flow module instance was found.
  ! Locals:
  Type(CommonAttribField_), Pointer :: C               ! The part of the attribute field common to all
                                                       ! attribute fields.
  Type(Flow_)                       :: Flow            !} Dummy variables for arguments in GetAttribKnownFlow.
  Type(Cloud_)                      :: Cloud           !}
  Type(Rain_)                       :: Rain            !}
  Type(Surface_)                    :: Surface         !}
  Type(Soil_)                       :: Soil            !}
  Type(Plant_)                      :: Plant           !}
  Integer                           :: iFlowModConvert !} Indices of the flow module instance selected
  Integer                           :: iFlowConvert    !} to perform any vertical coord conversions.
  Integer                           :: iFlowMod        !] Indices of the flow module instance selected
  Integer                           :: iFlow           !] to supply the attribute information.
  Type(ShortTime_)                  :: TravelTime      ! Dummy variable for TravelTime argument in WhichFlow.
  Type(Position_)                   :: Position        ! Position.
  Type(FlowMemory_)                 :: FlowMemory      ! Flow memory.
  Integer                           :: iX              !} Indices of grid point.
  Integer                           :: iY              !}
  Integer                           :: iZ              !}

  Error = .false.

  If (iAttribParam == A_Flow) Then
    C => FlowField%C
  Else If (iAttribParam == A_Cloud) Then
    C => CloudField%C
  Else If (iAttribParam == A_Rain) Then
    C => RainField%C
  Else
    Call Message('UNEXPECTED FATAL ERROR in GetAttribField', 4)
  End If

  If (Flows%iAttribs(iAttribParam) == 0) Then
    Error = .true.
    Return
  End If

  If (C%iZGrid < 1 .or. Grids%nZGrids < C%iZGrid .or. C%iHGrid < 1 .or. Grids%nHGrids < C%iHGrid) Then
    Call Message('UNEXPECTED FATAL ERROR in GetAttribField', 4)
  End If

  ! Attribute field at first time.

  C%Time1      = MetTime
  C%ShortTime1 = Time2ShortTime(MetTime)

  Do iX = 1, Grids%HGrids(C%iHGrid)%nX
  Do iY = 1, Grids%HGrids(C%iHGrid)%nY
  Do iZ = 1, Grids%ZGrids(C%iZGrid)%nZ

    Call ResetFlowMemory(Flows, FlowMemory)

    Position = X2Position(                       &
                 Coords,                         &
                 (/                              &
                   Grids%HGrids(C%iHGrid)%X(iX), &
                   Grids%HGrids(C%iHGrid)%Y(iY), &
                   Grids%ZGrids(C%iZGrid)%Z(iZ)  &
                 /),                             &
                 Grids%HGrids(C%iHGrid)%iHCoord, &
                 Grids%ZGrids(C%iZGrid)%iZCoord  &
               )

    ! Calculate best flow module instance for the convert attribute.
    Call WhichFlow(                              &
           Coords, Grids, Domains,               &
           Flows%iAttribs(A_Convert), iSubset,   &
           C%ShortTime1,                         &
           .true., TravelTime,                   &
           iFlowModConvert, iFlowConvert, Error, &
           Position,                             &
           Units, Mets, Flows,                   &
           FlowMemory                            &
         )
    If (Error) Return

    ! Calculate best flow module instance for the required attribute.
    Call WhichFlow(                               &
           Coords, Grids, Domains,                &
           Flows%iAttribs(iAttribParam), iSubset, &
           C%ShortTime1,                          &
           .true., TravelTime,                    &
           iFlowMod, iFlow, Error,                &
           Position,                              &
           Units, Mets, Flows,                    &
           FlowMemory                             &
         )
    If (Error) Return

    ! Get attribute info.
    If (iAttribParam == A_Flow) Then
      Call GetAttribKnownFlow(                    &
             iAttribParam,                        &
             Flows,                               &
             iFlowModConvert, iFlowConvert,       &
             iFlowMod,        iFlow,              &
             Coords, Grids,                       &
             .true., .true., .true.,              &
             C%ShortTime1, Position,              &
             FlowMemory,                          &
             FlowField%Flow(iX, iY, iZ, 1),       &
             Cloud,                               &
             Rain,                                &
             Surface,                             &
             Soil,                                &
             Plant,                               &
             FlowField%ProfileData(iX, iY, iZ, 1) &
           )
    Else If (iAttribParam == A_Cloud) Then
      Call GetAttribKnownFlow(                &
             iAttribParam,                    &
             Flows,                           &
             iFlowModConvert, iFlowConvert,   &
             iFlowMod,        iFlow,          &
             Coords, Grids,                   &
             .false., .false., .false.,       &
             C%ShortTime1, Position,          &
             FlowMemory,                      &
             Flow,                            &
             CloudField%Cloud(iX, iY, iZ, 1), &
             Rain,                            &
             Surface,                         &
             Soil,                            &
             Plant                            &
           )
    Else If (iAttribParam == A_Rain) Then
      Call GetAttribKnownFlow(              &
             iAttribParam,                  &
             Flows,                         &
             iFlowModConvert, iFlowConvert, &
             iFlowMod,        iFlow,        &
             Coords, Grids,                 &
             .false., .false., .false.,     &
             C%ShortTime1, Position,        &
             FlowMemory,                    &
             Flow,                          &
             Cloud,                         &
             RainField%Rain(iX, iY, iZ, 1), &
             Surface,                       &
             Soil,                          &
             Plant                          &
            )
    End If

  End Do
  End Do
  End Do

  ! Attribute field at second time.

  If (C%UseTwoTimes) Then

    If (Flows%C(1, 1)%P%FixedMet) Then
      Call Message('UNEXPECTED FATAL ERROR in GetAttribField', 4)
    End If

    C%Time2      = C%Time1 + C%Dt
    C%Time2      = TMin(C%Time2, Flows%OverallTValid)
    C%ShortTime2 = Time2ShortTime(C%Time2)

    Do iX = 1, Grids%HGrids(C%iHGrid)%nX
    Do iY = 1, Grids%HGrids(C%iHGrid)%nY
    Do iZ = 1, Grids%ZGrids(C%iZGrid)%nZ

      Call ResetFlowMemory(Flows, FlowMemory)

      Position = X2Position(                       &
                   Coords,                         &
                   (/                              &
                     Grids%HGrids(C%iHGrid)%X(iX), &
                     Grids%HGrids(C%iHGrid)%Y(iY), &
                     Grids%ZGrids(C%iZGrid)%Z(iZ)  &
                   /),                             &
                   Grids%HGrids(C%iHGrid)%iHCoord, &
                   Grids%ZGrids(C%iZGrid)%iZCoord  &
                 )

      ! Calculate best flow module instance for the convert attribute.
      Call WhichFlow(                              &
             Coords, Grids, Domains,               &
             Flows%iAttribs(A_Convert), iSubset,   &
             C%ShortTime2,                         &
             .true., TravelTime,                   &
             iFlowModConvert, iFlowConvert, Error, &
             Position,                             &
             Units, Mets, Flows,                   &
             FlowMemory                            &
           )
      If (Error) Return

      ! Calculate best flow module instance for the required attribute.
      Call WhichFlow(                               &
             Coords, Grids, Domains,                &
             Flows%iAttribs(iAttribParam), iSubset, &
             C%ShortTime2,                          &
             .true., TravelTime,                    &
             iFlowMod, iFlow, Error,                &
             Position,                              &
             Units, Mets, Flows,                    &
             FlowMemory                             &
           )
      If (Error) Return

      ! Get attribute info.
      If (iAttribParam == A_Flow) Then
        Call GetAttribKnownFlow(                    &
               iAttribParam,                        &
               Flows,                               &
               iFlowModConvert, iFlowConvert,       &
               iFlowMod,        iFlow,              &
               Coords, Grids,                       &
               .true., .true., .true.,              &
               C%ShortTime2, Position,              &
               FlowMemory,                          &
               FlowField%Flow(iX, iY, iZ, 2),       &
               Cloud,                               &
               Rain,                                &
               Surface,                             &
               Soil,                                &
               Plant,                               &
               FlowField%ProfileData(iX, iY, iZ, 2) &
             )
      Else If (iAttribParam == A_Cloud) Then
        Call GetAttribKnownFlow(                &
               iAttribParam,                    &
               Flows,                           &
               iFlowModConvert, iFlowConvert,   &
               iFlowMod,        iFlow,          &
               Coords, Grids,                   &
               .false., .false., .false.,       &
               C%ShortTime2, Position,          &
               FlowMemory,                      &
               Flow,                            &
               CloudField%Cloud(iX, iY, iZ, 2), &
               Rain,                            &
               Surface,                         &
               Soil,                            &
               Plant                            &
             )
      Else If (iAttribParam == A_Rain) Then
        Call GetAttribKnownFlow(              &
               iAttribParam,                  &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               iFlowMod,        iFlow,        &
               Coords, Grids,                 &
               .false., .false., .false.,     &
               C%ShortTime2, Position,        &
               FlowMemory,                    &
               Flow,                          &
               Cloud,                         &
               RainField%Rain(iX, iY, iZ, 2), &
               Surface,                       &
               Soil,                          &
               Plant                          &
             )
      End If

    End Do
    End Do
    End Do

  End If

End Subroutine GetAttribField

!-------------------------------------------------------------------------------------------------------------

Subroutine GetAttrib(                                   &
             iAttribParam,                              &
             Coords, Grids, Domains,                    &
             Moisture, Inhomog, Homog,                  &
             Time, AnyTravelTime, TravelTime, Position, &
             Units, Mets, Flows,                        &
             FlowMemory,                                &
             Flow, Cloud, Rain, Surface, Soil, Plant,   &
             ErrorCode,                                 &
             ProfileData                                &
           )
! Gets attribute information.

! Information on unrequired attributes is left unchanged.

  Implicit None
  ! Argument list:
  Integer,            Intent(In)             :: iAttribParam
  Type(Coords_),      Intent(In)             :: Coords
  Type(Grids_),       Intent(In)             :: Grids
  Type(Domains_),     Intent(In)             :: Domains
  Logical,            Intent(In)             :: Moisture
  Logical,            Intent(In)             :: Inhomog
  Logical,            Intent(In)             :: Homog
  Type(ShortTime_),   Intent(In)             :: Time
  Logical,            Intent(In)             :: AnyTravelTime
  Type(ShortTime_),   Intent(In)             :: TravelTime
  Type(Position_),    Intent(InOut)          :: Position
  Type(Units_),       Intent(InOut)          :: Units
  Type(Mets_),        Intent(InOut)          :: Mets
  Type(Flows_),       Intent(InOut)          :: Flows
  Type(FlowMemory_),  Intent(InOut)          :: FlowMemory
  Type(Flow_),        Intent(Out)            :: Flow
  Type(Cloud_),       Intent(Out)            :: Cloud
  Type(Rain_),        Intent(Out)            :: Rain
  Type(Surface_),     Intent(Out)            :: Surface
  Type(Soil_),        Intent(Out)            :: Soil
  Type(Plant_),       Intent(Out)            :: Plant
  Integer,            Intent(Out)            :: ErrorCode
  Type(ProfileData_), Intent(Out),  Optional :: ProfileData
  ! iAttribParam  :: Index of required attribute (in the parameter list in the common flow module).
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! Domains       :: Collection of domains.
  ! Moisture      :: Indicates Q is required.
  ! Inhomog       :: Indicates inhomogeneous quantities are required.
  ! Homog         :: Indicates homogeneous quantities are required.
  ! Time          :: Current time or time of the met to be used for supplying the
  !                  attribute information (these must be the same unless we have a fixed
  !                  met case; for fixed met cases the value is ignored).
  ! AnyTravelTime :: Indicates that choices of flow module instances are not to be
  !                  restricted by travel time.
  ! TravelTime    :: Travel time.
  ! Position      :: Position.
  ! Units         :: Collection of information on input/output unit numbers.
  ! Mets          :: Collection of met module instance states.
  ! Flows         :: Collection of flow module instance states.
  ! FlowMemory    :: Flow memory.
  ! Flow          :: Flow information.
  ! Cloud         :: Cloud information.
  ! Rain          :: Rain information.
  ! Surface       :: Surface information.
  ! Soil          :: Soil information.
  ! Plant         :: Plant information.
  ! ErrorCode     :: Error code. 0 = no error, 1 = no flow module for convert attribute, 2 = no flow module
  !                  for required attribute, 3 = required attribute not defined in collection of attributes.
  ! ProfileData   :: Data needed to construct idealised analytic mean flow and turbulence profiles.
  ! Locals:
  Type(ShortTime_) :: MetTime         ! Time of the met to be used for supplying the attribute information.
  Integer          :: iFlowModConvert !} Indices of the flow module instance selected
  Integer          :: iFlowConvert    !} to perform any vertical coord conversions.
  Integer          :: iFlowMod        !] Indices of the flow module instance selected
  Integer          :: iFlow           !] to supply the attribute information.
  Logical          :: Error           ! Indicates no suitable flow module instance was found.

  ErrorCode = 0

  If      (iAttribParam == A_Flow   ) Then
  Else If (iAttribParam == A_Cloud  ) Then
  Else If (iAttribParam == A_Rain   ) Then
  Else If (iAttribParam == A_Surface) Then
  Else If (iAttribParam == A_Soil   ) Then
  Else If (iAttribParam == A_Plant  ) Then
  Else
    Call Message('UNEXPECTED FATAL ERROR in GetAttrib', 4)
  End If

  If (Flows%iAttribs(iAttribParam) == 0) Then
    ErrorCode = 3
    Return
  End If

  If (Flows%C(1, 1)%P%FixedMet) Then
    MetTime = Flows%MetTime
  Else
    MetTime = Time
  End If

  ! Calculate best flow module instance for the convert attribute.
  Call WhichFlow(                               &
         Coords, Grids, Domains,                &
         Flows%iAttribs(A_Convert), 0, MetTime, &
         AnyTravelTime, TravelTime,             &
         iFlowModConvert, iFlowConvert, Error,  &
         Position,                              &
         Units, Mets, Flows,                    &
         FlowMemory                             &
       )
  If (Error) Then
    ErrorCode = 1
    Return
  End If

  ! Calculate best flow module instance for the required attribute.
  Call WhichFlow(                                  &
         Coords, Grids, Domains,                   &
         Flows%iAttribs(iAttribParam), 0, MetTime, &
         AnyTravelTime, TravelTime,                &
         iFlowMod, iFlow, Error,                   &
         Position,                                 &
         Units, Mets, Flows,                       &
         FlowMemory                                &
       )
  If (Error) Then
    ErrorCode = 2
    Return
  End If

  Call GetAttribKnownFlow(              &
         iAttribParam,                  &
         Flows,                         &
         iFlowModConvert, iFlowConvert, &
         iFlowMod, iFlow,               &
         Coords, Grids,                 &
         Moisture, Inhomog, Homog,      &
         MetTime, Position,             &
         FlowMemory,                    &
         Flow,                          &
         Cloud,                         &
         Rain,                          &
         Surface,                       &
         Soil,                          &
         Plant,                         &
         ProfileData                    &
       )

  If (iAttribParam == A_Flow) Then
    Flow%Backwards      = .False.
    Position%Topog      = Flow%Topog
    Position%PS         = Flow%PS
    Position%Rho        = Flow%Rho
    Position%TopogValid = .true.
    Position%PSValid    = .true.
    Position%RhoValid   = .true.
  End If

End Subroutine GetAttrib

!-------------------------------------------------------------------------------------------------------------

Subroutine GetAttribKnownFlow(              &
             iAttribParam,                  &
             Flows,                         &
             iFlowModConvert, iFlowConvert, &
             iFlowMod, iFlow,               &
             Coords, Grids,                 &
             Moisture, Inhomog, Homog,      &
             MetTime, Position,             &
             FlowMemory,                    &
             Flow, Cloud, Rain, Surface,    &
             Soil, Plant, ProfileData       &
           )
! Gets attribute information in situations where the flow module instances to be used for vertical coord
! conversions and for supplying the information are known.

! Information on unrequired attributes is left unchanged.

  Implicit None
  ! Argument list:
  Integer,            Intent(In)             :: iAttribParam
  Type(Flows_),       Intent(In)             :: Flows
  Integer,            Intent(In)             :: iFlowModConvert
  Integer,            Intent(In)             :: iFlowConvert
  Integer,            Intent(In)             :: iFlowMod
  Integer,            Intent(In)             :: iFlow
  Type(Coords_),      Intent(In)             :: Coords
  Type(Grids_),       Intent(In)             :: Grids
  Logical,            Intent(In)             :: Moisture
  Logical,            Intent(In)             :: Inhomog
  Logical,            Intent(In)             :: Homog
  Type(ShortTime_),   Intent(In)             :: MetTime
  Type(Position_),    Intent(InOut)          :: Position
  Type(FlowMemory_),  Intent(InOut)          :: FlowMemory
  Type(Flow_),        Intent(Out)            :: Flow
  Type(Cloud_),       Intent(Out)            :: Cloud
  Type(Rain_),        Intent(Out)            :: Rain
  Type(Surface_),     Intent(Out)            :: Surface
  Type(Soil_),        Intent(Out)            :: Soil
  Type(Plant_),       Intent(Out)            :: Plant
  Type(ProfileData_), Intent(Out),  Optional :: ProfileData
  ! iAttribParam    :: Index of required attribute (in the parameter list in the common flow module).
  ! Flows           :: Collection of flow module instance states.
  ! iFlowModConvert :} Indices of the flow module instance selected to perform any
  ! iFlowConvert    :} vertical coord conversions.
  ! iFlowMod        :] Indices of the flow module instance selected to supply the
  ! iFlow           :] attribute information.
  ! Coords          :: Collection of coord systems.
  ! Grids           :: Collection of grids.
  ! Moisture        :: Indicates Q is required.
  ! Inhomog         :: Indicates inhomogeneous quantities are required.
  ! Homog           :: Indicates homogeneous quantities are required.
  ! MetTime         :: Time of the met to be used for supplying the attribute information.
  ! Position        :: Position.
  ! FlowMemory      :: Flow memory.
  ! Flow            :: Flow information.
  ! Cloud           :: Cloud information.
  ! Rain            :: Rain information.
  ! Surface         :: Surface information.
  ! Soil            :: Soil information.
  ! Plant           :: Plant information.
  ! ProfileData     :: Data needed to construct idealised analytic mean flow and turbulence profiles.
  ! Locals:
  Integer :: nHCoords             !} Number and indices of the coords which are needed for supplying the
  Integer :: iHCoords(MaxHCoords) !} attribute information.
  Integer :: nZCoords             !}
  Integer :: iZCoords(MaxZCoords) !}
  Integer :: i                    ! Loop index.

  ! Flow.
  If (iAttribParam == A_Flow) Then

    If (.false.) Then

    ! Note the next two blocks are identical except for the form of the call to Get????Flow.

    Else If (Flows%nPrototypeFlows /= 0 .and. iFlowMod == Flows%PrototypeFlows(1)%C%iFlowMod) Then

      Call PrototypeFlowCoordIndices(     &
             Flows%PrototypeFlows(iFlow), &
             nHCoords, iHCoords,     &
             nZCoords, iZCoords      &
           )
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Do i = 1, nZCoords
        Call ConvertToZKnownFlow(             &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               Coords, Grids,                 &
               iZCoords(i),                   &
               MetTime, Position,             &
               FlowMemory                     &
             )
      End Do
      Call GetPrototypeFlow(                &
             Coords, Grids,            &
             Flows%PrototypeFlows(iFlow),   &
             Moisture, Inhomog, Homog, &
             MetTime, Position,        &
             Flow, ProfileData         &
           )

    Else If (Flows%nSingleSiteFlows /= 0 .and. iFlowMod == Flows%SingleSiteFlows(1)%C%iFlowMod) Then

      Call SingleSiteFlowCoordIndices(     &
             Flows%SingleSiteFlows(iFlow), &
             nHCoords, iHCoords,     &
             nZCoords, iZCoords      &
           )
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Do i = 1, nZCoords
        Call ConvertToZKnownFlow(             &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               Coords, Grids,                 &
               iZCoords(i),                   &
               MetTime, Position,             &
               FlowMemory                     &
             )
      End Do
      Call GetSingleSiteFlow(                &
             Coords, Grids,            &
             Flows%SingleSiteFlows(iFlow),   &
             Moisture, Inhomog, Homog, &
             MetTime, Position,        &
             Flow, ProfileData         &
           )

    Else If (Flows%nBuildingFlows /= 0 .and. iFlowMod == Flows%BuildingFlows(1)%C%iFlowMod) Then

      Call BuildingFlowCoordIndices(     &
             Flows%BuildingFlows(iFlow), &
             nHCoords, iHCoords,     &
             nZCoords, iZCoords      &
           )
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Do i = 1, nZCoords
        Call ConvertToZKnownFlow(             &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               Coords, Grids,                 &
               iZCoords(i),                   &
               MetTime, Position,             &
               FlowMemory                     &
             )
      End Do
      Call GetBuildingFlow(                &
             Coords, Grids,            &
             Flows%BuildingFlows(iFlow),   &
             Moisture, Inhomog, Homog, &
             MetTime, Position,        &
             Flow, ProfileData         &
           )

    Else If (Flows%nLINCOMFlows /= 0 .and. iFlowMod == Flows%LINCOMFlows(1)%C%iFlowMod) Then

      Call LINCOMFlowCoordIndices(     &
             Flows%LINCOMFlows(iFlow), &
             nHCoords, iHCoords,     &
             nZCoords, iZCoords      &
           )
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Do i = 1, nZCoords
        Call ConvertToZKnownFlow(             &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               Coords, Grids,                 &
               iZCoords(i),                   &
               MetTime, Position,             &
               FlowMemory                     &
             )
      End Do
      Call GetLINCOMFlow(                &
             Coords, Grids,            &
             Flows%LINCOMFlows(iFlow),   &
             Moisture, Inhomog, Homog, &
             MetTime, Position,        &
             Flow, ProfileData         &
           )


    Else If (Flows%nNWPFlows /= 0 .and. iFlowMod == Flows%NWPFlows(1)%C%iFlowMod) Then

      Call NWPFlowCoordIndices(     &
             Flows%NWPFlows(iFlow), &
             nHCoords, iHCoords,     &
             nZCoords, iZCoords      &
           )
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Do i = 1, nZCoords
        Call ConvertToZKnownFlow(             &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               Coords, Grids,                 &
               iZCoords(i),                   &
               MetTime, Position,             &
               FlowMemory                     &
             )
      End Do
      Call GetNWPFlow(                          &
             Coords, Grids,                      &
             Flows%NWPFlows(iFlow),             &
             Moisture, Inhomog, Homog,           &
             MetTime, Position,                  &
             FlowMemory%NWPFlowMemories(iFlow), &
             Flow, ProfileData                   &
           )


    Else
      Call Message('UNEXPECTED FATAL ERROR in GetAttribKnownFlow: flow module instance not found', 4)
    End If

  ! Cloud.
  Else If (iAttribParam == A_Cloud) Then

    If (.false.) Then

    ! Note the next two blocks are identical except for the form of the call to Get????Cloud.

    Else If (Flows%nSingleSiteFlows /= 0 .and. iFlowMod == Flows%SingleSiteFlows(1)%C%iFlowMod) Then

      Call SingleSiteCloudCoordIndices(    &
             Flows%SingleSiteFlows(iFlow), &
             nHCoords, iHCoords,     &
             nZCoords, iZCoords      &
           )
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Do i = 1, nZCoords
        Call ConvertToZKnownFlow(             &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               Coords, Grids,                 &
               iZCoords(i),                   &
               MetTime, Position,             &
               FlowMemory                     &
             )
      End Do
      Call GetSingleSiteCloud(             &
             Coords, Grids,          &
             Flows%SingleSiteFlows(iFlow), &
             MetTime, Position,      &
             Cloud                   &
           )


    Else If (Flows%nNWPFlows /= 0 .and. iFlowMod == Flows%NWPFlows(1)%C%iFlowMod) Then

      Call NWPCloudCoordIndices(    &
             Flows%NWPFlows(iFlow), &
             nHCoords, iHCoords,     &
             nZCoords, iZCoords      &
           )
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Do i = 1, nZCoords
        Call ConvertToZKnownFlow(             &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               Coords, Grids,                 &
               iZCoords(i),                   &
               MetTime, Position,             &
               FlowMemory                     &
             )
      End Do
      Call GetNWPCloud(                         &
             Coords, Grids,                      &
             Flows%NWPFlows(iFlow),             &
             MetTime, Position,                  &
             FlowMemory%NWPFlowMemories(iFlow), &
             Cloud                               &
           )


    Else
      Call Message('UNEXPECTED FATAL ERROR in GetAttribKnownFlow: flow module instance not found', 4)
    End If

  ! Rain.
  Else If (iAttribParam == A_Rain) Then

    If (.false.) Then

    ! Note the next two blocks are identical except for the form of the call to Get????Rain.

    Else If (Flows%nSingleSiteFlows /= 0 .and. iFlowMod == Flows%SingleSiteFlows(1)%C%iFlowMod) Then

      Call SingleSiteRainCoordIndices(     &
             Flows%SingleSiteFlows(iFlow), &
             nHCoords, iHCoords,     &
             nZCoords, iZCoords      &
           )
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Do i = 1, nZCoords
        Call ConvertToZKnownFlow(             &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               Coords, Grids,                 &
               iZCoords(i),                   &
               MetTime, Position,             &
               FlowMemory                     &
             )
      End Do
      Call GetSingleSiteRain(              &
             Coords, Grids,          &
             Flows%SingleSiteFlows(iFlow), &
             MetTime, Position,      &
             Rain                    &
           )

    Else If (Flows%nRadarFlows /= 0 .and. iFlowMod == Flows%RadarFlows(1)%C%iFlowMod) Then

      Call RadarRainCoordIndices(     &
             Flows%RadarFlows(iFlow), &
             nHCoords, iHCoords,     &
             nZCoords, iZCoords      &
           )
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Do i = 1, nZCoords
        Call ConvertToZKnownFlow(             &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               Coords, Grids,                 &
               iZCoords(i),                   &
               MetTime, Position,             &
               FlowMemory                     &
             )
      End Do
      Call GetRadarRain(              &
             Coords, Grids,          &
             Flows%RadarFlows(iFlow), &
             MetTime, Position,      &
             Rain                    &
           )


    Else If (Flows%nNWPFlows /= 0 .and. iFlowMod == Flows%NWPFlows(1)%C%iFlowMod) Then

      Call NWPRainCoordIndices(     &
             Flows%NWPFlows(iFlow), &
             nHCoords, iHCoords,     &
             nZCoords, iZCoords      &
           )
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Do i = 1, nZCoords
        Call ConvertToZKnownFlow(             &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               Coords, Grids,                 &
               iZCoords(i),                   &
               MetTime, Position,             &
               FlowMemory                     &
             )
      End Do
      Call GetNWPRain(                          &
             Coords, Grids,                      &
             Flows%NWPFlows(iFlow),             &
             MetTime, Position,                  &
             FlowMemory%NWPFlowMemories(iFlow), &
             Rain                                &
           )


    Else
      Call Message('UNEXPECTED FATAL ERROR in GetAttribKnownFlow: flow module instance not found', 4)
    End If

  ! Surface.
  Else If (iAttribParam == A_Surface) Then

    If (.false.) Then

    ! Note the next two blocks are identical except for the form of the call to Get????Surface.


    Else If (Flows%nNWPFlows /= 0 .and. iFlowMod == Flows%NWPFlows(1)%C%iFlowMod) Then

      Call NWPSurfaceCoordIndices(  &
             Flows%NWPFlows(iFlow), &
             nHCoords, iHCoords,     &
             nZCoords, iZCoords      &
           )
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Do i = 1, nZCoords
        Call ConvertToZKnownFlow(             &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               Coords, Grids,                 &
               iZCoords(i),                   &
               MetTime, Position,             &
               FlowMemory                     &
             )
      End Do
      Call GetNWPSurface(                       &
             Coords, Grids,                      &
             Flows%NWPFlows(iFlow),             &
             MetTime, Position,                  &
             FlowMemory%NWPFlowMemories(iFlow), &
             Surface                             &
           )


    Else
      Call Message('UNEXPECTED FATAL ERROR in GetAttribKnownFlow: flow module instance not found', 4)
    End If

  ! Soil.
  Else If (iAttribParam == A_Soil) Then

    If (.false.) Then

    ! Note the next two blocks are identical except for the form of the call to Get????Soil.


    Else If (Flows%nNWPFlows /= 0 .and. iFlowMod == Flows%NWPFlows(1)%C%iFlowMod) Then

      Call NWPSoilCoordIndices(     &
             Flows%NWPFlows(iFlow), &
             nHCoords, iHCoords,     &
             nZCoords, iZCoords      &
           )
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Do i = 1, nZCoords
        Call ConvertToZKnownFlow(             &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               Coords, Grids,                 &
               iZCoords(i),                   &
               MetTime, Position,             &
               FlowMemory                     &
             )
      End Do
      Call GetNWPSoil(                          &
             Coords, Grids,                      &
             Flows%NWPFlows(iFlow),             &
             MetTime, Position,                  &
             FlowMemory%NWPFlowMemories(iFlow), &
             Soil                                &
           )


    Else
      Call Message('UNEXPECTED FATAL ERROR in GetAttribKnownFlow: flow module instance not found', 4)
    End If

  ! Plant.
  Else If (iAttribParam == A_Plant) Then

    If (.false.) Then

    ! Note the next two blocks are identical except for the form of the call to Get????Plant.


    Else If (Flows%nNWPFlows /= 0 .and. iFlowMod == Flows%NWPFlows(1)%C%iFlowMod) Then

      Call NWPPlantCoordIndices(    &
             Flows%NWPFlows(iFlow), &
             nHCoords, iHCoords,     &
             nZCoords, iZCoords      &
           )
      Do i = 1, nHCoords
        Call ConvertToH(Coords, iHCoords(i), Position)
      End Do
      Do i = 1, nZCoords
        Call ConvertToZKnownFlow(             &
               Flows,                         &
               iFlowModConvert, iFlowConvert, &
               Coords, Grids,                 &
               iZCoords(i),                   &
               MetTime, Position,             &
               FlowMemory                     &
             )
      End Do
      Call GetNWPPlant(                         &
             Coords, Grids,                      &
             Flows%NWPFlows(iFlow),             &
             MetTime, Position,                  &
             FlowMemory%NWPFlowMemories(iFlow), &
             Plant                               &
           )


    Else
      Call Message('UNEXPECTED FATAL ERROR in GetAttribKnownFlow: flow module instance not found', 4)
    End If

  End If

End Subroutine GetAttribKnownFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine Reflect(                           &
             Coords, Grids, Domains,          &
             Time, AnyTravelTime, TravelTime, &
             Vel, OldPosition, Position, U,   &
             Units, Mets, Flows,              &
             FlowMemory,                      &
             Reflected,                       &
             ErrorCode                        &
           )
! Reflects position of particle or puff centroid.

  Implicit None
  ! Argument list:
  Type(Coords_),     Intent(In)    :: Coords
  Type(Grids_),      Intent(In)    :: Grids
  Type(Domains_),    Intent(In)    :: Domains
  Type(ShortTime_),  Intent(In)    :: Time
  Logical,           Intent(In)    :: AnyTravelTime
  Type(ShortTime_),  Intent(In)    :: TravelTime
  Logical,           Intent(In)    :: Vel
  Type(Position_),   Intent(In)    :: OldPosition
  Type(Position_),   Intent(InOut) :: Position
  Real(Std),         Intent(InOut) :: U(3)
  Type(Units_),      Intent(InOut) :: Units
  Type(Mets_),       Intent(InOut) :: Mets
  Type(Flows_),      Intent(InOut) :: Flows
  Type(FlowMemory_), Intent(InOut) :: FlowMemory
  Logical,           Intent(Out)   :: Reflected
  Integer,           Intent(Out)   :: ErrorCode
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! Domains       :: Collection of domains.
  ! Time          :: Current time or time of the met to be used for any vertical coord
  !                  conversions (these must be the same unless we have a fixed met
  !                  case; for fixed met cases the value is ignored).
  ! AnyTravelTime :: Indicates that choices of flow module instances are not to be
  !                  restricted by travel time.
  ! TravelTime    :: Travel time.
  ! Vel           :: Indicates a dispersion model with velocity memory is being used.
  ! OldPosition   :: Old particle position.
  ! Position      :: Current particle position.
  ! U(3)          :: Current particle velocity. Defined only if Vel is true.
  ! Units         :: Collection of information on input/output unit numbers.
  ! Mets          :: Collection of met module instance states.
  ! Flows         :: Collection of flow module instance states.
  ! FlowMemory    :: Flow memory.
  ! Reflected     :: Indicates a reflection has occurred.
  ! ErrorCode     :: Error code. 0 = no error, 1 = no flow module for convert
  !                  attribute, 2 = no flow module for flow attribute.
  ! Locals:
  Type(ShortTime_) :: MetTime         ! Time of the met to be used for any vertical
                                      ! coord conversions.
  Integer          :: iFlowModConvert !} Indices of the flow module instance selected
  Integer          :: iFlowConvert    !} to perform any vertical coord conversions.
  Integer          :: iFlowMod        !] Indices of the flow module instance selected
  Integer          :: iFlow           !] to perform the reflection.
  Logical          :: Error           ! Indicates no suitable flow module instance was
                                      ! found.

  ErrorCode = 0

  If (Flows%C(1, 1)%P%FixedMet) Then
    MetTime = Flows%MetTime
  Else
    MetTime = Time
  End If

  ! Calculate best flow module instance for the convert attribute.
  Call WhichFlow(                               &
         Coords, Grids, Domains,                &
         Flows%iAttribs(A_Convert), 0, MetTime, &
         AnyTravelTime, TravelTime,             &
         iFlowModConvert, iFlowConvert, Error,  &
         Position,                              &
         Units, Mets, Flows,                    &
         FlowMemory                             &
       )
  If (Error) Then
    ErrorCode = 1
    Return
  End If

  ! Calculate best flow module instance for the flow attribute.
  Call WhichFlow(                            &
         Coords, Grids, Domains,             &
         Flows%iAttribs(A_Flow), 0, MetTime, &
         AnyTravelTime, TravelTime,          &
         iFlowMod, iFlow, Error,             &
         Position,                           &
         Units, Mets, Flows,                 &
         FlowMemory                          &
       )
  If (Error) Then
    ErrorCode = 2
    Return
  End If

  Call ReflectKnownFlow(                         &
         Flows,                                  &
         iFlowModConvert, iFlowConvert,          &
         iFlowMod, iFlow,                        &
         Coords, Grids,                          &
         MetTime, Vel, OldPosition, Position, U, &
         FlowMemory,                             &
         Reflected                               &
       )

End Subroutine Reflect

!-------------------------------------------------------------------------------------------------------------

Subroutine ReflectKnownFlow(                         &
             Flows,                                  &
             iFlowModConvert, iFlowConvert,          &
             iFlowMod, iFlow,                        &
             Coords, Grids,                          &
             MetTime, Vel, OldPosition, Position, U, &
             FlowMemory,                             &
             Reflected                               &
           )
! Reflects position of particle or puff centroid in situations where the flow module
! instances to be used for vertical coord conversions and for the reflection are
! known.

  Implicit None
  ! Argument list:
  Type(Flows_),      Intent(In)    :: Flows
  Integer,           Intent(In)    :: iFlowModConvert
  Integer,           Intent(In)    :: iFlowConvert
  Integer,           Intent(In)    :: iFlowMod
  Integer,           Intent(In)    :: iFlow
  Type(Coords_),     Intent(In)    :: Coords
  Type(Grids_),      Intent(In)    :: Grids
  Type(ShortTime_),  Intent(In)    :: MetTime
  Logical,           Intent(In)    :: Vel
  Type(Position_),   Intent(In)    :: OldPosition
  Type(Position_),   Intent(InOut) :: Position
  Real(Std),         Intent(InOut) :: U(3)
  Type(FlowMemory_), Intent(InOut) :: FlowMemory
  Logical,           Intent(Out)   :: Reflected
  ! Flows           :: Collection of flow module instance states.
  ! iFlowModConvert :} Indices of the flow module instance selected to perform any
  ! iFlowConvert    :} vertical coord conversions.
  ! iFlowMod        :] Indices of the flow module instance selected to perform the
  ! iFlow           :] reflection.
  ! Coords          :: Collection of coord systems.
  ! Grids           :: Collection of grids.
  ! MetTime         :: Time of the met to be used for any vertical coord conversions.
  ! Vel             :: Indicates a dispersion model with velocity memory is being
  !                    used.
  ! OldPosition     :: Old particle position.
  ! Position        :: Current particle position.
  ! U(3)            :: Current particle velocity. Defined only if Vel is true.
  ! FlowMemory      :: Flow memory.
  ! Reflected       :: Indicates a reflection has occurred.
  ! Locals:
  Integer :: nHCoords             !} Number and indices of the coords which are needed
  Integer :: iHCoords(MaxHCoords) !} for performing the reflection.
  Integer :: nZCoords             !}
  Integer :: iZCoords(MaxZCoords) !}
  Integer :: i                    ! Loop index.

  If (.false.) Then

  Else If (                                   &
    Flows%nPrototypeFlows /= 0 .and.               &
    iFlowMod == Flows%PrototypeFlows(1)%C%iFlowMod &
  ) Then

    Call PrototypeReflectCoordIndices(  &
           Flows%PrototypeFlows(iFlow), &
           nHCoords, iHCoords,     &
           nZCoords, iZCoords      &
         )
    Do i = 1, nHCoords
      Call ConvertToH(Coords, iHCoords(i), Position)
    End Do
    Do i = 1, nZCoords
      Call ConvertToZKnownFlow(             &
             Flows,                         &
             iFlowModConvert, iFlowConvert, &
             Coords, Grids,                 &
             iZCoords(i),                   &
             MetTime, Position,             &
             FlowMemory                     &
           )
    End Do
    Call PrototypeReflect(                      &
           Coords, Flows%PrototypeFlows(iFlow), &
           Vel, OldPosition, Position, U,  &
           Reflected                       &
         )

  Else If (                                   &
    Flows%nSingleSiteFlows /= 0 .and.               &
    iFlowMod == Flows%SingleSiteFlows(1)%C%iFlowMod &
  ) Then

    Call SingleSiteReflectCoordIndices(  &
           Flows%SingleSiteFlows(iFlow), &
           nHCoords, iHCoords,     &
           nZCoords, iZCoords      &
         )
    Do i = 1, nHCoords
      Call ConvertToH(Coords, iHCoords(i), Position)
    End Do
    Do i = 1, nZCoords
      Call ConvertToZKnownFlow(             &
             Flows,                         &
             iFlowModConvert, iFlowConvert, &
             Coords, Grids,                 &
             iZCoords(i),                   &
             MetTime, Position,             &
             FlowMemory                     &
           )
    End Do
    Call SingleSiteReflect(                      &
           Coords, Flows%SingleSiteFlows(iFlow), &
           Vel, OldPosition, Position, U,  &
           Reflected                       &
         )

  Else If (                                   &
    Flows%nNWPFlows /= 0 .and.               &
    iFlowMod == Flows%NWPFlows(1)%C%iFlowMod &
  ) Then

    Call NWPReflectCoordIndices(  &
           Flows%NWPFlows(iFlow), &
           nHCoords, iHCoords,     &
           nZCoords, iZCoords      &
         )
    Do i = 1, nHCoords
      Call ConvertToH(Coords, iHCoords(i), Position)
    End Do
    Do i = 1, nZCoords
      Call ConvertToZKnownFlow(             &
             Flows,                         &
             iFlowModConvert, iFlowConvert, &
             Coords, Grids,                 &
             iZCoords(i),                   &
             MetTime, Position,             &
             FlowMemory                     &
           )
    End Do
    Call NWPReflect(                      &
           Coords, Flows%NWPFlows(iFlow), &
           Vel, OldPosition, Position, U,  &
           Reflected                       &
         )

  Else If (                                   &
    Flows%nBuildingFlows /= 0 .and.               &
    iFlowMod == Flows%BuildingFlows(1)%C%iFlowMod &
  ) Then

    Call BuildingReflectCoordIndices(  &
           Flows%BuildingFlows(iFlow), &
           nHCoords, iHCoords,     &
           nZCoords, iZCoords      &
         )
    Do i = 1, nHCoords
      Call ConvertToH(Coords, iHCoords(i), Position)
    End Do
    Do i = 1, nZCoords
      Call ConvertToZKnownFlow(             &
             Flows,                         &
             iFlowModConvert, iFlowConvert, &
             Coords, Grids,                 &
             iZCoords(i),                   &
             MetTime, Position,             &
             FlowMemory                     &
           )
    End Do
    Call BuildingReflect(                      &
           Coords, Flows%BuildingFlows(iFlow), &
           Vel, OldPosition, Position, U,  &
           Reflected                       &
         )

  Else If (                                   &
    Flows%nLINCOMFlows /= 0 .and.               &
    iFlowMod == Flows%LINCOMFlows(1)%C%iFlowMod &
  ) Then

    Call LINCOMReflectCoordIndices(  &
           Flows%LINCOMFlows(iFlow), &
           nHCoords, iHCoords,     &
           nZCoords, iZCoords      &
         )
    Do i = 1, nHCoords
      Call ConvertToH(Coords, iHCoords(i), Position)
    End Do
    Do i = 1, nZCoords
      Call ConvertToZKnownFlow(             &
             Flows,                         &
             iFlowModConvert, iFlowConvert, &
             Coords, Grids,                 &
             iZCoords(i),                   &
             MetTime, Position,             &
             FlowMemory                     &
           )
    End Do
    Call LINCOMReflect(                      &
           Coords, Flows%LINCOMFlows(iFlow), &
           Vel, OldPosition, Position, U,  &
           Reflected                       &
         )


  Else
    Call Message('UNEXPECTED FATAL ERROR in ReflectKnownFlow: flow module instance not found', 4)
  End If

End Subroutine ReflectKnownFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcdZdZ(                                    &
             Coords, Grids, Domains,                    &
             iZCoord1, iZCoord2,                        &
             Time, AnyTravelTime, TravelTime, Position, &
             Units, Mets, Flows,                        &
             FlowMemory,                                &
             dZdZ,                                      &
             ErrorCode                                  &
           )
! Calculates dZCoord1/dZCoord2 for two vertical coord systems ZCoord1 and dZCoord2.

  Implicit None
  ! Argument list:
  Type(Coords_),     Intent(In)    :: Coords
  Type(Grids_),      Intent(In)    :: Grids
  Type(Domains_),    Intent(In)    :: Domains
  Integer,           Intent(In)    :: iZCoord1
  Integer,           Intent(In)    :: iZCoord2
  Type(ShortTime_),  Intent(In)    :: Time
  Logical,           Intent(In)    :: AnyTravelTime
  Type(ShortTime_),  Intent(In)    :: TravelTime
  Type(Position_),   Intent(InOut) :: Position
  Type(Units_),      Intent(InOut) :: Units
  Type(Mets_),       Intent(InOut) :: Mets
  Type(Flows_),      Intent(InOut) :: Flows
  Type(FlowMemory_), Intent(InOut) :: FlowMemory
  Real(Std),         Intent(Out)   :: dZdZ
  Integer,           Intent(Out)   :: ErrorCode
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! Domains       :: Collection of domains.
  ! iZCoord1      :} Indices in Coords of coord systems.
  ! iZCoord2      :}
  ! Time          :: Current time or time of the met to be used for calculating dZdZ
  !                  (these must be the same unless we have a fixed met case; for
  !                  fixed met cases the value is ignored).
  ! AnyTravelTime :: Indicates that choices of flow module instances are not to be
  !                  restricted by travel time.
  ! TravelTime    :: Travel time.
  ! Position      :: Position.
  ! Units         :: Collection of information on input/output unit numbers.
  ! Mets          :: Collection of met module instance states.
  ! Flows         :: Collection of flow module instance states.
  ! FlowMemory    :: Flow memory.
  ! dZdZ          :: dZCoord1/dZCoord2.
  ! ErrorCode     :: Error code. 0 = no error, 1 = no flow module for convert
  !                  attribute, 2 = no flow module for flow attribute.
  ! Locals:
  Type(Flow_)    :: Flow    ! Flow information.
  Type(Cloud_)   :: Cloud   !} Dummy variables for arguments in GetAttrib.
  Type(Rain_)    :: Rain    !}
  Type(Surface_) :: Surface !}
  Type(Soil_)    :: Soil    !}
  Type(Plant_)   :: Plant   !}
  Real(Std)      :: Temp    ! Temporary variable.

  ! $$ this and the next few routines should have an if test for fixed met and replace time by mettime where
  ! appropriate

  Select Case (Coords%ZCoords(iZCoord1)%CoordType)

    Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)

      Select Case (Coords%ZCoords(iZCoord2)%CoordType)

        Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)

          Call CalcdZdZ_ZBased(                        &
            Coords, Grids, Domains,                    &
            iZCoord1, iZCoord2,                        &
            Time, AnyTravelTime, TravelTime, Position, &
            Units, Mets, Flows,                        &
            FlowMemory,                                &
            dZdZ,                                      &
            ErrorCode                                  &
          )

        Case (Z_P, Z_PAsZ, Z_PAsEta)

          Call CalcdZdZ_ZBased(                        &
            Coords, Grids, Domains,                    &
            iZCoord1, Flows%iZCoordMagl,               &
            Time, AnyTravelTime, TravelTime, Position, &
            Units, Mets, Flows,                        &
            FlowMemory,                                &
            Temp,                                      &
            ErrorCode                                  &
          )
          dZdZ = Temp
          If (.not.Position%RhoValid) Then
            Call GetAttrib(                              &
              A_Flow,                                    &
              Coords, Grids, Domains,                    &
              .false., .false., .false.,                 &
              Time, AnyTravelTime, TravelTime, Position, &
              Units, Mets, Flows,                        &
              FlowMemory,                                &
              Flow, Cloud, Rain,                         &
              Surface, Soil, Plant, ErrorCode            &
            )
          End If
          dZdZ = dZdZ / (Position%Rho * Gravity)
          Call CalcdZdZ_PBased(                        &
            Coords, Grids, Domains,                    &
            Flows%iZCoordPa, iZCoord2,                 &
            Time, AnyTravelTime, TravelTime, Position, &
            Units, Mets, Flows,                        &
            FlowMemory,                                &
            Temp,                                      &
            ErrorCode                                  &
          )
          dZdZ = dZdZ * Temp

        Case Default

          Call Message('UNEXPECTED FATAL ERROR in CalcdZdZ', 4)

      End Select

    Case (Z_P, Z_PAsZ, Z_PAsEta)

      Select Case (Coords%ZCoords(iZCoord2)%CoordType)

        Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)

          Call CalcdZdZ_PBased(                         &
             Coords, Grids, Domains,                    &
             iZCoord1, Flows%iZCoordPa,                 &
             Time, AnyTravelTime, TravelTime, Position, &
             Units, Mets, Flows,                        &
             FlowMemory,                                &
             Temp,                                      &
             ErrorCode                                  &
           )
          dZdZ = Temp
          If (.not.Position%RhoValid) Then
            Call GetAttrib(                              &
              A_Flow,                                    &
              Coords, Grids, Domains,                    &
              .false., .false., .false.,                 &
              Time, AnyTravelTime, TravelTime, Position, &
              Units, Mets, Flows,                        &
              FlowMemory,                                &
              Flow, Cloud, Rain,                         &
              Surface, Soil, Plant, ErrorCode            &
            )
          End If
          dZdZ = dZdZ * (Position%Rho * Gravity)
          Call CalcdZdZ_ZBased(                        &
            Coords, Grids, Domains,                    &
            Flows%iZCoordMagl, iZCoord2,               &
            Time, AnyTravelTime, TravelTime, Position, &
            Units, Mets, Flows,                        &
            FlowMemory,                                &
            Temp,                                      &
            ErrorCode                                  &
          )
          dZdZ = dZdZ * Temp

        Case (Z_P, Z_PAsZ, Z_PAsEta)

          Call CalcdZdZ_PBased(                         &
             Coords, Grids, Domains,                    &
             iZCoord1, iZCoord2,                        &
             Time, AnyTravelTime, TravelTime, Position, &
             Units, Mets, Flows,                        &
             FlowMemory,                                &
             dZdZ,                                      &
             ErrorCode                                  &
           )

        Case Default

          Call Message('UNEXPECTED FATAL ERROR in CalcdZdZ', 4)

      End Select

    Case Default

      Call Message('UNEXPECTED FATAL ERROR in CalcdZdZ', 4)

  End Select

End Subroutine CalcdZdZ

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcdZdZ_ZBased(                             &
             Coords, Grids, Domains,                    &
             iZCoord1, iZCoord2,                        &
             Time, AnyTravelTime, TravelTime, Position, &
             Units, Mets, Flows,                        &
             FlowMemory,                                &
             dZdZ,                                      &
             ErrorCode                                  &
           )
! Calculates dZCoord1/dZCoord2 for two height-based vertical coord systems ZCoord1 and
! dZCoord2.

  Implicit None
  ! Argument list:
  Type(Coords_),     Intent(In),   Target :: Coords
  Type(Grids_),      Intent(In)           :: Grids
  Type(Domains_),    Intent(In)           :: Domains
  Integer,           Intent(In)           :: iZCoord1
  Integer,           Intent(In)           :: iZCoord2
  Type(ShortTime_),  Intent(In)           :: Time
  Logical,           Intent(In)           :: AnyTravelTime
  Type(ShortTime_),  Intent(In)           :: TravelTime
  Type(Position_),   Intent(InOut)        :: Position
  Type(Units_),      Intent(InOut)        :: Units
  Type(Mets_),       Intent(InOut)        :: Mets
  Type(Flows_),      Intent(InOut)        :: Flows
  Type(FlowMemory_), Intent(InOut)        :: FlowMemory
  Real(Std),         Intent(Out)          :: dZdZ
  Integer,           Intent(Out)          :: ErrorCode
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! Domains       :: Collection of domains.
  ! iZCoord1      :} Indices in Coords of coord systems.
  ! iZCoord2      :}
  ! Time          :: Current time or time of the met to be used for calculating dZdZ
  !                  (these must be the same unless we have a fixed met case; for
  !                  fixed met cases the value is ignored).
  ! AnyTravelTime :: Indicates that choices of flow module instances are not to be
  !                  restricted by travel time.
  ! TravelTime    :: Travel time.
  ! Position      :: Position.
  ! Units         :: Collection of information on input/output unit numbers.
  ! Mets          :: Collection of met module instance states.
  ! Flows         :: Collection of flow module instance states.
  ! FlowMemory    :: Flow memory.
  ! dZdZ          :: dZCoord1/dZCoord2.
  ! ErrorCode     :: Error code. 0 = no error, 1 = no flow module for convert
  !                  attribute, 2 = no flow module for flow attribute.
  ! Locals:
  Type(ZCoord_), Pointer :: ZCoord1 !} Abbreviations for vertical coord systems.
  Type(ZCoord_), Pointer :: ZCoord2 !}
  Type(Flow_)            :: Flow    ! Flow information.
  Type(Cloud_)           :: Cloud   !} Dummy variables for arguments in GetAttrib.
  Type(Rain_)            :: Rain    !}
  Type(Surface_)         :: Surface !}
  Type(Soil_)            :: Soil !}
  Type(Plant_)           :: Plant !}

  ZCoord1 => Coords%ZCoords(iZCoord1)
  ZCoord2 => Coords%ZCoords(iZCoord2)

  If (ZCoord1%CoordType == Z_ZAsEta) Then

    Call ConvertToZ(                             &
      Coords, Grids, Domains,                    &
      iZCoord1,                                  &
      Time, AnyTravelTime, TravelTime, Position, &
      Units, Mets, Flows,                        &
      FlowMemory,                                &
      ErrorCode                                  &
    )
    If (.not.Position%TopogValid) Then
      Call GetAttrib(                              &
        A_Flow,                                    &
        Coords, Grids, Domains,                    &
        .false., .false., .false.,                 &
        Time, AnyTravelTime, TravelTime, Position, &
        Units, Mets, Flows,                        &
        FlowMemory,                                &
        Flow, Cloud, Rain,                         &
        Surface, Soil, Plant, ErrorCode            &
      )
    End If
    dZdZ = CalcdZdZZBased(ZCoord1, ZCoord2, Position%Z(iZCoord1), Position%Topog)

  Else If (ZCoord2%CoordType == Z_ZAsEta) Then

    Call ConvertToZ(                             &
      Coords, Grids, Domains,                    &
      iZCoord2,                                  &
      Time, AnyTravelTime, TravelTime, Position, &
      Units, Mets, Flows,                        &
      FlowMemory,                                &
      ErrorCode                                  &
    )
    If (.not.Position%TopogValid) Then
      Call GetAttrib(                              &
        A_Flow,                                    &
        Coords, Grids, Domains,                    &
        .false., .false., .false.,                 &
        Time, AnyTravelTime, TravelTime, Position, &
        Units, Mets, Flows,                        &
        FlowMemory,                                &
        Flow, Cloud, Rain,                         &
        Surface, Soil, Plant, ErrorCode            &
      )
    End If
    dZdZ = 1.0 / CalcdZdZZBased(ZCoord2, ZCoord1, Position%Z(iZCoord2), Position%Topog)

  Else

    dZdZ = CalcdZdZZBased(ZCoord1, ZCoord2)

  End If

End Subroutine CalcdZdZ_ZBased

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcdZdZ_PBased(                             &
             Coords, Grids, Domains,                    &
             iZCoord1, iZCoord2,                        &
             Time, AnyTravelTime, TravelTime, Position, &
             Units, Mets, Flows,                        &
             FlowMemory,                                &
             dZdZ,                                      &
             ErrorCode                                  &
           )
! Calculates dZCoord1/dZCoord2 for two pressure-based vertical coord systems ZCoord1
! and dZCoord2.

  Implicit None
  ! Argument list:
  Type(Coords_),     Intent(In),   Target :: Coords
  Type(Grids_),      Intent(In)           :: Grids
  Type(Domains_),    Intent(In)           :: Domains
  Integer,           Intent(In)           :: iZCoord1
  Integer,           Intent(In)           :: iZCoord2
  Type(ShortTime_),  Intent(In)           :: Time
  Logical,           Intent(In)           :: AnyTravelTime
  Type(ShortTime_),  Intent(In)           :: TravelTime
  Type(Position_),   Intent(InOut)        :: Position
  Type(Units_),      Intent(InOut)        :: Units
  Type(Mets_),       Intent(InOut)        :: Mets
  Type(Flows_),      Intent(InOut)        :: Flows
  Type(FlowMemory_), Intent(InOut)        :: FlowMemory
  Real(Std),         Intent(Out)          :: dZdZ
  Integer,           Intent(Out)          :: ErrorCode
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! Domains       :: Collection of domains.
  ! iZCoord1      :} Indices in Coords of coord systems.
  ! iZCoord2      :}
  ! Time          :: Current time or time of the met to be used for calculating dZdZ
  !                  (these must be the same unless we have a fixed met case; for
  !                  fixed met cases the value is ignored).
  ! AnyTravelTime :: Indicates that choices of flow module instances are not to be
  !                  restricted by travel time.
  ! TravelTime    :: Travel time.
  ! Position      :: Position.
  ! Units         :: Collection of information on input/output unit numbers.
  ! Mets          :: Collection of met module instance states.
  ! Flows         :: Collection of flow module instance states.
  ! FlowMemory    :: Flow memory.
  ! dZdZ          :: dZCoord1/dZCoord2.
  ! ErrorCode     :: Error code. 0 = no error, 1 = no flow module for convert
  !                  attribute, 2 = no flow module for flow attribute.
  ! Locals:
  Type(ZCoord_), Pointer :: ZCoord1 !} Abbreviations for vertical coord systems.
  Type(ZCoord_), Pointer :: ZCoord2 !}
  Type(Flow_)            :: Flow    ! Flow information.
  Type(Cloud_)           :: Cloud   !} Dummy variables for arguments in GetAttrib.
  Type(Rain_)            :: Rain    !}
  Type(Surface_)         :: Surface !}
  Type(Soil_)            :: Soil    !}
  Type(Plant_)           :: Plant   !}

  ZCoord1 => Coords%ZCoords(iZCoord1)
  ZCoord2 => Coords%ZCoords(iZCoord2)

  If (ZCoord1%CoordType == Z_ZAsEta) Then

    Call ConvertToZ(                             &
      Coords, Grids, Domains,                    &
      iZCoord1,                                  &
      Time, AnyTravelTime, TravelTime, Position, &
      Units, Mets, Flows,                        &
      FlowMemory,                                &
      ErrorCode                                  &
    )
    If (.not.Position%PSValid) Then
      Call GetAttrib(                              &
        A_Flow,                                    &
        Coords, Grids, Domains,                    &
        .false., .false., .false.,                 &
        Time, AnyTravelTime, TravelTime, Position, &
        Units, Mets, Flows,                        &
        FlowMemory,                                &
        Flow, Cloud, Rain,                         &
        Surface, Soil, Plant, ErrorCode            &
      )
    End If
    dZdZ = CalcdZdZPBased(ZCoord1, ZCoord2, Position%Z(iZCoord1), Position%PS)

  Else If (ZCoord2%CoordType == Z_ZAsEta) Then

    Call ConvertToZ(                             &
      Coords, Grids, Domains,                    &
      iZCoord2,                                  &
      Time, AnyTravelTime, TravelTime, Position, &
      Units, Mets, Flows,                        &
      FlowMemory,                                &
      ErrorCode                                  &
    )
    If (.not.Position%PSValid) Then
      Call GetAttrib(                              &
        A_Flow,                                    &
        Coords, Grids, Domains,                    &
        .false., .false., .false.,                 &
        Time, AnyTravelTime, TravelTime, Position, &
        Units, Mets, Flows,                        &
        FlowMemory,                                &
        Flow, Cloud, Rain,                         &
        Surface, Soil, Plant, ErrorCode            &
      )
    End If
    dZdZ = 1.0 / CalcdZdZPBased(ZCoord2, ZCoord1, Position%Z(iZCoord2), Position%PS)

  Else If (ZCoord1%CoordType == Z_PAsZ .or. ZCoord2%CoordType == Z_PAsZ) Then

    Call ConvertToZ(                             &
      Coords, Grids, Domains,                    &
      iZCoord1,                                  &
      Time, AnyTravelTime, TravelTime, Position, &
      Units, Mets, Flows,                        &
      FlowMemory,                                &
      ErrorCode                                  &
    )
    dZdZ = CalcdZdZPBased(ZCoord1, ZCoord2, Position%Z(iZCoord1))

  Else

    dZdZ = CalcdZdZPBased(ZCoord1, ZCoord2)

  End If

End Subroutine CalcdZdZ_PBased

!-------------------------------------------------------------------------------------------------------------

Subroutine XInTheDomain(                                &
             Coords, Grids, Domains,                    &
             iDomain,                                   &
             Time, AnyTravelTime, TravelTime, Position, &
             Units, Mets, Flows,                        &
             FlowMemory,                                &
             In, ErrorCode                              &
           )
! Finds out whether a point is in a domain (other than a flow module instance's
! domain).

! Note that a point is defined to lie within a domain if and only if it is in the
! domain when this is evaluated with any required coordinate conversions being
! performed by the preferred flow module instance for coordinate conversions (as
! determined by WhichFlow - this determination may of course make use of Time and
! TravelTime).

  Implicit None
  ! Argument list:
  Type(Coords_),     Intent(In)           :: Coords
  Type(Grids_),      Intent(In)           :: Grids
  Type(Domains_),    Intent(In),   Target :: Domains
  Integer,           Intent(In)           :: iDomain
  Type(ShortTime_),  Intent(In)           :: Time
  Logical,           Intent(In)           :: AnyTravelTime
  Type(ShortTime_),  Intent(In)           :: TravelTime
  Type(Position_),   Intent(InOut)        :: Position
  Type(Units_),      Intent(InOut)        :: Units
  Type(Mets_),       Intent(InOut)        :: Mets
  Type(Flows_),      Intent(InOut)        :: Flows
  Type(FlowMemory_), Intent(InOut)        :: FlowMemory
  Logical,           Intent(Out)          :: In
  Integer,           Intent(Out)          :: ErrorCode
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! Domains       :: Collection of domains.
  ! iDomain       :: Index of domain.
  ! Time          :: Current time or time of the met to be used for performing the
  !                  vertical coord conversion (these must be the same unless we have
  !                  a fixed met case; for fixed met cases the value is ignored).
  ! AnyTravelTime :: Indicates that choices of flow module instances are not to be
  !                  restricted by travel time.
  ! TravelTime    :: Travel time.
  ! Position      :: Position.
  ! Units         :: Collection of information on input/output unit numbers.
  ! Mets          :: Collection of met module instance states.
  ! Flows         :: Collection of flow module instance states.
  ! FlowMemory    :: Flow memory.
  ! In            :: Indicates whether the position is in the domain.
  ! ErrorCode     :: Error code. 0 = no error, 1 = no flow module for convert
  !                  attribute. In case of error, In is set to false.
  ! Locals:
  Type(Domain_), Pointer :: Domain ! Abbreviation for the domain.

  Domain => Domains%Domains(iDomain)

  If (.not. Domain%XUnbounded .or. .not. Domain%YUnbounded) Then
    Call ConvertToH(Coords, Domain%iHCoord, Position)
  End If

  In = HInDomain(Position, Domain, Coords)

  If (In) Then

    If (.not. Domain%ZUnbounded) Then
      Call ConvertToZ(                                  &
             Coords, Grids, Domains,                    &
             Domain%iZCoord,                            &
             Time, AnyTravelTime, TravelTime, Position, &
             Units, Mets, Flows,                        &
             FlowMemory,                                &
             ErrorCode                                  &
           )
      If (ErrorCode > 0) Then
        In = .false.
        Return
      End If
    End If

    In = ZInDomain(Position, Domain, Coords)

  End If

End Subroutine XInTheDomain

!-------------------------------------------------------------------------------------------------------------

Recursive Subroutine WhichFlow(                   &
                       Coords, Grids, Domains,    &
                       iAttrib, iSubset, MetTime, &
                       AnyTravelTime, TravelTime, &
                       iFlowMod, iFlow, Error,    &
                       Position,                  &
                       Units, Mets, Flows,        &
                       FlowMemory                 &
                     )
! Determines the preferred flow module instance to use for a given point and attribute.

! For attributes other than 'Convert', the preferred flow module instance for the Convert attribute should be
! determined first.

! The preferred flow module instance to use for a given point and attribute is defined as the highest priority
! (according to the priority order for the attribute) flow module instance which is 'suitable' for the point
! and attribute.

! A flow module instance is suitable for a point and attribute if and only if
!   (i)  the flow module instance is valid for the attribute at met time, and
!   (ii) the point, time and travel time are in the flow module instance's domain.
! Note 'suitability' depends on the time, travel time and met time as well as the point itself. [For non-fixed
! met, time = met time. For fixed met, time /= met time and the flow module instance's validity for the
! attribute depends on the met time but is, like the met time, unchanging (within a single case). The met time
! dependence arises because it may affect which flow module instances are valid for the attribute and because
! it is used for any coordinate conversions needed.] For non-fixed met no direct check is made on whether the
! time lies in the flow module instance's domain. This is because, for non-fixed met, the flow module validity
! is already restricted by the time limits of the domain.
!
! $$ Currently no check made for fixed met either. This assumes time unbounded domains with fixed met, but we
! could relax this in the future. For fixed met the domain should restrict the time and should be tested here,
! but the met time would not have to lie in the domain. See also CommonFlow.F90.

! 'A point being in a flow module instance's domain' is a slightly different concept from 'a point being in a
! domain'. The correct concept is described here.

! For flow module instances which are valid for the attribute in question but not for the Convert attribute at
! the met time, a point is in the flow module instance's domain if the point is in the domain when this is
! evaluated with any required coordinate conversions being performed by the preferred flow module instance for
! the given point and the Convert attribute.

! For flow module instances which are valid for the attribute in question and for the Convert attribute at the
! met time, the situation is more complex. Let A_1, ..., A_n be the flow module instances which are valid for
! the attribute in question and for the Convert attribute at the met time in order of priority (according to
! the priority order for the Convert attribute) with A_1 the highest priority. Then a point is defined to lie
! within the domain of A_1 if and only if it is within the domain according to A_1's conversion routines.
! Inductively, a point is defined to lie within the domain of A_m if either (i) the preferred flow module
! instance for coordinate conversions is A_p with p < m and the point lies within the domain according to
! A_p's conversion routines, or (ii) none of the A_i, i < m, is the preferred flow module instance for
! coordinate conversions and the point lies within the domain according to A_m's conversion routines. This is
! well defined because we can determine if any of A_i, i < m, is the preferred flow module instance for
! coordinate conversions once we know whether or not the point lies within the domains of A_i, i < m.

! For example, suppose there are just two flow module instances which are valid for the convert attribute, A
! and B say, with A a higher priority (according to the priority order for the Convert attribute) than B.
! Suppose A determines that the point lies outside the A domain, but B determines that it lies within both the
! A and B domains. B then becomes the preferred flow module instance for cooordinate conversions, but cannot
! then overrule A and say the point is in A (otherwise we would want to make A the preferred flow module
! instance for cooordinate conversions).

  Implicit None
  ! Argument list:
  Type(Coords_),     Intent(In)           :: Coords
  Type(Grids_),      Intent(In)           :: Grids
  Type(Domains_),    Intent(In),   Target :: Domains
  Integer,           Intent(In)           :: iAttrib
  Integer,           Intent(In)           :: iSubset
  Type(ShortTime_),  Intent(In)           :: MetTime
  Logical,           Intent(In)           :: AnyTravelTime
  Type(ShortTime_),  Intent(In)           :: TravelTime
  Integer,           Intent(Out)          :: iFlowMod
  Integer,           Intent(Out)          :: iFlow
  Logical,           Intent(Out)          :: Error
  Type(Position_),   Intent(InOut)        :: Position
  Type(Units_),      Intent(InOut)        :: Units
  Type(Mets_),       Intent(InOut)        :: Mets
  Type(Flows_),      Intent(InOut)        :: Flows
  Type(FlowMemory_), Intent(InOut)        :: FlowMemory
  ! Coords          :: Collection of coord systems.
  ! Grids           :: Collection of grids.
  ! Domains         :: Collection of domains.
  ! iAttrib         :: Index of attribute.
  ! iSubset         :: Index of subset giving a restricted choice of flow module
  !                    instances. 0 indicates no restriction.
  ! MetTime         :: Time of the met to be used for any coordinate conversions needed in determining the
  !                    preferred flow module instance for the attribute.
  ! AnyTravelTime   :: Indicates that choices of flow module instances are not to be
  !                    restricted by travel time.
  ! TravelTime      :: Travel time.
  ! iFlowMod        :} Indices of the preferred flow module instance for the attribute.
  ! iFlow           :}
  ! Error           :: Indicates no suitable flow module instance was found.
  ! Position        :: Position.
  ! Units           :: Collection of information on input/output unit numbers.
  ! Mets            :: Collection of met module instance states.
  ! Flows           :: Collection of flow module instance states.
  ! FlowMemory      :: FlowMemory.
  ! Locals:
  Type(Domain_), Pointer :: Domain             ! Abbreviation for domain.
  Integer                :: iHCoord            !} Indices of coord systems used by the domain (if bounded).
  Integer                :: iZCoord            !}
  Integer                :: i                  !] Loop indices.
  Integer                :: j                  !]
  Integer                :: k                  !]
  Logical                :: TravelTimeInDomain ! Indicates the domain is not excluded as a result of the
                                               ! travel time.
  Integer                :: iOrder             ! Index of ordered subset of flow module instances.

  ! For attributes other than convert, the preferred flow for the convert attribute should have been
  ! determined already.
# ifdef ExtraChecks
    If (iAttrib /= Flows%iAttribs(A_Convert) .and. .not. FlowMemory%iValid(Flows%iAttribs(A_Convert))) Then
      Call Message('UNEXPECTED FATAL ERROR in WhichFlow', 4)
    End If
# endif

  Error = .false.

  If (FlowMemory%iValid(iAttrib)) Then
    iFlowMod = FlowMemory%iFlowMod(iAttrib)
    iFlow    = FlowMemory%iFlow(iAttrib)
    Return
  End If

  iOrder = Flows%Attribs(iAttrib)%iOrder

  Do k = 1, Flows%Orders(iOrder)%nFlows

    i = Flows%Orders(iOrder)%iFlowMods(k)
    j = Flows%Orders(iOrder)%iFlows(k)

    ! Flow module instance valid (if not due for update)?

    If (.not. Flows%C(i, j)%P%DueForUpdate) Then
      If (.not. Flows%C(i, j)%P%ValidAttribParams(Flows%Attribs(iAttrib)%iAttribParam)) Cycle
    End If

    ! Flow module instance included in subset (if subset exists)?

    If (iSubset /= 0) Then
      If (.not. Flows%Subsets(iSubset)%Included(i, j)) Cycle
    End If


    If (.not. FlowMemory%InValid(i, j)) Then

      Domain => Domains%Domains(Flows%C(i, j)%P%iDomain)

      ! Travel time OK for flow module instance?

      If (AnyTravelTime) Then
        TravelTimeInDomain = .true.
      Else
        TravelTimeInDomain = SInDomain(TravelTime, Domain)
      End If

      If (.not.TravelTimeInDomain) Then
        FlowMemory%In     (i, j) = .false.
        FlowMemory%InValid(i, j) = .true.
        Cycle
      End If

      ! Horizontal position OK for flow module instance?

      If (.not. Domain%XUnbounded .or. .not. Domain%YUnbounded) Then
        iHCoord = Domain%iHCoord
        If (.not. Position%XYValid(iHCoord)) Then
          Call ConvertToH(Coords, iHCoord, Position)
        End If
        If (.not. HInDomain(Position, Domain, Coords)) Then
          FlowMemory%In     (i, j) = .false.
          FlowMemory%InValid(i, j) = .true.
          Cycle
        End If
      End If


      ! Case 1: We're trying to find a flow module instance to do a vertical coord conversion.
      If (iAttrib == Flows%iAttribs(A_Convert)) Then

        ! Vertical position OK for flow module instance? Flow module instance updated successfully (if due for
        ! update and needs updating to do vertical position test)?

        If (.not. Domain%ZUnbounded) Then

          iZCoord = Domain%iZCoord

          If (Position%ZValid(iZCoord)) Then

            If (.not. ZInDomain(Position, Domain, Coords)) Then
              FlowMemory%In     (i, j) = .false.
              FlowMemory%InValid(i, j) = .true.
              Cycle
            End If

          Else

            ! If due for update, try to update module.
            If (Flows%C(i, j)%P%DueForUpdate) Then

              Call SetFlowLock(i,j)

              ! Update flow module instance.
              Call UpdateFlow(                      &
                     Coords, Grids, Domains,        &
                     Flows%iCase, Flows%iMetCase,   &
                     i, j,                          &
                     ShortTime2Time(Flows%MetTime), & ! $$ rationalise times
                     Mets, Flows,                   &
                     Units                          &
                   )

              Call UnsetFlowLock(i,j)

              If (.not. Flows%C(i, j)%P%ValidAttribParams(Flows%Attribs(iAttrib)%iAttribParam)) Cycle

            End If

            Call ConvertToZKnownFlow( &
                   Flows,             &
                   i, j,              &
                   Coords, Grids,     &
                   iZCoord,           &
                   MetTime, Position, &
                   FlowMemory         &
                 )

            If (.not. ZInDomain(Position, Domain, Coords)) Then

              ! The Z coord calculated above is invalid if (i) it is calculated using a flow
              ! module instance which is not known (at the time of calculation) to be the
              ! correct flow module instance for the coord conversion, and (ii) the point
              ! isn't in the domain of the flow module instance used for the conversion. In
              ! this situation the flow memory of the flow module instance which attempted
              ! the conversion may also be invalid. Note however that FlowMemory%In(i, j)
              ! and FlowMemory%InValid(i, j) are still correct (see discussion above).
              Position%ZValid(iZCoord) = .false.
              Position%TopogValid      = .false.
              Position%PSValid         = .false.
              Position%RhoValid        = .false.
              Call ResetFlowMemorySpecificFlow(Flows, i, j, FlowMemory)

              FlowMemory%In     (i, j) = .false.
              FlowMemory%InValid(i, j) = .true.
              Cycle

            End If

          End If

        End If

      ! Case 2: We're trying to find a flow module instance for purposes other than vertical coord conversion.
      Else

        ! Vertical position OK for flow module instance?

        If (.not. Domain%ZUnbounded) Then
          iZCoord = Domain%iZCoord
          If (.not. Position%ZValid(iZCoord)) Then
            Call ConvertToZKnownFlow(                              &
                   Flows,                                          &
                   FlowMemory%iFlowMod(Flows%iAttribs(A_Convert)), &
                   FlowMemory%iFlow(Flows%iAttribs(A_Convert)),    &
                   Coords, Grids,                                  &
                   iZCoord,                                        &
                   MetTime, Position,                              &
                   FlowMemory                                      &
                 )
          End If
          If (.not. ZInDomain(Position, Domain, Coords)) Then
            FlowMemory%In     (i, j) = .false.
            FlowMemory%InValid(i, j) = .true.
            Cycle
          End If
        End If

      End If

      ! Flow module instance updated successfully (if due for update)?

      If (Flows%C(i, j)%P%DueForUpdate) Then

        Call SetFlowLock(i,j)

        ! Update flow module instance.
        Call UpdateFlow(                      &
               Coords, Grids, Domains,        &
               Flows%iCase, Flows%iMetCase,   &
               i, j,                          &
               ShortTime2Time(Flows%MetTime), & ! $$ rationalise times
               Mets, Flows,                   &
               Units                          &
             )

        Call UnsetFlowLock(i,j)

        If (.not. Flows%C(i, j)%P%ValidAttribParams(Flows%Attribs(iAttrib)%iAttribParam)) Cycle

      End If

      FlowMemory%In     (i, j) = .true.
      FlowMemory%InValid(i, j) = .true.

    End If

    If (FlowMemory%In(i, j)) Then
      iFlowMod = i
      iFlow    = j
      FlowMemory%iFlowMod(iAttrib) = iFlowMod
      FlowMemory%iFlow(iAttrib)    = iFlow
      FlowMemory%iValid(iAttrib)   = .true.
      Return
    End If

  End Do

  Error = .true.

End Subroutine WhichFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine ResetMetsFlows(Mets, Flows)
! Resets the state of the met and flow module instances for a new realisation.

  Implicit None
  ! Argument list:
  Type(Mets_),  Intent(InOut) :: Mets  !  Collection of met module instance states.
  Type(Flows_), Intent(InOut) :: Flows ! Collection of flow module instance states.

  Call ResetMets(Mets)
  Call ResetFlows(Flows)

End Subroutine ResetMetsFlows

!-------------------------------------------------------------------------------------------------------------

Subroutine ResetFlows(Flows)
! Resets the state of the flow module instances for a new realisation.

  Implicit None
  ! Argument list:
  Type(Flows_), Intent(InOut) :: Flows ! Collection of flow module instance states.
  ! Locals:
  Integer :: iFlowMod ! Flow module index.
  Integer :: iFlow    ! Flow module instance index.

  Do iFlowMod = 1, Flows%nFlowMods
  Do iFlow    = 1, Flows%nFlows(iFlowMod)

    Flows%C(iFlowMod, iFlow)%P%Valid                = .false.
    Flows%C(iFlowMod, iFlow)%P%ValidAttribParams(:) = .false.
    Flows%C(iFlowMod, iFlow)%P%ValidityUnimprovable = .false.
    Flows%C(iFlowMod, iFlow)%P%TValid               = InfPastTime()
    Flows%C(iFlowMod, iFlow)%P%DueForUpdate         = .false.

  End Do
  End Do

End Subroutine ResetFlows

!-------------------------------------------------------------------------------------------------------------

Subroutine CreateFlowMetLocks()
! Create locks for Flow Modules

  Implicit None

! Indices of Flow module
  Integer :: i, j

  Do i = 1, MaxFlowMods
    Do j = 1, MaxFlowsPerMod
!$      Call omp_init_lock(FlowLocks(i,j))
!$      Call omp_init_lock(MetLocks(i,j))
    End Do
  End Do

End Subroutine CreateFlowMetLocks

!-------------------------------------------------------------------------------------------------------------

Subroutine DestroyFlowMetLocks()
! Destroy locks for Flow Modules

  Implicit None

! Indices of Flow module
  Integer :: i, j

  Do i = 1, MaxFlowMods
    Do j = 1, MaxFlowsPerMod
!$      Call omp_destroy_lock(FlowLocks(i,j))
!$      Call omp_destroy_lock(MetLocks(i,j))
    End Do
  End Do

End Subroutine DestroyFlowMetLocks

!-------------------------------------------------------------------------------------------------------------

Subroutine SetFlowLock(i,j)
! Set lock for a particular Flow Module

  Implicit None

! Indices of Flow module
  Integer, Intent(In) :: i
  Integer, Intent(In) :: j

!$  Call omp_set_lock(FlowLocks(i,j))

End Subroutine SetFlowLock

!-------------------------------------------------------------------------------------------------------------

Subroutine UnsetFlowLock(i,j)
! Unset lock for a particular Flow Module

  Implicit None

! Indices of Flow module
  Integer, Intent(In) :: i
  Integer, Intent(In) :: j

!$  Call omp_unset_lock(FlowLocks(i,j))

End Subroutine UnsetFlowLock

!-------------------------------------------------------------------------------------------------------------

Subroutine SetMetLock(i,j)
! Set lock for a particular Met Module

  Implicit None

! Indices of Met module
  Integer, Intent(In) :: i
  Integer, Intent(In) :: j

!$  Call omp_set_lock(MetLocks(i,j))

End Subroutine SetMetLock

!-------------------------------------------------------------------------------------------------------------

Subroutine UnsetMetLock(i,j)
! Unset lock for a particular Flow Module

  Implicit None

! Indices of Met module
  Integer, Intent(In) :: i
  Integer, Intent(In) :: j

!$  Call omp_unset_lock(MetLocks(i,j))

End Subroutine UnsetMetLock

!-------------------------------------------------------------------------------------------------------------

End Module FlowsModule

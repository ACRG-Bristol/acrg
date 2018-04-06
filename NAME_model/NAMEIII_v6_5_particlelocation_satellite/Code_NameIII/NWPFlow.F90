! Module:  NWP Flow Module

Module NWPFlowModule

! This is a flow module designed to interpolate gridded NWP data and add turbulence
! information.

! It uses input from a single instance of the NWP met module.
!
! It does not use any input from other flow modules (and hence uses no update subset).
!
! It uses the topography and surface pressure information supplied by the NWP met
! module.
!
! It has attributes Update, Convert, Flow, Cloud and Rain.
!
! The coords stored in NWPFlow%C are
! Horizontal: 1: Coord system for horizontal grids used for met data.
!             2: Polar stereographic system at north pole of coord system 1. $$
!             3: Polar stereographic system at south pole of coord system 1. $$
! Vertical:   1: Coord system of vertical grids used for met data.
!             2: m agl.
!             3: Pressure (Pa).
! GetNWPFlow returns flow information in horizontal coord system 1, except near the
! poles where systems 2 or 3 are used as appropriate, and in vertical coord system 2.
!
! No grids are stored in NWPFlow%C.
!
! There are no restrictions on the domain. $$ NWPGrid? Upper boundary?
!
! It only works in a calendar time frame.
!
! The status of the optional flow module features is as follows (where ???? here =
! NWP):
!     ????FlowMemory_ data type           - used
!     SetUpCoordsEtc_????Flow routine     - used
!     SetUp????Flow_MetsCoordsEtc routine - used
!     Reset????FlowMemory routine         - used
!     ????ConvertCoordIndices routine     - used
!     ????ConvertToZ routine              - used (with flow memory argument)
!     ????FlowCoordIndices routine        - used
!     Get????Flow routine                 - used (with flow memory argument)
!     ????CloudCoordIndices routine       - used
!     Get????Cloud routine                - used (with flow memory argument)
!     ????RainCoordIndices routine        - used
!     Get????Rain routine                 - used (with flow memory argument)
!     ????SurfaceCoordIndices routine     - used
!     Get????Surface routine              - used (with flow memory argument)
!     ????SoilCoordIndices routine        - used
!     Get????Soil routine                 - used (with flow memory argument)
!     ????PlantCoordIndices routine       - used
!     Get????Plant routine                - used (with flow memory argument)
!     ????ReflectCoordIndices             - used
!     ????Reflect routine                 - used.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowAndFlowProfileModule
Use NWPMetModule
Use AncillaryMetModule
Use MetsModule
Use CommonFlowModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: NWPFlow_                   ! The state of an NWP flow module instance.
Public  :: NWPFlowMemory_             ! Information related to a particular space-time
                                      ! point of interest which an NWP flow module
                                      ! instance wishes to record.
Public  :: InitNWPFlow                ! Initialises the state of an NWP flow module
                                      ! instance.
Public  :: SetUpCoordsEtc_NWPFlow     ! Sets up Coords and Grids by adding any extra
                                      ! coords and grids which NWPFlow wants to
                                      ! define.
Public  :: SetUpNWPFlow_MetsCoordsEtc ! Sets up NWPFlow using information from
                                      ! EtaDefns, Coords, Grids and Mets.
Public  :: PrepareForUpdateNWPFlow    ! Prepares for updating an instance of the NWP flow module.
Public  :: NWPFlowReqs                ! Specifies what information the flow module
                                      ! instance wants from the other flow module
                                      ! instances.
Public  :: UpdateNWPFlow              ! Updates an instance of the NWP flow module.
Public  :: ResetNWPFlowMemory         ! Resets the flow memory of an NWP flow module
                                      ! instance.
Public  :: NWPConvertCoordIndices     ! Returns indices of coord systems in which
                                      ! positions need to be specified when using
                                      ! NWPConvertToZ.
Public  :: NWPConvertToZ              ! Converts vertical coords between coord systems
                                      ! using an NWP flow module instance.
Public  :: NWPFlowCoordIndices        ! Returns indices of coord systems in which
                                      ! positions need to be specified when using
                                      ! GetNWPFlow.
Public  :: GetNWPFlow                 ! Gets flow information from an NWP flow module
                                      ! instance.
Public  :: NWPCloudCoordIndices       ! Returns indices of coord systems in which
                                      ! positions need to be specified when using
                                      ! GetNWPCloud.
Public  :: GetNWPCloud                ! Gets cloud information from an NWP flow module
                                      ! instance.
Public  :: NWPRainCoordIndices        ! Returns indices of coord systems in which
                                      ! positions need to be specified when using
                                      ! GetNWPRain.
Public  :: GetNWPRain                 ! Gets rain information from an NWP flow module
                                      ! instance.
Public  :: NWPSurfaceCoordIndices     ! Returns indices of coord systems in which
                                      ! positions need to be specified when using
                                      ! GetNWPSurface.
Public  :: GetNWPSurface              ! Gets surface information from an NWP flow module
                                      ! instance.
Public  :: NWPSoilCoordIndices        ! Returns indices of coord systems in which
                                      ! positions need to be specified when using
                                      ! GetNWPSoil.
Public  :: GetNWPSoil                 ! Gets soil information from an NWP flow module
                                      ! instance.
Public  :: NWPPlantCoordIndices       ! Returns indices of coord systems in which
                                      ! positions need to be specified when using
                                      ! GetNWPPlant.
Public  :: GetNWPPlant                ! Gets plant information from an NWP flow module
                                      ! instance.
Public  :: NWPReflectCoordIndices     ! Returns indices of coord systems in which
                                      ! positions need to be specified when using
                                      ! NWPReflect.
Public  :: NWPReflect                 ! Reflects position of particle or puff centroid
                                      ! using an NWP flow module instance.

!-------------------------------------------------------------------------------------------------------------

Type :: NWPFlow_ ! The state of an NWP flow module instance.
  Type(CommonFlow_)            :: C                        ! Part of flow state common to all flow modules.
  Type(NWPMet_),       Pointer :: NWPMet                   !} Pointer to met modules containing required data
  Type(NWPMet_),       Pointer :: NWPMetSoilMoisture       !} corresponding to met attribute indices
  Type(NWPMet_),       Pointer :: NWPMetCanopyHeight       !} MA_NWPFlow, MA_NWPCloud, MA_NWPRain,
  Type(NWPMet_),       Pointer :: NWPMetCanopy             !} MA_NWPSoilMoisture or MA_AncillarySoilMoisture,
  Type(AncillaryMet_), Pointer :: AncillaryMetSoilMoisture !} MA_NWPCanopyHeight or MA_AncillaryCanopyHeight,
  Type(AncillaryMet_), Pointer :: AncillaryMetSoil         !} MA_NWPCanopy, MA_AncillarySoil,
  Type(AncillaryMet_), Pointer :: AncillaryMetLandUse      !} MA_AncillaryLandUse and MA_AncillaryLAI.
  Type(AncillaryMet_), Pointer :: AncillaryMetLAI          !}
  Type(AncillaryMet_), Pointer :: AncillaryMetCanopyHeight !}
  Integer                      :: iFlow                    !] Index in list of met modules in C of met module
  Integer                      :: iCloud                   !] for providing information corresponding to met
  Integer                      :: iRain                    !] attribute indices MA_NWPFlow, MA_NWPCloud,
  Integer                      :: iSoilMoisture            !] MA_NWPRain, MA_NWPSoilMoisture or
  Integer                      :: iSoil                    !] MA_AncillarySoilMoisture, MA_AncillarySoil,
  Integer                      :: iLandUse                 !] MA_AncillaryLandUse, MA_AncillaryLAI,
  Integer                      :: iLAI                     !] MA_NWPCanopyHeight or MA_AncillaryCanopyHeight
  Integer                      :: iCanopyHeight            !] and MA_NWPCanopy.
  Integer                      :: iCanopy                  !]
  Logical                      :: NWPSoilMoisture          ! Indicates soil moisture is to be obtained from an
                                                           ! NWP met module.
  Logical                      :: NWPCanopyHeight          ! Indicates canopy height is to be obtained from an
                                                           ! NWP met module.
  Logical                      :: UrbanCanopy              ! Indicates urban canopy effects are to be 
                                                           ! simulated.
  ! $$ more logical to split NWPMet into Flow, Cloud, Rain.
  ! $$ perhaps iFlow etc not needed here - could be local.
End Type NWPFlow_

!-------------------------------------------------------------------------------------------------------------

Type :: NWPFlowMemory_ ! Information related to a particular space-time point of
                       ! interest which an NWP flow module instance wishes to record.
                       ! The purpose of storing this information together is to
                       ! provide a memory to avoid the need to recalculate the
                       ! information.
  Private
  Type(HCoeffs_) :: HCoeffs    !} Interpolation coefficients for the main horizontal and
  Type(ZCoeffs_) :: ZCoeffs    !} vertical grids.
  Type(TCoeffs_) :: TCoeffs    ! Interpolation coefficients for interpolation in time.
  Type(TCoeffs_) :: TCoeffsEnd ! Interpolation coefficients for interpolation in time to a time at the end of
                               ! the NWP met time interval.
End Type NWPFlowMemory_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitNWPFlow(                              &
           FlowName,                               &
           MetModName, MetName, AncillaryMetNames, &
           DomainName,                             &
           FixedMet,                               &
           UpdateOnDemand,                         &
           UrbanCanopy                             &
         )                                         &
Result(NWPFlow)
! Initialises the state of an NWP flow module instance.

  Implicit None
  ! Argument list:
  Character(*),             Intent(In) :: FlowName             ! Name of flow module instance.
  Character(*),             Intent(In) :: MetModName           ! Name of used met module.
  Character(*),             Intent(In) :: MetName              ! Name of used met module instance.
  Character(MaxCharLength), Intent(In) :: AncillaryMetNames(:) ! Names of used ancillary met module instances.
                                                               ! $$ changed from * to support Intel compiler on Cray
  Character(*),             Intent(In) :: DomainName           ! Name of the domain of the flow module instance.
  Logical,                  Intent(In) :: FixedMet             ! Indicates that the met is fixed.
  Logical,                  Intent(In) :: UpdateOnDemand       ! Indicates the flow module instance is to be
                                                               ! updated using update-on-demand.
  Logical,                  Intent(In) :: UrbanCanopy          ! Indicates urban canopy effects are to be
                                                               ! simulated.
  ! Function result:
  Type(NWPFlow_) :: NWPFlow ! Initialised state of an NWP flow module instance.
  ! Locals:
  Integer                                               :: i             ! Loop index.
  Character(Max(Len(MetModName), Len('Ancillary Met'))) :: Temp(MaxMets) ! Temporary character array.

  Temp(1)                             = MetModName
  Temp(2:1 + Size(AncillaryMetNames)) = 'Ancillary Met'

  ! Met module instance used.
  If (.not.(MetModName .CIEq. 'NWP Met')) Then
    Call Message(                                                  &
           'Error in InitNWPFlow: met module name is given as ' // &
           Trim(MetModName)                                     // &
           ' instead of NWP Met',                                  &
           3                                                       &
         )
  End If

  If (.not. IsCalendar()) Then
    Call Message('FATAL ERROR in InitNWPFlow: The NWP flow module only works in a calendar time frame', 3)
  End If

  NWPFlow%C = InitCommonFlow(                                         &
                FlowModName    = 'NWP Flow',                          &
                FlowName       = FlowName,                            &
                nMets          = 1 + Size(AncillaryMetNames),         &
                MetModNames    = Temp(1:1 + Size(AncillaryMetNames)), &
                MetNames       = (/ MetName,    AncillaryMetNames /), &
                DomainName     = DomainName,                          &
                nAttribs       = 8,                                   &
                AttribNames    = (/                                   &
                                   'Update ',                         &
                                   'Convert',                         &
                                   'Flow   ',                         &
                                   'Cloud  ',                         &
                                   'Rain   ',                         &
                                   'Surface',                         &
                                   'Soil   ',                         &
                                   'Plant  '                          &
                                 /),                                  &
                FixedMet       = FixedMet,                            &
                UpdateOnDemand = UpdateOnDemand                       &
              )

  NWPFlow%UrbanCanopy = UrbanCanopy

End Function InitNWPFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpCoordsEtc_NWPFlow(NWPFlow, Coords, Grids)
! Sets up Coords and Grids by adding any extra coords and grids which NWPFlow wants
! to define.

  Implicit None
  ! Argument list:
  Type(NWPFlow_), Intent(In)    :: NWPFlow ! State of an NWP flow module instance.
  Type(Coords_),  Intent(InOut) :: Coords  ! Collection of coord systems.
  Type(Grids_),   Intent(InOut) :: Grids   ! Collection of grids.
  ! Locals:
  Type(HCoord_) :: HCoord  !} Coord systems to be added.
  Type(ZCoord_) :: ZCoordZ !}
  Type(ZCoord_) :: ZCoordP !}

  ! Construct coord systems to be added to Coords.
  HCoord  = HCoord_LatLong() ! $$ should be adding ps systems near poles.
  ZCoordZ = ZCoord_m_agl()
  ZCoordP = ZCoord_Pa()

  ! Add coord systems to Coords.
  Call AddHCoord(HCoord,  Coords)
  Call AddZCoord(ZCoordZ, Coords)
  Call AddZCoord(ZCoordP, Coords)

End Subroutine SetUpCoordsEtc_NWPFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpNWPFlow_MetsCoordsEtc(EtaDefns, Coords, Grids, Mets, NWPFlow)
! Sets up NWPFlow using information from EtaDefns, Coords, Grids and Mets.

  Implicit None
  ! Argument list:
  Type(EtaDefns_), Intent(In)    :: EtaDefns ! Collection of eta definitions.
  Type(Coords_),   Intent(In)    :: Coords   ! Collection of coord systems.
  Type(Grids_),    Intent(In)    :: Grids    ! Collection of grids.
  Type(Mets_),     Intent(In)    :: Mets     ! Set of met module instance states.
  Type(NWPFlow_),  Intent(InOut) :: NWPFlow  ! State of an NWP flow module instance.
  ! Locals:
  Integer                  :: iMetMod     ! Met module index.
  Integer                  :: iMet        ! Met module instance index.
  Character(MaxCharLength) :: HCoordName1 !} Names of coord systems to add to
  Character(MaxCharLength) :: HCoordName2 !} NWPFlow%C.
  Character(MaxCharLength) :: ZCoordName1 !}
  Character(MaxCharLength) :: ZCoordName2 !}
  Character(MaxCharLength) :: ZCoordName3 !}

  ! Find index of met module instance.
  Call FindMetIndex(                                      &
         NWPFlow%C%MetModNames(1), NWPFlow%C%MetNames(1), &
         Mets,                                            &
         iMetMod, iMet                                    &
       )

  ! Set up the names of the coord systems to add to NWPFlow%C.
  ! $$ should do for all met modules used or check coords (and grids?) the same
  HCoordName1 = Mets%C(iMetMod, iMet)%P%HCoordNames(1)
  HCoordName2 = 'Lat-Long'
  ZCoordName1 = Mets%C(iMetMod, iMet)%P%ZCoordNames(1)
  ZCoordName2 = 'm agl'
  ZCoordName3 = 'Pa'

  ! Add coord systems from met module instance to NWPFlow%C.
  Call AddCoordsToCommonFlow(                            &
         2, (/ HCoordName1, HCoordName2 /),              &
         3, (/ ZCoordName1, ZCoordName2, ZCoordName3 /), &
         NWPFlow%C                                       &
       )

  ! Check number of coord systems in NWPFlow.
  If (NWPFlow%C%nHCoords /= 2) Then
    Call Message(                                                &
           'Unexpected error in SetUpNWPFlow_MetsCoordsEtc: ' // &
           'error in numbering of horizontal coord systems '  // &
           'used by the NWP flow module instance "'           // &
           Trim(NWPFlow%C%FlowName)                           // &
           '"',                                                  &
           3                                                     &
         )
  End If
  If (NWPFlow%C%nZCoords /= 3) Then
    Call Message(                                                &
           'Unexpected error in SetUpNWPFlow_MetsCoordsEtc: ' // &
           'error in numbering of vertical coord systems '    // &
           'used by the NWP flow module instance "'           // &
           Trim(NWPFlow%C%FlowName)                           // &
           '"',                                                  &
           3                                                     &
         )
  End If

End Subroutine SetUpNWPFlow_MetsCoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForUpdateNWPFlow(  &
             Coords, Grids, Domains, &
             Mets,                   &
             iCase,                  &
             Time,                   &
             TValid, UpdateNow,      &
             NWPFlow,                &
             Units                   &
           )
! Prepares for updating an instance of the NWP flow module.

! This routine must set TValid and UpdateNow but must not alter the validity of the NWP flow module instance.
! $$ this comment should be with generic definition.

  Implicit None
  ! Argument list:
  Type(Coords_),  Intent(In)           :: Coords
  Type(Grids_),   Intent(In)           :: Grids
  Type(Domains_), Intent(In),   Target :: Domains
  Type(Mets_),    Intent(In)           :: Mets
  Integer,        Intent(In)           :: iCase
  Type(Time_),    Intent(In)           :: Time
  Type(Time_),    Intent(Out)          :: TValid
  Logical,        Intent(Out)          :: UpdateNow
  Type(NWPFlow_), Intent(InOut)        :: NWPFlow
  Type(Units_),   Intent(InOut)        :: Units
  ! Coords    :: Collection of coord systems.
  ! Grids     :: Collection of grids.
  ! Domains   :: Collection of domains.
  ! Mets      :: Set of met module instance states.
  ! iCase     :: Number of case.
  ! Time      :: Time for which the NWP flow module instance might be updated.
  ! TValid    :: Earliest time that the validity (overall or for any single attribute) of the flow module
  !              instance might change, except perhaps for a change caused by a change in the validity of the
  !              met and flow module instances acting as data sources, assuming the flow module instance is
  !              updated now. The value is that determined at the end of this routine (the actual time may be
  !              later).
  ! UpdateNow :: Indicates the flow module instance must be updated now (even if update-on-demand is
  !              specified). If set, TValid need not be set to any particular time.
  ! NWPFlow   :: State of an NWP flow module instance.
  ! Units     :: Collection of information on input/output unit numbers.
  ! Locals:
  Type(Domain_), Pointer :: Domain ! Abbreviation for domain.

  ! Abbreviations.
  Domain => Domains%Domains(NWPFlow%C%iDomain)

  ! Validity due to time domain.
  If (Time < StartTimeOfDomain(Domain)) Then
    TValid = StartTimeOfDomain(Domain)
  Else If (Time >= EndTimeOfDomain(Domain)) Then
    TValid = InfFutureTime()
  Else
    TValid = EndTimeOfDomain(Domain)
  End If
  UpdateNow = .false.

End Subroutine PrepareForUpdateNWPFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine NWPFlowReqs(                                                &
             NWPFlow, Time,                                            &
             WantFlowField,     WantCloudField,     WantRainField,     &
             iFlowUpdateSubset, iCloudUpdateSubset, iRainUpdateSubset, &
             FlowField,         CloudField,         RainField          &
           )
! Specifies what information the flow module instance wants from the other flow module
! instances.

  Implicit None
  ! Argument list:
  Type(NWPFlow_),    Intent(In)  :: NWPFlow
  Type(Time_),       Intent(In)  :: Time
  Logical,           Intent(Out) :: WantFlowField
  Logical,           Intent(Out) :: WantCloudField
  Logical,           Intent(Out) :: WantRainField
  Integer,           Intent(Out) :: iFlowUpdateSubset
  Integer,           Intent(Out) :: iCloudUpdateSubset
  Integer,           Intent(Out) :: iRainUpdateSubset
  Type(FlowField_),  Intent(Out) :: FlowField
  Type(CloudField_), Intent(Out) :: CloudField
  Type(RainField_),  Intent(Out) :: RainField
  ! NWPFlow            :: State of an NWP flow module instance.
  ! Time               :: Current time in non-fixed met cases and the time of the met
  !                       in fixed met cases.
  ! WantFlowField      :} Indicates the flow module instance wants information for
  ! WantCloudField     :} various attributes from other flow module instances.
  ! WantRainField      :}
  ! iFlowUpdateSubset  :] Index of update subsets for various attributes in the array
  ! iCloudUpdateSubset :] of update subset information in the common part of the flow
  ! iRainUpdateSubset  :] module instance.
  ! FlowField          :} Field of information for various attributes used here to
  ! CloudField         :} specify what information the flow module instance wants from
  ! RainField          :} the other flow module instances.

  ! Flow requirements.
  WantFlowField = .false.
  ! Avoid compiler warning.
  iFlowUpdateSubset = 0
  FlowField%C%Valid = .false.

  ! Cloud requirements.
  WantCloudField = .false.
  ! Avoid compiler warning.
  iCloudUpdateSubset = 0
  CloudField%C%Valid = .false.

  ! Rain requirements.
  WantRainField = .false.
  ! Avoid compiler warning.
  iRainUpdateSubset = 0
  RainField%C%Valid = .false.

End Subroutine NWPFlowReqs

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateNWPFlow(                      &
             Coords, Grids, Domains,           &
             Mets,                             &
             iCase,                            &
             Time,                             &
             FlowField, CloudField, RainField, &
             NWPFlow,                          &
             Units                             &
           )
! Updates an instance of the NWP flow module.

  Implicit None
  ! Argument list:
  Type(Coords_),     Intent(In)           :: Coords     ! Collection of coord systems.
  Type(Grids_),      Intent(In)           :: Grids      ! Collection of grids.
  Type(Domains_),    Intent(In),   Target :: Domains    ! Collection of domains.
  Type(Mets_),       Intent(In),   Target :: Mets       ! Set of met module instance states.
  Integer,           Intent(In)           :: iCase      ! Number of case.
  Type(Time_),       Intent(In)           :: Time       ! Time for which the NWP flow module instance is to be
                                                        ! updated.
  Type(FlowField_),  Intent(In)           :: FlowField  !} Field of information for various
  Type(CloudField_), Intent(In)           :: CloudField !} attributes needed by the flow
  Type(RainField_),  Intent(In)           :: RainField  !} module instance.
  Type(NWPFlow_),    Intent(InOut)        :: NWPFlow    ! State of an NWP flow module instance.
  Type(Units_),      Intent(InOut)        :: Units      ! Collection of information on input/output unit
                                                        ! numbers.
  ! Locals:
  Type(Domain_), Pointer :: Domain                                ! Abbreviation for domain.
  Logical                :: Valid                                 !} Validity due to time domain.
  Type(Time_)            :: TValid                                !}
  Logical                :: MValid(MaxMets)                       !] Validity of met modules.
  Type(Time_)            :: MTValid(MaxMets)                      !]
  Logical                :: MValidAttribs(MaxMetAttribs, MaxMets) ! Validity of attribs of met modules.
  Integer                :: i                                     ! Loop index.

  ! Abbreviations.
  Domain => Domains%Domains(NWPFlow%C%iDomain)

  ! Validity due to time domain.
  If (Time < StartTimeOfDomain(Domain)) Then
    Valid = .false.
    TValid = StartTimeOfDomain(Domain)
  Else If (Time >= EndTimeOfDomain(Domain)) Then
    Valid = .false.
    TValid = InfFutureTime()
  Else
    Valid = .true.
    TValid = EndTimeOfDomain(Domain)
  End If

  ! Validity of met modules.
  Do i = 1, NWPFlow%C%nMets
    Call MetValid(                                    &
           Mets,                                      &
           NWPFlow%C%iMetMods(i), NWPFlow%C%iMets(i), &
           Time,                                      &
           MValid(i), MTValid(i), MValidAttribs(:, i) &
         )
  End Do

  ! Decide on met module instance for each met attribute.
  NWPFlow%iFlow           = 0
  NWPFlow%iCloud          = 0
  NWPFlow%iRain           = 0
  NWPFlow%iSoilMoisture   = 0
  NWPFlow%iSoil           = 0
  NWPFlow%iLandUse        = 0
  NWPFlow%iLAI            = 0
  NWPFlow%iCanopyHeight   = 0
  NWPFlow%iCanopy         = 0
  NWPFlow%NWPSoilMoisture = .false.
  NWPFlow%NWPCanopyHeight = .false.
  Do i = 1, NWPFlow%C%nMets
    If (Mets%C(NWPFlow%C%iMetMods(i), NWPFlow%C%iMets(i))%P%MetModName == 'NWP Met') Then
      If (MValidAttribs(MA_NWPFlow,         i) .and. NWPFlow%iFlow         == 0) NWPFlow%iFlow         = i
      If (MValidAttribs(MA_NWPCloud,        i) .and. NWPFlow%iCloud        == 0) NWPFlow%iCloud        = i
      If (MValidAttribs(MA_NWPRain,         i) .and. NWPFlow%iRain         == 0) NWPFlow%iRain         = i
      If (MValidAttribs(MA_NWPSoilMoisture, i) .and. NWPFlow%iSoilMoisture == 0) NWPFlow%iSoilMoisture = i
      If (MValidAttribs(MA_NWPCanopyHeight, i) .and. NWPFlow%iCanopyHeight == 0) NWPFlow%iCanopyHeight = i
      If (MValidAttribs(MA_NWPCanopy,       i) .and. NWPFlow%iCanopy       == 0) NWPFlow%iCanopy       = i
    End If
  End Do
  If (NWPFlow%iSoilMoisture /= 0) NWPFlow%NWPSoilMoisture = .true.
  If (NWPFlow%iCanopyHeight /= 0) NWPFlow%NWPCanopyHeight = .true.
  Do i = 1, NWPFlow%C%nMets
    If (Mets%C(NWPFlow%C%iMetMods(i), NWPFlow%C%iMets(i))%P%MetModName == 'Ancillary Met') Then
      If (MValidAttribs(MA_AncillarySoil,         i) .and. NWPFlow%iSoil         == 0) &
        NWPFlow%iSoil         = i
      If (MValidAttribs(MA_AncillaryLandUse,      i) .and. NWPFlow%iLandUse      == 0) &
        NWPFlow%iLandUse      = i
      If (MValidAttribs(MA_AncillarySoilMoisture, i) .and. NWPFlow%iSoilMoisture == 0) &
        NWPFlow%iSoilMoisture = i
      If (MValidAttribs(MA_AncillaryLAI,          i) .and. NWPFlow%iLAI          == 0) &
        NWPFlow%iLAI          = i
      If (MValidAttribs(MA_AncillaryCanopyHeight, i) .and. NWPFlow%iCanopyHeight == 0) &
        NWPFlow%iCanopyHeight = i
    End If
  End Do

  ! Set pointers to met modules.
  If (NWPFlow%iFlow /= 0) Then
    NWPFlow%NWPMet => Mets%NWPMets(NWPFlow%C%iMets(NWPFlow%iFlow)) ! $$ cloud/rain assumed the same
  End If
  If (NWPFlow%iSoilMoisture /= 0) Then
    If (NWPFlow%NWPSoilMoisture) Then
      NWPFlow%NWPMetSoilMoisture => Mets%NWPMets(NWPFlow%C%iMets(NWPFlow%iSoilMoisture))
    Else
      NWPFlow%AncillaryMetSoilMoisture => Mets%AncillaryMets(NWPFlow%C%iMets(NWPFlow%iSoilMoisture))
    End If
  End If
  If (NWPFlow%iCanopyHeight /= 0) Then
    If (NWPFlow%NWPCanopyHeight) Then
      NWPFlow%NWPMetCanopyHeight => Mets%NWPMets(NWPFlow%C%iMets(NWPFlow%iCanopyHeight))
    Else
      NWPFlow%AncillaryMetCanopyHeight => Mets%AncillaryMets(NWPFlow%C%iMets(NWPFlow%iCanopyHeight))
    End If
  End If
  If (NWPFlow%iCanopy /= 0) Then
    NWPFlow%NWPMetCanopy => Mets%NWPMets(NWPFlow%C%iMets(NWPFlow%iCanopy))
  End If
  If (NWPFlow%iSoil /= 0) Then
    NWPFlow%AncillaryMetSoil => Mets%AncillaryMets(NWPFlow%C%iMets(NWPFlow%iSoil))
  End If
  If (NWPFlow%iLandUse /= 0) Then
    NWPFlow%AncillaryMetLandUse => Mets%AncillaryMets(NWPFlow%C%iMets(NWPFlow%iLandUse))
  End If
  If (NWPFlow%iLAI /= 0) Then
    NWPFlow%AncillaryMetLAI => Mets%AncillaryMets(NWPFlow%C%iMets(NWPFlow%iLAI))
  End If

  ! Overall validity.
  NWPFlow%C%ValidAttribParams(:) = .false.

  NWPFlow%C%ValidAttribParams(A_Convert) = Valid .and. NWPFlow%iFlow  /= 0
  NWPFlow%C%ValidAttribParams(A_Flow   ) = Valid .and. NWPFlow%iFlow  /= 0
  NWPFlow%C%ValidAttribParams(A_Cloud  ) = Valid .and. NWPFlow%iCloud /= 0
  NWPFlow%C%ValidAttribParams(A_Rain   ) = Valid .and. NWPFlow%iRain  /= 0
  NWPFlow%C%ValidAttribParams(A_Surface) = Valid                      .and. &
                                           NWPFlow%iLandUse      /= 0 .and. &
                                           NWPFlow%iSoilMoisture /= 0
  NWPFlow%C%ValidAttribParams(A_Soil   ) = Valid                      .and. &
                                           NWPFlow%iSoil         /= 0
  NWPFlow%C%ValidAttribParams(A_Plant  ) = Valid                      .and. &
                                           NWPFlow%iCanopyHeight /= 0 .and. &
                                           NWPFlow%iLAI          /= 0 .and. &
                                           NWPFlow%iCanopy       /= 0

  NWPFlow%C%Valid = NWPFlow%C%ValidAttribParams(A_Convert) .or. &
                    NWPFlow%C%ValidAttribParams(A_Flow   ) .or. &
                    NWPFlow%C%ValidAttribParams(A_Cloud  ) .or. &
                    NWPFlow%C%ValidAttribParams(A_Rain   ) .or. &
                    NWPFlow%C%ValidAttribParams(A_Surface) .or. &
                    NWPFlow%C%ValidAttribParams(A_Soil)    .or. &
                    NWPFlow%C%ValidAttribParams(A_Plant)

  NWPFlow%C%ValidityUnimprovable = NWPFlow%C%ValidAttribParams(A_Convert) .and. &
                                   NWPFlow%C%ValidAttribParams(A_Flow   ) .and. &
                                   NWPFlow%C%ValidAttribParams(A_Cloud  ) .and. &
                                   NWPFlow%C%ValidAttribParams(A_Rain   ) .and. &
                                   NWPFlow%C%ValidAttribParams(A_Surface) .and. &
                                   NWPFlow%C%ValidAttribParams(A_Soil   ) .and. &
                                   NWPFlow%C%ValidAttribParams(A_Plant  )
  NWPFlow%C%ValidityUnimprovable = NWPFlow%C%ValidityUnimprovable .or. NWPFlow%C%FixedMet

  NWPFlow%C%TValid = TValid

  If (NWPFlow%C%ValidAttribParams(A_Convert)) Then
    NWPFlow%C%TValid = TMin(NWPFlow%C%TValid, MTValid(NWPFlow%iFlow        ))
  End If
  If (NWPFlow%C%ValidAttribParams(A_Flow   )) Then
    NWPFlow%C%TValid = TMin(NWPFlow%C%TValid, MTValid(NWPFlow%iFlow        ))
  End If
  If (NWPFlow%C%ValidAttribParams(A_Cloud  )) Then
    NWPFlow%C%TValid = TMin(NWPFlow%C%TValid, MTValid(NWPFlow%iCloud       ))
  End If
  If (NWPFlow%C%ValidAttribParams(A_Rain   )) Then
    NWPFlow%C%TValid = TMin(NWPFlow%C%TValid, MTValid(NWPFlow%iRain        ))
  End If
  If (NWPFlow%C%ValidAttribParams(A_Surface)) Then
    NWPFlow%C%TValid = TMin(NWPFlow%C%TValid, MTValid(NWPFlow%iLandUse     ))
  End If
  If (NWPFlow%C%ValidAttribParams(A_Surface)) Then
    NWPFlow%C%TValid = TMin(NWPFlow%C%TValid, MTValid(NWPFlow%iSoilMoisture))
  End If
  If (NWPFlow%C%ValidAttribParams(A_Soil   )) Then
    NWPFlow%C%TValid = TMin(NWPFlow%C%TValid, MTValid(NWPFlow%iSoil        ))
  End If
  If (NWPFlow%C%ValidAttribParams(A_Plant  )) Then
    NWPFlow%C%TValid = TMin(NWPFlow%C%TValid, MTValid(NWPFlow%iLAI         ))
  End If
  If (NWPFlow%C%ValidAttribParams(A_Plant  )) Then
    NWPFlow%C%TValid = TMin(NWPFlow%C%TValid, MTValid(NWPFlow%iCanopyHeight))
  End If
  If (NWPFlow%C%ValidAttribParams(A_Plant  )) Then
    NWPFlow%C%TValid = TMin(NWPFlow%C%TValid, MTValid(NWPFlow%iCanopy      ))
  End If

  ! Added because Name II has Clay = Min(Clay, 1.0). Unclear if needed. $$
  If (NWPFlow%C%ValidAttribParams(A_Soil)) Then
    If (MaxVal(NWPFlow%AncillaryMetSoil%ClayMassFrac(:, :, :, :)) > 1.0) Then
      Call Message('WARNING: ClayMassFraction > 1', 2)
    End If
  End If

End Subroutine UpdateNWPFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine ResetNWPFlowMemory(NWPFlowMemory)
! Resets the flow memory of an NWP flow module instance.

  Implicit None
  ! Argument list:
  Type(NWPFlowMemory_), Intent(InOut) :: NWPFlowMemory ! Flow memory of NWP flow
                                                       ! module instance.

  ! Mark the interpolation coeeficients as not valid.
  NWPFlowMemory%HCoeffs%Valid    = .false.
  NWPFlowMemory%ZCoeffs%Valid    = .false.
  NWPFlowMemory%TCoeffs%Valid    = .false.
  NWPFlowMemory%TCoeffsEnd%Valid = .false.

End Subroutine ResetNWPFlowMemory

!-------------------------------------------------------------------------------------------------------------

Subroutine NWPConvertCoordIndices(NWPFlow, nHCoords, iHCoords)
! Returns indices of coord systems in which positions need to be specified when using
! NWPConvertToZ.

  Implicit None
  ! Argument list:
  Type(NWPFlow_), Intent(In)  :: NWPFlow              ! State of an NWP flow module
                                                      ! instance.
  Integer,        Intent(Out) :: nHCoords             !} Number and values of indices
  Integer,        Intent(Out) :: iHCoords(MaxHCoords) !} of horizontal coord systems.

  ! Horizontal coords needed:
  ! 1) Coord system for horizontal grids used for met data.
  nHCoords    = 1
  iHCoords(1) = NWPFlow%C%iHCoords(1)

End Subroutine NWPConvertCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine NWPConvertToZ(             &
             NWPFlow, Coords, Grids,  &
             iZCoord, Time, Position, &
             NWPFlowMemory            &
           )
! Converts vertical coords between coord systems using an NWP flow module instance.

  Implicit None
  ! Argument list:
  Type(NWPFlow_),       Intent(In)            :: NWPFlow
  Type(Coords_),        Intent(In),    Target :: Coords
  Type(Grids_),         Intent(In)            :: Grids
  Integer,              Intent(In)            :: iZCoord
  Type(ShortTime_),     Intent(In)            :: Time
  Type(Position_),      Intent(InOut)         :: Position
  Type(NWPFlowMemory_), Intent(InOut), Target :: NWPFlowMemory
  ! NWPFlow       :: State of an NWP flow module instance.
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! iZCoord       :: Index in Coords of coord system to convert to.
  ! Time          :: Current time in non-fixed met cases and the time of the met in
  !                  fixed met cases.
  ! Position      :: Coords of the point in various coord systems in Coords, with
  !                  flags to indicate whether the values are valid.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Locals:
  Type(ZCoord_), Pointer :: ZCoordIn  !} Abbreviations for the coord system to convert
  Type(ZCoord_), Pointer :: ZCoordOut !} from, the coord system to convert to, and two
  Type(ZCoord_), Pointer :: ZCoord1   !} intermediate coord systems.
  Type(ZCoord_), Pointer :: ZCoord2   !}
  Real(Std)              :: ZIn       !] Height in coord systems ZCoordIn, ZCoordOut,
  Real(Std)              :: ZOut      !] ZCoord1 and ZCoord2.
  Real(Std)              :: Z1        !]
  Real(Std)              :: Z2        !]
  Integer                :: iZUse     ! Index of the coord system used to convert
                                      ! from.
  Real(Std)              :: Topog     ! Topographic height above sea level (m).
  Real(Std)              :: PS        ! Surface pressure (Pa).
  Integer                :: i         ! Loop index.

  Do i = 1, Coords%nZCoords ! replace by efficient selection and store intermediate results $$
    If (Position%ZValid(i)) Then
      iZUse = i
      Exit
    End If
  End Do

  ZCoordIn  => Coords%ZCoords(iZUse)
  ZCoordOut => Coords%ZCoords(iZCoord)
  ZIn       =  Position%Z(iZUse)

  If (.not. (Position%TopogValid .and. Position%PSValid)) Then
    Call TopogAndPS(                              &
           Grids,                                 &
           NWPFlow%NWPMet,                        &
           Time,                                  &
           Position%XY(1, NWPFlow%C%iHCoords(1)), &
           Position%XY(2, NWPFlow%C%iHCoords(1)), &
           NWPFlowMemory, Topog, PS               &
         )
    Position%Topog      = Topog
    Position%PS         = PS
    Position%TopogValid = .true.
    Position%PSValid    = .true.
  End If
  Topog = Position%Topog
  PS    = Position%PS

  Select Case (ZCoordIn%CoordType)

    Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)

      Select Case (ZCoordOut%CoordType)

        Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)

          ZOut = ZBasedToZBased(ZCoordIn, ZCoordOut, Topog, ZIn)

        Case (Z_P, Z_PAsZ, Z_PAsEta)

          Select Case (Coords%ZCoords(NWPFlow%C%iZCoords(1))%CoordType) ! Met coord.
            Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)
              ZCoord1 => Coords%ZCoords(NWPFlow%C%iZCoords(1)) ! Met coord.
              ZCoord2 => Coords%ZCoords(NWPFlow%C%iZCoords(3)) ! Pressure.
              ! Note height-based eta to height-based eta conversions go via height
              ! when calculated with ZBasedToZBased. Hence this if-test.
              If (iZCoord /= NWPFlow%C%iZCoords(1)) Then
                Z1 = ZBasedToZBased(ZCoordIn, ZCoord1, Topog, ZIn)
              Else
                Z1 = ZIn
              End If
              Call EtaToP(                                  &
                     Coords, Grids,                         &
                     NWPFlow%NWPMet,                        &
                     Time,                                  &
                     Position%XY(1, NWPFlow%C%iHCoords(1)), &
                     Position%XY(2, NWPFlow%C%iHCoords(1)), &
                     Z1, Z2,                                &
                     NWPFlowMemory                          &
                   )
              ZOut = PBasedToPBased(ZCoord2, ZCoordOut, PS, Z2)
            Case (Z_P, Z_PAsZ, Z_PAsEta)
              ZCoord1 => Coords%ZCoords(NWPFlow%C%iZCoords(2)) ! Height above ground.
              ZCoord2 => Coords%ZCoords(NWPFlow%C%iZCoords(1)) ! Met coord.
              Z1 = ZBasedToZBased(ZCoordIn, ZCoord1, Topog, ZIn)
              Call ZToEta(                                  &
                     Coords, Grids,                         &
                     NWPFlow%NWPMet,                        &
                     Time,                                  &
                     Position%XY(1, NWPFlow%C%iHCoords(1)), &
                     Position%XY(2, NWPFlow%C%iHCoords(1)), &
                     Z1, Z2,                                &
                     NWPFlowMemory                          &
                   )
              ! Note pressure-based eta to pressure-based eta conversions go via
              ! pressure when calculated with PBasedToPBased. Hence this if-test.
              If (iZCoord /= NWPFlow%C%iZCoords(1)) Then
                ZOut = PBasedToPBased(ZCoord2, ZCoordOut, PS, Z2)
              Else
                ZOut = Z2
              End If
            Case Default
              Call Message('UNEXPECTED FATAL ERROR in NWPConvertToZ: unknown coordinate type', 4)
          End Select

        Case Default

          Call Message('UNEXPECTED FATAL ERROR in NWPConvertToZ: unknown coordinate type', 4)

      End Select

    Case (Z_P, Z_PAsZ, Z_PAsEta)

      Select Case (ZCoordOut%CoordType)

        Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)

          Select Case (Coords%ZCoords(NWPFlow%C%iZCoords(1))%CoordType) ! Met coord.
            Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)
              ZCoord1 => Coords%ZCoords(NWPFlow%C%iZCoords(3)) ! Pressure.
              ZCoord2 => Coords%ZCoords(NWPFlow%C%iZCoords(1)) ! Met coord.
              Z1 = PBasedToPBased(ZCoordIn, ZCoord1, PS, ZIn)
              Call PToEta(                                  &
                     Coords, Grids,                         &
                     NWPFlow%NWPMet,                        &
                     Time,                                  &
                     Position%XY(1, NWPFlow%C%iHCoords(1)), &
                     Position%XY(2, NWPFlow%C%iHCoords(1)), &
                     Z1, Z2,                                &
                     NWPFlowMemory                          &
                   )
              ! Note height-based eta to height-based eta conversions go via height
              ! when calculated with ZBasedToZBased. Hence this if-test.
              If (iZUse /= NWPFlow%C%iZCoords(1)) Then
                ZOut = ZBasedToZBased(ZCoord2, ZCoordOut, Topog, Z2)
              Else
                ZOut = Z2
              End If
            Case (Z_P, Z_PAsZ, Z_PAsEta)
              ZCoord1 => Coords%ZCoords(NWPFlow%C%iZCoords(1)) ! Met coord.
              ZCoord2 => Coords%ZCoords(NWPFlow%C%iZCoords(2)) ! Height above ground.
              ! Note pressure-based eta to pressure-based eta conversions go via pressure
              ! when calculated with PBasedToPBased. Hence this if-test.
              If (iZUse /= NWPFlow%C%iZCoords(1)) Then
                Z1 = PBasedToPBased(ZCoordIn, ZCoord1, PS, ZIn)
              Else
                Z1 = ZIn
              End If
              Call EtaToZ(                                  &
                     Coords, Grids,                         &
                     NWPFlow%NWPMet,                        &
                     Time,                                  &
                     Position%XY(1, NWPFlow%C%iHCoords(1)), &
                     Position%XY(2, NWPFlow%C%iHCoords(1)), &
                     Z1, Z2,                                &
                     NWPFlowMemory                          &
                   )
              ZOut = ZBasedToZBased(ZCoord2, ZCoordOut, Topog, Z2)
            Case Default
              Call Message('UNEXPECTED FATAL ERROR in NWPConvertToZ: unknown coordinate type', 4)
          End Select

        Case (Z_P, Z_PAsZ, Z_PAsEta)

          ZOut = PBasedToPBased(ZCoordIn, ZCoordOut, PS, ZIn)

        Case Default

          Call Message('UNEXPECTED FATAL ERROR in NWPConvertToZ: unknown coordinate type', 4)

      End Select

    Case Default

      Call Message('UNEXPECTED FATAL ERROR in NWPConvertToZ: unknown coordinate type', 4)

  End Select

  Position%Z(iZCoord) = ZOut

End Subroutine NWPConvertToZ

!-------------------------------------------------------------------------------------------------------------

Subroutine NWPFlowCoordIndices(NWPFlow, nHCoords, iHCoords, nZCoords, iZCoords)
! Returns indices of coord systems in which positions need to be specified when using GetNWPFlow.

  Implicit None
  ! Argument list:
  Type(NWPFlow_), Intent(In)  :: NWPFlow              ! State of an NWP flow module instance.
  Integer,        Intent(Out) :: nHCoords             !} Number and values of indices of horizontal and
  Integer,        Intent(Out) :: iHCoords(MaxHCoords) !} vertical coord systems.
  Integer,        Intent(Out) :: nZCoords             !}
  Integer,        Intent(Out) :: iZCoords(MaxHCoords) !}

  ! Horizontal coords needed:
  ! 1) Coord system of horizontal grids used for met data.
  nHCoords    = 1
  iHCoords(1) = NWPFlow%C%iHCoords(1)

  ! Vertical coords needed:
  ! 1) Coord system of vertical grids used for met data.
  ! 2) m agl.
  nZCoords    = 2
  iZCoords(1) = NWPFlow%C%iZCoords(1)
  iZCoords(2) = NWPFlow%C%iZCoords(2)

End Subroutine NWPFlowCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine GetNWPFlow(                 &
             Coords, Grids,            &
             NWPFlow,                  &
             Moisture, Inhomog, Homog, &
             Time, Position,           &
             NWPFlowMemory,            &
             Flow, ProfileData         &
           )
! Gets flow information from an NWP flow module instance.

  Implicit None
  ! Argument list:
  Type(Coords_),        Intent(In)             :: Coords
  Type(Grids_),         Intent(In)             :: Grids
  Type(NWPFlow_),       Intent(In)             :: NWPFlow
  Logical,              Intent(In)             :: Moisture
  Logical,              Intent(In)             :: Inhomog
  Logical,              Intent(In)             :: Homog
  Type(ShortTime_),     Intent(In)             :: Time
  Type(Position_),      Intent(In)             :: Position
  Type(NWPFlowMemory_), Intent(InOut)          :: NWPFlowMemory
  Type(Flow_),          Intent(Out)            :: Flow
  Type(ProfileData_),   Intent(Out),  Optional :: ProfileData
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! NWPFlow       :: State of an NWP flow module instance.
  ! Moisture      :: Indicates Q is required.
  ! Inhomog       :: Indicates inhomogeneous quantities are required.
  ! Homog         :: Indicates homogeneous quantities are required.
  ! Time          :: Current time in non-fixed met cases and the time of the met in
  !                  fixed met cases.
  ! Position      :: Coords of the point in various coord systems in Coords, with
  !                  flags to indicate whether the values are valid.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Flow          :: Flow information.
  ! ProfileData   :: Data needed to construct idealised analytic mean flow and
  !                  turbulence profiles.
  ! Locals:
  Real(Std)          :: X            !} Horizontal location in coord system of horizontal grids used for met
  Real(Std)          :: Y            !} data.
  Real(Std)          :: Eta          ! Height in coord system of vertical grids used for met data.
  Real(Std)          :: Z            ! Height (m agl).
  Type(ProfileData_) :: ProfileDataL ! Local copy of ProfileData.
  Logical            :: FastRun      ! Indicates fast, but less precise, calculation
                                     ! required.
  Real(Std)          :: TauMin       ! Minimum timescale.
  Real(Std)          :: ZMin         ! Minimum height for turbulence quantities.

  FastRun = .false. ! $$ This should be input (or option removed). But note Delta must vary smoothly.
  TauMin  = 20.0    ! $$ This should be input.
  ZMin    =  0.0    ! $$ This should be input (and should be non-zero).

  ! $$ use ps systems near poles.

  ! Coords of location.
  X   = Position%XY(1, NWPFlow%C%iHCoords(1))
  Y   = Position%XY(2, NWPFlow%C%iHCoords(1))
  Eta = Position%Z (   NWPFlow%C%iZCoords(1))
  Z   = Position%Z (   NWPFlow%C%iZCoords(2))

  ! Calculate the following parts of Flow and ProfileDataL:
  ! Flow: U, dUdX, T, Theta, dThetadX, Q, P, Rho, dRhodX, Topog and dTopogdX.
  ! ProfileDataL: Z0, Phi0, UStar, RecipLMO, H and WStar.
  Call FlowInfoFromMet(           &
         Coords, Grids,           &
         NWPFlow, NWPFlow%NWPMet, &
         Moisture, FastRun,       &
         Time, X, Y, Eta, Z,      &
         Flow, ProfileDataL,      &
         NWPFlowMemory            &
       )

  ! Set up other quantities required in Flow.
  Flow%dUdT(:)     = 0.0 ! $$ Also remember must be in chronological time direction - use IsBackwards
  Flow%ZS          = 0.0
  Flow%SmoothTopog = .true.

  ! Set up other quantities required in ProfileDataL.
  ProfileDataL%SigUUM     = NWPFlow%NWPMet%SigUUM
  ProfileDataL%TauUUM     = NWPFlow%NWPMet%TauUUM
  ProfileDataL%SigU2HPlus = NWPFlow%NWPMet%SigU2HPlus
  ProfileDataL%SigW2HPlus = NWPFlow%NWPMet%SigW2HPlus
  ProfileDataL%TauUHPlus  = NWPFlow%NWPMet%TauUHPlus
  ProfileDataL%TauWHPlus  = NWPFlow%NWPMet%TauWHPlus
  ProfileDataL%iHCoord    = NWPFlow%C%iHCoords(1)
  ProfileDataL%iZCoord    = NWPFlow%C%iZCoords(2)
  ProfileDataL%MeanFlow   = .false.
  ProfileDataL%Turb       = .true.

  ! 'Inhomogeneous turbulence quantities', 'inhomogeneous eddy diffusivities',
  ! 'homogeneous turbulence quantities', 'unresolved mesoscale motion quantities',
  ! 'boundary layer characteristics' and 'coord systems' parts of Flow.
  Call TurbProfiles(                 &
         Z           = Z,            &
         ProfileData = ProfileDataL, &
         Inhomog     = Inhomog,      &
         Homog       = Homog,        &
         UEqV        = .true.,       & ! Use sig u \= sig v $$
         TauMin      = TauMin,       &
         ZMin        = ZMin,         &
         Flow        = Flow          &
       )

  ! Time step limits.
  Flow%MaxDZ = Huge(Flow%MaxDZ)

  If (Inhomog) Then
    Flow%MaxDZ = Min(Flow%MaxDZ, Flow%H)
    If (Flow%dSigUUdX(3) /= 0.0) Then
      Flow%MaxDZ = Min(                                    &
                     Flow%MaxDZ,                           &
                     Abs(Flow%SigUU(3) / Flow%dSigUUdX(3)) &
                   )
    End If
    If (Flow%dTauUUdZ(3) /= 0.0) Then
      Flow%MaxDZ = Min(                               &
                     Flow%MaxDZ,                      &
                     Abs(Flow%TauUU(3) / Flow%dTauUUdZ(3)) &
                   )
    End If
    If (Flow%dKdX(3) /= 0.0) Then
      Flow%MaxDZ = Min(                            &
                     Flow%MaxDZ,                   &
                     Flow%K(3) / Abs(Flow%dKdX(3)) &
                   )
    End If
  End If

  ! Time step limits.
  Flow%MaxDT = Huge(Flow%MaxDT)
  
  ! 60s is the time step to be used when using UKV
  If (NWPFlow%C%FlowName(1:3).EQ.'UKV') Flow%MaxDT = 60.

  ! Puff size limits.
  Flow%DeltaI = Flow%H/4.0 ! $$

  ! ProfileData.
  If (Present(ProfileData)) ProfileData = ProfileDataL

End Subroutine GetNWPFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine NWPCloudCoordIndices(NWPFlow, nHCoords, iHCoords, nZCoords, iZCoords)
! Returns indices of coord systems in which positions need to be specified when using GetNWPCloud.

  Implicit None
  ! Argument list:
  Type(NWPFlow_), Intent(In)  :: NWPFlow              ! State of an NWP flow module instance.
  Integer,        Intent(Out) :: nHCoords             !} Number and values of indices of horizontal and
  Integer,        Intent(Out) :: iHCoords(MaxHCoords) !} vertical coord systems.
  Integer,        Intent(Out) :: nZCoords             !}
  Integer,        Intent(Out) :: iZCoords(MaxHCoords) !}

  ! Horizontal coords needed:
  ! 1) Coord system of horizontal grids used for met data.
  nHCoords    = 1
  iHCoords(1) = NWPFlow%C%iHCoords(1)

  ! Vertical coords needed:
  ! 1) Coord system of vertical grids used for met data.
  nZCoords    = 1
  iZCoords(1) = NWPFlow%C%iZCoords(1)

End Subroutine NWPCloudCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine GetNWPCloud(      &
             Coords, Grids,  &
             NWPFlow,        &
             Time, Position, &
             NWPFlowMemory,  &
             Cloud           &
           )
! Gets cloud information from an NWP flow module instance.

  Implicit None
  ! Argument list:
  Type(Coords_),        Intent(In)    :: Coords
  Type(Grids_),         Intent(In)    :: Grids
  Type (NWPFlow_),      Intent(In)    :: NWPFlow
  Type (ShortTime_),    Intent(In)    :: Time
  Type (Position_),     Intent(In)    :: Position
  Type(NWPFlowMemory_), Intent(InOut) :: NWPFlowMemory
  Type (Cloud_),        Intent(Out)   :: Cloud
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! NWPFlow       :: State of an NWP flow module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in fixed met cases.
  ! Position      :: Coords of the point in various coord systems in Coords, with flags to indicate whether
  !                  the values are valid.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Cloud         :: Cloud information.
  ! Locals:
  Real(Std) :: X   !} Horizontal location in coord system of horizontal grids used for met data.
  Real(Std) :: Y   !}
  Real(Std) :: Eta ! Height in coord system of vertical grids used for met data.

  X   = Position%XY(1, NWPFlow%C%iHCoords(1))
  Y   = Position%XY(2, NWPFlow%C%iHCoords(1))
  Eta = Position%Z (   NWPFlow%C%iZCoords(1))

  Call CloudInfoFromMet(  &
         Grids,           &
         NWPFlow%NWPMet,  &
         Time, X, Y, Eta, &
         Cloud,           &
         NWPFlowMemory    &
       )

  Cloud%iZCoord = NWPFlow%C%iZCoords(2) ! m agl.

End Subroutine GetNWPCloud

!-------------------------------------------------------------------------------------------------------------

Subroutine NWPRainCoordIndices(NWPFlow, nHCoords, iHCoords, nZCoords, iZCoords)
! Returns indices of coord systems in which positions need to be specified when using GetNWPRain.

  Implicit None
  ! Argument list:
  Type(NWPFlow_), Intent(In)  :: NWPFlow              ! State of an NWP flow module instance.
  Integer,        Intent(Out) :: nHCoords             !} Number and values of indices of horizontal and
  Integer,        Intent(Out) :: iHCoords(MaxHCoords) !} vertical coord systems.
  Integer,        Intent(Out) :: nZCoords             !}
  Integer,        Intent(Out) :: iZCoords(MaxHCoords) !}

  ! Horizontal coords needed:
  ! 1) Coord system of horizontal grids used for met data.
  nHCoords    = 1
  iHCoords(1) = NWPFlow%C%iHCoords(1)

  ! Vertical coords needed:
  ! 1) Coord system of vertical grids used for met data.
  nZCoords    = 1
  iZCoords(1) = NWPFlow%C%iZCoords(1)

End Subroutine NWPRainCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine GetNWPRain(       &
             Coords, Grids,  &
             NWPFlow,        &
             Time, Position, &
             NWPFlowMemory,  &
             Rain            &
           )
! Gets rain information from an NWP flow module instance.

  Implicit None
  ! Argument list:
  Type(Coords_),        Intent(In)    :: Coords
  Type(Grids_),         Intent(In)    :: Grids
  Type (NWPFlow_),      Intent(In)    :: NWPFlow
  Type (ShortTime_),    Intent(In)    :: Time
  Type (Position_),     Intent(In)    :: Position
  Type(NWPFlowMemory_), Intent(InOut) :: NWPFlowMemory
  Type (Rain_),         Intent(Out)   :: Rain
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! NWPFlow       :: State of an NWP flow module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in fixed met cases.
  ! Position      :: Coords of the point in various coord systems in Coords, with flags to indicate whether
  !                  the values are valid.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Rain          :: Rain information.
  ! Locals:
  Real(Std) :: X   !} Horizontal location in coord system of horizontal grids used for met data.
  Real(Std) :: Y   !}
  Real(Std) :: Eta ! Height in coord system of vertical grids used for met data.

  X   = Position%XY(1, NWPFlow%C%iHCoords(1))
  Y   = Position%XY(2, NWPFlow%C%iHCoords(1))
  Eta = Position%Z (   NWPFlow%C%iZCoords(1))

  Call RainInfoFromMet(   &
         Grids,           &
         NWPFlow%NWPMet,  &
         Time, X, Y, Eta, &
         Rain,            &
         NWPFlowMemory    &
       )

End Subroutine GetNWPRain

!-------------------------------------------------------------------------------------------------------------

Subroutine NWPSurfaceCoordIndices(NWPFlow, nHCoords, iHCoords, nZCoords, iZCoords)
! Returns indices of coord systems in which positions need to be specified when using GetNWPSurface.

  Implicit None
  ! Argument list:
  Type(NWPFlow_), Intent(In)  :: NWPFlow              ! State of an NWP flow module instance.
  Integer,        Intent(Out) :: nHCoords             !} Number and values of indices of horizontal and
  Integer,        Intent(Out) :: iHCoords(MaxHCoords) !} vertical coord systems.
  Integer,        Intent(Out) :: nZCoords             !}
  Integer,        Intent(Out) :: iZCoords(MaxHCoords) !}

  ! Horizontal coords needed:
  ! 1) Coord system of horizontal grids used for met data.
  ! 2) Coord system of horizontal grids used for ancillary data.
  If (NWPFlow%C%iHCoords(1) == NWPFlow%AncillaryMetLandUse%iHCoord) Then
    nHCoords    = 1
    iHCoords(1) = NWPFlow%C%iHCoords(1)
  Else
    nHCoords    = 2
    iHCoords(1) = NWPFlow%C%iHCoords(1)
    iHCoords(2) = NWPFlow%AncillaryMetLandUse%iHCoord
  End If
  ! $$ This needs a bit more thought. E.g. could set up the flow module only to supply ancil surface data -
  ! then NWPFlow%C%iHCoords(1) may not be defined.
  ! $$ First if-block may save a little time, but this may be insignificant. If so, simpler to just use the
  ! else block.

  ! Vertical coords needed: none.
  nZCoords    = 0
  iZCoords(1) = 0 ! avoids compiler warning

End Subroutine NWPSurfaceCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine GetNWPSurface(    &
             Coords, Grids,  &
             NWPFlow,        &
             Time, Position, &
             NWPFlowMemory,  &
             Surface         &
           )
! Gets surface information from an NWP flow module instance.

  Implicit None
  ! Argument list:
  Type(Coords_),        Intent(In)    :: Coords
  Type(Grids_),         Intent(In)    :: Grids
  Type (NWPFlow_),      Intent(In)    :: NWPFlow
  Type (ShortTime_),    Intent(In)    :: Time
  Type (Position_),     Intent(In)    :: Position
  Type(NWPFlowMemory_), Intent(InOut) :: NWPFlowMemory
  Type (Surface_),      Intent(Out)   :: Surface
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! NWPFlow       :: State of an NWP flow module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in fixed met cases.
  ! Position      :: Coords of the point in various coord systems in Coords, with flags to indicate whether
  !                  the values are valid.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Surface       :: Surface information.
  ! Locals:
  Real(Std) :: X      !} Horizontal location in coord system of horizontal grids used for met data.
  Real(Std) :: Y      !}
  Real(Std) :: XAncil !] Horizontal location in coord system of horizontal grids used for ancillary data.
  Real(Std) :: YAncil !]
  Real(Std) :: Eta    ! Height in coord system of vertical grids used for met data.

  X      = Position%XY(1, NWPFlow%C%iHCoords(1)              )
  Y      = Position%XY(2, NWPFlow%C%iHCoords(1)              )
  XAncil = Position%XY(1, NWPFlow%AncillaryMetLandUse%iHCoord)
  YAncil = Position%XY(2, NWPFlow%AncillaryMetLandUse%iHCoord)

  Call SurfaceInfoFromMet(                 &
         Grids,                            &
         NWPFlow%NWPMet,                   &
         NWPFlow%NWPMetSoilMoisture,       &
         NWPFlow%AncillaryMetSoilMoisture, &
         NWPFlow%AncillaryMetLandUse,      &
         NWPFlow%NWPSoilMoisture,          &
         Time, X, Y, XAncil, YAncil,       &
         Surface,                          &
         NWPFlowMemory                     &
       )

End Subroutine GetNWPSurface

!-------------------------------------------------------------------------------------------------------------

Subroutine NWPSoilCoordIndices(NWPFlow, nHCoords, iHCoords, nZCoords, iZCoords)
! Returns indices of coord systems in which positions need to be specified when using GetNWPSoil.

  Implicit None
  ! Argument list:
  Type(NWPFlow_), Intent(In)  :: NWPFlow              ! State of an NWP flow module instance.
  Integer,        Intent(Out) :: nHCoords             !} Number and values of indices of horizontal and
  Integer,        Intent(Out) :: iHCoords(MaxHCoords) !} vertical coord systems.
  Integer,        Intent(Out) :: nZCoords             !}
  Integer,        Intent(Out) :: iZCoords(MaxHCoords) !}

  ! Horizontal coords needed:
  ! 1) Coord system of horizontal grids used for met data.
  ! 2) Coord system of horizontal grids used for ancillary data.
  If (NWPFlow%C%iHCoords(1) == NWPFlow%AncillaryMetSoil%iHCoord) Then
    nHCoords    = 1
    iHCoords(1) = NWPFlow%C%iHCoords(1)
  Else
    nHCoords    = 2
    iHCoords(1) = NWPFlow%C%iHCoords(1)
    iHCoords(2) = NWPFlow%AncillaryMetSoil%iHCoord
  End If
  ! $$ This needs a bit more thought. E.g. could set up the flow module only to supply ancil soil data -
  ! then NWPFlow%C%iHCoords(1) may not be defined.
  ! $$ First if-block may save a little time, but this may be insignificant. If so, simpler to just use the
  ! else block.

  ! Vertical coords needed: none.
  nZCoords    = 0
  iZCoords(1) = 0 ! avoids compiler warning

End Subroutine NWPsoilCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine GetNWPSoil(       &
             Coords, Grids,  &
             NWPFlow,        &
             Time, Position, &
             NWPFlowMemory,  &
             Soil            &
           )
! Gets soil information from an NWP flow module instance.

  Implicit None
  ! Argument list:
  Type(Coords_),        Intent(In)    :: Coords
  Type(Grids_),         Intent(In)    :: Grids
  Type (NWPFlow_),      Intent(In)    :: NWPFlow
  Type (ShortTime_),    Intent(In)    :: Time
  Type (Position_),     Intent(In)    :: Position
  Type(NWPFlowMemory_), Intent(InOut) :: NWPFlowMemory
  Type (Soil_),         Intent(Out)   :: Soil
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! NWPFlow       :: State of an NWP flow module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in fixed met cases.
  ! Position      :: Coords of the point in various coord systems in Coords, with flags to indicate whether
  !                  the values are valid.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Soil          :: Soil information.
  ! Locals:
  Real(Std) :: X      !} Horizontal location in coord system of horizontal grids used for met data.
  Real(Std) :: Y      !}
  Real(Std) :: XAncil !] Horizontal location in coord system of horizontal grids used for ancillary data.
  Real(Std) :: YAncil !]
  Real(Std) :: Eta    ! Height in coord system of vertical grids used for met data.

  X      = Position%XY(1, NWPFlow%C%iHCoords(1))
  Y      = Position%XY(2, NWPFlow%C%iHCoords(1))
  XAncil = Position%XY(1, NWPFlow%AncillaryMetSoil%iHCoord)
  YAncil = Position%XY(2, NWPFlow%AncillaryMetSoil%iHCoord)

  Call SoilInfoFromMet(                    &
         Grids,                            &
         NWPFlow%NWPMet,                   &
         NWPFlow%AncillaryMetSoil,         &
         Time, X, Y, XAncil, YAncil,       &
         Soil,                             &
         NWPFlowMemory                     &
       )

End Subroutine GetNWPSoil

!-------------------------------------------------------------------------------------------------------------

Subroutine NWPPlantCoordIndices(NWPFlow, nHCoords, iHCoords, nZCoords, iZCoords)
! Returns indices of coord systems in which positions need to be specified when using GetNWPPlant.

  Implicit None
  ! Argument list:
  Type(NWPFlow_), Intent(In)  :: NWPFlow              ! State of an NWP flow module instance.
  Integer,        Intent(Out) :: nHCoords             !} Number and values of indices of horizontal and
  Integer,        Intent(Out) :: iHCoords(MaxHCoords) !} vertical coord systems.
  Integer,        Intent(Out) :: nZCoords             !}
  Integer,        Intent(Out) :: iZCoords(MaxHCoords) !}

  ! Horizontal coords needed:
  ! 1) Coord system of horizontal grids used for met data.
  ! 2) Coord system of horizontal grids used for ancillary data.
  If (NWPFlow%C%iHCoords(1) == NWPFlow%AncillaryMetLAI%iHCoord) Then
    nHCoords    = 1
    iHCoords(1) = NWPFlow%C%iHCoords(1)
  Else
    nHCoords    = 2
    iHCoords(1) = NWPFlow%C%iHCoords(1)
    iHCoords(2) = NWPFlow%AncillaryMetLAI%iHCoord
  End If
  ! $$ This needs a bit more thought. E.g. could set up the flow module only to supply ancil plant data -
  ! then NWPFlow%C%iHCoords(1) may not be defined.
  ! $$ First if-block may save a little time, but this may be insignificant. If so, simpler to just use the
  ! else block.

  ! Vertical coords needed: none.
  ! 1) Coord system of vertical grids used for met data.
  nZCoords    = 0
  iZCoords(1) = 0 ! avoids compiler warning

End Subroutine NWPPlantCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine GetNWPPlant(      &
             Coords, Grids,  &
             NWPFlow,        &
             Time, Position, &
             NWPFlowMemory,  &
             Plant           &
           )
! Gets plant information from an NWP flow module instance.

  Implicit None
  ! Argument list:
  Type(Coords_),        Intent(In)    :: Coords
  Type(Grids_),         Intent(In)    :: Grids
  Type (NWPFlow_),      Intent(In)    :: NWPFlow
  Type (ShortTime_),    Intent(In)    :: Time
  Type (Position_),     Intent(In)    :: Position
  Type(NWPFlowMemory_), Intent(InOut) :: NWPFlowMemory
  Type (Plant_),        Intent(Out)   :: Plant
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! NWPFlow       :: State of an NWP flow module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in fixed met cases.
  ! Position      :: Coords of the point in various coord systems in Coords, with flags to indicate whether
  !                  the values are valid.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Plant         :: Plant information.
  ! Locals:
  Real(Std) :: X      !} Horizontal location in coord system of horizontal grids used for met data.
  Real(Std) :: Y      !}
  Real(Std) :: XAncil !] Horizontal location in coord system of horizontal grids used for ancillary data.
  Real(Std) :: YAncil !]
  Real(Std) :: Eta    ! Height in coord system of vertical grids used for met data.

  X      = Position%XY(1, NWPFlow%C%iHCoords(1))
  Y      = Position%XY(2, NWPFlow%C%iHCoords(1))
  XAncil = Position%XY(1, NWPFlow%AncillaryMetLAI%iHCoord)
  YAncil = Position%XY(2, NWPFlow%AncillaryMetLAI%iHCoord)

  Call PlantInfoFromMet(                   &
         Grids,                            &
         NWPFlow%NWPMet,                   &
         NWPFlow%NWPMetCanopyHeight,       &
         NWPFlow%NWPMetCanopy,             &
         NWPFlow%AncillaryMetCanopyHeight, &
         NWPFlow%AncillaryMetLAI,          &
         NWPFlow%NWPCanopyHeight,          &
         Time, X, Y, XAncil, YAncil,       &
         Plant,                            &
         NWPFlowMemory                     &
       )

End Subroutine GetNWPPlant

!-------------------------------------------------------------------------------------------------------------

Subroutine NWPReflectCoordIndices(NWPFlow, nHCoords, iHCoords, nZCoords, iZCoords)
! Returns indices of coord systems in which positions need to be specified when using NWPReflect.

  Implicit None
  ! Argument list:
  Type(NWPFlow_), Intent(In)  :: NWPFlow              ! State of an NWP flow module instance.
  Integer,        Intent(Out) :: nHCoords             !} Number and values of indices of horizontal and
  Integer,        Intent(Out) :: iHCoords(MaxHCoords) !} vertical coord systems.
  Integer,        Intent(Out) :: nZCoords             !}
  Integer,        Intent(Out) :: iZCoords(MaxHCoords) !}

  ! Horizontal coords needed: none
  nHCoords    = 0
  iHCoords(1) = 0 ! avoids compiler warning

  ! Vertical coords needed:
  ! 1) m agl.
  nZCoords    = 1
  iZCoords(1) = NWPFlow%C%iZCoords(2)

End Subroutine NWPReflectCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine NWPReflect(                      &
             Coords, NWPFlow,               &
             Vel, OldPosition, Position, U, &
             Reflected                      &
           )
! Reflects position of particle or puff centroid using an NWP flow module instance.

  Implicit None
  ! Argument list:
  Type(Coords_),   Intent(In)    :: Coords
  Type(NWPFlow_),  Intent(In)    :: NWPFlow
  Logical,         Intent(In)    :: Vel
  Type(Position_), Intent(In)    :: OldPosition
  Type(Position_), Intent(InOut) :: Position
  Real(Std),       Intent(InOut) :: U(3)
  Logical,         Intent(Out)   :: Reflected
  ! Coords      :: Collection of coord systems.
  ! NWPFlow     :: State of an NWP flow module instance.
  ! Vel         :: Indicates a dispersion model with velocity memory is being used.
  ! OldPosition :: Old particle position.
  ! Position    :: Current particle position.
  ! U           :: Current particle velocity. Defined only if Vel is true.
  ! Reflected   :: Indicates a reflection has occurred.
  ! Locals:
  Integer :: iZCoordZ ! Index of m agl coord system.

  iZCoordZ = NWPFlow%C%iZCoords(2)

  Reflected = .false.

  ! Velocity memory.
  If (Vel) Then

    If (Position%Z(iZCoordZ) < 0.0) Then
      U(3)                               = - U(3)
      Position%Z(iZCoordZ)               = - Position%Z(iZCoordZ)
      Position%ZValid(1:Coords%nZCoords) = .false.
      Position%ZValid(iZCoordZ)          = .true.
      Position%RhoValid                  = .false.
      Reflected                          = .true.
    End If

  ! No velocity memory.
  Else

    If (Position%Z(iZCoordZ) < 0.0) Then
      Position%Z(iZCoordZ)               = - Position%Z(iZCoordZ)
      Position%ZValid(1:Coords%nZCoords) = .false.
      Position%ZValid(iZCoordZ)          = .true.
      Position%RhoValid                  = .false.
      Reflected                          = .true.
    End If

  End If

End Subroutine NWPReflect

!-------------------------------------------------------------------------------------------------------------

Subroutine FlowInfoFromMet(         &
             Coords, Grids,         &
             NWPFlow, M,            &
             Moisture, FastRun,     &
             Time, X, Y, Eta, ZAgl, &
             Flow, ProfileData,     &
             NWPFlowMemory          &
           )
! Extracts flow information at a particular location from the met fields.

  Implicit None
  ! Argument List:
  Type(Coords_),        Intent(In),    Target :: Coords
  Type(Grids_),         Intent(In),    Target :: Grids
  Type(NWPFlow_),       Intent(In)            :: NWPFlow ! $$ note parts of NWPFlow (M) are
                                                         ! also passed - rationalise
  Type(NWPMet_),        Intent(In)            :: M
  Logical,              Intent(In)            :: Moisture
  Logical,              Intent(In)            :: FastRun
  Type(ShortTime_),     Intent(In)            :: Time
  Real(Std),            Intent(In)            :: X
  Real(Std),            Intent(In)            :: Y
  Real(Std),            Intent(In)            :: Eta
  Real(Std),            Intent(In)            :: ZAgl
  Type(Flow_),          Intent(Out)           :: Flow
  Type(ProfileData_),   Intent(Out)           :: ProfileData
  Type(NWPFlowMemory_), Intent(InOut), Target :: NWPFlowMemory
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! M             :: State of an NWP met module instance.
  ! Moisture      :: Indicates Q is required.
  ! FastRun       :: Indicates fast, but less precise, calculation required.
  ! Time          :: Current time in non-fixed met cases and the time of the met in
  !                  fixed met cases.
  ! X             :} Horizontal location in coord system of horizontal grids used for met data.
  ! Y             :}
  ! Eta           :: Height in coord system of vertical grid used for met data.
  ! ZAgl          :: Height (m agl).
  ! Flow          :: Flow information.
  ! ProfileData   :: Data needed to construct idealised analytic mean flow and
  !                  turbulence profiles.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Locals:
  Real(Std)               :: HMax         ! Max metric coefficient.
  Real(Std)               :: dX           !} Grid separation in longitudinal and
  Real(Std)               :: dY           !} latitudinal direction in metres.
  Real(Std)               :: dZ           ! Height change between levels.
  Real(Std)               :: dZdX(2)      ! Horizontal height gradient (at constant
                                          ! Eta).
  Real(Std)               :: U0           ! X component of near surface wind.
  Real(Std)               :: V0           ! Y component of near surface wind.
  Real(Std)               :: RhoS         ! Surface density.
  Real(Std)               :: HeatFlux     ! Surface sensible heat flux.
  Real(Std)               :: ZAgl2        ! Height agl (m) at second level at point of interest.
  Type(HCoord_),  Pointer :: HCoord       !} Abbreviations for coords and grids.
  Type(ZCoord_),  Pointer :: ZCoord       !}
  Type(HGrid_),   Pointer :: HGrid        !}
  Type(HGrid_),   Pointer :: HGridU       !}
  Type(HGrid_),   Pointer :: HGridV       !}
  Type(ZGrid_),   Pointer :: ZGrid        !}
  Type(ZGrid_),   Pointer :: ZGridUV      !}
  Type(ZGrid_),   Pointer :: ZGridW       !}
  Type(HCoeffs_), Pointer ::  HCoeffs     !] Interpolation coefficients for various grids (including for
  Type(HCoeffs_)          :: dHCoeffsdX   !] obtaining derivatives of fields and for interpolating to a time
  Type(HCoeffs_)          :: dHCoeffsdY   !] at the end of the NWP met time interval). The pointers are
  Type(HCoeffs_)          ::  HCoeffsU    !] abbreviations for sets of interpolation coefficients in
  Type(HCoeffs_)          ::  HCoeffsV    !] FlowMemory.
  Type(ZCoeffs_), Pointer ::  ZCoeffs     !]
  Type(ZCoeffs_)          :: dZCoeffsdZ   !]
  Type(ZCoeffs_)          ::  ZCoeffsUV   !]
  Type(ZCoeffs_)          ::  ZCoeffsW    !]
  Type(TCoeffs_), Pointer ::  TCoeffs     !]
  Type(TCoeffs_), Pointer ::  TCoeffsEnd  !]
  Real(Std)               :: Temp(2)      ! Temporary array needed to avoid array-copying warnings on Linux
                                          ! Intel compiler with -C option.
                                          ! $$ review whether still needed.
  ! Locals for urban modifications: $$ Move into subroutine?
  Real(Std)      :: LandUseFracs(9)
  Real(Std)      :: CanopyHeight
  Real(Std)      :: LambdaF
  Real(Std)      :: LambdaP
  Real(Std)      :: D
  Real(Std)      :: CanopyFactor
  Real(Std)      :: CanopyFactor1
  Real(Std)      :: CanopyFactor2
  Real(Std)      :: LExp
  Real(Std)      :: ZHat
  Real(Std)      :: UStarG
  Real(Std)      :: UStarGByUVK
  Real(Std)      :: Z0G
  Type(ZCoeffs_) :: ZCoeffsUVCanopyTop ! $$ put in flow memory?
  Real(Std), Parameter :: CanopyDepthFactor = 3.0
  Integer        :: i

  ! Notes on derivatives: All differences are evaluated as 2 - 1, with the denominator
  ! dX, dY, dZ or dEta chosen consistently.

  ! 1) Set up abbreviations for coords, grids and interpolation coefficients.

  HCoord     => Coords%HCoords(M%iHCoord)
  ZCoord     => Coords%ZCoords(M%iZCoord)
  HGrid      => Grids%HGrids(M%iHGrid  )
  HGridU     => Grids%HGrids(M%iHGridU )
  HGridV     => Grids%HGrids(M%iHGridV )
  ZGrid      => Grids%ZGrids(M%iZGrid  )
  ZGridUV    => Grids%ZGrids(M%iZGridUV)
  ZGridW     => Grids%ZGrids(M%iZGridW )
  HCoeffs    => NWPFlowMemory%HCoeffs
  ZCoeffs    => NWPFlowMemory%ZCoeffs
  TCoeffs    => NWPFlowMemory%TCoeffs
  TCoeffsEnd => NWPFlowMemory%TCoeffsEnd

  ! 2) Calculate interpolation coefficients.

  ! Calculate temporal interpolation coefficients.
  If (.not. TCoeffs%Valid) Then
    Call GetTCoeffs(                                 &
           Time,                                     &
           M%OldTime, Int(ShortTime2RealTime(M%Dt)), &
           M%Old, M%New,                             &
           TCoeffs                                   &
         )
  End If

  ! Calculate time interpolation coefficients for interpolating to a time at the end of the NWP met time
  ! interval (note end here means in chronological time not run direction time - hence the use of
  ! IsBackwards). Note there's no need to check TCoeffs%Valid because of the code immediately above. Also
  ! there is no need to set TCoeffsEnd%Valid because this will be done through TCoeffsEnd = TCoeffs. The "If
  ! NextHeatFlux" is for efficiency as TCoeffsEnd is not used otherwise.
  If (M%MetDefn%NextHeatFlux) Then
    If (.not. TCoeffsEnd%Valid) Then
      TCoeffsEnd = TCoeffs
      If (IsBackwards()) Then
        TCoeffsEnd%T1     = 0.0
        TCoeffsEnd%T2     = 1.0
        TCoeffsEnd%iTNear = M%Old
      Else
        TCoeffsEnd%T1     = 1.0
        TCoeffsEnd%T2     = 0.0
        TCoeffsEnd%iTNear = M%New
      End If
    End If
  End If

  ! Calculate horizontal interpolation coefficients.
  If (.not. HCoeffs%Valid) Then
    Call GetHCoeffs(X, Y, HGrid, HCoeffs)
  End If
  dHCoeffsdX    =  HCoeffs
  dHCoeffsdY    =  HCoeffs
  dHCoeffsdX%X1 =  1.0
  dHCoeffsdX%X2 = -1.0
  dHCoeffsdY%Y1 =  1.0
  dHCoeffsdY%Y2 = -1.0
  Call GetHCoeffs(X, Y, HGridU, HCoeffsU)
  Call GetHCoeffs(X, Y, HGridV, HCoeffsV)

  ! Calculate vertical interpolation coefficients.
  If (.not. ZCoeffs%Valid) Then
    Call GetZCoeffs(Eta, ZGrid, ZCoeffs)
  End If
  dZCoeffsdZ    =  ZCoeffs
  dZCoeffsdZ%Z1 =  1.0
  dZCoeffsdZ%Z2 = -1.0
  Call GetZCoeffs(Eta, ZGridW, ZCoeffsW)

  ! Log vertical interpolation coefficients for horizontal wind velocity.
  Call InterpXYT(HCoeffs, TCoeffs, M%Z0,         ProfileData%Z0)
  Call InterpXYT(HCoeffs, TCoeffs, M%H,          ProfileData%H )
  Call InterpXYT(HCoeffs, TCoeffs, M%Z(:,:,2,:), ZAgl2         )
  Call GetZCoeffs(                                                                    &
         Eta, ZGridUV, ZCoeffsUV,                                                     &
         ProfileData%Z0, ProfileData%H, ZAgl2 / (ZGrid%Z(2) - ZGrid%Z(1)), ZGrid%Z(1) &
       )

  ! 3) Extract data from arrays.

  ! 3-D data: U, V, W, T, Q, P.
  Call InterpXYZT(HCoeffsU, ZCoeffsUV, TCoeffs, M%U, Flow%U(1))
  Call InterpXYZT(HCoeffsV, ZCoeffsUV, TCoeffs, M%V, Flow%U(2))
  Call InterpXYZT(HCoeffs,  ZCoeffsW,  TCoeffs, M%W, Flow%U(3))
  If (EulerianModelGlobal) Call InterpXYZT(HCoeffs,  ZCoeffsW,  TCoeffs, M%EtaDot, Flow%EtaDot) ! $$
  Call InterpXYZT(HCoeffs,  ZCoeffs,   TCoeffs, M%T, Flow%T   )
  If (Moisture) Then
    Call InterpXYZT(HCoeffs, ZCoeffs, TCoeffs, M%Q, Flow%Q)
  End If
  Call InterpXYZT(HCoeffs, ZCoeffs, TCoeffs, M%P, Flow%P)

  ! Block for urban modification. $$ Move into subroutine?

  ! Ensure canopy not used unless required.
  ProfileData%Canopy = .false.

  If (NWPFlow%iLandUse /= 0 .and. NWPFlow%UrbanCanopy) Then

    LandUseFracs(:) = NWPFlow%AncillaryMetLandUse%LandUseFracs(      &
                             HCoeffs%iXNear, HCoeffs%iYNear, :, 1) ! $$ check '1' OK in general

    If (LandUseFracs(6) > 0.05) Then ! $$ Urban tile is LandUseFracs(6)

      LambdaP      =   22.88 * LandUseFracs(6)**6 &  ! $$ could make this more efficient
                     - 59.47 * LandUseFracs(6)**5 &
                     + 57.75 * LandUseFracs(6)**4 &
                     - 25.11 * LandUseFracs(6)**3 &
                     +  4.33 * LandUseFracs(6)**2 &
                     +  0.19 * LandUseFracs(6)
      LambdaF      =   16.41 * LandUseFracs(6)**6 &
                     - 41.86 * LandUseFracs(6)**5 &
                     + 40.39 * LandUseFracs(6)**4 &
                     - 17.76 * LandUseFracs(6)**3 &
                     +  3.24 * LandUseFracs(6)**2 &
                     +  0.06 * LandUseFracs(6)
      CanopyHeight =   167.409 * LandUseFracs(6)**5 &
                     - 337.853 * LandUseFracs(6)**4 &
                     + 247.813 * LandUseFracs(6)**3 &
                     - 76.3678 * LandUseFracs(6)**2 &
                     + 11.4832 * LandUseFracs(6)    &
                     + 4.48226
      ! Traps for negative values - probably unnecessary but needs investigation $$
      If (LambdaP      <= 0.0) LambdaP      = 0.01
      If (LambdaF      <= 0.0) LambdaF      = 0.01
      If (CanopyHeight <= 0.0) CanopyHeight = 1.0

      ! Displacement height D.
      D = Sqrt(15.0 * LambdaF)
      If (D < 0.01) Then
        ! Use approximation D = D/2 - D^2/6 + D^3/24. 
        ! For D < 0.01 the approximation is accurate to 7 significant figures.
        D = D * (1.0/2.0 + D * (- 1.0/6.0 + D/24.0))
      Else
        D = 1.0 - (1.0 - Exp(- D)) / D 
      End If
      D = CanopyHeight * D

      Z0G = 0.1

      LExp = CanopyHeight / (9.6 * LambdaF)

      ProfileData%H = Max(ProfileData%H, CanopyDepthFactor * CanopyHeight)

      ! Set canopy elements of ProfileData. Note UAtHC, LExp, ZHat, UStarG, Z0G and LC only used if point is 
      ! inside canopy.
      ProfileData%Canopy  = .true.
      ProfileData%LambdaP = LambdaP
      ProfileData%D       = D
      ProfileData%HC      = CanopyHeight
      ProfileData%LExp    = LExp
      ProfileData%Z0G     = Z0G
      ProfileData%LC      = 100.0

      If (ZAgl < CanopyDepthFactor * CanopyHeight) Then

        If (ZAgl >= CanopyHeight) Then

          CanopyFactor  = Log((ZAgl                           - D) / ProfileData%Z0) /            &
                          Log((ZAgl                           + ProfileData%Z0) / ProfileData%Z0)
          CanopyFactor1 = Log((CanopyHeight                   - D) / ProfileData%Z0) /            &
                          Log((CanopyHeight                   + ProfileData%Z0) / ProfileData%Z0)
          CanopyFactor2 = Log((CanopyDepthFactor*CanopyHeight - D) / ProfileData%Z0) /            &
                          Log((CanopyDepthFactor*CanopyHeight + ProfileData%Z0) / ProfileData%Z0)
          CanopyFactor = CanopyFactor1 +                                          &
                         (CanopyFactor - CanopyFactor1) * (1.0 - CanopyFactor1) / &
                         (CanopyFactor2 - CanopyFactor1)
          Flow%U(1) = Flow%U(1) * CanopyFactor
          Flow%U(2) = Flow%U(2) * CanopyFactor

          ! Set other canopy elements of ProfileData in case of precision issues on the ZAgl >= CanopyHeight 
          ! test.
          ProfileData%UAtHC  = Sqrt(Flow%U(1)**2 + Flow%U(2)**2)
          ProfileData%ZHat   = 0.0
          ProfileData%UStarG = 0.0

        Else If (ZAgl < CanopyHeight) Then

          ! Get canopy top wind.
          Call GetZCoeffs(                                                                    &
                 (CanopyHeight/ZAgl2) * (ZGrid%Z(2) - ZGrid%Z(1)),                            &
                 ZGridUV, ZCoeffsUVCanopyTop,                                                 &
                 ProfileData%Z0, ProfileData%H, ZAgl2 / (ZGrid%Z(2) - ZGrid%Z(1)), ZGrid%Z(1) &
               )
          Call InterpXYZT(HCoeffsU, ZCoeffsUVCanopyTop, TCoeffs, M%U, Flow%U(1))
          Call InterpXYZT(HCoeffsV, ZCoeffsUVCanopyTop, TCoeffs, M%V, Flow%U(2))
          CanopyFactor = Log((CanopyHeight     - D) / ProfileData%Z0) /            &
                         Log((CanopyHeight     + ProfileData%Z0) / ProfileData%Z0)
          Flow%U(1) = Flow%U(1) * CanopyFactor
          Flow%U(2) = Flow%U(2) * CanopyFactor

          ! Solve LExp = (ZHat + Z0G) Log ((ZHat + Z0G) / Z0G).
          ZHat = Max(LExp, Z0G * Exp(1.0)) - Z0G ! Start larger than solution.
          Do i = 1, 3 ! $$ not sure really need 3 of these? Definitely need 1.
            ZHat = (LExp + ZHat + Z0G) / (Log((ZHat + Z0G)/Z0G) + 1.0) - Z0G
          End Do

          ! Limit ZHat to canopy height.
          ZHat = Min(ZHat, CanopyHeight)

          ! ZHat <= 0.0 part is just a precaution - shouldn't be needed with exact arithmetic.
          If (ZHat <= 0.0) Then
            UStarG = 0.0
          Else
            UStarGByUVK = Exp(- (CanopyHeight - ZHat) / LExp) / Log((ZHat + Z0G) / Z0G) 
            UStarG      = Sqrt(Flow%U(1)**2 + Flow%U(2)**2) * UStarGByUVK *VK
          End If

          If (ZAgl > ZHat) Then
            CanopyFactor = Exp(- (CanopyHeight - ZAgl) / LExp)
          Else
            CanopyFactor = UStarGByUVK * Log((ZAgl + Z0G) / Z0G)
          End If

          ! Set other canopy elements of ProfileData.
          ProfileData%UAtHC  = Sqrt(Flow%U(1)**2 + Flow%U(2)**2)
          ProfileData%ZHat   = ZHat
          ProfileData%UStarG = UStarG

          Flow%U(1) = Flow%U(1) * CanopyFactor
          Flow%U(2) = Flow%U(2) * CanopyFactor

        End If

      End If

    End If

  End If


  ! 3-D data: Rho gradients.
  If (FastRun) Then

    dZ =                                                                 &
      M%Z(HCoeffs%iXNear, HCoeffs%iYNear, ZCoeffs%iZ2, TCoeffs%iTNear) - &
      M%Z(HCoeffs%iXNear, HCoeffs%iYNear, ZCoeffs%iZ1, TCoeffs%iTNear)
    dZdX(1) =                                                            &
      M%Z(HCoeffs%iX2, HCoeffs%iYNear, ZCoeffs%iZNear, TCoeffs%iTNear) - &
      M%Z(HCoeffs%iX1, HCoeffs%iYNear, ZCoeffs%iZNear, TCoeffs%iTNear)
    dZdX(2) =                                                            &
      M%Z(HCoeffs%iXNear, HCoeffs%iY2, ZCoeffs%iZNear, TCoeffs%iTNear) - &
      M%Z(HCoeffs%iXNear, HCoeffs%iY1, ZCoeffs%iZNear, TCoeffs%iTNear)
    Flow%dRhodX(1) =                                                       &
      M%Rho(HCoeffs%iX2, HCoeffs%iYNear, ZCoeffs%iZNear, TCoeffs%iTNear) - &
      M%Rho(HCoeffs%iX1, HCoeffs%iYNear, ZCoeffs%iZNear, TCoeffs%iTNear)
    Flow%dRhodX(2) =                                                       &
      M%Rho(HCoeffs%iXNear, HCoeffs%iY2, ZCoeffs%iZNear, TCoeffs%iTNear) - &
      M%Rho(HCoeffs%iXNear, HCoeffs%iY1, ZCoeffs%iZNear, TCoeffs%iTNear)
    Flow%dRhodX(3) =                                                       &
      M%Rho(HCoeffs%iXNear, HCoeffs%iYNear, ZCoeffs%iZ2, TCoeffs%iTNear) - &
      M%Rho(HCoeffs%iXNear, HCoeffs%iYNear, ZCoeffs%iZ1, TCoeffs%iTNear)

  Else

    Call InterpXYZT(HCoeffs,    dZCoeffsdZ, TCoeffs, M%Z,   dZ            )
    Call InterpXYZT(dHCoeffsdX,    ZCoeffs, TCoeffs, M%Z,   dZdX(1)       )
    Call InterpXYZT(dHCoeffsdY,    ZCoeffs, TCoeffs, M%Z,   dZdX(2)       )
    Call InterpXYZT(dHCoeffsdX,    ZCoeffs, TCoeffs, M%Rho, Flow%dRhodX(1))
    Call InterpXYZT(dHCoeffsdY,    ZCoeffs, TCoeffs, M%Rho, Flow%dRhodX(2))
    Call InterpXYZT(HCoeffs,    dZCoeffsdZ, TCoeffs, M%Rho, Flow%dRhodX(3))

  End If

  ! 2-D data: U0, V0, T0, PS, Z0, FL, UStar, HeatFlux, H and Topog.
  If (FastRun) Then
    U0                = M%U         (HCoeffsU%iXNear, HCoeffsU%iYNear, 2, TCoeffs%iTNear)
    V0                = M%V         (HCoeffsV%iXNear, HCoeffsV%iYNear, 2, TCoeffs%iTNear)
    Flow%T0           = M%T         (HCoeffs%iXNear,  HCoeffs%iYNear,  2, TCoeffs%iTNear)
    Flow%PS           = M%P         (HCoeffs%iXNear,  HCoeffs%iYNear,  1, TCoeffs%iTNear)
    Flow%PSeaLevel    = M%PSeaLevel (HCoeffs%iXNear,  HCoeffs%iYNear,     TCoeffs%iTNear)
    ProfileData%Z0    = M%Z0        (HCoeffs%iXNear,  HCoeffs%iYNear,     TCoeffs%iTNear)
    Flow%FLPa         = M%FLPa      (HCoeffs%iXNear,  HCoeffs%iYNear,     TCoeffs%iTNear)
    If (M%MetDefn%NextHeatFlux) Then
      ProfileData%UStar = M%UStar   (HCoeffs%iXNear, HCoeffs%iYNear, TCoeffsEnd%iTNear)
      HeatFlux          = M%HeatFlux(HCoeffs%iXNear, HCoeffs%iYNear, TCoeffsEnd%iTNear)
    Else
      ProfileData%UStar = M%UStar   (HCoeffs%iXNear, HCoeffs%iYNear, TCoeffs%iTNear)
      HeatFlux          = M%HeatFlux(HCoeffs%iXNear, HCoeffs%iYNear, TCoeffs%iTNear)
    End If
    ProfileData%H     = M%H    (HCoeffs%iXNear, HCoeffs%iYNear, TCoeffs%iTNear)
    Flow%Topog        = M%Topog(HCoeffs%iXNear, HCoeffs%iYNear, TCoeffs%iTNear)
    Flow%dTopogdX(1)  = M%Topog(HCoeffs%iX2,    HCoeffs%iYNear, TCoeffs%iTNear) - &
                        M%Topog(HCoeffs%iX1,    HCoeffs%iYNear, TCoeffs%iTNear)
    Flow%dTopogdX(2)  = M%Topog(HCoeffs%iXNear, HCoeffs%iY2,    TCoeffs%iTNear) - &
                        M%Topog(HCoeffs%iXNear, HCoeffs%iY1,    TCoeffs%iTNear)
  Else
    Call InterpXYT( HCoeffsU,  TCoeffs,  M%U(:, :, 2, :), U0               )
    Call InterpXYT( HCoeffsV,  TCoeffs,  M%V(:, :, 2, :), V0               )
    Call InterpXYT( HCoeffs,   TCoeffs,  M%T(:, :, 2, :), Flow%T0          )
    Call InterpXYT( HCoeffs,   TCoeffs,  M%P(:, :, 1, :), Flow%PS          )
    Call InterpXYT( HCoeffs,   TCoeffs,  M%PSeaLevel,     Flow%PSeaLevel   )
    Call InterpXYT( HCoeffs,   TCoeffs,  M%Z0,            ProfileData%Z0   )
    Call InterpXYT( HCoeffs,   TCoeffs,  M%FLPa,          Flow%FLPa        )
    If (M%MetDefn%NextHeatFlux) Then
      Call InterpXYT( HCoeffs,   TCoeffsEnd, M%UStar,         ProfileData%UStar)
      Call InterpXYT( HCoeffs,   TCoeffsEnd, M%HeatFlux,      HeatFlux         )
    Else
      Call InterpXYT( HCoeffs,   TCoeffs, M%UStar,         ProfileData%UStar)
      Call InterpXYT( HCoeffs,   TCoeffs, M%HeatFlux,      HeatFlux         )
    End If
    Call InterpXYT( HCoeffs,   TCoeffs,  M%H,             ProfileData%H    )
    Call InterpXYT( HCoeffs,   TCoeffs,  M%Topog,         Flow%Topog       )
    Call InterpXYT(dHCoeffsdX, TCoeffs,  M%Topog,         Flow%dTopogdX(1) )
    Call InterpXYT(dHCoeffsdY, TCoeffs,  M%Topog,         Flow%dTopogdX(2) )
  End If

  ! Potential temperature and density.
  Flow%Theta = Flow%T * (PRef / Flow%P) ** (GasConstant/Cp)
  Flow%Rho   = Flow%P / (GasConstant * Flow%T)

  ! Surface wind direction, WT, WQ, Monin-Obukhov length and WStar.
  ProfileData%Phi0 = ATan2ZeroTest(V0, U0)
  RhoS = Flow%PS / (GasConstant * Flow%T0)
  Flow%WT = HeatFlux/(RhoS * Cp)
  Flow%WQ = 0 ! $$
  If (ProfileData%UStar == 0.0) Then
    ProfileData%RecipLMO = - Sign(200000.0, HeatFlux)
  Else
    ProfileData%RecipLMO = - (VK * Gravity * HeatFlux) /                  &
                           (RhoS * Cp * Flow%T0 * ProfileData%UStar ** 3)
  End If
  If (ProfileData%RecipLMO >= 0.0) Then
    ProfileData%WStar = 0.0
  Else
    ProfileData%WStar = (ProfileData%H * Gravity * HeatFlux / (RhoS * Cp * Flow%T0)) &
                        ** (1.0/3.0)
  End If

  ! Normalising gradients by dX, dY, dZ and correcting for change from X, Y, Eta to X,
  ! Y, Z coords.
  ! Note that (ignoring Y for simplicity) we have
  !     (i)  d / dZ |X  =  d / dEta |X  *  dEta / dZ |X
  !     (ii) d / dX |Z  =  d / dEta |X  *  dEta / dX |Z  +  d / dX |Eta * dX / dX |Z
  !                     =  d / dEta |X  *  dEta / dX |Z  +  d / dX |Eta.
  ! Putting Z in (ii) gives
  !                 0  =  dZ / dEta |X  *  dEta / dX |Z  +  dZ / dX |Eta
  ! and so
  !                 dEta / dX |Z  =  - [dZ / dX |Eta]  /  [dZ / dEta |X].
  ! Hence (ii) gives
  !   d / dX |Z  =  d / dX |Eta  -  d / dEta |X  *  [dZ / dX |Eta]  /  [dZ / dEta |X]
  !              =  d / dX |Eta  -  d / dZ |X  *  dZ / dX |Eta.

  Temp = (/ X, Y /)
  Call MetricCoeffs(HCoord, Temp, HMax, dX, dY)
  dX      = dX * HGrid%dX
  dY      = dY * HGrid%dY
  dZdX(1) = dZdX(1) / dX
  dZdX(2) = dZdX(2) / dY

  Flow%dTopogdX(1) = Flow%dTopogdX(1) / dX
  Flow%dTopogdX(2) = Flow%dTopogdX(2) / dY

  Flow%dRhodX(1) = Flow%dRhodX(1) / dX
  Flow%dRhodX(2) = Flow%dRhodX(2) / dY
  Flow%dRhodX(3) = Flow%dRhodX(3) / dZ
  Flow%dRhodX(1) = Flow%dRhodX(1) - Flow%dRhodX(3) * dZdX(1)
  Flow%dRhodX(2) = Flow%dRhodX(2) - Flow%dRhodX(3) * dZdX(2)

End Subroutine FlowInfoFromMet

!-------------------------------------------------------------------------------------------------------------

Subroutine CloudInfoFromMet(  &
             Grids,           &
             M,               &
             Time, X, Y, Eta, &
             Cloud,           &
             NWPFlowMemory    &
           )
! Extracts cloud information at a particular location from the met fields.

  Implicit None
  ! Argument List:
  Type(Grids_),         Intent(In),    Target :: Grids
  Type(NWPMet_),        Intent(In)            :: M
  Type(ShortTime_),     Intent(In)            :: Time
  Real(Std),            Intent(In)            :: X
  Real(Std),            Intent(In)            :: Y
  Real(Std),            Intent(In)            :: Eta
  Type(Cloud_),         Intent(Out)           :: Cloud
  Type(NWPFlowMemory_), Intent(InOut), Target :: NWPFlowMemory
  ! Grids         :: Collection of grids.
  ! M             :: State of an NWP met module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in
  !                  fixed met cases.
  ! X             :} Horizontal location in coord system of horizontal grids used for met data.
  ! Y             :}
  ! Eta           :: Height in coord system of vertical grid used for met data.
  ! Cloud         :: Cloud information.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Locals:
  Type(HGrid_),   Pointer :: HGrid      !} Abbreviations for grids and interpolation
  Type(ZGrid_),   Pointer :: ZGrid      !} coefficients.
  Type(HCoeffs_), Pointer :: HCoeffs    !}
  Type(ZCoeffs_), Pointer :: ZCoeffs    !}
  Type(TCoeffs_), Pointer :: TCoeffs    !}
  Type(TCoeffs_), Pointer :: TCoeffsEnd !}

  ! 1) Set up abbreviations for grids and interpolation coefficients.

  HGrid      => Grids%HGrids(M%iHGrid)
  ZGrid      => Grids%ZGrids(M%iZGrid)
  HCoeffs    => NWPFlowMemory%HCoeffs
  ZCoeffs    => NWPFlowMemory%ZCoeffs
  TCoeffs    => NWPFlowMemory%TCoeffs
  TCoeffsEnd => NWPFlowMemory%TCoeffsEnd

  ! 2) Calculate interpolation coefficients.

  ! Calculate temporal interpolation coefficients.
  If (.not. TCoeffs%Valid) Then
    Call GetTCoeffs(                                 &
           Time,                                     &
           M%OldTime, Int(ShortTime2RealTime(M%Dt)), &
           M%Old, M%New,                             &
           TCoeffs                                   &
         )
  End If

  ! Calculate time interpolation coefficients for interpolating to a time at the end of the NWP met time
  ! interval (note end here means in chronological time not run direction time - hence the use of
  ! IsBackwards). Note there's no need to check TCoeffs%Valid because of the code immediately above. Also
  ! there is no need to set TCoeffsEnd%Valid because this will be done through TCoeffsEnd = TCoeffs. The "If
  ! NextCloud" is for efficiency as TCoeffsEnd is not used otherwise.
  If (M%MetDefn%NextCloud) Then
    If (.not. TCoeffsEnd%Valid) Then
      TCoeffsEnd = TCoeffs
      If (IsBackwards()) Then
        TCoeffsEnd%T1     = 0.0
        TCoeffsEnd%T2     = 1.0
        TCoeffsEnd%iTNear = M%Old
      Else
        TCoeffsEnd%T1     = 1.0
        TCoeffsEnd%T2     = 0.0
        TCoeffsEnd%iTNear = M%New
      End If
    End If
  End If

  ! Calculate vertical interpolation coefficients.
  If (.not. ZCoeffs%Valid) Then
    Call GetZCoeffs(Eta, ZGrid, ZCoeffs)
  End If

  ! Calculate horizontal interpolation coefficients.
  If (.not. HCoeffs%Valid) Then
    Call GetHCoeffs(X, Y, HGrid, HCoeffs)
  End If

  ! 3) Extract data from arrays.

  ! 3-D data: Cloud3d, TotalOrDynCloudWater and TotalOrDynCloudIce.
   
  If (M%MetDefn%NextCloud) Then
    Call InterpXYZT(HCoeffs, ZCoeffs, TCoeffsEnd, M%Cloud3d,              Cloud%Cloud3d             )
    Call InterpXYZT(HCoeffs, ZCoeffs, TCoeffsEnd, M%TotalOrDynCloudWater, Cloud%TotalOrDynCloudWater)
    Call InterpXYZT(HCoeffs, ZCoeffs, TCoeffsEnd, M%TotalOrDynCloudIce,   Cloud%TotalOrDynCloudIce  )
  Else
    Call InterpXYZT(HCoeffs, ZCoeffs, TCoeffs, M%Cloud3d,              Cloud%Cloud3d             )
    Call InterpXYZT(HCoeffs, ZCoeffs, TCoeffs, M%TotalOrDynCloudWater, Cloud%TotalOrDynCloudWater)
    Call InterpXYZT(HCoeffs, ZCoeffs, TCoeffs, M%TotalOrDynCloudIce,   Cloud%TotalOrDynCloudIce  )
  End If


  ! 2-D data: ConCloud, ConCloudBase, ConCloudTop, Cloud, TotalOrDynCloudBase and
  ! TotalOrDynCloudTop.
  If (M%MetDefn%NextCloud) Then
    Call InterpXYT (HCoeffs, TCoeffsEnd, M%ConCloud,     Cloud%ConCloud    )
    Call InterpCXYT(HCoeffs, TCoeffsEnd, M%ConCloudBase, Cloud%ConCloudBase)
    Call InterpCXYT(HCoeffs, TCoeffsEnd, M%ConCloudTop,  Cloud%ConCloudTop )
    Call InterpCXYT(HCoeffs, TCoeffsEnd, M%ConCloudBasePa, Cloud%ConCloudBasePa)
    Call InterpCXYT(HCoeffs, TCoeffsEnd, M%ConCloudTopPa,  Cloud%ConCloudTopPa )

    Call InterpXYT (HCoeffs, TCoeffsEnd, M%Cloud,               Cloud%Cloud              )
    Call InterpCXYT(HCoeffs, TCoeffsEnd, M%TotalOrDynCloudBase, Cloud%TotalOrDynCloudBase)
    Call InterpCXYT(HCoeffs, TCoeffsEnd, M%TotalOrDynCloudTop,  Cloud%TotalOrDynCloudTop )
  Else
    Call InterpXYT (HCoeffs, TCoeffs, M%ConCloud,     Cloud%ConCloud    )
    Call InterpCXYT(HCoeffs, TCoeffs, M%ConCloudBase, Cloud%ConCloudBase)
    Call InterpCXYT(HCoeffs, TCoeffs, M%ConCloudTop,  Cloud%ConCloudTop )
    Call InterpCXYT(HCoeffs, TCoeffs, M%ConCloudBasePa, Cloud%ConCloudBasePa)
    Call InterpCXYT(HCoeffs, TCoeffs, M%ConCloudTopPa,  Cloud%ConCloudTopPa )

    Call InterpXYT (HCoeffs, TCoeffs, M%Cloud,               Cloud%Cloud              )
    Call InterpCXYT(HCoeffs, TCoeffs, M%TotalOrDynCloudBase, Cloud%TotalOrDynCloudBase)
    Call InterpCXYT(HCoeffs, TCoeffs, M%TotalOrDynCloudTop,  Cloud%TotalOrDynCloudTop )
  End If

  Cloud%TotalCloudFlag = M%TotalCloudFlag

End Subroutine CloudInfoFromMet

!-------------------------------------------------------------------------------------------------------------

Subroutine RainInfoFromMet(   &
             Grids,           &
             M,               &
             Time, X, Y, Eta, &
             Rain,            &
             NWPFlowMemory    &
           )
! Extracts rain information at a particular location from the met fields.

  Implicit None
  ! Argument List:
  Type(Grids_),         Intent(In),    Target :: Grids
  Type(NWPMet_),        Intent(In)            :: M
  Type(ShortTime_),     Intent(In)            :: Time
  Real(Std),            Intent(In)            :: X
  Real(Std),            Intent(In)            :: Y
  Real(Std),            Intent(In)            :: Eta
  Type(Rain_),          Intent(Out)           :: Rain
  Type(NWPFlowMemory_), Intent(InOut), Target :: NWPFlowMemory
  ! Grids         :: Collection of grids.
  ! M             :: State of an NWP met module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in
  !                  fixed met cases.
  ! X             :} Horizontal location in coord system of horizontal grids used for met data.
  ! Y             :}
  ! Eta           :: Height in coord system of vertical grid used for met data.
  ! Rain          :: Rain information.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Locals:
  Type(HGrid_),   Pointer :: HGrid      !} Abbreviations for grids and interpolation
  Type(ZGrid_),   Pointer :: ZGrid      !} coefficients.
  Type(HCoeffs_), Pointer :: HCoeffs    !}
  Type(ZCoeffs_), Pointer :: ZCoeffs    !}
  Type(TCoeffs_), Pointer :: TCoeffs    !}
  Type(TCoeffs_), Pointer :: TCoeffsEnd !}

  ! 1) Set up abbreviations for grids and interpolation coefficients.

  HGrid      => Grids%HGrids(M%iHGrid)
  ZGrid      => Grids%ZGrids(M%iZGrid)
  HCoeffs    => NWPFlowMemory%HCoeffs
  ZCoeffs    => NWPFlowMemory%ZCoeffs
  TCoeffs    => NWPFlowMemory%TCoeffs
  TCoeffsEnd => NWPFlowMemory%TCoeffsEnd

  ! 2) Calculate interpolation coefficients.

  ! Calculate temporal interpolation coefficients.
  If (.not. TCoeffs%Valid) Then
    Call GetTCoeffs(                                 &
           Time,                                     &
           M%OldTime, Int(ShortTime2RealTime(M%Dt)), &
           M%Old, M%New,                             &
           TCoeffs                                   &
         )
  End If

  ! Calculate time interpolation coefficients for interpolating to a time at the end of the NWP met time
  ! interval (note end here means in chronological time not run direction time - hence the use of
  ! IsBackwards). Note there's no need to check TCoeffs%Valid because of the code immediately above. Also
  ! there is no need to set TCoeffsEnd%Valid because this will be done through TCoeffsEnd = TCoeffs. The "If
  ! NextPrecip" is for efficiency as TCoeffsEnd is not used otherwise.
  If (M%MetDefn%NextPrecip) Then
    If (.not. TCoeffsEnd%Valid) Then
      TCoeffsEnd = TCoeffs
      If (IsBackwards()) Then
        TCoeffsEnd%T1     = 0.0
        TCoeffsEnd%T2     = 1.0
        TCoeffsEnd%iTNear = M%Old
      Else
        TCoeffsEnd%T1     = 1.0
        TCoeffsEnd%T2     = 0.0
        TCoeffsEnd%iTNear = M%New
      End If
    End If
  End If

  ! Calculate vertical interpolation coefficients.
  If (.not. ZCoeffs%Valid) Then
    Call GetZCoeffs(Eta, ZGrid, ZCoeffs)
  End If

  ! Calculate horizontal interpolation coefficients.
  If (.not. HCoeffs%Valid) Then
    Call GetHCoeffs(X, Y, HGrid, HCoeffs)
  End If

  ! 3) Extract data from arrays.

  ! 2-D data: DynPpt and ConPpt.
  If (M%MetDefn%NextPrecip) Then
    Call InterpXYT(HCoeffs, TCoeffsEnd, M%ConPpt, Rain%ConPpt)
    Call InterpXYT(HCoeffs, TCoeffsEnd, M%DynPpt, Rain%DynPpt)
  Else
    Call InterpXYT(HCoeffs, TCoeffs, M%ConPpt, Rain%ConPpt)
    Call InterpXYT(HCoeffs, TCoeffs, M%DynPpt, Rain%DynPpt)
  End If

End Subroutine RainInfoFromMet

!-------------------------------------------------------------------------------------------------------------

Subroutine SurfaceInfoFromMet(           &
             Grids,                      &
             NWPMet,                     &
             NWPMetSoilMoisture,         &
             AncillaryMetSoilMoisture,   &
             AncillaryMetLandUse,        &
             NWPSoilMoisture,            &
             Time, X, Y, XAncil, YAncil, &
             Surface,                    &
             NWPFlowMemory               &
           )
! Extracts surface information at a particular location from the met fields.

  Implicit None
  ! Argument List:
  Type(Grids_),         Intent(In),    Target :: Grids
  Type(NWPMet_),        Intent(In)            :: NWPMet
  Type(NWPMet_),        Intent(In)            :: NWPMetSoilMoisture
  Type(AncillaryMet_),  Intent(In)            :: AncillaryMetSoilMoisture ! $$ add comments
  Type(AncillaryMet_),  Intent(In)            :: AncillaryMetLandUse
  Logical,              Intent(In)            :: NWPSoilMoisture
  Type(ShortTime_),     Intent(In)            :: Time
  Real(Std),            Intent(In)            :: X
  Real(Std),            Intent(In)            :: Y
  Real(Std),            Intent(In)            :: XAncil
  Real(Std),            Intent(In)            :: YAncil
  Type(Surface_),       Intent(Out)           :: Surface
  Type(NWPFlowMemory_), Intent(InOut), Target :: NWPFlowMemory
  ! Grids         :: Collection of grids.
  ! M             :: State of an NWP met module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in
  !                  fixed met cases.
  ! X             :} Horizontal location in coord system of horizontal grids used for met data.
  ! Y             :}
  ! XAncil        :] Horizontal location in coord system of horizontal grids used for ancillary data.
  ! YAncil        :]
  ! Surface       :: Surface information.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Locals:
  Type(HGrid_),   Pointer :: HGrid        !} Abbreviations for grids and interpolation
  Type(HGrid_),   Pointer :: HGridAncil   !} coefficients.
  Type(HCoeffs_), Pointer :: HCoeffs      !}
  Type(HCoeffs_)          :: HCoeffsAncil !}
  Type(TCoeffs_), Pointer :: TCoeffs      !}
  Type(TCoeffs_), Pointer :: TCoeffsEnd   !}
  Integer                 :: i            ! Loop index.

  ! $$ Treatment of grids here should be more systematic. For the moment we assume HGrid same
  ! for all ancillaries
  ! data sources and use TCoeffsEnd for Ancillary Met, where we assume interpolation is not required.

  ! 1) Set up abbreviations for grids and interpolation coefficients.

  HGrid      => Grids%HGrids(NWPMet%iHGrid)
  HGridAncil => Grids%HGrids(AncillaryMetLandUse%iHGrid)   ! $$ Assumes all Ancillaries on same grid
  HCoeffs    => NWPFlowMemory%HCoeffs
  TCoeffs    => NWPFlowMemory%TCoeffs
  TCoeffsEnd => NWPFlowMemory%TCoeffsEnd

  ! 2) Calculate interpolation coefficients.

  ! Note time interpolation irrelevant for monthly data and infinite dt data. Monthly data is regarded as
  ! representative of the month and infinite dt data is valid for all times.
  ! $$ Could possibly do something for monthly and also support yearly and "yearly and monthly" options.
  ! $$ Currently TCoeffs and TCoeffsEnd same as for the NWPMet data - should be different for each ancil
  ! module (with traps for monthly and dt infinite).
  ! $$ see also "$$" comments in AncillaryMet.

  ! Calculate temporal interpolation coefficients.
  If (.not. TCoeffs%Valid) Then
    Call GetTCoeffs(                                           &
           Time,                                               &
           NWPMet%OldTime, Int(ShortTime2RealTime(NWPMet%Dt)), &
           NWPMet%Old, NWPMet%New,                             &
           TCoeffs                                             &
         )
  End If

  ! Calculate time interpolation coefficients for interpolating to a time at the end of the NWP met time
  ! interval (note "end" here means in chronological time not run direction time - hence the use of
  ! IsBackwards). Note there's no need to check TCoeffs%Valid because of the code immediately above. Also
  ! there is no need to set TCoeffsEnd%Valid because this will be done through TCoeffsEnd = TCoeffs.
  If (.not. TCoeffsEnd%Valid) Then
    TCoeffsEnd = TCoeffs
    If (IsBackwards()) Then
      TCoeffsEnd%T1     = 0.0
      TCoeffsEnd%T2     = 1.0
      TCoeffsEnd%iTNear = NWPMet%Old
    Else
      TCoeffsEnd%T1     = 1.0
      TCoeffsEnd%T2     = 0.0
      TCoeffsEnd%iTNear = NWPMet%New
    End If
  End If

  ! Calculate horizontal interpolation coefficients.
  If (.not. HCoeffs%Valid) Then
    Call GetHCoeffs(X, Y, HGrid, HCoeffs)
  End If

  Call GetHCoeffs(XAncil, YAncil, HGridAncil, HCoeffsAncil) ! $$ avoid recalculation if HGrid = HGridAncil

  ! 3) Extract data from arrays.

  ! 2-D data: Land use fraction (9 land use types), met data. Note we use the nearest point rather than
  ! interpolating to avoid problems caused by the convention that negative values indicate sea.
  Do i = 1, 9
    Surface%LandUseFracs(i) = AncillaryMetLandUse%LandUseFracs(                                &
                                HCoeffsAncil%iXNear, HCoeffsAncil%iYNear, i, TCoeffsEnd%iTNear &
                              )
  End Do
  Surface%LandFrac          = AncillaryMetLandUse%LandFrac(                                    &
                                HCoeffsAncil%iXNear, HCoeffsAncil%iYNear, 1, TCoeffsEnd%iTNear &
                              )
  If (NWPSoilMoisture) Then
    Surface%SoilMoisture    = NWPMetSoilMoisture%SoilMoisture(                                 &
                                HCoeffs%iXNear,      HCoeffs%iYNear,         TCoeffs%iTNear    &
                              )
  Else
    Surface%SoilMoisture    = AncillaryMetSoilMoisture%SoilMoisture(                           &
                                HCoeffsAncil%iXNear, HCoeffsAncil%iYNear, 1, TCoeffsEnd%iTNear &
                              )
  End If


  ! Nonsense - delete $$
  If (.false.) Then
    Surface%LandUseFracs(8)   = AncillaryMetLandUse%LandUseFracs(                                    &
                                  HCoeffs%iXNear, Max(1, HCoeffs%iYNear + 250), 8, TCoeffsEnd%iTNear &
                                )
    If (NWPSoilMoisture) Then
      Surface%SoilMoisture    = NWPMetSoilMoisture%SoilMoisture(                       &
                                  HCoeffs%iXNear, HCoeffs%iYNear,    TCoeffs%iTNear    &
                                )
    Else
      Surface%SoilMoisture    = AncillaryMetSoilMoisture%SoilMoisture(                 &
                                  HCoeffs%iXNear, Max(1, HCoeffs%iYNear + 250), 1, TCoeffsEnd%iTNear &
                                )
    End If
  End If

End Subroutine SurfaceInfoFromMet

!-------------------------------------------------------------------------------------------------------------

Subroutine SoilInfoFromMet(              &
             Grids,                      &
             NWPMet,                     &
             AncillaryMetSoil,           &
             Time, X, Y, XAncil, YAncil, &
             Soil,                       &
             NWPFlowMemory               &
           )
! Extracts soil information at a particular location from the met fields.

  Implicit None
  ! Argument List:
  Type(Grids_),         Intent(In),    Target :: Grids
  Type(NWPMet_),        Intent(In)            :: NWPMet
  Type(AncillaryMet_),  Intent(In)            :: AncillaryMetSoil
  Type(ShortTime_),     Intent(In)            :: Time
  Real(Std),            Intent(In)            :: X
  Real(Std),            Intent(In)            :: Y
  Real(Std),            Intent(In)            :: XAncil
  Real(Std),            Intent(In)            :: YAncil
  Type(Soil_),          Intent(Out)           :: Soil
  Type(NWPFlowMemory_), Intent(InOut), Target :: NWPFlowMemory
  ! Grids         :: Collection of grids.
  ! M             :: State of an NWP met module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in
  !                  fixed met cases.
  ! X             :} Horizontal location in coord system of horizontal grids used for met data.
  ! Y             :}
  ! XAncil        :] Horizontal location in coord system of horizontal grids used for ancillary data.
  ! YAncil        :]
  ! Soil          :: Soil information.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Locals:
  Type(HGrid_),   Pointer :: HGrid        !} Abbreviations for grids and interpolation
  Type(HGrid_),   Pointer :: HGridAncil   !} coefficients.
  Type(HCoeffs_), Pointer :: HCoeffs      !}
  Type(HCoeffs_)          :: HCoeffsAncil !}
  Type(TCoeffs_), Pointer :: TCoeffs      !}
  Type(TCoeffs_), Pointer :: TCoeffsEnd   !}
  Integer                 :: i            ! Loop index.

  ! $$ Treatment of grids here should be more systematic. For the moment we assume HGrid same
  ! for all ancillary
  ! data sources and use TCoeffsEnd for Ancillary Met, where we assume interpolation is not required.

  ! 1) Set up abbreviations for grids and interpolation coefficients.

  HGrid      => Grids%HGrids(NWPMet%iHGrid)
  HGridAncil => Grids%HGrids(AncillaryMetSoil%iHGrid)   ! $$ Assumes all Ancillaries on same grid
  HCoeffs    => NWPFlowMemory%HCoeffs
  TCoeffs    => NWPFlowMemory%TCoeffs
  TCoeffsEnd => NWPFlowMemory%TCoeffsEnd

  ! 2) Calculate interpolation coefficients.

  ! Note time interpolation irrelevant for monthly data and infinite dt data. Monthly data is regarded as
  ! representative of the month and infinite dt data is valid for all times.
  ! $$ Could possibly do something for monthly and also support yearly and "yearly and monthly" options.
  ! $$ Currently TCoeffs and TCoeffsEnd same as for the NWPMet data - should be different for each ancil
  ! module (with traps for monthly and dt infinite).
  ! $$ see also "$$" comments in AncillaryMet.

  ! Calculate temporal interpolation coefficients.
  If (.not. TCoeffs%Valid) Then
    Call GetTCoeffs(                                           &
           Time,                                               &
           NWPMet%OldTime, Int(ShortTime2RealTime(NWPMet%Dt)), &
           NWPMet%Old, NWPMet%New,                             &
           TCoeffs                                             &
         )
  End If

  ! Calculate time interpolation coefficients for interpolating to a time at the end of the NWP met time
  ! interval (note "end" here means in chronological time not run direction time - hence the use of
  ! IsBackwards). Note there's no need to check TCoeffs%Valid because of the code immediately above. Also
  ! there is no need to set TCoeffsEnd%Valid because this will be done through TCoeffsEnd = TCoeffs.
  If (.not. TCoeffsEnd%Valid) Then
    TCoeffsEnd = TCoeffs
    If (IsBackwards()) Then
      TCoeffsEnd%T1     = 0.0
      TCoeffsEnd%T2     = 1.0
      TCoeffsEnd%iTNear = NWPMet%Old
    Else
      TCoeffsEnd%T1     = 1.0
      TCoeffsEnd%T2     = 0.0
      TCoeffsEnd%iTNear = NWPMet%New
    End If
  End If

  ! Calculate horizontal interpolation coefficients.
  If (.not. HCoeffs%Valid) Then
    Call GetHCoeffs(X, Y, HGrid, HCoeffs)
  End If

    Call GetHCoeffs(XAncil, YAncil, HGridAncil, HCoeffsAncil)

  ! 3) Extract data from arrays.

  ! 2-D data: clay mass fraction and soil particle mass fractions (in particle size ranges).
  ! Note we use the nearest
  ! point rather than interpolating to avoid problems caused by the convention that negative values
  ! indicate sea.
  Soil%ClayFrac        = AncillaryMetSoil%ClayMassFrac(                                   &
                           HCoeffsAncil%iXNear, HCoeffsAncil%iYNear, 1, TCoeffsEnd%iTNear &
                         )
  Do i = 1, 6
    Soil%SoilFracs(i)  = AncillaryMetSoil%SoilMassFracs(                                  &
                           HCoeffsAncil%iXNear, HCoeffsAncil%iYNear, i, TCoeffsEnd%iTNear &
                         )
  End Do


  ! Nonsense - delete $$
  If (.false.) Then
    Soil%ClayFrac        = AncillaryMetSoil%ClayMassFrac(                                       &
                             HCoeffs%iXNear, Max(1, HCoeffs%iYNear + 250), 1, TCoeffsEnd%iTNear &
                           )
    Do i = 1, 6
      Soil%SoilFracs(i)  = AncillaryMetSoil%SoilMassFracs(                                      &
                             HCoeffs%iXNear, Max(1, HCoeffs%iYNear + 250), i, TCoeffsEnd%iTNear &
                           )
    End Do
  End If

End Subroutine SoilInfoFromMet

!-------------------------------------------------------------------------------------------------------------

Subroutine PlantInfoFromMet(             &
             Grids,                      &
             NWPMet,                     &
             NWPMetCanopyHeight,         &
             NWPMetCanopy,               &
             AncillaryMetCanopyHeight,   &
             AncillaryMetLAI,            &
             NWPCanopyHeight,            &
             Time, X, Y, XAncil, YAncil, &
             Plant,                      &
             NWPFlowMemory               &
           )
! Extracts plant information at a particular location from the met fields.

  Implicit None
  ! Argument List:
  Type(Grids_),         Intent(In),    Target :: Grids
  Type(NWPMet_),        Intent(In)            :: NWPMet
  Type(NWPMet_),        Intent(In)            :: NWPMetCanopyHeight
  Type(NWPMet_),        Intent(In)            :: NWPMetCanopy
  Type(AncillaryMet_),  Intent(In)            :: AncillaryMetCanopyHeight ! $$ add comments
  Type(AncillaryMet_),  Intent(In)            :: AncillaryMetLAI
  Logical,              Intent(In)            :: NWPCanopyHeight
  Type(ShortTime_),     Intent(In)            :: Time
  Real(Std),            Intent(In)            :: X
  Real(Std),            Intent(In)            :: Y
  Real(Std),            Intent(In)            :: XAncil
  Real(Std),            Intent(In)            :: YAncil
  Type(Plant_),         Intent(Out)           :: Plant
  Type(NWPFlowMemory_), Intent(InOut), Target :: NWPFlowMemory
  ! Grids         :: Collection of grids.
  ! M             :: State of an NWP met module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in
  !                  fixed met cases.
  ! X             :} Horizontal location in coord system of horizontal grids used for met data.
  ! Y             :}
  ! XAncil        :] Horizontal location in coord system of horizontal grids used for ancillary data.
  ! YAncil        :]
  ! Plant         :: Plant information.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Locals:
  Type(HGrid_),   Pointer :: HGrid        !} Abbreviations for grids and interpolation
  Type(HGrid_),   Pointer :: HGridAncil   !} coefficients.
  Type(HCoeffs_), Pointer :: HCoeffs      !}
  Type(HCoeffs_)          :: HCoeffsAncil !}
  Type(TCoeffs_), Pointer :: TCoeffs      !}
  Type(TCoeffs_), Pointer :: TCoeffsEnd   !}
  Integer                 :: i            ! Loop index.

  ! $$ Treatment of grids here should be more systematic. For the moment we assume HGrid same
  ! for all ancillary
  ! data sources and use TCoeffsEnd for Ancillary Met, where we assume interpolation is not required.

  ! 1) Set up abbreviations for grids and interpolation coefficients.

  HGrid      => Grids%HGrids(NWPMet%iHGrid)
  HGridAncil => Grids%HGrids(AncillaryMetLAI%iHGrid)   ! $$ Assumes all Ancillaries on same grid
  HCoeffs    => NWPFlowMemory%HCoeffs
  TCoeffs    => NWPFlowMemory%TCoeffs
  TCoeffsEnd => NWPFlowMemory%TCoeffsEnd

  ! 2) Calculate interpolation coefficients.

  ! Note time interpolation irrelevant for monthly data and infinite dt data. Monthly data is regarded as
  ! representative of the month and infinite dt data is valid for all times.
  ! $$ Could possibly do something for monthly and also support yearly and "yearly and monthly" options.
  ! $$ Currently TCoeffs and TCoeffsEnd same as for the NWPMet data - should be different for each ancil
  ! module (with traps for monthly and dt infinite).
  ! $$ see also "$$" comments in AncillaryMet.

  ! Calculate temporal interpolation coefficients.
  If (.not. TCoeffs%Valid) Then
    Call GetTCoeffs(                                           &
           Time,                                               &
           NWPMet%OldTime, Int(ShortTime2RealTime(NWPMet%Dt)), &
           NWPMet%Old, NWPMet%New,                             &
           TCoeffs                                             &
         )
  End If

  ! Calculate time interpolation coefficients for interpolating to a time at the end of the NWP met time
  ! interval (note "end" here means in chronological time not run direction time - hence the use of
  ! IsBackwards). Note there's no need to check TCoeffs%Valid because of the code immediately above. Also
  ! there is no need to set TCoeffsEnd%Valid because this will be done through TCoeffsEnd = TCoeffs.
  If (.not. TCoeffsEnd%Valid) Then
    TCoeffsEnd = TCoeffs
    If (IsBackwards()) Then
      TCoeffsEnd%T1     = 0.0
      TCoeffsEnd%T2     = 1.0
      TCoeffsEnd%iTNear = NWPMet%Old
    Else
      TCoeffsEnd%T1     = 1.0
      TCoeffsEnd%T2     = 0.0
      TCoeffsEnd%iTNear = NWPMet%New
    End If
  End If

  ! Calculate horizontal interpolation coefficients.
  If (.not. HCoeffs%Valid) Then
    Call GetHCoeffs(X, Y, HGrid, HCoeffs)
  End If

    Call GetHCoeffs(XAncil, YAncil, HGridAncil, HCoeffsAncil) ! $$ avoid recalculation if HGrid = HGridAncil

  ! 3) Extract data from arrays.

  ! 2-D data: Canopy height (5 plant types), leaf area index (5 plant types), canopy water (5 plant types)
  ! and stomatal conductance (5 plant types). Note we use the nearest point rather than interpolating
  ! to avoid problems caused by the convention that negative values indicate sea.
  Do i = 1, 5
    If (NWPCanopyHeight) Then
      Plant%CanopyHeight(i) = NWPMetCanopyHeight%CanopyHeight(                       &
                                HCoeffs%iXNear, HCoeffs%iYNear, i,   TCoeffs%iTNear  &
                              )
    Else
      Plant%CanopyHeight(i) = AncillaryMetCanopyHeight%CanopyHeight(                           &
                                HCoeffsAncil%iXNear, HCoeffsAncil%iYNear, i, TCoeffsEnd%iTNear &
                              )
    End If
    Plant%LAI(i)            = AncillaryMetLAI%LAI(                                             &
                                HCoeffsAncil%iXNear, HCoeffsAncil%iYNear, i, TCoeffsEnd%iTNear &
                              )
    Plant%CanopyWater(i)    = NWPMetCanopy%CanopyWater(                              &
                                HCoeffs%iXNear, HCoeffs%iYNear, i,   TCoeffs%iTNear  &
                              )
    Plant%StomataConduct(i) = NWPMetCanopy%StomataConduct(                           &
                                HCoeffs%iXNear, HCoeffs%iYNear, i,   TCoeffs%iTNear  &
                              )

    ! fix up stomatal conductance. 1.0E6 values are used for ice. NWP model may not agree with
    ! ancillaries (land use fractions) on location of ice.
    If (Plant%StomataConduct(i) > 1000.0) Plant%StomataConduct(i) = 0.01  ! $$ value of 0.01 arbitrary

  End Do


End Subroutine PlantInfoFromMet

!-------------------------------------------------------------------------------------------------------------

Subroutine EtaToZ(Coords, Grids, M, Time, X, Y, EtaIn, ZOut, NWPFlowMemory)
! Converts from the coord system of the vertical grid used for met data (called Eta
! here, but any coord system is supported) to height above ground.

  Implicit None
  ! Argument List:
  Type(Coords_),        Intent(In),    Target :: Coords
  Type(Grids_),         Intent(In),    Target :: Grids
  Type(NWPMet_),        Intent(In)            :: M
  Type(ShortTime_),     Intent(In)            :: Time
  Real,                 Intent(In)            :: X
  Real,                 Intent(In)            :: Y
  Real,                 Intent(In)            :: EtaIn
  Real,                 Intent(Out)           :: ZOut
  Type(NWPFlowMemory_), Intent(InOut), Target :: NWPFlowMemory
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! M             :: State of an NWP met module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in
  !                  fixed met cases.
  ! X             :} Horizontal location.
  ! Y             :}
  ! EtaIn         :: Height in coord system of the vertical grid used for met data.
  ! ZOut          :: Height (m agl).
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Locals:
  Real(Std)               :: Z1      ! Height at one level in from the edge of the
                                     ! grid.
  Real(Std)               :: T       ! Temperature at edge of grid.
  Type(ZCoord_),  Pointer :: ZCoord  !} Abbreviations for coords, grids and
  Type(HGrid_),   Pointer :: HGrid   !} interpolation coefficients.
  Type(ZGrid_),   Pointer :: ZGrid   !}
  Type(HCoeffs_), Pointer :: HCoeffs !}
  Type(ZCoeffs_), Pointer :: ZCoeffs !}
  Type(TCoeffs_), Pointer :: TCoeffs !}

  ! Set up abbreviations for coords, grids and interpolation coefficients.
  ZCoord  => Coords%ZCoords(M%iZCoord)
  HGrid   => Grids%HGrids(M%iHGrid)
  ZGrid   => Grids%ZGrids(M%iZGrid)
  HCoeffs => NWPFlowMemory%HCoeffs
  ZCoeffs => NWPFlowMemory%ZCoeffs
  TCoeffs => NWPFlowMemory%TCoeffs

  ! Calculate temporal interpolation coefficients.
  If (.not. TCoeffs%Valid) Then
    Call GetTCoeffs(                                 &
           Time,                                     &
           M%OldTime, Int(ShortTime2RealTime(M%Dt)), &
           M%Old, M%New,                             &
           TCoeffs                                   &
         )
  End If

  ! Calculate vertical interpolation coefficients.
  If (.not. ZCoeffs%Valid) Then
    Call GetZCoeffs(EtaIn, ZGrid, ZCoeffs)
  End If

  ! Calculate horizontal interpolation coefficients.
  If (.not. HCoeffs%Valid) Then
    Call GetHCoeffs(X, Y, HGrid, HCoeffs)
  End If

  ! Calculate ZOut (or the height at the grid edge if location is outside grid).
  Call InterpXYZT(HCoeffs, ZCoeffs, TCoeffs, M%Z, ZOut)

  ! Above and below ZGrid we assume constant temperature. Note that, above and below
  ! ZGrid, z-like eta is defined as having d eta / dz constant and equal to the grid
  ! edge value, and p-like eta is defined as proportional to pressure (i.e. defined as
  ! having d log eta / d p constant and equal to 1).
  If (ZCoeffs%ZOutside == 1) Then
    If (ZLike(ZCoord)) Then
      Call InterpXYT(HCoeffs, TCoeffs, M%Z(:, :, ZGrid%nZ - 1, :), Z1)
      ZOut = ZOut +                                       &
             (Z1 - ZOut) * (EtaIn - ZGrid%Z(ZGrid%nZ)) /  &
             (ZGrid%Z(ZGrid%nZ - 1) - ZGrid%Z(ZGrid%nZ))
    Else
      Call InterpXYT(HCoeffs, TCoeffs, M%T(:, :, ZGrid%nZ, :), T)
      ZOut = ZOut + Log(ZGrid%Z(ZGrid%nZ) / EtaIn) * GasConstant * T / Gravity
    End If
  Else If (ZCoeffs%ZOutside == -1) Then
    If (ZLike(ZCoord)) Then
      Call InterpXYT(HCoeffs, TCoeffs, M%Z(:, :, 2, :), Z1)
      ZOut = ZOut +                                                         &
             (Z1 - ZOut) * (EtaIn - ZGrid%Z(1)) / (ZGrid%Z(2) - ZGrid%Z(1))
    Else
      Call InterpXYT(HCoeffs, TCoeffs, M%T(:, :, 1, :), T)
      ZOut = ZOut + Log(ZGrid%Z(1) / EtaIn) * GasConstant * T / Gravity
    End If
  End If

End Subroutine EtaToZ

!-------------------------------------------------------------------------------------------------------------

Subroutine EtaToP(Coords, Grids, M, Time, X, Y, EtaIn, POut, NWPFlowMemory)
! Converts from the coord system of the vertical grid used for met data (called Eta
! here, but any coord system is supported) to pressure.

  Implicit None
  ! Argument List:
  Type(Coords_),        Intent(In),    Target :: Coords
  Type(Grids_),         Intent(In),    Target :: Grids
  Type(NWPMet_),        Intent(In)            :: M
  Type(ShortTime_),     Intent(In)            :: Time
  Real,                 Intent(In)            :: X
  Real,                 Intent(In)            :: Y
  Real,                 Intent(In)            :: EtaIn
  Real,                 Intent(Out)           :: POut
  Type(NWPFlowMemory_), Intent(InOut), Target :: NWPFlowMemory
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! M             :: State of an NWP met module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in
  !                  fixed met cases.
  ! X             :} Horizontal location.
  ! Y             :}
  ! EtaIn         :: Height in coord system of the vertical grid used for met data.
  ! POut          :: Pressure (Pa).
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Locals:
  Real(Std)               :: Z       ! Height at edge of grid.
  Real(Std)               :: Z1      ! Height at one level in from the edge of the
                                     ! grid.
  Real(Std)               :: T       ! Temperature at edge of grid.
  Type(ZCoord_),  Pointer :: ZCoord  !} Abbreviations for coords, grids and
  Type(HGrid_),   Pointer :: HGrid   !} interpolation coefficients.
  Type(ZGrid_),   Pointer :: ZGrid   !}
  Type(HCoeffs_), Pointer :: HCoeffs !}
  Type(ZCoeffs_), Pointer :: ZCoeffs !}
  Type(TCoeffs_), Pointer :: TCoeffs !}

  ! Set up abbreviations for coords, grids and interpolation coefficients.
  ZCoord  => Coords%ZCoords(M%iZCoord)
  HGrid   => Grids%HGrids(M%iHGrid)
  ZGrid   => Grids%ZGrids(M%iZGrid)
  HCoeffs => NWPFlowMemory%HCoeffs
  ZCoeffs => NWPFlowMemory%ZCoeffs
  TCoeffs => NWPFlowMemory%TCoeffs

  ! Calculate temporal interpolation coefficients.
  If (.not. TCoeffs%Valid) Then
    Call GetTCoeffs(                                 &
           Time,                                     &
           M%OldTime, Int(ShortTime2RealTime(M%Dt)), &
           M%Old, M%New,                             &
           TCoeffs                                   &
         )
  End If

  ! Calculate vertical interpolation coefficients.
  If (.not. ZCoeffs%Valid) Then
    Call GetZCoeffs(EtaIn, ZGrid, ZCoeffs)
  End If

  ! Calculate horizontal interpolation coefficients.
  If (.not. HCoeffs%Valid) Then
    Call GetHCoeffs(X, Y, HGrid, HCoeffs)
  End If

  ! Calculate POut (or the pressure at the grid edge if location is outside grid).
  Call InterpXYZT(HCoeffs, ZCoeffs, TCoeffs, M%P, POut)

  ! Above and below ZGrid we assume constant temperature. Note that, above and below
  ! ZGrid, z-like eta is defined as having d eta / dz constant and equal to the grid
  ! edge value, and p-like eta is defined as proportional to pressure (i.e. defined as
  ! having d log eta / d p constant and equal to 1).
  If (ZCoeffs%ZOutside == 1) Then
    If (.not.ZLike(ZCoord)) Then
      POut = POut * EtaIn / ZGrid%Z(ZGrid%nZ)
    Else
      Call InterpXYT(HCoeffs, TCoeffs, M%Z(:, :, ZGrid%nZ    , :), Z )
      Call InterpXYT(HCoeffs, TCoeffs, M%Z(:, :, ZGrid%nZ - 1, :), Z1)
      Call InterpXYT(HCoeffs, TCoeffs, M%T(:, :, ZGrid%nZ    , :), T )
      POut = POut * Exp(                                                         &
                      - Gravity *                                                &
                      (EtaIn - ZGrid%Z(ZGrid%nZ)) *                              &
                      ((Z - Z1) / (ZGrid%Z(ZGrid%nZ) - ZGrid%Z(ZGrid%nZ - 1))) / &
                      (GasConstant * T)                                          &
                    )
    End If
  Else If (ZCoeffs%ZOutside == -1) Then
    If (.not.ZLike(ZCoord)) Then
      POut = POut * EtaIn / ZGrid%Z(1)
    Else
      Call InterpXYT(HCoeffs, TCoeffs, M%Z(:, :, 1, :), Z )
      Call InterpXYT(HCoeffs, TCoeffs, M%Z(:, :, 2, :), Z1)
      Call InterpXYT(HCoeffs, TCoeffs, M%T(:, :, 1, :), T )
      POut = POut * Exp(                                       &
                      - Gravity *                              &
                      (EtaIn - ZGrid%Z(1)) *                   &
                      ((Z - Z1) / (ZGrid%Z(1) - ZGrid%Z(2))) / &
                      (GasConstant * T)                        &
                    )
    End If
  End If

End Subroutine EtaToP

!-------------------------------------------------------------------------------------------------------------

Subroutine ZToEta(Coords, Grids, M, Time, X, Y, ZIn, EtaOut, NWPFlowMemory)
! Converts from height above ground to the coord system of the vertical grid used for
! met data (called Eta here, but any coord system is supported).

  Implicit None
  ! Argument List:
  Type(Coords_),        Intent(In),    Target :: Coords
  Type(Grids_),         Intent(In),    Target :: Grids
  Type(NWPMet_),        Intent(In)            :: M
  Type(ShortTime_),     Intent(In)            :: Time
  Real,                 Intent(In)            :: X
  Real,                 Intent(In)            :: Y
  Real,                 Intent(In)            :: ZIn
  Real,                 Intent(Out)           :: EtaOut
  Type(NWPFlowMemory_), Intent(InOut), Target :: NWPFlowMemory
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! M             :: State of an NWP met module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in
  !                  fixed met cases.
  ! X             :} Horizontal location.
  ! Y             :}
  ! ZIn           :: Height (m agl).
  ! EtaOut        :: Height in coord system of the vertical grid used for met data.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Locals:
  Integer                 :: kMin    !} Indices of levels bracketing location.
  Integer                 :: kMax    !}
  Integer                 :: k       ! Index of a level between kMin and kMax.
  Real(Std)               :: ZMin    !} Height at levels bracketing location.
  Real(Std)               :: ZMax    !}
  Real(Std)               :: Z       ! Height at a level between ZMin and ZMax.
  Real(Std)               :: Z1      ! Height at one level in from the edge of the
                                     ! grid.
  Real(Std)               :: T       ! Temperature at edge of grid.
  Type(ZCoord_),  Pointer :: ZCoord  !} Abbreviations for coords, grids and
  Type(HGrid_),   Pointer :: HGrid   !} interpolation coefficients.
  Type(ZGrid_),   Pointer :: ZGrid   !}
  Type(HCoeffs_), Pointer :: HCoeffs !}
  Type(TCoeffs_), Pointer :: TCoeffs !}

  ! Set up abbreviations for coords, grids and interpolation coefficients.
  ZCoord  => Coords%ZCoords(M%iZCoord)
  HGrid   => Grids%HGrids(M%iHGrid)
  ZGrid   => Grids%ZGrids(M%iZGrid)
  HCoeffs => NWPFlowMemory%HCoeffs
  TCoeffs => NWPFlowMemory%TCoeffs

  ! Calculate temporal interpolation coefficients.
  If (.not. TCoeffs%Valid) Then
    Call GetTCoeffs(                                 &
           Time,                                     &
           M%OldTime, Int(ShortTime2RealTime(M%Dt)), &
           M%Old, M%New,                             &
           TCoeffs                                   &
         )
  End If

  ! Calculate horizontal interpolation coefficients.
  If (.not. HCoeffs%Valid) Then
    Call GetHCoeffs(X, Y, HGrid, HCoeffs)
  End If

  ! Find levels surrounding point.
  kMax = ZGrid%nZ + 1
  kMin = 0
  Do
    k = kMax - kMin
    If (k <= 1) Exit
    k = kMin + k/2
    Call InterpXYT(HCoeffs, TCoeffs, M%Z(:, :, k, :), Z)
    If ((ZIn - Z) * ZGrid%IncreasingHeight > 0.0) Then
      kMin = k
      ZMin = Z
    Else
      kMax = k
      ZMax = Z
    End If
  End Do

  ! Interpolate.
  ! Above and below ZGrid we assume constant temperature. Note that, above and below
  ! ZGrid, z-like eta is defined as having d eta / dz constant and equal to the grid
  ! edge value, and p-like eta is defined as proportional to pressure (i.e. defined as
  ! having d log eta / d p constant and equal to 1).
  If (kMin == ZGrid%nZ) Then
    If (ZLike(ZCoord)) Then
      Call InterpXYT(HCoeffs, TCoeffs, M%Z(:, :, kMin - 1, :), Z1)
      EtaOut = ZGrid%Z(kMin) +                                                  &
               (ZGrid%Z(kMin - 1) - ZGrid%Z(kMin)) * (ZIn - ZMin) / (Z1 - ZMin)
    Else
      Call InterpXYT(HCoeffs, TCoeffs, M%T(:, :, kMin, :), T)
      EtaOut = ZGrid%Z(kMin) * Exp(- Gravity * (ZIn - ZMin) / (GasConstant * T))
    End If
  Else If (kMax == 1) Then
    If (ZLike(ZCoord)) Then
      Call InterpXYT(HCoeffs, TCoeffs, M%Z(:, :, kMax + 1, :), Z1)
      EtaOut = ZGrid%Z(kMax) +                                                  &
               (ZGrid%Z(kMax + 1) - ZGrid%Z(kMax)) * (ZIn - ZMax) / (Z1 - ZMax)
    Else
      Call InterpXYT(HCoeffs, TCoeffs, M%T(:, :, kMax, :), T)
      EtaOut = ZGrid%Z(kMax) * Exp(- Gravity * (ZIn - ZMax) / (GasConstant * T))
    End If
  Else
    EtaOut = ZGrid%Z(kMin) +                                                &
             (ZGrid%Z(kMax) - ZGrid%Z(kMin)) * (ZIn - ZMin) / (ZMax - ZMin)
  End If

End Subroutine ZToEta

!-------------------------------------------------------------------------------------------------------------

Subroutine PToEta(Coords, Grids, M, Time, X, Y, PIn, EtaOut, NWPFlowMemory)
! Converts from pressure to the coord system of the vertical grid used for met data
! (called Eta here, but any coord system is supported).

  Implicit None
  ! Argument List:
  Type(Coords_),        Intent(In),    Target :: Coords
  Type(Grids_),         Intent(In),    Target :: Grids
  Type(NWPMet_),        Intent(In)            :: M
  Type(ShortTime_),     Intent(In)            :: Time
  Real,                 Intent(In)            :: X
  Real,                 Intent(In)            :: Y
  Real,                 Intent(In)            :: PIn
  Real,                 Intent(Out)           :: EtaOut
  Type(NWPFlowMemory_), Intent(InOut), Target :: NWPFlowMemory
  ! Coords        :: Collection of coord systems.
  ! Grids         :: Collection of grids.
  ! M             :: State of an NWP met module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in
  !                  fixed met cases.
  ! X             :} Horizontal location.
  ! Y             :}
  ! PIn           :: Pressure (Pa).
  ! EtaOut        :: Height in coord system of the vertical grid used for met data.
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Locals:
  Integer                 :: kMin    !} Indices of levels bracketing location.
  Integer                 :: kMax    !}
  Integer                 :: k       ! Index of a level between kMin and kMax.
  Real(Std)               :: PMin    !} Pressure at levels bracketing location.
  Real(Std)               :: PMax    !}
  Real(Std)               :: P       ! Pressure at a level between PMin and PMax.
  Real(Std)               :: Z       ! Height at the edge of the grid.
  Real(Std)               :: Z1      ! Height at one level in from the edge of the
                                     ! grid.
  Real(Std)               :: T       ! Temperature at edge of grid.
  Type(ZCoord_),  Pointer :: ZCoord  !} Abbreviations for coords, grids and
  Type(HGrid_),   Pointer :: HGrid   !} interpolation coefficients.
  Type(ZGrid_),   Pointer :: ZGrid   !}
  Type(HCoeffs_), Pointer :: HCoeffs !}
  Type(TCoeffs_), Pointer :: TCoeffs !}

  ! Set up abbreviations for coords, grids and interpolation coefficients.
  ZCoord  => Coords%ZCoords(M%iZCoord)
  HGrid   => Grids%HGrids(M%iHGrid)
  ZGrid   => Grids%ZGrids(M%iZGrid)
  HCoeffs => NWPFlowMemory%HCoeffs
  TCoeffs => NWPFlowMemory%TCoeffs

  ! Calculate temporal interpolation coefficients.
  If (.not. TCoeffs%Valid) Then
    Call GetTCoeffs(                                 &
           Time,                                     &
           M%OldTime, Int(ShortTime2RealTime(M%Dt)), &
           M%Old, M%New,                             &
           TCoeffs                                   &
         )
  End If

  ! Calculate horizontal interpolation coefficients.
  If (.not. HCoeffs%Valid) Then
    Call GetHCoeffs(X, Y, HGrid, HCoeffs)
  End If

  ! Find levels surrounding point.
  kMax = ZGrid%nZ + 1
  kMin = 0
  Do
    k = kMax - kMin
    If (k <= 1) Exit
    k = kMin + k/2
    Call InterpXYT(HCoeffs, TCoeffs, M%P(:, :, k, :), P)
    If ((PIn - P) * ZGrid%IncreasingHeight < 0.0) Then
      kMin = k
      PMin = P
    Else
      kMax = k
      PMax = P
    End If
  End Do

  ! Interpolate.
  ! Above and below ZGrid we assume constant temperature. Note that, above and below
  ! ZGrid, z-like eta is defined as having d eta / dz constant and equal to the grid
  ! edge value, and p-like eta is defined as proportional to pressure (i.e. defined as
  ! having d log eta / d p constant and equal to 1).
  If (kMin == ZGrid%nZ) Then
    If (.not.ZLike(ZCoord)) Then
      EtaOut = ZGrid%Z(kMin) * PIn / PMin
    Else
      Call InterpXYT(HCoeffs, TCoeffs, M%Z(:, :, kMin    , :), Z )
      Call InterpXYT(HCoeffs, TCoeffs, M%Z(:, :, kMin - 1, :), Z1)
      Call InterpXYT(HCoeffs, TCoeffs, M%T(:, :, kMin    , :), T )
      EtaOut = ZGrid%Z(kMin) +                                    &
               ((ZGrid%Z(kMin) - ZGrid%Z(kMin - 1)) / (Z - Z1)) * &
               Log(PMin / PIn) * GasConstant * T / Gravity
    End If
  Else If (kMax == 1) Then
    If (.not.ZLike(ZCoord)) Then
      EtaOut = ZGrid%Z(kMax) * PIn / PMax
    Else
      Call InterpXYT(HCoeffs, TCoeffs, M%Z(:, :, kMax    , :), Z )
      Call InterpXYT(HCoeffs, TCoeffs, M%Z(:, :, kMax + 1, :), Z1)
      Call InterpXYT(HCoeffs, TCoeffs, M%T(:, :, kMax    , :), T )
      EtaOut = ZGrid%Z(kMax) +                                    &
               ((ZGrid%Z(kMax) - ZGrid%Z(kMin + 1)) / (Z - Z1)) * &
               Log(PMax / PIn) * GasConstant * T / Gravity
    End If
  Else
    EtaOut = ZGrid%Z(kMin) +                                                &
             (ZGrid%Z(kMax) - ZGrid%Z(kMin)) * (PIn - PMin) / (PMax - PMin)
  End If

End Subroutine PToEta

!-------------------------------------------------------------------------------------------------------------

Subroutine TopogAndPS(Grids, M, Time, X, Y, NWPFlowMemory, Topog, PS)
! Calculate topographic height and surface pressure.

  Implicit None
  ! Argument List:
  Type(Grids_),         Intent(In), Target :: Grids
  Type(NWPMet_),        Intent(In)         :: M
  Type(ShortTime_),     Intent(In)         :: Time
  Real(Std),            Intent(In)         :: X
  Real(Std),            Intent(In)         :: Y
  Type(NWPFlowMemory_), Intent(In), Target :: NWPFlowMemory
  Real(Std),            Intent(Out)        :: Topog
  Real(Std),            Intent(Out)        :: PS
  ! Grids         :: Collection of grids.
  ! M             :: State of an NWP met module instance.
  ! Time          :: Current time in non-fixed met cases and the time of the met in
  !                  fixed met cases.
  ! X             :} Horizontal location.
  ! Y             :}
  ! NWPFlowMemory :: Flow memory of NWP flow module instance.
  ! Topog         :: Height of topography above sea level (m).
  ! PS            :: Surface pressure (Pa).
  ! Locals:
  Type(HGrid_),   Pointer :: HGrid   !} Abbreviations for grids and interpolation
  Type(HCoeffs_), Pointer :: HCoeffs !} coefficients.
  Type(TCoeffs_), Pointer :: TCoeffs !}

  ! Set up abbreviations for grids and interpolation coefficients.
  HGrid   => Grids%HGrids(M%iHGrid)
  HCoeffs => NWPFlowMemory%HCoeffs
  TCoeffs => NWPFlowMemory%TCoeffs

  ! Calculate temporal interpolation coefficients.
  If (.not. TCoeffs%Valid) Then
    Call GetTCoeffs(                                 &
           Time,                                     &
           M%OldTime, Int(ShortTime2RealTime(M%Dt)), &
           M%Old, M%New,                             &
           TCoeffs                                   &
         )
  End If

  ! Calculate horizontal interpolation coefficients.
  If (.not. HCoeffs%Valid) Then
    Call GetHCoeffs(X, Y, HGrid, HCoeffs)
  End If

  ! Calculate topographic height and surface pressure.
  Call InterpXYT(HCoeffs, TCoeffs, M%Topog,         Topog)
  Call InterpXYT(HCoeffs, TCoeffs, M%P(:, :, 1, :), PS   )

End Subroutine TopogAndPS

!-------------------------------------------------------------------------------------------------------------

End Module NWPFlowModule

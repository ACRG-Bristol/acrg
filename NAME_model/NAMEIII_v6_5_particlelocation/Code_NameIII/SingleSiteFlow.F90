! Module: Single Site Flow Module

Module SingleSiteFlowModule

! This is a flow module designed to use single site met data, construct mean flow
! profiles and add turbulence information.

! It uses input from a single instance of the single site met module.
!
! It does not use any input from other flow modules (and hence uses no update subset).
!
! It assumes no topography and takes the topographic height to be zero and the surface
! pressure to be the sea level pressure in the ICAO standard atmosphere.
!
! It has attributes Update, Convert, Flow, Cloud and Rain.
!
! The coords stored in SingleSiteFlow%C are
! Horizontal: 1: Cartesian coord system used in single site met module to define wind
!                directions.
! Vertical:   1: m agl.
! GetSingleSiteFlow returns flow information in horizontal coord system 1 and in
! vertical coord system 1.
!
! No grids are stored in SingleSiteFlow%C.
!
! The domain must exclude points not definable in horizontal coord system 1.
!
! It works in a calendar or relative time frame.
!
! The status of the optional flow module features is as follows (where ???? here =
! SingleSite):
!     ????FlowMemory_ data type           - not used
!     SetUpCoordsEtc_????Flow routine     - not used
!     SetUp????Flow_MetsCoordsEtc routine -     used
!     Reset????FlowMemory routine         - not used
!     ????ConvertCoordIndices routine     -     used
!     ????ConvertToZ routine              -     used (without flow memory argument)
!     ????FlowCoordIndices routine        -     used
!     Get????Flow routine                 -     used (without flow memory argument)
!     ????CloudCoordIndices routine       -     used
!     Get????Cloud routine                -     used (without flow memory argument)
!     ????RainCoordIndices routine        -     used
!     Get????Rain routine                 -     used (without flow memory argument)
!     ????ReflectCoordIndices             -     used
!     ????Reflect routine                 -     used.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowAndFlowProfileModule
Use SingleSiteMetModule
Use MetsModule
Use CommonFlowModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: SingleSiteFlow_                   ! The state of a single site flow module
                                             ! instance.
Public  :: InitSingleSiteFlow                ! Initialises the state of a single site
                                             ! flow module instance.
Public  :: SetUpSingleSiteFlow_MetsCoordsEtc ! Sets up SingleSiteFlow using
                                             ! information from EtaDefns, Coords,
                                             ! Grids and Mets.
Public  :: PrepareForUpdateSingleSiteFlow    ! Prepares for updating an instance of the single site flow
                                             ! module.
Public  :: SingleSiteFlowReqs                ! Specifies what information the flow
                                             ! module instance wants from the other
                                             ! flow module instances.
Public  :: UpdateSingleSiteFlow              ! Updates an instance of the single site
                                             ! flow module.
Public  :: SingleSiteConvertCoordIndices     ! Returns indices of coord systems in
                                             ! which positions need to be specified
                                             ! when using SingleSiteConvertToZ.
Public  :: SingleSiteConvertToZ              ! Converts vertical coords between coord
                                             ! systems using a single site flow module
                                             ! instance.
Public  :: SingleSiteFlowCoordIndices        ! Returns indices of coord systems in
                                             ! which positions need to be specified
                                             ! when using GetSingleSiteFlow.
Public  :: GetSingleSiteFlow                 ! Gets flow information from a single
                                             ! site flow module instance.
Public  :: SingleSiteCloudCoordIndices       ! Returns indices of coord systems in
                                             ! which positions need to be specified
                                             ! when using GetSingleSiteCloud.
Public  :: GetSingleSiteCloud                ! Gets cloud information from a single
                                             ! site flow module instance.
Public  :: SingleSiteRainCoordIndices        ! Returns indices of coord systems in
                                             ! which positions need to be specified
                                             ! when using GetSingleSiteRain.
Public  :: GetSingleSiteRain                 ! Gets rain information from a single
                                             ! site flow module instance.
Public  :: SingleSiteReflectCoordIndices     ! Returns indices of coord systems in
                                             ! which positions need to be specified
                                             ! when using SingleSiteReflect.
Public  :: SingleSiteReflect                 ! Reflects position of particle or puff
                                             ! centroid using a single site flow
                                             ! module instance.

!-------------------------------------------------------------------------------------------------------------

Type :: SingleSiteFlow_ ! The state of a single site flow module instance.
  Type (CommonFlow_)             :: C
  Type (SingleSiteMet_), Pointer :: SingleSiteMet
  ! C             :: Part of flow state common to all flow modules.
  ! SingleSiteMet :: Pointer to stored met data.
End Type SingleSiteFlow_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitSingleSiteFlow(    &
           FlowName,            &
           MetModName, MetName, &
           DomainName,          &
           FixedMet,            &
           UpdateOnDemand       &
         )                      &
Result(SingleSiteFlow)
! Initialises the state of a single site flow module instance.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: FlowName       ! Name of flow module instance.
  Character(*), Intent(In) :: MetModName     ! Name of used met module.
  Character(*), Intent(In) :: MetName        ! Name of used met module instance.
  Character(*), Intent(In) :: DomainName     ! Name of the domain of the flow module instance.
  Logical,      Intent(In) :: FixedMet       ! Indicates that the met is fixed.
  Logical,      Intent(In) :: UpdateOnDemand ! Indicates the flow module instance is to be updated using
                                             ! update-on-demand.
  ! Function result:
  Type(SingleSiteFlow_) :: SingleSiteFlow ! Initialised state of a single site flow module instance.

  ! Met module instance used.
  If (.not.(MetModName .CIEq. 'Single Site Met')) Then
    Call Message(                                                         &
           'Error in InitSingleSiteFlow: met module name is given as ' // &
           Trim(MetModName)                                            // &
           ' instead of Single Site Met',                                 &
           3                                                              &
         )
  End If

  SingleSiteFlow%C = InitCommonFlow(                        &
                       FlowModName    = 'Single Site Flow', &
                       FlowName       = FlowName,           &
                       nMets          = 1,                  &
                       MetModNames    = (/ MetModName /),   &
                       MetNames       = (/ MetName /),      &
                       DomainName     = DomainName,         &
                       nAttribs       = 5,                  &
                       AttribNames    = (/                  &
                                          'Update ',        &
                                          'Convert',        &
                                          'Flow   ',        &
                                          'Cloud  ',        &
                                          'Rain   '         &
                                        /),                 &
                       FixedMet       = FixedMet,           &
                       UpdateOnDemand = UpdateOnDemand      &
                     )

End Function InitSingleSiteFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpSingleSiteFlow_MetsCoordsEtc( &
             EtaDefns, Coords, Grids, Mets,   &
             SingleSiteFlow                   &
           )
! Sets up SingleSiteFlow using information from EtaDefns, Coords, Grids and Mets.

  Implicit None
  ! Argument list:
  Type(EtaDefns_),       Intent(In)           :: EtaDefns       ! Collection of eta
                                                                ! definitions.
  Type(Coords_),         Intent(In)           :: Coords         ! Collection of coord
                                                                ! systems.
  Type(Grids_),          Intent(In)           :: Grids          ! Collection of grids.
  Type(Mets_),           Intent(In),   Target :: Mets           ! Set of met module
                                                                ! instance states.
  Type(SingleSiteFlow_), Intent(InOut)        :: SingleSiteFlow ! State of a single
                                                                ! site flow module
                                                                ! instance.
  ! Locals:
  Integer :: iMetMod ! Met module index.
  Integer :: iMet    ! Met module instance index.

  ! Find index of met module instance.
  Call FindMetIndex(                                                    &
         SingleSiteFlow%C%MetModNames(1), SingleSiteFlow%C%MetNames(1), &
         Mets,                                                          &
         iMetMod, iMet                                                  &
       )

  ! Set pointer to met information.
  SingleSiteFlow%SingleSiteMet => Mets%SingleSiteMets(iMet)

  ! Add coord systems from met module instance to SingleSiteFlow%C.
  Call AddCoordsToCommonFlow(                             &
         1, (/ Mets%C(iMetMod, iMet)%P%HCoordNames(1) /), &
         1, (/ Mets%C(iMetMod, iMet)%P%ZCoordNames(1) /), &
         SingleSiteFlow%C                                 &
       )

  ! Check number of coord systems in SingleSiteFlow.
  If (SingleSiteFlow%C%nHCoords /= 1) Then
    Call Message(                                                       &
           'Unexpected error in SetUpSingleSiteFlow_MetsCoordsEtc: ' // &
           'error in numbering of horizontal coord systems '         // &
           'used by the single site flow module instance "'          // &
           Trim(SingleSiteFlow%C%FlowName)                           // &
           '"',                                                         &
           4                                                            &
         )
  End If
  If (SingleSiteFlow%C%nZCoords /= 1) Then
    Call Message(                                                       &
           'Unexpected error in SetUpSingleSiteFlow_MetsCoordsEtc: ' // &
           'error in numbering of vertical coord systems '           // &
           'used by the single site flow module instance "'          // &
           Trim(SingleSiteFlow%C%FlowName)                           // &
           '"',                                                         &
           4                                                            &
         )
  End If

End Subroutine SetUpSingleSiteFlow_MetsCoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForUpdateSingleSiteFlow( &
             Coords, Grids, Domains,       &
             Mets,                         &
             iCase,                        &
             Time,                         &
             TValid, UpdateNow,            &
             SingleSiteFlow,               &
             Units                         &
           )
! Prepares for updating an instance of the single site flow module.

! This routine must set TValid and UpdateNow but must not alter the validity of the single site flow module
! instance.
! $$ this comment should be with generic definition.

  Implicit None
  ! Argument list:
  Type(Coords_),         Intent(In)           :: Coords
  Type(Grids_),          Intent(In)           :: Grids
  Type(Domains_),        Intent(In),   Target :: Domains
  Type(Mets_),           Intent(In)           :: Mets
  Integer,               Intent(In)           :: iCase
  Type(Time_),           Intent(In)           :: Time
  Type(Time_),           Intent(Out)          :: TValid
  Logical,               Intent(Out)          :: UpdateNow
  Type(SingleSiteFlow_), Intent(InOut)        :: SingleSiteFlow
  Type(Units_),          Intent(InOut)        :: Units
  ! Coords         :: Collection of coord systems.
  ! Grids          :: Collection of grids.
  ! Domains        :: Collection of domains.
  ! Mets           :: Set of met module instance states.
  ! iCase          :: Number of case.
  ! Time           :: Time for which the single site flow module instance might be updated.
  ! TValid         :: Earliest time that the validity (overall or for any single attribute) of the flow module
  !                   instance might change, except perhaps for a change caused by a change in the validity of
  !                   the met and flow module instances acting as data sources, assuming the flow module
  !                   instance is updated now. The value is that determined at the end of this routine (the
  !                   actual time may be later).
  ! UpdateNow      :: Indicates the flow module instance must be updated now (even if update-on-demand is
  !                   specified). If set, TValid need not be set to any particular time.
  ! SingleSiteFlow :: State of a single site flow module instance.
  ! Units          :: Collection of information on input/output unit numbers.
  ! Locals:
  Type(Domain_), Pointer :: Domain ! Abbreviation for domain.

  ! Abbreviations.
  Domain => Domains%Domains(SingleSiteFlow%C%iDomain)

  ! Validity due to time domain.
  If (Time < StartTimeOfDomain(Domain)) Then
    TValid = StartTimeOfDomain(Domain)
  Else If (Time >= EndTimeOfDomain(Domain)) Then
    TValid = InfFutureTime()
  Else
    TValid = EndTimeOfDomain(Domain)
  End If
  UpdateNow = .false.

End Subroutine PrepareForUpdateSingleSiteFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine SingleSiteFlowReqs(                                         &
             SingleSiteFlow, Time,                                     &
             WantFlowField,     WantCloudField,     WantRainField,     &
             iFlowUpdateSubset, iCloudUpdateSubset, iRainUpdateSubset, &
             FlowField,         CloudField,         RainField          &
           )
! Specifies what information the flow module instance wants from the other flow module
! instances.

  Implicit None
  ! Argument list:
  Type(SingleSiteFlow_), Intent(In)  :: SingleSiteFlow
  Type(Time_),           Intent(In)  :: Time
  Logical,               Intent(Out) :: WantFlowField
  Logical,               Intent(Out) :: WantCloudField
  Logical,               Intent(Out) :: WantRainField
  Integer,               Intent(Out) :: iFlowUpdateSubset
  Integer,               Intent(Out) :: iCloudUpdateSubset
  Integer,               Intent(Out) :: iRainUpdateSubset
  Type(FlowField_),      Intent(Out) :: FlowField
  Type(CloudField_),     Intent(Out) :: CloudField
  Type(RainField_),      Intent(Out) :: RainField
  ! SingleSiteFlow     :: State of a single site flow module instance.
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

End Subroutine SingleSiteFlowReqs

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateSingleSiteFlow(               &
             Coords, Grids, Domains,           &
             Mets,                             &
             iCase,                            &
             Time,                             &
             FlowField, CloudField, RainField, &
             SingleSiteFlow,                   &
             Units                             &
           )
! Updates an instance of the single site flow module.

  Implicit None
  ! Argument list:
  Type(Coords_),         Intent(In)           :: Coords
  Type(Grids_),          Intent(In)           :: Grids
  Type(Domains_),        Intent(In),   Target :: Domains
  Type(Mets_),           Intent(In)           :: Mets
  Integer,               Intent(In)           :: iCase
  Type(Time_),           Intent(In)           :: Time
  Type(FlowField_),      Intent(In)           :: FlowField
  Type(CloudField_),     Intent(In)           :: CloudField
  Type(RainField_),      Intent(In)           :: RainField
  Type(SingleSiteFlow_), Intent(InOut)        :: SingleSiteFlow
  Type(Units_),          Intent(InOut)        :: Units
  ! Coords         :: Collection of coord systems.
  ! Grids          :: Collection of grids.
  ! Domains        :: Collection of domains.
  ! Mets           :: Set of met module instance states.
  ! iCase          :: Number of case.
  ! Time           :: Time for which the single site flow module instance is to be
  !                   updated.
  ! FlowField      :} Field of information for various attributes needed by the flow
  ! CloudField     :} module instance.
  ! RainField      :}
  ! SingleSiteFlow :: State of a single site flow module instance.
  ! Units          :: Collection of information on input/output unit numbers.
  ! Locals:
  Type(Domain_), Pointer :: Domain  ! Abbreviation for domain.
  Logical                :: Valid   !} Validity due to time domain.
  Type(Time_)            :: TValid  !}
  Logical                :: MValid  !] Validity due to met.
  Type(Time_)            :: MTValid !]

  ! Abbreviations.
  Domain => Domains%Domains(SingleSiteFlow%C%iDomain)

  ! Validity due to time domain.
  Valid = StartTimeOfDomain(Domain) <= Time .and. Time < EndTimeOfDomain(Domain)
  If (Valid) Then
    TValid = EndTimeOfDomain(Domain)
  Else
    TValid = StartTimeOfDomain(Domain)
    If (Time >= TValid) TValid = InfFutureTime()
  End If

  ! Validity due to met.
  Call MetValid(                                                  &
         Mets,                                                    &
         SingleSiteFlow%C%iMetMods(1), SingleSiteFlow%C%iMets(1), &
         Time,                                                    &
         MValid, MTValid                                          &
       )

  ! Overall validity.
  SingleSiteFlow%C%Valid                        = Valid .and. MValid
  SingleSiteFlow%C%ValidAttribParams(:)         = .false.
  SingleSiteFlow%C%ValidAttribParams(A_Convert) = SingleSiteFlow%C%Valid
  SingleSiteFlow%C%ValidAttribParams(A_Flow   ) = SingleSiteFlow%C%Valid
  SingleSiteFlow%C%ValidAttribParams(A_Cloud  ) = SingleSiteFlow%C%Valid
  SingleSiteFlow%C%ValidAttribParams(A_Rain   ) = SingleSiteFlow%C%Valid
  SingleSiteFlow%C%ValidityUnimprovable         = SingleSiteFlow%C%Valid .or. SingleSiteFlow%C%FixedMet
  SingleSiteFlow%C%TValid                       = TValid
  If (SingleSiteFlow%C%Valid) SingleSiteFlow%C%TValid = TMin(SingleSiteFlow%C%TValid, MTValid)

End Subroutine UpdateSingleSiteFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine SingleSiteConvertCoordIndices(SingleSiteFlow, nHCoords, iHCoords)
! Returns indices of coord systems in which positions need to be specified when using
! SingleSiteConvertToZ.

  Implicit None
  ! Argument list:
  Type(SingleSiteFlow_), Intent(In)  :: SingleSiteFlow       ! State of a single site
                                                             ! flow module instance.
  Integer,               Intent(Out) :: nHCoords             !} Number and values of
  Integer,               Intent(Out) :: iHCoords(MaxHCoords) !} indices of horizontal
                                                             !} coord systems.

  nHCoords    = 0
  iHCoords(1) = 0

End Subroutine SingleSiteConvertCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine SingleSiteConvertToZ(            &
             SingleSiteFlow, Coords, Grids, &
             iZCoord, Time, Position        &
           )
! Converts vertical coords between coord systems using a single site flow module
! instance.

  Implicit None
  ! Argument list:
  Type(SingleSiteFlow_), Intent(In)           :: SingleSiteFlow
  Type(Coords_),         Intent(In),   Target :: Coords
  Type(Grids_),          Intent(In)           :: Grids
  Integer,               Intent(In)           :: iZCoord
  Type(ShortTime_),      Intent(In)           :: Time
  Type(Position_),       Intent(InOut)        :: Position
  ! SingleSiteFlow :: State of a single site flow module instance.
  ! Coords         :: Collection of coord systems.
  ! Grids          :: Collection of grids.
  ! iZCoord        :: Index in Coords of coord system to convert to.
  ! Time           :: Current time in non-fixed met cases and the time of the met in
  !                   fixed met cases.
  ! Position       :: Coords of the point in various coord systems in Coords, with
  !                   flags to indicate whether the values are valid.
  ! Locals:
  Type(ZCoord_), Pointer :: ZCoordIn  !} Abreviations for the coord system used to
  Type(ZCoord_), Pointer :: ZCoordOut !} convert from, the coord system to convert to,
  Type(ZCoord_)          :: ZCoord1   !} and two intermediate coord systems.
  Type(ZCoord_)          :: ZCoord2   !}
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

  Topog = 0.0
  PS    = PAt0km

  Select Case (ZCoordIn%CoordType)

    Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)

      Select Case (ZCoordOut%CoordType)

        Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)

          ZOut = ZBasedToZBased(ZCoordIn, ZCoordOut, Topog, ZIn)

        Case (Z_P, Z_PAsZ, Z_PAsEta)

          ZCoord1 = ZCoord_m_asl()
          ZCoord2 = InitZCoord('PAsZ', Z_PAsZ, 1.0)
          Z1      = ZBasedToZBased    (ZCoordIn, ZCoord1,   Topog, ZIn)
          Z2      = AboveSeaToPAsZICAO(ZCoord1,  ZCoord2,          Z1 )
          ZOut    = PBasedToPBased    (ZCoord2,  ZCoordOut, PS,    Z2 )

        Case Default

          Call Message('Unexpected error in SingleSiteConvertToZ: unknown coordinate type', 4)

      End Select

    Case (Z_P, Z_PAsZ, Z_PAsEta)

      Select Case (ZCoordOut%CoordType)

        Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta)

          ZCoord1 = InitZCoord('Z_PAsZ', Z_PAsZ, 1.0)
          ZCoord2 = ZCoord_m_asl()
          Z1      = PBasedToPBased    (ZCoordIn, ZCoord1,   PS,    ZIn)
          Z2      = PAsZToAboveSeaICAO(ZCoord1,  ZCoord2,          Z1 )
          ZOut    = ZBasedToZBased    (ZCoord2,  ZCoordOut, Topog, Z2 )

        Case (Z_P, Z_PAsZ, Z_PAsEta)

          ZOut = PBasedToPBased(ZCoordIn, ZCoordOut, PS, ZIn)

        Case Default

          Call Message('Unexpected error in SingleSiteConvertToZ: unknown coordinate type', 4)

      End Select

    Case Default

      Call Message('Unexpected error in SingleSiteConvertToZ: unknown coordinate type', 4)

  End Select

  Position%Z(iZCoord) = ZOut

End Subroutine SingleSiteConvertToZ

!-------------------------------------------------------------------------------------------------------------

Subroutine SingleSiteFlowCoordIndices(              &
             SingleSiteFlow,                        &
             nHCoords, iHCoords, nZCoords, iZCoords &
           )
! Returns indices of coord systems in which positions need to be specified when using
! GetSingleSiteFlow.

  Implicit None
  ! Argument list:
  Type(SingleSiteFlow_), Intent(In)  :: SingleSiteFlow       ! State of a single site
                                                             ! flow module instance.
  Integer,               Intent(Out) :: nHCoords             !} Number and values of
  Integer,               Intent(Out) :: iHCoords(MaxHCoords) !} indices of horizontal
  Integer,               Intent(Out) :: nZCoords             !} and vertical coord
  Integer,               Intent(Out) :: iZCoords(MaxHCoords) !} systems.

  ! Horizontal coords needed:
  ! 1) Cartesian coord system used in single site met module to define wind
  !    directions.
  nHCoords    = 1
  iHCoords(1) = SingleSiteFlow%C%iHCoords(1)

  ! Vertical coords needed:
  ! 1) m agl.
  nZCoords    = 1
  iZCoords(1) = SingleSiteFlow%C%iZCoords(1)

End Subroutine SingleSiteFlowCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine GetSingleSiteFlow(          &
             Coords, Grids,            &
             SingleSiteFlow,           &
             Moisture, Inhomog, Homog, &
             Time, Position,           &
             Flow, ProfileData         &
           )
! Gets flow information from a single site flow module instance.

  Implicit None
  ! Argument list:
  Type(Coords_),         Intent(In)            :: Coords
  Type(Grids_),          Intent(In)            :: Grids
  Type(SingleSiteFlow_), Intent(In),  Target   :: SingleSiteFlow
  Logical,               Intent(In)            :: Moisture
  Logical,               Intent(In)            :: Inhomog
  Logical,               Intent(In)            :: Homog
  Type(ShortTime_),      Intent(In)            :: Time
  Type(Position_),       Intent(In)            :: Position
  Type(Flow_),           Intent(Out)           :: Flow
  Type(ProfileData_),    Intent(Out), Optional :: ProfileData
  ! Coords         :: Collection of coord systems.
  ! Grids          :: Collection of grids.
  ! SingleSiteFlow :: State of a single site flow module instance.
  ! Moisture       :: Indicates Q is required.
  ! Inhomog        :: Indicates inhomogeneous quantities are required.
  ! Homog          :: Indicates homogeneous quantities are required.
  ! Time           :: Current time in non-fixed met cases and the time of the met in fixed met cases.
  ! Position       :: Coords of the point in various coord systems in Coords, with flags to indicate whether
  !                   the values are valid.
  ! Flow           :: Flow information.
  ! ProfileData    :: Data needed to construct idealised analytic mean flow and turbulence profiles.
  ! Locals:
  Type(ProfileData_), Pointer :: Profile1    !} Abbreviations for profile data.
  Type(ProfileData_), Pointer :: Profile2    !}
  Type(ProfileData_)          :: Profile     ! Interpolated profile data.
  Type(TCoeffs_)              :: TCoeffs     ! Interpolation coefficients.
  Type(Flow_)                 :: Flow2       ! Flow information corresponding to ProfileData2.
  Real(Std)                   :: Z           ! Height (m agl).
  Real(Std)                   :: Speed       ! Wind speed.
  Real(Std)                   :: U           !} Variables used to interpolate stress and geostrophic wind as
  Real(Std)                   :: V           !} vectors.
  Real(Std)                   :: U1          !}
  Real(Std)                   :: V1          !}
  Real(Std)                   :: U2          !}
  Real(Std)                   :: V2          !}

  ! Coords of location.
  Z = Position%Z(SingleSiteFlow%C%iZCoords(1))

  ! Interpolation coefficients.
  Call GetTCoeffs(                                                     &
         Time    = Time,                                               &
         OldTime = Time2ShortTime(SingleSiteFlow%SingleSiteMet%Time1), &
         DtSec   = NInt(SingleSiteFlow%SingleSiteMet%ActualDt),        & ! $$ real?
         iOld    = 1,                                                  &
         iNew    = 2,                                                  &
         TCoeffs = TCoeffs                                             &
       )

  Profile1 => SingleSiteFlow%SingleSiteMet%ProfileData1
  Profile2 => SingleSiteFlow%SingleSiteMet%ProfileData2

  ! dUdT.
  Call MeanFlowProfiles(            & ! $$ speed up with thermodynamic switch
         Z              = Z,        &
         ProfileData    = Profile1, &
         Moisture       = .false.,  &
         HMinus         = .false.,  &
         Speed          = Speed,    &
         Flow           = Flow      &
       )
  Call MeanFlowProfiles(            &
         Z              = Z,        &
         ProfileData    = Profile2, &
         Moisture       = .false.,  &
         HMinus         = .false.,  &
         Speed          = Speed,    &
         Flow           = Flow2     &
       )
  Flow2%dUdT(:) = (Flow2%U(:) - Flow%U(:)) / SingleSiteFlow%SingleSiteMet%ActualDt ! $$ Fixed met etc

  ! Interpolate ProfileData.
  ! $$ Only some quantities interpolated. Other quantities? (these are currently fixed in SingleSiteMet.F90
  ! but may not always be.
  Profile = Profile1

  U1 = Profile1%UStar**2 * Cos(Profile1%Phi0)
  V1 = Profile1%UStar**2 * Sin(Profile1%Phi0)
  U2 = Profile2%UStar**2 * Cos(Profile2%Phi0)
  V2 = Profile2%UStar**2 * Sin(Profile2%Phi0)
  U = U1 * TCoeffs%T2 + U2 * TCoeffs%T1
  V = V1 * TCoeffs%T2 + V2 * TCoeffs%T1
  Profile%UStar = Sqrt(Sqrt(U**2 + V**2))
  Profile%Phi0  = ATan2ZeroTest(V, U)

  U1 = Profile1%UG * Cos(Profile1%PhiG)
  V1 = Profile1%UG * Sin(Profile1%PhiG)
  U2 = Profile2%UG * Cos(Profile2%PhiG)
  V2 = Profile2%UG * Sin(Profile2%PhiG)
  U = U1 * TCoeffs%T2 + U2 * TCoeffs%T1
  V = V1 * TCoeffs%T2 + V2 * TCoeffs%T1
  Profile%UG   = Sqrt(U**2 + V**2)
  Profile%PhiG = ATan2ZeroTest(V, U)

  Profile%DeltaPhi = Profile%PhiG - Profile%Phi0
  If (Profile%DeltaPhi >  Pi) Profile%DeltaPhi = Profile%DeltaPhi - 2.0*Pi
  If (Profile%DeltaPhi < -Pi) Profile%DeltaPhi = Profile%DeltaPhi + 2.0*Pi

  Profile%H          = Profile1%H          * TCoeffs%T2 + Profile2%H          * TCoeffs%T1
  Profile%T0         = Profile1%T0         * TCoeffs%T2 + Profile2%T0         * TCoeffs%T1
  Profile%WT         = Profile1%WT         * TCoeffs%T2 + Profile2%WT         * TCoeffs%T1
  Profile%THPlus     = Profile1%THPlus     * TCoeffs%T2 + Profile2%THPlus     * TCoeffs%T1
  Profile%dTdZHPlus  = Profile1%dTdZHPlus  * TCoeffs%T2 + Profile2%dTdZHPlus  * TCoeffs%T1
  Profile%Q0         = Profile1%Q0         * TCoeffs%T2 + Profile2%Q0         * TCoeffs%T1
  Profile%WQ         = Profile1%WQ         * TCoeffs%T2 + Profile2%WQ         * TCoeffs%T1
  Profile%RHHPlus    = Profile1%RHHPlus    * TCoeffs%T2 + Profile2%RHHPlus    * TCoeffs%T1
  Profile%dRHdZHPlus = Profile1%dRHdZHPlus * TCoeffs%T2 + Profile2%dRHdZHPlus * TCoeffs%T1
  Profile%PHPlus     = Profile1%PHPlus     * TCoeffs%T2 + Profile2%PHPlus     * TCoeffs%T1

  If (Profile%UStar /= 0.0) Then ! $ Better not to bother with these in Profile and calc as needed?
                                 ! Could also be usefully removed from Flow or substituted for.
                                 ! Could add functions for reciplmo etc to FlowAndFlowProfileModule
    Profile%RecipLMO = - VK * Profile%WT * Gravity / (Profile%T0 * Profile%UStar**3)
  Else If (Profile%WT > 0.0) Then
    Profile%RecipLMO = - 200000.0
  Else
    Profile%RecipLMO = 200000.0
  End If
  If (Profile%WT > 0.0) Then
    Profile%WStar = (Profile%H * Profile%WT * Gravity / Profile%T0) ** (1.0/3.0)
  Else
    Profile%WStar = 0.0
  End If
  If (Profile%UStar /= 0.0) Then
    Profile%TStar = - Profile%WT / Profile%UStar
  Else If (Profile%WT > 0.0) Then
    Profile%TStar = - 200000.0
  Else
    Profile%TStar = 200000.0
  End If
  If (Profile%UStar /= 0.0) Then
    Profile%QStar = - Profile%WQ / Profile%UStar
  Else If (Profile%WQ > 0.0) Then
    Profile%QStar = - 200000.0
  Else
    Profile%QStar = 200000.0
  End If

  ! 'Mean quantities', 'boundary layer characteristics' and 'coord systems' parts of Flow.
  Call MeanFlowProfiles(         &
         Z           = Z,        &
         ProfileData = Profile,  &
         Moisture    = Moisture, &
         HMinus      = .false.,  &
         Speed       = Speed,    &
         Flow        = Flow      &
       )

  ! dUdT.
  Flow%dUdT(:) = Flow2%dUdT(:)

  ! 'Inhomogeneous turbulence quantities', 'inhomogeneous eddy diffusivities', 'homogeneous turbulence
  ! quantities', 'unresolved mesoscale motion quantities', 'boundary layer characteristics' and 'coord systems' 
  ! parts of Flow.
  ! note taumin = 20.0 - should be input ? $$
  Call TurbProfiles(             &
         Z           = Z,        &
         ProfileData = Profile,  &
         Inhomog     = Inhomog,  &
         Homog       = Homog,    &
         UEqV        = .true.,   & ! $$
         TauMin      = 20.0,     &
         ZMin        = 0.0,      & ! $$
         Flow        = Flow      &
       )

  ! Topography.
  Flow%ZS          = 0.0
  Flow%Topog       = 0.0
  Flow%dTopogdX    = 0.0
  Flow%SmoothTopog = .true.

  ! Time step limits.
  Flow%MaxDZ = Huge(Flow%MaxDZ)
  If (Inhomog) Then
    Flow%MaxDZ = Min(Flow%MaxDZ, Flow%H)
    If (Flow%dSigUUdX(3) /= 0.0) Then
      Flow%MaxDZ = Min(Flow%MaxDZ, Abs(Flow%SigUU(3) / Flow%dSigUUdX(3)))
    End If
    If (Flow%dTauUUdZ(3) /= 0.0) Then
      Flow%MaxDZ = Min(Flow%MaxDZ, Abs(Flow%TauUU(3) / Flow%dTauUUdZ(3)))
    End If
    If (Flow%dKdX(3) /= 0.0) Then
      Flow%MaxDZ = Min(Flow%MaxDZ, Flow%K(3) / Abs(Flow%dKdX(3)))
    End If
  End If

  ! Puff size limits.
  Flow%DeltaI = Flow%H/4.0 ! $$
 ! Flow%DeltaI = Flow%H/40.0 ! $$

  ! ProfileData.
  If (Present(ProfileData)) Then
    ProfileData = Profile
    ProfileData%Canopy = .false.
  End If

End Subroutine GetSingleSiteFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine SingleSiteCloudCoordIndices(             &
             SingleSiteFlow,                        &
             nHCoords, iHCoords, nZCoords, iZCoords &
           )
! Returns indices of coord systems in which positions need to be specified when using
! GetSingleSiteCloud.

  Implicit None
  ! Argument list:
  Type(SingleSiteFlow_), Intent(In)  :: SingleSiteFlow       ! State of a single site
                                                             ! flow module instance.
  Integer,               Intent(Out) :: nHCoords             !} Number and values of
  Integer,               Intent(Out) :: iHCoords(MaxHCoords) !} indices of horizontal
  Integer,               Intent(Out) :: nZCoords             !} and vertical coord
  Integer,               Intent(Out) :: iZCoords(MaxHCoords) !} systems.

  ! Horizontal coords needed: none.
  nHCoords    = 0
  iHCoords(1) = 1 ! Avoids compiler warning.

  ! Vertical coords needed:
  ! 1) m agl.
  nZCoords    = 1
  iZCoords(1) = SingleSiteFlow%C%iZCoords(1)

End Subroutine SingleSiteCloudCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine GetSingleSiteCloud( &
             Coords, Grids,    &
             SingleSiteFlow,   &
             Time, Position,   &
             Cloud             &
           )
! Gets cloud information from a single site flow module instance.

  Implicit None
  ! Argument list:
  Type(Coords_),          Intent(In)  :: Coords
  Type(Grids_),           Intent(In)  :: Grids
  Type (SingleSiteFlow_), Intent(In)  :: SingleSiteFlow
  Type (ShortTime_),      Intent(In)  :: Time
  Type (Position_),       Intent(In)  :: Position
  Type (Cloud_),          Intent(Out) :: Cloud
  ! Coords         :: Collection of coord systems.
  ! Grids          :: Collection of grids.
  ! SingleSiteFlow :: State of a single site flow module instance.
  ! Time           :: Current time in non-fixed met cases and the time of the met in
  !                   fixed met cases.
  ! Position       :: Coords of the point in various coord systems in Coords, with
  !                   flags to indicate whether the values are valid.
  ! Cloud          :: Cloud information.
  ! Locals:
  Type(TCoeffs_) :: TCoeffs ! Interpolation coefficients.

  ! Interpolation coefficients.
  Call GetTCoeffs(                                                     &
         Time    = Time,                                               &
         OldTime = Time2ShortTime(SingleSiteFlow%SingleSiteMet%Time1), &
         DtSec   = NInt(SingleSiteFlow%SingleSiteMet%ActualDt),        & ! $$ real?
         iOld    = 1,                                                  &
         iNew    = 2,                                                  &
         TCoeffs = TCoeffs                                             &
       )

  Cloud%Cloud = SingleSiteFlow%SingleSiteMet%Cloud1 * TCoeffs%T2 + &
                SingleSiteFlow%SingleSiteMet%Cloud2 * TCoeffs%T1

  Cloud%Cloud3d              = 0.0
  Cloud%ConCloud             = 0.0
  Cloud%ConCloudBase         = -1.0
  Cloud%ConCloudTop          = -1.0
  Cloud%TotalOrDynCloudWater = 0.0
  Cloud%TotalOrDynCloudIce   = 0.0
  If (Cloud%Cloud == 0.0) Then
    Cloud%TotalOrDynCloudBase = -1.0
    Cloud%TotalOrDynCloudTop  = -1.0
  Else
    Cloud%TotalOrDynCloudBase = 2000.0
    Cloud%TotalOrDynCloudTop  = 6000.0
  End If

  Cloud%iZCoord        = SingleSiteFlow%C%iZCoords(1) ! m agl.
  Cloud%TotalCloudFlag = .true.

End Subroutine GetSingleSiteCloud

!-------------------------------------------------------------------------------------------------------------

Subroutine SingleSiteRainCoordIndices(              &
             SingleSiteFlow,                        &
             nHCoords, iHCoords, nZCoords, iZCoords &
           )
! Returns indices of coord systems in which positions need to be specified when using
! GetSingleSiteRain.

  Implicit None
  ! Argument list:
  Type(SingleSiteFlow_), Intent(In)  :: SingleSiteFlow       ! State of a single site
                                                             ! flow module instance.
  Integer,               Intent(Out) :: nHCoords             !} Number and values of
  Integer,               Intent(Out) :: iHCoords(MaxHCoords) !} indices of horizontal
  Integer,               Intent(Out) :: nZCoords             !} and vertical coord
  Integer,               Intent(Out) :: iZCoords(MaxHCoords) !} systems.

  ! Horizontal coords needed: none.
  nHCoords    = 0
  iHCoords(1) = 1 ! Avoids compiler warning.

  ! Vertical coords needed: none
  nZCoords    = 0
  iZCoords(1) = 1 ! Avoids compiler warning.

End Subroutine SingleSiteRainCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine GetSingleSiteRain( &
             Coords, Grids,   &
             SingleSiteFlow,  &
             Time, Position,  &
             Rain             &
           )
! Gets rain information from a single site flow module instance.

  Implicit None
  ! Argument list:
  Type(Coords_),          Intent(In)  :: Coords
  Type(Grids_),           Intent(In)  :: Grids
  Type (SingleSiteFlow_), Intent(In)  :: SingleSiteFlow
  Type (ShortTime_),      Intent(In)  :: Time
  Type (Position_),       Intent(In)  :: Position
  Type (Rain_),           Intent(Out) :: Rain
  ! Coords         :: Collection of coord systems.
  ! Grids          :: Collection of grids.
  ! SingleSiteFlow :: State of a single site flow module instance.
  ! Time           :: Current time in non-fixed met cases and the time of the met in
  !                   fixed met cases.
  ! Position       :: Coords of the point in various coord systems in Coords, with
  !                   flags to indicate whether the values are valid.
  ! Rain           :: Rain information.
  ! Locals:
  Type(TCoeffs_) :: TCoeffs ! Interpolation coefficients.

  ! Interpolation coefficients.
  Call GetTCoeffs(                                                     &
         Time    = Time,                                               &
         OldTime = Time2ShortTime(SingleSiteFlow%SingleSiteMet%Time1), &
         DtSec   = NInt(SingleSiteFlow%SingleSiteMet%ActualDt),        & ! $$ real?
         iOld    = 1,                                                  &
         iNew    = 2,                                                  &
         TCoeffs = TCoeffs                                             &
       )

  Rain%DynPpt = SingleSiteFlow%SingleSiteMet%Ppt1 * TCoeffs%T2 + &
                SingleSiteFlow%SingleSiteMet%Ppt2 * TCoeffs%T1

  Rain%ConPpt = 0.0

End Subroutine GetSingleSiteRain

!-------------------------------------------------------------------------------------------------------------

Subroutine SingleSiteReflectCoordIndices(                           &
             SingleSiteFlow, nHCoords, iHCoords, nZCoords, iZCoords &
           )
! Returns indices of coord systems in which positions need to be specified when using
! SingleSiteReflect.

  Implicit None
  ! Argument list:
  Type(SingleSiteFlow_), Intent(In)  :: SingleSiteFlow       ! State of a single site
                                                             ! flow module instance.
  Integer,               Intent(Out) :: nHCoords             !} Number and values of
  Integer,               Intent(Out) :: iHCoords(MaxHCoords) !} indices of horizontal
  Integer,               Intent(Out) :: nZCoords             !} and vertical coord
  Integer,               Intent(Out) :: iZCoords(MaxHCoords) !} systems.

  ! Horizontal coords needed: none
  nHCoords    = 0
  iHCoords(1) = 0 ! avoids compiler warning

  ! Vertical coords needed:
  ! 1) m agl.
  nZCoords    = 1
  iZCoords(1) = SingleSiteFlow%C%iZCoords(1)

End Subroutine SingleSiteReflectCoordIndices

!-------------------------------------------------------------------------------------------------------------

Subroutine SingleSiteReflect(               &
             Coords, SingleSiteFlow,        &
             Vel, OldPosition, Position, U, &
             Reflected                      &
           )
! Reflects position of particle or puff centroid using a single site flow module
! instance.

  Implicit None
  ! Argument list:
  Type(Coords_),         Intent(In)    :: Coords
  Type(SingleSiteFlow_), Intent(In)    :: SingleSiteFlow
  Logical,               Intent(In)    :: Vel
  Type(Position_),       Intent(In)    :: OldPosition
  Type(Position_),       Intent(InOut) :: Position
  Real(Std),             Intent(InOut) :: U(3)
  Logical,               Intent(Out)   :: Reflected
  ! Coords         :: Collection of coord systems.
  ! SingleSiteFlow :: State of a single site flow module instance.
  ! Vel            :: Indicates a dispersion model with velocity memory is being used.
  ! OldPosition    :: Old particle position.
  ! Position       :: Current particle position.
  ! U              :: Current particle velocity. Defined only if Vel is true.
  ! Reflected      :: Indicates a reflection has occurred.
  ! Locals:
  Integer :: iZCoordZ ! Index of m agl coord system.

  iZCoordZ = SingleSiteFlow%C%iZCoords(1)

  Reflected = .false.

  ! Velocity memory.
  If (Vel) Then

    If (Position%Z(iZCoordZ) < 0.0) Then
      U(3)                               = - U(3)
      Position%Z(iZCoordZ)               = - Position%Z(iZCoordZ)
      Position%ZValid(1:Coords%nZCoords) = .false.
      Position%ZValid(iZCoordZ)          = .true.
      Position%RhoValid                  = .false.
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

End Subroutine SingleSiteReflect

!-------------------------------------------------------------------------------------------------------------

End Module SingleSiteFlowModule



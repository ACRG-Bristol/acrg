! Module:  Flow and Flow Profile Module

Module FlowAndFlowProfileModule

! This module provides code for returning flow information from the flow modules to
! the dispersion calculation, for transfering information between the flow modules,
! and for use internally within the met and flow modules.

! This module uses the following conventions:
! Surface values are indicated by S, e.g. PS = surface pressure.
! Near surface (e.g. screen level) values are indicated by 0, e.g. T0 = near surface
!     temperature (but Z0 is roughness length and near surface height above ground is
!     ZT).
! Consistent with the above, surface layer wind direction could be written PhiS or
!     Phi0 - we adopt Phi0.
! Reference values are indicated by Ref, e.g. PRef = reference pressure for potential
!     temperature definition (PRef is defined in the global parameters module).

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: Flow_                 ! Flow properties at a particular location, plus some
                                 ! bulk boundary layer properties. Used to return flow
                                 ! information from the flow modules or, in
                                 ! conjunction with FlowField_, to transfer information
                                 ! between the flow modules.
Public  :: Cloud_                ! Cloud properties at a particular location. Used to
                                 ! return cloud information from the flow modules or,
                                 ! in conjunction with CloudField_, to transfer
                                 ! information between the flow modules.
Public  :: Rain_                 ! Rain properties at a particular location. Used to
                                 ! return rain information from the flow modules or,
                                 ! in conjunction with RainField_, to transfer
                                 ! information between the flow modules.
Public  :: Surface_              ! Surface properties at a particular location. Used to return surface
                                 ! information from the flow modules or, in conjunction with SurfaceField_, to
                                 ! transfer information between the flow modules.
Public  :: Soil_                 ! Soil properties at a particular location. Used to return soil
                                 ! information from the flow modules or, in conjunction with SoilField_, to
                                 ! transfer information between the flow modules.
Public  :: Plant_                ! Plant properties at a particular location. Used to return plant
                                 ! information from the flow modules or, in conjunction with PlantField_, to
                                 ! transfer information between the flow modules.
Public  :: ProfileData_          ! Data needed to contruct idealised analytic mean
                                 ! flow and turbulence profiles. Used internally
                                 ! within the flow modules, or to return flow
                                 ! information from the flow modules or, in
                                 ! conjunction with FlowField_, to transfer
                                 ! information between the flow modules.
Public  :: CommonAttribField_    ! The part of the attribute field common to all attribute fields. Used to
                                 ! transfer information between the flow modules.
Public  :: FlowField_            ! A flow field. Used to transfer information between the flow modules.
Public  :: CloudField_           ! A cloud field. Used to transfer information between the flow modules.
Public  :: RainField_            ! A rain field. Used to transfer information between the flow modules.
Public  :: ConvertFlow           ! Converts flow information between different coord systems.
Public  :: ReverseFlow           ! Changes the sign of the odd order velocity moments.
Public  :: MeanFlowFromFlowField ! Calculates mean flow properties at a particular
                                 ! location and time from an instance of FlowField_.
Public  :: TurbFromFlowField     ! Calculates turbulence properties at a particular
                                 ! location and time from an instance of FlowField_.
Public  :: MeanFlowProfiles      ! Calculates mean flow properties at a particular
                                 ! height from an instance of ProfileData_.
Public  :: TurbProfiles          ! Calculates turbulence properties at a particular
                                 ! height from an instance of ProfileData_.
Public  :: ModifyTurbByTerminalVelocity   ! Modifies turbulence values at a particular location
                                          ! by the sedimentation velocity to account for the
                                          ! trajectory-crossing effect and inertial effects.

!-------------------------------------------------------------------------------------------------------------

Type :: Flow_ ! Flow properties at a particular location, plus some bulk boundary layer properties. Used to
              ! return flow information from the flow modules or, in conjunction with FlowField_, to transfer
              ! information between the flow modules.

  ! Inhomogeneous turbulence quantities:
  Real(Std) :: SigUU(3)    ! Turbulent velocity covariances (currently diagonal).
  Real(Std) :: dSigUUdX(3) ! Divergence of SigUU.
  Real(Std) :: Sk          ! Vertical velocity skewness.
  Real(Std) :: dSkdZ       ! Vertical derivative of vertical velocity skewness.
  Real(Std) :: TauUU(3)    ! Lagrangian time scales (currently diagonal).
  Real(Std) :: dTauUUdZ(3) ! Z derivatives of Lagrangian timescales
  Real(Std) :: Eps         ! Dissipation rate per unit mass.

  ! Inhomogeneous eddy diffusivities:
  Real(Std) :: K(3)    ! Eddy diffusivity (currently diagonal).
  Real(Std) :: dKdX(3) ! Divergence of K.

  ! Homogeneous turbulence quantities:
  Real(Std) :: HSigUU(3) !} Homogeneous values of SigUU, TauUU, Eps and K.
  Real(Std) :: HTauUU(3) !}
  Real(Std) :: HEps      !}
  Real(Std) :: HK(3)     !}

  ! Turbulent layer information
  Integer   :: nTurbLayers                      ! Number of turbulent layers.
  Real(Std) :: SigW2TurbLayers(MaxTurbLayers)   !] Values of Sig_w^2 and tau_w
  Real(Std) :: TauWTurbLayers(MaxTurbLayers)    !] in turbulent layers.
  Real(Std) :: ZInterface(MaxTurbLayers - 1)    ! Heights of interfaces of turbulent layers.

  ! Unresolved mesoscale motion quantities:
  Real(Std) :: SigUUM ! Horizontal variance for unresolved mesoscale motions.
  Real(Std) :: TauUUM ! Horizontal Lagrangian time scale for unresolved mesoscale motions.

  ! Mean quantities:
  Real(Std) :: U(3)      ! Mean velocity.
  Real(Std) :: EtaDot    ! Eta dot.
  Real(Std) :: dUdT(3)   ! Mean rate of change of velocity.
  Real(Std) :: T         ! Temperature (K).
  Real(Std) :: Theta     ! Potential temperature.
  Real(Std) :: Q         ! Specific humidity.
  Real(Std) :: P         ! Pressure (Pa).
  Real(Std) :: Rho       ! Density.
  Real(Std) :: dRhodX(3) ! Density gradient.
  Real(Std) :: FLPa      ! Pressure of the vertical level nearer to the freezing point.

  ! Boundary layer characteristics:
  Real(Std) :: Z0        ! Roughness length.
  Real(Std) :: UStar     ! Friction velocity.
  Real(Std) :: RecipLMO  ! 1/Monin-Obukhov length (negative for HeatFlux > 0 and set to +/-200000 in calms).
  Real(Std) :: H         ! Boundary layer depth.
  Real(Std) :: WStar     ! Convective velocity scale if HeatFlux > 0, zero if HeatFlux <= 0.
  Real(Std) :: WT        ! Surface temperature flux.
  Real(Std) :: WQ        ! Surface specific humidity flux.
  Real(Std) :: T0        ! Near surface temperature (K).
  Real(Std) :: PS        ! Surface pressure (Pa).
  Real(Std) :: PSeaLevel ! Mean sea level pressure (Pa).

  ! Topography:
  Real(Std) :: ZS          ! Z coord value at the surface.
  Real(Std) :: Topog       ! Topographic height above sea level.
  Real(Std) :: dTopogdX(2) ! Horizontal gradient of topographic height above sea level.
  Logical   :: SmoothTopog ! Indicates topography is smooth (necessary for dTopogdX to be meaningful).

  ! Time step and puff size limits:
  Real(Std) :: MaxDT  ! Order of magnitude of desired maximum change in T over time-step. This is used to help
                      ! determine the time-step, although the dispersion scheme can choose to ignore it.
  Real(Std) :: MaxDZ  ! Order of magnitude of desired maximum change in Z over time-step. This is used to help
                      ! determine the time-step, although the dispersion scheme can choose to ignore it.
  Real(Std) :: DeltaI ! Limit on Delta due to inhomogeneity.

  ! Coord systems:
  Integer   :: iHCoord !} Indices in Coords (external to this module) of the coord systems in which the above
  Integer   :: iZCoord !} quantities are expressed. The horizontal coord system is used only to define the
                       !} orientation, with lengths expressed in metres. The vertical coord system must be
                       !} height based and use units of metres at the ground. These frames of reference are
                       !} also used to calculate the particle/puff evolution. H and ZS are expressed in the
                       !} vertical system defined by iZCoord (but P, Rho, Rho in dRhodX, and Topog use metres
                       !} as the length unit, with Topog being relative to sea level). The vertical coordinate
                       !} system must be height above sea if the topography isn't smooth.
  
  ! Logical flag for running backwards
  Logical   :: Backwards 
End Type Flow_

!-------------------------------------------------------------------------------------------------------------

Type :: Cloud_ ! Cloud properties at a particular location. Used to return cloud
               ! information from the flow modules or, in conjunction with
               ! CloudField_, to transfer information between the flow modules.
  Real(Std) :: Cloud3d              ! 3-d cloud amount (fraction).
  Real(Std) :: Cloud                ! Cloud amount (fraction).
  Real(Std) :: ConCloud             ! Convective cloud amount (fraction).
  Real(Std) :: ConCloudBase         !} Convective cloud base and top. Values < 0 
  Real(Std) :: ConCloudTop          !} indicate no cloud.
  Real(Std) :: ConCloudBasePa       !] Convective cloud base and top in Pa (for convection scheme).
  Real(Std) :: ConCloudTopPa        !] Values < 0 indicate no cloud.
  Real(Std) :: TotalOrDynCloudWater ! Total or dynamic cloud liquid water (kg/kg).
  Real(Std) :: TotalOrDynCloudIce   ! Total or dynamic cloud ice (kg/kg).
  Real(Std) :: TotalOrDynCloudBase  !} Total or dynamic cloud base and top. Values < 0 
  Real(Std) :: TotalOrDynCloudTop   !} indicate no cloud.
  Integer   :: iZCoord              ! Index in Coords (external to this module) of the
                                    ! vertical coord system in which the above quantities are
                                    ! expressed. The vertical coord system must be height
                                    ! based and use units of metres at the ground. This frame
                                    ! of reference are also used to calculate the
                                    ! particle/puff evolution.
  Logical   :: TotalCloudFlag       ! Flag denoting total (true) or dynamic (false) 
                                    ! cloud water/ice, cloud base/top

  ! concloudbase/top now < 0.0 for no cloud. $$ check this is understood where used
End Type Cloud_

!-------------------------------------------------------------------------------------------------------------

Type :: Rain_ ! Rain properties at a particular location. Used to return rain
              ! information from the flow modules or, in conjunction with RainField_,
              ! to transfer information between the flow modules.
  Real(Std) ConPpt ! Convective precipitation rate (mm/hr).
  Real(Std) DynPpt ! Dynamic precipitation rate (mm/hr).
End Type Rain_

!-------------------------------------------------------------------------------------------------------------

Type :: Surface_ ! Surface properties at a particular location. Used to return surface
                 ! information from the flow modules or, in conjunction with SurfaceField_,
                 ! to transfer information between the flow modules.
  Real(Std) :: LandUseFracs(9)        ! Land use fractions. $$ Moses types hardwired.
                                      ! Land use types:
                                      ! 1 = Broadleaf tree (Deciduous)
                                      ! 2 = Needleleaf tree (Coniferous)
                                      ! 3 = C3 Grass (temperate regions: temperate grasslands, prairie and
                                      !     crops)
                                      ! 4 = C4 Grass (tropical regions: savannah and other tropical
                                      !     grasslands)
                                      ! 5 = Shrubs
                                      ! 6 = Urban
                                      ! 7 = Inland water
                                      ! 8 = Bare soil
                                      ! 9 = Ice.
  Real(Std) :: SoilMoisture           ! Soil moisture of surface layer (kg/m^2).
  Real(Std) :: LandFrac               ! Land Fraction (0=sea, 1=land)
End Type Surface_

!-------------------------------------------------------------------------------------------------------------

Type :: Soil_ ! Soil properties at a particular location. Used to return soil
              ! information from the flow modules or, in conjunction with SoilField_,
              ! to transfer information between the flow modules.
  Real(Std) :: ClayFrac               ! Clay mass fraction.
  Real(Std) :: SoilFracs(6)           ! Soil particle mass fractions in particle size ranges. $$ 6 hardwired
  Real(Std) :: PreferentialSourceTerm ! Indication of soil erodability ($$ units??). $$ not yet used.
End Type Soil_

!-------------------------------------------------------------------------------------------------------------

Type :: Plant_ ! Plant properties at a particular location. Used to return plant
               ! information from the flow modules or, in conjunction with PlantField_,
               ! to transfer information between the flow modules.
  Real(Std) :: CanopyHeight(5)        ! Canopy height of 5 plant types (m).
  Real(Std) :: LAI(5)                 ! Leaf area index of 5 plant types.
  Real(Std) :: StomataConduct(5)      ! Stomatal conductance of 5 plant types (m/s).
  Real(Std) :: CanopyWater(5)         ! Canopy water of 5 plant types (kg/m^2)
End Type Plant_

!-------------------------------------------------------------------------------------------------------------

Type :: ProfileData_ ! Data needed to construct idealised analytic mean flow and
                     ! turbulence profiles. Used internally within the flow modules,
                     ! or to return flow information from the flow modules or, in
                     ! conjunction with FlowField_, to transfer information between
                     ! the flow modules.
  Real(Std) :: Z0              ! Roughness length.
  Real(Std) :: Phi0            ! Surface layer wind direction.
  Real(Std) :: PhiG            ! Geostrophic wind direction.
  Real(Std) :: DeltaPhi        ! Geostrophic wind direction minus surface layer wind
                               ! direction.
  Real(Std) :: UStar           ! Friction velocity.
  Real(Std) :: UG              ! Geostrophic wind speed.
  Real(Std) :: RecipLMO        ! 1/Monin-Obukhov length (negative for HeatFlux > 0 and
                               ! set to +/-200000 in calms).
  Real(Std) :: H               ! Boundary layer depth.
  Real(Std) :: WStar           ! Convective velocity scale if HeatFlux > 0, zero if
                               ! HeatFlux <= 0.
  Real(Std) :: SigUUM          ! Horizontal variance for unresolved mesoscale motions.
  Real(Std) :: TauUUM          ! Horizontal Lagrangian time scale for unresolved mesoscale motions.
  Real(Std) :: SigU2HPlus      !} Velocity variances above boundary layer (U here means
  Real(Std) :: SigW2HPlus      !} U and V).
  Real(Std) :: TauUHPlus       !] Lagrangian time scales above boundary layer (U here
  Real(Std) :: TauWHPlus       !] means U and V).
  Real(Std) :: ZT              ! Height of T0 and Q0 above ground.
  Real(Std) :: T0              ! Near surface temperature (K).
  Real(Std) :: WT              ! Surface temperature flux.
  Real(Std) :: TStar           ! Surface temperature scale (negative for HeatFlux > 0
                               ! and set to +/-200000 in calms).
  Real(Std) :: THPlus          ! Temperature just above the boundary layer (K).
  Real(Std) :: dTdZHPlus       ! d(temperature)/dz just above the boundary layer.
  Real(Std) :: Q0              ! Near surface specific humidity.
  Real(Std) :: WQ              ! Surface specific humidity flux.
  Real(Std) :: QStar           ! Surface specific humidity scale (negative for
                               ! LatentHeatFlux > 0 and set to +/-200000 in calms).
  Real(Std) :: RHHPlus         ! Relative humidity (percent) just above the boundary
                               ! layer.
  Real(Std) :: dRHdZHPlus      ! d(Relative Humidity (percent))/dz above the boundary
                               ! layer.
  Real(Std) :: PS              ! Surface pressure (Pa).
  Real(Std) :: PHPlus          ! Pressure just above the boundary layer (Pa).
  Logical   :: Canopy = .false. ! Indicates canopy effects are to be simulated.
  Real(Std) :: LambdaP         !}
  Real(Std) :: D               !}
  Real(Std) :: HC              !} For use with Urban wind profiles and turbulence
  Real(Std) :: UAtHC           !} $$ give full definitions.
  Real(Std) :: LExp            !} $$ If used with single site, need to think about where initialised.
  Real(Std) :: ZHat            !} $$ Update comments below on MeanFlow and Turb in relation to canopy.
  Real(Std) :: UStarG          !} $$ Useable with crop canopies?
  Real(Std) :: Z0G             !}
  Real(Std) :: LC              !}
  ! Coord systems:
  Integer   :: iHCoord    !} Indices in Coords (external to this module) of the
  Integer   :: iZCoord    !} coord systems in which the above quantities are
                          !} expressed. The horizontal coord system is used only
                          !} to define the wind directions. Wind directions are
                          !} the direction of the wind vector in radians measured
                          !} anti-clockwise from the x direction. The vertical
                          !} coord system must be height above the ground in
                          !} metres.
  ! Possible uses of the data:
  Logical   :: MeanFlow   ! Indicates that ProfileData can be used with the
                          ! MeanFlowProfiles routine to generate mean profiles in
                          ! the column above/below the location. The following
                          ! variables have to be set:
                          ! $$ List variables
                          ! $$ Note constraints on the physics variables (e.g. Z0
                          !    non-zero). H can be zero!
  Logical   :: Turb       ! Indicates that ProfileData can be used with the
                          ! TurbProfiles routine to generate turbulence profiles
                          ! in the column above/below the location. The following
                          ! variables have to be set:
                          !     Z0
                          !     Phi0
                          !     UStar
                          !     RecipLMO
                          !     H
                          !     WStar
                          !     SigUUM
                          !     TauUUM
                          !     SigU2HPlus
                          !     SigW2HPlus
                          !     TauUHPlus
                          !     TauWHPlus
                          !     iHCoord
                          !     iZCoord
                          !     MeanFlow
                          !     Turb
                          ! Note all the physics variables except Z0 can be zero.
End Type ProfileData_

!-------------------------------------------------------------------------------------------------------------

Type :: CommonAttribField_ ! The part of the attribute field common to all attribute fields. Used to transfer
                           ! information between the flow modules.
  Integer                  :: iHGrid
  Integer                  :: iZGrid
  Type(Time_)              :: Dt
  Logical                  :: UseTwoTimes
  Type(Time_)              :: Time1
  Type(Time_)              :: Time2
  Type(ShortTime_)         :: ShortTime1
  Type(ShortTime_)         :: ShortTime2
  Logical                  :: Valid
  ! iHGrid      :} Indices of grid on which attribute information is needed.
  ! iZGrid      :}
  ! Dt          :: Maximum time interval between Time1 and Time2.
  ! UseTwoTimes :: Indicates that information is needed at Time2 as well as Time1. Must be false for fixed
  !                met.
  ! Time1       :} Times at which attribute information is provided. Time1 is the time at which the data is
  ! Time2       :} requested and Time2 is Time1 + Dt unless data is not available as far ahead as this time.
  ! STime1      :] ShortTime_ versions of Time1 and Time2.
  ! STime2      :]
  ! Valid       :: Indicates the attribute information is valid.

  ! The flow module instance requiring attribute information from other flow modules must complete iHGrid,
  ! iZGrid, Dt, UseTwoTimes, and allocate the pointers or associate them with existing storage. This is done
  ! in the flow module instance's ????FlowReqs routine. When (if) the attribute information has been added,
  ! Valid is set to .true..

End Type CommonAttribField_

!-------------------------------------------------------------------------------------------------------------

Type :: FlowField_ ! A flow field. Used to transfer information between the flow modules.
  Type(CommonAttribField_)         :: C
  Type(Flow_),             Pointer :: Flow(:,:,:,:)
  Type(ProfileData_),      Pointer :: ProfileData(:,:,:,:)
  Logical                          :: NoHVariation
  ! C            :: The part of the attribute field common to all attribute fields.
  ! Flow         :: Array of flow information for TauMin = ZMin = 0 (i.e. with no lower limit on Lagrangian
  !                 time scales and no height below which the flow properties are held constant).
  ! ProfileData  :: Array of profile data which may be useable to generate profiles (depending on the values
  !                 of ProfileData%MeanProfile and ProfileData%TurbProfile).
  ! NoHVariation :: Indicates flow is uniform in the horizontal.

End Type FlowField_

!-------------------------------------------------------------------------------------------------------------

Type :: CloudField_ ! A cloud field. Used to transfer information between the flow modules.
  Type(CommonAttribField_)         :: C
  Type(Cloud_),            Pointer :: Cloud(:,:,:,:)
  ! C     :: The part of the attribute field common to all attribute fields.
  ! Cloud :: Array of cloud information.

End Type CloudField_

!-------------------------------------------------------------------------------------------------------------

Type :: RainField_ ! A rain field. Used to transfer information between the flow modules.
  Type(CommonAttribField_)         :: C
  Type(Rain_),             Pointer :: Rain(:,:,:,:)
  ! C    :: The part of the attribute field common to all attribute fields.
  ! Rain :: Array of rain information.

End Type RainField_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine ConvertFlow(Coords, iHCoord, iZCoord, Position, Flow)
! Converts flow information between different coord systems.

! $$ eventually will need to return an error code if e.g. ConvertToZ fails

  Implicit None
  ! Argument List:
  Type(Coords_),   Intent(In)    :: Coords   ! Collection of coord systems.
  Integer,         Intent(In)    :: iHCoord  !} Indices of desired horizontal and
  Integer,         Intent(In)    :: iZCoord  !} vertical coord systems.
  Type(Position_), Intent(InOut) :: Position ! Coords of location in various coord
                                             ! systems in Coords, with flags to
                                             ! indicate whether the values are valid.
  Type(Flow_),     Intent(InOut) :: Flow     ! Flow information.
  ! Locals:
  Real(Std) :: Angle ! Local rotation angle of the desired horizontal coord system
                     ! from Flow's current horizontal coord system (in radians, and
                     ! positive for anticlockwise rotations).
  Real(Std) :: U(2)  ! Copy of Flow%U(1:2).

  ! Currently only Mean flow converted in the horizontal $$

  If (iHCoord /= Flow%iHCoord) Then

    If (.not.Position%XYValid(Flow%iHCoord)) Then
      Call ConvertToH(Coords, iHCoord, Position)
    End If
    Call CalcHAngle(                                              &
           Coords%HCoords(Flow%iHCoord), Coords%HCoords(iHCoord), &
           Position%XY(:, Flow%iHCoord),                          &
           Angle                                                  &
         )
    U(1) = Flow%U(1)
    U(2) = Flow%U(2)
    Flow%U(1) = U(1)*Cos(Angle) + U(2)*Sin(Angle)
    Flow%U(2) = U(2)*Cos(Angle) - U(1)*Sin(Angle)

    Flow%iHCoord = iHCoord

  End If

  If (iZCoord /= Flow%iZCoord) Then

    Flow%iZCoord = iZCoord

  End If

End Subroutine ConvertFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine ReverseFlow(Flow)
! Changes the sign of the odd order velocity moments.

  Implicit None
  ! Argument List:
  Type(Flow_), Intent(InOut) :: Flow ! Flow information.

  Flow%Backwards = .Not. Flow%Backwards
  Flow%Sk        = -Flow%Sk
  Flow%dSkdZ     = -Flow%dSkdZ
  Flow%U(:)      = -Flow%U(:)

End Subroutine ReverseFlow

!-------------------------------------------------------------------------------------------------------------

Subroutine MeanFlowFromFlowField(                  &
             Coords, Grids, Time, X, Z, FlowField, &
             Moisture,                             &
             Speed, Flow,                          &
             HCoeffs, ZCoeffs, TCoeffs             &
           )
! Calculates mean flow properties at a particular location and time from an instance
! of FlowField_. The coord systems used for X, for the grid in FlowField, and within
! the Flow and ProfileData parts of FlowField, must all be the same. The same coord
! systems are used to return the flow properties. Returns the 'mean quantities',
! 'boundary layer characteristics', 'topography', 'time step and puff size limits' and
! 'coord systems' parts of Flow.

  Implicit None
  ! Argument List:
  Type(Coords_),    Intent(In)    :: Coords    ! Collection of coord systems.
  Type(Grids_),     Intent(In)    :: Grids     ! Collection of grids.
  Type(ShortTime_), Intent(In)    :: Time      ! Time.
  Real(Std),        Intent(In)    :: X(3)      ! Location.
  Real(Std),        Intent(In)    :: Z         ! Height above ground.
  Type(FlowField_), Intent(In)    :: FlowField ! FlowField.
  Logical,          Intent(In)    :: Moisture  ! Indicates Q is required.
  Real(Std),        Intent(Out)   :: Speed     ! Wind speed.
  Type(Flow_),      Intent(InOut) :: Flow      ! Flow information.
  Type(HCoeffs_),   Intent(InOut) :: HCoeffs   !} Interpolation coeeficients.
  Type(ZCoeffs_),   Intent(InOut) :: ZCoeffs   !}
  Type(TCoeffs_),   Intent(InOut) :: TCoeffs   !}
  ! Locals:
  Real(Std) :: dZdX(3) ! Rate of change of height above ground with respect to the
                       ! coords X.

  ! Need to convert coords used in Flow and ProfileData parts of FlowField before use.
  ! Routines to do this needed in this module, but called from Flows. $$

  ! Need to set and use NoHVariation for efficiency. $$

  ! The SABuilding block avoids the need to initialise certain parts of Coords, Grids
  ! and FlowField when used with the Stand-Alone Building Model.
# ifdef SABuilding

    Call MeanFlowProfiles(                    &
           Z, FlowField%ProfileData(1,1,1,1), &
           Moisture, .false.,                 &
           Speed, Flow                        &
         )

# else

    ! Single-column, single-time data.
    If (                                             &
      .not. FlowField%C%UseTwoTimes            .and. &
      Grids%HGrids(FlowField%C%iHGrid)%nX == 1 .and. &
      Grids%HGrids(FlowField%C%iHGrid)%nY == 1       &
    ) Then

      ! Valid profile data.
      ! $$ To use profile data require coords to be height above ground. add check
      If (FlowField%ProfileData(1,1,1,1)%MeanFlow) Then
        Call MeanFlowProfiles(                    &
               Z, FlowField%ProfileData(1,1,1,1), &
               Moisture, .false.,                 &
               Speed, Flow                        &
             )
        ! Calculate vertical velocity in right coord system.
        Call CalcdZdXZBased(                                           &
               Coords%ZCoords(FlowField%Flow(1,1,1,1)%iZCoord),        &
               Coords%ZCoords(FlowField%ProfileData(1,1,1,1)%iZCoord), &
               Z,                                                      &
               FlowField%Flow(1,1,1,1)%Topog,                          &
               FlowField%Flow(1,1,1,1)%dTopogdX,                       &
               dZdX                                                    &
        )
        Flow%U(3) = Flow%U(1)*dZdX(1) + Flow%U(2)*dZdX(2)

        ! $$ adjust coord system values

        Flow%ZS          = 0.0
        Flow%Topog       = FlowField%Flow(1,1,1,1)%Topog
        Flow%dTopogdX(:) = FlowField%Flow(1,1,1,1)%dTopogdX(:)
        Flow%SmoothTopog = FlowField%Flow(1,1,1,1)%SmoothTopog

        Flow%MaxDZ  = FlowField%Flow(1,1,1,1)%MaxDZ
        Flow%DeltaI = FlowField%Flow(1,1,1,1)%DeltaI

      ! Invalid profile data.
      Else
        ! $$

      End If

    ! Multi-column/time data.
    Else
      ! $$

    End If

# endif

End Subroutine MeanFlowFromFlowField

!-------------------------------------------------------------------------------------------------------------

Subroutine TurbFromFlowField(                      &
             Coords, Grids, Time, X, Z, FlowField, &
             Inhomog, Homog,                       &
             UEqV, TauMin, ZMin,                   &
             Flow,                                 &
             HCoeffs, ZCoeffs, TCoeffs             &
           )
! Calculates turbulence properties at a particular location and time from an instance
! of FlowField_. The coord systems used for X, for the grid in FlowField, and within
! the Flow and ProfileData parts of FlowField, must all be the same. The same coord
! systems are used to return the flow properties. Returns the 'inhomogeneous
! turbulence quantities', 'inhomogeneous eddy diffusivities', 'homogeneous turbulence
! quantities', 'unresolved mesoscale motion quantities', 'boundary layer characteristics' 
! (only Z0, UStar, RecipLMO, H and WStar), 'topography', 'time step and puff size limits' 
! and 'coord systems' parts of Flow.

  Implicit None
  ! Argument list:
  Type(Coords_),    Intent(In)    :: Coords      ! Collection of coord systems.
  Type(Grids_),     Intent(In)    :: Grids       ! Collection of grids.
  Type(ShortTime_), Intent(In)    :: Time        ! Time.
  Real(Std),        Intent(In)    :: X(3)        ! Location.
  Real(Std),        Intent(In)    :: Z           ! Height above ground.
  Type(FlowField_), Intent(In)    :: FlowField   ! FlowField.
  Logical,          Intent(In)    :: Inhomog     ! Inhomogeneous profiles needed.
  Logical,          Intent(In)    :: Homog       ! Homogeneous profiles needed.
  Logical,          Intent(In)    :: UEqV        ! Indicates SigU and TauU are set
                                                 ! equal to SigV and TauV.
  Real(Std),        Intent(In)    :: TauMin      ! Minimum Lagrangian time scale.
  Real(Std),        Intent(In)    :: ZMin        ! Minimum height used for calculating
                                                 ! turbulence statistics.
  Type(Flow_),      Intent(InOut) :: Flow        ! Flow information.
  Type(HCoeffs_),   Intent(InOut) :: HCoeffs     !} Interpolation coeeficients.
  Type(ZCoeffs_),   Intent(InOut) :: ZCoeffs     !}
  Type(TCoeffs_),   Intent(InOut) :: TCoeffs     !}

  ! The SABuilding block avoids the need to initialise certain parts of Coords, Grids
  ! and FlowField when used with the Stand-Alone Building Model.
# ifdef SABuilding

    Call TurbProfiles(X(3), FlowField%ProfileData(1,1,1,1), &
                      Inhomog, Homog, UEqV, TauMin, ZMin,   &
                      Flow)

# else

    ! $$ To use profile data require coords to be height above ground.
    ! need options for invalid profile data + horizontal interp
    Call TurbProfiles(X(3), FlowField%ProfileData(1,1,1,1), &
                      Inhomog, Homog, UEqV, TauMin, ZMin,   &
                      Flow)

    Flow%ZS          = 0.0
    Flow%Topog       = FlowField%Flow(1,1,1,1)%Topog
    Flow%dTopogdX(:) = FlowField%Flow(1,1,1,1)%dTopogdX(:)
    Flow%SmoothTopog = FlowField%Flow(1,1,1,1)%SmoothTopog

    Flow%MaxDZ  = FlowField%Flow(1,1,1,1)%MaxDZ
    Flow%DeltaI = FlowField%Flow(1,1,1,1)%DeltaI

# endif

End Subroutine TurbFromFlowField

!-------------------------------------------------------------------------------------------------------------

Subroutine MeanFlowProfiles( &
             Z, ProfileData, &
             Moisture,       &
             HMinus,         &
             Speed, Flow     &
           )
! Calculates mean flow properties at a particular height from an instance of
! ProfileData_. Returns the 'mean quantities', 'boundary layer characteristics' and
! 'coord systems' parts of Flow.

  Implicit None
  ! Argument List:
  Real(Std),          Intent(In)    :: Z           ! Height above ground.
  Type(ProfileData_), Intent(In)    :: ProfileData ! Data needed to construct
                                                   ! idealised analytic mean flow and
                                                   ! turbulence profiles.
  Logical,            Intent(In)    :: Moisture    ! Indicates Q is required.
  Logical,            Intent(In)    :: HMinus      ! Returns values using the 'below
                                                   ! H' formulae. This enables values
                                                   ! immediately below the boundary
                                                   ! layer top to be obtained without
                                                   ! relying on equality testing of
                                                   ! real numbers.
  Real(Std),          Intent(Out)   :: Speed       ! Wind Speed.
  Type(Flow_),        Intent(InOut) :: Flow        ! Flow information.
  ! Local parameters.
  Real(Std), Parameter :: AA = 0.7  !} Constants for stable profile.
  Real(Std), Parameter :: BB = 0.75 !}
  Real(Std), Parameter :: CC = 5.0  !}
  Real(Std), Parameter :: DD = 0.35 !}
  ! Locals:
  Real(Std) :: ZByL      ! (Z + Z0)/L.
  Real(Std) :: Z0ByL     ! Z0/L.
  Real(Std) :: ZTByL     ! (ZT + Z0)/L.
  Real(Std) :: PhiZ      ! Non-dimensional gradient at Z + Z0.
  Real(Std) :: X         !} 1/Phi at Z + Z0, Z0 and ZT + Z0.
  Real(Std) :: X0        !}
  Real(Std) :: XT        !}
  Real(Std) :: PsiZ      !] Integral of (Phi(Z/L) - 1.0)/Z at Z + Z0, Z0 and ZT + Z0.
  Real(Std) :: Psi0      !]
  Real(Std) :: PsiT      !]
  Real(Std) :: DPsi      ! PsiZ - Psi0 or PsiZ - PsiT.
  Real(Std) :: Phi       ! Wind direction.
  Real(Std) :: dPhidZ    ! d(wind direction)/dZ within the boundary layer.
  Real(Std) :: Profile   ! Factor reflecting profile shape.
  Real(Std) :: RH        ! Relative humidity (%).
  Real(Std) :: ThetaP    ! Potential temperature using surface pressure as reference
                         ! pressure.
  Real(Std) :: dThetaPdZ ! d(ThetaP)/dZ.
  Real(Std) :: ZAtLimit  !} Values of height above ground and pressure at which low
  Real(Std) :: PAtLimit  !} temperature fixup starts.

  ! Return neutral wind speed profile for use with the Stand-Alone Building Model.
# ifdef SABuilding
    Speed = (ALog((Z + ProfileData%Z0)/ProfileData%Z0))*ProfileData%UStar/VK
    Return
# endif

# ifdef ExtraChecks
    If (.not.ProfileData%MeanFlow) Then
      Call Message('UNEXPECTED FATAL ERROR in MeanFlowProfiles', 4)
    End If
# endif

  ! Check for ProfileData%H zero and ProfileData%DeltaPhi non-zero. The code can fail
  ! if this is so. (We don't check that inputs generally are sensible, but we make
  ! this check since the routine might be thought to work for this case.)
# ifdef ExtraChecks
    If (ProfileData%H == 0.0 .and. ProfileData%DeltaPhi /= 0.0) Then
      Call Message('UNEXPECTED FATAL ERROR in MeanFlowProfiles', 4)
    End If
# endif

  ! Z < H.
  If (Z < ProfileData%H .or. HMinus) Then

    ! (Z + Z0)/L, Z0/L and (ZT + Z0)/L.
    ZByL  = (Z + ProfileData%Z0) * ProfileData%RecipLMO
    Z0ByL = ProfileData%Z0 * ProfileData%RecipLMO
    ZTByL = (ProfileData%ZT + ProfileData%Z0) * ProfileData%RecipLMO

    ! Rate of change of wind direction (avoiding difficulties with ProfileData%H = 0).
    If (ProfileData%DeltaPhi == 0.0) Then
      dPhidZ = 0.0
    Else
      dPhidZ = ProfileData%DeltaPhi / ProfileData%H
    End If

    ! Stable Case.
    If (ProfileData%RecipLMO > 0.0) Then

      ! Wind profile.
      If (ProfileData%UStar == 0.0) Then ! Add calm flag to ProfileData to avoid real=0 tests $$
        Profile = 0.0
      Else
        PsiZ    = AA*ZByL  + BB*(ZByL  - CC/DD)*Exp(-DD*ZByL )
        Psi0    = AA*Z0ByL + BB*(Z0ByL - CC/DD)*Exp(-DD*Z0ByL)
        Profile = ALog((Z + ProfileData%Z0) / ProfileData%Z0) + PsiZ - Psi0
      End If
      Speed     = Profile * ProfileData%UStar / VK
      Phi       = ProfileData%Phi0 + dPhidZ * Z
      Flow%U(1) = Speed * Cos(Phi)
      Flow%U(2) = Speed * Sin(Phi)
      Flow%U(3) = 0.0

      ! Potential temperature and temperature profile.
      If (ProfileData%UStar == 0.0) Then
        Profile = 0.0 ! could do better? $$
      Else
        PsiT    = AA*ZTByL + BB*(ZTByL - CC/DD)*Exp(-DD*ZTByL)
        Profile = ALog((Z + ProfileData%Z0) / (ProfileData%ZT + ProfileData%Z0)) &
                  + PsiZ - PsiT
      End If
      ThetaP     = Profile * ProfileData%TStar / VK +             &
                   ProfileData%T0 + (Gravity/Cp) * ProfileData%ZT
      Flow%T     = ThetaP - (Gravity/Cp)*Z
      Flow%Theta = ThetaP * (PRef / ProfileData%PS) ** (GasConstant/Cp)

      ! Pressure and density profiles.
      Flow%P   = PRef * (Flow%T / Flow%Theta) ** (Cp/GasConstant)
      Flow%Rho = Flow%P / (GasConstant * Flow%T)

      ! Density gradient.
      If (ProfileData%UStar == 0.0) Then
        dThetaPdZ = 0.0 ! could do better? $$
      Else
        PhiZ = 1.0 + AA*ZByL + (BB*ZByL - BB*(ZByL - CC/DD)*DD*ZByL)*Exp(-DD*ZByL)
        dThetaPdZ = ProfileData%TStar * PhiZ / (VK * (Z + ProfileData%Z0))
      End If
      Flow%dRhodX(1:2) = 0.0
      Flow%dRhodX(3) = (Cp/GasConstant - 1.0) * (dThetaPdZ - Gravity/Cp) / Flow%T &
                       - (Cp/GasConstant) * dThetaPdZ / ThetaP
      Flow%dRhodX(3) = Flow%dRhodX(3) * Flow%Rho

      If (Moisture) Then

        ! Specific humidity.
        Flow%Q = Profile * ProfileData%QStar / VK
        Flow%Q = Flow%Q + ProfileData%Q0

        ! Limit relative humidity.
        RH = CalcRH(Flow%Q, Flow%T, Flow%P)
        If (RH <   0.0) RH =   0.0
        If (RH > 100.0) RH = 100.0
        Flow%Q = CalcQ(RH, Flow%T, Flow%P)

      End If

    ! Unstable case.
    Else

      ! Wind Profile.
      If (ProfileData%UStar == 0.0) Then
        Profile = 0.0
      Else
        X       = (1.0 - 16.0*ZByL )**0.25
        X0      = (1.0 - 16.0*Z0ByL)**0.25
        DPsi    = 2.0*(ATan(X) - ATan(X0)) -                                        &
                  ALog(((1.0 + X**2)*(1.0 + X)**2) / ((1.0 + X0**2)*(1.0 + X0)**2))
        Profile = ALog((Z + ProfileData%Z0) / ProfileData%Z0) + DPsi
      End If
      Speed     = Profile * ProfileData%UStar / VK
      Phi       = ProfileData%Phi0 + dPhidZ * Z
      Flow%U(1) = Speed * Cos(Phi)
      Flow%U(2) = Speed * Sin(Phi)
      Flow%U(3) = 0.0

      ! Potential temperature and temperature profile.
      If (ProfileData%UStar == 0.0) Then
        Profile = 0.0 ! $$ could do better if use free convective form (using WT not TStar)
      Else
        X       = (1.0 - 16.0*ZByL)**0.5
        XT      = (1.0 - 16.0*ZTByL)**0.5
        DPsi    = - ALog(((1.0 + X)**2) / ((1.0 + XT)**2))
        Profile = ALog((Z + ProfileData%Z0) / (ProfileData%ZT + ProfileData%Z0)) &
                  + DPsi
      End If
      ThetaP     = Profile * ProfileData%TStar / VK +             &
                   ProfileData%T0 + (Gravity/Cp) * ProfileData%ZT
      Flow%T     = ThetaP - (Gravity/Cp)*Z
      Flow%Theta = ThetaP * (PRef / ProfileData%PS) ** (GasConstant/Cp)

      ! Pressure and density profiles.
      Flow%P   = PRef * (Flow%T / Flow%Theta) ** (Cp/GasConstant)
      Flow%Rho = Flow%P / (GasConstant * Flow%T)

      ! Density gradient.
      If (ProfileData%UStar == 0.0) Then
        dThetaPdZ = 0.0 ! $$ could do better if use free convective form (using WT not TStar)
      Else
        dThetaPdZ = ProfileData%TStar / (VK * (Z + ProfileData%Z0) * X)
      End If
      Flow%dRhodX(1:2) = 0.0
      Flow%dRhodX(3) = (Cp/GasConstant - 1.0) * (dThetaPdZ - Gravity/Cp) / Flow%T &
                       - (Cp/GasConstant) * dThetaPdZ / ThetaP
      Flow%dRhodX(3) = Flow%dRhodX(3) * Flow%Rho

      If (Moisture) Then

        ! Specific humidity.
        Flow%Q = Profile * ProfileData%QStar / VK
        Flow%Q = Flow%Q + ProfileData%Q0

        ! Limit relative humidity.
        RH = CalcRH(Flow%Q, Flow%T, Flow%P)
        If (RH <   0.0) RH =   0.0
        If (RH > 100.0) RH = 100.0
        Flow%Q = CalcQ(RH, Flow%T, Flow%P)

      End If

    End If

  ! Z >= H.
  Else

    ! Wind Profile.
    Speed     = ProfileData%UG
    Flow%U(1) = Speed * Cos(ProfileData%PhiG)
    Flow%U(2) = Speed * Sin(ProfileData%PhiG)
    Flow%U(3) = 0.0

    ! Temperature profile.
    Flow%T = ProfileData%THPlus + ProfileData%dTdZHPlus * (Z - ProfileData%H)

    ! Normal case.
    If (Flow%T >= TAt11km .or. Flow%T >= ProfileData%THPlus) Then

      ! Pressure, density and potential temperature profiles.
      If (Abs(Flow%T / ProfileData%THPlus - 1.0) > 2.0*Sqrt(Epsilon(1.0_Std))) Then
        Flow%P = ProfileData%PHPlus *                              &
                 (Flow%T / ProfileData%THPlus) **                  &
                 (- Gravity/(GasConstant * ProfileData%dTdZHPlus))
      Else
        Flow%P = ProfileData%PHPlus *                              &
                 Exp(                                              &
                   - (Gravity/GasConstant) * (Z - ProfileData%H) / &
                   ProfileData%THPlus                              &
                 )
      End If
      Flow%Rho   = Flow%P / (GasConstant * Flow%T)
      Flow%Theta = Flow%T * (PRef / Flow%P) ** (GasConstant/Cp)

      ! Density gradient.
      Flow%dRhodX(1:2) = 0.0
      Flow%dRhodX(3) = - (Gravity * Flow%Rho) / Flow%P  &
                       - ProfileData%dTdZHPlus / Flow%T
      Flow%dRhodX(3) = Flow%dRhodX(3) * Flow%Rho

    ! Fix-up for low temperature.
    Else

      ! Temperature profile.
      Flow%T = Min(TAt11km, ProfileData%THPlus)

      ! Pressure, density and potential temperature profiles.
      If (Abs(Flow%T / ProfileData%THPlus - 1.0) > 2.0*Sqrt(Epsilon(1.0_Std))) Then
        PAtLimit = ProfileData%PHPlus *                              &
                   (Flow%T / ProfileData%THPlus) **                  &
                   (- Gravity/(GasConstant * ProfileData%dTdZHPlus))
        ZAtLimit = ProfileData%H +                                       &
                   (Flow%T - ProfileData%THPlus) / ProfileData%dTdZHPlus
        Flow%P   = PAtLimit *                                             &
                   Exp(- (Gravity/GasConstant) * (Z - ZAtLimit) / Flow%T)
      Else
        Flow%P = ProfileData%PHPlus *                                        &
                 Exp(- (Gravity/GasConstant) * (Z - ProfileData%H) / Flow%T)
      End If
      Flow%Rho   = Flow%P / (GasConstant * Flow%T)
      Flow%Theta = Flow%T * (PRef / Flow%P) ** (GasConstant/Cp)

      ! Density gradient.
      Flow%dRhodX(1:2) = 0.0
      Flow%dRhodX(3) = - Gravity * Flow%Rho**2 / Flow%P

    End If

    If (Moisture) Then

      ! Relative humidity.
      RH = ProfileData%RHHPlus + ProfileData%dRHdZHPlus * (Z - ProfileData%H)
      If (RH <   0.0) RH =   0.0
      If (RH > 100.0) RH = 100.0

      ! Specific humidity.
      Flow%Q = CalcQ(RH, Flow%T, Flow%P)

    End If

  End If

  ! Flow%dUdT.
  Flow%dUdT(:) = 0.0

  ! Boundary layer characteristics:
  Flow%Z0       = ProfileData%Z0
  Flow%UStar    = ProfileData%UStar
  Flow%RecipLMO = ProfileData%RecipLMO
  Flow%H        = ProfileData%H
  Flow%WStar    = ProfileData%WStar
  Flow%WT       = ProfileData%WT
  Flow%WQ       = ProfileData%WQ
  Flow%T0       = ProfileData%T0
  Flow%PS       = ProfileData%PS

  ! Coord systems:
  Flow%iHCoord = ProfileData%iHCoord
  Flow%iZCoord = ProfileData%iZCoord

End Subroutine MeanFlowProfiles

!-------------------------------------------------------------------------------------------------------------

Subroutine TurbProfiles(         &
             Z, ProfileData,     &
             Inhomog, Homog,     &
             UEqV, TauMin, ZMin, &
             Flow                &
           )
! Calculates turbulence properties at a particular height from an instance of
! ProfileData_. Returns the 'inhomogeneous turbulence quantities', 'inhomogeneous eddy
! diffusivities', 'homogeneous turbulence quantities', 'unresolved mesoscale motion quantities', 
! 'boundary layer characteristics' (only Z0, UStar, RecipLMO, H and WStar) and 'coord systems'
! parts of Flow.

  Implicit None
  ! Argument List:
  Real(Std),          Intent(In)    :: Z           ! Height above ground.
  Type(ProfileData_), Intent(In)    :: ProfileData ! Data needed to construct
                                                   ! idealised analytic mean flow and
                                                   ! turbulence profiles.
  Logical,            Intent(In)    :: Inhomog     ! Indicates inhomogeneous
                                                   ! quantities are required.
  Logical,            Intent(In)    :: Homog       ! Indicates homogeneous quantities
                                                   ! are required.
  Logical,            Intent(In)    :: UEqV        ! Indicates SigU and TauU are set
                                                   ! equal to SigV and TauV.
  Real(Std),          Intent(In)    :: TauMin      ! Minimum Lagrangian time scale.
  Real(Std),          Intent(In)    :: ZMin        ! Minimum height used for
                                                   ! calculating turbulence
                                                   ! statistics.
  Type(Flow_),        Intent(InOut) :: Flow        ! Flow information.
  ! Locals:
  Real(Std) :: Dummy1     !} Dummy variables needed in call to InhomogTurb.
  Real(Std) :: Dummy2     !}
  Real(Std) :: Dummy3     !}
  Real(Std) :: Dummy4     !}
  Real(Std) :: ZRightSide !] Values for height above ground which are definitely on
  Real(Std) :: ZWrongSide !] the same and on opposite sides of the boundary layer top
                          !] as Z.

# ifdef ExtraChecks
    If (.not.ProfileData%Turb) Then
      Call Message('UNEXPECTED FATAL ERROR in TurbProfiles', 4)
    End If
# endif

  If (Inhomog) Then

    Call InhomogTurb(                                               &
           Z,                                                       &
           ProfileData%Z0, ProfileData%UStar, ProfileData%RecipLMO, &
           ProfileData%H, ProfileData%WStar,                        &
           ProfileData%SigU2HPlus, ProfileData%SigW2HPlus,          &
           ProfileData%TauUHPlus, ProfileData%TauWHPlus,            &
           ProfileData%Canopy, ProfileData%LambdaP, ProfileData%D,  &
           ProfileData%HC, ProfileData%UAtHC, ProfileData%LExp,     &
           ProfileData%ZHat, ProfileData%UStarG, ProfileData%Z0G,   &
           ProfileData%LC,                                          &
           UEqV, TauMin, ZMin,                                      &
           Flow%SigUU(1), Flow%SigUU(2), Flow%SigUU(3),             &
           Flow%TauUU(1), Flow%TauUU(2), Flow%TauUU(3),             &
           Flow%K(1), Flow%K(2), Flow%K(3),                         &
           Flow%Eps, Flow%Sk,                                       &
           Dummy1, Dummy2, Flow%dSigUUdX(3),                        &
           Flow%dTauUUdZ(1), Flow%dTauUUdZ(2), Flow%dTauUUdZ(3),    &
           Dummy3, Dummy4, Flow%dKdX(3),                            &
           Flow%dSkdZ                                               &
         )

    Flow%dSigUUdX(1:2) = 0.0
    Flow%dKdX(1:2)     = 0.0

  End If

  If (Homog) Then

    Flow%nTurbLayers = 2  ! $$ change for turbulent layer scheme? Not currently clear if this module
                          ! should only support 2 layers (if so no changes needed here but changes needed
                          ! in flow modules) or should support multiple layers [probably for homog scheme
                          ! only] (in which case changes needed here).

    If (Z >= ProfileData%H) Then
      ZRightSide = 2.0*ProfileData%H
      ZWrongSide = 0.0
    Else
      ZRightSide = 0.0
      ZWrongSide = 2.0*ProfileData%H
    End If

    Call HomogTurb(                                        &
           ZWrongSide,                                     &
           ProfileData%UStar, ProfileData%RecipLMO,        &
           ProfileData%H, ProfileData%WStar,               &
           ProfileData%SigU2HPlus, ProfileData%SigW2HPlus, &
           ProfileData%TauUHPlus, ProfileData%TauWHPlus,   &
           UEqV, TauMin,                                   &
           Flow%HSigUU(1), Flow%HSigUU(2), Flow%HSigUU(3), &
           Flow%HTauUU(1), Flow%HTauUU(2), Flow%HTauUU(3), &
           Flow%HK(1), Flow%HK(2), Flow%HK(3),             &
           Flow%HEps                                       &
         )

    ! $$ Set first two elements of K3TurbLayers
    ! $$ will need to be changed for use with turbulent layer scheme?

    If (ZWrongSide >= ProfileData%H) Then
      Flow%SigW2TurbLayers(2) = Flow%HSigUU(3)
      Flow%TauWTurbLayers (2) = Flow%HTauUU(3)
    Else
      Flow%SigW2TurbLayers(1) = Flow%HSigUU(3)
      Flow%TauWTurbLayers (1) = Flow%HTauUU(3)
    EndIf

    Call HomogTurb(                                        &
           ZRightSide,                                     &
           ProfileData%UStar, ProfileData%RecipLMO,        &
           ProfileData%H, ProfileData%WStar,               &
           ProfileData%SigU2HPlus, ProfileData%SigW2HPlus, &
           ProfileData%TauUHPlus, ProfileData%TauWHPlus,   &
           UEqV, TauMin,                                   &
           Flow%HSigUU(1), Flow%HSigUU(2), Flow%HSigUU(3), &
           Flow%HTauUU(1), Flow%HTauUU(2), Flow%HTauUU(3), &
           Flow%HK(1), Flow%HK(2), Flow%HK(3),             &
           Flow%HEps                                       &
         )

    If (ZRightSide >= ProfileData%H) Then
      Flow%SigW2TurbLayers(2) = Flow%HSigUU(3)
      Flow%TauWTurbLayers (2) = Flow%HTauUU(3)
    Else
      Flow%SigW2TurbLayers(1) = Flow%HSigUU(3)
      Flow%TauWTurbLayers (1) = Flow%HTauUU(3)
    EndIf

    ! $$ Set element of ZInterface to be bl depth
    ! $$ will need to be changed for use with turbulent layer scheme?

    Flow%ZInterface(1) = ProfileData%H

    Flow%dTauUUdZ(1) = 0.0
    Flow%dTauUUdZ(2) = 0.0
    Flow%dTauUUdZ(3) = 0.0

  End If

  ! Convert from surface-layer-wind aligned coords:
  If (.not.UEqV) Then
    ! $$ use Phi0
  End If

  ! Unresolved mesoscale motions:
  Flow%SigUUM = ProfileData%SigUUM
  Flow%TauUUM = ProfileData%TauUUM

  ! Boundary layer characteristics:
  Flow%Z0       = ProfileData%Z0
  Flow%UStar    = ProfileData%UStar
  Flow%RecipLMO = ProfileData%RecipLMO
  Flow%H        = ProfileData%H
  Flow%WStar    = ProfileData%WStar

  ! Coord systems:
  Flow%iHCoord = ProfileData%iHCoord
  Flow%iZCoord = ProfileData%iZCoord

End Subroutine TurbProfiles

!-------------------------------------------------------------------------------------------------------------

Subroutine InhomogTurb(                                                  &
             Z,                                                          &
             Z0, UStar, RecipLMO, Zi, WStar,                             &
             SigU2HPlus, SigW2HPlus, TauUHPlus, TauWHPlus,               &
             Canopy, LambdaP, D, HC, UAtHC, LExp, ZHat, UStarG, Z0G, LC, &
             UEqV, TauMin, ZMin,                                         &
             SigU2, SigV2, SigW2, TauU, TauV, TauW, KX, KY, KZ, Eps, Sk, &
             dSigU2dZ, dSigV2dZ, dSigW2dZ, dTauUdZ, dTauVdZ, dTauWdZ,    &
             dKXdZ, dKYdZ, dKZdZ, dSkdZ                                  &
           )
! Calculates inhomogeneous velocity variances, Lagrangian timescales, eddy
! diffusivities, dissipation rate per unit mass, vertical velocity skewness and
! various derivatives. Here X and U are aligned with the mean flow.

  Implicit None
  ! Argument List:
  Real(Std), Intent(In)  :: Z          ! Height above ground.
  Real(Std), Intent(In)  :: Z0         ! Roughness length.
  Real(Std), Intent(In)  :: UStar      ! Friction velocity.
  Real(Std), Intent(In)  :: RecipLMO   ! 1/Monin-Obukhov length (set to +/-200000 in
                                       ! calms).
  Real(Std), Intent(In)  :: Zi         ! Boundary layer depth.
  Real(Std), Intent(In)  :: WStar      ! Convective velocity scale if HeatFlux > 0,
                                       ! zero if HeatFlux <= 0.
  Real(Std), Intent(In)  :: SigU2HPlus !} Velocity variances above boundary layer (U
  Real(Std), Intent(In)  :: SigW2HPlus !} here means U and V).
  Real(Std), Intent(In)  :: TauUHPlus  !] Lagrangian time scales above boundary layer
  Real(Std), Intent(In)  :: TauWHPlus  !] (U here means U and V).
  Logical,   Intent(In)  :: Canopy ! Canopy variables $$ define properly.
  Real(Std), Intent(In)  :: LambdaP
  Real(Std), Intent(In)  :: D
  Real(Std), Intent(In)  :: HC
  Real(Std), Intent(In)  :: UAtHC
  Real(Std), Intent(In)  :: LExp
  Real(Std), Intent(In)  :: ZHat
  Real(Std), Intent(In)  :: UStarG
  Real(Std), Intent(In)  :: Z0G
  Real(Std), Intent(In)  :: LC
  Logical,   Intent(In)  :: UEqV       ! Indicates SigU and TauU are set equal to SigV
                                       ! and TauV.
  Real(Std), Intent(In)  :: TauMin     ! Minimum Lagrangian time scale.
  Real(Std), Intent(In)  :: ZMin       ! Minimum height used for calculating turbulence
                                       ! statistics.
  Real(Std), Intent(Out) :: SigU2      !} Velocity variances.
  Real(Std), Intent(Out) :: SigV2      !}
  Real(Std), Intent(Out) :: SigW2      !}
  Real(Std), Intent(Out) :: TauU       !] Lagrangian time scales.
  Real(Std), Intent(Out) :: TauV       !]
  Real(Std), Intent(Out) :: TauW       !]
  Real(Std), Intent(Out) :: KX         !} Eddy diffusivities.
  Real(Std), Intent(Out) :: KY         !}
  Real(Std), Intent(Out) :: KZ         !}
  Real(Std), Intent(Out) :: Eps        ! Dissipation rate per unit mass.
  Real(Std), Intent(Out) :: Sk         ! Skewness of W.
  Real(Std), Intent(Out) :: dSigU2dZ   ! d(SigU2)/dZ.
  Real(Std), Intent(Out) :: dSigV2dZ   ! d(SigV2)/dZ.
  Real(Std), Intent(Out) :: dSigW2dZ   ! d(SigW2)/dZ.
  Real(Std), Intent(Out) :: dTauUdZ    ! d(TauU)/dZ.
  Real(Std), Intent(Out) :: dTauVdZ    ! d(TauV)/dZ.
  Real(Std), Intent(Out) :: dTauWdZ    ! d(TauW)/dZ.
  Real(Std), Intent(Out) :: dKXdZ      ! d(KX)/dZ.
  Real(Std), Intent(Out) :: dKYdZ      ! d(KY)/dZ.
  Real(Std), Intent(Out) :: dKZdZ      ! d(KZ)/dZ.
  Real(Std), Intent(Out) :: dSkdZ      ! d(Sk)/dZ.
  ! Locals:
  Real(Std) :: ZL             ! Local copy of Z, limited by ZMin.
  Real(Std) :: ZZ0            ! ZL + Z0.
  Real(Std) :: ZByZi          ! ZL/Zi.
  Real(Std) :: ZMD            ! ZL - D if canopy effects simulated; ZL + Z0 otherwise. 
                              ! $$ Not sure why need both ZZ0 and ZMD - ought to be able to rationalise this.
  Real(Std) :: C0             ! Constant for determining Lagrangian time scales.
  Real(Std) :: Sig2           ! Shape factor for sigma profiles.
  Real(Std) :: dSig2dZ        ! Shape factor for sigma gradient profiles.
  Real(Std) :: Tau            ! Shape factor for tau profiles.
  Real(Std) :: dTaudZ         ! Shape factor for tau gradient profiles.
  Real(Std) :: dEpsdZ         ! d(Eps)/dZ.
  Logical   :: InCanopy       ! Indicates canopy effects simulated and location is in the canopy.
  Real(Std) :: ExpDecay       ! Exponential factor for velocity decay in the canopy.
  Logical   :: UStarNotUStarG ! Indicates the dissipation rate in canopy uses UStar not UStarG.
  Real(Std) :: UC             ! Wind speed when location is in canopy.

  ! Return neutral boundary layer values with no turbulence above boundary layer for
  ! use with the Stand-Alone Building Model.
# ifdef SABuilding

    ! Calculate ZL.
    If (Z < ZMin) Then
      ZL = ZMin
    Else
      ZL = Z
    End If

    ! Calculate ZL + Z0.
    ZZ0 = ZL + Z0

    ! Case 1: Above boundary layer. We also include heights close to the boundary
    ! layer top and low values of UStar here. For these cases the full calculation
    ! will give small values of the sigma's and K's, but may involve overflows or
    ! divisions by zero. The criteria for being close to the boundary layer top or for
    ! UStar being small is somewhat arbitrary, and could be adjusted if necessary.
    If (ZL >= 0.9999 * Zi .or. UStar <= 0.0001) Then

      SigU2    = 0.0
      SigV2    = 0.0
      SigW2    = 0.0
      dSigU2dZ = 0.0
      dSigV2dZ = 0.0
      dSigW2dZ = 0.0
      Eps      = 0.0
      TauU     = 0.0
      TauV     = 0.0
      TauW     = 0.0
      dTauUdZ  = 0.0
      dTauVdZ  = 0.0
      dTauWdZ  = 0.0
      Sk       = 0.0
      dSkdZ    = 0.0
      KX       = 0.0
      KY       = 0.0
      KZ       = 0.0
      dKXdZ    = 0.0
      dKYdZ    = 0.0
      dKZdZ    = 0.0

    ! Case 2: In boundary layer.
    Else

      ZByZi = ZL / Zi

      Sig2  = UStar**2 * (1.0 - ZByZi)**1.5
      SigU2 = 6.25 * Sig2
      SigV2 = 4.0  * Sig2
      SigW2 = 1.69 * Sig2

      If (Z > ZMin) Then
        dSig2dZ = - 1.5 * UStar**2 * Sqrt(1.0 - ZByZi) / Zi
      Else
        dSig2dZ = 0.0
      End If
      dSigU2dZ = 6.25 * dSig2dZ
      dSigV2dZ = 4.0  * dSig2dZ
      dSigW2dZ = 1.69 * dSig2dZ

      Eps = UStar**3 * (1.0 - ZByZi) / (VK * ZZ0)

      Tau  = 2.0 * VK * ZZ0 * Sqrt(1.0 - ZByZi) / (5.0 * UStar)
      TauU = 6.25 * Tau
      TauV = 4.0  * Tau
      TauW = 1.69 * Tau

      If (Z > ZMin) Then
        dTaudZ  = Tau / ZZ0 - Tau * 0.5 / ((1.0 - ZByZi) * Zi)
      Else
        dTaudZ = 0.0
      End If
      dTauUdZ = 6.25 * dTaudZ
      dTauVdZ = 4.0  * dTaudZ
      dTauWdZ = 1.69 * dTaudZ

      Sk    = 0.0
      dSkdZ = 0.0

      KX    = SigU2 * TauU
      KY    = SigV2 * TauV
      KZ    = SigW2 * TauW
      dKXdZ = SigU2 * dTauUdZ + dSigU2dZ * TauU
      dKYdZ = SigV2 * dTauVdZ + dSigV2dZ * TauV
      dKZdZ = SigW2 * dTauWdZ + dSigW2dZ * TauW

    End If

    ! Impose minimum Lagrangian timescale TauMin.
    If (TauU < TauMin) Then
      TauU    = TauMin
      dTauUdZ = 0.0
      KX      = SigU2 * TauU
      dKXdZ   = SigU2 * dTauUdZ + dSigU2dZ * TauU
    End If
    If (TauV < TauMin) Then
      TauV    = TauMin
      dTauVdZ = 0.0
      KY      = SigV2 * TauV
      dKYdZ   = SigV2 * dTauVdZ + dSigV2dZ * TauV
    End If
    If (TauW < TauMin) Then
      TauW    = TauMin
      dTauWdZ = 0.0
      KZ      = SigW2 * TauW
      dKZdZ   = SigW2 * dTauWdZ + dSigW2dZ * TauW
    End If

    ! U statistics equal to V statistics.
    If (UEqV) Then
      SigU2    = SigV2
      TauU     = TauV
      dSigU2dZ = dSigV2dZ
      dTauUdZ  = dTauVdZ
      KX       = KY
      dKXdZ    = dKYdZ
    End If

    Return

# endif

  ! Calculate ZL.
  If (Z < ZMin) Then
    ZL = ZMin
  Else
    ZL = Z
  End If

  ! Calculate ZL + Z0.
  ZZ0 = ZL + Z0

  ! Calculate ZMD.
  If (Canopy) Then
    ZMD = ZL - D
  Else
    ZMD = ZL + Z0
  End If

  ! Calculate InCanopy.
  InCanopy = Canopy .and. ZL < HC

  ! Case 1: Above boundary layer. We also include heights close to the boundary layer
  ! top and low values of UStar and WStar here. For these cases the full calculation
  ! will agree with the above boundary layer values (except for convective SigU, SigV,
  ! TauU and TauV which may jump across the boundary layer top, and except for
  ! situations where the above boundary layer values are zero), but may involve
  ! overflows or divisions by zero. The criteria for being close to the boundary layer
  ! top or for UStar and WStar being small is somewhat arbitrary, and could be
  ! adjusted if necessary.
  If (ZL >= 0.9999 * Zi .or. WStar + UStar <= 0.0001) Then

    SigU2    = SigU2HPlus
    SigV2    = SigU2HPlus
    SigW2    = SigW2HPlus
    dSigU2dZ = 0.0
    dSigV2dZ = 0.0
    dSigW2dZ = 0.0
    If (TauWHPlus == 0.0) Then
      Eps    = 0.0
    Else
      Eps    = 2.0 * SigW2HPlus / (5.0 * TauWHPlus)
    End If
    TauU     = TauUHPlus
    TauV     = TauUHPlus
    TauW     = TauWHPlus
    dTauUdZ  = 0.0
    dTauVdZ  = 0.0
    dTauWdZ  = 0.0
    Sk       = 0.0
    dSkdZ    = 0.0
    KX       = SigU2 * TauU
    KY       = SigV2 * TauV
    KZ       = SigW2 * TauW
    dKXdZ    = SigU2 * dTauUdZ + dSigU2dZ * TauU
    dKYdZ    = SigV2 * dTauVdZ + dSigV2dZ * TauV
    dKZdZ    = SigW2 * dTauWdZ + dSigW2dZ * TauW

  ! Case 2: Stable conditions.
  Else If (RecipLMO >= 0.0) Then

    ZByZi = ZL / Zi

    Sig2  = UStar**2 * (1.0 - ZByZi)**1.5
    If (InCanopy) Then
      ExpDecay = Exp(- (HC - ZL) / LExp)
      Sig2 = Sig2 * Max(ExpDecay, UStarG / UStar)**2
    End If
    SigU2 = 6.25 * Sig2
    SigV2 = 4.0  * Sig2
    SigW2 = 1.69 * Sig2

    If (Z > ZMin) Then
      dSig2dZ = - 1.5 * UStar**2 * Sqrt(1.0 - ZByZi) / Zi
      If (InCanopy) Then
        If (ExpDecay > UStarG / UStar) Then
          dSig2dZ = dSig2dZ * ExpDecay**2 + Sig2 * 2.0 / LExp
        Else
          dSig2dZ = dSig2dZ * (UStarG / UStar)**2
        End If
      End If
    Else
      dSig2dZ = 0.0
    End If
    dSigU2dZ = 6.25 * dSig2dZ
    dSigV2dZ = 4.0  * dSig2dZ
    dSigW2dZ = 1.69 * dSig2dZ

    If (InCanopy) Then

      ! $$ Alternative idea - use normal formula but replace: 
      !    UStar -> UStar * ExpDecay } in Eps1
      !    ZMD   -> HC - D           }
      !    and
      !    UStar -> UStarG   } in Eps2.
      !    ZMD   -> ZL + Z0G }
      !    In convective case would scale WStar similarly?.
      ! Eps1 = (UStar * ExpDecay)**3 * (1.0 / (HC - D) + 4.0 * RecipLMO) * (1.0 - ZByZi) / VK
      ! Eps2 = UStarG**3 * (1.0 / (ZL + Z0G) + 4.0 * RecipLMO) * (1.0 - ZByZi) / VK
      ! If (Eps1 > Eps2) Then
      !   Eps = Eps1
      ! Else
      !   Eps = Eps2
      ! End If

      UStarNotUStarG = (UStar * ExpDecay)**3 / (HC - D) >= UStarG**3 / (ZL + Z0G)

      If (UStarNotUStarG) Then
        Eps = UStar**3 * (ExpDecay**3 / (HC - D) + 4.0 * RecipLMO) * (1.0 - ZByZi) / VK
      Else
        Eps = (UStarG**3 / (ZL + Z0G) + UStar**3 * 4.0 * RecipLMO) * (1.0 - ZByZi) / VK
      End If

      If (Z > ZMin) Then
        If (UStarNotUStarG) Then
          dEpsdZ = - UStar**3 * (ExpDecay**3 / (HC - D) + 4.0 * RecipLMO) / (VK * Zi)        &
                   + UStar**3 * (ExpDecay**3 / (HC - D)) * (1.0 - ZByZi) * 3.0 / (VK * LExp)
        Else
          dEpsdZ = - (UStarG**3 / (ZL + Z0G) + UStar**3 * 4.0 * RecipLMO) / (VK * Zi) &
                   - UStarG**3 * (1.0 - ZByZi) / (VK * (ZL + Z0G)**2)
        End If
      Else
        dEpsdZ = 0.0
      End If

    Else

      Eps = UStar**3 * (1.0 / ZMD + 4.0 * RecipLMO) * (1.0 - ZByZi) / VK

      If (Z > ZMin) Then
        dEpsdZ = - UStar**3 * (1.0 / ZMD + 4.0 * RecipLMO) / (VK * Zi) &
                 - UStar**3 * (1.0 - ZByZi) / (VK * ZZ0**2)
      Else
        dEpsdZ = 0.0
      End If

    End If

    Tau  = 2.0 * Sig2 / (5.0 * Eps)
    TauU = 6.25 * Tau
    TauV = 4.0  * Tau
    TauW = 1.69 * Tau

    If (Z > ZMin) Then
      dTaudZ = 2.0 * dSig2dZ / (5.0 * Eps) - 2.0 * Sig2 * dEpsdZ / (5.0 * Eps**2)
    Else
      dTaudZ = 0.0
    End If
    dTauUdZ = 6.25 * dTaudZ
    dTauVdZ = 4.0  * dTaudZ
    dTauWdZ = 1.69 * dTaudZ

    Sk    = 0.0
    dSkdZ = 0.0

    KX    = SigU2 * TauU
    KY    = SigV2 * TauV
    KZ    = SigW2 * TauW
    dKXdZ = SigU2 * dTauUdZ + dSigU2dZ * TauU
    dKYdZ = SigV2 * dTauVdZ + dSigV2dZ * TauV
    dKZdZ = SigW2 * dTauWdZ + dSigW2dZ * TauW

    ! Add variation due to spatial variations in the canopy (so-called "dispersive stresses").
    ! Note this addition isn't reflected in the vertical gradients which is OK I think for present 
    ! applications of the turbulence $$.
    If (InCanopy) Then
      If (ZL > ZHat) Then
        UC = UAtHC * ExpDecay
      Else
        UC = (UStarG / VK) * Log((Z + Z0G)/Z0G) 
      End If
      SigU2 = SigU2 + 0.5 * UC**2 * LambdaP
      SigV2 = SigV2 + 0.5 * UC**2 * LambdaP
      KX    = KX    + 0.5 * UC * LambdaP * LC
      KY    = KY    + 0.5 * UC * LambdaP * LC
      TauU  = KX / SigU2
      TauV  = KX / SigV2
    End If

    ! Limiting operations. Note that the sigma's and K's are directly related to the
    ! small and large time dispersion respectively. This is why we limit in terms of
    ! sigma's and K's, not sigma's and tau's. The latter can lead to K's which exceed
    ! both the original value and the above-boundary-layer value. Epsilon is limited
    ! separately and may well therefore be inconsistent with sigma's K's and tau's,
    ! but this is not a significant problem.
    ! (i) Limit by ZByZi * above boundary layer sigma's (keeping K's constant).
    If (SigU2 < SigU2HPlus * ZByZi) Then
      SigU2    = SigU2HPlus * ZByZi
      dSigU2dZ = SigU2HPlus / Zi
      TauU     = KX / SigU2
      dTauUdZ  = (dKXdZ - dSigU2dZ * TauU) / SigU2
    End If
    If (SigV2 < SigU2HPlus * ZByZi) Then
      SigV2    = SigU2HPlus * ZByZi
      dSigV2dZ = SigU2HPlus / Zi
      TauV     = KY / SigV2
      dTauVdZ  = (dKYdZ - dSigV2dZ * TauV) / SigV2
    End If
    If (SigW2 < SigW2HPlus * ZByZi) Then
      SigW2    = SigW2HPlus * ZByZi
      dSigW2dZ = SigW2HPlus / Zi
      TauW     = KZ / SigW2
      dTauWdZ  = (dKZdZ - dSigW2dZ * TauW) / SigW2
    End If
    ! (ii) Limit by ZByZi * above boundary layer K's (keeping sigma's constant).
    If (KX < SigU2HPlus * TauUHPlus * ZByZi) Then
      KX      = SigU2HPlus * TauUHPlus * ZByZi
      dKXdZ   = SigU2HPlus * TauUHPlus / Zi
      TauU    = KX / SigU2
      dTauUdZ  = (dKXdZ - dSigU2dZ * TauU) / SigU2
    End If
    If (KY < SigU2HPlus * TauUHPlus * ZByZi) Then
      KY      = SigU2HPlus * TauUHPlus * ZByZi
      dKYdZ   = SigU2HPlus * TauUHPlus / Zi
      TauV    = KY / SigV2
      dTauVdZ  = (dKYdZ - dSigV2dZ * TauV) / SigV2
    End If
    If (KZ < SigW2HPlus * TauWHPlus * ZByZi) Then
      KZ      = SigW2HPlus * TauWHPlus * ZByZi
      dKZdZ   = SigW2HPlus * TauWHPlus / Zi
      TauW    = KZ / SigW2
      dTauWdZ  = (dKZdZ - dSigW2dZ * TauW) / SigW2
    End If
    ! (iii) Limit by above boundary layer epsilon.
    If (TauWHPlus /= 0.0) Then
      Eps = Max(Eps, 2.0 * SigW2HPlus / (5.0 * TauWHPlus))
    End If

  ! Case 3: Unstable conditions.
  Else

    ZByZi = ZL / Zi

    SigU2 = 0.4  * WStar**2 +                                       &
            6.25 * UStar**2 * (1.0 - ZByZi)**1.5
    SigV2 = 0.4  * WStar**2 +                                       &
            4.0  * UStar**2 * (1.0 - ZByZi)**1.5
    SigW2 = 1.2  * WStar**2 * (1.0 - ZByZi) * (ZZ0/Zi)**(2.0/3.0) + &
            1.69 * UStar**2 * (1.0 - ZByZi)**1.5
    If (InCanopy) Then
      ExpDecay = Exp(- (HC - ZL) / LExp)
      SigU2 = SigU2 * Max(ExpDecay, UStarG / UStar)**2
      SigV2 = SigV2 * Max(ExpDecay, UStarG / UStar)**2
      SigW2 = SigW2 * Max(ExpDecay, UStarG / UStar)**2
    End If

    If (Z > ZMin) Then
      dSigU2dZ = - 6.25 * 1.5 * UStar**2 * Sqrt(1.0 - ZByZi) / Zi
      dSigV2dZ = - 4.0  * 1.5 * UStar**2 * Sqrt(1.0 - ZByZi) / Zi
      dSigW2dZ = (                                                        &
                   0.8 * WStar**2 * (1.0 - ZByZi) / (ZZ0/Zi)**(1.0/3.0) - &
                   1.2 * WStar**2 * (ZZ0/Zi)**(2.0/3.0) -                 &
                   1.69 * 1.5 * UStar**2 * Sqrt(1.0 - ZByZi)              &
                 ) / Zi
      If (InCanopy) Then
        If (ExpDecay > UStarG / UStar) Then
          dSigU2dZ = dSigU2dZ * ExpDecay**2 + SigU2 * 2.0 / LExp
          dSigV2dZ = dSigV2dZ * ExpDecay**2 + SigV2 * 2.0 / LExp
          dSigW2dZ = dSigW2dZ * ExpDecay**2 + SigW2 * 2.0 / LExp
        Else
          dSigU2dZ = dSigU2dZ * (UStarG / UStar)**2
          dSigV2dZ = dSigV2dZ * (UStarG / UStar)**2
          dSigW2dZ = dSigW2dZ * (UStarG / UStar)**2
        End If
      End If
    Else
      dSigU2dZ = 0.0
      dSigV2dZ = 0.0
      dSigW2dZ = 0.0
    End If

    If (InCanopy) Then

      UStarNotUStarG = (UStar * ExpDecay)**3 / (HC - D) >= UStarG**3 / (ZL + Z0G)

      If (UStarNotUStarG) Then
        Eps = WStar**3 * (1.5 - 1.2 * (ZZ0/(Zi + Z0))**(1.0/3.0)) / Zi + &
              (UStar * ExpDecay)**3 * (1.0 - ZByZi) / (VK * (HC - D))
      Else
        Eps = WStar**3 * (1.5 - 1.2 * (ZZ0/(Zi + Z0))**(1.0/3.0)) / Zi + &
              UStarG**3 * (1.0 - ZByZi) / (VK * (ZL + Z0G))
      End If

      If (Z > ZMin) Then
        If (UStarNotUStarG) Then
          dEpsdZ = - 0.4 * WStar**3 / ((ZZ0/(Zi + Z0))**(2.0/3.0) * Zi * (Zi + Z0))       &
                   - (UStar * ExpDecay)**3 / (VK * (HC - D) * Zi)                         &
                   + (UStar * ExpDecay)**3 * (1.0 - ZByZi) * 3.0 / (VK * (HC - D) * LExp)
        Else
          dEpsdZ = - 0.4 * WStar**3 / ((ZZ0/(Zi + Z0))**(2.0/3.0) * Zi * (Zi + Z0)) &
                   - UStarG**3 / (VK * (ZL + Z0G) * Zi)                             &
                   - UStarG**3 * (1.0 - ZByZi) / (VK * (ZL + Z0G)**2)
        End If
      Else
        dEpsdZ = 0.0
      End If

    Else

      Eps = WStar**3 * (1.5 - 1.2 * (ZZ0/(Zi + Z0))**(1.0/3.0)) / Zi + &
            UStar**3 * (1.0 - ZByZi) / (VK * ZMD)

      If (Z > ZMin) Then
        dEpsdZ = - 0.4 * WStar**3 / ((ZZ0/(Zi + Z0))**(2.0/3.0) * Zi * (Zi + Z0)) &
                 - UStar**3 / (VK * ZMD * Zi)                                     &
                 - UStar**3 * (1.0 - ZByZi) / (VK * ZMD**2)
      Else
        dEpsdZ = 0.0
      End If

    End If

    C0 = Max(Abs(RecipLMO), 1.0E-6)
    C0 = 2.1 * Log10(3.0 / C0) - 1.5
    C0 = Min(5.0, C0)
    C0 = Max(1.0, C0)

    TauU = 2.0 * SigU2 / (C0 * Eps)
    TauV = 2.0 * SigV2 / (C0 * Eps)
    TauW = 2.0 * SigW2 / (C0 * Eps)

    If (Z > ZMin) Then
      dTauUdZ = 2.0 * dSigU2dZ / (C0 * Eps) - 2.0 * SigU2 * dEpsdZ / (C0 * Eps**2)
      dTauVdZ = 2.0 * dSigV2dZ / (C0 * Eps) - 2.0 * SigV2 * dEpsdZ / (C0 * Eps**2)
      dTauWdZ = 2.0 * dSigW2dZ / (C0 * Eps) - 2.0 * SigW2 * dEpsdZ / (C0 * Eps**2)
    Else
      dTauUdZ = 0.0
      dTauVdZ = 0.0
      dTauWdZ = 0.0
    End If

    Sk    = 0.6
    dSkdZ = 0.0

    KX    = SigU2 * TauU
    KY    = SigV2 * TauV
    KZ    = SigW2 * TauW
    dKXdZ = SigU2 * dTauUdZ + dSigU2dZ * TauU
    dKYdZ = SigV2 * dTauVdZ + dSigV2dZ * TauV
    dKZdZ = SigW2 * dTauWdZ + dSigW2dZ * TauW

    ! Add variation due to spatial variations in the canopy (so-called "dispersive stresses").
    ! Note this addition isn't reflected in the vertical gradients which is OK I think for present 
    ! applications of the turbulence $$.
    If (InCanopy) Then
      If (ZL > ZHat) Then
        UC = UAtHC * ExpDecay
      Else
        UC = (UStarG / VK) * Log((Z + Z0G)/Z0G) 
      End If
      SigU2 = SigU2 + 0.5 * UC**2 * LambdaP
      SigV2 = SigV2 + 0.5 * UC**2 * LambdaP
      KX    = KX    + 0.5 * UC * LambdaP * LC
      KY    = KY    + 0.5 * UC * LambdaP * LC
      TauU  = KX / SigU2
      TauV  = KX / SigV2
    End If

    ! Limiting operations. Note that the sigma's and K's are directly related to the
    ! small and large time dispersion respectively. This is why we limit in terms of
    ! sigma's and K's, not sigma's and tau's. The latter can lead to K's which exceed
    ! both the original value and the above-boundary-layer value. Epsilon is limited
    ! separately and may well therefore be inconsistent with sigma's K's and tau's,
    ! but this is not a significant problem.
    ! (i) Limit by ZByZi * above boundary layer sigma's (keeping K's constant).
    If (SigU2 < SigU2HPlus * ZByZi) Then
      SigU2    = SigU2HPlus * ZByZi
      dSigU2dZ = SigU2HPlus / Zi
      TauU     = KX / SigU2
      dTauUdZ  = (dKXdZ - dSigU2dZ * TauU) / SigU2
    End If
    If (SigV2 < SigU2HPlus * ZByZi) Then
      SigV2    = SigU2HPlus * ZByZi
      dSigV2dZ = SigU2HPlus / Zi
      TauV     = KY / SigV2
      dTauVdZ  = (dKYdZ - dSigV2dZ * TauV) / SigV2
    End If
    If (SigW2 < SigW2HPlus * ZByZi) Then
      SigW2    = SigW2HPlus * ZByZi
      dSigW2dZ = SigW2HPlus / Zi
      TauW     = KZ / SigW2
      dTauWdZ  = (dKZdZ - dSigW2dZ * TauW) / SigW2
    End If
    ! (ii) Limit by ZByZi * above boundary layer K's (keeping sigma's constant).
    If (KX < SigU2HPlus * TauUHPlus * ZByZi) Then
      KX      = SigU2HPlus * TauUHPlus * ZByZi
      dKXdZ   = SigU2HPlus * TauUHPlus / Zi
      TauU    = KX / SigU2
      dTauUdZ  = (dKXdZ - dSigU2dZ * TauU) / SigU2
    End If
    If (KY < SigU2HPlus * TauUHPlus * ZByZi) Then
      KY      = SigU2HPlus * TauUHPlus * ZByZi
      dKYdZ   = SigU2HPlus * TauUHPlus / Zi
      TauV    = KY / SigV2
      dTauVdZ  = (dKYdZ - dSigV2dZ * TauV) / SigV2
    End If
    If (KZ < SigW2HPlus * TauWHPlus * ZByZi) Then
      KZ      = SigW2HPlus * TauWHPlus * ZByZi
      dKZdZ   = SigW2HPlus * TauWHPlus / Zi
      TauW    = KZ / SigW2
      dTauWdZ  = (dKZdZ - dSigW2dZ * TauW) / SigW2
    End If
    ! (iii) Limit by above boundary layer epsilon.
    If (TauWHPlus /= 0.0) Then
      Eps = Max(Eps, 2.0 * SigW2HPlus / (5.0 * TauWHPlus))
    End If

  End If

  ! Impose minimum Lagrangian timescale TauMin.
  If (TauU < TauMin) Then
    TauU    = TauMin
    dTauUdZ = 0.0
    KX      = SigU2 * TauU
    dKXdZ   = SigU2 * dTauUdZ + dSigU2dZ * TauU
  End If
  If (TauV < TauMin) Then
    TauV    = TauMin
    dTauVdZ = 0.0
    KY      = SigV2 * TauV
    dKYdZ   = SigV2 * dTauVdZ + dSigV2dZ * TauV
  End If
  If (TauW < TauMin) Then
    TauW    = TauMin
    dTauWdZ = 0.0
    KZ      = SigW2 * TauW
    dKZdZ   = SigW2 * dTauWdZ + dSigW2dZ * TauW
  End If

  ! U statistics equal to V statistics.
  If (UEqV) Then
    SigU2    = SigV2
    TauU     = TauV
    dSigU2dZ = dSigV2dZ
    dTauUdZ  = dTauVdZ
    KX       = KY
    dKXdZ    = dKYdZ
  End If

End Subroutine InhomogTurb

!-------------------------------------------------------------------------------------------------------------

Subroutine HomogTurb(                                               &
             Z,                                                     &
             UStar, RecipLMO, Zi, WStar,                            &
             SigU2HPlus, SigW2HPlus, TauUHPlus, TauWHPlus,          &
             UEqV, TauMin,                                          &
             SigU2, SigV2, SigW2, TauU, TauV, TauW, KX, KY, KZ, Eps &
           )
! Calculates homogeneous velocity variances, Lagrangian timescales, eddy diffusivities
! and dissipation rate per unit mass. Here X and U are aligned with the mean flow.

  Implicit None
  ! Argument List:
  Real(Std), Intent(In)  :: Z          ! Height above ground.
  Real(Std), Intent(In)  :: UStar      ! Friction velocity.
  Real(Std), Intent(In)  :: RecipLMO   ! 1/Monin-Obukhov length (set to +/-200000 in
                                       ! calms).
  Real(Std), Intent(In)  :: Zi         ! Boundary layer depth.
  Real(Std), Intent(In)  :: WStar      ! Convective velocity scale if HeatFlux > 0,
                                       ! zero if HeatFlux <= 0.
  Real(Std), Intent(In)  :: SigU2HPlus !} Velocity variances above boundary layer (U
  Real(Std), Intent(In)  :: SigW2HPlus !} here means U and V).
  Real(Std), Intent(In)  :: TauUHPlus  !] Lagrangian time scales above boundary layer
  Real(Std), Intent(In)  :: TauWHPlus  !] (U here means U and V).
  Logical,   Intent(In)  :: UEqV       ! Indicates SigU and TauU are set equal to SigV
                                       ! and TauV.
  Real(Std), Intent(In)  :: TauMin     ! Minimum Lagrangian time scale.
  Real(Std), Intent(Out) :: SigU2      !} Velocity variances.
  Real(Std), Intent(Out) :: SigV2      !}
  Real(Std), Intent(Out) :: SigW2      !}
  Real(Std), Intent(Out) :: TauU       !] Lagrangian time scales.
  Real(Std), Intent(Out) :: TauV       !]
  Real(Std), Intent(Out) :: TauW       !]
  Real(Std), Intent(Out) :: KX         !} Eddy diffusivities.
  Real(Std), Intent(Out) :: KY         !}
  Real(Std), Intent(Out) :: KZ         !}
  Real(Std), Intent(Out) :: Eps        ! Dissipation rate per unit mass.
  ! Locals:
  Real(Std) :: C0    ! Constant for determining Lagrangian time scales.
  Real(Std) :: ZiByL ! Zi/L.
  Real(Std) :: ZStar ! Temporary variable.
  Real(Std) :: Temp  ! Temporary variable.

  ! Case 1: Above boundary layer. We also include low values of UStar and WStar here.
  ! For these cases the full calculation will agree with the above boundary layer
  ! values, but may involve overflows or divisions by zero. The criteria for UStar and
  ! WStar being small is somewhat arbitrary, and could be adjusted if necessary.
  If (Z >= Zi .or. WStar + UStar <= 0.0001) Then

    SigU2 = SigU2HPlus
    SigV2 = SigU2HPlus
    SigW2 = SigW2HPlus
    TauU  = TauUHPlus
    TauV  = TauUHPlus
    TauW  = TauWHPlus
    KX    = SigU2 * TauU
    KY    = SigV2 * TauV
    KZ    = SigW2 * TauW
    If (TauWHPlus == 0.0) Then
      Eps    = 0.0
    Else
      Eps    = 2.0 * SigW2HPlus / (5.0 * TauWHPlus)
    End If

  ! Case 2: Stable conditions.
  Else If (RecipLMO >= 0.0) Then

    SigU2 = 6.25 * 0.4 * UStar**2
    SigV2 = 4.0  * 0.4 * UStar**2
    SigW2 = 1.69 * 0.4 * UStar**2

    ZiByL = Zi * RecipLMO
    If (ZiByL > 0.4 * Epsilon(1.0_Std) ** (2.0 / 9.0)) Then
      ! For large Zi/L, evaluate using exact expression.
      ZStar = 2.0 * Sqrt(ZiByL / (1.0 + 4.0 * ZiByL))
      Temp  = 8.0 + 6.0 / ZiByL -                                                    &
              3.0 / (ZStar * ZiByL) * Log((1.0 + ZStar)**2 * 4.0 * ZiByL / ZStar**2)
      Temp  = Temp / RecipLMO
    Else
      ! For small Zi/L, evaluate with error O(Zi * (Zi/L)**2). Exact result may have
      ! precision problems, and will fail for RecipLMO = 0.
      Temp = 32.0        * Zi         / (1.0 + 4.0 * ZiByL)    - &
             (96.0/5.0)  * Zi         / (1.0 + 4.0 * ZiByL)**2 - &
             (384.0/7.0) * Zi * ZiByL / (1.0 + 4.0 * ZiByL)**3
    End If
    TauU = 6.25 * VK * Temp / (120.0 * UStar)
    TauV = 4.0  * VK * Temp / (120.0 * UStar)
    TauW = 1.69 * VK * Temp / (120.0 * UStar)

    KX = SigU2 * TauU
    KY = SigV2 * TauV
    KZ = SigW2 * TauW

    Eps = 2.0 * SigW2 / (5.0 * TauW)

    ! Limiting operations. Note that the sigma's and K's are directly related to the
    ! small and large time dispersion respectively. This is why we limit in terms of
    ! sigma's and K's, not sigma's and tau's. The latter can lead to K's which exceed
    ! both the original value and the above-boundary-layer value. Epsilon is limited
    ! separately and may well therefore be inconsistent with sigma's K's and tau's,
    ! but this is not a significant problem.
    ! (i) Limit by half above boundary layer sigma's (keeping K's constant).
    If (SigU2 < SigU2HPlus * 0.5) Then
      SigU2 = SigU2HPlus * 0.5
      TauU  = KX / SigU2
    End If
    If (SigV2 < SigU2HPlus * 0.5) Then
      SigV2 = SigU2HPlus * 0.5
      TauV  = KY / SigV2
    End If
    If (SigW2 < SigW2HPlus * 0.5) Then
      SigW2 = SigW2HPlus * 0.5
      TauW  = KZ / SigW2
    End If
    ! (ii) Limit by half above boundary layer K's (keeping sigma's constant).
    If (KX < SigU2HPlus * TauUHPlus * 0.5) Then
      KX   = SigU2HPlus * TauUHPlus * 0.5
      TauU = KX / SigU2
    End If
    If (KY < SigU2HPlus * TauUHPlus * 0.5) Then
      KY   = SigU2HPlus * TauUHPlus * 0.5
      TauV = KY / SigV2
    End If
    If (KZ < SigW2HPlus * TauWHPlus * 0.5) Then
      KZ   = SigW2HPlus * TauWHPlus * 0.5
      TauW = KZ / SigW2
    End If
    ! (iii) Limit by above boundary layer epsilon.
    If (TauWHPlus /= 0.0) Then
      Eps = Max(Eps, 2.0 * SigW2HPlus / (5.0 * TauWHPlus))
    End If

  ! Case 3: Unstable conditions.
  Else

    SigU2 = 0.4  * WStar**2 + 6.25 * 0.4 * UStar**2
    SigV2 = 0.4  * WStar**2 + 4.0  * 0.4 * UStar**2
    SigW2 = 0.27 * WStar**2 + 1.69 * 0.4 * UStar**2

    C0 = Max(Min(Abs(RecipLMO), 1.0E6), 1.0E-6)
    C0 = 2.1 * Log10(3.0 / C0) - 1.5
    C0 = Min(5.0, C0)
    C0 = Max(1.0, C0)

    ! Trap for very small UStar to prevent overflow/divide by zero.
    If (UStar <= 0.000001) Then
      TauU = Zi / (C0 * WStar) * (3.125 * Log(5.0) - 3.5)
      TauV = Zi / (C0 * WStar) * (3.125 * Log(5.0) - 3.5)
      TauW = Zi / (C0 * WStar) * (335919.0 / 14336.0 - 114375.0 / 8192.0 * Log(5.0))
    Else
      TauU = Min(                                            &
               Zi / (C0 * WStar) * (3.125 * Log(5.0) - 3.5), &
               6.25 * 8.0 * VK * Zi / (15.0 * C0 * UStar)    &
             )
      TauV = Min(                                            &
               Zi / (C0 * WStar) * (3.125 * Log(5.0) - 3.5), &
               4.0  * 8.0 * VK * Zi / (15.0 * C0 * UStar)    &
             )
      TauW = Min(                                                                       &
               Zi / (C0 * WStar) * (335919.0 / 14336.0 - 114375.0 / 8192.0 * Log(5.0)), &
               1.69 * 8.0 * VK * Zi / (15.0 * C0 * UStar)                               &
             )
    End If

    KX = SigU2 * TauU
    KY = SigV2 * TauV
    KZ = SigW2 * TauW

    Eps = 2.0 * SigW2 / (C0 * TauW)

    ! Limiting operations. Note that the sigma's and K's are directly related to the
    ! small and large time dispersion respectively. This is why we limit in terms of
    ! sigma's and K's, not sigma's and tau's. The latter can lead to K's which exceed
    ! both the original value and the above-boundary-layer value. Epsilon is limited
    ! separately and may well therefore be inconsistent with sigma's K's and tau's,
    ! but this is not a significant problem.
    ! (i) Limit by half above boundary layer sigma's (keeping K's constant).
    If (SigU2 < SigU2HPlus * 0.5) Then
      SigU2 = SigU2HPlus * 0.5
      TauU  = KX / SigU2
    End If
    If (SigV2 < SigU2HPlus * 0.5) Then
      SigV2 = SigU2HPlus * 0.5
      TauV  = KY / SigV2
    End If
    If (SigW2 < SigW2HPlus * 0.5) Then
      SigW2 = SigW2HPlus * 0.5
      TauW  = KZ / SigW2
    End If
    ! (ii) Limit by half above boundary layer K's (keeping sigma's constant).
    If (KX < SigU2HPlus * TauUHPlus * 0.5) Then
      KX   = SigU2HPlus * TauUHPlus * 0.5
      TauU = KX / SigU2
    End If
    If (KY < SigU2HPlus * TauUHPlus * 0.5) Then
      KY   = SigU2HPlus * TauUHPlus * 0.5
      TauV = KY / SigV2
    End If
    If (KZ < SigW2HPlus * TauWHPlus * 0.5) Then
      KZ   = SigW2HPlus * TauWHPlus * 0.5
      TauW = KZ / SigW2
    End If
    ! (iii) Limit by above boundary layer epsilon.
    If (TauWHPlus /= 0.0) Then
      Eps = Max(Eps, 2.0 * SigW2HPlus / (5.0 * TauWHPlus))
    End If

  End If

  ! Impose minimum Lagrangian timescale TauMin.
  If (TauU < TauMin) Then
    TauU = TauMin
    KX   = SigU2 * TauU
  End If
  If (TauV < TauMin) Then
    TauV = TauMin
    KY   = SigV2 * TauV
  End If
  If (TauW < TauMin) Then
    TauW = TauMin
    KZ   = SigW2 * TauW
  End If

  ! U statistics equal to V statistics.
  If (UEqV) Then
    SigU2 = SigV2
    TauU  = TauV
    KX    = KY
  End If

End Subroutine HomogTurb

!-------------------------------------------------------------------------------------------------------------

Subroutine ModifyTurbByTerminalVelocity(WSed, Z, Inhomog, Flow)
! Modifies turbulence values by the sedimentation velocity to account for
! trajectory-crossing effects and inertial effects.

  Implicit None
  ! Argument list:
  Real(Std),    Intent(In)    :: WSed     ! Terminal velocity
  Real(Std),    Intent(In)    :: Z        ! Particle height
  Logical,      Intent(In)    :: Inhomog  ! Inhomogeneous flag
  Type(Flow_),  Intent(InOut) :: Flow     ! Flow
  ! Locals:
  Real(Std)                   :: TauMin   ! Minimum timescale
  Integer                     :: i        ! Loop variable

  TauMin  = 20.0    ! $$ This should be input.

  If (Inhomog) Then

    ! Tau derivatives
    Flow%dTauUUdZ(:) = Flow%dTauUUdZ(:) / Sqrt(1.0 + WSed**2/Flow%SigUU(3)) +          &
                       Flow%TauUU(:) * WSed**2 / (2.0 * Flow%SigUU(3)**2) *            &
                       Flow%dSigUUdX(3) / (1.0 + WSed**2 / Flow%SigUU(3))**1.5

    ! Trajectory-crossing effect modification
    Flow%TauUU(:)    = Flow%TauUU(:) / Sqrt(1.0 + WSed**2 / Flow%SigUU(3))


    ! Sig^2 derivatives
    Flow%dSigUUdX(3) = Flow%dSigUUdX(3) / (1.0 + WSed / (Gravity * Flow%TauUU(3)))     &
                       + WSed * Flow%SigUU(3) / (Gravity * Flow%TauUU(3)**2) *         &
                       Flow%dTauUUdZ(3) / (1.0 + WSed / (Gravity * Flow%TauUU(3)))**2

    ! Inertial effects modification
    Flow%SigUU(:) = Flow%SigUU(:) / (1.0 + WSed / (Gravity * Flow%TauUU(:)))

    ! K
    Flow%K(:)     = Flow%SigUU(:) * Flow%TauUU(:)
    Flow%dKdX(3)  = Flow%SigUU(3) * Flow%dTauUUdZ(3) + Flow%TauUU(3) * Flow%dSigUUdX(3)


    ! Impose minimum Lagrangian timescale TauMin.
    Do i = 1, 3
      If (Flow%TauUU(i) < TauMin) Then
        Flow%TauUU(i)    = TauMin
        Flow%dTauUUdZ(i) = 0.0
        Flow%K(i)        = Flow%SigUU(i) * Flow%TauUU(i)
        Flow%dKdX(i)     = Flow%dSigUUdX(i) * Flow%TauUU(i)
      End If
    End Do

  Else

    ! Trajectory-crossing effect modification
    Flow%HTauUU(:) = Flow%HTauUU(:) / Sqrt(1.0 + WSed**2 / Flow%HSigUU(3))

    ! Inertial effects modification
    Flow%HSigUU(:) = Flow%HSigUU(:) / (1.0 + WSed / (Gravity * Flow%HTauUU(:)))

    ! K
    Flow%HK(:) = Flow%HSigUU(:) * Flow%HTauUU(:)

    ! Adjust sig_w^2 and tau_w for interface scheme
    Flow%TauWTurbLayers(:)  = Max(Flow%TauWTurbLayers(:) / Sqrt(1.0 + WSed**2 / Flow%SigW2TurbLayers(:)), &
                                    TauMin)
    Flow%SigW2TurbLayers(:) = Flow%SigW2TurbLayers(:) / (1.0 + WSed / (Gravity * Flow%TauWTurbLayers(:)))


    ! Impose minimum Lagrangian timescale TauMin.
    Do i = 1, 3
      If (Flow%HTauUU(i) < TauMin) Then
        Flow%HTauUU(i) = TauMin
        Flow%HK(i)     = Flow%HSigUU(i) * Flow%HTauUU(i)
      End If
    End Do

  End If


End Subroutine ModifyTurbByTerminalVelocity

!-------------------------------------------------------------------------------------------------------------

End Module FlowAndFlowProfileModule

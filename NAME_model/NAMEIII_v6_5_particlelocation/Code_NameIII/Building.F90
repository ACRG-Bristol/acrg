! Module:  Building Module

module buildingmodule

  Use ServiceModule
  use FlowAndFlowProfileModule

! $$ replace write statements with calls to message once any issues concerning
!    the stand alone buildings model are sorted. Remove stop statements
! $$ add private/public statements

! *************************************************************************
! derived type defining dimensions of building and building effects region
! *************************************************************************

  type dimensions

! length scales in the x direction

    real :: lu     ! extent of upwind region
    real :: lur    ! extent of upwind recirculation region
    real :: lb     ! length of building
    real :: ldr    ! extent of downwind recirculation region
    real :: lw     ! extent of wake

! length scales in the y direction

    real :: wb     ! width of building
    real :: we     ! width of building effects region

! length scales in the z direction

    real :: hb     ! height of building
    real :: he     ! height of building effects region

  end type dimensions


! ****************************************
! derived type defining shape of building
! ****************************************

  type configuration

    real :: alpha
    real :: beta
    real :: gamma
    real :: zs

  end type configuration


! ******************************************
! derived type defining origin of l(lambda)
! ******************************************

  type origin_

    real :: y  ! y component of origin
    real :: z  ! z component of origin

  end type origin_


! ********************************************
! derived type defining gradient of l(lambda)
! ********************************************

  type gradient_

    real :: a  ! y component of gradient
    real :: b  ! z component of gradient

  end type gradient_


! ***********************************************
! derived type defining distance along l(lambda)
! ***********************************************

  type r_

    real :: o  ! outer boundary
    real :: i  ! inner boundary

  end type r_


! *********************************
! derived type defining turbulence
! *********************************

  type turbulence

    real :: sig2(3)          ! velocity variances
    real :: sig2_old(3)      ! velocity variances
    real :: dsig2_dx(3)      ! derivatives of velocity variances
    real :: tau(3)           ! correlation time scale
    real :: dTauUdX(3)       ! correlation time scale gradients
    real :: dTauVdX          ! this extra component is due to
                             ! not distinguishing between along
                             ! and cross wind directions in UK data.

  end type turbulence


! ************************************************
! derived type defining meteorological parameters
! ************************************************

  type meteorological_parameters

    real :: u_star           ! friction velocity
    real :: h                ! boundary layer depth
    real :: k                ! von Karman constant
    real :: z_0              ! roughness length

  end type meteorological_parameters

! ********************************
! derived type defining constants
! ********************************

  type constants_

    real :: pi               ! pi
    real :: z_tau_critical   ! tau(z)=constant if z < z_tau_critical
    real :: x_0              ! constant in function defining dividing line

  end type constants_

! *********************************
! derived type defining vortices
! *********************************

  type vortex_

    real    :: vorticity(4)            ! vorticity of vortex - conserved
    real    :: Vt_initial              ! maximum tangental velocity @ r_c
                                       ! set at trailing edge of building
    real    :: r_c_initial             ! radius of vortex core set at trailing
    real    :: growth_rate             ! growth of vortex as funct of x
    real    :: vortex_axis(2,4)        ! axis position of vortices
    real    :: building_phi            ! andgle of building
    real    :: building_trailing_edge  ! x position of building trailing edge

  end type vortex_


! ************************************
! derived type for overall flow state
! ************************************

  type flow_state

    type (dimensions)                :: cuboid
    type (configuration)             :: shape
    type (meteorological_parameters) :: met
    type (constants_)                :: constants
    type (vortex_)                   :: vortices
    type (FlowField_)                :: FlowField
    type (Coords_), Pointer          :: Coords
    type (Grids_),  Pointer          :: Grids

  end type flow_state


  contains

!**********************************************************************
!**********************************************************************
! Subroutines
!**********************************************************************
!**********************************************************************

! *******************************
! subroutine building_parameters
! *******************************

  subroutine building_parameters(flow_state_building)

  implicit none

! argument list:
  type (flow_state), intent (inout) :: flow_state_building

! locals:
  real :: tol

  tol = flow_state_building%cuboid%hb * epsilon(flow_state_building%cuboid%hb)


! upwind height of dividing line = alpha * Hb
  flow_state_building%shape%alpha = 1.5

! gamma is related to origin of radial fanning region.

  if ( abs(flow_state_building%cuboid%hb - flow_state_building%cuboid%wb) <= tol ) then
! cube
    flow_state_building%shape%gamma = 0.5
  elseif ( flow_state_building%cuboid%wb > flow_state_building%cuboid%hb ) then
! wide building
    flow_state_building%shape%gamma = 1.0      ! 0.5
    write(11,*)'Wide building; end regions (gamma*H_b) = ',                       &
                flow_state_building%shape%gamma * flow_state_building%cuboid%hb
  else
! tall building
!
! under this tall building formulation gamma = 2 for a cube. gamma = 4 is used as a limit for
! buildings where hb >= 3wb. Between the cube and this limit gamma is varied linearly.

    flow_state_building%shape%gamma = min(4.0, (2.0 +                                &
               (flow_state_building%cuboid%hb/flow_state_building%cuboid%wb - 1.0)))

    flow_state_building%shape%alpha = (flow_state_building%cuboid%hb           &
              + (flow_state_building%cuboid%hb/flow_state_building%shape%gamma)) &
                /flow_state_building%cuboid%hb
  endif

! upwind height of dividing line = alpha * Hb

 !   flow_state_building%shape%alpha = 1.5


! ADMS formula for dividing line evaluated at beta * Hb

  flow_state_building%shape%beta = 0.5

! downstream 'starting point' of vorticies

  flow_state_building%vortices%building_trailing_edge =                 &
                               flow_state_building%cuboid%lb/2.0

  flow_state_building%shape%zs = 0.0

! dimensions of building effects region

! upwind length

  flow_state_building%cuboid%lu = flow_state_building%cuboid%hb

! upwind recirculation region length

  flow_state_building%cuboid%lur = flow_state_building%cuboid%hb/4.0

! recirculation region length
! This uses Fackrell relationship. Ths relationship has limits of applicability
! 0.3 <= Lb/Hb <= 3.0. Out side of these limits we apply the end limit value.

  if ( flow_state_building%cuboid%lb/flow_state_building%cuboid%hb < 0.3 ) then

    flow_state_building%cuboid%ldr = (1.8 * flow_state_building%cuboid%wb)   &
                                    / (((0.3)**0.3) * (1.0 + (0.24         &
                                       * flow_state_building%cuboid%wb       &
                                        / flow_state_building%cuboid%hb)))

  elseif ( flow_state_building%cuboid%lb/flow_state_building%cuboid%hb > 3.0 ) then

    flow_state_building%cuboid%ldr = (1.8 * flow_state_building%cuboid%wb)   &
                                   / (((3.0)**0.3) * (1.0 + (0.24          &
                                      * flow_state_building%cuboid%wb        &
                                       /flow_state_building%cuboid%hb)))

  else

    flow_state_building%cuboid%ldr = (1.8 * flow_state_building%cuboid%wb)   &
                      / (((flow_state_building%cuboid%lb                   &
                         / flow_state_building%cuboid%hb)**0.3)              &
                      * (1.0 + (0.24 * flow_state_building%cuboid%wb     &
                           / flow_state_building%cuboid%hb)))

  end if

! height and width of building effects region

  if ( abs(flow_state_building%cuboid%hb - flow_state_building%cuboid%wb) <= tol ) then
! cube
    flow_state_building%cuboid%he = 5.0 * flow_state_building%cuboid%hb
    flow_state_building%cuboid%we = 6.0 * flow_state_building%cuboid%wb

  elseif ( flow_state_building%cuboid%wb > flow_state_building%cuboid%hb ) then
! wide building
    flow_state_building%cuboid%he = 5.0 * flow_state_building%cuboid%hb
!    flow_state_building%cuboid%we = flow_state_building%cuboid%wb +      &
!                                    4.0 * (2.0 * flow_state_building%shape%gamma) &
!                                        * flow_state_building%cuboid%wb
    flow_state_building%cuboid%we = flow_state_building%cuboid%wb +      &
                                    4.0 * (2.0 * flow_state_building%shape%gamma) &
                                        * flow_state_building%cuboid%hb

  else
! tall building
    flow_state_building%cuboid%he = flow_state_building%cuboid%hb +      &
                                    4.0 * (flow_state_building%cuboid%wb/2.0) &
                                        * flow_state_building%shape%gamma
!    flow_state_building%cuboid%we = 5.0 * flow_state_building%cuboid%wb   !  *12.0
    flow_state_building%cuboid%we = 4.0 * flow_state_building%cuboid%hb

  endif


! length of wake

  flow_state_building%cuboid%lw = 16.0 * flow_state_building%cuboid%hb


  end subroutine building_parameters


! *******************************
! subroutine meteorological_data
! *******************************

  subroutine meteorological_data (flow_state_building)

  implicit none

! argument list:
  type (flow_state), intent (out) :: flow_state_building

  Allocate(flow_state_building%FlowField%ProfileData(1,1,1,1))

! wind direction

  flow_state_building%FlowField%ProfileData(1,1,1,1)%Phi0 = 0.0

! friction velocity

  flow_state_building%FlowField%ProfileData(1,1,1,1)%UStar = 0.2
  flow_state_building%met%u_star = flow_state_building%FlowField%ProfileData(1,1,1,1)%UStar

! heat flux

  flow_state_building%FlowField%ProfileData(1,1,1,1)%WT = 0.0

! boundary layer depth

  flow_state_building%FlowField%ProfileData(1,1,1,1)%H = 400.0
  flow_state_building%met%h = flow_state_building%FlowField%ProfileData(1,1,1,1)%H

! roughness length

  flow_state_building%FlowField%ProfileData(1,1,1,1)%Z0 = 0.4
  flow_state_building%met%z_0 = flow_state_building%FlowField%ProfileData(1,1,1,1)%Z0

! Temperature

  flow_state_building%FlowField%ProfileData(1,1,1,1)%T0 = 288.0

! Height of temperature above ground

  flow_state_building%FlowField%ProfileData(1,1,1,1)%ZT = 1.2

! Reciprocal Monin-Obukhov length

  flow_state_building%FlowField%ProfileData(1,1,1,1)%RecipLMO = 0.0

! von Karman's constant

  flow_state_building%met%k = VK

  flow_state_building%FlowField%ProfileData(1,1,1,1)%MeanFlow = .true.
  flow_state_building%FlowField%ProfileData(1,1,1,1)%Turb     = .true.

  end subroutine meteorological_data


! *******************************
! subroutine calculate constants
! *******************************

  subroutine calculate_constants (flow_state_building, velocity_profile)

  implicit none

! argument list:
  type (flow_state), intent (inout) :: flow_state_building
  character(*),      intent (in)    :: velocity_profile

! locals:
  integer            :: i
  integer, parameter :: imax=10000
  real               :: a
  real               :: b
  real               :: c
  real               :: p
  real               :: q
  real               :: j0
  real               :: jp
  real               :: tol

! tau(z)=constant if z < z_tau_critical

  flow_state_building%constants%z_tau_critical = 1.0

! calculate x_0 using bisection method

  a = 0.5 * flow_state_building%shape%beta * flow_state_building%shape%beta
  b = 40.0 * flow_state_building%shape%beta * flow_state_building%shape%beta

  tol = 10.0 * epsilon(a)

  do i = 0, imax

    p = 0.5 * (a + b)
    q = 0.5 * (b - a)

    jp = d(p, flow_state_building)

    if ( jp /= 0.0 .and. q >= tol ) then

      j0 = d(a, flow_state_building)
      if ( j0 * jp > 0.0 ) then
        a = p
      else
        b = p
      endif

    elseif ( i /= imax ) then

      c = flow_state_building%cuboid%hb *                                        &
          u_mean(flow_state_building%shape%beta * flow_state_building%cuboid%hb, &
                 flow_state_building, velocity_profile)/                         &
          ( 8.0 * flow_state_building%met%k * flow_state_building%met%u_star )

      flow_state_building%constants%x_0 = 0.5 * flow_state_building%cuboid%lb +  &
                                          flow_state_building%cuboid%ldr - c * p

      exit

    else

      write (6,*) 'Error: bisection method failed'
      exit

    endif

  end do

  end subroutine calculate_constants


! **************
! function D(y)
! **************

  real function d(y, flow_state_building)

  implicit none

! argument list:
  real,              intent (in) :: y
  type (flow_state), intent (in) :: flow_state_building

! locals:
  real :: a1
  real :: a2
  real :: a

  if ( flow_state_building%cuboid%wb >= flow_state_building%cuboid%hb ) then

! wide building

    a1 = flow_state_building%shape%alpha * flow_state_building%shape%gamma *   &
         flow_state_building%cuboid%hb + 0.5 * flow_state_building%cuboid%wb - &
         flow_state_building%shape%gamma * flow_state_building%cuboid%hb

    a2 = sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) * &
         flow_state_building%shape%gamma * flow_state_building%cuboid%hb +               &
         0.5 * flow_state_building%cuboid%wb - flow_state_building%shape%gamma *         &
         flow_state_building%cuboid%hb

  else

! tall building

    a1 = flow_state_building%shape%alpha * flow_state_building%shape%gamma *   &
         0.5 * flow_state_building%cuboid%wb + flow_state_building%cuboid%hb - &
         flow_state_building%shape%gamma * 0.5 * flow_state_building%cuboid%wb

    a2 = sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) * &
         flow_state_building%shape%gamma * 0.5 * flow_state_building%cuboid%wb +         &
         flow_state_building%cuboid%hb - flow_state_building%shape%gamma * 0.5 *         &
         flow_state_building%cuboid%wb

  endif

  a = flow_state_building%shape%alpha * a1/ &
        ( sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) * a2 )

  d = y * y * exp( flow_state_building%shape%beta * flow_state_building%shape%beta/y ) -   &
      flow_state_building%shape%beta * 2.0**3.5 * flow_state_building%cuboid%wb/           &
      (5.0 * flow_state_building%constants%pi * (1.0 - a) * flow_state_building%cuboid%hb)

  end function d


! **********************************************************


! *******************************
! subroutine turbulence_building
! *******************************

  subroutine turbulence_building(x, flow_state_building, building_turbulence, velocity_profile, enhancment)

  implicit none

! argument list:
!  real,              intent (in)  :: x(:)
  real :: x(:)
  type (flow_state), intent (in)  :: flow_state_building
  type (turbulence), intent (out) :: building_turbulence
  character(*),      intent (in)  :: velocity_profile

  real, intent (out)          :: enhancment(3)

! locals:
  character(14) :: turbulence_enhancement
  real          :: u
  real          :: gamma
  real          :: d_gamma(3)
  real          :: dist_edge
  real          :: z

  real          :: dDdX(3)
!  real          :: dTauUdX(3)
  Type(ShortTime_) :: Time
  Type(Flow_)      :: Flow
  Type(HCoeffs_)   :: HCoeffs
  Type(ZCoeffs_)   :: ZCoeffs
  Type(TCoeffs_)   :: TCoeffs

  HCoeffs%Valid = .false. ! Could save for more efficient work with NWP data? $$
  ZCoeffs%Valid = .false.
  TCoeffs%Valid = .false.

!integer         :: i
!i = 1
!do i = 1,400
!x(1) = -1300.0
!x(2) = 0.0
!x(3) = real(i)

! covariance matrix and its derivatives

! gamma

! gamma and d_gamma set = 0 turns enhanced turbulence
! in building cavity off.
  enhancment = 0.0
  turbulence_enhancement = 'enhanced'

  if ( turbulence_enhancement == 'enhanced' ) then

    if ( x(1) >= 0.5 * flow_state_building%cuboid%lb .and. &
         x(1) <= (0.5 * flow_state_building%cuboid%lb +    &
          flow_state_building%cuboid%ldr) ) then

      gamma=max(gamma_prime(x,flow_state_building%cuboid, &
                turbulence_enhancement),0.0)

      enhancment(1) = 0.0
      enhancment(2) = 0.1
      enhancment(3) = 0.15

    else
      gamma=0.0
    endif

    if ( gamma > 0.0 ) then
      d_gamma=d_gamma_prime(x,flow_state_building%cuboid, &
                            turbulence_enhancement)
    else
      d_gamma=0.0
    endif

  elseif ( turbulence_enhancement == 'extra_enhanced' ) then
    gamma=max(gamma_prime(x,flow_state_building%cuboid, &
                          turbulence_enhancement),0.0)
    if ( gamma > 0.0 ) then
      d_gamma=d_gamma_prime(x,flow_state_building%cuboid, &
                            turbulence_enhancement)
    else
      d_gamma=0.0
    endif

  else
    gamma=0.0
    d_gamma=0.0
  endif

! mean velocity at 0.5 H_B

  u = u_mean(0.5 * flow_state_building%cuboid%hb,   &
             flow_state_building, velocity_profile)

! upwind velocity variances and derivatives

  call TurbFromFlowField(                                               &
         flow_state_building%Coords, flow_state_building%Grids,         &
         Time, (/0.0, 0.0, x(3)/), x(3), flow_state_building%FlowField, &
         .true., .false., .false., 0.0, 0.0,                            &
         Flow,                                                          &
         HCoeffs, ZCoeffs, TCoeffs                                      &
       )

! modified velocity variances


  building_turbulence%sig2(1) = Flow%SigUU(1) + gamma * 0.25 * u**2 * enhancment(1)

  building_turbulence%sig2(2) = Flow%SigUU(2) + gamma * 0.25 * u**2 * enhancment(2)

  building_turbulence%sig2(3) = Flow%SigUU(3) + gamma * 0.25 * u**2 * enhancment(3)

! modified derivatives

  building_turbulence%dsig2_dx(1) = d_gamma(1) * 0.25 * u**2 * enhancment(1)

  building_turbulence%dsig2_dx(2) = d_gamma(2) * 0.25 * u**2 * enhancment(2)

  building_turbulence%dsig2_dx(3) = Flow%dSigUUdX(3) + &
                                    d_gamma(3) * 0.25 * u**2 * enhancment(3)

!write(29,*)x(3),SigUU

!enddo
!stop

! correlation time scale tau

! time scale corresponding to vertical velocity component

    dDdX = 0.0

! directly above building

  if ( (abs(x(1)) <= flow_state_building%cuboid%lb*0.5) .and. &
       (abs(x(2)) <= flow_state_building%cuboid%wb*0.5) ) then

    z = max( (x(3)-flow_state_building%cuboid%hb), &
             flow_state_building%constants%z_tau_critical)

    dDdX(3) = 1.0

! down/upwind; below roof level; across building width

  elseif ( (abs(x(1)) > flow_state_building%cuboid%lb*0.5)  .and. &
           (abs(x(2)) <= flow_state_building%cuboid%wb*0.5) .and. &
           (x(3) <= flow_state_building%cuboid%hb) ) then

    if ( (abs(x(1)) - flow_state_building%cuboid%lb*0.5) < x(3) ) then

      z = max( abs(x(1)) - flow_state_building%cuboid%lb*0.5, &
               flow_state_building%constants%z_tau_critical)

      dDdX(1) = SIGN(1.0,x(1))

    else

      z = max( x(3), flow_state_building%constants%z_tau_critical)

      dDdX(3) = 1.0

    endif

! down/upwind; above roof level; across building width

  elseif ( (abs(x(1)) > flow_state_building%cuboid%lb*0.5)  .and. &
           (abs(x(2)) <= flow_state_building%cuboid%wb*0.5) .and. &
           (x(3) > flow_state_building%cuboid%hb) ) then

     dist_edge = SQRT((abs(x(1)) - flow_state_building%cuboid%lb*0.5)**2 &
                 + (x(3) - flow_state_building%cuboid%hb)**2)


     if ( dist_edge < x(3) ) then

       z = max(dist_edge, flow_state_building%constants%z_tau_critical)

       dDdX(1) = SIGN((abs(x(1)) - flow_state_building%cuboid%lb*0.5)/z,x(1))
       dDdX(3) = (x(3) - flow_state_building%cuboid%hb)/z

     else      !   if ( dist_edge >= x(3) ) then

       z = max(x(3), flow_state_building%constants%z_tau_critical)

       dDdX(3) = 1.0

     endif

! to side of building; below roof level; down/upwind of building

  elseif ( (abs(x(1)) > flow_state_building%cuboid%lb*0.5) .and. &
           (abs(x(2)) > flow_state_building%cuboid%wb*0.5) .and. &
           (x(3) <= flow_state_building%cuboid%hb) ) then

     dist_edge = SQRT((abs(x(2)) - flow_state_building%cuboid%wb*0.5)**2 &
                 + (abs(x(1)) - flow_state_building%cuboid%lb*0.5)**2)

     if ( dist_edge < x(3) ) then

       z = max(dist_edge, flow_state_building%constants%z_tau_critical)

       dDdX(1) = SIGN((abs(x(1)) - flow_state_building%cuboid%lb*0.5)/z,x(1))
       dDdX(2) = SIGN((abs(x(2)) - flow_state_building%cuboid%wb*0.5)/z,x(2))

     else      !   if ( dist_edge >= x(3) ) then

       z = max(x(3), flow_state_building%constants%z_tau_critical)

       dDdX(3) = 1.0

     endif

! to side of building; above roof level; down/upwind of building

  elseif ( (abs(x(1)) > flow_state_building%cuboid%lb*0.5) .and. &
           (abs(x(2)) > flow_state_building%cuboid%wb*0.5) .and. &
           (x(3) > flow_state_building%cuboid%hb) ) then

! this is distance to corner in this case
     dist_edge = SQRT( (abs(x(2)) - flow_state_building%cuboid%wb*0.5)**2 &
                   + (abs(x(1)) - flow_state_building%cuboid%lb*0.5)**2 &
                     + (x(3) - flow_state_building%cuboid%hb)**2)

     if ( dist_edge < x(3) ) then

       z = max(dist_edge, flow_state_building%constants%z_tau_critical)

       dDdX(1) = SIGN((abs(x(1)) - flow_state_building%cuboid%lb*0.5)/z,x(1))
       dDdX(2) = SIGN((abs(x(2)) - flow_state_building%cuboid%wb*0.5)/z,x(2))
       dDdX(3) = (x(3) - flow_state_building%cuboid%hb)/z

     else      !    if ( dist_edge >= x(3) ) then

       z = max(x(3), flow_state_building%constants%z_tau_critical)

       dDdX(3) = 1.0

     endif

! to side of building; below roof level; over length of building

  elseif ( (abs(x(1)) <= flow_state_building%cuboid%lb*0.5) .and. &
           (abs(x(2)) > flow_state_building%cuboid%wb*0.5)  .and. &
           (x(3) <= flow_state_building%cuboid%hb) ) then

     if ( abs(x(2)) - flow_state_building%cuboid%wb*0.5 < x(3) ) then

       z = max( abs(x(2)) - flow_state_building%cuboid%wb*0.5, &
                flow_state_building%constants%z_tau_critical)

       dDdX(2) = SIGN((abs(x(2)) - flow_state_building%cuboid%wb*0.5)/z,x(2))

     else

       z = max(x(3), flow_state_building%constants%z_tau_critical)

       dDdX(3) = 1.0

     endif

! to side of building; above roof level; over length of building

  elseif ( (abs(x(1)) <= flow_state_building%cuboid%lb*0.5) .and. &
           (abs(x(2)) > flow_state_building%cuboid%wb*0.5)  .and. &
           (x(3) > flow_state_building%cuboid%hb) ) then

     dist_edge = SQRT((abs(x(2)) - flow_state_building%cuboid%wb*0.5)**2 &
                 + (x(3) - flow_state_building%cuboid%hb)**2)

     if ( dist_edge < x(3) ) then

       z = max(dist_edge, flow_state_building%constants%z_tau_critical)

       dDdX(2) = SIGN((abs(x(2)) - flow_state_building%cuboid%wb*0.5)/z,x(2))
       dDdX(3) = (x(3) - flow_state_building%cuboid%hb)/z

     else

       z = max(x(3), flow_state_building%constants%z_tau_critical)

       dDdX(3) = 1.0

     endif

! everywhere else - N.B.there is of course no where else,
! if this happens then there is an error

  else

    write(6,*)'Error: Particle not within building effects domain',x

    stop

  endif

! chack if length scale limit has been used and if so then set
! gradients to zero.

  if (z < flow_state_building%constants%z_tau_critical) dDdX = 0.0


  call TurbFromFlowField(                                         &
         flow_state_building%Coords, flow_state_building%Grids,   &
         Time, (/0.0, 0.0, z/), z, flow_state_building%FlowField, &
         .true., .false., .false., 0.0, 0.0,                      &
         Flow,                                                    &
         HCoeffs, ZCoeffs, TCoeffs                                &
       )

! time scale corresponding to alongwind velocity component

  building_turbulence%tau(1) = Flow%TauUU(1)

! time scale corresponding to lateral velocity component

  building_turbulence%tau(2) = Flow%TauUU(2)

! time scale corresponding to vertical velocity component

  building_turbulence%tau(3) = Flow%TauUU(3)

! time scale gradients

  building_turbulence%dTauUdX = Flow%dTauUUdZ * dDdX

  building_turbulence%dTauVdX = Flow%dTauUUdZ(2) * dDdX(1)  ! this extra component is due to
                                                       ! not distinguishing between along
                                                       ! and cross wind directions in UK data.

  end subroutine turbulence_building


! ******************
! coefficient gamma
! ******************

  real function gamma_prime(x, cuboid, turbulence_enhancement)

  implicit none

! argument list:
  real,              intent(in) :: x(:)
  type (dimensions), intent(in) :: cuboid
  character(*),      intent(in) :: turbulence_enhancement

  if ( turbulence_enhancement == 'enhanced') then

    gamma_prime = 1.0 - ( x(3) * x(3)/(cuboid%hb * cuboid%hb) +              &
                        x(2) * x(2)/(0.25 * cuboid%wb * cuboid%wb) +         &
                        (x(1) - 0.5 * cuboid%lb) * (x(1) - 0.5 * cuboid%lb)/ &
                        (cuboid%ldr * cuboid%ldr) )

  else

    gamma_prime = ( 3.0 - ( x(3) * x(3)/(cuboid%hb * cuboid%hb)+             &
                          x(2) * x(2)/(0.25 * cuboid%wb * cuboid%wb) +       &
                          (x(1) - 0.5 * cuboid%lb) * (x(1) - 0.5*cuboid%lb)/ &
                          (cuboid%ldr * cuboid%ldr)))/3.0

  endif

  end function gamma_prime


! *********************
! derivatives of gamma
! *********************

  function d_gamma_prime (x, cuboid, turbulence_enhancement)

  implicit none

! argument list:
  real,              intent(in) :: x(3)
  type (dimensions), intent(in) :: cuboid
  character(*),      intent(in) :: turbulence_enhancement

! locals:
  real, dimension(3) :: d_gamma_prime

  if ( turbulence_enhancement == 'enhanced') then

    d_gamma_prime(1) = -2.0 * (x(1) - 0.5*cuboid%lb)/ &
                       (cuboid%ldr*cuboid%ldr)
    d_gamma_prime(2) = -2.0 * x(2)/(0.25 * cuboid%wb * cuboid%wb)
    d_gamma_prime(3) = -2.0 * x(3)/(cuboid%hb * cuboid%hb)

  else

    d_gamma_prime(1) = -2.0 * (x(1) - 0.5*cuboid%lb)/ &
                       (3.0 * cuboid%ldr * cuboid%ldr)
    d_gamma_prime(2) = -2.0 * x(2)/(0.75 * cuboid%wb * cuboid%wb)
    d_gamma_prime(3) = -2.0 * x(3)/(3.0 * cuboid%hb * cuboid%hb)

  endif

  end function d_gamma_prime


! ************************************
! x component of upwind mean velocity
! ************************************

  real function u_mean(z, flow_state_building, velocity_profile)

  implicit none

! argument list:
  real,              intent(in) :: z
  type (flow_state), intent(in) :: flow_state_building
  character(*),      intent(in) :: velocity_profile

! locals:
  real :: Speed
  real :: U(3)
  real :: T
  real :: Theta
  Type(ShortTime_) :: Time    ! Value irrelevant as buildings asks for single time
                              ! data only.
  Type(HCoeffs_)   :: HCoeffs
  Type(ZCoeffs_)   :: ZCoeffs
  Type(TCoeffs_)   :: TCoeffs
  Type(Flow_)      :: Flow

  HCoeffs%Valid = .false. ! Could save for more efficient work with NWP data? $$
  ZCoeffs%Valid = .false.
  TCoeffs%Valid = .false.

  if ( velocity_profile == 'constant' ) then

! constant velocity

    u_mean = 1.0

  elseif ( velocity_profile == 'linear' ) then

! linear velocity

    u_mean = flow_state_building%met%u_star * &
             (z + flow_state_building%met%z_0)/flow_state_building%met%z_0

  else

! velocity profile according to NAME-PPM

    call MeanFlowFromFlowField(                                   &
           flow_state_building%Coords, flow_state_building%Grids, &
           Time, (/ 0.0, 0.0, z/), z,                             &
           flow_state_building%FlowField,                         &
           .false.,                                               &
           Speed, Flow,                                           &
           HCoeffs, ZCoeffs, TCoeffs                              &
         )
    u_mean = Speed

  endif

  end function u_mean


! *********************
! uniform distribution
! *********************

  real function uniform(a, b)

  implicit none

! argument list:
  real, intent(in) :: a,b

! locals:
  real :: temp

#ifdef SABuilding
  call random_number(temp)
#endif
#ifndef SABuilding
  call GetRandomNumber(temp, 0) ! $$ in fact uniform is only used in the stand-alone building code
#endif

  uniform = (b - a) * temp + a

  end function uniform


!  **************************************************************
!  **************************************************************


! *****************************
! inverse of covariance matrix
! *****************************

  real function v_inverse(i, j, building_turbulence)

  implicit none

! argument list:
  integer,           intent(in) :: i,j
  type (turbulence), intent(in) :: building_turbulence

  if (i == j) then
    v_inverse = 1.0/building_turbulence%sig2(i)
  else
    v_inverse = 0.0
  endif

  end function v_inverse


! **************
! dV(i,j)/dx(j)
! **************

  real function dv_dx(i, j, building_turbulence)

  implicit none

! argument list:
  integer,           intent(in) :: i,j
  type (turbulence), intent(in) :: building_turbulence

  if (i == j) then
    dv_dx = building_turbulence%dsig2_dx(i)
  else
    dv_dx = 0.0
  endif

  end function dv_dx


! ********
! dV(i,j)
! ********

  real function dv(i, j, building_turbulence)

  implicit none

! argument list:
  integer,           intent(in) :: i,j
  type (turbulence), intent(in) :: building_turbulence

  if (i == j) then
    dv = building_turbulence%sig2(i) - building_turbulence%sig2_old(i)
  else
    dv = 0.0
  endif

  end function dv


! ************
! b^ij(x,u,t)
! ************

  real function b(i, building_turbulence)

  implicit none

! argument list:
  integer,           intent(in) :: i
  type (turbulence), intent(in) :: building_turbulence

  b = sqrt(2.0 * b_capital(i, i, building_turbulence))

  end function b


! ************
! B^ij(x,u,t)
! ************

  real function b_capital(i, j, building_turbulence)

  implicit none

! argument list:
  integer,           intent(in) :: i,j
  type (turbulence), intent(in) :: building_turbulence

  if (i == j) then
    b_capital = building_turbulence%sig2(i)/building_turbulence%tau(i)
  else
    b_capital = 0.0
  endif

  end function b_capital


! ***************************************************************
! ***************************************************************


! ****************
! subroutine mean
! ****************

  subroutine mean(x_old, mean_velocity, dt, flow_state_building, velocity_profile)

  implicit none

! argument list:
  real,              intent (in)  :: x_old(:)
  real,              intent (out) :: mean_velocity(:)
  real,              intent (in)  :: dt
  type (flow_state), intent (in)  :: flow_state_building
  character(*),      intent (in)  :: velocity_profile

! locals:
  type (r_)         :: r
  type (r_)         :: r_hat
  type (r_)         :: r_new
  type (origin_)    :: origin
  type (gradient_)  :: gradient
  character (len=5) :: region
  integer           :: sign_y
  real              :: x_hat(3)
  real              :: x(3)
  real              :: jacob
  real              :: r_p
  real              :: tol

  x = x_old

  x_hat = 0.0

! sign of y coordinate

  if ( x(2) < 0.0 ) then
    sign_y = -1
  else
    sign_y = 1
  endif

  x(2) = abs(x(2))

! check: is particle in recirculation region or building?

  if (x(1) >= -0.5*flow_state_building%cuboid%lb - flow_state_building%cuboid%lur .and. &
      x(1) <= 0.5*flow_state_building%cuboid%lb + flow_state_building%cuboid%ldr  .and. &
      x(3) <= xi_r_z(flow_state_building, x(1))                                   .and. &
      x(2) <= xi_r_y(flow_state_building, x(1))) then

    mean_velocity = 0.0

  else

! check: is particle in inner or outer region?

     if ( x(3) <= xi_d_z(flow_state_building, x(1), velocity_profile)  .and. &
          x(2) <= xi_d_y(flow_state_building, x(1), velocity_profile) ) then
       region = 'inner'
     else
       region = 'outer'
     endif

! calculate origin for l(lambda)

     call lambda (flow_state_building, x, origin, gradient)

     tol = flow_state_building%cuboid%hb * epsilon(x(1))

     if ( gradient%a <= tol .and. gradient%b <= tol ) then
       x(2) = real(sign_y) * x(2)
     else

! calculate R_O, R_I

       call boundary (flow_state_building, region, velocity_profile, x(1), origin, gradient, r)

! calculate upwind values of R_O, R_I

       x_hat(1) = -(flow_state_building%cuboid%lu + 0.5 * flow_state_building%cuboid%lb)

       call boundary (flow_state_building, region, velocity_profile, x_hat(1), origin, &
                      gradient, r_hat)

! calculate upwind P_y, P_z

       r_p = sqrt( gradient%a * gradient%a + gradient%b * gradient%b)

       call position (origin, gradient, r_p, r, r_hat, x_hat)

! calculate jacobian

       jacob = ( r_hat%o * r_hat%o - r_hat%i * r_hat%i )/( r%o *r%o - r%i * r%i )

! check: is Jacobian less than zero?

       if ( jacob < 0.0 ) then
         write (6,*) 'Error: Jacobian less than zero',jacob
         stop
       endif

! calculate x'

       x(1) = x(1) + u_mean(x_hat(3), flow_state_building, velocity_profile) * dt * jacob

! calculate R'

       call boundary (flow_state_building, region, velocity_profile, x(1), origin, &
                      gradient, r_new)

! calculate P'_y, P'_z

       call position (origin, gradient, r_p, r, r_new, x)

       x(2) = real(sign_y) * x(2)

     endif

! calculate components of mean velocity

     mean_velocity = (x - x_old)/dt

  endif

  end subroutine mean


! ******************
! subroutine lambda
! ******************

  subroutine lambda (flow_state_building, x, origin, gradient)

  implicit none

! argument list:
  type (flow_state), intent (in)  :: flow_state_building
  real,              intent (in)  :: x(:)
  type (origin_),    intent (out) :: origin
  type (gradient_),  intent (out) :: gradient

! locals:
  real           :: tol
  type (origin_) :: origin_prime

! origin of l(lambda)

  tol = 100.0 * flow_state_building%cuboid%hb * epsilon(x(1))

  if ( abs(flow_state_building%cuboid%hb - flow_state_building%cuboid%wb) <= tol ) then

! cube shaped building

    origin%y = 0.0
    origin%z = 0.0

  elseif ( flow_state_building%cuboid%hb > flow_state_building%cuboid%wb ) then

! tall and thin building

    if ( x(3) > flow_state_building%cuboid%hb -                                                 &
                0.5 * flow_state_building%shape%gamma * flow_state_building%cuboid%wb    .and.  &
         x(2) <= ( x(3) - flow_state_building%cuboid%hb )/flow_state_building%shape%gamma +     &
                 0.5 * flow_state_building%cuboid%wb ) then

      origin%y = 0.0
      origin%z = flow_state_building%cuboid%hb - &
                 0.5 * flow_state_building%shape%gamma * flow_state_building%cuboid%wb

    else

      call fanning_point(flow_state_building%cuboid%hb, 0.5 * flow_state_building%cuboid%wb, &
                         flow_state_building%shape%gamma, (/x(1), x(3), x(2)/), origin_prime)
      origin%y = origin_prime%z
      origin%z = origin_prime%y

    endif

  else

! short and wide building

    if ( x(2) > 0.5 * flow_state_building%cuboid%wb -                                             &
                flow_state_building%shape%gamma * flow_state_building%cuboid%hb             .and. &
         x(3) <= ( x(2) - 0.5 * flow_state_building%cuboid%wb )/flow_state_building%shape%gamma + &
                 flow_state_building%cuboid%hb ) then

      origin%y = 0.5 * flow_state_building%cuboid%wb - &
                 flow_state_building%shape%gamma * flow_state_building%cuboid%hb
      origin%z = 0.0

    else

      call fanning_point(0.5 * flow_state_building%cuboid%wb, flow_state_building%cuboid%hb, &
                         flow_state_building%shape%gamma, x, origin)

    endif

  endif

! gradient of l(lambda)

  gradient%a = x(2) - origin%y
  gradient%b = x(3) - origin%z

  if ( gradient%a < 0.0 .and. gradient%a >= -tol ) then
!    write (*,*) 'Warning: gradient of l(lambda) less than zero ',gradient
    gradient%a = 0.0
  endif

  if ( gradient%b < 0.0 .and. gradient%b >= -tol ) then
!    write (*,*) 'Warning: gradient of l(lambda) less than zero ',gradient
    gradient%b = 0.0
  endif

  if ( gradient%a < -tol .or. gradient%b < -tol ) then
    write (*,*) 'Error: gradient of l(lambda) less than zero ',gradient,x
    stop
  endif

  end subroutine lambda


! **********************
! subroutine boundary
! **********************

  subroutine boundary (flow_state_building, region, velocity_profile, x, origin, gradient, r)

  implicit none

! argument list:
  type (flow_state), intent (in)  :: flow_state_building
  character(*),      intent (in)  :: region
  character(*),      intent (in)  :: velocity_profile
  real,              intent (in)  :: x
  type (origin_),    intent (in)  :: origin
  type (gradient_),  intent (in)  :: gradient
  type (r_),         intent (out) :: r

! locals:
  real :: xi_o_y
  real :: xi_o_z
  real :: xi_i_y
  real :: xi_i_z
  real :: q_o_y
  real :: q_o_z
  real :: q_i_y
  real :: q_i_z

! determine functions describing inner and outer boundaries of region

  if ( region == 'outer' ) then

     xi_o_y = xi_e_y(flow_state_building)
     xi_o_z = xi_e_z(flow_state_building)

     xi_i_y = xi_d_y(flow_state_building, x, velocity_profile)
     xi_i_z = xi_d_z(flow_state_building, x, velocity_profile)

  else

     xi_o_y = xi_d_y(flow_state_building, x, velocity_profile)
     xi_o_z = xi_d_z(flow_state_building, x, velocity_profile)

     xi_i_y = xi_r_y(flow_state_building, x)
     xi_i_z = xi_r_z(flow_state_building, x)

  endif


! outer boundary

  if ( (xi_o_z - origin%z) * gradient%a - gradient%b * (xi_o_y - origin%y) >= 0.0 ) then

! l(lambda) intersects vertical boundary

    q_o_y = xi_o_y
    q_o_z = origin%z + gradient%b * (xi_o_y - origin%y)/gradient%a

  else

! l(lambda) intersects horizontal boundary

    q_o_y = origin%y + gradient%a * (xi_o_z - origin%z)/gradient%b
    q_o_z = xi_o_z

  endif


! inner boundary

  if ( flow_state_building%cuboid%wb >= flow_state_building%cuboid%hb ) then

! wide building and cube

    if ( gradient%b == 0.0 ) then

! l(lambda) coincides with z=0

      q_i_y = origin%y
      q_i_z = 0.0

    elseif ( gradient%b * (origin%y - xi_i_y) >= gradient%a * origin%z ) then

! l(lambda) intersects z=0

      q_i_y = origin%y - gradient%a * origin%z/gradient%b
      q_i_z = 0.0

    elseif ( (xi_i_z - origin%z) * gradient%a >= gradient%b * (xi_i_y - origin%y) ) then

! l(lambda) intersects vertical boundary

      q_i_y = xi_i_y
      q_i_z = origin%z + gradient%b * (xi_i_y - origin%y)/gradient%a

    else

! l(lambda) intersects horizontal boundary

      q_i_y = origin%y + gradient%a * (xi_i_z - origin%z)/gradient%b
      q_i_z = xi_i_z

    endif

  else

! tall building

    if ( gradient%a == 0.0 ) then

! l(lambda) coincides with y=0

      q_i_y = 0.0
      q_i_z = origin%z

    elseif ( gradient%a * (origin%z - xi_i_z) >= gradient%b * origin%y ) then

! l(lambda) intersects y=0

      q_i_y = 0.0
      q_i_z = origin%z - gradient%b * origin%y/gradient%a

    elseif ( (xi_i_z - origin%z) * gradient%a >= gradient%b * (xi_i_y - origin%y) ) then

! l(lambda) intersects vertical boundary

      q_i_y = xi_i_y
      q_i_z = origin%z + gradient%b * (xi_i_y - origin%y)/gradient%a

    else

! l(lambda) intersects horizontal boundary

      q_i_y = origin%y + gradient%a * (xi_i_z - origin%z)/gradient%b
      q_i_z = xi_i_z

    endif

  endif


! calculate R

  r%o = sqrt( (q_o_y - origin%y) * (q_o_y - origin%y) + &
              (q_o_z - origin%z) * (q_o_z - origin%z) )

  r%i = sqrt( (q_i_y - origin%y) * (q_i_y - origin%y) + &
              (q_i_z - origin%z) * (q_i_z - origin%z) )

  end subroutine boundary


! ********************
! subroutine position
! ********************

  subroutine position (origin, gradient, r, r_old, r_new, x)

  implicit none

! argument list:
  type (origin_),    intent (in)  :: origin
  type (gradient_),  intent (in)  :: gradient
  real,              intent (in)  :: r
  type (r_),         intent (in)  :: r_old
  type (r_),         intent (in)  :: r_new
  real,              intent (out) :: x(:)

! locals:
  real :: r_prime

! calculate R'

  r_prime = sqrt( (r_new%o * r_new%o - r_new%i * r_new%i) * (r * r - r_old%i * r_old%i)/  &
                  (r_old%o * r_old%o - r_old%i * r_old%i) + r_new%i * r_new%i )

! calculate P_y, P_z

  x(2) = origin%y + r_prime * sqrt( gradient%a * gradient%a/ &
                                    (gradient%a * gradient%a + gradient%b * gradient%b) )

  x(3) = origin%z + r_prime * sqrt( gradient%b * gradient%b/ &
                                    (gradient%a * gradient%a + gradient%b * gradient%b) )

  end subroutine position


! *************************
! subroutine fanning_point
! *************************

  subroutine fanning_point (w, h, gam, x, origin)

  implicit none

! argument list:
  real,           intent (in)  :: w
  real,           intent (in)  :: h
  real,           intent (in)  :: gam
  real,           intent (in)  :: x(:)
  type (origin_), intent (out) :: origin

! locals:
  real :: rho
  real :: cp
  real :: fp
  real :: s_b
  real :: c_b
  real :: s_g
  real :: c_g

! radius of circle

  rho = ( 1.0 - gam + sqrt(1.0 + gam * gam) ) * sqrt(1.0 + gam * gam) * (w - gam * h)/  &
        ( gam * ( 1.0 + gam - sqrt(1.0 + gam * gam) ) )

! lengths, sines and cosines

  cp = sqrt( (x(2) - rho) * (x(2) - rho) +                                &
             (x(3) + (1.0 + sqrt(1.0 + gam * gam)) * (w - gam * h)/gam) * &
             (x(3) + (1.0 + sqrt(1.0 + gam * gam)) * (w - gam * h)/gam) )
  fp = sqrt( cp * cp - rho * rho )

  s_b = rho/cp
  c_b = fp/cp

  s_g = (x(2) - rho)/cp
  c_g = (x(3) + (1.0 + sqrt(1.0 + gam * gam)) * (w - gam * h)/gam)/cp

! Y_0

  origin%y = rho * (1.0 - c_b * c_g + s_b * s_g)

! Z_0

  origin%z = rho * (s_b * c_g + c_b * s_g) - &
             (1.0 + sqrt(1.0 + gam * gam)) * (w - gam * h)/gam

  end subroutine fanning_point


! *****************************
! y component of function xi_e
! *****************************

  real function xi_e_y (flow_state_building)

  implicit none

! argument list:
  type (flow_state), intent(in) :: flow_state_building

  xi_e_y = 0.5 * flow_state_building%cuboid%we

  end function xi_e_y


! *****************************
! z component of function xi_e
! *****************************

  real function xi_e_z (flow_state_building)

  implicit none

! argument list:
  type (flow_state), intent(in) :: flow_state_building

  xi_e_z = flow_state_building%cuboid%he

  end function xi_e_z


! *****************************
! y component of function xi_d
! *****************************

  real function xi_d_y (flow_state_building, x, velocity_profile)

  implicit none

! argument list:
  type (flow_state), intent(in) :: flow_state_building
  real,              intent(in) :: x
  character(*),      intent(in) :: velocity_profile

! locals:
  real :: m
  real :: c
  real :: b
  real :: distance_oa

! wide building and cube

  if ( flow_state_building%cuboid%wb >= flow_state_building%cuboid%hb ) then

! xi_d_y upwind

    if ( x < -0.5 * flow_state_building%cuboid%lb ) then

      m = flow_state_building%shape%gamma * flow_state_building%cuboid%hb *                 &
          ( sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) - &
            flow_state_building%shape%alpha )/flow_state_building%cuboid%lu

      c = sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) * &
          flow_state_building%shape%gamma * flow_state_building%cuboid%hb

      xi_d_y = m * (x + 0.5 * flow_state_building%cuboid%lb) + c +                     &
               0.5 * flow_state_building%cuboid%wb - flow_state_building%shape%gamma * &
               flow_state_building%cuboid%hb

! xi_d_y above building and downwind recirculation region

    elseif ( x <= 0.5 * flow_state_building%cuboid%lb + flow_state_building%cuboid%ldr ) then

      xi_d_y = sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) * &
               flow_state_building%shape%gamma * flow_state_building%cuboid%hb +               &
               0.5 * flow_state_building%cuboid%wb - flow_state_building%shape%gamma *         &
               flow_state_building%cuboid%hb

! xi_d_y downwind

    else

      b = flow_state_building%shape%alpha * flow_state_building%cuboid%hb *       &
          ( flow_state_building%shape%alpha * flow_state_building%shape%gamma *   &
            flow_state_building%cuboid%hb + 0.5 * flow_state_building%cuboid%wb - &
            flow_state_building%shape%gamma * flow_state_building%cuboid%hb )

      distance_oa = 0.5 * flow_state_building%cuboid%wb - flow_state_building%shape%gamma * &
                    flow_state_building%cuboid%hb

      xi_d_y = distance_oa + sqrt( distance_oa * distance_oa + 4.0 * flow_state_building%shape%gamma * &
                                   b * jacob_downwind(x, flow_state_building, velocity_profile) )
      xi_d_y = 0.5 * xi_d_y

    endif

! tall building

  else

! xi_d_y upwind

    if ( x < -0.5 * flow_state_building%cuboid%lb ) then

      m = 0.5 * flow_state_building%cuboid%wb *                                             &
          ( sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) - &
            flow_state_building%shape%alpha )/flow_state_building%cuboid%lu

      c = sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) * &
          0.5 * flow_state_building%cuboid%wb

      xi_d_y = m * (x + 0.5 * flow_state_building%cuboid%lb) + c

! xi_d_y above building and downwind recirculation region

    elseif ( x <= 0.5 * flow_state_building%cuboid%lb + flow_state_building%cuboid%ldr ) then

      xi_d_y = sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) * &
               0.5 * flow_state_building%cuboid%wb

! xi_d_y downwind

    else

      b = flow_state_building%shape%alpha * flow_state_building%cuboid%wb *         &
          ( flow_state_building%shape%alpha * flow_state_building%shape%gamma *     &
            0.5 * flow_state_building%cuboid%wb + flow_state_building%cuboid%hb -   &
            flow_state_building%shape%gamma * 0.5 * flow_state_building%cuboid%wb )

      distance_oa = flow_state_building%cuboid%hb - flow_state_building%shape%gamma * 0.5 * &
                    flow_state_building%cuboid%wb

      xi_d_y = -distance_oa + sqrt( distance_oa * distance_oa + 2.0 * b *                      &
                                    flow_state_building%shape%gamma *                          &
                                    jacob_downwind(x, flow_state_building, velocity_profile) )
      xi_d_y = 0.5 * xi_d_y/flow_state_building%shape%gamma

    endif

  endif

  end function xi_d_y


! *****************************
! z component of function xi_d
! *****************************

  real function xi_d_z (flow_state_building, x, velocity_profile)

  implicit none

! argument list:
  type (flow_state), intent(in) :: flow_state_building
  real,              intent(in) :: x
  character(*),      intent(in) :: velocity_profile

! locals:
  real :: m
  real :: c
  real :: b
  real :: distance_oa

! wide building and cube

  if ( flow_state_building%cuboid%wb >= flow_state_building%cuboid%hb ) then

! xi_d_z upwind

    if ( x < -0.5 * flow_state_building%cuboid%lb ) then

      m = ( sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) - &
            flow_state_building%shape%alpha ) * flow_state_building%cuboid%hb/flow_state_building%cuboid%lu

      c = sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) * &
          flow_state_building%cuboid%hb

      xi_d_z = m * (x + 0.5 * flow_state_building%cuboid%lb) + c

! xi_d_z above building and downwind recirculation region

    elseif ( x <= 0.5 * flow_state_building%cuboid%lb + flow_state_building%cuboid%ldr ) then

      xi_d_z = sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) * &
               flow_state_building%cuboid%hb

! xi_d_z downwind

    else

      b = flow_state_building%shape%alpha * flow_state_building%cuboid%hb *       &
          ( flow_state_building%shape%alpha * flow_state_building%shape%gamma *   &
            flow_state_building%cuboid%hb + 0.5 * flow_state_building%cuboid%wb - &
            flow_state_building%shape%gamma * flow_state_building%cuboid%hb )

      distance_oa = 0.5 * flow_state_building%cuboid%wb - flow_state_building%shape%gamma * &
                    flow_state_building%cuboid%hb

      xi_d_z = -distance_oa + sqrt( distance_oa * distance_oa + 4.0 * flow_state_building%shape%gamma * &
                                    b *  jacob_downwind(x, flow_state_building, velocity_profile) )
      xi_d_z = 0.5 * xi_d_z/flow_state_building%shape%gamma

    endif

! tall building

  else

! xi_d_z upwind

    if ( x < -0.5 * flow_state_building%cuboid%lb ) then

      m = flow_state_building%shape%gamma * 0.5 * flow_state_building%cuboid%wb *           &
          ( sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) - &
            flow_state_building%shape%alpha )/flow_state_building%cuboid%lu

      c = sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) * &
          flow_state_building%shape%gamma * 0.5 * flow_state_building%cuboid%wb

      xi_d_z = m * (x + 0.5 * flow_state_building%cuboid%lb) + c + flow_state_building%cuboid%hb - &
               flow_state_building%shape%gamma * 0.5 * flow_state_building%cuboid%wb

! xi_d_z above building and downwind recirculation region

    elseif ( x <= 0.5 * flow_state_building%cuboid%lb + flow_state_building%cuboid%ldr ) then

      xi_d_z = sqrt(1.0 + flow_state_building%shape%alpha * flow_state_building%shape%alpha) * &
               flow_state_building%shape%gamma * 0.5 * flow_state_building%cuboid%wb +         &
               flow_state_building%cuboid%hb - flow_state_building%shape%gamma *               &
               0.5 * flow_state_building%cuboid%wb

! xi_d_z downwind

    else

      b = flow_state_building%shape%alpha * flow_state_building%cuboid%wb *         &
          ( flow_state_building%shape%alpha * flow_state_building%shape%gamma *     &
            0.5 * flow_state_building%cuboid%wb + flow_state_building%cuboid%hb -   &
            flow_state_building%shape%gamma * 0.5 * flow_state_building%cuboid%wb )

      distance_oa = flow_state_building%cuboid%hb - flow_state_building%shape%gamma * 0.5 * &
                    flow_state_building%cuboid%wb

      xi_d_z = distance_oa + sqrt( distance_oa * distance_oa + 2.0 * flow_state_building%shape%gamma * &
                                   b * jacob_downwind(x, flow_state_building, velocity_profile) )
      xi_d_z = 0.5 * xi_d_z

    endif

  endif

  end function xi_d_z


! *****************************
! y component of function xi_r
! *****************************

  real function xi_r_y (flow_state_building, x)

  implicit none

! argument list:
  type (flow_state), intent(in) :: flow_state_building
  real,              intent(in) :: x

  if ( x < -0.5*flow_state_building%cuboid%lb - flow_state_building%cuboid%lur ) then

    xi_r_y = 0.0

  elseif ( x <= -0.5*flow_state_building%cuboid%lb ) then

    xi_r_y = 0.5 * flow_state_building%cuboid%wb *                                       &
             (x + 0.5 * flow_state_building%cuboid%lb + flow_state_building%cuboid%lur)/ &
             flow_state_building%cuboid%lur

  elseif ( x <= 0.5*flow_state_building%cuboid%lb ) then

    xi_r_y = 0.5 * flow_state_building%cuboid%wb

  elseif ( x < flow_state_building%cuboid%ldr + 0.5*flow_state_building%cuboid%lb ) then

    xi_r_y = 0.5 * flow_state_building%cuboid%wb *                                       &
             (0.5 * flow_state_building%cuboid%lb + flow_state_building%cuboid%ldr - x)/ &
             flow_state_building%cuboid%ldr

  else

    xi_r_y = 0.0

  endif

  end function xi_r_y


! *****************************
! z component of function xi_r
! *****************************

  real function xi_r_z (flow_state_building,x)

  implicit none

! argument list:
  type (flow_state), intent(in) :: flow_state_building
  real,              intent(in) :: x

! locals
  real    :: a
  real    :: b
  real    :: c
  real    :: x_0
  real    :: z_max
  logical :: leading

! leading or trailing roof top seperation based on Fackrell(1982) ?.

  if (flow_state_building%cuboid%lb >= min(flow_state_building%cuboid%hb,            &
                                           0.5*flow_state_building%cuboid%wb)) then
    leading = .false.
!    write(11,*)'Trailing edge seperation'
  else
    leading = .true.
!    write(11,*)'Leading edge seperation'
  endif

  if (leading) then

! Leading edge separation.
! From Apsley(1988) based on Fackrell(1982) data)

    z_max = flow_state_building%cuboid%hb*(1                                     &
            + 0.7*(1 - exp(-(flow_state_building%cuboid%wb -                     &
            2.0*flow_state_building%cuboid%lb)/flow_state_building%cuboid%hb)))



    a     = z_max - flow_state_building%shape%zs

    c     = ( (z_max - flow_state_building%shape%zs)                                  &
            * (z_max - flow_state_building%shape%zs)                                  &
            - (flow_state_building%cuboid%hb - flow_state_building%shape%zs)          &
            * (flow_state_building%cuboid%hb - flow_state_building%shape%zs) )        &
            / ( (z_max - flow_state_building%shape%zs)                                &
            * (z_max - flow_state_building%shape%zs) )

    x_0   = ( (0.5 * flow_state_building%cuboid%lb + flow_state_building%cuboid%ldr)  &
             * sqrt(c) - 0.5 * flow_state_building%cuboid%lb )/(1.0 + sqrt(c))

    b     = 0.5 * flow_state_building%cuboid%lb + flow_state_building%cuboid%ldr

    if ( x < -0.5*flow_state_building%cuboid%lb - flow_state_building%cuboid%lur ) then

      xi_r_z = 0.0

    elseif ( x <= -0.5*flow_state_building%cuboid%lb ) then

      xi_r_z = sqrt( flow_state_building%cuboid%hb * flow_state_building%cuboid%hb *  &
               (x + 0.5 * flow_state_building%cuboid%lb                               &
               + flow_state_building%cuboid%lur) / flow_state_building%cuboid%lur)

    elseif ( x < 0.5*flow_state_building%cuboid%lb + flow_state_building%cuboid%ldr ) then

      xi_r_z = a * sqrt(1.0 - (x - x_0) * (x - x_0)/(b * b))

    else

      xi_r_z = 0.0

    endif

  else

! trailing edge separation

    if ( x < -0.5 * flow_state_building%cuboid%lb - flow_state_building%cuboid%lur ) then

      xi_r_z = 0.0

    elseif ( x <= -0.5*flow_state_building%cuboid%lb ) then

      xi_r_z = sqrt( flow_state_building%cuboid%hb * flow_state_building%cuboid%hb *   &
                     (x + 0.5 * flow_state_building%cuboid%lb                          &
                      + flow_state_building%cuboid%lur)/flow_state_building%cuboid%lur )

    elseif ( x <= 0.5*flow_state_building%cuboid%lb ) then

      xi_r_z = flow_state_building%cuboid%hb

    elseif ( x < 0.5 * flow_state_building%cuboid%lb + flow_state_building%cuboid%ldr ) then

      xi_r_z = sqrt( flow_state_building%cuboid%hb * flow_state_building%cuboid%hb * &
                     (flow_state_building%cuboid%ldr +                               &
                      0.5 * flow_state_building%cuboid%lb - x)/                      &
                      flow_state_building%cuboid%ldr)

    else

      xi_r_z = 0.0

    endif

  endif

  end function xi_r_z


! *****************************
! function jacob_downwind
! *****************************

  real function jacob_downwind(x, flow_state_building, velocity_profile)

  implicit none

! argument list:
  real,              intent(in) :: x
  type (flow_state), intent(in) :: flow_state_building
  character(*),      intent(in) :: velocity_profile

! locals:
  real :: z
  real :: xi
  real :: ratio

! variable xi

  z = flow_state_building%shape%beta * flow_state_building%cuboid%hb
  xi = z/lambda_z(flow_state_building, x, z, velocity_profile)

! ratio u(x,y,z)/u_mean(x,y,z)

  ratio = 1.0 - flow_state_building%cuboid%wb * flow_state_building%cuboid%hb *     &
                flow_state_building%cuboid%hb * g(xi) * f(flow_state_building, x)/  &
                ( 5.0 * flow_state_building%constants%pi *                          &
                  lambda_y(flow_state_building, x, z, velocity_profile) *           &
                  lambda_z(flow_state_building, x, z, velocity_profile) *           &
                  lambda_z(flow_state_building, x, z, velocity_profile) )

! Jacobian J_D

  if ( ratio < 0.0 ) then
     write (6,*) 'Error: Jacobian J_D less than zero',1.0/ratio
     stop
  endif

  jacob_downwind = 1.0/ratio

  end function jacob_downwind


! *****************************
! function g
! *****************************

  real function g(xi)

  implicit none

! argument list:
  real, intent(in) :: xi

  g = 0.5 * xi * exp(-0.25 * xi * xi)

  end function g


! *****************************
! function f
! *****************************

  real function f(flow_state_building, x)

  implicit none

! argument list:
  type (flow_state), intent (in) :: flow_state_building
  real,              intent (in) :: x

  f = sin( 0.5 * flow_state_building%constants%pi *                            &
           (1.0 + (x - 0.5 * flow_state_building%cuboid%lb -                   &
           flow_state_building%cuboid%ldr) /                                   &
           (flow_state_building%cuboid%lw - flow_state_building%cuboid%ldr)) )

  end function f


! *****************************
! function lambda_y
! *****************************

  real function lambda_y(flow_state_building, x, z, velocity_profile)

  implicit none

! argument list:
  type (flow_state), intent(in) :: flow_state_building
  real,              intent(in) :: x
  real,              intent(in) :: z
  character(*),      intent(in) :: velocity_profile

! locals:
  real :: dy

  dy = flow_state_building%met%k * flow_state_building%met%u_star * flow_state_building%cuboid%hb

  lambda_y = sqrt( dy * (x - flow_state_building%constants%x_0)/    &
                   u_mean(z, flow_state_building, velocity_profile) )

  end function lambda_y


! *****************************
! function lambda_z
! *****************************

  real function lambda_z(flow_state_building, x, z, velocity_profile)

  implicit none

! argument list:
  type (flow_state), intent(in) :: flow_state_building
  real,              intent(in) :: x
  real,              intent(in) :: z
  character(*),      intent(in) :: velocity_profile

! locals:
  real :: dz

  dz = 2.0 * flow_state_building%met%k*flow_state_building%met%u_star * flow_state_building%cuboid%hb

  lambda_z = sqrt( dz * (x - flow_state_building%constants%x_0)/   &
                   u_mean(z, flow_state_building, velocity_profile) )

  end function lambda_z


! ************************************************************
! ************************************************************

! **********************
! subroutine reflection
! **********************

  subroutine reflection (cuboid, x_old, x, u_prime, reflect)

  implicit none

! argument list:
  type (dimensions), intent (in)    :: cuboid
  real,              intent (inout) :: x_old(:)
  real,              intent (inout) :: x(:)
  real,              intent (inout) :: u_prime(:)
  logical,           intent (in)    :: reflect

! **********************************************
! check: does particle enter extended building?
! **********************************************

  if ( abs(x(1)) < 0.5*cuboid%lb .and. &
       abs(x(2)) < 0.5*cuboid%wb .and. &
       abs(x(3)) < cuboid%hb    ) then

     call reflection_building (cuboid, x_old, x, u_prime)

  else

     x_old=x

  endif

! ***********************
! reflection from ground
! ***********************

  if ( x(3) < 0.0 ) then
    x(3) = -x(3)
    u_prime(3) = -u_prime(3)
  endif

  if (reflect) then

! ****************************
! reflection from upwind edge
! ****************************

    if ( x(1) < -(cuboid%lu + 0.5*cuboid%lb) ) then
      x(1) = -2.0 * (cuboid%lu + 0.5*cuboid%lb)-x(1)
      u_prime(1) = -u_prime(1)
    endif

! **********************
! reflection from sides
! **********************

    if ( abs(x(2)) > 0.5*cuboid%we ) then
      if ( x(2) < 0.0 ) then
        x(2) = -cuboid%we-x(2)
      else
        x(2) = cuboid%we-x(2)
      endif
      u_prime(2) = -u_prime(2)
    endif

! ********************
! reflection from top
! ********************

    if ( x(3) > cuboid%he ) then
      x(3) = 2.0*cuboid%he-x(3)
      u_prime(3) = -u_prime(3)
    endif

  endif

  end subroutine reflection


! *******************************
! subroutine reflection_building
! *******************************

  subroutine reflection_building (cuboid, x_old, x, u_prime)

  implicit none

! argument list:
  type (dimensions), intent (in)    :: cuboid
  real,              intent (in)    :: x_old(:)
  real,              intent (inout) :: x(:)
  real,              intent (inout) :: u_prime(:)

!  locals:
  logical :: hit
  logical :: back
  logical :: left
  real    :: tol

  hit=.false.

  tol=cuboid%hb*epsilon(x(1))

!  *******************************
!  particle is upwind of building
!  *******************************

  if (x_old(1) < -0.5*cuboid%lb) then

     back=.false.

!  check: does particle hit front?

     call front (hit,back,cuboid,tol,x_old,x,u_prime)

!  check: does particle hit top?

     if ( (.not.hit) .and. (x_old(3) >= cuboid%hb) ) then
       call top (hit, cuboid, tol, x_old, x, u_prime)
     endif

!  check: does particle hit side?

     if (.not.hit) then
       if (x_old(2) < 0.0) then
         left = .true.
       else
         left = .false.
       endif
       call side (hit, left, cuboid, tol, x_old, x, u_prime)
     endif

!  ***************************************
!  particle is adjacent/above to building
!  ***************************************

  elseif (x_old(1) < 0.5*cuboid%lb) then

!  check: does particle hit top?

     if (x_old(3) >= cuboid%hb) then
       call top (hit,cuboid,tol,x_old,x,u_prime)
     endif

!  check: does particle hit side?

     if (.not.hit) then
       if (x_old(2) < 0.0) then
         left = .true.
       else
         left = .false.
       endif
       call side (hit,left,cuboid,tol,x_old,x,u_prime)
     endif

!  *********************************
!  particle is downwind of building
!  *********************************

  else

     back=.true.

!  check: does particle hit front?

     call front (hit,back,cuboid,tol,x_old,x,u_prime)

!  check: does particle hit top?

     if ( (.not.hit) .and. (x_old(3) >= cuboid%hb) ) then
       call top (hit,cuboid,tol,x_old,x,u_prime)
     endif

!  check: does particle hit side?

     if (.not.hit) then
       if (x_old(2) < 0.0) then
         left = .true.
       else
         left = .false.
       endif
       call side (hit,left,cuboid,tol,x_old,x,u_prime)
     endif

  endif

!  set x equal to x_old if absolute difference less than tolerance

  if (.not.hit) then
     x = x_old
  endif

  end subroutine reflection_building


!  *****************
!  subroutine front
!  *****************

  subroutine front (hit,back,cuboid,tol,x_old,x,u_prime)

  implicit none

!  argument list:
  logical,           intent (inout) :: hit
  logical,           intent (in)    :: back
  type (dimensions), intent (in)    :: cuboid
  real,              intent (in)    :: tol
  real,              intent (in)    :: x_old(:)
  real,              intent (inout) :: x(:)
  real,              intent (inout) :: u_prime(:)

!  locals:
  real :: p0
  real :: lambda
  real :: yp
  real :: zp

  if (back) then
    p0 =  0.5 * cuboid%lb
  else
    p0 = -0.5  *cuboid%lb
  endif

  if (abs(x(1) - x_old(1)) < tol) then
    hit = .false.
  else

!  calculate point of intersection with extrapolated face

    lambda = (p0 - x_old(1))/(x(1) - x_old(1))

    yp     = x_old(2) + lambda*(x(2) - x_old(2))
    zp     = x_old(3) + lambda*(x(3) - x_old(3))

!  check: does the intersection point lie on face?

    if ( (abs(yp) <= 0.5*cuboid%wb) .and. (abs(zp) <= cuboid%hb) ) then
      hit        = .true.
      x(1)       = 2.0*p0 - x(1)
      u_prime(1) = -u_prime(1)
    endif

  endif

  end subroutine front


!  ***************
!  subroutine top
!  ***************

  subroutine top (hit,cuboid,tol,x_old,x,u_prime)

  implicit none

!  argument list:
  logical,           intent (inout) :: hit
  type (dimensions), intent (in)    :: cuboid
  real,              intent (in)    :: tol
  real,              intent (in)    :: x_old(:)
  real,              intent (inout) :: x(:)
  real,              intent (inout) :: u_prime(:)

!  locals:
  real :: lambda
  real :: xp,yp

  if (abs(x(3)-x_old(3)) < tol) then
    hit = .false.
  else

!  calculate point of intersection with extrapolated face

    lambda = (cuboid%hb - x_old(3))/(x(3) - x_old(3))

    xp     = x_old(1) + lambda*(x(1) - x_old(1))
    yp     = x_old(2) + lambda*(x(2) - x_old(2))

!  check: does the intersection point lie on face?

    if ( (abs(xp) <= 0.5*cuboid%lb) .and. (abs(yp) <= 0.5*cuboid%wb) ) then
      hit        = .true.
      x(3)       = 2.0*cuboid%hb - x(3)
      u_prime(3) = -u_prime(3)
    endif

  endif

  end subroutine top


!  ****************
!  subroutine side
!  ****************

  subroutine side (hit,left,cuboid,tol,x_old,x,u_prime)

  implicit none

!  argument list:
  logical,           intent (inout) :: hit
  logical,           intent (in)    :: left
  type (dimensions), intent (in)    :: cuboid
  real,              intent (in)    :: tol
  real,              intent (in)    :: x_old(:)
  real,              intent (inout) :: x(:)
  real,              intent (inout) :: u_prime(:)

!  locals:
  real :: q0
  real :: lambda
  real :: xp,zp

  if (left) then
    q0 = -0.5*cuboid%wb
  else
    q0 = 0.5*cuboid%wb
  endif

  if (abs(x(2)-x_old(2)) < tol) then
    hit = .false.
  else

!  calculate point of intersection with extrapolated face

    lambda = (q0 - x_old(2))/(x(2) - x_old(2))

    xp     = x_old(1) + lambda*(x(1) - x_old(1))
    zp     = x_old(3) + lambda*(x(3) - x_old(3))

!  check: does the intersection point lie on face?

    if ( (abs(xp) <= 0.5*cuboid%lb) .and. (abs(zp) <= cuboid%hb) ) then
      hit        = .true.
      x(2)       = 2.0*q0 - x(2)
      u_prime(2) = -u_prime(2)
    endif

  endif

  end subroutine side


!  ************************
!  subroutine domain_check
!  ************************

  subroutine domain_check (x,cuboid,domain,finished)

  implicit none

!  argument list:
  real,              intent(in)    :: x(:)
  type (dimensions), intent(in)    :: cuboid
  logical,           intent(inout) :: domain
  logical,           intent(inout) :: finished

!  ******************************
!  check: is particle in domain?
!  ******************************

  if ( abs(x(1)) < 0.5*cuboid%lb .and. &
       abs(x(2)) < 0.5*cuboid%wb .and. &
       abs(x(3)) < cuboid%hb ) then

     domain=.false.

  elseif ( x(1)      >= -(0.5*cuboid%lb+cuboid%lu) .and. &
           x(1)      <= (0.5*cuboid%lb+cuboid%lw)  .and. &
           abs(x(2)) <= 0.5*cuboid%we              .and. &
           x(3)      <= cuboid%he                  .and. &
           x(3)      >= 0.0                       ) then

     domain=.true.

  elseif ( x(1) > (0.5*cuboid%lb+cuboid%lw) ) then

     finished=.true.

  endif

  end subroutine domain_check


!***************************************************************
!***************************************************************
!  subroutine vortex_parameters
!***************************************************************

  subroutine vortex_parameters(flow_state_building,velocity_profile)

  implicit none

! argument list:
  type (flow_state), intent(inout) :: flow_state_building

  character(*),      intent(in)    :: velocity_profile

! locals:
  type (turbulence) :: building_turbulence
  real              :: fz

  Type(ShortTime_) :: Time
  Type(Flow_)      :: Flow
  Type(HCoeffs_)   :: HCoeffs
  Type(ZCoeffs_)   :: ZCoeffs
  Type(TCoeffs_)   :: TCoeffs

  HCoeffs%Valid = .false. ! Could save for more efficient work with NWP data? $$
  ZCoeffs%Valid = .false.
  TCoeffs%Valid = .false.


! growth rate of vortex radius with downstream distance x(1)
! relating the vortex spread rate to vertical turbulent diffusion

  call TurbFromFlowField(                                          &
         flow_state_building%Coords, flow_state_building%Grids,    &
         Time,                                                     &
         (/0.0, 0.0, flow_state_building%cuboid%hb/), flow_state_building%cuboid%hb, &
         flow_state_building%FlowField,                            &
         .true., .false., .false., 0.0, 0.0,                       &
         Flow,                                                     &
         HCoeffs, ZCoeffs, TCoeffs                                 &
       )

!  flow_state_building%vortices%growth_rate =  0.7 * SQRT(SigUU(3))            &  ! 0.7
!                              / u_mean(flow_state_building%cuboid%hb,         &
!                              flow_state_building,velocity_profile)

  flow_state_building%vortices%growth_rate = 11.0 * SQRT(Flow%SigUU(3))       &
                              / u_mean(flow_state_building%cuboid%hb,  &
                              flow_state_building,velocity_profile)


! we set vortex size at downstream edge of building
! to 1/10 of building height (we are only considering
! a cube at the moment so actual building length scale
! not critical.

  flow_state_building%vortices%r_c_initial = 0.5*min(flow_state_building%cuboid%hb,  &
                                                 0.8*flow_state_building%cuboid%wb,  &   ! 0.8
                                                 flow_state_building%cuboid%lb)



! we set maximum vortex velocity at downstream edge of building
! i.e., Vt which occurs at r_c
! make Vt dependant on building angle


  flow_state_building%vortices%Vt_initial =                               &
                    sin(2.0 * flow_state_building%vortices%building_phi)  &
                     * u_mean(flow_state_building%cuboid%hb,              &
                       flow_state_building,velocity_profile)/4.7  +       & ! 4.7 decided on  6/6/02.
                        sin(flow_state_building%constants%pi / 30.0)      & ! fixed vortex strength
                        * u_mean(flow_state_building%cuboid%hb,           & ! for zero degree downwash
                        flow_state_building,velocity_profile)/4.5           ! 4.5

! then we calculate the total vorticity

  flow_state_building%vortices%vorticity =                                   &
                        (2.0 * flow_state_building%vortices%Vt_initial     &
                          / flow_state_building%vortices%r_c_initial)      &
                           * flow_state_building%constants%pi              &
                           * flow_state_building%vortices%r_c_initial**2

  flow_state_building%vortices%vorticity(1) = -flow_state_building%vortices%vorticity(1)
  flow_state_building%vortices%vorticity(4) = -flow_state_building%vortices%vorticity(4)

! this then becomes our concerved quantity

  write (6,*)  'Vorticies On'
  write (11,*) 'Vortex Parameters'

  write (11,*) 'X position of building trailing edge  ',               &
                flow_state_building%vortices%building_trailing_edge

  write (11,*) 'Vortices radius r_c at building trailing edge  ',      &
                flow_state_building%vortices%r_c_initial

  write (11,*) 'Vortices max Velocity Vt at building trailing edge  ', &
                flow_state_building%vortices%Vt_initial

  write (11,*) 'Vortex growth rate  ',flow_state_building%vortices%growth_rate

end subroutine vortex_parameters


!***********************************************************
!  subroutine vortex_axis_position
!
! NOTE: called only once at the moment because axial position
! is constant of all downstream positions.
!
!***********************************************************

subroutine vortex_axis_position(flow_state_building)

  implicit none

  type (flow_state) :: flow_state_building

!
! these positional data could be represented as functions
! of x. Also for building orientations other than zero the
! lateral start position should relate to the
! most lateral point of the building on each side.
!
! A relationship can be derived to predict the effect of each
! vortex on the others. This is represented by a velocity
! that advects the core of each vortex. This will result in the
! downward movement of the 'real' vortices and ten their
! lateral movement away from the buildign centre line as
! the effect from the reflected vortices becomes more significant.
!V = gamma/(2 pi b) where V is the advecting velocity,
!gamma the vortex circulation and b the vortex pair seperation.
!

! primary 'real' voticies

! LATERAL POSITION MUST BE CALCULATED USING PROJECTED FRONTAL DIMENSIONS.
! This is VERY important!!!!!!!!

  flow_state_building%vortices%vortex_axis(1,1) = -                                          &
       ( flow_state_building%cuboid%wb**2 + flow_state_building%cuboid%lb**2 )**0.5 *        &
       cos((atan(flow_state_building%cuboid%lb/flow_state_building%cuboid%wb)) -             &
       (flow_state_building%vortices%building_phi)) * 0.5 +                                  &
       (0.25 * flow_state_building%vortices%r_c_initial *                                    &
       sin(flow_state_building%vortices%building_phi*2.0))


! EXTRA HEIGHT ABOVE ROOF LEVEL BASED ON VORTEX SIZE.
! If done as multiple of Hb then could end up very high above building for tall buildings.
  flow_state_building%vortices%vortex_axis(2,1) =  flow_state_building%cuboid%hb  +   &
                                                   0.15 * flow_state_building%vortices%r_c_initial  ! z pos

  flow_state_building%vortices%vortex_axis(1,2) =  - flow_state_building%vortices%vortex_axis(1,1)  ! y pos
  flow_state_building%vortices%vortex_axis(2,2) =  flow_state_building%vortices%vortex_axis(2,1)    ! z pos


! reflection vortices
  flow_state_building%vortices%vortex_axis(1,3) =  flow_state_building%vortices%vortex_axis(1,1)
  flow_state_building%vortices%vortex_axis(2,3) = -flow_state_building%vortices%vortex_axis(2,1)
  flow_state_building%vortices%vortex_axis(1,4) =  flow_state_building%vortices%vortex_axis(1,2)
  flow_state_building%vortices%vortex_axis(2,4) = -flow_state_building%vortices%vortex_axis(2,2)


  write (11,*) 'Vortex (1) axis position (y,z) =  ',              &
                flow_state_building%vortices%vortex_axis(1,1),    &
                flow_state_building%vortices%vortex_axis(2,1)

end subroutine vortex_axis_position


!*************************************************************************
!*************************************************************************
!
!  vortex wake main routine
!
!  Owner: Matthew Hort
!  Date : Nov 2001
!
!*************************************************************************
!*************************************************************************

subroutine building_vortices(x,flow_state_building,u_vortex)

 implicit none

 type (flow_state),  intent(in)    :: flow_state_building

 real,  intent(in)    :: x(3)               ! particle coordinates
 real,  intent(out)   :: u_vortex(3)        ! total vortex induced velocity

! local
 integer              :: j                  ! loop counter
 real                 :: x_vort(3)          ! vortex axis coordinates
 real                 :: r                  ! radius
 real                 :: vortex_u_r(3,4)    ! velocity components @ r
 real                 :: vortex_axis(2,4)    ! velocity components @ r

! initialise arrays
 vortex_u_r = 0.0
 u_vortex   = 0.0

! particle must be downwind of building trailing edge
! to be affected by vortices
 if (x(1) > flow_state_building%vortices%building_trailing_edge) then

!
! adjust vortex axis height as function of downstream distance.
! this may help glc especially with tall buildings.
! The configuration used here is far from finalised. MH 27/02/03
!
   vortex_axis(1,1) = flow_state_building%vortices%vortex_axis(1,1)
   vortex_axis(1,2) = flow_state_building%vortices%vortex_axis(1,2)
   vortex_axis(1,3) = flow_state_building%vortices%vortex_axis(1,3)
   vortex_axis(1,4) = flow_state_building%vortices%vortex_axis(1,4)

! cube and wide building case
   If (flow_state_building%cuboid%hb < flow_state_building%cuboid%wb + 0.1) Then ! 0.1 is a tolerance

     vortex_axis(2,1) = flow_state_building%vortices%vortex_axis(2,1)
     vortex_axis(2,2) = flow_state_building%vortices%vortex_axis(2,2)
     vortex_axis(2,3) = flow_state_building%vortices%vortex_axis(2,3)
     vortex_axis(2,4) = flow_state_building%vortices%vortex_axis(2,4)

   Else   ! tall buildings only

     vortex_axis(2,1) =    &
         max((flow_state_building%cuboid%hb + (0.15 * flow_state_building%vortices%r_c_initial)    &
         - (x(1)-flow_state_building%vortices%building_trailing_edge)**0.5),flow_state_building%cuboid%wb)

     vortex_axis(2,2) =  vortex_axis(2,1)
     vortex_axis(2,3) = -vortex_axis(2,1)
     vortex_axis(2,4) = -vortex_axis(2,1)

   EndIf

! loop through effect from each of the 4 vortices

     vortices_loop: do j = 1,4

! radius of particle from vortex axis

       r = sqrt( (x(2) - vortex_axis(1,j))**2 + (x(3) - vortex_axis(2,j))**2 )

       x_vort(1) = x(1)
       x_vort(2) = x(2) - vortex_axis(1,j)
       x_vort(3) = x(3) - vortex_axis(2,j)


! velocity components from each vortex

       call vortex_velocity(flow_state_building,x_vort,r,j,vortex_u_r)

     end do vortices_loop

! advection of particle based on sum of all velocities

   call particle_velocity(vortex_u_r,u_vortex)

 else

   u_vortex = 0.0

 endif

end subroutine building_vortices


!***********************************************************
!  subroutine vortex_velocity
!***********************************************************

subroutine vortex_velocity(flow_state_building,x_vort,r,j,vortex_u_r)

  implicit none

  type (flow_state),  intent(in) :: flow_state_building

  integer, intent(in)  :: j                   ! loop counter
  real,    intent(in)  :: x_vort(3)           ! vortex axis coordinates
  real,    intent(in)  :: r                   ! radius
  real,    intent(out) :: vortex_u_r(3,4)     ! velocity components @ r

! local
  real    :: local_core_radius
  real    :: local_Vt
  real    :: vt_r
  integer :: i

! radius of core at downstream position x(1)

  local_core_radius = flow_state_building%vortices%r_c_initial * (1.0 + &
                    ((x_vort(1) - flow_state_building%vortices%building_trailing_edge)/       &
                    flow_state_building%vortices%building_trailing_edge)**(0.5)  &
                    * flow_state_building%vortices%growth_rate)

! maximum tangental velocity

  local_Vt = (flow_state_building%vortices%vorticity(j)      &
              / (2.0 * flow_state_building%constants%pi      &
                * local_core_radius))

! tangental velocity at radius r

  if (abs(r) <= local_core_radius) then

    vt_r = local_Vt * r / local_core_radius

  else

    vt_r = local_Vt * local_core_radius / r

  endif

! velocity components

  vortex_u_r(2,j) = -vt_r * x_vort(3) / sqrt(x_vort(2)**2 + x_vort(3)**2)

  vortex_u_r(3,j) =  vt_r * x_vort(2) / sqrt(x_vort(2)**2 + x_vort(3)**2)

!write(36,*)local_Vt,local_core_radius,r,vt_r

end subroutine vortex_velocity


!***********************************************************
!  subroutine particle_velocity
!***********************************************************

subroutine particle_velocity(vortex_u_r,u_vortex)

  implicit none

!  type (flow_state) :: flow_state_building

  real, intent(in)    :: vortex_u_r(3,4)     ! velocity components @ r
  real, intent(out)   :: u_vortex(3)         ! total vortex induced velocity


  u_vortex(1) = 0.0

  u_vortex(2) = + vortex_u_r(2,1) + vortex_u_r(2,2)   &
                + vortex_u_r(2,3) + vortex_u_r(2,4)

  u_vortex(3) = + vortex_u_r(3,1) + vortex_u_r(3,2)   &
                + vortex_u_r(3,3) + vortex_u_r(3,4)

end subroutine particle_velocity


!***********************************************************

end module buildingmodule

!***********************************************************
!***********************************************************


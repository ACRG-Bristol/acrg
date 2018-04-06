! Module:  Semi Lagrangian

Module SemiLagrangianModule

! This module provides the semi-Lagrangian advection scheme

! **********************************************************
! $$ issue with filter in mass conservation routine - check!
! **********************************************************

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use SpeciesModule, Only: Con ! $$ used only to switch concentration precision for chemistry. Remove if it turns out its not needed.

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: sl_init
Public :: sl_advect

   integer, Save                             :: Ni, Nj, Nk
   integer, Save                             :: istrt, iend
   integer, Save                             :: jstrt, jend
   integer, Save                             :: i_p_start, i_p_end
   integer, Save                             :: j_p_start, j_p_end
   integer, Save                             :: halo

   real(4), Save, allocatable, dimension(:)     :: eta_levels, deta
   real(4), Save, allocatable, dimension(:)     :: deta_ratio_m, deta_ratio_p
   real(4), Save, allocatable, dimension(:)     :: x_pnts, y_pnts
   real(4), Save, allocatable, dimension(:)     :: cs_y_p, sn_y_p
   real(4), Save, allocatable, dimension(:,:,:) :: volume
   real(4), Save                                :: dlambda, dphi
   real(4), Save                                :: Lx, Ly, xo, yo

   logical, Save                                :: global

!-------------------------------------------------------------------------------------------------------------

   Contains

!-------------------------------------------------------------------------------------------------------------

   subroutine sl_init(Nx, Ny, Nz, ModelLevels, xmax, ymax, y0, dx, dy, global_flag)

   implicit none

   integer, intent(in) :: Nx, Ny, Nz
   real(4), intent(in) :: ModelLevels(:)
   real(4), intent(in) :: xmax, ymax, y0
   real(4), intent(in) :: dx, dy
   logical, intent(in) :: global_flag

   real(4) :: x, y
   integer :: i, j, k
   real(4) :: temp1

! Initialise array dimensions used for interpoloation

!  Ni = number of points East-West (Ni+1 equiv 1 for periodicity)
!  Nj = number of points North-South 1 is the South pole and Nj
!                                      is the North pole.
!  Nk = number of model levels
!  Nt = number of tracers - no longer global (module) variable
! 

   Ni = Nx
   Nj = Ny
   Nk = Nz

   Lx = xmax
   Ly = ymax
   dlambda = dx
   dphi    = dy
   xo = 0.0
   yo = y0

   i_p_start = 1
   i_p_end   = Ni
   j_p_start = 1
   j_p_end   = Nj

! halo size >= 2 for cubic lagrange

   halo = 2

   istrt = 1 - halo
   iend = Ni + halo
   jstrt = 1 - halo
   jend = Nj + halo

! Global flag
  global = global_flag  

! set Vertical levels for "p" points and grid spacing

  allocate(eta_levels(Nk), deta(Nk))

  eta_levels(:) = ModelLevels(:)
  deta(1) = eta_levels(1)
  do i = 2, Nk
    deta(i) = eta_levels(i) - eta_levels(i-1)
  end do

  allocate(deta_ratio_m(Nk),deta_ratio_p(Nk))

   do k = 2, Nk-2
     temp1 = eta_levels(k+1) - eta_levels(k)
     deta_ratio_p(k) = ( eta_levels(k+2) - eta_levels(k+1) ) / temp1
     deta_ratio_m(k) = ( eta_levels(k  ) - eta_levels(k-1) ) / temp1
   enddo
   deta_ratio_p(1) = 1.0
   deta_ratio_m(1) = 1.0
   deta_ratio_p(Nk-1:Nk) = 1.0
   deta_ratio_m(Nk-1:Nk) = 1.0

! Set up array of grid points
! i.e. x_pnts(i) = (i-1)*dlambda

   allocate( x_pnts(Ni), y_pnts(Nj) )
   allocate( cs_y_p(Nj), sn_y_p(Nj) )
   allocate( volume(Ni,Nj,Nk) )

   do j = 1, Nj
      y_pnts(j) = yo + (j-1)*dphi
      cs_y_p(j) = cos(y_pnts(j))
      sn_y_p(j) = sin(y_pnts(j))
   enddo

   do i = 1, Ni
      x_pnts(i) = xo + (i-1)*dlambda
   enddo

   !
   ! set the volumes around the rho-points
   ! correctly (I am just puting something
   !
   do i = 1, Ni
      do j = 1, Nj
         do k =1, Nk
             volume(i,j,k) = 1.0
            !volume(i,j,k) = dlambda*dphi* r**2 * dr
         enddo
      enddo
   enddo
   end subroutine sl_init

!-------------------------------------------------------------------------------------------------------------

   subroutine sl_advect(fld_n, fld_np1, rho_n, rho_np1, source_np1,  &
                        u_n, v_n, w_n, u_np1, v_np1, w_np1, timestep, Nt)
   !
   ! this routine solves the advection equation:
   !
   ! Df/Dt = S1(n) + S2(n+1)
   ! f(n+1) - f(n)_d = dt*S1(n)_d + dt*S2(n+1)
   ! f(n+1) = [f(n)+dt*S1(n)]_d + dt*S2(n+1)
   ! usually fld_n = f(n) + dt*S1(n) and source_np1 = dt*S2(n+1)
   !
   ! On input fld_n is the field(s) to be advected
   ! u_n, v_n are the collocated angular velocities (u/r v/r) and w_n is etadot
   !          (also collocated to a "p" point)
   ! u_np1, v_np1, w_np1 are the similarly defined values at time level n+1
   !
   ! on output fld_np1 contains the result after advection.

   implicit none
   real(Con), dimension(istrt:,jstrt:,:,:), intent(IN)    :: fld_n
   real(4), dimension(istrt:,jstrt:,:,:), intent(IN)    :: source_np1
   real(4), dimension(istrt:,jstrt:,:  ), intent(IN)    :: rho_n
   real(4), dimension(istrt:,jstrt:,:  ), intent(INOUT) :: rho_np1
   real(4), dimension(istrt:,jstrt:,:  ), intent(IN)    :: u_n, v_n, w_n
   real(4), dimension(istrt:,jstrt:,:  ), intent(IN)    :: u_np1, v_np1, w_np1
   real(Con), dimension(istrt:,jstrt:,:,:), intent(OUT)   :: fld_np1
   real(4),                                        intent(IN)    :: timestep
   integer,                                        intent(IN)    :: Nt

   integer,   dimension(Ni,Nj,Nk) :: i_pd, j_pd, k_pd
   real(4), dimension(Ni,Nj,Nk) :: wght_px, wght_py, wght_pz

   REAL(4), DIMENSION(Nt)                   :: fld_min
   INTEGER, DIMENSION(Nt)                     :: fld_switches
   real(4), allocatable, dimension(:,:,:)   :: f_temp1
   real(4), allocatable, dimension(:,:,:,:) :: f_temp2, f_temp3, f_temp4

   real(4) :: rho_min, rho_min2(1)
   integer   :: kk, ii(1)
   integer   :: use_mass_conservation_tracers = 1
   integer   :: use_mass_conservation_rho = 0

   fld_min = 0.0    ! set minimum value required for fld_np1(:,:,:,kk)
                    ! here default zero for all fld_np1(:,:,:,kk=1:nt)
   rho_min = 0.0    ! min value for dry density

!-------------------------------------------------------------------
! 1.  compute departure of rho-points (i.e., centre of boxes)
!-------------------------------------------------------------------

   call depart_pnt(i_pd, j_pd, k_pd, wght_px, wght_py, wght_pz,                &
                   u_n, v_n, w_n, u_np1, v_np1, w_np1, timestep)

!-------------------------------------------------------------------
! 2.  interpolate fld_n to obtain fld_np1 = (fld_n)_d
!-------------------------------------------------------------------

   call cubic_lagrange(fld_n, fld_np1, i_pd, j_pd, k_pd,                       &
                       wght_px, wght_py, wght_pz, Nt)

   ! add source terms [ note source_np1 is dt*S2(n+1) ]

   fld_np1 = fld_np1 + source_np1

!-------------------------------------------------------------------
! 3.  enforce conservation
!-------------------------------------------------------------------

   If ( use_mass_conservation_rho == 1 ) then  ! enforce mass conservation of dry-density

      allocate( f_temp1(istrt:iend,jstrt:jend,Nk  ),  &
                f_temp2(istrt:iend,jstrt:jend,Nk,1),  &
                f_temp3(istrt:iend,jstrt:jend,Nk,1),  &
                f_temp4(istrt:iend,jstrt:jend,Nk,1)   )

      f_temp1 = 1.0
      f_temp2(:,:,:,1) = rho_n
      f_temp3(:,:,:,1) = rho_np1
      f_temp4(:,:,:,1) = 0.0
      rho_min2(1) = rho_min
      ii(1) = 1

  !    Next call commented out because not used with use_mass_conservation_rho = 0
  !    and hard to fix up to deal with possible precision changes of the concentration arrays
  !    (as controlled by "Con"). Note rho_np1 isn't changed by sl_advect so is in practice 
  !    Intent(In).
  !    call mass_conservation(f_temp1,f_temp1,f_temp2,f_temp3,f_temp4,rho_min2,ii,1)
      rho_np1 = f_temp3(:,:,:,1)

      deallocate(f_temp1, f_temp2, f_temp3, f_temp4)

   Else  ! just make sure  that  rho_np1 >= rho_min
         ! assumption is already   rho_n >= rho_min

      rho_np1 = max(rho_np1, rho_min)

   Endif

   IF ( use_mass_conservation_tracers == 1 ) THEN

       fld_switches = 1  ! 1 switch on or 0 switch off conservation option for
                         ! fld_np1(:,:,:,kk)
                         ! here it is set to on for all fld_np1(:,:,:,kk=1:nt)

       call mass_conservation(rho_n, rho_np1, fld_n, fld_np1, source_np1,   &
                              fld_min, fld_switches,  Nt                    )

    ELSE ! simply insure fld_np1 >= fld_min

      DO kk = 1, Nt
         fld_np1(:,:,:,kk) = max(fld_np1(:,:,:,kk),fld_min(kk))
      END DO

   ENDIF

!-------------------------------------------------------------------
! 4.  fill halo at the end of step
!-------------------------------------------------------------------

   if ( global ) call fill_haloes(fld_np1,.false.)

 end subroutine sl_advect

!-------------------------------------------------------------------------------------------------------------

   subroutine cubic_lagrange(flds_in, flds_out,                                &
                             i_pd, j_pd, k_pd,                                 &
                             wght_px, wght_py, wght_pz, Nt)

   implicit none

   real(Con), dimension(istrt:,jstrt:,:,:) :: flds_in
   real(Con), dimension(istrt:,jstrt:,:,:) :: flds_out
   integer, dimension(Ni,Nj,Nk)          :: i_pd, j_pd, k_pd
   real(4), dimension(Ni,Nj,Nk)          :: wght_px, wght_py, wght_pz
   integer                               :: Nt

   real(4),    parameter                 :: eps = 1.e-10
   real(4),    parameter                 :: one_6th = 1.0/6.0
   real(4),    parameter                 :: half = 0.5
   real(4),    parameter                 :: one_minus_eps = 1.0-eps
   integer                               :: i, j, k, ii, jj, kk
   integer                               :: ipd, jpd, kpd
   real(4),    dimension(-1:2)           :: Ax, Ay, Az
   integer                               :: k_start, k_end
   real(4)                               :: Sx, Sy, Sz
   real(4),    dimension(Nt,-1:2)        :: fi, fim1, fip1, fip2
   real(4),    dimension(Nt)             :: fk, fkm1, fkp1, fkp2
   real(4)                               :: Szm, Szp

!! !$omp parallel do private(i,j, sx,sy,sz, Ax, Ay, Az, fim1,fi,fip1,fip2,            &
!! !$omp& jkm1, fk, fkp1, fkp2, ii,  k_start,k_end)
   do k = 1, Nk
      do j = j_p_start, j_p_end
         do i = i_p_start, i_p_end

! We can now perform all the required interpolations

            sx = wght_px(i,j,k)
            sy = wght_py(i,j,k)
            sz = wght_pz(i,j,k)

            ipd = i_pd(i,j,k)
            jpd = j_pd(i,j,k)
            kpd = k_pd(i,j,k)

! First calculate interpolation coefs. (uniform mesh)

            Ax(-1) = -Sx*(Sx-1.0)*(Sx-2.0)      * one_6th
            Ax(0)  = (Sx+1.0)*(Sx-1.0)*(Sx-2.0) * half
            Ax(1)  = -Sx*(Sx+1.0)*(Sx-2.0)      * half
            Ax(2)  =  Sx*(Sx+1.0)*(Sx-1.0)      * one_6th

            Ay(-1) = -Sy*(Sy-1.0)*(Sy-2.0)      * one_6th
            Ay(0)  = (Sy+1.0)*(Sy-1.0)*(Sy-2.0) * half
            Ay(1)  = -Sy*(Sy+1.0)*(Sy-2.0)      * half
            Ay(2)  =  Sy*(Sy+1.0)*(Sy-1.0)      * one_6th

            if( kpd == 1 ) then
                  Az(-1)     = 0.0
                  Az(2)      = 0.0
                  Az(0)      = 1.0 - Sz
                  Az(1)      = Sz

                  fim1(:,-1) = 0.0
                  fi(:,-1)   = 0.0
                  fip1(:,-1) = 0.0
                  fip2(:,-1) = 0.0
                  fim1(:,2)  = 0.0
                  fi(:,2)    = 0.0
                  fip1(:,2)  = 0.0
                  fip2(:,2)  = 0.0
                  k_start    = 0
                  k_end      = 1
               elseif( kpd == Nk ) then
                  Az(-1)     = Sz
                  Az(0)      = 1.0 - Sz
                  Az(1)      = 0.0
                  Az(2)      = 0.0

                  fim1(:,1)  = 0.0
                  fi(:,1)    = 0.0
                  fip1(:,1)  = 0.0
                  fip2(:,1)  = 0.0
                  fim1(:,2)  = 0.0
                  fi(:,2)    = 0.0
                  fip1(:,2)  = 0.0
                  fip2(:,2)  = 0.0
                  k_start    = -1
                  k_end      = 0
               elseif( kpd == Nk-1 ) then
                  Az(-1)     = 0.
                  Az(0)      = 1.0 - Sz
                  Az(1)      = Sz
                  Az(2)      = 0.0

                  fim1(:,1)  = 0.0
                  fi(:,1)    = 0.0
                  fip1(:,1)  = 0.0
                  fip2(:,1)  = 0.0
                  fim1(:,2)  = 0.0
                  fi(:,2)    = 0.0
                  fip1(:,2)  = 0.0
                  fip2(:,2)  = 0.0
                  k_start    = 0
                  k_end      = 1
               else

                  ! this for equal vertical mesh

                  !Az(-1)  =- Sz*(Sz-1.0)*(Sz-2.0)      * one_6th
                  !Az(0)   = (Sz+1.0)*(Sz-1.0)*(Sz-2.0) * half
                  !Az(1)   =- Sz*(Sz+1.0)*(Sz-2.0)      * half
                  !Az(2)   = Sz*(Sz+1.0)*(Sz-1.0)       * one_6th

                  ! this for variable vertical mesh

                 Szm = deta_ratio_m(k)
                 Szp = deta_ratio_p(k)

                 Az(-1)  = Sz*(Sz-1.0)*(1.0+Szp-Sz) / ( Szm*(Szm+1.0)*(1.0+Szp+Szm) )
                 Az(0)   = -(Szm+Sz)*(Sz-1.0)*(1.0+Szp-Sz) / ( Szm*(1.0+Szp) )
                 Az(1)   =  (Szm+Sz)*Sz*(1.0+Szp-Sz) / ( Szp*(1.0+Szm) )
                 Az(2)   = (Szm+Sz)*Sz*(Sz-1.0) / ( Szp*(Szp+1.0)*(1.0+Szp+Szm) )

                 k_start = -1
                 k_end   = 2

               endif

            do ii = k_start, k_end

               fim1(:,ii) = Ay(-1)*flds_in(ipd-1,jpd-1,kpd+ii,:) +                   &
                            Ay(0) *flds_in(ipd-1,jpd  ,kpd+ii,:) +                   &
                            Ay(1) *flds_in(ipd-1,jpd+1,kpd+ii,:) +                   &
                            Ay(2) *flds_in(ipd-1,jpd+2,kpd+ii,:)

               fi(:,ii)   = Ay(-1)*flds_in(ipd,jpd-1,kpd+ii,:)   +                   &
                            Ay(0) *flds_in(ipd,jpd  ,kpd+ii,:)   +                   &
                            Ay(1) *flds_in(ipd,jpd+1,kpd+ii,:)   +                   &
                            Ay(2) *flds_in(ipd,jpd+2,kpd+ii,:)

               fip1(:,ii) = Ay(-1)*flds_in(ipd+1,jpd-1,kpd+ii,:) +                   &
                            Ay(0) *flds_in(ipd+1,jpd  ,kpd+ii,:) +                   &
                            Ay(1) *flds_in(ipd+1,jpd+1,kpd+ii,:) +                   &
                            Ay(2) *flds_in(ipd+1,jpd+2,kpd+ii,:)

               fip2(:,ii) = Ay(-1)*flds_in(ipd+2,jpd-1,kpd+ii,:) +                   &
                            Ay(0) *flds_in(ipd+2,jpd  ,kpd+ii,:) +                   &
                            Ay(1) *flds_in(ipd+2,jpd+1,kpd+ii,:) +                   &
                            Ay(2) *flds_in(ipd+2,jpd+2,kpd+ii,:)

               enddo

            fkm1(:) = Ax(-1)*fim1(:,-1) + Ax(0)*fi(:,-1) + Ax(1)*fip1(:,-1)    &
                                        + Ax(2)*fip2(:,-1)
            fk(:)   = Ax(-1)*fim1(:,0)  + Ax(0)*fi(:,0)  + Ax(1)*fip1(:,0)     &
                                        + Ax(2)*fip2(:,0)
            fkp1(:) = Ax(-1)*fim1(:,1)  + Ax(0)*fi(:,1)  + Ax(1)*fip1(:,1)     &
                                        + Ax(2)*fip2(:,1)
            fkp2(:) = Ax(-1)*fim1(:,2)  + Ax(0)*fi(:,2)  + Ax(1)*fip1(:,2)     &
                                        + Ax(2)*fip2(:,2)

            flds_out(i,j,k,:) = Az(-1)*fkm1 + Az(0)*fk + Az(1)*fkp1            &
                                            + Az(2)*fkp2

            enddo
         enddo
      enddo
!! !$omp end parallel do

   end subroutine cubic_lagrange

!-------------------------------------------------------------------------------------------------------------

   subroutine depart_pnt(i_pd, j_pd, k_pd,                                     &
                         wght_px, wght_py, wght_pz,                            &
                         u_n, v_n, w_n, u_np1, v_np1, w_np1, timestep)

   implicit none
   integer, parameter                                        :: Num_its = 1
   integer, dimension(Ni,Nj,Nk),                 intent(OUT) :: i_pd, j_pd, k_pd
   real(4), dimension(Ni,Nj,Nk),                 intent(OUT) :: wght_px, wght_py, wght_pz
   real(4), dimension(istrt:iend,jstrt:jend,Nk), intent(IN)  :: u_n, v_n, w_n,            &
                                                                 u_np1, v_np1, w_np1
   real(4),                                      intent(IN)  :: timestep

   real(4),    dimension(Ni,Nj,Nk)                   :: x_pd, y_pd, z_pd
   real(4),    dimension(Ni,Nj,Nk)                   :: u_int_d, v_int_d, w_int_d
   real(4)                                           :: Sx, Sy, Sz
   real(4)                                           :: x_d, y_d, z_d
   real(4)                                           :: u_a, v_a, w_a
   real(4)                                           :: u_d, v_d, w_d
   real(4)                                           :: rm_11, rm_12
   real(4)                                           :: snx_d, csx_d, sny_d, csy_d
   real(4)                                           :: gam, q_33, Csalpha
   integer                                           :: i, j, k, it

! Initial guess at departure point

!! !$omp parallel do private(i,j,u_a,v_a,w_a,x_d,y_d,z_d)
         do k = 1, Nk
            do j = j_p_start, j_p_end
               do i = i_p_start, i_p_end

! Arrival points

                  u_a    = u_np1(i,j,k)
                  v_a    = v_np1(i,j,k)
                  w_a    = w_np1(i,j,k)

! Update departure points

                  x_d = -timestep*u_a
                  y_d = -timestep*v_a
                  z_d = sqrt(1.0 - x_d**2 - y_d**2)

                  x_pd(i,j,k) = x_pnts(i) + ATAN2(x_d,z_d*cs_y_p(j)-y_d*sn_y_p(j))
                  y_pd(i,j,k) = ASIN(y_d*cs_y_p(j)+z_d*sn_y_p(j))
                  z_pd(i,j,k) = eta_levels(k) - timestep*w_a
                  z_pd(i,j,k) = max(eta_levels(1),min(eta_levels(Nk),z_pd(i,j,k)))

               enddo
            enddo
         enddo
!! !$omp end parallel do

   call locate_dps(x_pd,y_pd,z_pd, i_pd, j_pd, k_pd, wght_px, wght_py, wght_pz)

   do it = 1, Num_its

      call interpolate_wind(u_int_d, v_int_d, w_int_d, u_n, v_n, w_n,          &
                            i_pd, j_pd, k_pd, wght_px, wght_py, wght_pz)

!! !$omp parallel do
!! !private(i,j,u_a,v_a,w_a,u_d,v_d,w_d,x_d,y_d,z_d,rm_11,rm_12, &
!! !$omp& gam,q_33,Csalpha,snx_d,csx_d,sny_d,csy_d)
         do k = 1, Nk
            do j = j_p_start, j_p_end
               do i = i_p_start, i_p_end

                  snx_d = SIN(x_pd(i,j,k)-x_pnts(i))
                  csx_d = COS(x_pd(i,j,k)-x_pnts(i))
                  sny_d = SIN(y_pd(i,j,k))
                  csy_d = COS(y_pd(i,j,k))

                  q_33    = sny_d*sn_y_p(j)+csy_d*cs_y_p(j)*csx_d
                  gam     = (1.0 + q_33)
                  Csalpha = 1.0/gam
                  gam     = 0.5*gam*timestep

! Components of rotation matrix
                  rm_11 = (csy_d*cs_y_p(j) + (1.0+sn_y_p(j)*sny_d)             &
                          *csx_d)*Csalpha
                  rm_12 =-snx_d*(sn_y_p(j) + sny_d)*Csalpha

! Arrival points

                  u_a    = u_np1(i,j,k)
                  v_a    = v_np1(i,j,k)
                  w_a    = w_np1(i,j,k)

! Rotated velocities at departure points

                  u_d = rm_11*u_int_d(i,j,k) + rm_12*v_int_d(i,j,k)
                  v_d =-rm_12*u_int_d(i,j,k) + rm_11*v_int_d(i,j,k)
                  w_d = w_int_d(i,j,k)

! Update departure points

                  x_d = -0.5*gam*(u_a + u_d)
                  y_d = -0.5*gam*(v_a + v_d)
                  z_d = sqrt(1.0 - x_d**2 - y_d**2)

                  x_pd(i,j,k) = x_pnts(i) + ATAN2(x_d,z_d*cs_y_p(j)-y_d*sn_y_p(j))
                  y_pd(i,j,k) = ASIN(y_d*cs_y_p(j)+z_d*sn_y_p(j))
                  z_pd(i,j,k) = eta_levels(k) - 0.5*timestep*(w_a + w_d)
                  z_pd(i,j,k) = max(eta_levels(1),min(eta_levels(Nk),z_pd(i,j,k)))
                  
               enddo
            enddo
         enddo
!! !$omp end parallel do

    call locate_dps(x_pd,y_pd,z_pd, i_pd, j_pd, k_pd, wght_px, wght_py, wght_pz)

    enddo

   end subroutine depart_pnt

!-------------------------------------------------------------------------------------------------------------

   subroutine fill_haloes(fld,flag)

   implicit none
   real(Con), dimension(istrt:,jstrt:,:,:) :: fld
   logical                               :: flag
   real(4)                               :: isgn
   integer                               :: i, j, k

   isgn = 1.0
   if( flag) isgn = -1.0

   fld(Ni+1:Ni+halo,:,:,:) = fld(1:halo,:,:,:)
   fld(1-halo:0,:,:,:)     = fld(Ni+1-halo:Ni,:,:,:)

   k = Ni/2 - halo + 1
   do i = 1-halo, Ni+halo
      do j = 1, halo
         fld(i,Nj+j,:,:) = isgn*fld(k,Nj-j,:,:)
         fld(i,1-j,:,:)  = isgn*fld(k,j+1,:,:)
         enddo
      k = k + 1
      if( k > Ni ) k = 1
      enddo

   end subroutine fill_haloes

!-------------------------------------------------------------------------------------------------------------

   subroutine interpolate_wind(u_d,v_d,w_d, u,v,w,                             &
                               i_d,j_d,k_d, sx,sy,sz)

   implicit none
   real(4), dimension(Ni,Nj,Nk),  intent(OUT)   :: u_d, v_d, w_d
   real(4), dimension(Ni,Nj,Nk),  intent(IN)    :: sx, sy, sz
   real(4), dimension(istrt:iend,jstrt:jend,Nk),                               &
                                  intent(IN)    :: u, v, w
   integer, dimension(Ni,Nj,Nk),  intent(IN)    :: i_d, j_d, k_d

   real(4)                                    :: r0, s0, t0
   real(4)                                    :: r1, s1, t1
   integer                                      :: ii, jj, kk
   integer                                      :: i, j, k

!! !$omp parallel do private(ii,jj, i,j,k, r0,r1, s0,s1, t0,t1)
   do kk = 1, Nk
      do jj = j_p_start, j_p_end
         do ii = i_p_start, i_p_end

            i = i_d(ii,jj,kk)
            j = j_d(ii,jj,kk)
            k = k_d(ii,jj,kk)

            r0 = sx(ii,jj,kk)
            s0 = sy(ii,jj,kk)
            t0 = sz(ii,jj,kk)

            r1 = 1.0 - r0
            s1 = 1.0 - s0
            t1 = 1.0 - t0

            u_d(ii,jj,kk) = t0*( s0*( r0*u(i,j,k)     + r1*u(i+1,j,k) )           &
                                +s1*( r0*u(i,j+1,k)   + r1*u(i+1,j+1,k) ) )       &
                           +t1*( s0*( r0*u(i,j,k+1)   + r1*u(i+1,j,k+1) )         &
                                +s1*( r0*u(i,j+1,k+1) + r1*u(i+1,j+1,k+1) ) )

            v_d(ii,jj,kk) = t0*( s0*( r0*v(i,j,k)     + r1*v(i+1,j,k) )           &
                                +s1*( r0*v(i,j+1,k)   + r1*v(i+1,j+1,k) ) )       &
                           +t1*( s0*( r0*v(i,j,k+1)   + r1*v(i+1,j,k+1) )         &
                                +s1*( r0*v(i,j+1,k+1) + r1*v(i+1,j+1,k+1) ) )

            w_d(ii,jj,kk) = t0*( s0*( r0*w(i,j,k)     + r1*w(i+1,j,k) )           &
                                +s1*( r0*w(i,j+1,k)   + r1*w(i+1,j+1,k) ) )       &
                           +t1*( s0*( r0*w(i,j,k+1)   + r1*w(i+1,j,k+1) )         &
                                +s1*( r0*w(i,j+1,k+1) + r1*w(i+1,j+1,k+1) ) )

            enddo
         enddo
      enddo
!! !$omp end parallel do

   end subroutine interpolate_wind

!-------------------------------------------------------------------------------------------------------------

   subroutine  locate_dps(x_d,y_d,z_d, i_d, j_d, k_d,                          &
                          wght_x, wght_y, wght_z)

   implicit none
   real(4), dimension(Ni,Nj,Nk),    intent(INOUT) :: x_d, y_d, z_d
   real(4), dimension(Ni,Nj,Nk),    intent(OUT)   :: wght_x, wght_y, wght_z
   integer,   dimension(Ni,Nj,Nk),    intent(OUT)   :: i_d, j_d, k_d

   integer                                        :: i, j, k
   integer                                        :: k_lower, k_upper, k_mid

! do horizontal grid first since it's uniform
! NOTE: we don't need to check crossing the pole since
!       the asin in the departure points can only return
!       a value between -pi/2 and +pi/2.

!! !$omp parallel do private(i,j)
   do k = 1, Nk
      do j = j_p_start, j_p_end
         do i = i_p_start, i_p_end

            if ( global ) then
              if ( x_d(i,j,k) < xo     ) x_d(i,j,k) = x_d(i,j,k) + Lx
              if ( x_d(i,j,k) >= xo+Lx ) x_d(i,j,k) = x_d(i,j,k) - Lx
            else
              if ( x_d(i,j,k) < xo     ) x_d(i,j,k) = xo
              if ( x_d(i,j,k) >= xo+Lx ) x_d(i,j,k) = xo + Lx
            end if
            i_d(i,j,k) = 1 + floor( (x_d(i,j,k) - xo)/dlambda )
            if ( i_d(i,j,k) > i_p_end ) i_d(i,j,k) = i_p_end
            if ( i_d(i,j,k) < i_p_start ) i_d(i,j,k) = i_p_start

            if ( .not.global ) then
              if ( y_d(i,j,k) < yo     ) y_d(i,j,k) = yo 
              if ( y_d(i,j,k) >= yo+Ly ) y_d(i,j,k) = yo + Ly
            end if
            j_d(i,j,k) = 1 + floor( (y_d(i,j,k) - yo)/dphi )
            if ( j_d(i,j,k) > j_p_end ) j_d(i,j,k) = j_p_end
            if ( j_d(i,j,k) < j_p_start ) j_d(i,j,k) = j_p_start

            wght_x(i,j,k) = (x_d(i,j,k) - x_pnts(i_d(i,j,k)))/dlambda
            wght_y(i,j,k) = (y_d(i,j,k) - y_pnts(j_d(i,j,k)))/dphi

         enddo
      enddo
   enddo

!! !$omp end parallel do

!! !$omp parallel do private(k_lower,k_upper,k_mid,i,j)
  do k = 1, Nk
    do j = j_p_start, j_p_end
        do i = i_p_start, i_p_end
           k_mid = Nk/2
           if ( z_d(i,j,k) > eta_levels(Nk - 1) ) z_d(i,j,k) = eta_levels(Nk-1)
           if ( z_d(i,j,k) < eta_levels(1) ) z_d(i,j,k) = eta_levels(1)
           if ( z_d(i,j,k) <= eta_levels(k_mid) ) then
              k_lower = 1
              k_upper = k_mid
           else
              k_lower = k_mid
              k_upper = Nk-1
           endif

           do
              k_mid = (k_lower + k_upper)/2
              if( k_mid == k_lower ) exit

              if( z_d(i,j,k) < eta_levels(k_mid) ) THEN
                 k_upper = k_mid
              else
                 k_lower = k_mid
              endIF
           enddo
           k_d(i,j,k) = k_mid
           wght_z(i,j,k) = (z_d(i,j,k) - eta_levels(k_d(i,j,k)))/( eta_levels(k_upper) - eta_levels(k_lower) )
        enddo
      enddo
    enddo
!! !$omp end parallel do

   end subroutine locate_dps

!-------------------------------------------------------------------------------------------------------------

 SUBROUTINE mass_conservation(rho_n, rho_np1, qs_n, qs_np1, qs_s,  &
                              qsmin, qswitches, number_qs          )

! ----------------------------------------------------------------------
!  This routine assumes rho_species = rho * q
!
!  therefore, it should work for both:
!
!            (a) rho = rho_dry & q = mixing_ratios, and
!            (b) rho = rho_wet & q = specific humidity
!
!=========================================================================
!
! The quantity we want to conserve is:
!
! (1) SUM{ volume(i,j,k)*rho_np1(i,j,k)*qs_np1(i,j,k) } = C (known value)
!
! where volume(i,j,k) is the volume of the box surrounding a rho-point
! or rho(i,j,k) is the centre of volume(i,j,k), and,
!
! (2) C = SUM(  volume(i,j,k) *   rho_n(i,j,k) *  qs_n(i,j,k) )  +
!         SUM(  volume(i,j,k) * rho_np1(i,j,k) *  qs_s(i,j,k) )  
!
!     Note that the source_term qs_s are at (n+1) and refers to Physics2
!     the sources due to physics1 is contained within qs_n
!
! ----------------------------------------------------------------------

       implicit none
       INTEGER,                                    INTENT(IN)    :: number_qs
       REAL(4),   DIMENSION(istrt:, jstrt:, :),    INTENT(IN)    :: rho_n, rho_np1
       REAL(Con), DIMENSION(istrt:, jstrt:, :, :), INTENT(IN)    :: qs_n
       REAL(4),   DIMENSION(istrt:, jstrt:, :, :), INTENT(IN)    :: qs_s
       REAL(Con), DIMENSION(istrt:, jstrt:, :, :), INTENT(INOUT) :: qs_np1
       REAL(4),   DIMENSION(:),                    INTENT(IN)    :: qsmin
       INTEGER,   DIMENSION(:),                    INTENT(IN)    :: qswitches

! locals

       REAL(4), DIMENSION(Ni,Nj,Nk)      :: psi_n, psi_np1
       REAL(4), DIMENSION(number_qs,3)   :: global_sums
       REAL(4), DIMENSION(number_qs+1)   :: global_sums2
       REAL(4), DIMENSION(number_qs)     :: qs_max

       REAL(4), DIMENSION(Ni,Nj,Nk,number_qs)                  :: lap
       REAL(4), DIMENSION(istrt:iend,jstrt:jend,Nk,number_qs)  :: qs_lap
       REAL(4), DIMENSION(Ni,Nj,Nk)                            :: lap_rho

       REAL(4)                :: dfdxp,dfdxm,dfdyp,dfdym,dfdzp,dfdzm,delsq,gamma
       REAL(4), DIMENSION(Nk) :: dz2, dz3

       REAL(4)          :: mass_deficit, one_over_total_lap, temp1, temp2
       INTEGER            :: i,j,k,kk, im, ip, jm, jp
       REAL(4)          :: small_tol = 1.0E-30
       Integer            :: filter_data = 0

! ----------------------------------------------------------------------
! Section 1 : make qs(n+1) > qsmin
! ----------------------------------------------------------------------

      DO kk = 1, number_qs
         qs_np1(:,:,:,kk) = max(qs_np1(:,:,:,kk),qsmin(kk))
      END DO
      qs_lap = qs_np1

! ----------------------------------------------------------------------
! Section 2 : filter data if noisy (optional)
! ----------------------------------------------------------------------

      If ( filter_data == 1)  Then

        ! filter qs_np1 (apply 1-2-1 filter in all directions)
        ! these filtered data are used only in computing "lap" as
        ! a measure of smoothness -- this is to avoid giving too much
        ! weight for noisy data with insignificant magnitudes

      DO k = 1, Nk
        DO j = jstrt, jend
           DO i = istrt+1, iend-1
              qs_lap(i,j,k,:) = 0.50 * qs_lap(i  ,j,k,:) + &
                                0.25 * qs_lap(i-1,j,k,:) + &
                                0.25 * qs_lap(i+1,j,k,:)
           END DO
         END DO
      END DO
      DO k = 1, Nk
        DO j = jstrt+1, jend-1
           DO i = istrt, iend
              qs_lap(i,j,k,:) = 0.50 * qs_lap(i,j  ,k,:) + &
                                0.25 * qs_lap(i,j-1,k,:) + &
                                0.25 * qs_lap(i,j+1,k,:)
           END DO
         END DO
      END DO
      DO k = 1+1, Nk-1
        DO j = jstrt, jend
           DO i = istrt, iend
              qs_lap(i,j,k,:) = 0.50 * qs_lap(i,j,k  ,:) + &
                                0.25 * qs_lap(i,j,k-1,:) + &
                                0.25 * qs_lap(i,j,k+1,:)
           END DO
         END DO
      END DO

      Endif

!-------------------------------------------------------------
! section 3: Compute the 3D array  psi = rho * volume
!-------------------------------------------------------------

      DO k = 1, Nk
        DO j = 1, Nj
          DO i = 1, Ni
               psi_n(i,j,k)   =   rho_n(i,j,k) * volume(i,j,k)
               psi_np1(i,j,k) = rho_np1(i,j,k) * volume(i,j,k)
            END DO
         END DO
       END DO

! -----------------------------------------------------------------------------
!  section 4:  Calculate the laplacian (lap) of qs. This lap is computed like
!              Cartesian (and this consistent with interpolation).
!              Here constant mesh in x and y directions (horizontally) is assumed
! -----------------------------------------------------------------------------

       DO k = 2 , Nk - 1
          dz2(k) = 0.5 * ( eta_levels(k+1) - eta_levels(k-1) )
       END DO
       DO k = 1 , Nk - 1
          dz3(k) = eta_levels(k+1) - eta_levels(k)
       END DO

       lap = 0.0

       DO kk = 1, number_qs
        IF ( qswitches(kk) == 1 ) THEN

         k = 1
         DO j = 1, Nj
            DO i = 1, Ni

               dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))
               dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))

               dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))
               dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))

               delsq  = (dfdxp-dfdxm) +  (dfdyp-dfdym)

               lap(i,j,k,kk) = abs(delsq) * abs(qs_np1(i,j,k,kk)) * psi_np1(i,j,k)

           END DO
         END DO

         DO k = 2, Nk - 1
           DO j = 1, Nj
             DO i = 1, Ni

               dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))
               dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))

               dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))
               dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))

               dfdzp  = (qs_lap(i,j,k+1,kk) - qs_lap(i,j,k  ,kk))/dz3(k  )
               dfdzm  = (qs_lap(i,j,k  ,kk) - qs_lap(i,j,k-1,kk))/dz3(k-1)

               delsq  = (dfdxp-dfdxm) + (dfdyp-dfdym) +  &
                        (dz3(k)*dz3(k-1)) * (dfdzp-dfdzm)/dz2(k)

               lap(i,j,k,kk) = abs(delsq) * abs(qs_np1(i,j,k,kk)) * psi_np1(i,j,k)

             END DO
           END DO
         END DO


         k = Nk
         DO j = 1, Nj
             DO i = 1, Ni

               dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))
               dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))

               dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))
               dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))

               delsq  = (dfdxp-dfdxm) + (dfdyp-dfdym)

               lap(i,j,k,kk) = abs(delsq) * abs(qs_np1(i,j,k,kk)) * psi_np1(i,j,k)

           END DO
         END DO

        END IF
       END DO

! ----------------------------------------------------------------------
!  section 5:  Calculate the global sums in (1) & (2) and total lap
! ----------------------------------------------------------------------

       DO kk = 1, number_qs

         global_sums(kk,1) = 0.
         global_sums(kk,2) = 0.
         global_sums(kk,3) = 0.

         IF ( qswitches(kk) == 1 ) THEN
           DO k = 1, Nk
             DO j = 1, Nj
                 DO i = 1, Ni

                 global_sums(kk,1) =  global_sums(kk,1)                &
                                   +    psi_n(i,j,k) * qs_n(i,j,k,kk)  &
                                   +  psi_np1(i,j,k) * qs_s(i,j,k,kk)

                 global_sums(kk,2) =  global_sums(kk,2)                &
                                   +  psi_np1(i,j,k) * qs_np1(i,j,k,kk)

                 global_sums(kk,3) =  global_sums(kk,3) +  lap(i,j,k,kk)
               END DO
             END DO
           END DO
         END IF

       END DO

! ----------------------------------------------------------------------
!  section 6:  Modify qs to achieve mass conservation
! ----------------------------------------------------------------------

       DO kk = 1, number_qs
         IF ( qswitches(kk) == 1 .AND. abs(global_sums(kk,3)) > small_tol ) THEN
             mass_deficit = global_sums(kk,1) -  global_sums(kk,2)
             one_over_total_lap = 1.0/global_sums(kk,3)
             DO k = 1, Nk
               DO j = 1, Nj
                 DO i = 1, Ni
                   gamma = lap(i,j,k,kk) * mass_deficit * one_over_total_lap
                   qs_np1(i,j,k,kk) = qs_np1(i,j,k,kk) + gamma/psi_np1(i,j,k)
                 END DO
               END DO
             END DO
         END IF
      END DO

! ----------------------------------------------------------------------
!  section 7: Apply a second correction to force q >= qmin
!             while maintaining mass conservation
! ----------------------------------------------------------------------

       ! clip the field points that are below qsmin

       Do kk = 1, number_qs
         qs_np1(:,:,:,kk) = max(qs_np1(:,:,:,kk),qsmin(kk))
       END DO

       ! compute new masses at (n+1) after cliping

       DO kk = 1, number_qs
         global_sums2(kk) = 0.0
         DO k = 1, Nk
           DO j = 1, Nj
             DO i = 1, Ni
                global_sums2(kk) = global_sums2(kk)             &
                                 + psi_np1(i,j,k) * qs_np1(i,j,k,kk)
             END DO
           END DO
         END DO
       END DO

       ! this sum doesn't invlove qs (it's the sum above with q=1)

       global_sums2(number_qs+1) = 0.0
       DO k = 1, Nk
         DO j = 1, Nj
             DO i = 1, Ni
                global_sums2(number_qs+1) = global_sums2(number_qs+1)  &
                                           + psi_np1(i,j,k)
            END DO
          END DO
       END DO

       DO kk = 1, number_qs
        mass_deficit = global_sums2(kk)  - qsmin(kk)*global_sums2(number_qs+1)
        IF ( qswitches(kk) == 1 .AND. abs(mass_deficit) > small_tol ) THEN
          temp1 = ( global_sums(kk,1) - qsmin(kk)*global_sums2(number_qs+1) ) / &
                    mass_deficit
          temp2  = qsmin(kk)*(1.0 - temp1)
          qs_np1(:,:,:,kk) = temp1 * qs_np1(:,:,:,kk) + temp2
        END IF
       END DO

      RETURN
END SUBROUTINE mass_conservation

!-------------------------------------------------------------------------------------------------------------

End Module SemiLagrangianModule

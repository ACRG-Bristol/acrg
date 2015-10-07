Module transd_inv

contains


SUBROUTINE hbtdmcmc(beta,k, x, h_agg,y, z, plon, plat, regions_v, &
pdf_param1, pdf_param2, lon,lat, h_v, sigma_yt, sigma_ys, sigma_measure, &
R_indices, timeindex_nonzero, sigma_clon, sigma_clat, tau, deltatime, distance, &
pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, sigma_yt_hparams, sigma_ys_hparams, tau_hparams, &
stepsize, stepsize_pdf_p1, stepsize_pdf_p2,  stepsize_sigma_yt, stepsize_sigma_ys, stepsize_tau, stepsize_y, &
x_pdf_all, pdf_param1_pdf, pdf_param2_pdf, sigma_yt_pdf, sigma_ys_pdf, tau_pdf, &
rho, nu, rho_hparams, nu_hparams, stepsize_rho, stepsize_nu, rho_pdf, nu_pdf, rjmcmc, &
lonmin, lonmax, latmin, latmax, kmin, sigma_bd, burn_in, nIt, nmeasuretotal, nsub, nIC, nit_sub,  &
nmeasure, nmeasuremax, ydim1, ydim2, Ngrid, nlon, kICmax, nlat, nbeta, kmax, numsites, nIC1, &
k_out, x_out, regions_out, plon_out, plat_out, n0T_out, y_out, &
sigma_yt_out, sigma_ys_out, tau_out, pdf_param1_out, pdf_param2_out, rho_out, nu_out, &
accept, reject, accept_birth, reject_birth, accept_death, reject_death, &
accept_move, reject_move, accept_swap, reject_swap, accept_sigma_yt, reject_sigma_yt, &
accept_sigma_ys, reject_sigma_ys, accept_tau, reject_tau, accept_rho, reject_rho, accept_y, reject_y)

IMPLICIT NONE
!! INPUTS !!!!!!!!!!!!!!!!!!!
! Dimensions
INTEGER nbeta
INTEGER kmax
INTEGER kICmax
INTEGER nmeasure
INTEGER nmeasuretotal
INTEGER numsites
INTEGER nmeasuremax
INTEGER Ngrid
INTEGER nlon
INTEGER nlat
INTEGER nit_sub
INTEGER nIt
INTEGER burn_in
INTEGER nsub
INTEGER kmin
INTEGER nIC
INTEGER ydim1
INTEGER ydim2
INTEGER nIC1
! Single Variables
REAL lonmin
REAL lonmax
REAL latmin
REAL latmax
REAL sigma_bd
REAL sigma_clon
REAL sigma_clat
REAL stepsize_sigma_yt
REAL stepsize_sigma_ys
REAL stepsize_y
REAL deltatime
REAL stepsize_tau
REAL stepsize_rho
REAL stepsize_nu
INTEGER tau_pdf
INTEGER sigma_yt_pdf
INTEGER sigma_ys_pdf
INTEGER pdf_param1_pdf
INTEGER pdf_param2_pdf
INTEGER rho_pdf
INTEGER nu_pdf
INTEGER rjmcmc
! Input arrays
INTEGER x_pdf_all(nIC1)
REAL stepsize_pdf_p1(nIC1)
REAL stepsize_pdf_p2(nIC1)
REAL stepsize(nIC1)
REAL beta(nbeta)
INTEGER k(nbeta)
REAL z(nmeasure)
REAL x(kICmax, nbeta)               
REAL pdf_param1(nIC1,nbeta)            
REAL pdf_param2(nIC1,nbeta)                
REAL h_agg(nmeasuremax,kICmax,nbeta)
REAL y(nmeasuremax, nbeta) 
REAL sigma_yt(ydim2,nbeta)
REAL sigma_ys(numsites,nbeta)
REAL tau(nbeta)
INTEGER R_indices(ydim1,ydim2)
INTEGER timeindex_nonzero(nmeasure)
REAL sigma_measure(nmeasure)
REAL plon(kmax,nbeta)
REAL plat(kmax,nbeta)
INTEGER regions_v(Ngrid,nbeta)
REAL lon(nlon)
REAL lat(nlat)
REAL h_v(nmeasuremax, Ngrid)
REAL pdf_p2_hparam1(nIC1)
REAL pdf_p2_hparam2(nIC1)
REAL pdf_p1_hparam1(nIC1)
REAL pdf_p1_hparam2(nIC1)
REAL sigma_yt_hparams(2)
REAL sigma_ys_hparams(2)
REAL tau_hparams(2)
REAL rho_hparams(2)
REAL nu_hparams(2)
REAL rho(nbeta)
REAL nu(nbeta)
REAL distance(numsites,numsites)
! Outputs
INTEGER k_out(nit_sub)
REAL x_out(kICmax,nit_sub)              
REAL pdf_param1_out(nIC1,nit_sub), pdf_param2_out(nIC1,nit_sub)    
REAL plon_out(kmax,nit_sub)
REAL plat_out(kmax,nit_sub)
INTEGER regions_out(Ngrid,nit_sub)
REAL sigma_yt_out(ydim2, nit_sub)
REAL sigma_ys_out(numsites,nit_sub)
REAL n0T_out(nit_sub)
REAL tau_out(nit_sub)
REAL rho_out(nit_sub)
REAL nu_out(nit_sub)
REAL y_out(nmeasuremax,nit_sub)
INTEGER accept(nIC1)
INTEGER reject(nIC1) 
INTEGER accept_birth, reject_birth
INTEGER accept_death, reject_death, accept_move, reject_move, accept_swap
INTEGER accept_sigma_yt, reject_sigma_yt, accept_tau, reject_tau
INTEGER accept_sigma_ys, reject_sigma_ys, accept_y, reject_y, reject_swap
INTEGER accept_rho, reject_rho
! INTERMEDIATE VARIABLES
INTEGER it, ibeta, remain_it, pair1,pair2, ib, it_sub, remain, kIC, ii,jj       !remain_dim
INTEGER remain_swap, ti, remain_hparam
REAL u, u1,u2, u3, randomu,pT_chain, beta1,beta2, q_small
INTEGER k_it(nit_sub)
REAL x_it(kICmax,nit_sub)                             
REAL pdf_param1_it(nIC1,nit_sub), pdf_param2_it(nIC1,nit_sub)     
REAL plon_it(kmax,nit_sub)               
REAL plat_it(kmax,nit_sub)
INTEGER regions_it(Ngrid,nit_sub)
REAL sigma_yt_it(ydim2,nit_sub)
REAL n0T_it(nit_sub), tau_it(nit_sub), sigma_ys_it(numsites,nit_sub)
REAL rho_it(nit_sub), nu_it(nit_sub)
REAL y_it(nmeasuremax,nit_sub)           
REAL detval(nbeta), detval_T(nbeta), detval_S(nbeta), detval_U(nbeta)
REAL Rinv(nmeasuremax,nmeasuremax,nbeta)
REAL Tinv(nmeasuretotal,nmeasuretotal,nbeta)
REAL Sinv(numsites,numsites,nbeta)
REAL Qinv(nmeasuretotal,nmeasuretotal,nbeta)
REAL Uinv(numsites,numsites,nbeta)
REAL U_temp(numsites,numsites)
REAL Rinv_temp(nmeasuremax,nmeasuremax)
REAL Sinv_temp(numsites,numsites)
REAL Tinv_temp(nmeasuretotal,nmeasuretotal)
REAL n0(nmeasuremax, nbeta), m0(nmeasure) 
REAL n0T(nbeta), m0T(nbeta)
REAL y_small(nmeasure)
REAL sigma_ys_ap(numsites), sigma_yt_ap(ydim2)
REAL y_vart(nmeasuretotal), y_vart_inv(nmeasuretotal), autocorr_vec(nmeasuretotal)
REAL y_vars_inv(numsites), y_vars(numsites)
REAL ga
real (kind =8)  alpha, arg_dum, k_arg_dum
real (kind =8) arg(numsites,numsites)
integer (kind = 4) nb, ize, ncalc
! SUBROUTINE INPUTS
REAL betaib, n0Tib, detvalib, tauib, m0Tib, detval_Sib, detval_Tib
REAL rhoib, nuib, detval_Uib
REAL xib(kICmax), plonib(kmax), platib(kmax), n0ib(nmeasuremax)           
REAL pdf_param1ib(nIC1), pdf_param2ib(nIC1), yib(nmeasuremax)
REAL h_aggib(nmeasuremax,kICmax)
REAL sigma_ysib(numsites), sigma_ytib(ydim2)
INTEGER kib
INTEGER regions_vib(Ngrid)
REAL Tinvib(nmeasuretotal,nmeasuretotal), Rinvib(nmeasuremax,nmeasuremax)
REAL Qinvib(nmeasuretotal,nmeasuretotal), Sinvib(numsites,numsites), Uinvib(numsites,numsites)
! SUBROUTINE OUTPUTS
REAL n0Tib1, detvalib1, tauib1, m0Tib1, detval_Sib1, detval_Tib1
REAL rhoib1, nuib1, detval_Uib1                     
REAL xib1(kICmax), plonib1(kmax), platib1(kmax), n0ib1(nmeasuremax)      
REAL pdf_param1ib1(nIC1), pdf_param2ib1(nIC1)
REAL h_aggib1(nmeasuremax,kICmax), yib1(nmeasuremax)
REAL sigma_ysib1(numsites), sigma_ytib1(ydim2)
INTEGER kib1, rejectib1, acceptib1, acceptyib1, rejectyib1
INTEGER acceptxib1(nIC1), rejectxib1(nIC1)
INTEGER regions_vib1(Ngrid)
REAL Tinvib1(nmeasuretotal, nmeasuretotal), Rinvib1(nmeasuremax,nmeasuremax)
REAL Qinvib1(nmeasuretotal, nmeasuretotal), Sinvib1(numsites,numsites), Uinvib1(numsites,numsites)
! OTHER INTERMEDIATES


!f2py intent(in) beta,k, x, h_agg,y, z, plon, plat, regions_v 
!f2py intent(in) pdf_param1, pdf_param2, lon,lat, h_v, sigma_yt, sigma_ys, sigma_measure
!f2py intent(in) R_indices, timeindex_nonzero, sigma_clon, sigma_clat, tau, deltatime, distance
!f2py intent(in) pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2
!f2py intent(in) sigma_yt_hparams, sigma_ys_hparams, tau_hparams, stepsize, stepsize_pdf_p1
!f2py intent(in) stepsize_pdf_p2, stepsize_sigma_yt, stepsize_sigma_ys, stepsize_tau, stepsize_y
!f2py intent(in) x_pdf_all, pdf_param1_pdf, pdf_param2_pdf, sigma_yt_pdf, sigma_ys_pdf, tau_pdf
!f2py intent(in) lonmin, lonmax, latmin, latmax, kmin, sigma_bd, burn_in, nIt, nsub, nit_sub, nIC
!f2py intent(in) nIC1, nit_sub, nIC, nbeta, kmax, kICmax, rjmcmc
!f2py intent(in) nmeasure, nmeasuretotal, nmeasuremax, numsites, Ngrid, nlon,nlat, ydim1, ydim2
!f2py intent(in) nu, rho, nu_hparams, rho_hparams, stepsize_rho, stepsize_nu, rho_pdf, nu_pdf
!f2py intent(out) k_out, x_out, regions_out, plon_out, plat_out, n0T_out, y_out
!f2py intent(out) sigma_yt_out, sigma_ys_out, tau_out, pdf_param1_out, pdf_param2_out
!f2py intent(out) nu_out, rho_out, accept_rho, reject_rho
!f2py intent(out) accept, reject, accept_birth, reject_birth, accept_death, reject_death
!f2py intent(out) accept_move, reject_move, accept_swap, reject_swap, accept_sigma_yt, reject_sigma_yt
!f2py intent(out) accept_sigma_ys, reject_sigma_ys, accept_tau, reject_tau, accept_y, reject_y

!!!!!!call!f2py intent(out) k_out, x_out, regions_out, plon_out, plat_out, n0T_out, y_out
!!!!!!calls!f2py intent(out) sigma_yt_out, sigma_ys_out, tau_out, pdf_param1_out, pdf_param2_out
!!f intent(hide),depend(beta) :: nbeta=shape(beta,0)       !f2py intent(in) nit_sub, nIC, nbeta, kmax, kICmax
!!!f2 intent(hide),depend(plon) :: kmax=shape(plon,0)      !f2py intent(in) nmeasure, nmeasuretotal, nmeasuremax, numsites, Ngrid, nlon,nlat, ydim1, ydim2
!!!!fy intent(hide),depend(x) :: kICmax=shape(x,0)
!!!!fy intent(hide),depend(y) :: nmeasuremax=shape(y,0)
!!!fy intent(hide),depend(z) :: nmeasure=shape(z,0)
!!!fy intent(hide),depend(space_corr) :: numsites=shape(space_corr,0)
!!!f2 integer intent(hide),depend(h_v) :: Ngrid=shape(h_v,1)
!!!f2 integer intent(hide),depend(lon) :: nlon=shape(lon,0)
!!!f2 integer intent(hide),depend(lat) :: nlat=shape(lat,0)
!!!f2 integer intent(hide),depend(R_indices) :: ydim1=shape(R_indices,0), ydim2=shape(R_indices,1)

 ! call OMP_SET_NUM_THREADS(nbeta)

accept=0
accept_birth=0
accept_death=0
accept_move=0
accept_sigma_yt=0
accept_sigma_ys=0
accept_tau=0
accept_rho=0
accept_y=0
accept_swap=0

reject=0
reject_birth=0
reject_death=0
reject_move=0
reject_sigma_yt=0
reject_sigma_ys=0
reject_tau=0
reject_rho=0
reject_y=0
reject_swap=0
it_sub=1

sigma_ys_ap = sigma_ys(:,1)
sigma_yt_ap = sigma_yt(:,1)

 !! call init_random_seed()

! ******** Set up covariance matrices ************************
! Should be same for all beta
Sinv(:,:,:)=0.
Qinv(:,:,:) = 0.
Uinv(:,:,:)=0.

do ibeta=1,nbeta
!  ibeta=1

  do jj=1,ydim2   
      y_vart(R_indices(:,jj)) = sigma_yt(jj,ibeta)    ! Provided ydim2 isn't too big then should be fine
  enddo  

  y_vart_inv = 1./y_vart

  q_small = exp(-1.*deltatime/tau(ibeta))

  !Qinv(:,:,ibeta) = 0.
  Qinv(1,1,ibeta) = 1.
  Qinv(nmeasuretotal,nmeasuretotal,ibeta)=1.
  Qinv(1,2,ibeta) = q_small*(-1.)
  Qinv(nmeasuretotal, nmeasuretotal-1,ibeta)=q_small*(-1.)

  do ii=2, nmeasuretotal-1
          Qinv(ii,ii,ibeta) = 1 + q_small**2
          Qinv(ii,ii+1,ibeta) = q_small*(-1.)
          Qinv(ii,ii-1,ibeta) = q_small*(-1.)
  enddo
               
  Qinv(:,:,ibeta) = Qinv(:,:,ibeta)/(1-q_small**2)

  do ti=1,nmeasuretotal
        autocorr_vec = y_vart_inv(ti)*y_vart_inv*Qinv(:,ti,ibeta)
        Tinv(:,ti,ibeta) = autocorr_vec      
  enddo

  if (numsites .GT. 1) then

     y_vars = sigma_ys(:,ibeta)
     y_vars_inv = 1./sigma_ys(:,ibeta) 

     nb = int(nu(ibeta))
     alpha = nu(ibeta) - nb
     nb = nb+1
     arg = sqrt(2*nu(ibeta))*distance/rho(ibeta)
     ize = 2
     ga= gamma(nu(ibeta))


     do ii = 1,numsites 
         do jj = 1,numsites
              arg_dum = arg(ii,jj)
           !   if (arg_dum==0.) then
           !        U_temp(ii,jj) = 1.
           !   else
           !        call rkbesl(arg_dum, alpha, nb, ize, k_arg_dum, ncalc)
           !        k_arg_dum = k_arg_dum*exp(-1*arg(ii,jj))
           !        U_temp(ii,jj) = 1/(2**(nu(ibeta)-1)*ga)*(sqrt(2*nu(ibeta))*distance(ii,jj)/rho(ibeta))**nu(ibeta)*k_arg_dum
                   U_temp(ii,jj) = exp(-1.*distance(ii,jj)/rho(ibeta))
           !   endif
         enddo
     enddo
 
       
     Uinv(:,:,ibeta) = Ainv(U_temp)

     detval_U(ibeta) = det(U_temp)
 
     do ii=1,numsites
           Sinv(:,ii,ibeta) = y_vars_inv(ii)*y_vars_inv*Uinv(:,ii,ibeta)
           !Sinv(ii,ii,ibeta) = space_vec
         !  Sinv(ii,ii,ibeta) = 1./y_vars(ii)**2      ! Ignore spatital correlation for now

     enddo
   

     Sinv_temp = Sinv(:,:,ibeta)
     Tinv_temp = Tinv(:,:,ibeta)

     call kron(Rinv_temp, Sinv_temp, Tinv_temp, numsites, nmeasuretotal, nmeasuremax)
     !call kron(Rinv_temp, Sinv(:,:,ibeta), Tinv(:,:,ibeta))
     !Rinv_temp(:,:)=1.

     Rinv(:,:,ibeta) = Rinv_temp

     ! This is wrong as currently stands - only the case for uncorrelated spatial component
    ! detval_S(ibeta) = sum(alog(y_vars))*nmeasuretotal 

     detval_S(ibeta) = (sum(alog(y_vars)) + detval_U(ibeta))*nmeasuretotal
   ! detval_S(ibeta)  = det(Sinv_temp)

  else 

     Rinv(:,:,ibeta) = Tinv(:,:,ibeta)
     detval_S(ibeta) = 0.
     detval_U(ibeta) = 0.

  endif

  detval_T(ibeta) = (sum(alog(y_vart))+(alog(1-q_small**2)*((nmeasuretotal-1)/2.)))*numsites   ! DETERMINANT IN LOG SPACE 
  
  detval(ibeta) = detval_S(ibeta) + detval_T(ibeta)

 

  !n0 = matmul(H_agg,x_agg) - y   ! but y = matmul(H_agg,x_agg) so n0(:)=0 and n0T=0.
  !C = matmul(n0,Rinv)
  !n0T= dot_product(n0,C)
  n0(:,ibeta)=0.
  n0T(ibeta)=0.

  y_small = y(timeindex_nonzero, ibeta)
  m0 = z - y_small
  m0T(ibeta)=sum((m0/sigma_measure)**2)

enddo    ! beta loop 


! MCMC loop
!###############################################################
do it=1,(nIt+burn_in)
   
   call random_number(u)


   ! DO REVERSIBLE JUMP
 !  if (rjmcmc .EQ. 1) then
 !     remain_it = FLOOR(7*u) + 1  
 !   ! FIXED-DIMENSION CASE
 !  else if (rjmcmc .NE. 1) then
 !     remain_it = FLOOR(2*u) + 4 
 !  else
 !       WRITE(*,*) 'rjmcmc not specified correctly'
 !       CALL ABORT
 !  endif

   ! DO REVERSIBLE JUMP
   if (rjmcmc .EQ. 1 .and. numsites .EQ. 1) then
       remain_it = FLOOR(5*u) + 1  
   else if (rjmcmc .EQ. 1 .and. numsites .GT. 1) then
       remain_it = FLOOR(7*u) + 1    
 
   ! FIXED-DIMENSION CASE
   else if (rjmcmc .NE. 1 .and. numsites .EQ. 1) then
       remain_it = FLOOR(2*u) + 4 
   else if (rjmcmc .NE. 1 .and. numsites .GT. 1) then
      remain_it = FLOOR(4*u) + 4 

   ! ELSE  ABORT - HAS TO BE ONE OF THEM!!
   else
        WRITE(*,*) 'Either rjmcmc or numsites not specified correctly'
        CALL ABORT
   endif

   !remain_it = FLOOR(3*u) + 1 
   !remain_it = 8

!$OMP PARALLEL DO DEFAULT(SHARED) private(ibeta, betaib,kib,xib,pdf_param1ib,pdf_param2ib, plonib, platib), &
!$OMP& private(regions_vib,h_aggib,n0ib,n0Tib, sigma_ytib, sigma_ysib, detval_Tib, detval_Sib, detvalib), &
!$OMP& private(tauib, Qinvib, Rinvib, Sinvib, Tinvib, yib, m0Tib, kIC), &
!$OMP& private(kib1,xib1,pdf_param1ib1,pdf_param2ib1, plonib1, platib1), &
!$OMP& private(regions_vib1,h_aggib1,n0ib1,n0Tib1, sigma_ytib1, sigma_ysib1), &
!$OMP& private(detval_Tib1, detval_Sib1, detvalib1, acceptxib1, rejectxib1, acceptib1, rejectib1), &
!$OMP& private(tauib1, Qinvib1, Rinvib1, Sinvib1, Tinvib1, yib1, m0Tib1, u, acceptyib1, rejectyib1), &
!$OMP& private(rhoib, nuib, Uinvib, detval_Uib, rhoib1, nuib1, Uinvib1, detval_Uib1), &
!$OMP& shared(x,n0,n0T, k, pdf_param1, pdf_param2, h_agg, plon,plat, regions_v)
  ! do ibeta=1,nbeta

       ibeta=1
       betaib = beta(ibeta)
       kib = k(ibeta)
       xib  = x(:,ibeta)
       pdf_param1ib = pdf_param1(:,ibeta)
       pdf_param2ib = pdf_param2(:,ibeta)

       plonib = plon(:,ibeta)
       platib = plat(:,ibeta)
       regions_vib = regions_v(:,ibeta)
       h_aggib = h_agg(:,:,ibeta)
       n0ib = n0(:,ibeta)
       n0Tib = n0T(ibeta)

       sigma_ytib = sigma_yt(:,ibeta)
       sigma_ysib = sigma_ys(:,ibeta)

       detval_Uib = detval_U(ibeta)
       detval_Tib = detval_T(ibeta)
       detval_Sib = detval_S(ibeta)
       detvalib = detval(ibeta)

       tauib = tau(ibeta)
       rhoib = rho(ibeta)
       nuib = nu(ibeta)

       Qinvib = Qinv(:,:,ibeta)
       Rinvib = Rinv(:,:,ibeta)
       Sinvib = Sinv(:,:,ibeta)
       Tinvib = Tinv(:,:,ibeta)
       Uinvib = Uinv(:,:,ibeta)

       yib = y(:,ibeta)

       m0Tib = m0T(ibeta)

       kIC = kib+nIC
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Y UPDATE EVERY ITERATION
       call  y_update(betaib, yib, z, n0ib, n0Tib, m0Tib, Rinvib, sigma_measure, accept_y, reject_y, &
                    stepsize_y, timeindex_nonzero, nmeasuremax, nmeasure, it, burn_in,  &
                    n0Tib1, m0Tib1, yib1, n0ib1, acceptyib1, rejectyib1)

            n0T(ibeta) = n0Tib1
            m0T(ibeta) = m0Tib1
            n0(:,ibeta) = n0ib1
            y(:,ibeta) = yib1

            n0Tib = n0Tib1       
            n0ib = n0ib1
            yib = yib1

               
            if (betaib .EQ. 1.) then 
               accept_y = acceptyib1
               reject_y = rejectyib1
            endif
 
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (remain_it .EQ. 4) then              ! X UPDATE

            call x_update(betaib,kib, xib, pdf_param1ib,pdf_param2ib, &
                         h_aggib,n0ib,n0Tib, Rinvib, Sinvib, Tinvib,stepsize, &
                         pdf_p1_hparam1, pdf_p1_hparam2, stepsize_pdf_p1, pdf_param1_pdf, &
                         pdf_p2_hparam1, pdf_p2_hparam2, stepsize_pdf_p2, pdf_param2_pdf, &
                         accept,reject,x_pdf_all, it, burn_in, nIC, kICmax, nmeasuremax, nIC1, &
                         numsites, nmeasuretotal, &
                         xib1, pdf_param1ib1, pdf_param2ib1, n0ib1, n0Tib1, acceptxib1, rejectxib1)

            x(:,ibeta) = xib1
            n0(:,ibeta) = n0ib1
            n0T(ibeta) = n0Tib1 
            pdf_param1(:,ibeta) = pdf_param1ib1
            pdf_param2(:,ibeta) = pdf_param2ib1
            
            if (betaib .EQ. 1.) then 
               accept = acceptxib1
               reject = rejectxib1
            endif
            
       elseif (remain_it .EQ. 1) then       ! BIRTH
               
               call birth(betaib,kib, xib, h_aggib,yib,n0ib,n0Tib,Rinvib, plonib, platib, regions_vib, lon,lat, & 
                          h_v, pdf_param1ib(nIC1), pdf_param2ib(nIC1), x_pdf_all(nIC1), Sinvib, Tinvib, &
                          lonmin, lonmax, latmin,latmax, sigma_bd, &
                          accept_birth, reject_birth,it,burn_in,nIC,kICmax,kmax,nmeasuremax, &
                          numsites, nmeasuretotal, Ngrid,nlon,nlat, &
                          kib1, xib1, h_aggib1, n0ib1, n0Tib1, regions_vib1, plonib1, platib1, acceptib1, rejectib1)

                k(ibeta) = kib1
                x(:,ibeta) = xib1
                plon(:,ibeta) = plonib1
                plat(:,ibeta) = platib1
                regions_v(:,ibeta) = regions_vib1
                h_agg(:,:,ibeta) = h_aggib1
                n0(:,ibeta) = n0ib1
                n0T(ibeta) = n0Tib1


               if (betaib .EQ. 1.) then 
                   accept_birth=acceptib1
                   reject_birth=rejectib1
               endif

           elseif (remain_it .EQ. 2) then    ! DEATH

               call death(betaib,kib, xib, h_aggib, yib, n0ib,n0Tib,Rinvib, Sinvib, Tinvib, plonib, platib, regions_vib, &
                          lon,lat,h_v, pdf_param1ib(nIC1), pdf_param2ib(nIC1), x_pdf_all(nIC1), sigma_bd, &
                          accept_death, reject_death, it, burn_in,nIC, kICmax, kmin, kmax, nmeasuremax, &
                          numsites, nmeasuretotal, Ngrid,nlon,nlat, &
                          kib1, xib1, h_aggib1,n0ib1, n0Tib1, regions_vib1, plonib1, platib1, acceptib1, rejectib1)
                     
               k(ibeta) = kib1
               x(:,ibeta) = xib1
               plon(:,ibeta) = plonib1
               plat(:,ibeta) = platib1
               regions_v(:,ibeta) = regions_vib1
               h_agg(:,:,ibeta) = h_aggib1
               n0(:,ibeta) = n0ib1
               n0T(ibeta) = n0Tib1
      
               if (betaib .EQ. 1.) then 
                   accept_death=acceptib1
                   reject_death=rejectib1
               endif

           elseif (remain_it .EQ. 3) then    ! MOVE
               
               call move(betaib,kib, xib, h_aggib, yib ,n0ib,n0Tib,Rinvib, Sinvib, Tinvib, plonib, platib, regions_vib, & 
                         lon,lat,h_v, lonmin, lonmax, latmin,latmax, sigma_clon, sigma_clat, accept_move, reject_move, it, &
                         burn_in, nIC, kICmax, kIC, kmax, nmeasuremax, numsites, nmeasuretotal, Ngrid,nlon,nlat, &
                         h_aggib1, n0ib1, n0Tib1, regions_vib1, plonib1, platib1, acceptib1, rejectib1)
               
               plon(:,ibeta) = plonib1
               plat(:,ibeta) = platib1
               regions_v(:,ibeta) = regions_vib1
               h_agg(:,:,ibeta) = h_aggib1
               n0(:,ibeta) = n0ib1
               n0T(ibeta) = n0Tib1
               
                !do jj=1,ydim2   
                !     y_vart(R_indices(:,jj)) = sigma_ytib   ! Provided ydim2 isn't too big then should be fine
                !enddo  

               !n0T2 = sum((n0ib1/y_vart)**2)

               !write(*,*) n0T2, n0Tib1
               
               if (betaib .EQ. 1.) then 
                  accept_move=acceptib1
                  reject_move=rejectib1
               endif

          elseif (remain_it .EQ. 5) then  ! HYPERPARAMETERS UPDATE
             

           !   call random_number(u3)
              
             ! Only one site
           !  if (numsites .EQ. 1) then
           !     remain_hparam = FLOOR(2*u3) + 1  
           !  ! MULTIPLE SITE CASE
           !  else if (numsites .GT. 1) then
           !     remain_hparam = FLOOR(4*u3) + 1 
           !  endif
             
              ! SIGMA_YT UPDATE

             ! if (remain_hparam .EQ. 1) then

              call sigma_yt_update(betaib, sigma_ytib, sigma_yt_ap, &
                   detvalib, detval_Tib, detval_Sib, sigma_yt_hparams, stepsize_sigma_yt, sigma_yt_pdf, R_indices, &
                   Tinvib, Rinvib, Qinvib, Sinvib, deltatime, tauib, n0ib,n0Tib, &
                   accept_sigma_yt, reject_sigma_yt, it, burn_in, nmeasuretotal, nmeasuremax, numsites, ydim1, ydim2, &
                   n0Tib1, acceptib1, rejectib1, sigma_ytib1, detvalib1, detval_Tib1, Rinvib1, Tinvib1) 


           
              sigma_yt(:,ibeta) = sigma_ytib1  
              n0T(ibeta) = n0Tib1
              detval(ibeta) = detvalib1
              detval_T(ibeta) = detval_Tib1
              Tinv(:,:,ibeta)=Tinvib1 
              Rinv(:,:,ibeta)=Rinvib1 
        
              n0Tib=n0Tib1
              sigma_ytib=sigma_ytib1  
              detvalib=detvalib1
              detval_Tib=detval_Tib1
              Tinvib=Tinvib1
              Rinvib=Rinvib1

              if (betaib .EQ. 1.) then 
               accept_sigma_yt = acceptib1
               reject_sigma_yt = rejectib1
              endif
             
             
             ! elseif (remain_hparam .EQ. 2) then

              ! TAU UPDATE 

             ! if (rjmcmc .eq. 5) then
              call tau_update(betaib, tauib, sigma_ytib, R_indices,  &
                   detvalib, detval_Tib, detval_Sib, tau_hparams(1), tau_hparams(2), stepsize_tau, tau_pdf,  &
                   Rinvib, Tinvib, Qinvib, Sinvib, deltatime, n0ib, n0Tib, &
                   accept_tau, reject_tau, it, burn_in, nmeasuretotal, nmeasuremax, numsites, ydim1,ydim2, &
                   n0Tib1, acceptib1, rejectib1, tauib1, Rinvib1, Tinvib1, Qinvib1, detvalib1, detval_Tib1)

              n0T(ibeta) = n0Tib1
              tau(ibeta) = tauib1
              Tinv(:,:,ibeta)=Tinvib1 
              Rinv(:,:,ibeta)=Rinvib1 
              Qinv(:,:,ibeta)=Qinvib1 
              detval(ibeta) = detvalib1
              detval_T(ibeta) = detval_Tib1

              !n0Tib = n0Tib1
              !tauib = tauib1
              !Tinvib=Tinvib1 
              !Rinvib=Rinvib1 
              !Qinvib=Qinvib1 
              !detvalib = detvalib1
              !detval_Tib = detval_Tib1

              if (betaib .EQ. 1.) then 
               accept_tau = acceptib1
               reject_tau = rejectib1
              endif

  
                   
              elseif (remain_it .EQ. 6) then  ! SIGMA_YS UPDATE
             
                   call sigma_ys_update(betaib, sigma_ysib, sigma_ys_ap, &
                        detvalib, detval_Sib, detval_Tib, detval_Uib, &
                        sigma_ys_hparams, stepsize_sigma_ys, sigma_ys_pdf, &
                        Sinvib, Rinvib, Tinvib, Uinvib, n0ib, n0Tib, &
                        accept_sigma_ys, reject_sigma_ys, it, burn_in, nmeasuretotal, nmeasuremax, numsites, &
                        n0Tib1, acceptib1, rejectib1, sigma_ysib1, detvalib1, detval_Sib1, Rinvib1, Sinvib1) 


                    sigma_ys(:,ibeta) = sigma_ysib1
                    n0T(ibeta) = n0Tib1
                    detval(ibeta) = detvalib1
                    detval_S(ibeta) = detval_Sib1
                    Sinv(:,:,ibeta)=Sinvib1 
                    Rinv(:,:,ibeta)=Rinvib1 
        
                    !sigma_ysib = sigma_ysib1
                    !n0Tib = n0Tib1
                    !detvalib = detvalib1
                    !detval_Sib = detval_Sib1
                    !Sinvib=Sinvib1 
                    !Rinvib=Rinvib1 


                    if (betaib .EQ. 1.) then 
                       accept_sigma_ys = acceptib1
                       reject_sigma_ys = rejectib1
                    endif
                
                 
               elseif (remain_it .EQ. 7) then ! RHO + NU UPDATE
                    

                   call rho_nu_update(betaib, rhoib, nuib, sigma_ysib,  &
                        detvalib, detval_Tib, detval_Sib, detval_Uib, &
                        rho_hparams(1), rho_hparams(2), stepsize_rho, rho_pdf,  &
                        nu_hparams(1), nu_hparams(2), stepsize_nu, nu_pdf,  &
                        Rinvib, Tinvib, Sinvib, Uinvib, distance, &
                        n0ib, n0Tib, accept_rho, reject_rho, it, burn_in, nmeasuretotal, nmeasuremax, numsites, &
                        n0Tib1, acceptib1, rejectib1, rhoib1, nuib1, Rinvib1, Sinvib1, &
                        Uinvib1, detvalib1, detval_Sib1, detval_Uib1) 


                    n0T(ibeta) = n0Tib1
                    rho(ibeta) = rhoib1
                    nu(ibeta) = nuib1
                    Rinv(:,:,ibeta)=Rinvib1 
                    Sinv(:,:,ibeta)=Sinvib1 
                    Uinv(:,:,ibeta)=Uinvib1 
                    detval(ibeta) = detvalib1
                    detval_S(ibeta) = detval_Sib1
                    detval_U(ibeta) = detval_Uib1

                    if (betaib .EQ. 1.) then 
                        accept_rho = acceptib1
                        reject_rho = rejectib1
                    endif
                  

              ! endif   !remian_hparam

          endif     ! remain_it

           
  ! enddo    ! beta loop
!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Store xit and swap betas
  
   IF (it .GT. burn_in/2) THEN      ! Begin swaps after half of burn-in time
        remain_swap = modulo(it,2)
        IF (remain_swap .EQ. 1) THEN
            call random_number(u1)   
            pair1 = FLOOR(nbeta*u1) + 1
            call random_number(u2)  
            pair2 = FLOOR(nbeta*u2) + 1
            if (pair1 .EQ. pair2) THEN
                if (pair2 .EQ. 1) then
                    pair2=pair1+1
                else 
                    pair2=pair1-1
                endif
            endif

            beta1=beta(pair1)*1.
            beta2=beta(pair2)*1.
            !pT_chain = (beta2-beta1)*(n0T(pair2)/2.-n0T(pair1)/2.+detval(pair2)-detval(pair1))  ! detvals should be inverse determinants so signs this way round


           pT_chain = (beta2-beta1)*(n0T(pair2)/2.-n0T(pair1)/2.+ m0T(pair2)/2.-m0T(pair1)/2.+detval(pair2)-detval(pair1))
          ! pT_chain = (beta2-beta1)*(m0T(pair2)/2.-m0T(pair1)/2.) 

            !write(*,*) pT_chain, m0T(pair2), m0T(pair1), n0T(pair2), n0T(pair1)
            !write(*,*) m0T(pari2), m0T(pair1)
           ! write(*,*) beta2, n0T(pair2), detval(pair2)
           ! stop
   
            call random_number(randomu)
            if (alog(randomu) .LE. pT_chain) then
           !     beta(pair2)=beta1*1.
           !     beta(pair1)=beta2*1.
                accept_swap=accept_swap+1
            else
                reject_swap=reject_swap+1
            endif      ! pT_chain if      
         ENDIF      ! reamin_it =0 if
   ENDIF          ! it > burn_in/2
    


   IF (it .GT. burn_in) THEN     
        remain = modulo(it,nsub)          ! nsub typically = 100
        if (remain .EQ. 0) then
            do ib=1,nbeta
               if (beta(ib) .EQ. 1.) then
                  x_it(:,it_sub)=x(:,ib)
                  plon_it(:,it_sub)=plon(:,ib)
                  plat_it(:,it_sub)=plat(:,ib)
                  k_it(it_sub)=k(ib)
                  regions_it(:,it_sub)=regions_v(:,ib)
                  sigma_yt_it(:,it_sub)=sigma_yt(:,ib)
                  sigma_ys_it(:,it_sub)=sigma_ys(:,ib)
                  tau_it(it_sub) = tau(ib)  
                  rho_it(it_sub) = rho(ib)  
                  nu_it(it_sub) = nu(ib)  
                  y_it(:,it_sub) = y(:,ib)
                  n0T_it(it_sub)=n0T(ib)
                  pdf_param1_it(:,it_sub)=pdf_param1(:,ib)   
                  pdf_param2_it(:,it_sub)=pdf_param2(:,ib)  
                  it_sub=it_sub+1
               endif
            enddo
        endif
   ENDIF           ! it >= burn_in


   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
enddo  ! It loop    ! END OF MCMC LOOP

! OUTPUTS
k_out = k_it
x_out=x_it
regions_out = regions_it
plon_out = plon_it
plat_out=plat_it
sigma_yt_out=sigma_yt_it
sigma_ys_out=sigma_ys_it
tau_out = tau_it
rho_out=rho_it
nu_out=nu_it
n0T_out=n0T_it
pdf_param1_out=pdf_param1_it
pdf_param2_out=pdf_param2_it
y_out=y_it
END SUBROUTINE hbtdmcmc



SUBROUTINE x_update(beta,k, x, pdf_param1_all,pdf_param2_all,  &
h_agg,n0,n0T,Rinv, Sinv, Tinv, stepsize, &
pdf_p1_hparam1_all, pdf_p1_hparam2_all, stepsize_pdf_p1, pdf_param1_pdf, &
pdf_p2_hparam1_all, pdf_p2_hparam2_all, stepsize_pdf_p2, pdf_param2_pdf, &
accept, reject, x_pdf_all, it, burn_in, nIC, kICmax, nmeasuremax, nIC1, &
numsites, nmeasuretotal, &
x_out, pdf_param1_out, pdf_param2_out, n0_out, n0T_out, accept_out, reject_out) 

Implicit none 
INTEGER nmeasuremax, it, burn_in, k, nIC, kICmax, nIC1
INTEGER nmeasuretotal, numsites
REAL beta, n0T, n0T_out 
REAL x(kICmax) 
REAL x_out(kICmax) 
REAL h_agg(nmeasuremax,kICmax)   
REAL n0(nmeasuremax) 
REAL n0_out(nmeasuremax) 
REAL Rinv(nmeasuremax, nmeasuremax)
REAL Sinv(numsites, numsites)
REAL Tinv(nmeasuretotal, nmeasuretotal)
REAL dy(nmeasuremax)
REAL n1(nmeasuremax) 
INTEGER x_pdf_all(nIC1)
REAL pdf_param1_all(nIC1), pdf_param2_all(nIC1) 
REAL pdf_param1        
REAL pdf_param2
REAL stepsize(nIC1)
INTEGER accept(nIC1), reject(nIC1)
INTEGER accept_out(nIC1), reject_out(nIC1)
INTEGER xi, x_pdf
REAL dx, n1T, pT, randomu, p0,p1, u
REAL C(nmeasuremax)
INTEGER pdf_param1_pdf, pdf_param2_pdf
REAL pdf_p1_hparam1_all(nIC1), pdf_p1_hparam2_all(nIC1) 
REAL pdf_p2_hparam1_all(nIC1), pdf_p2_hparam2_all(nIC1) 
REAL pdf_param1_out(nIC1), pdf_param2_out(nIC1)
REAL pdf_p2_hparam1, pdf_p2_hparam2, pdf_p1_hparam1, pdf_p1_hparam2
REAL stepsize_pdf_p1(nIC1), stepsize_pdf_p2(nIC1)
REAL stepsize_pdf_p10, stepsize_pdf_p20, stepsize0
INTEGER ki, xx
INTEGER elem(2)
REAL p0_temp, p1_temp
REAL dpdf_param1, pdf_param1_new, dpdf_param2, pdf_param2_new
REAL p0_pdf_param2, p1_pdf_param2, p0_pdf_param1, p1_pdf_param1
REAL x_new(kICmax)
REAL n1_arr(nmeasuretotal,numsites)
REAL C_arr(nmeasuretotal,numsites)   

!INPUTS:
! ARRAYS :: beta, nbeta, pdf_param1, pdf_param2, x, n0, h_agg, sigma_y, n0T, k
! INTEGERS :: nbeta, accept, reject, nMeasure, x_pdf, k_max, it, burn_in    REAL h_agg(nMeasure,kmax,nbeta) 
! REAL ::  xmin, xmax

! OUTPUTS:
! x_out, n0_out, n0T_out, accept_out, reject_out

   ! CHANGE OF EMISSIONS VALUES
 ! call random_number(u)   
 ! elem(1)= FLOOR((nIC)*u)+1 
 ! elem(2) = FLOOR(k*u)+1+nIC

 ! do xx=1,2
 ! xi=elem(xx)

 do xi = 1, k+nIC

  !xi=FLOOR((nIC+k)*u)+1
  if (xi .LE. nIC) then
     x_pdf = x_pdf_all(xi)
     pdf_param1 = pdf_param1_all(xi)
     pdf_param2 = pdf_param2_all(xi)
     pdf_p1_hparam1 = pdf_p1_hparam1_all(xi)
     pdf_p1_hparam2 = pdf_p1_hparam2_all(xi)
     pdf_p2_hparam1 = pdf_p2_hparam1_all(xi)
     pdf_p2_hparam2 = pdf_p2_hparam2_all(xi)
     stepsize0 = stepsize(xi)
     stepsize_pdf_p10 = stepsize_pdf_p1(xi)
     stepsize_pdf_p20 = stepsize_pdf_p2(xi)
  else if (xi .GT. nIC) then
     x_pdf = x_pdf_all(nIC1)
     pdf_param1 = pdf_param1_all(nIC1)
     pdf_param2 = pdf_param2_all(nIC1)
     pdf_p1_hparam1 = pdf_p1_hparam1_all(nIC1)
     pdf_p1_hparam2 = pdf_p1_hparam2_all(nIC1) 
     pdf_p2_hparam1 = pdf_p2_hparam1_all(nIC1)
     pdf_p2_hparam2 = pdf_p2_hparam2_all(nIC1)
     stepsize0 = stepsize(nIC1)
     stepsize_pdf_p10 = stepsize_pdf_p1(nIC1)
     stepsize_pdf_p20 = stepsize_pdf_p2(nIC1) 
  endif
  !stepsize_param2=pdf_param2*stepsize_pdf_p2
  !dx = random_normal()*stepsize
  !dpdf_param1 = random_normal()*stepsize_pdf_p1*pdf_param1
  !dpdf_param2 = random_normal()*stepsize_pdf_p2*pdf_param2
  !pdf_param1_new = pdf_param1 + dpdf_param1
  !pdf_param2_new = pdf_param2 + dpdf_param2

  dx = random_normal()*stepsize0
  dpdf_param1 = random_normal()*stepsize_pdf_p10
  dpdf_param2 = random_normal()*stepsize_pdf_p20
  pdf_param1_new = pdf_param1 + dpdf_param1
  pdf_param2_new = pdf_param2 + dpdf_param2

  x_new = x
  x_new(xi) =x(xi)+dx
  p0=0.
  p1=0.
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (x_pdf .EQ. 1) THEN  ! 1=UNIFORM
     if (x(xi)+dx .GT. pdf_param1 .and. x(xi)+dx .LT. pdf_param2) THEN     ! pdf_param1 = xmin pdf_param2=xmax          
          !Only accept if x is positive
          dy=h_agg(:,xi)*dx  
          n1=n0+dy
          !n1T=sum((n1/sigma_y)**2)

          C = matmul(n1,Rinv)
          n1T= dot_product(n1,C)
                  
          pT=(-0.5)*(n1T - n0T)*beta
 
          call random_number(randomu)
          if (alog(randomu) .le. pT) THEN
             !ACCEPT
             x(xi)=x(xi)+dx       
             n0=n1
             n0T=n1T
              if (xi .LE. nIC) then
                if (beta .EQ. 1. .and. it .GT. burn_in) accept(xi) = accept(xi) + 1
             else if (xi .GT. nIC) then      
                if (beta .EQ. 1. .and. it .GT. burn_in) accept(nIC1) = accept(nIC1) + 1
             endif
             !if (beta.EQ. 1. .and. it .GT. burn_in)  accept = accept + 1
          else
                 !REJECT
                 if (beta .EQ. 1. .and. it .GT. burn_in) then 
                    if (xi .LE. nIC) then 
                     reject(xi) = reject(xi) + 1
                    else
                     reject(nIC1) = reject(nIC1) + 1  
                    endif
                 endif
                 !if (beta .EQ. 1. .and. it .GT. burn_in) reject = reject + 1
          endif           ! randomu condition  
      else 
          !Reject if x is negative
          if (beta .EQ. 1. .and. it .GT. burn_in) then 
                 if (xi .LE. nIC) then 
                     reject(xi) = reject(xi) + 1
                 else
                     reject(nIC1) = reject(nIC1) + 1  
                 endif
          endif
          !if (beta .EQ. 1. .and. it .GT. burn_in) reject = reject + 1

      endif           ! xmin, xmax condition     

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
  else if (x_pdf .GE. 2) THEN    ! 2=GAUSSIAN  3=LOGNORMAL

         dy=h_agg(:,xi)*dx 
         n1=n0+dy

         if (numsites .EQ. 1) then
            C = matmul(n1,Rinv)
            n1T= dot_product(n1,C)
         else
            n1_arr = reshape(n1, (/nmeasuretotal, numsites/))
            C_arr = matmul(matmul(Tinv,n1_arr),Sinv)
            C = reshape(C_arr,(/nmeasuremax/))
            n1T = dot_product(n1,C)
         endif   ! numsites

! n1_arr = reshape(n1, (/nmeasuretotal, numsites/))
!call sgemm ('N', 'N', nmeasuretotal, numsites, nmeasuretotal, aa, Tinv_current, nmeasuretotal, n1_arr, nmeasuretotal, bb, dum3, nmeasuretotal)
!call sgemm ('N', 'N', nmeasuretotal, numsites, numsites, aa, dum3, nmeasuretotal, Sinv_current, numsites, bb, C_arr, nmeasuretotal)
!C = reshape(C_arr,(/nmeasuremax/))


         !n1T=sum((n1/sigma_y)**2)
                 
         call calc_pdf(pdf_param1,pdf_p1_hparam1,pdf_p1_hparam2,pdf_param1_pdf, p0_pdf_param1) 
         call  calc_pdf(pdf_param1_new, pdf_p1_hparam1,pdf_p1_hparam2,pdf_param1_pdf, p1_pdf_param1)

         call calc_pdf(pdf_param2,pdf_p2_hparam1,pdf_p2_hparam2,pdf_param2_pdf, p0_pdf_param2) 
         call  calc_pdf(pdf_param2_new, pdf_p2_hparam1,pdf_p2_hparam2,pdf_param2_pdf, p1_pdf_param2)

         if (xi .LE. nIC) then
           
            !do ki =1,nIC
               call calc_pdf(x(xi),pdf_param1,pdf_param2,x_pdf, p0_temp)           ! Will apply whatever the PDF
               call calc_pdf(x_new(xi),pdf_param1_new,pdf_param2_new,x_pdf, p1_temp)        
               p0=p0+p0_temp
               p1=p1+p1_temp
            !enddo

         else if (xi .GT. nIC) then

            do ki =nIC+1,k+nIC
               call calc_pdf(x(ki),pdf_param1,pdf_param2,x_pdf, p0_temp)           ! Will apply whatever the PDF
               call calc_pdf(x_new(ki),pdf_param1_new,pdf_param2_new,x_pdf, p1_temp)        
               p0=p0+p0_temp
               p1=p1+p1_temp
            enddo
         endif


         !call calc_pdf(x(xi),pdf_param1,pdf_param2,x_pdf, p0)           ! Will apply whatever the PDF
         !call calc_pdf(x(xi)+dx,pdf_param1,pdf_param2,x_pdf, p1)
        
         ! pT = p1-p0-0.5*(n1T - n0T)*beta  ! All p1s already in log space

          pT = p1+p1_pdf_param1+p1_pdf_param2-p0-p0_pdf_param1-p0_pdf_param2 -0.5*(n1T - n0T)*beta ! All p1s already in log space


         if (pdf_param1_pdf .eq. 1) then
             if (pdf_param1_new .lt. pdf_p1_hparam1) pT = -1.e20
             if (pdf_param1_new .gt. pdf_p1_hparam2) pT = -1.e20
         endif

         if (pdf_param2_pdf .eq. 1) then
             if (pdf_param2_new .lt. pdf_p2_hparam1) pT = -1.e20
             if (pdf_param2_new .gt. pdf_p2_hparam2) pT = -1.e20
         endif

         call random_number(randomu)
         if (alog(randomu) .LE. pT) THEN

             !ACCEPT
             x(xi)=x(xi)+dx
             n0=n1
             n0T=n1T
             if (xi .LE. nIC) then
                pdf_param1_all(xi)=pdf_param1_new
                pdf_param2_all(xi)=pdf_param2_new
                if (beta .EQ. 1. .and. it .GT. burn_in) accept(xi) = accept(xi) + 1
             else if (xi .GT. nIC) then
                pdf_param1_all(nIC1)=pdf_param1_new
                pdf_param2_all(nIC1)=pdf_param2_new
                if (beta .EQ. 1. .and. it .GT. burn_in) accept(nIC1) = accept(nIC1) + 1
             endif



             !if (xi .LE. nIC) then
             !   pdf_param1_all(1)=pdf_param1_new
             !   pdf_param2_all(1)=pdf_param2_new
             !else if (xi .GT. nIC) then
             !   pdf_param1_all(2)=pdf_param1_new
             !   pdf_param2_all(2)=pdf_param2_new
             !endif
             
             !if (beta .EQ. 1. .and. it .GT. burn_in) accept = accept + 1
                                
         else
             !REJECT
             if (beta .EQ. 1. .and. it .GT. burn_in) then 
                 if (xi .LE. nIC) then 
                     reject(xi) = reject(xi) + 1
                 else
                     reject(nIC1) = reject(nIC1) + 1  
                 endif
             endif
             !if (beta .EQ. 1. .and. it .GT. burn_in) reject = reject + 1
                              
         endif   ! randomu condition

  endif   ! x_pdf condition
  enddo   ! xi loop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

x_out=x
n0_out=n0
n0T_out=n0T
pdf_param1_out=pdf_param1_all
pdf_param2_out=pdf_param2_all
accept_out=accept
reject_out=reject
END SUBROUTINE x_update



SUBROUTINE birth(beta,k, x, h_agg,y,n0,n0T,Rinv, plon, plat, regions_v, lon,lat, & 
h_v,pdf_param1, pdf_param2, x_pdf, Sinv, Tinv, &
lonmin, lonmax, latmin,latmax, sigma_bd, &
accept_birth, reject_birth, it, burn_in, nIC, kICmax, kmax, nmeasuremax, &
numsites, nmeasuretotal, Ngrid,nlon,nlat, &
k_out, x_out, h_agg_out, n0_out, n0T_out, regions_v_out, plon_out, plat_out, accept_out, reject_out)

IMPLICIT NONE
! Dimensions
INTEGER nmeasuremax, Ngrid, nlon,nlat, k,kmax, nIC, kICmax
INTEGER nmeasuretotal, numsites
REAL lonmin,lonmax,latmin,latmax,sigma_bd
! Single Variables
INTEGER accept_birth, reject_birth, it, burn_in, x_pdf 
REAL beta, n0T
! Input arrays
REAL x(kICmax) 
REAL pdf_param1 
REAL pdf_param2
REAL h_agg(nmeasuremax,kICmax)  
REAL y(nmeasuremax) 
REAL n0(nmeasuremax) 
REAL Rinv(nmeasuremax,nmeasuremax)
REAL Tinv(nmeasuretotal,nmeasuretotal)
REAL Sinv(numsites,numsites)
REAL plon(kmax)
REAL plat(kmax)
INTEGER regions_v(Ngrid)
REAL lon(nlon)
REAL lat(nlat)
REAL h_v(nmeasuremax, Ngrid)
! Outputs
INTEGER k_out
REAL x_out(kICmax) 
REAL h_agg_out(nmeasuremax,kICmax)
REAL n0_out(nmeasuremax) 
REAL n0T_out
INTEGER regions_v_out(Ngrid)
REAL plon_out(kmax)
REAL plat_out(kmax)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER ri, rib, k1, errstat,jj, kIC      
REAL u, u2, plon_new, plat_new, c_new, x_new, n1Tb, pT_birth, randomu
REAL mu, sigma
INTEGER regions_v1b(Ngrid)       
REAL n1b(nmeasuremax)  
REAL C(nmeasuremax)          
REAL n1_arr(nmeasuretotal,numsites)
REAL C_arr(nmeasuretotal,numsites)  
! Allocatable arrays
REAL, DIMENSION(:),ALLOCATABLE :: plon1b, plat1b, x1b
REAL, DIMENSION(:,:), ALLOCATABLE :: h_agg2
REAL,PARAMETER     :: pi = 3.14159265 


k1=k+1
kIC=k1+nIC

if (k1 .LT. kmax) THEN            
  
  allocate(plon1b(k1))     
  allocate(plat1b(k1))     
  allocate(h_agg2(nmeasuremax, kIC))
  allocate(x1b(kIC))               

  ! 1. Select new cell location - needs to be different to current locations
  call random_number(u)
  plon_new = lonmin+(lonmax-lonmin)*u
  call random_number(u2)
  plat_new = latmin+(latmax-latmin)*u2
  
  plon1b(1:k) = plon(1:k)
  plon1b(k1) = plon_new
  plat1b(1:k) = plat(1:k)
  plat1b(k1) = plat_new
                                  
  ! 2. Recalculate voronoi cells

  call closest_grid(regions_v1b,lon,lat, plon1b(1:k1),plat1b(1:k1), k1, nlon,nlat)

   h_agg2(:,:)=0.

   if (nIC .GT. 0) then 
       h_agg2(:,1:nIC)=h_agg(:,1:nIC)
   endif
       
    do ri=1,k1
       do jj=1,Ngrid
         if (regions_v1b(jj) .EQ. ri) then
            h_agg2(:,ri+nIC) = h_agg2(:,ri+nIC)+h_v(:,jj)
         endif
      enddo
   enddo
 
 !#######################################################################
                    
  call closest_grid_sng(rib,plon_new,plat_new, plon(1:k),plat(1:k), k) 
                                  
  c_new = x(rib+nIC)
  x_new = random_normal()*sigma_bd+c_new
  
  if (x_pdf .EQ. 1) THEN  ! 1=UNIFORM  
                
      if (x_new .GT. pdf_param1 .and. x_new .LT. pdf_param2) THEN
                       
          x1b(1:k+nIC)=x(1:k+nIC)
          x1b(kIC) = x_new

          n1b=matmul(h_agg2,x1b)-y                    
          !n1Tb=sum((n1b/sigma_y)**2)
          C = matmul(n1b,Rinv)
          n1Tb= dot_product(n1b,C)

          pT_birth = alog(sqrt(2*pi)*sigma_bd/(pdf_param2-pdf_param1)) &
              + ((x_new-c_new)**2)/2./sigma_bd**2 - 0.5*(n1Tb - n0T)*beta
       

         call random_number(randomu)             
         if (alog(randomu) .LE. pT_birth) THEN
           !ACCEPT
           k=k1
           x(:)=0.
           x(1:kIC)=x1b(1:kIC)
           h_agg(:,:)=0.
           h_agg(:,1:kIC)=h_agg2(:,1:kIC)
           n0=n1b
           n0T=n1Tb
           regions_v(:)=regions_v1b
           plon(:)=0.
           plat(:)=0.
           plon(1:k1)=plon1b(1:k1)
           plat(1:k1)=plat1b(1:k1)
           if (beta .EQ. 1. .and. it .GT. burn_in) accept_birth=accept_birth+1    
         else 
           !REJECT
           if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1  
         endif

      else
          !REJECT if x_new is negative
          if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1
      endif   

  else if (x_pdf .GE. 2) THEN     ! 2=GAUSSIAN 3=LOGNORMAL 

          x1b(1:k+nIC)=x(1:k+nIC)
          x1b(kIC) = x_new
          
          n1b=matmul(h_agg2,x1b)-y                                         
          !n1Tb=sum((n1b/sigma_y)**2)

          if (numsites .EQ. 1) then
             C = matmul(n1b,Rinv)
             n1Tb= dot_product(n1b,C)
          else
             n1_arr = reshape(n1b, (/nmeasuretotal, numsites/))
             C_arr = matmul(matmul(Tinv,n1_arr),Sinv)
             C = reshape(C_arr,(/nmeasuremax/))
             n1Tb = dot_product(n1b,C)
          endif    ! numsites

          if (x_pdf .EQ. 2) THEN   !2=GAUSSIAN

              pT_birth = alog(sigma_bd/pdf_param2) + ((x_new-c_new)**2)/2./sigma_bd**2 - &
                         ((x_new-pdf_param1)**2)/2./pdf_param2**2 &
                         - 0.5*(n1Tb - n0T)*beta 

          else if (x_pdf .EQ. 3) THEN
 
              mu = alog(pdf_param1) - 0.5*alog(1. + pdf_param2**2/pdf_param1**2)
              sigma=sqrt(alog((pdf_param2/pdf_param1)**2 + 1.))

               pT_birth = alog(sigma_bd/sigma/x_new) + ((x_new-c_new)**2)/2./sigma_bd**2 - &
                         ((alog(x_new)-mu)**2)/2./sigma**2 &
                          - 0.5*(n1Tb - n0T)*beta             

         else
             stop
         endif
       
         call random_number(randomu)             
         if (alog(randomu) .LE. pT_birth) THEN
           !ACCEPT
           k=k1
           x(:)=0.
           x(1:kIC)=x1b(1:kIC)
           h_agg(:,:)=0.
           h_agg(:,1:kIC)=h_agg2(:,1:kIC)
           n0=n1b
           n0T=n1Tb
           regions_v(:)=regions_v1b
           plon(:)=0.
           plat(:)=0.
           plon(1:k1)=plon1b(1:k1)
           plat(1:k1)=plat1b(1:k1)
           if (beta .EQ. 1. .and. it .GT. burn_in) accept_birth=accept_birth+1    
         else 
           !REJECT
           if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1  
         endif

  endif       ! x_pdf

        
else
    !REJECT if k1 > kmax
    if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1 
endif


!! Deallocate arrays in each loop
if(allocated(plon1b))  deallocate(plon1b,stat=errstat)
  if (errstat /= 0) stop
if(allocated(plat1b))  deallocate(plat1b,stat=errstat)
  if (errstat /= 0) stop
if(allocated(h_agg2))  deallocate(h_agg2,stat=errstat)
  if (errstat /= 0) stop
if(allocated(x1b))  deallocate(x1b,stat=errstat)
  if (errstat /= 0) stop

k_out=k
x_out=x
h_agg_out=h_agg
n0_out=n0
n0T_out=n0T
regions_v_out=regions_v
plon_out=plon
plat_out=plat
accept_out=accept_birth
reject_out=reject_birth


END SUBROUTINE birth


SUBROUTINE death(beta,k, x, h_agg,y,n0,n0T,Rinv, Sinv, Tinv, plon, plat, regions_v, lon,lat, & 
h_v, pdf_param1, pdf_param2, x_pdf, &
sigma_bd, accept_death, reject_death, &
it, burn_in, nIC, kICmax, kmin, kmax, nmeasuremax, numsites, nmeasuretotal, Ngrid,nlon,nlat, &
k_out, x_out, h_agg_out, n0_out, n0T_out, regions_v_out, &
plon_out, plat_out, accept_out, reject_out)



IMPLICIT NONE
! Dimensions
INTEGER kmax,nmeasuremax,Ngrid,nlon,nlat, accept_death, reject_death, it, burn_in, kmin, k, nIC, kICmax
INTEGER numsites, nmeasuretotal
REAL sigma_bd, beta, n0T 
! Input arrays
REAL x(kICmax) 
REAL pdf_param1
REAL pdf_param2
REAL h_agg(nmeasuremax,kICmax) 
REAL y(nmeasuremax) 
REAL n0(nmeasuremax) 
REAL Rinv(nmeasuremax,nmeasuremax)
REAL Sinv(numsites,numsites)
REAL Tinv(nmeasuretotal,nmeasuretotal)
REAL plon(kmax)
REAL plat(kmax)
INTEGER regions_v(Ngrid)
REAL lon(nlon)
REAL lat(nlat)
REAL h_v(nmeasuremax, Ngrid)
INTEGER x_pdf
! Outputs
INTEGER k_out 
REAL n0T_out
REAL x_out(kICmax)
REAL h_agg_out(nmeasuremax,kICmax)
REAL n0_out(nmeasuremax) 
INTEGER regions_v_out(Ngrid)
REAL plon_out(kmax)
REAL plat_out(kmax)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER ri, rid, k1d, jj, ci_rm, errstat, kIC     
REAL u, plon_rm, plat_rm, x_cell, x_rm, n1Td, pT_death, randomu
REAL mu,sigma
INTEGER regions_v1d(Ngrid)      
REAL n1d(nmeasuremax)  
REAL C(nmeasuremax)           
REAL n1_arr(nmeasuretotal,numsites)
REAL C_arr(nmeasuretotal,numsites)  
! Allocatable arrays
REAL, DIMENSION(:),ALLOCATABLE :: plon1d, plat1d, x1d 
REAL, DIMENSION(:,:), ALLOCATABLE :: h_agg2d

REAL,PARAMETER     :: pi = 3.14159265 

!DEATH
k1d=k-1

kIC=k1d+nIC

if (k1d .GE. kmin) THEN            

  allocate(plon1d(k1d))     
  allocate(plat1d(k1d))     
  allocate(h_agg2d(nmeasuremax, kIC))
  allocate(x1d(kIC))  
           
  ! 1. Select new cell location - needs to be different to current locations
   
  call random_number(u)   
  ci_rm = FLOOR((k1d+1)*u) + 1
  
  plon_rm = plon(ci_rm)
  plat_rm = plat(ci_rm)
  x_rm = x(ci_rm+nIC)
  

  IF (nIC .GT. 0) THEN
     x1d(1:nIC) = x(1:nIC)
  ENDIF

  IF (ci_rm .EQ. 1) THEN
      plon1d(1:k1d) = plon(2:(k1d+1))
      plat1d(1:k1d) = plat(2:(k1d+1))
      x1d(nIC+1:kIC) = x(nIC+2:(kIC+1))
  
  ELSEIF (ci_rm .EQ. (k1d+1)) THEN
      plon1d(1:k1d) = plon(1:k1d)
      plat1d(1:k1d) = plat(1:k1d)
      x1d(nIC+1:kIC) = x(nIC+1:kIC)

  ELSE   
      plon1d(1:(ci_rm-1)) = plon(1:(ci_rm-1))                   
      plon1d(ci_rm:k1d) = plon((ci_rm+1):(k1d+1))
      plat1d(1:(ci_rm-1)) = plat(1:(ci_rm-1))                   
      plat1d(ci_rm:k1d) = plat((ci_rm+1):(k1d+1)) 

      x1d(nIC+1:(ci_rm+nIC-1)) = x(nIC+1:(ci_rm+nIC-1)) 
      x1d(ci_rm+nIC:kIC) = x((ci_rm+nIC+1):(kIC+1))

  ENDIF
                                            
  ! 2. Recalculate voronoi cells

  call closest_grid(regions_v1d,lon,lat, plon1d,plat1d, k1d, nlon,nlat)

   h_agg2d(:,:)=0.
  
   if (nIC .GT. 0) then
       h_agg2d(:,1:nIC)=h_agg(:,1:nIC)
   endif

   do ri=1,k1d
      do jj=1,Ngrid
         if (regions_v1d(jj) .EQ. ri) then
            h_agg2d(:,ri+nIC) = h_agg2d(:,ri+nIC)+h_v(:,jj)
         endif
      enddo
   enddo

 !#######################################################################
                    
  call closest_grid_sng(rid,plon_rm,plat_rm, plon1d,plat1d, k1d) 
                                  
  x_cell = x1d(rid+nIC)
                    
  n1d=matmul(h_agg2d,x1d)-y
  !n1Td=sum((n1d/sigma_y)**2) 

  if (numsites .EQ. 1) then
     C = matmul(n1d,Rinv)
     n1Td= dot_product(n1d,C)
  else
     n1_arr = reshape(n1d, (/nmeasuretotal, numsites/))
     C_arr = matmul(matmul(Tinv,n1_arr),Sinv)
     C = reshape(C_arr,(/nmeasuremax/))
     n1Td = dot_product(n1d,C)
  endif       ! numsites

                          
  ! ACCEPTANCE PROBABILITY  

  IF (x_pdf .EQ. 1) THEN  ! 1 = UNIFORM
     pT_death = alog((pdf_param2-pdf_param1)/sqrt(2.*pi)/sigma_bd) &
      - ((x_cell-x_rm)**2)/2./sigma_bd**2 -0.5*(n1Td - n0T)*beta

  ELSE IF (x_pdf .EQ. 2) THEN   ! 2=GAUSSIAN

      pT_death = alog(pdf_param2/sigma_bd) - ((x_cell-x_rm)**2)/2./sigma_bd**2 + &
                         ((x_rm-pdf_param1)**2)/2./pdf_param2**2 &
                          - 0.5*(n1Td - n0T)*beta


  ELSE IF (x_pdf .EQ. 3) THEN   ! 3=LOGNORMAL

      mu = alog(pdf_param1) - 0.5*alog(1. + pdf_param2**2/pdf_param1**2)
      sigma=sqrt(alog((pdf_param2/pdf_param1)**2 + 1.))
  
      pT_death = alog(sigma*x_rm/sigma_bd) - ((x_cell-x_rm)**2)/2./sigma_bd**2 + &
                         ((alog(x_rm)-mu)**2)/2./sigma**2 &
                          - 0.5*(n1Td - n0T)*beta

   
  ELSE
     ! write(*,*) 'There is an error with x_pdf in death loop - exiting!'
      stop
  ENDIF 

  !write(*,*) pT_death                             
  call random_number(randomu)             
  if (alog(randomu) .LE. pT_death) THEN
          !ACCEPT
           k=k1d
           x(:)=0.
           x(1:kIC)=x1d
           h_agg(:,:)=0.
           h_agg(:,1:kIC)=h_agg2d
           n0=n1d
           n0T=n1Td
           regions_v(:)=regions_v1d
           plon(:)=0.
           plat(:)=0.
           plon(1:k1d)=plon1d
           plat(1:k1d)=plat1d
           if (beta .EQ. 1. .and. it .GT. burn_in) accept_death=accept_death+1   

       else 
           !REJECT
           if (beta .EQ. 1. .and. it .GT. burn_in) reject_death=reject_death+1  
       endif
        
else
    !REJECT if k1d < kmin
    if (beta .EQ. 1. .and. it .GT. burn_in) reject_death=reject_death+1 
endif

!! Deallocate arrays in each loop
if(allocated(plon1d))  deallocate(plon1d,stat=errstat)
  if (errstat /= 0) stop
if(allocated(plat1d))  deallocate(plat1d,stat=errstat)
  if (errstat /= 0) stop
if(allocated(h_agg2d))  deallocate(h_agg2d,stat=errstat)
  if (errstat /= 0) stop
if(allocated(x1d))  deallocate(x1d,stat=errstat)
  if (errstat /= 0) stop


k_out=k
x_out=x
h_agg_out=h_agg
n0_out=n0
n0T_out=n0T
regions_v_out=regions_v
plon_out=plon
plat_out=plat
accept_out=accept_death
reject_out=reject_death


END SUBROUTINE death


SUBROUTINE move(beta,k, x, h_agg, y,n0,n0T,Rinv, Sinv, Tinv, plon, plat, regions_v, lon,lat, & 
h_v, lonmin, lonmax, latmin,latmax, sigma_clon, sigma_clat, accept_move, reject_move, it, &
burn_in, nIC, kICmax, kIC, kmax, nmeasuremax, numsites, nmeasuretotal, Ngrid,nlon,nlat, &
h_agg_out, n0_out, n0T_out, regions_v_out, plon_out, plat_out, accept_out, reject_out)



IMPLICIT NONE
! Dimensions
INTEGER kmax, nmeasuremax, Ngrid, nlon,nlat, accept_move, reject_move, it, burn_in,k, nIC, kIC, kICmax
INTEGER numsites, nmeasuretotal
! Single Variables
REAL lonmin,lonmax,latmin,latmax, sigma_clon, sigma_clat, beta 
! Input arrays
REAL x(kICmax) 
REAL h_agg(nmeasuremax,kICmax)  
REAL y(nmeasuremax) 
REAL n0(nmeasuremax) 
REAL n0T
REAL Rinv(nmeasuremax,nmeasuremax)
REAL Sinv(numsites,numsites)
REAL Tinv(nmeasuretotal,nmeasuretotal)
REAL plon(kmax)
REAL plat(kmax)
INTEGER regions_v(Ngrid)
REAL lon(nlon)
REAL lat(nlat)
REAL h_v(nmeasuremax, Ngrid)
REAL plon1m(k)
REAL plat1m(k)
REAL x1m(kIC)
REAL h_agg2m(nmeasuremax,kIC)
! Outputs
REAL h_agg_out(nmeasuremax,kICmax)
REAL n0_out(nmeasuremax) 
REAL n0T_out
INTEGER regions_v_out(Ngrid)
REAL plon_out(kmax)
REAL plat_out(kmax)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER  ri, k1, jj, ci_mv    
REAL u, n1Tm, pT_move, randomu
INTEGER regions_v1m(Ngrid)      
REAL n1m(nmeasuremax) 
REAL C(nmeasuremax)
REAL n1_arr(nmeasuretotal,numsites)
REAL C_arr(nmeasuretotal,numsites)  
! Allocatable arrays
! None

REAL,PARAMETER     :: pi = 3.14159265 

   !MOVE
   k1=k
           
   ! 1. Select new cell location - needs to be different to current locations

   x1m=x(1:kIC)               
   plon1m = plon(1:k1)
   plat1m = plat(1:k1)
   call random_number(u)   
   ci_mv = FLOOR((k1)*u) + 1
                
   ! 2. Move voronoi cell to new random location
   ! Location based on gaussian PDF about current position
    
   plon1m(ci_mv) = random_normal()*sigma_clon+plon1m(ci_mv)
   plat1m(ci_mv) = random_normal()*sigma_clat+plat1m(ci_mv)             
         
   ! Need to reject if outside of lon/lat range.
   IF (plon1m(ci_mv) .GT. lonmax .OR. plon1m(ci_mv) .LT. lonmin) THEN
       if (beta .EQ. 1. .and. it .GT. burn_in) reject_move=reject_move+1
                         
   ELSEIF (plat1m(ci_mv) .GT. latmax .OR. plat1m(ci_mv) .LT. latmin) THEN           
       if (beta .EQ. 1. .and. it .GT. burn_in) reject_move=reject_move+1 

   ELSE    
                                              
      ! 2. Recalculate voronoi cells

      call closest_grid(regions_v1m,lon,lat, plon1m,plat1m, k1, nlon,nlat)

      h_agg2m(:,:)=0.
      
      if (nIC .GT. 0) then 
          h_agg2m(:,1:nIC)=h_agg(:,1:nIC)
      endif

      do ri=1,k1
         do jj=1,Ngrid
            if (regions_v1m(jj) .EQ. ri) then
                h_agg2m(:,ri+nIC) = h_agg2m(:,ri+nIC)+h_v(:,jj)
            endif
         enddo
      enddo

 !#######################################################################
     
     n1m=matmul(h_agg2m,x1m)-y
     !n1Tm=sum((n1m/sigma_y)**2)
    
     if (numsites .EQ. 1) then
        C = matmul(n1m,Rinv)
        n1Tm= dot_product(n1m,C)
     else
        n1_arr = reshape(n1m, (/nmeasuretotal, numsites/))
        C_arr = matmul(matmul(Tinv,n1_arr),Sinv)
        C = reshape(C_arr,(/nmeasuremax/))
        n1Tm = dot_product(n1m,C)
     endif    !numsites
                                     
     ! ACCEPTANCE PROBABILITY   
     pT_move = (n1Tm - n0T)*(-0.5)*beta  

     !write(*,*) pT_move, n1Tm, n0T          
                    
     call random_number(randomu)             
     if (alog(randomu) .LE. pT_move) THEN
           !ACCEPT
           h_agg(:,:)=0.
           h_agg(:,1:kIC)=h_agg2m
           n0=n1m
           n0T=n1Tm
           regions_v=regions_v1m
           plon(:)=0.
           plat(:)=0.
           plon(1:k1)=plon1m
           plat(1:k1)=plat1m
           if (beta .EQ. 1. .and. it .GT. burn_in) accept_move=accept_move+1   

      else 
           !REJECT
           if (beta .EQ. 1. .and. it .GT. burn_in) reject_move=reject_move+1  
      endif
        
   
  ENDIF    ! plon & lat within range
     

h_agg_out=h_agg
n0_out=n0
n0T_out=n0T
regions_v_out=regions_v
plon_out=plon
plat_out=plat
accept_out=accept_move
reject_out=reject_move

END SUBROUTINE move


SUBROUTINE sigma_yt_update(beta, sigma_yt_current, sigma_yt_ap, &
detval_current, detval_T_current, detval_S, sigma_yt_hparams, stepsize_sigma_yt, sigma_yt_pdf, R_indices, &
Tinv_current, Rinv_current, Qinv, Sinv, deltatime, tau, &
n0,n0T, accept, reject, it, burn_in, nmeasuretotal, nmeasuremax, numsites, dim1, dim2, &
n0T_out, accept_out, reject_out, sigma_yt_out, detval_out, detval_T_out, Rinv_out, Tinv_out) 

IMPLICIT NONE

! Dimensions
INTEGER nmeasuretotal, accept, reject, it, burn_in, dim1, dim2
INTEGER nmeasuremax, numsites
! Single Variables
REAL beta 
! Input arrays
REAL n0(nmeasuremax) 
REAL n0T
REAL detval_current, detval_T_current, detval_S
REAL sigma_yt_current(dim2)
REAL sigma_yt_ap(dim2)
INTEGER R_indices(dim1,dim2)
REAL stepsize_sigma_yt
REAL sigma_yt_hparams(2)
INTEGER sigma_yt_pdf
REAL Tinv_current(nmeasuretotal,nmeasuretotal)
REAL Sinv(numsites,numsites)
REAL Rinv_current(nmeasuremax,nmeasuremax)
REAL Qinv(nmeasuretotal,nmeasuretotal)
REAL deltatime
REAL tau
! Outputs
REAL n0T_out, detval_out, detval_T_out
REAL sigma_yt_out(dim2)
INTEGER accept_out, reject_out
REAL Tinv_out(nmeasuretotal,nmeasuretotal)
REAL Rinv_out(nmeasuremax,nmeasuremax)
! Intermediate variables
INTEGER  yi, jj, ti
REAL randomu, dsigma_yt, sigma_yt_new, u
REAL p0_sigma_yt, p1_sigma_yt, n1T, detval_new, detval_T_new, pT, q_small
REAL y_vart_new(nmeasuretotal), y_vart_new_inv(nmeasuretotal)
REAL Tinv_new(nmeasuretotal,nmeasuretotal), Rinv_new(nmeasuremax,nmeasuremax) 
REAL C(nmeasuremax), autocorr_vec(nmeasuretotal)
REAL n0_arr(nmeasuretotal,numsites)
REAL C_arr(nmeasuretotal,numsites)  


 !do yi=1,dim2
    call random_number(u)   
    yi = FLOOR(dim2*u)+1
		
    ! Generate new value of sigma_y
    dsigma_yt = random_normal()*stepsize_sigma_yt*sigma_yt_ap(yi)
    sigma_yt_new = sigma_yt_current(yi) + dsigma_yt       

    call calc_pdf(sigma_yt_current(yi), sigma_yt_hparams(1), sigma_yt_hparams(2), sigma_yt_pdf, p0_sigma_yt)
    call calc_pdf(sigma_yt_new, sigma_yt_hparams(1), sigma_yt_hparams(2), sigma_yt_pdf, p1_sigma_yt)

	
    do jj=1,dim2   
        y_vart_new(R_indices(:,jj)) = sigma_yt_current(jj)    ! Provided dim2 isn't too big then should be fine
    enddo  

    ! Change one element of sigma_y
    y_vart_new(R_indices(:,yi)) = sigma_yt_new  ! R_indices = array of indices to which each sigma_y applies 

    y_vart_new_inv = 1./y_vart_new

    do ti=1,nmeasuretotal
        autocorr_vec = y_vart_new_inv(ti)*y_vart_new_inv*Qinv(:,ti)
        Tinv_new(:,ti) = autocorr_vec      
    enddo  

    if (numsites .GT. 1) then
        call kron(Rinv_new, Sinv, Tinv_new, numsites, nmeasuretotal,nmeasuremax)
    else
        Rinv_new = Tinv_new
    endif

    q_small = exp(-1.*deltatime/tau)
    !detval_T_new = (sum(alog(y_vart_new**2))+(alog(1-q_small**2)*(nmeasuretotal-1)))*numsites   ! DETERMINANT IN LOG SPACE - need to actually take the log

    detval_T_new = (sum(alog(y_vart_new))+(alog(1.-q_small**2)*((nmeasuretotal-1.)/2.)))*numsites   ! DETERMINANT IN LOG SPACE log of sqrt of determinant
    !detval_T_new = det(T_new)*numsites
    detval_new = detval_S + detval_T_new


    if (numsites .EQ. 1) then
       C = matmul(n0,Rinv_new)
       n1T = dot_product(n0,C)
    else
       n0_arr = reshape(n0, (/nmeasuretotal, numsites/))
       C_arr = matmul(matmul(Tinv_new,n0_arr),Sinv)
       C = reshape(C_arr,(/nmeasuremax/))
       n1T = dot_product(n0,C)
    endif

    !n1T=sum((n0/sigma_y_new)**2)
    !detval_new = sum(alog(sigma_y_new))           ! This is actually the inverse of the determinant

    ! Compute P1/P0 	
    pT= p1_sigma_yt-p0_sigma_yt - detval_new*beta + detval_current*beta - 0.5*(n1T - n0T)*beta       !*beta      ! -detval_new becasue it's inverse of true determinant

    !if (pT .NE. 0.) then
    !write(*,*) n1T, n0T, sum(Rinv_current), sum(Rinv_new)
    !endif

    if (sigma_yt_pdf .eq. 1) then
       if (sigma_yt_new .lt. sigma_yt_hparams(1)) pT = -1.e20
       if (sigma_yt_new .gt. sigma_yt_hparams(2)) pT = -1.e20
    endif

    call random_number(randomu)     ! Generates uniformly distributed random number
    
    if(alog(randomu) .le. pT) then      
       !ACCEPT	
       sigma_yt_current(yi) = sigma_yt_new
       !y_vart_current = y_vart_new
       detval_current = detval_new
       detval_T_current = detval_T_new
       Rinv_current = Rinv_new
       Tinv_current = Tinv_new
       n0T=n1T
       if(beta .eq. 1. .and. it .gt. burn_in) accept=accept + 1
    else
       !;REJECT					
       if(beta .eq. 1. .and. it .gt. burn_in) reject=reject + 1
    endif

n0T_out=n0T
sigma_yt_out=sigma_yt_current
!y_vart_out=y_vart_current
accept_out=accept
reject_out=reject
detval_out=detval_current
detval_T_out=detval_T_current
Rinv_out=Rinv_current
Tinv_out=Tinv_current
END SUBROUTINE sigma_yt_update



SUBROUTINE sigma_ys_update(beta, sigma_ys_current, sigma_ys_ap, &
detval_current, detval_S_current, detval_T,  detval_U, &
sigma_ys_hparams, stepsize_sigma_ys, sigma_ys_pdf, &
Sinv_current, Rinv_current, Tinv, Uinv, &
n0,n0T, accept, reject, it, burn_in, nmeasuretotal, nmeasuremax, numsites, &
n0T_out, accept_out, reject_out, sigma_ys_out, detval_out, detval_S_out, Rinv_out, Sinv_out) 

IMPLICIT NONE

! Dimensions
INTEGER nmeasuretotal, nmeasuremax, accept, reject, it, burn_in, numsites
! Input arrays
REAL beta
REAL n0(nmeasuremax) 
REAL n0T
REAL detval_current, detval_S_current, detval_T, detval_U
REAL sigma_ys_current(numsites)
REAL sigma_ys_ap(numsites)
REAL stepsize_sigma_ys
REAL sigma_ys_hparams(2)
INTEGER sigma_ys_pdf
REAL Tinv(nmeasuretotal,nmeasuretotal)
REAL Sinv_current(numsites,numsites)
REAL Rinv_current(nmeasuremax,nmeasuremax)
REAL Uinv(numsites,numsites)
! Outputs
REAL n0T_out, detval_out, detval_S_out
REAL sigma_ys_out(numsites)
INTEGER accept_out, reject_out
REAL Sinv_out(numsites,numsites)
REAL Rinv_out(nmeasuremax,nmeasuremax)
! Intermediate variables
INTEGER  ii, yi
REAL randomu, dsigma_ys, sigma_ys_new, u
REAL p0_sigma_ys, p1_sigma_ys, n1T, detval_new, detval_S_new, pT
REAL y_vars_new(numsites), y_vars_inv_new(numsites)
REAL Sinv_new(numsites,numsites)
REAL Rinv_new(nmeasuremax,nmeasuremax) 
REAL C(nmeasuremax)
REAL n0_arr(nmeasuretotal,numsites)
REAL C_arr(nmeasuretotal,numsites)  

    ! )nly one value for sigma_ys
    call random_number(u)   
    yi = FLOOR(numsites*u)+1
    ! Generate new value of sigma_y
    dsigma_ys = random_normal()*stepsize_sigma_ys*sigma_ys_ap(yi)
    sigma_ys_new = sigma_ys_current(yi) + dsigma_ys       

    call calc_pdf(sigma_ys_current(yi), sigma_ys_hparams(1), sigma_ys_hparams(2), sigma_ys_pdf, p0_sigma_ys)
    call calc_pdf(sigma_ys_new, sigma_ys_hparams(1), sigma_ys_hparams(2), sigma_ys_pdf, p1_sigma_ys)

    !S_new = S_current*(sigma_ys_new/sigma_ys_current)**2
    y_vars_new(:) = sigma_ys_current
    y_vars_new(yi) = sigma_ys_new
   
    y_vars_inv_new(:) = 1./sigma_ys_current
    y_vars_inv_new(yi) = 1./sigma_ys_new 

    Sinv_new(:,:)=0.

      
    do ii=1,numsites
        Sinv_new(:,ii) = y_vars_inv_new(ii)*y_vars_inv_new*Uinv(:,ii)
       ! Sinv_new(ii,ii) = 1./y_vars_new(ii)**2      ! Ignore spatital correlation for now        
    enddo

    call kron(Rinv_new, Sinv_new, Tinv, numsites, nmeasuretotal, nmeasuremax)     ! Is there a way to speed up if Sinv is diagonal?

    ! DETERMINANT IN LOG SPACE - need to actually take the log
    !detval_T_new = det(T_new)*numsites
    
    ! detval_S_new = sum(alog(y_vars_new))*nmeasuretotal      !!!! DIAGONAL  !!

    !detval_S_new = det(Sinv_new)*nmeasuretotal*(-1)
    detval_S_new = (sum(alog(y_vars_new)) + detval_U)*nmeasuretotal


    detval_new = detval_S_new + detval_T

    if (numsites .EQ. 1) then
       C = matmul(n0,Rinv_new)
       n1T = dot_product(n0,C)
    else
       n0_arr = reshape(n0, (/nmeasuretotal, numsites/))
       C_arr = matmul(matmul(Tinv,n0_arr),Sinv_new)
       C = reshape(C_arr,(/nmeasuremax/))
       n1T = dot_product(n0,C)
    endif

    ! Compute P1/P0 	
    pT= p1_sigma_ys-p0_sigma_ys - detval_new*beta + detval_current*beta - 0.5*(n1T - n0T)*beta       !*beta      ! -detval_new becasue it's inverse of true determinant

    !if (pT .NE. 0.) then
    !write(*,*) n1T, n0T, sum(Rinv_current), sum(Rinv_new)
    !endif

    if (sigma_ys_pdf .eq. 1) then
       if (sigma_ys_new .lt. sigma_ys_hparams(1)) pT = -1.e20
       if (sigma_ys_new .gt. sigma_ys_hparams(2)) pT = -1.e20
    endif

    call random_number(randomu)     ! Generates uniformly distributed random number
    
    if(alog(randomu) .le. pT) then      
       !ACCEPT	
       sigma_ys_current(yi) = sigma_ys_new
       detval_current = detval_new
       detval_S_current = detval_S_new
       Rinv_current = Rinv_new
       Sinv_current = Sinv_new
       n0T=n1T
       if(beta .eq. 1. .and. it .gt. burn_in) accept=accept + 1
    else
       !;REJECT					
       if(beta .eq. 1. .and. it .gt. burn_in) reject=reject + 1
    endif

n0T_out=n0T
sigma_ys_out=sigma_ys_current
accept_out=accept
reject_out=reject
detval_out=detval_current
detval_S_out=detval_S_current
Rinv_out=Rinv_current
Sinv_out=Sinv_current
END SUBROUTINE sigma_ys_update

 
SUBROUTINE tau_update(beta, tau_current, sigma_yt, R_indices,  &
detval_current, detval_T_current, detval_S, tau_hparam1, tau_hparam2, stepsize_tau, tau_pdf,  &
Rinv_current, Tinv_current, Qinv_current, Sinv, deltatime, &
n0, n0T, accept, reject, it, burn_in, nmeasuretotal, nmeasuremax, numsites, dim1,dim2, &
n0T_out, accept_out, reject_out, tau_out, Rinv_out, Tinv_out, Qinv_out, detval_out, detval_T_out) 

IMPLICIT NONE

! Dimensions
INTEGER nmeasuretotal, accept, reject, it, burn_in, dim1, dim2
INTEGER nmeasuremax, numsites
! Single Variables
REAL beta 
REAL tau_current
REAL stepsize_tau
REAL tau_hparam1, tau_hparam2
INTEGER tau_pdf
REAL n0T
REAL detval_current, detval_T_current, detval_S
REAL deltatime
! Input arrays
REAL n0(nmeasuremax) 
REAL Rinv_current(nmeasuremax,nmeasuremax)
REAL Tinv_current(nmeasuretotal,nmeasuretotal)
REAL Qinv_current(nmeasuretotal,nmeasuretotal)
REAL Sinv(numsites,numsites)
INTEGER R_indices(dim1,dim2)
REAL sigma_yt(dim2)
! Outputs
REAL n0T_out, detval_out, tau_out, detval_T_out
INTEGER accept_out, reject_out
REAL Rinv_out(nmeasuremax,nmeasuremax), Tinv_out(nmeasuretotal,nmeasuretotal)
REAL Qinv_out(nmeasuretotal, nmeasuretotal)
! Intermediate variables
INTEGER  ii, ti, jj
REAL randomu, dtau, tau_new
REAL p0_tau, p1_tau, n1T, detval_new, pT
REAL q_small, detval_T_new
REAL Rinv_new(nmeasuremax, nmeasuremax)
REAL Tinv_new(nmeasuretotal, nmeasuretotal)
REAL Qinv_new(nmeasuretotal, nmeasuretotal)
REAL C(nmeasuremax), autocorr_vec(nmeasuretotal)
REAL y_vart(nmeasuretotal), y_vart_inv(nmeasuretotal)
REAL n0_arr(nmeasuretotal,numsites)
REAL C_arr(nmeasuretotal,numsites)  



! Generate new value of tau
dtau = random_normal()*stepsize_tau
tau_new = tau_current + dtau 

!if(tau_new .lt. 0.) then
!    ! Skip if it generates a tau less than or equal to 0 or matrix will blow-up
!else
if (tau_new .GT. tau_hparam1 .and. tau_new .LT. tau_hparam2) THEN
!if(tau_new .gt. 0.) then
! Compute P1 for new value of tau
   call calc_pdf(tau_current, tau_hparam1, tau_hparam2, tau_pdf, p0_tau)
   call calc_pdf(tau_new, tau_hparam1, tau_hparam2, tau_pdf, p1_tau)
    		
    do jj=1,dim2   
        y_vart(R_indices(:,jj)) = sigma_yt(jj)    ! Provided dim2 isn't too big then should be fine
    enddo  


    y_vart_inv = 1./y_vart

    q_small = exp(-1.*deltatime/tau_new)

    Qinv_new(:,:) = 0.
    Qinv_new(1,1) = 1.
    Qinv_new(nmeasuretotal,nmeasuretotal)=1.
    Qinv_new(1,2) = q_small*(-1.)
    Qinv_new(nmeasuretotal, nmeasuretotal-1)=q_small*(-1.)

    do ii=2, nmeasuretotal-1
             Qinv_new(ii,ii) = 1. + q_small**2
             Qinv_new(ii,ii+1) = q_small*(-1.)
             Qinv_new(ii,ii-1) = q_small*(-1.)
    enddo
               
    Qinv_new = Qinv_new/(1.-q_small**2)
    do ti=1,nmeasuretotal
        autocorr_vec = y_vart_inv(ti)*y_vart_inv*Qinv_new(:,ti)
        Tinv_new(:,ti) = autocorr_vec      
    enddo
	
    if (numsites .GT. 1) then
        call kron(Rinv_new, Sinv, Tinv_new, numsites, nmeasuretotal, nmeasuremax)
    else
        Rinv_new = Tinv_new
    endif

    !detval_T_new = (sum(alog(y_vart**2))+(alog(1-q_small**2)*(nmeasuretotal-1)))*numsites   ! DETERMINANT IN LOG SPACE - need to actually take the log

     detval_T_new = (sum(alog(y_vart))+(alog(1-q_small**2)*((nmeasuretotal-1)/2.)))*numsites   ! DETERMINANT IN LOG SPACE log of sqrt of determinant

    detval_new = detval_S + detval_T_new


    if (numsites .EQ. 1) then
       C = matmul(n0,Rinv_new)
       n1T = dot_product(n0,C)
    else
       n0_arr = reshape(n0, (/nmeasuretotal, numsites/))
       C_arr = matmul(matmul(Tinv_new,n0_arr),Sinv)
       C = reshape(C_arr,(/nmeasuremax/))
       n1T = dot_product(n0,C)
    endif
	
    ! compute P1/P0 
    pT=p1_tau - p0_tau - detval_new*beta + detval_current*beta - 0.5*(n1T - n0T)*beta    ! This way round or should detvals be opposite? 


    !write(*,*) pT, p1_tau, p0_tau, detval_new, detval_current, n1T, n0T

    call random_number(randomu)        ! Generates uniformly distributed random number
 
    if(alog(randomu) .le. pT) then      
       !;ACCEPT
       tau_current = tau_new
       Qinv_current = Qinv_new
       Tinv_current = Tinv_new
       Rinv_current= Rinv_new
       detval_current = detval_new
       detval_T_current = detval_T_new
       n0T=n1T
       if(beta .eq. 1. .and. it .gt. burn_in) accept=accept + 1
    else 
       !;REJECT	
        if(beta .eq. 1. .and. it .gt. burn_in) reject=reject + 1
    endif
else
   if(beta .eq. 1. .and. it .gt. burn_in) reject=reject + 1
endif

tau_out=tau_current
detval_out = detval_current
detval_T_out = detval_T_current
Qinv_out=Qinv_current
Tinv_out=Tinv_current
Rinv_out=Rinv_current
n0T_out=n0T
accept_out=accept
reject_out=reject

END SUBROUTINE tau_update

SUBROUTINE rho_nu_update(beta, rho_current, nu_current, sigma_ys,  &
detval_current, detval_T, detval_S_current, detval_U_current, &
rho_hparam1, rho_hparam2, stepsize_rho, rho_pdf,  &
nu_hparam1, nu_hparam2, stepsize_nu, nu_pdf,  &
Rinv_current, Tinv, Sinv_current, Uinv_current, distance, &
n0, n0T, accept, reject, it, burn_in, nmeasuretotal, nmeasuremax, numsites, &
n0T_out, accept_out, reject_out, rho_out, nu_out, Rinv_out, Sinv_out, Uinv_out, detval_out, detval_S_out, detval_U_out) 

IMPLICIT NONE

! Dimensions
INTEGER nmeasuretotal, accept, reject, it, burn_in
INTEGER nmeasuremax, numsites
! Single Variables
REAL beta 
REAL rho_current, nu_current
REAL stepsize_rho, stepsize_nu
REAL rho_hparam1, rho_hparam2, nu_hparam1, nu_hparam2
INTEGER rho_pdf, nu_pdf
REAL n0T
REAL detval_current, detval_T, detval_S_current, detval_U_current
! Input arrays
REAL n0(nmeasuremax) 
REAL Rinv_current(nmeasuremax,nmeasuremax)
REAL Tinv(nmeasuretotal,nmeasuretotal)
REAL Uinv_current(numsites,numsites)
REAL Sinv_current(numsites,numsites)
REAL sigma_ys(numsites)
REAL distance(numsites,numsites)
! Outputs
REAL n0T_out, detval_out, rho_out, nu_out, detval_S_out, detval_U_out
INTEGER accept_out, reject_out
REAL Rinv_out(nmeasuremax,nmeasuremax), Sinv_out(numsites,numsites)
REAL Uinv_out(numsites, numsites)
! Intermediate variables
INTEGER  ii, ti, jj
REAL randomu, n1T, pT
REAL drho, rho_new, p0_rho, p1_rho
REAL dnu, nu_new, p0_nu, p1_nu
REAL detval_S_new, detval_U_new, detval_new
REAL Rinv_new(nmeasuremax, nmeasuremax)
REAL Sinv_new(numsites, numsites)
REAL Uinv_new(numsites, numsites), U_new(numsites,numsites)
REAL C(nmeasuremax), autocorr_vec(numsites)
REAL y_vars_inv(numsites)
REAL ga
REAL n0_arr(nmeasuretotal,numsites)
REAL C_arr(nmeasuretotal,numsites)  
real (kind =8)  alpha, arg_dum, k_arg_dum
real (kind =8) arg(numsites,numsites)
integer ( kind = 4 ) nb, ize, ncalc

! Generate new value of rho
drho = random_normal()*stepsize_rho
rho_new = rho_current + drho 

! Generate new value of nu
dnu = random_normal()*stepsize_nu
nu_new = nu_current + dnu 

U_new(:,:)=0.

if (rho_new .GT. rho_hparam1 .and. rho_new .LT. rho_hparam2 &
    .and. nu_new .GE. nu_hparam1 .and. nu_new .LT. nu_hparam2) THEN

! Compute P1 for new value of rho
   call calc_pdf(rho_current, rho_hparam1, rho_hparam2, rho_pdf, p0_rho)
   call calc_pdf(rho_new, rho_hparam1, rho_hparam2, rho_pdf, p1_rho)

   call calc_pdf(nu_current, nu_hparam1, nu_hparam2, nu_pdf, p0_nu)
   call calc_pdf(nu_new, nu_hparam1, nu_hparam2, nu_pdf, p1_nu)
    		
    y_vars_inv = 1./sigma_ys
    
    nb = int(nu_new)
    alpha = nu_new - nb
    nb = nb+1
    arg = sqrt(2*nu_new)*distance/rho_new
    ize = 2
    ga= gamma(nu_new)
    

    do ii = 1,numsites 
         do jj = 1,numsites
              arg_dum = arg(ii,jj)
            !  if (arg_dum==0.) then
            !       U_new(ii,jj) = 1.
            !  else
            !       call rkbesl(arg_dum, alpha, nb, ize, k_arg_dum, ncalc)
            !       k_arg_dum = k_arg_dum*exp(-1*arg(ii,jj))
            !       U_new(ii,jj) = 1/(2**(nu_new-1)*ga)*(sqrt(2*nu_new)*distance(ii,jj)/rho_new)**nu_new*k_arg_dum
                   U_new(ii,jj) = exp(-1.*distance(ii,jj)/rho_new)
            !  endif
         enddo
    enddo

       
    Uinv_new = Ainv(U_new)

    do ti=1,numsites
        autocorr_vec = y_vars_inv(ti)*y_vars_inv*Uinv_new(:,ti)
        Sinv_new(:,ti) = autocorr_vec      
    enddo
	
   
    call kron(Rinv_new, Sinv_new, Tinv, numsites, nmeasuretotal, nmeasuremax)
    
    detval_U_new = det(U_new)
   
    detval_S_new = (sum(alog(sigma_ys)) + detval_U_new)*nmeasuretotal
   
    detval_new = detval_S_new + detval_T
		
    if (numsites .EQ. 1) then		
       C = matmul(n0,Rinv_new)
       n1T = dot_product(n0,C)
    else
       n0_arr = reshape(n0, (/nmeasuretotal, numsites/))
       C_arr = matmul(matmul(Tinv,n0_arr),Sinv_new)
       C = reshape(C_arr,(/nmeasuremax/))
       n1T = dot_product(n0,C)
    endif	
	
    ! compute P1/P0 
    pT= p1_rho - p0_rho  + p1_nu - p0_nu - detval_new*beta + detval_current*beta - 0.5*(n1T - n0T)*beta    ! This way round or should detvals be opposite? 


    call random_number(randomu)        ! Generates uniformly distributed random number
 
    if(alog(randomu) .le. pT) then      
       !;ACCEPT
       rho_current = rho_new
       nu_current = nu_new
       Uinv_current = Uinv_new
       Sinv_current = Sinv_new
       Rinv_current = Rinv_new
       detval_current = detval_new
       detval_S_current = detval_S_new
       detval_U_current = detval_U_new
       n0T=n1T
       if(beta .eq. 1. .and. it .gt. burn_in) accept=accept + 1
    else 
       !;REJECT	
        if(beta .eq. 1. .and. it .gt. burn_in) reject=reject + 1
    endif
else
   if(beta .eq. 1. .and. it .gt. burn_in) reject=reject + 1
endif

rho_out=rho_current
nu_out = nu_current
detval_out = detval_current
detval_S_out = detval_S_current
detval_U_out=detval_U_current
Uinv_out=Uinv_current
Sinv_out=Sinv_current
Rinv_out=Rinv_current
n0T_out=n0T
accept_out=accept
reject_out=reject

END SUBROUTINE rho_nu_update





SUBROUTINE y_update(beta, y_current, z, n0, n0T, m0T, Rinv, sigma_measure, accept_y, reject_y, &
stepsize_y, timeindex_nonzero, nmeasuremax, nmeasure, it, burn_in,  &
n0T_out, m0T_out, y_out, n0_out, accept_y_out, reject_y_out)

Implicit None
! Dimensions
INTEGER nmeasure, nmeasuremax, it, burn_in
! Inputs
REAL beta 
REAL y_current(nmeasuremax)
REAL z(nmeasure)
REAL n0(nmeasuremax)
REAL Rinv(nmeasuremax,nmeasuremax)
REAL sigma_measure(nmeasure)
REAL n0T, m0T
REAL stepsize_y
INTEGER timeindex_nonzero(nmeasure)
INTEGER accept_y, reject_y
! Outputs
REAL n0T_out, m0T_out
INTEGER accept_y_out, reject_y_out
REAL y_out(nmeasuremax), n0_out(nmeasuremax)
! Intermediates
INTEGER  yi
REAL dymod, m1T, n1T
REAL dum, dumb,dumc,dumd
REAL pT, randomu
REAL y_new(nmeasuremax), n1(nmeasuremax)
REAL y_small(nmeasure), m1(nmeasure)
REAL dumx(nmeasuremax), dumy(nmeasuremax)

do yi=1,nmeasuremax

   y_new = y_current
	
   ! Generate new value of y
   dymod = random_normal()*stepsize_y*y_current(yi)
   y_new(yi) = y_current(yi) + dymod

   y_small = y_new(timeindex_nonzero)
   !y_small = y_current(timeindex_nonzero)

   m1 = z - y_small
 
   !do ii=1,nmeasure
   !   dum2(ii) = m1(ii)*Dinv(ii)*m1(ii)
   !enddo

   m1T = sum((m1/sigma_measure)**2)

   !m0T = sum(Dinv*m0**2)
   !m1T = sum(Dinv*m1**2)

   !m1T = sum(dum2)

   n1 = n0
   n1(yi) = n0(yi) - dymod

   dumx = n0(yi)*Rinv(:,yi)
   dumy = n1(yi)*Rinv(:,yi)

   dum  = dot_product(n0,dumx)  
   dumb = dot_product(n1,dumy)
   dumc = dot_product(n0,Rinv(:,yi))
   dumd = dot_product(n1,Rinv(:,yi))

   n1T = n0T - dum + dumb - n0(yi)*dumc + n1(yi)*dumd - n1(yi)*n1(yi)*Rinv(yi,yi) + n0(yi)*n0(yi)*Rinv(yi,yi)

   ! Compute P1/P0 
   pT= -0.5*(m1T - m0T)*beta -0.5*(n1T - n0T)*beta 

   call random_number(randomu)       !Generates uniformly distributed random number

   if(alog(randomu) .le. pT) then      
      !;ACCEPT
      y_current(yi) = y_new(yi)
      n0=n1
      n0T=n1T
      m0T = m1T
      if(it .gt. burn_in) accept_y=accept_y+ 1
   else
      !;REJECT
      if(it .gt. burn_in) reject_y=reject_y+ 1
   endif
enddo 

accept_y_out=accept_y
reject_y_out=reject_y
y_out = y_current
n0_out = n0
n0T_out =n0T
m0T_out = m0T
End subroutine y_update

FUNCTION random_normal() RESULT(fn_val)
IMPLICIT NONE
! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

REAL :: fn_val

!     Local variables
REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
            r1 = 0.27597, r2 = 0.27846, half =0.5, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
RETURN

END FUNCTION random_normal

SUBROUTINE closest_grid(region_v,lon,lat, plon,plat, np, nlon,nlat)
implicit none

integer :: region_v(nlon*nlat)
real    :: lon(nlon), lat(nlat), plon(np), plat(np)
integer :: lai, loi, la,lo,pi
integer :: nlon,nlat,np
real    :: maxdist, dist
    lai=1
    do la =1,nlat
    !do lo =1,nlon
        loi=1
        do lo=1, nlon
            maxdist=1.e6
            do pi=1,np 
                dist=(lat(la) - plat(pi))*(lat(la) - plat(pi)) + (lon(lo) - plon(pi))*(lon(lo) - plon(pi))  
                ! Can't work out which way round??????????????
                if (dist .LT. maxdist) THEN
                    !!region(lai, loi)=pi
                    region_v(loi+(lai-1)*nlon)=pi
                    !region_v(lai+(loi-1)*nlat)=pi
                    maxdist=dist
                endif
            enddo
            loi=loi+1
        enddo
        lai=lai+1
    enddo

END SUBROUTINE closest_grid


SUBROUTINE closest_grid_sng(region,lon,lat, plon,plat, np)
implicit none

  integer :: region, pi, np
  real    :: plon(np), plat(np)
  real    :: lon, lat, maxdist, dist

  maxdist=1.e6
  do pi=1,np 
     dist=(lat - plat(pi))*(lat - plat(pi)) + (lon - plon(pi))*(lon - plon(pi))
     if (dist .LT. maxdist) THEN  
        region=pi
        maxdist=dist
     endif
  enddo

END SUBROUTINE closest_grid_sng

subroutine calc_pdf(x,pdf_param1,pdf_param2,pdf,p1)     ! Subroutine for calculating P1 for specified PDF

 Implicit none
 ! Args
 integer,intent(in)   :: pdf
 real,intent(in)      :: x, pdf_param1, pdf_param2
 real,intent(out)     :: p1
 integer              :: positive
 real                 :: mu, sigma
 ! Locals
real,parameter        :: pi = 3.14159265        ! a fixed constant

!! ALL PDFS ARE IN LOG-SPACE
!  PDF .eq. 'UNIFORM'
  If(pdf .eq. 1) then                 ! Uniform PDF = 1	

     if (x .lt. pdf_param1 .OR. x .gt. pdf_param2) then
        p1 = -1.e20
     else
        p1= alog(1./(pdf_param2 - pdf_param1))
     endif
     positive=0
  Endif

!  PDF .eq. 'GAUSSIAN'
  If(pdf .eq. 2) then                ! Gaussian PDF = 2
     p1=alog(1./(pdf_param2*sqrt(2.*pi)))+(-1.*(x - pdf_param1)**2 / 2./(pdf_param2**2))
     positive=0
  Endif

!  PDF .eq. 'LOGNORMAL'
  If(pdf .eq. 3) then                    ! Lognormal PDF = 3

    mu = alog(pdf_param1) - 0.5*alog(1. + pdf_param2**2/pdf_param1**2)
    sigma=sqrt(alog((pdf_param2/pdf_param1)**2 + 1.))


     p1=alog(1./(x*sigma*sqrt(2.*pi)))+( -1.*(alog(x) - mu)**2 / 2./(sigma**2))
     positive=1
  Endif

!  PDF .eq. 'EXPONENTIAL'
  If(pdf .eq. 4) then                    ! Exponential PDF = 4
     p1=alog(pdf_param1)+(-1.*pdf_param1*x)
     positive=1
  Endif


  if((x .le. 0.) .and. (positive .eq. 1)) p1 = -1.e20       ! Since log(0) isn't defined need to make p1 a very large negative number to ensure proposal is rejected

end subroutine calc_pdf

subroutine kron(K, A, B, x, y, z)
! Calculate the Kronecker product of two matrices
    IMPLICIT NONE
    REAL, INTENT(IN)  :: A(x,x), B(y,y)
    REAL, INTENT(INOUT) :: K(z,z)
    INTEGER, INTENT(IN) :: x,y,z
    INTEGER :: I, J, MA, NA, MB, NB
    MA = UBOUND(A, 1)
    NA = UBOUND(A, 2)
    MB = UBOUND(B, 1)
    NB = UBOUND(B, 2)

    IF (SIZE(K,1) /= MA*MB .OR. SIZE(K,2) /= NA*NB) THEN
        WRITE(*,*) 'K has invalid size'
        CALL ABORT
    END IF
    FORALL(I=1:MA, J=1:NA)
        K(MB*(I-1)+1:MB*I,NB*(J-1)+1:NB*J) = A(I,J)*B
    END FORALL
END subroutine kron

subroutine init_random_seed()
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed

  call random_seed(size = n)
  allocate(seed(n))

  call system_clock(count=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)

  deallocate(seed)
end subroutine

function Ainv(A) 
  real, dimension(:,:), intent(in) :: A
  real, dimension(size(A,1),size(A,2)) :: Ainv

  real, dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external SGETRF
  external SGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call SGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'In Ainv: Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call SGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function Ainv


function det(A)
  !Returns the log of the "square root of the determinant" of a postive definite matrix calculated from the Cholesky decomposition
  ! Uses LAPACK routines
  real, dimension(:,:), intent(in) :: A
  real, dimension(size(A,1),size(A,2)) :: Adum
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  real :: det
  integer :: n, info,j

  ! External procedures defined in LAPACK
  external SPOTRF

  ! Store A in Adum to prevent it from being overwritten by LAPACK
  Adum = A
  n = size(A,1)

  call SPOTRF('U',n,Adum,n,info)
  det = 0

  do j = 1,n
  det = det + alog(Adum(j,j))
  enddo

end function det

subroutine rkbesl ( x, alpha, nb, ize, k_arg, ncalc )

!*****************************************************************************80
!
!! RKBESL calculates K Bessel function with non-integer orders.
!
!  Discussion:
!
!    This routine calculates modified Bessel functions of the second
!    kind, K SUB(N+ALPHA) (X), for non-negative argument X, and
!    non-negative order N+ALPHA, with or without exponential scaling.
!
!    This program is based on a program written by J. B. Campbell
!    that computes values of the Bessel functions K of real
!    argument and real order.  Modifications include the addition
!    of non-scaled functions, parameterization of machine
!    dependencies, and the use of more accurate approximations
!    for SINH and SIN.
!
!    In case of an error, NCALC .NE. NB, and not all K's are
!    calculated to the desired accuracy.
!
!    NCALC < -1:  An argument is out of range. For example,
!    NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
!    XMAX.  In this case, the B-vector is not calculated,
!    and NCALC is set to MIN0(NB,0)-2  so that NCALC .NE. NB.
!
!    NCALC = -1:  Either  K(ALPHA,X) .GE. XINF  or
!    K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) .GE. XINF.  In this case,
!    the B-vector is not calculated.  Note that again
!    NCALC .NE. NB.
!
!    0 < NCALC < NB: Not all requested function values could
!    be calculated accurately.  BK(I) contains correct function
!    values for I <= NCALC, and contains the ratios
!    K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!    22 January 2014 A. Ganesan - added bk(1:nb) = 0 (ln 396) after  'ncalc = min ( nb, 0 ) - 2'
!    16 January 2014 A. Ganesan - only the final order of bk is output k_arg = bk(nb)
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    JB Campbell,
!    On Temme's Algorithm for the Modified Bessel Functions of the
!    Third Kind,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 581-586.
!
!    JB Campbell,
!    A FORTRAN IV Subroutine for the Modified Bessel Functions of
!    the Third Kind of Real Order and Real Argument,
!    Report NRC/ERB-925,
!    National Research Council, Canada.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the non-negative argument for which
!    K's or exponentially scaled K's (K*EXP(X))
!    are to be calculated.  If K's are to be calculated,
!    X must not be greater than XMAX.
!
!    Input, real ( kind = 8 ) ALPHA, the fractional part of order for which
!    K's or exponentially scaled K's (K*EXP(X)) are to be calculated.
!    0 <= ALPHA < 1.0.
!
!    Input, integer ( kind = 4 ) NB, the number of functions to be calculated, 
!    NB .GT. 0.  The first function calculated is of order ALPHA, and the
!    last is of order (NB - 1 + ALPHA).
!
!    Input, integer ( kind = 4 ) IZE, scaling option.
!    1, unscaled functions are to calculated,
!    2, exponentially scaled functions are to be calculated.
!
!    Output, real ( kind = 8 ) BK(NB), the results.  If the routine
!    terminates normally, with NCALC = NB, the vector BK contains the
!    functions K(ALPHA,X), ... , K(NB-1+ALPHA,X), or the corresponding
!    exponentially scaled functions.
!    If (0 < NCALC < NB), BK(I) contains correct function
!    values for I <= NCALC, and contains the ratios
!    K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
!
!    Output, integer ( kind = 4 ) NCALC, error indicator.  If NCALC = NB, then 
!    all the requested values were calculated to the desired accuracy.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) blpha
  real ( kind = 8 ) bk(1)
  real ( kind = 8 ) bk1
  real ( kind = 8 ) bk2
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) dm
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) d3
  real ( kind = 8 ) enu
  real ( kind = 8 ) estf(7)
  real ( kind = 8 ) estm(6)
  real ( kind = 8 ) ex
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iend
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) ize
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mplus1
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) ncalc
  real ( kind = 8 ) p(8)
  real ( kind = 8 ) p0
  real ( kind = 8 ) q(7)
  real ( kind = 8 ) q0
  real ( kind = 8 ) r(5)
  real ( kind = 8 ) ratio
  real ( kind = 8 ) s(4)
  real ( kind = 8 ) sqxmin
  real ( kind = 8 ) t(6)
  real ( kind = 8 ) tinyx
  real ( kind = 8 ) twonu
  real ( kind = 8 ) twox
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) wminf
  real ( kind = 8 ) x
  real ( kind = 8 ) xinf
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) x2by4
  real ( kind = 8 ) k_arg
!
!  Mathematical constants
!    A = LOG(2.D0) - Euler's constant
!    D = SQRT(2.D0/PI)
!
  data tinyx / 1.0d-10/
  data a / 0.11593151565841244881d0/
  data d /0.797884560802865364d0/
!
!  Machine dependent parameters
!
  data sqxmin / 1.49d-154 /
  data xinf / 1.79d+308 /
  data xmin / 2.23d-308 /
  data xmax / 705.342d0 /
!
!  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA
!                                         + Euler's constant
!  Coefficients converted from hex to decimal and modified
!  by W. J. Cody, 2/26/82
!
!  R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA)
!  T    - Approximation for SINH(Y)/Y
!
  data p/ 0.805629875690432845d00,    0.204045500205365151d02, &
          0.157705605106676174d03,    0.536671116469207504d03, &
          0.900382759291288778d03,    0.730923886650660393d03, &
          0.229299301509425145d03,    0.822467033424113231d00/
  data q/ 0.294601986247850434d02,    0.277577868510221208d03, &
          0.120670325591027438d04,    0.276291444159791519d04, &
          0.344374050506564618d04,    0.221063190113378647d04, &
          0.572267338359892221d03/
  data r/-0.48672575865218401848d+0,  0.13079485869097804016d+2, &
         -0.10196490580880537526d+3,  0.34765409106507813131d+3, &
          0.34958981245219347820d-3/
  data s/-0.25579105509976461286d+2,  0.21257260432226544008d+3, &
         -0.61069018684944109624d+3,  0.42269668805777760407d+3/
  data t/ 0.16125990452916363814d-9, 0.25051878502858255354d-7, &
          0.27557319615147964774d-5, 0.19841269840928373686d-3, &
          0.83333333333334751799d-2, 0.16666666666666666446d+0/
  data estm / 5.20583d1, 5.7607d0, 2.7782d0, 1.44303d1, 1.853004d2, &
            9.3715d0/
  data estf / 4.18341d1, 7.1075d0, 6.4306d0, 4.25110d1, 1.35633d0, &
            8.45096d1, 2.0d1/

  ex = x
  enu = alpha
  ncalc = min ( nb, 0 ) - 2
  !bk(1:nb) = 0

  if ( 0 < nb .and. &
    ( 0.0D+00 <= enu .and. enu < 1.0D+00 ) .and. &
    ( 1 <= ize .and. ize <= 2 ) .and. &
    ( ize /= 1 .or. ex <= xmax ) .and. &
    0.0D+00 < ex )  then

    k = 0
    if ( enu < sqxmin ) then
      enu = 0.0D+00
    end if

    if ( 0.5D+00 < enu ) then
      k = 1
      enu = enu - 1.0D+00
    end if

    twonu = enu + enu
    iend = nb + k - 1
    c = enu * enu
    d3 = -c

    if ( ex <= 1.0D+00 ) then
!
!  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA,
!                 Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA.
!
      d1 = 0.0D+00
      d2 = p(1)
      t1 = 1.0D+00
      t2 = q(1)

      do i = 2, 7, 2
        d1 = c * d1 + p(i)
        d2 = c * d2 + p(i+1)
        t1 = c * t1 + q(i)
        t2 = c * t2 + q(i+1)
      end do

      d1 = enu * d1
      t1 = enu * t1
      f1 = log ( ex )
      f0 = a + enu * ( p(8) &
         - enu * ( d1 + d2 ) / ( t1 + t2 ) ) - f1
      q0 = exp ( -enu * ( a - enu * &
         ( p(8) + enu * ( d1 - d2 ) / ( t1 - t2 ) ) - f1 ) )
      f1 = enu * f0
      p0 = exp ( f1 )
!
!  Calculation of F0.
!
      d1 = r(5)
      t1 = 1.0D+00
      do i = 1, 4
        d1 = c * d1 + r(i)
        t1 = c * t1 + s(i)
      end do

      if ( abs ( f1 ) <= 0.5D+00 ) then
        f1 = f1 * f1
        d2 = 0.0D+00
        do i = 1, 6
          d2 = f1 * d2 + t(i)
        end do
        d2 = f0 + f0 * f1 * d2
      else
        d2 = sinh ( f1 ) / enu
      end if

      f0 = d2 - enu * d1 / ( t1 * p0 )
!
!  X <= 1.0E-10.
!
!  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X).
!
      if ( ex <= tinyx ) then

        bk(1) = f0 + ex * f0

        if ( ize == 1 ) then
          bk(1) = bk(1) - ex * bk(1)
        end if

        ratio = p0 / f0
        c = ex * xinf
!
!  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X),
!  1/2 <= ALPHA.
!
        if ( k /= 0 ) then

          ncalc = -1

          if ( c / ratio <= bk(1) ) then
!	  k_arg = bk(nb)
            return
          end if

          bk(1) = ratio * bk(1) / ex
          twonu = twonu + 2.0D+00
          ratio = twonu

        end if

        ncalc = 1

        if ( nb == 1 ) then
          k_arg = bk(nb)
          return
        end if
!
!  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),  L  =  1, 2, ... , NB-1.
!
        ncalc = -1
        do i = 2, nb
          if ( c <= ratio ) then
!	  k_arg = bk(nb)
            return
          end if
          bk(i) = ratio / ex
          twonu = twonu + 2.0D+00
          ratio = twonu
        end do

        ncalc = 1
        j = ncalc + 1

        do i = j, nb
          if ( xinf / bk(i) <= bk(ncalc) ) then
            k_arg = bk(nb)
            return
          end if
          bk(i) = bk(ncalc) * bk(i)
          ncalc = i
        end do
        k_arg = bk(nb)
        return
!
!  1.0E-10 < X <= 1.0.
!
      else

        c = 1.0D+00
        x2by4 = ex * ex / 4.0D+00
        p0 = 0.5D+00 * p0
        q0 = 0.5D+00 * q0
        d1 = - 1.0D+00
        d2 = 0.0D+00
        bk1 = 0.0D+00
        bk2 = 0.0D+00
        f1 = f0
        f2 = p0

  100       continue

        d1 = d1 + 2.0D+00
        d2 = d2 + 1.0D+00
        d3 = d1 + d3
        c = x2by4 * c / d2
        f0 = ( d2 * f0 + p0 + q0 ) / d3
        p0 = p0 / ( d2 - enu )
        q0 = q0 / ( d2 + enu )
        t1 = c * f0
        t2 = c * ( p0 - d2 * f0 )
        bk1 = bk1 + t1
        bk2 = bk2 + t2

        if ( epsilon ( t1 ) < abs ( t1 / ( f1 + bk1 ) ) .or. &
             epsilon ( t2 ) < abs ( t2 / ( f2 + bk2 ) ) )  then
          go to 100
        end if

        bk1 = f1 + bk1
        bk2 = 2.0D+00 * ( f2 + bk2 ) / ex

        if ( ize == 2 ) then
          d1 = exp ( ex )
          bk1 = bk1 * d1
          bk2 = bk2 * d1
        end if

        wminf = estf(1) * ex + estf(2)

      end if
!
!  1/EPS < X.
!
    else if ( 1.0D+00 < epsilon ( ex ) * ex ) then

      ncalc = nb
      bk1 = 1.0D+00 / ( d * sqrt ( ex ) )
      do i = 1, nb
        bk(i) = bk1
      end do
      k_arg = bk(nb)
      return

    else
!
!  1 < X.
!
      twox = ex + ex
      blpha = 0.0D+00
      ratio = 0.0D+00

      if ( ex <= 4.0D+00 ) then
!
!  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 <= X <= 4.0.
!
        d2 = aint ( estm(1) / ex + estm(2) )
        m = int ( d2 )
        d1 = d2 + d2
        d2 = d2 - 0.5D+00
        d2 = d2 * d2
        do i = 2, m
          d1 = d1 - 2.0D+00
          d2 = d2 - d1
          ratio = ( d3 + d2 ) / ( twox + d1 - ratio )
        end do
!
!  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
!  recurrence and K(ALPHA,X) from the Wronskian.
!
        d2 = aint ( estm(3) * ex + estm(4) )
        m = int ( d2 )
        c = abs ( enu )
        d3 = c + c
        d1 = d3 - 1.0D+00
        f1 = xmin
        f0 = ( 2.0D+00 * ( c + d2 ) / ex &
           + 0.5D+00 * ex / ( c + d2 + 1.0D+00 ) ) * xmin

        do i = 3, m
          d2 = d2 - 1.0D+00
          f2 = ( d3 + d2 + d2 ) * f0
          blpha = ( 1.0D+00 + d1 / d2 ) * ( f2 + blpha )
          f2 = f2 / ex + f1
          f1 = f0
          f0 = f2
        end do

        f1 = ( d3 + 2.0D+00 ) * f0 / ex + f1
        d1 = 0.0D+00
        t1 = 1.0D+00
        do i = 1, 7
          d1 = c * d1 + p(i)
          t1 = c * t1 + q(i)
        end do

        p0 = exp ( c * ( a + c * ( p(8) &
           - c * d1 / t1 ) - log ( ex ) ) ) / ex
        f2 = ( c + 0.5D+00 - ratio ) * f1 / ex
        bk1 = p0 + ( d3 * f0 - f2 + f0 + blpha ) &
          / ( f2 + f1 + f0 ) * p0

        if ( ize == 1 ) then
          bk1 = bk1 * exp ( - ex )
        end if

        wminf = estf(3) * ex + estf(4)

      else
!
!  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward
!  recurrence, for 4 < X.
!
        dm = aint ( estm(5) / ex + estm(6) )
        m = int ( dm )
        d2 = dm - 0.5D+00
        d2 = d2 * d2
        d1 = dm + dm

        do i = 2, m
          dm = dm - 1.0D+00
          d1 = d1 - 2.0D+00
          d2 = d2 - d1
          ratio = ( d3 + d2 ) / ( twox + d1 - ratio )
          blpha = ( ratio + ratio * blpha ) / dm
        end do

        bk1 = 1.0D+00 / ( ( d + d * blpha ) * sqrt ( ex ) )

        if ( ize == 1 ) then
          bk1 = bk1 * exp ( - ex )
        end if

        wminf = estf(5) * ( ex - abs ( ex - estf(7) ) ) + estf(6)

      end if
!
!  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and
!  K(ALPHA+1,X)/K(ALPHA,X).
!
      bk2 = bk1 + bk1 * ( enu + 0.5D+00 - ratio ) / ex

    end if
!
!  Calculation of 'NCALC', K(ALPHA+I,X), I  =  0, 1, ... , NCALC-1,
!  K(ALPHA+I,X)/K(ALPHA+I-1,X), I  =  NCALC, NCALC+1, ... , NB-1.
!
    ncalc = nb
    bk(1) = bk1

    if ( iend == 0 ) then
      k_arg = bk(nb)
      return
    end if

    j = 2 - k

    if ( 0 < j ) then
      bk(j) = bk2
    end if

    if ( iend == 1 ) then
      k_arg = bk(nb)
      return
    end if

    m = min ( int ( wminf - enu ), iend )

    do i = 2, m

      t1 = bk1
      bk1 = bk2
      twonu = twonu + 2.0D+00

      if ( ex < 1.0D+00 ) then

        if ( ( xinf / twonu ) * ex <= bk1 ) then
          exit
        end if

      else

        if ( xinf / twonu <= bk1 / ex ) then
          exit
        end if

      end if

      bk2 = twonu / ex * bk1 + t1
      itemp = i
      j = j + 1

      if ( 0 < j ) then
        bk(j) = bk2
      end if

    end do

    m = itemp

    if ( m == iend ) then
      k_arg = bk(nb)
      return
    end if

    ratio = bk2 / bk1
    mplus1 = m + 1
    ncalc = -1

    do i = mplus1, iend

      twonu = twonu + 2.0D+00
      ratio = twonu / ex + 1.0D+00 / ratio
      j = j + 1

      if ( 1 < j ) then
        bk(j) = ratio
      else
        if ( xinf / ratio <= bk2 ) then
           k_arg = bk(nb)
          return
        end if
        bk2 = ratio * bk2
      end if

    end do

    ncalc = max ( mplus1 - k, 1 )

    if ( ncalc == 1 ) then
      bk(1) = bk2
    end if

    if ( nb == 1 ) then
      k_arg = bk(nb)
      return
    end if

    j = ncalc + 1

    do i = j, nb
      if ( xinf / bk(i) <= bk(ncalc) ) then
        k_arg = bk(nb)
        return
      end if
      bk(i) = bk(ncalc) * bk(i)
      ncalc = i
    end do

  end if
  k_arg = bk(nb)
  return
end subroutine rkbesl


End module transd_inv

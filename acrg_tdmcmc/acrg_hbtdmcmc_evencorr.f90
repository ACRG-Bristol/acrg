
Module transd_evencorr

contains


SUBROUTINE hbtdmcmc(beta,k, x, h_agg,y,n0, plon, plat, regions_v, &
pdf_param1, pdf_param2, lon,lat, h_v, sigma_model, sigma_measure, error_structure, &
R_indices, sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y, sigma_model_pdf, &
tau, tau_hparams, stepsize_tau, tau_pdf, deltatime, &
y_hparam1, y_hparam2, stepsize_y, y_pdf, timeindex_zero, &
sigma_clon, sigma_clat, rjmcmc, para_temp, nmeasure_site, nsite_max, & 
lonmin, lonmax, latmin,latmax, sigma_bd, kmin, x_pdf_all, burn_in, &
pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, pdf_param1_pdf, pdf_param2_pdf, &
stepsize, stepsize_pdf_p1,stepsize_pdf_p2, nIt, nsub, nit_sub, nIC, &
nbeta, kmax, kICmax, nmeasure, Ngrid, nlon,nlat, ydim1, ydim2, numsites,nIC1, nzero,  &
k_out, x_out, regions_out, plon_out, plat_out, sigma_y_out, sigma_model_out, n0T_out, &
pdf_param1_out, pdf_param2_out, tau_out, y_it,accept, reject, &
accept_birth, reject_birth, accept_death, reject_death, accept_move, reject_move, &
accept_swap, reject_swap, accept_sigma_y, reject_sigma_y, accept_tau, reject_tau, &
accept_y, reject_y,tot_acc_x, tot_acc_p1, tot_acc_p2, tot_acc_sigma_y, tot_acc_tau, &
accept_all, reject_all, accept_birth_all, reject_birth_all, &
accept_death_all, reject_death_all, accept_move_all, reject_move_all, &
accept_sigma_y_all, reject_sigma_y_all, accept_tau_all, reject_tau_all)


IMPLICIT NONE

!! INPUTS !!!!!!!!!!!!!!!!!!!
! Dimensions
INTEGER nbeta
INTEGER kmax
INTEGER kICmax
INTEGER nmeasure
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
INTEGER numsites
INTEGER nzero
! Single Variables
REAL lonmin
REAL lonmax
REAL latmin
REAL latmax
REAL sigma_bd
REAL sigma_clon
REAL sigma_clat
REAL stepsize_sigma_y(ydim2)
REAL stepsize_tau
INTEGER sigma_model_pdf
INTEGER pdf_param1_pdf
INTEGER pdf_param2_pdf
INTEGER tau_pdf
INTEGER rjmcmc
INTEGER para_temp
INTEGER nsite_max
INTEGER y_pdf
REAL stepsize_y
REAL y_hparam1
REAL y_hparam2
! Input arrays
REAL stepsize(nIC1)
REAL stepsize_pdf_p1(nIC1)
REAL stepsize_pdf_p2(nIC1)
INTEGER x_pdf_all(nIC1)
REAL beta(nbeta)
INTEGER k(nbeta)
REAL x(kICmax, nbeta)               
REAL pdf_param1(kICmax,nbeta)            
REAL pdf_param2(kICmax,nbeta)                
REAL pdf_p1_hparam1(nIC1)
REAL pdf_p1_hparam2(nIC1)
REAL pdf_p2_hparam1(nIC1)
REAL pdf_p2_hparam2(nIC1)
REAL h_agg(nmeasure,kICmax,nbeta)
REAL y(nmeasure) 
REAL n0(nmeasure, nbeta) 
REAL sigma_model(ydim2, nbeta)
INTEGER R_indices(ydim1,ydim2)
REAL sigma_measure(nmeasure)
REAL error_structure(nmeasure) 
REAl sigma_model_hparam1(ydim2)
REAl sigma_model_hparam2(ydim2)
REAL plon(kmax,nbeta)
REAL plat(kmax,nbeta)
INTEGER regions_v(Ngrid,nbeta)
REAL lon(nlon)
REAL lat(nlat)
REAL h_v(nmeasure, Ngrid)
REAL tau(numsites,nbeta)
REAL tau_hparams(2)
REAL deltatime
INTEGER nmeasure_site(numsites)
INTEGER timeindex_zero(nzero)
! Outputs
INTEGER k_out(nit_sub)
REAL x_out(kICmax,nit_sub)  
REAL pdf_param1_out(kICmax,nit_sub)             
REAL pdf_param2_out(kICmax,nit_sub)   
REAL plon_out(kmax,nit_sub)
REAL plat_out(kmax,nit_sub)
INTEGER regions_out(Ngrid,nit_sub)
REAL sigma_y_out(nmeasure, nit_sub)
REAL sigma_model_out(ydim2, nit_sub)
REAL n0T_out(nit_sub)
REAL tau_out(numsites,nit_sub)
INTEGER accept(nIC1)
INTEGER reject(nIC1)
INTEGER accept_all(nIC1,nbeta)
INTEGER reject_all(nIC1,nbeta)
REAL tot_acc_x(nIC1)
REAL tot_acc_p1(nIC1)
REAL tot_acc_p2(nIC1)
REAL tot_acc_sigma_y(ydim2)
REAL tot_acc_tau
INTEGER accept_birth_all(nbeta), reject_birth_all(nbeta)
INTEGER accept_death_all(nbeta), reject_death_all(nbeta)
INTEGER accept_move_all(nbeta), reject_move_all(nbeta)
INTEGER accept_sigma_y_all(nbeta), reject_sigma_y_all(nbeta)
INTEGER accept_tau_all(nbeta), reject_tau_all(nbeta)
INTEGER accept_birth, reject_birth
INTEGER accept_death, reject_death, accept_move, reject_move, accept_swap
INTEGER accept_sigma_y, reject_sigma_y, accept_tau, reject_tau, reject_swap
INTEGER accept_y, reject_y
REAL y_it(nmeasure,nit_sub)
! INTERMEDIATE VARIABLES
INTEGER it, ibeta, remain_it, pair1,pair2, ib, it_sub, remain, kIC       !remain_dim
INTEGER remain_swap
REAL u, u1,u2, randomu,pT_chain, beta1,beta2
INTEGER k_it(nit_sub)
REAL x_it(kICmax,nit_sub)                             
REAL pdf_param1_it(kICmax,nit_sub)   
REAL pdf_param2_it(kICmax,nit_sub)   
REAL plon_it(kmax,nit_sub)               
REAL plat_it(kmax,nit_sub)              
REAL sigma_model_ap(ydim2)            
INTEGER regions_it(Ngrid,nit_sub)
REAL sigma_y_it(nmeasure,nit_sub)
REAL sigma_model_it(ydim2, nit_sub) 
REAL n0T_it(nit_sub)
REAL tau_it(numsites,nit_sub)
REAL Rinv(nmeasure,nmeasure,nbeta), Qinv(nmeasure,nmeasure,nbeta)
REAL detval(nbeta), detval_Q(nbeta), n0T(nbeta)
REAL sigma_y(nmeasure,nbeta)
REAL detval_Q_block(numsites,nbeta)
INTEGER acc_h_batch(nIC1)
INTEGER rej_h_batch(nIC1)
INTEGER accept_batch(nIC1)
INTEGER reject_batch(nIC1)
INTEGER acc_y_batch(ydim2)
INTEGER rej_y_batch(ydim2)
INTEGER acc_tau_batch
INTEGER rej_tau_batch
! INTERMEDIATE TEMP VARIABLES
INTEGER ti,jj
REAL Rinv_temp(nmeasure,nmeasure), Qinv_temp(nmeasure,nmeasure)
!REAL Q_temp(nmeasure,nmeasure)
REAL detval_temp, detval_Q_temp, n0T_temp
REAL n0_temp(nmeasure), sigma_y_temp(nmeasure), y_error(nmeasure)
REAL sigma_yinv(nmeasure), C(nmeasure)
! SUBROUTINE INPUTS
REAL betaib, n0Tib, tauib(numsites)
REAL xib(kICmax), plonib(kmax), platib(kmax), n0ib(nmeasure)           
REAL pdf_param1ib(kICmax), pdf_param2ib(kICmax)
REAL h_aggib(nmeasure,kICmax)
REAL sigma_yib(nmeasure), sigma_modelib(ydim2)
REAL Rinvib(nmeasure,nmeasure), Qinvib(nmeasure,nmeasure)
INTEGER kib
INTEGER regions_vib(Ngrid)
! SUBROUTINE OUTPUTS
REAL n0Tib1, tauib1(numsites)                     
REAL xib1(kICmax), plonib1(kmax), platib1(kmax), n0ib1(nmeasure)      
REAL pdf_param1ib1(kICmax), pdf_param2ib1(kICmax)
REAL h_aggib1(nmeasure,kICmax)
REAL sigma_yib1(nmeasure), sigma_modelib1(ydim2)
INTEGER kib1, rejectib1, acceptib1, reject_yib1, accept_yib1
INTEGER acceptxib1(nIC1), rejectxib1(nIC1), acc_bxib1(nIC1), rej_bxib1(nIC1)
INTEGER acc_byib1(ydim2), rej_byib1(ydim2)
INTEGER acc_tauib1, rej_tauib1
INTEGER regions_vib1(Ngrid)
REAL detvalib, detvalib1, detval_Qib, detval_Qib1
REAL detval_Q_blockib(numsites), detval_Q_blockib1(numsites)
REAL Rinvib1(nmeasure,nmeasure), Qinvib1(nmeasure,nmeasure)
REAL y_out(nmeasure)
REAL stepsize_sig_ib1(ydim2)
REAL stepsize_tau_ib1
REAL stepsize_ib1(nIC1)     
REAL stepsize_p1_ib1(nIC1)
INTEGER acc_prob_p1_ib1(nIC1)
REAL stepsize_p2_ib1(nIC1)
INTEGER rej_prob_p1_ib1(nIC1)
! BLOCK INTERMEDIATES
INTEGER cum_nmeas, si, ii
!REAL Q_block_temp(nsite_max,nsite_max)
!REAL Q_block_inv_temp(nsite_max,nsite_max)
REAL detval_Q_block_temp(numsites)
!REAL n0T_block(numsites)
REAL aa, bb, q_small


! OTHER INTERMEDIATES

!f2py intent(in) beta,k, x, h_agg,y,n0, plon, plat, regions_v
!f2py intent(in) pdf_param1, pdf_param2, lon,lat, h_v, sigma_clon, sigma_clat  
!f2py intent(in) tau, tau_hparams, stepsize_tau, tau_pdf, deltatime, error_structure 
!f2py intent(in) sigma_model, sigma_measure, rjmcmc, para_temp, nmeasure_site, nsite_max
!f2py intent(in) R_indices, sigma_model_hparam1, sigma_model_hparam2 ,stepsize_sigma_y
!f2py intent(in) sigma_model_pdf, lonmin, lonmax, latmin,latmax, sigma_bd, kmin, x_pdf_all, burn_in
!f2py intent(in) pdf_p2_hparam1, pdf_p2_hparam2, stepsize_pdf_p2, pdf_param2_pdf
!f2py intent(in) pdf_p1_hparam1, pdf_p1_hparam2, stepsize_pdf_p1, pdf_param1_pdf
!f2py intent(in) stepsize, nIt,nsub,nit_sub, nIC1, numsites
!f2py intent(in) nIC, nbeta, kmax, kICmax, nmeasure, Ngrid, nlon,nlat, ydim1, ydim2
!f2py intent(in) y_hparam1,y_hparam2,y_pdf,stepsize_y,timeindex_zero, nzero
!f2py intent(out) k_out, x_out, regions_out, plon_out, plat_out, sigma_y_out, n0T_out
!f2py intent(out) pdf_param2_out, pdf_param1_out, sigma_model_out
!f2py intent(out) accept, reject, accept_swap, accept_birth, accept_death, accept_move
!f2py intent(out) reject_birth, reject_death, reject_move, accept_sigma_y, reject_sigma_y
!f2py intent(out) reject_swap, accept_tau, reject_tau, tau_out, y_it, accept_y, reject_y
!f2py intent(out) tot_acc_x, tot_acc_p1, tot_acc_p2, tot_acc_sigma_y, tot_acc_tau
!f2py intent(out) accept_all, reject_all, accept_birth_all, reject_birth_all
!f2py intent(out) accept_death_all, reject_death_all, accept_move_all, reject_move_all
!f2py intent(out) accept_sigma_y_all, reject_sigma_y_all, accept_tau_all, reject_tau_all


  !call OMP_SET_NUM_THREADS(nbeta)      ! Uncomment for Parallel Tempering

 call init_random_seed()          ! Ensure random number generation starts from new point each time program is run
                                  ! Random seed only needs to be called once in a program.  
				    ! It is used to generate the first random number. 
				   ! Any random numbers generated after the first will use the previous random number as its seed.
aa = 1
bb = 0
accept(:)=0
accept_birth=0
accept_death=0
accept_move=0
reject(:)=0
reject_birth=0
reject_death=0
reject_move=0
accept_sigma_y=0
reject_sigma_y=0
accept_swap=0
reject_swap=0
accept_tau=0
reject_tau=0
it_sub=1

accept_all(:,:)=0
accept_birth_all(:)=0
accept_death_all(:)=0
accept_move_all(:)=0
reject_all(:,:)=0
reject_birth_all(:)=0
reject_death_all(:)=0
reject_move_all(:)=0
accept_sigma_y_all(:)=0
reject_sigma_y_all(:)=0
accept_tau_all(:)=0
reject_tau_all(:)=0

acc_h_batch(:)=0.
rej_h_batch(:)=0.
accept_batch(:)=0
reject_batch(:) = 0
acc_y_batch(:)=0
rej_y_batch(:)=0
acc_tau_batch=0
rej_tau_batch=0


 cum_nmeas=0
 Qinv_temp(:,:)=0.
sigma_model_ap=sigma_model(:,1)

  do jj=1,ydim2   
        !y_error(R_indices(:,jj)) = sigma_model_ap(jj)
        y_error(R_indices(:,jj)) = sigma_model_ap(jj)*error_structure(R_indices(:,jj))    ! Provided dim2 isn't too big then should be fine
  enddo  

  sigma_y_temp=sqrt(y_error**2 + sigma_measure**2)   
  sigma_yinv = 1./sigma_y_temp

  !Q_temp = exp((-1.)*deltatime/tau(1,1))
  q_small = exp((-1.)*deltatime/tau(1,1))   ! Assumes all tau start from the same value

  do si=1,numsites

    !Q_block_temp(:,:)=0.
    !Q_block_inv_temp(:,:)=0.

    !Q_block_temp(1:nmeasure_site(si),1:nmeasure_site(si)) = &
    !Q_temp(cum_nmeas+1:cum_nmeas+nmeasure_site(si), cum_nmeas+1:cum_nmeas+nmeasure_site(si))
    
    !Q_block_inv_temp(1:nmeasure_site(si),1:nmeasure_site(si)) = Ainv(Q_block_temp(1:nmeasure_site(si),1:nmeasure_site(si)))

    !detval_Q_block_temp(si) = det(Q_block_temp(1:nmeasure_site(si),1:nmeasure_site(si)))

    !Qinv_temp(cum_nmeas+1:cum_nmeas+nmeasure_site(si), cum_nmeas+1:cum_nmeas+nmeasure_site(si)) = &
    !Q_block_inv_temp(1:nmeasure_site(si),1:nmeasure_site(si))

    !cum_nmeas = cum_nmeas + nmeasure_site(si)

    Qinv_temp(cum_nmeas+1,cum_nmeas+1) = 1.
    Qinv_temp(cum_nmeas+nmeasure_site(si),cum_nmeas+nmeasure_site(si)) = 1.
    Qinv_temp(cum_nmeas+1,cum_nmeas+2) = q_small*(-1.)
    Qinv_temp(cum_nmeas+nmeasure_site(si),cum_nmeas+nmeasure_site(si)-1) = q_small*(-1.)

    do ii=2, nmeasure_site(si)-1
         Qinv_temp(cum_nmeas+ii,cum_nmeas+ii) = 1. + q_small**2
         Qinv_temp(cum_nmeas+ii,cum_nmeas+ii+1) = q_small*(-1.)
         Qinv_temp(cum_nmeas+ii,cum_nmeas+ii-1) = q_small*(-1.)
    enddo
   
    detval_Q_block_temp(si) = alog(1-q_small**2)*((nmeasure_site(si)-1)/2.)   ! DETERMINANT IN LOG SPACE log of sqrt of determinant
   
    cum_nmeas = cum_nmeas + nmeasure_site(si)

  enddo

  Qinv_temp = Qinv_temp/(1.-q_small**2)   ! Assumes all taus start at smae value


  !Qinv_temp = Ainv(Q_temp)
  !detval_Q_temp = det(Q_temp)

  !write(*,*) detval_Q_block
  !stop

  detval_Q_temp = sum(detval_Q_block_temp)

  do ti=1,nmeasure
       Rinv_temp(:,ti) = sigma_yinv(ti)*sigma_yinv*Qinv_temp(:,ti)
  enddo
	
  detval_temp =  sum(alog(sigma_y_temp)) + detval_Q_temp  ! DETERMINANT IN LOG SPACE log of sqrt of determinant

  n0_temp = n0(:,1)

!  cum_nmeas = 0
!  do si=1, numsites
!
!  C_temp(cum_nmeas+1:cum_nmeas+nmeasure_site(si)) = matmul(n0_temp(cum_nmeas+1:cum_nmeas+nmeasure_site(si)), &
!           Rinv_temp(cum_nmeas+1:cum_nmeas+nmeasure_site(si), cum_nmeas+1:cum_nmeas+nmeasure_site(si)))
!
!  n0T_block(si) = dot_product(n0_temp(cum_nmeas+1:cum_nmeas+nmeasure_site(si)),C_temp(cum_nmeas+1:cum_nmeas+nmeasure_site(si)))
!
!  cum_nmeas = cum_nmeas + nmeasure_site(si)
!  enddo


  C = matmul(n0_temp,Rinv_temp)
  n0T_temp = dot_product(n0_temp,C)

  !call sgemm ('N', 'N', nmeasure, 1, nmeasure, aa, Rinv_temp, nmeasure, n0_temp, nmeasure, bb, C, nmeasure)
  !n0T_temp = dot_product(n0_temp,C)

  !write(*,*) n0T_temp, n0T_temp2
  !stop

  do ibeta = 1,nbeta
     sigma_y(:,ibeta) = sigma_y_temp
     Rinv(:,:,ibeta)=Rinv_temp
     Qinv(:,:,ibeta)=Qinv_temp
     detval(ibeta)=detval_temp
     detval_Q(ibeta)=detval_Q_temp
     detval_Q_block(:,ibeta)=detval_Q_block_temp
     n0T(ibeta) = n0T_temp
  enddo


 !write(*,*) detval_Q_temp, detval_temp
 !stop

! MCMC loop
!###############################################################
do it=1,(nIt+burn_in)
   
  ! call random_number(u) 

   if (rjmcmc .EQ. 1) then
       remain_it= modulo(it,7)+1
       !remain_it = FLOOR(7*u) + 1    ! Choose random number between 1 and 7 to choose what to update
   else
  
       remain_it= modulo(it,4)+1
       !if (modulo(it,8) .EQ. 5) then
       !    remain_it=3
       !elseif (modulo(it,8) .GE. 6) then
       !    remain_it=2
       !else
       !    remain_it=1                ! Weight probability in favour of x_update. Severely deweight tau update since it's slow
       !endif
       !remain_it = FLOOR(3*u) + 1    ! Choose random number between 1 and 2 - no reversible jump.
   endif


!$OMP PARALLEL DO DEFAULT(SHARED) private(ibeta, betaib,kib,xib,pdf_param1ib,pdf_param2ib, plonib, platib), &
!$OMP& private(regions_vib,h_aggib,n0ib,n0Tib,kIC,xib1,n0ib1,n0Tib1,acceptib1,rejectib1,regions_vib1), &
!$OMP& private(u, kib1, h_aggib1, plonib1, platib1, acceptxib1, rejectxib1), &
!$OMP& private(sigma_yib, sigma_modelib, sigma_yib1, sigma_modelib1, accept_yib1, reject_yib1), &
!$OMP& private(pdf_param1ib1, pdf_param2ib1, detvalib, detvalib1, detval_Qib, detval_Qib1), &
!$OMP& private(Rinvib, Qinvib, tauib, Rinvib1, Qinvib1, tauib1), &
!$OMP& private(detval_Q_blockib, detval_Q_blockib1), & 
!$OMP& private(stepsize_ib1, acc_bxib1, rej_bxib1),&
!$OMP& private(stepsize_p1_ib1, acc_prob_p1_ib1,rej_prob_p1_ib1),&
!$OMP& private(stepsize_p2_ib1, stepsize_sig_ib1, acc_byib1, rej_byib1),& 
!$OMP& private(stepsize_tau_ib1, acc_tauib1, rej_tauib1),&             
!$OMP& shared(x,n0,n0T, k, pdf_param1, pdf_param2, h_agg, plon,plat, regions_v)
   do ibeta=1,nbeta

     if (para_temp .EQ. 1 .or. ibeta .EQ. 1) then 
      
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

       sigma_yib = sigma_y(:,ibeta)
       sigma_modelib = sigma_model(:,ibeta)
       Rinvib = Rinv(:,:,ibeta)
       Qinvib = Qinv(:,:,ibeta)

       detvalib = detval(ibeta)
       detval_Qib = detval_Q(ibeta)
       detval_Q_blockib = detval_Q_block(:,ibeta)
       tauib = tau(:,ibeta)

       kIC = kib+nIC
       
       
       if (remain_it .EQ. 1) then              ! X UPDATE


            call x_hparam_update(betaib,kib, xib, pdf_param1ib,pdf_param2ib, &
                 pdf_p1_hparam1, pdf_p1_hparam2, stepsize_pdf_p1, pdf_param1_pdf, &
                 pdf_p2_hparam1, pdf_p2_hparam2, stepsize_pdf_p2, pdf_param2_pdf, &
                 acc_h_batch, rej_h_batch, x_pdf_all, it, burn_in, nIC, kICmax, nIC1, &
                 pdf_param1ib1, pdf_param2ib1, stepsize_p1_ib1, stepsize_p2_ib1, acc_prob_p1_ib1, rej_prob_p1_ib1) 

            call x_update(betaib,kib, xib, pdf_param1ib1,pdf_param2ib1, &
                         h_aggib,n0ib,n0Tib,Rinvib, stepsize, &
                         accept_batch, reject_batch, x_pdf_all, it, burn_in, nIC, kICmax, nmeasure, nIC1, &
                         xib1, n0ib1, n0Tib1, acceptxib1, rejectxib1, stepsize_ib1, acc_bxib1, rej_bxib1)


            x(:,ibeta) = xib1
            n0(:,ibeta) = n0ib1
            n0T(ibeta) = n0Tib1 
            pdf_param1(:,ibeta) = pdf_param1ib1
            pdf_param2(:,ibeta) = pdf_param2ib1
            
            accept_all(:,ibeta) = accept_all(:,ibeta) + acceptxib1
            reject_all(:,ibeta) = reject_all(:,ibeta) + rejectxib1
            if (betaib .EQ. 1.) then 
               accept(:) = accept(:) + acceptxib1
               reject(:) = reject(:) + rejectxib1
               stepsize=stepsize_ib1
               stepsize_pdf_p1=stepsize_p1_ib1
               stepsize_pdf_p2=stepsize_p2_ib1
               acc_h_batch=acc_prob_p1_ib1
               rej_h_batch=rej_prob_p1_ib1
               accept_batch=acc_bxib1
               reject_batch=rej_bxib1
            endif
            
           

       elseif (remain_it .EQ. 5) then       ! BIRTH
               
               call birth(betaib,kib, xib, h_aggib,y,n0ib,n0Tib,Rinvib, plonib, platib, regions_vib, lon,lat, & 
                          h_v, pdf_param1ib, pdf_param2ib, x_pdf_all(nIC1), &
                          sigma_bd,it,burn_in,nIC,kICmax,kmax, &
                          nmeasure,Ngrid,nlon,nlat, &
                          kib1, xib1, h_aggib1, n0ib1, n0Tib1, regions_vib1, plonib1, platib1, acceptib1, rejectib1, &
                          pdf_param1ib1, pdf_param2ib1)

                k(ibeta) = kib1
                x(:,ibeta) = xib1
                plon(:,ibeta) = plonib1
                plat(:,ibeta) = platib1
                regions_v(:,ibeta) = regions_vib1
                h_agg(:,:,ibeta) = h_aggib1
                n0(:,ibeta) = n0ib1
                n0T(ibeta) = n0Tib1
                pdf_param1(:,ibeta) = pdf_param1ib1
                pdf_param2(:,ibeta) = pdf_param2ib1

                accept_birth_all(ibeta)= accept_birth_all(ibeta) + acceptib1
                reject_birth_all(ibeta)= reject_birth_all(ibeta) + rejectib1
                if (betaib .EQ. 1.) then 
                   accept_birth= accept_birth + acceptib1
                   reject_birth= reject_birth + rejectib1
                endif

           elseif (remain_it .EQ. 6) then    ! DEATH

               call death(betaib,kib, xib, h_aggib, y,n0ib,n0Tib,Rinvib, plonib, platib, regions_vib, lon,lat, & 
                          h_v, pdf_param1ib, pdf_param2ib, x_pdf_all(nIC1), sigma_bd, &
                          it, burn_in,nIC, kICmax, kmin, kmax, nmeasure, &
                          Ngrid,nlon,nlat, &
                          kib1, xib1, h_aggib1,n0ib1, n0Tib1, regions_vib1, plonib1, platib1, acceptib1, rejectib1,&
                          pdf_param1ib1, pdf_param2ib1)
                     
               k(ibeta) = kib1
               x(:,ibeta) = xib1
               plon(:,ibeta) = plonib1
               plat(:,ibeta) = platib1
               regions_v(:,ibeta) = regions_vib1
               h_agg(:,:,ibeta) = h_aggib1
               n0(:,ibeta) = n0ib1
               n0T(ibeta) = n0Tib1
               pdf_param1(:,ibeta) = pdf_param1ib1
               pdf_param2(:,ibeta) = pdf_param2ib1
          
               accept_death_all(ibeta)= accept_death_all(ibeta) + acceptib1
               reject_death_all(ibeta)= reject_death_all(ibeta) + rejectib1
               if (betaib .EQ. 1.) then 
                   accept_death= accept_death + acceptib1
                   reject_death= reject_death + rejectib1
               endif

           elseif (remain_it .EQ. 7) then    ! MOVE
                
               
               call move(betaib,kib, xib, h_aggib, y,n0ib,n0Tib,Rinvib, plonib, platib, regions_vib, lon,lat, & 
                         h_v, lonmin, lonmax, latmin,latmax, sigma_clon, sigma_clat, it, &
                         burn_in, nIC, kICmax, kIC, kmax, nmeasure, Ngrid,nlon,nlat, &
                         h_aggib1, n0ib1, n0Tib1, regions_vib1, plonib1, platib1, acceptib1, rejectib1)
               
               plon(:,ibeta) = plonib1
               plat(:,ibeta) = platib1
               regions_v(:,ibeta) = regions_vib1
               h_agg(:,:,ibeta) = h_aggib1
               n0(:,ibeta) = n0ib1
               n0T(ibeta) = n0Tib1
               
               accept_move_all(ibeta)= accept_move_all(ibeta) + acceptib1
               reject_move_all(ibeta)= reject_move_all(ibeta) + rejectib1
               if (betaib .EQ. 1.) then 
                  accept_move= accept_move + acceptib1
                  reject_move= reject_move + rejectib1
               endif

          elseif (remain_it .EQ. 2) then  ! SIGMA_Y UPDATE
             
              call sigma_y_update(betaib, sigma_modelib, sigma_measure, sigma_yib, error_structure, detvalib, &
                 detval_Qib, Rinvib, Qinvib, sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y, &
                 sigma_model_pdf, R_indices, &
                 n0ib,n0Tib, acc_y_batch, rej_y_batch, it, burn_in, nmeasure, ydim1, ydim2, &
                 n0Tib1, accept_yib1, reject_yib1, sigma_yib1, sigma_modelib1, Rinvib1, detvalib1, &
                 stepsize_sig_ib1, acc_byib1, rej_byib1) 
                      
              n0T(ibeta) = n0Tib1
              sigma_y(:,ibeta) = sigma_yib1
              sigma_model(:,ibeta) = sigma_modelib1
              Rinv(:,:,ibeta) = Rinvib1 
              detval(ibeta) = detvalib1 

              accept_sigma_y_all(ibeta) = accept_sigma_y_all(ibeta) + accept_yib1
              reject_sigma_y_all(ibeta) = reject_sigma_y_all(ibeta) + reject_yib1
              if (betaib .EQ. 1.) then 
               accept_sigma_y = accept_sigma_y + accept_yib1
               reject_sigma_y = reject_sigma_y + reject_yib1
               stepsize_sigma_y=stepsize_sig_ib1
               acc_y_batch=acc_byib1
               rej_y_batch=rej_byib1
              endif


           elseif (remain_it .EQ. 3) then  ! TAU_UPDATE

            ! if (betaib .EQ. 1.) then 
              call tau_update(betaib, tauib, sigma_yib, detvalib, detval_Qib, detval_Q_blockib, &
                              tau_hparams(1), tau_hparams(2), stepsize_tau, tau_pdf,  &
                              Rinvib, Qinvib, deltatime, nmeasure_site, n0ib, n0Tib, &
                              acc_tau_batch, rej_tau_batch, it, burn_in, nmeasure, numsites, &
                              n0Tib1, acceptib1, rejectib1, tauib1, Rinvib1, Qinvib1, detvalib1, &
                              detval_Qib1, detval_Q_blockib1, stepsize_tau_ib1, acc_tauib1, rej_tauib1) 

              n0T(ibeta) = n0Tib1
              tau(:,ibeta) = tauib1
              Rinv(:,:,ibeta) = Rinvib1 
              Qinv(:,:,ibeta) = Qinvib1 
              detval(ibeta) = detvalib1 
              detval_Q(ibeta) = detval_Qib1 
              detval_Q_block(:,ibeta) = detval_Q_blockib1 

              accept_tau_all(ibeta) = accept_tau_all(ibeta) + acceptib1
              reject_tau_all(ibeta) = reject_tau_all(ibeta) + rejectib1
              if (betaib .EQ. 1.) then 
                accept_tau = accept_tau + acceptib1
                reject_tau = reject_tau + rejectib1
                stepsize_tau=stepsize_tau_ib1
                acc_tau_batch=acc_tauib1
                rej_tau_batch=rej_tauib1
              endif

           elseif (remain_it .EQ. 4) then  ! Y_UPDATE

              if (timeindex_zero(1) .NE. -999) then

                if (betaib .EQ. 1.) then
                  call y_update(betaib, y, n0ib, n0Tib, &
                            Rinvib, y_hparam1, y_hparam2, timeindex_zero, &
                            stepsize_y, y_pdf, accept_y, reject_y, it, burn_in, nmeasure, nzero, &
                            n0ib1, n0Tib1, y_out, acceptib1, rejectib1) 

             
                 n0T(ibeta) = n0Tib1
                 n0(:,ibeta) = n0ib1
                 y = y_out 
             
                 !if (betaib .EQ. 1.) then 
                   accept_y = acceptib1
                   reject_y = rejectib1
                 endif
              endif            ! timeindex_zero .NE. -999

           endif     ! remain_it
        
       endif     !para_temp .EQ. 1 .or. ibeta .EQ. 1)    

   enddo    ! beta loop  UNCOMMENT IF DOING PT
!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Store xit and swap betas
  
   IF (it .GT. burn_in/2) THEN      ! Begin swaps after half of burn-in time
        remain_swap = modulo(it,2)
        IF (remain_swap .EQ. 1) THEN
          if (para_temp .EQ. 1) then
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
            pT_chain = (beta2-beta1)*(n0T(pair2)/2.-n0T(pair1)/2.+detval(pair2)-detval(pair1))  ! detvals should be inverse determinants so signs this way round
            call random_number(randomu)
            if (alog(randomu) .LE. pT_chain) then
                beta(pair2)=beta1*1.                     ! UNCOMMENT IF DOING PT
                beta(pair1)=beta2*1.                     ! UNCOMMENT IF DOING PT
                accept_swap=accept_swap+1
            else
                reject_swap=reject_swap+1
            endif      ! pT_chain if      
           else
                reject_swap=reject_swap+1
           endif   ! para_temp=1  

         ENDIF      ! reamin_it =0 if
   ENDIF          ! it > burn_in/2
    


   IF (it .GT. burn_in) THEN     
        remain = modulo(it,nsub)          ! nsub typically = 100
        if (remain .EQ. 0) then
      !            ib=1
            do ib=1,nbeta                             ! UNCOMMENT IF DOING PT
               if (beta(ib) .EQ. 1.) then             ! UNCOMMENT IF DOING PT
                  x_it(:,it_sub)=x(:,ib)
                  plon_it(:,it_sub)=plon(:,ib)
                  plat_it(:,it_sub)=plat(:,ib)
                  k_it(it_sub)=k(ib)
                  regions_it(:,it_sub)=regions_v(:,ib)
                  sigma_y_it(:,it_sub)=sigma_y(:,ib)
                  sigma_model_it(:,it_sub)=sigma_model(:,ib)
                  n0T_it(it_sub)=n0T(ib)
                  tau_it(:,it_sub) = tau(:,ib)
                  pdf_param1_it(:,it_sub)=pdf_param1(:,ib)
                  pdf_param2_it(:,it_sub)=pdf_param2(:,ib) 
                  y_it(:,it_sub)=y      
                  it_sub=it_sub+1
               endif                                  ! UNCOMMENT IF DOING PT
            enddo                                     ! UNCOMMENT IF DOING PT
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
sigma_y_out=sigma_y_it
sigma_model_out=sigma_model_it
n0T_out=n0T_it
tau_out=tau_it
pdf_param1_out=pdf_param1_it
pdf_param2_out=pdf_param2_it
tot_acc_sigma_y=stepsize_sigma_y
tot_acc_x=stepsize
tot_acc_p1=stepsize_pdf_p1
tot_acc_p2=stepsize_pdf_p2
tot_acc_tau=stepsize_tau
END SUBROUTINE hbtdmcmc



SUBROUTINE x_hparam_update(beta,k, x, pdf_param1_all,pdf_param2_all, &
pdf_p1_hparam1_all, pdf_p1_hparam2_all, stepsize_pdf_p1, pdf_param1_pdf, &
pdf_p2_hparam1_all, pdf_p2_hparam2_all, stepsize_pdf_p2, pdf_param2_pdf, &
accept_batch, reject_batch, x_pdf_all, it, burn_in, nIC, kICmax, nIC1, &
pdf_param1_out, pdf_param2_out, stepsize_p1_out, stepsize_p2_out, accept_batch_out, reject_batch_out) 

Implicit none
INTEGER it, burn_in, k, nIC, kICmax, nIC1
REAL av_acc,beta
REAL x(kICmax) 
INTEGER x_pdf_all(nIC1)
INTEGER accept_batch(nIC1), reject_batch(nIC1)
REAL pdf_param1_all(kICmax), pdf_param2_all(kICmax)
REAL pdf_p1_hparam1_all(nIC1), pdf_p1_hparam2_all(nIC1) 
REAL  pdf_p2_hparam1_all(nIC1), pdf_p2_hparam2_all(nIC1)
REAL pdf_param1_out(kICmax), pdf_param2_out(kICmax)
REAL stepsize_p1_out(nIC1), stepsize_p2_out(nIC1)
INTEGER accept_batch_out(nIC1), reject_batch_out(nIC1)
INTEGER elem(5)
REAL u(5)
REAL pdf_param1        
REAL pdf_param2
REAL accep_prob(nIC1)
REAL stepsize_pdf_p2(nIC1), stepsize_pdf_p1(nIC1)
INTEGER xi, x_pdf, xx
REAL pT, randomu, p0,p1
INTEGER pdf_param2_pdf, pdf_param1_pdf
REAL pdf_p2_hparam1, pdf_p2_hparam2, pdf_p1_hparam1, pdf_p1_hparam2
REAL dpdf_param2, pdf_param2_new, dpdf_param1, pdf_param1_new
REAL p0_pdf_param2, p1_pdf_param2, p0_pdf_param1, p1_pdf_param1
REAL stepsize_pdf_p10, stepsize_pdf_p20

   call random_number(u)   
    
  elem = FLOOR(k*u)+1+nIC

  do xx=1,nIC+5
  
  if (xx .LE. nIC) then
     xi = xx
  else
     xi=elem(xx-nIC)
  endif

  if (xi .LE. nIC) then
     x_pdf = x_pdf_all(xi)
     pdf_param1 = pdf_param1_all(xi)
     pdf_param2 = pdf_param2_all(xi)
     pdf_p1_hparam1 = pdf_p1_hparam1_all(xi)
     pdf_p1_hparam2 = pdf_p1_hparam2_all(xi)
     pdf_p2_hparam1 = pdf_p2_hparam1_all(xi)
     pdf_p2_hparam2 = pdf_p2_hparam2_all(xi)
     stepsize_pdf_p10 = stepsize_pdf_p1(xi)
     stepsize_pdf_p20 = stepsize_pdf_p2(xi)
  else if (xi .GT. nIC) then
     x_pdf = x_pdf_all(nIC1)
     pdf_param1 = pdf_param1_all(xi)
     pdf_param2 = pdf_param2_all(xi)
     pdf_p1_hparam1 = pdf_p1_hparam1_all(nIC1)
     pdf_p1_hparam2 = pdf_p1_hparam2_all(nIC1)
     pdf_p2_hparam1 = pdf_p2_hparam1_all(nIC1)
     pdf_p2_hparam2 = pdf_p2_hparam2_all(nIC1) 
     stepsize_pdf_p10 = stepsize_pdf_p1(nIC1)
     stepsize_pdf_p20 = stepsize_pdf_p2(nIC1)
  endif
  dpdf_param1 = random_normal()*stepsize_pdf_p10
  pdf_param1_new = pdf_param1 + dpdf_param1

  dpdf_param2 = random_normal()*stepsize_pdf_p20
  pdf_param2_new = pdf_param2 + dpdf_param2
  

         call calc_pdf(pdf_param1,pdf_p1_hparam1,pdf_p1_hparam2,pdf_param1_pdf, p0_pdf_param1) 
         call  calc_pdf(pdf_param1_new, pdf_p1_hparam1,pdf_p1_hparam2,pdf_param1_pdf, p1_pdf_param1)

         call calc_pdf(pdf_param2,pdf_p2_hparam1,pdf_p2_hparam2,pdf_param2_pdf, p0_pdf_param2) 
         call  calc_pdf(pdf_param2_new, pdf_p2_hparam1,pdf_p2_hparam2,pdf_param2_pdf, p1_pdf_param2)

      
         call calc_pdf(x(xi),pdf_param1,pdf_param2,x_pdf, p0)           ! Will apply whatever the PDF
         call calc_pdf(x(xi),pdf_param1_new,pdf_param2_new,x_pdf, p1)        
               
         pT = p1+p1_pdf_param1+p1_pdf_param2-p0-p0_pdf_param1-p0_pdf_param2  ! All p1s already in log space
           ! Likelihood values will remain unchanged since x doesn't change

         if (pdf_param2_pdf .eq. 1) then
             if (pdf_param2_new .lt. pdf_p2_hparam1) pT = -1.e20
             if (pdf_param2_new .gt. pdf_p2_hparam2) pT = -1.e20
         endif
         if (pdf_param1_pdf .eq. 1) then
             if (pdf_param1_new .lt. pdf_p1_hparam1) pT = -1.e20
             if (pdf_param1_new .gt. pdf_p1_hparam2) pT = -1.e20
         endif

         call random_number(randomu)
         if (alog(randomu) .LE. pT) THEN
             !ACCEPT   
             pdf_param1_all(xi)=pdf_param1_new
             pdf_param2_all(xi)=pdf_param2_new    
              !if(it .le. burn_in) then
              if(beta .eq. 1. .and. it .le. burn_in) then  
                  if (xi .LE. nIC) then
                    accept_batch(xi) = accept_batch(xi) + 1
                 else if (xi .GT. nIC) then
                    accept_batch(nIC1) = accept_batch(nIC1) + 1
                 endif  ! xi lt nIC
              endif  ! it lt burn_in
         else
              if(beta .eq. 1. .and. it .le. burn_in) then  
              !if(it .le. burn_in) then                 
                 if (xi .LE. nIC) then 
                     reject_batch(xi) = reject_batch(xi) + 1
                 else
                     reject_batch(nIC1) = reject_batch(nIC1) + 1  
                 endif   ! xi lt nIC
              endif ! it lt burn_in
                 
         endif   ! randomu condition

         if(beta .eq. 1.) then
         if(it .le. burn_in .and. modulo(it,560) .eq. 0) then
             if (it .gt. 100) then
             if (xi .LE. nIC) then
                 if(accept_batch(xi)+reject_batch(xi) .gt. 0) then
                     accep_prob(xi) = real(accept_batch(xi))/(accept_batch(xi)+reject_batch(xi))
                     av_acc = max(0.01,1.0/sqrt(real(it)/500.)) !1.0/sqrt(real(accept(xi)+reject(xi)))
                     if(accep_prob(xi) .lt. 0.2) then
                         stepsize_pdf_p1(xi) = exp(alog(stepsize_pdf_p1(xi)) - av_acc)
                         stepsize_pdf_p2(xi) = exp(alog(stepsize_pdf_p2(xi)) - av_acc)
                     endif
                     if(accep_prob(xi) .gt. 0.6) then
                         stepsize_pdf_p1(xi) = exp(alog(stepsize_pdf_p1(xi)) + av_acc)
                         stepsize_pdf_p2(xi) = exp(alog(stepsize_pdf_p2(xi)) + av_acc)
                     endif
                     accept_batch(xi) = 0
                     reject_batch(xi) = 0
                 endif
             else if (xi .GT. nIC) then
                 if(accept_batch(nIC1)+reject_batch(nIC1) .gt. 0) then    
                      accep_prob(nIC1) = real(accept_batch(nIC1))/(accept_batch(nIC1)+reject_batch(nIC1))
                      av_acc = max(0.01,1.0/sqrt(real(it)/500.)) !1.0/sqrt(real(accept(nIC1)+reject(nIC1)))
                      if(accep_prob(nIC1) .lt. 0.2) then
                         stepsize_pdf_p1(nIC1) = exp(alog(stepsize_pdf_p1(nIC1)) - av_acc)
                         stepsize_pdf_p2(nIC1) = exp(alog(stepsize_pdf_p2(nIC1)) - av_acc)
                      endif
                      if(accep_prob(nIC1) .gt. 0.6) then
                         stepsize_pdf_p1(nIC1) = exp(alog(stepsize_pdf_p1(nIC1)) + av_acc)
                         stepsize_pdf_p2(nIC1) = exp(alog(stepsize_pdf_p2(nIC1)) + av_acc)
                      endif
                      accept_batch(nIC1) = 0
                      reject_batch(nIC1) = 0
                 endif
             endif
             endif

          endif
          endif
  enddo

pdf_param1_out=pdf_param1_all
pdf_param2_out=pdf_param2_all
stepsize_p1_out=stepsize_pdf_p1
stepsize_p2_out=stepsize_pdf_p2
accept_batch_out=accept_batch
reject_batch_out=reject_batch
END SUBROUTINE x_hparam_update

SUBROUTINE x_update(beta,k, x, pdf_param1_all,pdf_param2_all,  &
h_agg,n0,n0T,Rinv, stepsize, &
accept_batch, reject_batch, x_pdf_all, it, burn_in, nIC, kICmax, nmeasure, nIC1, &
x_out, n0_out, n0T_out, accept_out, reject_out, stepsize_out, acc_batch_out, rej_batch_out) 


Implicit none 
INTEGER nmeasure, it, burn_in, k, nIC, kICmax, nIC1
REAL beta, n0T, n0T_out 
REAL av_acc
REAL x(kICmax) 
REAL x_out(kICmax) 
REAL h_agg(nmeasure,kICmax)   
REAL n0(nmeasure) 
REAL n0_out(nmeasure) 
REAL Rinv(nmeasure,nmeasure)
REAL dy(nmeasure)
REAL n1(nmeasure) 
INTEGER x_pdf_all(nIC1)
REAL pdf_param1_all(kICmax), pdf_param2_all(kICmax)
INTEGER elem(5), elem2(5)
REAL u(5),u2(5)
REAL pdf_param1        
REAL pdf_param2
REAL stepsize(nIC1)
REAL accep_prob(nIC1)
INTEGER accept(nIC1), reject(nIC1)
INTEGER accept_batch(nIC1)
INTEGER reject_batch(nIC1)
INTEGER accept_out(nIC1), reject_out(nIC1)
REAL stepsize_out(nIC1)
REAL C(nmeasure)
INTEGER xi, x_pdf,xx
REAL dx, n1T, pT, randomu, p0,p1
REAL x_new(kICmax)
REAL stepsize0
INTEGER acc_batch_out(nIC1)
INTEGER rej_batch_out(nIC1)

accept=0
reject=0

   ! CHANGE OF EMISSIONS VALUES
  call random_number(u)   
  call random_number(u2)  
  elem = FLOOR(k*u)+1+nIC
  elem2 = FLOOR(nIC*u2)+1
  
  !do xx=1,nIC+5
  do xx=1,10
  
  if (xx .LE. 5) then
  !if (xx .LE. nIC) then
     !xi = xx
     xi=elem2(xx)
  else
     !xi=elem(xx-nIC)
     xi=elem(xx-5)
  endif

  if (xi .LE. nIC) then
     x_pdf = x_pdf_all(xi)
     pdf_param1 = pdf_param1_all(xi)
     pdf_param2 = pdf_param2_all(xi)
     stepsize0 = stepsize(xi)

  else if (xi .GT. nIC) then
     x_pdf = x_pdf_all(nIC1)
     pdf_param1 = pdf_param1_all(xi)
     pdf_param2 = pdf_param2_all(xi)
     stepsize0 = stepsize(nIC1)
    
  endif
  dx = random_normal()*stepsize0

  x_new = x
  x_new(xi) =x(xi)+dx
  p0=0.
  p1=0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
         dy=h_agg(:,xi)*dx 
         n1=n0+dy

         C = matmul(n1,Rinv)
         n1T= dot_product(n1,C)
          
         ! hyperparams are fixed below to be a single number - will only apply when x is a scaling of the prior
                  
         call calc_pdf(x(xi),pdf_param1,pdf_param2,x_pdf, p0)           ! Will apply whatever the PDF
         call calc_pdf(x_new(xi),pdf_param1,pdf_param2,x_pdf, p1)        
       
         pT = p1-p0 -0.5*(n1T - n0T)*beta ! All p1s already in log space

         if (x_pdf .eq. 1) then
             if (x_new(xi) .lt. pdf_param1) pT = -1.e20
             if (x_new(xi) .gt. pdf_param2) pT = -1.e20
         endif

         call random_number(randomu)
         if (alog(randomu) .LE. pT) THEN

             !ACCEPT
             x(xi)=x(xi)+dx
             n0=n1
             n0T=n1T
             if (xi .LE. nIC) then
                !if (beta .EQ. 1. .and. it .GT. burn_in) accept(xi) = accept(xi) + 1
                if (it .GT. burn_in) accept(xi) = accept(xi) + 1
             else if (xi .GT. nIC) then
                !if (beta .EQ. 1. .and. it .GT. burn_in) accept(nIC1) = accept(nIC1) + 1
                if (it .GT. burn_in) accept(nIC1) = accept(nIC1) + 1
             endif
             
             !if (beta .EQ. 1. .and. it .GT. burn_in) accept = accept + 1
                                
         else
             !REJECT
             !if (beta .EQ. 1. .and. it .GT. burn_in) then
             if (it .GT. burn_in) then  
                 if (xi .LE. nIC) then 
                     reject(xi) = reject(xi) + 1
                 else
                     reject(nIC1) = reject(nIC1) + 1  
                 endif
             endif                    
         endif   ! randomu condition

         ! Stepsize tuning
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if(beta .eq. 1. .and. it .le. burn_in) then
             if (alog(randomu) .LE. pT) THEN

                 if (xi .LE. nIC) then
                     accept_batch(xi) = accept_batch(xi) + 1
                 else if (xi .GT. nIC) then
                    accept_batch(nIC1) = accept_batch(nIC1) + 1
                 endif
             else
                 if (xi .LE. nIC) then 
                     reject_batch(xi) = reject_batch(xi) + 1
                 else
                     reject_batch(nIC1) = reject_batch(nIC1) + 1  
                 endif

             endif

         endif

         if(beta .eq. 1. .and. it .le. burn_in .and. modulo(it,560) .eq. 0) then
             if (it .gt. 100) then
             if (xi .LE. nIC) then
                 if(accept_batch(xi)+reject_batch(xi) .gt. 0) then
                     accep_prob(xi) = real(accept_batch(xi))/(accept_batch(xi)+reject_batch(xi))
                     av_acc = max(0.01,1.0/sqrt(real(it/500))) !1.0/sqrt(real(accept(xi)+reject(xi)))
                     if(accep_prob(xi) .lt. 0.2) stepsize(xi) = exp(alog(stepsize(xi)) - av_acc)
                     if(accep_prob(xi) .gt. 0.6) stepsize(xi) = exp(alog(stepsize(xi)) + av_acc)
                     accept_batch(xi) = 0
                     reject_batch(xi) = 0
                 endif
             else if (xi .GT. nIC) then
                 if(accept_batch(nIC1)+reject_batch(nIC1) .gt. 0) then    
                      accep_prob(nIC1) = real(accept_batch(nIC1))/(accept_batch(nIC1)+reject_batch(nIC1))
                      av_acc = max(0.01,1.0/sqrt(real(it/500))) !1.0/sqrt(real(accept(nIC1)+reject(nIC1)))
                      if(accep_prob(nIC1) .lt. 0.2) stepsize(nIC1) = exp(alog(stepsize(nIC1)) - av_acc)
                      if(accep_prob(nIC1) .gt. 0.6) stepsize(nIC1) = exp(alog(stepsize(nIC1)) + av_acc)
                      accept_batch(nIC1) = 0
                      reject_batch(nIC1) = 0
                 endif
             endif
             endif

         endif
         ! End of stepsize tuning
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

x_out=x
n0_out=n0
n0T_out=n0T
accept_out=accept
reject_out=reject
stepsize_out=stepsize
acc_batch_out = accept_batch
rej_batch_out = reject_batch
END SUBROUTINE x_update



SUBROUTINE birth(beta,k, x, h_agg,y,n0,n0T,Rinv, plon, plat, regions_v, lon,lat, & 
h_v,pdf_param1, pdf_param2, x_pdf, sigma_bd, &
it, burn_in, nIC, kICmax, kmax, nmeasure, Ngrid,nlon,nlat, &
k_out, x_out, h_agg_out, n0_out, n0T_out, regions_v_out, plon_out, plat_out, accept_out, reject_out, &
pdf_param1_out, pdf_param2_out)

IMPLICIT NONE
! Dimensions
INTEGER nmeasure, Ngrid, nlon,nlat, k,kmax, nIC, kICmax
!REAL lonmin,lonmax,latmin,latmax
REAL sigma_bd
! Single Variables
INTEGER accept_birth, reject_birth, it, burn_in, x_pdf 
REAL beta, n0T
! Input arrays
REAL x(kICmax) 
REAL pdf_param1(kICmax) 
REAL pdf_param2(kICmax)
REAL h_agg(nmeasure,kICmax)  
REAL y(nmeasure) 
REAL n0(nmeasure) 
REAL Rinv(nmeasure,nmeasure)
REAL plon(kmax)
REAL plat(kmax)
INTEGER regions_v(Ngrid)
REAL lon(nlon)
REAL lat(nlat)
REAL h_v(nmeasure, Ngrid)
! Outputs
INTEGER k_out
REAL x_out(kICmax) 
REAL pdf_param1_out(kICmax), pdf_param2_out(kICmax)
REAL h_agg_out(nmeasure,kICmax)
REAL n0_out(nmeasure) 
REAL n0T_out
INTEGER regions_v_out(Ngrid)
REAL plon_out(kmax)
REAL plat_out(kmax)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER ri, rib, k1, errstat,jj, kIC      
REAL u, u2, plon_new, plat_new, c_new, x_new, n1Tb, pT_birth, randomu
REAL mu, sigma, pdf_param1_new, pdf_param2_new
INTEGER regions_v1b(Ngrid)       
REAL n1b(nmeasure), C(nmeasure)     
INTEGER ilon,ilat, reject_stat,zi       
! Allocatable arrays
REAL, DIMENSION(:),ALLOCATABLE :: plon1b, plat1b, x1b
REAL, DIMENSION(:,:), ALLOCATABLE :: h_agg2
REAL,PARAMETER     :: pi = 3.14159265 

accept_birth=0
reject_birth=0

k1=k+1
kIC=k1+nIC

if (k1 .LT. kmax) THEN            
  
  allocate(plon1b(k1))     
  allocate(plat1b(k1))     
  allocate(h_agg2(nmeasure, kIC))
  allocate(x1b(kIC))               

  ! 1. Select new cell location - needs to be different to current locations
  !call random_number(u)
  !plon_new = lonmin+(lonmax-lonmin)*u
  !call random_number(u2)
  !plat_new = latmin+(latmax-latmin)*u2
  
  ! 1. Select new nucleus location - different from current locations and on hte underlying grid
  call random_number(u) 
  call random_number(u2)
  ilon=FLOOR(nlon*u)+1   ! Choose random number between 1 and nlon
  ilat=FLOOR(nlat*u2)+1   ! Choose random number between 1 and nlat
        
  plon_new = lon(ilon)
  plat_new = lat(ilat)
  reject_stat=0
  do zi=1,k
     if (plon(zi) .EQ. plon_new .and. plat(zi) .eq. plat_new) then
        reject_stat=1  
     endif
  enddo

  if (reject_stat .ne. 1) then

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

  pdf_param1_new = pdf_param1(rib+nIC)
  pdf_param2_new = pdf_param2(rib+nIC)
  
  if (x_pdf .EQ. 1) THEN  ! 1=UNIFORM  
                
      if (x_new .GT. pdf_param1_new .and. x_new .LT. pdf_param2_new) THEN
                       
          x1b(1:k+nIC)=x(1:k+nIC)
          x1b(kIC) = x_new

          n1b=matmul(h_agg2,x1b)-y 
          C = matmul(n1b,Rinv)
          n1Tb= dot_product(n1b,C)                   
          !n1Tb=sum((n1b/nmeasure)**2)

          pT_birth = alog(sqrt(2*pi)*sigma_bd/(pdf_param2_new-pdf_param1_new)) &
              + ((x_new-c_new)**2)/2./sigma_bd**2 - 0.5*(n1Tb - n0T)*beta
       
         call random_number(randomu)             
         if (alog(randomu) .LE. pT_birth) THEN
           !ACCEPT
           k=k1
           x(:)=0.
           x(1:kIC)=x1b(1:kIC)
           pdf_param1(kIC)=pdf_param1_new
           pdf_param2(kIC)=pdf_param2_new
           h_agg(:,:)=0.
           h_agg(:,1:kIC)=h_agg2(:,1:kIC)
           n0=n1b
           n0T=n1Tb
           regions_v(:)=regions_v1b
           plon(:)=0.
           plat(:)=0.
           plon(1:k1)=plon1b(1:k1)
           plat(1:k1)=plat1b(1:k1)
           !if (beta .EQ. 1. .and. it .GT. burn_in) accept_birth=accept_birth+1 
           if (it .GT. burn_in) accept_birth=accept_birth+1       
         else 
           !REJECT
           !if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1  
           if (it .GT. burn_in) reject_birth=reject_birth+1  
         endif

      else
          !REJECT if x_new is negative
          !if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1
          if (it .GT. burn_in) reject_birth=reject_birth+1
      endif   

  else if (x_pdf .GE. 2) THEN     ! 2=GAUSSIAN 3=LOGNORMAL 

          x1b(1:k+nIC)=x(1:k+nIC)
          x1b(kIC) = x_new
          
          n1b=matmul(h_agg2,x1b)-y  
          C = matmul(n1b,Rinv)
          n1Tb= dot_product(n1b,C)                                       
          !n1Tb=sum((n1b/sigma_y)**2)

          if (x_pdf .EQ. 2) THEN   !2=GAUSSIAN

              pT_birth = alog(sigma_bd/pdf_param2_new) + ((x_new-c_new)**2)/2./sigma_bd**2 - &
                         ((x_new-pdf_param1_new)**2)/2./pdf_param2_new**2 &
                         - 0.5*(n1Tb - n0T)*beta  

          else if (x_pdf .EQ. 3) THEN
 
              mu = alog(pdf_param1_new) - 0.5*alog(1. + pdf_param2_new**2/pdf_param1_new**2)
              sigma=sqrt(alog((pdf_param2_new/pdf_param1_new)**2 + 1.))

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
           pdf_param1(kIC)=pdf_param1_new
           pdf_param2(kIC)=pdf_param2_new
           h_agg(:,:)=0.
           h_agg(:,1:kIC)=h_agg2(:,1:kIC)
           n0=n1b
           n0T=n1Tb
           regions_v(:)=regions_v1b
           plon(:)=0.
           plat(:)=0.
           plon(1:k1)=plon1b(1:k1)
           plat(1:k1)=plat1b(1:k1)
          !if (beta .EQ. 1. .and. it .GT. burn_in) accept_birth=accept_birth+1    
           if (it .GT. burn_in) accept_birth=accept_birth+1    
         else 
           !REJECT
           !if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1  
           if (it .GT. burn_in) reject_birth=reject_birth+1  
         endif

  endif       ! x_pdf

 else
    !REJECT if plon_new and plat_new are the same as any other point
    !if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1
    if (it .GT. burn_in) reject_birth=reject_birth+1  
 endif
        
else
    !REJECT if k1 > kmax
    !if (beta .EQ. 1. .and. it .GT. burn_in) reject_birth=reject_birth+1
    if (it .GT. burn_in) reject_birth=reject_birth+1  
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
pdf_param1_out=pdf_param1
pdf_param2_out=pdf_param2

END SUBROUTINE birth


SUBROUTINE death(beta,k, x, h_agg,y,n0,n0T,Rinv, plon, plat, regions_v, lon,lat, & 
h_v, pdf_param1, pdf_param2, x_pdf, sigma_bd, &
it, burn_in, nIC, kICmax, kmin, kmax, nmeasure, Ngrid,nlon,nlat, &
k_out, x_out, h_agg_out, n0_out, n0T_out, regions_v_out, &
plon_out, plat_out, accept_out, reject_out, pdf_param1_out, pdf_param2_out)



IMPLICIT NONE
! Dimensions
INTEGER kmax,nmeasure,Ngrid,nlon,nlat, accept_death, reject_death, it, burn_in, kmin, k, nIC, kICmax
REAL sigma_bd, beta, n0T 
! Input arrays
REAL x(kICmax) 
REAL pdf_param1(kICmax)
REAL pdf_param2(kICmax)
REAL h_agg(nmeasure,kICmax) 
REAL y(nmeasure) 
REAL n0(nmeasure) 
REAL Rinv(nmeasure,nmeasure)
REAL plon(kmax)
REAL plat(kmax)
INTEGER regions_v(Ngrid)
REAL lon(nlon)
REAL lat(nlat)
REAL h_v(nmeasure, Ngrid)
INTEGER x_pdf
! Outputs
INTEGER k_out 
REAL n0T_out
REAL x_out(kICmax)
REAL pdf_param1_out(kICmax), pdf_param2_out(kICmax)
REAL h_agg_out(nmeasure,kICmax)
REAL n0_out(nmeasure) 
INTEGER regions_v_out(Ngrid)
REAL plon_out(kmax)
REAL plat_out(kmax)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER ri, rid, k1d, jj, ci_rm, errstat, kIC     
REAL u, plon_rm, plat_rm, x_cell, x_rm, n1Td, pT_death, randomu
REAL mu,sigma, pdf_param1_rm, pdf_param2_rm
INTEGER regions_v1d(Ngrid)      
REAL n1d(nmeasure), C(nmeasure)             
! Allocatable arrays
REAL, DIMENSION(:),ALLOCATABLE :: plon1d, plat1d, x1d, pdf_param1d, pdf_param2d 
REAL, DIMENSION(:,:), ALLOCATABLE :: h_agg2d

REAL,PARAMETER     :: pi = 3.14159265 

accept_death=0
reject_death=0

!DEATH
k1d=k-1

kIC=k1d+nIC

if (k1d .GE. kmin) THEN            

  allocate(plon1d(k1d))     
  allocate(plat1d(k1d))     
  allocate(h_agg2d(nmeasure, kIC))
  allocate(x1d(kIC))  
  allocate(pdf_param1d(kIC)) 
  allocate(pdf_param2d(kIC))  
           
  ! 1. Select new cell location - needs to be different to current locations
   
  call random_number(u)   
  ci_rm = FLOOR((k1d+1)*u) + 1
  
  plon_rm = plon(ci_rm)
  plat_rm = plat(ci_rm)
  x_rm = x(ci_rm+nIC)
  pdf_param1_rm = pdf_param1(ci_rm+nIC)  
  pdf_param2_rm = pdf_param2(ci_rm+nIC) 
  

  IF (nIC .GT. 0) THEN
     x1d(1:nIC) = x(1:nIC)
     pdf_param1d(1:nIC) = pdf_param1(1:nIC)
     pdf_param2d(1:nIC) = pdf_param2(1:nIC)
  ENDIF

  IF (ci_rm .EQ. 1) THEN
      plon1d(1:k1d) = plon(2:(k1d+1))
      plat1d(1:k1d) = plat(2:(k1d+1))
      x1d(nIC+1:kIC) = x(nIC+2:(kIC+1))
      pdf_param1d(nIC+1:kIC) = pdf_param1(nIC+2:(kIC+1))  
      pdf_param2d(nIC+1:kIC) = pdf_param2(nIC+2:(kIC+1))  
  
  ELSEIF (ci_rm .EQ. (k1d+1)) THEN
      plon1d(1:k1d) = plon(1:k1d)
      plat1d(1:k1d) = plat(1:k1d)
      x1d(nIC+1:kIC) = x(nIC+1:kIC)
      pdf_param1d(nIC+1:kIC) = pdf_param1(nIC+1:kIC)  
      pdf_param2d(nIC+1:kIC) = pdf_param2(nIC+1:kIC)  

  ELSE   
      plon1d(1:(ci_rm-1)) = plon(1:(ci_rm-1))                   
      plon1d(ci_rm:k1d) = plon((ci_rm+1):(k1d+1))
      plat1d(1:(ci_rm-1)) = plat(1:(ci_rm-1))                   
      plat1d(ci_rm:k1d) = plat((ci_rm+1):(k1d+1)) 

      x1d(nIC+1:(ci_rm+nIC-1)) = x(nIC+1:(ci_rm+nIC-1)) 
      x1d(ci_rm+nIC:kIC) = x((ci_rm+nIC+1):(kIC+1))
      pdf_param1d(nIC+1:(ci_rm+nIC-1)) = pdf_param1(nIC+1:(ci_rm+nIC-1)) 
      pdf_param1d(ci_rm+nIC:kIC) = pdf_param1((ci_rm+nIC+1):(kIC+1))
      pdf_param2d(nIC+1:(ci_rm+nIC-1)) = pdf_param2(nIC+1:(ci_rm+nIC-1)) 
      pdf_param2d(ci_rm+nIC:kIC) = pdf_param2((ci_rm+nIC+1):(kIC+1))

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
  C = matmul(n1d,Rinv)
  n1Td= dot_product(n1d,C)
  !n1Td=sum((n1d/sigma_y)**2)                           
  ! ACCEPTANCE PROBABILITY  

  IF (x_pdf .EQ. 1) THEN  ! 1 = UNIFORM
     pT_death = alog((pdf_param2_rm-pdf_param1_rm)/sqrt(2.*pi)/sigma_bd) &
      - ((x_cell-x_rm)**2)/2./sigma_bd**2 -0.5*(n1Td - n0T)*beta

  ELSE IF (x_pdf .EQ. 2) THEN   ! 2=GAUSSIAN

      pT_death = alog(pdf_param2_rm/sigma_bd) - ((x_cell-x_rm)**2)/2./sigma_bd**2 + &
                         ((x_rm-pdf_param1_rm)**2)/2./pdf_param2_rm**2 &
                          - 0.5*(n1Td - n0T)*beta


  ELSE IF (x_pdf .EQ. 3) THEN   ! 3=LOGNORMAL

      mu = alog(pdf_param1_rm) - 0.5*alog(1. + pdf_param2_rm**2/pdf_param1_rm**2)
      sigma=sqrt(alog((pdf_param2_rm/pdf_param1_rm)**2 + 1.))
     
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
           pdf_param1(:)=0.
           pdf_param1(1:kIC)=pdf_param1d
           pdf_param2(:)=0.
           pdf_param2(1:kIC)=pdf_param2d
           h_agg(:,:)=0.
           h_agg(:,1:kIC)=h_agg2d
           n0=n1d
           n0T=n1Td
           regions_v(:)=regions_v1d
           plon(:)=0.
           plat(:)=0.
           plon(1:k1d)=plon1d
           plat(1:k1d)=plat1d
           !if (beta .EQ. 1. .and. it .GT. burn_in) accept_death=accept_death+1   
           if (it .GT. burn_in) accept_death=accept_death+1   

       else 
           !REJECT
           !if (beta .EQ. 1. .and. it .GT. burn_in) reject_death=reject_death+1
           if (it .GT. burn_in) reject_death=reject_death+1    
       endif
        
else
    !REJECT if k1d < kmin
    !if (beta .EQ. 1. .and. it .GT. burn_in) reject_death=reject_death+1 
    if (it .GT. burn_in) reject_death=reject_death+1 
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
if(allocated(pdf_param1d))  deallocate(pdf_param1d,stat=errstat)
  if (errstat /= 0) stop
if(allocated(pdf_param2d))  deallocate(pdf_param2d,stat=errstat)
  if (errstat /= 0) stop


k_out=k
x_out=x
pdf_param1_out=pdf_param1
pdf_param2_out=pdf_param2
h_agg_out=h_agg
n0_out=n0
n0T_out=n0T
regions_v_out=regions_v
plon_out=plon
plat_out=plat
accept_out=accept_death
reject_out=reject_death


END SUBROUTINE death


SUBROUTINE move(beta,k, x, h_agg, y,n0,n0T,Rinv, plon, plat, regions_v, lon,lat, & 
h_v, lonmin, lonmax, latmin,latmax, sigma_clon, sigma_clat, it, &
burn_in, nIC, kICmax, kIC, kmax, nmeasure, Ngrid,nlon,nlat, &
h_agg_out, n0_out, n0T_out, regions_v_out, plon_out, plat_out, accept_out, reject_out)



IMPLICIT NONE
! Dimensions
INTEGER kmax, nmeasure, Ngrid, nlon,nlat, accept_move, reject_move, it, burn_in,k, nIC, kIC, kICmax
! Single Variables
REAL lonmin,lonmax,latmin,latmax, sigma_clon, sigma_clat, beta 
! Input arrays
REAL x(kICmax) 
REAL h_agg(nmeasure,kICmax)  
REAL y(nmeasure) 
REAL n0(nmeasure) 
REAL n0T
REAL Rinv(nmeasure,nmeasure)
REAL plon(kmax)
REAL plat(kmax)
INTEGER regions_v(Ngrid)
REAL lon(nlon)
REAL lat(nlat)
REAL h_v(nmeasure, Ngrid)
REAL plon1m(k)
REAL plat1m(k)
REAL x1m(kIC)
REAL h_agg2m(nmeasure,kIC)
! Outputs
REAL h_agg_out(nmeasure,kICmax)
REAL n0_out(nmeasure) 
REAL n0T_out
INTEGER regions_v_out(Ngrid)
REAL plon_out(kmax)
REAL plat_out(kmax)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER  ri, k1, jj, ci_mv    
REAL u, n1Tm, pT_move, randomu
INTEGER regions_v1m(Ngrid)      
REAL n1m(nmeasure), C(nmeasure) 
INTEGER reject_stat, zi
! Allocatable arrays
! None
REAL,PARAMETER     :: pi = 3.14159265 

accept_move=0
reject_move=0

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
    
   !plon1m(ci_mv) = random_normal()*sigma_clon+plon1m(ci_mv)
   !plat1m(ci_mv) = random_normal()*sigma_clat+plat1m(ci_mv)

   ! Move onto centre of underlying grid 
   plon1m(ci_mv) = FLOOR(random_normal()*sigma_clon)*(lon(2)-lon(1))+plon1m(ci_mv)
   plat1m(ci_mv) = FLOOR(random_normal()*sigma_clat)*(lat(2)-lat(1))+plat1m(ci_mv)
            
   reject_stat=0
   do zi=1,k
     if (plon(zi) .EQ. plon1m(ci_mv) .and. plat(zi) .eq. plat1m(ci_mv)) then
        reject_stat=1  
     endif
  enddo    

   if (reject_stat .ne. 1) then

   ! Need to reject if outside of lon/lat range.
   IF (plon1m(ci_mv) .GT. lonmax .OR. plon1m(ci_mv) .LT. lonmin) THEN
       !if (beta .EQ. 1. .and. it .GT. burn_in) reject_move=reject_move+1
       if (it .GT. burn_in) reject_move=reject_move+1
                         
   ELSEIF (plat1m(ci_mv) .GT. latmax .OR. plat1m(ci_mv) .LT. latmin) THEN           
       !if (beta .EQ. 1. .and. it .GT. burn_in) reject_move=reject_move+1 
       if (it .GT. burn_in) reject_move=reject_move+1 

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
     C = matmul(n1m,Rinv)
     n1Tm= dot_product(n1m,C)
     !n1Tm=sum((n1m/sigma_y)**2)
                                     
     ! ACCEPTANCE PROBABILITY   
     pT_move = (n1Tm - n0T)*(-0.5)*beta  
                              
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
           !if (beta .EQ. 1. .and. it .GT. burn_in) accept_move=accept_move+1   
           if (it .GT. burn_in) accept_move=accept_move+1  

      else 
           !REJECT
           !if (beta .EQ. 1. .and. it .GT. burn_in) reject_move=reject_move+1  
           if (it .GT. burn_in) reject_move=reject_move+1  
      endif
        
   
  ENDIF    ! plon & lat within range

  else   ! lon_new, lat_new on same location as another point
      !REJECT
       !if (beta .EQ. 1. .and. it .GT. burn_in) reject_move=reject_move+1  
       if (it .GT. burn_in) reject_move=reject_move+1  
endif   

h_agg_out=h_agg
n0_out=n0
n0T_out=n0T
regions_v_out=regions_v
plon_out=plon
plat_out=plat
accept_out=accept_move
reject_out=reject_move

END SUBROUTINE move

SUBROUTINE sigma_y_update(beta, sigma_model_current, sigma_measure, sigma_y_current, error_structure, &
detval_current, detval_Q, Rinv_current, Qinv, sigma_model_hparam1, sigma_model_hparam2, &
stepsize_sigma_y, sigma_model_pdf, R_indices,n0,n0T, accept_batch, reject_batch, &
it, burn_in, nmeasure, dim1, dim2, &
n0T_out, accept_out, reject_out, sigma_y_out, sigma_model_out, Rinv_out, detval_out, &
stepsize_sig_out, accept_batch_out, reject_batch_out) 


IMPLICIT NONE

! Dimensions
INTEGER nmeasure, accept, reject, it, burn_in, dim1, dim2
! Single Variables
REAL beta 
! Input arrays
REAL n0(nmeasure) 
REAL n0T
REAL detval_current, detval_Q
REAL sigma_measure(nmeasure)
REAL sigma_model_current(dim2)
REAl sigma_y_current(nmeasure)
REAL error_structure(nmeasure)
INTEGER R_indices(dim1,dim2)
REAL Rinv_current(nmeasure,nmeasure), Qinv(nmeasure,nmeasure)
REAL stepsize_sigma_y(dim2)
REAL sigma_model_hparam1(dim2)
REAL sigma_model_hparam2(dim2)
INTEGER sigma_model_pdf
! Outputs
REAL n0T_out, detval_out
REAL sigma_y_out(nmeasure)
REAL sigma_model_out(dim2)
REAL Rinv_out(nmeasure,nmeasure)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER  yi, jj, ti,ii
REAL randomu, dsigma_y, sigma_model_new, u, av_acc
REAL p0_sigma_y, p1_sigma_y, n1T, detval_new, pT
REAL y_error_new(nmeasure), autocorr_vec(nmeasure)
REAL sigma_y_new(nmeasure), sigma_yinv_new(nmeasure)
REAL Rinv_new(nmeasure,nmeasure), C(nmeasure)
REAL stepsize_sig_out(dim2)
INTEGER accept_batch(dim2), reject_batch(dim2)
INTEGER accept_batch_out(dim2), reject_batch_out(dim2)
REAL accep_prob

accept=0
reject=0

 !do yi=1,dim2
    call random_number(u)   
    yi = FLOOR(dim2*u)+1
		
    ! Generate new value of sigma_y
    dsigma_y = random_normal()*stepsize_sigma_y(yi)
    sigma_model_new = sigma_model_current(yi) + dsigma_y       

    call calc_pdf(sigma_model_current(yi), sigma_model_hparam1(yi), sigma_model_hparam2(yi), sigma_model_pdf, p0_sigma_y)
    call calc_pdf(sigma_model_new, sigma_model_hparam1(yi), sigma_model_hparam2(yi), sigma_model_pdf, p1_sigma_y)
	
    do jj=1,dim2   
        !y_error_new(R_indices(:,jj)) = sigma_model_current(jj)    ! Provided dim2 isn't too big then should be fine
        y_error_new(R_indices(:,jj)) = sigma_model_current(jj)*error_structure(R_indices(:,jj)) 
    enddo  

    ! Change one element of sigma_y
    !y_error_new(R_indices(:,yi)) = sigma_model_new  ! R_indices = array of indices to which each sigma_y applies 
    y_error_new(R_indices(:,yi)) = sigma_model_new*error_structure(R_indices(:,yi))
    sigma_y_new=sqrt(y_error_new**2 + sigma_measure**2)   
    sigma_yinv_new = 1./sigma_y_new

     do ti=1,nmeasure
        autocorr_vec = sigma_yinv_new(ti)*sigma_yinv_new*Qinv(:,ti)
        Rinv_new(:,ti) = autocorr_vec      
    enddo  


    ! caluclate new determinant by scaling
    detval_new = sum(alog(sigma_y_new)) + detval_Q

 
    C = matmul(n0,Rinv_new)
    n1T= dot_product(n0,C)
    

    ! Compute P1/P0 	
    pT= p1_sigma_y-p0_sigma_y - detval_new*beta + detval_current*beta - 0.5*(n1T - n0T)*beta       !*beta      ! -detval_new becasue it's inverse of true determinant
    if (sigma_model_pdf .eq. 1) then
       if (sigma_model_new .lt. sigma_model_hparam1(yi)) pT = -1.e20
       if (sigma_model_new .gt. sigma_model_hparam2(yi)) pT = -1.e20
    endif
  
    call random_number(randomu)     ! Generates uniformly distributed random number
    
    if(alog(randomu) .le. pT) then      
       !ACCEPT	
       sigma_model_current(yi) = sigma_model_new
       sigma_y_current = sigma_y_new
       detval_current = detval_new
       Rinv_current = Rinv_new
       n0T=n1T
       if(beta .eq. 1. .and. it .gt. burn_in) accept=accept + 1
       if(beta .eq. 1. .and. it .le. burn_in) accept_batch(yi)=accept_batch(yi) + 1
    else
       !;REJECT					
       if(beta .eq. 1. .and. it .gt. burn_in) reject=reject + 1
       if(beta .eq. 1. .and. it .le. burn_in) reject_batch(yi)=reject_batch(yi) + 1
    endif
  
     if(beta .eq. 1. .and. it .le. burn_in .and. modulo(it,560) .eq. 1) then
       if (it .gt. 100) then
        av_acc =  max(0.01,1.0/sqrt(real(it)/500.))
        do ii=1,dim2 
          if(accept_batch(ii)+reject_batch(ii) .gt. 0) then          
             accep_prob = real(accept_batch(ii))/(accept_batch(ii)+reject_batch(ii))
             if(accep_prob .lt. 0.2) stepsize_sigma_y(ii) = exp(alog(stepsize_sigma_y(ii)) - av_acc)
             if(accep_prob .gt. 0.6) stepsize_sigma_y(ii) = exp(alog(stepsize_sigma_y(ii)) + av_acc)
             accept_batch(ii) = 0
             reject_batch(ii) = 0
          endif
        enddo
       endif   ! it .gt. 1
    endif

  !enddo   ! yi loop

n0T_out=n0T
sigma_y_out=sigma_y_current
sigma_model_out=sigma_model_current
detval_out=detval_current
Rinv_out=Rinv_current
accept_out=accept
reject_out=reject
stepsize_sig_out=stepsize_sigma_y
accept_batch_out=accept_batch
reject_batch_out=reject_batch
END SUBROUTINE sigma_y_update


SUBROUTINE tau_update(beta, tau_current, sigma_y,  &
detval_current, detval_Q_current, detval_Q_block, tau_hparam1, tau_hparam2,  &
stepsize_tau, tau_pdf, Rinv_current, Qinv_current, deltatime, nmeasure_site, &
n0, n0T, accept_batch, reject_batch, it, burn_in, nmeasure, numsites,  &
n0T_out, accept_out, reject_out, tau_out, Rinv_out, Qinv_out, detval_out, &
detval_Q_out, detval_Q_block_out, &
stepsize_tau_out, accept_batch_out,reject_batch_out) 

IMPLICIT NONE

! Dimensions
INTEGER nmeasure, accept, reject, it, burn_in, numsites
! Single Variables
REAL beta 
REAL tau_current(numsites)
REAL stepsize_tau
REAL tau_hparam1, tau_hparam2
INTEGER tau_pdf
REAL n0T
REAL detval_current, detval_Q_current 
REAL detval_Q_block(numsites)
! Input arrays
REAL n0(nmeasure), sigma_y(nmeasure) 
REAL Rinv_current(nmeasure,nmeasure), Qinv_current(nmeasure,nmeasure)
REAL deltatime
INTEGER nmeasure_site(numsites)
! Outputs
REAL n0T_out, detval_out, tau_out(numsites), detval_Q_out
INTEGER accept_out, reject_out
REAL Rinv_out(nmeasure,nmeasure) 
REAL Qinv_out(nmeasure, nmeasure)
REAL detval_Q_block_out(numsites)
! Intermediate variables
INTEGER ti, yi, ii
REAL randomu, dtau, tau_new(numsites), u
REAL p0_tau, p1_tau, n1T, detval_new, pT
REAL detval_Q_new
REAL Rinv_new(nmeasure, nmeasure)
REAL Qinv_new(nmeasure, nmeasure)
REAL C(nmeasure), autocorr_vec(nmeasure)
REAL sigma_yinv(nmeasure)
INTEGER cum_nmeas
!REAL Q_block_new(nsite_max,nsite_max)
!REAL Q_block_inv_new(nsite_max,nsite_max)
REAL detval_Q_block_new(numsites)
REAL q_small
REAL stepsize_tau_out
INTEGER accept_batch, reject_batch
INTEGER accept_batch_out, reject_batch_out
REAL accep_prob, av_acc

Qinv_new=Qinv_current
detval_Q_block_new=detval_Q_block
tau_new=tau_current

accept=0
reject=0

 call random_number(u)   
 yi = FLOOR(numsites*u)+1
! Generate new value of tau
dtau = random_normal()*stepsize_tau
tau_new(yi) = tau_current(yi) + dtau 

!if (tau_new(yi) .GT. tau_hparam1 .and. tau_new(yi) .LT. tau_hparam2) THEN

! Compute P1 for new value of tau
   call calc_pdf(tau_current(yi), tau_hparam1, tau_hparam2, tau_pdf, p0_tau)
   call calc_pdf(tau_new(yi), tau_hparam1, tau_hparam2, tau_pdf, p1_tau)
    		
    sigma_yinv = 1./sigma_y

    if (yi .EQ. 1) THEN
        cum_nmeas = 0
    elseif (yi .EQ. 2) THEN
        cum_nmeas = nmeasure_site(1)
    else
        cum_nmeas = sum(nmeasure_site(1:yi-1))
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       q_small = exp(-1.*deltatime/tau_new(yi))    ! deltatime = single number; q_small=single number

       Qinv_new(cum_nmeas+1,cum_nmeas+1) = 1.
       Qinv_new(cum_nmeas+nmeasure_site(yi),cum_nmeas+nmeasure_site(yi)) = 1.
       Qinv_new(cum_nmeas+1,cum_nmeas+2) = q_small*(-1.)
       Qinv_new(cum_nmeas+nmeasure_site(yi),cum_nmeas+nmeasure_site(yi)-1) = q_small*(-1.)

       do ii=2, nmeasure_site(yi)-1
           Qinv_new(cum_nmeas+ii,cum_nmeas+ii) = 1. + q_small**2
           Qinv_new(cum_nmeas+ii,cum_nmeas+ii+1) = q_small*(-1.)
           Qinv_new(cum_nmeas+ii,cum_nmeas+ii-1) = q_small*(-1.)
       enddo
   
       detval_Q_block_new(yi) = alog(1-q_small**2)*((nmeasure_site(yi)-1)/2.)   ! DETERMINANT IN LOG SPACE log of sqrt of determinant
   
       !cum_nmeas = cum_nmeas + nmeasure_site(yi)
     
       Qinv_new(cum_nmeas+1:cum_nmeas+nmeasure_site(yi), cum_nmeas+1:cum_nmeas+nmeasure_site(yi)) = &
       Qinv_new(cum_nmeas+1:cum_nmeas+nmeasure_site(yi), cum_nmeas+1:cum_nmeas+nmeasure_site(yi))/(1.-q_small**2)
       detval_Q_new = sum(detval_Q_block_new)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    do ti=1,nmeasure
        autocorr_vec = sigma_yinv(ti)*sigma_yinv*Qinv_new(:,ti)
        Rinv_new(:,ti) = autocorr_vec      
    enddo

    detval_new =  sum(alog(sigma_y)) + detval_Q_new  ! DETERMINANT IN LOG SPACE log of sqrt of determinant

    C = matmul(n0,Rinv_new)
    n1T = dot_product(n0,C)
		
    ! compute P1/P0 
    pT=p1_tau - p0_tau - detval_new*beta + detval_current*beta - 0.5*(n1T - n0T)*beta    ! This way round or should detvals be opposite? 
     if (tau_pdf .eq. 1) then
       if (tau_new(yi) .lt. tau_hparam1) pT = -1.e20
       if (tau_new(yi) .gt. tau_hparam2) pT = -1.e20
    endif

    call random_number(randomu)        ! Generates uniformly distributed random number
 
    if(alog(randomu) .le. pT) then      
       !;ACCEPT
       tau_current = tau_new
       Qinv_current = Qinv_new
       Rinv_current= Rinv_new
       detval_current = detval_new
       detval_Q_current = detval_Q_new
       detval_Q_block = detval_Q_block_new
       n0T=n1T
       if(beta .eq. 1. .and. it .gt. burn_in) accept=accept + 1
       if(beta .eq. 1. .and. it .le. burn_in) accept_batch=accept_batch + 1
    else 
       !;REJECT	
        if(beta .eq. 1. .and. it .gt. burn_in) reject=reject + 1
        if(beta .eq. 1. .and. it .le. burn_in) reject_batch=reject_batch + 1
    endif


   if(beta .eq. 1 .and. it .le. burn_in .and. modulo(it,560) .eq. 2) then
       if (it .gt. 100) then
        av_acc =  max(0.01,1.0/sqrt(real(it)/500.))       
          if(accept_batch+reject_batch .gt. 0) then          
             accep_prob = real(accept_batch)/(accept_batch+reject_batch)
             if(accep_prob .lt. 0.2) stepsize_tau = exp(alog(stepsize_tau) - av_acc)
             if(accep_prob .gt. 0.6) stepsize_tau = exp(alog(stepsize_tau) + av_acc)
             accept_batch = 0
             reject_batch = 0
          endif
       endif   ! it .gt. 2
    endif   ! beta=1 etc.

tau_out=tau_current
Qinv_out=Qinv_current
Rinv_out=Rinv_current
detval_out = detval_current
detval_Q_out = detval_Q_current
detval_Q_block_out=detval_Q_block
n0T_out=n0T
accept_out=accept
reject_out=reject
stepsize_tau_out=stepsize_tau
accept_batch_out=accept_batch
reject_batch_out=reject_batch
END SUBROUTINE tau_update


SUBROUTINE y_update(beta, y_current, n0, n0T, &
Rinv, y_hparam1, y_hparam2, timeindex_zero, &
stepsize_y, y_pdf, accept, reject, it, burn_in, nmeasure, nzero, &
n0_out, n0T_out, y_out, accept_out, reject_out) 

IMPLICIT NONE

! Dimensions
INTEGER nmeasure, accept, reject, it, burn_in, nzero
! Single Variables
REAL beta 
! Input arrays
INTEGER timeindex_zero(nzero)
REAL y_current(nmeasure)
REAL n0(nmeasure) 
REAL n0T
REAL Rinv(nmeasure,nmeasure)
REAL stepsize_y
REAL y_hparam1
REAL y_hparam2
INTEGER y_pdf
! Outputs
REAL n0T_out, n0_out(nmeasure)
REAL y_out(nmeasure)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER  yi, xi
REAL randomu, dy, y_new
REAL p0_y, p1_y, n1T, pT
REAL n1(nmeasure) 
REAL dum, dumb,dumc,dumd
REAL dumx(nmeasure), dumy(nmeasure)

 do yi=1,nzero
		
    ! Generate new value of sigma_y
    n1=n0
    xi = timeindex_zero(yi)
    dy = random_normal()*stepsize_y
    y_new = y_current(xi) + dy       

    call calc_pdf(y_current(xi), y_hparam1, y_hparam2, y_pdf, p0_y)
    call calc_pdf(y_new, y_hparam1, y_hparam2, y_pdf, p1_y)
	
    !n1=matmul(h_agg2,x1)-y 
    n1(xi)=n0(xi)-dy

    dumx = n0(xi)*Rinv(:,xi)
    dumy = n1(xi)*Rinv(:,xi)

   dum  = dot_product(n0,dumx)  
   dumb = dot_product(n1,dumy)
   dumc = dot_product(n0,Rinv(:,xi))
   dumd = dot_product(n1,Rinv(:,xi))

   n1T = n0T - dum + dumb - n0(xi)*dumc + n1(xi)*dumd - n1(xi)*n1(xi)*Rinv(xi,xi) + n0(xi)*n0(xi)*Rinv(xi,xi)

    !C = matmul(n1,Rinv)
    !n1T= dot_product(n1,C)
    
    ! Compute P1/P0 	
    pT= p1_y-p0_y - 0.5*(n1T - n0T)*beta       !*beta      

    if (y_pdf .eq. 1) then
       if (y_new .lt. y_hparam1) pT = -1.e20
       if (y_new .gt. y_hparam2) pT = -1.e20
    endif
  
    call random_number(randomu)     ! Generates uniformly distributed random number
    
    if(alog(randomu) .le. pT) then      
       !ACCEPT	
       y_current(xi) = y_new
       n0=n1
       n0T=n1T
       if(beta .eq. 1. .and. it .gt. burn_in) accept=accept + 1
    else
       !;REJECT					
       if(beta .eq. 1. .and. it .gt. burn_in) reject=reject + 1
    endif

  enddo   ! yi loop

n0_out=n0
n0T_out=n0T
y_out=y_current
accept_out=accept
reject_out=reject
END SUBROUTINE y_update


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

subroutine init_random_seed()
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid, t(2), s, getpid
            integer(8) :: count, tms
           
            call random_seed(size = n)
            allocate(seed(n))

            ! First try if the OS provides a random number generator
            open(unit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(count)
               if (count /= 0) then
                  t = transfer(count, t)
               else
                  call date_and_time(values=dt)
                  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
                  t = transfer(tms, t)
               end if
               s = ieor(t(1), t(2))
               pid = getpid() + 1099279 ! Add a prime
               s = ieor(s, pid)
               if (n >= 3) then
                  seed(1) = t(1) + 36269
                  seed(2) = t(2) + 72551
                  seed(3) = pid
                  if (n > 3) then
                     seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                  end if
               else
                  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
               end if
            end if
            call random_seed(put=seed)
end subroutine init_random_seed

End Module transd_evencorr


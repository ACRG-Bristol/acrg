
SUBROUTINE hbtdmcmc(beta,k, x, h_agg,y,n0,n0T, plon, plat, regions_v, &
pdf_param1, pdf_param2, lon,lat, h_v, sigma_y, sigma_model, sigma_measure, &
R_indices, sigma_model_hparams, stepsize_sigma_y, sigma_model_pdf, &
sigma_clon, sigma_clat, & 
lonmin, lonmax, latmin,latmax, sigma_bd, kmin, x_pdf_all, burn_in, &
pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, pdf_param1_pdf, pdf_param2_pdf, &
stepsize, stepsize_pdf_p1,stepsize_pdf_p2, nIt, nsub, nit_sub, nIC, &
nbeta, kmax, kICmax, nmeasure, Ngrid, nlon,nlat, ydim1, ydim2, nIC1, &
k_out, x_out, regions_out, plon_out, plat_out, sigma_y_out, n0T_out, &
pdf_param1_out, pdf_param2_out, accept, reject, &
accept_birth, reject_birth, accept_death, reject_death, accept_move, reject_move, &
accept_swap, accept_sigma_y, reject_sigma_y)



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
! Single Variables
REAL lonmin
REAL lonmax
REAL latmin
REAL latmax
REAL sigma_bd
REAL sigma_clon
REAL sigma_clat
REAL stepsize_sigma_y
INTEGER sigma_model_pdf
INTEGER pdf_param1_pdf
INTEGER pdf_param2_pdf
! Input arrays
REAL stepsize(nIC1)
REAL stepsize_pdf_p1(nIC1)
REAL stepsize_pdf_p2(nIC1)
INTEGER x_pdf_all(nIC1)
REAL beta(nbeta)
INTEGER k(nbeta)
REAL x(kICmax, nbeta)               
REAL pdf_param1(nIC1,nbeta)            
REAL pdf_param2(nIC1,nbeta)                
REAL pdf_p1_hparam1(nIC1)
REAL pdf_p1_hparam2(nIC1)
REAL pdf_p2_hparam1(nIC1)
REAL pdf_p2_hparam2(nIC1)
REAL h_agg(nmeasure,kICmax,nbeta)
REAL y(nmeasure) 
REAL n0(nmeasure, nbeta) 
REAL n0T(nbeta)
REAL sigma_y(nmeasure,nbeta)
REAL sigma_model(ydim2, nbeta)
INTEGER R_indices(ydim1,ydim2)
REAL sigma_measure(nmeasure)
REAl sigma_model_hparams(2)
REAL plon(kmax,nbeta)
REAL plat(kmax,nbeta)
INTEGER regions_v(Ngrid,nbeta)
REAL lon(nlon)
REAL lat(nlat)
REAL h_v(nmeasure, Ngrid)
! Outputs
INTEGER k_out(nit_sub)
REAL x_out(kICmax,nit_sub)  
REAL pdf_param1_out(nIC1,nit_sub)             
REAL pdf_param2_out(nIC1,nit_sub)   
REAL plon_out(kmax,nit_sub)
REAL plat_out(kmax,nit_sub)
INTEGER regions_out(Ngrid,nit_sub)
REAL sigma_y_out(nmeasure, nit_sub)
REAL n0T_out(nit_sub)
INTEGER accept(nIC1)
INTEGER reject(nIC1)
INTEGER accept_birth, reject_birth
INTEGER accept_death, reject_death, accept_move, reject_move, accept_swap
INTEGER accept_sigma_y, reject_sigma_y
! INTERMEDIATE VARIABLES
INTEGER it, ibeta, remain_it, pair1,pair2, ib, it_sub, remain, kIC       !remain_dim
INTEGER remain_swap
REAL u, u1,u2, randomu,pT_chain, beta1,beta2
INTEGER k_it(nit_sub)
REAL x_it(kICmax,nit_sub)                             
REAL pdf_param1_it(nIC1,nit_sub)   
REAL pdf_param2_it(nIC1,nit_sub)   
REAL plon_it(kmax,nit_sub)               
REAL plat_it(kmax,nit_sub)              
REAL sigma_model_ap(ydim2)            
INTEGER regions_it(Ngrid,nit_sub)
REAL sigma_y_it(nmeasure,nit_sub)
REAL n0T_it(nit_sub)
! SUBROUTINE INPUTS
REAL betaib, n0Tib
REAL xib(kICmax), plonib(kmax), platib(kmax), n0ib(nmeasure)           
REAL pdf_param1ib(nIC1), pdf_param2ib(nIC1)
REAL h_aggib(nmeasure,kICmax)
REAL sigma_yib(nmeasure), sigma_modelib(ydim2)
INTEGER kib
INTEGER regions_vib(Ngrid)
! SUBROUTINE OUTPUTS
REAL n0Tib1                     
REAL xib1(kICmax), plonib1(kmax), platib1(kmax), n0ib1(nmeasure)      
REAL pdf_param1ib1(nIC1), pdf_param2ib1(nIC1)
REAL h_aggib1(nmeasure,kICmax)
REAL sigma_yib1(nmeasure), sigma_modelib1(ydim2)
INTEGER kib1, rejectib1, acceptib1, reject_yib1, accept_yib1
INTEGER acceptxib1(nIC1), rejectxib1(nIC1)
INTEGER regions_vib1(Ngrid)
REAL detval(nbeta)
REAL detvalib, detvalib1

! OTHER INTERMEDIATES

!f2py intent(in) beta,k, x, h_agg,y,n0,n0T, plon, plat, regions_v
!f2py intent(in) pdf_param1, pdf_param2, lon,lat, h_v, sigma_clon, sigma_clat  
!f2py intent(in) sigma_y, sigma_model, sigma_measure
!f2py intent(in) R_indices, sigma_model_hparams, stepsize_sigma_y, sigma_model_pdf
!f2py intent(in) lonmin, lonmax, latmin,latmax, sigma_bd, kmin, x_pdf_all, burn_in
!f2py intent(in) pdf_p2_hparam1, pdf_p2_hparam2, stepsize_pdf_p2, pdf_param2_pdf
!f2py intent(in) pdf_p1_hparam1, pdf_p1_hparam2, stepsize_pdf_p1, pdf_param1_pdf
!f2py intent(in) stepsize, nIt,nsub,nit_sub, nIC1
!f2py intent(in) nIC, nbeta, kmax, kICmax, nmeasure, Ngrid, nlon,nlat, ydim1, ydim2
!f2py intent(out) k_out, x_out, regions_out, plon_out, plat_out, sigma_y_out, n0T_out
!f2py intent(out) pdf_param2_out, pdf_param1_out
!f2py intent(out) accept, reject, accept_swap, accept_birth, accept_death, accept_move
!f2py intent(out) reject_birth, reject_death, reject_move, accept_sigma_y, reject_sigma_y

  call OMP_SET_NUM_THREADS(nbeta)

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
it_sub=1


sigma_model_ap=sigma_model(:,1)

detval(:) = sum(alog(sigma_y(:,1)))



! MCMC loop
!###############################################################
do it=1,(nIt+burn_in)
   
   call random_number(u) 
   remain_it = FLOOR(5*u) + 1    ! Choose random number between 1 and 5 to chose what to update

!$OMP PARALLEL DO DEFAULT(SHARED) private(ibeta, betaib,kib,xib,pdf_param1ib,pdf_param2ib, plonib, platib), &
!$OMP& private(regions_vib,h_aggib,n0ib,n0Tib,kIC,xib1,n0ib1,n0Tib1,acceptib1,rejectib1,regions_vib1), &
!$OMP& private(u, kib1, h_aggib1, plonib1, platib1, acceptxib1, rejectxib1), &
!$OMP& private(sigma_yib, sigma_modelib, sigma_yib1, sigma_modelib1, accept_yib1, reject_yib1), &
!$OMP& private(pdf_param1ib1, pdf_param2ib1, detvalib, detvalib1),&
!$OMP& shared(x,n0,n0T, k, pdf_param1, pdf_param2, h_agg, plon,plat, regions_v)
   do ibeta=1,nbeta

       !ibeta=1
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
       detvalib = detval(ibeta)

       kIC = kib+nIC
       
       
       if (remain_it .EQ. 1) then              ! X UPDATE

            call x_update(betaib,kib, xib, pdf_param1ib,pdf_param2ib, &
                         h_aggib,n0ib,n0Tib,sigma_yib, stepsize, &
                         pdf_p1_hparam1, pdf_p1_hparam2, stepsize_pdf_p1, pdf_param1_pdf, &
                         pdf_p2_hparam1, pdf_p2_hparam2, stepsize_pdf_p2, pdf_param2_pdf, &
                         accept,reject,x_pdf_all, it, burn_in, nIC, kICmax, nmeasure, nIC1, &
                         xib1, pdf_param1ib1, pdf_param2ib1, n0ib1, n0Tib1, acceptxib1, rejectxib1)


            x(:,ibeta) = xib1
            n0(:,ibeta) = n0ib1
            n0T(ibeta) = n0Tib1 
            pdf_param1(:,ibeta) = pdf_param1ib1
            pdf_param2(:,ibeta) = pdf_param2ib1
            
            if (betaib .EQ. 1.) then 
               accept(:) = acceptxib1
               reject(:) = rejectxib1
            endif
            
           

       elseif (remain_it .EQ. 5) then       ! BIRTH
               
               call birth(betaib,kib, xib, h_aggib,y,n0ib,n0Tib,sigma_yib, plonib, platib, regions_vib, lon,lat, & 
                          h_v, pdf_param1ib(nIC1), pdf_param2ib(nIC1), x_pdf_all(nIC1), &
                          lonmin, lonmax, &
                          latmin,latmax, sigma_bd,accept_birth, reject_birth,it,burn_in,nIC,kICmax,kmax, &
                          nmeasure,Ngrid,nlon,nlat, &
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

           elseif (remain_it .EQ. 3) then    ! DEATH

               call death(betaib,kib, xib, h_aggib, y,n0ib,n0Tib,sigma_yib, plonib, platib, regions_vib, lon,lat, & 
                          h_v, pdf_param1ib(nIC1), pdf_param2ib(nIC1), x_pdf_all(nIC1), sigma_bd, &
                          accept_death, reject_death, it, burn_in,nIC, kICmax, kmin, kmax, nmeasure, &
                          Ngrid,nlon,nlat, &
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

           elseif (remain_it .EQ. 4) then    ! MOVE
                
               
               call move(betaib,kib, xib, h_aggib, y,n0ib,n0Tib,sigma_yib, plonib, platib, regions_vib, lon,lat, & 
                         h_v, lonmin, lonmax, latmin,latmax, sigma_clon, sigma_clat, accept_move, reject_move, it, &
                         burn_in, nIC, kICmax, kIC, kmax, nmeasure, Ngrid,nlon,nlat, &
                         h_aggib1, n0ib1, n0Tib1, regions_vib1, plonib1, platib1, acceptib1, rejectib1)
               
               plon(:,ibeta) = plonib1
               plat(:,ibeta) = platib1
               regions_v(:,ibeta) = regions_vib1
               h_agg(:,:,ibeta) = h_aggib1
               n0(:,ibeta) = n0ib1
               n0T(ibeta) = n0Tib1
               
               
               if (betaib .EQ. 1.) then 
                  accept_move=acceptib1
                  reject_move=rejectib1
               endif

          elseif (remain_it .EQ. 2) then  ! SIGMA_Y UPDATE
             
              call sigma_y_update(betaib, sigma_modelib, sigma_model_ap, sigma_measure, sigma_yib, detvalib, &
                 sigma_model_hparams, stepsize_sigma_y, sigma_model_pdf, R_indices, &
                 n0ib,n0Tib, accept_sigma_y, reject_sigma_y, it, burn_in, nmeasure, ydim1, ydim2, &
                 n0Tib1, accept_yib1, reject_yib1, sigma_yib1, sigma_modelib1, detvalib1) 

              sigma_y(:,ibeta) = sigma_yib1
              sigma_model(:,ibeta) = sigma_modelib1
              n0T(ibeta) = n0Tib1
              detval(ibeta) = detvalib1 

              if (betaib .EQ. 1.) then 
               accept_sigma_y = accept_yib1
               reject_sigma_y = reject_yib1
              endif

           endif     ! remain_it
           
   enddo    ! beta loop
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
            pT_chain = (beta2-beta1)*(n0T(pair2)/2.-n0T(pair1)/2.+detval(pair2)-detval(pair1))  ! detvals should be inverse determinants so signs this way round
            call random_number(randomu)
            if (alog(randomu) .LE. pT_chain) then
                beta(pair2)=beta1*1.
                beta(pair1)=beta2*1.
                accept_swap=accept_swap+1
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
                  sigma_y_it(:,it_sub)=sigma_y(:,ib)
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
sigma_y_out=sigma_y_it
n0T_out=n0T_it
pdf_param1_out=pdf_param1_it
pdf_param2_out=pdf_param2_it

END SUBROUTINE hbtdmcmc



SUBROUTINE x_update(beta,k, x, pdf_param1_all,pdf_param2_all,  &
h_agg,n0,n0T,sigma_y, stepsize, &
pdf_p1_hparam1_all, pdf_p1_hparam2_all, stepsize_pdf_p1, pdf_param1_pdf, &
pdf_p2_hparam1_all, pdf_p2_hparam2_all, stepsize_pdf_p2, pdf_param2_pdf, &
accept, reject, x_pdf_all, it, burn_in, nIC, kICmax, nmeasure, nIC1, &
x_out, pdf_param1_out, pdf_param2_out, n0_out, n0T_out, accept_out, reject_out) 

Implicit none 
INTEGER nmeasure, it, burn_in, k, nIC, kICmax, nIC1
REAL beta, n0T, n0T_out 
REAL x(kICmax) 
REAL x_out(kICmax) 
REAL h_agg(nmeasure,kICmax)   
REAL n0(nmeasure) 
REAL n0_out(nmeasure) 
REAL sigma_y(nmeasure)
REAL dy(nmeasure)
REAL n1(nmeasure) 

INTEGER x_pdf_all(nIC1)
REAL pdf_param1_all(nIC1), pdf_param2_all(nIC1)
REAL pdf_p1_hparam1_all(nIC1), pdf_p1_hparam2_all(nIC1) 
REAL  pdf_p2_hparam1_all(nIC1), pdf_p2_hparam2_all(nIC1)
REAL pdf_param1_out(nIC1), pdf_param2_out(nIC1)
INTEGER, DIMENSION(2) :: elem
REAL pdf_param1        
REAL pdf_param2
REAL stepsize(nIC1)
INTEGER accept(nIC1), reject(nIC1)
INTEGER accept_out(nIC1), reject_out(nIC1)
REAL stepsize_pdf_p2(nIC1), stepsize_pdf_p1(nIC1)
INTEGER xi, x_pdf, ki,xx
REAL dx, n1T, pT, randomu, random_normal, p0,p1, u
INTEGER pdf_param2_pdf, pdf_param1_pdf
REAL pdf_p2_hparam1, pdf_p2_hparam2, pdf_p1_hparam1, pdf_p1_hparam2
REAL p0_temp, p1_temp
REAL dpdf_param2, pdf_param2_new, dpdf_param1, pdf_param1_new
REAL p0_pdf_param2, p1_pdf_param2, p0_pdf_param1, p1_pdf_param1
REAL x_new(kICmax)
REAL stepsize0, stepsize_pdf_p10, stepsize_pdf_p20

   ! CHANGE OF EMISSIONS VALUES
  call random_number(u)   
  elem(1)= FLOOR((nIC)*u)+1 
  elem(2) = FLOOR(k*u)+1+nIC
  !xi=FLOOR((nIC+k)*u)+1

  do xx=1,2
  xi=elem(xx)
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
  dx = random_normal()*stepsize0
  dpdf_param1 = random_normal()*stepsize_pdf_p10
  pdf_param1_new = pdf_param1 + dpdf_param1

  dpdf_param2 = random_normal()*stepsize_pdf_p20
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
          n1T=sum((n1/sigma_y)**2)
                  
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

         n1T=sum((n1/sigma_y)**2)
          
         ! hyperparams are fixed below to be a single number - will only apply when x is a scaling of the prior

         call calc_pdf(pdf_param1,pdf_p1_hparam1,pdf_p1_hparam2,pdf_param1_pdf, p0_pdf_param1) 
         call  calc_pdf(pdf_param1_new, pdf_p1_hparam1,pdf_p1_hparam2,pdf_param1_pdf, p1_pdf_param1)

         call calc_pdf(pdf_param2,pdf_p2_hparam1,pdf_p2_hparam2,pdf_param2_pdf, p0_pdf_param2) 
         call  calc_pdf(pdf_param2_new, pdf_p2_hparam1,pdf_p2_hparam2,pdf_param2_pdf, p1_pdf_param2)

         if (xi .LE. nIC) then
           
            do ki =1,nIC
               call calc_pdf(x(ki),pdf_param1,pdf_param2,x_pdf, p0_temp)           ! Will apply whatever the PDF
               call calc_pdf(x_new(ki),pdf_param1_new,pdf_param2_new,x_pdf, p1_temp)        
               p0=p0+p0_temp
               p1=p1+p1_temp
            enddo

         else if (xi .GT. nIC) then

            do ki =nIC+1,k+nIC
               call calc_pdf(x(ki),pdf_param1,pdf_param2,x_pdf, p0_temp)           ! Will apply whatever the PDF
               call calc_pdf(x_new(ki),pdf_param1_new,pdf_param2_new,x_pdf, p1_temp)        
               p0=p0+p0_temp
               p1=p1+p1_temp
            enddo
         endif


          pT = p1+p1_pdf_param2+p1_pdf_param1-p0-p0_pdf_param2-p0_pdf_param1 -0.5*(n1T - n0T)*beta ! All p1s already in log space

         if (pdf_param2_pdf .eq. 1) then
             if (pdf_param2_new .lt. pdf_p2_hparam1) pT = -1.e20
             if (pdf_param2_new .gt. pdf_p2_hparam2) pT = -1.e20
         else if (pdf_param1_pdf .eq. 1) then
             if (pdf_param1_new .lt. pdf_p1_hparam1) pT = -1.e20
             if (pdf_param1_new .gt. pdf_p1_hparam2) pT = -1.e20
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
                              
         endif   ! randomu condition

  endif   ! x_pdf condition
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

x_out=x
n0_out=n0
n0T_out=n0T
pdf_param1_out=pdf_param1_all
pdf_param2_out=pdf_param2_all
accept_out=accept
reject_out=reject
END SUBROUTINE x_update



SUBROUTINE birth(beta,k, x, h_agg,y,n0,n0T,sigma_y, plon, plat, regions_v, lon,lat, & 
h_v,pdf_param1, pdf_param2, x_pdf,  &
lonmin, lonmax, latmin,latmax, sigma_bd, &
accept_birth, reject_birth, it, burn_in, nIC, kICmax, kmax, nmeasure, Ngrid,nlon,nlat, &
k_out, x_out, h_agg_out, n0_out, n0T_out, regions_v_out, plon_out, plat_out, accept_out, reject_out)

IMPLICIT NONE
! Dimensions
INTEGER nmeasure, Ngrid, nlon,nlat, k,kmax, nIC, kICmax
REAL lonmin,lonmax,latmin,latmax,sigma_bd
! Single Variables
INTEGER accept_birth, reject_birth, it, burn_in, x_pdf 
REAL beta, n0T
! Input arrays
REAL x(kICmax) 
REAL pdf_param1 
REAL pdf_param2
REAL h_agg(nmeasure,kICmax)  
REAL y(nmeasure) 
REAL n0(nmeasure) 
REAL sigma_y(nmeasure)
REAL plon(kmax)
REAL plat(kmax)
INTEGER regions_v(Ngrid)
REAL lon(nlon)
REAL lat(nlat)
REAL h_v(nmeasure, Ngrid)
! Outputs
INTEGER k_out
REAL x_out(kICmax) 
REAL h_agg_out(nmeasure,kICmax)
REAL n0_out(nmeasure) 
REAL n0T_out
INTEGER regions_v_out(Ngrid)
REAL plon_out(kmax)
REAL plat_out(kmax)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER ri, rib, k1, errstat,jj, kIC      
REAL u, u2, plon_new, plat_new, c_new, x_new, n1Tb, pT_birth, randomu, random_normal
REAL mu, sigma
INTEGER regions_v1b(Ngrid)       
REAL n1b(nmeasure)            
! Allocatable arrays
REAL, DIMENSION(:),ALLOCATABLE :: plon1b, plat1b, x1b
REAL, DIMENSION(:,:), ALLOCATABLE :: h_agg2
REAL,PARAMETER     :: pi = 3.14159265 


k1=k+1
kIC=k1+nIC

if (k1 .LT. kmax) THEN            
  
  allocate(plon1b(k1))     
  allocate(plat1b(k1))     
  allocate(h_agg2(nmeasure, kIC))
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
          n1Tb=sum((n1b/sigma_y)**2)

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
          n1Tb=sum((n1b/sigma_y)**2)

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


SUBROUTINE death(beta,k, x, h_agg,y,n0,n0T,sigma_y, plon, plat, regions_v, lon,lat, & 
h_v, pdf_param1, pdf_param2, x_pdf, &
sigma_bd, accept_death, reject_death, &
it, burn_in, nIC, kICmax, kmin, kmax, nmeasure, Ngrid,nlon,nlat, &
k_out, x_out, h_agg_out, n0_out, n0T_out, regions_v_out, &
plon_out, plat_out, accept_out, reject_out)



IMPLICIT NONE
! Dimensions
INTEGER kmax,nmeasure,Ngrid,nlon,nlat, accept_death, reject_death, it, burn_in, kmin, k, nIC, kICmax
REAL sigma_bd, beta, n0T 
! Input arrays
REAL x(kICmax) 
REAL pdf_param1
REAL pdf_param2
REAL h_agg(nmeasure,kICmax) 
REAL y(nmeasure) 
REAL n0(nmeasure) 
REAL sigma_y(nmeasure)
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
REAL h_agg_out(nmeasure,kICmax)
REAL n0_out(nmeasure) 
INTEGER regions_v_out(Ngrid)
REAL plon_out(kmax)
REAL plat_out(kmax)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER ri, rid, k1d, jj, ci_rm, errstat, kIC     
REAL u, plon_rm, plat_rm, x_cell, x_rm, n1Td, pT_death, randomu
REAL mu,sigma
INTEGER regions_v1d(Ngrid)      
REAL n1d(nmeasure)             
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
  allocate(h_agg2d(nmeasure, kIC))
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
  n1Td=sum((n1d/sigma_y)**2)                           
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


SUBROUTINE move(beta,k, x, h_agg, y,n0,n0T,sigma_y, plon, plat, regions_v, lon,lat, & 
h_v, lonmin, lonmax, latmin,latmax, sigma_clon, sigma_clat, accept_move, reject_move, it, &
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
REAL sigma_y(nmeasure)
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
REAL u, n1Tm, pT_move, randomu, random_normal
INTEGER regions_v1m(Ngrid)      
REAL n1m(nmeasure) 
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
     n1Tm=sum((n1m/sigma_y)**2)
                                     
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

SUBROUTINE sigma_y_update(beta, sigma_model_current, sigma_model_ap, sigma_measure, sigma_y_current, &
detval_current,sigma_model_hparams, stepsize_sigma_y, sigma_model_pdf, R_indices, &
n0,n0T, accept, reject, it, burn_in, nmeasure, dim1, dim2, &
n0T_out, accept_out, reject_out, sigma_y_out, sigma_model_out, detval_out) 

IMPLICIT NONE

! Dimensions
INTEGER nmeasure, accept, reject, it, burn_in, dim1, dim2
! Single Variables
REAL beta 
! Input arrays
REAL n0(nmeasure) 
REAL n0T
REAL detval_current
REAL sigma_measure(nmeasure)
REAL sigma_model_current(dim2)
REAL sigma_model_ap(dim2)
REAl sigma_y_current(nmeasure)
INTEGER R_indices(dim1,dim2)
REAL stepsize_sigma_y
REAL sigma_model_hparams(2)
INTEGER sigma_model_pdf
! Outputs
REAL n0T_out, detval_out
REAL sigma_y_out(nmeasure)
REAL sigma_model_out(dim2)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER  yi, jj
REAL randomu, random_normal, dsigma_y, sigma_model_new, u
REAL p0_sigma_y, p1_sigma_y, n1T, detval_new, pT
REAL y_error_new(nmeasure)
REAL sigma_y_new(nmeasure)



 !do yi=1,dim2
    call random_number(u)   
    yi = FLOOR(dim2*u)+1
		
    ! Generate new value of sigma_y
    dsigma_y = random_normal()*stepsize_sigma_y*sigma_model_ap(yi)
    sigma_model_new = sigma_model_current(yi) + dsigma_y       

    call calc_pdf(sigma_model_current(yi), sigma_model_hparams(1), sigma_model_hparams(2), sigma_model_pdf, p0_sigma_y)
    call calc_pdf(sigma_model_new, sigma_model_hparams(1), sigma_model_hparams(2), sigma_model_pdf, p1_sigma_y)
	
    do jj=1,dim2   
        y_error_new(R_indices(:,jj)) = sigma_model_current(jj)    ! Provided dim2 isn't too big then should be fine
    enddo  

    ! Change one element of sigma_y
    y_error_new(R_indices(:,yi)) = sigma_model_new  ! R_indices = array of indices to which each sigma_y applies 
   			   
    sigma_y_new=sqrt(y_error_new**2 + sigma_measure**2)   

    n1T=sum((n0/sigma_y_new)**2)

    detval_new = sum(alog(sigma_y_new))           ! This is actually the inverse of the determinant

    ! Compute P1/P0 	
    pT= p1_sigma_y-p0_sigma_y - detval_new*beta + detval_current*beta - 0.5*(n1T - n0T)*beta       !*beta      ! -detval_new becasue it's inverse of true determinant
    if (sigma_model_pdf .eq. 1) then
       if (sigma_model_new .lt. sigma_model_hparams(1)) pT = -1.e20
       if (sigma_model_new .gt. sigma_model_hparams(2)) pT = -1.e20
    endif
  
    call random_number(randomu)     ! Generates uniformly distributed random number
    
    if(alog(randomu) .le. pT) then      
       !ACCEPT	
       sigma_model_current(yi) = sigma_model_new
       sigma_y_current = sigma_y_new
       detval_current = detval_new
       n0T=n1T
       if(beta .eq. 1. .and. it .gt. burn_in) accept=accept + 1
    else
       !;REJECT					
       if(beta .eq. 1. .and. it .gt. burn_in) reject=reject + 1
    endif

n0T_out=n0T
sigma_y_out=sigma_y_current
sigma_model_out=sigma_model_current
accept_out=accept
reject_out=reject
detval_out=detval_current
END SUBROUTINE sigma_y_update




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


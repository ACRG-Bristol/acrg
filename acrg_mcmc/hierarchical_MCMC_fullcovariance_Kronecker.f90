!  Purpose:
!  Hierarchical Bayesian inversion using Markov Chain Monte Carlo 
!  
!  Inputs: see below inputs to be stored in .nc file
!   
!  Outputs:
!  Must create a folder called 'outputs' in the current directory where a file 'output_file.nc' will be created (./outputs/output_file.nc) 
!
!  Requires:
!  Netcdf file *input_file.nc*, stored in folder called 'inputs' in current directory (./inputs/input_file.nc)
!   
!  Netcdf & LAPACK libraries
!   
!  Modules Used: mod_MCMC.f90
!   
!  
!  Written by: A.Ganesan, University of Bristol, July 18 2014

!  Calculates z = Cy + e1 where diag(e1) = D (measurement uncertainty, nugget term) (assumed Gaussian)
!  y = Hx + e2 where corr(e2) = the Kronecker product of spatial (S) and temporal (T) correlation matrices (assumed Gaussian)
!  Samples through y to find values for missing observations
!  Solves for separate variances for time and space and modifies determinant to be det(SXT)=det(S)^nmeasuretotal*det(T)^numsites
!  Exploits properties of the kronecker product such that the matrix-vector multiplication C = Rinv*n0 can be rewritten as
!  	Rinv*n0 = (Ainv kron Binv)*n0 = Y = Binv*N0_arr*Ainv' where N0_arr = reshape(n0,numsites,nmeasuretotal) and C = reshape(Y,nmeasuretotal*numsites,1)
!  This code never only requires explicit computation of the Kronecker product, inverse and determinant once in the code
!  This code is ideally suited to deal with in situ data where the fraction of missing observations is small, as it requires model output for all times, even when no observations exist.

Program hierarchical_MCMC

use mod_MCMC ! Module contains pdf subroutine, randomn function and random_seed subroutine, inverse and (log sqrt) determinant fuctions, Kronecker product (not used in this version), gamma function, bessel
use netcdf

Implicit none
 
!********************* INPUTS [dimensions] ********************
! Dimensions:
! STATESIZE  		= # of elements in emissions state vector, x_ap
! NMEASURE   		= # of elements in measurement vector, z 
! DIM1 			= maximum # of measurements that correspond to an element of the temporal model-measurement uncertainty (e.g., 100 if 100 observations out of nmeasuretotal correspond to element 1, and 80 correspond to element 2, etc)
! DIM2 			= # of elements of temporal model-measurement uncertainty estimated in the inversion (e.g., daytime and nighttime uncertainties over the month = 2 elements)
! DIM3 			= # of correlation parameters (e.g., tau or nu or rho) to be estimated for each of the spatial and temporal matrices (currently only value allowed is 1)
! NMEASURETOTAL		= maximum possible number of measurements during the time period of interest (e.g., 372 measurements in July if estimated every 2 hours). 
! NUMSITES		= number of measurement sites in the inversion
!			  NMEASUREMAX = NMEASURETOTAL*NUMSITES is computed in the code
! DIM4			= number of gridcells in inversion domain

! Variables:
! NIT        		= Number of iterations [scalar]
! BURN_IN   	 	= Length of burn-in chain [scalar]
! z          		= Vector of observations [nmeasure]
! H         	 	= sensitvity matrix [nmeasuremax x statesize]
! D         	 	= measurement uncertainty (nugget) matrix, assuming uncorrelated observations [nmeasure x nmeasure]

! X_AP       		= Prior parameter estimate vector [statesize]
! PDF_PARAM1 		= Prior of primary parameter describing PDF (e.g. Mean) [statesize]
! PDF_PARAM2 		= Prior of secondary parameter describing PDF (e.g. Standard Deviation) [statesize]
! SIGMA_Y 		= Prior of model-measurement uncertainty in time [dim2]
! SIGMA_YS 		= Prior of model-measurement uncertainty in space [currently fixed at numsites]
! TAU 			= Prior of autocorrelation timescale [dim3]
! NU 			= Prior of Matern shape parameter for spatial autocorrelation [dim3]
! RHO 			= Prior of Matern rate parameter for spatial autocorrelation [dim3]

! STEPSIZE   		= Size of proposal distriubtion for x [statesize]
! STEPSIZE_PDF_PARAM1 	= Size of proposal distriubtion for pdf_param1 [statesize]
! STEPSIZE_PDF_PARAM2 	= Size of proposal distriubtion for pdf_param2 [statesize]
! STEPSIZE_SIGMA_Y 	= Size of proposal distribution for sigma_y [dim2]
! STEPSIZE_SIGMA_YS 	= Size of proposal distribution for sigma_ys [numsites]
! STEPSIZE_TAU 		= Size of proposal distribution for tau [dim3]
! STEPSIZE_NU 		= Size of proposal distribution for nu [dim3]
! STEPSIZE_RHO 		= Size of proposal distribution for rho [dim3]
! STEPSIZE_Y 		= Size of proposal distribution for y [nmeasuremax]

! PDF type parameters ( 0 = LOGNORMAL, 1 = EXPONENTIAL, 2 = GAUSSIAN, 3 = UNIFORM) [statesize]

! X_PDF      		= Type of PDF associated with each element of x[statesize]
! PDF_PARAM1_PDF 	= Type of PDF associated with each element of pdf_param1 [statesize]
! PDF_PARAM2_PDF 	= Type of PDF associated with each element of pdf_param2 [statesize]
! SIGMA_Y_PDF		= Type of PDF associated with each element of sigma_y [dim2]
! SIGMA_YS_PDF		= Type of PDF associated with each element of sigma_ys [numsites]
! TAU_PDF   		= Type of PDF associated with each element of tau [dim3]
! NU_PDF   		= Type of PDF associated with each element of nu [dim3]
! RHO_PDF   		= Type of PDF associated with each element of rho [dim3]

! PDF_HYPERPARAM1 	= Parameter governing PDF of pdf_param1 (e.g., standard deviation in pdf_param1 "hyper-PDF") [statesize]
! PDF_HYPERPARAM2 	= Parameter governing PDF of pdf_param2 (e.g., uncertainty in pdf_param2 "hyper-PDF") [statesize]
! SIGMA_Y_HYPERPARAM	= Parameter governing PDF of sigma_y (e.g., uncertainty in sigma_y "hyper-PDF") [dim2]
! SIGMA_YS_HYPERPARAM	= Parameter governing PDF of sigma_ys (e.g., uncertainty in sigma_ys "hyper-PDF") [numsites]
! TAU_HYPERPARAM1	= Parameter governing PDF of tau [dim3]
! TAU_HYPERPARAM2	= Parameter governing PDF of tau [dim3]
! NU_HYPERPARAM1	= Parameter1 governing PDF of nu [dim3]
! RHO_HYPERPARAM1	= Parameter1 governing PDF of rho [dim3]
! NU_HYPERPARAM2	= Parameter2 governing PDF of nu [dim3]
! RHO_HYPERPARAM2	= Parameter2 governing PDF of rho [dim3]

! T_INDICES 		= Array containing indices of measurements that correspond to each element of sigma_y (e.g., 1st column corresponds to the indices governed by first element of sigma_y, 2nd column  corresponds to indices governed by the second element of sigma_y, and so forth) [dim1 x dim2]

! DELTATIME 		= Array of 'delta_t' in same units as TAU to calculate temporal autocorrelation (e.g., if tau is in days, measurements that occur 3 hours apart would have a value of 0.125 days, 6 hours apart would have a value of 0.5, etc.) [nmeasuretotal x nmeasuretotal]

! DISTANCE 		= Array of Euclidean distance between sites in same units as rho to calculate spatial autocorrelation  [numsites x numsites]

! DATENUMBER		= Vector of serial date number (or equivalent) for the time period (e.g., July 2012 estimated every 3 hours would have 248 datenumbers). This is not used in the MCMC routine but is saved into the outputs for the convenience of analysis. Ignore if desired.  [nmeasuretotal]

! TIMEINDEX_NONZERO	= Vector of indices corresponding to when a measurement was avaiable during time period (1:nmeasuremax).  [nmeasure]

! N_OBS			= Vector of number observations per site (not used in analysis)  [numsites]

! KRON_FLAG             = Flag for orientation of y,H and kronecker product. Value = 0 if spatially stacked, = 1 if temporally stacked. [Scalar]

! NUMTHREADS		= Number of threads to use for multithreaded applications

! SITESPRESENT		= Index corresponding to which sites are present using whatever numbering system the user has implemented (not used in analysis but important for tracking which sites are used). For example if Mace Head = 1, Ridge Hill = 2, Tacolneston = 3, Angus = 4 then if only MHD and TAC are available for the analysis then this would be [1 3] [numsites]

! DISTRIBUTION		= Percentage of emissions from grid cell contained in larger region that has been aggregated for the inversion (not used in the analysis) [dim4]

!! *************** OUTPUTS [dimensions] ***************************
! X_IT       		= Array containing value of x at each iteration (nIt x statesize) 
! PDF_PARAM1_IT		= Array containing value of pdf_param1 at each iteration (nIt x statesize) 
! PDF_PARAM2_IT		= Array containing value of pdf_param2 at each iteration (nIt x statesize) 
! SIGMA_Y_IT		= Array containing value of sigma_y at each iteration (nIt x dim2) 
! SIGMA_YS_IT		= Array containing value of sigma_ys at each iteration (nIt x numsites) 
! TAU_IT		= Array containing value of tau at each iteration (nIt x dim3) 
! NU_IT			= Array containing value of nu at each iteration (nIt x dim3) 
! RHO_IT		= Array containing value of rho at each iteration (nIt x dim3) 
! Y_IT			= Array containing value of y at each iteration (nIt x nmeasuremax) 
! ACCEPTANCE     	= Vector of acceptance ratio for each element of x/pdf_param1/pdf_param2 [statesize] 
! ACCEPTANCE_SIGMAY	= Vector of acceptance ratio for each element of sigma_y [dim2] 
! ACCEPTANCE_SIGMAYS	= Vector of acceptance ratio for each element of sigma_ys [numsites] 
! ACCEPTANCE_TAU	= Vector of acceptance ratio for each element of tau [dim3] 
! ACCEPTANCE_NU		= Vector of acceptance ratio for each element of nu and rho [dim3] 
! ACCEPTANCE_Y		= Vector of acceptance ratio for each element of y [nmeasuremax] 
! H			= Same as input - restored for convenience of analyzing output
! Z			= Same as input - restored for convenience of analyzing output
! DATENUMBER		= Same as input - restored for convenience of analyzing output
! TIMEINDEX_NONZERO	= Same as input - restored for convenience of analyzing output
! X_AP			= Same as input - restored for convenience of analyzing output
! N_OBS			= Same as input - restored for convenience of analyzing output
! sitespresent		= Same as input - restored for convenience of analyzing output
! deltatime		= Same as input - restored for convenience of analyzing output
! distance		= Same as input - restored for convenience of analyzing output
! T_indices		= Same as input - restored for convenience of analyzing output
! DISTRIBUTION		= Same as input - restored for convenience of analyzing output
! D			= Same as input - restored for convenience of analyzing output

! Variable declarations

! dimensions
  integer                           :: 	nIt, statesize, nmeasure, nmeasuretotal, nmeasuremax, numsites, burn_in, dim1, dim2, dim3, kron_flag, numthreads, dim4
                                                                               
! counters                                       
  integer                           :: 	ii, jj, xi, yi, It, status, errstat, count_0, count_1, count_max, tstart, tend, &
					reject, reject_sigma_y, reject_tau , reject_nu, reject_y, reject_sigma_ys, omp_get_num_threads

! input and output variable and dimension IDs         
  integer                           :: 	inputID, outputID, burn_inID, nItID, statesizeID, nMeasureID, x_apID, x_pdfID, pdf_param1ID, pdf_param2ID, stepsizeID, & 
                                       	zID, HID, x_itID, y_itID, acceptID, statesizeDID, nMeasureDID, nItDID, stepsizepdfparam1ID, stepsizepdfparam2ID, &
				  	pdf_hyperparam1ID, pdf_hyperparam2ID, pdf_param1_pdfID, pdf_param2_pdfID, sigmayID, TindicesID, sigmayhyperparamID, &
				       	stepsizesigmayID, sigmaypdfID, dim1ID, dim2ID, dim3ID, pdf_param1_itID, pdf_param2_itID, sigma_y_itID, dim1DID, dim2DID, dim3DID, &
     				        acceptsigmayID, acceptyID, deltatimeID, distanceID, tauID, taupdfID, tau_itID, accepttauID, stepsizetauID, tauhyperparam1ID, tauhyperparam2ID, &
				        HoID, zoID, nmeasuretotalID, datenumberID, timeindexnonzeroID, nmeasuretotalDID, datenumberoID, timeindexnonzerooID, &
					stepsizenuID, stepsizerhoID, nuID, rhoID, nuhyperparam1ID, rhohyperparam1ID, nuhyperparam2ID, rhohyperparam2ID, nupdfID, rhopdfID, & 						acceptnuID, nu_itID, rho_itID, numsitesID, DID, stepsizeyID, nmeasuremaxDID, xapoID, nobsID, nobsoID, numsitesDID, &
					sigmaysID, sigmays_itID, acceptsigmaysID, stepsizesigmaysID, sigmayshyperparamID, sigmayspdfID, kron_flagID, & 				 						numthreadsID, sitespresentID,sitespresentoID, deltatimeoID, distanceoID, TindicesoID, dim4ID, distributionID, dim4DID, &
					distributionoID, DoID                           

  integer		            :: 	tau_pdf, nu_pdf, rho_pdf 

  integer,allocatable, dimension(:) :: 	x_pdf, pdf_param1_pdf, pdf_param2_pdf, sigma_y_pdf, sigma_ys_pdf, timeindex_nonzero, nobs, sitespresent, num_T_indices  

  real                              :: 	t0, t1,count_rate, n0T, n1T, m0T, m1T, pT, randomu, &
                                       	p1, p1_x, p1_pdf_param1, p1_pdf_param2, p1_sigma_y, p1_tau, p0_tau, p0_nu, p0_rho, p1_nu, p1_rho, p1_sigma_ys, &
				       	dx, dpdf_param1, dpdf_param2, dsigma_y, dtau,dnu, drho, dymod, dsigma_ys, &
  				       	tau, tau_current, tau_new, stepsize_tau, tau_hyperparam1,tau_hyperparam2, accept_tau, accept_nu, &
					nu, nu_current, nu_new, stepsize_nu, nu_hyperparam1, nu_hyperparam2, &
					rho, rho_current, rho_new, stepsize_rho, rho_hyperparam1, rho_hyperparam2, &
					detval_current,detval_new, detval, detval2, detval_T_current, detval_T_new, detval_S_current, detval_S_new, ga, &
					aa, bb, &			
					dum, dumb, dumc, dumd

  real,allocatable, dimension(:)    :: 	stepsize, stepsize_pdf_param1, stepsize_pdf_param2, stepsize_sigma_y, stepsize_sigma_ys, stepsize_y, &
					p0_x,p0_pdf_param1, p0_pdf_param2, p0_sigma_y, p0_sigma_ys, &
					x , x_ap, z, y, dy, y_small, &
				        pdf_param1, pdf_param2, pdf_param1_current, pdf_param1_new, pdf_param2_current, pdf_param2_new, &
				       	pdf_hyperparam1, pdf_hyperparam2, y_current, y_new, &
                                       	y_error_t, dumx, dumy, dum2, &
                                       	reject_vector, reject_sigma_y_vector, reject_sigma_ys_vector, reject_y_vector, accept_vector, accept_sigma_y_vector, &
					accept_sigma_ys_vector, accept_y_vector, &
				       	sigma_y, sigma_y_current, sigma_y_new, sigma_y_hyperparam, autocorr_vec, &
					sigma_ys, sigma_ys_current, sigma_ys_new, sigma_ys_hyperparam, &
					datenumber, distribution, &
					n0, n1, m0, m1, C, C2

  real,allocatable, dimension(:,:)  :: 	H, D, Dinv, deltatime, distance, x_it, pdf_param1_it, pdf_param2_it, sigma_y_it, tau_it, nu_it, rho_it, y_it, sigma_ys_it, &
					Rinv_current, T_indices, T_current, T_new, S_current, S_new, &
					Sinv_current, Sinv_new, Tinv_current, Tinv_new, n0_arr,n1_arr, C_arr, dum3

  real (kind =8) nu_double, rho_double, alpha, ga_double, arg_dum, k_arg_dum

  real (kind =8), allocatable, dimension(:,:)   :: distance_double, arg

  integer ( kind = 4 ) nb, ize, ncalc

 call system_clock(count_0,count_rate,count_max)
 t0=count_0/count_rate

 ! *************** Open Netcdf input and output files *****************
 STATUS=NF90_OPEN('./inputs/input_file.nc', NF90_NOWRITE, ncid=inputID)
 STATUS = NF90_CREATE ('./outputs/output_file.nc', cmode=0, ncid=outputID)

 ! ************** Read input file *************************

 STATUS = NF90_INQ_VARID (inputID, 'burn_in', burn_inID)
 STATUS = NF90_INQ_VARID (inputID, 'nit', nItID)
 STATUS = NF90_INQ_VARID (inputID, 'statesize', statesizeID)
 STATUS = NF90_INQ_VARID (inputID, 'nmeasure', nMeasureID)
 STATUS = NF90_INQ_VARID (inputID, 'dim1', dim1ID)
 STATUS = NF90_INQ_VARID (inputID, 'dim2', dim2ID)
 STATUS = NF90_INQ_VARID (inputID, 'dim3', dim3ID)
 STATUS = NF90_INQ_VARID (inputID, 'nmeasuretotal', nmeasuretotalID)
 STATUS = NF90_INQ_VARID (inputID, 'numsites', numsitesID)
 STATUS = NF90_INQ_VARID (inputID, 'dim4', dim4ID)

 STATUS = NF90_GET_VAR (inputID, burn_inID, burn_in)
 STATUS = NF90_GET_VAR (inputID, nItID, nIt)
 STATUS = NF90_GET_VAR (inputID, statesizeID, statesize)
 STATUS = NF90_GET_VAR (inputID, nMeasureID, nMeasure)
 STATUS = NF90_GET_VAR (inputID, dim1ID, dim1)
 STATUS = NF90_GET_VAR (inputID, dim2ID, dim2)
 STATUS = NF90_GET_VAR (inputID, dim3ID, dim3)
 STATUS = NF90_GET_VAR (inputID, nmeasuretotalID, nmeasuretotal)
 STATUS = NF90_GET_VAR (inputID, numsitesID, numsites)
 STATUS = NF90_GET_VAR (inputID, dim4ID, dim4)

 nmeasuremax = nmeasuretotal*numsites

 !!***************** Allocate space for each variable ***************************

 allocate(z(nMeasure),stat=errstat)            
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate z'
     	stop
 endif

 allocate(y(nmeasuremax),stat=errstat)            
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate y'
     	stop
 endif

 allocate(nobs(numsites),stat=errstat)            
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate nobs'
     	stop
 endif

 allocate(sitespresent(numsites),stat=errstat)            
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate sitespresent'
     	stop
 endif

 allocate(y_current(nmeasuremax),stat=errstat)            
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate y_current'
     	stop
 endif

 allocate(y_new(nmeasuremax),stat=errstat)            
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate y_new'
     	stop
 endif

 allocate(y_small(nmeasure),stat=errstat)            
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate y_small'
     	stop
 endif

 allocate(n0(nmeasuremax),stat=errstat)           
 if (errstat /= 0) then
    	write (*,*) 'ERROR: could not allocate n0'
     	stop
  endif

 allocate(n1(nmeasuremax),stat=errstat)         
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate n1'
     	stop
 endif

 allocate(m0(nmeasure),stat=errstat)           
 if (errstat /= 0) then
    	write (*,*) 'ERROR: could not allocate m0'
     	stop
  endif

 allocate(m1(nmeasure),stat=errstat)         
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate m1'
     	stop
 endif

 allocate(dy(nmeasuremax),stat=errstat)            
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate dy'
     	stop
 endif

 allocate(H(nmeasuremax,statesize),stat=errstat)           
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate H'
     	stop
 endif

 allocate(D(nMeasure,nMeasure),stat=errstat)           
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate D'
     	stop
 endif

 allocate(Dinv(nMeasure,nMeasure),stat=errstat)           
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate Dinv'
     	stop
 endif

 allocate(T_indices(dim1,dim2),stat=errstat)             
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate T_indices'
     	stop
 endif

 allocate(deltatime(nmeasuretotal,nmeasuretotal),stat=errstat)           
 if (errstat /= 0) then
   	 write (*,*) 'ERROR: could not allocate deltatime'
    	 stop
 endif

 allocate(distance(numsites,numsites),stat=errstat)           
 if (errstat /= 0) then
   	 write (*,*) 'ERROR: could not allocate distance'
    	 stop
 endif

 allocate(distance_double(numsites,numsites),stat=errstat)           
 if (errstat /= 0) then
   	 write (*,*) 'ERROR: could not allocate distance_double'
    	 stop
 endif

 allocate(arg(numsites,numsites),stat=errstat)           
 if (errstat /= 0) then
   	 write (*,*) 'ERROR: could not allocate arg'
    	 stop
 endif

 allocate(x_ap(statesize),stat=errstat)             
 if (errstat /= 0) then
   	write (*,*) 'ERROR: could not allocate x_ap'
    	stop
 endif

 allocate(x(statesize),stat=errstat)             
 if (errstat /= 0) then
 	write (*,*) 'ERROR: could not allocate x'
     	stop
 endif

 allocate(pdf_param1(statesize),stat=errstat)             
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate pdf_param1'
    	stop
 endif

 allocate(pdf_param1_current(statesize),stat=errstat)            
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate pdf_param1_current'
     	stop
 endif

 allocate(pdf_param1_new(statesize),stat=errstat)            
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate pdf_param1_new'
     	stop
 endif

 allocate(pdf_param2(statesize),stat=errstat)             
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate pdf_param2'
     	stop
 endif

 allocate(pdf_param2_current(statesize),stat=errstat)            
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate pdf_param2_current'
     	stop
 endif

 allocate(pdf_param2_new(statesize),stat=errstat)            
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate pdf_param1_new'
     	stop
 endif

 allocate(sigma_y(dim2),stat=errstat)             
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate sigma_y'
     	stop
 endif

 allocate(sigma_y_current(dim2),stat=errstat)             
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate sigma_y_current'
     	stop
 endif

 allocate(sigma_y_new(dim2),stat=errstat)           
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate sigma_y_new'
     	stop
 endif

 allocate(sigma_ys(numsites),stat=errstat)             
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate sigma_ys'
     	stop
 endif

 allocate(sigma_ys_current(numsites),stat=errstat)             
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate sigma_ys_current'
     	stop
 endif

 allocate(sigma_ys_new(numsites),stat=errstat)           
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate sigma_ys_new'
     	stop
 endif

 allocate(stepsize(statesize),stat=errstat)            
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate stepsize'
     	stop
 endif

 allocate(stepsize_pdf_param1(statesize),stat=errstat)             
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate stepsize_pdf_param1'
     	stop
 endif

 allocate(stepsize_pdf_param2(statesize),stat=errstat)             
 if (errstat /= 0) then
     	write (*,*) 'ERROR: could not allocate stepsize_pdf_param2'
     	stop
 endif

allocate(stepsize_sigma_y(dim2),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate stepsize_sigma_y'
     stop
  endif

allocate(stepsize_sigma_ys(numsites),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate stepsize_sigma_ys'
     stop
  endif

allocate(stepsize_y(nmeasuremax),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate stepsize_y'
     stop
  endif

allocate(x_pdf(statesize),stat=errstat)	
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate x_pdf'
     stop
  endif


allocate(pdf_param1_pdf(statesize),stat=errstat)	
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate pdf_param1_pdf'
     stop
  endif

allocate(pdf_param2_pdf(statesize),stat=errstat)	
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate pdf_param2_pdf'
     stop
  endif

allocate(sigma_y_pdf(dim2),stat=errstat)	
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate sigma_y_pdf'
     stop
  endif

allocate(sigma_ys_pdf(numsites),stat=errstat)	
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate sigma_ys_pdf'
     stop
  endif

allocate(pdf_hyperparam1(statesize),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate pdf_hyperparam1'
     stop
  endif

allocate(pdf_hyperparam2(statesize),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate pdf_hyperparam2'
     stop
  endif

allocate(sigma_y_hyperparam(dim2),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate sigma_y_hyperparam1'
     stop
  endif

allocate(sigma_ys_hyperparam(numsites),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate sigma_ys_hyperparam1'
     stop
  endif

allocate(p0_x(statesize),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate p0_x'
     stop
  endif

allocate(p0_pdf_param1(statesize),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate p0_pdf_param1'
     stop
  endif

allocate(p0_pdf_param2(statesize),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate p0_pdf_param2'
     stop
  endif

allocate(p0_sigma_y(dim2),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate p0_sigma_y'
     stop
  endif

allocate(p0_sigma_ys(numsites),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate p0_sigma_ys'
     stop
  endif


allocate(y_error_t(nmeasuretotal),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate y_error_t array'
     stop
  endif

allocate(dumx(nmeasuremax),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate dumx'
     stop
  endif

allocate(dumy(nmeasuremax),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate dumy'
     stop
  endif

allocate(dum2(nmeasure),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate dum2'
     stop
  endif

allocate(C(nmeasuremax),stat=errstat)           
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate C'
     stop
  endif

allocate(C2(nmeasure),stat=errstat)           
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate C2'
     stop
  endif

allocate(T_current(nmeasuretotal,nmeasuretotal),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate T_current'
     stop
  endif

allocate(T_new(nmeasuretotal,nmeasuretotal),stat=errstat)         
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate T_new'
     stop
  endif

allocate(Tinv_current(nmeasuretotal,nmeasuretotal),stat=errstat)         
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate Tinv_current'
     stop
  endif

allocate(Tinv_new(nmeasuretotal,nmeasuretotal),stat=errstat)         
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate Tinv_new'
     stop
  endif

allocate(S_current(numsites,numsites),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate S_current'
     stop
  endif

allocate(S_new(numsites,numsites),stat=errstat)         
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate S_new'
     stop
  endif

allocate(Sinv_current(numsites,numsites),stat=errstat)         
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate Sinv_current'
     stop
  endif

allocate(Sinv_new(numsites,numsites),stat=errstat)         
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate Sinv_new'
     stop
  endif

allocate(Rinv_current(nmeasuremax,nmeasuremax),stat=errstat)        
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate Rinv_current'
     stop
  endif

allocate(autocorr_vec(nmeasuretotal),stat=errstat)          
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate autocorr_vec'
     stop
  endif

allocate(x_it(nIt,statesize),stat=errstat)           
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate x_it'
     stop
  endif

allocate(pdf_param1_it(nIt,statesize),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate pdf_param1_it'
     stop
  endif

allocate(pdf_param2_it(nIt,statesize),stat=errstat)           
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate pdf_param2_it'
     stop
  endif

allocate(sigma_y_it(nIt,dim2),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate sigma_y_it'
     stop
  endif

allocate(tau_it(nIT,dim3),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate tau_it'
     stop
  endif

allocate(nu_it(nIT,dim3),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate nu_it'
     stop
  endif

allocate(rho_it(nIT,dim3),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate rho_it'
     stop
  endif

allocate(sigma_ys_it(nIT,numsites),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate sigma_ys_it'
     stop
  endif

allocate(y_it(nIT,nmeasuremax),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate y_it'
     stop
  endif

allocate(reject_vector(statesize),stat=errstat)         
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate reject_vector'
     stop
  endif

allocate(reject_sigma_y_vector(dim2),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate reject_sigma_y_vector'
     stop
  endif

allocate(reject_sigma_ys_vector(numsites),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate reject_sigma_ys_vector'
     stop
  endif

allocate(reject_y_vector(nmeasuremax),stat=errstat)            
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate reject_y_vector'
     stop
  endif

allocate(accept_vector(statesize),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate accept_vector'
     stop
  endif

allocate(accept_sigma_y_vector(dim2),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate accept_sigma_y_vector'
     stop
  endif

allocate(accept_sigma_ys_vector(numsites),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate accept_sigma_ys_vector'
     stop
  endif

allocate(accept_y_vector(nmeasuremax),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate accept_y_vector'
     stop
  endif

allocate(datenumber(nmeasuretotal),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate datenumber'
     stop
  endif

allocate(timeindex_nonzero(nmeasure),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate timeindex_nonzero'
     stop
  endif

allocate(distribution(dim4),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate distribution'
     stop
  endif

allocate(n0_arr(nmeasuretotal,numsites),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate n0_arr'
     stop
  endif

allocate(n1_arr(nmeasuretotal,numsites),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate n1_arr'
     stop
  endif

allocate(C_arr(nmeasuretotal,numsites),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate C_arr'
     stop
  endif

allocate(dum3(nmeasuretotal,nmeasuretotal),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate dum3'
     stop
  endif

allocate(num_T_indices(dim2),stat=errstat)             
  if (errstat /= 0) then
     write (*,*) 'ERROR: could not allocate num_T_indices'
     stop
  endif


! **************** Read in variables from input netcdf file *********************

STATUS = NF90_INQ_VARID (inputID, 'x_ap', x_apID)
STATUS = NF90_INQ_VARID (inputID, 'stepsize', stepsizeID)
STATUS = NF90_INQ_VARID (inputID, 'stepsize_pdf_param1', stepsizepdfparam1ID)
STATUS = NF90_INQ_VARID (inputID, 'stepsize_pdf_param2', stepsizepdfparam2ID)
STATUS = NF90_INQ_VARID (inputID, 'stepsize_sigma_y', stepsizesigmayID)
STATUS = NF90_INQ_VARID (inputID, 'stepsize_tau', stepsizetauID)
STATUS = NF90_INQ_VARID (inputID, 'stepsize_nu', stepsizenuID)
STATUS = NF90_INQ_VARID (inputID, 'stepsize_rho', stepsizerhoID)
STATUS = NF90_INQ_VARID (inputID, 'stepsize_y', stepsizeyID)
STATUS = NF90_INQ_VARID (inputID, 'stepsize_sigma_ys', stepsizesigmaysID)
STATUS = NF90_INQ_VARID (inputID, 'pdf_param1', pdf_param1ID)
STATUS = NF90_INQ_VARID (inputID, 'pdf_param2', pdf_param2ID)
STATUS = NF90_INQ_VARID (inputID, 'sigma_y', sigmayID)
STATUS = NF90_INQ_VARID (inputID, 'sigma_ys', sigmaysID)
STATUS = NF90_INQ_VARID (inputID, 'tau', tauID)
STATUS = NF90_INQ_VARID (inputID, 'nu', nuID)
STATUS = NF90_INQ_VARID (inputID, 'rho', rhoID)
STATUS = NF90_INQ_VARID (inputID, 'pdf_hyperparam1', pdf_hyperparam1ID)
STATUS = NF90_INQ_VARID (inputID, 'pdf_hyperparam2', pdf_hyperparam2ID)
STATUS = NF90_INQ_VARID (inputID, 'sigma_y_hyperparam', sigmayhyperparamID)
STATUS = NF90_INQ_VARID (inputID, 'tau_hyperparam1', tauhyperparam1ID)
STATUS = NF90_INQ_VARID (inputID, 'tau_hyperparam2', tauhyperparam2ID)
STATUS = NF90_INQ_VARID (inputID, 'sigma_ys_hyperparam', sigmayshyperparamID)
STATUS = NF90_INQ_VARID (inputID, 'nu_hyperparam1', nuhyperparam1ID)
STATUS = NF90_INQ_VARID (inputID, 'rho_hyperparam1', rhohyperparam1ID)
STATUS = NF90_INQ_VARID (inputID, 'nu_hyperparam2', nuhyperparam2ID)
STATUS = NF90_INQ_VARID (inputID, 'rho_hyperparam2', rhohyperparam2ID)
STATUS = NF90_INQ_VARID (inputID, 'z', zID)
STATUS = NF90_INQ_VARID (inputID, 'nobs', nobsID)
STATUS = NF90_INQ_VARID (inputID, 'H', HID)
STATUS = NF90_INQ_VARID (inputID, 'D', DID)
STATUS = NF90_INQ_VARID (inputID, 'x_pdf', x_pdfID)
STATUS = NF90_INQ_VARID (inputID, 'pdf_param1_pdf', pdf_param1_pdfID)
STATUS = NF90_INQ_VARID (inputID, 'pdf_param2_pdf', pdf_param2_pdfID)
STATUS = NF90_INQ_VARID (inputID, 'tau_pdf', taupdfID)
STATUS = NF90_INQ_VARID (inputID, 'nu_pdf', nupdfID)
STATUS = NF90_INQ_VARID (inputID, 'rho_pdf', rhopdfID)
STATUS = NF90_INQ_VARID (inputID, 'sigma_y_pdf', sigmaypdfID)
STATUS = NF90_INQ_VARID (inputID, 'sigma_ys_pdf', sigmayspdfID)
STATUS = NF90_INQ_VARID (inputID, 'T_indices', TindicesID)
STATUS = NF90_INQ_VARID (inputID, 'deltatime', deltatimeID)
STATUS = NF90_INQ_VARID (inputID, 'distance', distanceID)
STATUS = NF90_INQ_VARID (inputID, 'datenumber', datenumberID)
STATUS = NF90_INQ_VARID (inputID, 'timeindex_nonzero', timeindexnonzeroID)
STATUS = NF90_INQ_VARID (inputID, 'kron_flag', kron_flagID)
STATUS = NF90_INQ_VARID (inputID, 'numthreads', numthreadsID)
STATUS = NF90_INQ_VARID (inputID, 'sitespresent', sitespresentID)
STATUS = NF90_INQ_VARID (inputID, 'distribution', distributionID)

STATUS = NF90_GET_VAR (inputID, x_apID, x_ap)
STATUS = NF90_GET_VAR (inputID, stepsizeID, stepsize)
STATUS = NF90_GET_VAR (inputID, stepsizepdfparam1ID, stepsize_pdf_param1)
STATUS = NF90_GET_VAR (inputID, stepsizepdfparam2ID, stepsize_pdf_param2)
STATUS = NF90_GET_VAR (inputID, stepsizesigmayID, stepsize_sigma_y)
STATUS = NF90_GET_VAR (inputID, stepsizetauID, stepsize_tau)
STATUS = NF90_GET_VAR (inputID, stepsizenuID, stepsize_nu)
STATUS = NF90_GET_VAR (inputID, stepsizerhoID, stepsize_rho)
STATUS = NF90_GET_VAR (inputID, stepsizeyID, stepsize_y)
STATUS = NF90_GET_VAR (inputID, stepsizesigmaysID, stepsize_sigma_ys)
STATUS = NF90_GET_VAR (inputID, pdf_param1ID, pdf_param1)
STATUS = NF90_GET_VAR (inputID, pdf_param2ID, pdf_param2)
STATUS = NF90_GET_VAR (inputID, sigmayID, sigma_y)
STATUS = NF90_GET_VAR (inputID, sigmaysID, sigma_ys)
STATUS = NF90_GET_VAR (inputID, tauID, tau)
STATUS = NF90_GET_VAR (inputID, nuID, nu)
STATUS = NF90_GET_VAR (inputID, rhoID, rho)
STATUS = NF90_GET_VAR (inputID, pdf_hyperparam1ID, pdf_hyperparam1)
STATUS = NF90_GET_VAR (inputID, pdf_hyperparam2ID, pdf_hyperparam2)
STATUS = NF90_GET_VAR (inputID, sigmayhyperparamID, sigma_y_hyperparam)
STATUS = NF90_GET_VAR (inputID, tauhyperparam1ID, tau_hyperparam1)
STATUS = NF90_GET_VAR (inputID, tauhyperparam2ID, tau_hyperparam2)
STATUS = NF90_GET_VAR (inputID, nuhyperparam1ID, nu_hyperparam1)
STATUS = NF90_GET_VAR (inputID, rhohyperparam1ID, rho_hyperparam1)
STATUS = NF90_GET_VAR (inputID, nuhyperparam2ID, nu_hyperparam2)
STATUS = NF90_GET_VAR (inputID, rhohyperparam2ID, rho_hyperparam2)
STATUS = NF90_GET_VAR (inputID, sigmayshyperparamID, sigma_ys_hyperparam)
STATUS = NF90_GET_VAR (inputID, zID, z)
STATUS = NF90_GET_VAR (inputID, nobsID, nobs)
STATUS = NF90_GET_VAR (inputID, HID, H)
STATUS = NF90_GET_VAR (inputID, DID, D)
STATUS = NF90_GET_VAR (inputID, x_pdfID, x_pdf)
STATUS = NF90_GET_VAR (inputID, pdf_param1_pdfID, pdf_param1_pdf)
STATUS = NF90_GET_VAR (inputID, pdf_param2_pdfID, pdf_param2_pdf)
STATUS = NF90_GET_VAR (inputID, taupdfID, tau_pdf)
STATUS = NF90_GET_VAR (inputID, nupdfID, nu_pdf)
STATUS = NF90_GET_VAR (inputID, rhopdfID, rho_pdf)
STATUS = NF90_GET_VAR (inputID, sigmaypdfID, sigma_y_pdf)
STATUS = NF90_GET_VAR (inputID, sigmayspdfID, sigma_ys_pdf)
STATUS = NF90_GET_VAR (inputID, TindicesID, T_indices)
STATUS = NF90_GET_VAR (inputID, deltatimeID, deltatime)
STATUS = NF90_GET_VAR (inputID, distanceID, distance)
STATUS = NF90_GET_VAR (inputID, datenumberID, datenumber)
STATUS = NF90_GET_VAR (inputID, timeindexnonzeroID, timeindex_nonzero)
STATUS = NF90_GET_VAR (inputID, kron_flagID, kron_flag)
STATUS = NF90_GET_VAR (inputID, numthreadsID, numthreads)
STATUS = NF90_GET_VAR (inputID, sitespresentID, sitespresent)
STATUS = NF90_GET_VAR (inputID, distributionID, distribution)

STATUS = nf90_close(inputID)

!****************** Hierarchical Bayesian inversion********************* 
! Set number of threads for multi-threaded routines
 call OMP_SET_NUM_THREADS(numthreads)
 call MKL_SET_NUM_THREADS(numthreads)

 ! Set first value of each parameter from prior values
 x=x_ap  
 pdf_param1_current = pdf_param1
 pdf_param2_current = pdf_param2 
 sigma_y_current = sigma_y
 sigma_ys_current = sigma_ys
 tau_current = tau
 nu_current = nu
 rho_current = rho

 ! Create autoregressive AR(1) temporal autocorrelation matrix

 do ii=1,dim2
 	y_error_t(T_indices(:,ii)) = sigma_y_current(ii)
 enddo

!$OMP PARALLEL DO private(ii, autocorr_vec) shared(T_current)
 do ii=1,nmeasuretotal
 	autocorr_vec = y_error_t(ii)*y_error_t*exp(-1*(deltatime(:,ii)/tau_current))
	T_current(ii,:) = autocorr_vec
	T_current(:,ii) = autocorr_vec      
 enddo
!$OMP END PARALLEL DO

 ! Compute inverse of T

 Tinv_current = Ainv(T_current)

! find number of elements corresponding to each column of T_indices in order to scale determinant

 do ii = 1,dim2
	num_T_indices(ii) =maxloc(abs(T_indices(:,ii)), 1)
 enddo	

 ! Create Matern covariance for spatial correlation

 nu_double = nu_current
 rho_double = rho_current
 distance_double = distance
 nb = int(nu_double)
 alpha = nu_double - nb
 nb = nb+1
 ize = 2
 
 arg = sqrt(2*nu_double)*distance_double/rho_double

 call gamma(nu_double, ga_double)
 ga = ga_double

!$OMP PARALLEL DO private(ii,jj, arg_dum, k_arg_dum) shared(S_current) 
	do ii = 1,numsites 
		do jj = 1,numsites
			arg_dum = arg(ii,jj)
			if (arg_dum==0) then
				S_current(ii,jj) = 1*sigma_ys_current(ii)**2
			else
 				call rkbesl(arg_dum, alpha, nb, ize, k_arg_dum, ncalc)
	 			k_arg_dum = k_arg_dum*exp(-1*arg(ii,jj))
				S_current(ii,jj) = sigma_ys_current(ii)*sigma_ys_current(jj)*(1/(2**(nu_current-1)*ga)*(sqrt(2*nu_current)*distance(ii,jj)/rho_current)**nu_current*k_arg_dum)
			endif
		enddo
	enddo
!$OMP END PARALLEL DO

 ! Compute inverse of S
 Sinv_current = Ainv(S_current)

 ! Calulate R using the Kronecker product
 if (kron_flag==0) then
    call kron(Rinv_current, Sinv_current, Tinv_current)
 elseif (kron_flag==1) then
    call kron(Rinv_current, Tinv_current, Sinv_current)
 endif

 ! Constants for Intel fortran matrix multiplication routines
 aa = 1
 bb = 0

 ! Compute (log sqrt) determinant of R using properties of the Kronecker product
 detval_S_current = det(S_current)*nmeasuretotal
 detval_T_current = det(T_current)*numsites
 detval_current = detval_S_current + detval_T_current

 ! Calclulate the inverse of the measurement uncertainty nugget term (assumed diagonal)
 Dinv = Ainv(D)

 ! Intial value of model simulated mole fractions
 y_current = matmul(H,x)

 ! Subsample at measurement times
 y_small = y_current(timeindex_nonzero)

 ! Initial value of the difference between observations and model mole fractions at observation times
 m0 = z - y_small

 ! Compute (z-y)^T*Dinv*(z-y)
 call sgemm ('N', 'N', nmeasure, 1, nmeasure, aa, Dinv, nmeasure, m0, nmeasure, bb, C2, nmeasure)
 
 m0T = dot_product(m0,C2)

 ! Compute difference between sampled value of model mole fractions and model derived ones based on value of x
 n0 = matmul(H,x) - y_current   

 ! Compute (Hx-y)^T*Rinv*(Hx-y) by reshaping into an array and using properties of the Kronecker product
 n0_arr = reshape(n0, (/nmeasuretotal, numsites/))
 C_arr = matmul(matmul(Tinv_current,n0_arr),Sinv_current)
 C = reshape(C_arr,(/nmeasuremax/))
									
 n0T = dot_product(n0,C)


 ! Compute P0 for starting state   
 do xi=1,statesize	
   	call pdf_calc(x(xi), pdf_param1(xi), pdf_param2(xi), x_pdf(xi), p1)
  	p0_x(xi) = p1

   	call pdf_calc(pdf_param1_current(xi), pdf_param1(xi), pdf_hyperparam1(xi), pdf_param1_pdf(xi), p1)
   	p0_pdf_param1(xi) = p1

   	call pdf_calc(pdf_param2_current(xi), pdf_param2(xi), pdf_hyperparam2(xi), pdf_param2_pdf(xi), p1)
   	p0_pdf_param2(xi) = p1
 enddo 

 do yi = 1,dim2
 	call pdf_calc(sigma_y_current(yi), sigma_y(yi), sigma_y_hyperparam(yi), sigma_y_pdf(yi), p1)
   	p0_sigma_y(yi) = p1
 enddo

 do yi = 1,numsites
 	call pdf_calc(sigma_ys_current(yi), sigma_ys(yi), sigma_ys_hyperparam(yi), sigma_ys_pdf(yi), p1)
   	p0_sigma_ys(yi) = p1
 enddo

 call pdf_calc(tau_current, tau_hyperparam1, tau_hyperparam2, tau_pdf, p1)
 p0_tau = p1

 call pdf_calc(nu_current, nu_hyperparam1, nu_hyperparam2, nu_pdf, p1)
 p0_nu = p1

 call pdf_calc(rho_current, rho_hyperparam1, rho_hyperparam2, rho_pdf, p1)
 p0_rho = p1

! Initialize arrays that will store values and acceptance ratios
 x_it(1,:)=x
 pdf_param1_it(1,:) = pdf_param1
 pdf_param2_it(1,:) = pdf_param2
 sigma_y_it(1,:) = sigma_y
 sigma_ys_it(1,:) = sigma_ys
 tau_it(1,1) = tau
 nu_it(1,1) = nu
 rho_it(1,1) = rho
 y_it(1,:) = y

 reject=0
 reject_vector(:)=0
 reject_sigma_y=0
 reject_sigma_y_vector(:)=0
 reject_sigma_ys = 0
 reject_sigma_ys_vector(:)=0
 reject_tau=0
 reject_nu=0
 reject_y = 0
 reject_y_vector(:) = 0

 call init_random_seed()              ! Ensure random number generation starts from new point each time program is run
                                      ! Random seed only needs to be called once in a program.  
				      ! It is used to generate the first random number. 
				      ! Any random numbers generated after the first will use the previous random number as its seed.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MCMC MAIN LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 do It=1,nIt+burn_in
			

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!X LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
	do xi=1,statesize
		! Calculate new values of x, pdf_param1, pdf_param2 using random number generator
    		dpdf_param1 = random_normal()
       		dpdf_param1 = dpdf_param1*stepsize_pdf_param1(xi)
       		pdf_param1_new(xi) = pdf_param1_current(xi) + dpdf_param1

       		dpdf_param2 = random_normal()
       		dpdf_param2 = dpdf_param2*stepsize_pdf_param2(xi)
       		pdf_param2_new(xi) = pdf_param2_current(xi) + dpdf_param2

       		dx = random_normal()  		  
       		dx = dx*stepsize(xi)              
       		dy=H(:,xi)*dx

		! Compute P1 for each parameter
    		call pdf_calc(pdf_param1_new(xi), pdf_param1(xi), pdf_hyperparam1(xi), pdf_param1_pdf(xi), p1)
    		p1_pdf_param1 = p1

    		call pdf_calc(pdf_param2_new(xi), pdf_param2(xi), pdf_hyperparam2(xi), pdf_param2_pdf(xi), p1)
    		p1_pdf_param2 = p1

    		call pdf_calc((x(xi)+dx), pdf_param1_new(xi), pdf_param2_new(xi), x_pdf(xi), p1)  
    		p1_x = p1
	
		! Compute n1T = (y-Hx)^T*Rinv*(y-HX) for new value of x
      		n1=(n0+dy)

		n1_arr = reshape(n1, (/nmeasuretotal, numsites/))
		call sgemm ('N', 'N', nmeasuretotal, numsites, nmeasuretotal, aa, Tinv_current, nmeasuretotal, n1_arr, nmeasuretotal, bb, dum3, nmeasuretotal)
		call sgemm ('N', 'N', nmeasuretotal, numsites, numsites, aa, dum3, nmeasuretotal, Sinv_current, numsites, bb, C_arr, nmeasuretotal)
		C = reshape(C_arr,(/nmeasuremax/))

		n1T = dot_product(n1,C)
	
		! Compute P1/P0
     		pT=alog(p1_x*p1_pdf_param1*p1_pdf_param2/(p0_x(xi)*p0_pdf_param1(xi)*p0_pdf_param2(xi)))-0.5*(n1T - n0T)
      
      		call random_number(randomu)        	 ! Generates uniformly distributed random number
      
      		if(alog(randomu) .le. pT) then      
        		!;ACCEPT
         		p0_x(xi)=p1_x
	 		p0_pdf_param1(xi) = p1_pdf_param1
	 		p0_pdf_param2(xi) = p1_pdf_param2
         	
         		x(xi)=x(xi) + dx         		
         		pdf_param1_current(xi) = pdf_param1_new(xi)
         		pdf_param2_current(xi) = pdf_param2_new(xi)

         		n0=n1
         		n0T=n1T

			n0_arr = n1_arr
      		else 
         		!;REJECT					
         		if(it .gt. burn_in) reject=reject+ 1		
         		if(it .gt. burn_in) reject_vector(xi)=reject_vector(xi)+ 1
      		endif
	enddo                       

    	if(it .gt. burn_in) then
        	if(mod(It+1,100) .eq. 0) then           ! Print acceptance ratio for every 1000th iteration to see if everything looks ok
           		call system_clock(count_1,count_rate,count_max)
           		t1=count_1/count_rate
           		write(*,*) t1-t0, It - burn_in, 'x acceptance: ', 1. - real(reject)/real(It-burn_in)/real(statesize)   
        	endif 
        	x_it(it-burn_in,:)=x        		! Record x value for It > burn_in
        	pdf_param1_it(it-burn_in,:) = pdf_param1_current
        	pdf_param2_it(it-burn_in,:) = pdf_param2_current  

    	endif


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SIGMA_Y LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do yi=1,dim2

		sigma_y_new = sigma_y_current
	
		! Generate new value of sigma_y
    		dsigma_y = random_normal()
    		dsigma_y = dsigma_y*stepsize_sigma_y(yi)
    		sigma_y_new(yi) = sigma_y_current(yi) + dsigma_y

		! Calculate P1 for new value of sigma_y
   		call pdf_calc(sigma_y_new(yi), sigma_y(yi), sigma_y_hyperparam(yi), sigma_y_pdf(yi), p1)
    		p1_sigma_y = p1

		! Calculate new value of inverse by scaling old value 
		Tinv_new = Tinv_current
		Tinv_new(T_indices(:,yi),:) = 1/(sigma_y_new(yi)/sigma_y_current(yi))*Tinv_current(T_indices(:,yi),:)		
		Tinv_new(:,T_indices(:,yi)) = 1/(sigma_y_new(yi)/sigma_y_current(yi))*Tinv_new(:,T_indices(:,yi))

		! calculate new value of n1T = (y-Hx)^T * Rinv * (y-Hx) 	 

		call sgemm ('N', 'N', nmeasuretotal, numsites, nmeasuretotal, aa, Tinv_new, nmeasuretotal, n0_arr, nmeasuretotal, bb, dum3, nmeasuretotal)
		call sgemm ('N', 'N', nmeasuretotal, numsites, numsites, aa, dum3, nmeasuretotal, Sinv_current, numsites, bb, C_arr, nmeasuretotal)
		C = reshape(C_arr,(/nmeasuremax/))

		n1T = dot_product(n0,C)

		! caluclate new determinant by scaling
		detval_T_new = (num_T_indices(yi)*alog(sigma_y_new(yi)/sigma_y_current(yi)))*numsites + detval_T_current

		detval_new = detval_S_current + detval_T_new	
		
		! Compute P1/P0 
		pT=alog(p1_sigma_y/p0_sigma_y(yi)) - detval_new + detval_current - 0.5*(n1T - n0T)            

     		call random_number(randomu)        	 ! Generates uniformly distributed random number
      
      		if(alog(randomu) .le. pT) then      
        		 !;ACCEPT
         		p0_sigma_y(yi)=p1_sigma_y	
        		sigma_y_current(yi) = sigma_y_new(yi)
			Tinv_current = Tinv_new
			detval_current = detval_new
			detval_T_current = detval_T_new
         		n0T=n1T
		else
        		 !;REJECT					
         		if(it .gt. burn_in) reject_sigma_y=reject_sigma_y+ 1		
         		if(it .gt. burn_in) reject_sigma_y_vector(yi)=reject_sigma_y_vector(yi)+ 1
      		endif

	enddo 
         

    	if(it .gt. burn_in) then
       		if(mod(It+1,100) .eq. 0) then              ! Print sigma_y acceptance ratio for every 1000th iteration
           	call system_clock(count_1,count_rate,count_max)
           	t1=count_1/count_rate
           	write(*,*) t1-t0, It - burn_in, 'sigma_y acceptance: ', 1. - real(reject_sigma_y)/real(It-burn_in)/real(dim2)   
       		endif 
        	sigma_y_it(it-burn_in,:)=sigma_y_current   ! Record x value for It > burn_in
       
    	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SIGMA_YS LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	do yi=1,numsites

		sigma_ys_new = sigma_ys_current
 
		! Generate new value of sigma_ys
    		dsigma_ys = random_normal()
    		dsigma_ys = dsigma_ys*stepsize_sigma_ys(yi)
    		sigma_ys_new(yi) = sigma_ys_current(yi) + dsigma_ys

		! Calculate P1 for new value of sigma_y
   		call pdf_calc(sigma_ys_new(yi), sigma_ys(yi), sigma_ys_hyperparam(yi), sigma_ys_pdf(yi), p1)
    		p1_sigma_ys = p1

		! Calculate new inverse by scaling based on value of sigma_ys
		Sinv_new = Sinv_current
		Sinv_new(yi,:) = 1/(sigma_ys_new(yi)/sigma_ys_current(yi))*Sinv_current(yi,:)
		Sinv_new(:,yi) = 1/(sigma_ys_new(yi)/sigma_ys_current(yi))*Sinv_new(:,yi)

		! calculate new value of n1T = (y-Hx)^T * Rinv * (y-Hx) 	 

		call sgemm ('N', 'N', nmeasuretotal, numsites, nmeasuretotal, aa, Tinv_current, nmeasuretotal, n0_arr, nmeasuretotal, bb, dum3, nmeasuretotal)
		call sgemm ('N', 'N', nmeasuretotal, numsites, numsites, aa, dum3, nmeasuretotal, Sinv_new, numsites, bb, C_arr, nmeasuretotal)
 		C = reshape(C_arr,(/nmeasuremax/))

		n1T = dot_product(n0,C)

		! calculate new determinant based on scaling

		detval_S_new = (alog(sigma_ys_new(yi)/sigma_ys_current(yi)))*nmeasuretotal + detval_S_current

		detval_new = detval_T_current + detval_S_new 

		! Compute P1/P0 

		pT=alog(p1_sigma_ys/p0_sigma_ys(yi)) - detval_new + detval_current - 0.5*(n1T - n0T)

     		call random_number(randomu)        	 ! Generates uniformly distributed random number
      
      		if(alog(randomu) .le. pT) then      
        		 !;ACCEPT
         		p0_sigma_ys(yi) = p1_sigma_ys	
        		sigma_ys_current(yi) = sigma_ys_new(yi)
			Sinv_current = Sinv_new
			S_current = S_new
			detval_current = detval_new
			detval_S_current = detval_S_new
         		n0T=n1T
		else
        		 !;REJECT					
         		if(it .gt. burn_in) reject_sigma_ys=reject_sigma_ys+ 1		
         		if(it .gt. burn_in) reject_sigma_ys_vector(yi)=reject_sigma_ys_vector(yi)+ 1	
      		endif
        enddo  
		

    	if(it .gt. burn_in) then
       		if(mod(It+1,100) .eq. 0) then              ! Print sigma_ys acceptance ratio for every 1000th iteration
           	call system_clock(count_1,count_rate,count_max)
           	t1=count_1/count_rate
           	write(*,*) t1-t0, It - burn_in, 'sigma_ys acceptance: ', 1. - real(reject_sigma_ys)/real(It-burn_in)/real(numsites) 
       		endif 
        	sigma_ys_it(it-burn_in,:)=sigma_ys_current   ! Record x value for It > burn_in
       
    	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TAU LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Generate new value of tau
	dtau = random_normal()
    	dtau = dtau*stepsize_tau
    	tau_new = tau_current + dtau 

 	if(tau_new .le. 0) then
	! Skip if it generates a tau less than or equal to 0 or matrix will blow-up
		if(it .gt. burn_in) reject_tau=reject_tau+ 1	
 	else

		! Compute P1 for new value of tau
    		call pdf_calc(tau_new, tau_hyperparam1, tau_hyperparam2, tau_pdf, p1)
    		p1_tau = p1

		! create new T matrix
    		do ii=1,dim2
    			y_error_t(T_indices(:,ii)) = sigma_y_current(ii)
    		enddo

!$OMP PARALLEL DO private(ii, autocorr_vec) shared(T_new)
		 do ii=1,nmeasuretotal
		 	autocorr_vec = y_error_t(ii)*y_error_t*exp(-1*(deltatime(:,ii)/tau_new))
			T_new(ii,:) = autocorr_vec
			T_new(:,ii) = autocorr_vec   
 		enddo
!$OMP END PARALLEL DO

	 	! compute new inverse
		Tinv_new = Ainv(T_new)	

		! calculate new value of n1T = (y-Hx)^T * Rinv * (y-Hx) 

		call sgemm ('N', 'N', nmeasuretotal, numsites, nmeasuretotal, aa, Tinv_new, nmeasuretotal, n0_arr, nmeasuretotal, bb, dum3, nmeasuretotal)
		call sgemm ('N', 'N', nmeasuretotal, numsites, numsites, aa, dum3, nmeasuretotal, Sinv_current, numsites, bb, C_arr, nmeasuretotal)
 		C = reshape(C_arr,(/nmeasuremax/))

		n1T = dot_product(n0,C)

		! calculate new determinant
		
		detval_T_new = det(T_new)*numsites
		detval_new = detval_S_current + detval_T_new

		! compute P1/P0 
   		pT=alog(p1_tau/p0_tau) - detval_new + detval_current - 0.5*(n1T - n0T)

    		call random_number(randomu)        	 ! Generates uniformly distributed random number
 
     		if(alog(randomu) .le. pT) then      
         		!;ACCEPT
         		p0_tau=p1_tau	
         		tau_current = tau_new
			Tinv_current = Tinv_new
			detval_T_current = detval_T_new
			detval_current = detval_new
         		n0T=n1T
      		else 
        		!;REJECT					
         		if(it .gt. burn_in) reject_tau=reject_tau+ 1		
      		endif
 	endif

	if(it .gt. burn_in) then
        	if(mod(It+1,100) .eq. 0) then              ! Print tau acceptance ratio for every 1000th iteration
           		call system_clock(count_1,count_rate,count_max)
           		t1=count_1/count_rate
           		write(*,*) t1-t0, It - burn_in, 'tau acceptance: ', 1. - real(reject_tau)/real(It-burn_in)  
        	endif 
        	tau_it(it-burn_in,:)=tau_current        ! Record x value for It > burn_in
       
    	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NU and RHO LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Generate new value of nu and rho
	dnu = random_normal()
    	dnu = dnu*stepsize_nu
    	nu_new = nu_current + dnu 

	drho = random_normal()
    	drho = drho*stepsize_rho
    	rho_new = rho_current + drho

 	if(nu_new .le. 0 .OR. rho_new .le. 0 ) then
	! Skip if it generates a nu or rho less than or equal to 0 or matrix will blow-up
		if(it .gt. burn_in) reject_nu=reject_nu+ 1	
 	else

		! Compute P1 for new value of rho and nu
    		call pdf_calc(nu_new, nu_hyperparam1, nu_hyperparam2, nu_pdf, p1)
    		p1_nu = p1

    		call pdf_calc(rho_new, rho_hyperparam1, rho_hyperparam2, rho_pdf, p1)
    		p1_rho = p1

		! Create new S matrix

 		nu_double = nu_new
 		rho_double = rho_new
 		nb = int(nu_double)
 		alpha = nu_double - nb
 		nb = nb+1
 
 		arg = sqrt(2*nu_double)*distance_double/rho_double

 		ize = 2

 		call gamma(nu_double, ga_double)
		ga = ga_double

!$OMP PARALLEL DO private(ii,jj, arg_dum, k_arg_dum) shared(S_new) 
		do ii = 1,numsites 
			do jj = 1,numsites
				arg_dum = arg(ii,jj)
				if (arg_dum==0) then
					S_new(ii,jj) = 1*sigma_ys_current(ii)**2
				else
	 				call rkbesl(arg_dum, alpha, nb, ize, k_arg_dum, ncalc)
		 			k_arg_dum = k_arg_dum*exp(-1*arg(ii,jj))
					S_new(ii,jj) = sigma_ys_current(ii)*sigma_ys_current(jj)*(1/(2**(nu_new-1)*ga)*(sqrt(2*nu_new)*distance(ii,jj)/rho_new)**nu_new*k_arg_dum)
				endif
			enddo
		enddo
!$OMP END PARALLEL DO

		! compute new inverse

		Sinv_new = Ainv(S_new)	

		! calculate new value of n1T = (y-Hx)^T * Rinv * (y-Hx) 

		call sgemm ('N', 'N', nmeasuretotal, numsites, nmeasuretotal, aa, Tinv_current, nmeasuretotal, n0_arr, nmeasuretotal, bb, dum3, nmeasuretotal)
		call sgemm ('N', 'N', nmeasuretotal, numsites, numsites, aa, dum3, nmeasuretotal, Sinv_new, numsites, bb, C_arr, nmeasuretotal)
 		C = reshape(C_arr,(/nmeasuremax/))

		n1T = dot_product(n0,C)

		! calculate new determinant

   		detval_S_new = det(S_new)*nmeasuretotal
   		detval_new = detval_T_current + detval_S_new	

		! compute P1/P0 
   		pT=alog(p1_nu*p1_rho/p0_nu/p0_rho) - detval_new + detval_current - 0.5*(n1T - n0T)

    		call random_number(randomu)        	 ! Generates uniformly distributed random number
 
     		if(alog(randomu) .le. pT) then      
         		!;ACCEPT
         		p0_nu=p1_nu	
         		p0_rho=p1_rho
         		nu_current = nu_new
         		rho_current = rho_new
			Sinv_current = Sinv_new
			S_current = S_new
			detval_S_current = detval_S_new
			detval_current = detval_new
         		n0T=n1T
      		else 
        		!;REJECT					
         		if(it .gt. burn_in) reject_nu=reject_nu+ 1		
      		endif
 	endif

	if(it .gt. burn_in) then
        	if(mod(It+1,100) .eq. 0) then              ! Print tau acceptance ratio for every 1000th iteration
           		call system_clock(count_1,count_rate,count_max)
           		t1=count_1/count_rate
           		write(*,*) t1-t0, It - burn_in, 'nu_rho acceptance: ', 1. - real(reject_nu)/real(It-burn_in)  
        	endif 
        	nu_it(it-burn_in,:)=nu_current        ! Record x value for It > burn_in
		rho_it(it-burn_in,:)=rho_current
       
    	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Y LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! create full covariance matrix R based on current state of S and T

                if (kron_flag==0) then
                    call kron(Rinv_current, Sinv_current, Tinv_current)
                elseif (kron_flag==1) then
                    call kron(Rinv_current, Tinv_current, Sinv_current)
                endif

	do yi=1,nmeasuremax

		y_new = y_current
	
		! Generate new value of y
   		dymod = random_normal()
    		dymod = dymod*stepsize_y(yi)
    		y_new(yi) = y_current(yi) + dymod
	
		! subsample at observation times
		y_small = y_new(timeindex_nonzero)

		! calculate (z-y)^T*Dinv*(z-y)
		m1 = z - y_small

!$OMP PARALLEL DO private(ii) shared(dum2) 
		do ii=1,nmeasure		
  			dum2(ii) = m1(ii)*Dinv(ii,ii)*m1(ii)
		end do
!$OMP END PARALLEL DO
		m1T = sum(dum2)
		
		! compute (Hx-y)^T*Rinv*(Hx-y) based on current value of y
		n1 = n0
		n1(yi) = n0(yi) - dymod

		dumx = n0(yi)*Rinv_current(yi,:)
		dumy = n1(yi)*Rinv_current(yi,:)


		dum  = dot_product(n0,dumx)  
		dumb = dot_product(n1,dumy)
		dumc = dot_product(n0,Rinv_current(yi,:))
		dumd = dot_product(n1,Rinv_current(yi,:))

		n1T = n0T - dum + dumb - n0(yi)*dumc + n1(yi)*dumd - n1(yi)*n1(yi)*Rinv_current(yi,yi) + n0(yi)*n0(yi)*Rinv_current(yi,yi)


		! Compute P1/P0 
		pT= -0.5*(m1T - m0T) -0.5*(n1T - n0T) 


     		call random_number(randomu)        	 ! Generates uniformly distributed random number
      
      		if(alog(randomu) .le. pT) then      
        		 !;ACCEPT
        		y_current(yi) = y_new(yi)
         		n0=n1
        		n0T=n1T
			m0T = m1T
      		else 
        		 !;REJECT					
         		if(it .gt. burn_in) reject_y=reject_y+ 1		
         		if(it .gt. burn_in) reject_y_vector(yi)=reject_y_vector(yi)+ 1
      		endif
	enddo 
         

    	if(it .gt. burn_in) then
       		if(mod(It+1,100) .eq. 0) then              ! Print y acceptance ratio for every 1000th iteration
           	call system_clock(count_1,count_rate,count_max)
           	t1=count_1/count_rate
           	write(*,*) t1-t0, It - burn_in, 'y acceptance: ', 1. - real(reject_y)/real(It-burn_in)/real(nmeasuremax)   
       		endif 
        	y_it(it-burn_in,:)=y_current   ! Record x value for It > burn_in
       
    	endif
    
 end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END MCMC MAIN LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute acceptance ratios
 accept_vector = 1. - real(reject_vector)/real(nIt)

 accept_sigma_y_vector = 1. - real(reject_sigma_y_vector)/real(nIt)

 accept_sigma_ys_vector = 1. - real(reject_sigma_ys_vector)/real(nIt)

 accept_tau = 1. - real(reject_tau)/real(nIt)

 accept_nu = 1. - real(reject_nu)/real(nIt)

 accept_y_vector = 1. - real(reject_y_vector)/real(nIt)

 write(*,*) 'Finished HB_MCMC'

  ! ********************************** Write output to netcdf **************************************

  STATUS = NF90_DEF_DIM (outputID, 'nit', nIt, nItDID)
  STATUS = NF90_DEF_DIM (outputID, 'statesize', statesize, statesizeDID)
  STATUS = NF90_DEF_DIM (outputID, 'nmeasure', nMeasure, nMeasureDID)
  STATUS = NF90_DEF_DIM (outputID, 'dim1', dim1, dim1DID)
  STATUS = NF90_DEF_DIM (outputID, 'dim2', dim2, dim2DID)
  STATUS = NF90_DEF_DIM (outputID, 'dim3', dim3, dim3DID)
  STATUS = NF90_DEF_DIM (outputID, 'nmeasuretotal', nmeasuretotal, nmeasuretotalDID)
  STATUS = NF90_DEF_DIM (outputID, 'nmeasuremax', nmeasuremax, nmeasuremaxDID)
  STATUS = NF90_DEF_DIM (outputID, 'numsites', numsites, numsitesDID)
  STATUS = NF90_DEF_DIM (outputID, 'dim4', dim4, dim4DID)


  STATUS = NF90_DEF_VAR (outputID, 'x_it', NF90_DOUBLE, (/ nItDID, statesizeDID/), x_itID)
  STATUS = NF90_DEF_VAR (outputID, 'pdf_param1_it', NF90_DOUBLE, (/ nItDID, statesizeDID/), pdf_param1_itID)
  STATUS = NF90_DEF_VAR (outputID, 'pdf_param2_it', NF90_DOUBLE, (/ nItDID, statesizeDID/), pdf_param2_itID)
  STATUS = NF90_DEF_VAR (outputID, 'sigma_y_it', NF90_DOUBLE, (/ nItDID, dim2DID/), sigma_y_itID)
  STATUS = NF90_DEF_VAR (outputID, 'sigma_ys_it', NF90_DOUBLE, (/ nItDID, numsitesDID/), sigmays_itID)
  STATUS = NF90_DEF_VAR (outputID, 'tau_it', NF90_DOUBLE, nItDID, tau_itID)
  STATUS = NF90_DEF_VAR (outputID, 'nu_it', NF90_DOUBLE, nItDID, nu_itID)
  STATUS = NF90_DEF_VAR (outputID, 'rho_it', NF90_DOUBLE, nItDID, rho_itID)
  STATUS = NF90_DEF_VAR (outputID, 'y_it', NF90_DOUBLE, (/ nItDID, nmeasuremaxDID/), y_itID)
  STATUS = NF90_DEF_VAR (outputID, 'acceptance', NF90_DOUBLE, statesizeDID, acceptID)
  STATUS = NF90_DEF_VAR (outputID, 'acceptance_sigmay', NF90_DOUBLE, dim2DID, acceptsigmayID)
  STATUS = NF90_DEF_VAR (outputID, 'acceptance_sigmays', NF90_DOUBLE, numsitesDID , acceptsigmaysID)
  STATUS = NF90_DEF_VAR (outputID, 'acceptance_tau', NF90_DOUBLE, dim3DID , accepttauID)
  STATUS = NF90_DEF_VAR (outputID, 'acceptance_nu', NF90_DOUBLE, dim3DID , acceptnuID)
  STATUS = NF90_DEF_VAR (outputID, 'acceptance_y', NF90_DOUBLE, nmeasuremaxDID, acceptyID)
  STATUS = NF90_DEF_VAR (outputID, 'H', NF90_DOUBLE, (/ nMeasuremaxDID, statesizeDID/), HoID) ! include for convenience
  STATUS = NF90_DEF_VAR (outputID, 'z', NF90_DOUBLE, nMeasureDID, zoID) ! include for convenience
  STATUS = NF90_DEF_VAR (outputID, 'datenumber', NF90_DOUBLE, nmeasuretotalDID, datenumberoID) ! include for convenience
  STATUS = NF90_DEF_VAR (outputID, 'timeindex_nonzero', NF90_DOUBLE, nMeasureDID, timeindexnonzerooID) ! include for convenience
  STATUS = NF90_DEF_VAR (outputID, 'x_ap', NF90_DOUBLE, statesizeDID, xapoID) ! include for convenience
  STATUS = NF90_DEF_VAR (outputID, 'nobs', NF90_DOUBLE, numsitesDID, nobsoID) ! include for convenience
  STATUS = NF90_DEF_VAR (outputID, 'sitespresent', NF90_DOUBLE, numsitesDID, sitespresentoID) ! include for convenience
  STATUS = NF90_DEF_VAR (outputID, 'deltatime', NF90_DOUBLE, (/ nmeasuretotalDID, nmeasuretotalDID/), deltatimeoID) ! include for convenience
  STATUS = NF90_DEF_VAR (outputID, 'distance', NF90_DOUBLE, (/ numsitesDID, numsitesDID/), distanceoID) ! include for convenience
  STATUS = NF90_DEF_VAR (outputID, 'T_indices', NF90_DOUBLE, (/ dim1DID, dim2DID/), TindicesoID) ! include for convenience
  STATUS = NF90_DEF_VAR (outputID, 'distribution', NF90_DOUBLE, dim4DID, distributionoID) ! include for convenience
  STATUS = NF90_DEF_VAR (outputID, 'D', NF90_DOUBLE, (/ nMeasureDID, nMeasureDID/), DoID) ! include for convenience
  STATUS = NF90_ENDDEF (outputID)

  STATUS = NF90_PUT_VAR (outputID, x_itID, x_it)
  STATUS = NF90_PUT_VAR (outputID, pdf_param1_itID, pdf_param1_it)
  STATUS = NF90_PUT_VAR (outputID, pdf_param2_itID, pdf_param2_it)
  STATUS = NF90_PUT_VAR (outputID, sigma_y_itID, sigma_y_it)
  STATUS = NF90_PUT_VAR (outputID, sigmays_itID, sigma_ys_it)
  STATUS = NF90_PUT_VAR (outputID, tau_itID, tau_it)
  STATUS = NF90_PUT_VAR (outputID, nu_itID, nu_it)
  STATUS = NF90_PUT_VAR (outputID, rho_itID, rho_it)
  STATUS = NF90_PUT_VAR (outputID, y_itID, y_it)
  STATUS = NF90_PUT_VAR (outputID, acceptID, accept_vector)
  STATUS = NF90_PUT_VAR (outputID, acceptsigmayID, accept_sigma_y_vector)
  STATUS = NF90_PUT_VAR (outputID, acceptsigmaysID, accept_sigma_ys_vector)
  STATUS = NF90_PUT_VAR (outputID, accepttauID, accept_tau)
  STATUS = NF90_PUT_VAR (outputID, acceptnuID, accept_nu)
  STATUS = NF90_PUT_VAR (outputID, acceptyID, accept_y_vector)
  STATUS = NF90_PUT_VAR (outputID, HoID, H)
  STATUS = NF90_PUT_VAR (outputID, zoID, z)
  STATUS = NF90_PUT_VAR (outputID, datenumberoID, datenumber)
  STATUS = NF90_PUT_VAR (outputID, timeindexnonzerooID, timeindex_nonzero)
  STATUS = NF90_PUT_VAR (outputID, xapoID, x_ap)
  STATUS = NF90_PUT_VAR (outputID, nobsoID, nobs)
  STATUS = NF90_PUT_VAR (outputID, sitespresentoID, sitespresent)
  STATUS = NF90_PUT_VAR (outputID, deltatimeoID, deltatime)
  STATUS = NF90_PUT_VAR (outputID, distanceoID, distance)
  STATUS = NF90_PUT_VAR (outputID, TindicesoID, T_indices)
  STATUS = NF90_PUT_VAR (outputID, distributionoID, distribution)
  STATUS = NF90_PUT_VAR (outputID, DoID, D)
  STATUS = NF90_CLOSE(outputID)
  
!!************************************************************************************************************

! All arrays that are not explicitly saved will be deallocated upon exit


 call system_clock(count_1,count_rate,count_max)
 t1=count_1/count_rate
 write(*,*) 'Total time elapsed (s):', t1-t0

!!*****************************************************************************
!!*****************************************************************************  
end program hierarchical_MCMC


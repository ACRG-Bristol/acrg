# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 16:30:21 2016

tdmcmc inputs.py

Inputs file for doing a transdimensional inversion

Inputs are fed into template_tdmcmc.py
In turn this script runs the fortran mcmc loop that does all the work

Please note: The default is not to perform parallel tempering and the rjmcmc 
is performed using a single chain on a single processor. If you think parallel
tempering is required this is easily changed. Speak to me (ML).

After the function definition at the top of the script the next section
is the inputs. This contains all the basic stuff you will probably want and 
have to change or tune.

Further changes to stepsizes and parameter uncertainties can be applied 
further down the script once the dimensions have been set. You will almost 
certainly need to tune your stepsizes appropriately for different parameter 
types (i.e. baseline and emissions). You will also want to adjust the 
uncertainties differently (emissions uncertainty > baseline uncertainty). 

@author: ml12574
"""
import numpy as np
import run_tdmcmc as acrg_tdmcmc
import acrg_name as name
#############################################################
#parser = argparse.ArgumentParser(description='This is a demo script by Mark.')
#parser.add_argument("start", help="Start date string yyyy-mm-dd")                  
#parser.add_argument("end", help="End date sting yyyy-mm-dd")
#args = parser.parse_args()
#
#start_date = args.start
#end_date = args.end
  
##############################################################
# Variables and definitions that need to be set
  
sites=['MHD', 'TAC']
meas_period = ['2H', '2H']   # Frequency to read in measurements 
av_period=['24H', '24H']   # Frequency to average footprints and measuerements for inversion

species='CH4'
start_date = '2014-03-01'
end_date = '2014-04-01'
domain='EUROPE'
network="test"

fp_basis_case = 'transd'
bc_basis_case = 'NESW'

################################################################
# SET OUTPUT DIRECTORY
output_dir="/path/to/output/directory/"

#######################################################
# DO YOU WANT TO DO REVERSIBLE JUMP OR NOT?????
rjmcmc=1           # 1 = do reversible jump; any other number = don't

# DO YOU WANT CORRELATED MEASUREMENTS OR NOT??
inv_type = 'uncorrelated'    # Options are 'correlated', 'evencorr', 'corr'

bl_period = 7      # No. of days for which each sigma_model value applies
kmin=4             # Minimum number of regions
kmax=400           # Maximum number of regions
k_ap = 50         # Starting number of regions
nIt=25000          # of iterations
burn_in=25000      # of discarded burn-in iterations 
nsub=10         # nsub=100=store every 100th iteration)
   
############################################################
# LOCALNESS THRESHOLD
threshold=0.4         # Fraction above which air is considered "local"
                      # Change to 1 if you don't want to separate high local times
       
#####################################################
# PARALLEL TEMPERING PARAMETERS          
nbeta=2            # Number of parallel chains - needs to be defined even if no Parallel tempering
beta= np.array((1.,1./2.)) # Values of beta for tempered chains 

#################################################################################
# DEFINE FORM OF PDFs FOR EMISSIONS, EMISSIONS UNCERTAINTIES AND MODEL UNCERTAINTY
x_pdf0=3        # 1 = UNIFORM, 2=GAUSSIAN, 3=LOGNORMAL  1st term for fixed terms, 2nd for variable
pdf_param1_pdf0 = 1
pdf_param2_pdf0 = 1
sigma_model_pdf = 1   

##########################################################################
# Hyperparameters of emissions parameters
pdf_param10=1.         # Mean of lognormal or normal distribution
pdf_param20=0.4        # Std. of lognormal or normal distribution

pdf_p1_hparam10=0.8    # Lower bound of uniform distribution of pdf_param1
pdf_p1_hparam20=1.2    # Upper bound of uniform distribution of pdf_param1
pdf_p2_hparam10=0.05   # Lower bound of uniform distribution of pdf_param2
pdf_p2_hparam20=0.5    # Upper bound of uniform distribution of pdf_param2

######################################################################
# Model-Measurement starting value and uncertainty parameters
sigma_model_ap = 20.   # Initial starting value of sigma_model (in same units as y, e.g ppb)
sigma_model_hparams = np.array([0.1*sigma_model_ap,5.*sigma_model_ap]) # upper and lower bounds of uniform dist. 
####################################################################
# DEFINE STEPSIZE FOR PROPOSAL DISTRIBUTIONS 
# TO TUNE INDIVIDUAL ELEMENTS SEE LINE 
stepsize=0.5   # Stepsize for proposal distirbution of x_update 
stepsize_pdf_p1=0.1   # Stepsize for proposal distirbution fof pdf_param1 
stepsize_pdf_p2=0.1   # Stepsize for proposal distirbution of pdf_param2 
stepsize_sigma_y=0.5  # Stepsize for proposal distribution of sigma_model

stepsize_clon = 8.    # Stepsize for longitude for move
stepsize_clat = 5.    # Stepsize for latitude for move
stepsize_bd=2.        # Stepsize for change in x during birth step

################################################
# TUNING OF INDIVIDUAL PARAMETER STEPSIZES AND UNCERTAINTIES

if fp_basis_case in('transd', 'INTEM'):
    nfixed = 8

if bc_basis_case in('NESW'):
    nBC = 4

if 'GOSAT' in(sites):
    nBias = 1           # Change as needed 
else: 
    nBias=0

nIC=nfixed+nBC+nBias
nIC1=nIC+1
kICmax=kmax+nIC    
pdf_param1 = np.zeros((kICmax,nbeta))
pdf_param2 = np.zeros((kICmax,nbeta))

if fp_basis_case in("INTEM"):
    basis_func = name.name.basis(domain = domain, basis_case = fp_basis_case)
    k_ap=len(np.unique(basis_func.basis.values))

##########################################################
# TUNE INDIVIDUAL STEPSIZES AND HYPERPARAMETER VALUES
"""
You'll need to run the inversion once first and then check the output file
for the acceptance ratios to see which elements need their stepsizes altered.
You should be able to get away with a few thousand iterations to get a good 
grasp on acceptance ratios. 
Find ratio in output file by doing:
ratio=post_mcmc.accepts/(post_mcmc.accepts+post_mcmc.rejects)
As ever aim for 25-50%  (Roberts et al. 1997) for each element.
The births deaths and moves are not so intuitive so maybe don't worry about
those so much.
"""
stepsize_all=np.zeros((nIC1))+stepsize
stepsize_pdf_p1_all=np.zeros((nIC1))+(stepsize_pdf_p1*pdf_param10)
stepsize_pdf_p2_all=np.zeros((nIC1))+(stepsize_pdf_p2*pdf_param20)
stepsize_all[:nBC]=stepsize_all[:nBC]/200.
stepsize_all[1:3]=stepsize_all[1:3]*70.
stepsize_pdf_p1_all[:nBC]=stepsize_pdf_p1_all[:nBC]/10.
stepsize_pdf_p2_all[:nBC]=stepsize_pdf_p2_all[:nBC]/10.

stepsize_all[-1]=stepsize_all[-1]/2.
stepsize_pdf_p1_all[-1]=stepsize_pdf_p1_all[-1]
stepsize_pdf_p2_all[-1]=stepsize_pdf_p2_all[-1]

pdf_param1[:,:]=pdf_param10
pdf_param2[:,:]=pdf_param20

pdf_p1_hparam1=np.zeros((nIC1))+pdf_p1_hparam10
pdf_p1_hparam2=np.zeros((nIC1))+pdf_p1_hparam20

pdf_p2_hparam1=np.zeros((nIC1))+pdf_p2_hparam10
pdf_p2_hparam2=np.zeros((nIC1))+pdf_p2_hparam20

x_pdf = np.zeros((nIC1), dtype=np.int)+x_pdf0
x_pdf[:nBC]=2
pdf_param1_pdf=np.zeros((nIC1), dtype=np.int)+pdf_param1_pdf0
pdf_param2_pdf=np.zeros((nIC1), dtype=np.int)+pdf_param2_pdf0

pdf_param2[:nBC,:]=0.05
pdf_p2_hparam1[:nBC]=0.01
pdf_p2_hparam2[:nBC]=0.1
pdf_p2_hparam1[-1]=0.2
pdf_p2_hparam2[-1]=2.
pdf_param2[nIC:,:]=1.

post_mcmc=acrg_tdmcmc.run_tdmcmc(sites,meas_period,av_period,species,start_date,end_date, 
    domain,network,fp_basis_case ,bc_basis_case,rjmcmc,bl_period,kmin,kmax,
    k_ap,nIt,burn_in,nsub,threshold,nbeta,beta,sigma_model_pdf,sigma_model_ap, 
    sigma_model_hparams,stepsize_sigma_y,stepsize_clon,stepsize_clat,
    stepsize_bd,stepsize_all,stepsize_pdf_p1_all,stepsize_pdf_p2_all,
    pdf_param1,pdf_param2,pdf_p1_hparam1,pdf_p1_hparam2,pdf_p2_hparam1,
    pdf_p2_hparam2,x_pdf ,pdf_param1_pdf,pdf_param2_pdf, inv_type, output_dir)
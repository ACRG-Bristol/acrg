# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 16:30:21 2016

tdmcmc_inputs.py

Run file for doing a transdimensional inversion

Inputs are fed into template_tdmcmc.py
In turn this script runs the fortran mcmc loop that does all the work

Please note: The default is currently to perform parallel tempering, 
using cpus = nbeta parameter.

All inputs are defined within a configuration file. By default this is called
"param.ini" and should be within the acrg_tdmcmc directory. This 
can be changed within the script (or by feeding inputs from the command line).

Three parameters can be specified externally if running on the command line:
    - start_date (positional - first argument supplied)
    - end_date (positional - second argument supplied)
    - configuration file (using -c flag)
e.g. 
>> python tdmcmc_inputs.py 2015-01-01 2015-04-01 -c param.ini

Further changes to stepsizes and parameter uncertainties can be applied within 
this script.
You will almost certainly need to tune your stepsizes appropriately for 
different parameter types (i.e. baseline and emissions). You will also want 
to adjust the uncertainties differently (emissions uncertainty > baseline uncertainty). 

@author: ml12574 (updated by rt17603)
"""

import os
import argparse
import glob
import shutil
import numpy as np
import datetime as dt
import acrg_name as name
from acrg_tdmcmc import run_tdmcmc 
from acrg_tdmcmc import tdmcmc_post_process as process
import tdmcmc_config
#import acrg_config as configread

acrg_path = os.getenv("ACRG_PATH")
data_path = os.getenv("DATA_PATH")

config_file = 'param.ini'
config_path = os.path.join(acrg_path,"acrg_tdmcmc")
config_file = os.path.join(config_path,config_file)

#############################################################
parser = argparse.ArgumentParser(description='Running mcmc script')
parser.add_argument("start", help="Start date string yyyy-mm-dd",nargs="?")                  
parser.add_argument("end", help="End date sting yyyy-mm-dd",nargs="?")
parser.add_argument("-c","--config",help="Configuration filename",default=config_file)
args = parser.parse_args()

start_date = args.start
end_date = args.end
config_file = args.config

#############################################################

verbose = True
if verbose:
    print '\n---------------\n'
    print 'Input configuration file: {0}'.format(config_file)

# Extract parameters from configuration file
param = tdmcmc_config.all_mcmc_param(config_file)

# Variables and definitions that need to be set

sites = param['sites']
meas_period = param['meas_period'] # Frequency to read in measurements 
av_period = param['av_period']     # Frequency to average footprints and measuerements for inversion

species = param['species']
if not start_date:
    start_date = param['start_date']
    if not start_date:
        raise Exception("Start date of observations must be specified either on the command line or within the configuration file")
if not end_date:
    end_date = param['end_date']
    if not end_date:
        raise Exception("End date of observations must be specified either on the command line or within the configuration file")
domain = param['domain']
network = param['network']

fp_basis_case = param['fp_basis_case']
bc_basis_case = param['bc_basis_case']

if verbose:
    print 'Measurement details: sites - {0}, species - {1}, domain - {2}, network - {3}'.format(sites,species,domain,network)
    print 'Date range: {0} - {1}'.format(start_date,end_date)
    print 'Basis case for footprint: {0}'.format(fp_basis_case)
    print 'Basis case for boundary conditions: {0}\n'.format(bc_basis_case)
    print '\n---------------\n'

################################################################
# SET OUTPUT DIRECTORY AND OUTPUT DETAILS
output_dir = param['output_dir']
unique_copy = param['unique_copy']

if output_dir == "/path/to/output/directory/":
    raise Exception("Please set output directory (output_dir) parameter within configuration file: {0}".format(config_file))
elif output_dir.startswith("$ACRG_PATH"):
    output_dir = output_dir.replace("$ACRG_PATH",acrg_path)
elif output_dir.startswith("$DATA_PATH"):
    output_dir = output_dir.replace("$DATA_PATH",data_path)

#######################################################
# DO YOU WANT TO DO REVERSIBLE JUMP OR NOT?????
reversible_jump = param['reversible_jump']         # True = do reversible jump; False = don't
parallel_tempering = param['parallel_tempering']   # True = do parallel tempering

# DO YOU WANT CORRELATED MEASUREMENTS OR NOT??
inv_type = param['inv_type']  # Options are 'uncorrelated', 'evencorr', 'corr'

# Have to re-create value explicitly to ensure any 'is' statements return True
#if inv_type == 'uncorrelated':
#    inv_type = 'uncorrelated'
#elif inv_type == 'evencorr':
#    inv_type = 'evencorr'
#elif inv_type == 'corr':
#    inv_type = 'corr'

kmin = param['kmin']          # Minimum number of regions
kmax = param['kmax']          # Maximum number of regions
k_ap = param['k_ap']          # Starting number of regions
nIt = param['nIt']            # of iterations
burn_in = param['burn_in']    # of discarded burn-in iterations 
nsub = param['nsub']          # nsub=100=store every 100th iteration)

if verbose:
    print 'Inversion type: {0}'.format(inv_type)
    print 'Regions in trans-dimesional grid - minimum allowed: {0}, maximum allowed: {1}, starting value: {2}'.format(kmin,kmax,k_ap)
    print 'Burn-in iterations: {0}'.format(burn_in)
    print 'Number of iterations to run: {0} (nsub = {1}, {2} iterations will be saved)\n'.format(nIt,nsub,nIt/nsub)
    print '\n---------------\n' 
    
############################################################
# FILTERS
"""
# Include a list of filters you want to include here to filter the data
# Purpose is to remove observation times with potential biases

Options are:
"daily_median": What it says
 "daytime": Only between 11:00 - 15:00 inclusive 
 "nighttime": Only b/w 23:00 - 03:00 inclusive
 "noon": Only 12:00 fp and obs used
 "pblh_gt_500":
 "pblh_gt_250": 
 "local_influence": Only keep times when localness is low
 "six_hr_mean":
 "ferry_loc": GAUGE-FERRY specific - Used to filter out dodgy ferry locations
 "ferry_mf": GAUGE-FERRY specific - Used to filter out dodg ferry 
 "ferry_fp_zero": GAUGE-FERRY specific 

"""
filters = param['filters']
       
#####################################################
# PARALLEL TEMPERING PARAMETERS     
# PARALLEL TEMPERING PARAMETERS     
if parallel_tempering == False:          
    nbeta=2            # Number of parallel chains - needs to be defined even if no Parallel tempering
    beta= np.array((1.,1./2.)) # Values of beta for tempered chains 
elif parallel_tempering == True: 
    nbeta = param['nbeta']    # Number of parallel chains - needs to be defined even if no Parallel tempering
    series = np.linspace(0.,1.,num=nbeta)
    beta_1 = np.exp(0 + (np.log(250)-0)*series)
    beta= 1./beta_1


#################################################################################
# DEFINE FORM OF PDFs FOR EMISSIONS, EMISSIONS UNCERTAINTIES AND MODEL UNCERTAINTY
x_pdf0 = param['x_pdf0']   # 1 = UNIFORM, 2=GAUSSIAN, 3=LOGNORMAL  1st term for fixed terms, 2nd for variable
pdf_param1_pdf0 = param['pdf_param1_pdf0']
pdf_param2_pdf0 = param['pdf_param2_pdf0']
sigma_model_pdf = param['sigma_model_pdf'] 

##########################################################################
# Hyperparameters of emissions parameters
pdf_param10 = param['pdf_param10']         # Mean of lognormal or normal distribution
pdf_param20 = param['pdf_param20']         # Std. of lognormal or normal distribution

pdf_p1_hparam10 = param['pdf_p1_hparam10']  # Lower bound of uniform distribution of pdf_param1
pdf_p1_hparam20 = param['pdf_p1_hparam20']  # Upper bound of uniform distribution of pdf_param1
pdf_p2_hparam10 = param['pdf_p2_hparam10']  # Lower bound of uniform distribution of pdf_param2
pdf_p2_hparam20 = param['pdf_p2_hparam20']  # Upper bound of uniform distribution of pdf_param2

######################################################################
# Model-Measurement starting value and uncertainty parameters
sigma_model_ap = param['sigma_model_ap']   # Initial starting value of sigma_model (in same units as y, e.g ppb)
sigma_model_hparams = param['sigma_model_hparams'] # upper and lower bounds of uniform dist. - percentages of sigma_model_ap

bl_period = param['bl_period']   # No. of days for which each sigma_model value applies 
bl_split = param['bl_split']     # Set to true if want to split sigma_model values by BL depth rather than days
levels = param['levels']         # Banding of bL depths to solve for different sigma_model
if levels == 0:                  # e.g. levels=[0.,500.,1000.,10000.] Set if bl_split = True
    levels=None         

####################################################################
# DEFINE STEPSIZE FOR PROPOSAL DISTRIBUTIONS 
# TO TUNE INDIVIDUAL ELEMENTS SEE LINE 
stepsize = param['stepsize']                 # Stepsize for proposal distirbution of x_update 
stepsize_pdf_p1 = param['stepsize_pdf_p1']   # Stepsize for proposal distirbution fof pdf_param1 
stepsize_pdf_p2 = param['stepsize_pdf_p2']   # Stepsize for proposal distirbution of pdf_param2 
stepsize_sigma_y = param['stepsize_sigma_y'] # Stepsize for proposal distribution of sigma_model

stepsize_clon = param['stepsize_clon']       # Stepsize for longitude for move
stepsize_clat = param['stepsize_clat']       # Stepsize for latitude for move
stepsize_bd = param['stepsize_bd']           # Stepsize for change in x during birth step

################################################
# TAU
#"Only need if inv_type = ('evencorr', 'correlated'):"
tau_ap = param['tau_ap']
tau_hparams = param['tau_hparams']
stepsize_tau = param['stepsize_tau'] 
tau_pdf = param['tau_pdf'] 

# TUNING OF INDIVIDUAL PARAMETER STEPSIZES AND UNCERTAINTIES
################################################
# SET DIMENSIONS of nBC and nIC
f_list=glob.glob(data_path + "/NAME/basis_functions/" 
                    + domain + "/" + fp_basis_case + 
                    "_" + domain + "_*.nc") 

if len(f_list) > 0:
    ds = process.open_ds(f_list[0]) 
    if fp_basis_case in('sub-transd'):
        nfixed = len(np.unique(ds.basis))-1
    else:
        nfixed = len(np.unique(ds.basis))
else:
    raise LookupError("No file exists for that fp_basis_case and domain")

f_list2=glob.glob(data_path + "/NAME/bc_basis_functions/" 
                    + domain + "/" + bc_basis_case + 
                    "_" + domain + "_*.nc") 
if len(f_list2) > 0:                    
    ds2 = process.open_ds(f_list2[0]) 
    nBC = len(ds2.region)
else:
    raise LookupError("No file exists for that bc_basis_case and domain")

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
Stepsize tuning occurs automatically during the burn-in period. 
However, some manual tuning may be required to start with so step sizes 
are at least of the right order of magnitude.

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

post_mcmc=run_tdmcmc.run_tdmcmc(sites, meas_period, av_period, species, start_date, end_date,  
    domain, network, fp_basis_case, bc_basis_case, reversible_jump, parallel_tempering,    
    bl_period, kmin, kmax, k_ap, nIt, burn_in ,nsub,    
    nbeta, beta, sigma_model_pdf, sigma_model_ap,     
    sigma_model_hparams, stepsize_sigma_y, stepsize_clon, stepsize_clat,    
    stepsize_bd, stepsize_all, stepsize_pdf_p1_all, stepsize_pdf_p2_all,    
    pdf_param1, pdf_param2, pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1,    
    pdf_p2_hparam2, x_pdf, pdf_param1_pdf, pdf_param2_pdf, inv_type,     
    output_dir,filters=filters,bl_split=bl_split, bl_levels=levels,
    tau_ap=tau_ap, tau_hparams=tau_hparams, stepsize_tau=stepsize_tau, tau_pdf=tau_pdf)

if unique_copy:
    # Create date-stamped sub-directory of the form:
    #   Output_"sites"_"species"_"start_date"_"creation_dt" e.g. Output_MHD-TAC_CH4_20080101_20171110T12-00-00
    now = dt.datetime.now().replace(microsecond=0)
    site_str = '-'.join(sites)
    start_date = start_date.replace('-','')
    creation_dt = dt.datetime.strftime(now,'%Y%m%dT%H-%M-%S')
    sub_dir = 'Output_{0}_{1}_{2}_{3}'.format(site_str, species, start_date, creation_dt)
    datestamp_output_dir = os.path.join(output_dir,sub_dir)
    
    print 'Creating new unique directory for output: {0}'.format(datestamp_output_dir)
    os.mkdir(datestamp_output_dir)    
    
    # Copy config file to output sub-directory
    shutil.copy(config_file,datestamp_output_dir)
    
    # Write output from MCMC code again but this time to the subdirectory
    fname=os.path.join(datestamp_output_dir,
                            "output_" + network + "_" + species + "_" + start_date + ".nc")
    for key in post_mcmc.keys():
        post_mcmc[key].encoding['zlib'] = True
    post_mcmc.to_netcdf(path=fname, mode='w')
else:
    shutil.copy(config_file,output_dir)


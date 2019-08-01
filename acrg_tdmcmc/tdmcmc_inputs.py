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
from __future__ import print_function
from __future__ import division
from builtins import str

from past.utils import old_div
import os
import argparse
import glob
import shutil
import sys
import numpy as np
import pandas as pd
import datetime as dt
import acrg_name as name
from acrg_tdmcmc import run_tdmcmc 
from acrg_tdmcmc import tdmcmc_post_process as process
import acrg_tdmcmc.tdmcmc_config as tdmcmc_config
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
parser.add_argument("-c","--config",help="Name (including path) of configuration file",default=config_file)
parser.add_argument("-g","--group",help="Group to apply to run. All members of this group will be saved to the same sub-folder named after the group.")
parser.add_argument("-r","--regenerate",action="store_true",help="Regenerate configuration file from template and exit. This will overwrite any pre-existing config file.")
args = parser.parse_args()

regenerate_then_exit = args.regenerate
start_date = args.start
end_date = args.end
config_file = args.config
group = args.group

if regenerate_then_exit:
    tdmcmc_config.regenerate_tdmcmc_config(config_file)
    print("Note: new configuration file includes all parameters including optional ones.")
    sys.exit("Exiting.")

#############################################################

verbose = True
if verbose:
    print('\n---------------\n')
    print('Input configuration file: {0}'.format(config_file))

# Extract parameters from configuration file
param = tdmcmc_config.all_mcmc_param(config_file)

# Variables and definitions that need to be set

sites = param['sites']
meas_period = param['meas_period'] # Frequency to read in measurements 
av_period = param['av_period']     # Frequency to average footprints and measuerements for inversion

max_level = param["max_level"] # Only relevant for satellite data

if "site_modifier" in param and param["site_modifier"] is not None:
    site_modifier = param["site_modifier"]
else:
    site_modifier = {}

species = param['species']
if not start_date:
    start_date = str(param['start_date'])
    if not start_date:
        raise Exception("Start date of observations must be specified either on the command line or within the configuration file")
else:
    start_date = str(start_date)
if not end_date:
    end_date = str(param['end_date'])
    if not end_date:
        raise Exception("End date of observations must be specified either on the command line or within the configuration file")
else:
    end_date = str(end_date)
    
domain = param['domain']
network = param['network']

fp_basis_case = param['fp_basis_case']
bc_basis_case = param['bc_basis_case']

if verbose:
    print('Measurement details: sites - {0}, species - {1}, domain - {2}, network - {3}'.format(sites,species,domain,network))
    print('Date range: {0} - {1}'.format(start_date,end_date))
    print('Basis case for footprint: {0}'.format(fp_basis_case))
    print('Basis case for boundary conditions: {0}\n'.format(bc_basis_case))
    print('\n---------------\n')

################################################################
# SET OUTPUT DIRECTORY AND OUTPUT DETAILS
output_dir = param['output_dir']

data_dir = param['data_dir']
fp_dir = param['fp_dir']
flux_dir = param['flux_dir']
basis_dir = param['basis_dir']
bc_basis_dir = param['bc_basis_dir']
bc_dir = param['bc_dir']


unique_copy = param['unique_copy']

if output_dir == "/path/to/output/directory/":
    raise Exception("Please set output directory (output_dir) parameter within configuration file: {0}".format(config_file))
elif output_dir.startswith("$ACRG_PATH"):
    output_dir = output_dir.replace("$ACRG_PATH",acrg_path)
elif output_dir.startswith("$DATA_PATH"):
    output_dir = output_dir.replace("$DATA_PATH",data_path)

if not os.path.isdir(output_dir):
    raise Exception("Output directory: {} does not exist.".format(output_dir))
    
if group:
    output_dir = os.path.join(output_dir,group)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
if data_dir is not None and data_dir.startswith("$ACRG_PATH"):
    data_dir = data_dir.replace("$ACRG_PATH",acrg_path)
    print(("data_directory is ", data_dir))

if fp_dir is not None and fp_dir.startswith("$ACRG_PATH"):
    fp_dir = fp_dir.replace("$ACRG_PATH",acrg_path)
    print(("fp_directory is ", fp_dir))

if flux_dir is not None and flux_dir.startswith("$ACRG_PATH"):
    flux_dir = flux_dir.replace("$ACRG_PATH",acrg_path)
    print(("flux_directory is ", flux_dir))

if bc_dir is not None and bc_dir.startswith("$ACRG_PATH"):
    bc_dir = bc_dir.replace("$ACRG_PATH",acrg_path)
    print(("bc_directory is ", bc_dir))    

if basis_dir is not None and basis_dir.startswith("$ACRG_PATH"):
    basis_dir = basis_dir.replace("$ACRG_PATH",acrg_path)
    print(("basis_directory is ", basis_dir))

if bc_basis_dir is not None and bc_basis_dir.startswith("$ACRG_PATH"):
    bc_basis_dir = bc_basis_dir.replace("$ACRG_PATH",acrg_path)
    print(("bc_basis_directory is ", bc_basis_dir))


################################################
# CALCULATE PRIOR UNCERTAINTY
# This switches modes to allow calculation prior uncertainty reduction
# This is as alternative to completing a posterior estimation based on the input measurements
# **This parameter should be set to True only in those circumstances.**

if "prior_uncertainty" in param:
    prior_uncertainty = param["prior_uncertainty"]
    if prior_uncertainty == True:
        print("WARNING: Run will calculate prior uncertainty reduction by setting dmf to a large number so there is no improvement based on measurements.")
        print("Output of this run should NOT be used as an improved posterior output.")
else:
    prior_uncertainty = False


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
    print('Inversion type: {0}'.format(inv_type))
    print('Regions in trans-dimensional grid - minimum allowed: {0}, maximum allowed: {1}, starting value: {2}'.format(kmin,kmax,k_ap))
    print('Burn-in iterations: {0}'.format(burn_in))
    print('Number of iterations to run: {0} (nsub = {1}, {2} iterations will be saved)\n'.format(nIt,nsub,old_div(nIt,nsub)))
    print('\n---------------\n') 
    
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
# nfixed - number of fixed areas within the basis function (nfixed = number from basis function -1 for trans-d model)
# nBC - number of fixed areas within boundary condition basis function
# nBIAS - seems to be specific to GOSAT?
# nIC - a combination of all the above dimesions

f_list=glob.glob(data_path + "/LPDM/basis_functions/" 
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

f_list2=glob.glob(data_path + "/LPDM/bc_basis_functions/" 
                    + domain + "/" + bc_basis_case + 
                    "_" + domain + "_*.nc") 
if len(f_list2) > 0:                    
    ds2 = process.open_ds(f_list2[0]) 
    nBC_basis = len(ds2.region)
else:
    raise LookupError("No file exists for that bc_basis_case and domain")

for site in sites:
    if 'GOSAT' in site and len(sites) > 1: # Will need to update to base on platform rather than searching for "GOSAT"
        nBias = 1
        break
else: 
    nBias = 0

## Define nBC based on time, so each month is scaled individually
pd_start=pd.to_datetime(start_date)
pd_end=pd.to_datetime(end_date)

# Calculate number of months in inversion period
if pd_end.day == 1:
    nmonths = pd_end.to_period('M') - pd_start.to_period('M')   
else:
    nmonths = pd_end.to_period('M') - pd_start.to_period('M')+1
    

nBC = nBC_basis*nmonths   # No. of bc_basis functions x nmonths

nBC+=nBias
nIC=nfixed+nBC
#nIC=nfixed+nBC+nBias
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
    output_dir,fp_dir=fp_dir, flux_dir = flux_dir, data_dir=data_dir, basis_dir=basis_dir, bc_basis_dir=bc_basis_dir, bc_dir = bc_dir,
    filters=filters,bl_split=bl_split, bl_levels=levels,
    tau_ap=tau_ap, tau_hparams=tau_hparams, stepsize_tau=stepsize_tau, tau_pdf=tau_pdf,
    max_level=max_level, site_modifier=site_modifier,prior_uncertainty=prior_uncertainty)

if unique_copy:
    shutil.copy(config_file,output_dir)
    # Create date-stamped sub-directory of the form:
    #   Output_"sites"_"species"_"start_date"_"creation_dt" e.g. Output_MHD-TAC_CH4_20080101_20171110T12-00-00
    now = dt.datetime.now().replace(microsecond=0)
    site_str = '-'.join(sites)
    #start_date_d = start_date.replace('-','')
    creation_dt = dt.datetime.strftime(now,'%Y-%m-%dT%H-%M-%S')
    sub_dir = 'Output_{0}_{1}_{2}_{3}'.format(site_str, species, start_date, creation_dt)
    datestamp_output_dir = os.path.join(output_dir,sub_dir)
    
    print('Creating new unique directory for output: {0}'.format(datestamp_output_dir))
    os.mkdir(datestamp_output_dir)    
    
    # Copy config file to output sub-directory
    shutil.copy(config_file,datestamp_output_dir)
    
    network_w = network.split('/')[-1]
    # Write output from MCMC code again but this time to the subdirectory
    fname = os.path.join(datestamp_output_dir,
                            "output_{network}_{species}_{date}.nc".format(network=network_w,species=species,date=start_date))
    for key in list(post_mcmc.keys()):
        post_mcmc[key].encoding['zlib'] = True
    post_mcmc.to_netcdf(path=fname, mode='w')
    
    if prior_uncertainty:
        pu_string = "This run is to calculate the prior uncertainty reduction and does NOT provide an improved posterior output based on measurements.\n"
        pu_string += "Note: Measurement uncertainity was set to an artifically high value (mf*10000) to achieve this."
        pu_fname = os.path.join(datestamp_output_dir,"README_uncertainty_reduction.txt")
        pu_warn_file = open(pu_fname,'w')
        pu_warn_file.write(pu_string)
        pu_warn_file.close()
else:
    shutil.copy(config_file,output_dir)


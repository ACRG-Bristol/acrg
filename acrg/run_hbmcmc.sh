#!/bin/bash
#SBATCH --partition=compute,short
#SBATCH --mem=64GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=sathbmcmc1
#SBTACH --output=sathbmcmc1
#SBATCH --time=24:00:00


## missing octreal, novreal,
#module load lang/python/anaconda/3.7-2019.10

source activate /user/work/ef17148/oldstuff/ef17148/.conda/envs/acrg

python /user/work/ef17148/acrg/acrg/hbmcmc/run_hbmcmc.py 2016-06-01 2016-07-01 -c /user/work/ef17148/acrg/acrg/hbmcmc/hbmcmc_input_GOSAT-brazil_v1.ini
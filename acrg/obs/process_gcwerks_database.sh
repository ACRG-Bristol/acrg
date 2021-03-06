#!/bin/sh
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=00:05:00
#PBS -N gcwerks_db

source ~/.bashrc

conda activate $ACRG_CONDA_ENV

# Find ACRG path
acrg_path=`python -c 'from acrg.config.paths import Paths; print(Paths.acrg)'`

# Run database processing script
python -c "import acrg.obs as acrg_obs; acrg_obs.utils.obs_database()" > $acrg_path/tmp/process_gcwerks_database

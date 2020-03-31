#!/bin/sh
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 1-20
# #PBS -o process_gcwerks_^array_index^_stdout
# #PBS -e process_gcwerks_^array_index^_stderr
#PBS -l walltime=00:02:00
#PBS -N process_gcwerks

/bin/bash
conda activate matt

acrg_path=`python -c 'from acrg_config.paths import paths; print(paths.acrg)'`

# mkdir $acrg_path/tmp
# cd ~/code/test/test$PBS_ARRAY_INDEX

python -c "import acrg_obs.process_gcwerks as process; process.array_job($PBS_ARRAY_INDEX)" > $acrg_path/tmp/process_gcwerks_$PBS_ARRAY_INDEX

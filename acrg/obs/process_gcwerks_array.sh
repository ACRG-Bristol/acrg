#!/bin/sh
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=00:05:00
#PBS -N process_gcwerks

# run this qsub line to call this array shell script to process all AGAGE data
# qsub -V -J 1-40 -k oe -j oe -Wsandbox=PRIVATE process_gcwerks_array.sh

/bin/bash
source ~/.bashrc

conda activate $ACRG_CONDA_ENV

acrg_path=`python -c 'from acrg.config.paths import Paths; print(Paths.acrg)'`

python -c "import acrg.obs.process_gcwerks as process; process.array_job($PBS_ARRAY_INDEX)" > $acrg_path/tmp/process_gcwerks_$PBS_ARRAY_INDEX

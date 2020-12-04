#!/bin/sh

source ~/.bashrc

conda activate $ACRG_CONDA_ENV

# Find ACRG path
acrg_path=`python -c 'from acrg_config.paths import paths; print(paths.acrg)'`

# Get files from DAGAGE2
$acrg_path/acrg_obs/process_gcwerks_dagage2.sh

# Run array processing job
jobout=`qsub -V -J 1-40 -k oe -j oe -Wsandbox=PRIVATE $acrg_path/acrg_obs/process_gcwerks_array.sh`
jobid=${jobout%.*}

echo "Launched array processing job $jobout"
echo "Launching db update after job $jobid"

# Launch database job
qsub -V  -k oe -j oe -Wsandbox=PRIVATE -W depend=afterany:$jobid $acrg_path/acrg_obs/process_gcwerks_database.sh

echo "Done"
#!/bin/bash
#SBATCH --job-name=slurm_script_test
#SBATCH --output=/user/work/bm13805/slurm_script_test/logs/test_%a.out
#SBATCH --error=/user/work/bm13805/slurm_script_test/logs/test_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=40gb
#SBATCH --account=chem007981
#SBATCH --array=1-48

# Set up Python environment
module purge
module load lang/python/anaconda
module load tools/git
source /user/work/bm13805/.venv_inversions_standard/bin/activate

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
    mkdir -p /user/work/bm13805/slurm_script_test/logs
    pip list >> /user/work/bm13805/slurm_script_test/logs/packages.txt
fi

inversions_path=$(pip list | grep "openghg-inversion" | awk '{ print $3 }')
git_branch=$(git -C $inversions_path status | awk 'NR==1{ print $3 }')
git_commit=$(git -C $inversions_path log --oneline -n 1 $git_branch | awk '{ print $1 }')
echo "Using commit $git_commit on branch $git_branch in repo $inversions_path" >> /user/work/bm13805/slurm_script_test/logs/git_info.txt

# Specify the path to the config file
config=/user/work/bm13805/slurm_script_test/config.txt
nkwargs=5

# Extract start and end dates for the current $SLURM_ARRAY_TASK_ID
start=$(awk -v array_task_id=$SLURM_ARRAY_TASK_ID -F'\t' '$1==array_task_id {print $2}' $config)
end=$(awk -v array_task_id=$SLURM_ARRAY_TASK_ID -F'\t' '$1==array_task_id {print $3}' $config)

# Keyword arg handling
if [ $nkwargs -gt 3 ]; then
    # Extract keyword args for the current $SLURM_ARRAY_TASK_ID
    kwargs=$(awk -v array_task_id=$SLURM_ARRAY_TASK_ID -F'\t' '$1==array_task_id {print $4}' $config)

    # Extract output dir name corresponding to keyword args
    child_out_dir=$(awk -v array_task_id=$SLURM_ARRAY_TASK_ID -F'\t' '$1==array_task_id {print $5}' $config)
fi


# Print a message with current $SLURM_ARRAY_TASK_ID, start date, and end date
echo "Running array task ${SLURM_ARRAY_TASK_ID}: start date ${start}, end date is ${end}."

# Check if there are kwargs
if [ $nkwargs -le 3 ]; then
    python ${inversions_path}/openghg_inversions/hbmcmc/run_hbmcmc.py "${start}" "${end}" -c /user/work/bm13805/ini_files/test.ini --output-path=/user/work/bm13805/slurm_script_test
else
    echo "Keyword args: ${kwargs}"
    echo "Output directory name: ${child_out_dir}"
    python ${inversions_path}/openghg_inversions/hbmcmc/run_hbmcmc.py "${start}" "${end}" -c /user/work/bm13805/ini_files/test.ini --kwargs="${kwargs}" --output-path="/user/work/bm13805/slurm_script_test/${child_out_dir}"
fi

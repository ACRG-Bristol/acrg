#!/bin/bash
#SBATCH --partition=compute,short
#SBATCH --mem=64GB
#SBATCH --ntasks=1
#SBATCH --array=2-12
#SBATCH --cpus-per-task=1
#SBATCH --job-name=bc_test
#SBATCH --output=bc_test_e%a.out
#SBATCH --time=24:00:00
#SBATCH --account=chem007981
## chem007981 SEMT030444

source activate /user/work/ef17148/oldstuff/ef17148/.conda/envs/acrg_new

month=$SLURM_ARRAY_TASK_ID
year=2017

# Format the month to always have two digits
formatted_month=$(printf "%02d" $month)

# Define the start date as the first day of the current month
start_date="$year-$formatted_month-01"

# Define the end date as the first day of the next month
# If the month is December (12), set the end date to December 31
if [ "$formatted_month" -eq "12" ]; then
    end_date="$year-12-31"
else
    next_month=$(printf "%02d" $(($month + 1)))
    end_date="$year-$next_month-01"
fi    

echo $start_date $end_date

python /user/work/ef17148/acrg/acrg/hbmcmc/run_hbmcmc.py $start_date $end_date -c /user/work/ef17148/acrg/acrg/bc_inversions/load_file_true.ini

date

python /user/work/ef17148/acrg/acrg/hbmcmc/run_hbmcmc.py $start_date $end_date -c /user/work/ef17148/acrg/acrg/bc_inversions/load_file_pred.ini

date

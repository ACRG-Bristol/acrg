#!/bin/sh
#SBATCH --job-name=sf6-bc
#SBATCH --output=sf6-bc.out
#SBATCH --error=sf6-bc.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=30gb
#SBATCH --account=chem007981

# Set up Python environment
module purge
module load lang/python/anaconda
source /user/work/bm13805/.venv_inversions_standard/bin/activate


python baseline.py "sf6" 2018 5 "1e-12" -o "/user/work/bm13805/paris_flat" -u "/user/work/bm13805/paris_runs/sf6" --standardise --bc-input "mhd_clear_sector" --store "user"

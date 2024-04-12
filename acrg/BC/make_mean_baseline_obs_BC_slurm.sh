#!/bin/sh
#SBATCH --job-name=make-bc
#SBATCH --output=/user/work/qq24644/my_paris/make_BC/make-bc.out
#SBATCH --error=/user/work/qq24644/my_paris/make_BC/make-bc.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=30gb
#SBATCH --account=geog011430

# Set up Python environment
module purge
module load lang/python/anaconda
source /user/home/qq24644/openghg_inv/bin/activate

# List of species to process
species=("hfc143a2" "c3f8" "hfc125" "hfc227ea" "hfc32" "sf6" "clch2ch2cl" "hfo1234yf" "hfo1234zee" "hcfo1233zde")

# Process each species
for spec in "${species[@]}"; do
    echo "Processing $spec"
    python /user/home/qq24644/acrg/acrg/BC/make_mean_baseline_obs_BC.py "$spec" 2018 5 "1e-12" \
    -o "/user/home/qq24644/work/my_paris/make_BC" \
    -u "/user/home/qq24644/work/my_paris/make_BC" \
    --standardise --bc-input "mhd_clean_sector" \
    --output-store "paris_openghg_store" \
    --obs-store "obs_nir_2024_01_25_store" \
    --fp-store "paris_openghg_store"
done
#!/bin/bash

air_shared_folder="/data/shared/"
snowy_shared_folder="/shared_data/snowy/shared/"


# Air to Snowy, one way
##############################################

folders="NAME/basis_functions
         NAME/boundary_conditions
         NAME/emissions
         NAME/fp_netcdf"

for f in $folders
do
    dirsync -svc "$air_shared_folder$f" "$snowy_shared_folder$f"
done


# Air to Snowy, two way
##############################################

folders="obs/"

for f in $folders
do
    dirsync -svc2 "$air_shared_folder$f" "$snowy_shared_folder$f"
done

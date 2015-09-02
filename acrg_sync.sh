#!/bin/bash

air_shared_folder="/data/shared/"
snowy_shared_folder="/shared_data/snowy/shared/"

# Air to Snowy, one way
##############################################

folders="NAME/basis_functions
         NAME/bc
         NAME/emissions
         NAME/fp
         NAME/bc_basis_functions
         NAME/countries"

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

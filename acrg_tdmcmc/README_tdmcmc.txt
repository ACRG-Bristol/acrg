README file for acrg_tdmcmc

Files in this directory:

acrg_hbtdmcmc.f90  - Fortran script for running transdimensional reversible jump mcmc
tdmcmc_template.py - Python template script to create right inputs for acrg_hbtdmcmc.f90

tdmcmc_post_process.py - Python script containing routines for plotting tdmcmmc output and generating country totals
post_process_template.py - Template file to show how you might call functions contained in tdmcmc_post_process.py


***** Templates should be copied to a local area for you to edit *****************
****** They exist for reference only, at some point you have to actually do some work yourself!!!!! *************************

*** Compiling with f2py ***

The file hbtdmcmc.so was created with:

f2py -c -m hbtdmcmc --f90flags='-fopenmp' -lgomp acrg_hbtdmcmc.f90

This compiles with gfortran, with openmp. To compile for serial just do:

f2py -c -m hbtdmcmc acrg_hbtdmcmc.f90

To compile using intel compiler use: 

f2py -c -m hbtdmcmc --fcompiler=intelem acrg_hbtdmcmc.f90

*********** WARNING - INTELEM DOES NOT WORK WITH OPENMP - it compiles but won't run. ***************


*** Running from python ***

import module_name

Where module_name is the name of module_name.so

When calling the trasdimensional mcmc function follow what's in the template file, but if using a different module name make sure you change that.

************** If using openMP then you have to run from ipython not Spyder *********************************




README file for acrg_tdmcmc

Files in this directory:

acrg_hbtdmcmc_uncorr.f90  - Fortran script for running uncorrelated transdimensional, hierarchical reversible jump mcmc
acrg_hbtdmcmc_evencorr.f90  - Fortran script for running temporally correlated transdimensional, hierarchical reversible jump mcmc. Assumes every data point is evenly spaced in time.
acrg_hbtdmcmc_corr.f90  - Fortran script for running temporally correlated transdimensional, hierarchical reversible jump mcmc. Suitable for any separation of data points.

tdmcmc_inputs.py - Python template script to create right inputs for template_tdmcmc.f90

run_tdmcmc.py - Python script that reads in all specified inputs and calls fortran subroutine, writes ouput file and returns dataset of outputs 
                   - Other than major changes you shouldn't have to touch this on a day-to-day basis

tdmcmc_post_process.py - Python script containing routines for plotting tdmcmmc output and generating country totals
post_process_inputs.py - Template file to show how you might call functions contained in tdmcmc_post_process.py


***** Templates and inputs should be copied to a local area for you to edit *****************
****** They exist for reference only, at some point you have to actually do some work yourself!!!!! *************************

In order for the run_tdmcmc.py script to work you need to have compiled both Fortran subroutines using f2py.

It is suggested that each fortran code is compiled twice - once for use with parallel tempering and one without.
The difference to make to the fortran before compiling is simply the inclusion (or exclusion) of the line:

call OMP_SET_NUM_THREADS(nbeta)      ! Uncomment for Parallel Tempering

This needs to be included for parallel tempering but should be commented out otherwise.

By default the .so files must be called:
tdmcmc_uncorr.so, tdmcmc_uncorr_pt.so, tdmcmc_evencorr.so, tdmcmc_evencorr_pt.so, tdmcmc_corr.so, tdmcmc_corr_pt.so

("_pt" - for parallel tempering versions i.e. with OMP_SET_NUM_THREADS line included).

But you don't have to stick with these defaults.

Known issues
^^^^^^^^^^^^

If you get a segmentation fault for an OMP run, then you will need to set the memory allocation using:

export KMP_STACKSIZE=512m

Sometimes this doesn't seem to work either, in which case also try

export OMP_STACKSIZE=512m
export GOMP_STACKSIZE=512m

This appears to be foolproof, but not sure which of these commands is key. Stick them in your .bashrc so you don't have to type them every time.

**************** Module load and setting up environments ******************************************

For most servers (except snowy and air) it is often the case that you will need to load at additional 
libraries to be able to compile your code.

This normally involves running `module load` command but the exact command depend on the system being used.

In addition, an appropriate acrg environment should be set up using the acrg_environment.yml file 
contained within the repository.

These are some of the specific modules needed if running on Bristol systems

*Bristol - Blue Crystal 4 (15/07/2020) *

Anaconda:
 $ module load languages/anaconda3/3.7
Intel (if using intelem as fcompiler):
 $ module load intel/2017.01


*Bristol - Blue Pebble 1 (15/07/2020) *

Anaconda:
 $ module load lang/intel-parallel-studio-xe/2019.u3
Intel (if using intelem as fcompiler):
 $ module load lang/python/anaconda/3.7-2019.10


*Other*

For other systems use whichever modules are available (anaconda using python 3 and intel if using).
To check all available modules you can often run:
$ module avail

*ALL - ACRG environment*

To create (if needed) and activate the acrg environment these lines should be run *after* the 
appropriate modules have been loaded:
 $ conda env create -f ${ACRG_PATH}/acrg_environment.yml
 $ source activate acrg

***************** Compiling uncorrelated version with f2py ****************************************

gfortran
^^^^^^^^

One way to compile for is to use `gfortran`. The file tdmcmc_uncorr.so is created with::

 $ f2py -c -m tdmcmc_uncorr acrg_hbtdmcmc_uncorr.f90

For parallel tempering, to use the gfortran compiler for multiple processors using OMP do::

 $ f2py -c -m tdmcmc_uncorr_pt --f90flags='-fopenmp' -lgomp acrg_hbtdmcmc_uncorr.f90

intelem
^^^^^^^

However, if possible, it is better to compile with the `intel` compiler as this should increase the speed of the run.

The file tdmcmc_uncorr.so is created with::

 $ f2py -c -m tdmcmc_uncorr --fcompiler=intelem acrg_hbtdmcmc_uncorr.f90

For parallel tempering (for python 3.*) run the following::

 $ f2py -c -m tdmcmc_uncorr_pt --f90flags='-qopenmp' -liomp5 --fcompiler=intelem acrg_hbtdmcmc_uncorr.f90


***************** Compiling evencorr hierarchical version with f2py ****************************************

See the previous section, just include different input and output filenames e.g. to compile for single thread using the intel compiler::

 $ f2py -c -m tdmcmc_evencorr --fcompiler=intelem acrg_hbtdmcmc_evencorr.f90


***************** Compiling correlated hierarchical version with f2py ****************************************


$ f2py -L/usr/lib64 -llapack -c -m tdmcmc_corr_s acrg_hbtdmcmc_corr.f90

$ f2py -L/usr/lib64 -llapack -c -m tdmcmc_corr_pt --f90flags='-fopenmp' -lgomp acrg_hbtdmcmc_corr.f90


***************** Compiling correlated hierarchical version with ifort ****************************************

For the compilation of the correlated code in this way there are extra intel libraries that need to be linked.
How to link to these libraries can depend on the system being used.


Note: This should be faster than gfortran, my (Mark's) tests have shown run time improvements of 
~2x to 4x faster.

*Bristol - snowy*

For the Bristol snowy server use this version:

 $ f2py -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -c -m tdmcmc_corr_pt --fcompiler=intelem --f90flags='-fast -qopenmp' -liomp5 acrg_hbtdmcmc_corr.f90

Make sure you've got the following line in your .bashrc to set up the correct intel environments:

source /opt/intel/bin/compilervars.sh intel64


*Bristol - Blue Crystal 4*

 $ f2py -lmkl_avx2 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -c -m tdmcmc_corr_pt 
--fcompiler=intelem --f90flags='-fast -qopenmp' -liomp5 acrg_hbtdmcmc_corr.f90


*Bristol - Blue Pebble 1*

 $ f2py -lmkl_avx512 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -c -m tdmcmc_corr_pt 
 --fcompiler=intelem --f90flags='-fast -qopenmp' -liomp5 acrg_hbtdmcmc_corr.f90


*Other*

Note: the main difference between the BC4 and BP1 compilations is the -lmkl_av* flags needed. 
Either of these flags may work on a new system but adding a -lmkl_def flag instead should be a 
catch all term if these aren't appropriate (this should use the default):

 $ f2py -lmkl_def -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -c -m tdmcmc_corr_pt 
 --fcompiler=intelem --f90flags='-fast -qopenmp' -liomp5 acrg_hbtdmcmc_corr.f90


** Previous compilations using this method ** 

This was the previous snowy compilation in python 2.7::

 $ f2py -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -c -m tdmcmc_corr_pt --fcompiler=intelem --f90flags='-fast -openmp' -liomp5 acrg_hbtdmcmc_corr.f90


To get this to work on air you'll need to copy the following into a terminal, or put in your bashrc:

export LD_PRELOAD=/opt/intel/mkl/lib/intel64/libmkl_core.so:/opt/intel/mkl/lib/intel64/libmkl_sequential.so

For some reason it works on snowy without this. Not sure about non-Bristol servers...


********************************* Running from python *************************************************

import module_name

Where module_name is the name of module_name.so

When calling the trasdimensional mcmc function follow what's in the template file, but if using a different module name make sure you change that.

************** If using openMP then you have to run from ipython not Spyder *********************************




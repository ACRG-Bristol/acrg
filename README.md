# __ACRG Repository __
Shared Python code for the Atmospheric Chemistry Research Group (University of Bristol). Dealing with observations, NAME output and inversion code.

## __Getting Started__

### Prerequisites

- python 3
- git 
- anaconda (recommended)
- fortran compiler (optional)

The ACRG code is written primarily using python 3. It is recommended to use Anaconda and an environment file is provided 
to setup a compatible python environment. The TDMCMC section of the repository contains FORTRAN code that must be compiled by the user to be used. Git is required both for obtaining the code and for stamping code outputs with version numbers.

### Installing

1. Clone the repository
2. Setup paths:
    - copy acrg_config/templates/paths_default.yaml to acrg_config/paths.yaml
    - Overwrite the default folder paths in acrg_config/paths.yaml with your system specific values

3. Setup the python environment: 
    - conda env create -f acrg_environment.yml
    - conda activate acrg
    - Note: creating the environment may take several minutes

4. The ACRG repository is ready to use! Start by running the tests (see below)

### Running tests

To ensure the code has been setup properly, and to aid in development and maintenance, the repository comes with a suite of tests to be run. These should be run when you first install the code to ensure it is setup correctly, and frequently during code development work to prevent adding bugs to the code.

From the /tests/ folder run the terminal command 'pytest'. The complete test suite may take a few minutes to run. Run with 'pytest -v' for more information if any tests fail.

## __Credits__
The code is developed and maintained by the ACRG Modelling team (http://www.bristol.ac.uk/chemistry/research/acrg/)

## __License__
The code currently has no license. This should change to enable sharing and further users. 


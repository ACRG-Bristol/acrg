Read me file for Name III

Components:
----------

Name III proper - Main model.

Preprocessor - Preprocesses some of the Name III fortran before compilation.

Lincom - Linear flow model (optionally invoked by Name III as separate executable).

ECMWF GRIB API software - software for reading and writing GRIB files (optionally 
                          linked into the Name III executable).

Name to Adms - Converts certain Name III output files to Adms met input files.

IDL Graphics - Plots results.

Python - Python utilities for processing and plotting output results.

Fortran to Html Converter - Converts fortran to Html.

ReadMe.txt

NameIII.dsw } Files for Microsoft Developer Studio / Compaq Visual Fortran.
*.dsp       }

NameIII.sln } Files for Microsoft Visual Studio / Intel Visual Fortran.
*.vfproj    }


Folder Structure:
----------------

The source code, together with any compilation scripts and change information is kept 
in folders beginning with 'Code_'.

Separate 'Code_' folders are used for Name III, Preprocessor, Lincom, Name to Adms,
IDL Graphics, Grib API, Python and Fortran to Html.

Name III, Preprocessor, Lincom and Name to Adms are compiled to produce executables.

GRIB API is compiled to produce a library which is optionally linked into the NAME III 
executable.

IDL Graphics, Python and Fortran to Html are not compiled but interpreted at run time.
 
When an executable is produced it is copied to ...\Executables_Win or 
.../Executables_Linux. It can be executed from there or from its original location. 
Alternatively it can be copied somewhere else and executed from there. 

...\NotesOnUsedRoutines contains notes on code which is used but which is not included 
in the Name III Library.

...\Runs and ...\Examples contain various example scripts and input files and are 
convenient folders for such files, but any folder can be used. ...\Resources contains
various resources such as met data for use in runs.


Compilers:
---------

Compilers used are Compaq Visual Fortran Standard Edition 6.6-1684-47B6E (on Windows) 
and Intel Fortran Compiler XE for Linux, Version 12.0.4.191 (on Linux).

On Windows machines without the compiler installed, the file FQWIN.HLP is needed as 
well as the exe file.


Compiler directives for fpp for use with Name III, Preprocessor and Name to Adms:
--------------------------------------------------------------------------------

The following variables can be defined to control the compilation.

ExtraChecks: Used to specify that extra checks are made in the code, leading to
             increased reliability and bug identification at the expense of 
             increased run time.

CompaqPCCompiler: Used when the Compaq Fortran complier is used on PCs.

IntelLinCompiler: Used when the Intel Fortran compiler is used on Linux PCs.

StdIs64Bit: Used to specify that precision for standard reals is 64-bit.

PosIs64Bit: Used to specify that precision for position coords is 64-bit (if not set,
            precision is the same as for standard reals).

SABuilding: Used to configure the FlowAndFlowProfile Module for use in the Stand-Alone
            Building Model. Avoids use of Time Module and removes some non-neutral
            stability calculations.

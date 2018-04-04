Read me file for Name III

Components: 
----------  

Global Parameters
Error And Message 
String
System
Screen
Unit        
Time
Maths  
Physics           
Headed File
Error And Message II 
Coordinate System
Grid And Domain
Sort
Service
Flow And Flow Profile
Single Site
Building
Common Met            
Prototype Met   } Met Modules
Single Site Met }
NWP Met         }
Radar Met       }
Ancillary Met   }
Mets
Common Flow            
Prototype Flow   } Flow Modules
Single Site Flow }
NWP Flow         }
Building Flow    }
Radar Flow       }
Lincom Flow      }
Flows
Particle Size Distribution
Species                         
Source            
Fluctuation
Output
Plume Rise
Particle          
Puff
ChemistryScheme     
Chemistry
Case
Restart
Input       
Timers
OpenMP
IOThread
CompilerInfo
PhysicalUnits

Main Name III

LinuxIntelDebug           } Scripts for compiling code on Linux with Intel Fortran compiler.
LinuxIntelRelease         } LinuxIntelRelease_GribAPI includes the ECMWF GRIB API software.
LinuxIntelRelease_GribAPI }

Makefile        } Make file for compiling code on Linux PCs with Intel compiler using 
                } the Linux utility make.
makedepend.pl   } Perl script for automatically generating module dependencies
                } for make
NameIII_sources } List of source files, needed by makedepend.pl

CatCode     } Linux and MS-DOS scripts for concatenating the Fortran code.
CatCode.bat }

Changes.txt } Changes file.

ReadMe.txt

Note Name III can also invoke { Lincom (as separate exe file)
                              { ECMWF GRIB API (as linked library).


Compilation: 
-----------  

Name III is compiled to produce an executable which is copied to ...\Executables_Win 
or .../Executables_Linux.

Preprocessor.exe is used to preprocess the Mets and Flows modules, converting the 
.P90 files to .F90 files. All the .F90 files are preprocessed with the Fortran 
pre-processor (fpp) before compilation.


1) On Windows PCs with Compaq compiler:

The compilation is performed in the developer studio and the compiler options are set 
up there. 

In compiling Name III, the preprocessing of the Mets and Flows modules by 
Preprocessor.exe is part of the build process in the developer studio - custom build 
steps are defined for Mets.P90 and Flows.P90 with outputs Mets.F90 and Flows.F90 (note 
no dependencies need to be specified for these steps as there are no other files, such 
as include files, involved). 

Post build steps are defined to automatically copy the executables to 
...\Executables_Win.

Name III is compiled as a quick-win application and preprocessing by fpp is invoked by 
using compiler option /fpp. The Compaq version of fpp doesn't by default do macro 
expansion outside of pre-processor lines (which is a good thing). The output of the 
preprocessing stage can be kept with the option /keep.

Various fpp symbols are defined by using the compiler option /define:SymbolName.
Symbols defined are: 
   Debug version:   CompaqPCCompiler, UseConvert, ExtraChecks  
   Release version: CompaqPCCompiler, UseConvert

The stack size is increased from the default using the linker option /stack: 
   Debug version:   /stack:0x10000000 
   Release version: /stack:0x10000000
(0x10000000 is 256 MB in hexadecimal). 
 
No other compiler options are used (other than defaults as set by the developer 
studio).
 

1a) Intel Fortran on Windows:

Various fpp symbols are defined by using the compiler option /define:SymbolName.
Symbols defined are: 
   Debug version:   IntelWin, UseConvert, ExtraChecks  
   Release version: IntelWin, UseConvert


2) Intel Fortran on Linux PCs:

The compilation is performed by a compilation script and the compiler options are set 
up there.  

The script automatically copies the executables to ...\Executables_Linux.

The preprocessing of the Mets and Flows modules by Preprocessor.exe is carried out by
the script. 

The code is preprocessed with the Fortran pre-processor (fpp) - this is automatic 
provided the source code ends in .F90 (as opposed to .f90). The option -Wp,-macro=no 
is needed to ensure fpp doesn't do macro expansion outside of pre-processor lines (the 
-Wp bit passes what follows to fpp; the -macro=no bit seems undocumented by Intel but 
is what worked with the NAG compiler). The output of the pre-processing stage can be 
obtained by using the option -F or -P (but without any options other than those 
related to fpp - no compilation is performed). 

Various fpp symbols are defined by using the compiler option -DSymbolName.
Symbols defined are:
   LinuxIntelDebug:   IntelLinCompiler, UseConvert, ExtraChecks 
   LinuxIntelRelease: IntelLinCompiler, UseConvert

For other compiler options see the scripts.

The makefile Makefile can also be invoked by typing "make".


3) Intel Fortran on Macs:

Not guaranteed to work without some effort, but instructions for Intel Fortran on 
Linux PCs should work.


4) Sun compilers on Sun or Intel machines:

Not guaranteed to work without some effort, but the following instructions should 
work.

Use makefile Makefile (after editing it). 
The fpp symbol 'sun' is defined automatically by the compiler.


Compiling (parallel) NAME with make
-----------------------------------

The standard compile scripts will only compile the code in serial mode. 

To compile the parallel version of the code use 'make all' ('make' or 'make nameiii' 
will also work; 'make clean' removes object files from previous compilations). The 
Makefile automatically preprocesses the code, resolves all dependencies, determines 
the correct architecture and constructs the correct name of the executable which is 
moved to ../Executables_.... The name of the default executable is 'nameiii.exe' but 
the following labels are added to the base of the filename:

'GribAPI', if compiled with USEGRIBAPI=true
'Debug', if compiled with COMPILERMODE=Debugging
'_64bit', if compiled on a 64-bit system (automatically detected)
'_par', if compiled in parallel mode (MODEL=parallel)

The following options can be passed to make, either my modifying the Makefile
or by specifying them in the form OPTION=<value> as a command line parameter:

MAKE: name of the make executable
COMPILEROPTIONS: the compiler and system to be used, so far only IntelLinux,
IntelMac, IntelSun and IBMUnix are supported
COMPILEROPTIONSFILE: If COMPILEROPTIONS is set to Custom, all options are
loaded from this file
COMPILERMODE: can be Optimized or Debugging
PROFILER: If this is set to gprof, the code is compiled with support for the
gprof profiler
MODEL: if this is set to serial, the code will be compiled in serial mode, if
it is parallel (default), the code will be compiled with OpenMP support.
USETIMERS: Set this to true to use the timer module. Note that this is
currently only supported together with MODEL=parallel
USEGRIBAPI: Set this to true to use the ECMWF GRIB API
OPENMPVERSION: Selects the version of OpenMP to use (available options are
2.5, 3.0, or ifortautodetect, which only works with the Intel compiler)

IT IS ESSENTIAL to run 'make clean' when changing one of the above options or 
when moving to a different architecture.

To use the ifort 11.0 compiler on 64-bit systems, run the following command:
. /opt/intel/Compiler/11.0/083/bin/ifortvars.sh intel64. Note that the
standard compiler (intel 9.0) does not link the OpenMP library properly and should 
not be used.

The following compiler options are used in the Makefile:

Debugging mode: -Wp,-macro=no -C -w -extend_source -nbs -Vaxlib -auto 
-DIntelLinCompiler -DExtraChecks -heap-arrays -sox -traceback -g

Optimized mode: -Wp,-macro=no -O3 -w -extend_source -nbs -Vaxlib -auto 
-DIntelLinCompiler -heap-arrays -sox

If compiled with OpenMP: -openmp (and optionally -DUseOpenMP3_0)
If compiled with v11 of Intel compiler (or later): -DUseConvert
If compiled with the timer module: -DUseTimers

-V             print compiler version information
-Wp,-macro=no  pass -macro=no to (fortran) preprocessor, whichs turns off macro 
               expansion during preprocessing
-C             Enable runtime checks
-O3            Optimization level 3
-w             Disable all warning messages
-extend_source Allow up to 132 columns in fortran source code
-nbs           Treat backslash as normal parameter
-Vaxlib        Link with portability library (needed for SYSTEM commands and +U77). 
               Probably not needed for ifort 10 and newer
-auto          Causes all local, non-SAVEd variables to be allocated to the run-time 
               stack.
-openmp        Compile with OpenMP support
-DIntelLinCompiler  Define IntelLinCompiler (for use in #ifdef)
-DExtraChecks       Define ExtraChecks (for use in #ifdef)
-DUseTimers         Use the timer module, so far this is only supported in 
                    combination with -openmp.
-DUseOpenMP3_0  Use extra parallelisation functionality available at OpenMP 3.0
                (default is to use routines supported at version 2.5).
-DUseConvert    Use convert statement for endian conversion when reading unformatted
                binary files.
-heap-arrays 	Move some of the processing onto the heap. This was necessary due to a 
                memory fault occurring in the parallel version on 20/01/2010 which
                is probably related to the issue that Helen W found on the 32bit 
                machines.



Running:
-------
 
When running Name III, there is full control over the locations and (at least aspects 
of) the names of most files and folders. The files and folders can be described with 
full paths or from the current folder. This is done through the specification of the 
executable file, the command line arguments, user input in response to queries, and 
file and folder references in the input files. There are a few files where there is no 
control over their location and/or name (e.g. Lincom input/output). The current folder 
is used for these. 

Name III produces some output in quick-win windows on Windows PCs with the Compaq 
compiler. This includes output to unit 6. When run on other compilers the unit 6 
output is treated as usual, some of the rest of this output (namely the text output) 
is written to compiler determined files (usually of the form 'fort.n'), and the 
remainder (namely the graphical output) is lost.

If Lincom is invoked, the unit 6 output from Lincom is piped to LincomOut.txt in the 
current folder.

When executed from within the developer studio, the working directory is set to 
...\Runs and the immediately triggered executable is in ...\Debug or ...\Release not 
...\Executables_Win.

When running on Linux the following commands are needed
export F_UFMTENDIAN="big;little:100"  (if not using the -DUseConvert compiler option)
ulimit -s unlimited                   (this command sets the stack size)

Running parallel NAME:
----------------------

The standard compiler on 64-bit systems is intel fortran 9.0. This means that 
even if the code has been compiled with intel fortran 11.0, it will not find the 
(dynamically linked) OpenMP library in the fortran 9.0 library path. To set the 
correct library path either set it manually or run the following command to
set up the correct environment variables for ifort 11.0:

. /opt/intel/Compiler/11.0/083/bin/ifortvars.sh intel64

processes can be bound to individual cores by setting the environment variable 
KMP_AFFINITY:

export KMP_AFFINITY=verbose,compact
see http://www.intel.com/software/products/compilers/docs/flin/main_for/mergedprojects/optaps_for/common/optaps_openmp_thread_affinity.htm

It is recommended (and sometimes essential) to set the OMP_STACKSIZE environment
variable when running with multiple threads under OpenMP, e.g.,

 export OMP_STACKSIZE='4m'  ... on Linux systems
 setenv OMP_STACKSIZE='4m'  ... on Windows systems
 
The OMP_STACKSIZE variable defines the stack size that is allocated for any
additional threads that are created. The default, if not set explicitly, is
4 MB ('4m') and this is usually sufficient when running NAME in parallel.
However it may sometimes be necessary to set a larger value than this.
The NAME run will fail with a 'memory fault' if the stack size limit is
insufficient for the run.

OpenMP references:

[1] Introduction to OpenMP tutorial from Lawrence Livermore National Lab: 
https://computing.llnl.gov/tutorials/openMP/
[2] "An Introduction to OpenMP" by Ruud van der Pas, IWOMP 2005: 
http://www.nic.uoregon.edu/iwomp2005/iwomp2005_tutorial_openmp_rvdp.pdf
[3] "Parallel Programming in Fortran 95 using OpenMP" by Miguel Hermanns: 
http://www.openmp.org/presentations/miguel/F95_OpenMPv1_v2.pdf

Command line arguments:
----------------------

These are case insensitive with the exception of file names and parts of file names
which are case sensitive when the operating system is.

-InputFile=???? 
   Here ???? is the primary input file.
   If this is the first command line option the '-InputFile=' can be omitted.
   If not specified the user is asked for the primary input file name.

-LogFolder=???? 
   Here ???? is the log folder (for keeping the log file in).
   If not specified the log folder is the folder containing the primary input file.

-RestartFolder=????
   Here ???? is the restart folder (for keeping the restart files in).
   If not specified the restart folder is the folder containing the primary input 
   file. 
   This can be specified whether or not the run is restartable (although it will have
   no effect if the run is not restartable).

-Restart=????
   If present this indicates that the run is to restart from a restart file. Otherwise
   the run starts from the beginning.
   Here ???? is the part of the restart file name between 'Restart_' and '.dat' (this
   part characterises the point in the run at which the restart file was created).
   If =???? is omitted the latest restart file is used.
   If starting from the beginning (i.e. -Restart absent) then the 'restart catalogue 
   validity files' must not be present. If restarting from the latest restart file 
   (i.e. -Restart present but =???? absent) and there are no restart files, the run 
   starts at the beginning.

-ClosePromptLevel=? 
   Here ? indicates the error level at or above which the user is prompted to 
   close the window. Below this value the window closes automatically. The error
   levels are as follows:
       0 indicates no errors or warnings,
       1 indicates some warnings but no errors, 
       2 indicates some non-fatal errors but no fatal errors,
       3 indicates fatal errors but no unexpected fatal errors,
       4 indicates unexpected fatal errors.
   Unexpected errors are those which require the code or system to be changed to 
   correct (i.e. not due to errors in the user input).
   Note this does not apply to errors which occur before the command line 
   processing is complete which always produce the standard quick win termination 
   message box or to untrapped fatal errors which are handled directly by the 
   compiler.
   If not specified, the value is 0.


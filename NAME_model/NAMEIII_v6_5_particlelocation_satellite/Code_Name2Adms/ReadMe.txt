Read me file for Name to Adms

Components: 
----------  

Global Parameters    } identical to components of Name III
ErrorAndMessage      } 
String               }
System               }
Screen               }
Unit                 }    
Time                 }
Maths                }
Physics              }          
Headed File          }
Error And Message II }
Coordinate System    }
Grid And Domain      }
Sort                 }
Physical Units       }
Service              }
Main Name to Adms 

LinuxIntelDebug   } Script for compiling code on Linux PCs with Intel compiler.
LinuxIntelRelease } 

Name2Adms.dsp    } Project files for use with Compaq developer studio on Windows PCs.
Name2Adms.plg    }

Name2Adms.vfproj } Project file for use with Visual Studio 2008 on Windows PCs.

Changes.txt } Changes file.

ReadMe.txt  


Compilation: 
-----------  

Name to Adms is compiled to produce an executable which is copied to 
...\Executables_Win or .../Executables_Linux.

All the .F90 files are preprocessed with the Fortran pre-processor (fpp) before 
compilation.


1) On Windows PCs with Compaq compiler:

The compilation is performed in the developer studio and the compiler options are set 
up there. 

Post build steps are defined to automatically copy the executables to 
...\Executables_Win.

Name to Adms is compiled as a quick-win application and preprocessing by fpp is 
invoked by using compiler option /fpp. The Compaq version of fpp doesn't by default do 
macro expansion outside of pre-processor lines (which is a good thing). The output of 
the preprocessing stage can be kept with the option /keep.

Various fpp symbols are defined by using the compiler option /define:SymbolName.
Symbols defined are: 
   Debug version:   CompaqPCCompiler, ExtraChecks  
   Release version: CompaqPCCompiler, ExtraChecks

No other compiler options are used (other than defaults as set by the developer 
studio).
 

2) Intel Fortran on Linux PCs:

The compilation is performed by a compilation script and the compiler options are set 
up there.  

The script automatically copies the executables to ...\Executables_Linux.

The code is preprocessed with the Fortran pre-processor (fpp) - this is automatic 
provided the source code ends in .F90 (as opposed to .f90). The option -Wp,-macro=no 
is needed to ensure fpp doesn't do macro expansion outside of pre-processor lines (the 
-Wp bit passes what follows to fpp; the -macro=no bit seems undocumented by Intel but 
is what worked with the NAG compiler). The output of the pre-processing stage can be 
obtained by using the option -F or -P (but without any options other than those 
related to fpp - no compilation is performed). 

Various fpp symbols are defined by using the compiler option -DSymbolName.
Symbols defined are:
   LinuxIntelDebug:   IntelLinCompiler, Linux64Intel, ExtraChecks 
   LinuxIntelRelease: IntelLinCompiler, Linux64Intel, ExtraChecks 

For other compiler options see the scripts.


Running:
-------

When executed from within the developer studio, the working directory is set to 
...\Runs and the immediately triggered executable is in ...\Debug or ...\Release not 
...\Executables_Win.

 
Command line arguments:
----------------------

First argument: input file name 

Second argument: output file name

# Makefile for NAME

#######################################################################
#
#  Makefile for Name III
#
#######################################################################
#
#  Usage:
#
#     'make', 'make all' or 'make nameiii' will compile name and 
#     copy the executable to the executables directory.
#   
#     Compiler options can be specified by umcommenting the appropriate option
#     in the OPTIONS section below.
#     Alternatively these options can be set as command line options, 
#     in which case they will override the default options below. 
#     For example
#
#        make COMPILEROPTIONS=IntelLinux COMPILERMODE=Debugging PROFILER=gprof
#
#     will compile with the options for the intel compiler on a linux
#     architecture and the debugging and gprof options for the compiler.
#
#        make COMPILEROPTIONS=IntelLinux COMPILERMODE=Optimized
#
#     will compile with options for the intel compiler on Linux,
#     using optimized compiler options.
#
#     During development work, Compiler flags can also be specified in
#     the separate file COMPILEROPTIONSFILE (Default: Makefile_cflags).
#     This file is loaded if COMPILEROPTIONS equals Custom.
#
#
#  Eike Mueller, Atmospheric Dispersion Group, MetOffice 26 Feb 2010 
#
#######################################################################



#######################################
# OPTIONS BEGIN
#######################################

# Compiler executable

# make executable
MAKE=gmake

# Compiler options

#COMPILEROPTIONS=IBMUnix
COMPILEROPTIONS=IntelLinux
#COMPILEROPTIONS=IntelMac
#COMPILEROPTIONS=IntelSun
#COMPILEROPTIONS=IntelCLE
#COMPILEROPTIONS=CrayCLE
#COMPILEROPTIONS=Custom

# File with custom compiler flags

COMPILEROPTIONSFILE=Makefile_cflags 

# Compiler mode

#COMPILERMODE=Debugging
COMPILERMODE=Optimized

# Profiler
#PROFILER=gprof
PROFILER=

# vTune
#VTUNE=true
VTUNE=false

# Parallel/Serial implementation
#MODEL=serial
MODEL=parallel

# Use the timer Module? (Currently only supported for MODEL=parallel)
USETIMERS=true
#USETIMERS=false

# Compile with GRIB_API support?
#USEGRIBAPI=true
USEGRIBAPI=false

# Compile with NetCDF support?
#USENETCDF=true
USENETCDF=false

# Which OpenMP version do we use?
#
# 2.5             : use OpenMP 2.5
# 3.0             : use OpenMP 3.0
# ifortautodetect : this automatically detects the version of the 
#                   Intel Fortran compiler that is used and then set
#                   uses OpenMP 3.0 for version 11 and OpenMP 2.5 for
#                   all other version. Only works together with the Intel
#                   compiler.
OPENMPVERSION=ifortautodetect


#######################################
# OPTIONS END
#######################################

#######################################
# FILE EXTENSIONS
#######################################

F90_EXT=F90
P90_EXT=P90
OBJ_EXT=o
EXE_EXT=exe

# XL Unix Compiler
# ----------------
ifeq ($(COMPILEROPTIONS),IBMUnix)
  FC=xlf90_r
# Note: Macro expansion NOT suppressed
# Display version information?
# Debugging options
  ifeq ($(COMPILERMODE),Debugging)
    CFLAGS = -WF,-DIBMUnixCompiler -WF,-DExtraChecks \
			  		 -qnoescape -qcheck=all -w -g -qfree -qreport
  endif
  ifeq ($(COMPILERMODE),Optimized)
    CFLAGS = -WF,-DIBMUnixCompiler -WF,-DExtraChecks \
			  		 -O3 -qnoescape -w -qfree -qreport
	endif
  LFLAGS=$(CFLAGS)
endif

# Intel Linux Compiler
# --------------------
ifeq ($(COMPILEROPTIONS),IntelLinux)
# Intel compiler on Linux and Mac
  FC=ifort
# Debugging options
  ifeq ($(COMPILERMODE),Debugging)
    CFLAGS = -Wp,-macro=no -C -check noarg_temp_created -w -extend_source -nbs -Vaxlib -auto \
             -DIntelLinCompiler -DExtraChecks -heap-arrays -sox -traceback                   \
						 -g
	endif
# Optimized options
  ifeq ($(COMPILERMODE),Optimized)
    CFLAGS = -Wp,-macro=no -O3 -w -extend_source -nbs -Vaxlib -auto \
             -DIntelLinCompiler -heap-arrays -sox
  endif
  LFLAGS=$(CFLAGS)
endif

# Intel CLE Compiler
# --------------------
# -std95        add to see non-standard code constructs - useful to refactor code at some point in future!
# -Vaxlib       removed as this is now deprecated
# -heap-arrays  removed
# -fp-model precise -dynamic -shared-intel -mcmodel=medium are options used with the UM

ifeq ($(COMPILEROPTIONS),IntelCLE)
# Intel compiler on Cray Linux Environment
  FC=ftn
# Debugging options
  ifeq ($(COMPILERMODE),Debugging)
    CFLAGS = -Wp,-macro=no -C -check noarg_temp_created -extend_source -nbs -auto -no-vec \
             -DIntelLinCompiler -DExtraChecks -sox -traceback -g                          \
             -fp-model precise -dynamic -shared-intel -mcmodel=medium
  endif
# Optimized options
  ifeq ($(COMPILERMODE),Optimized)
    CFLAGS = -Wp,-macro=no -O3 -w -extend_source -nbs -auto           \
             -DIntelLinCompiler -sox                                  \
	     -fp-model precise -dynamic -shared-intel -mcmodel=medium
  endif
  LFLAGS=$(CFLAGS)
endif

# Cray CLE Compiler
# -----------------
ifeq ($(COMPILEROPTIONS),CrayCLE)
# Cray compiler on Cray Linux Environment
  FC=ftn
# Debugging options
  ifeq ($(COMPILERMODE),Debugging)
    CFLAGS = -ea -em -en -dF -g -R bc -N 132 -DCrayCLECompiler -DExtraChecks
  endif
# Optimized options
  ifeq ($(COMPILERMODE),Optimized)
    CFLAGS = -ea -em -dF -O3 -N 132 -DCrayCLECompiler
  endif
  LFLAGS=$(CFLAGS)
endif

# Intel Mac Compiler
# ------------------
ifeq ($(COMPILEROPTIONS),IntelMac)
# Intel compiler on Linux and Mac
  FC=ifort
# Debugging options
  ifeq ($(COMPILERMODE),Debugging)
    CFLAGS = -Wp,-macro=no -DIntelLinCompiler -DExtraChecks \
             -m64 -extend-source -g -C -ftrapuv -fpe0 -convert big_endian -nbs
  endif
# Optimized options
	ifeq ($(COMPILERMODE),Optimized)
    CFLAGS = -Wp,-macro=no -DIntelLinCompiler \
             -fast -no-prec-div -m64 -fpe0 -convert big_endian -nbs
  endif    
  LFLAGS := $(CFLAGS)
endif

# SUN Intel hardware and compiler
# -------------------------------
# Note: preprocessor relies on 'sun' symbol being defined.

ifeq ($(COMPILEROPTIONS),IntelSun)
  FC=f95
  ifeq ($(COMPILERMODE),Debugging)
# Debugging options
    CFLAGS = -g -fnonstd -DExtraChecks -e -C -xfilebyteorder=big16:%all
  endif
	ifeq ($(COMPILERMODE),Optimized)
# Optimized options
    CFLAGS = -fast -fnonstd -e -xfilebyteorder=big16:%all
  endif    
  LFLAGS := $(CFLAGS)
endif

ifeq ($(COMPILEROPTIONS),Custom)
  -include $(COMPILEROPTIONSFILE)
endif

# Profile with gprof
# -------------------------------
ifeq ($(PROFILER),gprof)
  CFLAGS += -pg
  LFLAGS += -pg
endif

# Optimized with vtune
# -------------------------------
ifeq ($(VTUNE),true)
  CFLAGS += -g
  LFLAGS += -g
endif

# OpenMP version
# -------------------------------
ifeq ($(OPENMPVERSION), ifortautodetect)
  ifeq ($(FC), ifort)
    $(shell ifort -V 2> ifortversion.tmp)
    IFORTVERSION=$(shell grep Version ifortversion.tmp | sed 's/.*Version \([0-9]*\).*/\1/')
    USEOPENMPTHREE := $(shell test $(IFORTVERSION) -ge 11 && echo 1)
    ifeq ($(USEOPENMPTHREE), 1)
      # --- Currently an issue with using OpenMP 3.0 with v12.0 of the Intel compiler.
      #OPENMPVERSION=3.0
      OPENMPVERSION=2.5
    else
      OPENMPVERSION=2.5
    endif
  else
    OPENMPVERSION=2.5
  endif
endif

ifeq ($(OPENMPVERSION), 3.0)
  CFLAGS += -DUseOpenMP3_0
endif

# Use 'convert' option when opening unformatted binary data files
# -------------------------------
USECONVERT=false
ifeq ($(FC), ifort)
  $(shell ifort -V 2> ifortversion.tmp)
  IFORTVERSION=$(shell grep Version ifortversion.tmp | sed 's/.*Version \([0-9]*\).*/\1/')
  IFORTUSESCONVERT := $(shell test $(IFORTVERSION) -ge 11 && echo 1)
  ifeq ($(IFORTUSESCONVERT), 1)
    USECONVERT=true
  endif
endif
ifeq ($(FC), ftn)
  USECONVERT=true
endif

ifeq ($(USECONVERT),true)
  CFLAGS += -DUseConvert
endif

#######################################
# ARCHITECTURE DEPENDENT SETTINGS
#######################################

SYSTEM=$(shell uname)
ifeq ($(SYSTEM), Linux)
  ARCHNAME=$(shell uname -p)
  ifeq ($(ARCHNAME),x86_64)
    ARCHITECTURE=64bit
  else
    ARCHITECTURE=32bit
  endif
  PREPROCESSORDIR=../Executables_Linux/
  EXEDIR=../Executables_Linux/
  SHAREDLIBDIR=../SharedLibraries_Linux/
endif

ifeq ($(SYSTEM), AIX)
  ARCHITECTURE=64bit
  PREPROCESSORDIR=../Executables_Unix/
  EXEDIR=../Executables_Unix/
endif

# Parallelisation with OpenMP
# -------------------------------
# Note that Cray compiler on CLE has a default of OpenMP being on, so treat this case separately.

ifeq ($(COMPILEROPTIONS),CrayCLE)
  # Cray compiler on Cray Linux Environment
  ifeq ($(MODEL),serial)
    CFLAGS += -xomp
    LFLAGS += -xomp
  endif
else
  ifeq ($(MODEL),parallel)
    ifeq ($(SYSTEM), Linux)  
      CFLAGS += -openmp
      LFLAGS += -openmp
    endif
    ifeq ($(SYSTEM), AIX)  
      CFLAGS += -qsmp=omp
      LFLAGS += -qsmp=omp
    endif 
  endif
endif

# Use Timer Module
# -------------------------------
ifeq ($(USETIMERS),true)
  ifeq ($(SYSTEM), Linux)  
    CFLAGS += -DUseTimers
    LFLAGS += -DUseTimers 
  endif
  ifeq ($(SYSTEM), AIX)  
    CFLAGS += -WF,-DUseTimers=UseTimers
    LFLAGS += -WF,-DUseTimers=UseTimers
  endif
endif


# Construct name of executable and preprocessor

ifeq ($(COMPILERMODE),Debugging)
  EXL_DEBUG=Debug
else
  EXL_DEBUG=
endif

ifeq ($(ARCHITECTURE),64bit)
	EXL_ARCH=_64bit
else
	EXL_ARCH=
endif

ifeq ($(MODEL),parallel)
	EXL_PAR=_par
else
	EXL_PAR=
endif

ifeq ($(USEGRIBAPI),true)
	EXL_GRIBAPI=GribAPI
else
	EXL_GRIBAPI=
endif

# Compile with GRIB_API support
# -------------------------------
ifeq ($(USEGRIBAPI),true)
  CFLAGS += -DGRIBsupport -I$(SHAREDLIBDIR)include${EXL_ARCH} -fPIC
  LFLAGS += -L$(SHAREDLIBDIR)lib -lgrib_api_f90$(EXL_ARCH) -lgrib_api$(EXL_ARCH)
endif

ifeq ($(USENETCDF),true)
        EXL_NETCDF=NetCDF
else
        EXL_NETCDF=
endif

# Compile with NetCDF support
# -------------------------------
ifeq ($(USENETCDF),true)
  CFLAGS += -DNetCDFsupport -I/usr/local/sci/include/netcdf_fortran/ifort_composerxe
  LFLAGS += -L/usr/local/sci/lib/netcdf_fortran/ifort_composerxe  -lnetcdff
endif


#######################################
# COMPILER INFO
#######################################

ifeq ($(COMPILEROPTIONS),IBMUnix)
  CINFOFLAGS+=-WF,-DCOMPILETIME=\""$(shell date)"\"
  CINFOFLAGS+=-WF,-DCOMPILERVERSION=\""$(shell $(FC) -qversion | sed 's/,//g')"\"
else
  ifeq ($(COMPILEROPTIONS),CrayCLE)
    CINFOFLAGS+=-DCOMPILERVERSION=\""$(shell $(FC) -V 2>&1)"\"
    CINFOFLAGS+=-DCOMPILETIME=\""$(shell date)"\"
  else
    CINFOFLAGS+=-DCOMPILERVERSION=\""$(shell $(FC) --version)"\"
    CINFOFLAGS+=-DCOMPILETIME=\""$(shell date)"\"
  endif
endif


PREPROCESSOR=$(PREPROCESSORDIR)Preprocessor$(EXL_ARCH).exe
EXECUTABLE=$(EXEDIR)nameiii$(EXL_GRIBAPI)$(EXL_NETCDF)$(EXL_DEBUG)$(EXL_ARCH)$(EXL_PAR).$(EXE_EXT)

#######################################
# RULES
#######################################


# Default rule for generating object files, i.e.
# file.$(F90_EXT) is compiled to generate an object file
# file.$(OBJ_EXT)

%.$(OBJ_EXT) : %.$(F90_EXT)
	$(FC) $(CFLAGS) -c $< -o $@

# Preprocessing of P90 files
%.$(F90_EXT) : %.$(P90_EXT)
	$(PREPROCESSOR) $< $@

# Source files, $(SOURCES) the list of all source files defined
# in the file NameIII_sources
include NameIII_sources

# Object files
OBJS=$(SOURCES:$(F90_EXT)=$(OBJ_EXT))

# Construct F90 files which are based on a P90 file. To do this,
# replace F90 -> P90 in every source filename (step 1) and then use
# wildcard to check whether the P90 file exists in the current directory 
# (step 2). Finally replace back P90 -> F90 (step 3)
TMP1:=$(SOURCES:$(F90_EXT)=$(P90_EXT))             # (step 1)
TMP2=$(foreach FILE,$(TMP1),$(wildcard $(FILE)))   # (step 2)
PSOURCES=$(TMP2:$(P90_EXT)=$(F90_EXT))             # (step 3)

#######################################
# TARGETS
#######################################

.PHONY : all
all: nameiii

-include Makefile_dependencies

.PHONY : depend
depend:
	@echo $(SOURCES) > Makefile_sources.tmp
	@$(MAKE) $(PSOURCES)
	perl makedepend.pl
	rm Makefile_sources.tmp

.PHONY : nameiii
nameiii:
	@echo "+-----------------------+"
	@echo "! Makefile for Name III !"
	@echo "+-----------------------+"
	@echo
	@echo "Compiler Options : $(COMPILEROPTIONS)"
	@echo "Compilermode     : $(COMPILERMODE)"
	@echo "Architecture     : $(ARCHNAME)"
	@echo "Preprocessor     : $(PREPROCESSOR)"
	$(MAKE) depend
	$(MAKE) nameiii.$(EXE_EXT)
	cp nameiii.$(EXE_EXT) $(EXECUTABLE)

nameiii.$(EXE_EXT): CompilerInfo.$(OBJ_EXT) $(OBJS)
	 $(FC) -V -o nameiii.$(EXE_EXT) $(LFLAGS) $(OBJS) CompilerInfo.$(OBJ_EXT)

CompilerInfo.$(OBJ_EXT): FORCE
	$(FC) $(CINFOFLAGS) -c CompilerInfo.$(F90_EXT) -o CompilerInfo.$(OBJ_EXT)

FORCE:


# Clean up
.PHONY : clean
clean:
	rm -rf *.$(OBJ_EXT) *.mod *.$(EXE_EXT) ifortversion.tmp Makefile_dependencies


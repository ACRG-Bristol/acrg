# Microsoft Developer Studio Project File - Name="NameIII" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) QuickWin Application" 0x0107

CFG=NameIII - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "NameIII.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "NameIII.mak" CFG="NameIII - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "NameIII - Win32 Release" (based on "Win32 (x86) QuickWin Application")
!MESSAGE "NameIII - Win32 Debug" (based on "Win32 (x86) QuickWin Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "NameIII - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /libs:qwin /nologo /warn:nofileopt
# ADD F90 /compile_only /define:"CompaqPCCompiler" /define:"UseConvert" /fpp /libs:qwin /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /entry:"WinMainCRTStartup" /subsystem:windows /machine:I386 /nodefaultlib:"dfconsol.lib"
# ADD LINK32 kernel32.lib /nologo /stack:0x10000000 /entry:"WinMainCRTStartup" /subsystem:windows /machine:I386 /nodefaultlib:"dfconsol.lib"
# SUBTRACT LINK32 /pdb:none
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=copy .\Release\NameIII.exe ..\Executables_Win\NameIII.exe
# End Special Build Tool

!ELSEIF  "$(CFG)" == "NameIII - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /libs:qwin /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /define:"CompaqPCCompiler" /define:"ExtraChecks" /define:"UseConvert" /fpp /libs:qwin /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /entry:"WinMainCRTStartup" /subsystem:windows /debug /machine:I386 /nodefaultlib:"dfconsol.lib" /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /stack:0x10000000 /entry:"WinMainCRTStartup" /subsystem:windows /incremental:no /debug /machine:I386 /nodefaultlib:"dfconsol.lib" /pdbtype:sept
# SUBTRACT LINK32 /pdb:none
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=copy .\Debug\NameIII.exe ..\Executables_Win\NameIIIDebug.exe
# End Special Build Tool

!ENDIF 

# Begin Target

# Name "NameIII - Win32 Release"
# Name "NameIII - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\AncillaryMet.F90
DEP_F90_ANCIL=\
	".\Debug\CommonMetModule.mod"\
	".\Debug\ServiceModule.mod"\
	
NODEP_F90_ANCIL=\
	".\Debug\grib_api.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Building.F90
DEP_F90_BUILD=\
	".\Debug\FlowAndFlowProfileModule.mod"\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BuildingFlow.F90
DEP_F90_BUILDI=\
	".\Debug\BuildingModule.mod"\
	".\Debug\CommonFlowModule.mod"\
	".\Debug\FlowAndFlowProfileModule.mod"\
	".\Debug\MetsModule.mod"\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Case.F90
DEP_F90_CASE_=\
	".\Debug\ChemistryModule.mod"\
	".\Debug\EulerianInterfaceModule.mod"\
	".\Debug\FlowAndFlowProfileModule.mod"\
	".\Debug\FlowsModule.mod"\
	".\Debug\NWPMetModule.mod"\
	".\Debug\OpenMPModule.mod"\
	".\Debug\OutputModule.mod"\
	".\Debug\ParticleModule.mod"\
	".\Debug\PuffModule.mod"\
	".\Debug\SemiLagrangianModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\SizeDistModule.mod"\
	".\Debug\SourceModule.mod"\
	".\Debug\SpeciesModule.mod"\
	".\Debug\TimerModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Chemistry.F90
DEP_F90_CHEMI=\
	".\Debug\FlowsModule.mod"\
	".\Debug\NameChemistrySchemeModule.mod"\
	".\Debug\OpenMPModule.mod"\
	".\Debug\OutputModule.mod"\
	".\Debug\ParticleModule.mod"\
	".\Debug\PuffModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\SpeciesModule.mod"\
	".\Debug\TimerModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\ChemistryScheme.F90
DEP_F90_CHEMIS=\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\CommonFlow.F90
DEP_F90_COMMO=\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\CommonMet.F90
DEP_F90_COMMON=\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\CompilerInfo.F90
# End Source File
# Begin Source File

SOURCE=.\CoordinateSystem.F90
DEP_F90_COORD=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\GlobalParametersModule.mod"\
	".\Debug\MathsModule.mod"\
	".\Debug\PhysicsModule.mod"\
	".\Debug\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\ErrorAndMessage.F90
DEP_F90_ERROR=\
	".\Debug\GlobalParametersModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\ErrorAndMessageII.F90
DEP_F90_ERRORA=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\GlobalParametersModule.mod"\
	".\Debug\HeadedFileModule.mod"\
	".\Debug\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\EulerianInterface.F90
DEP_F90_EULER=\
	".\Debug\FlowsModule.mod"\
	".\Debug\NWPMetModule.mod"\
	".\Debug\SemiLagrangianModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\SpeciesModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\FlowAndFlowProfile.F90
DEP_F90_FLOWA=\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Flows.F90
DEP_F90_FLOWS=\
	".\Debug\BuildingFlowModule.mod"\
	".\Debug\CommonFlowModule.mod"\
	".\Debug\FlowAndFlowProfileModule.mod"\
	".\Debug\LINCOMFlowModule.mod"\
	".\Debug\MetsModule.mod"\
	".\Debug\NWPFlowModule.mod"\
	".\Debug\PrototypeFlowModule.mod"\
	".\Debug\RadarFlowModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\SingleSiteFlowModule.mod"\
	
NODEP_F90_FLOWS=\
	".\Debug\omp_lib.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Fluctuation.F90
DEP_F90_FLUCT=\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\GlobalParameters.F90
# End Source File
# Begin Source File

SOURCE=.\GridAndDomain.F90
DEP_F90_GRIDA=\
	".\Debug\CoordinateSystemModule.mod"\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\GlobalParametersModule.mod"\
	".\Debug\MathsModule.mod"\
	".\Debug\PhysicsModule.mod"\
	".\Debug\StringModule.mod"\
	".\Debug\TimeModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\HeadedFile.F90
DEP_F90_HEADE=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\GlobalParametersModule.mod"\
	".\Debug\StringModule.mod"\
	".\Debug\TimeModule.mod"\
	".\Debug\UnitModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Input.F90
DEP_F90_INPUT=\
	".\Debug\AncillaryMetModule.mod"\
	".\Debug\CaseModule.mod"\
	".\Debug\ChemistryModule.mod"\
	".\Debug\FlowsModule.mod"\
	".\Debug\NWPMetModule.mod"\
	".\Debug\OpenMPModule.mod"\
	".\Debug\OutputModule.mod"\
	".\Debug\ParticleModule.mod"\
	".\Debug\PuffModule.mod"\
	".\Debug\RestartModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\SizeDistModule.mod"\
	".\Debug\SourceModule.mod"\
	".\Debug\SpeciesModule.mod"\
	".\Debug\TimerModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\IOThread.F90
DEP_F90_IOTHR=\
	".\Debug\CommonMetModule.mod"\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\MetsModule.mod"\
	".\Debug\NWPMetModule.mod"\
	".\Debug\OpenMPModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\TimerModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\LincomFlow.F90
DEP_F90_LINCO=\
	".\Debug\CommonFlowModule.mod"\
	".\Debug\CommonMetModule.mod"\
	".\Debug\FlowAndFlowProfileModule.mod"\
	".\Debug\MetsModule.mod"\
	".\Debug\NWPMetModule.mod"\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\MainNameIII.F90
DEP_F90_MAINN=\
	".\Debug\CaseModule.mod"\
	".\Debug\ChemistryModule.mod"\
	".\Debug\CompilerInfoModule.mod"\
	".\Debug\FlowsModule.mod"\
	".\Debug\InputModule.mod"\
	".\Debug\IOThreadModule.mod"\
	".\Debug\NWPMetModule.mod"\
	".\Debug\OpenMPModule.mod"\
	".\Debug\OutputModule.mod"\
	".\Debug\ParticleModule.mod"\
	".\Debug\PuffModule.mod"\
	".\Debug\RestartModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\SizeDistModule.mod"\
	".\Debug\SourceModule.mod"\
	".\Debug\SpeciesModule.mod"\
	".\Debug\TimerModule.mod"\
	
NODEP_F90_MAINN=\
	".\Debug\omp_lib.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Maths.F90
DEP_F90_MATHS=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\GlobalParametersModule.mod"\
	".\Debug\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Mets.F90
DEP_F90_METS_=\
	".\Debug\AncillaryMetModule.mod"\
	".\Debug\CommonMetModule.mod"\
	".\Debug\NWPMetModule.mod"\
	".\Debug\OpenMPModule.mod"\
	".\Debug\PrototypeMetModule.mod"\
	".\Debug\RadarMetModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\SingleSiteMetModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NWPFlow.F90
DEP_F90_NWPFL=\
	".\Debug\AncillaryMetModule.mod"\
	".\Debug\CommonFlowModule.mod"\
	".\Debug\FlowAndFlowProfileModule.mod"\
	".\Debug\MetsModule.mod"\
	".\Debug\NWPMetModule.mod"\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NWPMet.F90
DEP_F90_NWPME=\
	".\Debug\CommonMetModule.mod"\
	".\Debug\OpenMPModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\TimerModule.mod"\
	
NODEP_F90_NWPME=\
	".\Debug\grib_api.mod"\
	".\Debug\omp_lib.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\OpenMP.F90
DEP_F90_OPENM=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\StringModule.mod"\
	".\Debug\TimeModule.mod"\
	
NODEP_F90_OPENM=\
	".\Debug\omp_lib.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Output.F90
DEP_F90_OUTPU=\
	".\Debug\FlowsModule.mod"\
	".\Debug\FluctuationModule.mod"\
	".\Debug\OpenMPModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\SizeDistModule.mod"\
	".\Debug\SourceModule.mod"\
	".\Debug\SpeciesModule.mod"\
	".\Debug\TimerModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Particle.F90
DEP_F90_PARTI=\
	".\Debug\FlowsModule.mod"\
	".\Debug\OutputModule.mod"\
	".\Debug\PlumeRiseModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\SizeDistModule.mod"\
	".\Debug\SourceModule.mod"\
	".\Debug\SpeciesModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\ParticleSizeDistribution.F90
DEP_F90_PARTIC=\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\PhysicalUnits.F90
DEP_F90_PHYSI=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\GlobalParametersModule.mod"\
	".\Debug\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Physics.F90
DEP_F90_PHYSIC=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\GlobalParametersModule.mod"\
	".\Debug\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\PlumeRise.F90
DEP_F90_PLUME=\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\PrototypeFlow.F90
DEP_F90_PROTO=\
	".\Debug\CommonFlowModule.mod"\
	".\Debug\FlowAndFlowProfileModule.mod"\
	".\Debug\MetsModule.mod"\
	".\Debug\PrototypeMetModule.mod"\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\PrototypeMet.F90
DEP_F90_PROTOT=\
	".\Debug\CommonMetModule.mod"\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Puff.F90
DEP_F90_PUFF_=\
	".\Debug\FlowsModule.mod"\
	".\Debug\OutputModule.mod"\
	".\Debug\ParticleModule.mod"\
	".\Debug\PlumeRiseModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\SourceModule.mod"\
	".\Debug\SpeciesModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\RadarFlow.F90
DEP_F90_RADAR=\
	".\Debug\CommonFlowModule.mod"\
	".\Debug\FlowAndFlowProfileModule.mod"\
	".\Debug\MetsModule.mod"\
	".\Debug\RadarMetModule.mod"\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\RadarMet.F90
DEP_F90_RADARM=\
	".\Debug\CommonMetModule.mod"\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Restart.F90
DEP_F90_RESTA=\
	".\Debug\CaseModule.mod"\
	".\Debug\ChemistryModule.mod"\
	".\Debug\OutputModule.mod"\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Screen.F90
DEP_F90_SCREE=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\GlobalParametersModule.mod"\
	".\Debug\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\SemiLagrangian.F90
DEP_F90_SEMIL=\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Service.F90
DEP_F90_SERVI=\
	".\Debug\CoordinateSystemModule.mod"\
	".\Debug\ErrorAndMessageIIModule.mod"\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\GlobalParametersModule.mod"\
	".\Debug\GridAndDomainModule.mod"\
	".\Debug\HeadedFileModule.mod"\
	".\Debug\MathsModule.mod"\
	".\Debug\PhysicalUnitsModule.mod"\
	".\Debug\PhysicsModule.mod"\
	".\Debug\ScreenModule.mod"\
	".\Debug\SortModule.mod"\
	".\Debug\StringModule.mod"\
	".\Debug\SystemModule.mod"\
	".\Debug\TimeModule.mod"\
	".\Debug\UnitModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\SingleSite.F90
DEP_F90_SINGL=\
	".\Debug\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\SingleSiteFlow.F90
DEP_F90_SINGLE=\
	".\Debug\CommonFlowModule.mod"\
	".\Debug\FlowAndFlowProfileModule.mod"\
	".\Debug\MetsModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\SingleSiteMetModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\SingleSiteMet.F90
DEP_F90_SINGLES=\
	".\Debug\CommonMetModule.mod"\
	".\Debug\FlowAndFlowProfileModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\SingleSiteModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Sort.F90
DEP_F90_SORT_=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\GlobalParametersModule.mod"\
	
NODEP_F90_SORT_=\
	".\Debug\Iflport.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Source.F90
DEP_F90_SOURC=\
	".\Debug\FlowsModule.mod"\
	".\Debug\ServiceModule.mod"\
	".\Debug\SizeDistModule.mod"\
	".\Debug\SpeciesModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Species.F90
DEP_F90_SPECI=\
	".\Debug\ServiceModule.mod"\
	".\Debug\SizeDistModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\String.F90
DEP_F90_STRIN=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\GlobalParametersModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\System.F90
DEP_F90_SYSTE=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\StringModule.mod"\
	
NODEP_F90_SYSTE=\
	".\Debug\Iflport.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Time.F90
DEP_F90_TIME_=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\GlobalParametersModule.mod"\
	".\Debug\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Timers.F90
DEP_F90_TIMER=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\UnitModule.mod"\
	
NODEP_F90_TIMER=\
	".\Debug\omp_lib.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Unit.F90
DEP_F90_UNIT_=\
	".\Debug\ErrorAndMessageModule.mod"\
	".\Debug\GlobalParametersModule.mod"\
	".\Debug\StringModule.mod"\
	
NODEP_F90_UNIT_=\
	".\Debug\Iflport.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\CatCode
# End Source File
# Begin Source File

SOURCE=.\Changes.txt
# End Source File
# Begin Source File

SOURCE=.\Flows.P90

!IF  "$(CFG)" == "NameIII - Win32 Release"

# Begin Custom Build
InputPath=.\Flows.P90

"Flows.F90" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	..\Executables_Win\Preprocessor.exe Flows.P90 Flows.F90

# End Custom Build

!ELSEIF  "$(CFG)" == "NameIII - Win32 Debug"

# Begin Custom Build
InputPath=.\Flows.P90

"Flows.F90" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	..\Executables_Win\Preprocessor.exe Flows.P90 Flows.F90

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\LinuxIntelDebug
# End Source File
# Begin Source File

SOURCE=.\LinuxIntelRelease
# End Source File
# Begin Source File

SOURCE=.\LinuxIntelRelease_GRIBEX
# End Source File
# Begin Source File

SOURCE=.\Makefile
# End Source File
# Begin Source File

SOURCE=.\Mets.P90

!IF  "$(CFG)" == "NameIII - Win32 Release"

# Begin Custom Build
InputPath=.\Mets.P90

"Mets.F90" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	..\Executables_Win\Preprocessor.exe Mets.P90 Mets.F90

# End Custom Build

!ELSEIF  "$(CFG)" == "NameIII - Win32 Debug"

# Begin Custom Build
InputPath=.\Mets.P90

"Mets.F90" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	..\Executables_Win\Preprocessor.exe Mets.P90 Mets.F90

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\ReadMe.txt
# End Source File
# End Target
# End Project

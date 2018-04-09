# Microsoft Developer Studio Project File - Name="Name2Adms" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) QuickWin Application" 0x0107

CFG=Name2Adms - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Name2Adms.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Name2Adms.mak" CFG="Name2Adms - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Name2Adms - Win32 Release" (based on "Win32 (x86) QuickWin Application")
!MESSAGE "Name2Adms - Win32 Debug" (based on "Win32 (x86) QuickWin Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Name2Adms - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /libs:qwin /nologo /warn:nofileopt
# ADD F90 /compile_only /define:"CompaqPCCompiler" /define:"ExtraChecks" /fpp /libs:qwin /nologo /warn:nofileopt
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
# ADD LINK32 kernel32.lib /nologo /entry:"WinMainCRTStartup" /subsystem:windows /machine:I386 /nodefaultlib:"dfconsol.lib"
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=copy .\Release\Name2Adms.exe ..\Executables_Win\Name2Adms.exe
# End Special Build Tool

!ELSEIF  "$(CFG)" == "Name2Adms - Win32 Debug"

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
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /define:"CompaqPCCompiler" /define:"ExtraChecks" /fpp /libs:qwin /nologo /traceback /warn:argument_checking /warn:nofileopt
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
# ADD LINK32 kernel32.lib /nologo /entry:"WinMainCRTStartup" /subsystem:windows /incremental:no /debug /machine:I386 /nodefaultlib:"dfconsol.lib" /pdbtype:sept
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=copy .\Debug\Name2Adms.exe ..\Executables_Win\Name2AdmsDebug.exe
# End Special Build Tool

!ENDIF 

# Begin Target

# Name "Name2Adms - Win32 Release"
# Name "Name2Adms - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\Code_NameIII\CoordinateSystem.F90
NODEP_F90_COORD=\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\GlobalParametersModule.mod"\
	".\Release\MathsModule.mod"\
	".\Release\PhysicsModule.mod"\
	".\Release\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\ErrorAndMessage.F90
NODEP_F90_ERROR=\
	".\Release\GlobalParametersModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\ErrorAndMessageII.F90
NODEP_F90_ERRORA=\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\GlobalParametersModule.mod"\
	".\Release\HeadedFileModule.mod"\
	".\Release\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\GlobalParameters.F90
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\GridAndDomain.F90
NODEP_F90_GRIDA=\
	".\Release\CoordinateSystemModule.mod"\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\GlobalParametersModule.mod"\
	".\Release\MathsModule.mod"\
	".\Release\PhysicsModule.mod"\
	".\Release\StringModule.mod"\
	".\Release\TimeModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\HeadedFile.F90
NODEP_F90_HEADE=\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\GlobalParametersModule.mod"\
	".\Release\StringModule.mod"\
	".\Release\TimeModule.mod"\
	".\Release\UnitModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\MainName2Adms.F90
NODEP_F90_MAINN=\
	".\Release\ServiceModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\Maths.F90
NODEP_F90_MATHS=\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\GlobalParametersModule.mod"\
	".\Release\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\PhysicalUnits.F90
NODEP_F90_PHYSI=\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\GlobalParametersModule.mod"\
	".\Release\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\Physics.F90
NODEP_F90_PHYSIC=\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\GlobalParametersModule.mod"\
	".\Release\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\Screen.F90
NODEP_F90_SCREE=\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\GlobalParametersModule.mod"\
	".\Release\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\Service.F90
NODEP_F90_SERVI=\
	".\Release\CoordinateSystemModule.mod"\
	".\Release\ErrorAndMessageIIModule.mod"\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\GlobalParametersModule.mod"\
	".\Release\GridAndDomainModule.mod"\
	".\Release\HeadedFileModule.mod"\
	".\Release\MathsModule.mod"\
	".\Release\PhysicalUnitsModule.mod"\
	".\Release\PhysicsModule.mod"\
	".\Release\ScreenModule.mod"\
	".\Release\SortModule.mod"\
	".\Release\StringModule.mod"\
	".\Release\SystemModule.mod"\
	".\Release\TimeModule.mod"\
	".\Release\UnitModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\Sort.F90
NODEP_F90_SORT_=\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\GlobalParametersModule.mod"\
	".\Release\IFLPort.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\String.F90
NODEP_F90_STRIN=\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\GlobalParametersModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\System.F90
NODEP_F90_SYSTE=\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\IFLPort.mod"\
	".\Release\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\Time.F90
NODEP_F90_TIME_=\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\GlobalParametersModule.mod"\
	".\Release\StringModule.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Code_NameIII\Unit.F90
NODEP_F90_UNIT_=\
	".\Release\ErrorAndMessageModule.mod"\
	".\Release\GlobalParametersModule.mod"\
	".\Release\IFLPort.mod"\
	".\Release\StringModule.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\Changes.txt
# End Source File
# Begin Source File

SOURCE=.\LinuxIntelDebug
# End Source File
# Begin Source File

SOURCE=.\LinuxIntelRelease
# End Source File
# Begin Source File

SOURCE=.\ReadMe.txt
# End Source File
# End Target
# End Project

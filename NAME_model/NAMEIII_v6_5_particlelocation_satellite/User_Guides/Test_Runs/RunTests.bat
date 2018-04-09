echo off
rem # ******************************************************************************
rem #
rem # Project: NAME III documentation --> User Guides --> Testing Runs
rem #
rem # File:    batch file for executing NAME test runs on Windows systems
rem #
rem # Author:  Andrew Jones, Atmospheric Dispersion, UK Met Office
rem #
rem # Date:    31/01/2012
rem #
rem # ******************************************************************************

rem # Get the test directory

set testdir=%CD%

rem #
rem # Unzip the sample NWP met files
rem #

set metdir=%testdir%\Sample_NWP_Met_Data

echo Unzipping sample NWP met files (if these are still gzipped) ...
for %%f in ("%metdir%\MO*.gz") do "C:\Program Files\WINZIP\WinZip32.exe" -min -e -o "%%f" "%metdir%"
echo ... done

rem #
rem # NAME run using single-site met data
rem #

echo Performing NAME test run using single-site met data ...
..\..\Executables_Win\NameIII.exe Example_SingleSiteMet.txt -ClosePromptLevel=0
echo  ... done

rem # Set the NAME output directory

set workdir=%testdir%\Output_SingleSiteMet\

rem #
rem # Generate IDL graphics for single-site met run
rem #

echo Generating IDL graphics for NAME run based on single-site met data ...
call .\genIDLplots.bat "%workdir%"
echo ... done

rem #
rem # NAME run using NWP met data
rem #

echo Performing NAME test run using NWP met data ...
..\..\Executables_Win\NameIII.exe Example_NWPMet.txt -ClosePromptLevel=0
echo ... done

rem # Set the NAME output directory

set workdir=%testdir%\Output_NWPMet\

rem #
rem # Generate IDL graphics for NWP met run
rem #

echo Generating IDL graphics for NAME run based on NWP met data ...
call .\genIDLplots.bat "%workdir%"
echo ... done

echo Script completed
echo on


echo off
rem # ******************************************************************************
rem #
rem # Project: NAME III documentation --> User Guides --> Testing Runs
rem #
rem # File:    batch file for calling IDL plotting procedure on Windows systems
rem #
rem # Author:  Andrew Jones, Atmospheric Dispersion, UK Met Office
rem #
rem # Date:    31/01/2012
rem #
rem # ******************************************************************************

rem # Get the NAME output directory

set WORKDIR=%1

rem # Set the main testing directory and IDL graphics directory

cd %WORKDIR%
cd ..
set TESTDIR=%CD%

cd ..\..\Code_IDLGraphics
set GRAPHICSDIR="%CD%"

cd %TESTDIR%

rem # Run IDL and call plotting procedure in genidlplots.pro

echo .compile "%TESTDIR%\genidlplots.pro"    > genplots
echo genidlplots, %GRAPHICSDIR%, %WORKDIR%  >> genplots
echo exit, status=0                         >> genplots

echo Calling IDL ...
echo on
"C:\Program Files\ITT\IDL70\idlde\idlde.exe" -batch "%TESTDIR%\genplots"
echo off

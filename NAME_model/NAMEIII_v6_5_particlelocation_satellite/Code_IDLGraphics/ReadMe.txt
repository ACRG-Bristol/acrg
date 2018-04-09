Read me file for PV-Wave Graphics.

Note not all the routines work with Name III. 


Although we didn't generally record PV-Wave changes prior to 4.3 (these are kept with 
NAME II) we note here the following changes in usage (but not extensions of 
functionality) introduced on 6/7/06 (not for NAME II format output):
selcase value needed for plotseries changed from e.g. 'C1' to '1' (with a default of 1)
selcase introduced for plotfield (same values as for plotseries)
selgrid value needed for plotfield changed from e.g. 'grid1_C1_T' to 'Fields_grid1' 


Example input file for plotting fields:
selspecies=['SO2']
selfield=['Air Concentration']
sellevel=['Z = 50.00000 m agl']
SelTimeAverage=['No time averaging']
selgrid='Fields_grid1'
selcase='C'
seldatavalue='D = 50.00000'
selname='Req 1'
plotfield,'../Input/Output',selgrid,/exact,$
 selspecies=selspecies,selfield=selfield,sellevel=sellevel,$
 plotpmsl=0,plotuv=0,group='time',$
 multiplot=[0,2,2],plotheader=1,plotlegend=1,/genanim,$
 projection=1,ncontours=11,$
 namever=3,timeframe='absolute',SelTimeAverage=SelTimeAverage,$
 selcase=selcase,seldatavalue=seldatavalue,selname=selname

This requres Name III output to be produced with: 
  Separate File = 'T'
  Across = 'TZ' (D also needed for quantities depending on a data grid - must not
                 be a 'floating' D-Grid; can also include D if no D dependence)
  Output Format options to include I & A but not 2
  A T-grid, a structured regular H-grid, and no S-grid 
  
Note that selecting fields by specifying the 'selname' option will override the other
variable specifications (although these may also be required to ensure that headers,
etc. are correctly plotted in the graphics).

We can also add '2' to 'Output Format' and invoke pv-wave as for Name II. Here the Output 
Group must begin 'Fields_'. 

Example input file for plotting time series:
 selspecies=['SO2']
 selfield=['Air Concentration']
 sellevel=['Z = 50.00000 m agl']
 sellocation=['Receptor']
 selgrid='TimeSeries'
 selcase='1'
 tssplit='NONE'
 plotseries,'../NameIII/Output',selgrid,log=0,$
  selspecies=selspecies,selfield=selfield,sellocation=sellocation,$
  sellevel=sellevel,tssplit=tssplit,$
  multiplot=[0,1,1],plotheader=1,plotlegend=1,genanim=0,$
  namever=3,timeframe='absolute',selcase=selcase

This requres Name III output to be produced with: 
  Separate File = 'XY' or blank 
  Across = 'XYZ' (D also needed for quantities depending on a data grid - must not
                  be a 'floating' D-Grid; can also include D if no D dependence)
  Output Format options to include A & Z but not I or 2
  A T-grid, an unstructured H-grid of named points, and no S-grid 

Note that the TimeAverage option is not supported for time series yet
(so, e.g., if two series exist in the same file which have the same parameters except
for time-averaging information, only the first time series would be plotted).

We can also add '2' to 'Output Format' and invoke pv-wave as for Name II. Here the Output 
Group must begin 'Time_series_'. 

The plottraj routine also supports NAME III format output. The output format is
specified by the namever option (defaults to namever=2). All other options and
features for plottraj are the same.

 needs example call of plottraj ?? $$

 need more examples - e.g. for Name II format output $$

 keep these comments aligned with those describing NameII option at start of Output.F90 $$
 
All options are described in NAME documentation and listed in header of plotfield.pro. 
The user should note the directory path which should indicate the location of the 
NAMEIII output files.

Specific NAME III options are: 
namever =       3 for normal format output from NAME III 
                2 for NAME output and NAME II format output from NAME III.
timeframe =     'absolute' for absolute time frame
                'relative' for relative time frame.
selcase =       the case identifier in NAME III.
                This is normally 1 unless multiple cases are being run.
                PVwave routines will default to 1 if not entered.
seldatavalue =  'D = ?' where ? is a data value (probability or percentile).
selname =       Name of field requirement (can be used to specify fields which
                are not uniquely specified by the other recognised parameters).
SelTimeAverage = Time averaging/integrating information for a field.

Options with different usage between NAME and NAME III

selgrid =       for NAME III this is the 'Output Group' label as defined in
                the input file under  Output Requirements - Fields:
                e.g. TimeSeries.

Note not all Name III output can be processed by the PV Wave Graphics.

PV-Wave output goes in the same directory as the NAME III output.


Running:
-------

On PCs:

cd ..\Code_PVWaveGraphics
c:\vni\wave\bin\bin.i386nt\wave.exe PVWaveInputFile


----------------------------------------------------------------------------------
IDL Virtual Machine Read Me
----------------------------------------------------------------------------------

The IDL Virtual Machine (VM) is a program which runs pre-compiled IDL code. Unlike 
IDL the user is not required to have a license to run the code. Programs in this
folder which begin 'vm_' contain virtual machine programs intended to mirror 
their non-pre-compiled IDL counterparts. Currently only plotfield.pro, plottraj.pro
and plotseries.pro have VM counterparts. 
This file contains a short description of how to compile and run the IDL VM code for ADG.

Compiling VM code
~~~~~~~~~~~~~~~~~
For ease of compilation all necessary programs should be stored in one folder.
To compile:
cd to folder containing programs, then at the IDL command line type:
.FULL_RESET_SESSION                         ;this clears IDL's memory
.compile vm_plotfield                       ;compiles vm_plotfield
resolve_all                                 ;compiles subroutines used by plotfield
save, /routines, filename = <filename.sav>  ;saves compiled program and subroutine
                                            ;replace the name between <> with the full
                                            ;file path

Running VM code
~~~~~~~~~~~~~~~
Requires an input file which contains the inputs fed to plotfield, plottraj or
plotseries. 

Example input file contents (note that comments can be added before the first
line) and PLOTFIELD should be replaced with PLOTTRAJ or PLOTSERIES as relevant.
Also note that no quote marks or brackets are required. Multiple fields, levels 
etc. should be separated by a comma.

--
;;PLOTFIELD INPUTS
datadir=/home/h02/apsl/name_runs/volcano_run
selgrid=Fields_grid1
selfield=Air Concentration
sellevel= Z = 7.500000 FL
seltimeaverage=6hr 0min average
multiplot=0,2,2
exact=1
plotuv=0
namever=3
plotheader=1
plotlegend=1
ct_blue2red=0
genanim=0
gengif=0
--

To run on Linux:
(a) cd to the directory where the IDL code is stored then type:
     idl -vm="vm_plotfield.sav"
    program will prompt for input file

(b) Use a src file with contents as written below
  
-- 
#!/bin/ksh
echo `date`

WORKDIR=/home/h02/apsl/name_runs/volcano_run
OPTFILE=${WORKDIR}/plotfield_input.txt

GRAPHDIR=/home/h03/apdg/NameIII/NameIIIDevVersion/Version5_4a/Code_IDLGraphics

cd $GRAPHDIR
idl -vm="vm_plotfield.sav" -arg ${OPTFILE}
--

To run on a PC:
Locate the vm_plotfield.sav file and double click to run
Program will then prompt for input file

Note that IDL VM only requires the following files to run:
MO_Master_B.jpg          fieldtemplate.sav      vm_plotseries.sav
fieldheadtemplate.sav    readplotfieldopts.sav  vm_plottraj.sav
fieldheadtemplateV2.sav  vm_plotfield.sav



Read me file for Runs
Date:   24/8/07
Author: Dave Thomson

...\Runs contains various example scripts and input files and is a 
convenient folder for such files, but any folder can be used. ...\Resources contains
various resources such as met data for use in runs. 

($$ give more details / move some of this to a readme in resources?
This file (and the folders Resources and Runs) needs tidying up.)
                            
WinNameIII.bat - Script for running the code on Windows PCs. Double clicking launches
                 NameIII.

MetRestore.bat - Script for restoring met on Windows PCs.

In1.txt } Headed input files for running 4 test cases. The first three 
In2.txt } cases use prototype-met-module met files Met1.dat, Met2.dat and 
In3.txt } Met3.dat. The fourth case uses the NWP-met-module met files 
In4.txt } HP200101190000.REGH2001 and HP200101190300.REGH2001, the 
          NWP-met-module topography file Topogrh, and a headed input file 
          UMMetDefnRH2001.txt which contains the met definition for the 
          NWP met module. The HP... and Topogrh files are in the NAME III 
          library and will need to be copied to the Met folder. 
In5.txt     } Headed input files for some other test cases involving the 
In6.txt     } single site met and flow modules, buildings, plume rise,
In7.txt     } deposition and decay, and fast run with large time steps.
In8.txt     } In9.txt needs Nimrod data which is also in the NAME III 
In9.txt     } library and needs to be copied to the Nimrod folder.
In4fast.txt }
plus other input files

Note: NWP Met (including associated topography data) and Nimrod data is 
      unformatted and so needs to be ftp'ed as binary between PC and Unix 
      platforms.

Met1.dat } Example met files for the prototype met module.  
Met2.dat } 
Met3.dat }

MetDemo.met - Example met file for the single site met module.

Kincaid met files.rtf } Kincaid met files and notes on these.
kincaidfullrep.met    }
kincaidfullreped.met  }

NotesOnUMMetFiles.txt } Headed input files containing met definitions for 
UMMetDefnGH2001.txt   } using UM data in the NWP met module, together  
UMMetDefnRH2001.txt   } with notes on these and on the UM met files 
UMMetDefnMH2001.txt   } themselves.
UMMetDefnGUM5.txt     }
UMMetDefnRUM5.txt     }
UMMetDefnMUM5.txt     }

Some simple example runs - input files, scripts to run PV-Wave on Windows and Linux,
and plots of the results:

Note the PV-Wave graphics output may not include plots of all available fields. Other
fields could be plotted by editing the wave script. The GIS plots cannot currently be 
produced without extra software not included with Name III and the current Met Office 
GIS software requires that the run is done with 'Name II format output' requested.

Input file:  Example_NWPMet.txt             
             Example_NWPMet_PlumeRise.txt
             Example_SingleSiteMet.txt
Wave script: WinWaveExample.bat LinuxWaveExample.sh
Wave output: .\Output\Plot_NWPMet.ps
             .\Output\Plot_NWPMet_PlumeRise.ps
             .\Output\Plot_SingleSiteMet.ps
                       
Input file:  Example_NWPMet_ShortRange.txt
Wave script: WinWaveExample_NWPMet_ShortRange.bat 
Wave output: Plot_NWPMet_ShortRange.pdf
GIS output:  AirConcentration_Boundarylayer_TRACER_20070216{1900|2000|2100|2200}.png





#######################################################################
#######################################################################
##
##  Main program for generating example plots
##
#######################################################################
##
##
## Eike Mueller, Atmospheric DispersionGroup, Met Office July 2011
##
#######################################################################
#######################################################################

import sys
import os
import getopt
import re
import glob
import Field2d
import math
import numpy as np
import datetime
from matplotlib import pyplot as plt
import matplotlib
from mpl_toolkits.basemap import Basemap

#######################################################################
#######################################################################
##
##  M A I N
##
#######################################################################
#######################################################################
def Main(DataDir,
         GridName,
         PlotDir,
         FieldName,
         FieldLabel):

  # Find all files with a given grid name
  NAMEOutputFiles = glob.glob(DataDir+'/Fields_'+GridName+'*.txt')
  TMin = datetime.datetime(3000,1,1)
  TMax = datetime.datetime(   1,1,1)
  nT = len(NAMEOutputFiles)
  for Filename in NAMEOutputFiles:
    R = re.match('.*([0-9]{4})([0-9]{2})([0-9]{2})([0-9]{2})([0-9]{2})\.txt',Filename)
    Year = int(R.group(1))
    Month = int(R.group(2))
    Day = int(R.group(3))
    Hour = int(R.group(4))
    Minute = int(R.group(5))
    T = datetime.datetime(Year,Month,Day,Hour,Minute)
    if (T > TMax):
      TMax = T
    if (T < TMin):
      TMin = T
  if (nT > 1):
    dT = (TMax - TMin)/(nT-1)
  else:
    dT = datetime.timedelta(hours=1)
  T0 = TMin

  if (nT > 0):
    # Read fields from file    
    (PlotRange, Thresholds, FieldList) = ReadFields(DataDir,
                                                    GridName,
                                                    FieldName,
                                                    T0, dT, nT)
    # Plot fields
    PlotFields(PlotDir, 
               FieldLabel,
               PlotRange,
               Thresholds,
               T0, dT, nT,
               FieldList)
    
#######################################################################
# Read a set of fields from disk.
#######################################################################
# Read nT data files, which are assumed to be of the form
# Fields_grid<FieldIndex>_C<CaseIndex>_T<iT>_<TimeStamp>.txt
# The first timestep to be plotted is given by T0 and the time increment
# by dT.
#######################################################################
def ReadFields(OutputDir,
               GridName,
               FieldName,
               T0,
               dT,
               nT,
               CaseIndex=1):
  # This defines the number of orders of magnitude that are plotted, 
  # e.g. if the maximal value of a field is 1.45E-9, and nDecades is 4,
  # then only concentrations in the range [1E-12, 1E-8] will be plotted
  nDecades = 4
  T = T0
  FieldList = []
  # The range of concentrations found in the files is stored in these 
  # variables.
  MinField = 1.0E12
  MaxField = -1.0
  # Loop over all files and extract field data into Field2d object
  for iT in range(1,nT+1):
    # Construct filename
    TimeStamp = T.strftime('%Y%m%d%H%M')    
    Filename = OutputDir+'Fields_'+GridName+'_C'+str(CaseIndex)+'_T'+str(iT)+'_'+TimeStamp+'.txt'
    print 'Reading field \''+FieldName+'\' from file \''+Filename+'\''
    # Read field data from file
    Field = Field2d.Field2d(Filename,Name=FieldName)
    FieldList.append(Field)
    # Extract minimal and maximal value of these fields (ignore 0.0 s)
    maxData = Field.Data.max()
    minData = (np.ma.masked_less_equal(Field.Data, 0.0)).min()
    if (maxData > MaxField):
      MaxField = maxData
    if (minData < MinField):
      MinField = minData
    T += dT
  # Determine plotting boundaries from latitude/longitude arrays
  LatMin = min(Field.Lats)
  LatMax = max(Field.Lats)
  LonMin = min(Field.Lons)
  LonMax = max(Field.Lons)

  # Work out upper and lower boundaries of the logarithmic plot range
  UpperThreshold = np.floor(np.log10(MaxField))+1
  LowerThreshold = UpperThreshold - nDecades
  Thresholds = {'Lower': LowerThreshold, 
                'Upper': UpperThreshold}
  PlotRange = {'LatMin': LatMin, 
              'LatMax': LatMax,
              'LonMin': LonMin,
              'LonMax': LonMax}
  return (PlotRange, Thresholds, FieldList)
      
#######################################################################
# Plot fields to file
#######################################################################
def PlotFields(PlotDir,
               FieldName,
               PlotRange,
               Thresholds,               
               T0, dT, nT,
               FieldList):
   
  norm_graysteps=matplotlib.colors.BoundaryNorm(np.arange(Thresholds['Lower'],Thresholds['Upper']+0.5,0.5),ncolors=256, clip = False)

  # set up map
  m = Basemap(resolution='h',projection='merc',\
              llcrnrlon=PlotRange['LonMin'], urcrnrlon=PlotRange['LonMax'],
              llcrnrlat=PlotRange['LatMin'], urcrnrlat=PlotRange['LatMax'])
  
  T = T0
  iT = 1
  # Now loop over all timesteps and plot the fields on a logarithmic scale
  for Field in FieldList:

    plt.clf()
    
    # Extract source location
    # Field.FileHeader['ReleaseLocation'] is of the form 'Lon E/W Lat N/S', e.g. '5.0 W 45.0 N'
    R = re.match('\s*([0-9\.]+)\s*([EW])\s*([0-9\.]+)\s*([NS])\s*', Field.FileHeader['ReleaseLocation'])
    if R:
      SourceLon = float(R.group(1))
      SourceLat = float(R.group(3))
      if (R.group(2) == 'W'):
         SourceLon *= -1.0
      if (R.group(4) == 'W'):
         SourceLat *= -1.0

    # Mark and label the source location
    Lats = [ SourceLat ]
    Lons = [ SourceLon ]    
    Locations = ['Source']
    X,Y = m(Lons,Lats)
    m.plot(X,Y,'ro',markersize=2.0)
    for name,xpt,ypt in zip(Locations,X,Y):
      plt.text(xpt+10,ypt+10,name,fontdict=dict(size=12,color='black'))
      
    # draw parallels, meridians coastlines etc.
    m.drawcoastlines()
#    m.drawparallels(np.arange(-90.,90.,1.),labels=[1,1,1,1])
#    m.drawmeridians(np.arange(-180.,180.,1.),labels=[1,1,1,1])
    ColorLand='#dcc2b0'
    ColorWater='#e9f2ff'
    m.fillcontinents(color=ColorLand,lake_color=ColorWater,zorder=-1)
    m.drawmapboundary(fill_color=ColorWater)
    m.drawcountries()

    
    X,Y = m(np.outer(Field.Lons,np.ones(len(Field.Lats))),
            np.outer(np.ones(len(Field.Lons)),Field.Lats) )
    # Plot field data
    PlotField = m.pcolor(X,Y,np.ma.masked_outside(np.log10(Field.Data),Thresholds['Lower'],Thresholds['Upper']),
                         cmap=plt.get_cmap('jet'),
                         norm=norm_graysteps,
                         shading='flat')
    cbar = plt.colorbar(orientation='horizontal',
                        ticks=np.arange(Thresholds['Lower'],Thresholds['Upper']+1,1),
                        shrink=0.6,
                        format='$10^{%d}$',
                        pad=0.04,
                        fraction=0.07)
    cbar.set_label(Field.ColumnHeader['Name']+' ['+Field.ColumnHeader['Units']+']')
    # Save plot to output file
    TimeStamp = Field.ColumnHeader['T']
    plt.title('Run name: '+Field.FileHeader['RunName']+'\nValid at: '+TimeStamp)
    iTstr = str(iT)
    while (len(iTstr) < 5):
      iTstr = '0'+iTstr        
    OutFilename = PlotDir+FieldName+'_'+'T'+iTstr+'.png'
    print 'Plotting file \''+OutFilename+'\''
    plt.savefig(OutFilename)
    T += dT
    iT += 1

#######################################################################
# Call main program
#######################################################################
      
if (__name__ == '__main__'):
   argc = len(sys.argv)

   optlist, args = getopt.getopt(sys.argv[1:], '', ['datadir=', 'plotdir=', 'grid=','fieldname=','label='])
   DataDir = None
   GridName = None
   PlotDir = None
   FieldName = None
   FieldLabel = None
   for Key, Value in optlist:
     if (Key == '--datadir'):
       DataDir = Value
     if (Key == '--grid'):
       GridName = Value
     if (Key == '--plotdir'):
       PlotDir = Value
     if (Key == '--fieldname'):
       FieldName = Value
     if (Key == '--label'):
       FieldLabel = Value
   if ( (DataDir == None) or 
        (GridName == None) or 
        (PlotDir == None) or
        (FieldName == None) or 
        (FieldLabel == None) ):
     print 'Usage: python '+sys.argv[0]+' --datadir=<datadir>  --plotdir=<plotdir>  --grid=<grid>  --fieldname=<fieldname>  -label=<label>'
     sys.exit(-1)
   else:
     Main(DataDir, GridName, PlotDir, FieldName, FieldLabel)

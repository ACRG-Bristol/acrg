import re
import numpy as np

#######################################################################
#######################################################################
##
##  Class representing a two dimensional horizontal fiels
##
#######################################################################
##
##  This class represents a two dimensional field in a NAME 3 formatted
##  output file.
##  In addition to the actual data (in self.data) and
##  the latitude and longitude arrays (in self.Lats and self.Lons) it 
##  also stores metadata on the field in the dictionaries FileHeader
##  (header of the file the data has been read from) and ColumnHeader
##  (header of the column in the file).
##
##  Data is read from the file when the object is instantiated.
##  
##
##  Eike Mueller, Atmospheric Dispersion Group, Met Office July 2011
##
#######################################################################
#######################################################################

class Field2d():

#-------------------------------------------------  
# Constructor. 
#
#   Create a new Field2d object by reading the 
#   column specified by the keywords from the 
#   file Filename.
#
#-------------------------------------------------  
  def __init__(self,
               Filename,                   # Name of file to read
               SpeciesCategory='',         # Species category, e.g. 'Tracer'
               Name='',                    # Name of field
               Quantity='',                # Quantity, e.g. 'Air Concentration'
               Species='',                 # Species, e.g. 'Volcanic Ash'
               Units='',                   # Units, e.g. 'g / m^3'
               Source='',                  # Source, e.g. 'All sources'
               EnsembleAvgInfo='',         # Ensemble avg. informattion, e.g. 'No ensemble averaging'
               TimeAvgInfo='',             # Time avg. information, e.g. '3hr 0min average'
               HorizontalAvgInfo='',       # Horizontal avg. information, e.g. 'No horizontal averaging'
               VerticalAvgInfo='',         # Vertical avg. information, e.g. 'No vertical averaging'
               ProbPerc='',                # Probability and percentiles
               ProbPercEnsemble='',        # Probabilities and percentiles over ensemble
               ProbPercTime='',            # Probabilities and percentiles over time
               T='',                       # Time of validity string, e.g. '23/05/2011 12:00 UTC'
               Z=''                        # Vertical position string, e.g. 'Z = 50.00000 m agl'
              ):
     
  # Initialise fields
    self.verbose=False
    self.Filename = Filename
    self.ColumnHeader = {'SpeciesCategory' : SpeciesCategory,
                         'Name' : Name, 
                         'Quantity' : Quantity,
                         'Species' : Species,
                         'Units' : Units,
                         'Source' : Source,
                         'EnsembleAvgInfo' : EnsembleAvgInfo,
                         'TimeAvgInfo' : TimeAvgInfo,
                         'HorizontalAvgInfo' : HorizontalAvgInfo,
                         'VerticalAvgInfo' : VerticalAvgInfo,
                         'ProbPerc' : ProbPerc,
                         'ProbPercEnsemble' : ProbPercEnsemble,
                         'ProbPercTime' : ProbPercTime,
                         'T' : T,
                         'Z' : Z}
   # Read Header and data
    self.FileHeader = None
    self.X = None
    self.Y = None
    self.Data = None
    self.ReadFileHeader()
    self.ReadColumnData()
   
#-------------------------------------------------  
# Read file header from file
#-------------------------------------------------  
  def ReadFileHeader(self):
    self.FileHeader = {}
    File = None
    try:
      File = open(self.Filename,'r')
    except IOError:
      print 'Error opening \''+self.Filename+'\' for reading.'
    if (self.verbose):
        print 'Reading file header from \''+self.Filename+'\''
    KeywordList = {'Run name':'RunName',
                   'Run time':'RunTime',
                   'Met data':'MetData',
                   'Start of release':'StartOfRelease',
                   'End of Release':'EndOfRelease',
                   'Source strength':'SourceStrength',
                   'Release location':'ReleaseLocation',
                   'Release height':'ReleaseHeight',
                   'Run duration':'RunDuration',
                   'X grid origin':'XGridOrigin',
                   'Y grid origin':'YGridOrigin',
                   'X grid size':'XGridSize',
                   'Y grid size':'YGridSize',
                   'X grid resolution':'XGridResolution',
                   'Y grid resolution':'YGridResolution'
                   }
    if File:
      for Line in File:
        for Keyword in KeywordList.keys():
            m = re.match('\s*'+Keyword+'\s*:\s*(.*)',Line)
            if m:
                self.FileHeader[KeywordList[Keyword]] = m.group(1).strip()
                if (self.verbose):
                    print '   '+Keyword+' = '+self.FileHeader[KeywordList[Keyword]]
                
        if re.match('\s*Fields:\s*',Line):
            break
      File.close()

#-------------------------------------------------  
# Read data column from file
#-------------------------------------------------  
  def ReadColumnData(self):

    def ColMatches(A,B):
       return ((A.strip().lstrip() == B.strip()) or (B.strip() == ''))
    File = None
    try:
      File = open(self.Filename,'r')
    except IOError:
      print 'Error opening \''+self.Filename+'\' for reading.'
    if File:
      for Line in File:
        if re.match('\s*Fields:\s*',Line):
            break
      ColumnHeads = []
      for i in range(0,15):
        ColumnHeads.append(File.next().split(','))
      nColumns = len(ColumnHeads[0])
      iColumn = -1
      ColumnFound = False
  # Find first column that matches all search fields
      for i in range(0,nColumns):        
        if (   ColMatches(ColumnHeads[0][i],self.ColumnHeader['SpeciesCategory'])
          and  ColMatches(ColumnHeads[1][i],self.ColumnHeader['Name'])
          and  ColMatches(ColumnHeads[2][i],self.ColumnHeader['Quantity'])
          and  ColMatches(ColumnHeads[3][i],self.ColumnHeader['Species'])
          and  ColMatches(ColumnHeads[4][i],self.ColumnHeader['Units'])
          and  ColMatches(ColumnHeads[5][i],self.ColumnHeader['Source'])
          and  ColMatches(ColumnHeads[6][i],self.ColumnHeader['EnsembleAvgInfo'])
          and  ColMatches(ColumnHeads[7][i],self.ColumnHeader['TimeAvgInfo'])
          and  ColMatches(ColumnHeads[8][i],self.ColumnHeader['HorizontalAvgInfo'])
          and  ColMatches(ColumnHeads[9][i],self.ColumnHeader['VerticalAvgInfo'])
          and  ColMatches(ColumnHeads[10][i],self.ColumnHeader['ProbPerc'])
          and  ColMatches(ColumnHeads[11][i],self.ColumnHeader['ProbPercEnsemble'])
          and  ColMatches(ColumnHeads[12][i],self.ColumnHeader['ProbPercTime'])
          and  ColMatches(ColumnHeads[13][i],self.ColumnHeader['T'])
          and  ColMatches(ColumnHeads[14][i],self.ColumnHeader['Z'])):
          iColumn = i
          ColumnFound = True
          self.ColumnHeader['SpeciesCategory'] = ColumnHeads[0][i].strip().lstrip()
          self.ColumnHeader['Name'] = ColumnHeads[1][i].strip().lstrip()
          self.ColumnHeader['Quantity'] = ColumnHeads[2][i].strip().lstrip()
          self.ColumnHeader['Species'] = ColumnHeads[3][i].strip().lstrip()
          self.ColumnHeader['Units'] = ColumnHeads[4][i].strip().lstrip()
          self.ColumnHeader['Source'] = ColumnHeads[5][i].strip().lstrip()
          self.ColumnHeader['EnsembleAvgInfo'] = ColumnHeads[6][i].strip().lstrip()
          self.ColumnHeader['TimeAvgInfo'] = ColumnHeads[7][i].strip().lstrip()
          self.ColumnHeader['HorizontalAvgInfo'] = ColumnHeads[8][i].strip().lstrip()
          self.ColumnHeader['VerticalAvgInfo'] = ColumnHeads[9][i].strip().lstrip()
          self.ColumnHeader['ProbPerc'] = ColumnHeads[10][i].strip().lstrip()
          self.ColumnHeader['ProbPercEnsemble'] = ColumnHeads[11][i].strip().lstrip()
          self.ColumnHeader['ProbPercTime'] = ColumnHeads[12][i].strip().lstrip()
          self.ColumnHeader['T'] = ColumnHeads[13][i].strip().lstrip()
          self.ColumnHeader['Z'] = ColumnHeads[14][i].strip().lstrip()
          break
      if (ColumnFound):
        RegEx = '^'+2*'\s*([0-9]+)\s*,'+(nColumns-3)*'\s*([0-9\.\+\-Ee]+)\s*,'+'\s*$'
    
        nX = int(self.FileHeader['XGridSize'])
        nY = int(self.FileHeader['YGridSize'])
        XOrigin = float(self.FileHeader['XGridOrigin'])
        YOrigin = float(self.FileHeader['YGridOrigin'])
        XResolution = float(self.FileHeader['XGridResolution'])
        YResolution = float(self.FileHeader['YGridResolution'])
    
        self.Lons = np.zeros((nX+1))
        self.Lats = np.zeros((nY+1))
        for iX in range(0,nX+1):
          self.Lons[iX] = XOrigin+iX*XResolution
        for iY in range(0,nY+1):
          self.Lats[iY] = YOrigin+iY*YResolution
        self.Data = np.zeros((nX,nY))
        for Line in File:
          if (not re.match('\s*X Index.*',Line)):
            m = re.match(RegEx,Line)
            if m:
              iX = int(m.group(1))
              iY = int(m.group(2))
              self.Data[iX-1,iY-1] = float(m.group(iColumn+1))
      else:
        self.Lons = None
        self.Lats = None
        self.Data = None
        raise Field2dError('Can not find data column.')
      File.close()

#######################################################################
##
##  Class for exceptions of Field 2d class
##
#######################################################################

class Field2dError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

# -*- coding: utf-8 -*-
"""
Created on Tue May 31 09:52:07 2016

@author: bs15965
"""
import subprocess

#Arguments for R script
usr='bs15965' # change to your username
start="01/06/2006"  #format dd/mm/yyyy
end="02/06/2006"    #format dd/mm/yyyy
tz = 'UTC'   #time zone
ntimes = '1' #number of equal periods to split the date range into
ghgName = 'ch4'    # choice of "ch4", "co2", "n2o"
proj = "OSGB" # only option is OSGB at the he moment but should be "OSGB" and "LonLat" soon

unitType = "mol" #choice of "mol", "g"
unitSIprefix = "nano" #choice of "kilo", "none", "milli", "micro", "nano", "pico"
sectorList = 'None' # default 1:10

# Define command and arguments
command = 'Rscript'
path2script = '/home/'+usr+'/acrg/ukghg.R'    #assumes repository is in your home directory
# Variable number of args in a list
args = [start,end,tz,ntimes,ghgName,proj,unitType,unitSIprefix ,sectorList]

# Build subprocess command
cmd = [command, path2script] + args

# runs the R script and outputs ncdf files to the working directory
# there currently seems to be no option to put them elsewhere so will need to move
x = subprocess.check_call(cmd)

sector='total'  #choose which of the sector ncdf files you want to import

path='/home/'+usr+'/acrg/uk_flux_'+sector+'_'+ghgName+'_OSGB.nc' #path to ncdf file
with xray.open_dataset(path) as ds:
            ds.load()

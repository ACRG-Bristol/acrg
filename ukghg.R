# Fetch command line arguments
myArgs <- commandArgs(trailingOnly = TRUE)
library(ukghg)
startDate <- as.POSIXct(strptime(myArgs[1], "%d/%m/%Y"), tz = myArgs[3])
endDate <- as.POSIXct(strptime(myArgs[2], "%d/%m/%Y"), tz = myArgs[3])
nTimes <- as.integer(myArgs[4])
# create a sequence of timestamps
datect <- seq(startDate, endDate, length = nTimes)
# calculate fluxes for these times
myFlux <- calcFlux(myArgs[5], datect, myArgs[6], myArgs[7], myArgs[8],sectorList = ifelse(myArgs[9]=='None',1:10,myArgs[9]))

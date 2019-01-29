##########################################################################################
#
# The purpose of this script is to examine drying-ratios in the events detected by
# ../find_ARs.percentiles.pro. This script uses composites of the WRF synoptic fields at
# corresponding to event times. The composites are created by ../compositeEvents.pro.
#
##########################################################################################
import os,fnmatch,netCDF4,numpy as np
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

# Directory containing event composites.
dir = '/data/dswales/NA-CORDEX/ARdet/composites/'

# What composites to look at?
# WRF horizontal resolution (50 or 25)
res         = 50
# Event strength threshold (percentile)
ptile       = '97.00'
# Event persistence threshold (hours)
persistence = 24
# Model ID
modelID = 'mpi'

# Get filenames
fileH_Composite = find('composite.raw.events.'+str(res)+'km.'+ptile+'ptile.'+str(persistence)+'hrs.'+modelID+'.nc', dir)
fileF_Composite = find('composite.raw.events.future.'+str(res)+'km.'+ptile+'ptile.'+str(persistence)+'hrs.'+modelID+'.nc', dir)

# Read in data
dataIN = netCDF4.Dataset(fileH_Composite[0],'r')
year   = dataIN.variables['year'][:]
month  = dataIN.variables['month'][:]
day    = dataIN.variables['day'][:]
hour   = dataIN.variables['hour'][:]
lon    = dataIN.variables['lon'][:]
lat    = dataIN.variables['lat'][:]
nlon   = len(lon[:,0])
nlat   = len(lat[0,:])
ntime  = len(year)

for iTime in range(0,ntime-1):
    ivtU = dataIN.variables['ivtU'][iTime,:,:]
    ivtV = dataIN.variables['ivtV'][iTime,:,:]


##########################################################################################
# END PROGRAM
##########################################################################################

##########################################################################################
#
# The purpose of this script is to examine drying-ratios in the events detected by
# ../find_ARs.percentiles.pro. This script uses composites of the WRF synoptic fields at
# corresponding to event times. The composites are created by ../compositeEvents.pro.
#
##########################################################################################
import os,fnmatch,netCDF4,numpy as np
from math import radians, cos, sin, asin, sqrt
##########################################################################################
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result
##########################################################################################
def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    # Radius of earth in kilometers is 6371
    km = 6371* c
    return km
##########################################################################################
# Configuration
##########################################################################################

# Assign some name for configuration, do this to avoid overwriting previous results.
configID = "Take1"

# What composites to look at?
# WRF horizontal resolution (50 or 25)
res         = 50
# Event strength threshold (percentile)
ptile       = '96.00'
# Event persistence threshold (hours)
persistence = 24
# Model ID
modelID = 'mpi'

# Where to output results?
dirOUT  = "/data/dswales/NA-CORDEX/ARdet/dryingRatios/"
fileOUT = dirOUT+"dryingRatios."+configID+"."+str(res)+'km.'+ptile+'ptile.'+str(persistence)+'hrs.'+modelID+".nc"

# Where to compute drying ratios?
# 1) Define a line ([x0,y0],[xF,yF]), number of points (nDRs), and a distance downstream
#    (d_dnstream; in Degrees) from that line
x0         = -123.50
xF         = -122.75
y0         =  42.50
yF         =  47.50
nDRs       = 20
d_dnstream = 0.5
lon_upstream  = [i for i in np.arange(x0,xF+((xF-x0)/(nDRs-1)),((xF-x0)/(nDRs-1)))]
lat_upstream  = [i for i in np.arange(y0,yF+((yF-y0)/(nDRs-1)),((yF-y0)/(nDRs-1)))]
x0d   = x0 + d_dnstream
xFd   = xF + d_dnstream
lon_dnstream  = [i for i in np.arange(x0d,xFd+((xFd-x0d)/(nDRs-1)),((xFd-x0d)/(nDRs-1)))]
lat_dnstream  = [i for i in np.arange(y0,yF+((yF-y0)/(nDRs-1)),((yF-y0)/(nDRs-1)))]
# OR 2) Define a set of points (if you just want to manually enter a few points)
#lon_upstream = [-110.1, -120.7]
#lat_upstream = [  45.0,   46.0]
#lon_dnstream = [-109.1, -119.7]
#lat_dnstream = [  45.5,   46.5]
#nDRs         = len(lon_upstream)

# Directory containing event composites.
dir = '/data/dswales/NA-CORDEX/ARdet/composites/'

##########################################################################################
# i) Read in metadata.
#  - Read in grid and determine WRF gridpoint nearest to requested drying-ratio calculation
##########################################################################################
# Get filenames
fileH_Composite = find('composite.raw.events.'+str(res)+'km.'+ptile+'ptile.'+str(persistence)+'hrs.'+modelID+'.nc', dir)
fileF_Composite = find('composite.raw.events.future.'+str(res)+'km.'+ptile+'ptile.'+str(persistence)+'hrs.'+modelID+'.nc', dir)

# Open file
dataIN = netCDF4.Dataset(fileH_Composite[0],'r')
lon    = dataIN.variables['lon'][:]
lat    = dataIN.variables['lat'][:]
yearH  = dataIN.variables['year'][:]
monthH = dataIN.variables['month'][:]
dayH   = dataIN.variables['day'][:]
hourH  = dataIN.variables['hour'][:]
ntimeH = len(yearH)
nlon   = len(lon[:,0])
nlat   = len(lat[0,:])
dlonlat_upstream = np.zeros([nlon,nlat])
dlonlat_dnstream = np.zeros([nlon,nlat])
loni_upstream    = np.zeros(nDRs)
lati_upstream    = np.zeros(nDRs)
loni_dnstream    = np.zeros(nDRs)
lati_dnstream    = np.zeros(nDRs)
for il in range(0,nDRs):
    for ij in range(0,nlon):
        for ik in range(0,nlat):
            dlonlat_upstream[ij,ik] = haversine(lon[ij,ik], lat[ij,ik], lon_upstream[il], lat_upstream[il])
            dlonlat_dnstream[ij,ik] = haversine(lon[ij,ik], lat[ij,ik], lon_dnstream[il], lat_dnstream[il])
    # Store indices.
    kup               = dlonlat_upstream.argmin()
    ncol              = dlonlat_upstream.shape[1]
    loni_upstream[il] = kup/ncol
    lati_upstream[il] = kup%ncol
    kdn               = dlonlat_dnstream.argmin()
    ncol              = dlonlat_dnstream.shape[1]
    loni_dnstream[il] = kdn/ncol
    lati_dnstream[il] = kdn%ncol

##########################################################################################
# ii) Create output file
##########################################################################################
dataOUT = netCDF4.Dataset(fileOUT,'w',format='NETCDF4_CLASSIC')
dataOUT.createDimension("case",nDRs)
dataOUT.createDimension("timeH",ntimeH)
yearH_out     = dataOUT.createVariable("yearH",         "i4", ("timeH"))
monthH_out    = dataOUT.createVariable("monthH",        "i4", ("timeH"))
dayH_out      = dataOUT.createVariable("dayH",          "i4", ("timeH"))
hourH_out     = dataOUT.createVariable("hourH",         "i4", ("timeH"))
ivt_upstreamH = dataOUT.createVariable("IVT_upstreamH", "f4", ("case","timeH"))
ivt_dnstreamH = dataOUT.createVariable("IVT_dnstreamH", "f4", ("case","timeH"))
z0k_upstreamH = dataOUT.createVariable("z0k_upstreamH", "f4", ("case","timeH"))
z0k_dnstreamH = dataOUT.createVariable("z0k_dnstreamH", "f4", ("case","timeH"))

##########################################################################################
# iii) Now that we know the indices of the points we are interested in, readin in data at
#      these points, compute drying-ratio, and write output to netCDF.
##########################################################################################
# i) Historical period.

for iTime in range(0,ntimeH):
    yearH_out[iTime]  = dataIN.variables['year' ][iTime]
    monthH_out[iTime] = dataIN.variables['month'][iTime]
    dayH_out[iTime]   = dataIN.variables['day'  ][iTime]
    hourH_out[iTime]  = dataIN.variables['hour' ][iTime]
    for il in range(0,nDRs):
        ivt_upstreamH[il,iTime] = dataIN.variables['ivt'][iTime,lati_upstream[il],lon_upstream[il]]
        z0k_upstreamH[il,iTime] = dataIN.variables['z0k'][iTime,lati_upstream[il],lon_upstream[il]]
        ivt_dnstreamH[il,iTime] = dataIN.variables['ivt'][iTime,lati_dnstream[il],lon_dnstream[il]]
        z0k_dnstreamH[il,iTime] = dataIN.variables['z0k'][iTime,lati_dnstream[il],lon_dnstream[il]]

# ii) Future period (only GCMs)
if (modelID != 'erain'):
    dataIN = netCDF4.Dataset(fileF_Composite[0],'r')
    yearF  = dataIN.variables['year'][:]
    monthF = dataIN.variables['month'][:]
    dayF   = dataIN.variables['day'][:]
    hourF  = dataIN.variables['hour'][:]
    ntimeF = len(yearF)
    # Add future fields to outputs
    dataOUT.createDimension("timeF",ntimeF)
    yearF_out     = dataOUT.createVariable("yearF",         "i4", ("timeF"))
    monthF_out    = dataOUT.createVariable("monthF",        "i4", ("timeF"))
    dayF_out      = dataOUT.createVariable("dayF",          "i4", ("timeF"))
    hourF_out     = dataOUT.createVariable("hourF",         "i4", ("timeF"))
    ivt_upstreamF = dataOUT.createVariable("IVT_upstreamF", "f4", ("case","timeF"))
    ivt_dnstreamF = dataOUT.createVariable("IVT_dnstreamF", "f4", ("case","timeF"))
    z0k_upstreamF = dataOUT.createVariable("z0k_upstreamF", "f4", ("case","timeF"))
    z0k_dnstreamF = dataOUT.createVariable("z0k_dnstreamF", "f4", ("case","timeF"))
    for iTime in range(0,ntimeF):
        yearF_out[iTime]  = dataIN.variables['year' ][iTime]
        monthF_out[iTime] = dataIN.variables['month'][iTime]
        dayF_out[iTime]   = dataIN.variables['day'  ][iTime]
        hourF_out[iTime]  = dataIN.variables['hour' ][iTime]
        for il in range(0,nDRs):
            ivt_upstreamF[il,iTime] = dataIN.variables['ivt'][iTime,lati_upstream[il],lon_upstream[il]]
            z0k_upstreamF[il,iTime] = dataIN.variables['z0k'][iTime,lati_upstream[il],lon_upstream[il]]
            ivt_dnstreamF[il,iTime] = dataIN.variables['ivt'][iTime,lati_dnstream[il],lon_dnstream[il]]
            z0k_dnstreamF[il,iTime] = dataIN.variables['z0k'][iTime,lati_dnstream[il],lon_dnstream[il]]

    
# Close output file
dataOUT.close()
##########################################################################################
# END PROGRAM
##########################################################################################

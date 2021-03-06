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

# Output location
dirOUT  = "/data/dswales/NA-CORDEX/ARdet/dryingRatios/"

# What composites to look at?
# WRF horizontal resolution (50 or 25)
res         = 50
# Event strength threshold (percentile)
ptile       = ["96.00","97.00","98.00","99.00"]
nps         = len(ptile)
# Event persistence threshold (hours)
persistence = 36
# Model ID
#modelID     = ["gfdl","hadgem","mpi"]
modelID = ["mpi"]
nmodels     = len(modelID)

# PDF Configuration
nBins = 21   # Number of drying-ratio/precipitation bins
xMax  = 0.6  # Drying-ratio (max)
xMin  = 0.0  # Drying-ratio (min)
yMax  = 1.   # Precipitation (max)
yMin  = 0.   # Precipitation (min)

# Where to compute drying ratios?
# 1) Define a line ([x0,y0],[xF,yF]), number of points (nDRs), and a distance downstream
#    (d_dnstream; in Degrees) from that line
#x0         = -123.25
#xF         = -122.50
#y0         =  44.50
#yF         =  47.50
#nDRs       = 20
#d_dnstream = 2.0
## Assign some name for configuration, do this to avoid overwriting previous results.
#configID = "Take1"
x0         = -123.50
xF         = -122.50
y0         =  44.50
yF         =  47.50
nDRs       = 50
d_dnstream = 2.0
# Assign some name for configuration, do this to avoid overwriting previous results.
configID = "Take4"
x0d        = x0 + d_dnstream
xFd        = xF + d_dnstream
# Compute array of starting/ending points
lon_upstream  = [i for i in np.arange(x0, xF, ((xF-x0)/(nDRs)))]
lat_upstream  = [i for i in np.arange(y0, yF, ((yF-y0)/(nDRs)))]
lon_dnstream  = [i for i in np.arange(x0d,xFd,((xFd-x0d)/(nDRs)))]
lat_dnstream  = [i for i in np.arange(y0, yF, ((yF-y0)/(nDRs)))]

# OR 2) Define a set of points (if you just want to manually enter a few points)
#lon_upstream = [-110.1, -120.7]
#lat_upstream = [  45.0,   46.0]
#lon_dnstream = [-109.1, -119.7]
#lat_dnstream = [  45.5,   46.5]
#nDRs         = len(lon_upstream)

# Directory containing event composites.
dir = '/data/dswales/NA-CORDEX/ARdet/composites/'

for imod in range(0,nmodels):
    print(modelID[imod])
    for ip in range(0,nps):
        # Where to output results?
        fileOUT = dirOUT+"dryingRatios."+configID+"."+str(res)+'km.'+ptile[ip]+'ptile.'+str(persistence)+'hrs.'+modelID[imod]+".nc"

        ##########################################################################################
        # i) Read in metadata.
        #  - Read in grid and determine WRF gridpoint nearest to requested drying-ratio calculation
        ##########################################################################################
        # Get filenames
        fileH_precip = find('composite.raw.precip.events.'+str(res)+'km.'+ptile[ip]+'ptile.'+str(persistence)+'hrs.'+modelID[imod]+'.nc', dir)
        fileH_ivt    = find('composite.raw.ivt.events.'+str(res)+'km.'+ptile[ip]+'ptile.'+str(persistence)+'hrs.'+modelID[imod]+'.nc', dir)
        fileH_z0k    = find('composite.raw.z0k.events.'+str(res)+'km.'+ptile[ip]+'ptile.'+str(persistence)+'hrs.'+modelID[imod]+'.nc', dir)
        fileF_precip = find('composite.raw.precip.events.future.'+str(res)+'km.'+ptile[ip]+'ptile.'+str(persistence)+'hrs.'+modelID[imod]+'.nc', dir)
        fileF_ivt    = find('composite.raw.ivt.events.future.'+str(res)+'km.'+ptile[ip]+'ptile.'+str(persistence)+'hrs.'+modelID[imod]+'.nc', dir)
        fileF_z0k    = find('composite.raw.z0k.events.future.'+str(res)+'km.'+ptile[ip]+'ptile.'+str(persistence)+'hrs.'+modelID[imod]+'.nc', dir)
        
        # Open file
        dataIN1 = netCDF4.Dataset(fileH_precip[0],'r')
        dataIN2 = netCDF4.Dataset(fileH_ivt[0],'r')
        dataIN3 = netCDF4.Dataset(fileH_z0k[0],'r')
        lon    = dataIN1.variables['lon'][:]
        lat    = dataIN1.variables['lat'][:]
        yearH  = dataIN1.variables['year'][:]
        monthH = dataIN1.variables['month'][:]
        dayH   = dataIN1.variables['day'][:]
        hourH  = dataIN1.variables['hour'][:]
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
            # Also store indices for points neighboring the transect from [upstrem,downstream]
            if (il==0):
                dyTransect = int(lati_dnstream[il])-int(lati_upstream[il])+1
                dxTransect = int(loni_upstream[il])-int(loni_dnstream[il])+1
                nPtsAlongTransect = dxTransect*dyTransect
                loni_transect = np.zeros([nDRs,nPtsAlongTransect])
                lati_transect = np.zeros([nDRs,nPtsAlongTransect])
            count = 0
            for ix in range(0,dxTransect):
                for iy in range(0,dyTransect):
                    loni_transect[il,count] = loni_upstream[il]-ix
                    lati_transect[il,count] = lati_upstream[il]+iy
                    count = count+1

        ##########################################################################################
        # ii) Create output file
        ##########################################################################################
        dataOUT = netCDF4.Dataset(fileOUT,'w',format='NETCDF4_CLASSIC')
        dataOUT.createDimension("case",nDRs)
        dataOUT.createDimension("timeH",ntimeH)
        dataOUT.createDimension("xTransect",nPtsAlongTransect)
        dataOUT.createDimension("nXpdf",nBins-1)
        dataOUT.createDimension("nYpdf",nBins-1)
        dataOUT.createDimension("pdfBins",nBins)
        lonup_out        = dataOUT.createVariable("lon_upstream",     "f4", ("case"                    ))
        londn_out        = dataOUT.createVariable("lon_dnstream",     "f4", ("case"                    ))
        latup_out        = dataOUT.createVariable("lat_upstream",     "f4", ("case"                    ))
        latdn_out        = dataOUT.createVariable("lat_dnstream",     "f4", ("case"                    ))
        yearH_out        = dataOUT.createVariable("yearH",            "i4", (                   "timeH"))
        monthH_out       = dataOUT.createVariable("monthH",           "i4", (                   "timeH"))
        dayH_out         = dataOUT.createVariable("dayH",             "i4", (                   "timeH"))
        hourH_out        = dataOUT.createVariable("hourH",            "i4", (                   "timeH"))
        precip_transectH = dataOUT.createVariable("precip_transectH", "f4", ("case","xTransect","timeH"))
        ivt_transectH    = dataOUT.createVariable("ivt_transectH",    "f4", ("case","xTransect","timeH"))
        z0k_transectH    = dataOUT.createVariable("z0k_transectH",    "f4", ("case","xTransect","timeH"))
        drH              = dataOUT.createVariable("drH",              "f4", ("case",            "timeH"))
        xpdf_center      = dataOUT.createVariable("xpdf_center",      "f4", ("nXpdf"                   ))
        ypdf_center      = dataOUT.createVariable("ypdf_center",      "f4", (        "nYpdf"           ))
        xpdf_edge        = dataOUT.createVariable("xpdf_edge",        "f4", ("pdfBins"                 ))
        ypdf_edge        = dataOUT.createVariable("ypdf_edge",        "f4", (        "pdfBins"         ))
        pdfH             = dataOUT.createVariable("pdfH",             "f4", ("nYpdf","nXpdf"           ))
        lonup_out[:]     = lon_upstream[:]
        latup_out[:]     = lat_upstream[:]
        londn_out[:]     = lon_dnstream[:]
        latdn_out[:]     = lat_dnstream[:]
            
        ##########################################################################################
        # iii) Now that we know the indices of the points we are interested in, readin in data at
        #      these points, compute drying-ratio, and write output to netCDF.
        ##########################################################################################
        # i) Historical period.
        yearH_out[:]  = dataIN1.variables['year' ][:]
        monthH_out[:] = dataIN1.variables['month'][:]
        dayH_out[:]   = dataIN1.variables['day'  ][:]
        hourH_out[:]  = dataIN1.variables['hour' ][:]
        for il in range(0,nDRs):
            for im in range(0,nPtsAlongTransect):
                ivt_transectH[il,im,:]    = dataIN2.variables['ivt'    ][:,lati_transect[il,im],loni_transect[il,im]]
                precip_transectH[il,im,:] = dataIN1.variables['precip' ][:,lati_transect[il,im],loni_transect[il,im]]
                z0k_transectH[il,im,:]    = dataIN3.variables['z0k'    ][:,lati_transect[il,im],loni_transect[il,im]]
        drHa     = 1-ivt_transectH[:,0,:]/ivt_transectH[:,nPtsAlongTransect-1,:]
        drH[:,:] = drHa

        # ii) Future period (only GCMs)
        if (modelID[imod] != 'erain'):
            # Read in future file metadata
            dataIN1 = netCDF4.Dataset(fileF_precip[0],'r')
            dataIN2 = netCDF4.Dataset(fileF_ivt[0],'r')
            dataIN3 = netCDF4.Dataset(fileF_z0k[0],'r')
            yearF  = dataIN1.variables['year'][:]
            monthF = dataIN1.variables['month'][:]
            dayF   = dataIN1.variables['day'][:]
            hourF  = dataIN1.variables['hour'][:]
            ntimeF = len(yearF)
            
            # Add future fields to output file
            dataOUT.createDimension("timeF",ntimeF)
            yearF_out        = dataOUT.createVariable("yearF",            "i4", (                   "timeF"))
            monthF_out       = dataOUT.createVariable("monthF",           "i4", (                   "timeF"))
            dayF_out         = dataOUT.createVariable("dayF",             "i4", (                   "timeF"))
            hourF_out        = dataOUT.createVariable("hourF",            "i4", (                   "timeF"))
            precip_transectF = dataOUT.createVariable("precip_transectF", "f4", ("case","xTransect","timeF"))
            ivt_transectF    = dataOUT.createVariable("ivt_transectF",    "f4", ("case","xTransect","timeF"))
            z0k_transectF    = dataOUT.createVariable("z0k_transectF",    "f4", ("case","xTransect","timeF"))
            drF              = dataOUT.createVariable("drF",              "f4", ("case",            "timeF"))
            pdfF             = dataOUT.createVariable("pdfF",             "f4", ("nYpdf","nXpdf"           ))
                
            # Read in data
            yearF_out[:]  = dataIN1.variables['year' ][:]
            monthF_out[:] = dataIN1.variables['month'][:]
            dayF_out[:]   = dataIN1.variables['day'  ][:]
            hourF_out[:]  = dataIN1.variables['hour' ][:]
            for il in range(0,nDRs):
                for im in range(0,nPtsAlongTransect):
                    ivt_transectF[il,im,:]    = dataIN2.variables['ivt'    ][:,lati_transect[il,im],loni_transect[il,im]]
                    precip_transectF[il,im,:] = dataIN1.variables['precip' ][:,lati_transect[il,im],loni_transect[il,im]]
                    z0k_transectF[il,im,:]    = dataIN3.variables['z0k'    ][:,lati_transect[il,im],loni_transect[il,im]]            
            drFa = 1-ivt_transectF[:,0,:]/ivt_transectF[:,nPtsAlongTransect-1,:]
            drF[:,:]=drFa
    
        ##########################################################################################
        # iv) Compute and output joint-PDFs of statistics.
        ##########################################################################################
        # Which DRs have precipitation observed along the transect?
        precipTotTransH = np.sum(precip_transectH,axis=1)
        dryTransectsH = np.argwhere((precipTotTransH != 0) & (drHa <= 0))
        wetTransectsH = np.argwhere((precipTotTransH != 0) & (drHa > 0))
        nWetH = len(wetTransectsH)
        nDryH = len(dryTransectsH)
        if (modelID[imod] != 'erain'):
            precipTotTransF = np.sum(precip_transectF,axis=1)
            dryTransectsF = np.argwhere((precipTotTransF != 0) & (drFa <= 0))
            wetTransectsF = np.argwhere((precipTotTransF != 0) & (drFa > 0))
            nWetF = len(wetTransectsF)
            nDryF = len(dryTransectsF)

        # Compute histogram bin-boundaries
        dX   = (xMax-xMin)/(float(nBins)-1)
        dY   = (yMax-yMin)/(float(nBins)-1)
        xaxis = [i for i in np.arange(xMin, xMax+dX, dX)]
        yaxis = [i for i in np.arange(yMin, yMax+dY, dY)]
        
        # Historical 
        [pdf2DH,xi1,yi1]=np.histogram2d(drHa[wetTransectsH[:,0],wetTransectsH[:,1]], \
                                        precipTotTransH[wetTransectsH[:,0],wetTransectsH[:,1]], \
                                        bins = [xaxis,yaxis])
        pdfH[:,:] = pdf2DH#/sum(pdf2DH)
        
        # Future
        if (modelID[imod] != 'erain'):
            [pdf2DF,xi2,yi2]=np.histogram2d(drFa[wetTransectsF[:,0],wetTransectsF[:,1]], \
                                            precipTotTransF[wetTransectsF[:,0],wetTransectsF[:,1]], \
                                            bins = [xaxis,yaxis])
        pdfF[:,:] = pdf2DF#/sum(pdf2DF)
    
        # Store bin centers/edges
        xpdf_center[:] = (xi1[1:nBins]+xi1[0:nBins-1])/2.
        ypdf_center[:] = (yi1[1:nBins]+yi1[0:nBins-1])/2.
        xpdf_edge[:]   = xi1
        ypdf_edge[:]   = yi1

        # Close output file
        dataOUT.close()
        
    
##########################################################################################
# END PROGRAM
##########################################################################################

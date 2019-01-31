##########################################################################################
##########################################################################################
import os,netCDF4,numpy as np
import matplotlib.pyplot as plt
from pylab import figure

# What event configuration to look at
configID    = "Take1"
res         = 50
ptile       = '96.00'
persistence = 24
modelID     = 'gfdl'

# Input file
fileIN ='/data/dswales/NA-CORDEX/ARdet/dryingRatios/dryingRatios.'+configID+'.'+\
    str(res)+'km.'+ptile+'ptile.'+str(persistence)+'hrs.'+modelID+'.nc'

# Read in data
dataIN  = netCDF4.Dataset(fileIN,'r')
yearH   = dataIN.variables['yearH'][:]
monthH  = dataIN.variables['monthH'][:]
dayH    = dataIN.variables['dayH'][:]
hourH   = dataIN.variables['hourH'][:]
precipH = dataIN.variables['precip_transectH'][:]
ivtH    = dataIN.variables['ivt_transectH'][:]
z0kH    = dataIN.variables['z0k_transectH'][:]
ntransectH = len(precipH[0,:,0])
if (modelID != 'erain'):
    yearF   = dataIN.variables['yearF'][:]
    monthF  = dataIN.variables['monthF'][:]
    dayF    = dataIN.variables['dayF'][:]
    hourF   = dataIN.variables['hourF'][:]
    precipF = dataIN.variables['precip_transectF'][:]
    ivtF    = dataIN.variables['ivt_transectF'][:]
    z0kF    = dataIN.variables['z0k_transectF'][:]
    ntransectF = len(precipF[0,:,0])

# Compute drying-ratios
drH = 1-ivtH[:,0,:]/ivtH[:,ntransectH-1,:]
if (modelID != 'erain'):
    drF = 1-ivtF[:,0,:]/ivtF[:,ntransectF-1,:]

# Which DRs have precipitation observed along the transect?
precipTotTransH = np.sum(precipH,axis=1)
dryTransectsH = np.argwhere((precipTotTransH != 0) & (drH <= 0))
wetTransectsH = np.argwhere((precipTotTransH != 0) & (drH > 0))
nWetH = len(wetTransectsH)
nDryH = len(dryTransectsH)
if (modelID != 'erain'):
    precipTotTransF = np.sum(precipF,axis=1)
    dryTransectsF = np.argwhere((precipTotTransF != 0) & (drF <= 0))
    wetTransectsF = np.argwhere((precipTotTransF != 0) & (drF > 0))
    nWetH = len(wetTransectsF)
    nDryH = len(dryTransectsF)

# Make plot(s)
fig = figure(0, figsize=(16,10), dpi=80, facecolor='w', edgecolor='k')
#
if(modelID == 'erain'): plt.subplot(211)
if(modelID != 'erain'): plt.subplot(221)
plt.plot(drH[wetTransectsH[:,0],wetTransectsH[:,1]],precipTotTransH[wetTransectsH[:,0],wetTransectsH[:,1]],'*')
plt.axis([-1,1,0,100])
plt.ylabel('Precipitation (mm)')
plt.title('Precipitation w/ DR>0 (westerly-flow regime) (historical)')
#
if(modelID == 'erain'): plt.subplot(212)
if(modelID != 'erain'): plt.subplot(223)
plt.plot(drH[dryTransectsH[:,0],dryTransectsH[:,1]],precipTotTransH[dryTransectsH[:,0],dryTransectsH[:,1]],'*')
plt.axis([-1,1,0,100])
plt.ylabel('Precipitation (mm)')
plt.xlabel('Drying-Ratio')
plt.title('Precipitation w/ DR<0 (easterly-flow regime?) (historical)')
if (modelID != 'erain'):
    plt.subplot(222)
    plt.plot(drF[wetTransectsF[:,0],wetTransectsF[:,1]],precipTotTransF[wetTransectsF[:,0],wetTransectsF[:,1]],'r*')
    plt.axis([-1,1,0,100])
    plt.ylabel('Precipitation (mm)')
    plt.title('Precipitation w/ DR>0 (westerly-flow regime) (2070-2100)')
    #
    plt.subplot(224)
    plt.plot(drF[dryTransectsF[:,0],dryTransectsF[:,1]],precipTotTransF[dryTransectsF[:,0],dryTransectsF[:,1]],'r*')
    plt.axis([-1,1,0,100])
    plt.ylabel('Precipitation (mm)')
    plt.xlabel('Drying-Ratio')
    plt.title('Precipitation w/ DR<0 (easterly-flow regime?) (2070-2100)')
#
plt.show()

##########################################################################################
# END PROGRAM
##########################################################################################

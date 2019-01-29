; ##################################################################### 
;
; The purpose of this program is to find the nearest WRF gridpoints
; for a set on longitude and latitude points, loni and lati.
;
; These gridpoints are stored in a text file for later use by
; pullOutSlab.pro
;
; ##################################################################### 
@find_nearestWRFpt.pro

loni = [-127.00, -126.25, -125.50, -124.75, -124.60, -124.40, -124.20,$
        -124.00, -124.00, -124.10, -124.15, -124.00, -124.15, -124.25,$
        -124.25, -124.25, -124.25, -124.00, -124.25, -124.00, -124.00,$
        -124.00, -123.75, -123.50, -123.00, -122.50, -122.00, -122.00,$
        -121.50, -121.00, -120.50, -120.50, -119.00, -118.00, -117.00,$
        -117.00, -117.00, -116.50]
lati = 50-indgen(n_elements(loni))*0.5
npts = n_elements(loni)

; WRF grid resolution
res = 50

; Configuration ID
configID = 'coast1.'+string(res,format='(i2)')

; Read in WRF grid from NA-CORDEX experiments
fileIN = '/Projects/HydroMet/dswales/NA-CORDEX/WRF_erain_'+string(res,format='(i2)')+'km/wrfout_198001.nc'
fileID = ncdf_open(fileIN)
ncdf_varget,fileID,ncdf_varid(fileID,'XLONG'),lon
ncdf_varget,fileID,ncdf_varid(fileID,'XLAT'),lat
ncdf_close,fileID

xi = intarr(npts)
yi = intarr(npts)
for ij=0,npts-1 do begin
   find_nearestWRFpt,lon,lat,loni(ij),lati(ij),xii,yii
   xi(ij) = xii
   yi(ij) = yii
end

; plot
cgmap_set,/continents,/usa,limit=[30,-130,55,-105],/grid,/isotropic
for ij=0,npts-1 do begin
   cgplot,lon(xi(ij),yi(ij)),lat(xi(ij),yi(ij)),psym=2,/overplot,color='red'
end
saveimage,'plots/coastalSlabXY.'+configID+'.png'

openw,101,'data/coastalSlabXY.'+configID+'.txt'
printf,101,n_elements(xi)
printf,101,xi
printf,101,yi
close,101

; END PROGRAM
end

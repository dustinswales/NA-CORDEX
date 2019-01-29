; #####################################################################
;
; The purpose of this program is to pull out the nearest WRF
; gridpoint for a given location.
;
; Name    Type    Size    Dscription
;
; Inputs
; lon     float   [:,:]   gridded longitude array.
; lat     float   [:,:]   gridded latitude array.
; xi      float   [1]     Longitude of gridpoint.
; yi      float   [1]     Latitude of gridpoint.
;
; Outputs
; xii     int     [1]     Index into 'lon' for nearest longitude.
; yii     int     [1]     Index into 'lat' for nearest latitude.
;
; #####################################################################
pro find_nearestWRFpt,lon,lat,xi,yi,xii,yii

; Compute distance from each gridpoint and find the smallest distance.
d = fltarr(n_elements(lon(*,0)),n_elements(lat(0,*)))
for ij=0,n_elements(lon(*,0))-1 do begin
   for ik=0,n_elements(lat(0,*))-1 do begin
      d(ij,ik) = MAP_2POINTS(xi, yi, lon(ij,ik), lat(ij,ik),/MILES)
   end
end
xyi=array_indices(lon,where(d eq min(d)))
xii = xyi(0)
yii = xyi(1)

; END PROGRAM
end

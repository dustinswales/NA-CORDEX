; This program is to play with drying ratio locations.

; Plotting setup
y0 = 30                         ; min lat
x0 = -125                       ; min lon
y1 = 50                         ; max lat
x1 = -105                       ; max lon
nlevels=32                      ; # of levels in plot
c_colors=floor((32.0/nlevels)*(indgen(nlevels)))


; Read in CFSR topography data
file='/data/jscott/gtopo30/topo.global.30sec.nc'
fileID=ncdf_open(file)
ncdf_varget,fileID,ncdf_varid(fileID,'lat'),topo_lat,stride=20
ncdf_varget,fileID,ncdf_varid(fileID,'lon'),topo_lon,stride=20
ncdf_varget,fileID,ncdf_varid(fileID,'elev'),topo_z,stride=[20,20]
ncdf_close,fileID

; Subset topography
xxx=where(topo_lon ge x0+360 and topo_lon le x1+360)
yyy=where(topo_lat ge y0 and topo_lat le y1)
topo_lonS=topo_lon(xxx)
topo_latS=topo_lat(yyy)
topo_zS=fltarr(n_elements(xxx),n_elements(yyy))
for ij=0,n_elements(xxx)-1 do begin
   for ik=0,n_elements(yyy)-1 do begin
      topo_zS(ij,ik)=topo_z(xxx(ij),yyy(ik))
   end
end

; Read in Jamies topography color table
openr,101,'/home/jscott/ncl/colormaps/land.usgs.rgb'
cin=intarr(3,32)
readf,101,cin
close,101
colortable=intarr(3,256)
colortable(*,0:31)=cin
colortable(*,32)=0
colortable(*,0)=255
r=reform(colortable(0,*)) & g=reform(colortable(1,*)) & b=reform(colortable(2,*))
tvlct,r,g,b

; 1) Define a line ([xS,yS],[xF,yF]), number of points (nDRs), and a distance downstream
;    (d_dnstream; in Degrees) from that line
xS         = -123.50
xF         = -122.50
yS         =  44.50
yF         =  47.50
nDRs       = 50
d_dnstream = 2.0
xSd        = xS + d_dnstream
xFd        = xF + d_dnstream
;
lonup = xS+findgen(nDRs)*((xF-xS)/(nDRs-1))
londn = xSd+findgen(nDRs)*((xFd-xSd)/(nDRs-1))
latup = yS+findgen(nDRs)*((yF-yS)/(nDRs-1))
latdn = yS+findgen(nDRs)*((yF-yS)/(nDRs-1))


; plot
margin=0.10
position=aspect(1.0,margin=margin)+margin/2.0
!p.background=0
map_set,/continents,/us,limit=[40,-125,50,-115],/grid,color=32,$
  latlab=x0,lonlab=y0,pos=position,/isotropic,/noborder
topo_zs(where(topo_zs lt 0 and topo_zs gt -500))=0.0
cgcontour,smooth(topo_zs(0:*,0:*),0),topo_lons(0:*,0:*),topo_lats(0:*,0:*),$
          /overplot,/fill,c_colors=c_colors,nlevels=nlevels
map_set,/continents,/us,limit=[40,-125,50,-115],/grid,color=32,/isotropic,/noborder,$
  latlab=x0,lonlab=y0,pos=position,/noerase,charsize=2,mlinethick=2
cgplot,lonup,latup,psym=2,/overplot
cgplot,londn,latdn,psym=2,/overplot
range=[min(topo_zs(where(topo_zs ne -500))),max(topo_zs)]
cgcolorbar,position=[position(0),margin/2.0,position(2),margin],bottom=0,ncolors=max(c_colors),color=32,$
  range=range


; END PROGRAM
end

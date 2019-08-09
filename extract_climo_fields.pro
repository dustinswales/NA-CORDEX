; #####################################################################
;
; The purpose of this program is to compute the climo for fields from
; the WRF NA-CORDEX. Currently precipitation and lifting condensation
; height (zLCL) are the only fields computed.
;
; #####################################################################
pro extract_climo_fields,modelID,yearS,yearF
  
; Physical constants
rv      = 461.50
rd      = 287.04
lv      = 2.501e6
cp      = 1005.0
cv      = 718.0
epsilon = rd/rv
  
; 50 or 25km data?
res = '50'
; Cool-season months only (JFM.OND)
cool=1

if (yearS lt 2010 and yearF le 2010) then period='clim' else period='future'

; WRF data
expID = 'WRF_'+modelID+'_'+res+'km'
dir   = '/Projects/HydroMet/dswales/NA-CORDEX/'+expID+'/'

; Output location
dirOUT = '/Projects/HydroMet/dswales/NA-CORDEX/analysis/climo/'

; Get files to read in
files = file_search(dir,'wrfout_*.nc')
files = files(where(strpos(files,'mshf') eq -1))
nfiles=n_elements(files)
mon = intarr(nfiles)
yy  = intarr(nfiles)
for ij=0,nfiles-1 do begin
   mon(ij) = fix(strmid(files(ij),strlen(dir)+11,2))
   yy(ij)  = fix(strmid(files(ij),strlen(dir)+7,4))
end
cool_mon = where((mon le 3 or mon ge 10))
if cool then begin
   files=files(cool_mon)
   mon = mon(cool_mon)
   yy  = yy(cool_mon)
endif
nfiles=n_elements(files)

; Subset in time
files = files(where(yy ge yearS and yy le yearF))
nfiles=n_elements(files)

; Read in grid description.
fileID=ncdf_open(files(0))
ncdf_varget,fileID,ncdf_varid(fileID,'XLONG'),lon
ncdf_varget,fileID,ncdf_varid(fileID,'XLAT'),lat
nlon = n_elements(lon(*,0))
nlat = n_elements(lon(0,*))
ncdf_close,fileID

; Read in each latitude slice at a time, save time/memory
print,'Readin in data...'
var_temp1=fltarr(nlon,nlat)
var_temp1a=fltarr(nlon,nlat)
t0 = 0
for ij=0,nfiles-1 do begin
   print,ij,'of',nfiles
   fileID=ncdf_open(files(ij))
   ncdf_varget,fileID,ncdf_varid(fileID,'Day'),dayA
   ncdf_varget,fileID,ncdf_varid(fileID,'Hour'),hourA
   ncdf_varget,fileID,ncdf_varid(fileID,'Month'),monthA
   ncdf_varget,fileID,ncdf_varid(fileID,'Year'),yearA
   nt = n_elements(dayA)
   ;
   var_temp2 = fltarr(nlon,nlat,nt)
   ncdf_varget,fileID,ncdf_varid(fileID,'RAINNC'),rainncA
   ncdf_varget,fileID,ncdf_varid(fileID,'RAINC'),raincA
   rainncA = rainncA+raincA
   var_temp2(*,*,0)=var_temp1
   var_temp2(*,*,1:nt-1)=rainncA(*,*,1:nt-1)-rainncA(*,*,0:nt-2)
   var_temp2(where(var_temp2 lt 0)) = 0
   var_temp1 = var_temp2(*,*,nt-1)
   ;
   ncdf_varget,fileID,ncdf_varid(fileID,'IVTU'), ivtU
   ncdf_varget,fileID,ncdf_varid(fileID,'IVTV'), ivtV
   ivtA = sqrt(ivtU*ivtU+ivtV*ivtV)
   ncdf_varget,fileID,ncdf_varid(fileID,'Z0K'), z0ka
   ncdf_varget,fileID,ncdf_varid(fileID,'Z500'), z500a
   ncdf_varget,fileID,ncdf_varid(fileID,'SNOWE'), swea
   ;
   ncdf_varget,fileID,ncdf_varid(fileID,'PSFC'),sfcP
   ncdf_varget,fileID,ncdf_varid(fileID,'Q2m'),q2m
   ncdf_varget,fileID,ncdf_varid(fileID,'T2m'),t2m
   ; Compute saturation vapor pressure (es) and vapor pressure (pe)
   es = 610.94*exp(lv*(1./273.15 - 1./t2m)/rv)
   pe = sfcP/((epsilon/q2m)+1)
   pe(where(pe lt 0))=0.0
   ; Compute dew-point temperature
   Td = 1./((1./273.15) - (rv/lv)*alog(pe/610.94))
   ; Compute heigt lifting condensation level
   zLCLa=(t2m-Td)/8.
   zLCLa(where(zLCLa lt 0))=0.0
   ncdf_close,fileID
   rainnca=reform(temporary(rainnca))
         
   ; Store data
   rainnc = fltarr(nlon,nlat,nt)
   rainnc(*,*,1:nt-1) = rainncA(*,*,1:nt-1)-rainncA(*,*,0:nt-2)
   ivt          = ivta
   swe          = swea
   z500         = z500a
   zLCL         = zLCLa
   day          = dayA
   month        = monthA
   year         = yearA
   hour         = hourA

   ; Write output

   ; 1) Precipitation
   if (ij eq 0) then begin
      if (strcmp(period,'clim')) then $
         fileOUT = 'climo.rainnc.'+expID+'.3hr.nc'
      if (strcmp(period,'future')) then $
         fileOUT = 'climo.future.rainnc.'+expID+'.3hr.nc'
      
      fileID1  = ncdf_create(dirOUT+fileOUT,/clobber)
      dimID1  = ncdf_dimdef(fileID1,'lon',nlon)
      dimID2  = ncdf_dimdef(fileID1,'lat',nlat)
      dimID4  = ncdf_dimdef(fileID1,'ntime',/unlimited)
      varID1  = ncdf_vardef(fileID1,'lon',    [dimID1,dimID2],/float)
      varID2  = ncdf_vardef(fileID1,'lat',    [dimID1,dimID2],/float)
      varID3  = ncdf_vardef(fileID1,'rainnc', [dimID1,dimID2,dimID4],/float)
      ncdf_control,fileID1,/endef
      ncdf_varput,fileID1,varID1,lon
      ncdf_varput,fileID1,varID2,lat
   endif
   ncdf_varput,fileID1,varID3,rainnc,offset=[0,0,t0]

   ; 2) IVT
   if (ij eq 0) then begin
      if (strcmp(period,'clim')) then $
         fileOUT = 'climo.ivt.'+expID+'.3hr.nc'
      if (strcmp(period,'future')) then $
         fileOUT = 'climo.future.ivt.'+expID+'.3hr.nc'
      
      fileID2  = ncdf_create(dirOUT+fileOUT,/clobber)
      dimID1  = ncdf_dimdef(fileID2,'lon',nlon)
      dimID2  = ncdf_dimdef(fileID2,'lat',nlat)
      dimID4  = ncdf_dimdef(fileID2,'ntime',/unlimited)
      varID1  = ncdf_vardef(fileID2,'lon', [dimID1,dimID2],       /float)
      varID2  = ncdf_vardef(fileID2,'lat', [dimID1,dimID2],       /float)
      varID3  = ncdf_vardef(fileID2,'ivt', [dimID1,dimID2,dimID4],/float)
      ncdf_control,fileID2,/endef
      ncdf_varput,fileID2,varID1,lon
      ncdf_varput,fileID2,varID2,lat
   endif
   ncdf_varput,fileID2,varID3,ivt,offset=[0,0,t0]

   ; 3) zLCL
   if (ij eq 0) then begin
      if (strcmp(period,'clim')) then $
         fileOUT = 'climo.zLCL.'+expID+'.3hr.nc'
      if (strcmp(period,'future')) then $
         fileOUT = 'climo.future.zLCL.'+expID+'.3hr.nc'
      
      fileid3  = ncdf_create(dirOUT+fileOUT,/clobber)
      dimID1  = ncdf_dimdef(fileid3,'lon',nlon)
      dimID2  = ncdf_dimdef(fileid3,'lat',nlat)
      dimID4  = ncdf_dimdef(fileid3,'ntime',/unlimited)
      varID1  = ncdf_vardef(fileid3,'lon',  [dimID1,dimID2],       /float)
      varID2  = ncdf_vardef(fileid3,'lat',  [dimID1,dimID2],       /float)
      varID3  = ncdf_vardef(fileid3,'zLCL', [dimID1,dimID2,dimID4],/float)
      ncdf_control,fileid3,/endef
      ncdf_varput,fileid3,varID1,lon
      ncdf_varput,fileid3,varID2,lat
   endif
   ncdf_varput,fileid3,varID3,zLCL,offset=[0,0,t0]

   ; 4) z500
   if (ij eq 0) then begin
      if (strcmp(period,'clim')) then $
         fileOUT = 'climo.z500.'+expID+'.3hr.nc'
      if (strcmp(period,'future')) then $
         fileOUT = 'climo.future.z500.'+expID+'.3hr.nc'

      fileID4  = ncdf_create(dirOUT+fileOUT,/clobber)
      dimID1  = ncdf_dimdef(fileID4,'lon',nlon)
      dimID2  = ncdf_dimdef(fileID4,'lat',nlat)
      dimID4  = ncdf_dimdef(fileID4,'ntime',/unlimited)
      varID1  = ncdf_vardef(fileID4,'lon',  [dimID1,dimID2],       /float)
      varID2  = ncdf_vardef(fileID4,'lat',  [dimID1,dimID2],       /float)
      varID3  = ncdf_vardef(fileID4,'z500', [dimID1,dimID2,dimID4],/float)
      ncdf_control,fileID4,/endef
      ncdf_varput,fileID4,varID1,lon
      ncdf_varput,fileID4,varID2,lat
   endif
   ncdf_varput,fileID4,varID3,z500,offset=[0,0,t0]

   ; 5) swe
   if (ij eq 0) then begin
      if (strcmp(period,'clim')) then $
         fileOUT = 'climo.swe.'+expID+'.3hr.nc'
      if (strcmp(period,'future')) then $
         fileOUT = 'climo.future.swe.'+expID+'.3hr.nc'
      
      fileID5  = ncdf_create(dirOUT+fileOUT,/clobber)
      dimID1  = ncdf_dimdef(fileID5,'lon',nlon)
      dimID2  = ncdf_dimdef(fileID5,'lat',nlat)
      dimID4  = ncdf_dimdef(fileID5,'ntime',/unlimited)
      varID1  = ncdf_vardef(fileID5,'lon',  [dimID1,dimID2],       /float)
      varID2  = ncdf_vardef(fileID5,'lat',  [dimID1,dimID2],       /float)
      varID3  = ncdf_vardef(fileID5,'swe',  [dimID1,dimID2,dimID4],/float)
      ncdf_control,fileID5,/endef
      ncdf_varput,fileID5,varID1,lon
      ncdf_varput,fileID5,varID2,lat
   endif
   ncdf_varput,fileID5,varID3,swe,offset=[0,0,t0]

   t0 = t0+n_elements(hour)
end

; Close netCDF files
ncdf_close,fileID1
ncdf_close,fileID2
ncdf_close,fileID3
ncdf_close,fileID4
ncdf_close,fileID5

print,'Finished'
; #####################################################################
; END PROGRAM
; #####################################################################
end

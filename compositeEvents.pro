; ######################################################################### 
;
; The purpose of this program is to pull out the WRF data for the 
; events found in find_ARs.pro
; Name                Type    Size    Dscription
;
; Inputs
; modelID             string  [1]     Shortname for model (e.g. 'gfdl')
; persistenceThresh   int     [1]     Persistence threshold (hours)
; ivtThresh           float   [1]     Strength threshold (percentile)
; dirOUT              string  [1]     Output directory location.
; dirIN               string  [1]     Input directory
; future              int     [1]     Set to true(1) if compositing
;                                     over future epoch.
;
; Output is written to netCDF and stored in dirOUT.
;
; ######################################################################### 
pro compositeEvents,dirIN,dirOUT,modelID,ivtThresh,res,persistenceThresh,future

doPrecip = 0
doIVT    = 1
doLCL    = 1
doz0k    = 0
doSWE    = 1
doZ500   = 1
  
; WRF data
expID = 'WRF_'+modelID+'_'+res+'km'
dir   = '/Projects/HydroMet/dswales/NA-CORDEX/'+expID+'/'

; For composites, how many hours after event to include?
lingerTime = 12

; Physical constants
rv      = 461.50
rd      = 287.04
lv      = 2.501e6
cp      = 1005.0
cv      = 718.0
epsilon = rd/rv

; Read in file containing detected AR events.
fprefix = 'events.'
if (future) then fprefix = 'events.future.'
fileIN = fprefix+res+'km.'+string(ivtThresh,format='(f5.2)')+'ptile.'+string(persistenceThresh,format='(i2.2)')+'hrs.'+modelID+'.nc'
fileID = ncdf_open(dirIN+fileIN)
ncdf_varget,fileID,ncdf_varid(fileID,'year_start'),yearS
ncdf_varget,fileID,ncdf_varid(fileID,'year_end'),yearF
ncdf_varget,fileID,ncdf_varid(fileID,'month_start'),monthS
ncdf_varget,fileID,ncdf_varid(fileID,'month_end'),monthF
ncdf_varget,fileID,ncdf_varid(fileID,'day_start'),dayS
ncdf_varget,fileID,ncdf_varid(fileID,'day_end'),dayF
ncdf_varget,fileID,ncdf_varid(fileID,'hour_start'),hourS
ncdf_varget,fileID,ncdf_varid(fileID,'hour_end'),hourF
ncdf_varget,fileID,ncdf_varid(fileID,'length'),event_length
ncdf_close,fileID
nEvents = n_elements(yearS)

; If requested, extend composite time from end of event.
if (lingerTime gt 0) then begin
   for ij=0,nEvents-1 do begin
      ti = julday(monthF(ij),dayF(ij),yearF(ij),hourF(ij))
      caldat,ti+lingerTime/24.,m,d,y,h
      yearF(ij)=y
      monthF(ij)=m
      dayF(ij)=d
      hourF(ij)=h
      event_length(ij)=event_length(ij)+lingerTime
   end
end

; Read in WRF data for each event and make plots.
init=1
for ij=0,nEvents-1 do begin
   print,string(monthS(ij),format='(i2.2)')+'/'+string(dayS(ij),format='(i2.2)')+'/'+string(yearS(ij),format='(i4)')+':'+string(hourS(ij),format='(i2.2)')+$
         'Z - '+string(monthF(ij),format='(i2.2)')+'/'+string(dayF(ij),format='(i2.2)')+'/'+string(yearF(ij),format='(i4)')+':'+string(hourF(ij),format='(i2.2)')+'Z'
   f0 = 'wrfout_'+string(yearS(ij),format='(i4)')+string(monthS(ij),format='(i2.2)')+'.nc'
   f1 = 'wrfout_'+string(yearF(ij),format='(i4)')+string(monthF(ij),format='(i2.2)')+'.nc'
   if (strcmp(f0,f1)) then begin
      ;print,"data in single file"
      fileID = ncdf_open(dir+f0)
      
      ; Determine what time slice to read in?
      ncdf_varget,fileID,ncdf_varid(fileID,'Day'),day
      ncdf_varget,fileID,ncdf_varid(fileID,'Hour'),hour
      ncdf_varget,fileID,ncdf_varid(fileID,'Month'),month
      ncdf_varget,fileID,ncdf_varid(fileID,'Year'),year
      if (ij eq 0) then begin
         ncdf_varget,fileID,ncdf_varid(fileID,'XLONG'),lon
         ncdf_varget,fileID,ncdf_varid(fileID,'XLAT'),lat
         nlon = n_elements(lon(*,0))
         nlat = n_elements(lat(0,*))
      endif
      t0 = where(Year eq YearS(ij) and Month eq MonthS(ij) and Day eq DayS(ij) and Hour eq HourS(ij))
      tF = where(Year eq YearF(ij) and Month eq MonthF(ij) and Day eq DayF(ij) and Hour eq HourF(ij))
      if (t0(0) eq -1 or n_elements(t0) gt 1 or tF(0) eq -1 or n_elements(tF) gt 1) then begin
         print,'ERROR: something went wrong finding start/end time in WRF dataset... Skip event'
         goto, badEvent
      endif
                                ; Read in fields
      if (doPrecip) then begin
         ncdf_varget,fileID,ncdf_varid(fileID,'RAINNC'),var_temp2,offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
         ncdf_varget,fileID,ncdf_varid(fileID,'RAINC'),var_temp2a,offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
         var_temp2 = var_temp2 + var_temp2a
         nt        = n_elements(var_temp2(0,0,*))
         rainnc    = fltarr(nlon,nlat,nt)
         rainnc(*,*,1:nt-1)=var_temp2(*,*,1:nt-1)-var_temp2(*,*,0:nt-2)
         rainnc(where(rainnc lt 0)) = 0
      endif
      if (doIVT)  then ncdf_varget,fileID,ncdf_varid(fileID,'IVTU'), ivtU, offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
      if (doIVT)  then ncdf_varget,fileID,ncdf_varid(fileID,'IVTV'), ivtV, offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
      if (doz0k)  then ncdf_varget,fileID,ncdf_varid(fileID,'Z0K'), z0k,   offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
      if (doz500) then ncdf_varget,fileID,ncdf_varid(fileID,'Z500'), z500, offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
      if (doswe)  then ncdf_varget,fileID,ncdf_varid(fileID,'SNOWE'), swe, offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
      if (doLCL)  then begin
         ncdf_varget,fileID,ncdf_varid(fileID,'PSFC'),sfcP,offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
         ncdf_varget,fileID,ncdf_varid(fileID,'Q2m'),  q2m,offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
         ncdf_varget,fileID,ncdf_varid(fileID,'T2m'),  t2m,offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
  
         ; Compute saturation vapor pressure (es) and vapor pressure (pe)
         es = 610.94*exp(lv*(1./273.15 - 1./t2m)/rv)
         pe = sfcP/((epsilon/q2m)+1)
         pe(where(pe lt 0))=0.0
         ; Compute dew-point temperature
         Td = 1./((1./273.15) - (rv/lv)*alog(pe/610.94))
         ; Compute heigt lifting condensation level
         zLCL=(t2m-Td)/8.
         zLCL(where(zLCL lt 0))=0.0
      endif
      
      ; Store stuff for plots
      yy = year(t0:t0+event_length(ij)/3-1)
      mm = month(t0:t0+event_length(ij)/3-1)
      dd = day(t0:t0+event_length(ij)/3-1)
      hh = hour(t0:t0+event_length(ij)/3-1)
      
      ; Close file
      ncdf_close,fileID
      
   endif
   if (not strcmp(f0,f1)) then begin
      ;print,"data NOT in single file"

      ; File 1
      fileID = ncdf_open(dir+f0)
      
      ; Determine what time slice to read in?
      ncdf_varget,fileID,ncdf_varid(fileID,'Day'),day
      ncdf_varget,fileID,ncdf_varid(fileID,'Hour'),hour
      ncdf_varget,fileID,ncdf_varid(fileID,'Month'),month
      ncdf_varget,fileID,ncdf_varid(fileID,'Year'),year
      ncdf_varget,fileID,ncdf_varid(fileID,'XLONG'),lon
      ncdf_varget,fileID,ncdf_varid(fileID,'XLAT'),lat
      nlon = n_elements(lon(*,0))
      nlat = n_elements(lat(0,*))
      t0 = where(Year eq YearS(ij) and Month eq MonthS(ij) and Day eq DayS(ij) and Hour eq HourS(ij))
      if (t0(0) eq -1) then begin
         print,'ERROR: something went wrong?'
         stop
      endif

      ; Read in fields
      nt0 = n_elements(day)-t0
      if (doPrecip) then begin
         ncdf_varget,fileID,ncdf_varid(fileID,'RAINNC'),rainnca,offset=[0,0,t0],count=[nlon,nlat,nt0]
         ncdf_varget,fileID,ncdf_varid(fileID,'RAINC'), rainca, offset=[0,0,t0],count=[nlon,nlat,nt0]
         rainnca = rainnca + rainca
         nt   = n_elements(rainnca(0,0,*))
         var_temp2a = fltarr(nlon,nlat,nt)
         if (nt gt 1) then begin
            var_temp2a(*,*,1:nt-1)=rainnca(*,*,1:nt-1)-rainnca(*,*,0:nt-2)
            var_temp2a(where(var_temp2a lt 0)) = 0
         endif
      endif
      if (doIVT)  then ncdf_varget,fileID,ncdf_varid(fileID,'IVTU'), ivtUa, offset=[0,0,t0],count=[nlon,nlat,nt0]
      if (doIVT)  then ncdf_varget,fileID,ncdf_varid(fileID,'IVTV'), ivtVa, offset=[0,0,t0],count=[nlon,nlat,nt0]
      if (doz0k)  then ncdf_varget,fileID,ncdf_varid(fileID,'Z0K'), z0ka,   offset=[0,0,t0],count=[nlon,nlat,nt0]
      if (doz500) then ncdf_varget,fileID,ncdf_varid(fileID,'Z500'), z500a, offset=[0,0,t0],count=[nlon,nlat,nt0]
      if (doswe)  then ncdf_varget,fileID,ncdf_varid(fileID,'SNOWE'), swea, offset=[0,0,t0],count=[nlon,nlat,nt0]
      if (doLCL)  then begin
         ncdf_varget,fileID,ncdf_varid(fileID,'PSFC'),sfcPa,offset=[0,0,t0],count=[nlon,nlat,nt0]
         ncdf_varget,fileID,ncdf_varid(fileID,'Q2m'), q2ma, offset=[0,0,t0],count=[nlon,nlat,nt0]
         ncdf_varget,fileID,ncdf_varid(fileID,'T2m') ,t2ma, offset=[0,0,t0],count=[nlon,nlat,nt0]
      
         ; Compute saturation vapor pressure (es) and vapor pressure (pe)
         es = 610.94*exp(lv*(1./273.15 - 1./t2ma)/rv)
         pe = sfcPa/((epsilon/q2ma)+1)
         pe(where(pe lt 0))=0.0
         ; Compute dew-point temperature
         Td = 1./((1./273.15) - (rv/lv)*alog(pe/610.94))
         ; Compute heigt lifting condensation level
         zLCLa=(t2ma-Td)/8.
         zLCLa(where(zLCLa lt 0))=0.0
      endif
      
      yya = year(t0:t0+nt0-1)
      mma = month(t0:t0+nt0-1)
      dda = day(t0:t0+nt0-1)
      hha = hour(t0:t0+nt0-1)
      
      ; Close file 1
      ncdf_close,fileID

      ; File 2
      fileID = ncdf_open(dir+f1)
      ncdf_varget,fileID,ncdf_varid(fileID,'Day'),day
      ncdf_varget,fileID,ncdf_varid(fileID,'Hour'),hour
      ncdf_varget,fileID,ncdf_varid(fileID,'Month'),month
      ncdf_varget,fileID,ncdf_varid(fileID,'Year'),year
      nt1 = event_length(ij)-nt0
      if (doPrecip) then begin
         ncdf_varget,fileID,ncdf_varid(fileID,'RAINNC'),rainncb,count=[nlon,nlat,nt1]
         ncdf_varget,fileID,ncdf_varid(fileID,'RAINC'), raincb, count=[nlon,nlat,nt1]
         rainncb = rainncb + raincb
         nt   = n_elements(rainncb(0,0,*))
         var_temp2b = fltarr(nlon,nlat,nt)
         var_temp2b(*,*,1:nt-1)=rainncB(*,*,1:nt-1)-rainncB(*,*,0:nt-2)
         var_temp2b(where(var_temp2b lt 0)) = 0
      endif
      if (doIVT)  then ncdf_varget,fileID,ncdf_varid(fileID,'IVTU'), ivtUb, count=[nlon,nlat,nt1]
      if (doIVT)  then ncdf_varget,fileID,ncdf_varid(fileID,'IVTV'), ivtVb, count=[nlon,nlat,nt1]
      if (doz0k)  then ncdf_varget,fileID,ncdf_varid(fileID,'Z0K'), z0kb,   count=[nlon,nlat,nt1]
      if (doz500) then ncdf_varget,fileID,ncdf_varid(fileID,'Z500'), z500b, count=[nlon,nlat,nt1] 
      if (doswe)  then ncdf_varget,fileID,ncdf_varid(fileID,'SNOWE'), sweb, count=[nlon,nlat,nt1]
      if (doLCL)  then begin
         ncdf_varget,fileID,ncdf_varid(fileID,'PSFC'),sfcPb,count=[nlon,nlat,nt1]
         ncdf_varget,fileID,ncdf_varid(fileID,'Q2m'), q2mb, count=[nlon,nlat,nt1]
         ncdf_varget,fileID,ncdf_varid(fileID,'T2m'), t2mb, count=[nlon,nlat,nt1]

         ; Compute saturation vapor pressure (es) and vapor pressure (pe)
         es = 610.94*exp(lv*(1./273.15 - 1./t2mb)/rv)
         pe = sfcPb/((epsilon/q2mb)+1)
         pe(where(pe lt 0))=0.0
         ; Compute dew-point temperature
         Td = 1./((1./273.15) - (rv/lv)*alog(pe/610.94))
         ; Compute heigt lifting condensation level
         zLCLb=(t2mb-Td)/8.
         zLCLb(where(zLCLb lt 0))=0.0
      endif
      
      yyb = year(0:nt1-1)
      mmb = month(0:nt1-1)
      ddb = day(0:nt1-1)
      hhb = hour(0:nt1-1)

      ; Close file 2
      ncdf_close,fileID
      
      ; Combine
      if (doLCL)    then zLCL   = [[[zLCLa]],[[zLCLb]]]
      if (doIVT)    then ivtU   = [[[ivtUa]],[[ivtUb]]]
      if (doIVT)    then ivtV   = [[[ivtVa]],[[ivtVb]]]
      if (doz0k)    then z0k    = [[[z0ka]],[[z0kb]]]
      if (doz500)   then z500   = [[[z500a]],[[z500b]]]
      if (doswe)    then swe    = [[[swea]],[[sweb]]]
      if (doPrecip) then rainnc = [[[var_temp2A]],[[var_temp2B]]]
      yy   = [yya,yyb]
      mm   = [mma,mmb]
      dd   = [dda,ddb]
      hh   = [hha,hhb]
   endif

   ; Store data from event(s).
   if (init) then begin
      event_hour   = hh 
      event_day    = dd
      event_month  = mm
      event_year   = yy
      if (doPrecip) then event_precip = rainnc
      if (doIVT)    then event_ivtU   = ivtU
      if (doIVT)    then event_ivtV   = ivtV
      if (doz0k)    then event_z0k    = z0k
      if (doLCL)    then event_zLCL   = zLCL
      if (doZ500)   then event_z500   = z500
      if (doSWE)    then event_swe    = swe
   endif
   if (not init) then begin
      event_hour   = [event_hour, hh]
      event_day    = [event_day,  dd]
      event_month  = [event_month,mm]
      event_year   = [event_year, yy]
      if (doPrecip) then event_precip = [[[event_precip]],[[rainnc]]]
      if (doIVT)    then event_ivtU   = [[[event_ivtU]],  [[ivtU]]]
      if (doIVT)    then event_ivtV   = [[[event_ivtV]],  [[ivtV]]]
      if (doz0k)    then event_z0k    = [[[event_z0k]],   [[z0k]]]
      if (doLCL)    then event_zLCL   = [[[event_zLCL]],  [[zLCL]]]
      if (doZ500)   then event_z500   = [[[event_z500]],  [[z500]]]
      if (doSWE)    then event_swe    = [[[event_swe]],   [[swe]]]
   endif
   init = 0

   badEvent:
end

if (doz0k) then begin
   bground = where(event_z0k lt 0)
   if(bground(0) ne -1) then event_z0k(bground)=0
endif

if (doIVT) then  event_ivt = sqrt(event_ivtU*event_ivtU + event_ivtV*event_ivtV)

; ####################################################################################
; Write output
; ####################################################################################
print,'Writing to output...'

; Precip
if (doPrecip) then begin
   fileOUT_precip = 'composite.raw.precip.'+fileIN
   fileID_precip  = ncdf_create(dirOUT+fileOUT_precip,/clobber)
   dimID1  = ncdf_dimdef(fileID_precip,'lon',nlon)
   dimID2  = ncdf_dimdef(fileID_precip,'lat',nlat)
   dimID3  = ncdf_dimdef(fileID_precip,'ntime',n_elements(event_day))
   varID1  = ncdf_vardef(fileID_precip,'lon',  [dimID1,dimID2       ],/float)
   varID2  = ncdf_vardef(fileID_precip,'lat',  [dimID1,dimID2       ],/float)
   varID3  = ncdf_vardef(fileID_precip,'year', [              dimID3],/long)
   varID4  = ncdf_vardef(fileID_precip,'month',[              dimID3],/long)
   varID5  = ncdf_vardef(fileID_precip,'day',  [              dimID3],/long)
   varID9  = ncdf_vardef(fileID_precip,'hour', [              dimID3],/long)
   varID6  = ncdf_vardef(fileID_precip,'precip', [dimID1,dimID2,dimID3],/float)
   ncdf_control,fileID_precip,/endef
   ncdf_varput,fileID_precip,varID1,lon
   ncdf_varput,fileID_precip,varID2,lat
   ncdf_varput,fileID_precip,varID3,event_year
   ncdf_varput,fileID_precip,varID4,event_month
   ncdf_varput,fileID_precip,varID5,event_day
   ncdf_varput,fileID_precip,varID9,event_hour
   ncdf_varput,fileID_precip,varID6,event_precip
   ncdf_close,fileID_precip
endif

; SWE
if (doSwe) then begin
   fileOUT_swe = 'composite.raw.swe.'+fileIN
   fileID_swe  = ncdf_create(dirOUT+fileOUT_swe,/clobber)
   dimID1      = ncdf_dimdef(fileID_swe,'lon',nlon)
   dimID2      = ncdf_dimdef(fileID_swe,'lat',nlat)
   dimID3      = ncdf_dimdef(fileID_swe,'ntime',n_elements(event_day))
   varID1      = ncdf_vardef(fileID_swe,'lon',  [dimID1,dimID2       ],/float)
   varID2      = ncdf_vardef(fileID_swe,'lat',  [dimID1,dimID2       ],/float)
   varID3      = ncdf_vardef(fileID_swe,'year', [              dimID3],/long)
   varID4      = ncdf_vardef(fileID_swe,'month',[              dimID3],/long)
   varID5      = ncdf_vardef(fileID_swe,'day',  [              dimID3],/long)
   varID9      = ncdf_vardef(fileID_swe,'hour', [              dimID3],/long)
   varID6      = ncdf_vardef(fileID_swe,'swe',  [dimID1,dimID2,dimID3],/float)
   ncdf_control,fileID_swe,/endef
   ncdf_varput,fileID_swe,varID1,lon
   ncdf_varput,fileID_swe,varID2,lat
   ncdf_varput,fileID_swe,varID3,event_year
   ncdf_varput,fileID_swe,varID4,event_month
   ncdf_varput,fileID_swe,varID5,event_day
   ncdf_varput,fileID_swe,varID9,event_hour
   ncdf_varput,fileID_swe,varID6,event_swe
   ncdf_close,fileID_swe
endif

; z500
if (doZ500) then begin
   fileOUT_z500 = 'composite.raw.z500.'+fileIN
   fileID_z500  = ncdf_create(dirOUT+fileOUT_z500,/clobber)
   dimID1       = ncdf_dimdef(fileID_z500,'lon',nlon)
   dimID2       = ncdf_dimdef(fileID_z500,'lat',nlat)
   dimID3       = ncdf_dimdef(fileID_z500,'ntime',n_elements(event_day))
   varID1       = ncdf_vardef(fileID_z500,'lon',  [dimID1,dimID2       ],/float)
   varID2       = ncdf_vardef(fileID_z500,'lat',  [dimID1,dimID2       ],/float)
   varID3       = ncdf_vardef(fileID_z500,'year', [              dimID3],/long)
   varID4       = ncdf_vardef(fileID_z500,'month',[              dimID3],/long)
   varID5       = ncdf_vardef(fileID_z500,'day',  [              dimID3],/long)
   varID9       = ncdf_vardef(fileID_z500,'hour', [              dimID3],/long)
   varID6       = ncdf_vardef(fileID_z500,'z500', [dimID1,dimID2,dimID3],/float)
   ncdf_control,fileID_z500,/endef
   ncdf_varput,fileID_z500,varID1,lon
   ncdf_varput,fileID_z500,varID2,lat
   ncdf_varput,fileID_z500,varID3,event_year
   ncdf_varput,fileID_z500,varID4,event_month
   ncdf_varput,fileID_z500,varID5,event_day
   ncdf_varput,fileID_z500,varID9,event_hour
   ncdf_varput,fileID_z500,varID6,event_z500
   ncdf_close,fileID_z500
endif

; IVT
if (doIVT) then begin
   fileOUT_ivt = 'composite.raw.ivt.'+fileIN
   fileID_ivt  = ncdf_create(dirOUT+fileOUT_ivt,/clobber)
   dimID1      = ncdf_dimdef(fileID_ivt,'lon',nlon)
   dimID2      = ncdf_dimdef(fileID_ivt,'lat',nlat)
   dimID3      = ncdf_dimdef(fileID_ivt,'ntime',n_elements(event_day))
   varID1      = ncdf_vardef(fileID_ivt,'lon',  [dimID1,dimID2       ],/float)
   varID2      = ncdf_vardef(fileID_ivt,'lat',  [dimID1,dimID2       ],/float)
   varID3      = ncdf_vardef(fileID_ivt,'year', [              dimID3],/long)
   varID4      = ncdf_vardef(fileID_ivt,'month',[              dimID3],/long)
   varID5      = ncdf_vardef(fileID_ivt,'day',  [              dimID3],/long)
   varID9      = ncdf_vardef(fileID_ivt,'hour', [              dimID3],/long)
   varID7      = ncdf_vardef(fileID_ivt,'ivt',  [dimID1,dimID2,dimID3],/float)
   ncdf_control,fileID_ivt,/endef
   ncdf_varput,fileID_ivt,varID1,lon
   ncdf_varput,fileID_ivt,varID2,lat
   ncdf_varput,fileID_ivt,varID3,event_year
   ncdf_varput,fileID_ivt,varID4,event_month
   ncdf_varput,fileID_ivt,varID5,event_day
   ncdf_varput,fileID_ivt,varID9,event_hour
   ncdf_varput,fileID_ivt,varID7,event_ivt
   ncdf_close,fileID_ivt
endif

; Z0k
if (doz0k) then begin
   fileOUT_z0k = 'composite.raw.z0k.'+fileIN
   fileID_z0k  = ncdf_create(dirOUT+fileOUT_z0k,/clobber)
   dimID1      = ncdf_dimdef(fileID_z0k,'lon',nlon)
   dimID2      = ncdf_dimdef(fileID_z0k,'lat',nlat)
   dimID3      = ncdf_dimdef(fileID_z0k,'ntime',n_elements(event_day))
   varID1      = ncdf_vardef(fileID_z0k,'lon',  [dimID1,dimID2       ],/float)
   varID2      = ncdf_vardef(fileID_z0k,'lat',  [dimID1,dimID2       ],/float)
   varID3      = ncdf_vardef(fileID_z0k,'year', [              dimID3],/long)
   varID4      = ncdf_vardef(fileID_z0k,'month',[              dimID3],/long)
   varID5      = ncdf_vardef(fileID_z0k,'day',  [              dimID3],/long)
   varID9      = ncdf_vardef(fileID_z0k,'hour', [              dimID3],/long)
   varID10     = ncdf_vardef(fileID_z0k,'z0k',  [dimID1,dimID2,dimID3],/float)
   ncdf_control,fileID_z0k,/endef
   ncdf_varput,fileID_z0k,varID1,lon
   ncdf_varput,fileID_z0k,varID2,lat
   ncdf_varput,fileID_z0k,varID3,event_year
   ncdf_varput,fileID_z0k,varID4,event_month
   ncdf_varput,fileID_z0k,varID5,event_day
   ncdf_varput,fileID_z0k,varID9,event_hour
   ncdf_varput,fileID_z0k,varID10,event_z0k
   ncdf_close,fileID_z0k
endif

; zLCL
if (doLCL) then begin
   fileOUT_zLCL = 'composite.raw.zLCL.'+fileIN
   fileID_zLCL  = ncdf_create(dirOUT+fileOUT_zLCL,/clobber)
   dimID1       = ncdf_dimdef(fileID_zLCL,'lon',nlon)
   dimID2       = ncdf_dimdef(fileID_zLCL,'lat',nlat)
   dimID3       = ncdf_dimdef(fileID_zLCL,'ntime',n_elements(event_day))
   varID1       = ncdf_vardef(fileID_zLCL,'lon',  [dimID1,dimID2       ],/float)
   varID2       = ncdf_vardef(fileID_zLCL,'lat',  [dimID1,dimID2       ],/float)
   varID3       = ncdf_vardef(fileID_zLCL,'year', [              dimID3],/long)
   varID4       = ncdf_vardef(fileID_zLCL,'month',[              dimID3],/long)
   varID5       = ncdf_vardef(fileID_zLCL,'day',  [              dimID3],/long)
   varID9       = ncdf_vardef(fileID_zLCL,'hour', [              dimID3],/long)
   varID10      = ncdf_vardef(fileID_zLCL,'zLCL', [dimID1,dimID2,dimID3],/float)
   ncdf_control,fileID_zLCL,/endef
   ncdf_varput,fileID_zLCL,varID1,lon
   ncdf_varput,fileID_zLCL,varID2,lat
   ncdf_varput,fileID_zLCL,varID3,event_year
   ncdf_varput,fileID_zLCL,varID4,event_month
   ncdf_varput,fileID_zLCL,varID5,event_day
   ncdf_varput,fileID_zLCL,varID9,event_hour
   ncdf_varput,fileID_zLCL,varID10,event_zLCL
   ncdf_close,fileID_zLCL
endif

; END PROGRAM
end

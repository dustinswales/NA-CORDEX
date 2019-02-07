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
pro compositeEvents,dirIN,dirOUT,modelID,ivtThresh,persistenceThresh,future,fileOUT_precip,fileOUT_ivt,fileOUT_z0k
; Output directory
; dirOUT = '/data/dswales/NA-CORDEX/ARdet/'
; IVT threshold (standardized)
;ivtThresh = 99.
; How long are AR events? (in hours)
;persistenceThresh = 24
; Which RCM
;modelID='gfdl'
  
; 50 or 25km data?
res = '50'

; WRF data
expID = 'WRF_'+modelID+'_'+res+'km'
dir   = '/Projects/HydroMet/dswales/NA-CORDEX/'+expID+'/'

; For composites, how many hours after event to include?
lingerTime = 12

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
      ncdf_varget,fileID,ncdf_varid(fileID,'IVTU'),ivtU,offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
      ncdf_varget,fileID,ncdf_varid(fileID,'IVTV'),ivtV,offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
      ncdf_varget,fileID,ncdf_varid(fileID,'Z0K'),z0k,offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
      ncdf_varget,fileID,ncdf_varid(fileID,'RAINNC'),var_temp2,offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]
      nt   = n_elements(ivtU(0,0,*))
      rainnc = fltarr(nlon,nlat,nt)
      rainnc(*,*,1:nt-1)=var_temp2(*,*,1:nt-1)-var_temp2(*,*,0:nt-2)
      rainnc(where(rainnc lt 0)) = 0
      
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
      ncdf_varget,fileID,ncdf_varid(fileID,'IVTU'),ivtUa,offset=[0,0,t0],count=[nlon,nlat,nt0]
      ncdf_varget,fileID,ncdf_varid(fileID,'IVTV'),ivtVa,offset=[0,0,t0],count=[nlon,nlat,nt0]
      ncdf_varget,fileID,ncdf_varid(fileID,'Z0K'),z0ka,offset=[0,0,t0],count=[nlon,nlat,nt0]
      ncdf_varget,fileID,ncdf_varid(fileID,'RAINNC'),rainnca,offset=[0,0,t0],count=[nlon,nlat,nt0]
      nt   = n_elements(ivtUa(0,0,*))
      var_temp2a = fltarr(nlon,nlat,nt)
      if (nt gt 1) then begin
         var_temp2a(*,*,1:nt-1)=rainnca(*,*,1:nt-1)-rainnca(*,*,0:nt-2)
         var_temp2a(where(var_temp2a lt 0)) = 0
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
      ncdf_varget,fileID,ncdf_varid(fileID,'IVTU'),ivtUb,count=[nlon,nlat,nt1]
      ncdf_varget,fileID,ncdf_varid(fileID,'IVTV'),ivtVb,count=[nlon,nlat,nt1]
      ncdf_varget,fileID,ncdf_varid(fileID,'Z0K'),z0kb,count=[nlon,nlat,nt1]
      ncdf_varget,fileID,ncdf_varid(fileID,'RAINNC'),rainncb,count=[nlon,nlat,nt1]
      nt   = n_elements(ivtUb(0,0,*))
      var_temp2b = fltarr(nlon,nlat,nt)
      var_temp2b(*,*,1:nt-1)=rainncB(*,*,1:nt-1)-rainncB(*,*,0:nt-2)
      var_temp2b(where(var_temp2b lt 0)) = 0
      yyb = year(0:nt1-1)
      mmb = month(0:nt1-1)
      ddb = day(0:nt1-1)
      hhb = hour(0:nt1-1)

      ; Close file 2
      ncdf_close,fileID

      ; Combine
      ivtU = [[[ivtUa]],[[ivtUb]]]
      ivtV = [[[ivtVa]],[[ivtVb]]]
      z0k = [[[z0ka]],[[z0kb]]]
      rainnc = [[[var_temp2A]],[[var_temp2B]]]
      yy   = [yya,yyb]
      mm   = [mma,mmb]
      dd   = [dda,ddb]
      hh   = [hha,hhb]
   endif
   nt = n_elements(ivtU(0,0,*))

   ; Store data from event(s).
   if (init) then begin
      event_hour   = hh 
      event_day    = dd
      event_month  = mm
      event_year   = yy
      event_precip = rainnc
      event_ivtU   = ivtU
      event_ivtV   = ivtV
      event_z0k    = z0k
   endif
   if (not init) then begin
      event_hour   = [event_hour, hh]
      event_day    = [event_day,  dd]
      event_month  = [event_month,mm]
      event_year   = [event_year, yy]
      event_precip = [[[event_precip]],[[rainnc]]]
      event_ivtU   = [[[event_ivtU]],[[ivtU]]]
      event_ivtV   = [[[event_ivtV]],[[ivtV]]]
      event_z0k    = [[[event_z0k]],[[z0k]]]
   endif
   init = 0

   badEvent:
end

bground = where(event_z0k lt 0)
if(bground(0) ne -1) then event_z0k(bground)=0

event_ivt = sqrt(event_ivtU*event_ivtU + event_ivtV*event_ivtV)

; Write output
print,'Writing to output...'
fileOUT_precip = 'composite.raw.precip.'+fileIN
fileOUT_ivt    = 'composite.raw.ivt.'+fileIN
fileOUT_z0k    = 'composite.raw.z0k.'+fileIN
; Precip file
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
; IVT file
fileID_ivt  = ncdf_create(dirOUT+fileOUT_ivt,/clobber)
dimID1  = ncdf_dimdef(fileID_ivt,'lon',nlon)
dimID2  = ncdf_dimdef(fileID_ivt,'lat',nlat)
dimID3  = ncdf_dimdef(fileID_ivt,'ntime',n_elements(event_day))
varID1  = ncdf_vardef(fileID_ivt,'lon',  [dimID1,dimID2       ],/float)
varID2  = ncdf_vardef(fileID_ivt,'lat',  [dimID1,dimID2       ],/float)
varID3  = ncdf_vardef(fileID_ivt,'year', [              dimID3],/long)
varID4  = ncdf_vardef(fileID_ivt,'month',[              dimID3],/long)
varID5  = ncdf_vardef(fileID_ivt,'day',  [              dimID3],/long)
varID9  = ncdf_vardef(fileID_ivt,'hour', [              dimID3],/long)
varID7  = ncdf_vardef(fileID_ivt,'ivt',  [dimID1,dimID2,dimID3],/float)
ncdf_control,fileID_ivt,/endef
ncdf_varput,fileID_ivt,varID1,lon
ncdf_varput,fileID_ivt,varID2,lat
ncdf_varput,fileID_ivt,varID3,event_year
ncdf_varput,fileID_ivt,varID4,event_month
ncdf_varput,fileID_ivt,varID5,event_day
ncdf_varput,fileID_ivt,varID9,event_hour
ncdf_varput,fileID_ivt,varID7,event_ivt
ncdf_close,fileID_ivt
; Z0k file
fileID_z0k  = ncdf_create(dirOUT+fileOUT_z0k,/clobber)
dimID1  = ncdf_dimdef(fileID_z0k,'lon',nlon)
dimID2  = ncdf_dimdef(fileID_z0k,'lat',nlat)
dimID3  = ncdf_dimdef(fileID_z0k,'ntime',n_elements(event_day))
varID1  = ncdf_vardef(fileID_z0k,'lon',  [dimID1,dimID2       ],/float)
varID2  = ncdf_vardef(fileID_z0k,'lat',  [dimID1,dimID2       ],/float)
varID3  = ncdf_vardef(fileID_z0k,'year', [              dimID3],/long)
varID4  = ncdf_vardef(fileID_z0k,'month',[              dimID3],/long)
varID5  = ncdf_vardef(fileID_z0k,'day',  [              dimID3],/long)
varID9  = ncdf_vardef(fileID_z0k,'hour', [              dimID3],/long)
varID10 = ncdf_vardef(fileID_z0k,'z0k',  [dimID1,dimID2,dimID3],/float)
ncdf_control,fileID_z0k,/endef
ncdf_varput,fileID_z0k,varID1,lon
ncdf_varput,fileID_z0k,varID2,lat
ncdf_varput,fileID_z0k,varID3,event_year
ncdf_varput,fileID_z0k,varID4,event_month
ncdf_varput,fileID_z0k,varID5,event_day
ncdf_varput,fileID_z0k,varID9,event_hour
ncdf_varput,fileID_z0k,varID10,event_z0k
ncdf_close,fileID_z0k

; END PROGRAM
end

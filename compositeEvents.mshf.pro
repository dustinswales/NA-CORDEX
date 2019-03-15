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
   f0 = 'wrfout_'+string(yearS(ij),format='(i4)')+string(monthS(ij),format='(i2.2)')+'.mshf.nc'
   f1 = 'wrfout_'+string(yearF(ij),format='(i4)')+string(monthF(ij),format='(i2.2)')+'.mshf.nc'
   if (strcmp(f0,f1)) then begin
      ;print,"data in single file"
      fileID = ncdf_open(dir+f0)
      
      ; Determine what time slice to read in?
      ncdf_varget,fileID,ncdf_varid(fileID,'day'),day
      ncdf_varget,fileID,ncdf_varid(fileID,'hour'),hour
      ncdf_varget,fileID,ncdf_varid(fileID,'month'),month
      ncdf_varget,fileID,ncdf_varid(fileID,'year'),year
      if (ij eq 0) then begin
         ncdf_varget,fileID,ncdf_varid(fileID,'lon'),lon
         ncdf_varget,fileID,ncdf_varid(fileID,'lat'),lat
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
      ncdf_varget,fileID,ncdf_varid(fileID,'mshf'),mshf,offset=[0,0,t0],count=[nlon,nlat,event_length(ij)/3.]

      
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
      ncdf_varget,fileID,ncdf_varid(fileID,'day'),day
      ncdf_varget,fileID,ncdf_varid(fileID,'hour'),hour
      ncdf_varget,fileID,ncdf_varid(fileID,'month'),month
      ncdf_varget,fileID,ncdf_varid(fileID,'year'),year
      ncdf_varget,fileID,ncdf_varid(fileID,'lon'),lon
      ncdf_varget,fileID,ncdf_varid(fileID,'lat'),lat
      nlon = n_elements(lon(*,0))
      nlat = n_elements(lat(0,*))
      t0 = where(Year eq YearS(ij) and Month eq MonthS(ij) and Day eq DayS(ij) and Hour eq HourS(ij))
      if (t0(0) eq -1) then begin
         print,'ERROR: something went wrong?'
         stop
      endif

      ; Read in fields
      nt0 = n_elements(day)-t0
      ncdf_varget,fileID,ncdf_varid(fileID,'mshf'),mshfa,offset=[0,0,t0],count=[nlon,nlat,nt0]
      nt = n_elements(mshfa(0,0,*))
      
      ; Store data
      yya = year(t0:t0+nt0-1)
      mma = month(t0:t0+nt0-1)
      dda = day(t0:t0+nt0-1)
      hha = hour(t0:t0+nt0-1)
      
      ; Close file 1
      ncdf_close,fileID

      ; File 2
      fileID = ncdf_open(dir+f1)
      ncdf_varget,fileID,ncdf_varid(fileID,'day'),day
      ncdf_varget,fileID,ncdf_varid(fileID,'hour'),hour
      ncdf_varget,fileID,ncdf_varid(fileID,'month'),month
      ncdf_varget,fileID,ncdf_varid(fileID,'year'),year
      nt1 = event_length(ij)-nt0
      ncdf_varget,fileID,ncdf_varid(fileID,'mshf'),mshfb,count=[nlon,nlat,nt1]
      yyb = year(0:nt1-1)
      mmb = month(0:nt1-1)
      ddb = day(0:nt1-1)
      hhb = hour(0:nt1-1)
      
      ; Close file 2
      ncdf_close,fileID

      ; Combine
      mshf = [[[mshfa]],[[mshfb]]]
      yy   = [yya,yyb]
      mm   = [mma,mmb]
      dd   = [dda,ddb]
      hh   = [hha,hhb]
   endif
   nt = n_elements(mshf(0,0,*))

   ; Store data from event(s).
   if (init) then begin
      event_hour   = hh 
      event_day    = dd
      event_month  = mm
      event_year   = yy
      event_mshf   = mshf
   endif
   if (not init) then begin
      event_hour  = [event_hour, hh]
      event_day   = [event_day,  dd]
      event_month = [event_month,mm]
      event_year  = [event_year, yy]
      event_mshf  = [[[event_mshf]],[[mshf]]]
   endif
   init = 0
   badEvent:
end

; Write output
print,'Writing to output...'
fileOUT_mshf = 'composite.raw.mshf.'+fileIN

; Mshf file
fileID_mshf  = ncdf_create(dirOUT+fileOUT_mshf,/clobber)
dimID1  = ncdf_dimdef(fileID_mshf,'lon',nlon)
dimID2  = ncdf_dimdef(fileID_mshf,'lat',nlat)
dimID3  = ncdf_dimdef(fileID_mshf,'ntime',n_elements(event_day))
varID1  = ncdf_vardef(fileID_mshf,'lon',  [dimID1,dimID2       ],/float)
varID2  = ncdf_vardef(fileID_mshf,'lat',  [dimID1,dimID2       ],/float)
varID3  = ncdf_vardef(fileID_mshf,'year', [              dimID3],/long)
varID4  = ncdf_vardef(fileID_mshf,'month',[              dimID3],/long)
varID5  = ncdf_vardef(fileID_mshf,'day',  [              dimID3],/long)
varID9  = ncdf_vardef(fileID_mshf,'hour', [              dimID3],/long)
varID6  = ncdf_vardef(fileID_mshf,'mshf', [dimID1,dimID2,dimID3],/float)
ncdf_control,fileID_mshf,/endef
ncdf_varput,fileID_mshf,varID1,lon
ncdf_varput,fileID_mshf,varID2,lat
ncdf_varput,fileID_mshf,varID3,event_year
ncdf_varput,fileID_mshf,varID4,event_month
ncdf_varput,fileID_mshf,varID5,event_day
ncdf_varput,fileID_mshf,varID9,event_hour
ncdf_varput,fileID_mshf,varID6,event_mshf
ncdf_close,fileID_mshf

; END PROGRAM
end

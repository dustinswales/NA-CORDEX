
doIVT  = 1
doSWE  = 1
doz500 = 1
dozLCL = 1

; WRF RCM
modelID = ['gfdl','hadgem','mpi','erain']
nmod    = n_elements(modelID)

; Persistence threshold
phold = [24,36]
np1   = n_elements(phold)
ptile = [0.95,0.97];,0.98,0.99]
np2   = n_elements(ptile)

; Directory w/ climotology
dirClimo = '/Projects/HydroMet/dswales/NA-CORDEX/analysis/climo/'
; Output directory
dirOUT="/data/dswales/NA-CORDEX/ARdet/composites/"

; ############################################################
; Compute precipitaion climo @ each gridpoint
; ############################################################
for imod=0,nmod-1 do begin
   
   ; ############################################################
   ; Read in all precipitation data.
   ; ############################################################
   filePrecipH = dirClimo+'climo.rainnc.WRF_'        + modelID(imod) + '_50km.3hr.nc'
   filePrecipF = dirClimo+'climo.future.rainnc.WRF_' + modelID(imod) + '_50km.3hr.nc'
   if (strcmp(modelID(imod),'erain')) then filePrecipF=filePrecipH

   ; Read in raw precipitation data (historical). This file contains all
   ; precipitation data at each point, from ONDJFM.
   fileID = ncdf_open(filePrecipH)
   ncdf_varget,fileID,ncdf_varid(fileID,'lon'),lon
   ncdf_varget,fileID,ncdf_varid(fileID,'lat'),lat
   ncdf_varget,fileID,ncdf_varid(fileID,'rainnc'),rainnc
   ncdf_close,fileID
   nlon=n_elements(lon(*,0))
   nlat=n_elements(lon(0,*))
   
   ; Read in raw precipitation data (future). This file contains all
   ; precipitation data at each point, from ONDJFM.   
   fileID = ncdf_open(filePrecipF)
   ncdf_varget,fileID,ncdf_varid(fileID,'rainnc'),rainncf
   ncdf_close,fileID
   ;
   test1 = where(finite(rainnc) eq 0)
   if (test1(0) ne -1) then rainnc(test1)=0
   test2 = where(finite(rainncf) eq 0)
   if (test2(0) ne -1) then rainncf(test2)=0

   ; Read in other fields (if requested)
   if (doIVT) then begin
      fileIVTH = dirClimo+'climo.ivt.WRF_'        + modelID(imod) + '_50km.3hr.nc'
      fileIVTF = dirClimo+'climo.future.ivt.WRF_' + modelID(imod) + '_50km.3hr.nc'
      if (strcmp(modelID(imod),'erain')) then fileIVTF=fileIVTH
      fileID = ncdf_open(fileIVTH)
      ncdf_varget,fileID,ncdf_varid(fileID,'ivt'),ivt
      ncdf_close,fileID
      fileID = ncdf_open(fileIVTF)
      ncdf_varget,fileID,ncdf_varid(fileID,'ivt'),ivtf
      ncdf_close,fileID
      ;
      test1 = where(finite(ivt) eq 0)
      if (test1(0) ne -1) then ivt(test1)=0
      test2 = where(finite(ivtf) eq 0)
      if (test2(0) ne -1) then ivtf(test2)=0
      test3 = where(ivt gt 2000)
      if (test3(0) ne -1) then ivt(test3)=0
      test4 = where(ivtf gt 2000)
      if (test4(0) ne -1) then ivtf(test4)=0
   endif
   ;
   if (doZLCL) then begin
      filezLCLH = dirClimo+'climo.zLCL.WRF_'        + modelID(imod) + '_50km.3hr.nc'
      filezLCLF = dirClimo+'climo.future.zLCL.WRF_' + modelID(imod) + '_50km.3hr.nc'
      if (strcmp(modelID(imod),'erain')) then filezLCLF=filezLCLH
      fileID = ncdf_open(filezLCLH)
      ncdf_varget,fileID,ncdf_varid(fileID,'zLCL'),zLCL
      ncdf_close,fileID
      fileID = ncdf_open(filezLCLF)
      ncdf_varget,fileID,ncdf_varid(fileID,'zLCL'),zLCLf
      ncdf_close,fileID
      ;
      test1 = where(finite(zLCL) eq 0)
      if (test1(0) ne -1) then zLCL(test1)=0
      test2 = where(finite(zLCLf) eq 0)
      if (test2(0) ne -1) then zLCLf(test2)=0
   endif
   ;
   if (doZ500) then begin
      filez500H = dirClimo+'climo.z500.WRF_'        + modelID(imod) + '_50km.3hr.nc'
      filez500F = dirClimo+'climo.future.z500.WRF_' + modelID(imod) + '_50km.3hr.nc'
      if (strcmp(modelID(imod),'erain')) then filez500F=filez500H
      fileID = ncdf_open(filez500H)
      ncdf_varget,fileID,ncdf_varid(fileID,'z500'),z500
      ncdf_close,fileID
      fileID = ncdf_open(filez500F)
      ncdf_varget,fileID,ncdf_varid(fileID,'z500'),z500f
      ncdf_close,fileID
      test1 = where(finite(z500) eq 0)
      if (test1(0) ne -1) then z500(test1)=0
      test2 = where(finite(z500f) eq 0)
      if (test2(0) ne -1) then z500f(test2)=0
   endif
   ;
   if (doSWE) then begin
      filesweH = dirClimo+'climo.swe.WRF_'        + modelID(imod) + '_50km.3hr.nc'
      filesweF = dirClimo+'climo.future.swe.WRF_' + modelID(imod) + '_50km.3hr.nc'
      if (strcmp(modelID(imod),'erain')) then filesweF=filesweH
      fileID = ncdf_open(filesweH)
      ncdf_varget,fileID,ncdf_varid(fileID,'swe'),swe
      ncdf_close,fileID
      fileID = ncdf_open(filesweF)
      ncdf_varget,fileID,ncdf_varid(fileID,'swe'),swef
      ncdf_close,fileID
      test1 = where(finite(swe) eq 0)
      if (test1(0) ne -1) then swe(test1)=0
      test2 = where(finite(swef) eq 0)
      if (test2(0) ne -1) then swef(test2)=0
   endif
   
   ; ############################################################
   ; Read in event composites (precipitation)
   ; ############################################################
   for ip1=0,np1-1 do begin
      for ip2=0,np2-1 do begin
         print,modelID(imod),phold(ip1),ptile(ip2)*100.
         ; Get filenames for requested case
         ; Precipitation (default)
         fileComposite_precip_hist = file_search(dirOUT,'composite.raw.precip.events.50km.'         +$
                                                 string(100*ptile(ip2),format='(f5.2)')+'ptile.'    +$
                                                 string(phold(ip1),format='(i2.2)')+'hrs.'          +$
                                                 modelID(imod)+'.nc')
         fileComposite_precip_future = file_search(dirOUT,'composite.raw.precip.events.future.50km.'+$
                                                   string(100*ptile(ip2),format='(f5.2)')+'ptile.'  +$
                                                   string(phold(ip1),format='(i2.2)')+'hrs.'        +$
                                                   modelID(imod)+'.nc')
         if (strcmp(modelID(imod),'erain')) then begin
            fileComposite_precip_future = fileComposite_precip_hist
         endif
         ; historical
         fileID_precip = ncdf_open(fileComposite_precip_hist)
         ncdf_varget,fileID_precip,ncdf_varid(fileID_precip,'year'),   event_year_hist
         ncdf_varget,fileID_precip,ncdf_varid(fileID_precip,'month'),  event_month_hist
         ncdf_varget,fileID_precip,ncdf_varid(fileID_precip,'day'),    event_day_hist
         ncdf_varget,fileID_precip,ncdf_varid(fileID_precip,'precip'), event_precip_hist
         ncdf_close,fileID_precip
         ntime_event_hist = n_elements(event_day_hist)

         ; future
         fileID_precip = ncdf_open(fileComposite_precip_future)
         ncdf_varget,fileID_precip,ncdf_varid(fileID_precip,'year'),   event_year_future
         ncdf_varget,fileID_precip,ncdf_varid(fileID_precip,'month'),  event_month_future
         ncdf_varget,fileID_precip,ncdf_varid(fileID_precip,'day'),    event_day_future
         ncdf_varget,fileID_precip,ncdf_varid(fileID_precip,'precip'), event_precip_future
         ncdf_close,fileID_precip
         ntime_event_future = n_elements(event_day_future)
         
         ; IVT (optional)
         if (doIVT) then begin
            fileComposite_ivt_hist = file_search(dirOUT,'composite.raw.ivt.events.50km.'            +$
                                                 string(100*ptile(ip2),format='(f5.2)')+'ptile.'    +$
                                                 string(phold(ip1),format='(i2.2)')+'hrs.'          +$
                                                 modelID(imod)+'.nc')
            fileComposite_ivt_future = file_search(dirOUT,'composite.raw.ivt.events.future.50km.'   +$
                                                   string(100*ptile(ip2),format='(f5.2)')+'ptile.'  +$
                                                   string(phold(ip1),format='(i2.2)')+'hrs.'        +$
                                                   modelID(imod)+'.nc')
            if (strcmp(modelID(imod),'erain')) then begin
               fileComposite_ivt_future = fileComposite_ivt_hist
            endif
            ; historical
            fileID_ivt = ncdf_open(fileComposite_ivt_hist)
            ncdf_varget,fileID_ivt,ncdf_varid(fileID_ivt,'ivt'), event_ivt_hist
            ncdf_close,fileID_ivt

            ; future
            fileID_ivt = ncdf_open(fileComposite_ivt_future)
            ncdf_varget,fileID_ivt,ncdf_varid(fileID_ivt,'ivt'), event_ivt_future
            ncdf_close,fileID_ivt
         endif
         
         ; SWE (optional)
         if (doSWE) then begin
            fileComposite_swe_hist = file_search(dirOUT,'composite.raw.swe.events.50km.'            +$
                                                 string(100*ptile(ip2),format='(f5.2)')+'ptile.'    +$
                                                 string(phold(ip1),format='(i2.2)')+'hrs.'          +$
                                                 modelID(imod)+'.nc')
            fileComposite_swe_future = file_search(dirOUT,'composite.raw.swe.events.future.50km.'   +$
                                                   string(100*ptile(ip2),format='(f5.2)')+'ptile.'  +$
                                                   string(phold(ip1),format='(i2.2)')+'hrs.'        +$
                                                   modelID(imod)+'.nc')
            if (strcmp(modelID(imod),'erain')) then begin
               fileComposite_swe_future = fileComposite_swe_hist
            endif
            ; historical
            fileID_swe = ncdf_open(fileComposite_swe_hist)
            ncdf_varget,fileID_swe,ncdf_varid(fileID_swe,'swe'), event_swe_hist
            ncdf_close,fileID_swe

            ; future
            fileID_swe = ncdf_open(fileComposite_swe_future)
            ncdf_varget,fileID_swe,ncdf_varid(fileID_swe,'swe'), event_swe_future
            ncdf_close,fileID_swe            
         endif
         
         ; zLCL (optional)
         if (dozLCL) then begin
            fileComposite_zLCL_hist = file_search(dirOUT,'composite.raw.zLCL.events.50km.'          +$
                                                 string(100*ptile(ip2),format='(f5.2)')+'ptile.'    +$
                                                 string(phold(ip1),format='(i2.2)')+'hrs.'          +$
                                                 modelID(imod)+'.nc')
            fileComposite_zLCL_future = file_search(dirOUT,'composite.raw.zLCL.events.future.50km.' +$
                                                   string(100*ptile(ip2),format='(f5.2)')+'ptile.'  +$
                                                   string(phold(ip1),format='(i2.2)')+'hrs.'        +$
                                                   modelID(imod)+'.nc')
            if (strcmp(modelID(imod),'erain')) then begin
               fileComposite_zLCL_future = fileComposite_zLCL_hist
            endif
            ; historical
            fileID_zLCL = ncdf_open(fileComposite_zLCL_hist)
            ncdf_varget,fileID_zLCL,ncdf_varid(fileID_zLCL,'zLCL'), event_zLCL_hist
            ncdf_close,fileID_zLCL

            ; future
            fileID_zLCL = ncdf_open(fileComposite_zLCL_future)
            ncdf_varget,fileID_zLCL,ncdf_varid(fileID_zLCL,'zLCL'), event_zLCL_future
            ncdf_close,fileID_zLCL            
         endif
         
         ; z500 (optional)
         if (doz500) then begin
            fileComposite_z500_hist = file_search(dirOUT,'composite.raw.z500.events.50km.'          +$
                                                 string(100*ptile(ip2),format='(f5.2)')+'ptile.'    +$
                                                 string(phold(ip1),format='(i2.2)')+'hrs.'          +$
                                                 modelID(imod)+'.nc')
            fileComposite_z500_future = file_search(dirOUT,'composite.raw.z500.events.future.50km.' +$
                                                   string(100*ptile(ip2),format='(f5.2)')+'ptile.'  +$
                                                   string(phold(ip1),format='(i2.2)')+'hrs.'        +$
                                                   modelID(imod)+'.nc')
            if (strcmp(modelID(imod),'erain')) then begin
               fileComposite_z500_future = fileComposite_z500_hist
            endif
            ; historical
            fileID_z500 = ncdf_open(fileComposite_z500_hist)
            ncdf_varget,fileID_z500,ncdf_varid(fileID_z500,'z500'), event_z500_hist
            ncdf_close,fileID_z500

            ; future
            fileID_z500 = ncdf_open(fileComposite_z500_future)
            ncdf_varget,fileID_z500,ncdf_varid(fileID_z500,'z500'), event_z500_future
            ncdf_close,fileID_z500
         endif
         
         ; Composites over entire time period(s)
         ntime_hist   = n_elements(rainnc(0,0,*))
         ntime_future = n_elements(rainncf(0,0,*))
         precip_hist   = total(rainnc,3)  / ntime_hist
         precip_future = total(rainncf,3) / ntime_future
         ivt_hist      = total(ivt,3)     / ntime_hist
         ivt_future    = total(ivtf,3)    / ntime_future
         swe_hist      = total(swe,3)     / ntime_hist
         swe_future    = total(swef,3)    / ntime_future
         z500_hist     = total(z500,3)    / ntime_hist
         z500_future   = total(z500f,3)   / ntime_future
         zLCL_hist     = total(zLCL,3)    / ntime_hist
         zLCL_future   = total(zLCLf,3)   / ntime_future

         ; Composites for non-event times
         nevent_precip_hist   = (total(rainnc,3)  - total(event_precip_hist,3))/  (ntime_hist   - ntime_event_hist   )
         nevent_precip_future = (total(rainncf,3) - total(event_precip_future,3))/(ntime_future - ntime_event_future )
         nevent_ivt_hist      = (total(ivt,3)     - total(event_ivt_hist,3))/     (ntime_hist   - ntime_event_hist   )
         nevent_ivt_future    = (total(ivtf,3)    - total(event_ivt_future,3))/   (ntime_future - ntime_event_future )
         nevent_swe_hist      = (total(swe,3)     - total(event_swe_hist,3))/     (ntime_hist   - ntime_event_hist   )
         nevent_swe_future    = (total(swef,3)    - total(event_swe_future,3))/   (ntime_future - ntime_event_future )
         nevent_z500_hist     = (total(z500,3)    - total(event_z500_hist,3))/    (ntime_hist   - ntime_event_hist   )
         nevent_z500_future   = (total(z500f,3)   - total(event_z500_future,3))/  (ntime_future - ntime_event_future )
         nevent_zLCL_hist     = (total(zLCL,3)    - total(event_zLCL_hist,3))/    (ntime_hist   - ntime_event_hist   )
         nevent_zLCL_future   = (total(zLCLf,3)   - total(event_zLCL_future,3))/  (ntime_future - ntime_event_future )
         
         ; Write output
         print,'Writing to output...'
         fileOUT = dirOUT+'composite.frequency.'+string(phold(ip1),format='(i2.2)')+'hrs.'+$
                   string(100.*ptile(ip2),format='(f5.2)')+'ptile.'+modelID(imod)+'.nc'
         fileID  = ncdf_create(fileOUT,/clobber)
         dimID1  = ncdf_dimdef(fileID,'lon',nlon)
         dimID2  = ncdf_dimdef(fileID,'lat',nlat)
         varID1  = ncdf_vardef(fileID,'lon',                 [dimID1,dimID2],/float)
         varID2  = ncdf_vardef(fileID,'lat',                 [dimID1,dimID2],/float)
         varID3 = ncdf_vardef(fileID,'ntime_event_hist',    1,               /long)
         varID4 = ncdf_vardef(fileID,'ntime_event_future',  1,               /long)
         varID5 = ncdf_vardef(fileID,'ntime_hist',          1,               /long)
         varID6 = ncdf_vardef(fileID,'ntime_future',        1,               /long)
         ;
         varID7  = ncdf_vardef(fileID,'total_precip_hist',   [dimID1,dimID2],/float)
         varID8  = ncdf_vardef(fileID,'total_precip_future', [dimID1,dimID2],/float)
         varID9  = ncdf_vardef(fileID,'total_ivt_hist',      [dimID1,dimID2],/float)
         varID10 = ncdf_vardef(fileID,'total_ivt_future',    [dimID1,dimID2],/float)
         varID11 = ncdf_vardef(fileID,'total_z500_hist',     [dimID1,dimID2],/float)
         varID12 = ncdf_vardef(fileID,'total_z500_future',   [dimID1,dimID2],/float)
         varID13 = ncdf_vardef(fileID,'total_zLCL_hist',     [dimID1,dimID2],/float)
         varID14 = ncdf_vardef(fileID,'total_zLCL_future',   [dimID1,dimID2],/float)
         varID15 = ncdf_vardef(fileID,'total_swe_hist',      [dimID1,dimID2],/float)
         varID16 = ncdf_vardef(fileID,'total_swe_future',    [dimID1,dimID2],/float)
         ;
         varID17 = ncdf_vardef(fileID,'event_precip_hist',   [dimID1,dimID2],/float)
         varID18 = ncdf_vardef(fileID,'event_precip_future', [dimID1,dimID2],/float)         
         varID19 = ncdf_vardef(fileID,'event_ivt_hist',      [dimID1,dimID2],/float)
         varID20 = ncdf_vardef(fileID,'event_ivt_future',    [dimID1,dimID2],/float)
         varID21 = ncdf_vardef(fileID,'event_swe_hist',      [dimID1,dimID2],/float)
         varID22 = ncdf_vardef(fileID,'event_swe_future',    [dimID1,dimID2],/float)
         varID23 = ncdf_vardef(fileID,'event_z500_hist',     [dimID1,dimID2],/float)
         varID24 = ncdf_vardef(fileID,'event_z500_future',   [dimID1,dimID2],/float)
         varID25 = ncdf_vardef(fileID,'event_zLCL_hist',     [dimID1,dimID2],/float)
         varID26 = ncdf_vardef(fileID,'event_zLCL_future',   [dimID1,dimID2],/float)
         ;
         varID29 = ncdf_vardef(fileID,'nevent_precip_hist',  [dimID1,dimID2],/float)
         varID30 = ncdf_vardef(fileID,'nevent_precip_future',[dimID1,dimID2],/float)         
         varID31 = ncdf_vardef(fileID,'nevent_ivt_hist',     [dimID1,dimID2],/float)
         varID32 = ncdf_vardef(fileID,'nevent_ivt_future',   [dimID1,dimID2],/float)
         varID33 = ncdf_vardef(fileID,'nevent_swe_hist',     [dimID1,dimID2],/float)
         varID34 = ncdf_vardef(fileID,'nevent_swe_future',   [dimID1,dimID2],/float)
         varID35 = ncdf_vardef(fileID,'nevent_z500_hist',    [dimID1,dimID2],/float)
         varID36 = ncdf_vardef(fileID,'nevent_z500_future',  [dimID1,dimID2],/float)
         varID37 = ncdf_vardef(fileID,'nevent_zLCL_hist',    [dimID1,dimID2],/float)
         varID38 = ncdf_vardef(fileID,'nevent_zLCL_future',  [dimID1,dimID2],/float)
         ;
         ncdf_control,fileID,/endef
         ncdf_varput,fileID,varID1,lon
         ncdf_varput,fileID,varID2,lat
         ncdf_varput,fileID,varID3,ntime_event_hist
         ncdf_varput,fileID,varID4,ntime_event_future
         ncdf_varput,fileID,varID5,ntime_hist
         ncdf_varput,fileID,varID6,ntime_future
         ;
         ncdf_varput,fileID,varID7,precip_hist
         ncdf_varput,fileID,varID8,precip_future
         ncdf_varput,fileID,varID9,ivt_hist
         ncdf_varput,fileID,varID10,ivt_future
         ncdf_varput,fileID,varID11,z500_hist
         ncdf_varput,fileID,varID12,z500_future
         ncdf_varput,fileID,varID13,zLCL_hist
         ncdf_varput,fileID,varID14,zLCL_future
         ncdf_varput,fileID,varID15,swe_hist
         ncdf_varput,fileID,varID16,swe_future
         ;
         ncdf_varput,fileID,varID17,total(event_precip_hist,3)   / ntime_event_hist
         ncdf_varput,fileID,varID18,total(event_precip_future,3) / ntime_event_future        
         ncdf_varput,fileID,varID19,total(event_ivt_hist,3)      / ntime_event_hist
         ncdf_varput,fileID,varID20,total(event_ivt_future,3)    / ntime_event_future
         ncdf_varput,fileID,varID21,total(event_swe_hist,3)      / ntime_event_hist
         ncdf_varput,fileID,varID22,total(event_swe_future,3)    / ntime_event_future
         ncdf_varput,fileID,varID23,total(event_z500_hist,3)     / ntime_event_hist
         ncdf_varput,fileID,varID24,total(event_z500_future,3)   / ntime_event_future
         ncdf_varput,fileID,varID25,total(event_zLCL_hist,3)     / ntime_event_hist
         ncdf_varput,fileID,varID26,total(event_zLCL_future,3)   / ntime_event_future
         ;
         ncdf_varput,fileID,varID29,nevent_precip_hist
         ncdf_varput,fileID,varID30,nevent_precip_future
         ncdf_varput,fileID,varID31,nevent_ivt_hist
         ncdf_varput,fileID,varID32,nevent_ivt_future
         ncdf_varput,fileID,varID33,nevent_swe_hist
         ncdf_varput,fileID,varID34,nevent_swe_future
         ncdf_varput,fileID,varID35,nevent_z500_hist
         ncdf_varput,fileID,varID36,nevent_z500_future
         ncdf_varput,fileID,varID37,nevent_zLCL_hist
         ncdf_varput,fileID,varID38,nevent_zLCL_future
         ;
         ncdf_close,fileID
      end                       ; Model
   end                          ; Persistence
end                             ; Strength



; END PROGRAM
end

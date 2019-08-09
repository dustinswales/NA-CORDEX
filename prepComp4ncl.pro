
doIVT  = 1
doSWE  = 1
doz500 = 1
dozLCL = 1

; WRF RCM
modelID = ['gfdl','hadgem','erain','mpi']
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
; Climo @ each gridpoint
; ############################################################
for imod=0,nmod-1 do begin
   ; ############################################################
   ; Read in raw precipitation data.
   ; ############################################################
   fileClimoH = dirClimo+'climo.rainnc.WRF_'        + modelID(imod) + '_50km.3hr.nc'
   fileClimoF = dirClimo+'climo.future.rainnc.WRF_' + modelID(imod) + '_50km.3hr.nc'
   if (strcmp(modelID(imod),'erain')) then fileClimoF=fileClimoH

   ; Read in raw precipitation data (historical). This file contains all
   ; precipitation data at each point, from ONDJFM.
   fileID = ncdf_open(fileClimoH)
   ncdf_varget,fileID,ncdf_varid(fileID,'lon'),lon
   ncdf_varget,fileID,ncdf_varid(fileID,'lat'),lat
   ncdf_varget,fileID,ncdf_varid(fileID,'rainnc_precip'),rainnc_precip
   ncdf_close,fileID
   nlon=n_elements(lon(*,0))
   nlat=n_elements(lon(0,*))
   
   ; Read in raw precipitation data (future). This file contains all
   ; precipitation data at each point, from ONDJFM.   
   fileID = ncdf_open(fileClimoF)
   ncdf_varget,fileID,ncdf_varid(fileID,'rainnc_precip'),rainncf_precip
   ncdf_close,fileID
   
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

         ; future
         fileID_precip = ncdf_open(fileComposite_precip_future)
         ncdf_varget,fileID_precip,ncdf_varid(fileID_precip,'year'),   event_year_future
         ncdf_varget,fileID_precip,ncdf_varid(fileID_precip,'month'),  event_month_future
         ncdf_varget,fileID_precip,ncdf_varid(fileID_precip,'day'),    event_day_future
         ncdf_varget,fileID_precip,ncdf_varid(fileID_precip,'precip'), event_precip_future
         ncdf_close,fileID_precip
         ntime = n_elements(event_day_hist)
         
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
         
         ; ###########################################################################
         ; At each point, how do these precipitation events measure-up statistically?   
         ; ###########################################################################
         nPrecipDays_hist    = fltarr(nlon,nlat)
         nEvent_all_hist     = fltarr(nlon,nlat)
         nEvent_event_hist   = fltarr(nlon,nlat)
         rainnc_total_hist   = fltarr(nlon,nlat)
         event_total_hist    = fltarr(nlon,nlat)
         frac_event_hist     = fltarr(nlon,nlat)
         nPrecipDays_future  = fltarr(nlon,nlat)
         nEvent_all_future   = fltarr(nlon,nlat)
         nEvent_event_future = fltarr(nlon,nlat)
         rainnc_total_future = fltarr(nlon,nlat)
         event_total_future  = fltarr(nlon,nlat)
         frac_event_future   = fltarr(nlon,nlat)

         for ij=0,nlon-1 do begin
            for ik=0,nlat-1 do begin
               ; ###########################
               ; Historical 
               ; ###########################
               ; Does this gridpoint observe precipitation?
               Wet = where(rainnc_precip(ij,ik,*) gt 0)
               nWet = n_elements(Wet)

               ; If so, look at events and see how they relate to statistically
               ; large precipitation events
               if (Wet(0) ne -1 and nWet gt 1) then begin
                  nPrecipDays_hist(ij,ik)   = nWet
                  ; Create CDF for precipitation, use percentile to composite
                  y=histogram(rainnc_precip(ij,ik,Wet),location=x,min=0,max=max(rainnc_precip(ij,ik,Wet)),nbins=200)
                  yCDF=total(y,/cumulative)
               
                  ; Compute frequencies
                  nEvent_all_hist(ij,ik)   = n_elements(where(rainnc_precip(ij,ik,Wet)   gt interpol(x,yCDF/max(yCDF),ptile(ip2))))
                  nEvent_event_hist(ij,ik) = n_elements(where(event_precip_hist(ij,ik,*) gt interpol(x,yCDF/max(yCDF),ptile(ip2))))
                  frac_event_hist(ij,ik)   = nEvent_event_hist(ij,ik)/nEvent_all_hist(ij,ik)
                  rainnc_total_hist(ij,ik) = total(rainnc_precip(ij,ik,Wet))
                  event_total_hist(ij,ik)  = total(event_precip_hist(ij,ik,*))
                  ; Do the same aggreation (event vs total) for IVT, SWE, zLCL and z500
                  
               endif

               ; ###########################
               ; Future
               ; ###########################
               ; Does this gridpoint observe precipitation?
               Wet = where(rainncf_precip(ij,ik,*) gt 0)
               nWet = n_elements(Wet)

               ; If so, look at events and see how they relate to statistically
               ; large precipitation events
               if (Wet(0) ne -1 and nWet gt 1) then begin
                  nPrecipDays_future(ij,ik)   = nWet
                  ; Create CDF for precipitation, use percentile to composite
                  y=histogram(rainncf_precip(ij,ik,Wet),location=x,min=0,max=max(rainncf_precip(ij,ik,Wet)),nbins=200)
                  yCDF=total(y,/cumulative)
            
                  ; Compute frequencies
                  nEvent_all_future(ij,ik)   = n_elements(where(rainncf_precip(ij,ik,Wet)    gt interpol(x,yCDF/max(yCDF),ptile(ip2))))
                  nEvent_event_future(ij,ik) = n_elements(where(event_precip_future(ij,ik,*) gt interpol(x,yCDF/max(yCDF),ptile(ip2))))
                  frac_event_future(ij,ik)   = nEvent_event_future(ij,ik)/nEvent_all_future(ij,ik)
                  rainnc_total_future(ij,ik) = total(rainncf_precip(ij,ik,Wet))
                  event_total_future(ij,ik)  = total(event_precip_future(ij,ik,*))
               endif
            end
         end

         ; Write output
         print,'Writing to output...'
         fileOUT = dirOUT+'composite.frequency.'+string(phold(ip1),format='(i2.2)')+'hrs.'+$
                   string(100.*ptile(ip2),format='(f5.2)')+'ptile.'+modelID(imod)+'.nc'
         fileID  = ncdf_create(fileOUT,/clobber)
         dimID1  = ncdf_dimdef(fileID,'lon',nlon)
         dimID2  = ncdf_dimdef(fileID,'lat',nlat)
         dimID3  = ncdf_dimdef(fileID,'time',ntime)
         varID1  = ncdf_vardef(fileID,'lon',                 [dimID1,dimID2],/float)
         varID2  = ncdf_vardef(fileID,'lat',                 [dimID1,dimID2],/float)
         varID3  = ncdf_vardef(fileID,'frac_event_hist',     [dimID1,dimID2],/float)
         varID4  = ncdf_vardef(fileID,'event_precip_hist',   [dimID1,dimID2],/float)
         varID5  = ncdf_vardef(fileID,'total_precip_hist',   [dimID1,dimID2],/float)
         varID6  = ncdf_vardef(fileID,'frac_event_future',   [dimID1,dimID2],/float)
         varID7  = ncdf_vardef(fileID,'event_precip_future', [dimID1,dimID2],/float)
         varID8  = ncdf_vardef(fileID,'total_precip_future', [dimID1,dimID2],/float)
         varID9  = ncdf_vardef(fileID,'event_ivt_hist',      [dimID1,dimID2],/float)
         varID10 = ncdf_vardef(fileID,'event_ivt_future',    [dimID1,dimID2],/float)
         varID11 = ncdf_vardef(fileID,'event_swe_hist',      [dimID1,dimID2],/float)
         varID12 = ncdf_vardef(fileID,'event_swe_future',    [dimID1,dimID2],/float)
         varID13 = ncdf_vardef(fileID,'event_z500_hist',     [dimID1,dimID2],/float)
         varID14 = ncdf_vardef(fileID,'event_z500_future',   [dimID1,dimID2],/float)
         varID15 = ncdf_vardef(fileID,'event_zLCL_hist',     [dimID1,dimID2],/float)
         varID16 = ncdf_vardef(fileID,'event_zLCL_future',   [dimID1,dimID2],/float)         
         ncdf_control,fileID,/endef
         ncdf_varput,fileID,varID1,lon
         ncdf_varput,fileID,varID2,lat
         ncdf_varput,fileID,varID3,frac_event_hist
         ncdf_varput,fileID,varID4,event_total_hist
         ncdf_varput,fileID,varID5,rainnc_total_hist
         ncdf_varput,fileID,varID6,frac_event_future
         ncdf_varput,fileID,varID7,event_total_future
         ncdf_varput,fileID,varID8,rainnc_total_future
         ncdf_varput,fileID,varID9,total(event_ivt_hist,3)
         ncdf_varput,fileID,varID10,total(event_ivt_future,3)
         ncdf_varput,fileID,varID11,total(event_swe_hist,3)
         ncdf_varput,fileID,varID12,total(event_swe_future,3)
         ncdf_varput,fileID,varID13,total(event_z500_hist,3)
         ncdf_varput,fileID,varID14,total(event_z500_future,3)
         ncdf_varput,fileID,varID15,total(event_zLCL_hist,3)
         ncdf_varput,fileID,varID16,total(event_zLCL_future,3)

         
         ncdf_close,fileID
      end                       ; Model
   end                          ; Persistence
end                             ; Strength



; END PROGRAM
end

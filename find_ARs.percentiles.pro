; #########################################################################
; The purpose of this program is to detect anomalously high IVT events
; at the coast. Events are flagged when anomalously high IVT is
; observed (ivtThresh) for more an extended period of time
; (persistenceThresh).
;
;
; Name                Type    Size    Dscription
;
; Inputs
; modelID             string  [1]     Shortname for model (e.g. 'gfdl')
; persistenceThresh   int     [1]     Persistence threshold (hours)
; ivtThresh           float   [1]     Strength threshold (percentile)
; yearS               int     [1]     Starting year for event detection.
; yearF               int     [1]     Ending year for event detection.
; dirOUT              string  [1]     Output directory location.

; Output is written to netCDF and stored in dirOUT.
;
; #########################################################################
pro find_ARs,modelID,persistenceThresh,ivtThresh,yearS,yearF,dirOUT,nevents
writeOutput=1
if (strcmp(dirOUT,'')) then begin
   writeOutput=0
endif

; IVT threshold 
;ivtThresh=99
; How long are AR events? (in hours)
;persistenceThresh = 24
; Which RCM
;modelID='hadgem'
; Temporal configuration
;yearS=1980
;yearF=2010;
; Output directory
;dirOUT = 'data/events/'

; 50 or 25km data?
res = '50'
dt = 3                          ; time interval (hours)

; Make plots of stats?
mkplots=0

if (yearS lt 2010 and yearF le 2010) then period='clim' else period='future'

; Read in WRF grid from NA-CORDEX experiments
fileIN = '/Projects/HydroMet/dswales/NA-CORDEX/WRF_erain_'+string(res,format='(i2)')+'km/wrfout_198001.nc'
fileID = ncdf_open(fileIN)
ncdf_varget,fileID,ncdf_varid(fileID,'XLONG'),lon
ncdf_varget,fileID,ncdf_varid(fileID,'XLAT'),lat
ncdf_close,fileID

; Read in slab indices
openr,102,'data/coastalSlabXY.coast1.'+res+'.txt'
readf,102,npts
xi = intarr(npts)
yi = intarr(npts)
readf,102,xi
readf,102,yi
close,102

; What are the lon,lats for the slab points?
lon_slab = fltarr(npts)
lat_slab = fltarr(npts)
for ij=0,npts-1 do begin
   lon_slab(ij) = lon(xi(ij),yi(ij))
   lat_slab(ij) = lat(xi(ij),yi(ij))
end

; Read in IVT data
file = 'data/coastalSlabXY.coast1.'+res+'.'+modelID+'.nc'
file = '/data/dswales/NA-CORDEX/ARdet/slabs/coastalSlabXY.coast1.'+res+'.'+modelID+'.nc'
fileID = ncdf_open(file)
ncdf_varget,fileID,ncdf_varid(fileID,'year'),year
ncdf_varget,fileID,ncdf_varid(fileID,'month'),month
ncdf_varget,fileID,ncdf_varid(fileID,'day'),day
ncdf_varget,fileID,ncdf_varid(fileID,'hour'),hour
ncdf_varget,fileID,ncdf_varid(fileID,'ivtU'),ivtU
ncdf_varget,fileID,ncdf_varid(fileID,'ivtV'),ivtV
ncdf_close,fileID

; Only do stats for cool-season months
cool_season = where(year ge yearS and year le yearF and (month le 3 or month ge 10))
nTimes = n_elements(cool_season)

; Compute IVT and some stats
ivt = sqrt(ivtU(*,cool_season)*ivtU(*,cool_season)+ivtV(*,cool_season)*ivtV(*,cool_season))
ivtThreshA = [50,90,95,ivtThresh]
nThresh   = n_elements(ivtThreshA)
ivtThreshSlab = fltarr(npts,nThresh)
for ij=0,npts-1 do begin
   y=histogram(ivt(ij,*),location=x,min=0,max=1000,nbins=99)
   yCDF=total(y,/cumulative)/total(y)
   for ik=0,nThresh-1 do ivtThreshSlab(ij,ik)=interpol(x,yCDF,ivtThreshA(ik)/100.) 
end

; What does the climatology look like, make PDFs of IVT at each point.
if (mkplots) then begin
   window,xsize=1000,ysize=800
   !p.multi=[0,1,2,0,0]
   cgloadct,50,ncolors=npts
   ; IVT stats
   cgplot,lat_slab,ivtThreshSlab(*,0),xtitle='Latitude',ytitle='IVT',xrange=[50,30],yrange=[0,800],thick=3
   if (nThresh gt 1) then begin
      for ik=1,nThresh-1 do begin
         cgplot,lat_slab,ivtThreshSlab(*,ik),/overplot,color='grey',thick=2
      end
   end

   ; PDF plot
   for ij=0,npts-1 do begin
      y = histogram(ivt(ij,*),location=x,min=10,max=1000,nbins=200)
      color='dodger_blue'
      if(lat_slab(ij) lt 40) then color='green'
      if (ij eq 0) then cgplot,x,smooth(y/total(y),3),psym=-3,yrange=[0,.05],xrange=[0,800],color=ij,$;color,$
                               ytitle='Frequency',xtitle="IVT" 
      if (ij ne 0) then cgplot,x,smooth(y/total(y),3),psym=-3,/overplot,color=ij;color
   end
   saveimage,'plots/stats.coastalSlabXY.coast1.'+res+'.'+modelID+'.png'
endif

; #######################################################################
; BEGIN AR detection.
; Generally this is based on "How often do we observe IVT' persistence
; at the coast?"
;
; #######################################################################

; First, make mask for when anamolously high IVT occuring at the coast
iflag = intarr(npts,nTimes)
for ij=0,npts-1 do begin
   bigJuice = where(ivt(ij,*) gt ivtThreshSlab(ij,nThresh-1))
   if (bigJuice(0) ne -1) then iflag(ij,bigJuice) = 1
end

; Now, identify these events in the data record.
nevents    = 0
eventCount = 0
init = 1
for ij=0,nTimes-1 do begin
   ; Are we in the midst of an event? If so, make note
   if (total(iflag(*,ij)) gt 0) then eventCount = eventCount+1

   ; Has an event just ended? If so, save start/stop times.
   if (eventCount ne 0 and total(iflag(*,ij)) eq 0) then begin
      event_len = julday(month(cool_season(ij-1)), day(cool_season(ij-1)),$
                         year(cool_season(ij-1)),hour(cool_season(ij-1)))-$
                  julday(month(cool_season(ij-eventCount)), day(cool_season(ij-eventCount)),$
                         year(cool_season(ij-eventCount)),hour(cool_season(ij-eventCount)))
      if (event_len lt 10) then begin ; Just to avoid discontinuities in dataset (e.g. Mar-Oct data gap)
         ; Is this event long enough?
         if (eventCount gt persistenceThresh/dt) then begin
            ; Archive start/end times for event
            if (not init) then begin
               event_year   = [[event_year], [year(cool_season(ij-eventCount)), year(cool_season(ij-1))]]
               event_month  = [[event_month],[month(cool_season(ij-eventCount)),month(cool_season(ij-1))]]
               event_day    = [[event_day],  [day(cool_season(ij-eventCount)),  day(cool_season(ij-1))]]
               event_hour   = [[event_hour], [hour(cool_season(ij-eventCount)), hour(cool_season(ij-1))]]
               event_length = [event_length, event_len]
               event_lat    = [[event_lat],[mean(lat_slab(where(iflag(*,ij-eventCount) gt 0))),mean(lat_slab(where(iflag(*,ij-1) gt 0)))]]
            endif
            if (init) then begin
               event_year   = [year(cool_season(ij-eventCount)), year(cool_season(ij-1))]
               event_month  = [month(cool_season(ij-eventCount)),month(cool_season(ij-1))]
               event_day    = [day(cool_season(ij-eventCount)),  day(cool_season(ij-1))]
               event_hour   = [hour(cool_season(ij-eventCount)), hour(cool_season(ij-1))]
               event_length = event_len
               event_lat    = [mean(lat_slab(where(iflag(*,ij-eventCount) gt 0))),mean(lat_slab(where(iflag(*,ij-1) gt 0)))]
               init = 0
            endif
            nevents = nevents+1
         endif
      endif
      eventCount = 0
   endif
end

; #######################################################################
; Print some diagnostics
; #######################################################################
;event_dlat = event_lat(0,*)-event_lat(1,*)
; Quasi stationary events(?)
;event_qs   = where(event_dlat lt 3)
; "Sweeping" events
;event_sw   = where(event_dlat ge 3)

;
if (nevents gt 0) then begin
   a=julday(event_month,event_day,event_year,event_hour)
   printf,101,FORMAT='(A6,I16,F16.2,2I16,F16.2)',modelID,persistenceThresh,ivtThresh,nevents,ceil(total(a(1,*)-a(0,*))),(total(a(1,*)-a(0,*)))/54.75 
   print,FORMAT='(A6,I16,F16.2,2I16,F16.2)',modelID,persistenceThresh,ivtThresh,nevents,ceil(total(a(1,*)-a(0,*))),(total(a(1,*)-a(0,*)))/54.75
endif
if (nevents eq 0) then begin
   printf,101,FORMAT='(A6,I16,F16.2,2I16,F16.2)',modelID,persistenceThresh,ivtThresh,0,0,0.
   print,FORMAT='(A6,I16,F16.2,2I16,F16.2)',modelID,persistenceThresh,ivtThresh,0,0,0.
endif

; #######################################################################
; Write to output
; #######################################################################
if (nevents gt 0 and writeOutput) then begin
   if (strcmp(period,'clim')) then $
      fileOUT = 'events.'+res+'km.'+string(ivtThresh,format='(f5.2)')+'ptile.'+string(persistenceThresh,format='(i2.2)')+'hrs.'+modelID+'.nc'
   if (strcmp(period,'future')) then $
      fileOUT = 'events.future.'+res+'km.'+string(ivtThresh,format='(f5.2)')+'ptile.'+string(persistenceThresh,format='(i2.2)')+'hrs.'+modelID+'.nc'
   fileID  = ncdf_create(dirOUT+fileOUT,/clobber)
   dimID1  = ncdf_dimdef(fileID,'nEvents',nevents)
   varID1  = ncdf_vardef(fileID,'year_start', [dimID1],/long)
   varID2  = ncdf_vardef(fileID,'year_end',   [dimID1],/long)
   varID3  = ncdf_vardef(fileID,'month_start',[dimID1],/long)
   varID4  = ncdf_vardef(fileID,'month_end',  [dimID1],/long)
   varID5  = ncdf_vardef(fileID,'day_start',  [dimID1],/long)
   varID6  = ncdf_vardef(fileID,'day_end',    [dimID1],/long)
   varID7  = ncdf_vardef(fileID,'hour_start', [dimID1],/long)
   varID8  = ncdf_vardef(fileID,'hour_end',   [dimID1],/long)
   varID9  = ncdf_vardef(fileID,'lat_start',  [dimID1],/float)
   varID10 = ncdf_vardef(fileID,'lat_end',    [dimID1],/float)
   varID11 = ncdf_vardef(fileID,'length',     [dimID1],/long)
   ncdf_control,fileID,/endef
   ncdf_varput,fileID,varID1,reform(event_year(0,*))
   ncdf_varput,fileID,varID2,reform(event_year(1,*))
   ncdf_varput,fileID,varID3,reform(event_month(0,*))
   ncdf_varput,fileID,varID4,reform(event_month(1,*))
   ncdf_varput,fileID,varID5,reform(event_day(0,*))
   ncdf_varput,fileID,varID6,reform(event_day(1,*))
   ncdf_varput,fileID,varID7,reform(event_hour(0,*))
   ncdf_varput,fileID,varID8,reform(event_hour(1,*))
   ncdf_varput,fileID,varID9,reform(event_lat(0,*))
   ncdf_varput,fileID,varID10,reform(event_lat(1,*))
   ncdf_varput,fileID,varID11,event_length*24.
   ncdf_close,fileID
endif


; Plots
; a) Domain
;cgmap_set,/continents,/usa,limit=[30,-140,55,-110],/grid,/isotropic
;for ij=0,npts-1 do begin
;   cgplot,lon(xi(ij),yi(ij)),lat(xi(ij),yi(ij)),psym=2,/overplot,color='red'
;end

; END PROGRAM
end

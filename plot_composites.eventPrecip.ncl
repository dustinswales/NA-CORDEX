; #########################################################################
;
;#########################################################################
;
; #########################################################################
; BEGIN
; #########################################################################
begin

; Which variable to plot? Precip is default
doPrecip = 0
doIVT  = 01
doZ500 = 0
doZLCL = 0
doSWE  = 0

; Which RCM?
modelID = (/"erain","hadgem","gfdl","mpi"/)
nmods   = dimsizes(modelID)

; Switch for future epoch composites
future=01

; Extreme parameters
;ptile       = (/"95.00","97.00"/)
ptile = "95.00"
np1 = dimsizes(ptile)
;persistence = (/"24","36"/)
persistence = "24"
np2 = dimsizes(persistence)

; Composite location
dirIN="/data/dswales/NA-CORDEX/ARdet/composites/"

; Output location
dirOUT = "/home/dswales/Projects/NA-CORDEX/scripts/ARdet/plots/composites/"

; #########################################################################
; Plot configuration
; #########################################################################
wks_type = "png"
nLevels  = 40.0

; What varaible to plot? (default is precip)
if (doPrecip .eq. 1) then
   varName  = "precip"
   varRange = (/1.,10./)                         
   ctable   = "NCV_bright"
   units    = "mm"
end if
if (doIVT .eq. 1) then
   varName  = "ivt"
   varRange = (/100.,400./)                       
   ctable   = "NCV_bright"
   units    = "Kg/m/s"
end if
;
if (doZ500 .eq. 1) then
   varName  = "z500"
   varRange = (/5600.,6200./)                       
   ctable   = "NCV_bright"
   units    = "m"
end if
;
if (doZLCL .eq. 1) then
   varName  = "zLCL"
   varRange = (/0.,2000./)                       
   ctable   = "NCV_bright"
   units    = "m"
end if
;
if (doSWE .eq. 1) then
   varName  = "SWE"
   varRange = (/0.,20./)                       
   ctable   = "NCV_bright"
   units    = "m"
end if

; Variable spacing.
dVarRange = ((varRange(1)-varRange(0))/nLevels)

; Load colortables for plots
ctableD   = "NCV_blue_red"				; Sequential
ctable    = "NCV_bright"				; Divergent
cmap      = read_colormap_file(ctable)
cmap(0,:) = namedcolor2rgba("white")
ncmap     = dimsizes(cmap(:,0))                      ; # colors
cmapD     = read_colormap_file(ctableD)
ncmapD    = dimsizes(cmapD(:,0))                      ; # colors

; Configuration for plot.
res                         = True
res@gsnDraw                 = False                    ; don't draw
res@gsnFrame                = False                    ; don't advance frame
res@cnInfoLabelOn           = False                    ; turn off cn info label
res@cnFillOn                = True                     ; turn on color
res@lbLabelBarOn            = False                    ; turn off individual cb's
res@cnLinesOn               = False                    ; Turn off lines
res@cnLineLabelsOn          = False                    ; Turn off line labels
res@tmYRLabelsOn            = False                    ; no right labels
res@tmYROn                  = False                    ; no right tickmarks
res@tmXTOn                  = False                    ; no top tickmarks
res@cnLevelSelectionMode    = "ManualLevels"           ; Manually set levels
res@mpOutlineBoundarySets   = "GeophysicalAndUSStates" ; Plot boundary set
res@gsnPaperOrientation     = "landscape"              ; orientation
res@mpMaxLatF               = 75.                      ; specify the plot domain
res@mpMinLatF               = 5.                       ;                         
res@mpMinLonF               = 205.                     ;
res@mpMaxLonF               = 270.                     ;
res@mpCenterLonF            = 235.                     ;
res@cnMissingValFillPattern = 17                       ; set the missing value fill pattern to 17
res@tmXBLabelFontHeightF    = 0.015                    ; resize tick labels
res@tmYLLabelFontHeightF    = 0.015		       ;

; Load variable specific configuration
;
resC                        = res
resC@cnFillPalette          = cmap(:ncmap-3,:)         ; skip last two colors
resC@cnMinLevelValF         = varRange(0)              ; Lower range
resC@cnMaxLevelValF         = varRange(1)              ; Uppermost level
resC@cnLevelSpacingF        = dVarRange                ; level spacing


; Loop over all thresholds
do ip1=0,np1-1
do ip2=0,np2-1

; #########################################################################
; Loop over all models and create plots.
; #########################################################################
plot    = new(nmods*3,graphic)
if (future .eq. 0) then
   fileOUT = dirOUT+"composite."+varName+".events.historical.MultiModel."+persistence(ip2)+"hrs."+ptile(ip1)+"ptile"+".v3"
end if
if (future .eq. 1) then
   fileOUT = dirOUT+"composite."+varName+".events.future.MultiModel."+persistence(ip2)+"hrs."+ptile(ip1)+"ptile"+".v3"
end if
print(fileOUT+"."+wks_type)

wks     = gsn_open_wks(wks_type,fileOUT); send graphics to PNG file
do ij=0,nmods-1
   ; ##########################################   
   ; Read in data
   ; ##########################################
   fileIN = dirIN+"composite.frequency."+persistence(ip2)+"hrs."+ptile(ip1)+"ptile."+modelID(ij)+".nc"
   fileID = addfile(fileIN,"r")
   lon                   = fileID->lon(5::,5::)
   lat                   = fileID->lat(5::,5::)
   nTime_event_hist      = fileID->ntime_event_hist
   nTime_hist            = fileID->ntime_hist
   nTime_event_future    = fileID->ntime_event_future
   nTime_future          = fileID->ntime_future
   if (doPrecip .eq. 1) then
      precip_total_hist     = fileID->total_precip_hist(5::,5::)
      precip_event_hist     = fileID->event_precip_hist(5::,5::)
      precip_nonevent_hist  = fileID->nevent_precip_hist(5::,5::)
      precip_total_future   = fileID->total_precip_future(5::,5::)
      precip_event_future   = fileID->event_precip_future(5::,5::)
      precip_nonevent_future= fileID->nevent_precip_future(5::,5::)
   end if
   if (doIVT .eq. 1) then
      ivt_total_hist       = fileID->total_ivt_hist(5::,5::)
      ivt_event_hist       = fileID->event_ivt_hist(5::,5::)
      ivt_nonevent_hist    = fileID->nevent_ivt_hist(5::,5::)
      ivt_total_future     = fileID->total_ivt_future(5::,5::)
      ivt_event_future     = fileID->event_ivt_future(5::,5::)
      ivt_nonevent_future  = fileID->nevent_ivt_future(5::,5::)
      ;
      ;ivt_total_hist@_FillValue = -999.
      ;ivt_total_hist = where(abs(ivt_total_hist).gt.1e20, ivt_total_hist@_FillValue, ivt_total_hist)
      ;
      ;ivt_nonevent_hist@_FillValue = -999.
      ;ivt_nonevent_hist = where(abs(ivt_nonevent_hist).gt.1e20, ivt_nonevent_hist@_FillValue, ivt_nonevent_hist)
      ;
      ;ivt_total_future@_FillValue = -999.
      ;ivt_total_future = where(abs(ivt_total_future).gt.1e20, ivt_total_future@_FillValue, ivt_total_future)
      ;
      ;ivt_nonevent_future@_FillValue = -999.
      ;ivt_nonevent_future = where(abs(ivt_nonevent_future).gt.1e20, ivt_nonevent_future@_FillValue, ivt_nonevent_future)
   end if
   if (doZ500 .eq. 1) then
      z500_total_hist       = fileID->total_z500_hist(5::,5::)
      z500_event_hist       = fileID->event_z500_hist(5::,5::)
      z500_nonevent_hist    = fileID->nevent_z500_hist(5::,5::)
      z500_total_future     = fileID->total_z500_future(5::,5::)
      z500_event_future     = fileID->event_z500_future(5::,5::)
      z500_nonevent_future  = fileID->nevent_z500_future(5::,5::)
   end if
   if (doZLCL .eq. 1) then
      zLCL_total_hist       = fileID->total_zLCL_hist(5::,5::)
      zLCL_event_hist       = fileID->event_zLCL_hist(5::,5::)
      zLCL_nonevent_hist    = fileID->nevent_zLCL_hist(5::,5::)
      zLCL_total_future     = fileID->total_zLCL_future(5::,5::)
      zLCL_event_future     = fileID->event_zLCL_future(5::,5::)
      zLCL_nonevent_future  = fileID->nevent_zLCL_future(5::,5::)
   end if
   if (doSWE .eq. 1) then
      swe_total_hist       = fileID->total_swe_hist(5::,5::)
      swe_event_hist       = fileID->event_swe_hist(5::,5::)
      swe_nonevent_hist    = fileID->nevent_swe_hist(5::,5::)
      swe_total_future     = fileID->total_swe_future(5::,5::)
      swe_event_future     = fileID->event_swe_future(5::,5::)
      swe_nonevent_future  = fileID->nevent_swe_future(5::,5::)
   end if

   ; ##########################################
   ; Add attributes
   ; ##########################################
   if (doPrecip .eq. 1) then
      precip_total_hist  = precip_total_hist*8.	; Express as mean daily precipitation	
      precip_total_hist@lat2d = lat
      precip_total_hist@lon2d = lon
      precip_total_hist@units = " "
      ;
      precip_event_hist  = precip_event_hist*8.	; Express as mean daily precipitation
      precip_event_hist@lat2d = lat
      precip_event_hist@lon2d = lon
      precip_event_hist@units = " "
      ;
      precip_nonevent_hist  = precip_nonevent_hist*8.	; Express as mean daily precipitation
      precip_nonevent_hist@lat2d = lat
      precip_nonevent_hist@lon2d = lon
      precip_nonevent_hist@units = " "
      ;
      precip_total_future = precip_total_future*8.	; Express as mean daily precipitation
      precip_total_future@lat2d = lat
      precip_total_future@lon2d = lon
      precip_total_future@units = " "
      ;
      precip_event_future = precip_event_future*8.	; Express as mean daily precipitation	
      precip_event_future@lat2d = lat
      precip_event_future@lon2d = lon
      precip_event_future@units = " "
      ;
      precip_nonevent_future  = precip_nonevent_future*8.	; Express as mean daily precipitation
      precip_nonevent_future@lat2d = lat
      precip_nonevent_future@lon2d = lon
      precip_nonevent_future@units = " "
   end if
   ;
   if (doIVT .eq. 1) then
      ivt_total_hist  = ivt_total_hist
      ivt_total_hist@lat2d = lat
      ivt_total_hist@lon2d = lon
      ivt_total_hist@units = " "
      ;
      ivt_event_hist  = ivt_event_hist
      ivt_event_hist@lat2d = lat
      ivt_event_hist@lon2d = lon
      ivt_event_hist@units = " "
      ;
      ivt_total_future = ivt_total_future
      ivt_total_future@lat2d = lat
      ivt_total_future@lon2d = lon
      ivt_total_future@units = " "
      ;
      ivt_event_future = ivt_event_future
      ivt_event_future@lat2d = lat
      ivt_event_future@lon2d = lon
      ivt_event_future@units = " "
      ;
      ivt_nonevent_hist       = ivt_nonevent_hist
      ivt_nonevent_hist@lat2d = lat
      ivt_nonevent_hist@lon2d = lon
      ivt_nonevent_hist@units = " "
      ;
      ivt_nonevent_future       = ivt_nonevent_future
      ivt_nonevent_future@lat2d = lat
      ivt_nonevent_future@lon2d = lon
      ivt_nonevent_future@units = " "
   end if
   ;
   if (doZ500 .eq. 1) then
      z500_total_hist  = z500_total_hist	
      z500_total_hist@lat2d = lat
      z500_total_hist@lon2d = lon
      z500_total_hist@units = " "
      ;
      z500_event_hist  = z500_event_hist
      z500_event_hist@lat2d = lat
      z500_event_hist@lon2d = lon
      z500_event_hist@units = " "
      ;
      z500_total_future = z500_total_future
      z500_total_future@lat2d = lat
      z500_total_future@lon2d = lon
      z500_total_future@units = " "
      ;
      z500_event_future = z500_event_future
      z500_event_future@lat2d = lat
      z500_event_future@lon2d = lon
      z500_event_future@units = " "
      ;
      z500_nonevent_hist       = z500_nonevent_hist
      z500_nonevent_hist@lat2d = lat
      z500_nonevent_hist@lon2d = lon
      z500_nonevent_hist@units = " "
      ;
      z500_nonevent_future       = z500_nonevent_future
      z500_nonevent_future@lat2d = lat
      z500_nonevent_future@lon2d = lon
      z500_nonevent_future@units = " "
   end if
   ;
   if (doZLCL .eq. 1) then
      zLCL_total_hist  = zLCL_total_hist	
      zLCL_total_hist@lat2d = lat
      zLCL_total_hist@lon2d = lon
      zLCL_total_hist@units = " "
      ;
      zLCL_event_hist  = zLCL_event_hist
      zLCL_event_hist@lat2d = lat
      zLCL_event_hist@lon2d = lon
      zLCL_event_hist@units = " "
      ;
      zLCL_total_future = zLCL_total_future
      zLCL_total_future@lat2d = lat
      zLCL_total_future@lon2d = lon
      zLCL_total_future@units = " "
      ;
      zLCL_event_future = zLCL_event_future
      zLCL_event_future@lat2d = lat
      zLCL_event_future@lon2d = lon
      zLCL_event_future@units = " "
      ;
      zLCL_nonevent_hist       = zLCL_nonevent_hist
      zLCL_nonevent_hist@lat2d = lat
      zLCL_nonevent_hist@lon2d = lon
      zLCL_nonevent_hist@units = " "
      ;
      zLCL_nonevent_future       = zLCL_nonevent_future
      zLCL_nonevent_future@lat2d = lat
      zLCL_nonevent_future@lon2d = lon
      zLCL_nonevent_future@units = " "
   end if
   ;
   if (doSWE .eq. 1) then
      swe_total_hist  = swe_total_hist	
      swe_total_hist@lat2d = lat
      swe_total_hist@lon2d = lon
      swe_total_hist@units = " "
      ;
      swe_event_hist  = swe_event_hist
      swe_event_hist@lat2d = lat
      swe_event_hist@lon2d = lon
      swe_event_hist@units = " "
      ;
      swe_total_future = swe_total_future
      swe_total_future@lat2d = lat
      swe_total_future@lon2d = lon
      swe_total_future@units = " "
      ;
      swe_event_future = swe_event_future
      swe_event_future@lat2d = lat
      swe_event_future@lon2d = lon
      swe_event_future@units = " "
      ;
      swe_nonevent_hist       = swe_nonevent_hist
      swe_nonevent_hist@lat2d = lat
      swe_nonevent_hist@lon2d = lon
      swe_nonevent_hist@units = " "
      ;
      swe_nonevent_future       = swe_nonevent_future
      swe_nonevent_future@lat2d = lat
      swe_nonevent_future@lon2d = lon
      swe_nonevent_future@units = " "
   end if
   
   ; ##########################################
   ; Make plots
   ; ##########################################

   ; Which variable are we plotting? (default is precip)
   if (doPrecip .eq. 1) then
      var_total_hist      = precip_total_hist
      var_event_hist      = precip_event_hist
      var_nonevent_hist   = precip_nonevent_hist
      var_total_future    = precip_total_future
      var_event_future    = precip_event_future
      var_nonevent_future = precip_nonevent_future
   end if
   if (doIVT .eq. 1) then
      var_total_hist      = ivt_total_hist
      var_event_hist      = ivt_event_hist
      var_nonevent_hist   = ivt_nonevent_hist
      var_total_future    = ivt_total_future
      var_event_future    = ivt_event_future
      var_nonevent_future = ivt_nonevent_future
   end if
   ;
   if (doZ500 .eq. 1) then
      var_total_hist      = z500_total_hist
      var_event_hist      = z500_event_hist
      var_nonevent_hist   = z500_nonevent_hist
      var_total_future    = z500_total_future
      var_event_future    = z500_event_future
      var_nonevent_future = z500_nonevent_future
   end if
   ;
   if (doZLCL .eq. 1) then
      var_total_hist      = zLCL_total_hist
      var_event_hist      = zLCL_event_hist
      var_nonevent_hist   = zLCL_nonevent_hist
      var_total_future    = zLCL_total_future
      var_event_future    = zLCL_event_future
      var_nonevent_future = zLCL_nonevent_future
   end if
   ;
   if (doSWE .eq. 1) then
      var_total_hist      = swe_total_hist
      var_event_hist      = swe_event_hist
      var_nonevent_hist   = swe_nonevent_hist
      var_total_future    = swe_total_future
      var_event_future    = swe_event_future
      var_nonevent_future = swe_nonevent_future
   end if   

   resC@tmYLLabelsOn = False
   resC@tmXBLabelsOn = False
   resC@tmYLLabelsOn = False
   resC@tmXBLabelsOn = True
   if (ij .eq. 0) then
      resC@tmYLLabelsOn = True
      resC@tmYLLabelsOn = True
   end if
   resC@gsnCenterString   = ""
   if (future .eq. 0) then
      resC@tiMainString = "Cummulative ("+modelID(ij)+")"
      plot(ij)          = gsn_csm_contour_map(wks,var_total_hist,resC)
      resC@tiMainString = "Event-Days ("+modelID(ij)+")"
      plot(ij+nmods)    = gsn_csm_contour_map(wks,var_event_hist,resC)
      resC@tmXBLabelsOn = True
      resC@tiMainString = "non-Event Days ("+modelID(ij)+")"
      plot(ij+nmods*2)  = gsn_csm_contour_map(wks,var_nonevent_hist,resC)
   end if
   if (future .eq. 1) then
      if (ij .eq. 0) then
      	 var_total_future    = 0
	 var_event_future    = 0
	 var_nonevent_future = 0
      end if 
      resC@tiMainString = "Cummulative ("+modelID(ij)+")"
      plot(ij)          = gsn_csm_contour_map(wks,var_total_future,resC)
      resC@tiMainString = "Event-Days ("+modelID(ij)+")"
      plot(ij+nmods)    = gsn_csm_contour_map(wks,var_event_future,resC)
      resC@tmXBLabelsOn = True
      resC@tiMainString = "non-Event-Days ("+modelID(ij)+")"
      plot(ij+nmods*2)  = gsn_csm_contour_map(wks,var_nonevent_future,resC)
   end if
end do

; #########################################################################
; Make panel plot (1 upper row, 2 lower rows)
; #########################################################################
; Panel configuration
res_L                   = True
res_L@gsnFrame          = False
res_L@gsnPanelLabelBar  = True
if (doPrecip .eq. 1) then
   res_L@lbTitleString     = "Precipitation (mm/day)"
end if
if (doIVT .eq. 1) then
   res_L@lbTitleString     = "IVT (kg/m/s)"
end if
if (doZ500 .eq. 1) then
   res_L@lbTitleString     = "Geopotential height @ 500hPa (m)"
end if
if (doZLCL .eq. 1) then
   res_L@lbTitleString     = "Lifting Condensation Level (m)"
end if
if (doSWE .eq. 1) then
   res_L@lbTitleString     = "Snow Water Equivalent (m)"
end if
res_L@lbTitleFontHeightF= .01
res_L@pmLabelBarOrthogonalPosF = -0.015              ; move whole thing down a tad

; Panel upper and lower
gsn_panel(wks,plot,(/3,nmods/),res_L)
frame(wks)


end do	; Loop over duration
end do 	; Loop over persistence
; #########################################################################
; END PROGRAM
; #########################################################################
end
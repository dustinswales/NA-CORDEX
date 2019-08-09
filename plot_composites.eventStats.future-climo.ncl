; #########################################################################
;
;#########################################################################
;
; #########################################################################
; BEGIN
; #########################################################################
begin

; Which RCM?
modelID = (/"hadgem","mpi","gfdl"/)
nmods   = dimsizes(modelID)
; Extreme parameters
ptile       = (/"97.00","95.00"/)
np1 = dimsizes(ptile)
persistence = (/"24","36"/)
np2 = dimsizes(persistence)

; Composite location
dirIN="/data/dswales/NA-CORDEX/ARdet/composites/"

; Output location
dirOUT = "/home/dswales/Projects/NA-CORDEX/scripts/ARdet/plots/composites/"

; #########################################################################
; Plot configuration
; #########################################################################
wks_type = "png"
varName  = "precip"
varRangeCumPrecip = (/-100.,100./)                 
varRangeFrac = (/0.,1./)                    
varRangedFrac = (/-50,50/)                    
ctable   = "NCV_bright"
units    = "mm"

; Variable spacing.
nLevels       = 40.0
dVarCumPrecip = ((varRangeCumPrecip(1)-varRangeCumPrecip(0))/nLevels)
dVarFrac      = ((varRangeFrac(1)-varRangeFrac(0))/nLevels)
dVardFrac     = ((varRangedFrac(1)-varRangedFrac(0))/nLevels)

; Load colortables for plots
ctableD      = "NCV_blue_red"
ctable       = "NCV_bright"
cmap      = read_colormap_file(ctable)
cmap(0,:) = namedcolor2rgba("white")
ncmap     = dimsizes(cmap(:,0))                      ; # colors
cmapD     = read_colormap_file(ctableD)
cmapD = cmapD(::-1,:) 
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
res@mpMinLatF               = 5.                      ;                         
res@mpMinLonF               = 200.                     ;
res@mpMaxLonF               = 270.                     ;
res@mpCenterLonF            = 235.                     ;
res@cnMissingValFillPattern = 17                       ; set the missing value fill pattern to 17
res@tmXBLabelFontHeightF    = 0.015                    ; resize tick labels
res@tmYLLabelFontHeightF    = 0.015		       ;

; Load variable specific configuration
res@cnFillPalette           = cmap(:ncmap-3,:)         ; skip last two colors
res@cnMinLevelValF          = varRangeFrac(0)          ; Lower range
res@cnMaxLevelValF          = varRangeFrac(1)          ; Uppermost level
res@cnLevelSpacingF         = dVarFrac                 ; level spacing
;
resC                        = res
resC@cnFillPalette          = cmapD(:ncmap-3,:)         ; skip last two colors
resC@cnMinLevelValF         = varRangeCumPrecip(0)     ; Lower range
resC@cnMaxLevelValF         = varRangeCumPrecip(1)     ; Uppermost level
resC@cnLevelSpacingF        = dVarCumPrecip            ; level spacing
;
resD                        = res
resD@cnFillPalette          = cmapD(:ncmapD-3,:)       ; skip last two colors
resD@cnMinLevelValF         = varRangedFrac(0)         ; Lower range
resD@cnMaxLevelValF         = varRangedFrac(1)         ; Uppermost level
resD@cnLevelSpacingF        = dVardFrac                ; level spacing


; Loop over all thresholds
do ip1=0,np1-1
do ip2=0,np2-1

; #########################################################################
; Loop over all models and create plots.
; #########################################################################
plot1    = new(nmods,graphic)
plot2    = new(nmods*2,graphic)
fileOUT = dirOUT+"composite.eventStats.future-historical.MultiModel."+persistence(ip2)+"hrs."+ptile(ip1)+"ptile"
print(fileOUT+"."+wks_type)
wks     = gsn_open_wks(wks_type,fileOUT); send graphics to PNG file
do ij=0,nmods-1
   ; ##########################################   
   ; Read in data
   ; ##########################################
   ;fileIN = "data/composites/composite.frequency."+modelID(ij)+".nc"
   fileIN = dirIN+"composite.frequency."+persistence(ip2)+"hrs."+ptile(ip1)+"ptile."+modelID(ij)+".nc"
   fileID = addfile(fileIN,"r")
   lon                     = fileID->lon(5::,5::)
   lat                     = fileID->lat(5::,5::)
   precip_total_hist       = fileID->total_precip_hist(5::,5::)
   precip_event_hist       = fileID->event_precip_hist(5::,5::)
   precip_eventfrac_hist   = fileID->frac_event_hist(5::,5::)
   precip_total_future     = fileID->total_precip_future(5::,5::)
   precip_event_future     = fileID->event_precip_future(5::,5::)
   precip_eventfrac_future = fileID->frac_event_future(5::,5::)
   
   ; ##########################################
   ; Compute event fraction precipitation
   ; ##########################################
   event_frac_hist = 100.*precip_event_hist/where(precip_total_hist .eq. 0,-999.,precip_total_hist)
   event_frac_hist@lat2d = lat
   event_frac_hist@lon2d = lon
   event_frac_hist@units = " "
   ;
   event_frac_future = 100.*precip_event_future/where(precip_total_future .eq. 0,-999.,precip_total_future)
   event_frac_future@lat2d = lat
   event_frac_future@lon2d = lon
   event_frac_future@units = " "

   ; ##########################################
   ; Compute changes (future-historical)
   ; ##########################################
   dTotalPrecip = (precip_total_future/(30.*6.))-(precip_total_hist/(31.*6.))
   dTotalPrecip@lat2d = lat
   dTotalPrecip@lon2d = lon
   dTotalPrecip@units = " "
   ;
   dEventFrac = event_frac_future - event_frac_hist
   dEventFrac@lat2d = lat
   dEventFrac@lon2d = lon
   dEventFrac@units = " "
   ;
   dExtremeFrac = 100.*(precip_eventfrac_future - precip_eventfrac_hist)
   dExtremeFrac@lat2d = lat
   dExtremeFrac@lon2d = lon
   dEventFrac@units   = " "   

   ; ##########################################
   ; Make plots
   ; ##########################################
   resD@tmYLLabelsOn = False
   resD@tmXBLabelsOn = False
   resC@tmYLLabelsOn = False
   resC@tmXBLabelsOn = True
   if (ij .eq. 0) then
      resD@tmYLLabelsOn = True
      resC@tmYLLabelsOn = True
   end if
   resC@gsnCenterString   = modelID(ij)
   plot1(ij)        = gsn_csm_contour_map(wks,dTotalPrecip,resC)
   plot2(ij)        = gsn_csm_contour_map(wks,dEventFrac,resD)
   resD@tmXBLabelsOn = True
   plot2(ij+nmods)  = gsn_csm_contour_map(wks,dExtremeFrac,resD)

end do

; #########################################################################
; Make panel plot (1 upper row, 2 lower row)
; #########################################################################
; Upper panel configuration
res_U                   = True
res_U@gsnFrame          = False
res_U@gsnPanelTop       = 0.95
res_U@gsnPanelBottom    = 0.63
res_U@gsnPanelLabelBar  = True
res_U@lbTitleString     = "(mm/month)"
res_U@lbTitleFontHeightF= .01 
res_U@pmLabelBarOrthogonalPosF = -0.015              ; move whole thing down a tad

; Lower panel configuration
res_L                   = True
res_L@gsnFrame          = False
res_L@gsnPanelTop       = 0.60
res_L@gsnPanelBottom    = 0.05
res_L@gsnPanelLabelBar  = True
res_L@lbTitleString     = "Frequency (%)"
res_L@lbTitleFontHeightF= .01
res_L@pmLabelBarOrthogonalPosF = -0.015              ; move whole thing down a tad

; Panel upper and lower
gsn_panel(wks,plot1,(/1,nmods/),res_U)
gsn_panel(wks,plot2,(/2,nmods/),res_L)
frame(wks)



end do	; Loop over duration
end do 	; Loop over persistence
; #########################################################################
; END PROGRAM
; #########################################################################
end
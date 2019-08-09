; #########################################################################
;
;#########################################################################
;
; #########################################################################
; BEGIN
; #########################################################################
begin

; Which RCM?
modelID = (/"erain","hadgem","mpi","gfdl"/)
nmods   = dimsizes(modelID)

; Switch for future epoch composites
future=0

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
varRangeCumPrecip = (/10.,400./)                 
varRangeFrac = (/0.,1./)                    
varRangedFrac = (/-0.3,0.3/)                    
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
res@mpMinLonF               = 205.                     ;
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
resC@cnFillPalette          = cmap(:ncmap-3,:)         ; skip last two colors
resC@cnMinLevelValF         = varRangeCumPrecip(0)     ; Lower range
resC@cnMaxLevelValF         = varRangeCumPrecip(1)     ; Uppermost level
resC@cnLevelSpacingF        = dVarCumPrecip            ; level spacing


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
   precip_total_hist     = fileID->total_precip_hist(5::,5::)
   precip_event_hist     = fileID->event_precip_hist(5::,5::)
   precip_total_future   = fileID->total_precip_future(5::,5::)
   precip_event_future   = fileID->event_precip_future(5::,5::)
   
   ; ##########################################
   ; Add attributes
   ; ##########################################
   precip_total_hist  = precip_total_hist/(31.*6.) ; Express as mean monthly precip	
   precip_total_hist@lat2d = lat
   precip_total_hist@lon2d = lon
   precip_total_hist@units = " "
   ;
   precip_event_hist  = precip_event_hist/(31.*6.) ; Express as mean monthly precip
   precip_event_hist@lat2d = lat
   precip_event_hist@lon2d = lon
   precip_event_hist@units = " "
   ;
   precip_total_future = precip_total_future/(31.*6.) ; Express as mean monthly precip	
   precip_total_future@lat2d = lat
   precip_total_future@lon2d = lon
   precip_total_future@units = " "
   ;
   precip_event_future = precip_event_future/(31.*6.) ; Express as mean monthly precip	
   precip_event_future@lat2d = lat
   precip_event_future@lon2d = lon
   precip_event_future@units = " "

   ; ##########################################
   ; Compute non-event precipitation
   ; ##########################################
   precip_nonevent_hist = precip_total_hist-precip_event_hist
   precip_nonevent_hist@lat2d = lat
   precip_nonevent_hist@lon2d = lon
   precip_nonevent_hist@units = " "
   ;
   precip_nonevent_future = precip_total_future-precip_event_future
   precip_nonevent_future@lat2d = lat
   precip_nonevent_future@lon2d = lon
   precip_nonevent_future@units = " "
   
   ; ##########################################
   ; Make plots
   ; ##########################################
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
      resC@tiMainString= modelID(ij)
      plot(ij)          = gsn_csm_contour_map(wks,precip_total_hist,resC)
      plot(ij+nmods)    = gsn_csm_contour_map(wks,precip_event_hist,resC)
      resC@tmXBLabelsOn = True
      plot(ij+nmods*2)  = gsn_csm_contour_map(wks,precip_nonevent_hist,resC)
   end if
   if (future .eq. 1) then
      if (ij .eq. 0) then
      	 precip_total_future    = 0
	 precip_event_future    = 0
	 precip_nonevent_future = 0
      end if 
      resC@tiMainString= modelID(ij)
      plot(ij)          = gsn_csm_contour_map(wks,precip_total_future,resC)
      plot(ij+nmods)    = gsn_csm_contour_map(wks,precip_event_future,resC)
      resC@tmXBLabelsOn = True
      plot(ij+nmods*2)  = gsn_csm_contour_map(wks,precip_nonevent_future,resC)
   end if
end do

; #########################################################################
; Make panel plot (1 upper row, 2 lower rows)
; #########################################################################
; Panel configuration
res_L                   = True
res_L@gsnFrame          = False
res_L@gsnPanelLabelBar  = True
res_L@lbTitleString     = "Precipitation (mm/month)"
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
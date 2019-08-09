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
varRangeCumPrecip = (/0.,400./)
varRangedPrecip   = (/-100.,100./)                
ctable   = "NCV_bright"
units    = "mm"

; Variable spacing.
nLevels       = 40.0
dVarCumPrecip = ((varRangeCumPrecip(1)-varRangeCumPrecip(0))/nLevels)
dVardPrecip   = ((varRangedPrecip(1)-varRangedPrecip(0))/nLevels)

; Load colortables for plots
ctableD      = "NCV_blue_red"
ctable       = "NCV_bright"
cmap      = read_colormap_file(ctable)
cmap(0,:) = namedcolor2rgba("white")
ncmap     = dimsizes(cmap(:,0))                      ; # colors
cmapD     = read_colormap_file(ctableD)
cmapD     = cmapD(::-1,:) 
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
resD                        = res
resD@cnFillPalette          = cmapD                    ; divergent colormap
resD@cnMinLevelValF         = varRangedPrecip(0)       ; Lower range
resD@cnMaxLevelValF         = varRangedPrecip(1)       ; Uppermost level
resD@cnLevelSpacingF        = dVardPrecip              ; level spacing

; Loop over all thresholds
do ip1=0,np1-1
do ip2=0,np2-1

; #########################################################################
; Loop over all models and create plots.
; #########################################################################
plot    = new(nmods*3,graphic)
fileOUT = dirOUT+"composite."+varName+".events.future-historical.MultiModel."+persistence(ip2)+"hrs."+ptile(ip1)+"ptile"+".v3"
print(fileOUT+"."+wks_type)

wks     = gsn_open_wks(wks_type,fileOUT); send graphics to PNG file
do ij=0,nmods-1
   ; ##########################################   
   ; Read in data
   ; ##########################################
   fileIN = dirIN+"composite.frequency."+persistence(ip2)+"hrs."+ptile(ip1)+"ptile."+modelID(ij)+".nc"
   fileID = addfile(fileIN,"r")
   lon                   = fileID->lon(5:,5:)
   lat                   = fileID->lat(5:,5:)
   precip_total_hist     = fileID->total_precip_hist(5:,5:)
   precip_event_hist     = fileID->event_precip_hist(5:,5:)
   precip_total_future   = fileID->total_precip_future(5:,5:)
   precip_event_future   = fileID->event_precip_future(5:,5:)
   
   ; ##########################################
   ; Display precipitation totals as mean monthly values
   ; ##########################################
   precip_total_hist   = precip_total_hist/(31.*6.)	
   precip_event_hist   = precip_event_hist/(31.*6.)
   precip_total_future = precip_total_future/(31.*6.)	
   precip_event_future = precip_event_future/(31.*6.)	

   ; ##########################################
   ; Compute non-event precipitation
   ; ##########################################
   precip_nonevent_hist   = precip_total_hist   - precip_event_hist
   precip_nonevent_future = precip_total_future - precip_event_future

   ; ##########################################
   ; Compute future-historical
   ; ##########################################
   dprecip_total   = precip_total_future    - precip_total_hist
   dprecip_event   = precip_event_future    - precip_event_hist
   dprecip_noevent = precip_nonevent_future - precip_nonevent_hist
   
   ; ##########################################
   ; Add attributes
   ; ##########################################
   dprecip_total@lat2d = lat
   dprecip_total@lon2d = lon
   dprecip_total@units = " "
   ;
   dprecip_event@lat2d = lat
   dprecip_event@lon2d = lon
   dprecip_event@units = " "
   ;
   dprecip_noevent@lat2d = lat
   dprecip_noevent@lon2d = lon
   dprecip_noevent@units = " "

   ; ##########################################
   ; Make plots
   ; ##########################################
   resD@tmYLLabelsOn = False
   resD@tmXBLabelsOn = False
   resD@tmYLLabelsOn = False
   resD@tmXBLabelsOn = True
   resD@gsnCenterString   = ""
   if (ij .eq. 0) then
      resD@tmYLLabelsOn = True
      resD@tmYLLabelsOn = True
   end if
   resD@tiMainString = modelID(ij)
   plot(ij)         = gsn_csm_contour_map(wks,dprecip_total,resD)
   resD@tiMainString = ""
   plot(ij+nmods)   = gsn_csm_contour_map(wks,dprecip_event,resD)
   resD@tmXBLabelsOn = True
   plot(ij+nmods*2) = gsn_csm_contour_map(wks,dprecip_noevent,resD)
 
end do

; #########################################################################
; Make panel plot (1 upper row, 2 lower rows)
; #########################################################################
; Upper panel configuration
res_U                   = True
res_U@gsnFrame          = False
res_U@gsnPanelTop       = 0.9999
res_U@gsnPanelBottom    = 0.60
res_U@gsnPanelLabelBar  = False
res_U@lbTitleString     = "(mm)"
res_U@lbTitleFontHeightF= .01 
res_U@pmLabelBarOrthogonalPosF = -0.015              ; move whole thing down a tad

; Lower panel configuration
res_L                   = True
res_L@gsnFrame          = False
;res_L@gsnPanelTop       = 0.60
;res_L@gsnPanelBottom    = 0.1
res_L@gsnPanelLabelBar  = True
res_L@lbTitleString     = "Precipitation (mm)"
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
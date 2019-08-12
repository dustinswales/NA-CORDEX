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
doIVT    = 01

; Which RCM?
modelID = (/"erain","hadgem","gfdl","mpi"/)
nmods   = dimsizes(modelID)

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

if (doPrecip .eq. 1) then
   varName  = "precip"	
   varRange   = (/-10.,10./)                
   ctable   = "NCV_bright"
   units    = "mm"
end if
if (doIVT .eq. 1) then
   varName  = "ivt"
   varRange = (/-200.,200./)                       
   ctable   = "NCV_bright"
   units    = "Kg/m/s"
end if

; Variable spacing.
dVardPrecip   = ((varRange(1)-varRange(0))/nLevels)

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
res@mpMinLonF               = 205.                     ;
res@mpMaxLonF               = 270.                     ;
res@mpCenterLonF            = 235.                     ;
res@cnMissingValFillPattern = 17                       ; set the missing value fill pattern to 17
res@tmXBLabelFontHeightF    = 0.015                    ; resize tick labels
res@tmYLLabelFontHeightF    = 0.015		       ;

; Load variable specific configuration
resD                        = res
resD@cnFillPalette          = cmapD                    ; divergent colormap
resD@cnMinLevelValF         = varRange(0)       ; Lower range
resD@cnMaxLevelValF         = varRange(1)       ; Uppermost level
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
do ij=1,nmods-1
   ; ##########################################   
   ; Read in data
   ; ##########################################
   fileIN = dirIN+"composite.frequency."+persistence(ip2)+"hrs."+ptile(ip1)+"ptile."+modelID(ij)+".nc"
   fileID = addfile(fileIN,"r")
   lon                   = fileID->lon(5:,5:)
   lat                   = fileID->lat(5:,5:)
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
      
      ; Display precipitation totals as daily values
      precip_total_hist      = precip_total_hist*8.
      precip_event_hist      = precip_event_hist*8.
      precip_nonevent_hist   = precip_nonevent_hist*8.
      precip_total_future    = precip_total_future*8.
      precip_event_future    = precip_event_future*8.
      precip_nonevent_future = precip_nonevent_future*8.
   end if
   if (doIVT .eq. 1) then
      ivt_total_hist       = fileID->total_ivt_hist(5::,5::)
      ivt_event_hist       = fileID->event_ivt_hist(5::,5::)
      ivt_nonevent_hist    = fileID->nevent_ivt_hist(5::,5::)
      ivt_total_future     = fileID->total_ivt_future(5::,5::)
      ivt_event_future     = fileID->event_ivt_future(5::,5::)
      ivt_nonevent_future  = fileID->nevent_ivt_future(5::,5::)
   end if

   ; ##########################################
   ; Compute future-historical
   ; ##########################################
   if (doPrecip .eq. 1) then
      dprecip_total   = precip_total_future    - precip_total_hist
      dprecip_event   = precip_event_future    - precip_event_hist
      dprecip_noevent = precip_nonevent_future - precip_nonevent_hist
      if (ij .eq. 0) then
      	 dprecip_total   = 0
      	 dprecip_event   = 0
      	 dprecip_noevent = 0
       end if
   
       ; Add attributes
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

       ;
       dVar_total   = dprecip_total
       dVar_event   = dprecip_event
       dVar_noevent = dprecip_noevent
   end if
   ;
   if (doIVT .eq. 1) then
      divt_total   = ivt_total_future    - ivt_total_hist
      divt_event   = ivt_event_future    - ivt_event_hist
      divt_noevent = ivt_nonevent_future - ivt_nonevent_hist
      if (ij .eq. 0) then
      	 divt_total   = 0
      	 divt_event   = 0
      	 divt_noevent = 0
       end if
   
       ; Add attributes
       divt_total@lat2d = lat
       divt_total@lon2d = lon
       divt_total@units = " "
       ;
       divt_event@lat2d = lat
       divt_event@lon2d = lon
       divt_event@units = " "
       ;
       divt_noevent@lat2d = lat
       divt_noevent@lon2d = lon
       divt_noevent@units = " "

       ;
       dVar_total   = divt_total
       dVar_event   = divt_event
       dVar_noevent = divt_noevent
   end if

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
   resD@tiMainString = "Cummulative ("+modelID(ij)+")"
   plot(ij)         = gsn_csm_contour_map(wks,dVar_total,resD)
   resD@tiMainString = ""
   resD@tiMainString = "Event-Days ("+modelID(ij)+")"
   plot(ij+nmods)   = gsn_csm_contour_map(wks,dVar_event,resD)
   resD@tmXBLabelsOn = True
   resD@tiMainString = "non-Event Days ("+modelID(ij)+")"
   plot(ij+nmods*2) = gsn_csm_contour_map(wks,dVar_noevent,resD)
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
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
doIVT    = 0
doZ500   = 0
doSWE    = 1

; Which RCM?
;modelID = (/"erain","hadgem","gfdl","mpi"/)
modelID = (/"erain","hadgem"/)
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
   varRange = (/-50.,50./)                
   units    = "mm"
   ctable   = "NCV_blue_red"
end if
if (doIVT .eq. 1) then
   varName  = "ivt"
   varRange = (/10.,100./)                       
   units    = "Kg/m/s"
;   ctable   = "NCV_bright"
   ctable   = "CBR_wet"
end if
;
if (doZ500 .eq. 1) then
   varName  = "z500"
   varRange = (/50,200/)                       
   units    = "m"
   ctable   = "NCV_bright"
end if
;
if (doSWE .eq. 1) then
   varName  = "SWE"
   varRange = (/-15,15/)                       
   ctable   = "NCV_blue_red"
   units    = "cm"
end if

; Variable spacing.
dVardPrecip   = ((varRange(1)-varRange(0))/nLevels)

; Load colortables for plots
cmap      = read_colormap_file(ctable)
;cmap(0,:) = namedcolor2rgba("white")
ncmap     = dimsizes(cmap(:,0))                      ; # colors
if (doPrecip .eq. 1) then
   cmap      = cmap(::-1,:)
end if

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
resD                        = res
resD@cnFillPalette          = cmap                     ; divergent colormap
resD@cnMinLevelValF         = varRange(0)              ; Lower range
resD@cnMaxLevelValF         = varRange(1)              ; Uppermost level
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
      
      ; Display precipitation totals as monthly mean values
      precip_total_hist  = precip_total_hist*nTime_hist(0)	/ (8.*30.) ; Express as mean monthly totals
      precip_total_hist@lat2d = lat
      precip_total_hist@lon2d = lon
      precip_total_hist@units = " "
      ;
      precip_event_hist  = precip_event_hist*nTime_event_hist(0) / (8.*30.) ; Express as mean monthly totals
      precip_event_hist@lat2d = lat
      precip_event_hist@lon2d = lon
      precip_event_hist@units = " "
      ;
      precip_nonevent_hist  = precip_nonevent_hist*(nTime_hist(0) - nTime_event_hist(0)) / (8.*30.) ; Express as mean monthly totals
      precip_nonevent_hist@lat2d = lat
      precip_nonevent_hist@lon2d = lon
      precip_nonevent_hist@units = " "
      ;
      precip_total_future = precip_total_future*nTime_future(0) / (8.*30.) ; Express as mean monthly totals
      precip_total_future@lat2d = lat
      precip_total_future@lon2d = lon
      precip_total_future@units = " "
      ;
      precip_event_future = precip_event_future*nTime_event_future(0) / (8.*30.) ; Express as mean monthly totals
      precip_event_future@lat2d = lat
      precip_event_future@lon2d = lon
      precip_event_future@units = " "
      ;
      precip_nonevent_future  = precip_nonevent_future*(nTime_future(0) - nTime_event_future(0)) / (8.*30.) ; Express as mean monthly totals
      precip_nonevent_future@lat2d = lat
      precip_nonevent_future@lon2d = lon
      precip_nonevent_future@units = " "
   end if
   if (doIVT .eq. 1) then
      ivt_total_hist       = fileID->total_ivt_hist(5::,5::)
      ivt_event_hist       = fileID->event_ivt_hist(5::,5::)
      ivt_nonevent_hist    = fileID->nevent_ivt_hist(5::,5::)
      ivt_total_future     = fileID->total_ivt_future(5::,5::)
      ivt_event_future     = fileID->event_ivt_future(5::,5::)
      ivt_nonevent_future  = fileID->nevent_ivt_future(5::,5::)
   end if
   if (doZ500 .eq. 1) then
      z500_total_hist       = fileID->total_z500_hist(5::,5::)
      z500_event_hist       = fileID->event_z500_hist(5::,5::)
      z500_nonevent_hist    = fileID->nevent_z500_hist(5::,5::)
      z500_total_future     = fileID->total_z500_future(5::,5::)
      z500_event_future     = fileID->event_z500_future(5::,5::)
      z500_nonevent_future  = fileID->nevent_z500_future(5::,5::)
   end if
   if (doSWE .eq. 1) then
      swe_total_hist       = fileID->total_swe_hist(5::,5::)
      swe_event_hist       = fileID->event_swe_hist(5::,5::)
      swe_nonevent_hist    = fileID->nevent_swe_hist(5::,5::)
      swe_total_future     = fileID->total_swe_future(5::,5::)
      swe_event_future     = fileID->event_swe_future(5::,5::)
      swe_nonevent_future  = fileID->nevent_swe_future(5::,5::)

      ; Express in cm
      swe_total_hist  = swe_total_hist/10. ; Express in cm.	
      swe_total_hist@lat2d = lat
      swe_total_hist@lon2d = lon
      swe_total_hist@units = " "
      ;
      swe_event_hist  = swe_event_hist/10. ; Express in cm.
      swe_event_hist@lat2d = lat
      swe_event_hist@lon2d = lon
      swe_event_hist@units = " "
      ;
      swe_total_future = swe_total_future/10. ; Express in cm.
      swe_total_future@lat2d = lat
      swe_total_future@lon2d = lon
      swe_total_future@units = " "
      ;
      swe_event_future = swe_event_future/10. ; Express in cm.
      swe_event_future@lat2d = lat
      swe_event_future@lon2d = lon
      swe_event_future@units = " "
      ;
      swe_nonevent_hist       = swe_nonevent_hist/10. ; Express in cm.
      swe_nonevent_hist@lat2d = lat
      swe_nonevent_hist@lon2d = lon
      swe_nonevent_hist@units = " "
      ;
      swe_nonevent_future       = swe_nonevent_future/10. ; Express in cm.
      swe_nonevent_future@lat2d = lat
      swe_nonevent_future@lon2d = lon
      swe_nonevent_future@units = " "
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
   ;
   if (doZ500 .eq. 1) then
      dz500_total   = z500_total_future    - z500_total_hist
      dz500_event   = z500_event_future    - z500_event_hist
      dz500_noevent = z500_nonevent_future - z500_nonevent_hist
      if (ij .eq. 0) then
      	 dz500_total   = 0
      	 dz500_event   = 0
      	 dz500_noevent = 0
       end if
   
       ; Add attributes
       dz500_total@lat2d = lat
       dz500_total@lon2d = lon
       dz500_total@units = " "
       ;
       dz500_event@lat2d = lat
       dz500_event@lon2d = lon
       dz500_event@units = " "
       ;
       dz500_noevent@lat2d = lat
       dz500_noevent@lon2d = lon
       dz500_noevent@units = " "
       ;
       dVar_total   = dz500_total
       dVar_event   = dz500_event
       dVar_noevent = dz500_noevent
   end if
   ;
   if (doSWE .eq. 1) then
      dswe_total   = swe_total_future    - swe_total_hist
      dswe_event   = swe_event_future    - swe_event_hist
      dswe_noevent = swe_nonevent_future - swe_nonevent_hist
      if (ij .eq. 0) then
      	 dswe_total   = 0
      	 dswe_event   = 0
      	 dswe_noevent = 0
       end if
   
       ; Add attributes
       dswe_total@lat2d = lat
       dswe_total@lon2d = lon
       dswe_total@units = " "
       ;
       dswe_event@lat2d = lat
       dswe_event@lon2d = lon
       dswe_event@units = " "
       ;
       dswe_noevent@lat2d = lat
       dswe_noevent@lon2d = lon
       dswe_noevent@units = " "
       ;
       dVar_total   = dswe_total
       dVar_event   = dswe_event
       dVar_noevent = dswe_noevent
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
   res_L@lbTitleString     = "Precipitation (mm/month)"
end if
if (doIVT .eq. 1) then
   res_L@lbTitleString     = "IVT (kg/m/s)"
end if
if (doZ500 .eq. 1) then
   res_L@lbTitleString     = "Geopotential height @ 500hPa (m)"
end if
if (doSWE .eq. 1) then
   res_L@lbTitleString     = "Snow Water Equivalent (cm)"
end if
res_L@lbTitleFontHeightF= .01
res_L@pmLabelBarOrthogonalPosF = -0.0075              ; move whole thing down a tad

; Panel upper and lower
gsn_panel(wks,plot,(/3,nmods/),res_L)
frame(wks)

end do	; Loop over duration
end do 	; Loop over persistence
; #########################################################################
; END PROGRAM
; #########################################################################
end
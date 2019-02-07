; #########################################################################
; #########################################################################

; #########################################################################
; BEGIN PROGRAM
; #########################################################################
begin
;
; What composites to look at?
;
; WRF horizontal resolution (50 or 25)
resolution  = 50
; Event strength threshold (percentile)
ptile       = "97.00"
; Event persistence threshold (hours)
persistence = 24
; Model ID
modelID = "hadgem"
configID = "Take1"

; Input directory
dirIN = "/data/dswales/NA-CORDEX/ARdet/dryingRatios/"

; Input file
fileIN = dirIN+"dryingRatios."+configID+"."+tostring(resolution)+"km."+ptile+"ptile."+tostring(persistence)+"hrs."+modelID+".nc"

; Output file
dirOUT = "/home/dswales/Projects/NA-CORDEX/scripts/ARdet/dryingRatioAnalysis/"
fileOUT = "jointHistograms.dryingRatios."+configID+"."+tostring(resolution)+"km."+ptile+"ptile."+tostring(persistence)+"hrs."+modelID

; #########################################################################
; Read in data
; #########################################################################
fileID    = addfile(fileIN,"r")
xBins     = fileID->xpdf_center
yBins     = fileID->ypdf_center
xEdge     = fileID->xpdf_edge
yEdge     = fileID->ypdf_edge
pdfFuture = fileID->pdfF
pdfHist   = fileID->pdfH
nxBins    = dimsizes(xBins)
nyBins    = dimsizes(yBins)

; Normalize PDFs, express in percent
pdfFuture = 100.*pdfFuture/sum(pdfFuture)
pdfHist   = 100.*pdfHist/sum(pdfHist)

; Compute relative difference (not used kinda ugly...)
dVarRel = 100.*(1-pdfHist/where(pdfFuture .eq. 0,999.,pdfFuture))
dVarRel = where(pdfFuture .eq. 0., 0., dVarRel)
eps = 0.001
dVarRel = where(pdfHist .lt. eps .and. pdfFuture .lt. eps, 0., dVarRel)

; Plot range
nLevels  = 20.0
varRange =(/0,3./)
dvarRange =(/-2,2/)
dvarRangeRelative = (/-100.,100./)
dVar     = ((varRange(1)-varRange(0))/nLevels)
ddVar    = ((dvarRange(1)-dvarRange(0))/nLevels)
ddVarRelative   = ((dvarRangeRelative(1)-dvarRangeRelative(0))/nLevels)

; #########################################################################
; Plot setup
; #########################################################################
plot  = new(2,graphic)
plotD = new(1,graphic)
wks_type = "png"
wks  = gsn_open_wks(wks_type,dirOUT+fileOUT); send graphics to PNG file

; Read in colortables
cmap = read_colormap_file("CBR_wet")
;do ij=0,1
;   cmap(ij,:) = namedcolor2rgba("white")
;end do

cmapD = read_colormap_file("NCV_blu_red")
;cmapD = cmapD(::-1,:) ; reverse the color map
do ij=120,135
   cmapD(ij,:) = 1
end do

; Common plot configuration
res                      = True
res@gsnDraw              = False          ; don't draw
res@gsnFrame             = False          ; don't advance frame
res@cnInfoLabelOn        = False          ; turn off cn info label   
res@cnFillOn             = True           ; Turn on contour fill
res@lbLabelBarOn         = False          ; turn off individual cb's
res@cnFillMode           = "RasterFill"   ; Turn on raster fill
res@cnLineLabelsOn       = False          ; Turn off line labels
res@cnLinesOn            = False          ; Turn off contour lines
res@cnLevelSelectionMode = "ManualLevels" ;
res@tmYROn               = False          ; Turn off Y tick marks on right side of plot
res@tmXTOn               = False          ; Turn off X tick marks on top of plot
res@tmXBMode             = "Explicit"     ;
res@tmXBMinorOn          = False          ;
res@tmXBValues           = xEdge          ;
res@tmXBLabels           = xEdge          ;
res@sfXArray             = res@tmXBValues ;
res@tmYLMode             = "Explicit"     ;
res@tmYLMinorOn          = False          ;
res@tmYLValues           = yEdge          ;
res@tmYLLabels           = yEdge          ;
res@sfYArray             = res@tmYLValues ;
res@tiXAxisString        = ""
res@tiYAxisString        = "" 
res@tmXBLabelFontHeightF        = 0.0075         ; Resize tick labels
res@tmYLLabelFontHeightF        = 0.0075	 ; Resize tick labels
res@tmXBLabelAngleF             = 45		 ; Tilt x-axis labels
res@tiXAxisFontHeightF          = 0.015
res@tiYAxisFontHeightF          = 0.015

; Configuration for absolute plots
resA                      = res
resA@cnFillPalette        = cmap           ;
resA@cnMinLevelValF       = varRange(0)    ; Lower range
resA@cnMaxLevelValF       = varRange(1)    ; Uppermost level
resA@cnLevelSpacingF      = dVar           ;

; Configuration for difference (absolute) plots
resDA                      = res
resDA@cnFillPalette        = cmapD           ;
resDA@cnMinLevelValF       = dvarRange(0)   ; Lower range
resDA@cnMaxLevelValF       = dvarRange(1)   ; Uppermost level
resDA@cnLevelSpacingF      = ddVar
;resDA@lbLabelBarOn         = True          ; turn on individual cb's
;resDA@lbTitleString        = "Change in Frequency (Absolute) (%)"
;resDA@lbTitleFontHeightF   = .01

; Configuration for difference (relative) plots
resDR                      = res
resDR@cnFillPalette        = cmapD           ;
resDR@cnMinLevelValF       = dvarRangeRelative(0)   ; Lower range
resDR@cnMaxLevelValF       = dvarRangeRelative(1)   ; Uppermost level
resDR@cnLevelSpacingF      = ddVarRelative        
;resDR@lbLabelBarOn         = True          ; turn on individual cb's
;resDR@lbTitleString        = "Change in Frequency (Relative) (%)"
;resDR@lbTitleFontHeightF   = .01

; #########################################################################
; Make panel plot (2 upper, 1 lower)
; #########################################################################
; Upper panel configuration
res_U                   = True
res_U@gsnFrame          = False
res_U@gsnPanelTop       = 0.90
res_U@gsnPanelBottom    = 0.50
res_U@gsnPanelLabelBar  = True
res_U@lbTitleString     = "Frequency of Occurence(%)"
res_U@lbTitleFontHeightF= .01 

; Lower panel configuration
res_L                   = True
res_L@gsnFrame          = False
res_L@gsnPanelTop       = 0.48
res_L@gsnPanelLabelBar  = True
res_L@lbTitleString     = "Change in Frequency (Absolute) (%)"
res_L@lbTitleFontHeightF= .01

; Make plots
; 1) Historical (Upper left)
resA@tiXAxisString    = "Drying-ratio (1)"
resA@tiYAxisString    = "Precipitation along transect (mm)"
plot(0) = gsn_csm_contour(wks,transpose(pdfHist),resA)
; 2) Future (Upper right)
plot(1) = gsn_csm_contour(wks,transpose(pdfFuture),resA)
; 3) Difference (absolute) (Lower left)
resDA@tiXAxisString    = "Drying-ratio (1)"
resDA@tiYAxisString    = "Precipitation along transect (mm)"
plotD(0) = gsn_csm_contour(wks,transpose(pdfFuture-pdfHist),resDA)
;; 4) Difference (relative) (Lower right)
;resDR@tiXAxisString    = "Drying-ratio (1)"
;resDR@tiYAxisString    = "Precipitation along transect (mm)"
;plotD(1) = gsn_csm_contour(wks,transpose(dVarRel),resDR)

; Panel upper and lower
gsn_panel(wks,plot,(/1,2/),res_U)
gsn_panel(wks,plotD,(/1,1/),res_L)
frame(wks)

; #########################################################################
; END PROGRAM
; #########################################################################
end
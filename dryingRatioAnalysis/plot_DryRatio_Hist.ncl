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
modelID = "mpi"
configID = "Take4"

; Input directory
dirIN = "/data/dswales/NA-CORDEX/ARdet/dryingRatios/"

; Input file
fileIN = dirIN+"dryingRatios."+configID+"."+tostring(resolution)+"km."+ptile+"ptile."+tostring(persistence)+"hrs."+modelID+".nc"

; Output file
dirOUT = "/home/dswales/Projects/NA-CORDEX/scripts/ARdet/dryingRatioAnalysis/"
fileOUT = "histograms.dryingRatios."+configID+"."+tostring(resolution)+"km."+ptile+"ptile."+tostring(persistence)+"hrs."+modelID

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

; Collapse joint-PDFs into 1D
pdfF = dim_sum_n(pdfFuture,1)
pdfH = dim_sum_n(pdfHist,1)

; #########################################################################
; Plot
; #########################################################################
wks_type = "png"
wks  = gsn_open_wks(wks_type,dirOUT+fileOUT); send graphics to PNG file

; Plot configuration
res                = True
res@xyLineThicknesses = (/  2.0,   2.0, 2.0, 1.0/)          
res@xyLineColors      = (/"blue","red", "black"/)
res@xyDashPatterns    = (/0.0,0.0,0.0,1.0/) 
res@tmYROn            = False          ; Turn off Y tick marks on right side of plot
res@tmXTOn            = False          ; Turn off X tick marks on top of plot
res@tiXAxisString     = "Drying-Ratio (1)"
res@tiYAxisString     = "Frequency (%)"
res@trYMinF           = -3         
res@trYMaxF           = 15         

; Put data to be plotted into array (nLine,datalength)
data = new((/4,nxBins/),float)
data(0,:) = pdfH
data(1,:) = pdfF
data(2,:) = pdfF-pdfH
aline =  new(nxBins,float)
aline(:)=0.0
data(3,:) = aline

; make plot
plot  = gsn_csm_xy (wks,xBins,data,res) ; create plot

; #########################################################################
; END PROGRAM
; #########################################################################
end
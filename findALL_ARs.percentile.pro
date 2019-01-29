; #########################################################################
;
; This code is a wrapper to call find_ARs.percentiles.
;
; #########################################################################
@find_ARs.percentiles.pro

; Output location
dirOUT = '/data/dswales/NA-CORDEX/ARdet/events/'

; Intensity percentiles to examine
ptiles = [80+indgen(19),99+findgen(9)/10.,99.9+findgen(10)/100.]
nptiles = n_elements(ptiles)

; Temporal range
time_range=[1980,2010]
;time_range=[2070,2099]
if (time_range(0) lt 2010 and time_range(1) le 2010) then period='historical' else period='future'

; Persistence threshold
persistence=36

; Models
models = ['gfdl','mpi','hadgem']
if (strcmp(period,'historical')) then models=[models,'erain']
nmods = n_elements(models)

; Create output text file
fileOUT = 'data/event_frequency.'+string(persistence,format='(i2)')+'hr.'+period+'.txt'
openw,101,fileOUT
printf,101,'  Model        Persistence      Percentile       nEvents        nEventDays   EventFrequency'
;
for ij=0,nmods-1 do begin
   for ik=0,nptiles-1 do begin
      find_ARs,models(ij),persistence,ptiles(ik),time_range(0),time_range(1),dirOUT
   end
end
close,101


; END PROGRAM
end

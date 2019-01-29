; #######################################################################
;
; The purpose of this program is to look at drying ratios for detected
; events.
;
; #######################################################################

res         = 50
ptile       = 97
persistence = 24

; Data directory
dir = '/home/dswales/Projects/NA-CORDEX/scripts/ARdet/data/composites/' 
files = file_search(dir,'*raw*'+string(res,format='(i2)')+'km.'+string(ptile,format='(f5.2)')+$
                    'ptile.'+string(persistence,format='(i2)')+'hrs*')
nfiles = n_elements(files)


; #######################################################################
; END PROGRAM
; #######################################################################
end

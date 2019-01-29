; #########################################################################
;
; This code is a wrapper to call compositeEvents.pro.
;
; #########################################################################
@compositeEvents.pro

dirOUT = '/data/dswales/NA-CORDEX/ARdet/composites/'
dirIN  = '/data/dswales/NA-CORDEX/ARdet/events/'

; RCMs
models = ['erain','mpi','gfdl','hadgem']
models = ['gfdl','hadgem']
nmods  = n_elements(models)

; IVT threshold (percentile)
ivtThresh = [96.00,97.00,98.00,99.00]
nithresh  = n_elements(ivtThresh)

; Persistence threshold
pthresh  = [24,36]
npthresh = n_elements(pthresh)

for ij=0,nmods-1 do begin
   for ik=0,nIthresh-1 do begin
      for il=0,npthresh-1 do begin
         print,'########################################################'
         print,models(ij),ivtThresh(ik),pThresh(il)
         compositeEvents,dirIN,dirOUT,models(ij),ivtThresh(ik),pThresh(il),0
         if (ij ne 0) then compositeEvents,dirOUT,models(ij),ivtThresh(ik),pThresh(il),1
      end
   end
end


; END PROGRAM
end

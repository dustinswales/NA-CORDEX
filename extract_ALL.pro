@extract_climo_fields

; WRF RCM
modelID = ['erain','gfdl','hadgem','mpi']
nmod    = n_elements(modelID)

for ij=0,nmod-1 do begin
   extract_climo_fields,modelID(ij),1980,2010
   if (ij gt 0) then extract_climo_fields,modelID(ij),2070,2100
end



; END PROGRAM
end

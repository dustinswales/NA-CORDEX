; ##################################################################### 
;
; The purpose of this program is to extract WRF IVT data at grid
; points defined in build_costalSlab.pro
;
; Name    Type    Size    Dscription
;
; Inputs
; res     string  [1]     String ID for WRF resolution ('50' or '25')
; modelID string  [1]     Shortname for WRF-GCM (e.g. 'gfdl')
; dirOUT  string  [1]     Output directory location.
;
; Output is written to netCDF and storred in dirOUT
;
; #####################################################################
pro pullOutSlab,res,modelID,dirOUT
  ;res = '50'
  ;modelID = 'gfdl'

  ; Output location
  ;dirOUT = '/data/dswales/NA-CORDEX/ARdet/slabs/'

  ; Input data location
  expID = 'WRF_'+modelID+'_'+res+'km'
  dir   ='/Projects/HydroMet/dswales/NA-CORDEX/'+expID+'/'

  ; What points to pull out?
  openr,101,'data/coastalSlabXY.coast1.'+res+'.txt'
  readf,101,npts
  xi = intarr(npts)
  yi = intarr(npts)
  readf,101,xi
  readf,101,yi
  close,101
  
  ; Get file names
  files  = file_search(dir,'wrfout_*.nc')
  nfiles = n_elements(files)
 
  ; Loop over all files
  for ij=0,nfiles-1 do begin
     print,ij+1,' of ',nfiles
     fileID=ncdf_open(files(ij))
     ncdf_varget,fileID,ncdf_varid(fileID,'Year'),year1
     ncdf_varget,fileID,ncdf_varid(fileID,'Month'),month1
     ncdf_varget,fileID,ncdf_varid(fileID,'Day'),day1
     ncdf_varget,fileID,ncdf_varid(fileID,'Hour'),hour1
     ntime = n_elements(day1)
     ivtUtemp = fltarr(npts,ntime)
     ivtVtemp = fltarr(npts,ntime)
     Ttemp    = fltarr(npts,ntime,10)
     Qtemp    = fltarr(npts,ntime,10)
     Utemp    = fltarr(npts,ntime,10)
     Vtemp    = fltarr(npts,ntime,10)
     for ik=0,npts-1 do begin
        ncdf_varget,fileID,ncdf_varid(fileID,'IVTU'),ivtUtemp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        ncdf_varget,fileID,ncdf_varid(fileID,'IVTV'),ivtVtemp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        ivtUtemp(ik,*) = ivtUtemp1
        ivtVtemp(ik,*) = ivtVtemp1
        ;
        ncdf_varget,fileID,ncdf_varid(fileID,'T1000'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Ttemp(ik,*,0) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'T950'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Ttemp(ik,*,1) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'T900'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Ttemp(ik,*,2) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'T850'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Ttemp(ik,*,3) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'T800'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Ttemp(ik,*,4) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'T750'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Ttemp(ik,*,5) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'T700'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Ttemp(ik,*,6) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'T600'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Ttemp(ik,*,7) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'T500'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Ttemp(ik,*,8) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'T250'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Ttemp(ik,*,9) = temp1
        ;
        ncdf_varget,fileID,ncdf_varid(fileID,'Q1000'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Qtemp(ik,*,0) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'Q950'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Qtemp(ik,*,1) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'Q900'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Qtemp(ik,*,2) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'Q850'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Qtemp(ik,*,3) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'Q800'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Qtemp(ik,*,4) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'Q750'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Qtemp(ik,*,5) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'Q700'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Qtemp(ik,*,6) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'Q600'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Qtemp(ik,*,7) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'Q500'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Qtemp(ik,*,8) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'Q250'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Qtemp(ik,*,9) = temp1
        ;
        ncdf_varget,fileID,ncdf_varid(fileID,'U1000'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Utemp(ik,*,0) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'U950'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Utemp(ik,*,1) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'U900'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Utemp(ik,*,2) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'U850'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Utemp(ik,*,3) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'U800'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Utemp(ik,*,4) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'U750'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Utemp(ik,*,5) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'U700'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Utemp(ik,*,6) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'U600'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Utemp(ik,*,7) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'U500'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Utemp(ik,*,8) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'U250'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Utemp(ik,*,9) = temp1
        
        ;
        ncdf_varget,fileID,ncdf_varid(fileID,'V1000'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Vtemp(ik,*,0) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'V950'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Vtemp(ik,*,1) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'V900'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Vtemp(ik,*,2) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'V850'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Vtemp(ik,*,3) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'V800'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Vtemp(ik,*,4) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'V750'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Vtemp(ik,*,5) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'V700'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Vtemp(ik,*,6) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'V600'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Vtemp(ik,*,7) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'V500'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Vtemp(ik,*,8) = temp1
        ncdf_varget,fileID,ncdf_varid(fileID,'V250'),temp1,offset=[xi(ik),yi(ik),0],count=[1,1,ntime]
        Vtemp(ik,*,9) = temp1
     end
     ncdf_close,fileID
     
     ; Guard against bad IVT data
     testU = where(abs(ivtUtemp) gt 10000)
     if (testU(0) ne -1) then ivtUtemp(testU)=0.0
     testV = where(abs(ivtVtemp) gt 10000)
     if (testV(0) ne -1) then ivtVtemp(testV)=0.0

     if (ij eq 0) then begin
        ivtU  = ivtUtemp
        ivtV  = ivtVtemp
        year  = year1
        month = month1
        day   = day1
        hour  = hour1
        U     = Utemp
        V     = Vtemp
        Q     = Qtemp
        T     = Ttemp
     endif
     if (ij ne 0) then begin
        ivtU  = [[ivtU],[ivtUtemp]]
        ivtV  = [[ivtV],[ivtVtemp]]
        U     = [[U],[Utemp]]
        V     = [[V],[Vtemp]]
        Q     = [[Q],[Qtemp]]
        T     = [[T],[Ttemp]]        
        year  = [year,year1]
        month = [month,month1]
        day   = [day,day1]
        hour  = [hour,hour1]
     endif
  end
  
  fileID  = ncdf_create(dirOUT+'coastalSlabXY.coast1.'+res+'.'+modelID+'.nc',/clobber)
  dimID1  = ncdf_dimdef(fileID,'nPts',npts)
  dimID2  = ncdf_dimdef(fileID,'nTime',n_elements(day))
  dimID3  = ncdf_dimdef(fileID,'lev', 10)
  varID1  = ncdf_vardef(fileID,'year', [       dimID2       ],/long)
  varID2  = ncdf_vardef(fileID,'month',[       dimID2       ],/long)
  varID3  = ncdf_vardef(fileID,'day',  [       dimID2       ],/long)
  varID4  = ncdf_vardef(fileID,'hour', [       dimID2       ],/long)
  varID11 = ncdf_vardef(fileID,'lev',  [              dimID3],/long)
  varID5  = ncdf_vardef(fileID,'ivtU', [dimID1,dimID2       ],/float)
  varID6  = ncdf_vardef(fileID,'ivtV', [dimID1,dimID2       ],/float)
  varID7  = ncdf_vardef(fileID,'U',    [dimID1,dimID2,dimID3],/float)
  varID8  = ncdf_vardef(fileID,'V',    [dimID1,dimID2,dimID3],/float)
  varID9  = ncdf_vardef(fileID,'T',    [dimID1,dimID2,dimID3],/float)
  varID10 = ncdf_vardef(fileID,'Q',    [dimID1,dimID2,dimID3],/float)
  ncdf_control,fileID,/endef
  ncdf_varput,fileID,varID1,year
  ncdf_varput,fileID,varID2,month
  ncdf_varput,fileID,varID3,day
  ncdf_varput,fileID,varID4,hour
  ncdf_varput,fileID,varID11,[1000,950,900,850,800,750,700,600,500,250]
  ncdf_varput,fileID,varID5,ivtU
  ncdf_varput,fileID,varID6,ivtV
  ncdf_varput,fileID,varID7,U
  ncdf_varput,fileID,varID8,V
  ncdf_varput,fileID,varID9,T
  ncdf_varput,fileID,varID10,Q
  ncdf_close,fileID


; END PROGRAM
end
  

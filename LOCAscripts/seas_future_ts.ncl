load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; seas_future_ts.ncl 
; EA Burakowski
; 2018-05-08
; 
; seas_future_ts.ncl calculates the annual average times series 
;	seasonal temperature (min & max) and precipitation (total) for
;	rcp4.5 and rcp8.5 (2006-2099), using CMIP5 LOCAlized Constructed Analogs (LOCA) 
;	daily downscaled data (16th degree resolution). 
;	Inputs include daily LOCA data for three variables:
;		(1) tasmin (surface minimum air temperature)
;		(2) tasmax (surface maximum air temperature)
;		(3) precipitation (daily total)
;	for 29 CMIP5 models (see below, mdl), 2006-2099
;
; Output includes:
;   .nc files for seasonal time series in each variable.
;   Abbreviations for the seasons are as follows:
; 	- DJF = winter
;	- MAM = spring
;	- JJA = summer
;	- SON = fall
;
; Use plotSEAS_future_ts.ncl and plotSEAS_individ.ncl to plot up the .nc files
;	

;======================================================================
;                         The main code
;======================================================================

begin

;---- LOCA constants

  dir	= "/net/nfs/merrimack/raid/Northeast_US_Downscaling_cmip5/"
  var	= (/"pr","tasmax","tasmin"/)
  mdl	= (/"ACCESS1-0","ACCESS1-3",\
	    "bcc-csm1-1-m","bcc-csm1-1",\
	    "CanESM2",\
	    "CCSM4","CESM1-BGC","CESM1-CAM5",\
	    "CMCC-CM","CMCC-CMS",\
	    "CNRM-CM5",\
	    "CSIRO-Mk3-6-0",\
            "FGOALS-g2",\
	    "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
	    "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES",\
	    "inmcm4",\
	    "IPSL-CM5A-LR","IPSL-CM5A-MR",\
	    "MIROC5","MIROC-ESM-CHEM","MIROC-ESM",\
	    "MPI-ESM-LR","MPI-ESM-MR",\
	    "MRI-CGCM3",\
	    "NorESM1-M"/)

  ;--- Annoying stuff to write model names as strings to netCDF4 
  model	= ispan(0,dimsizes(mdl)-1,1)
  model!0 = "model"
  modelname = new(dimsizes(mdl),"string")
  do i = 0,dimsizes(mdl)-1  
    modelname(i)	= str_join((/tostring(model(i)), mdl(i)/)," - ")
  end do

  modelname!0 = "model"

  ;--- Scenario
  scn	= (/"rcp45","rcp85"/)

  ;--- Time periods
  fyr   	= ispan(2006,2099,1)
  fyr!0 	= "year"
  fyr&year 	= fyr
  mon		= ispan(1,12,1)

;---- Output directories

  odir	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/"	


;---- Seasonal Total and Averages
;
;      Loop over models and years to calculate annual average
;	tasmax, tasmin, pr for each year in record
;
;----

;---- Load one future file to get the lat/lon dims

  tfil 	= addfile(dir+"pr_day_ACCESS1-0_rcp45_r1i1p1_20810101-20811231.16th.nc","r")
  tmp1	= tfil->$var(0)$
  tvar  = tmp1(lat|:,lon|:,time|:)
  dims	= dimsizes(tvar)
  lat	= tvar&lat
  lon	= tvar&lon

printVarSummary(tvar)

;--- Loop over var, mod, scen daily files to calculate annual average
  do ivar = 0,dimsizes(var)-1
    print("----  Working on "+var(ivar)+"  -----")
      do iscn = 0,dimsizes(scn)-1	
        print("---- Working on "+scn(iscn)+" ----")

	;--- Initialize "+scn(iscn)+" seasonal trends
	  fmam_ts	= new((/dims(0),dims(1),dimsizes(fyr),dimsizes(mdl)/),float,-9999)
	
	;--- Assign coordinates
	  fmam_ts!0	= "lat"
	  fmam_ts&lat	= lat
	  fmam_ts!1	= "lon"
	  fmam_ts&lon	= lon
	  fmam_ts!2	= "year"
	  fmam_ts&year	= fyr
	  fmam_ts!3	= "model"
	  fmam_ts&model = model
	
	;--- Initialize other seasons based on spring (fmam_ts)  
	  fjja_ts	= fmam_ts
	  fson_ts	= fmam_ts
	  fdjf_ts	= fmam_ts
	  
	;--- loop over models  
          do imdl = 0,dimsizes(mdl)-1
	  print("----  Working on "+mdl(imdl)+"  -----")

      	  ;--- If/else statement for Precip vs. Tmin/Tmax
            if ivar.eq.0

   ;-----------------------------
   ;
   ;   PRECIP CALCS
   ;
   ;-----------------------------

	    ;--- First, the easy seasons (MAM, JJA, SON)
            ;--- Loop over years and calculate season total precip
	      do iyr = 0,dimsizes(fyr)-1
	  
	      ;--- if it's first year of the analysis, skip winter (@_FillValue)
	      if iyr.eq.0
                fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_"+scn(iscn)+"_r*_"+fyr(iyr)+"0101*.nc")
	        a		= addfile(dir+fil,"r")
                v		= a->$var(ivar)$
	    
	      ;--- convert kg m-2 s-2 to mm/day
	        vmm1	= v
	        vmm1	= v*86400
	        vmm	= vmm1(lat|:, lon|:, time|:)
  
              ;--- identify months and seasons in file
	        time	= a->time
	        utc_date= cd_calendar(time,0)
	        mm	= tointeger(utc_date(:,1))  ; really only need month 
	        mam	= ind(mm.ge.3 .and. mm.le.5)
	        jja	= ind(mm.ge.6 .and. mm.le.8)
	        son	= ind(mm.ge.9 .and. mm.le.11)
	
	      ;--- calculate the seasonal total, rearrange dims
	        fmam_ts(:,:,iyr,imdl)	= dim_sum_n_Wrap(vmm(:,:,mam),2)
	        fjja_ts(:,:,iyr,imdl)	= dim_sum_n_Wrap(vmm(:,:,jja),2)
	        fson_ts(:,:,iyr,imdl)	= dim_sum_n_Wrap(vmm(:,:,son),2)
	        fdjf_ts(:,:,iyr,imdl)	= fdjf_ts@_FillValue

	      ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
              delete([/fil,a,v,vmm1,vmm,time,utc_date,mm,mam,jja,son/])
	  
	     else  ; if it's not the first year
	    
	    ;--- 
	        fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_"+scn(iscn)+"_r*_"+fyr(iyr)+"0101*.nc")
	        a		= addfile(dir+fil,"r")
                v		= a->$var(ivar)$
	    
	      ;--- convert kg m-2 s-2 to mm/day
	        vmm1	= v
	        vmm1	= v*86400
	    
	      ;--- identify months and seasons in file
	        time	= a->time
	        utc_date	= cd_calendar(time,0)
	        mm	= tointeger(utc_date(:,1))  ; really only need month 
	        mam	= ind(mm.ge.3 .and. mm.le.5)
	        jja	= ind(mm.ge.6 .and. mm.le.8)
	        son	= ind(mm.ge.9 .and. mm.le.11)
	        jf	= ind(mm.le.2) ; jan & feb only
	      
	      ;--- December of previous year
	        fil_dec         = systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_"+scn(iscn)+"_r*_"+fyr(iyr-1)+"0101*.nc")
                a_dec           = addfile(dir+fil_dec,"r")
                v_dec           = a_dec->$var(ivar)$
                vmm1_dec        = v_dec
                vmm1_dec        = v_dec*86400
                time_dec        = a_dec->time
                utc_dec         = cd_calendar(time_dec,0)
                mm_dec          = tointeger(utc_dec(:,1))
                dec             = ind(mm_dec.eq.12)

              ;--- combine december of previous year (iyr-1) with jan-feb of current year (iyr)
                djf1            = array_append_record(vmm1_dec(dec,:,:),vmm1(jf,:,:),0)

              ;--- now rearrange dims
                vmm             = vmm1(lat|:, lon|:, time|:)
                djf             = djf1(lat|:, lon|:, time|:)

              ;--- Calculate seasonal sums
                fmam_ts(:,:,iyr,imdl)   = dim_sum_n_Wrap(vmm(:,:,mam),2)
                fjja_ts(:,:,iyr,imdl)   = dim_sum_n_Wrap(vmm(:,:,jja),2)
                fson_ts(:,:,iyr,imdl)   = dim_sum_n_Wrap(vmm(:,:,son),2)
                fdjf_ts(:,:,iyr,imdl)   = dim_sum_n_Wrap(djf,2)
		
              ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
                delete([/fil,a,v,vmm1,vmm,time,utc_date,mm,mam,jja,son,jf/])
                delete([/fil_dec,a_dec,v_dec,vmm1_dec,time_dec,utc_dec,mm_dec,dec,djf,djf1/])

	    end if   ; first year loop	  
    
          end do   ; "+scn(iscn)+" years
	
	else   ; pr/temp
	
   ;-----------------------------
   ;
   ;   TEMPERATURE CALCS
   ;
   ;-----------------------------   
	    
          ;--- Loop over years and calculate season total precip
	    do iyr = 0,dimsizes(fyr)-1
	  
	      ;--- if it's first year of the analysis, skip winter (@_FillValue)
	      if iyr.eq.0
              fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_"+scn(iscn)+"_r*_"+fyr(iyr)+"0101*.nc")
	      a		= addfile(dir+fil,"r")
              v		= a->$var(ivar)$
	    
	    ;--- convert K to C
	      vC1	= v
	      vC1	= v-273.15
	      vC	= vC1(lat|:, lon|:, time|:)
  
            ;--- identify months and seasons in file
	      time	= a->time
	      utc_date	= cd_calendar(time,0)
	      mm	= tointeger(utc_date(:,1))  ; really only need month 
	      mam	= ind(mm.ge.3 .and. mm.le.5)
	      jja	= ind(mm.ge.6 .and. mm.le.8)
	      son	= ind(mm.ge.9 .and. mm.le.11)
	
	    ;--- calculate the seasonal total, rearrange dims
	      fmam_ts(:,:,iyr,imdl)	= dim_avg_n_Wrap(vC(:,:,mam),2)
	      fjja_ts(:,:,iyr,imdl)	= dim_avg_n_Wrap(vC(:,:,jja),2)
	      fson_ts(:,:,iyr,imdl)	= dim_avg_n_Wrap(vC(:,:,son),2)
	      fdjf_ts(:,:,iyr,imdl)	= fdjf_ts@_FillValue

	    ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
              delete([/fil,a,v,vC1,vC,time,utc_date,mm,mam,jja,son/])
	  
	  else  ; if it's not the first year
	    
	    ;--- 
	      fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_"+scn(iscn)+"_r*_"+fyr(iyr)+"0101*.nc")
	      a		= addfile(dir+fil,"r")
              v		= a->$var(ivar)$
	    
	    ;--- convert K to C
	      vC1	= v
	      vC1	= v-273.15
	    
	    ;--- identify months and seasons in file
	      time	= a->time
	      utc_date	= cd_calendar(time,0)
	      mm	= tointeger(utc_date(:,1))  ; really only need month 
	      mam	= ind(mm.ge.3 .and. mm.le.5)
	      jja	= ind(mm.ge.6 .and. mm.le.8)
	      son	= ind(mm.ge.9 .and. mm.le.11)
	      jf	= ind(mm.le.2) ; jan & feb only
	      
	    ;--- December of previous year
	      fil_dec   = systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_"+scn(iscn)+"_r*_"+fyr(iyr-1)+"0101*.nc")
              a_dec     = addfile(dir+fil_dec,"r")
              v_dec     = a_dec->$var(ivar)$
              v1_dec    = v_dec
              vC1_dec   = v1_dec-273.15
              time_dec  = a_dec->time
              utc_dec   = cd_calendar(time_dec,0)
              mm_dec    = tointeger(utc_dec(:,1))
              dec       = ind(mm_dec.eq.12)

            ;--- combine december of previous year (iyr-1) with jan-feb of current year (iyr)
              djf1             = array_append_record(vC1_dec(dec,:,:),vC1(jf,:,:),0)
printVarSummary(djf1)

            ;--- now rearrange dims
              vC             = vC1(lat|:, lon|:, time|:)
              djf            = djf1(lat|:, lon|:, time|:)
printVarSummary(djf)
            
	    ;--- Calculate seasonal average
              fmam_ts(:,:,iyr,imdl)   = dim_avg_n_Wrap(vC(:,:,mam),2)
              fjja_ts(:,:,iyr,imdl)   = dim_avg_n_Wrap(vC(:,:,jja),2)
              fson_ts(:,:,iyr,imdl)   = dim_avg_n_Wrap(vC(:,:,son),2)
              fdjf_ts(:,:,iyr,imdl)   = dim_avg_n(djf,2)

            ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
              delete([/fil,a,v,vC1,vC,time,utc_date,mm,mam,jja,son,jf/])
              delete([/fil_dec,a_dec,v_dec,v1_dec,vC1_dec/])
	      delete([/time_dec,utc_dec,mm_dec,dec,djf,djf1/])

	    end if   ; first year

          end do   ; "+scn(iscn)+" years
	  
        end if ; pr/temp 
	   
      end do ; models 

  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ;
  ;--- create .nc file of seasonal trends
  ;
  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ;--- assign coordinate vars for climatological averages
    
      ;-- Time series
      fmam_ts@short_name 	= scn(iscn)+": LOCA Climatological spring time series of "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fmam_ts@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) Spring "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+" time series"
      if ivar.eq.0
        fmam_ts@units		= "mm"
      else
        fmam_ts@units 		= "degrees C"
      end if

      fjja_ts@short_name 	= scn(iscn)+": LOCA Climatological summer time series of "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fjja_ts@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) Summer "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+" time series"
      if ivar.eq.0
        fjja_ts@units		= "mm"
      else
        fjja_ts@units 		= "degrees C"
      end if
      
      fson_ts@short_name 	= scn(iscn)+": LOCA Climatological fall time series of "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fson_ts@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) Fall "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+" time series"
      if ivar.eq.0
        fson_ts@units		= "mm"
      else
        fson_ts@units 		= "degrees C"
      end if
      
      fdjf_ts@short_name 	= scn(iscn)+": LOCA Climatological winter time series of "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fdjf_ts@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) Winter "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+" time series"
      if ivar.eq.0
        fdjf_ts@units		= "mm"
      else
        fdjf_ts@units 		= "degrees C"
      end if
      
      
    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")

    ;--- create file name
      fn 		= "SEAS_"+var(ivar)+"_"+scn(iscn)+"_"+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+"_timeseries"
print("Filename = "+fn)

    ;--- remove old file, if it exists
      system("/bin/rm -f "+odir+fn+".nc")
      
    ;--- create new file   
      ncdf		= addfile(odir+fn+".nc","c")	
      fAtt		= True
      fAtt@title	= fn
      fAtt@orig_fil	= var(ivar)+"_day_LOCA_CMIP5_models_"+scn(iscn)+"_"+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fAtt@Conventions	= "COARDS/CF-1.0"
      fAtt@creation_date= systemfunc("date")
      fAtt@author	= "Elizabeth Burakowski (elizabeth.burakowski@unh.edu)"
      
      fAtt@description 	= scn(iscn)+": Seasonal "+var(ivar)+" time series calculated from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      
    ;--- file attributes, size
      fileattdef(ncdf,fAtt)
      dimNames		= (/"lat", "lon", "year", "model"/)
      dimSizes		= (/ dims(0), dims(1), dimsizes(fyr), dimsizes(model) /)
      dimUnlim          = (/ False, False, False, False /)
      filedimdef(ncdf,dimNames,dimSizes,dimUnlim)

    ;--- Define file variables
      filevardef(ncdf,"model",typeof(model),getvardims(model))
      filevardef(ncdf,"lat",typeof(lat),getvardims(lat))
      filevardef(ncdf,"lon",typeof(lon),getvardims(lon))
      filevardef(ncdf,"MAMTIMESERIES",typeof(fmam_ts),getvardims(fmam_ts))
      filevardef(ncdf,"JJATIMESERIES",typeof(fjja_ts),getvardims(fjja_ts))
      filevardef(ncdf,"SONTIMESERIES",typeof(fson_ts),getvardims(fson_ts))
      filevardef(ncdf,"DJFTIMESERIES",typeof(fdjf_ts),getvardims(fdjf_ts))
      filevardef(ncdf,"Year",typeof(fyr),getvardims(fyr))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      
    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"MAMTIMESERIES",fmam_ts)
      filevarattdef(ncdf,"JJATIMESERIES",fjja_ts)
      filevarattdef(ncdf,"SONTIMESERIES",fson_ts)
      filevarattdef(ncdf,"DJFTIMESERIES",fdjf_ts)
      filevarattdef(ncdf,"Year",fyr)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model		= (/model/)
      ncdf->lat			= (/lat/)
      ncdf->lon	 		= (/lon/)
      ncdf->MAMTIMESERIES	= (/fmam_ts/)
      ncdf->JJATIMESERIES	= (/fjja_ts/)
      ncdf->SONTIMESERIES	= (/fson_ts/)
      ncdf->DJFTIMESERIES	= (/fdjf_ts/)
      ncdf->Year		= (/fyr/)
      ncdf->ModelName		= (/modelname/)
      
    ;--- delete vars
      delete([/fmam_ts,fjja_ts,fson_ts,fdjf_ts/])
    
    end do   ; scn  
     
  end do  ; variables

end


 

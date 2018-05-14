load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; seas_hist_ts.ncl 
; EA Burakowski
; 2018-05-10
; 
; seas_hist_ts.ncl calculates the seasonally averaged time series in 
;	temperature (min & max) and precipitation (total) for 
;	CMIP5 LOCAlized Constructed Analogs (LOCA) daily downscaled data. 
;	Inputs include daily LOCA data for three variables:
;		(1) tasmin (surface minimum air temperature)
;		(2) tasmax (surface maximum air temperature)
;		(3) precipitation (daily total)
;	for 29 CMIP5 models (see below, mdl), 1980-2005
;
; Output includes:
;   .nc files for seasonal 25-yr time series in each variable.
;   Abbreviations for the seasons are as follows:
; 	- DJF = winter
;	- MAM = spring
;	- JJA = summer
;	- SON = fall
;
; Use plotSEAS.ncl and plotSEAS_trends.ncl to plot up the .nc files
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

  ;--- Write model names as strings to netCDF4 
  model	= ispan(0,dimsizes(mdl)-1,1)
  model!0 = "model"
  modelname = new(dimsizes(mdl),"string")
  do i = 0,dimsizes(mdl)-1  
    modelname(i)	= str_join((/tostring(model(i)), mdl(i)/)," - ")
  end do

  modelname!0 = "model"

  hyr   = ispan(1980,2005,1)
  hyr!0	= "year"
  hyr&year = hyr
  mon	= ispan(1,12,1)

;---- Output directories

  odir	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/"	


;---- Seasonal Total and Averages
;
;      Loop over models and years to calculate seasonal average
;	(tasmax, tasmin) and total (pr) for each year in record
;
;----

;---- Load one historical file to get the lat/lon dims

  tfil 	= addfile(dir+"pr_day_ACCESS1-0_historical_r1i1p1_19800101-19801231.16th.nc","r")
  tmp1	= tfil->$var(0)$
  tvar  = tmp1(lat|:,lon|:,time|:)
  dims	= dimsizes(tvar)
  lat	= tvar&lat
  lon	= tvar&lon

printVarSummary(tvar)

;--- Loop over var, mod, scen daily files to calculate annual average
  do ivar = 0,dimsizes(var)-1
    print("----  Working on "+var(ivar)+"  -----")
      ;--- Initialize historical seasonal time sereis
        hmam_ts       = new((/dims(0),dims(1),dimsizes(hyr),dimsizes(mdl)/),float,-9999)

      ;--- Assign coordinates
        hmam_ts!0     = "lat"
        hmam_ts&lat   = lat
        hmam_ts!1     = "lon"
        hmam_ts&lon   = lon
        hmam_ts!2     = "year"
        hmam_ts&year  = hyr
        hmam_ts!3     = "model"
        hmam_ts&model = model

      ;--- Initialize other seasons based on spring (hmam_ts)
        hjja_ts       = hmam_ts
        hson_ts       = hmam_ts
        hdjf_ts       = hmam_ts
	  
	;--- Loop over models  
        do imdl = 0,dimsizes(mdl)-1
	  print("----  Working on "+mdl(imdl)+"  -----")

        ;--- If/else statement for Precip vs. Tmin/Tmax
        if ivar.eq.0

   ;-----------------------------
   ;
   ;   PRECIP CALCS
   ;
   ;----------------------------- 
	    
          ;--- Loop over years and calculate season total precip
	  do iyr = 0,dimsizes(hyr)-1
	  
	    ;--- if it's first year of the analysis, skip winter (@_FillValue)
	    if iyr.eq.0
              fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_historical_r*_"+hyr(iyr)+"0101*.nc")
	      a		= addfile(dir+fil,"r")
              v		= a->$var(ivar)$
	    
	    ;--- convert kg m-2 s-2 to mm/day
	      vmm1	= v
	      vmm1	= v*86400
	      vmm	= vmm1(lat|:, lon|:, time|:)
  
            ;--- identify months and seasons in file
	      time	= a->time
	      utc_date	= cd_calendar(time,0)
	      mm	= tointeger(utc_date(:,1))  ; really only need month 
	      mam	= ind(mm.ge.3 .and. mm.le.5)
	      jja	= ind(mm.ge.6 .and. mm.le.8)
	      son	= ind(mm.ge.9 .and. mm.le.11)
	
	    ;--- calculate the seasonal total, rearrange dims
	      hmam_ts(:,:,iyr,imdl)	= dim_sum_n_Wrap(vmm(:,:,mam),2)
	      hjja_ts(:,:,iyr,imdl)	= dim_sum_n_Wrap(vmm(:,:,jja),2)
	      hson_ts(:,:,iyr,imdl)	= dim_sum_n_Wrap(vmm(:,:,son),2)
	      hdjf_ts(:,:,iyr,imdl)	= hdjf_ts@_FillValue

	    ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
              delete([/fil,a,v,vmm1,vmm,time,utc_date,mm,mam,jja,son/])
	  
	  else  ; if it's not the first year
	    
	    ;--- 
	      fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_historical_r*_"+hyr(iyr)+"0101*.nc")
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
	      fil_dec         = systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_historical_r*_"+hyr(iyr-1)+"0101*.nc")
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
              hmam_ts(:,:,iyr,imdl)   = dim_sum_n_Wrap(vmm(:,:,mam),2)
              hjja_ts(:,:,iyr,imdl)   = dim_sum_n_Wrap(vmm(:,:,jja),2)
              hson_ts(:,:,iyr,imdl)   = dim_sum_n_Wrap(vmm(:,:,son),2)
              hdjf_ts(:,:,iyr,imdl)   = dim_sum_n_Wrap(djf,2)

            ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
              delete([/fil,a,v,vmm1,vmm,time,utc_date,mm,mam,jja,son,jf/])
              delete([/fil_dec,a_dec,v_dec,vmm1_dec,time_dec,utc_dec,mm_dec,dec,djf,djf1/])

	    end if   ; first year loop	  
    
          end do   ; historical years
	
	else   ; pr/temp
	
   ;-----------------------------
   ;
   ;   TEMPERATURE CALCS
   ;
   ;-----------------------------   
	    
          ;--- Loop over years and calculate season total precip
	  do iyr = 0,dimsizes(hyr)-1
	  
	    ;--- if it's first year of the analysis, skip winter (@_FillValue)
	    if iyr.eq.0
              fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_historical_r*_"+hyr(iyr)+"0101*.nc")
	      a		= addfile(dir+fil,"r")
              v		= a->$var(ivar)$
	    
	    ;--- convert K to C
	      vC1	= v
	      vC1	= v - 273.15
	      vC	= vC1(lat|:, lon|:, time|:)
  
            ;--- identify months and seasons in file
	      time	= a->time
	      utc_date	= cd_calendar(time,0)
	      mm	= tointeger(utc_date(:,1))  ; really only need month 
	      mam	= ind(mm.ge.3 .and. mm.le.5)
	      jja	= ind(mm.ge.6 .and. mm.le.8)
	      son	= ind(mm.ge.9 .and. mm.le.11)
	
	    ;--- calculate the seasonal total, rearrange dims
	      hmam_ts(:,:,iyr,imdl)	= dim_avg_n_Wrap(vC(:,:,mam),2)
	      hjja_ts(:,:,iyr,imdl)	= dim_avg_n_Wrap(vC(:,:,jja),2)
	      hson_ts(:,:,iyr,imdl)	= dim_avg_n_Wrap(vC(:,:,son),2)
	      hdjf_ts(:,:,iyr,imdl)	= hdjf_ts@_FillValue

	    ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
              delete([/fil,a,v,vC1,vC,time,utc_date,mm,mam,jja,son/])
	  
	  else  ; if it's not the first year
	    
	    ;--- 
	      fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_historical_r*_"+hyr(iyr)+"0101*.nc")
	      a		= addfile(dir+fil,"r")
              v		= a->$var(ivar)$
	    
	    ;--- convert K to C
	      vC1	= v
	      vC1	= v - 273.15
	    
	    ;--- identify months and seasons in file
	      time	= a->time
	      utc_date	= cd_calendar(time,0)
	      mm	= tointeger(utc_date(:,1))  ; really only need month 
	      mam	= ind(mm.ge.3 .and. mm.le.5)
	      jja	= ind(mm.ge.6 .and. mm.le.8)
	      son	= ind(mm.ge.9 .and. mm.le.11)
	      jf	= ind(mm.le.2) ; jan & feb only
	      
	    ;--- December of previous year
	      fil_dec         = systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_historical_r*_"+hyr(iyr-1)+"0101*.nc")
              a_dec           = addfile(dir+fil_dec,"r")
              v_dec           = a_dec->$var(ivar)$
              vC1_dec         = v_dec
              vC1_dec         = v_dec - 273.15
              time_dec        = a_dec->time
              utc_dec         = cd_calendar(time_dec,0)
              mm_dec          = tointeger(utc_dec(:,1))
              dec             = ind(mm_dec.eq.12)

            ;--- combine december of previous year (iyr-1) with jan-feb of current year (iyr)
              djf1            = array_append_record(vC1_dec(dec,:,:),vC1(jf,:,:),0)

            ;--- now rearrange dims
              vC              = vC1(lat|:, lon|:, time|:)
              djf             = djf1(lat|:, lon|:, time|:)

            ;--- Calculate seasonal averages
              hmam_ts(:,:,iyr,imdl)   = dim_avg_n_Wrap(vC(:,:,mam),2)
              hjja_ts(:,:,iyr,imdl)   = dim_avg_n_Wrap(vC(:,:,jja),2)
              hson_ts(:,:,iyr,imdl)   = dim_avg_n_Wrap(vC(:,:,son),2)
              hdjf_ts(:,:,iyr,imdl)   = dim_avg_n_Wrap(djf,2)

            ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
              delete([/fil,a,v,vC1,vC,time,utc_date,mm,mam,jja,son,jf/])
              delete([/fil_dec,a_dec,v_dec,vC1_dec,time_dec,utc_dec,mm_dec,dec,djf,djf1/])

	    end if   ; first year

          end do   ; historical years
	  
        end if ; pr/temp 

      end do ; models 

  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ;
  ;--- create .nc file of spring trends
  ;
  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ;--- assign coordinate vars for climatological averages
      hmam_ts@short_name 	= "LOCA spring time series of "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hmam_ts@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) Spring "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+" time series"
      if ivar.eq.0
        hmam_ts@units		= "mm"
      else
        hmam_ts@units 		= "degrees C"
      end if

    ;--- assign coordinate vars for climatological averages
      hjja_ts@short_name 	= "LOCA spring time series of "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hjja_ts@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) Spring "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+" time series"
      if ivar.eq.0
        hjja_ts@units		= "mm"
      else
        hjja_ts@units 		= "degrees C"
      end if
    
    ;--- assign coordinate vars for climatological averages
      hson_ts@short_name 	= "LOCA spring time series of "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hson_ts@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) Spring "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+" time series"
      if ivar.eq.0
        hson_ts@units		= "mm"
      else
        hson_ts@units 		= "degrees C"
      end if

    ;--- assign coordinate vars for climatological averages
      hdjf_ts@short_name 	= "LOCA spring time series of "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hdjf_ts@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) Spring "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+" time series"
      if ivar.eq.0
        hdjf_ts@units		= "mm"
      else
        hdjf_ts@units 		= "degrees C"
      end if
      
    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")

    ;--- create file name
      fn 		= "SEAS_"+var(ivar)+"_historical_"+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+"_timeseries"
print("Filename = "+fn)

    ;--- remove old file, if it exists
      system("/bin/rm -f "+odir+fn+".nc")
      
    ;--- create new file   
      ncdf		= addfile(odir+fn+".nc","c")	
      fAtt		= True
      fAtt@title	= fn
      fAtt@orig_fil	= var(ivar)+"_day_LOCA_CMIP5_models_historical_"+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      fAtt@Conventions	= "COARDS/CF-1.0"
      fAtt@creation_date= systemfunc("date")
      fAtt@author	= "Elizabeth Burakowski (elizabeth.burakowski@unh.edu)"
      
      fAtt@description 	= "Seasonal "+var(ivar)+" time series from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      
    ;--- file attributes, size
      fileattdef(ncdf,fAtt)
      dimNames		= (/"lat", "lon", "year", "model"/)
      dimSizes		= (/ dims(0), dims(1), dimsizes(hyr), dimsizes(model) /)
      dimUnlim          = (/ False, False, False, False /)
      filedimdef(ncdf,dimNames,dimSizes,dimUnlim)

    ;--- Define file variables
      filevardef(ncdf,"model",typeof(model),getvardims(model))
      filevardef(ncdf,"lat",typeof(lat),getvardims(lat))
      filevardef(ncdf,"lon",typeof(lon),getvardims(lon))
      filevardef(ncdf,"MAMTIMESERIES",typeof(hmam_ts),getvardims(hmam_ts))
      filevardef(ncdf,"JJATIMESERIES",typeof(hjja_ts),getvardims(hjja_ts))
      filevardef(ncdf,"SONTIMESERIES",typeof(hson_ts),getvardims(hson_ts))
      filevardef(ncdf,"DJFTIMESERIES",typeof(hdjf_ts),getvardims(hdjf_ts))
      filevardef(ncdf,"Year",typeof(hyr),getvardims(hyr))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      
    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"MAMTIMESERIES",hmam_ts)
      filevarattdef(ncdf,"JJATIMESERIES",hjja_ts)
      filevarattdef(ncdf,"SONTIMESERIES",hson_ts)
      filevarattdef(ncdf,"DJFTIMESERIES",hdjf_ts)
      filevarattdef(ncdf,"Year",hyr)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model		= (/model/)
      ncdf->lat			= (/lat/)
      ncdf->lon	 		= (/lon/)
      ncdf->MAMTIMESERIES	= (/hmam_ts/)
      ncdf->JJATIMESERIES	= (/hjja_ts/)
      ncdf->SONTIMESERIES	= (/hson_ts/)
      ncdf->DJFTIMESERIES	= (/hdjf_ts/)
      ncdf->ModelName		= (/modelname/)
      
    ;--- delete vars
      delete([/hmam_ts,hjja_ts,hson_ts,hdjf_ts/])
      	     
  end do  ; variables

end


 

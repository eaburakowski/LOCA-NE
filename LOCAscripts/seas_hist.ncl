load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; seas_hist.ncl 
; EA Burakowski
; 2017-07-17
; 
; seas_hist.ncl calculates the climatological average & linear trend in 
;	seasonal temperature (min & max) and precipitation (total) for 
;	CMIP5 LOCAlized Constructed Analogs (LOCA) daily downscaled data. 
;	Inputs include daily LOCA data for three variables:
;		(1) tasmin (surface minimum air temperature)
;		(2) tasmax (surface maximum air temperature)
;		(3) precipitation (daily total)
;	for 29 CMIP5 models (see below, mdl), 1980-2005
;
; Output includes:
;   .nc files for seasonal 30-yr climatological mean and trend in each variable.
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

  ;--- Annoying stuff to write model names as strings to netCDF4 
  model	= ispan(0,dimsizes(mdl)-1,1)
  model!0 = "model"
  modelname = new(dimsizes(mdl),"string")
  do i = 0,dimsizes(mdl)-1  
    modelname(i)	= str_join((/tostring(model(i)), mdl(i)/)," - ")
  end do

  modelname!0 = "model"

  scn	= (/"historical"/)
  hyr   = ispan(1980,2005,1)
  mon	= ispan(1,12,1)
  early = ind(hyr.ge.2010 .and. hyr.le.2039)
  mid   = ind(hyr.ge.2040 .and. hyr.le.2069)
  late  = ind(hyr.ge.2070 .and. hyr.le.2099)

;---- Output directories

  odir	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/"	


;---- Seasonal Total and Averages
;
;      Loop over models and years to calculate annual average
;	tasmax, tasmin, pr for each year in record
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
	
	;--- Initialize historical seasonal trends
	  hmam_cl	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
	  hmam_tr 	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
	  hjja_cl	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
	  hjja_tr 	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
	  hson_cl	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
	  hson_tr 	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
	  hdjf_cl	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
	  hdjf_tr	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
          
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
	 
	  ;--- intialize hmam (historical spring variable time series)
	    hmam 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
            copy_VarCoords_1(tvar,hmam)
   	    hmam!2	= "year"
	    hmam&year	= hyr
	    
	  ;--- intialize hjja (historical summer variable time series)
	    hjja 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
            copy_VarCoords_1(tvar,hjja)
   	    hjja!2	= "year"
	    hjja&year	= hyr
	   
	  ;--- intialize hson (historical fall variable time series)
	    hson 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
            copy_VarCoords_1(tvar,hson)
   	    hson!2	= "year"
	    hson&year	= hyr
	    
	  ;--- intialize hdjf (historical fall variable time series)
	    hdjf 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
            copy_VarCoords_1(tvar,hdjf)
   	    hdjf!2	= "year"
	    hdjf&year	= hyr  
	    
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
	      hmam(:,:,iyr)	= dim_sum_n_Wrap(vmm(:,:,mam),2)
	      hjja(:,:,iyr)	= dim_sum_n_Wrap(vmm(:,:,jja),2)
	      hson(:,:,iyr)	= dim_sum_n_Wrap(vmm(:,:,son),2)
	      hdjf(:,:,iyr)	= hdjf@_FillValue

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
              hmam(:,:,iyr)   = dim_sum_n_Wrap(vmm(:,:,mam),2)
              hjja(:,:,iyr)   = dim_sum_n_Wrap(vmm(:,:,jja),2)
              hson(:,:,iyr)   = dim_sum_n_Wrap(vmm(:,:,son),2)
              hdjf(:,:,iyr)   = dim_sum_n_Wrap(djf,2)

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
	 
	  ;--- intialize hmam (historical spring variable time series)
	    hmam 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
            copy_VarCoords_1(tvar,hmam)
   	    hmam!2	= "year"
	    hmam&year	= hyr
	    
	  ;--- intialize hjja (historical summer variable time series)
	    hjja 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
            copy_VarCoords_1(tvar,hjja)
   	    hjja!2	= "year"
	    hjja&year	= hyr
	   
	  ;--- intialize hson (historical fall variable time series)
	    hson 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
            copy_VarCoords_1(tvar,hson)
   	    hson!2	= "year"
	    hson&year	= hyr
	  
	  ;--- intialize hdjf (historical fall variable time series)
	    hdjf 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
            copy_VarCoords_1(tvar,hdjf)
   	    hdjf!2	= "year"
	    hdjf&year	= hyr  
	    
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
	      hmam(:,:,iyr)	= dim_avg_n_Wrap(vC(:,:,mam),2)
	      hjja(:,:,iyr)	= dim_avg_n_Wrap(vC(:,:,jja),2)
	      hson(:,:,iyr)	= dim_avg_n_Wrap(vC(:,:,son),2)
	      hdjf(:,:,iyr)	= hdjf@_FillValue

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
              djf            = djf1(lat|:, lon|:, time|:)

            ;--- Calculate seasonal sums
              hmam(:,:,iyr)   = dim_avg_n_Wrap(vC(:,:,mam),2)
              hjja(:,:,iyr)   = dim_avg_n_Wrap(vC(:,:,jja),2)
              hson(:,:,iyr)   = dim_avg_n_Wrap(vC(:,:,son),2)
              hdjf(:,:,iyr)   = dim_avg_n_Wrap(djf,2)

            ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
              delete([/fil,a,v,vC1,vC,time,utc_date,mm,mam,jja,son,jf/])
              delete([/fil_dec,a_dec,v_dec,vC1_dec,time_dec,utc_dec,mm_dec,dec,djf,djf1/])

	    end if   ; first year

          end do   ; historical years
	  
        end if ; pr/temp 

	;--- calculate annual climatological mean, historical
   	  hmam_cl(:,:,imdl)	= dim_avg_n(hmam,2)
	  hjja_cl(:,:,imdl)	= dim_avg_n(hjja,2)
	  hson_cl(:,:,imdl)	= dim_avg_n(hson,2)
	  hdjf_cl(:,:,imdl)	= dim_avg_n(hdjf,2)
 
        ;--- calculate annual trend, historical
          hmam_tr(:,:,imdl)	= regCoef_n(hyr,hmam,0,2)
	  hjja_tr(:,:,imdl)	= regCoef_n(hyr,hjja,0,2)
	  hson_tr(:,:,imdl)	= regCoef_n(hyr,hson,0,2)
	  hdjf_tr(:,:,imdl)	= regCoef_n(hyr,hdjf,0,2) ; double check this...
	   
      end do ; models 

  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ;
  ;--- create .nc file of spring trends
  ;
  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ;--- assign coordinate vars for climatological averages
      copy_VarCoords_1(tvar,hmam_cl)
      hmam_cl!2			= "model"
      hmam_cl&model		= model
      hmam_cl@short_name 	= "LOCA Climatological mean "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hmam_cl@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Spring "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      if ivar.eq.0
        hmam_cl@units		= "mm"
      else
        hmam_cl@units 		= "degrees C"
      end if

    ;--- assign coordinate vars for trends
      copy_VarCoords_1(tvar,hmam_tr)
      hmam_tr!2			= "model"
      hmam_tr&model		= model
      hmam_tr@short_name 	= "LOCA Spring "+var(ivar)+" trend, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hmam_tr@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) Spring Temperature ("+var(ivar)+"), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      if ivar.eq.0
        hmam_tr@units		= "mm year-1"
      else
        hmam_tr@units 		= "degrees C per year"
      end if

    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")

    ;--- create file name
      fn 		= "MAM_"+var(ivar)+"_historical_"+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
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
      
      fAtt@description 	= "Spring "+var(ivar)+" climatological mean and trend calculated from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      
    ;--- file attributes, size
      fileattdef(ncdf,fAtt)
      dimNames		= (/"lat", "lon", "model"/)
      dimSizes		= (/ dims(0), dims(1), dimsizes(model) /)
      dimUnlim          = (/ False, False, False /)
      filedimdef(ncdf,dimNames,dimSizes,dimUnlim)

    ;--- Define file variables
      filevardef(ncdf,"model",typeof(model),getvardims(model))
      filevardef(ncdf,"lat",typeof(lat),getvardims(lat))
      filevardef(ncdf,"lon",typeof(lon),getvardims(lon))
      filevardef(ncdf,"MAMCLIM",typeof(hmam_cl),getvardims(hmam_cl))
      filevardef(ncdf,"MAMTREND",typeof(hmam_tr),getvardims(hmam_tr))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      
    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"MAMCLIM",hmam_cl)
      filevarattdef(ncdf,"MAMTREND",hmam_tr)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model	= (/model/)
      ncdf->lat		= (/lat/)
      ncdf->lon	 	= (/lon/)
      ncdf->MAMCLIM	= (/hmam_cl/)
      ncdf->MAMTREND	= (/hmam_tr/)
      ncdf->ModelName	= (/modelname/)
      
  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ;
  ;--- create .nc file of summer trends
  ;
  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
   ;--- assign coordinate vars for climatological averages   
    copy_VarCoords_1(tvar,hjja_cl)
      hjja_cl!2			= "model"
      hjja_cl&model		= model
      hjja_cl@short_name 	= "LOCA Climatological mean "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hjja_cl@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Summer "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      if ivar.eq.0
        hjja_cl@units		= "mm"
      else
        hjja_cl@units 		= "degrees C"
      end if

    ;--- assign coordinate vars for trends
      copy_VarCoords_1(tvar,hjja_tr)
      hjja_tr!2			= "model"
      hjja_tr&model		= model
      hjja_tr@short_name 	= "LOCA Summer "+var(ivar)+" trend, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hjja_tr@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) Summer Temperature ("+var(ivar)+"), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      if ivar.eq.0 
        hjja_tr@units		= "mm year-1"
      else  
        hjja_tr@units 		= "degrees C per year"
      end if

    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")

    ;--- create file name
      fn 		= "JJA_"+var(ivar)+"_historical_"+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
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
      
      fAtt@description 	= "Summer "+var(ivar)+" climatological mean and trend calculated from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      
    ;--- file attributes, size
      fileattdef(ncdf,fAtt)
      dimNames		= (/"lat", "lon", "model"/)
      dimSizes		= (/ dims(0), dims(1), dimsizes(model) /)
      dimUnlim          = (/ False, False, False /)
      filedimdef(ncdf,dimNames,dimSizes,dimUnlim)

    ;--- Define file variables
      filevardef(ncdf,"model",typeof(model),getvardims(model))
      filevardef(ncdf,"lat",typeof(lat),getvardims(lat))
      filevardef(ncdf,"lon",typeof(lon),getvardims(lon))
      filevardef(ncdf,"JJACLIM",typeof(hjja_cl),getvardims(hjja_cl))
      filevardef(ncdf,"JJATREND",typeof(hjja_tr),getvardims(hjja_tr))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      
    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"JJACLIM",hjja_cl)
      filevarattdef(ncdf,"JJATREND",hjja_tr)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model	= (/model/)
      ncdf->lat		= (/lat/)
      ncdf->lon	 	= (/lon/)
      ncdf->JJACLIM	= (/hjja_cl/)
      ncdf->JJATREND	= (/hjja_tr/)
      ncdf->ModelName	= (/modelname/)
      
  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ;
  ;--- create .nc file of Fall trends
  ;
  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
   ;--- assign coordinate vars for climatological averages   
    copy_VarCoords_1(tvar,hson_cl)
      hson_cl!2			= "model"
      hson_cl&model		= model
      hson_cl@short_name 	= "LOCA Climatological mean "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hson_cl@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Fall "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      if ivar.eq.0
        hson_cl@units		= "mm"
      else      
        hson_cl@units 		= "degrees C"
      end if

    ;--- assign coordinate vars for trends
      copy_VarCoords_1(tvar,hson_tr)
      hson_tr!2			= "model"
      hson_tr&model		= model
      hson_tr@short_name 	= "LOCA Fall "+var(ivar)+" trend, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hson_tr@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) Fall Temperature ("+var(ivar)+"), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hson_tr@units 		= "degrees C per year"

    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")

    ;--- create file name
      fn 		= "SON_"+var(ivar)+"_historical_"+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
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
      
      fAtt@description 	= "Fall "+var(ivar)+" climatological mean and trend calculated from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      
    ;--- file attributes, size
      fileattdef(ncdf,fAtt)
      dimNames		= (/"lat", "lon", "model"/)
      dimSizes		= (/ dims(0), dims(1), dimsizes(model) /)
      dimUnlim          = (/ False, False, False /)
      filedimdef(ncdf,dimNames,dimSizes,dimUnlim)

    ;--- Define file variables
      filevardef(ncdf,"model",typeof(model),getvardims(model))
      filevardef(ncdf,"lat",typeof(lat),getvardims(lat))
      filevardef(ncdf,"lon",typeof(lon),getvardims(lon))
      filevardef(ncdf,"SONCLIM",typeof(hson_cl),getvardims(hson_cl))
      filevardef(ncdf,"SONTREND",typeof(hson_tr),getvardims(hson_tr))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      
    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"SONCLIM",hson_cl)
      filevarattdef(ncdf,"SONTREND",hson_tr)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model	= (/model/)
      ncdf->lat		= (/lat/)
      ncdf->lon	 	= (/lon/)
      ncdf->SONCLIM	= (/hson_cl/)
      ncdf->SONTREND	= (/hson_tr/)
      ncdf->ModelName	= (/modelname/)


;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ;
  ;--- create .nc file of Winter trends
  ;
  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
   ;--- assign coordinate vars for climatological averages   
    copy_VarCoords_1(tvar,hdjf_cl)
      hdjf_cl!2			= "model"
      hdjf_cl&model		= model
      hdjf_cl@short_name 	= "LOCA Climatological mean "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hdjf_cl@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Winter "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      if ivar.eq.0
        hdjf_cl@units		= "mm"
      else      
        hdjf_cl@units 		= "degrees C"
      end if

    ;--- assign coordinate vars for trends
      copy_VarCoords_1(tvar,hdjf_tr)
      hdjf_tr!2			= "model"
      hdjf_tr&model		= model
      hdjf_tr@short_name 	= "LOCA Winter "+var(ivar)+" trend, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hdjf_tr@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) Winter Temperature ("+var(ivar)+"), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hdjf_tr@units 		= "degrees C per year"

    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")

    ;--- create file name
      fn 		= "DJF_"+var(ivar)+"_historical_"+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
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
      
      fAtt@description 	= "Winter "+var(ivar)+" climatological mean and trend calculated from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      
    ;--- file attributes, size
      fileattdef(ncdf,fAtt)
      dimNames		= (/"lat", "lon", "model"/)
      dimSizes		= (/ dims(0), dims(1), dimsizes(model) /)
      dimUnlim          = (/ False, False, False /)
      filedimdef(ncdf,dimNames,dimSizes,dimUnlim)

    ;--- Define file variables
      filevardef(ncdf,"model",typeof(model),getvardims(model))
      filevardef(ncdf,"lat",typeof(lat),getvardims(lat))
      filevardef(ncdf,"lon",typeof(lon),getvardims(lon))
      filevardef(ncdf,"DJFCLIM",typeof(hdjf_cl),getvardims(hdjf_cl))
      filevardef(ncdf,"DJFTREND",typeof(hdjf_tr),getvardims(hdjf_tr))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      
    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"DJFCLIM",hdjf_cl)
      filevarattdef(ncdf,"DJFTREND",hdjf_tr)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model	= (/model/)
      ncdf->lat		= (/lat/)
      ncdf->lon	 	= (/lon/)
      ncdf->DJFCLIM	= (/hdjf_cl/)
      ncdf->DJFTREND	= (/hdjf_tr/)
      ncdf->ModelName	= (/modelname/)
      
    ;--- delete vars
      delete([/hmam_cl,hmam_tr,hjja_cl,hjja_tr,hson_cl,hson_tr/])
      	     
  end do  ; variables

end


 

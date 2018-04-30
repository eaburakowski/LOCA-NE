load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; seas_future.ncl 
; EA Burakowski
; 2017-04-26
; 
; seas_future.ncl calculates the climatological average & linear trend in 
;	seasonal temperature (min & max) and precipitation (total) fori
;	rcp4.5 and rcp8.5 for early (2010-2039), mid (2040-2069) and late 
;	(2070-2099) century, using CMIP5 LOCAlized Constructed Analogs (LOCA) 
;	daily downscaled data (16th degree resolution). 
;	Inputs include daily LOCA data for three variables:
;		(1) tasmin (surface minimum air temperature)
;		(2) tasmax (surface maximum air temperature)
;		(3) precipitation (daily total)
;	for 29 CMIP5 models (see below, mdl), 2006-2099
;
; Output includes:
;   .nc files for seasonal 30-yr climatological mean and trend in each variable.
;   Abbreviations for the seasons are as follows:
; 	- DJF = winter
;	- MAM = spring
;	- JJA = summer
;	- SON = fall
;
; Use plotSEAS_future.ncl and plotSEAS_individ.ncl to plot up the .nc files
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
  fyr   = ispan(2006,2099,1)
  mon	= ispan(1,12,1)
  early = ind(fyr.ge.2010 .and. fyr.le.2039)
  mid   = ind(fyr.ge.2040 .and. fyr.le.2069)
  late  = ind(fyr.ge.2070 .and. fyr.le.2099)

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
	  fmam_early	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
	  fmam_mid 	= fmam_early
	  fmam_late	= fmam_early
	
	  fjja_early	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
	  fjja_mid 	= fjja_early
	  fjja_late	= fjja_early
	  
	  fson_early	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
	  fson_mid 	= fson_early
	  fson_late	= fson_early
	  
	  fdjf_early	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
	  fdjf_mid 	= fdjf_early
	  fdjf_late	= fdjf_early	

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
	 
	    ;--- intialize fmam (future spring variable time series)
	      fmam 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
              copy_VarCoords_1(tvar,fmam)
   	      fmam!2	= "year"
	      fmam&year	= fyr
	    
	    ;--- intialize fjja (future summer variable time series)
	      fjja 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
              copy_VarCoords_1(tvar,fjja)
   	      fjja!2	= "year"
	      fjja&year	= fyr
	   
	    ;--- intialize fson (future fall variable time series)
	      fson 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
              copy_VarCoords_1(tvar,fson)
   	      fson!2	= "year"
	      fson&year	= fyr
	    
	    ;--- intialize fdjf (future fall variable time series)
	      fdjf 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
              copy_VarCoords_1(tvar,fdjf)
   	      fdjf!2	= "year"
	      fdjf&year	= fyr  
	    
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
	        fmam(:,:,iyr)	= dim_sum_n_Wrap(vmm(:,:,mam),2)
	        fjja(:,:,iyr)	= dim_sum_n_Wrap(vmm(:,:,jja),2)
	        fson(:,:,iyr)	= dim_sum_n_Wrap(vmm(:,:,son),2)
	        fdjf(:,:,iyr)	= fdjf@_FillValue

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
                djf             = array_append_record(vmm1_dec(dec,:,:),vmm1(jf,:,:),0)

              ;--- now rearrange dims
                vmm             = vmm1(lat|:, lon|:, time|:)
                djf1            = djf(lat|:, lon|:, time|:)

              ;--- Calculate seasonal sums
                fmam(:,:,iyr)   = dim_sum_n_Wrap(vmm(:,:,mam),2)
                fjja(:,:,iyr)   = dim_sum_n_Wrap(vmm(:,:,jja),2)
                fson(:,:,iyr)   = dim_sum_n_Wrap(vmm(:,:,son),2)
                fdjf(:,:,iyr)   = dim_sum_n_Wrap(djf1,2)
		
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
	 
	  ;--- intialize fmam (future spring variable time series)
	    fmam 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
            copy_VarCoords_1(tvar,fmam)
   	    fmam!2	= "year"
	    fmam&year	= fyr
	    
	  ;--- intialize fjja (future summer variable time series)
	    fjja 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
            copy_VarCoords_1(tvar,fjja)
   	    fjja!2	= "year"
	    fjja&year	= fyr
	   
	  ;--- intialize fson (future fall variable time series)
	    fson 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
            copy_VarCoords_1(tvar,fson)
   	    fson!2	= "year"
	    fson&year	= fyr
	  
	  ;--- intialize fdjf (future fall variable time series)
	    fdjf 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
            copy_VarCoords_1(tvar,fdjf)
   	    fdjf!2	= "year"
	    fdjf&year	= fyr  
	    
          ;--- Loop over years and calculate season total precip
	    do iyr = 0,dimsizes(fyr)-1
	  
	      ;--- if it's first year of the analysis, skip winter (@_FillValue)
	      if iyr.eq.0
              fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_"+scn(iscn)+"_r*_"+fyr(iyr)+"0101*.nc")
	      a		= addfile(dir+fil,"r")
              v		= a->$var(ivar)$
printVarSummary(v)	    
	    ;--- convert K to C
	      vC1	= v
	      vC1	= v-273.15
	      vC	= vC1(lat|:, lon|:, time|:)
printVarSummary(vC)
print("Max vC: "+max(vC))
print("Min vC: "+min(vC))
  
            ;--- identify months and seasons in file
	      time	= a->time
	      utc_date	= cd_calendar(time,0)
	      mm	= tointeger(utc_date(:,1))  ; really only need month 
	      mam	= ind(mm.ge.3 .and. mm.le.5)
	      jja	= ind(mm.ge.6 .and. mm.le.8)
	      son	= ind(mm.ge.9 .and. mm.le.11)
	
	    ;--- calculate the seasonal total, rearrange dims
	      fmam(:,:,iyr)	= dim_avg_n_Wrap(vC(:,:,mam),2)
	      fjja(:,:,iyr)	= dim_avg_n_Wrap(vC(:,:,jja),2)
	      fson(:,:,iyr)	= dim_avg_n_Wrap(vC(:,:,son),2)
	      fdjf(:,:,iyr)	= fdjf@_FillValue

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
              fmam(:,:,iyr)   = dim_avg_n_Wrap(vC(:,:,mam),2)
              fjja(:,:,iyr)   = dim_avg_n_Wrap(vC(:,:,jja),2)
              fson(:,:,iyr)   = dim_avg_n_Wrap(vC(:,:,son),2)
              fdjf(:,:,iyr)   = dim_avg_n(djf,2)
print("Max DJF "+var(ivar)+" "+fyr(iyr)+" = "+max(fdjf(:,:,iyr)))
print("Min DJF "+var(ivar)+" "+fyr(iyr)+" = "+min(fdjf(:,:,iyr)))

            ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
              delete([/fil,a,v,vC1,vC,time,utc_date,mm,mam,jja,son,jf/])
              delete([/fil_dec,a_dec,v_dec,v1_dec,vC1_dec/])
	      delete([/time_dec,utc_dec,mm_dec,dec,djf,djf1/])

	    end if   ; first year

          end do   ; "+scn(iscn)+" years
	  
        end if ; pr/temp 

	;--- calculate seasonal climatological mean, future
   	  fmam_early(:,:,imdl)	= dim_avg_n(fmam(:,:,early),2)
	  fmam_mid(:,:,imdl)	= dim_avg_n(fmam(:,:,mid),2)
	  fmam_late(:,:,imdl)	= dim_avg_n(fmam(:,:,late),2)
	  
	  fjja_early(:,:,imdl)	= dim_avg_n(fjja(:,:,early),2)
	  fjja_mid(:,:,imdl)	= dim_avg_n(fjja(:,:,mid),2)
	  fjja_late(:,:,imdl)	= dim_avg_n(fjja(:,:,late),2)
	  
	  fson_early(:,:,imdl)	= dim_avg_n(fson(:,:,early),2)
	  fson_mid(:,:,imdl)	= dim_avg_n(fson(:,:,mid),2)
	  fson_late(:,:,imdl)	= dim_avg_n(fson(:,:,late),2)
	  
	  fdjf_early(:,:,imdl)	= dim_avg_n(fdjf(:,:,early),2)
	  fdjf_mid(:,:,imdl)	= dim_avg_n(fdjf(:,:,mid),2)
	  fdjf_late(:,:,imdl)	= dim_avg_n(fdjf(:,:,late),2)

print("Max DJF EARLY: "+max(fdjf_early))	   
      end do ; models 

  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ;
  ;--- create .nc file of spring trends
  ;
  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ;--- assign coordinate vars for climatological averages
    
      ;-- Early (2010-2039)    
      copy_VarCoords_1(tvar,fmam_early)
      fmam_early!2		= "model"
      fmam_early&model		= model
      fmam_early@short_name 	= scn(iscn)+": LOCA Climatological mean "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fmam_early@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Spring "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      if ivar.eq.0
        fmam_early@units		= "mm"
      else
        fmam_early@units 		= "degrees C"
      end if
      
      ;-- Mid (2040-2069)    
      copy_VarCoords_1(tvar,fmam_mid)
      fmam_mid!2		= "model"
      fmam_mid&model		= model
      fmam_mid@short_name 	= scn(iscn)+": LOCA Climatological mean "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fmam_mid@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Spring "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      if ivar.eq.0
        fmam_mid@units		= "mm"
      else
        fmam_mid@units 		= "degrees C"
      end if
      
      ;-- Late (2070-2099)    
      copy_VarCoords_1(tvar,fmam_late)
      fmam_late!2		= "model"
      fmam_late&model		= model
      fmam_late@short_name 	= scn(iscn)+": LOCA Climatological mean "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fmam_late@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Spring "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      if ivar.eq.0
        fmam_late@units		= "mm"
      else
        fmam_late@units 	= "degrees C"
      end if

    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")

    ;--- create file name
      fn 		= "MAM_"+var(ivar)+"_"+scn(iscn)+"_"+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+"_clim"
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
      
      fAtt@description 	= scn(iscn)+": Spring "+var(ivar)+" climatological mean and trend calculated from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      
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
      filevardef(ncdf,"MAMCLIM_EARLY",typeof(fmam_early),getvardims(fmam_early))
      filevardef(ncdf,"MAMCLIM_MID",typeof(fmam_mid),getvardims(fmam_mid))
      filevardef(ncdf,"MAMCLIM_LATE",typeof(fmam_late),getvardims(fmam_late))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      
    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"MAMCLIM_EARLY",fmam_early)
      filevarattdef(ncdf,"MAMCLIM_MID",fmam_mid)
      filevarattdef(ncdf,"MAMCLIM_LATE",fmam_late)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model		= (/model/)
      ncdf->lat			= (/lat/)
      ncdf->lon	 		= (/lon/)
      ncdf->MAMCLIM_EARLY	= (/fmam_early/)
      ncdf->MAMCLIM_MID		= (/fmam_mid/)
      ncdf->MAMCLIM_LATE	= (/fmam_late/)
      ncdf->ModelName		= (/modelname/)
      
    ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ;
  ;--- create .nc file of Summer trends
  ;
  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ;-- Early (2010-2039)    
      copy_VarCoords_1(tvar,fjja_early)
      fjja_early!2		= "model"
      fjja_early&model		= model
      fjja_early@short_name 	= scn(iscn)+": LOCA Climatological mean "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fjja_early@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Summer "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      if ivar.eq.0
        fjja_early@units		= "mm"
      else
        fjja_early@units 		= "degrees C"
      end if
      
      ;-- Mid (2040-2069)    
      copy_VarCoords_1(tvar,fjja_mid)
      fjja_mid!2		= "model"
      fjja_mid&model		= model
      fjja_mid@short_name 	= scn(iscn)+": LOCA Climatological mean "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fjja_mid@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Spring "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      if ivar.eq.0
        fjja_mid@units		= "mm"
      else
        fjja_mid@units 		= "degrees C"
      end if
      
      ;-- Late (2070-2099)    
      copy_VarCoords_1(tvar,fjja_late)
      fjja_late!2		= "model"
      fjja_late&model		= model
      fjja_late@short_name 	= scn(iscn)+": LOCA Climatological mean "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fjja_late@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Summer "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      if ivar.eq.0
        fjja_late@units		= "mm"
      else
        fjja_late@units 	= "degrees C"
      end if

    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")

    ;--- create file name
      fn 		= "JJA_"+var(ivar)+"_"+scn(iscn)+"_"+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+"_clim"
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
      
      fAtt@description 	= scn(iscn)+": Summer "+var(ivar)+" climatological mean and trend calculated from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      
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
      filevardef(ncdf,"JJACLIM_EARLY",typeof(fjja_early),getvardims(fjja_early))
      filevardef(ncdf,"JJACLIM_MID",typeof(fjja_mid),getvardims(fjja_mid))
      filevardef(ncdf,"JJACLIM_LATE",typeof(fjja_late),getvardims(fjja_late))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      
    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"JJACLIM_EARLY",fjja_early)
      filevarattdef(ncdf,"JJACLIM_MID",fjja_mid)
      filevarattdef(ncdf,"JJACLIM_LATE",fjja_late)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model		= (/model/)
      ncdf->lat			= (/lat/)
      ncdf->lon	 		= (/lon/)
      ncdf->JJACLIM_EARLY	= (/fjja_early/)
      ncdf->JJACLIM_MID		= (/fjja_mid/)
      ncdf->JJACLIM_LATE	= (/fjja_late/)
      ncdf->ModelName		= (/modelname/)
      
    ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ;
  ;--- create .nc file of Fall trends
  ;
  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ;-- Early (2010-2039)    
      copy_VarCoords_1(tvar,fson_early)
      fson_early!2		= "model"
      fson_early&model		= model
      fson_early@short_name 	= scn(iscn)+": LOCA Climatological mean "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fson_early@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Fall "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      if ivar.eq.0
        fson_early@units		= "mm"
      else
        fson_early@units 		= "degrees C"
      end if
      
      ;-- Mid (2040-2069)    
      copy_VarCoords_1(tvar,fson_mid)
      fson_mid!2		= "model"
      fson_mid&model		= model
      fson_mid@short_name 	= scn(iscn)+": LOCA Climatological mean "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fson_mid@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Fall "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      if ivar.eq.0
        fson_mid@units		= "mm"
      else
        fson_mid@units 		= "degrees C"
      end if
      
      ;-- Late (2070-2099)    
      copy_VarCoords_1(tvar,fson_late)
      fson_late!2		= "model"
      fson_late&model		= model
      fson_late@short_name 	= scn(iscn)+": LOCA Climatological mean "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fson_late@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Fall "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      if ivar.eq.0
        fson_late@units		= "mm"
      else
        fson_late@units 	= "degrees C"
      end if

    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")

    ;--- create file name
      fn 		= "SON_"+var(ivar)+"_"+scn(iscn)+"_"+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+"_clim"
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
      
      fAtt@description 	= scn(iscn)+": Fall "+var(ivar)+" climatological mean and trend calculated from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      
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
      filevardef(ncdf,"SONCLIM_EARLY",typeof(fson_early),getvardims(fson_early))
      filevardef(ncdf,"SONCLIM_MID",typeof(fson_mid),getvardims(fson_mid))
      filevardef(ncdf,"SONCLIM_LATE",typeof(fson_late),getvardims(fson_late))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      
    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"SONCLIM_EARLY",fson_early)
      filevarattdef(ncdf,"SONCLIM_MID",fson_mid)
      filevarattdef(ncdf,"SONCLIM_LATE",fson_late)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model		= (/model/)
      ncdf->lat			= (/lat/)
      ncdf->lon	 		= (/lon/)
      ncdf->SONCLIM_EARLY	= (/fson_early/)
      ncdf->SONCLIM_MID		= (/fson_mid/)
      ncdf->SONCLIM_LATE	= (/fson_late/)
      ncdf->ModelName		= (/modelname/)


  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ;
  ;--- create .nc file of Winter trends
  ;
  ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ;-- Early (2010-2039)    
      copy_VarCoords_1(tvar,fdjf_early)
      fdjf_early!2		= "model"
      fdjf_early&model		= model
      fdjf_early@short_name 	= scn(iscn)+": LOCA Climatological mean "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fdjf_early@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Winter "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      if ivar.eq.0
        fdjf_early@units		= "mm"
      else
        fdjf_early@units 		= "degrees C"
      end if
      
      ;-- Mid (2040-2069)    
      copy_VarCoords_1(tvar,fdjf_mid)
      fdjf_mid!2		= "model"
      fdjf_mid&model		= model
      fdjf_mid@short_name 	= scn(iscn)+": LOCA Climatological mean "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fdjf_mid@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Winter "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      if ivar.eq.0
        fdjf_mid@units		= "mm"
      else
        fdjf_mid@units 		= "degrees C"
      end if
      
      ;-- Late (2070-2099)    
      copy_VarCoords_1(tvar,fdjf_late)
      fdjf_late!2		= "model"
      fdjf_late&model		= model
      fdjf_late@short_name 	= scn(iscn)+": LOCA Climatological mean "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fdjf_late@long_name	= scn(iscn)+": Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean Winter "+var(ivar)+", "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      if ivar.eq.0
        fdjf_late@units		= "mm"
      else
        fdjf_late@units 	= "degrees C"
      end if

    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")

    ;--- create file name
      fn 		= "DJF_"+var(ivar)+"_"+scn(iscn)+"_"+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+"_clim"
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
      
      fAtt@description 	= scn(iscn)+": Winter "+var(ivar)+" climatological mean and trend calculated from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      
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
      filevardef(ncdf,"DJFCLIM_EARLY",typeof(fdjf_early),getvardims(fdjf_early))
      filevardef(ncdf,"DJFCLIM_MID",typeof(fdjf_mid),getvardims(fdjf_mid))
      filevardef(ncdf,"DJFCLIM_LATE",typeof(fdjf_late),getvardims(fdjf_late))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      
    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"DJFCLIM_EARLY",fdjf_early)
      filevarattdef(ncdf,"DJFCLIM_MID",fdjf_mid)
      filevarattdef(ncdf,"DJFCLIM_LATE",fdjf_late)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model		= (/model/)
      ncdf->lat			= (/lat/)
      ncdf->lon	 		= (/lon/)
      ncdf->DJFCLIM_EARLY	= (/fdjf_early/)
      ncdf->DJFCLIM_MID		= (/fdjf_mid/)
      ncdf->DJFCLIM_LATE	= (/fdjf_late/)
      ncdf->ModelName		= (/modelname/)
      
    ;--- delete vars
      delete([/fmam_early,fmam_mid,fmam_late,fjja_early,fjja_mid,fjja_late/])
      delete([/fson_early,fson_mid,fson_late,fdjf_early,fdjf_mid,fdjf_late/])
    
    end do   ; scn  
     
  end do  ; variables

end


 

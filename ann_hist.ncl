load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; ann_hist.ncl 
; EA Burakowski
; 2017-07-17
; 
; ann_hist.ncl calculates the climatological average & linear trend in 
;	annual temperature (min & max) and precipitation (total) for 
;	CMIP5 LOCAlized Constructed Analogs (LOCA) daily downscaled data. 
;	Inputs include daily LOCA data for three variables:
;		(1) tasmin (surface minimum air temperature)
;		(2) tasmax (surface maximum air temperature)
;		(3) precipitation (daily total)
;	for 29 CMIP5 models (see below, mdl), 1980-2005
;
; Output includes:
;   .nc files for 30-yr climatological mean and annual trend in each variable:
; 	- ann_pr_historical_1980-2005.nc
;	- ann_tasmax_historical_1980-2005.nc
;	- ann_tasmin_historical_1980-2005.nc
;	- ann_pr_historical_1980-2005_trend.nc
;	- ann_tasmax_historical_1980-2005_trend.nc
;	- ann_tasmin_historical_1980-2005_trend.nc
;
; Use plotANN.ncl and plotANN_trends.ncl to plot up the .nc files
;	

;======================================================================
;                         The main code
;======================================================================

begin

;---- LOCA directory

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

;---- Output directories

  odir	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/"	


;---- Annual Averages
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
  do ivar = 1,dimsizes(var)-1
    print("----  Working on "+var(ivar)+"  -----")
	
	;--- Initialize historical annual trend 
	  hann_cl	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
	  hann_tr 	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
         
        do imdl = 0,dimsizes(mdl)-1
	  print("----  Working on "+mdl(imdl)+"  -----")

      ;--- If/else statement for Precip vs. Tmin/Tmax
        if ivar.eq.0

   ;-----------------------------
   ;
   ;   PRECIP CALCS
   ;
   ;-----------------------------
	    
	  ;--- intialize hann (historical annual variable time series)
	    hann 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
            copy_VarCoords_1(tvar,hann)
   	    hann!2	= "year"
	    hann&year	= hyr

	  do iyr = 0,dimsizes(hyr)-1
            fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_historical_r*_"+hyr(iyr)+"0101*.nc")
	    a		= addfile(dir+fil,"r")
            v		= a->$var(ivar)$
	    
	  ;--- convert kg m-2 s-2 to mm/day
	    vmm		= v
	    vmm		= v*86400
  
  print("Maximum precip (mm/day) = "+max(vmm))
  print("Minimum precip (mm/day) = "+min(vmm))  
	
	  ;--- calculate the annual total, rearrange dims
	    hann(:,:,iyr)	= dim_sum_n_Wrap(vmm(lat|:,lon|:,time|:),2)

	  ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
            delete([/fil,a,v,vmm/])

          end do   ; historical years
	
	else   ; pr/temp
	
   ;-----------------------------
   ;
   ;   TEMPERATURE CALCS
   ;
   ;-----------------------------   
  
	  ;--- intialize hann (historical annual variable time series)
	    hann 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
            copy_VarCoords_1(tvar,hann)
   	    hann!2	= "year"
	    hann&year	= hyr

	  do iyr = 0,dimsizes(hyr)-1
            fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_historical_r*_"+hyr(iyr)+"0101*.nc")
	    a		= addfile(dir+fil,"r")
            v		= a->$var(ivar)$
	    
	  ;--- Convert K to C
            vC		= v
	    vC		= v - 273.15
	
	  ;--- calculate the annual average, rearrange dims
	    hann(:,:,iyr)	= dim_avg_n_Wrap(vC(lat|:,lon|:,time|:),2)

	  ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
            delete([/fil,a,v,vC/])

          end do   ; historical years
	  
        end if ; pr/temp 

	;--- calculate annual climatological mean, historical
   	  hann_cl(:,:,imdl)	= dim_avg_n(hann,2)
 
        ;--- calculate annual trend, historical
          hann_tr(:,:,imdl)	= regCoef_n(hyr,hann,0,2)
	   
      end do ; models 

  ;--- create .nc file of annual trends
    
    ;--- assign coordinate vars for climatological averages
      copy_VarCoords_1(tvar,hann_cl)
      hann_cl!2			= "model"
      hann_cl&model		= model
      hann_cl@short_name 	= "LOCA Climatological mean "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hann_cl@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hann_cl@units 		= "degrees C"
      
      
      
    ;--- assign coordinate vars for trends
      copy_VarCoords_1(tvar,hann_tr)
      hann_tr!2			= "model"
      hann_tr&model		= model
      hann_tr@short_name 	= "LOCA Annual "+var(ivar)+" trend, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hann_tr@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) Annual Temperature ("+var(ivar)+"), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      hann_tr@units 		= "degrees C per year"


    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")

    ;--- create file name
      fn 		= "ANN_"+var(ivar)+"_historical_"+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
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
      
      fAtt@description 	= "Annual "+var(ivar)+" climatological mean and trend calculated from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      
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
      filevardef(ncdf,"ANNCLIM",typeof(hann_cl),getvardims(hann_cl))
      filevardef(ncdf,"ANNTREND",typeof(hann_tr),getvardims(hann_tr))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      
    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"ANNCLIM",hann_cl)
      filevarattdef(ncdf,"ANNTREND",hann_tr)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model	= (/model/)
      ncdf->lat		= (/lat/)
      ncdf->lon	 	= (/lon/)
      ncdf->ANNCLIM	= (/hann_cl/)
      ncdf->ANNTREND	= (/hann_tr/)
      ncdf->ModelName	= (/modelname/)

    ;--- delete vars
      delete([/hann_cl,hann_tr/])	     
  end do  ; variables






end


 

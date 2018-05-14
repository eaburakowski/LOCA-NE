
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"


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

  ;--- Time periods
  hyr	= ispan(1980,2005,1)
  hyr!0	= "year"
  hyr&year = hyr
 
;---- Output directories

  odir	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/"	

;---- Annual Averages
;
;      Loop over models and scenarios to calculate annual average
;	tasmax, tasmin, tasmean (daily calc from tasmin and max)
;
;----

;---- Load one future file to get the lat/lon dims

  tfil 	= addfile(dir+"pr_day_ACCESS1-0_historical_r1i1p1_19810101-19811231.16th.nc","r")
  tmp1	= tfil->$var(0)$
  tvar  = tmp1(lat|:,lon|:,time|:)
  dims	= dimsizes(tvar)
  lat	= tvar&lat
  lon	= tvar&lon

printVarSummary(tvar)

;--- Loop over var, mod, scen daily files to calculate annual average
  do ivar = 0,dimsizes(var)-1
    print("----  Working on "+var(ivar)+"  -----")
	
	;--- Initialize future annual trend (1980-2005), and climatologies (early, mid, late)
	  ann_ts 	= new((/dims(0),dims(1),dimsizes(hyr),dimsizes(mdl)/),float,-9999)
	  ann_ts!0	= "lat"
	  ann_ts&lat	= lat
	  ann_ts!1	= "lon"
	  ann_ts&lon	= lon
	  ann_ts!2	= "year"
	  ann_ts&year	= hyr
	  ann_ts!3	= "model"
	  ann_ts&model	= model

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

          do iyr = 0,dimsizes(hyr)-1
            fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_historical_r*_"+hyr(iyr)+"0101*.nc")
            a 		= addfile(dir+fil,"r")
            v		= a->$var(ivar)$
	    vmm		= v
	    vmm		= v*86400

          ;--- calculate the annual average for full time series (1980-2005), rearrange dims
            ann_ts(:,:,iyr,imdl)       = dim_sum_n_Wrap(vmm(lat|:,lon|:,time|:),2)
	  
	  ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
 	    delete([/fil,a,v,vmm/])

	  end do  ; future years
	
	else   ; pr/temp
	
   ;-----------------------------
   ;
   ;   TEMPERATURE CALCS
   ;
   ;-----------------------------   

          do iyr = 0,dimsizes(hyr)-1
            fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_historical_r*_"+hyr(iyr)+"0101*.nc")
            a 		= addfile(dir+fil,"r")
            v		= a->$var(ivar)$
	    
	  ;--- convert K to C
	    vC		= v
            vC		= v-273.15
            
          ;--- calculate the annual average, rearrange dims
            ann_ts(:,:,iyr,imdl)       = dim_avg_n_Wrap(vC(lat|:,lon|:,time|:),2)	

	  ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
 	    delete([/fil,a,v,vC/])

	  end do  ; future years

        end if ; pr/temp		   
	       
       end do ; models 

  ;--- create .nc file of annual time series

    ;--- assign coordinate vars for trends
      ann_ts@short_name 	= "LOCA Annual "+var(ivar)+" time series, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      ann_ts@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) Annual ("+var(ivar)+"), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      if ivar.eq.0
        ann_ts@units 		= "mm"
      else
        ann_ts@units		= "degrees C"
      end if
 
     ;--- Filename 
      fn 		= "ANN_"+var(ivar)+"_historical_"+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+"_timeseries"
print("Filename = "+fn)
    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")
 
    ;--- remove old file
      system("/bin/rm -f "+odir+fn+".nc")
      
    ;--- create new file   
      ncdf		= addfile(odir+fn+".nc","c")	
      fAtt		= True
      fAtt@title	= fn
      fAtt@orig_fil	= var(ivar)+"_day_LOCA_CMIP5_models_historical_"+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      fAtt@Conventions	= "COARDS/CF-1.0"
      fAtt@creation_date= systemfunc("date")
      fAtt@author	= "Elizabeth Burakowski (elizabeth.burakowski@unh.edu)"
      
      fAtt@description 	= "Annual "+var(ivar)+" time series calculated from daily LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      
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
      filevardef(ncdf,"ANNTIMESERIES",typeof(ann_ts),getvardims(ann_ts))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      filevardef(ncdf,"Year",typeof(hyr),getvardims(hyr))

    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"ANNTIMESERIES",ann_ts)
      filevarattdef(ncdf,"ModelName",modelname)
      filevarattdef(ncdf,"Year",hyr)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model		= (/model/)
      ncdf->lat			= (/lat/)
      ncdf->lon	 		= (/lon/)
      ncdf->ANNTIMESERIES	= (/ann_ts/)
      ncdf->ModelName 		= (/modelname/)
      ncdf->Year		= (/hyr/)
     
    ;--- delete variables for next loop
      delete([/ann_ts/])
    
  end do  ; variables


end


 


 

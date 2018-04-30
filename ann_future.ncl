
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

  ;--- Scenario
  scn	= (/"rcp45","rcp85"/)

  ;--- Time periods
  fyr	= ispan(2006,2099,1)
  early	= ind(fyr.ge.2010 .and. fyr.le.2039)
  mid 	= ind(fyr.ge.2040 .and. fyr.le.2069)
  late  = ind(fyr.ge.2070 .and. fyr.le.2099)
 
;---- Output directories

  odir	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/"	

;---- Annual Averages
;
;      Loop over models and scenarios to calculate annual average
;	tasmax, tasmin, tasmean (daily calc from tasmin and max)
;
;----

;---- Load one future file to get the lat/lon dims

  tfil 	= addfile(dir+"pr_day_ACCESS1-0_rcp85_r1i1p1_20810101-20811231.16th.nc","r")
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
        print("----  Working on "+scn(iscn)+"  -----")
	
	;--- Initialize future annual trend (2006-2099), and climatologies (early, mid, late)
	  ann_tr 	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
printVarSummary(ann_tr)
	  ann_2010_2039 = ann_tr
	  ann_2040_2069 = ann_tr
	  ann_2070_2099 = ann_tr

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
           
          ;--- intialize fann (future annual variable time series)
	    fann	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
	    copy_VarCoords_1(tvar,fann)
	    fann!2	= "year"
	    fann&year	= fyr

          do iyr = 0,dimsizes(fyr)-1
            fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_"+scn(iscn)+"_r*_"+fyr(iyr)+"0101*.nc")
            a 		= addfile(dir+fil,"r")
            v		= a->$var(ivar)$
	    vmm		= v
	    vmm		= v*86400

          ;--- calculate the annual average for full time series (2006-2099), rearrange dims
            fann(:,:,iyr)       = dim_sum_n_Wrap(vmm(lat|:,lon|:,time|:),2)
	  
	  ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
 	    delete([/fil,a,v,vmm/])

	  end do  ; future years
	
	else   ; pr/temp
	
   ;-----------------------------
   ;
   ;   TEMPERATURE CALCS
   ;
   ;-----------------------------   

          ;--- intialize fann (future annual variable time series)
	    fann	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
	    copy_VarCoords_1(tvar,fann)
	    fann!2	= "year"
	    fann&year	= fyr

          do iyr = 0,dimsizes(fyr)-1
            fil		= systemfunc("cd "+dir+"; ls "+var(ivar)+"_day_"+mdl(imdl)+"_"+scn(iscn)+"_r*_"+fyr(iyr)+"0101*.nc")
            a 		= addfile(dir+fil,"r")
            v		= a->$var(ivar)$
	    
	  ;--- convert K to C
	    vC		= v
            vC		= v-273.15
            
          ;--- calculate the annual average, rearrange dims
            fann(:,:,iyr)       = dim_avg_n_Wrap(vC(lat|:,lon|:,time|:),2)	

	  ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
 	    delete([/fil,a,v,vC/])

	  end do  ; future years

        end if ; pr/temp		   
      
        ;-- Calculate annual trend in future
        
	;--- calculate annual trend (linear) for imdl
	    ann_tr(:,:,imdl)	= regCoef_n(fyr,fann,0,2)
	    
	;--- calculate early (2010-2039) climatology
 	    ann_2010_2039(:,:,imdl) = dim_avg_n_Wrap(fann(:,:,early),2)
	    ann_2040_2069(:,:,imdl) = dim_avg_n_Wrap(fann(:,:,mid),2)
	    ann_2070_2099(:,:,imdl) = dim_avg_n_Wrap(fann(:,:,late),2)	

        ;--- delete fann for next loop
            delete([/fann/])
	       
       end do ; models 

  ;--- create .nc file of annual trends and climatology (early, mid, late)

    ;--- assign coordinate vars for trends
      copy_VarCoords_1(tvar,ann_tr)
      ann_tr!2			= "model"
      ann_tr&model		= model
      ann_tr@short_name 	= "LOCA Annual "+var(ivar)+" trend, "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      ann_tr@long_name		= "Statistically downscaled Localized Constructed Analogs (LOCA) Annual ("+var(ivar)+"), "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      ann_tr@units 		= "degrees C per year"
     
    ;--- assign coordinate vars for climatology (early, mid, late)
     ;--- Early (2010-2039) 
      copy_VarCoords_1(tvar,ann_2010_2039)
      ann_2010_2039!2		= "model"
      ann_2010_2039&model       = model
      ann_2010_2039@short_name  = "LOCA Annual "+var(ivar)+" Climatology, 2010-2039"
      ann_2010_2039@long_name   = "Statistically downscaled Localized Constructed Analogs (LOCA) Annual ("+var(ivar)+"), 2010-2039"
      ann_2010_2039@units       = "degrees C"
      
     ;--- Mid (2040-2069) 
      copy_VarCoords_1(tvar,ann_2040_2069)
      ann_2040_2069!2			= "model"
      ann_2040_2069&model             	= model
      ann_2040_2069@short_name       	= "LOCA Annual "+var(ivar)+" Climatology, 2040-2069"
      ann_2040_2069@long_name         	= "Statistically downscaled Localized Constructed Analogs (LOCA) Annual ("+var(ivar)+"), 2040-2069"
      ann_2040_2069@units             	= "degrees C"
      
     ;-- Late (2070-2099) 
      copy_VarCoords_1(tvar,ann_2070_2099)
      ann_2070_2099!2		= "model"
      ann_2070_2099&model	= model
      ann_2070_2099@short_name  = "LOCA Annual "+var(ivar)+" Climatology, 2070-2099"
      ann_2070_2099@long_name   = "Statistically downscaled Localized Constructed Analogs (LOCA) Annual ("+var(ivar)+"), 2070-2099"
      ann_2070_2099@units       = "degrees C"
    
     ;--- Filename 
      fn 		= "ANN_"+var(ivar)+"_"+scn(iscn)+"_"+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+"_trend_clim"
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
      fAtt@orig_fil	= var(ivar)+"_day_LOCA_CMIP5_models_"+scn(iscn)+"_"+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fAtt@Conventions	= "COARDS/CF-1.0"
      fAtt@creation_date= systemfunc("date")
      fAtt@author	= "Elizabeth Burakowski (elizabeth.burakowski@unh.edu)"
      
      fAtt@description 	= "Annual "+var(ivar)+" trend and climatology (early, mid, late) calculated from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      
    ;--- file attributes, size
      fileattdef(ncdf,fAtt)
      dimNames		= (/"lat", "lon", "model"/)
      dimSizes		= (/ dims(0), dims(1), dimsizes(model) /)
      dimUnlim          = (/ False, False, False /)
      filedimdef(ncdf,dimNames,dimSizes,dimUnlim)
print(typeof(model))
print(getvardims(model))

    ;--- Define file variables
      filevardef(ncdf,"model",typeof(model),getvardims(model))
      filevardef(ncdf,"lat",typeof(lat),getvardims(lat))
      filevardef(ncdf,"lon",typeof(lon),getvardims(lon))
      filevardef(ncdf,"ANNTREND",typeof(ann_tr),getvardims(ann_tr))
      filevardef(ncdf,"ANNCLIM_EARLY",typeof(ann_2010_2039),getvardims(ann_2010_2039))
      filevardef(ncdf,"ANNCLIM_MID",typeof(ann_2040_2069),getvardims(ann_2040_2069))
      filevardef(ncdf,"ANNCLIM_LATE",typeof(ann_2070_2099),getvardims(ann_2070_2099))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))

    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"ANNTREND",ann_tr)
      filevarattdef(ncdf,"ANNCLIM_EARLY",ann_2010_2039)
      filevarattdef(ncdf,"ANNCLIM_MID",ann_2040_2069)
      filevarattdef(ncdf,"ANNCLIM_LATE",ann_2070_2099)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model	= (/model/)
      ncdf->lat		= (/lat/)
      ncdf->lon	 	= (/lon/)
      ncdf->ANNTREND	= (/ann_tr/)
      ncdf->ANNCLIM_EARLY = (/ann_2010_2039/)
      ncdf->ANNCLIM_MID	= (/ann_2040_2069/)
      ncdf->ANNCLIM_LATE= (/ann_2070_2099/)
      ncdf->ModelName 	= (/modelname/)
     
    ;--- delete variables for next loop
      delete([/ann_tr,ann_2010_2039,ann_2040_2069,ann_2070_2099/])
    
    end do ; scenario	     
  end do  ; variables


end


 


 

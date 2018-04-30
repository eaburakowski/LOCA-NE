load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; daysExtremeTemp_hist.ncl 
; EA Burakowski
; 2018-04-26
; 
; daysExtremeTemp_hist.ncl calculates the climatological average
;	days per year:  
;	- tmax > 90F (32.2C) 
;	- tmax > 95F (35C)
;	- tmin < 32F (0C)
;	- tmin < 0F (-17.8C)
;	- tmax on hottest day of year 
;	- tmin on coldest day of year
;	
;	... from CMIP5 LOCAlized Constructed Analogs (LOCA) daily 
;	downscaled data (1980-2005).  
;	Inputs include daily LOCA data for variables:
;		(1) tasmax (surface maximum air temperature)
;	 	(2) taxmin (surface minimum air temperature)
;	for 29 CMIP5 models (see below, mdl), 1980-2005
;
; Output includes:
;   .nc file for 30-yr climatological mean and annual trend:
;	- ExtremeTemp_historical_1980-2005.nc
;
; Use plotExtremeTemp.ncl and plotExtremeTemp_trends.ncl to plot up the .nc files
;	

;======================================================================
;                         The main code
;======================================================================

begin

;---- LOCA directory

  dir	= "/net/nfs/merrimack/raid/Northeast_US_Downscaling_cmip5/"
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


;---- Annual Calculations
;
;      Loop over models and years to calculate annual 
;	extreme temperature metrics
;
;----

;---- Load one historical file to get the lat/lon dims

  tfil 	= addfile(dir+"tasmin_day_ACCESS1-0_historical_r1i1p1_19800101-19801231.16th.nc","r")
  tmp1	= tfil->tasmin
  tvar  = tmp1(lat|:,lon|:,time|:)
  dims	= dimsizes(tvar)
  lat	= tvar&lat
  lon	= tvar&lon

printVarSummary(tvar)

;--- Loop model, scenario daily files to calculate annual extreme temp metrics

;--- Initialize historical climatology and trend (reuse first variable)
days90F_cl	= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
days90F_tr 	= days90F_cl
days95F_cl	= days90F_cl
days95F_tr	= days90F_cl
days32F_cl	= days90F_cl
days32F_tr	= days90F_cl
days0F_cl	= days90F_cl
days0F_tr	= days90F_cl
txHottest_cl	= days90F_cl
txHottest_tr	= days90F_cl
txColdest_cl	= days90F_cl
txColdest_tr	= days90F_cl
tnColdest_cl	= days90F_cl
tnColdest_tr	= days90F_cl
tnHottest_cl	= days90F_cl
tnHottest_tr	= days90F_cl

  ;--- Loop over models
  do imdl = 0,dimsizes(mdl)-28
    print("----  Working on "+mdl(imdl)+"  -----")  
  
  ;--- intialize hann (historical annual variable time series)
    days90F 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
    copy_VarCoords_1(tvar,days90F)
    days90F!2	= "year"
    days90F&year= hyr

    days95F 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
    copy_VarCoords_1(tvar,days95F)
    days95F!2	= "year"
    days95F&year= hyr

    days32F 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
    copy_VarCoords_1(tvar,days32F)
    days32F!2	= "year"
    days32F&year= hyr

    days0F 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
    copy_VarCoords_1(tvar,days0F)
    days0F!2	= "year"
    days0F&year= hyr

    txHottest 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
    copy_VarCoords_1(tvar,txHottest)
    txHottest!2	= "year"
    txHottest&year= hyr

    tnColdest 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
    copy_VarCoords_1(tvar,tnColdest)
    tnColdest!2	= "year"
    tnColdest&year= hyr
    
    tnHottest 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
    copy_VarCoords_1(tvar,tnHottest)
    tnHottest!2	= "year"
    tnHottest&year= hyr

    txColdest 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
    copy_VarCoords_1(tvar,txColdest)
    txColdest!2	= "year"
    txColdest&year= hyr

    ;--- Loop over years
    do iyr = 0,dimsizes(hyr)-1
      tn_fil	= systemfunc("cd "+dir+"; ls tasmin_day_"+mdl(imdl)+"_historical_r*_"+hyr(iyr)+"0101*.nc")
      a		= addfile(dir+tn_fil,"r")
      tn	= a->tasmin

      tx_fil	= systemfunc("cd "+dir+"; ls tasmax_day_"+mdl(imdl)+"_historical_r*_"+hyr(iyr)+"0101*.nc")
      b		= addfile(dir+tx_fil,"r")
      tx	= b->tasmax

      ;--- Convert K to C
      txC	= tx
      txC	= tx - 273.15
      tnC	= tn
      tnC	= tn - 273.15
      
      ;--- Calculate number of days above and below thresholds
      days90F(:,:,iyr)	= dim_num_n(tx.ge.32.2,0)
      days95F(:,:,iyr)	= dim_num_n(tx.ge.35,0)
      days32F(:,:,iyr)	= dim_num_n(tn.le.0,0)
      days0F(:,:,iyr)	= dim_num_n(tn.le.-17.8,0)
      
      ;--- Calculate hottest and coldest day of year
      txHottest(:,:,iyr)= dim_max_n(tx,0)
      txColdest(:,:,iyr)= dim_min_n(tx,0)
      
      ;--- Calculate hottest and coldest night of the year
      tnHottest(:,:,iyr)= dim_max_n(tn,0)
      txColdest(:,:,iyr)= dim_min_n(tn,0)

      ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
      delete([/tx_fil,a,tx,tn_fil,b,tn,txC,tnC/])

    end do   ; historical years

    ;--- calculate annual climatological mean, historical
    days90F_cl(:,:,imdl)		= dim_avg_n(days90F,2)
    days95F_cl(:,:,imdl)		= dim_avg_n(days95F,2)
    days32F_cl(:,:,imdl)		= dim_avg_n(days32F,2)
    days0F_cl(:,:,imdl)		= dim_avg_n(days0F,2)
    txHottest_cl(:,:,imdl)	= dim_avg_n(txHottest,2)
    txColdest_cl(:,:,imdl)	= dim_avg_n(txColdest,2)
    tnHottest_cl(:,:,imdl)	= dim_avg_n(tnHottest,2)
    tnColdest_cl(:,:,imdl)	= dim_avg_n(tnColdest,2)

    ;--- calculate trend, historical
    days90F_tr(:,:,imdl)		= regCoef_n(hyr,days90F,0,2)
    days95F_tr(:,:,imdl)		= regCoef_n(hyr,days95F,0,2)
    days32F_tr(:,:,imdl)		= regCoef_n(hyr,days32F,0,2)
    days0F_tr(:,:,imdl)		= regCoef_n(hyr,days0F,0,2)
    txHottest_tr(:,:,imdl)	= regCoef_n(hyr,txHottest,0,2)
    txColdest_tr(:,:,imdl)	= regCoef_n(hyr,txColdest,0,2)
    tnHottest_tr(:,:,imdl)	= regCoef_n(hyr,tnHottest,0,2)
    tnColdest_tr(:,:,imdl)	= regCoef_n(hyr,tnColdest,0,2)
  
  end do ; models 

  ;--- create .nc file of annual trends
    
    ;--- assign coordinate vars for climatological averages
      copy_VarCoords_1(tvar,days90F_cl)
      days90F_cl!2		= "model"
      days90F_cl&model		= model
      days90F_cl@short_name 	= "LOCA Climatological mean days per year > 90F (32.2C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days90F_cl@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 90F (32.2C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days90F_cl@units 		= "days" 
      
      copy_VarCoords_1(tvar,days95F_cl)
      days95F_cl!2		= "model"
      days95F_cl&model		= model
      days95F_cl@short_name 	= "LOCA Climatological mean days per year > 95F (35C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days95F_cl@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 95F (35C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days95F_cl@units 		= "days"
      
      copy_VarCoords_1(tvar,days32F_cl)
      days32F_cl!2		= "model"
      days32F_cl&model		= model
      days32F_cl@short_name 	= "LOCA Climatological mean days per year < 32F (0C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days32F_cl@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year < 32F (0C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days32F_cl@units 		= "days"
      
      copy_VarCoords_1(tvar,days0F_cl)
      days0F_cl!2		= "model"
      days0F_cl&model		= model
      days0F_cl@short_name 	= "LOCA Climatological mean days per year < 0F (-17.8C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days0F_cl@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year < 0F (-17.8C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days0F_cl@units 		= "days"
      
      copy_VarCoords(tvar,txHottest_cl)
      txHottest_cl!2		= "model"
      txHottest_cl&model	= model
      txHottest_cl@short_name	= "LOCA Climatological mean hottest day of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      txHottest_cl@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean hottest day of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      txHottest_cl@units 	= "degrees C"
      
      copy_VarCoords(tvar,txColdest_cl)
      txColdest_cl!2		= "model"
      txColdest_cl&model	= model
      txColdest_cl@short_name	= "LOCA Climatological mean coldest day of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      txColdest_cl@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean coldest day of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      txColdest_cl@units 	= "degrees C"
      
      copy_VarCoords(tvar,tnColdest_cl)
      tnColdest_cl!2		= "model"
      tnColdest_cl&model	= model
      tnColdest_cl@short_name	= "LOCA Climatological mean coldest night of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      tnColdest_cl@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean coldest night of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      tnColdest_cl@units 	= "degrees C"
      
      copy_VarCoords(tvar,tnHottest_cl)
      tnHottest_cl!2		= "model"
      tnHottest_cl&model	= model
      tnHottest_cl@short_name	= "LOCA Climatological mean Hottest night of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      tnHottest_cl@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean hottest night of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      tnHottest_cl@units 	= "degrees C"
      
      
    ;--- assign coordinate vars for historical trends
      copy_VarCoords_1(tvar,days90F_tr)
      days90F_tr!2		= "model"
      days90F_tr&model		= model
      days90F_tr@short_name 	= "LOCA trend days per year > 90F (32.2C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days90F_tr@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) trend days per year > 90F (32.2C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days90F_tr@units 		= "days" 
      
      copy_VarCoords_1(tvar,days95F_tr)
      days95F_tr!2		= "model"
      days95F_tr&model		= model
      days95F_tr@short_name 	= "LOCA trend days per year > 95F (35C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days95F_tr@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) trend days per year > 95F (35C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days95F_tr@units 		= "days"
      
      copy_VarCoords_1(tvar,days32F_tr)
      days32F_tr!2		= "model"
      days32F_tr&model		= model
      days32F_tr@short_name 	= "LOCA trend days per year < 32F (0C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days32F_tr@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) trend days per year < 32F (0C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days32F_tr@units 		= "days"
      
      copy_VarCoords_1(tvar,days0F_tr)
      days0F_tr!2		= "model"
      days0F_tr&model		= model
      days0F_tr@short_name 	= "LOCA trend days per year < 0F (-17.8C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days0F_tr@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) trend days per year < 0F (-17.8C), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      days0F_tr@units 		= "days"
      
      copy_VarCoords(tvar,txHottest_tr)
      txHottest_tr!2		= "model"
      txHottest_tr&model	= model
      txHottest_tr@short_name= "LOCA trend hottest day of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      txHottest_tr@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) trend hottest day of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      txHottest_tr@units 	= "degrees C"
      
      copy_VarCoords(tvar,txColdest_tr)
      txColdest_tr!2		= "model"
      txColdest_tr&model	= model
      txColdest_tr@short_name= "LOCA trend coldest day of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      txColdest_tr@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) trend coldest day of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      txColdest_tr@units 	= "degrees C"
      
      copy_VarCoords(tvar,tnColdest_tr)
      tnColdest_tr!2		= "model"
      tnColdest_tr&model	= model
      tnColdest_tr@short_name= "LOCA trend coldest night of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      tnColdest_tr@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) trend coldest night of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      tnColdest_tr@units 	= "degrees C"
      
      copy_VarCoords(tvar,tnHottest_tr)
      tnHottest_tr!2		= "model"
      tnHottest_tr&model	= model
      tnHottest_tr@short_name= "LOCA trend Hottest night of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      tnHottest_tr@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) trend hottest night of year, "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      tnHottest_tr@units 	= "degrees C"
      
      
    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")

    ;--- create file name
      fn 		= "ExtremeTemp_historical_"+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
print("Filename = "+fn)

    ;--- remove old file, if it exists
      system("/bin/rm -f "+odir+fn+".nc")
      
    ;--- create new file   
      ncdf		= addfile(odir+fn+".nc","c")	
      fAtt		= True
      fAtt@title	= fn
      fAtt@orig_fil	= "tasmax_day_LOCA_CMIP5_models_historical_"+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      fAtt@Conventions	= "COARDS/CF-1.0"
      fAtt@creation_date= systemfunc("date")
      fAtt@author	= "Elizabeth Burakowski (elizabeth.burakowski@unh.edu)"
      
      fAtt@description 	= "Extreme Temperature climatological means and trends calculated from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
      
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
      filevardef(ncdf,"DAYS90FCLIM",typeof(days90F_cl),getvardims(days90F_cl))
      filevardef(ncdf,"DAYS95FCLIM",typeof(days95F_cl),getvardims(days95F_cl))
      filevardef(ncdf,"DAYS32FCLIM",typeof(days32F_cl),getvardims(days32F_cl))
      filevardef(ncdf,"DAYS0FCLIM",typeof(days0F_cl),getvardims(days0F_cl))
      filevardef(ncdf,"HOTTESTDAYCLIM",typeof(txHottest_cl),getvardims(txHottest_cl))
      filevardef(ncdf,"COLDESTDAYCLIM",typeof(txColdest_cl),getvardims(txColdest_cl))
      filevardef(ncdf,"HOTTESTNIGHTCLIM",typeof(tnHottest_cl),getvardims(tnHottest_cl))
      filevardef(ncdf,"COLDESTNIGHTCLIM",typeof(tnColdest_cl),getvardims(tnColdest_cl))
      filevardef(ncdf,"DAYS90FTREND",typeof(days90F_tr),getvardims(days90F_tr))
      filevardef(ncdf,"DAYS95FTREND",typeof(days95F_tr),getvardims(days95F_tr))
      filevardef(ncdf,"DAYS32FTREND",typeof(days32F_tr),getvardims(days32F_tr))
      filevardef(ncdf,"DAYS0FTREND",typeof(days0F_tr),getvardims(days0F_tr))
      filevardef(ncdf,"HOTTESTDAYTREND",typeof(txHottest_tr),getvardims(txHottest_tr))
      filevardef(ncdf,"COLDESTDAYTREND",typeof(txColdest_tr),getvardims(txColdest_tr))
      filevardef(ncdf,"HOTTESTNIGHTTREND",typeof(tnHottest_tr),getvardims(tnHottest_tr))
      filevardef(ncdf,"COLDESTNIGHTTREND",typeof(tnColdest_tr),getvardims(tnColdest_tr))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      
    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"DAYS90FCLIM",days90F_cl)
      filevarattdef(ncdf,"DAYS95FCLIM",days95F_cl)
      filevarattdef(ncdf,"DAYS32FCLIM",days32F_cl)
      filevarattdef(ncdf,"DAYS0FCLIM",days0F_cl)
      filevarattdef(ncdf,"HOTTESTDAYCLIM",txHottest_cl)
      filevarattdef(ncdf,"COLDESTDAYCLIM",txColdest_cl)
      filevarattdef(ncdf,"HOTTESTNIGHTCLIM",tnHottest_cl)
      filevarattdef(ncdf,"COLDESTNIGHTCLIM",tnColdest_cl)
      filevarattdef(ncdf,"DAYS90FTREND",days90F_tr)
      filevarattdef(ncdf,"DAYS95FTREND",days95F_tr)
      filevarattdef(ncdf,"DAYS32FTREND",days32F_tr)
      filevarattdef(ncdf,"DAYS0FTREND",days0F_tr)
      filevarattdef(ncdf,"HOTTESTDAYTREND",txHottest_tr)
      filevarattdef(ncdf,"COLDESTDAYTREND",txColdest_tr)
      filevarattdef(ncdf,"HOTTESTNIGHTTREND",tnHottest_tr)
      filevarattdef(ncdf,"COLDESTNIGHTTREND",tnColdest_tr)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model		= (/model/)
      ncdf->lat			= (/lat/)
      ncdf->lon	 		= (/lon/)
      ncdf->DAYS95FCLIM		= (/days95F_cl/)
      ncdf->DAYS90FCLIM		= (/days90F_cl/)
      ncdf->DAYS32FCLIM		= (/days32F_cl/)
      ncdf->DAYS0FCLIM		= (/days0F_cl/)
      ncdf->HOTTESTDAYCLIM	= (/txHottest_cl/)
      ncdf->COLDESTDAYCLIM	= (/txColdest_cl/)
      ncdf->HOTTESTNIGHTCLIM	= (/tnHottest_cl/)
      ncdf->COLDESTNIGHTCLIM	= (/tnColdest_cl/)
      ncdf->DAYS95FTREND	= (/days95F_tr/)
      ncdf->DAYS90FTREND	= (/days90F_tr/)
      ncdf->DAYS32FTREND	= (/days32F_tr/)
      ncdf->DAYS0FTREND		= (/days0F_tr/)
      ncdf->HOTTESTDAYTREND	= (/txHottest_tr/)
      ncdf->COLDESTDAYTREND	= (/txColdest_tr/)
      ncdf->HOTTESTNIGHTTREND	= (/tnHottest_tr/)
      ncdf->COLDESTNIGHTTREND	= (/tnColdest_tr/)
      ncdf->ModelName	= (/modelname/)
	     
end


 

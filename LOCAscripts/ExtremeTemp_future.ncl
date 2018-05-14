load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; daysExtremeTemp_future.ncl 
; EA Burakowski
; 2018-04-30
; 
; daysExtremeTemp_future.ncl calculates the climatological average
;	days per year:  
;	- tmax > 90F (32.2C) 
;	- tmax > 95F (35C)
;	- tmin < 32F (0C)
;	- tmin < 0F (-17.8C)
;	- tmax on hottest day of year 
;	- tmin on coldest day of year
;	
;	... from CMIP5 LOCAlized Constructed Analogs (LOCA) daily 
;	downscaled data (2005-2099).  
;	Inputs include daily LOCA data for variables:
;		(1) tasmax (surface maximum air temperature)
;	 	(2) taxmin (surface minimum air temperature)
;	for 29 CMIP5 models (see below, mdl), 2006-2099
;
; Output includes:
;   .nc file for 30-yr climatological mean and annual trend:
;	- ExtremeTemp_"+scn(iscn)+"2006-2099.nc
;
; Use plotExtremeTemp_future.ncl to plot up the .nc files
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

;---- Scenarios
  scn	= (/"rcp45","rcp85"/)

;---- Time Periods  
  fyr   = ispan(2006,2099,1)
  early	= ind(fyr.ge.2010 .and. fyr.le.2039)
  mid 	= ind(fyr.ge.2040 .and. fyr.le.2069)
  late	= ind(fyr.ge.2070 .and. fyr.le.2099)

;---- Output directories

  odir	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/"	


;---- Annual Calculations
;
;      Loop over models and years to calculate annual 
;	extreme temperature metrics
;
;----

;---- Load one future file to get the lat/lon dims

  tfil 	= addfile(dir+"tasmin_day_ACCESS1-0_rcp45_r1i1p1_20810101-20811231.16th.nc","r")
  tmp1	= tfil->tasmin
  tvar  = tmp1(lat|:,lon|:,time|:)
  dims	= dimsizes(tvar)
  lat	= tvar&lat
  lon	= tvar&lon

printVarSummary(tvar)

;--- Loop model, scenario daily files to calculate annual extreme temp metrics
do iscn = 0,dimsizes(iscn)
;--- Initialize future climatology (reuse first variable)
days90F_early		= new((/dims(0),dims(1),dimsizes(mdl)/),float,-9999)
days90F_early!0		= "lat"
days90F_early&lat	= lat
days90F_early!1		= "lon"
days90F_early&lon	= lon

;-- Reuse first variable (days90F) format
days90F_mid		= days90F_early
days90F_late		= days90F_early
days95F_early		= days90F_early
days95F_mid		= days90F_early
days95F_late		= days90F_early
days32F_early		= days90F_early
days32F_mid		= days90F_early
days32F_late		= days90F_early
days0F_early		= days90F_early
days0F_mid		= days90F_early
days0F_late		= days90F_early
txHottest_early		= days90F_early
txHottest_mid		= days90F_early
txHottest_late		= days90F_early
txColdest_early		= days90F_early
txColdest_mid		= days90F_early
txColdest_late		= days90F_early
tnHottest_early		= days90F_early
tnHottest_mid		= days90F_early
tnHottest_late		= days90F_early
tnColdest_early		= days90F_early
tnColdest_mid		= days90F_early
tnColdest_late		= days90F_early

  ;--- Loop over models
  do imdl = 0,dimsizes(mdl)-1
    print("----  Working on "+mdl(imdl)+"  -----")  
  
  ;--- intialize hann (future variable time series)
    days90F 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
    copy_VarCoords_1(tvar,days90F)
    days90F!2	= "year"
    days90F&year= fyr

    days95F 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
    copy_VarCoords_1(tvar,days95F)
    days95F!2	= "year"
    days95F&year= fyr

    days32F 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
    copy_VarCoords_1(tvar,days32F)
    days32F!2	= "year"
    days32F&year= fyr

    days0F 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
    copy_VarCoords_1(tvar,days0F)
    days0F!2	= "year"
    days0F&year= fyr

    txHottest 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
    copy_VarCoords_1(tvar,txHottest)
    txHottest!2	= "year"
    txHottest&year= fyr

    tnColdest 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
    copy_VarCoords_1(tvar,tnColdest)
    tnColdest!2	= "year"
    tnColdest&year= fyr
    
    tnHottest 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
    copy_VarCoords_1(tvar,tnHottest)
    tnHottest!2	= "year"
    tnHottest&year= fyr

    txColdest 	= new((/dims(0),dims(1),dimsizes(fyr)/),float,-9999)
    copy_VarCoords_1(tvar,txColdest)
    txColdest!2	= "year"
    txColdest&year= fyr

    ;--- Loop over years
    do iyr = 0,dimsizes(fyr)-1
      tn_fil	= systemfunc("cd "+dir+"; ls tasmin_day_"+mdl(imdl)+"_"+scn(iscn)+"*_"+fyr(iyr)+"0101*.nc")
      a		= addfile(dir+tn_fil,"r")
      tn	= a->tasmin

      tx_fil	= systemfunc("cd "+dir+"; ls tasmax_day_"+mdl(imdl)+"_"+scn(iscn)+"*_"+fyr(iyr)+"0101*.nc")
      b		= addfile(dir+tx_fil,"r")
      tx	= b->tasmax

      ;--- Convert K to C
      txC	= tx
      txC	= tx - 273.15
      tnC	= tn
      tnC	= tn - 273.15
      
      ;--- Calculate number of days above and below thresholds
      days90F(:,:,iyr)	= dim_num_n(txC.ge.32.2,0)
      days95F(:,:,iyr)	= dim_num_n(txC.ge.35,0)
      days32F(:,:,iyr)	= dim_num_n(tnC.le.0,0)
      days0F(:,:,iyr)	= dim_num_n(tnC.le.-17.8,0)
      
      ;--- Calculate hottest and coldest day of year
      txHottest(:,:,iyr)= dim_max_n(txC,0)
      txColdest(:,:,iyr)= dim_min_n(txC,0)
      
      ;--- Calculate hottest and coldest night of the year
      tnHottest(:,:,iyr)= dim_max_n(tnC,0)
      tnColdest(:,:,iyr)= dim_min_n(tnC,0)

      ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
      delete([/tx_fil,a,tx,tn_fil,b,tn,txC,tnC/])

    end do   ; historical years

    ;--- calculate annual climatological mean, future early, mid and late century
    days90F_early(:,:,imdl)	= dim_avg_n(days90F(:,:,early),2)
    days90F_mid(:,:,imdl)	= dim_avg_n(days90F(:,:,mid),2)
    days90F_late(:,:,imdl)	= dim_avg_n(days90F(:,:,late),2)
    days95F_early(:,:,imdl)	= dim_avg_n(days95F(:,:,early),2)
    days95F_mid(:,:,imdl)	= dim_avg_n(days95F(:,:,mid),2)
    days95F_late(:,:,imdl)	= dim_avg_n(days95F(:,:,late),2)
    days32F_early(:,:,imdl)	= dim_avg_n(days32F(:,:,early),2)
    days32F_mid(:,:,imdl)	= dim_avg_n(days32F(:,:,mid),2)
    days32F_late(:,:,imdl)	= dim_avg_n(days32F(:,:,late),2)
    days0F_early(:,:,imdl)	= dim_avg_n(days0F(:,:,early),2)
    days0F_mid(:,:,imdl)	= dim_avg_n(days0F(:,:,mid),2)
    days0F_late(:,:,imdl)	= dim_avg_n(days0F(:,:,late),2)
    txHottest_early(:,:,imdl)	= dim_avg_n(txHottest(:,:,early),2)
    txHottest_mid(:,:,imdl)	= dim_avg_n(txHottest(:,:,mid),2)
    txHottest_late(:,:,imdl)	= dim_avg_n(txHottest(:,:,late),2)
    txColdest_early(:,:,imdl)	= dim_avg_n(txColdest(:,:,early),2)
    txColdest_mid(:,:,imdl)	= dim_avg_n(txColdest(:,:,mid),2)
    txColdest_late(:,:,imdl)	= dim_avg_n(txColdest(:,:,late),2)
    tnHottest_early(:,:,imdl)	= dim_avg_n(tnHottest(:,:,early),2)
    tnHottest_mid(:,:,imdl)	= dim_avg_n(tnHottest(:,:,mid),2)
    tnHottest_late(:,:,imdl)	= dim_avg_n(tnHottest(:,:,late),2)
    tnColdest_early(:,:,imdl)	= dim_avg_n(tnColdest(:,:,early),2)
    tnColdest_mid(:,:,imdl)	= dim_avg_n(tnColdest(:,:,mid),2)
    tnColdest_late(:,:,imdl)	= dim_avg_n(tnColdest(:,:,late),2)
 print("Max TnColdest_early = "+max(tnColdest_early)) 
 print("Min TnColdest_early = "+min(tnColdest_early))
 print("Max TxHottest_late = "+max(txHottest_late))

  end do ; models 

  ;--- create .nc file of annual trends
    
    ;--- assign coordinate vars for climatological averages
      copy_VarCoords_1(tvar,days90F_early)
      days90F_early!2		= "model"
      days90F_early&model	= model
      days90F_early@short_name 	= "LOCA Climatological mean days per year > 90F (32.2C), 2010-2039"
      days90F_early@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 90F (32.2C), 2010-2039"
      days90F_early@units 	= "days" 
      
      copy_VarCoords_1(tvar,days90F_mid)
      days90F_mid!2		= "model"
      days90F_mid&model		= model
      days90F_mid@short_name 	= "LOCA Climatological mean days per year > 90F (32.2C), 2040-2069"
      days90F_mid@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 90F (32.2C), 2040-2069"
      days90F_mid@units 	= "days" 
      
      copy_VarCoords_1(tvar,days90F_late)
      days90F_late!2		= "model"
      days90F_late&model	= model
      days90F_late@short_name 	= "LOCA Climatological mean days per year > 90F (32.2C), 2070-2099"
      days90F_late@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 90F (32.2C), 2070-2099"
      days90F_late@units 	= "days" 
      
      copy_VarCoords_1(tvar,days95F_early)
      days95F_early!2		= "model"
      days95F_early&model	= model
      days95F_early@short_name 	= "LOCA Climatological mean days per year > 95F (35C), 2010-2039"
      days95F_early@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 95F (32.2C), 2010-2039"
      days95F_early@units 	= "days" 
      
      copy_VarCoords_1(tvar,days95F_mid)
      days95F_mid!2		= "model"
      days95F_mid&model		= model
      days95F_mid@short_name 	= "LOCA Climatological mean days per year > 95F (35C), 2040-2069"
      days95F_mid@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 95F (32.2C), 2040-2069"
      days95F_mid@units 	= "days" 
      
      copy_VarCoords_1(tvar,days95F_late)
      days95F_late!2		= "model"
      days95F_late&model	= model
      days95F_late@short_name 	= "LOCA Climatological mean days per year > 95F (35C), 2070-2099"
      days95F_late@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 95F (32.2C), 2070-2099"
      days95F_late@units 	= "days" 
      
      copy_VarCoords_1(tvar,days32F_early)
      days32F_early!2		= "model"
      days32F_early&model	= model
      days32F_early@short_name 	= "LOCA Climatological mean days per year < 32F (0C), 2010-2039"
      days32F_early@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 32F (32.2C), 2010-2039"
      days32F_early@units 	= "days" 
      
      copy_VarCoords_1(tvar,days32F_mid)
      days32F_mid!2		= "model"
      days32F_mid&model		= model
      days32F_mid@short_name 	= "LOCA Climatological mean days per year < 32F (0C), 2040-2069"
      days32F_mid@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 32F (32.2C), 2040-2069"
      days32F_mid@units 	= "days" 
      
      copy_VarCoords_1(tvar,days32F_late)
      days32F_late!2		= "model"
      days32F_late&model	= model
      days32F_late@short_name 	= "LOCA Climatological mean days per year < 32F (0C), 2070-2099"
      days32F_late@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 32F (32.2C), 2070-2099"
      days32F_late@units 	= "days" 
      
      copy_VarCoords_1(tvar,days0F_early)
      days0F_early!2		= "model"
      days0F_early&model	= model
      days0F_early@short_name 	= "LOCA Climatological mean days per year < 0F (-17.8C), 2010-2039"
      days0F_early@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 0F (32.2C), 2010-2039"
      days0F_early@units 	= "days" 
      
      copy_VarCoords_1(tvar,days0F_mid)
      days0F_mid!2		= "model"
      days0F_mid&model		= model
      days0F_mid@short_name 	= "LOCA Climatological mean days per year < 0F (-17.8C), 2040-2069"
      days0F_mid@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 0F (32.2C), 2040-2069"
      days0F_mid@units 	= "days" 
      
      copy_VarCoords_1(tvar,days0F_late)
      days0F_late!2		= "model"
      days0F_late&model		= model
      days0F_late@short_name 	= "LOCA Climatological mean days per year < 0F (-17.8C), 2070-2099"
      days0F_late@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual days per year > 0F (32.2C), 2070-2099"
      days0F_late@units 	= "days" 
      
      copy_VarCoords_1(tvar,txHottest_early)
      txHottest_early!2		= "model"
      txHottest_early&model	= model
      txHottest_early@short_name= "LOCA Climatological mean hottest day of year, 2010-2039"
      txHottest_early@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual hottest day of year, 2010-2039"
      txHottest_early@units 	= "degrees C" 
      
      copy_VarCoords_1(tvar,txHottest_mid)
      txHottest_mid!2		= "model"
      txHottest_mid&model	= model
      txHottest_mid@short_name 	= "LOCA Climatological mean hottest day of year, 2040-2069"
      txHottest_mid@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual hottest day of year, 2040-2069"
      txHottest_mid@units 	= "degrees C" 
      
      copy_VarCoords_1(tvar,txHottest_late)
      txHottest_late!2		= "model"
      txHottest_late&model	= model
      txHottest_late@short_name = "LOCA Climatological mean hottest day of year, 2070-2099"
      txHottest_late@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual hottest day of year, 2070-2099"
      txHottest_late@units 	= "degrees C" 
      
      copy_VarCoords_1(tvar,txColdest_early)
      txColdest_early!2		= "model"
      txColdest_early&model	= model
      txColdest_early@short_name= "LOCA Climatological mean coldest day of year, 2010-2039"
      txColdest_early@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual coldest day of year, 2010-2039"
      txColdest_early@units 	= "degrees C" 
      
      copy_VarCoords_1(tvar,txColdest_mid)
      txColdest_mid!2		= "model"
      txColdest_mid&model	= model
      txColdest_mid@short_name 	= "LOCA Climatological mean coldest day of year, 2040-2069"
      txColdest_mid@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual coldest day of year, 2040-2069"
      txColdest_mid@units 	= "degrees C" 
      
      copy_VarCoords_1(tvar,txColdest_late)
      txColdest_late!2		= "model"
      txColdest_late&model	= model
      txColdest_late@short_name = "LOCA Climatological mean coldest day of year, 2070-2099"
      txColdest_late@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual coldest day of year, 2070-2099"
      txColdest_late@units 	= "degrees C"
      
      copy_VarCoords_1(tvar,tnColdest_early)
      tnColdest_early!2		= "model"
      tnColdest_early&model	= model
      tnColdest_early@short_name= "LOCA Climatological mean coldest night of year, 2010-2039"
      tnColdest_early@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual coldest night of year, 2010-2039"
      tnColdest_early@units 	= "degrees C" 
      
      copy_VarCoords_1(tvar,tnColdest_mid)
      tnColdest_mid!2		= "model"
      tnColdest_mid&model	= model
      tnColdest_mid@short_name 	= "LOCA Climatological mean coldest night of year, 2040-2069"
      tnColdest_mid@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual coldest night of year, 2040-2069"
      tnColdest_mid@units 	= "degrees C" 
      
      copy_VarCoords_1(tvar,tnColdest_late)
      tnColdest_late!2		= "model"
      tnColdest_late&model	= model
      tnColdest_late@short_name = "LOCA Climatological mean coldest night of year, 2070-2099"
      tnColdest_late@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual coldest night of year, 2070-2099"
      tnColdest_late@units 	= "degrees C"
      
      copy_VarCoords_1(tvar,tnHottest_early)
      tnHottest_early!2		= "model"
      tnHottest_early&model	= model
      tnHottest_early@short_name= "LOCA Climatological mean Hottest night of year, 2010-2039"
      tnHottest_early@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual Hottest night of year, 2010-2039"
      tnHottest_early@units 	= "degrees C" 
      
      copy_VarCoords_1(tvar,tnHottest_mid)
      tnHottest_mid!2		= "model"
      tnHottest_mid&model	= model
      tnHottest_mid@short_name 	= "LOCA Climatological mean Hottest night of year, 2040-2069"
      tnHottest_mid@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual Hottest night of year, 2040-2069"
      tnHottest_mid@units 	= "degrees C" 
      
      copy_VarCoords_1(tvar,tnHottest_late)
      tnHottest_late!2		= "model"
      tnHottest_late&model	= model
      tnHottest_late@short_name = "LOCA Climatological mean Hottest night of year, 2070-2099"
      tnHottest_late@long_name	= "Statistically downscaled Localized Constructed Analogs (LOCA) climatological mean annual Hottest night of year, 2070-2099"
      tnHottest_late@units 	= "degrees C"
      
    ;--- Write variables to nc4 file
      setfileoption("nc","FileStructure","Advanced")
      setfileoption("nc","Format","NetCDF4")

    ;--- create file name
      fn 		= "ExtremeTemp_"+scn(iscn)+"_"+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
print("Filename = "+fn)

    ;--- remove old file, if it exists
      system("/bin/rm -f "+odir+fn+".nc")
      
    ;--- create new file   
      ncdf		= addfile(odir+fn+".nc","c")	
      fAtt		= True
      fAtt@title	= fn
      fAtt@orig_fil	= "tasmax_day_LOCA_CMIP5_models_"+scn(iscn)+"_"+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      fAtt@Conventions	= "COARDS/CF-1.0"
      fAtt@creation_date= systemfunc("date")
      fAtt@author	= "Elizabeth Burakowski (elizabeth.burakowski@unh.edu)"
      
      fAtt@description 	= scn(iscn)+": Extreme Temperature climatological means calculated from LOCA statistically downscaled (16th degree) data (Pierce et al. 2014), "+fyr(0)+"-"+fyr(dimsizes(fyr)-1)+""
      
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
      filevardef(ncdf,"DAYS90FCLIM_EARLY",typeof(days90F_early),getvardims(days90F_early))
      filevardef(ncdf,"DAYS90FCLIM_MID",typeof(days90F_mid),getvardims(days90F_mid))
      filevardef(ncdf,"DAYS90FCLIM_LATE",typeof(days90F_late),getvardims(days90F_late))
      filevardef(ncdf,"DAYS95FCLIM_EARLY",typeof(days95F_early),getvardims(days95F_early))
      filevardef(ncdf,"DAYS95FCLIM_MID",typeof(days95F_mid),getvardims(days95F_mid))
      filevardef(ncdf,"DAYS95FCLIM_LATE",typeof(days95F_late),getvardims(days95F_late))
      filevardef(ncdf,"DAYS32FCLIM_EARLY",typeof(days32F_early),getvardims(days32F_early))
      filevardef(ncdf,"DAYS32FCLIM_MID",typeof(days32F_mid),getvardims(days32F_mid))
      filevardef(ncdf,"DAYS32FCLIM_LATE",typeof(days32F_late),getvardims(days32F_late))
      filevardef(ncdf,"DAYS0FCLIM_EARLY",typeof(days0F_early),getvardims(days0F_early))
      filevardef(ncdf,"DAYS0FCLIM_MID",typeof(days0F_mid),getvardims(days0F_mid))
      filevardef(ncdf,"DAYS0FCLIM_LATE",typeof(days0F_late),getvardims(days0F_late))
      filevardef(ncdf,"HOTTESTDAYCLIM_EARLY",typeof(txHottest_early),getvardims(txHottest_early))
      filevardef(ncdf,"HOTTESTDAYCLIM_MID",typeof(txHottest_mid),getvardims(txHottest_mid))
      filevardef(ncdf,"HOTTESTDAYCLIM_LATE",typeof(txHottest_late),getvardims(txHottest_late))
      filevardef(ncdf,"COLDESTDAYCLIM_EARLY",typeof(txColdest_early),getvardims(txColdest_early))
      filevardef(ncdf,"COLDESTDAYCLIM_MID",typeof(txColdest_mid),getvardims(txColdest_mid))
      filevardef(ncdf,"COLDESTDAYCLIM_LATE",typeof(txColdest_late),getvardims(txColdest_late))
      filevardef(ncdf,"HOTTESTNIGHTCLIM_EARLY",typeof(tnHottest_early),getvardims(tnHottest_early))
      filevardef(ncdf,"HOTTESTNIGHTCLIM_MID",typeof(tnHottest_mid),getvardims(tnHottest_mid))
      filevardef(ncdf,"HOTTESTNIGHTCLIM_LATE",typeof(tnHottest_late),getvardims(tnHottest_late))
      filevardef(ncdf,"COLDESTNIGHTCLIM_EARLY",typeof(tnColdest_early),getvardims(tnColdest_early))
      filevardef(ncdf,"COLDESTNIGHTCLIM_MID",typeof(tnColdest_mid),getvardims(tnColdest_mid))
      filevardef(ncdf,"COLDESTNIGHTCLIM_LATE",typeof(tnColdest_late),getvardims(tnColdest_late))
      filevardef(ncdf,"ModelName",typeof(modelname),getvardims(modelname))
      
    ;--- Define file attributes
      filevarattdef(ncdf,"model",model)
      filevarattdef(ncdf,"lat",lat)
      filevarattdef(ncdf,"lon",lon)
      filevarattdef(ncdf,"DAYS90FCLIM_EARLY",days90F_early)
      filevarattdef(ncdf,"DAYS90FCLIM_MID",days90F_mid)
      filevarattdef(ncdf,"DAYS90FCLIM_LATE",days90F_late)
      filevarattdef(ncdf,"DAYS95FCLIM_EARLY",days95F_early)
      filevarattdef(ncdf,"DAYS95FCLIM_MID",days95F_mid)
      filevarattdef(ncdf,"DAYS95FCLIM_LATE",days95F_late)
      filevarattdef(ncdf,"DAYS32FCLIM_EARLY",days32F_early)
      filevarattdef(ncdf,"DAYS32FCLIM_MID",days32F_mid)
      filevarattdef(ncdf,"DAYS32FCLIM_LATE",days32F_late)
      filevarattdef(ncdf,"DAYS0FCLIM_EARLY",days0F_early)
      filevarattdef(ncdf,"DAYS0FCLIM_MID",days0F_mid)
      filevarattdef(ncdf,"DAYS0FCLIM_LATE",days0F_late)
      filevarattdef(ncdf,"HOTTESTDAYCLIM_EARLY",txHottest_early)
      filevarattdef(ncdf,"HOTTESTDAYCLIM_MID",txHottest_mid)
      filevarattdef(ncdf,"HOTTESTDAYCLIM_LATE",txHottest_late)
      filevarattdef(ncdf,"COLDESTDAYCLIM_EARLY",txColdest_early)
      filevarattdef(ncdf,"COLDESTDAYCLIM_MID",txColdest_mid)
      filevarattdef(ncdf,"COLDESTDAYCLIM_LATE",txColdest_late)
      filevarattdef(ncdf,"HOTTESTNIGHTCLIM_EARLY",tnHottest_early)
      filevarattdef(ncdf,"HOTTESTNIGHTCLIM_MID",tnHottest_mid)
      filevarattdef(ncdf,"HOTTESTNIGHTCLIM_LATE",tnHottest_late)
      filevarattdef(ncdf,"COLDESTNIGHTCLIM_EARLY",tnColdest_early)
      filevarattdef(ncdf,"COLDESTNIGHTCLIM_MID",tnColdest_mid)
      filevarattdef(ncdf,"COLDESTNIGHTCLIM_LATE",tnColdest_late)
      filevarattdef(ncdf,"ModelName",modelname)

      setfileoption(ncdf,"DefineMode",False)
      
    ;--- write variable to file
      ncdf->model			= (/model/)
      ncdf->lat				= (/lat/)
      ncdf->lon	 			= (/lon/)
      ncdf->DAYS90FCLIM_EARLY		= (/days90F_early/)
      ncdf->DAYS90FCLIM_MID		= (/days90F_mid/)
      ncdf->DAYS90FCLIM_LATE		= (/days90F_late/)
      ncdf->DAYS95FCLIM_EARLY		= (/days95F_early/)
      ncdf->DAYS95FCLIM_MID		= (/days95F_mid/)
      ncdf->DAYS95FCLIM_LATE		= (/days95F_late/)
      ncdf->DAYS32FCLIM_EARLY		= (/days32F_early/)
      ncdf->DAYS32FCLIM_MID		= (/days32F_mid/)
      ncdf->DAYS32FCLIM_LATE		= (/days32F_late/)
      ncdf->DAYS0FCLIM_EARLY		= (/days0F_early/)
      ncdf->DAYS0FCLIM_MID		= (/days0F_mid/)
      ncdf->DAYS0FCLIM_LATE		= (/days0F_late/)
      ncdf->HOTTESTDAYCLIM_EARLY	= (/txHottest_early/)
      ncdf->HOTTESTDAYCLIM_MID		= (/txHottest_mid/)
      ncdf->HOTTESTDAYCLIM_LATE		= (/txHottest_late/)
      ncdf->COLDESTDAYCLIM_EARLY	= (/txColdest_early/)
      ncdf->COLDESTDAYCLIM_MID		= (/txColdest_mid/)
      ncdf->COLDESTDAYCLIM_LATE		= (/txColdest_late/)
      ncdf->HOTTESTNIGHTCLIM_EARLY	= (/tnHottest_early/)
      ncdf->HOTTESTNIGHTCLIM_MID	= (/tnHottest_mid/)
      ncdf->HOTTESTNIGHTCLIM_LATE	= (/tnHottest_late/)
      ncdf->COLDESTNIGHTCLIM_EARLY	= (/tnColdest_early/)
      ncdf->COLDESTNIGHTCLIM_MID	= (/tnColdest_mid/)
      ncdf->COLDESTNIGHTCLIM_LATE	= (/tnColdest_late/)
      ncdf->ModelName	= (/modelname/)

end do    ; scenario 	    
 
end

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; obs_ann_hist.ncl 
; EA Burakowski
; 2017-08-07
; 
; obs_ann_hist.ncl 
;	(1) selects a regional subset upon which to 
;	(2) calculate annual sum (pr) and mean (tmax, tmin), then
;	(3) calculate the climatological average & linear trend in 
;	annual temperature (min & max) and precipitation (total) for 
;	a hybrid NLDAS2-PRISM generated observational daily data. 
;	Inputs include daily NLDAS2-PRISM data for three variables:
;		(1) tasmin (surface minimum air temperature)
;		(2) tasmax (surface maximum air temperature)
;		(3) precipitation (daily total)
;		over the period 1980-2005.  
;	
;	More info on the hybrid NLDAS2-PRISM
;	data can be found in:
;	 	Abatzoglou et al. 2013, Int. J. Climatology,33: 121-131. 
;		doi: 10.1002/joc.3413. 
;	
; Output includes:
;   .nc files for 30-yr climatological mean and annual trend in each variable:
; 	- obs_ann_pr_historical_1980-2005.nc
;	- obs_ann_tasmax_historical_1980-2005.nc
;	- obs_ann_tasmin_historical_1980-2005.nc
;	- obs_ann_pr_historical_1980-2005_trend.nc
;	- obs_ann_tasmax_historical_1980-2005_trend.nc
;	- obs_ann_tasmin_historical_1980-2005_trend.nc
;
; Use plotANN_clim.ncl and plotANN_trends.ncl to plot up the .nc files
;	

;======================================================================
;                         The main code
;======================================================================

begin

;---- METDATA directory & constants

  dir	= "/net/nfs/merrimack/raid2/data/METDATA/"
  ovar	= (/"pr","tmmx","tmmn"/)
  var	= (/"pr","tasmax","tasmin"/)
  hyr   = ispan(1980,2005,1)
  
;---- LOCA directory & constants
  ldir 	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/clim/"
  lfil	= addfile(ldir+"ANN_pr_historical_1980-2005.nc","r")
  tvar	= lfil->ANNCLIM
  llat 	= tvar&lat
  llon	= tvar&lon
  llon2	= llon-360  ; convert to degrees_west (-)
  
;---- (1) select regional subset bounds (NE US below)
  latmin		= min(llat) 
  latmax		= max(llat)
  lonmin		= min(llon2)
  lonmax		= max(llon2)

;---- Output directories

  odir	= "/net/home/eos/ean2/METDATA_processed/"	

;---- Add one met file to extract dims
  a	= addfile(dir+ovar(0)+"/pr_2009.nc","r")
  pr1	= a->precipitation_amount(:,{latmin:latmax},{lonmin:lonmax})
  pr	= pr1(lat|:,lon|:,day|:)  ; rearrange for copy_VarCoords_1 below
  dims 	= dimsizes(pr)  ; lat, lon, day 
  lat	= pr&lat
  lon	= pr&lon  

;---- Annual Averages ---------------------------------------
;
;      Loop over variable and then years to calculate the 
;	annual average/total for each year, then the 
;	climatological mean (1979-2005).
;
;------------------------------------------------------------

;---- Loop over vars (use obs var abbreviations = ovar)

  do ivar = 0,dimsizes(ovar)-1
    print("----- Working on variable "+ovar(ivar)+" -----")
   
   ;-- initialize empty variables to hold climatological mean (cl) and trend
    hann_cl	= new((/dims(0),dims(1)/),float,-9999)
    hann_tr	= new((/dims(0),dims(1)/),float,-9999)
    
   ;-- if/else statement for precip vs. tmin/tmax
    if ivar.eq.0
    
    ;-- Initialize empty array to hold annual totals
     hann 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
     copy_VarCoords_1(pr,hann)
     hann!2	= "year"
     hann&year	= hyr
	
	;-- Loop over years to calculate annual total pr for each year
	  do iyr = 1,dimsizes(hyr)-1
            fil		= systemfunc("cd "+dir+"/"+ovar(ivar)+" ; ls *.nc")
	    a		= addfile(dir+"/"+ovar(ivar)+"/"+fil(iyr),"r")
            v		= a->precipitation_amount(:,{latmin:latmax},{lonmin:lonmax})
	    delete(v@missing_value)
	    v		= where(v.lt.0,0,v)	; change negative values to '0' 
	
	  ;--- (2) calculate the annual total for each year, rearrange dims
	    hann(:,:,iyr)	= dim_sum_n_Wrap(v(lat|:,lon|:,day|:),2)

	  ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
            delete([/fil,a,v/])

          end do   ; historical years
	
    else   ; pr/temp
	
   ;-----------------------------
   ;
   ;   TEMPERATURE CALCS
   ;
   ;-----------------------------   
    ;-- initialize empty array to hold annual means
     hann 	= new((/dims(0),dims(1),dimsizes(hyr)/),float,-9999)
     copy_VarCoords_1(pr,hann)
     hann!2	= "year"
     hann&year	= hyr
	
	;-- Loop over years to calculate annual total pr for each year
	  do iyr = 1,dimsizes(hyr)-1 	; start at 1 because METDATA starts in 1979
            fil		= systemfunc("cd "+dir+"/"+ovar(ivar)+" ; ls *.nc")
	    a		= addfile(dir+"/"+ovar(ivar)+"/"+fil(iyr),"r")
            v		= a->air_temperature(:,{latmin:latmax},{lonmin:lonmax})
	
	  ;--- (2) calculate the annual total for each year, rearrange dims
	    hann(:,:,iyr)	= dim_avg_n_Wrap(v(lat|:,lon|:,day|:),2)

	  ;--- delete temporary vars (leap yr will trip up ->LH/RH assignment mismatch)
            delete([/fil,a,v/])

          end do   ; historical years
	  
    end if ; pr/temp 

 ;--- (3) calculate annual climatological mean, historical
   hann_cl(:,:)	= dim_avg_n(hann,2)
   
print("Max 30-yr clim: "+max(hann_cl))
print("Min 30-yr clim: "+min(hann_cl))
 
 ;--- (3) calculate annual trend, historical
   hann_tr(:,:)	= regCoef_n(hyr,hann,0,2)

print("Max trend = "+max(hann_tr))
print("Min trend = "+min(hann_tr))

 ;--- create .nc file of annual climatological averages and trends
    
 ;--- assign coordinate vars for climatological averages
   hann_cl!0		= "lat"
   hann_cl&lat		= lat
   hann_cl!1		= "lon"
   hann_cl&lon		= lon
   hann_cl@short_name 	= "NLDAS2-PRISM Climatological mean "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
   hann_cl@long_name	= "NLDAS2-PRISM Climatological Mean "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
   if ivar.eq.0
     hann_cl@units 	= "mm"
   else 
     hann_cl@units	= "degrees C"
   end if
 

printVarSummary(hann_cl)

 ;--- assign coordinate vars for trends
   hann_tr!0		= "lat"
   hann_tr&lat		= lat
   hann_tr!1		= "lon"
   hann_tr&lon		= lon
   hann_tr@short_name 	= "NLDAS2-PRISM Climatological mean "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
   hann_tr@long_name	= "NLDAS2-PRISM Climatological mean "+var(ivar)+", "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
   if ivar.eq.0
     hann_tr@units 	= "mm year-1"
   else
     hann_tr@units	= "degrees C year-1"
   end if

printVarSummary(hann_tr)

 ;--- Write variables to nc4 file
   setfileoption("nc","FileStructure","Advanced")
   setfileoption("nc","Format","NetCDF4")

 ;--- create file name
   fn 		= "Obs_ANN_"+var(ivar)+"_historical_"+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""
print("Filename = "+fn)

 ;--- remove old file, if it exists
   system("/bin/rm -f "+odir+fn+".nc")
      
 ;--- create new file   
   ncdf			= addfile(odir+fn+".nc","c")	
   fAtt			= True
   fAtt@title		= fn
   fAtt@orig_fil	= ovar(ivar)+"_"+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+".nc"
   fAtt@Conventions	= "COARDS/CF-1.0"
   fAtt@creation_date	= systemfunc("date")
   fAtt@author		= "Elizabeth Burakowski (elizabeth.burakowski@unh.edu)"
      
   fAtt@description 	= "Annual "+var(ivar)+" climatological mean and trend calculated from Abatzoglou et al. 2013 NLDAS2-PRISM hybrid data (Abatzoglou 2013), "+hyr(0)+"-"+hyr(dimsizes(hyr)-1)+""

      
 ;--- file attributes, size
   fileattdef(ncdf,fAtt)
   dimNames		= (/"lat", "lon"/)
   dimSizes		= (/ dims(0), dims(1) /)
   dimUnlim          = (/ False, False /)
   filedimdef(ncdf,dimNames,dimSizes,dimUnlim)

 ;--- Define file variables
   filevardef(ncdf,"lat",typeof(lat),getvardims(lat))
   filevardef(ncdf,"lon",typeof(lon),getvardims(lon))
   filevardef(ncdf,"ANNCLIM",typeof(hann_cl),getvardims(hann_cl))
   filevardef(ncdf,"ANNTREND",typeof(hann_tr),getvardims(hann_tr))

print("Defined file variables")
      
 ;--- Define file attributes
   filevarattdef(ncdf,"lat",lat)
   filevarattdef(ncdf,"lon",lon)
   filevarattdef(ncdf,"ANNCLIM",hann_cl)
   filevarattdef(ncdf,"ANNTREND",hann_tr)

   setfileoption(ncdf,"DefineMode",False)
      
 ;--- write variable to file
   ncdf->lat		= (/lat/)
   ncdf->lon	 	= (/lon/)
   ncdf->ANNCLIM	= (/hann_cl/)
   ncdf->ANNTREND	= (/hann_tr/)

print("Wrote variables to file")

 ;--- delete vars
   delete([/hann_cl,hann_tr/])	     

print("Made it to end of loop")
end do  ; variables


end


 

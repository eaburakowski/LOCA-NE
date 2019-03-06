load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; writeLOCA_GreatBayDaily_csv.ncl 
; EA Burakowski
; 2019-03-04
; 
; writeLOCA_GreatBayDaily_csv.ncl writes LOCA statistically downscaled time series 
; for historical (1980-2005) and future (2006-2099) scenarios 
; (rcp45, rcp85) for select stations in the Great Bay Watershed:
;	
;	Concord NH	43.2049  -71.5027
;	Durham NH	43.1432	 -70.9501
;	Epping NH	43.0302  -71.0838
;	Lawrence MA	42.6993	 -71.1659
;	Sanford ME	43.4569  -70.7802
; 
; User can choose from the following climate time series:
;	- Daily Precipitation (pr)
; 	- Daily Maximum Temperature (tasmax)
;	- Daily Minimum Temperature (tasmin)
;
; Output includes:
;   .csv files of complete 1980-2099 daily time series for selected variables.
;
;======================================================================
;                         The main code
;======================================================================

begin

;--- CMIP5 model names
  ovars   = (/"Year,Month,Day,ACCESS1-0,ACCESS1-3,bcc-csm1-1-m,bcc-csm1-1,CanESM2,CCSM4,CESM1-BGC,CESM1-CAM5,CMCC-CM,CMCC-CMS,CNRM-CM5,CSIRO-Mk3-6-0,FGOALS-g2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,HadGEM2-AO,HadGEM2-CC,HadGEM2-ES,inmcm4,IPSL-CM5A-LR,IPSL-CM5A-MR,MIROC5,MIROC-ESM-CHEM,MIROC-ESM,MPI-ESM-LR,MPI-ESM-MR,MRI-CGCM3,NorESM1-M"/)
  
  mdl   = (/"ACCESS1-0","ACCESS1-3",\
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

;---- LOCA .nc directories
  dir   	= "/net/nfs/merrimack/raid/Northeast_US_Downscaling_cmip5/"
  odir		= "/net/home/eos/ean2/LOCA_processed/LOCA_csv/csv/"
  var		= (/"pr","tasmax","tasmin"/)
  scn		= (/"rcp45","rcp85"/)

;---- Years
  hyr 	= ispan(1980,2005,1)
  fyr 	= ispan(2006,2099,1)
  yr	= array_append_record(hyr,fyr,0)

;-- Read in GreatBayLOCAstations.csv
  gb_dir           = "/net/home/eos/ean2/LOCA_processed/LOCA_csv/GreatBay/"
  filcsv        = "GreatBayLOCAstations.csv"
  delim         = ","
  vcsv          = asciiread(gb_dir+filcsv,-1,"string")
  latcsv        = tofloat(str_get_field(vcsv(1:),3,delim))
  loncsv        = tofloat(str_get_field(vcsv(1:),4,delim))
  loncsv        = where(loncsv.lt.0,loncsv+360,loncsv)
  stationnames  = str_get_field(vcsv(1:),1,delim)

;--- Change directory to csv directory
  cd = systemfunc("cd "+odir)
	
;--- Loop over vars to make one file per variable
  do ivar = 0,dimsizes(var)-1
    print("----  Working on "+var(ivar)+"  -----")
      
      ;--- Loop over scenarios
      do iscn = 0,dimsizes(scn)-1	    
        print("----  Working on "+scn(iscn)+"  -----")

	;--- Loop over stations
	do istation = 0,dimsizes(latcsv)-6
	print("--- Now Working on "+stationnames(istation)+"  -----")

	;--- Create new variable to hold all models
	  ts_GB	= new((/44195,dimsizes(mdl)/),float,-9999)

	;--- Loop over models
        do imodel = 0, dimsizes(mdl)-1
          print("---- Working on "+mdl(imodel)+" -----")

	  ;--- Add historical LOCA daily files
	    hfils	= systemfunc("cd "+dir+" ; ls "+var(ivar)+"_day_"+mdl(imodel)+"_historical_*.nc")
	    a		= addfiles(dir+hfils,"r")
	    hist	= a[:]->$var(ivar)$(:,{latcsv(istation)},{loncsv(istation)})

	    ;--- unit conversions
	      if ivar.eq.0 ; precip
	        hist_mmC	= hist
	   	hist_mmC	= hist*86400
	      else	   ; temperature
		hist_mmC	= hist
	        hist_mmC	= hist - 273.15
	      end if  

	  ;--- Add future LOCA daily files
	    ffils	= systemfunc("cd "+dir+" ; ls "+var(ivar)+"_day_"+mdl(imodel)+"_"+scn(iscn)+"*.nc")
	    b		= addfiles(dir+ffils,"r")
  	    future	= b[:]->$var(ivar)$(:,{latcsv(istation)},{loncsv(istation)})
	  
	    ;--- unit conversions
              if ivar.eq.0 ; precip
                future_mmC	= future
                future_mmC      = future*86400
              else         ; temperature
                future_mmC      = future
                future_mmC      = future - 273.15 
              end if

	  ;--- Append future to historical
	    ts		= array_append_record(hist_mmC,future_mmC,0)

	  ;--- One model is one day short. Find it and replace missing day with a zero.
	    if imodel.eq.3 .and. iscn.eq.1
	      ts_GB(0:dimsizes(ts)-1,imodel) 	= ts
	      ts_GB(dimsizes(ts),imodel)	= 0
	    else
	      ts_GB(:,imodel)			= ts
	    end if

	  ;--- Generate year, month, day from time
	    if imodel.eq.3 .and. iscn.eq.1
	      print("do not generate time variable")
	    else
	      t		= ts&time
	      utc	= cd_calendar(t,0)
	      year	= tointeger(utc(:,0))
	      month	= tointeger(utc(:,1))
	      day	= tointeger(utc(:,2))
	    end if

	  ;--- Clear unnecessary variables
	    delete([/hfils,a,hist,hist_mmC,ffils,b,future,future_mmC,ts/])

       end do ; models

          ;----------------------------------------------------------------------
          ;     Write data to table
          ;----------------------------------------------------------------------
	    header		= (/ovars/)
	    hlist		= [/header/]
	    alist		= [/year,month,day,ts_GB(:,0),ts_GB(:,1),ts_GB(:,2),ts_GB(:,3),ts_GB(:,4),ts_GB(:,4),ts_GB(:,5),\
				    ts_GB(:,6),ts_GB(:,7),ts_GB(:,8),ts_GB(:,9),ts_GB(:,10),\
				    ts_GB(:,11),ts_GB(:,12),ts_GB(:,13),ts_GB(:,14),ts_GB(:,15),\
				    ts_GB(:,16),ts_GB(:,17),ts_GB(:,18),ts_GB(:,19),ts_GB(:,20),\
				    ts_GB(:,21),ts_GB(:,22),ts_GB(:,23),ts_GB(:,24),ts_GB(:,25),ts_GB(:,26),ts_GB(:,27),ts_GB(:,28)/]
	    write_table(stationnames(istation)+"_"+var(ivar)+"_"+scn(iscn)+"_1980-2099_daily.csv","w",hlist,"%s")
	    write_table(stationnames(istation)+"_"+var(ivar)+"_"+scn(iscn)+"_1980-2099_daily.csv","a",alist,"%i,%i,%i,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f,%2.3f")
         
      ;--- Delete unnecessary variables
      delete([/ts_GB,year,month,day,alist/])

      end do ; stations 
    end do ; scenario
 end do ; variable

end

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; writeLOCA_csv.ncl 
; EA Burakowski
; 2018-05-02
; 
; writeLOCA_csv.ncl writes LOCA statistically downscaled time series 
; for historical (1980-2005) and future (2006-2099) scenarios 
; (rcp45, rcp85) for individual stations specified in LOCA_latlon.csv. 
; User can choose from the following climate indicator time series:
; 	- Annual Maximum Temperature (ANN_tasmax)
;	- Annual Minimum Temperature (ANN_tasmin)
;	- Annual Precipitation (ANN_pr)
;	- Seasonal Maximum (tasmax), minimum (tasmin) Temperature or Precipitation
;	  Seasons designated by DJF, MAM, JJA, and/or SON
;	    ie: DJF_tasmax is winter maximum temperature for Dec - Feb
;	- Extreme Temperature metrics, including:
;	    - Days per year > 90 F (DAYS90F)
;	    - Days per year > 95 F (DAYS95F)
;	    - Days per year < 32 F (DAYS32F)
;	    - Days per year < 0F (DAYS0F)
;	    - Temperature on hottest day of year (HOTTESTDAY) from tasmax
;	    - Temperature on coldest day of year (COLDESTDAY) from tasmax
;	    - Temperature on hottest night of year (HOTTESTNIGHT) from tasmin
;	    - Temperature on coldest night of year (COLDESTNIGHT) from tasmin
;	- Extreme Precipitation metrics, including:
;	    - TBD
;
; Output includes:
;   .csv files of complete 1980-2099 time series for stations indicated in 
;   input file LOCAstations.csv. The LOCAstations.csv is formatted as follows:
;	- Column 1: City
;	- Column 2: State
;	- Column 3: Lat (N)
;	- Column 4: Lon (- west)
;	- GHCN station id (optional), enter as -9999 if blank
;
;======================================================================
;                         The main code
;======================================================================

begin

;---- LOCA .nc directories
  hdir 		= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/historical/timeseries/"
  fdir		= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/future/timeseries/"
  odir		= "/net/home/eos/ean2/LOCA_processes/LOCA_csv/csv/"
  filprefix	= (/"ANN_pr","ANN_tasmax","ANN_tasmin",\
		    "SEAS_pr","SEAS_tasmax","SEAS_tasmin",\
		    "ExtremeTemp_"/)
  seas_var	= (/"NaN", "NaN", "NaN",\
		    "pr", "tasmax", "tasmin"/)
  scn		= (/"rcp45","rcp85"/)

;---- Years
  hyr 	= ispan(1980,2005,1)
  fyr 	= ispan(2006,2099,1)
  yr	= array_append_record(hyr,fyr,0)
		
;---- Load one historical climatology (1980-2005) file (created using
;	ann_hist.ncl) to get the lat/lon dims:
	
  hfil 	= addfile(hdir+filprefix(0)+"_historical_1980-2005_timeseries.nc","r")
  hist	= hfil->ANNTIMESERIES
  dims	= dimsizes(hist)
  lat	= hist&lat
  lon	= hist&lon
  
  ;-- delete 'hist', no longer needed.
  delete(hist)

;-- Read in LOCAstations.csv
  dir  		= "/net/home/eos/ean2/LOCA_processed/LOCA_csv/"
  filcsv	= "LOCAstations.csv"
  delim		= ","
  vcsv		= asciiread(dir+filcsv,-1,"string")
  latcsv	= tofloat(str_get_field(vcsv(1:),3,delim))
  loncsv	= tofloat(str_get_field(vcsv(1:),4,delim))
  loncsv	= where(loncsv.lt.0,loncsv+360,loncsv)
  stationnames	= str_get_field(vcsv(1:),1,delim)

;--- Change directory to csv directory
  cd = systemfunc("cd "+odir)
	
;--- Loop over vars to make one plot per variable
  do ifil = 0,dimsizes(filprefix)-2
    print("----  Working on "+filprefix(ifil)+"  -----")
    
      ;--- if statement for annual, seasonal, extreme temp
      if ifil.le.2	; annual
      
        do iscn = 0, dimsizes(scn)-1
          print("---- Working on "+scn(iscn)+" -----")

	  ;--- loop over stations
	  do istation = 0,dimsizes(latcsv)-1
	  print("---- Working on "+stationnames(istation)+"  ------")

	  ;--- Add historical LOCA time series file, calculate ensemble mean, max, and min
	    a		= addfile(hdir+filprefix(ifil)+"_historical_1980-2005_timeseries.nc","r")
	    hnames	= getfilevarnames(a)
	    
	    hist	= a->$hnames(3)$({latcsv(istation)},{loncsv(istation)},:,:)
	    hist_ens 	= dim_avg_n_Wrap(hist,1)
	    hist_max	= dim_max_n_Wrap(hist,1)
	    hist_min	= dim_min_n_Wrap(hist,1)
	  
	  ;--- Add future LOCA time series file, calculate ensemble mean, max, and min
	  ;    Climatology files created using ann_future.ncl
	    b		= addfile(fdir+filprefix(ifil)+"_"+scn(iscn)+"_2006-2099_timeseries.nc","r")
	    fnames	= getfilevarnames(b)
	    
  	    future	= b->$fnames(3)$({latcsv(istation)},{loncsv(istation)},:,:)	 
	    future_ens	= dim_avg_n_Wrap(future,1)
	    future_max	= dim_max_n_Wrap(future,1)
	    future_min	= dim_min_n_Wrap(future,1)

	  ;--- Combine historical and future into single time series (1980-2005)
	    ts_ens	= array_append_record(hist_ens,future_ens,0)
	    ts_max	= array_append_record(hist_max,future_max,0)
	    ts_min	= array_append_record(hist_min,future_min,0)

          ;----------------------------------------------------------------------
          ;     Write data to table
          ;----------------------------------------------------------------------
	    ovars		= (/"Year,Ensemble_mean,Ensemble_min,Ensemble_max"/)
	    header		= (/ovars/)
	    hlist		= [/header/]
	    alist		= [/yr,ts_ens,ts_min,ts_max/]
	    write_table(stationnames(istation)+"_"+filprefix(ifil)+"_"+scn(iscn)+"_1980-2099_timeseries.csv","w",hlist,"%s")
	    write_table(stationnames(istation)+"_"+filprefix(ifil)+"_"+scn(iscn)+"_1980-2099_timeseries.csv","a",alist,"%i,%2.2f,%2.2f,%2.2f")
          
	  ;--- Clear temporary variables
	  delete([/a,hnames,hist,hist_ens,hist_max,hist_min/])
          delete([/b,fnames,future,future_ens,future_max,future_min/])
	  delete([/ts_ens,ts_max,ts_min,ovars,header,hlist,alist/])

	  end do   ; stations
          
        end do; scenario
	
      ;--------------------------------------------
      else if ifil.ge.3 .and. ifil.le.5	;  seasonal
      ;--------------------------------------------
      
        do iscn = 0, dimsizes(scn)-1
          print("---- Working on "+scn(iscn)+" -----")

	  ;--- loop over stations
	  do istation = 0,dimsizes(latcsv)-1
	  print("---- Working on "+stationnames(istation)+"  ------")

	  ;--- Add historical LOCA time series file, calculate ensemble mean, max, and min
	    a		= addfile(hdir+filprefix(ifil)+"_historical_1980-2005_timeseries.nc","r")
	    hnames	= getfilevarnames(a)
	    histmam	= a->$hnames(3)$({latcsv(istation)},{loncsv(istation)},:,:)
	    histmam_ens = dim_avg_n_Wrap(histmam,1)
	    histmam_max	= dim_max_n_Wrap(histmam,1)
	    histmam_min	= dim_min_n_Wrap(histmam,1)
	    
	    histjja	= a->$hnames(4)$({latcsv(istation)},{loncsv(istation)},:,:)
	    histjja_ens = dim_avg_n_Wrap(histjja,1)
	    histjja_max	= dim_max_n_Wrap(histjja,1)
	    histjja_min	= dim_min_n_Wrap(histjja,1)
	    
	    histson	= a->$hnames(5)$({latcsv(istation)},{loncsv(istation)},:,:)
	    histson_ens = dim_avg_n_Wrap(histson,1)
	    histson_max	= dim_max_n_Wrap(histson,1)
	    histson_min	= dim_min_n_Wrap(histson,1)
	    
	    histdjf	= a->$hnames(6)$({latcsv(istation)},{loncsv(istation)},:,:)
	    histdjf_ens = dim_avg_n_Wrap(histdjf,1)
	    histdjf_max	= dim_max_n_Wrap(histdjf,1)
	    histdjf_min	= dim_min_n_Wrap(histdjf,1)
	  
	  ;--- Add future LOCA time series file, calculate ensemble mean, max, and 
	  ;    min for each season
	    b			= addfile(fdir+filprefix(ifil)+"_"+scn(iscn)+"_2006-2099_timeseries.nc","r")
	    fnames		= getfilevarnames(b)
	    
	    futuremam		= b->$fnames(3)$({latcsv(istation)},{loncsv(istation)},:,:)
	    futuremam_ens 	= dim_avg_n_Wrap(futuremam,1)
	    futuremam_max	= dim_max_n_Wrap(futuremam,1)
	    futuremam_min	= dim_min_n_Wrap(futuremam,1)
	    
	    futurejja		= b->$fnames(4)$({latcsv(istation)},{loncsv(istation)},:,:)
	    futurejja_ens 	= dim_avg_n_Wrap(futurejja,1)
	    futurejja_max	= dim_max_n_Wrap(futurejja,1)
	    futurejja_min	= dim_min_n_Wrap(futurejja,1)
	    
	    futureson		= b->$fnames(5)$({latcsv(istation)},{loncsv(istation)},:,:)
	    futureson_ens 	= dim_avg_n_Wrap(futureson,1)
	    futureson_max	= dim_max_n_Wrap(futureson,1)
	    futureson_min	= dim_min_n_Wrap(futureson,1)
	    
	    futuredjf		= b->$fnames(6)$({latcsv(istation)},{loncsv(istation)},:,:)
	    futuredjf_ens 	= dim_avg_n_Wrap(futuredjf,1)
	    futuredjf_max	= dim_max_n_Wrap(futuredjf,1)
	    futuredjf_min	= dim_min_n_Wrap(futuredjf,1)

	  ;--- Combine historical and future into single time series (1980-2099) for each season
	    tsmam_ens	= array_append_record(histmam_ens,futuremam_ens,0)
	    tsmam_max	= array_append_record(histmam_max,futuremam_max,0)
	    tsmam_min	= array_append_record(histmam_min,futuremam_min,0)
	    
	    tsjja_ens	= array_append_record(histjja_ens,futurejja_ens,0)
	    tsjja_max	= array_append_record(histjja_max,futurejja_max,0)
	    tsjja_min	= array_append_record(histjja_min,futurejja_min,0)
	    
	    tsson_ens	= array_append_record(histson_ens,futureson_ens,0)
	    tsson_max	= array_append_record(histson_max,futureson_max,0)
	    tsson_min	= array_append_record(histson_min,futureson_min,0)
	    
	    tsdjf_ens	= array_append_record(histdjf_ens,futuredjf_ens,0)
	    tsdjf_max	= array_append_record(histdjf_max,futuredjf_max,0)
	    tsdjf_min	= array_append_record(histdjf_min,futuredjf_min,0)

          ;----------------------------------------------------------------------
          ;     Write data to table, separate table for each season
          ;----------------------------------------------------------------------
	    mamovars		= (/"Year,Ensemble_mean,Ensemble_min,Ensemble_max"/)
	    mamheader		= (/mamovars/)
	    mamhlist		= [/mamheader/]
	    mamalist		= [/yr,tsmam_ens,tsmam_min,tsmam_max/]
	    write_table(stationnames(istation)+"_Spring_"+seas_var(ifil)+"_"+scn(iscn)+"_1980-2099_timeseries.csv","w",mamhlist,"%s")
	    write_table(stationnames(istation)+"_Spring_"+seas_var(ifil)+"_"+scn(iscn)+"_1980-2099_timeseries.csv","a",mamalist,"%i,%2.2f,%2.2f,%2.2f")
	    
	    jjaovars		= (/"Year,Ensemble_mean,Ensemble_min,Ensemble_max"/)
	    jjaheader		= (/jjaovars/)
	    jjahlist		= [/jjaheader/]
	    jjaalist		= [/yr,tsjja_ens,tsjja_min,tsjja_max/]
	    write_table(stationnames(istation)+"_Summer_"+seas_var(ifil)+"_"+scn(iscn)+"_1980-2099_timeseries.csv","w",jjahlist,"%s")
	    write_table(stationnames(istation)+"_Summer_"+seas_var(ifil)+"_"+scn(iscn)+"_1980-2099_timeseries.csv","a",jjaalist,"%i,%2.2f,%2.2f,%2.2f")
	    
	    sonovars		= (/"Year,Ensemble_mean,Ensemble_min,Ensemble_max"/)
	    sonheader		= (/sonovars/)
	    sonhlist		= [/sonheader/]
	    sonalist		= [/yr,tsson_ens,tsson_min,tsson_max/]
	    write_table(stationnames(istation)+"_Fall_"+seas_var(ifil)+"_"+scn(iscn)+"_1980-2099_timeseries.csv","w",sonhlist,"%s")
	    write_table(stationnames(istation)+"_Fall_"+seas_var(ifil)+"_"+scn(iscn)+"_1980-2099_timeseries.csv","a",sonalist,"%i,%2.2f,%2.2f,%2.2f")
	    
	    djfovars		= (/"Year,Ensemble_mean,Ensemble_min,Ensemble_max"/)
	    djfheader		= (/djfovars/)
	    djfhlist		= [/djfheader/]
	    djfalist		= [/yr,tsdjf_ens,tsdjf_min,tsdjf_max/]
	    write_table(stationnames(istation)+"_Winter_"+seas_var(ifil)+"_"+scn(iscn)+"_1980-2099_timeseries.csv","w",djfhlist,"%s")
	    write_table(stationnames(istation)+"_Winter_"+seas_var(ifil)+"_"+scn(iscn)+"_1980-2099_timeseries.csv","a",djfalist,"%i,%2.2f,%2.2f,%2.2f")

	  ;--- Clear temporary variables
	  delete([/a,hnames,histmam,histmam_ens,histmam_max,histmam_min/])
	  delete([/histjja,histjja_ens,histjja_max,histjja_min/])
          delete([/histson,histson_ens,histson_max,histson_min/])
          delete([/histdjf,histdjf_ens,histdjf_max,histdjf_min/])
          delete([/b,fnames,futuremam,futuremam_ens,futuremam_max,futuremam_min/])
	  delete([/futurejja,futurejja_ens,futurejja_max,futurejja_min/])
          delete([/futureson,futureson_ens,futureson_max,futureson_min/])
          delete([/futuredjf,futuredjf_ens,futuredjf_max,futuredjf_min/])
	  delete([/tsmam_ens,tsmam_max,tsmam_min,tsjja_ens,tsjja_max,tsjja_min/])
	  delete([/tsson_ens,tsson_max,tsson_min,tsdjf_ens,tsdjf_max,tsdjf_min/])
          delete([/mamovars,mamheader,mamhlist,mamalist,jjaovars,jjaheader,jjahlist,jjaalist/])
          delete([/sonovars,sonheader,sonhlist,sonalist,djfovars,djfheader,djfhlist,djfalist/])
   
          end do   ; stations
          
        end do; scenario
	
      ;--------------------------------
      else if ifil.eq.6  ; Extreme Temp
      ;--------------------------------
      
        do iscn = 0, dimsizes(scn)-1
          print("---- Working on "+scn(iscn)+" -----")

	  ;--- loop over stations
	  do istation = 0,dimsizes(latcsv)-1
	  print("---- Working on "+stationnames(istation)+"  ------")

	  ;--- Add historical LOCA time series file, calculate ensemble mean, max, and min
	    a			= addfile(hdir+filprefix(ifil)+"_historical_1980-2005_timeseries.nc","r")
	    hnames		= getfilevarnames(a)
	    
	    hdays90		= a->$hnames(3)$({latcsv(istation)},{loncsv(istation)},:,:)
	    hdays90_ens 	= dim_avg_n_Wrap(hdays90,1)
	    hdays90_max		= dim_max_n_Wrap(hdays90,1)
	    hdays90_min		= dim_min_n_Wrap(hdays90,1)
	    
	    hdays95		= a->$hnames(4)$({latcsv(istation)},{loncsv(istation)},:,:)
	    hdays95_ens 	= dim_avg_n_Wrap(hdays95,1)
	    hdays95_max		= dim_max_n_Wrap(hdays95,1)
	    hdays95_min		= dim_min_n_Wrap(hdays95,1)
	    
	    hdays32		= a->$hnames(5)$({latcsv(istation)},{loncsv(istation)},:,:)
	    hdays32_ens 	= dim_avg_n_Wrap(hdays32,1)
	    hdays32_max		= dim_max_n_Wrap(hdays32,1)
	    hdays32_min		= dim_min_n_Wrap(hdays32,1)
	    
	    hdays0		= a->$hnames(6)$({latcsv(istation)},{loncsv(istation)},:,:)
	    hdays0_ens		= dim_avg_n_Wrap(hdays0,1)
	    hdays0_max		= dim_max_n_Wrap(hdays0,1)
	    hdays0_min		= dim_min_n_Wrap(hdays0,1)
	  
	  ;--- Add future LOCA time series file, calculate ensemble mean, max, and 
	  ;    min for each season
	    b			= addfile(fdir+filprefix(ifil)+"_"+scn(iscn)+"_2006-2099_timeseries.nc","r")
	    fnames		= getfilevarnames(b)
	    
	    fdays90		= b->$fnames(3)$({latcsv(istation)},{loncsv(istation)},:,:)
	    fdays90_ens 	= dim_avg_n_Wrap(fdays90,1)
	    fdays90_max		= dim_max_n_Wrap(fdays90,1)
	    fdays90_min		= dim_min_n_Wrap(fdays90,1)
	    
	    fdays95		= b->$fnames(4)$({latcsv(istation)},{loncsv(istation)},:,:)
	    fdays95_ens 	= dim_avg_n_Wrap(futurejja,1)
	    fdays95_max		= dim_max_n_Wrap(futurejja,1)
	    fdays95_min		= dim_min_n_Wrap(futurejja,1)
	    
	    fdays32		= b->$fnames(5)$({latcsv(istation)},{loncsv(istation)},:,:)
	    fdays32_ens 	= dim_avg_n_Wrap(fdays32,1)
	    fdays32_max		= dim_max_n_Wrap(fdays32,1)
	    fdays32_min		= dim_min_n_Wrap(fdays32,1)
	    
	    fdays0		= b->$fnames(6)$({latcsv(istation)},{loncsv(istation)},:,:)
	    fdays0_ens 		= dim_avg_n_Wrap(fdays0,1)
	    fdays0_max		= dim_max_n_Wrap(fdays0,1)
	    fdays0_min		= dim_min_n_Wrap(fdays0,1)

	  ;--- Combine historical and future into single time series (1980-2099) for each season
	    tsdays90_ens	= array_append_record(hdays90_ens,fdays90_ens,0)
	    tsdays90_max	= array_append_record(hdays90_max,fdays90_max,0)
	    tsdays90_min	= array_append_record(hdays90_min,fdays90_min,0)
	    
	    tsdays95_ens	= array_append_record(hdays95_ens,fdays95_ens,0)
	    tsdays95_max	= array_append_record(hdays95_max,fdays95_max,0)
	    tsdays95_min	= array_append_record(hdays95_min,fdays95_min,0)
	    
	    tsdays32_ens	= array_append_record(hdays32_ens,fdays32_ens,0)
	    tsdays32_max	= array_append_record(hdays32_max,fdays32_max,0)
	    tsdays32_min	= array_append_record(hdays32_min,fdays32_min,0)
	    
	    tsdays0_ens		= array_append_record(hdays0_ens,fdays0_ens,0)
	    tsdays0_max		= array_append_record(hdays0_max,fdays0_max,0)
	    tsdays0_min		= array_append_record(hdays0_min,fdays0_min,0)

          ;----------------------------------------------------------------------
          ;     Write data to table, separate tabel for each season
          ;----------------------------------------------------------------------
	    days90ovars		= (/"Year,Ensemble_mean,Ensemble_min,Ensemble_max"/)
	    days90header	= (/days90ovars/)
	    days90hlist		= [/days90header/]
	    days90alist		= [/yr,tsdays90_ens,tsdays90_min,tsdays90_max/]
	    write_table(stationnames(istation)+"_"+scn(iscn)+"_Days90F_1980-2099_timeseries.csv","w",days90hlist,"%s")
	    write_table(stationnames(istation)+"_"+scn(iscn)+"_Days90F_1980-2099_timeseries.csv","a",days90alist,"%i,%2.1f,%2.1f,%2.1f")
	    
	    days95ovars		= (/"Year,Ensemble_mean,Ensemble_min,Ensemble_max"/)
	    days95header	= (/days95ovars/)
	    days95hlist		= [/days95header/]
	    days95alist		= [/yr,tsdays95_ens,tsdays95_min,tsdays95_max/]
	    write_table(stationnames(istation)+"_"+scn(iscn)+"_Days95F_1980-2099_timeseries.csv","w",days95hlist,"%s")
	    write_table(stationnames(istation)+"_"+scn(iscn)+"_Days95F_1980-2099_timeseries.csv","a",days95alist,"%i,%2.1f,%2.1f,%2.1f")
	    
	    days32ovars		= (/"Year,Ensemble_mean,Ensemble_min,Ensemble_max"/)
	    days32header	= (/days32ovars/)
	    days32hlist		= [/days32header/]
	    days32alist		= [/yr,tsdays32_ens,tsdays32_min,tsdays32_max/]
	    write_table(stationnames(istation)+"_"+scn(iscn)+"_Days32F_1980-2099_timeseries.csv","w",days32hlist,"%s")
	    write_table(stationnames(istation)+"_"+scn(iscn)+"_Days32F_1980-2099_timeseries.csv","a",days32alist,"%i,%2.1f,%2.1f,%2.1f")
	    
	    days0ovars		= (/"Year,Ensemble_mean,Ensemble_min,Ensemble_max"/)
	    days0header		= (/days0ovars/)
	    days0hlist		= [/days0header/]
	    days0alist		= [/yr,tsdays0_ens,tsdays0_min,tsdays0_max/]
	    write_table(stationnames(istation)+"_"+scn(iscn)+"_Days0F_1980-2099_timeseries.csv","w",days0hlist,"%s")
	    write_table(stationnames(istation)+"_"+scn(iscn)+"_Days0F_1980-2099_timeseries.csv","a",days0alist,"%i,%2.1f,%2.1f,%2.1f")

	    ;--- Clear temporary variables
            delete([/a,hnames,histdays90,histdays90_ens,histdays90_max,histdays90_min/])
            delete([/histdays95,histdays95_ens,histdays95_max,histdays95_min/])
            delete([/histdays32,histdays32_ens,histdays32_max,histdays32_min/])
            delete([/histdays0,histdays0_ens,histdays0_max,histdays0_min/])
            delete([/b,fnames,futuredays90,futuredays90_ens,futuredays90_max,futuredays90_min/])
            delete([/futuredays95,futuredays95_ens,futuredays95_max,futuredays95_min/])
            delete([/futuredays32,futuredays32_ens,futuredays32_max,futuredays32_min/])
            delete([/futuredays0,futuredays0_ens,futuredays0_max,futuredays0_min/])
            delete([/tsdays90_ens,tsdays90_max,tsdays90_min,tsdays95_ens,tsdays95_max,tsdays95_min/])
            delete([/tsdays32_ens,tsdays32_max,tsdays32_min,tsdays0_ens,tsdays0_max,tsdays0_min/])
            delete([/days90ovars,days90header,days90hlist,days90alist,days95ovars,days95header,days95hlist,days95alist/])
            delete([/days32ovars,days32header,days32hlist,days32alist,days0ovars,days0header,days0hlist,days0alist/])	
     
          end do   ; stations
          
        end do; scenario

     end if     ; if ifil.le.2
     end if	; if ifil.ge.3 .and. ifil.le.5
     end if	; if ifil.ge.6 
	
   end do  ; filprefix

end


 

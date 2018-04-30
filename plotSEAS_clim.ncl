load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; plotSEAS_clims.ncl 
; EA Burakowski
; 2018-04-29
; 
; plotSEAS_clim.ncl plots climatological 30-yr mean seasonal minimum temperature 
;	(tmin), maximum temperature (tmax), and total precipitation derived from 29 
;	individual CMIP5 LOCAlized Constructed Analogs (LOCA) daily downscaled 
;	data, 1980-2005.  Also includes ensemble mean climatology. Seasonal .nc files
;	created using seas_hist.ncl. 
;
; Output includes:
;   .png files for annual clim in each variable:
; 	- $SEAS$_pr_historical_1980-2005_clim.png
;	- $SEAS$_tasmax_historical_1980-2005_clim.png
;	- $SEAS$_tasmin_historical_1980-2005_clim.png
;
; Seasons ($SEAS$) include:
;	- MAM (March, April, May)
;	- JJA (June, July, August)
;	- SON (September, October, November)
;	- DJF (December, January, February
;	

;======================================================================
;                         The main code
;======================================================================

begin

;---- LOCA .nc directory

  dir	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/clim/"
  var	= (/"pr","tasmax","tasmin"/)
  varname = (/"Total Precip., mm","Tmax, C","Tmin, C"/)
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
	    "NorESM1-M","NLDAS2-PRISM"/)
  seas	= (/"MAM","JJA","SON","DJF"/)   
	   
;---- Observed NLDAS2-PRISM directory
  obdir = "/net/home/eos/ean2/METDATA_processed/"
 
;---- Output directories
  odir	= "/net/home/eos/ean2/LOCA_processed/LOCA_Figures/clim/"	

;---- Load one historical clim file to get the lat/lon dims
  hfil 	= addfile(dir+"MAM_"+var(0)+"_historical_1980-2005.nc","r")
  clim	= hfil->MAMCLIM
printVarSummary(clim)
  dims	= dimsizes(clim)
  lat	= clim&lat
  lon	= clim&lon

;-- Regional domain boundaries  
  minlat = min(lat)
  maxlat = max(lat)
  minlon = min(lon)
  maxlon = max(lon)
  bndadd = 0.20


;--- Loop over vars to make one plot per variable
  do ivar = 0,dimsizes(var)-1
    print("----  Working on "+var(ivar)+"  -----")
      do iseas = 0,dimsizes(seas)
	
	;---- Create new graphic variable to hold plot contents
  	  plot 		= new(30,graphic)
	  
	;--- Add historical LOCA annual clim file 
	  a		= addfile(dir+seas(iseas)+"_"+var(ivar)+"_historical_1980-2005.nc","r")
	  vNames	= getfilevarnames(a)
	  clim		= a->vNames(3)
printVarSummary(clim)
	  
	  
	;--- Add observed NLDAS2-PRISM annual clim file
          ofil 		= addfile(obdir+"Obs_ANN_"+var(ivar)+"_historical_1980-2005.nc","r")
          obsK		= ofil->ANNCLIM
	  obs		= obsK
	  obs		= obsK - 273.15    ; convert to C
printVarSummary(obs)
 
	;--- Regrid obs to LOCA grid
	
	 ;---Specify interpolation to be used
	  if ivar.eq.0
            method      = "conserve"	; use conservative for precip (D. Shea, pers. comm)
	  else
	    method	= "bilinear"	; use bilinear for temperature (D. Shea, pers. comm.)
	  end if
          
	  wgtFileDir  = "/net/home/eos/ean2/NCLscripts/METDATAscripts/esmf_regrid/"
          wgtFileName = "NLDAS2-PRISM_4kmM2_to_LOCA."+method+"_wgt.nc"
          wgtFilePath =  wgtFileDir+wgtFileName 
          print("wgtFilePath="+wgtFilePath)

         ;-- regrid with weights
          obs_rg     	= ESMF_regrid_with_weights(obs,wgtFilePath,False); 

printVarSummary(obs_rg)  ; check dims of regridded output

	;--- Conform dims of regridded obs to clim  
	  obs3d		= conform_dims(dimsizes(clim),obs_rg,(/0,1/))
	  
printVarSummary(obs3d)
	  
	;--- Calculate modeled - obs
	  modobs_diff	= clim
	  modobs_diff	= clim - obs3d

printVarSummary(modobs_diff)
print("max mod-obs = "+max(modobs_diff))
print("min mod-obs = "+min(modobs_diff))


;--- Regrid obs to WRF grid
	
	 ;---Specify interpolation to be used
	  if ivar.eq.0
            method      = "conserve"	; use conservative for precip (D. Shea, pers. comm)
	  else
	    method	= "bilinear"	; use bilinear for temperature (D. Shea, pers. comm.)
	  end if
          
	  wgtFileDir  = "/net/home/eos/ean2/NCLscripts/METDATAscripts/esmf_regrid/"
          wgtFileName = "NLDAS2-PRISM_4kmM2_to_WRF3km."+method+"_wgt.nc"
          wgtFilePath =  wgtFileDir+wgtFileName 
          print("wgtFilePath="+wgtFilePath)

         ;-- regrid with weights
          obs_rg     	= ESMF_regrid_with_weights(obs,wgtFilePath,False); 

printVarSummary(obs_rg)  ; check dims of regridded output

	;--- Conform dims of regridded obs to clim  
	  obs3d		= conform_dims(dimsizes(clim),obs_rg,(/0,1/))
	  
printVarSummary(obs3d)
	  
	;--- Calculate modeled - obs
	  wrfobs_diff	= clim
	  wrfobs_diff	= clim - obs3d

printVarSummary(modobs_diff)
print("max mod-obs = "+max(modobs_diff))
print("min mod-obs = "+min(modobs_diff))

    ;----------------------------------------------------------------------
    ; Climatology Plotting options section
    ;----------------------------------------------------------------------

        pltType = "png"                                 ; plot type
        pltDir = "/net/home/eos/ean2/LOCA_processed/LOCA_Figures/clim/"                           
        pltName = "ANN_"+var(ivar)+"Clim_1980-2005"    ; plot name 
        pltPath = pltDir+pltName                        ; plot path

     wks = gsn_open_wks(pltType,pltPath)		; create workstation for ps file
	
	if ivar.eq.0
          gsn_define_colormap(wks,"CBR_wet")      ; precip color table 
	else
	  gsn_define_colormap(wks,"WhiteBlueGreenYellowRed") ; temperature color table
	end if
        res = True
	
	res@gsnMaximize		 = True
        res@gsnDraw              = False
        res@gsnFrame             = False

        res@cnFillOn             = True              ; color plot desired
        res@cnLinesOn            = False             ; turn off contour lines
        res@cnLineLabelsOn       = False             ; turn off contour labels
        res@cnInfoLabelOn        = False             ; turn off info label (top labels of indvid. plots)
        res@cnFillMode           = "RasterFill"      ; turn raster on
        res@cnLevelSelectionMode = "ManualLevels"    ; Set contour levels manually
	if ivar.eq.0
	  res@cnMinLevelValF	 = 0
	  maxclim		 = max(clim)
	  maxobs		 = max(obs)
	  res@cnMaxLevelValF	 = round(max((/maxclim,maxobs/)),3)
	  res@cnLevelSpacingF	 = 100		     ; contour for prcp
	else
	  res@cnMinLevelValF 	 = round(max((/maxclim,maxobs/)),3)*-1 
          res@cnMaxLevelValF     = round(max((/maxclim,maxobs/)),3)  
	  res@cnLevelSpacingF	 = 2		     	
	end if
 	res@lbLabelBarOn	 = False
	res@mpOutlineOn		 = True
	res@mpDataBaseVersion	 = "MediumRes"
	res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
	res@mpProjection	 = "CylindricalEquidistant"
       
        res@mpLimitMode           = "LatLon"           ; required
        res@mpMinLatF             = minlat-bndadd
        res@mpMaxLatF             = maxlat+bndadd
        res@mpMinLonF             = minlon-bndadd
        res@mpMaxLonF             = maxlon+bndadd
        res@mpCenterLonF          = (minlon + maxlon)*0.5
        res@mpCenterLatF          = (minlat + maxlat)*0.5
	res@pmTickMarkDisplayMode = "Always"
	res@tmXBOn		  = False
	res@tmYROn		  = False
	
	res@gsnLeftString         = ""                ; Turn off left subtitle
        res@gsnRightString        = ""                ; Turn off right subtitle
        res@gsnMajorLatSpacing    = 1
        res@gsnMajorLonSpacing    = 1
        res@gsnMinorLonSpacing    = 1
	
        res@gsnAddCyclic          = False             ; regional grid (changes central meridian)/xwo
        
	;--- Loop over model clims and put plot into graphic variable
	do iplot = 0,dimsizes(mdl)-2
	  plot(iplot) = gsn_csm_contour_map(wks,clim(:,:,iplot),res)
	end do

        ;--- Add plot of hybrid NLDAS-PRISM obs data
	plot(29) = gsn_csm_contour_map(wks,obs_rg,res)

	delete(res@cnLevelSpacingF)
		
	;------------------------------------------------------------
	; Panel Plotting options section
	;------------------------------------------------------------

        pres                    	= True
        pres@gsnMaximize        	= True
        pres@gsnPanelLabelBar   	= True
	
	pres@lbOrientation		= "Horizontal"
	pres@lbLabelFontHeightF		= 0.02
	if ivar.eq.0
	  pres@lbTitleString		= "mm/year, 1980-2005"  ; precip title label
	else
	pres@lbTitleString		= "~F34~0~F~ C, 1980-2005"  ; temperature title label
	end if
	pres@lbTitlePosition		= "Bottom"
	pres@lbTitleFontHeightF		=  0.02
	pres@lbTitleOffsetF		= 0.25
	
	pres@pmLabelBarWidthF		= 0.4
	pres@pmLabelBarOrthogonalPosF 	= -0.05
	pres@pmLabelBarParallelPosF	= -0.005
	
	pres@gsnPanelBottom		= 0.05
	
	pres@gsnPanelYWhiteSpace	= 5
	pres@gsnPanelXWhiteSpace	= 5
	pres@txString           	= varname(ivar)
	pres@txPosXF			= 0.1

	pres@gsnPanelFigureStrings 		= mdl 
	pres@gsnPanelFigureStringsFontHeightF 	= 0.007
	pres@gsnPanelFigureStringsPerimOn 	= False
	
	;--- Now panel plot
	gsn_panel(wks,plot,(/5,6/),pres)

        draw(plot)
        frame(wks)  

  delete(res)
  
    ;----------------------------------------------------------------------
    ; Model-Obs Climatology Plotting options section
    ;----------------------------------------------------------------------

        pltType = "png"                                 ; plot type
        pltDir = "/net/home/eos/ean2/LOCA_processed/LOCA_Figures/clim/"                           
        pltName = "ANN_"+var(ivar)+"Clim_1980-2005_diff"    ; plot name 
        pltPath = pltDir+pltName                        ; plot path

     wks = gsn_open_wks(pltType,pltPath)		; create workstation for ps file
	
	if ivar.eq.0
          gsn_define_colormap(wks,"precip_diff_12lev")      ; precip color table 
	else
	  gsn_define_colormap(wks,"temp_diff_18lev") ; temperature color table
	end if
        res = True
	
	res@gsnMaximize		 = True
        res@gsnDraw              = False
        res@gsnFrame             = False

        res@cnFillOn             = True              ; color plot desired
        res@cnLinesOn            = False             ; turn off contour lines
        res@cnLineLabelsOn       = False             ; turn off contour labels
        res@cnInfoLabelOn        = False             ; turn off info label (top labels of indvid. plots)
        res@cnFillMode           = "RasterFill"      ; turn raster on
        res@cnLevelSelectionMode = "ManualLevels"    ; Set contour levels manually
	if ivar.eq.0
	  res@cnMaxLevelValF	 = 350		; force pr diff symmetrical
          res@cnMinLevelValF	 = -350
          res@cnLevelSpacingF	 = 50
	else
	  maxclim		 	 = max(modobs_diff)
 	  minclim			 = min(modobs_diff)
	  res@cnMaxLevelValF	 = round(maxclim,3)
	  res@cnMinLevelValF	 = round(maxclim*-1,3)
	  res@cnLevelSpacingF	 = 0.2		     	
	end if
 	res@lbLabelBarOn	 = False
	res@mpOutlineOn		 = True
	res@mpDataBaseVersion	 = "MediumRes"
	res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
	res@mpProjection	 = "CylindricalEquidistant"
       
        res@mpLimitMode           = "LatLon"           ; required
        res@mpMinLatF             = minlat-bndadd
        res@mpMaxLatF             = maxlat+bndadd
        res@mpMinLonF             = minlon-bndadd
        res@mpMaxLonF             = maxlon+bndadd
        res@mpCenterLonF          = (minlon + maxlon)*0.5
        res@mpCenterLatF          = (minlat + maxlat)*0.5
	res@pmTickMarkDisplayMode = "Always"
	res@tmXBOn		  = False
	res@tmYROn		  = False
	
	res@gsnLeftString         = ""                ; Turn off left subtitle
        res@gsnRightString        = ""                ; Turn off right subtitle
        res@gsnMajorLatSpacing    = 1
        res@gsnMajorLonSpacing    = 1
        res@gsnMinorLonSpacing    = 1
	
        res@gsnAddCyclic          = False             ; regional grid (changes central meridian)/xwo
        
	;--- Loop over model clims and put plot into graphic variable
	do iplot = 0,dimsizes(mdl)-2
	  plot(iplot) = gsn_csm_contour_map(wks,modobs_diff(:,:,iplot),res)
	end do

	delete(res@cnLevelSpacingF)
		
	;------------------------------------------------------------
	; Panel Plotting options section
	;------------------------------------------------------------

        pres                    	= True
        pres@gsnMaximize        	= True
        pres@gsnPanelLabelBar   	= True
	
	pres@lbOrientation		= "Horizontal"
	pres@lbLabelFontHeightF		= 0.02
	if ivar.eq.0
	  pres@lbTitleString		= "MOD-OBS in mm/year, 1980-2005"  ; precip title label
	else
	pres@lbTitleString		= "MOD-OBS in ~F34~0~F~ C, 1980-2005"  ; temperature title label
	end if
	pres@lbTitlePosition		= "Bottom"
	pres@lbTitleFontHeightF		= 0.02
	pres@lbTitleOffsetF		= 0.25
	
	pres@pmLabelBarWidthF		= 0.4
	pres@pmLabelBarOrthogonalPosF 	= -0.05
	pres@pmLabelBarParallelPosF	= -0.005
	
	pres@gsnPanelBottom		= 0.05
	
	pres@gsnPanelYWhiteSpace	= 5
	pres@gsnPanelXWhiteSpace	= 5
	pres@txString           	= varname(ivar)
	pres@txPosXF			= 0.1

	pres@gsnPanelFigureStrings 		= mdl 
	pres@gsnPanelFigureStringsFontHeightF 	= 0.007
	pres@gsnPanelFigureStringsPerimOn 	= False
	
	;--- Now panel plot
	gsn_panel(wks,plot,(/5,6/),pres)

        draw(plot)
        frame(wks) 


  end do  ; variables

end


 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; plotANN_future.ncl 
; EA Burakowski
; 2018-04-20
; 
; plotANN_future.ncl plots the climatological 30-yr mean annual minimum temperature 
;	(tmin), maximum temperature (tmax), and total precipitation derived from 29 
;	individual CMIP5 LOCAlized Constructed Analogs (LOCA) daily downscaled 
;	data for early (2010-2039), mid (2040-2069), and late (2070-2099) century.  
;	Panel plots are ensemble mean climatology. 
;
; Output includes:
;   .png files for annual clim in each variable:
; 	- ANN_pr_rcp45_1980_2099.png
;	- ANN_pr_rcp85_1980_2099.png
;	- ANN_tasmax_rcp45_1980_2099.png
;	- ANN_tasmax_rcp85_1980_2099.png
;	- ANN_tasmin_rcp45_1980_2099.png
;	- ANN_tasmin_rcp85_1980_2099.png

;======================================================================
;                         The main code
;======================================================================

begin

;---- LOCA .nc directory
  hdir 	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/historical/"
  fdir	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/future/"
  var	= (/"pr","tasmax","tasmin"/)
  scn	= (/"rcp45","rcp85"/)
  		
;---- Load one historical climatology (1980-2005) file (created using
;	ann_hist.ncl) to get the lat/lon dims:
	
  hfil 	= addfile(hdir+"ANN_"+var(0)+"_historical_1980-2005.nc","r")
  hist	= hfil->ANNCLIM
  dims	= dimsizes(hist)
  lat	= hist&lat
  lon	= hist&lon

;-- Regional domain boundaries  
  minlat = min(lat)
  maxlat = max(lat)
  minlon = min(lon)
  maxlon = max(lon)
  bndadd = 0.20


;--- Loop over vars to make one plot per variable
  do ivar = 0,dimsizes(var)-1
    print("----  Working on "+var(ivar)+"  -----")
      do iscn = 0, dimsizes(scn)-1
        print("---- Working on "+scn(iscn)+" -----")
	
	;---- Create new graphic variable to hold plot contents
  	  plot 	= new(3,graphic)
	  
	;--- Add historical LOCA annual clim file, calculate ensemble mean
	  a		= addfile(hdir+"ANN_"+var(ivar)+"_historical_1980-2005.nc","r")
	  hist		= a->ANNCLIM
	  hist_ens 	= dim_avg_n_Wrap(hist,2)
	  
	;--- Add future LOCA annual climatology file, calculate ensemble mean
	;    Climatology files created using ann_future.ncl
	
	  b		= addfile(fdir+"ANN_"+var(ivar)+"_"+scn(iscn)+"_2006-2099_trend_clim.nc","r")
	  early		= b->ANNCLIM_EARLY
	  early_diff	= hist
	  early_diff	= early - hist
	  early_ens	= dim_avg_n_Wrap(early_diff,2)
          mid		= b->ANNCLIM_MID
 	  mid_diff	= hist
	  mid_diff	= mid - hist
	  mid_ens	= dim_avg_n_Wrap(mid_diff,2)
	  late		= b->ANNCLIM_LATE
	  late_diff	= hist
	  late_diff	= late - hist
	  late_ens	= dim_avg_n_Wrap(late_diff,2)
print("Mean "+var(ivar)+" "+scn(iscn)+"Early Ens = "+avg(early_ens))
print("Mean "+var(ivar)+" "+scn(iscn)+"Mid Ens = "+avg(mid_ens))
print("Mean "+var(ivar)+" "+scn(iscn)+"Late Ens = "+avg(late_ens))
	 
    ;----------------------------------------------------------------------
    ; Climatology Plotting options section
    ;----------------------------------------------------------------------

        pltType = "png"                                 ; plot type
        pltDir = "/net/home/eos/ean2/LOCA_processed/LOCA_Figures/LOCA_panel/"                     
        pltName = "ANN_"+var(ivar)+"_"+scn(iscn)+"_1980-2099_clim"    ; plot name 
        pltPath = pltDir+pltName                        ; plot path

     wks = gsn_open_wks(pltType,pltPath)		; create workstation for ps file
	
	if ivar.eq.0
          gsn_define_colormap(wks,"CBR_wet")      ; precip color table 
	else
	  gsn_define_colormap(wks,"MPL_YlOrRd") ; temperature color table
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
	  res@cnMaxLevelValF	 = 250
	  res@cnLevelSpacingF	 = 25		     ; contour for prcp
	else
	  res@cnMinLevelValF 	 = 0
          res@cnMaxLevelValF     = 7.0 
	  res@cnLevelSpacingF	 = 0.25		     	
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
        
	;--- Early (2010-2039)
	plot(0)	= gsn_csm_contour_map(wks,early_ens,res)

	;--- Mid (2040-2069)
	plot(1) = gsn_csm_contour_map(wks,mid_ens,res)

	;--- Late (2070-2099)
	plot(2) = gsn_csm_contour_map(wks,late_ens,res)

	delete(res@cnLevelSpacingF)	; not sure if this needs to be done...
		
	;------------------------------------------------------------
	; Panel Plotting options section
	;------------------------------------------------------------

        pres                    	= True
        pres@gsnMaximize        	= True
        pres@gsnPanelLabelBar   	= True
	
	pres@lbOrientation		= "Horizontal"
	pres@lbLabelFontHeightF		= 0.02
	if ivar.eq.0
	  pres@lbTitleString		= "Change in Annual Precipitation (mm)"  ; precip title label
	else
	pres@lbTitleString		= "Change in Annual Temperature (~F34~0~F~ C)"  ; temperature title label
	end if
	pres@lbTitlePosition		= "Bottom"
	pres@lbTitleFontHeightF		=  0.02
	pres@lbTitleOffsetF		= 0.25
	pres@lbBoxEndCapStyle		="TriangleBothEnds"
	
	pres@pmLabelBarWidthF		= 0.75
	pres@pmLabelBarOrthogonalPosF 	= -0.05 ; 
	pres@pmLabelBarParallelPosF	= 0 ; changed this one from -0.005	
 	pres@gsnPanelBottom		= 0.05
	
	pres@gsnPanelYWhiteSpace	= 5
	pres@gsnPanelXWhiteSpace	= 5
	pres@txPosXF			= 0.1

	period				= (/"Early, 2010-2039","Mid, 2040-2069","Late, 2070-2099"/)
	pres@gsnPanelFigureStrings 		= period
	pres@gsnPanelFigureStringsFontHeightF 	= 0.013
	pres@gsnPanelFigureStringsPerimOn 	= False
	
	;--- Now panel plot
	gsn_panel(wks,plot,(/1,3/),pres)

        draw(plot)

  	delete(res)
  
      end do; scenario
   end do  ; variables

end


 

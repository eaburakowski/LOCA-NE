load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; plotSEAS_future.ncl 
; EA Burakowski
; 2018-04-20
; 
; plotSEAS_future.ncl plots the climatological 30-yr mean seasonal minimum temperature 
;	(tmin), maximum temperature (tmax), and total precipitation derived from 29 
;	individual CMIP5 LOCAlized Constructed Analogs (LOCA) daily downscaled 
;	data for early (2010-2039), mid (2040-2069), and late (2070-2099) century.  
;	Panel plots are ensemble mean climatology. 
;
; Output includes:
;   Three-panel plot .png files for seasonal ($SEAS$ : MAM, JJA, SON, DJF) 
;	clim in each variable:
; 	- $SEAS$_pr_rcp45_1980_2099.png
;	- $SEAS$_pr_rcp85_1980_2099.png

;======================================================================
;                         The main code
;======================================================================

begin

;---- LOCA .nc directory
  hdir 	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/historical/"
  fdir	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/future/"
  var	= (/"pr","tasmax","tasmin"/)
  scn	= (/"rcp45","rcp85"/)
  seas	= (/"MAM","JJA","SON","DJF"/)
  		
;---- Load one historical climatology (1980-2005) file (created using
;	ann_hist.ncl) to get the lat/lon dims:
	
  hfil 	= addfile(hdir+"MAM_"+var(0)+"_historical_1980-2005.nc","r")
  hist	= hfil->MAMCLIM
printVarSummary(hist)
  dims	= dimsizes(hist)
  lat	= hist&lat
  lon	= hist&lon

;-- Regional domain boundaries  
  minlat = min(lat)
  maxlat = max(lat)
  minlon = min(lon)
  maxlon = max(lon)
  bndadd = 0.20

;-- Contour levels and interval for precip (season specific)
;   Order is spring, summer, fall, winter.
pr_MinLev	= (/-5, -5, -5, -5/)
pr_MaxLev	= (/10, 10, 10, 10/)
pr_LevSpace	= (/ 1,  1,  1,  1/)

;--- Loop over vars to make one plot per variable
  do ivar = 0,dimsizes(var)-1
    print("----  Working on "+var(ivar)+"  -----")
      do iscn = 0, dimsizes(scn)-1
        print("---- Working on "+scn(iscn)+" -----")
	do iseas = 0,dimsizes(seas)-1

	;---- Create new graphic variable to hold plot contents
  	  plot 	= new(3,graphic)
	  
	;--- Add historical LOCA seasonal clim file, calculate ensemble mean
	  a		= addfile(hdir+seas(iseas)+"_"+var(ivar)+"_historical_1980-2005.nc","r")
	  vNames_hist	= getfilevarnames(a)
	  hist		= a->$vNames_hist(3)$
printVarSummary(hist)
print("Max "+var(ivar)+" Historical = "+max(hist))
	  hist_ens 	= dim_avg_n_Wrap(hist,2)
	  
	;--- Add future LOCA seasonal climatology file, calculate ensemble mean
	;    Climatology files created using seas_future.ncl
	
	  b		= addfile(fdir+seas(iseas)+"_"+var(ivar)+"_"+scn(iscn)+"_2006-2099_clim.nc","r")
	  vNames_future	= getfilevarnames(b)
	  early		= b->$vNames_future(3)$	; assumes all files have var names in same order (they should!)
	  early_diff	= hist
	  early_diff	= early - hist
	  early_ens	= dim_avg_n_Wrap(early_diff,2)
printVarSummary(early_ens)
print("Max "+var(ivar)+" "+scn(iscn)+" Early Ens Avg = "+max(early_ens))
          mid		= b->$vNames_future(4)$
 	  mid_diff	= hist
	  mid_diff	= mid - hist
	  mid_ens	= dim_avg_n_Wrap(mid_diff,2)
print("Max "+var(ivar)+" "+scn(iscn)+" Mid Ens Avg = "+max(mid_ens))
	  late		= b->$vNames_future(5)$
	  late_diff	= hist
	  late_diff	= late - hist
	  late_ens	= dim_avg_n_Wrap(late_diff,2)
print("Max "+var(ivar)+" "+scn(iscn)+" Late Ensemble avg = "+max(late_ens))	 
    ;----------------------------------------------------------------------
    ; Climatology Plotting options section
    ;----------------------------------------------------------------------

        pltType = "png"                                 ; plot type
        pltDir = "/net/home/eos/ean2/LOCA_processed/LOCA_Figures/LOCA_panel/"                     
        pltName = seas(iseas)+"_"+var(ivar)+"_"+scn(iscn)+"_1980-2099_clim"    ; plot name 
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
        res@cnLevelSelectionMode = "ManualLevels"    ; Set contour levels manuallyi

	;-- Sets of min and max contour intervals by season
	if ivar.eq.0
	  res@cnMinLevelValF	 = pr_MinLev(iseas)
	  res@cnMaxLevelValF	 = pr_MaxLev(iseas)
	  res@cnLevelSpacingF	 = pr_LevSpace(iseas)		     ; contour for prcp
	else
	  res@cnMinLevelValF 	 = 0
          res@cnMaxLevelValF     = 9
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

	;-- Season and scenario labels
	seas_label 			= (/"Spring","Summer","Fall","Winter"/)
	scen_label			= (/"RCP4.5","RCP8.5"/)
	
        ;-- Variable (pr, tasmax, tasmin) specific lables
	if ivar.eq.0
	  pres@lbTitleString		= scen_label(iscn)+" : Change in "+seas(iseas)+" Precipitation (mm), relative to 1980-2005"  ; precip title label
	else
	  pres@lbTitleString		= scen_label(iscn)+" : Change in "+seas(iseas)+" Temperature (~F34~0~F~ C)"  ; temperature title label
	end if
	
	pres@lbTitlePosition		= "Bottom"
	pres@lbTitleFontHeightF		=  0.02
	pres@lbTitleOffsetF		= 0.25
	pres@lbBoxEndCapStyle		="TriangleBothEnds"
	
	pres@pmLabelBarWidthF		= 0.5
	pres@pmLabelBarOrthogonalPosF 	= -0.05 ; 
;	pres@pmLabelBarParallelPosF	= 0 ; changed this one from -0.005	
 	pres@gsnPanelBottom		= 0.05
	
	pres@gsnPanelYWhiteSpace	= 5
	pres@gsnPanelXWhiteSpace	= 5
	pres@txPosXF			= 0.1

	period				= (/"Early, 2010-2039","Mid, 2040-2069","Late, 2070-2099"/)
	pres@gsnPanelFigureStrings 		= period
	pres@gsnPanelFigureStringsFontHeightF 	= 0.01
	pres@gsnPanelFigureStringsPerimOn 	= False
	
	;--- Now panel plot
	gsn_panel(wks,plot,(/1,3/),pres)

        draw(plot)

  	delete(res)
	
	end do  ; season  
      end do; scenario
   end do  ; variables

end

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
  scn	= (/"rcp45","rcp85"/)
  metrics = (/"DAYS90F","DAYS95F","DAYS32F","DAYS0F","HOTTESTDAY","COLDESTDAY","HOTTESTNIGHT","COLDESTNIGHT"/)
  		
;---- Load one historical climatology (1980-2005) file (created using
;	ExtremeTemp_hist.ncl) to get the lat/lon dims:
	
  hfil 	= addfile(hdir+"ExtremeTemp_historical_1980-2005.nc","r")
  hist	= hfil->DAYS90FCLIM
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
  do iscn = 0,dimsizes(scn)-1
    print("----  Working on "+scn(iscn)+"  -----")
      a 	= addfile(hdir+"ExtremeTemp_historical_1980-2005.nc","r")
      hvNames 	= getfilevarnames(a)
      b 	= addfile(fdir+"ExtremeTemp_"+scn(iscn)+"_2006-2099.nc","r")
      fvNames 	= getfilevarnames(b)
 
      ;--- Loop over extreme temp metrics
      do imetric = 0, dimsizes(metrics)-1
        print("---- Working on "+metrics(imetric)+" -----")
	
	;---- Create new graphic variable to hold plot contents
  	  plot 	= new(3,graphic)
	  
	;--- Add historical LOCA annual clim file, calculate ensemble mean
	  hist		= a->$hvNames(imetric+3)$
	  hist_ens 	= dim_avg_n_Wrap(hist,2)
	  
	;--- Add future LOCA annual climatology file, calculate ensemble mean
	;    Climatology files created using ann_future.ncl
	  early		= b->$fvNames((3*imetric)+3)$
	  early_diff	= hist
	  early_diff	= early - hist
	  early_ens	= dim_avg_n_Wrap(early_diff,2)
          mid		= b->$fvNames((3*imetric)+4)$
 	  mid_diff	= hist
	  mid_diff	= mid - hist
	  mid_ens	= dim_avg_n_Wrap(mid_diff,2)
	  late		= b->$fvNames((3*imetric)+5)$
	  late_diff	= hist
	  late_diff	= late - hist
	  late_ens	= dim_avg_n_Wrap(late_diff,2)
print("Mean historical "+hvNames(imetric+3)+" = "+avg(hist_ens))
print("Mean "+fvNames((imetric*3)+3)+" "+scn(iscn)+" Early Ens = "+avg(early_ens))
print("Mean "+fvNames((imetric*3)+4)+" "+scn(iscn)+" Mid Ens = "+avg(mid_ens))
print("Mean "+fvNames((imetric*3)+5)+" "+scn(iscn)+" Late Ens = "+avg(late_ens))

          ;--- delete temporary vars
	  delete([/hist,hist_ens,early,early_diff,early_ens,mid,mid_diff,mid_ens/])
	  delete([/late,late_diff,late_ens/])
    end do   ; metrics
  end do   ; scn	 
end
  exit
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


 

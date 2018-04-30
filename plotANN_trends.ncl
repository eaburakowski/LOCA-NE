load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; plotANN_trends.ncl 
; EA Burakowski
; 2017-08-01
; 
; plotANN_trends.ncl plots the annual trend in minimum temperature (tmin), 
;	maximum temperature (tmax), and total precipitation derived from 29 
;	individual CMIP5 LOCAlized Constructed Analogs (LOCA) daily downscaled 
;	data, 1980-2005.  Also includes ensemble mean trend. 
;
; Output includes:
;   .png files for annual trend in each variable:
; 	- ANN_pr_historical_1980-2005_trend.png
;	- ANN_tasmax_historical_1980-2005_trend.png
;	- ANN_tasmin_historical_1980-2005_trend.png
;	

;======================================================================
;                         The main code
;======================================================================

begin

;---- LOCA .nc directory

  dir	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/trends/"
  var	= (/"pr","tasmax","tasmin"/)
  varname = (/"Total Precip., mm/dec","Tmax, C/dec","Tmin, C/dec"/)
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
  

;---- Output directories

  odir	= "/net/home/eos/ean2/LOCA_processed/LOCA_nc/"	


;---- Annual Averages
;
;      Loop over models to calculate annual average
;	tasmax, tasmin, pr
;
;----

;---- Load one historical trend file to get the lat/lon dims

  hfil 	= addfile(dir+"ANN_"+var(0)+"_historical_1980-2005_trend.nc","r")
  trnd	= hfil->ANNTREND
printVarSummary(trnd)
  dims	= dimsizes(trnd)
  lat	= trnd&lat
  lon	= trnd&lon

  ;-- Regional domain boundaries  
    minlat = min(lat)
    maxlat = max(lat)
    minlon = min(lon)
    maxlon = max(lon)
    bndadd = 0.20


;--- Loop over vars to make one plot per variable
  do ivar = 0,dimsizes(var)-1
    print("----  Working on "+var(ivar)+"  -----")
	
	;---- Create new graphic variable to hold plot contents
  	  plot 	= new(30,graphic)
	  
	;--- Add historical annual trend file 
	  a		= addfile(dir+"ANN_"+var(ivar)+"_historical_1980-2005_trend.nc","r")
	  trend		= a->ANNTREND
printVarSummary(trend)
	;--- Convert units from degC/yr to degC/decade
	  degCdec	= trend
	  degCdec	= trend*10

print("max trend degC/dec = "+max(degCdec))
print("min trend degC/dec = "+min(degCdec))


	;----------------------------------------------------------------------
	; Plotting options section
	;----------------------------------------------------------------------

       pltType = "png"                                 ; plot type
       pltDir = "/net/home/eos/ean2/LOCA_processed/LOCA_Figures/"                           
       pltName = "ANN_"+var(ivar)+"Trend_1980-2005"    ; plot name 
       pltPath = pltDir+pltName                        ; plot path

   wks = gsn_open_wks(pltType,pltPath)		; create workstation for ps file
	
	if ivar.eq.0
          gsn_define_colormap(wks,"MPL_PuOr")      ; define color table 
	else
	  gsn_define_colormap(wks,"BlueWhiteOrangeRed")
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
        res@cnMinLevelValF       = round(max(degCdec),3)*-1        ; minimum contour, degC/decade
        res@cnMaxLevelValF       = round(max(degCdec),3)        ; maximum contour, decC/decade 
	if ivar.eq.0
	  res@cnLevelSpacingF	= 20		     ; contour if prcp
        else
	  res@cnLevelSpacingF	 = 0.1		     ; countour interval
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
        
	;--- Loop over model trends and put plot into graphic variable
	do iplot = 0,dimsizes(mdl)-1
	  plot(iplot) = gsn_csm_contour_map(wks,degCdec(:,:,iplot),res)
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
	  pres@lbTitleString		= "mm/decade, 1980-2005"
	else
	pres@lbTitleString		= "~F34~0~F~ C/decade, 1980-2005"
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

  end do  ; variables

end


 

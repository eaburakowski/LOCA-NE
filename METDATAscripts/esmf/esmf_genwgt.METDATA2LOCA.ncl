;======================================================================
; ESMF Regrid:   NLDAS2-PRISM (rectilinear) to LOCA (rectilinear)
;======================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;======================================================================
;                         The main code
;======================================================================

begin

;---Specify interpolation to be used
    method = "bilinear"

;---Input file (Source)
    srcFileDir  = "/net/home/eos/ean2/METDATA_processed/"
    srcFileName = "Obs_ANN_pr_historical_1980-2005.nc"       ; any METDATA file

;---Retrieve a sample variable (used in plot)
    sfile       = addfile(srcFileDir+srcFileName,"r")
    srcVarName  = "ANNCLIM"	; [y | 161] x [x | 225]
    var_in      = sfile->$srcVarName$
printVarSummary(var_in)

;---File that contains the destination file coordinates
    dstFileDir  = "/net/home/eos/ean2/LOCA_processed/LOCA_nc/clim/"
    dstFileName = "ANN_pr_historical_1980-2005.nc"         ; any LOCA processed file

;---open LOCA file that contains the LOCA grid coordinates
    dfile       = addfile(dstFileDir+dstFileName,"r")
    dstvar	= dfile->ANNCLIM(:,:,0)

;---Retrieve coordinates of destination grid from a WRF file
    lat = dfile->lat
    lon = dfile->lon

;---Convert NLDAS2-PRISM to SCRIP convention file
    srcGridName	       = "NLDAS2-PRISM_SCRIP.nc"
    dstGridName	       = "LOCA_SCRIP.nc"

    wgtFileName        = "NLDAS2-PRISM2Rect.nc"

    Opt                = True
    Opt@GridMask       = where(.not.ismissing(var_in),1,0)
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    Opt@Title          = "NLDAS2-PRISM (4km X 4km)"
    
    rectilinear_to_SCRIP(srcGridName,var_in&lat,var_in&lon,Opt)

    delete(Opt)

;---Convert LOCA grid to SCRIP convention file
    Opt		       = True
    Opt@GridMask       = where(.not.ismissing(dstvar),1,0)    
    Opt@ForceOverwrite = True
    Opt@Title	       = "LOCA 16th degree grid"
 
    rectilinear_to_SCRIP(dstGridName,lat,lon,Opt)

    delete(Opt)

;---Generate weights
    wgtFileName        = "NLDAS2-PRISM_4kmM2_to_LOCA."+method+"_wgt.nc"  ; *** This is the wgt file
    Opt		       = True
    Opt@InterpMethod   = method
    Opt@ForceOverwrite = True
  ;;Opt@PrintTimings   = True

    ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtFileName,Opt)
  
    delete(Opt)

end

# LOCA-NE
This repository consists of NCL and MATLAB code for processing, analyzing, and plotting LOcalized Constructed Analogs (LOCA; Pierce et al. 2014) statistically downscaled, gridded climate data.  The gridded climate data includes daily maximum surface air temperature (tasmax), surface air minimum temperature (tasmin), and daily precipitation (pr) at 1/16th degree resolution over the Northeastern United States domain for 29 Coupled Model Intercomparison Project 5 (CMIP5) models.  Simulations include historical (1980-2005), lower emissions rcp4.5 (2006-2099), and higher emissions rcp8.5 (2006-2099). 

This is a repository in progress and updates can be expected throughout Summer 2018.  Code will not run "as-is", but will require customization (i.e., paths, sub-folders). 

Below is the list of scripts and general order of usage in my own processing. 

  ann_hist.ncl
  - calculates the historical (1980-2005) climatological average & linear trend in annual temperature (min & max) and 
  precipitation (total).  Inputs include daily LOCA data for three variables:
     (1) tasmin (surface minimum air temperature)
     (2) tasmax (surface maximum air temperature)
     (3) precipitation (daily total)

  Output includes network Common Data Format (netCDF, or .nc) files for climatological mean and annual trend in each variable 
  for each model:
     - ann_pr_historical_1980-2005.nc
     - ann_tasmax_historical_1980-2005.nc
     - ann_tasmin_historical_1980-2005.nc
     
  ann_future.ncl
  
  seas_hist.ncl
  
  seas_future.ncl
  
  ExtremeTemp_hist.ncl
  
  ExtremeTemp_future.ncl





References: 

Pierce, DW, DR Cayan, and B Thrasher. 2014. Statistical Downscaling Using Localized Constructed Analogs (LOCA). J Hydromet. 15: 2558-2585. doi: 10.1175/JHM-D-14-0082.1. 

# LOCA-NE
This repository consists of NCL and MATLAB code for processing, analyzing, and plotting LOcalized Constructed Analogs (LOCA; Pierce et al. 2014) statistically downscaled, gridded climate data.  The gridded climate data includes daily maximum surface air temperature (tasmax), surface air minimum temperature (tasmin), and daily precipitation (pr) at 1/16th degree resolution over the Northeastern United States domain for 29 Coupled Model Intercomparison Project 5 (CMIP5) models.  Simulations include historical (1980-2005), lower emissions rcp4.5 (2006-2099), and higher emissions rcp8.5 (2006-2099). 

This is a repository in progress and updates can be expected throughout Summer 2018.  Code will not run "as-is", but will require customization (i.e., paths, sub-folders). 

Below is the list of scripts and general order of usage in my own processing. 

  ann_hist.ncl
  - calculates the historical (1980-2005) climatological average & linear trend in annual temperature (min & max) and 
  precipitation (total).  
  - Inputs include daily LOCA data for three variables:
     - tasmin (surface minimum air temperature)
     - tasmax (surface maximum air temperature)
     - precipitation (daily total)
  - Output includes network Common Data Format (netCDF, or .nc) files for climatological mean and annual trend in each
  variable for each model:
     - ann_pr_historical_1980-2005.nc
     - ann_tasmax_historical_1980-2005.nc
     - ann_tasmin_historical_1980-2005.nc
     
  ann_future.ncl
  - calculates the future projected (2006-2099) climatological average annual temperature (min & max) and 
    precipitation (total) for early (2010-2039), mid (2040-2069), and late (2070-2099) century.  
  - Inputs include daily LOCA data for three variables:
     - tasmin (surface minimum air temperature)
     - tasmax (surface maximum air temperature)
     - precipitation (daily total)
  - Output includes network Common Data Format (netCDF, or .nc) files for climatological mean and annual trend in each
  variable for each model:
     - ann_pr_historical_1980-2005.nc
     - ann_tasmax_historical_1980-2005.nc
     - ann_tasmin_historical_1980-2005.nc
  
  seas_hist.ncl
  - calculates the future projected (2006-2099) climatological average seasonal temperature (min & max) and 
    precipitation (total) for early (2010-2039), mid (2040-2069), and late (2070-2099) century.  Seasons include spring 
    (March, April, May; MAM), summer (June, July, and August; JJA), fall (September, October, and November; SON) and winter 
    (December, January, February; DJF). Winter years are associated with the January of a given year (i.e., Winter 2018 is 
    December 2017, January and Feburary 2018).
  - Inputs include daily LOCA data for three variables:
     - tasmin (surface minimum air temperature)
     - tasmax (surface maximum air temperature)
     - precipitation (daily total)
  - Output includes network Common Data Format (netCDF, or .nc) files for climatological mean and seasonal ($SEAS$ = MAM, JJA, 
  SON, or DJF) trend in each variable for each model:
     - $SEAS$_pr_historical_1980-2005.nc
     - $SEAS$_tasmax_historical_1980-2005.nc
     - $SEAS$_tasmin_historical_1980-2005.nc
     
  seas_future.ncl
  
  ExtremeTemp_hist.ncl
  
  ExtremeTemp_future.ncl





References: 

Pierce, DW, DR Cayan, and B Thrasher. 2014. Statistical Downscaling Using Localized Constructed Analogs (LOCA). J Hydromet. 15: 2558-2585. doi: 10.1175/JHM-D-14-0082.1. 

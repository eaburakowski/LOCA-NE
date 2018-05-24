#!/bin/bash
# compute daily tasmean from tasmax and tasmin
# 
# EA Burakowski
# 2018-05-22
# uses nco ncks and ncap2:
# http://nco.sourceforge.net/nco.html 

# define file path
p="/net/nfs/merrimack/raid/Northeast_US_Downscaling_cmip5/"

# list all the tasmax files (will use the suffix to match with tasmin pair)
FILES=${p}tasmax*.nc

# loop over tasmax files and use tasmin pair to calculate tmean in new file
for infiles in ${FILES}
  do
  
  # separate the base filename from file path
  bn=$(basename $infiles)
  echo $bn

  # pull out suffix from tasmax filename
  suff=${bn#*_}
  echo $suff

  # use the suffix to identify the tasmin and tasmax pair
  txFile="tasmax_${suff}"
  tnFile="tasmin_${suff}"
  
  echo $txFile
  echo $tnFile

  # if the input tasmin and tasmax files do not exist, exit
  [ ! -e ${p}${txFile} ] && echo "${p}${txFile} does not exist!" && exit
  [ ! -e ${p}${tnFile} ] && echo "${p}${tnFile} does not exist!" && exit

  # create new filename for tasmean using suffix from input files
  ofil="tasmean_${suff}"
  echo $ofil

  # copy the tasmax file, using the new file name for tasmean
  cp ${p}${txFile} ${p}${ofil}

  # execute nco 'kitchen sink' command to append tasmin to ofil
  ncks -A -v tasmin $p$tnFile $p$ofil

  # ofil now has tasmax and tasmin in the same file. Use ncap2 to calculate average, or tasmean
  ncap2 -O -s 'tasmean=(tasmax+tasmin)/2' ${p}${ofil} ${p}${ofil} 

  # delete tasmin and tasmax variables to reduce file size
  ncks -C -O -x -v tasmin,tasmax ${p}${ofil} ${p}${ofil}

done



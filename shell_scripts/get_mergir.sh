#!/bin/bash

SCRIPT_PATH="/home/disk/eos15/jnug/cold-point-overshoot/python_scripts"
OUT_PATH="/home/disk/eos15/jnug/MERGIR"
LINKS_PATH="/home/disk/eos15/jnug/MERGIR/mergir_links_lists"

MONTHS="DJF"
YEAR=2009
declare -a RegionArr=("AMZ" "IOS" "SPC1" "SPC2" "ECP") # DJF
# declare -a RegionArr=("AFR" "IOE" "WPC" "ECP") # JJA

for region in "${RegionArr[@]}"; do
    url_list=$LINKS_PATH/subset_$MONTHS$YEAR"_"$region".txt"
    subdir_out=$OUT_PATH/$region"_"$MONTHS$YEAR
    mkdir -p $subdir_out
    cd $subdir_out
    
    # get all of the individual files into a subdirectory
    wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies --content-disposition -i $url_list
    
    # concat into one file 
    out_file=$OUT_PATH/MERGIR_Tb_4km_$MONTHS$YEAR"_"$region".nc4"
    cdo cat $subdir_out/*.nc4.nc4 $out_file
    
done

#!/bin/bash

module load cdo

declare -a SeasonArr=("DJF") #"JJA" )
FILE_PATH="/work/bb1153/b380887/big_obs_climo"
OUT_PATH="/scratch/b/b380887"

for season in "${SeasonArr[@]}"; do
    if [ "$season" == "DJF" ] ; then
        declare -a RegionArr=(SPC1 SPC2 ECP) #AMZ IOS SPC1 SPC2 ECP)
    elif [ "$season" == "JJA" ] ; then
        declare -a RegionArr=(AFR WPC ECP IOE)
    fi
    echo $season
    
    for region in "${RegionArr[@]}"; do
        echo $region
        echo $out_file
        out_file=$OUT_PATH/MERGIR_Tb_4km_$season"2007-2010_"$region".nc4"
        cdo cat $FILE_PATH/$season/MERGIR_Tb_4km_$season"*"$region".nc4" $out_file
    done
done

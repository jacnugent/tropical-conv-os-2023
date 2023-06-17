#!/bin/bash

set -evx # verbose messages and crash message

module load cdo

declare -a SeasonArr=("DJF" "JJA" )
FILE_PATH="/work/bb1153/b380887/big_obs_climo"
OUT_PATH="/scratch/b/b380887"

for season in "${SeasonArr[@]}"; do
    if [ "$season" == "DJF" ] ; then
        declare -a RegionArr=(AMZ SPC1 SPC2 IOS ECP)
    elif [ "$season" == "JJA" ] ; then
        declare -a RegionArr=(AFR WPC IOE ECP)
    fi
    echo $season
    
    for region in "${RegionArr[@]}"; do
        echo $region
        out_file_tb=$OUT_PATH/MERGIR_Tb_4km_$season"2007-2010_"$region".nc4"
        out_file_cpt=$OUT_PATH/ERA5_cpT_reindexed_$season"2007-2010_"$region".nc"
        cdo cat $FILE_PATH/$season/MERGIR_Tb_4km_$season"*"$region".nc4" $out_file_tb
        cdo cat $FILE_PATH/$season/ERA5_cpT_reindexed_$season"*"$region".nc" $out_file_cpt
    done
done

#!/bin/bash

module load python3
source activate /home/b/b380887/.conda/envs/d2env

SEASON="DJF"

SCRIPT_PATH="/home/b/b380887/cold-point-overshoot/python_scripts"
FILE_PATH="/work/bb1153/b380887/big_obs_climo/"$SEASON
PICKLE_PATH="/home/b/b380887/cold-point-overshoot/pickle_files/climo_tb_cp_hists/"
# OUT_PATH="/scratch/b/b380887/"
OUT_PATH=$FILE_PATH

declare -a YearArr=(2007 2008 2009 2010)

if [ "$SEASON" == "DJF" ] ; then
    declare -a RegionArr=(SPC) #AMZ IOS SPC1 SPC2 SPC ECP)
elif [ "$SEASON" == "JJA" ] ; then
    declare -a RegionArr=(AFR WPC ECP IOE)
fi
    

for year in "${YearArr[@]}"; do
    echo $year
    for region in "${RegionArr[@]}"; do
        if [ "$region" != "SPC" ] ; then
            # only get regridded files for SPC1 & SPC2
            echo "Getting regridded cold point temperature file..."
            python $SCRIPT_PATH/cold_point_reindex.py -y $year -r $region -m $SEASON -f $FILE_PATH -o $OUT_PATH
        fi
        echo "Get the histogram dictionary..."
        python $SCRIPT_PATH/biv_hist.py -y $year -r $region -m $SEASON -f $FILE_PATH -p $PICKLE_PATH
        echo "...done with "$region
    done
done
    

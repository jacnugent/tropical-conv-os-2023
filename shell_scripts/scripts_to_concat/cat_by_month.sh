#!/bin/bash

set -evx # verbose messages and crash message

module load python3
source activate /home/b/b380887/.conda/envs/d2env

YEAR=2010
SEASON="JJA"

SCRIPT_PATH="/home/b/b380887/cold-point-overshoot/python_scripts"
FILE_PATH="/scratch/b/b380887/ERA5_ml/"
OUT_PATH="/work/bb1153/b380887/big_obs_climo/"$SEASON

if [ "$SEASON" == "DJF" ] ; then
    declare -a RegionArr=(SPC1 SPC2) #AMZ IOS SPC1 SPC2 ECP)
    months="Dec Jan Feb"
elif [ "$SEASON" == "JJA" ] ; then
    declare -a RegionArr=(AFR IOE WPC ECP)
    months="Jun Jul Aug"
fi

python $SCRIPT_PATH/cat_files.py -y $YEAR -f $FILE_PATH -o $OUT_PATH -p "ERA5_T_0.25deg_ml_12-20km" -e ".nc" -r "ECP" -m $months 

for region in "${RegionArr[@]}"; do
    echo $region":"
    python $SCRIPT_PATH/cat_files.py -y $YEAR -f $FILE_PATH -o $OUT_PATH -p "ERA5_T_0.25deg_ml_12-20km" -e ".nc" -r $region -m $months 
    python $SCRIPT_PATH/cat_files.py -y $YEAR -f $FILE_PATH -o $OUT_PATH -p "ERA5_zg_0.25deg_ml"  -e ".nc" -r $region -m $months 
done

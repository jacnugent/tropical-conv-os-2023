#!/bin/bash

module load python3
source activate /home/b/b380887/.conda/envs/d2env

YEAR=2008
SEASON="JJA"

MIN_HEIGHT=12000
MAX_HEIGHT=25080

SCRIPT_PATH="/home/b/b380887/cold-point-overshoot/python_scripts"
FILE_PATH="/scratch/b/b380887/DARDAR_v3/"$SEASON
OUT_PATH="/scratch/b/b380887/DARDAR_v3"


# set coord arrays to the correct season
if [ "$SEASON" == "DJF" ] ; then
    declare -A CoordDict
    CoordDict["AMZ"]="-72 -47 -30 0"
    CoordDict["SPC1"]="165 180 -20 -5" # W half
    CoordDict["SPC2"]="-180 -145 -20 -5" # E half
    CoordDict["ECP"]="-150 -100 0 15"
    CoordDict["IOS"]="50 100 -15 0"
else
    declare -A CoordDict
    CoordDict["AFR"]="-7 35 0 18"
    CoordDict["WPC"]="130 180 0 15"
    CoordDict["ECP"]="-150 -100 0 15"
    CoordDict["IOE"]="53 95 -12 6"
fi

for region in "${!CoordDict[@]}"; do
    echo $region":"
    coords_str="${CoordDict[$region]}"
    coords=($coords_str)
    lon0=${coords[0]}
    lon1=${coords[1]}
    lat0=${coords[2]}
    lat1=${coords[3]}
    
    python3 $SCRIPT_PATH/process_dardar_v3.py -y $YEAR -lt $lat0 $lat1 -ln $lon0 $lon1 -zmin $MIN_HEIGHT -zmax $MAX_HEIGHT -r $region -m $SEASON -fp $FILE_PATH -op $OUT_PATH

done

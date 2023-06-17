#!/bin/bash

module load python3

YEAR=2007
SEASON="JJA"

SCRIPT_PATH="/home/b/b380887/cold-point-overshoot/python_scripts"
OUT_PATH="/scratch/b/b380887/ERA5_ml/tq_regional"
WORK_PATH="/work/bb1153/b380887/big_obs_climo"
# GT_PATH=$OUT_PATH
GT_PATH="/scratch/b/b380887/ERA5_ml"

# set month/coord arrays to the correct season
if [ "$SEASON" == "DJF" ] ; then
    declare -a MonthArr=(Dec Jan Feb)
    months="Dec Jan Feb"
    declare -a RegionArr=(AMZ IOS SPC1 SPC2 ECP)
    declare -A CoordDict
    CoordDict["AMZ"]="-72,-47,-30,0"
    CoordDict["SPC1"]="165,180,-20,-5" # W half
    CoordDict["SPC2"]="-180,-145,-20,-5" # E half
    CoordDict["ECP"]="-150,-100,0,15"
    CoordDict["IOS"]="50,100,-15,0"
else
    declare -a MonthArr=(Jun Jul Aug)
    months="Jun Jul Aug"
    declare -a RegionArr=(AFR IOE WPC ECP)
    declare -A CoordDict
    CoordDict["AFR"]="-7,35,0,18"
    CoordDict["WPC"]="130,180,0,15"
    CoordDict["ECP"]="-150,-100,0,15"
    CoordDict["IOE"]="53,95,-12,6"
fi


# process renalysis data for all months/regions in that season
for month in "${MonthArr[@]}"; do
    if [ "$month" == "Dec" ] ; then
        # year=expr $YEAR - 1
        year=$(($YEAR - 1))
    else
        year=$YEAR
    fi
    echo $month $year
    temp_gt=$GT_PATH/ERA5_T_0.25deg_ml_12-20km_$month$year"_GT.nc"
    # tq_gt=$OUT_PATH/tq_$month$year"_GT.grib"
    tq_gt=$GT_PATH/tq_$month$year"_GT.grib"
    zlnsp_gt=$OUT_PATH/zlnsp_$month$year"_GT.grib"
    
    for region in "${!CoordDict[@]}"; do
        echo $region":"
        # echo "Subsetting global tropics files..."
        coords="${CoordDict[$region]}"
        temp_out=$OUT_PATH/ERA5_T_0.25deg_ml_12-20km_$month$year"_"$region".nc"
        tq_out=$OUT_PATH/tq_$month$year"_"$region".grib"
        zlnsp_out=$OUT_PATH/zlnsp_$month$year"_"$region".grib"
        
        # if [ ! -f $temp_out ] ; then
        #     cdo sellonlatbox,$coords $temp_gt $temp_out
        # fi
        if [ ! -f $tq_out ] ; then
            cdo sellonlatbox,$coords $tq_gt $tq_out
        fi
        # if [ ! -f $zlnsp_out ] ; then
        #     cdo sellonlatbox,$coords $zlnsp_gt $zlnsp_out
        # fi
        
#         echo "Computing geopotential..."
#         zout=$OUT_PATH/z_out_$month$year"_"$region".grib"
#         if [ ! -f $zout ] ; then
#             python3 $SCRIPT_PATH/compute_geopotential_on_ml.py -o $zout $tq_out $zlnsp_out
#         fi
        
#         echo "Converting to netcdf..."
#         zout_nc=$OUT_PATH/ERA5_zg_0.25deg_ml_$month$year"_"$region".nc"
#         if [ ! -f $zout_nc ] ; then
#             grib_to_netcdf -o $zout_nc $zout
#         fi

        echo "Converting tq file to netcdf..."
        tq_out_nc=$OUT_PATH/ERA5_tq_0.25deg_ml_$month$year"_"$region".nc"
        if [ ! -f $tq_out_nc ] ; then
            grib_to_netcdf -o $tq_out_nc $tq_out
        fi
        echo "..."$region" done"

    done
    
done

# concatenate files by month for each region
# source deactivate
source activate /home/b/b380887/.conda/envs/d2env
echo "Concatenating monthly files into seasonal files..."
# out_path_szn=$WORK_PATH/$SEASON
out_path_szn=$OUT_PATH
for region in "${RegionArr[@]}"; do
    echo $region":"
    # python $SCRIPT_PATH/cat_files.py -y $YEAR -f $OUT_PATH -o $out_path_szn -p "ERA5_T_0.25deg_ml_12-20km" -e ".nc" -r $region -m $months 
    # python $SCRIPT_PATH/cat_files.py -y $YEAR -f $OUT_PATH -o $out_path_szn -p "ERA5_zg_0.25deg_ml"  -e ".nc" -r $region -m $months 
    python $SCRIPT_PATH/cat_files.py -y $YEAR -f $OUT_PATH -o $out_path_szn -p "ERA5_tq-t_only_0.25deg_ml"  -e ".nc" -r $region -m $months -v "t"    
    python $SCRIPT_PATH/cat_files.py -y $YEAR -f $OUT_PATH -o $out_path_szn -p "ERA5_tq-q_only_0.25deg_ml"  -e ".nc" -r $region -m $months -v "q" 
done


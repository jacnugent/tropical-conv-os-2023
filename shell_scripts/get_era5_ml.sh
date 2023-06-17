#!/bin/bash

module load python3

SCRIPT_PATH="/home/b/b380887/cold-point-overshoot/python_scripts"
OUT_PATH="/scratch/b/b380887/ERA5_ml/"

COORDS="30 -180 -30 180"
REGION="GT"
YEAR=2010
# declare -a MonthArr=("06" "07" "08")
# declare -a MonthArr=("12" "01" "02")
declare -a MonthArr=("02") # just Feb for 2010

source activate /home/b/b380887/.conda/envs/era5env

for month in "${MonthArr[@]}"; do
    if [ "$month" == "12" ] ; then
        # year=expr $YEAR - 1
        year=$(($YEAR - 1))
    else
        year=$YEAR
    fi
    if [[ "$month" == @("04"|"06"|"09"|"11") ]]; then
        endday="30"
    elif [ "$month" == "02" ]; then
        endday="28"
    else
        endday="31"
    fi
    date1=$year"-"$month"-01"
    date2=$year"-"$month"-"$endday
    echo "$date1 $date2"
    
    # get string name of the month
    if [ "$month" == "12" ]; then
        monthstr="Dec"
    elif [ "$month" == "01" ]; then
        monthstr="Jan"
    elif [ "$month" == "02" ]; then
        monthstr="Feb"
    elif [ "$month" == "06" ]; then
        monthstr="Jun"
    elif [ "$month" == "07" ]; then
        monthstr="Jul"
    elif [ "$month" == "08" ]; then
        monthstr="Aug"
    else
        monthstr=$month
    fi
    
    # only download it if it doesn't exist already
    # if [[ ! -f $OUT_PATH/tq_zlnsp_$monthstr$year"_GT.grib" ]] ; then
    #     python $SCRIPT_PATH/get_era5_climo_ml.py -s $date1 -e $date2 -c $COORDS -o $OUT_PATH"/" -r $REGION
    # fi
    
    # TEMP - just get tq for JJA 2007/2009 and DJF 2009/2010
    python $SCRIPT_PATH/get_era5_climo_ml.py -s $date1 -e $date2 -c $COORDS -o $OUT_PATH/ -r $REGION -v temp_and_qv
    
#     # split the files
#     comb_in=$GT_PATH/tq_zlnsp_$month"_GT.grib"
#     tq_out=$OUT_PATH/tq_$month$year"_GT.grib"
#     zlnsp_out=$OUT_PATH/zlnsp_$month$year"_GT.grib"
#     echo "getting tq:"
#     cdo selname,t,q $comb_in $tq_out
#     echo "getting zlnsp:"
#     cdo selname,z,lnsp $comb_in $zlnsp_out
    
#     # remove global combined file (if everything worked)
#     if [[ -f $tq_out && -f $zlnsp_out ]] ; then
#         rm $comb_in
#     fi
    
done

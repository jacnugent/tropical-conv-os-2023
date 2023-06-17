#!/bin/bash

module load lftp

HOST=ftp.icare.univ-lille1.fr
# USER=
# PASSWORD=

YEAR=2008
SEASON="JJA"

DAR_PATH=/SPACEBORNE/CLOUDSAT/DARDAR-CLOUD.v3.10
OUT_PATH=/scratch/b/b380887/DARDAR_v3
declare -a DayArr=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31)

if [ "$SEASON" == "DJF" ] ; then
    declare -a MonthArr=(12 01 02)
elif [ "$SEASON" == "JJA" ] ; then
    declare -a MonthArr=(06 07 08)
fi

{
    echo user $USER $PASSWORD
    echo !mkdir -p $OUT_PATH/$SEASON
    echo lcd $OUT_PATH/$SEASON
    
    for month in ${MonthArr[@]}; do
    
        if [ "$month" == "12" ] ; then
            year=$(($YEAR - 1))
        else
            year=$YEAR
        fi
        
        for day in ${DayArr[@]}; do
            subdir=$DAR_PATH/$year/$year"_"$month"_"$day
            echo cd $subdir
            echo !echo $subdir
            echo mget *.nc
        done
        
    done

} | lftp $HOST

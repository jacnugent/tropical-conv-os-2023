#!/bin/bash

module load python3
source activate /home/b/b380887/.conda/envs/d2env

declare -a YearArr=(2007 2008 2009 2010) 
declare -a SeasonArr=(DJF JJA)

SCRIPT_PATH="/home/b/b380887/cold-point-overshoot/python_scripts"
FILE_PATH="/work/bb1153/b380887/big_obs_climo"
PICKLE_PATH="/home/b/b380887/cold-point-overshoot/pickle_files/binned_by_tb/"
PLOT_PATH="/home/b/b380887/cold-point-overshoot/plots/obs_paper/standard_plots/"
mkdir -p $PICKLE_PATH
mkdir -p $PLOT_PATH

for season in "${SeasonArr[@]}"; do
    echo $season
    if [ "$season" == "DJF" ] ; then
        declare -a RegionArr=(AMZ IOS SPC1 SPC2 ECP)
    elif [ "$season" == "JJA" ] ; then
        declare -a RegionArr=(AFR WPC ECP IOE)
    fi
    
    for region in "${RegionArr[@]}"; do
        # echo "Regridding data onto cold point-relative levels..."
        # python $SCRIPT_PATH/regrid_data_cp.py -y $YEAR -r $region -m $season -f $FILE_PATH -o $FILE_PATH/$season
        # echo "Binning by brightness temperature and saving standard plots..."
        # python $SCRIPT_PATH/bin_obs_overshoot.py -y $YEAR -r $region -m $season -f $FILE_PATH/$season -p $PICKLE_PATH -s $PLOT_PATH
        # echo "...done with "$region"."
        
        # TEMP - year loop
        for year in "${YearArr[@]}"; do
            echo $year
            echo "Binning by brightness temperature and saving standard plots..."
            python $SCRIPT_PATH/bin_obs_overshoot.py -y $year -r $region -m $season -f $FILE_PATH/$season -p $PICKLE_PATH"Tb_minus_cpT/" -s $PLOT_PATH"Tb_minus_cpT/" 
        done
        
        echo "...done with "$region"."
    done
    
done

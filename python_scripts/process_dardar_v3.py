""" 
process_dardar_v3.py

Module to process DARDAR-CLOUD version 3 global files
and subset for a given region.

** must specify the year!**

------------------------------------------------------
usage: process_dardar_v3.py [-h] -y YEAR [-ln LONS [LONS ...]] [-lt LATS [LATS ...]] [-m MONTHS] [-zmin MIN_HEIGHT] [-zmax MAX_HEIGHT] [-fp FILE_PATH] [-op OUT_PATH] [-r REGION]

optional arguments:
  -h, --help            show this help message and exit
  -y YEAR, --year YEAR  Year
  -ln LONS [LONS ...], --lons LONS [LONS ...]
                        Longitudes (-180, 180) in the form "west east", e.g., 143 153
  -lt LATS [LATS ...], --lats LATS [LATS ...]
                        Latitudes in the form "south north", e.g., -5 5
  -m MONTHS, --months MONTHS
                        Month(s) name for the output file
  -zmin MIN_HEIGHT, --min_height MIN_HEIGHT
                        Minimum height in m
  -zmax MAX_HEIGHT, --max_height MAX_HEIGHT
                        Maximum height in m
  -fp FILE_PATH, --file_path FILE_PATH
                        Input path for files
  -op OUT_PATH, --out_path OUT_PATH
                        Output path to save processed files
  -r REGION, --region REGION
                        Region name for the output file
----------------------------------------------------

1. To get the DARDAR files, first edit the script
"get_dardar_v3_via_ftp.sh" to specify the year, 
which downloads the files for some months into FILE_PATH 
using ftp. (or modify for different months, etc.)
2. Run this script in the command line to process
the files.
3. Remove the global .nc files from FILE_PATH (they're huge). 
4. Repeat as needed.

Default is to get only heights 12-22 km.

Code to process modified from Adam Sokol's notebook (https://github.com/jacnugent/dyamond-ms/blob/master/
python_scripts/old/iwp_hists_for_dyamond_From_Adam.ipynb)
"""

import sys
import glob
import argparse

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from datetime import datetime, timedelta


FILE_PATH = "/home/disk/eos15/jnug/DARDAR_v3/"
OUT_PATH = "home/disk/eos15/jnug/DARDAR_v3/"
LATS = [-5, 5]
LONS = [143, 153]
REGION = "TWP"
MONTHS = "Jan"
MIN_HEIGHT = 12000
MAX_HEIGHT = 22000


def in_region(file_name, latlims, lonlims):
    """ 
    Give a DARDAR file path and lat/lon lims. Function tells you if the file has data
    from within those limits
    """    
    ds = xr.open_dataset(file_name)
    lat = ds.latitude.values
    lon = ds.longitude.values
    in_box = (lat>=latlims[0]) & (lat<=latlims[1]) & (lon>=lonlims[0]) & (lon<=lonlims[1])
    if np.nansum(in_box)>0:
        return True
    else:
        return False
    

def get_iwp_ds(file_name, latlims, lonlims, zmin, zmax):
    """
    Makes a dataset with only certain variables for the DARDAR file if profiles are within 
    the lat/lon limits. Also converts the time to datetime
    """
    ds_all = xr.open_dataset(file_name)
    varlist=['latitude','longitude','time','iwc','effective_radius','instrument_flag','ln_effective_radius_error','ln_iwc_error','land_water_mask','day_night_flag']
    ds = ds_all[varlist]
    
    # subset for region
    ds_sub = ds.where((ds.latitude>=latlims[0])&(ds.latitude<=latlims[1])&(ds.longitude>=lonlims[0])&(ds.longitude<=lonlims[1]),drop=True)
    
    # add IWP variable (BEFORE subsetting heights - so you get total-column IWP)
    ds_sub['iwp'] = (('time'), np.nansum(ds_sub.iwc.values*60,axis=1))
    ds_sub['iwp'].attrs = {'units': 'kg m-2',
                           'long_name': 'Calculated Ice Water Path (total column)'
                          }

    # subset for height limits
    ds_sub = ds_sub.where((ds_sub.height<=zmax)&(ds_sub.height>=zmin), drop=True)
    
    # change lat & lon to be coordinates
    ds_sub = ds_sub.rename({"latitude": "lat", "longitude": "lon"})
    ds_sub = ds_sub.assign_coords({"lat": ds_sub.lat, "lon": ds_sub.lon})

    return ds_sub


def parse_args():
    """ Parse command-line arguments
    """
    parser = argparse.ArgumentParser()
    
    # required
    parser.add_argument("-y", "--year", help="Year", type=int, required=True)
    
    # optional
    parser.add_argument("-ln", "--lons", nargs="+", help="Longitudes (-180, 180) in the form \"west east\", e.g., 143 153", type=float)
    parser.add_argument("-lt", "--lats", nargs="+", help="Latitudes in the form \"south north\", e.g., -5 5", type=float)
    parser.add_argument("-m", "--months", default=MONTHS, help="Month(s) name for the output file")
    parser.add_argument("-zmin", "--min_height", default=MIN_HEIGHT, help="Minimum height in m", type=float)
    parser.add_argument("-zmax", "--max_height", default=MAX_HEIGHT, help="Maximum height in m", type=float)
    parser.add_argument("-fp", "--file_path", default=FILE_PATH, help="Input path for files")
    parser.add_argument("-op", "--out_path", help="Output path to save processed files")
    parser.add_argument("-r", "--region", default=REGION, help="Region name for the output file")

    args = parser.parse_args()

    return args


def main():
    """ Subset lat/lon and save as netcdf. Prints progress.
    """   
    args = parse_args()
    print(args)

    lats = args.lats
    lons = args.lons
    region = args.region
    year = args.year
    months = args.months
    min_height = args.min_height
    max_height = args.max_height
    file_path = args.file_path
    out_path = args.out_path
    
    # if coords are not specified, use the defaults
    if lons is None:
        lons = LONS
    if lats is None:
        lats = LATS

    # if an out_path is not specified, use "year" subdirectory inside of the file_path
    if out_path is None:
        out_path = file_path + str(year)

    # default is that all target files are stored in file_path; if not,
    # assumes they are stored in subdirectories for each day, YYYY_MM_DD
    if len(sorted(glob.glob(file_path + "/DARDAR-CLOUD*.nc"))) > 0:
        files = sorted(glob.glob(file_path + "/DARDAR-CLOUD*.nc"))
    else:
        try:
            month_num = datetime.strptime(months, "%b").strftime("%m")
        except:
            month_num = datetime.strptime(months, "%B").strftime("%m")
        files = glob.glob(file_path + "/{y}_{m}_??/DARDAR-CLOUD*.nc".format(y=year, m=month_num))
    print("{} total files".format(len(files)))

    reg_files = [file for file in files if in_region(file, lats, lons)]
    print("{} files in the region".format(len(reg_files)))
    
    ds_list = [get_iwp_ds(file, latlims=lats, lonlims=lons, zmin=min_height, zmax=max_height) for file in reg_files]
    print("Subset each .nc global file for select variables and that region")
    
    ds_cat = xr.concat(ds_list, dim="time").sortby("time")
    print("Concatenated the list of datasets into one (by time)")

    out_name = out_path + "/DARDAR-v3_iwc_{m}{y}_{r}.nc".format(m=months, y=year, r=region)
    ds_cat.to_netcdf(out_name)
    print("Dataset saved to", out_name)


# run in the command line
if __name__ == "__main__":
    main()
    

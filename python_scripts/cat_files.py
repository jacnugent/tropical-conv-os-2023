"""
cat_files.py

--------------------------------------------------------------
usage: cat_files.py [-h] -y YEAR -fp FILE_PATH -r REGIONS [REGIONS ...] -m MONTHS
                    [MONTHS ...] [-op OUT_PATH] [-cl] [-np] [-t]

optional arguments:
  -h, --help            show this help message and exit
  -y YEAR, --year YEAR  Year
  -f FILE_PATH, --file_path FILE_PATH
                        Input path for files
  -r REGIONS [REGIONS ...], --regions REGIONS [REGIONS ...]
                        Region abbreviation(s)
  -m MONTHS [MONTHS ...], --months MONTHS [MONTHS ...]
                        Month abbreviations
  -o OUT_PATH, --out_path OUT_PATH
                        Output path to save processed files
  -c, --cat_lon        Concatenate all regions by longitude instead of time
                        (default=False); must input region names in order from west to
                        east
  -n, --no_print       Do not print the list(s) of files that will be concatenated
  -t, --time_dim_only   The only file dimension is time (so will concatenate "longitude"
                        by the time dim
  -p, --file_pre        Beginning of file(s) to concatenate
  -e, --file_ext        Extension of file(s) to concatenate
  -v, --variable_name   Name of the variable to take from the dataset (if you only want
                         to do one)
----------------------------------------------------------------

Concatenates monthly files into one season or concatenates split
regional files by longitude (e.g., you have a W + E half). 
Name convention is that files must include three-letter
month abbreviation, year, and then region abbreviation as follows:
    
    *_Jun2009_AFR.[ext] for June 2009 in an Africa/"AFR" region 
    with file extension .[ext]
    
Uses xarray (sometimes cdo cat doesn't work with the time coordinate
in these files).
Default is to output to the same directory the files were in.
"""
import argparse
import glob

import xarray as xr


def parse_args():
    """ Parse command-line arguments
    """
    parser = argparse.ArgumentParser()
    
    # required
    parser.add_argument("-y", "--year", help="Year", type=int, required=True)
    parser.add_argument("-f", "--file_path", help="Input path for files", required=True)
    parser.add_argument("-r", "--regions", nargs="+", help="Region abbreviation(s)", required=True)
    parser.add_argument("-m", "--months", nargs="+", help="Month abbreviations", required=True)

    
    # optional
    parser.add_argument("-o", "--out_path", help="Output path to save processed files")
    parser.add_argument('-c', "--cat_lon", help="Concatenate all regions by longitude instead of time (default=False); must input region names in order from west to east", action='store_true')
    parser.add_argument("-n", "--no_print", help="Do not print the list(s) of files that will be concatenated", action='store_true')
    parser.add_argument("-t", "--time_dim_only", help="The only file dimension is time (so will concatenate \"longitude\" by the time dim", action='store_true')
    parser.add_argument("-p", "--file_pre", help="Beginning of file(s) to concatenate")
    parser.add_argument("-e", "--file_ext", help="Beginning of file(s) to concatenate")
    parser.add_argument("-v", "--variable_name", help="Name of the variable to take from the dataset (if you only want to do one")

    return parser.parse_args()


def cat_lon(regions, month, year, file_pre, file_ext, file_path, out_path, time_dim_only, no_print, return_ds=False,
            variable_name=None):
    """ Concatenate files for one month and region (e.g., SPC1 + SPC2) by longitude
    """
    ds_list = []
    for i, region in enumerate(regions):
        if file_pre is None:
            pre = "*"
        else:
            pre = file_pre
        if file_ext is None:
            ext = ".*"
        else:
            ext = file_ext
        filename = glob.glob(file_path + "{p}_{m}{y}_{r}{e}".format(p=pre, m=month, y=year, r=regions, e=ext))[0]
        if not no_print:
            print(filename)
        if variable_name is not None:
            ds_list.append(xr.open_dataset(filename))
        else:
            ds_list.append(xr.open_dataset(filename)[variable_name])

        # assumes file naming convention is the same for all
        if i == 0:
            if file_pre is None:
                pre_list = filename.split("/")[-1].split("_")[:-2]
                file_pre = [pre_list[i] + "_" + pre_list[i+1] for i in range(len(pre_list) - 1)][0] 
            if file_ext is None:
                file_ext = filename.split(".")[-1]
                
    # do the concatenating
    if time_dim_only:
        ds_cat = xr.concat(ds_list, dim="time").sortby("time")
    else:
        ds_cat = xr.concat(ds_list, dim="lon")
    cat_region = region[:-1] # cut last character (ordinal of the region list)
    if variable_name is not None:
        file_pre = file_pre + "-{v}_only".format(v=variable_name)
    out_name = file_pre + "_{m}{y}_{c}.".format(m=months[0], y=year, c=cat_region) + file_ext

    if out_path[-1] != "/":
        out_path = out_path + "/"
    ds_cat.to_netcdf(out_path + out_name)
    print("Dataset saved to", out_path + out_name)
    
    if return_ds:
        return ds_cat

    
def cat_time(region, months, year, file_pre, file_ext, file_path, out_path,  time_dim_only, no_print, return_ds=False,
             variable_name=None):
    """ Concatenate files for one season and region by time
    """    
    ds_list = []
    for j, month in enumerate(months):
        if month == "Dec" and j == 0:
            file_year = year - 1
        else:
            file_year = year
        if file_pre is None:
            pre = "*"
        else:
            pre = file_pre
        if file_ext is None:
            ext = ".*"
        else:
            ext = file_ext
        filename = glob.glob(file_path + "{p}_{m}{y}_{r}{e}".format(p=pre, m=month, y=file_year, r=region, e=ext))[0]            
        if not no_print:
            print(filename)
        if variable_name is None:
            ds_list.append(xr.open_dataset(filename))
        else:
            ds_list.append(xr.open_dataset(filename)[variable_name])

    if file_pre is None:
        pre_list = filename.split("/")[-1].split("_")[:-2]
        file_pre = [pre_list[i] + "_" + pre_list[i+1] for i in range(len(pre_list) - 1)][0] 
    if file_ext is None:
        file_ext = filename.split(".")[-1]

    season = ""
    for month in months:
        season += month[0]
    ds_cat = xr.concat(ds_list, dim="time")
    if variable_name is not None:
        file_pre = file_pre + "-{v}_only".format(v=variable_name)
    out_name = file_pre + "_{s}{y}_{r}".format(s=season, y=year, r=region) + file_ext

    if out_path[-1] != "/":
        out_path = out_path + "/"
    ds_cat.to_netcdf(out_path + out_name)
    print("Dataset saved to", out_path + out_name)
        
    if return_ds:
        return ds_cat
    
    
def main():
    """ 
    Concatenate files (monthly --> season or W/E --> single region by lon). Saves the concatenated file.
    """   
    args = parse_args()
    print(args)

    year = args.year
    file_path = args.file_path
    regions = args.regions
    months = args.months
    out_path = args.out_path
    cat_lon = args.cat_lon
    no_print = args.no_print
    time_dim_only = args.time_dim_only
    file_pre = args.file_pre
    file_ext = args.file_ext
    variable_name = args.variable_name
    
    if out_path is None:
        out_path = file_path
    if file_path[-1] != "/":
        file_path = file_path + "/"
    
    if cat_lon:
        cat_lon(regions, months[0], year, file_pre, file_ext, file_path, out_path, 
                time_dim_only, no_print, variable_name=variable_name)
    else:
        for i, region in enumerate(regions):
            cat_time(region, months, year, file_pre, file_ext, file_path, out_path, 
                     time_dim_only, no_print, variable_name=variable_name)
        
    
if __name__ == "__main__":
    main()   

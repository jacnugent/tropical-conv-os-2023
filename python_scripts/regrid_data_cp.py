"""
regrid_data_cp.py

Regrid MERGIR brightness temperature (geostationary) and ERA5 cold point 
height/temperature (lat-lon) onto the 1-D DARDAR grid (sun-synchronous).
Gets the ice water content, effective radius, and instrument flag at levels
relative to the cold point (+/- 1000 m, +/- 500m, at cold point) and saves as
a netcdf.

---------------

usage: regrid_data_cp.py [-h] -y YEAR -r REGION -m MONTHS [-f FILE_PATH]
                         [-o OUT_PATH] [-p]

optional arguments:
  -h, --help            show this help message and exit
  -y YEAR, --year YEAR  Year
  -r REGION, --region REGION
                        Region abbreviation
  -m MONTHS, --months MONTHS
                        Month/season abbreviation (e.g., DJF)
  -f FILE_PATH, --file_path FILE_PATH
                        Input path for files
  -o OUT_PATH, --out_path OUT_PATH
                        Output path to save files
  -p, --dask_progress   Print dask progress bars (default=False)
  
-------------

"""
import argparse 
import dask
import xarray as xr
import numpy as np
from dask.diagnostics import ProgressBar


FILE_PATH = "/work/bb1153/b380887/big_obs_climo/"
OUT_PATH = "/scratch/b/b380887/temp_for_dardar/"
DZ_DARDAR = 60. # vertical grid spacing for DARDAR (m)


def parse_args():
    """ Parse command-line arguments
    """
    parser = argparse.ArgumentParser()
    
    # required
    parser.add_argument("-y", "--year", help="Year", type=int, required=True)
    parser.add_argument("-r", "--region", help="Region abbreviation", required=True)
    parser.add_argument("-m", "--months", help="Month/season abbreviation (e.g., DJF)", required=True)
    
    # optional
    parser.add_argument("-f", "--file_path", help="Input path for files", default=FILE_PATH)
    parser.add_argument("-o", "--out_path", help="Output path to save files", default=OUT_PATH)
    parser.add_argument("-p", "--dask_progress", help="Print dask progress bars (default=False)", action='store_true')

    return parser.parse_args()

        
def format_string(string, months, year, region):
    """ Function to format file paths (as strings)
    """
    string_fmt = string.format(m=months, y=year, r=region)
    return string_fmt


def get_files(months, year, region, file_path=FILE_PATH, out_path=OUT_PATH):
    """ 
    Open the files for all variables needed (ERA5 temp/height,
    DARDAR iwc, and MERGIR brightness temperature).
    ** change file name conventions here if necessary**
    """
    tb = xr.open_dataset(file_path + format_string("{m}/MERGIR_Tb_4km_{m}{y}_{r}.nc4", months, year, region))["Tb"]
    ds_iwc = xr.open_dataset(file_path + format_string("{m}/DARDAR-v3_iwc_{m}{y}_{r}.nc", months, year, region))
    temp_era5 = xr.open_dataset(file_path + format_string("{m}/ERA5_T_0.25deg_ml_12-20km_{m}{y}_{r}.nc", months, year, region))["t"]
    z_era5 = xr.open_dataset(file_path + format_string("{m}/ERA5_zg_0.25deg_ml_{m}{y}_{r}.nc", months, year, region))["z"].sortby("time")/9.81
    
    # check the times worked
    assert(np.array_equal(temp_era5.time.values, z_era5.time.values))

    # convert longitude from 0-360 to -180-180 if needed (on ERA5 height; changed in z computation)
    if (z_era5.longitude.min() > 180) or (z_era5.longitude.max() > 180):
        z_era5.coords['longitude'] = (z_era5.coords['longitude'] + 180) % 360 - 180
    
    # if one is -180 and the other is 180 at some point:
    if not (np.array_equal(temp_era5.longitude.values, z_era5.longitude.values)):
        inds_diff = np.where((temp_era5.longitude.values - z_era5.longitude.values) > 0)[0]
        if len(inds_diff) > 1:
            raise Exception("ERA5 longitude values are different for temperature & height!")
        else:
            ind = inds_diff[0]
            if abs(temp_era5.longitude.values[ind]) == 180 and abs(z_era5.longitude.values[ind]) == 180:
                new_lon = z_era5.longitude.values
                new_lon[ind] = temp_era5.longitude[ind]
                z_era5 = z_era5.assign_coords({"longitude": new_lon})
                assert(np.array_equal(temp_era5.longitude.values, z_era5.longitude.values))
            else:
                raise Exception("ERA5 longitude values are different for temperature & height!")

        
    # flip the ERA5 coordinates to be increasing if they're not
    if temp_era5.latitude[0] > temp_era5.latitude[-1]:
        temp_era5 = temp_era5.isel(latitude=slice(None, None, -1))
    if temp_era5.longitude[0] > temp_era5.longitude[-1]:
        temp_era5 = temp_era5.isel(longitude=slice(None, None, -1))
    if z_era5.latitude[0] > z_era5.latitude[-1]:
        z_era5 = z_era5.isel(latitude=slice(None, None, -1))
    if z_era5.longitude[0] > z_era5.longitude[-1]:
        z_era5 = z_era5.isel(longitude=slice(None, None, -1))

    return tb, ds_iwc, temp_era5, z_era5

    
    
def regrid_tb_cpz(iwc, tb, cpz):
    """ 
    Regrid brightness temperature & cold point height data arrays
    onto the DARDAR IWC 1-D grid. Returns the regridded data arrays.
    """
    iwc2d = iwc.isel(height=0)

    if ((tb.lon[-1] - tb.lon[0]) < 0):
        if tb.lon.max() > 180:
            # if you have 360 --> 0 coords, change it so they're all increasing
            cpz.coords['longitude'] = ((cpz.coords['longitude'] - 180) % 360) -180
            tb.coords['lon'] = ((tb.coords['lon'] - 180) % 360) - 180
            iwc2d.coords['lon'] = ((iwc2d.coords['lon'] - 180) % 360) - 180
            iwc.coords['lon'] = ((iwc.coords['lon'] - 180) % 360) - 180
        else:
            # if you have 180 --> -180 coords, change it so they're all increasing
            cpz.coords['longitude'] = (cpz.coords['longitude'] % 360)
            tb.coords['lon'] = (tb.coords['lon'] % 360)
            iwc2d.coords['lon'] = (iwc2d.coords['lon'] % 360)
            iwc.coords['lon'] = (iwc.coords['lon'] % 360)
    
    tb_dar = tb.sel(time=iwc.time, lat=iwc2d.lat, lon=iwc2d.lon, method="nearest")
    cpz_dar = cpz.sel(time=iwc.time, latitude=iwc2d.lat, longitude=iwc2d.lon, method="nearest")

    return tb_dar, cpz_dar


def regrid_tb_only(iwc, tb):
    """ 
    Regrid brightness temperature data array
    onto the DARDAR IWC 1-D grid. Returns the regridded data array.
    """
    iwc2d = iwc.isel(height=0)
    if (iwc2d.lon[-1] - iwc2d.lon[0]) > 0:
        tb_dar = tb.sel(time=iwc.time, lat=iwc2d.lat, lon=iwc2d.lon, method="nearest")
    else:
        tb.coords['lon'] = (tb.coords['lon'] % 360)
        iwc2d.coords['lon'] = (iwc2d.coords['lon'] % 360)
        tb_dar = tb.sel(time=iwc.time, lat=iwc2d.lat, lon=iwc2d.lon, method="nearest")

    return tb_dar


def get_cp_relative_dict(ds_iwc, cpz_dar, dz=DZ_DARDAR):
    """ 
    Gets data arrays of regridded IWC, effective radius, & instrument flag at each 
    cold point-relative level (at cp and +/- 1000m) and saves as a dataset. 
    Returns three dictionaries with
    - Keys: 
        ["iwc_a1000", "iwc_a500", "iwc_cp", "iwc_b500", "iwc_b1000"]
        ["Re_a1000", "Re_a500", "Re_cp", "Re_b500", "Re_b1000"]
        ["iflag_a1000", "iflag_a500", "iflag_cp", "iflag_b500", "iflag_b1000"]
    - Values:
        [(cp+1000m), (cp+500m), (cp), (cp-500m), (cp-1000m)],
        etc.

    * cpz_dar = the data array of ERA5 cold point heights at each
    DARDAR time (get it from interp_dardar).
    * iwc_rg = the data array of regridded IWC (at all possible cp-relative heights)
    """
    iwc = ds_iwc["iwc"]
    Re = ds_iwc["effective_radius"] # TEMP
    flag = ds_iwc["instrument_flag"] # TEMP
    
    iwc_keys = ["iwc_a1000", "iwc_a500", "iwc_cp", "iwc_b500", "iwc_b1000"]
    Re_keys = ["Re_a1000", "Re_a500", "Re_cp", "Re_b500", "Re_b1000"] # TEMP
    iflag_keys = ["iflag_a1000", "iflag_a500", "iflag_cp", "iflag_b500", "iflag_b1000"] # TEMP
    
    dar_heights = iwc.height
    cpz_dar = cpz_dar.assign_coords({"time": iwc.time})

    # get number of indices for a difference of 1000 and 500 m
    # make it negative because DARDAR is ordered top-down
    dind_1000 = -1*int(np.round(1000/dz))
    dind_500 = -1*int(np.round(500/dz))
    
    # get the smallest height you can have with the data given; might need to
    # cut out the lowest heights in some cases 
    cp_lower_thresh = dar_heights.min().values - dz*dind_1000
    
    # cut out the lowest values by setting nans to an absurdly high value so
    # argmin won't pick them (throws all nan slice encountered error if you leave them in)
    if cpz_dar.min() < cp_lower_thresh:
        cpz_dar = cpz_dar.where(cpz_dar >= cp_lower_thresh, drop=True).fillna(99999)

    cp_inds_dar = np.abs(cpz_dar - dar_heights).argmin(dim="height") 
    print("got the cp inds")

    # if the cp height +1000m index is above the dardar max height
    # (which can happen if there was no ice there anywhere), 
    # we need to drop it! (set to nan) - so get the min height
    # that it could realistically be (so the indices don't get set 
    # back to negative)
    a1000_min = cpz_dar.min().values - dz*dind_1000
    a500_min = cpz_dar.min().values - dz*dind_500
    print("got the min values - so you can drop if the indices flipped")
    
    # rename the heights so they don't conflict when you make the dataset
    # and do the checks for +1000 and +500
    iwc_a1000 = iwc.isel(height=(cp_inds_dar + dind_1000))
    Re_a1000 = Re.isel(height=(cp_inds_dar + dind_1000))
    iflag_a1000 = flag.isel(height=(cp_inds_dar + dind_1000))
    print("+1000 m done")
    iwc_a500 = iwc.isel(height=(cp_inds_dar + dind_500))
    Re_a500 = Re.isel(height=(cp_inds_dar + dind_500))
    iflag_a500 = flag.isel(height=(cp_inds_dar + dind_500))
    print("+500 m done")
    iwc_cp = iwc.isel(height=cp_inds_dar)  
    Re_cp = Re.isel(height=cp_inds_dar)  
    iflag_cp = flag.isel(height=cp_inds_dar)  
    print("cold point done")
    iwc_b500 = iwc.isel(height=(cp_inds_dar - dind_500))
    Re_b500 = Re.isel(height=(cp_inds_dar - dind_500))
    iflag_b500 = flag.isel(height=(cp_inds_dar - dind_500))
    print("-500 m done")
    iwc_b1000 = iwc.isel(height=(cp_inds_dar - dind_1000))
    Re_b1000 = Re.isel(height=(cp_inds_dar - dind_1000))
    iflag_b1000 = flag.isel(height=(cp_inds_dar - dind_1000))
    print("-1000 m done")
        
    iwc_a1000_chunked = iwc_a1000.chunk({"time": 1000})
    iwc_a1000 = iwc_a1000_chunked.where(iwc_a1000_chunked.height >= a1000_min).compute()
    Re_a1000_chunked = Re_a1000.chunk({"time": 1000})
    Re_a1000 = Re_a1000_chunked.where(Re_a1000_chunked.height >= a1000_min).compute()
    iflag_a1000_chunked = iflag_a1000.chunk({"time": 1000})
    iflag_a1000 = iflag_a1000_chunked.where(iflag_a1000_chunked.height >= a1000_min).compute()
    print("got +1000m where it's bigger than min")
    
    iwc_a500_chunked = iwc_a500.chunk({"time": 1000})
    iwc_a500 = iwc_a500_chunked.where(iwc_a500_chunked.height >= a500_min).compute()
    Re_a500_chunked = Re_a500.chunk({"time": 1000})
    Re_a500 = Re_a500_chunked.where(Re_a500_chunked.height >= a500_min).compute()
    iflag_a500_chunked = iflag_a500.chunk({"time": 1000})
    iflag_a500 = iflag_a500_chunked.where(iflag_a500_chunked.height >= a500_min).compute()
    print("got +500 where it's bigger than min")
    
    print("iwc at levels is done; need to rename before saving")
    iwc_cp = iwc_cp.rename({"height": "height_cp"})
    iwc_b500 = iwc_b500.rename({"height": "height_b500"})
    iwc_b1000 = iwc_b1000.rename({"height": "height_b1000"})
    iwc_a1000 = iwc_a1000.rename({"height": "height_a1000"})
    iwc_a500 = iwc_a500.rename({"height": "height_a500"})
    
    Re_cp = Re_cp.rename({"height": "height_cp"})
    Re_b500 = Re_b500.rename({"height": "height_b500"})
    Re_b1000 = Re_b1000.rename({"height": "height_b1000"})
    Re_a1000 = Re_a1000.rename({"height": "height_a1000"})
    Re_a500 = Re_a500.rename({"height": "height_a500"})
    
    iflag_cp = iflag_cp.rename({"height": "height_cp"})
    iflag_b500 = iflag_b500.rename({"height": "height_b500"})
    iflag_b1000 = iflag_b1000.rename({"height": "height_b1000"})
    iflag_a1000 = iflag_a1000.rename({"height": "height_a1000"})
    iflag_a500 = iflag_a500.rename({"height": "height_a500"})
    
    iwc_da_list = [iwc_a1000, iwc_a500, iwc_cp, iwc_b500, iwc_b1000]
    Re_da_list = [Re_a1000, Re_a500, Re_cp, Re_b500, Re_b1000]
    iflag_da_list = [iflag_a1000, iflag_a500, iflag_cp, iflag_b500, iflag_b1000]
    
    cp_rel_dict_list = [dict(zip(iwc_keys, iwc_da_list)), dict(zip(Re_keys, Re_da_list)), dict(zip(iflag_keys, iflag_da_list))]
    return cp_rel_dict_list


def save_cp_relative(months, year, region, cp_rel_dict_list, out_path, return_ds=False):
    """ 
    Save a dataset of the cold point-relative IWCs (at cold point
    and +/- 1000m and 500m). Returns the dataset if requested.
    """
    out_file_name_cp = out_path + format_string("DARDAR_cp_relative_iwc_{m}{y}_{r}.nc", months, year, region)
    iwc_cp_ds = xr.Dataset(cp_rel_dict_list[0])
    Re_cp_ds = xr.Dataset(cp_rel_dict_list[1])
    iflag_cp_ds = xr.Dataset(cp_rel_dict_list[2])
    
    ds_all = xr.merge([iwc_cp_ds, Re_cp_ds, iflag_cp_ds])
    ds_all.to_netcdf(out_file_name_cp)
    
    if return_ds:
        return ds_all


def compute_cp_rel_dardar(months, year, region, ds_iwc, cpz_dar, dz=DZ_DARDAR, out_path=OUT_PATH, return_dict=False):
    """ 
    Runs get_cp_relative_dict and save_cp_relative, which
    does the following:
        Takes the IWC file (Dataset), gets data arrays at the 
        regridded cold point-relative levels (cp and +/- 1000m, 500m) and 
        saves as a netcdf. 
    """
    cp_rel_dict_list = get_cp_relative_dict(ds_iwc, cpz_dar, dz)
    print("Dictionary done")
    save_cp_relative(months, year, region, cp_rel_dict_list, out_path, return_ds=False)
    print("Saved!")
    
    if return_dict:
        return cp_rel_dict_list
    

def get_era5_cp(temp, z, return_temp=False):
    """ 
    Return data arrays with cold point height (& temperature, if
    return_temp=True) for ERA5. 
    """      
    cp_inds = temp.argmin(dim="level")
    cpT = temp.isel(level=cp_inds)
    cpz = z.sel(level=temp.level).isel(level=cp_inds)
    
    if return_temp:
        return cpz, cpT
    else: 
        return cpz
    
    
    
def main(dz=DZ_DARDAR):
    """ 
    To run in the command line! Saves two netcdf files to out_path.
    
    Regrids the MERGIR Tb and ERA5 cold point height files, subsets
    the native IWC file at regridded cold point-relative levels, and saves 
    resulting dataset to a netcdf file.
    """
    args = parse_args()
    print(args)

    year = args.year
    region = args.region
    months = args.months
    file_path = args.file_path
    out_path = args.out_path
    dask_progress = args.dask_progress

    if file_path[-1] != "/":
        file_path = file_path + "/"
    if out_path[-1] != "/":
        out_path = out_path + "/"

    if dask_progress:
        pbar = ProgressBar()
        pbar.register()
        
    # get files, cold point heights, & regrid to match IWC 1D grid
    tb, ds_iwc, temp_era5, z_era5 = get_files(months, year, region, file_path, out_path)
    cpz = get_era5_cp(temp_era5, z_era5, return_temp=False)
    tb_dar, cpz_dar = regrid_tb_cpz(ds_iwc["iwc"], tb, cpz)
    
    # get & save the dictionary
    compute_cp_rel_dardar(months, year, region, ds_iwc, cpz_dar, dz, out_path)


if __name__ == "__main__":
    main()

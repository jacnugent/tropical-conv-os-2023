"""
interp_dardar_cp.py

Regrids DARDAR IWC to a (regular) lat-lon
grid at ERA5 cold point-relative levels 
(+/- 1000m and 500m). Generates 2 netcdf files:
    * Regridded DARDAR IWC (at all possible cold point (+/- 1km) levels
    * Regridded DARDAR IWC at five levels: 
      cold point and cold point +/- 1000 m and 500 m

Assumes you want to regrid/interpolate it to match
the MERGIR brightness temperature ("tb") grid,
but should work for any! 

Ignores any cold point values that are below (~12 km + 1km) so
that (cold point - 1000m) will always be in the dataset (this excludes
a very small percentage of profiles in only some seasons/regions, 
e.g., ~0.3% for a 10x10 region in South America in JJA 2015).

Assumes files are named like this (so change if 
they're not): 
    * ERA5_T_0.25deg_ml_12-20km_Feb2015_TWP.nc
    * ERA5_zg_0.25deg_ml_Feb2015_TWP.nc
    * MERGIR_tb_4km_Feb2015_TWP.nc4
    * DARDAR-v3_iwc_Feb2015_TWP.nc

Input months, year, and region when you run, e.g.:
    python3 interp_dardar_cp.py Feb 2015 TWP
        or
    interp_dardar.interp_dardar(month="Feb", year=2015, region="TWP")
    interp_dardar.compute_cp_rel_dardar(month="Feb", year=2015, region="TWP")

If running in the command line, can specify alternate file and out paths:
    python3 interp_dardar_cp.py Feb 2015 TWP (file_path) (out_path)

"""
import os.path
import sys
import gc
import dask
import glob

import xarray as xr
import pandas as pd
import numpy as np
import datetime as dt

from scipy.interpolate import griddata
from dask.diagnostics import ProgressBar
from numba import jit
from pathlib import Path


FILE_PATH = "/work/bb1153/b380887/10x10/obs_climo/"
OUT_PATH = "/scratch/b/b380887/temp_for_dardar/"
DZ_DARDAR = 60. # vertical grid spacing for DARDAR (m)
IWC_MIN = 1e-7 # minimum IWC value to count as ice (in kg/m3)
DASK_PROGRESS = False


def format_string(string, months, year, region):
    """ Function to format file paths (as strings)
    """
    string_fmt = string.format(m=months, y=year, r=region)
    return string_fmt


def closest_time(times, options=None):
    """ 
    Convert the dardar times to the closest half-hour time
    on that day/month/year. 
    Possible times are hours 0-23 and mins 0 and 30
    """
    if options is None:
        options = [
            dt.time(hour=h, minute=m) 
            for h in np.arange(24)
            for m in [0, 30]
        ]
        
    closest_times = [[]]*len(times)
    
    for i, time in enumerate(times):
        time_dt = pd.to_datetime(time)
        options_dt = [dt.datetime(year=time_dt.year, month=time_dt.month,
                                  day=time_dt.day, hour=option.hour,
                                  minute=option.minute) 
                      for option in options]
    
        closest_times[i] = min(options_dt, key=lambda x: abs(x - time_dt))
    

    return closest_times


def get_era5_cp(temp, z, iwc=None, return_temp=False):
    """ 
    Return data arrays with cold point height (& temperature, if
    return_temp=True) for ERA5. If you pass in a data array for iwc, returns the 
    cold points only at times when you have any DARDAR data.
    """      
    cp_inds = temp.argmin(dim="level")
    cpT = temp.isel(level=cp_inds)
    cpz = z.isel(level=cp_inds)
    
    if iwc is not None:
        # at the hours when you have DARDAR output
        hr_min = iwc["time.hour"].min().values
        hr_max = iwc["time.hour"].max().values
        cpz_hd = cpz.where((cpz["time.hour"] < hr_max+1) & \
                           (cpz["time.hour"] >= hr_min), drop=True
                          )
        
        # get list of the time in the ERA5 dataset that's closest to the dardar time
        hrly_options = [dt.time(hour=h, minute=0) for h in np.arange(hr_min, hr_max+1)]
        dar_times_hrly = closest_time(iwc.time.values, options=hrly_options)
    
        # get the cold point heights at each dardar time
        cpz_dar = cpz_hd.sel(time=sorted(list(set(dar_times_hrly))))

    if iwc is None:
        return cpz, cpT
    else: 
        return cpz


def get_files(months, year, region, get_era5, file_path=FILE_PATH, out_path=OUT_PATH):
    """ 
    Open the files for all variables needed (ERA5 temp/height,
    DARDAR iwc, and MERGIR brightness temperature). 
    
    ** change file name conventions here if necessary**
    
    """
    tb = xr.open_dataset(file_path + format_string("MERGIR_tb_4km_{m}{y}_{r}.nc4", months, year, region))["Tb"]
    iwc = xr.open_dataset(file_path + format_string("DARDAR-v3_iwc_{m}{y}_{r}.nc", months, year, region))["iwc"]
    
    # for DARDAR: check that it doesn't include the first day of the next month!
    # check by making sure first and last data point have the same month
    if iwc["time.month"][-1] != iwc["time.month"][0]:
        month_num = int(dt.datetime.strptime(months, "%b").strftime("%m"))
        iwc = iwc.where(iwc["time.month"]==month_num, drop=True)
        
    if get_era5:
        temp_era5 = xr.open_dataset(file_path + format_string("ERA5_T_0.25deg_ml_12-20km_{m}{y}_{r}.nc", months, year, region))["t"]

        # have to sort ERA5 height by time (gets jumbled with the compute/to_netcdf)
        z_era5 = xr.open_dataset(file_path + format_string("ERA5_zg_0.25deg_ml_{m}{y}_{r}.nc", months, year, region))["z"].sel(level=temp_era5.level).sortby("time")/9.81
        assert(np.array_equal(temp_era5.time.values, z_era5.time.values))

        # convert longitude from 0-360 to -180-180 if needed (on ERA5 height; changed in z computation)
        if (z_era5.longitude.min() > 180) or (z_era5.longitude.max() > 180):
            z_era5.coords['longitude'] = (z_era5.coords['longitude'] + 180) % 360 - 180
        assert(np.array_equal(temp_era5.longitude.values, z_era5.longitude.values))

        # flip the ERA5 coordinates to be increasing if they're not
        if temp_era5.latitude[0] > temp_era5.latitude[-1]:
            temp_era5 = temp_era5.isel(latitude=slice(None, None, -1))
        if temp_era5.longitude[0] > temp_era5.longitude[-1]:
            temp_era5 = temp_era5.isel(longitude=slice(None, None, -1))
        if z_era5.latitude[0] > z_era5.latitude[-1]:
            z_era5 = z_era5.isel(latitude=slice(None, None, -1))
        if z_era5.longitude[0] > z_era5.longitude[-1]:
            z_era5 = z_era5.isel(longitude=slice(None, None, -1))

        return tb, iwc, temp_era5, z_era5
    else:
        return tb, iwc


def get_cp_height_slice(cpz_dar, dz):
    """
    Get the height range of (min ERA5 cold point - 1km) to
    the max height in DARDAR subset(here, ~22km). 
    The -dz at the end ensures the slice takes ALL values (right index 
    not inclusive). 
    Returns a slice object.
    """
    cpz_min = cpz_dar.min().values - 1000
    
    return slice(None, cpz_min-dz)


def get_time_dict_coords(iwc):
    """ 
    For each IWC time step, get the time in the MERGIR dataset (half hourly) 
    that's closest to that time step.
    
    Returns: 
        * dictionary mapping unique DARDAR times (matching the MERGIR dataset)
          to an index of that time 
          (e.g., {"2020-01-01T05:30": 0, "2020-01-06T07:00": 1})
        * list (of length len(iwc.time)) of indices of the unique times
          (e.g., [0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 2])
    """
    dar_times_30m = closest_time(iwc.time.values)
    
    # get the unique times
    dar_times_30m_unique = sorted(list(set(dar_times_30m)))
    
    # get the dictionary and coordinate list
    time_dict = dict(zip(dar_times_30m_unique, np.arange(len(dar_times_30m_unique))))
    time_coords = [time_dict[time] for time in dar_times_30m]
    
    return time_dict, time_coords


@jit(nopython=True)
def populate_target_arr(arr, target_arr, target_lats, target_lons, arr_lats, arr_lons):
    """
    Populates the target array with iwc values at the closest lat and lon
    in the target grid to the original grid. Uses numba to speed up the for loop.
    """
    for i in range(np.shape(arr)[0]):
        lat_ind = (np.abs(target_lats - arr_lats[i])).argmin()
        lon_ind = (np.abs(target_lons - arr_lons[i])).argmin()
        target_arr[i, :, lat_ind, lon_ind] = arr[i, :]
        
    return target_arr
    

def get_target_da(iwc_nz, tb, tc):
    """ 
    Returns the interpolated IWC at some time that matches
    a time step in the MERGIR dataset. Takes the mean over time, 
    so if a cell has multiple IWC values within that time step
    (e.g., 7:09, 7:15, 7:20 for 7:30), it averages them.
    
    tc is the time coordinate index for that time.
    """
    # get only the values in that time coordinate group
    iwc_nz_tc = iwc_nz.where(iwc_nz.tc == tc, drop=True)
    print("# time steps in this tc group:", len(iwc_nz_tc.time))
    del iwc_nz
    gc.collect()
    
    # lats and lons of the target grid
    target_lats = tb.lat.values
    target_lons = tb.lon.values
    
    # --- if there are > 1200 time steps, split it up! (may throw memory error otherwise...) ---
    # (this is faster if you don't loop over a list of the two arrays)
    tc_time_length = len(iwc_nz_tc.time)
    if tc_time_length > 1200:
        if tc_time_length > 1800:
            print("Splitting up large tc group into 3: {a}, {b}, and {c}".format(a=600, b=600, c=tc_time_length-1200))
        else:
            print("Splitting up large tc group into 2: {a} and {b}".format(a=600, b=tc_time_length-600))
        
        # -- part 1 --
        iwc_nz_tc1 = iwc_nz_tc.isel(time=slice(None, 600))
        target_arr1 = np.zeros((len(iwc_nz_tc1.time), len(iwc_nz_tc1.height),\
                               len(tb.lat), len(tb.lon)))
        dar_lats1 = iwc_nz_tc1.isel(height=0).lat.values
        dar_lons1 = iwc_nz_tc1.isel(height=0).lon.values
        target_arr1 = populate_target_arr(iwc_nz_tc1.values, target_arr1, target_lats,\
                                          target_lons, dar_lats1, dar_lons1)
        target_da1_chunked = xr.DataArray(
            target_arr1, 
            dims=["time", "height", "lat", "lon"],
            coords={"time": iwc_nz_tc1.time, "height": iwc_nz_tc1.height,
                    "lat": tb.lat, "lon": tb.lon}
        ).chunk({"lat": 20, "lon": 20})
        del target_arr1 
        gc.collect() 
        target_da_mean1 = target_da1_chunked.where(target_da1_chunked > 0).mean(dim="time").compute()
        
        # -- part 2 --
        if tc_time_length > 1800:
            iwc_nz_tc2 = iwc_nz_tc.isel(time=slice(600, 1800))
        else:
            iwc_nz_tc2 = iwc_nz_tc.isel(time=slice(600, None))
        target_arr2 = np.zeros((len(iwc_nz_tc2.time), len(iwc_nz_tc2.height),\
                               len(tb.lat), len(tb.lon)))
        dar_lats2 = iwc_nz_tc2.isel(height=0).lat.values
        dar_lons2 = iwc_nz_tc2.isel(height=0).lon.values
        target_arr2 = populate_target_arr(iwc_nz_tc2.values, target_arr2, target_lats,\
                                          target_lons, dar_lats2, dar_lons2)
        target_da2_chunked = xr.DataArray(
            target_arr2, 
            dims=["time", "height", "lat", "lon"],
            coords={"time": iwc_nz_tc2.time, "height": iwc_nz_tc2.height,
                    "lat": tb.lat, "lon": tb.lon}
        ).chunk(chunks={"lat": 20, "lon": 20})
        del target_arr2
        gc.collect()
        target_da_mean2 = target_da2_chunked.where(target_da2_chunked > 0).mean(dim="time").compute()
        
        # -- part 3, if needed --
        if tc_time_length > 1800:
            iwc_nz_tc3 = iwc_nz_tc.isel(time=slice(1800, None))
            target_arr3 = np.zeros((len(iwc_nz_tc3.time), len(iwc_nz_tc3.height),\
                                   len(tb.lat), len(tb.lon)))
            dar_lats3 = iwc_nz_tc3.isel(height=0).lat.values
            dar_lons3 = iwc_nz_tc3.isel(height=0).lon.values
            target_arr3 = populate_target_arr(iwc_nz_tc3.values, target_arr3, target_lats,\
                                              target_lons, dar_lats3, dar_lons3)
            target_da3_chunked = xr.DataArray(
                target_arr3, 
                dims=["time", "height", "lat", "lon"],
                coords={"time": iwc_nz_tc3.time, "height": iwc_nz_tc3.height,
                        "lat": tb.lat, "lon": tb.lon}
            ).chunk(chunks={"lat": 20, "lon": 20})
            del target_arr3
            gc.collect()
            target_da_mean3 = target_da3_chunked.where(target_da3_chunked > 0).mean(dim="time").compute()
        
        # now take the mean together
        if tc_time_length > 1800:
            target_da_cat = xr.concat([target_da_mean1, target_da_mean2, target_da_mean3], dim="time").chunk({"lat": 20, "lon": 20})
        else:
            target_da_cat = xr.concat([target_da_mean1, target_da_mean2], dim="time").chunk({"lat": 20, "lon": 20})
        target_da_mean = target_da_cat.where(target_da_cat > 0).mean(dim="time").compute()
        del target_da_mean1
        del target_da_mean2
        if tc_time_length > 1800:
            del target_da_mean3
        gc.collect()
 
    # --- if there are <=1200 time steps ---
    else:
        target_arr = np.zeros((len(iwc_nz_tc.time), len(iwc_nz_tc.height),\
                               len(tb.lat), len(tb.lon)))
        
        # get the target array with the numba function
        dar_lats = iwc_nz_tc.isel(height=0).lat.values
        dar_lons = iwc_nz_tc.isel(height=0).lon.values
        target_arr = populate_target_arr(iwc_nz_tc.values, target_arr, target_lats,\
                                         target_lons, dar_lats, dar_lons)
        
        # make into an xarray
        target_da = xr.DataArray(
            target_arr, 
            dims=["time", "height", "lat", "lon"],
            coords={"time": iwc_nz_tc.time, "height": iwc_nz_tc.height,
                    "lat": tb.lat, "lon": tb.lon}
        )

        # target_da_chunked = target_da.chunk({"lat": 10, "lon": 10})
        target_da_chunked = target_da.chunk({"lat": 20, "lon": 20})
        del target_da
        gc.collect()

        # before you take mean, set zeros (empty values left in the array) to nans
        target_da_mean = target_da_chunked.where(target_da_chunked > 0).mean(dim="time").compute()

    return target_da_mean


def save_target_nc(iwc_nz, tb, tc, months, year, region, time_dict, out_path,# return_da=False,
                   print_filename=True):
    """ 
    Runs get_target_da and saves the resulting regridded 
    data array as a netcdf (with variable name "iwc"). Returns the target
    data array if return_da=True. Prints the output file name if print_filename=True
    (default).
    """
    target_da = get_target_da(iwc_nz, tb, tc=tc)

    if tc < 10:
        tc_str = "00" + str(tc)
    elif tc < 100:
        tc_str = "0" + str(tc)
    else:
        tc_str = str(tc)
     
    # don't need to hang on to the dataset
    # but try = make it a variable so we can delete it once saved??
    out_file_name = "tc_{t}_{m}_{y}_{r}.nc".format(t=tc_str, m=months, y=year, r=region)
    xr.Dataset({"iwc": target_da}, attrs = {"time_coord": tc,
                                            "time": str(list(time_dict.keys())[tc])}
              ).to_netcdf(out_path + out_file_name)
    print("saved the target netcdf successfully")
    if print_filename:
        print("saved to", (out_path + out_file_name), "\n")


def interp_dardar(months=None, year=None, region=None, dz=DZ_DARDAR, 
                  iwc_min=IWC_MIN, file_path=FILE_PATH, out_path=OUT_PATH,
                  cpz_dar=None): #return_cpz_dar=False):
    """ 
    Runs all of the necessary functions to interpolate a DARDAR
    (time, height) file to a regular grid for some set of months in 
    one year & region. Returns the cpz_dar data array if return_cpz_dar=True
    (need this when running the main() function). If save_cpz_dar = True (default),
    saves the DARDAR cold points (from cpz_dar) as a netcdf (variable name "cpz_dar").
    """
    if months is None or year is None or region is None:
        raise Exception("Must input months, year, AND region")
                        
    # read in the files (only need ERA5 if you need to calculate the cold point heights)
    if cpz_dar is None:
        tb, iwc_full, temp_era5, z_era5 = get_files(months, year, region, True, file_path, out_path)
        cpz_dar = get_era5_cp(temp_era5, z_era5, iwc_full)
        print("Saving the ERA5 cold point heights file")
        out_name_cpz_dar = format_string("ERA5_cpz_DARDAR_times_{m}{y}_{r}.nc", months, year, region)
        cpz_dar_ds = xr.Dataset({"cpz_dar": cpz_dar})
        cpz_dar_ds.to_netcdf(out_path + out_name_cpz_dar)
    else:
        tb, iwc_full = get_files(months, year, region, False, file_path, out_path)
        if type(cpz_dar) == xr.core.dataset.Dataset:
            cpz_dar_ds = cpz_dar
            cpz_dar = cpz_dar_ds["cpz_dar"] # get the data array only
        save_cpz_dar = False # override this argument if it was True
                        
    # get iwc at the possible cold point heights (+/- 1000 m)
    height_slice = get_cp_height_slice(cpz_dar, dz)
    iwc = iwc_full.sel(height=height_slice)
    iwc_full = None
    
    # get the time coords dictionary & assign the coords to the iwc time dim
    time_dict, time_coords = get_time_dict_coords(iwc)
    iwc_with_tc = iwc.assign_coords(tc=("time", time_coords))
    iwc = None
    
    # restrict to only the times that have ice according to iwc_min threshold
    # and get the list of time coordinate indices for what's left
    iwc_nz = iwc_with_tc.where(iwc_with_tc > iwc_min, drop=True)
    tc_groups = list(set(iwc_nz.tc.values))
    print("Total number of time coords groups: {}".format(len(tc_groups)))
    
    # loop through the groups, getting the regridded dataset (all possible cp-relative
    # levels) & saves it to a netcdf
    da_list = [[]]*len(tc_groups)
    for i, tc in enumerate(tc_groups):
        # # skip the ones you already did!
        # if (months == "Jun") and (tc <= 107):
        #     print("skipping {}; already done".format(tc))
        #     continue
        # else:
        print("Regridding IWC files for time coords group {}...".format(tc))
        save_target_nc(iwc_nz, tb, tc, months, year, region, time_dict, out_path)
        
    # if return_cpz_dar:
    try:
        return cpz_dar_ds
    except:
        return cpz_dar
    

def concat_save_rg_files(months, year, region, out_path=OUT_PATH, 
                         return_iwc_rg=False, iwc_min=IWC_MIN):
    """ 
    Concatenate the regridded files for each time group into one 
    data array and save it to a netcdf (variable name: "iwc").
    """
    # read in the files and concatenate into one regridded file
    group_files = list(sorted(glob.glob(out_path + "*tc*{m}*{y}*{r}*.nc".format(m=months, y=year, r=region))))
    da_list = [open_rg_add_time(file) for file in group_files]
    da_cat = xr.concat(da_list, dim="time")
    del da_list
    gc.collect()
    
    # save the regridded data array  
    da_cat.attrs = {
        "iwc_min": iwc_min,
        "units": "kg/m3",
        "long_name": "regridded ice water content"
    }
    out_file_name = out_path + format_string("DARDAR_regridded_iwc_{m}{y}_{r}.nc", months, year, region)
    iwc_rg = xr.Dataset({"iwc": da_cat})
    iwc_rg.to_netcdf(out_file_name)
        
    if return_iwc_rg:
        return iwc_rg


def open_rg_add_time(file_name):
    """ 
    Open a regridded IWC file and give it a 1-D time dimension with the
    time from the attributes as coordinate. Returns that data array.
    """
    ds = xr.open_dataset(file_name)
    time = pd.to_datetime(ds.time)
    da_exp = ds["iwc"].expand_dims({"time": [time]})
    return da_exp


def get_cp_relative_dict(iwc_rg, cpz_dar, dz=DZ_DARDAR):
    """ 
    Gets data arrays of regridded IWC at each cold point-relative level 
    (at cp and +/- 1000m) and saves as a dataset. Returns a dictionary with
    - Keys: 
        ["iwc_a1000", "iwc_a500", "iwc_cp", "iwc_b500", "iwc_b1000"]
    - Values:
        [(cp+1000m), (cp+500m), (cp), (cp-500m), (cp-1000m)]

    * cpz_dar = the data array of ERA5 cold point heights at each
    DARDAR time (get it from interp_dardar).
    * iwc_rg = the data array of regridded IWC (at all possible cp-relative heights)
    """
    keys = ["iwc_a1000", "iwc_a500", "iwc_cp", "iwc_b500", "iwc_b1000"]
    dar_heights = iwc_rg.height

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
     
    print("reindexing...")
    # reindex the cold point heights to match dardar
    cpz_ri = cpz_dar.rename({"latitude": "lat", "longitude": "lon"}).reindex_like(iwc_rg, method="nearest")
    print("reindexing done; getting indices...") 
    cp_inds_dar = np.abs(cpz_ri - dar_heights).argmin(dim="height") 
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
    iwc_a1000 = iwc_rg.isel(height=(cp_inds_dar + dind_1000))
    print("+1000 m done")
    iwc_a500 = iwc_rg.isel(height=(cp_inds_dar + dind_500))
    print("+500 m done")
    iwc_cp = iwc_rg.isel(height=cp_inds_dar)   
    print("cold point done")
    iwc_b500 = iwc_rg.isel(height=(cp_inds_dar - dind_500))
    print("-500 m done")
    iwc_b1000 = iwc_rg.isel(height=(cp_inds_dar - dind_1000))
    print("-1000 m done")
    
    iwc_a1000_chunked = iwc_a1000.chunk({"time": 6})
    iwc_a1000 = iwc_a1000_chunked.where(iwc_a1000_chunked.height >= a1000_min).compute()
    print("got +1000m where it's bigger than min")
    
    iwc_a500_chunked = iwc_a500.chunk({"time": 6})
    iwc_a500 = iwc_a500_chunked.where(iwc_a500_chunked.height >= a500_min).compute()
    print("got +500 where it's bigger than min")
    
    print("iwc at levels is done; need to rename before saving")
    iwc_cp = iwc_cp.rename({"height": "height_cp"})
    iwc_b500 = iwc_b500.rename({"height": "height_b500"})
    iwc_b1000 = iwc_b1000.rename({"height": "height_b1000"})
    iwc_a1000 = iwc_a1000.rename({"height": "height_a1000"})
    iwc_a500 = iwc_a500.rename({"height": "height_a500"})

    iwc_da_list = [iwc_a1000, iwc_a500, iwc_cp, iwc_b500, iwc_b1000]

    return dict(zip(keys, iwc_da_list))


def save_cp_relative(months, year, region, cp_rel_dict, out_path, return_ds=False):
    """ 
    Save a dataset of the cold point-relative IWCs (at cold point
    and +/- 1000m and 500m). Returns the dataset if requested.
    """
    out_file_name_cp = out_path + format_string("DARDAR_cp_relative_iwc_{m}{y}_{r}.nc", months, year, region)
    iwc_cp_ds = xr.Dataset(cp_rel_dict)
    iwc_cp_ds.to_netcdf(out_file_name_cp)
    
    if return_ds:
        return iwc_cp_ds


def compute_cp_rel_dardar(months, year, region, iwc_rg, cpz_dar, dz=DZ_DARDAR, out_path=OUT_PATH):
    """ 
    Runs get_cp_relative_dict and save_cp_relative, which
    does the following:
        Takes the regridded IWC file, gets data arrays at the 
        cold point-relative levels (cp and +/- 1000m, 500m) and 
        saves as a netcdf. 
    """
    cp_rel_dict = get_cp_relative_dict(iwc_rg, cpz_dar, dz)
    print("dictionary done")
    try:
        save_cp_relative(months, year, region, cp_rel_dict, out_path, return_ds=False)
        print("Saved!")
    except:
        return cp_rel_dict


def main(dz=DZ_DARDAR, iwc_min=IWC_MIN, file_path=FILE_PATH, out_path=OUT_PATH, dask_progress=DASK_PROGRESS):
    """ 
    To run in the command line! Saves two netcdf files to out_path.
    
    Regrids the file, saves it, subsets
    the regridded file at cold point-relative levels, and saves that too.
    """
    months = str(sys.argv[1])
    year = str(sys.argv[2])
    region = str(sys.argv[3])
    
    if dask_progress:
        pbar = ProgressBar()
        pbar.register()
        
    # change the file paths if those are input
    if len(sys.argv) == 6:
        file_path = str(sys.argv[4])
        out_path = str(sys.argv[5])
        
    # read in the cpz_dar file, or make it if it doesn't exist yet,
    # and interpolate the files & save as regridded files for each group
    out_name_cpz_dar = format_string("ERA5_cpz_DARDAR_times_{m}{y}_{r}.nc", months, year, region)
    if os.path.isfile(out_path + out_name_cpz_dar):
        cpz_dar = xr.open_dataset(out_path + out_name_cpz_dar)
        interp_dardar(months, year, region, dz, iwc_min, file_path, out_path,
                             cpz_dar=cpz_dar)
    else:
        cpz_dar = interp_dardar(months, year, region, dz, iwc_min, file_path, out_path,
                                cpz_dar=None) 

    # read in the files and concatenate into one regridded file
    iwc_rg = concat_save_rg_files(months, year, region, out_path, return_iwc_rg=True)
    
    # # read in the files you've already made - for debugging, if needed
    # cpz_dar = xr.open_dataset(out_path + "ERA5_cpz_DARDAR_times_{m}{y}_{r}.nc".format(m=months, y=year, r=region))
    # iwc_rg = xr.open_dataset(out_path + "DARDAR_regridded_iwc_{m}{y}_{r}.nc".format(m=months, y=year, r=region))
    
    # ---- get cold point-relative regridded IWCs and save ----
    print("Getting the cp-relative dict!")
    compute_cp_rel_dardar(months, year, region, iwc_rg["iwc"], cpz_dar["cpz_dar"], dz, out_path)
    
    
# to run in the command line
if __name__ == "__main__":
    main()

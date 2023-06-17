"""
biv_hist.py

Script to calculate bivariate histogram (cold point temperature
and brightness temperature). Saves the histogram as a dictionary.

------------------
usage: biv_hist.py [-h] -y YEAR -f FILE_PATH -r REGION -m MONTHS
                   [-p SAVE_PICKLE_PATH] [-c CHUNK_SIZE]

optional arguments:
  -h, --help            show this help message and exit
  -y YEAR, --year YEAR  Year
  -f FILE_PATH, --file_path FILE_PATH
                        Input path for files
  -r REGION, --region REGION
                        Region abbreviation
  -m MONTHS, --months MONTHS
                        Months abbreviation (e.g., DJF)
  -p SAVE_PICKLE_PATH, --save_pickle_path SAVE_PICKLE_PATH
                        Output path to save pickle files (default =
                        file_path)
  -c CHUNK_SIZE, --chunk_size CHUNK_SIZE
                        Chunk size (number elements) to compute the
                        histogram
------------------       

Assumes the following file name convention:
* cold point: ERA5_cpT_reindexed_{m}{y}_{r}.nc
* brightness temp: MERGIR_Tb_4km_{m}{y}_{r}.nc4

Saves the dictionary as
    Tb-cpT_hist_dict_{m}{y}_{r}.pickle

Default bins are (175, 176, ..., 205) K for cold point
and (175, 180, 185, ..., 330) K for brightness temperature.
"""
import argparse
import dask
import gc
import pickle
import xarray as xr
import numpy as np


# bins for the bivariate histogram
CPT_BINS = np.arange(175, 206, 1)
TB_BINS = np.arange(175, 331, 5)


def parse_args():
    """ Parse command-line arguments
    """
    parser = argparse.ArgumentParser()
    
    # required
    parser.add_argument("-y", "--year", help="Year", type=int, required=True)
    parser.add_argument("-f", "--file_path", help="Input path for files", required=True)
    parser.add_argument("-r", "--region", help="Region abbreviation", required=True)
    parser.add_argument("-m", "--months", help="Months abbreviation (e.g., DJF)", required=True)

    # optional
    parser.add_argument("-p", "--save_pickle_path", help="Output path to save pickle files (default = file_path)")
    parser.add_argument("-c", "--chunk_size", help="Chunk size (number elements) to compute the histogram", default="5e6")

    return parser.parse_args()


def get_dask_arrays(cpT, tb, chunk_size):
    """ Convert DataArrays into Dask arrays. Takes ~4 min
    """
    cpT_da = dask.array.from_array(cpT.values.flatten(), chunks=chunk_size)
    tb_da = dask.array.from_array(tb.values.flatten(), chunks=chunk_size)

    return cpT_da, tb_da


def compute_histogram(cpT_da, tb_da, tb_count, tb_bins=TB_BINS, cpT_bins=CPT_BINS):
    """ 
    Compute the bivariate (cold point temp & brightness temp) histogram.
    
    Returns the dict with the histogram counts, x & y bins and edges, and 
    length of non-nan brightness temperature elements (need to get from histogram
    counts to frequencies).
    """
    binned_stat, xedges, yedges = dask.array.histogram2d(cpT_da, tb_da, bins=(cpT_bins, tb_bins))
    hist_computed = binned_stat.compute()
    hist_dict = {
        "hist_computed": hist_computed, 
        "xedges": xedges, 
        "yedges": yedges, 
        "tb_bins": tb_bins, 
        "cpT_bins": cpT_bins,
        "nan_len": int(tb_count),
    }
    return hist_dict


def main():
    """ docstring
    """
    args = parse_args()
    print(args)

    year = args.year
    file_path = args.file_path
    region = args.region
    months = args.months
    save_pickle_path = args.save_pickle_path
    chunk_size = float(args.chunk_size)
    
    if file_path[-1] != "/":
        file_path = file_path + "/"
    if save_pickle_path is None:
        save_pickle_path = file_path
    else:
        if save_pickle_path[-1] != "/":
            save_pickle_path = save_pickle_path + "/"
    
    if region == "SPC":
        cpT_file_1 = "ERA5_cpT_reindexed_{m}{y}_{r}1.nc".format(m=months, y=year, r=region)
        cpT_file_2 = "ERA5_cpT_reindexed_{m}{y}_{r}2.nc".format(m=months, y=year, r=region)
        cpT_1 = xr.open_dataset(file_path + cpT_file_1)["t"]
        cpT_2 = xr.open_dataset(file_path + cpT_file_2)["t"]
        cpT = xr.concat([cpT_1, cpT_2], dim="lon")
        print("SPC1 & SPC2 cpT joined along longitude")
        tb_file_1 = "MERGIR_Tb_4km_{m}{y}_{r}1.nc4".format(m=months, y=year, r=region)
        tb_file_2 = "MERGIR_Tb_4km_{m}{y}_{r}2.nc4".format(m=months, y=year, r=region)
        tb1 = xr.open_dataset(file_path + tb_file_1)["Tb"]
        tb2 = xr.open_dataset(file_path + tb_file_2)["Tb"]
        tb = xr.concat([tb1, tb2], dim="lon")
        print("SPC1 & SPC2 Tb joined along longitude")
        tb_count = tb.count().values
        del tb1
        del tb2
        del cpT_1
        del cpT_2
        gc.collect()
    else:
        cpT_file = "ERA5_cpT_reindexed_{m}{y}_{r}.nc".format(m=months, y=year, r=region)
        tb_file = "MERGIR_Tb_4km_{m}{y}_{r}.nc4".format(m=months, y=year, r=region)
        tb = xr.open_dataset(file_path + tb_file)["Tb"]
        tb_count = tb.count().values
        cpT = xr.open_dataset(file_path + cpT_file)["t"]
    
    # reindex the cold point temperature to half hourly (match Tb)
    cpT_30m = cpT.reindex({"time": tb.time}, method="nearest")
    if region == "SPC":
        del cpT
        gc.collect()
    
    # compute & save histogram dictionary
    dict_file_name = "Tb-cpT_hist_dict_{m}{y}_{r}.pickle".format(m=months, y=year, r=region)
    cpT_da, tb_da = get_dask_arrays(cpT_30m, tb, chunk_size)
    if region == "SPC":
        del cpT_30m
        gc.collect()
    hist_dict = compute_histogram(cpT_da, tb_da, tb_count)
    with open(save_pickle_path + dict_file_name, "wb") as handle:
        pickle.dump(hist_dict, handle)
    print("Dictionary saved to " + save_pickle_path + dict_file_name)
    

if __name__ == "__main__":
    main()      

"""
bin_obs_overshoot.py

Module to make plots of IWC binned by brightness
temperature for DARDAR observations. Companion
to bin_overshoot.py.
"""
import argparse
import pickle
import dask
import sys
sys.path.append("/home/b/b380887/cold-point-overshoot/python_scripts/")

import xarray as xr
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import bin_overshoot as bin_os
import regrid_data_cp as rdc

from dask.diagnostics import ProgressBar


# Based on SHIELD median cold point/pressure levels
REF_DENSITIES = [
    0.180602, # -1000 m
    0.164667, # -500 m
    0.151853, # at cold point
    0.136337, # +500 m
    0.122172  # +1000 m
]

# Defaults bins for Tb and (Tb-cpT)
TB_BINS = np.arange(186, 301, 2) # new
DIFF_BINS = np.arange(-20, 101, 2)

QI_MIN = 1e-6 # kg/kg
IWC_MIN = 1e-7 # kg/m3
FILE_PATH = "/work/bb1153/b380887/big_obs_climo/"
PICKLE_DIR = "/home/b/b380887/cold-point-overshoot/pickle_files/binned_by_tb/"


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
    parser.add_argument("-n", "--no_save", help="Do not save dictionaries (as pickle files) or plots", action='store_true')
    parser.add_argument("-p", "--save_pickle_path", help="Output path to save pickle files (default = file_path)")
    parser.add_argument("-s", "--save_fig_path", help="Output path to save figures (default = file_path)")
    parser.add_argument("-d", "--bin_by_diffs", help="Bin by (Tb - cpT) instead of Tb", action='store_true')

    return parser.parse_args()


def get_data(month, year, region, file_path, sortby_time=False, bin_by_diffs=False):
    """ 
    Returns datasets of seasonal or monthly IWC and Tb for that
    year/region. If month = "JJA" or "DJF", tries to open those files first,
    then reads in/concatenates the individual monthly files if it fails.
    
    If sortby_time=True, sorts the IWC dataset by the time dimension 
    (need for some Jul/Aug 2009 regions).
    If bin_by_diffs=True, also returns the reindexed cold point temp (to match Tb grid)
    """
    # --- IWC ---
    ds_all = xr.open_dataset(file_path + "DARDAR_cp_relative_iwc_{m}{y}_{r}.nc".format(m=month, y=year, r=region))
    if sortby_time:
        ds_all = ds_all.sortby("time")
        
    # --- Tb ---
    tb = xr.open_dataset(
        file_path + "MERGIR_Tb_4km_{m}{y}_{r}.nc4".format(m=month, y=year, r=region)
    )["Tb"]
    
    # --- regridded cpT ---
    if bin_by_diffs:
        cpT_ri = xr.open_dataset(
            file_path + "ERA5_cpT_reindexed_{m}{y}_{r}.nc".format(m=month, y=year, r=region)
        )["t"]

    if bin_by_diffs:
        return [ds_all, tb, cpT_ri]
    else:
        return [ds_all, tb]


def count_dardar_retrievals(month, year, region, file_path, dar_list_ds=None):
    """ Return the total number of time steps/columns, i.e. retrievals, for DARDAR
    """
    if dar_list_ds is None:
        da = xr.open_dataset(file_path + "DARDAR-v3_iwc_{m}{y}_{r}.nc".format(m=month, y=year, r=region))["iwc"]
        dar_count = len(da.time)
    else:
        dar_count = 0
        for ds in dar_list_ds:
            dar_count += len(ds["iwc"].time)

    return dar_count


def calc_avg_cold_point(month, year, region, file_path, return_cpT=False):
    """ 
    Returns the time- and area-mean cold point temperature for that
    region/year/month. If return_cpT=True, also returns the timexlatxlon cold
    point temperature.
    """
    temp = xr.open_dataset(file_path + "ERA5_T_0.25deg_ml_12-20km_{m}{y}_{r}.nc".format(m=month, y=year, r=region))["t"]
    cpT_inds = temp.argmin(dim="level")
    cpT = temp.isel(level=cpT_inds)
    cpT_avg = cpT.mean().values
    
    if return_cpT:
        return [cpT, cpT_avg]
    else:
        return cpT_avg
    
        
def get_qi_levs_dict(ds_all, convert_iwc=True, ref_densities=REF_DENSITIES, iflag=None):
    """
    Converts iwc --> qi (if convert_iwc=True) and makes a dictionary with keys of the offset
    index as a string ("-2" for 1 km below, "0" for at cold point, etc.)
    and values of the qi variable at that level. Takes only the values where
    instrument_flag = iflag if iflag value is passed in (default=None). 
    If convert_iwc=False, keeps variable as IWC (i.e., "reference densities" are set to 1)
    """
    # check that it's bottom--> up and flip if not
    if convert_iwc:
        if ref_densities[4] < ref_densities[0]:
            ref_densities = ref_densities[::-1]
    else:
        ref_densities = [1]*len(ref_densities)
    
    if iflag is not None:
        qi_a1000 = ds_all["iwc_a1000"].where(ds_all["iflag_a1000"] == iflag)/ref_densities[4]
        qi_a500 = ds_all["iwc_a500"].where(ds_all["iflag_a500"] == iflag)/ref_densities[3]
        qi_cp = ds_all["iwc_cp"].where(ds_all["iflag_cp"] == iflag)/ref_densities[2]
        qi_b500 = ds_all["iwc_b500"].where(ds_all["iflag_b500"] == iflag)/ref_densities[1]
        qi_b1000 = ds_all["iwc_b1000"].where(ds_all["iflag_b1000"] == iflag)/ref_densities[0]
    else:
        qi_a1000 = ds_all["iwc_a1000"]/ref_densities[4]
        qi_a500 = ds_all["iwc_a500"]/ref_densities[3]
        qi_cp = ds_all["iwc_cp"]/ref_densities[2]
        qi_b500 = ds_all["iwc_b500"]/ref_densities[1]
        qi_b1000 = ds_all["iwc_b1000"]/ref_densities[0]
        
    qi_levs_dict = {
        "-2": qi_b1000,
        "-1": qi_b500,
        "0": qi_cp,
        "1": qi_a500,
        "2": qi_a1000
    }
    
    return qi_levs_dict


def get_Re_levs_dict(ds_all, iflag=3):
    """
    Makes a dictionary with keys of the offset
    index as a string ("-2" for 1 km below, "0" for at cold point, etc.)
    and values of the effective radius variable at that level (when the instrument
    flag has value iflag, default=3/both).
    """
    re_levs_dict = {
        "-2": ds_all["Re_b1000"].where(ds_all["iflag_b1000"] == iflag),
        "-1": ds_all["Re_b500"].where(ds_all["iflag_b500"] == iflag),
        "0": ds_all["Re_cp"].where(ds_all["iflag_cp"] == iflag),
        "1": ds_all["Re_a500"].where(ds_all["iflag_a500"] == iflag),
        "2": ds_all["Re_a1000"].where(ds_all["iflag_a1000"] == iflag)
    }
    
    return re_levs_dict


def get_tb_dar(tb, ds_all):
    """ Returns brightness temp at DARDAR times only
    """
    return tb.sel(time=ds_all.time)
    

def bin_var_by_tb(ds_all, tb_dar, varname, bins=TB_BINS, save_dicts=False, pickle_dir=PICKLE_DIR,
                 month=None, year=None, region=None, qi_min=QI_MIN, ref_densities=REF_DENSITIES, iflag=None):
    """ 
    Bins the qi or effective radius at each cold point-relative level by the brightness temperature.
    Defaults to bins=TB_BINS if bins not provided and then returns the bins.
    Returns list of dictionaries and the bins: [bin_means_dict, bin_counts_dict, bins].
    """
    if varname == "qi":
        var_levs_dict = get_qi_levs_dict(ds_all, ref_densities, iflag)
    elif varname == "Re":
        var_levs_dict = get_Re_levs_dict(ds_all, iflag)
       
    var_bin_means = {}
    var_bin_counts = {}

    for key in var_levs_dict.keys():
        if varname == "qi":
            var_lev = var_levs_dict[key].where(var_levs_dict[key] > qi_min)
        else:
            var_lev = var_levs_dict[key]

        # if there are no values at that level
        if var_lev.count().values == 0:
            var_bin_counts[key] = [np.nan]*(len(bins) - 1)
            var_bin_means[key] = [np.nan]*(len(bins) - 1)
        else:
            counts = bin_os.bin_single_level(var_lev, tb_dar, bins=bins, statistic="count")
            var_bin_counts[key] = np.where(counts > 0, counts, np.nan)
            means = bin_os.bin_single_level(var_lev, tb_dar, bins=bins, statistic="mean")
            var_bin_means[key] = np.where(counts > 0, means, np.nan)
            
    if save_dicts:
        if pickle_dir is None or month is None or year is None or region is None:
            raise Exception("Must input pickle_dir, month (or season), year, AND region if you want to save the dictionaries")
        if iflag is not None:
            mean_name = pickle_dir + "{v}_iflag{f}_30min_bin_means_{m}{y}_{r}.pickle".format(v=varname, f=iflag, m=month, y=year, r=region)
            count_name = pickle_dir + "{v}_iflag{f}_30min_bin_counts_{m}{y}_{r}.pickle".format(v=varname, f=iflag, m=month, y=year, r=region)
        else:
            mean_name = pickle_dir + "{v}_30min_bin_means_{m}{y}_{r}.pickle".format(v=varname, m=month, y=year, r=region)
            count_name = pickle_dir + "{v}_30min_bin_counts_{m}{y}_{r}.pickle".format(v=varname, m=month, y=year, r=region)
        with open(mean_name, 'wb') as handle:
            pickle.dump(var_bin_means, handle)
        with open(count_name, 'wb') as handle:
            pickle.dump(var_bin_counts, handle)

    return [var_bin_means, var_bin_counts, bins]


def bin_var_by_diffs(ds_all, diffs, varname, bins=DIFF_BINS, save_dicts=False, pickle_dir=PICKLE_DIR,
                      month=None, year=None, region=None, qi_min=QI_MIN, ref_densities=REF_DENSITIES, iwc_min=IWC_MIN, iflag=None):
    """ 
    Bins the qi, iwc, or effective radius at each cold point-relative level by 
    (brightness temperature - cold point temperature). 
    Input diffs is the (brightness temp - regridded cold point temperature).
    Defaults to bins=DIFF_BINS if bins not provided and then returns the bins.
    Returns list of dictionaries and the bins: [bin_means_dict, bin_counts_dict, bins].
    """
    if varname == "qi":
        var_levs_dict = get_qi_levs_dict(ds_all, convert_iwc=True, ref_densities=ref_densities,
                                         iflag=iflag)
    elif varname == "iwc":
        var_levs_dict = get_qi_levs_dict(ds_all, convert_iwc=False, iflag=iflag)
    elif varname == "Re":
        var_levs_dict = get_Re_levs_dict(ds_all, iflag)
           
    var_bin_means = {}
    var_bin_counts = {}

    for key in var_levs_dict.keys():
        if varname == "qi":
            var_lev = var_levs_dict[key].where(var_levs_dict[key] > qi_min)
        elif varname == "iwc":
            var_lev = var_levs_dict[key].where(var_levs_dict[key] > iwc_min)
        else:
            var_lev = var_levs_dict[key]

        # if there are no values at that level
        if var_lev.count().values == 0:
            var_bin_counts[key] = [np.nan]*(len(bins) - 1)
            var_bin_means[key] = [np.nan]*(len(bins) - 1)
        else:
            counts = bin_os.bin_single_level(var_lev, diffs, bins=bins, statistic="count")
            var_bin_counts[key] = np.where(counts > 0, counts, np.nan)
            means = bin_os.bin_single_level(var_lev, diffs, bins=bins, statistic="mean")
            var_bin_means[key] = np.where(counts > 0, means, np.nan)
            
    if save_dicts:
        if pickle_dir is None or month is None or year is None or region is None:
            raise Exception("Must input pickle_dir, month (or season), year, AND region if you want to save the dictionaries")
        if iflag is not None:
            mean_name = pickle_dir + "{v}_iflag{f}_30min_bin_means_{m}{y}_{r}.pickle".format(v=varname, f=iflag, m=month, y=year, r=region)
            count_name = pickle_dir + "{v}_iflag{f}_30min_bin_counts_{m}{y}_{r}.pickle".format(v=varname, f=iflag, m=month, y=year, r=region)
        else:
            mean_name = pickle_dir + "{v}_30min_bin_means_{m}{y}_{r}.pickle".format(v=varname, m=month, y=year, r=region)
            count_name = pickle_dir + "{v}_30min_bin_counts_{m}{y}_{r}.pickle".format(v=varname, m=month, y=year, r=region)
        with open(mean_name, 'wb') as handle:
            pickle.dump(var_bin_means, handle)
        with open(count_name, 'wb') as handle:
            pickle.dump(var_bin_counts, handle)

    return [var_bin_means, var_bin_counts, bins]


def load_saved_dicts(varname, month, region, year, pickle_dir=PICKLE_DIR, iflag=None):
    """ Returns saved bin mean and bin count dictionaries
    """
    if iflag is not None:
        mean_name = pickle_dir + "{v}_iflag{f}_30min_bin_means_{m}{y}_{r}.pickle".format(v=varname, f=iflag, m=month, y=year, r=region)
        count_name = pickle_dir + "{v}_iflag{f}_30min_bin_counts_{m}{y}_{r}.pickle".format(v=varname, f=iflag, m=month, y=year, r=region)
    else:
        mean_name = pickle_dir + "{v}_30min_bin_means_{m}{y}_{r}.pickle".format(v=varname, m=month, y=year, r=region)
        count_name = pickle_dir + "{v}_30min_bin_counts_{m}{y}_{r}.pickle".format(v=varname, m=month, y=year, r=region)
            
    with open(mean_name, 'rb') as handle:
        var_bin_means = pickle.load(handle)
    with open(count_name, 'rb') as handle:
        var_bin_counts = pickle.load(handle)
        
    return [var_bin_means, var_bin_counts]
    

def plot_tb_hist(tb, month, region, year, bins=TB_BINS, dar_times_only=False, fsize=13, 
                    tsize=14, figsize=(7, 4), ylim=(0, 0.035),
                    gridlines=True, annotate_total=True, save=False, save_dir=None, dar_file_path=FILE_PATH, dar_list_ds=None, return_hist=False):
    """ 
    Plot a brightness temperature histogram for MERGIR. Annotate with # of total DARDAR
    retrievals if annotate_total=True.  If dar_times_only is True, will plot only Tb values
    conditioned on times when you have DARDAR retrievals; otherwise, will plot all Tb values for 
    that month (or season)/year/region.
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # don't need dask if it's conditional on dardar retrieval times
    if dar_times_only:
        n_pts = len(tb.values.flatten())
        hist_comp, bin_edges = np.histogram(tb.values.flatten(), bins=bins, range=[bins[0], bins[-1]])
    else:
        n_pts = tb.size
        tb_chunked = tb.chunk({"time": 8})
        tbhist, bin_edges = dask.array.histogram(tb_chunked, bins=bins, range=[bins[0], bins[-1]])
        with ProgressBar():
            hist_comp = tbhist.compute()
    bin_means = 0.5*(bin_edges[:-1] + bin_edges[1:])
    ax.bar(bin_means, hist_comp/n_pts, width=2, edgecolor="w")
        
    ax.set_ylabel("Frequency", fontsize=fsize)
    ax.set_xlabel("Brightness temperature (K)", fontsize=fsize)
    ax.set_xlim((bins[0], bins[-1]))
    ax.set_ylim(ylim)
    ax.tick_params(axis="both", labelsize=fsize-1)
    
    if annotate_total:
        ndar = count_dardar_retrievals(month, year, region, file_path=dar_file_path, dar_list_ds=dar_list_ds)
        annotation = "Total DARDAR retrievals: {:.1e}".format(ndar)
        ax.annotate(annotation, xy=(0.05, 0.9), xycoords="axes fraction", fontsize=fsize-1)
    if gridlines:
        ax.grid(color="grey", linestyle=":")
    
    # if values are conditioned on dardar or not
    if dar_times_only:
        ax.set_title(
            "MERGIR {m} {y} ({r}): frequency of Tb bins\ntimes of DARDAR retrievals only".format(m=month, y=year, r=region),
            fontsize=tsize
        )
    else:
        ax.set_title(
            "MERGIR {m} {y} ({r}): frequency of Tb bins\nall times".format(m=month, y=year, r=region),
            fontsize=tsize
        )
        
    if save:
        if save_dir is None:
            raise Exception("Must provide save_dir if you want to save the plot")
        if dar_times_only:
            plt.savefig(save_dir + "MERGIR_Tb_bin_frequency_dar_times_only_{m}{y}_{r}.png".format(m=month, y=year, r=region))
        else:
            plt.savefig(save_dir + "MERGIR_Tb_bin_frequency_{m}{y}_{r}.png".format(m=month, y=year, r=region))
    
    plt.show()
    
    if return_hist:
        return hist_comp, bins, bin_edges, n_pts
    
    
def plot_diffs_hist(diffs, month, region, year, bins=DIFF_BINS, dar_times_only=False, fsize=13, 
                    tsize=14, figsize=(7, 4), ylim=(0, 0.035),
                    gridlines=True, annotate_total=True, save=False, save_dir=None, dar_file_path=FILE_PATH, dar_list_ds=None, return_hist=False):
    """ 
    Plot a histogram of (brightness temperature - cold point temperature). Annotate with # of total DARDAR
    retrievals if annotate_total=True.  If dar_times_only is True, will plot only diff values
    conditioned on times when you have DARDAR retrievals; otherwise, will plot all diff values for 
    that month (or season)/year/region.
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # don't need dask if it's conditional on dardar retrieval times
    if dar_times_only:
        n_pts = len(diffs.values.flatten())
        hist_comp, bin_edges = np.histogram(diffs.values.flatten(), bins=bins, range=[bins[0], bins[-1]])
    else:
        n_pts = diffs.size
        diffs_chunked = diff.chunk({"time": 8})
        diffhist, bin_edges = dask.array.histogram(diffs_chunked, bins=bins, range=[bins[0], bins[-1]])
        with ProgressBar():
            hist_comp = diffhist.compute()
    bin_means = 0.5*(bin_edges[:-1] + bin_edges[1:])
    ax.bar(bin_means, hist_comp/n_pts, width=2, edgecolor="w")
        
    ax.set_ylabel("Frequency", fontsize=fsize)
    ax.set_xlabel("$T_b - CP_T$ (K)", fontsize=fsize)
    ax.set_xlim((bins[0], bins[-1]))
    ax.set_ylim(ylim)
    ax.tick_params(axis="both", labelsize=fsize-1)
    
    if annotate_total:
        ndar = count_dardar_retrievals(month, year, region, file_path=dar_file_path, dar_list_ds=dar_list_ds)
        annotation = "Total DARDAR retrievals: {:.1e}".format(ndar)
        ax.annotate(annotation, xy=(0.05, 0.9), xycoords="axes fraction", fontsize=fsize-1)
    if gridlines:
        ax.grid(color="grey", linestyle=":")
    
    # if values are conditioned on dardar or not
    if dar_times_only:
        ax.set_title(
            "{m} {y} ({r}): frequency of ($T_b - CP_T$) bins\ntimes of DARDAR retrievals only".format(m=month, y=year, r=region),
            fontsize=tsize
        )
    else:
        ax.set_title(
            "{m} {y} ({r}): frequency of ($T_b - CP_T$) bins\nall times".format(m=month, y=year, r=region),
            fontsize=tsize
        )
        
    if save:
        if save_dir is None:
            raise Exception("Must provide save_dir if you want to save the plot")
        if dar_times_only:
            plt.savefig(save_dir + "Tb-cpT_bin_frequency_dar_times_only_{m}{y}_{r}.png".format(m=month, y=year, r=region))
        else:
            plt.savefig(save_dir + "Tb-cpT_bin_frequency_{m}{y}_{r}.png".format(m=month, y=year, r=region))
    
    plt.show()
    
    if return_hist:
        return hist_comp, bins, bin_edges, n_pts
    
    
def plot_binned_by_tb(bin_dict, tb_dar, statistic, month, region, year, bins=TB_BINS,
                      stat_lims=None, fsize=26, tsize=30, figsize=[13, 6],
                       colormap="Spectral", cbar_ext="max", lognorm=True, save=False,
                      save_dir=None, varname="qi", save_extr=""): 
    """
    Plot DARDAR qi values binned by MERGIR brightness temperatures at levels relative to the cold 
    point. Variable "statistic" is "count" or "mean"
    """
    figsize = [16, 6] # TEMP
    if varname == "qi":
        long_name = "cloud ice"
        units = "kg/kg"
    elif varname == "reff":
        long_name = "eff. radius"
        units = "$\mu m$"
        
    if stat_lims is None:
        if statistic == "mean":
            if varname == "reff":
                stat_lims = (0, 100)
            elif varname == "qi":
                stat_lims = (1e-6, 5e-4)
        elif statistic == "count":
            stat_lims = (1, 1000)

    
    ind_offsets = [-2, -1, 0, 1, 2]  
    offset_labs = ["~1000 m below", 
           "~500 m below", 
           "At cold point", 
           "~500 m above", 
           "~1000 m above"
          ]
    
    # PLOT!
    fig, ax1 = plt.subplots(figsize=figsize)
    plt.subplots_adjust(hspace=0.5)
    bin_mean_values = (bins[:-1] + bins[1:])/2
            
    # turn the dict values into arrays
    bin_dict_arr = np.array(list(bin_dict.values()))
    if varname == "reff":
        bin_dict_arr = bin_dict_arr*1e6 # convert from m --> um
        

    if lognorm:
        pcm = ax1.pcolormesh(bin_mean_values, ind_offsets, bin_dict_arr, cmap=colormap,
                             norm=mcolors.LogNorm(vmin=stat_lims[0], vmax=stat_lims[1])
                            )
    else:
        pcm = ax1.pcolormesh(bin_mean_values, ind_offsets, bin_dict_arr, cmap=colormap,
                             vmin=stat_lims[0], vmax=stat_lims[1]
                            )

    ax1.set_xlabel("Brightness temperature (K)", fontsize=fsize)
    ax1.set_yticks(ind_offsets)
    ax1.set_yticklabels(offset_labs)
    ax1.tick_params(axis="y", labelsize=fsize, length=0)
    ax1.set_xticks(bins[::4])        
    ax1.tick_params(axis="x", rotation=45, labelsize=fsize-1)
    ax1.xaxis.set_minor_locator(mticker.MultipleLocator(2))
    ax1.tick_params(which='minor', length=10)
    ax1.tick_params(which='major', length=14)
    
    plt.suptitle("{m} {y} ({r})".format(m=month, y=year, r=region), fontsize=tsize)

    cb = plt.colorbar(pcm, ax=ax1, fraction=0.046, pad=0.04, extend=cbar_ext)
    cb.ax.tick_params(labelsize=fsize)
    if statistic == "count":
        cb.set_label("Bin counts", fontsize=fsize)
    elif statistic == "mean":
        cb.set_label("Bin-mean {n} ({u})".format(n=long_name, u=units), fontsize=fsize)

    # add horizontal borders
    for ind in ind_offsets[:-1]:
        gap = 0.5
        ax1.axhline(ind+0.5, color="k", linewidth=2)


    if save:
        if save_dir is None:
            raise Exception("Must provide save_dir if you want to save the plot")           
        plt.savefig(save_dir + "DARDAR_{v}_binned_by_Tb_{s}_{m}{y}_{r}{e}.png".format(v=varname, 
                                                                                      s=statistic,
                                                                                  m=month,
                                                                                  y=year, 
                                                                                  r=region,
                                                                                      e=save_extr,
                                                                                 ),
                    dpi=300,
                    bbox_inches="tight"
                   )

    plt.show()
    
    
def plot_binned_by_diffs(bin_dict, diffs, statistic, month, region, year, bins=DIFF_BINS,
                      stat_lims=None, fsize=26, tsize=30, figsize=[13, 6],
                       colormap="Spectral", cbar_ext="max", lognorm=True, save=False,
                      save_dir=None, varname="qi", save_extr=""): 
    """
    Plot DARDAR qi values binned by (Tb - cold point temp) at levels relative to the cold 
    point. Variable "statistic" is "count" or "mean"
    """
    figsize = [16, 6] # TEMP
    if varname == "qi":
        long_name = "cloud ice"
        units = "kg/kg"
    elif varname == "reff":
        long_name = "eff. radius"
        units = "$\mu m$"
        
    if stat_lims is None:
        if statistic == "mean":
            if varname == "reff":
                stat_lims = (0, 100)
            elif varname == "qi":
                stat_lims = (1e-6, 5e-4)
        elif statistic == "count":
            stat_lims = (1, 1000)

    
    ind_offsets = [-2, -1, 0, 1, 2]  
    offset_labs = ["~1000 m below", 
           "~500 m below", 
           "At cold point", 
           "~500 m above", 
           "~1000 m above"
          ]
    
    # PLOT!
    fig, ax1 = plt.subplots(figsize=figsize)
    plt.subplots_adjust(hspace=0.5)
    bin_mean_values = (bins[:-1] + bins[1:])/2
            
    # turn the dict values into arrays
    bin_dict_arr = np.array(list(bin_dict.values()))
    if varname == "reff":
        bin_dict_arr = bin_dict_arr*1e6 # convert from m --> um
        

    if lognorm:
        pcm = ax1.pcolormesh(bin_mean_values, ind_offsets, bin_dict_arr, cmap=colormap,
                             norm=mcolors.LogNorm(vmin=stat_lims[0], vmax=stat_lims[1])
                            )
    else:
        pcm = ax1.pcolormesh(bin_mean_values, ind_offsets, bin_dict_arr, cmap=colormap,
                             vmin=stat_lims[0], vmax=stat_lims[1]
                            )

    ax1.set_xlabel("$T_b - CP_T$ (K)", fontsize=fsize)
    ax1.set_yticks(ind_offsets)
    ax1.set_yticklabels(offset_labs)
    ax1.tick_params(axis="y", labelsize=fsize, length=0)
    ax1.set_xticks(bins[::4])        
    ax1.tick_params(axis="x", rotation=45, labelsize=fsize-1)
    ax1.xaxis.set_minor_locator(mticker.MultipleLocator(2))
    ax1.tick_params(which='minor', length=10)
    ax1.tick_params(which='major', length=14)
    
    plt.suptitle("{m} {y} ({r})".format(m=month, y=year, r=region), fontsize=tsize)

    cb = plt.colorbar(pcm, ax=ax1, fraction=0.046, pad=0.04, extend=cbar_ext)
    cb.ax.tick_params(labelsize=fsize)
    if statistic == "count":
        cb.set_label("Bin counts", fontsize=fsize)
    elif statistic == "mean":
        cb.set_label("Bin-mean {n} ({u})".format(n=long_name, u=units), fontsize=fsize)

    # add horizontal borders
    for ind in ind_offsets[:-1]:
        gap = 0.5
        ax1.axhline(ind+0.5, color="k", linewidth=2)
    
    # add zero line
    ax1.axvline(0, color="k", linestyle="--", linewidth=3)


    if save:
        if save_dir is None:
            raise Exception("Must provide save_dir if you want to save the plot")           
        plt.savefig(save_dir + "DARDAR_{v}_binned_by_Tb-cpT_{s}_{m}{y}_{r}{e}.png".format(v=varname, 
                                                                                      s=statistic,
                                                                                  m=month,
                                                                                  y=year, 
                                                                                  r=region,
                                                                                      e=save_extr,
                                                                                 ),
                    dpi=300,
                    bbox_inches="tight"
                   )

    plt.show()
    
    
def get_tb_and_cp_ri(month, year, region, dar_times_only=False, 
                     ds_all=None, file_path=FILE_PATH):
    """ 
    Returns brightness temperature coarsened to hourly and cold point 
    temperature reindexed to match hourly tb grid.
    """
    if dar_times_only and ds_all is None:
        raise Exception("Must provide ds_all to get tb and cp at times of retrievals")
    
    if month == "JJA" or month == "DJF":
        tb_30m = xr.open_dataset(
            file_path + "{m}/MERGIR_Tb_4km_{m}{y}_{r}.nc4".format(m=month, y=year, r=region)
        )["Tb"]
    else:
        tb_30m = xr.open_dataset(
            file_path + "MERGIR_Tb_4km_{m}{y}_{r}.nc4".format(m=month, y=year, r=region)
        )["Tb"]
    if dar_times_only:
        tb_30m = tb_30m.sel(time=ds_all.time)
    tb = tb_30m.resample(time="1h").mean()

    cpT, _ = calc_avg_cold_point(month, year, region, return_cpT=True, file_path=file_path)
    cpT_ri = cpT.rename({"latitude": "lat", "longitude": "lon"}).reindex_like(tb, method="nearest")
    
    return [tb, cpT_ri]
    

def joint_tb_cp_hist(tb, cpT, month, year, region, plot_type, offset=None, file_path=FILE_PATH, 
                     dar_times_only=False, cpT_bins=np.arange(175, 206, 1), 
                     tb_bins=np.arange(176, 241, 2), levels=np.arange(-6, -1), figsize=(6, 4), 
                     fsize=14, tsize=16, save=False, save_dir=None, Tb_limit=240,
                    save_hist=True, pickle_dir=PICKLE_DIR):
    """ 
    Calculate and plot joint histogram of brightness temperature and cold point temperature.
    Can do pcolormesh (plot_type = "pcolormesh") or filled contours (plot_type = "contourf").
    Input tb and cpT as None to load in a previously-saved histogram.
    """
    # read in histogram if it's there
    if tb is None or cpT is None:
        out_file = pickle_dir + "OBS_tb-cp_joint_hist_{m}{y}_{r}.pickle".format(m=month, y=year, r=region)
        with open(out_file, 'rb') as handle:
            hist_dict = pickle.load(handle)
        xbin_means = hist_dict["xbin_means"]
        ybin_means = hist_dict["ybin_means"]
        hist_normed_nonzero = hist_dict["hist_normed_nonzero"]
    
    # calculate it yourself if it;s not
    else:  
        # take only Tb values below a certain threshold
        if Tb_limit is not None:
            tb = tb.where(tb < Tb_limit)
        # use dask to compute the histogram
        cpT_da = dask.array.from_array(cpT.values.ravel(), chunks="auto")
        tb_da = dask.array.from_array(tb.values.ravel(), chunks="auto")
        binned_stat, xedges, yedges = dask.array.histogram2d(cpT_da, tb_da, bins=(cpT_bins, tb_bins))
        with ProgressBar():
            hist_computed = binned_stat.compute()

        # normalize the histogram
        nan_len = tb.count().values
        hist_normed = hist_computed/(nan_len) # normalized bin count
        xbin_means, ybin_means = (xedges[:-1]+xedges[1:])/2, (yedges[:-1]+yedges[1:])/2
        hist_normed_nonzero = np.where(hist_normed > 0, hist_normed, np.nan)
    
        # save the histogram
        if save_hist:
            hist_dict = {
                "hist_normed": hist_normed,
                "hist_normed_nonzero": hist_normed_nonzero,
                "xbin_means": xbin_means,
                "ybin_means": ybin_means
            }
            out_file = pickle_dir + "OBS_tb-cp_joint_hist_{m}{y}_{r}.pickle".format(m=month, y=year, r=region)
            with open(out_file, 'wb') as handle:
                pickle.dump(hist_dict, handle)
        
        
    # plot
    fig, ax = plt.subplots(figsize=figsize)

    if plot_type == "pcolormesh":
        csn = ax.pcolormesh(xbin_means, ybin_means, np.log10(hist_normed_nonzero.T),  
                          vmin=levels[0], vmax=levels[-1])
        cb = plt.colorbar(csn, ax=ax, extend="both")
    elif plot_type == "contourf":
        ax.grid(linestyle=":", color="gray")
        if levels is None:
            csn = ax.contourf(xbin_means, ybin_means, np.log10(hist_normed_nonzero.T), extend='both')
        else:
            csn = ax.contourf(xbin_means, ybin_means, np.log10(hist_normed_nonzero.T), extend='both', levels=levels, 
                              vmin=levels[0], vmax=levels[-1])
        ax.contour(csn, colors='k', linestyles='solid', linewidths=1) 
        cb = plt.colorbar(csn, ax=ax)
    
    cb.set_label("log$_{10}$(PDF)", fontsize=fsize-1)
    ax.set_xlabel("Cold Point Temperature (K)", fontsize=fsize)
    ax.set_ylabel("Brightness Temperature (K)", fontsize=fsize)
    cb.ax.tick_params(axis="y", labelsize=fsize-2)
    ax.tick_params(axis="both", labelsize=fsize-2)
    ax.set_ylim((tb_bins[0], tb_bins[-1]))
    ax.set_xlim((cpT_bins[0], cpT_bins[-1]))
    
    if dar_times_only:
        ax.set_title("{m} {y} ({r})\ntimes of DARDAR retrievals only".format(m=month, y=year, r=region), fontsize=tsize)
    else:
        ax.set_title("{m} {y} ({r})\nall times".format(m=month, y=year, r=region), fontsize=tsize)
    
    # offset line (Tb = cp + offset) and 1-1 Tb=cp line
    ax.plot(cpT_bins, cpT_bins, color="k", linestyle="--")
    if offset is not None:
        ax.plot(cpT_bins, [x+offset for x in cpT_bins], color="k")
    
    # save
    if save:
        plt.savefig(save_dir + "OBS_Tb_cpT_joint_hists_{m}{y}_{r}.png".format(m=month, y=year, r=region), dpi=300, bbox_inches="tight")
    
    plt.show()
    
    
def main():
    """ Bin qi by Tb, save dictionaries, and make "default" versions of the plots
    """
    args = parse_args()
    print(args)

    year = args.year
    file_path = args.file_path
    region = args.region
    months = args.months
    no_save = args.no_save
    save_fig_path = args.save_fig_path
    save_pickle_path = args.save_pickle_path
    bin_by_diffs = args.bin_by_diffs
    
    # make sure output path names end in "/"
    if file_path[-1] != "/":
        file_path = file_path + "/"
    
    # default to file_path if not specified
    if save_fig_path is None:
        save_fig_path = file_path
    else:
        if save_fig_path[-1] != "/":
            save_fig_path = save_fig_path + "/"
    
    # default to file_path if not specified
    if save_pickle_path is None:
        save_pickle_path = file_path
    else:
        if save_pickle_path[-1] != "/":
            save_pickle_path = save_pickle_path + "/"
            
    if no_save:
        save_dicts = False
        save_plots = False
    else:
        save_dicts = True
        save_plots = True
    

    # get the regridded data
    if bin_by_diffs:
        ds_all, tb, cpT_ri = get_data(months, year, region, file_path=file_path, bin_by_diffs=True)
        tb_dar = tb.sel(time=ds_all["iwc_cp"].time, lat=ds_all["iwc_cp"].lat,
                         lon=ds_all["iwc_cp"].lon, method="nearest")
        cpT_ri_dar = cpT_ri.sel(time=ds_all["iwc_cp"].time, lat=ds_all["iwc_cp"].lat, 
                                lon=ds_all["iwc_cp"].lon, method="nearest")
        
        
        diffs_vals = tb_dar.values - cpT_ri_dar.values
        diffs_dar = xr.DataArray(diffs_vals, dims=["time"], coords={"time": ds_all.time})
    else:
        ds_all, tb = get_data(months, year, region, file_path=file_path)
        tb_dar = tb.sel(
            time=ds_all["iwc_cp"].time, lat=ds_all["iwc_cp"].lat, lon=ds_all["iwc_cp"].lon, 
            method="nearest"
        )
    
        
    # bin qi and re at cold point-relative levels by Tb or (Tb-cpT)
    if bin_by_diffs:        
        output_qi = bin_var_by_diffs(ds_all, diffs_dar, "qi", save_dicts=save_dicts, 
                                     pickle_dir=save_pickle_path,
                                   month=months, year=year, region=region)
        output_re = bin_var_by_diffs(ds_all, diffs_dar, "Re", save_dicts=save_dicts, 
                                     pickle_dir=save_pickle_path,
                                   month=months, year=year, region=region, iflag=3)
        output_qi3 = bin_var_by_diffs(ds_all, diffs_dar, "qi", save_dicts=save_dicts, 
                                      pickle_dir=save_pickle_path,
                                   month=months, year=year, region=region, iflag=3)
    else:
        output_qi = bin_var_by_tb(ds_all, tb_dar, "qi", save_dicts=save_dicts, pickle_dir=save_pickle_path,
                                   month=months, year=year, region=region)
        output_re = bin_var_by_tb(ds_all, tb_dar, "Re", save_dicts=save_dicts, pickle_dir=save_pickle_path,
                                   month=months, year=year, region=region, iflag=3)
        output_qi3 = bin_var_by_tb(ds_all, tb_dar, "qi", save_dicts=save_dicts, pickle_dir=save_pickle_path,
                                   month=months, year=year, region=region, iflag=3)
    
    # load in the dictionaries (to check that they saved properly)
    qi_means, qi_counts = load_saved_dicts("qi", months, region, year, save_pickle_path)
    re_means, re_counts = load_saved_dicts("Re", months, region, year, save_pickle_path, iflag=3)
    qi3_means, qi3_counts = load_saved_dicts("qi", months, region, year, save_pickle_path, iflag=3)

    if bin_by_diffs:
        plot_binned_by_diffs(qi_means, diffs_dar, "mean", months, region, year,
                              save=save_plots, save_dir=save_fig_path)
        plot_binned_by_diffs(qi_counts, diffs_dar, "count", months, region, year,
                              save=save_plots, save_dir=save_fig_path)
        plot_binned_by_diffs(qi3_means, diffs_dar, "mean", months, region, year, save_extr="_both_instr",
                              save=save_plots, save_dir=save_fig_path, stat_lims=(3e-5, 1e-3))
        plot_binned_by_diffs(qi3_counts, diffs_dar, "count", months, region, year, save_extr="_both_instr",
                              save=save_plots, save_dir=save_fig_path)
        plot_binned_by_diffs(re_means, diffs_dar, "mean", months, region, year, save_extr="_both_instr",
                          save=save_plots, save_dir=save_fig_path, varname="reff", lognorm=False)

    else:
        plot_binned_by_tb(qi_means, tb_dar, "mean", months, region, year,
                              save=save_plots, save_dir=save_fig_path)
        plot_binned_by_tb(qi_counts, tb_dar, "count", months, region, year,
                              save=save_plots, save_dir=save_fig_path)
        plot_binned_by_tb(qi3_means, tb_dar, "mean", months, region, year, save_extr="_both_instr",
                              save=save_plots, save_dir=save_fig_path, stat_lims=(3e-5, 1e-3))
        plot_binned_by_tb(qi3_counts, tb_dar, "count", months, region, year, save_extr="_both_instr",
                              save=save_plots, save_dir=save_fig_path)
        plot_binned_by_tb(re_means, tb_dar, "mean", months, region, year, save_extr="_both_instr",
                          save=save_plots, save_dir=save_fig_path, varname="reff", lognorm=False)

    
if __name__ == "__main__":
    main()            
                        

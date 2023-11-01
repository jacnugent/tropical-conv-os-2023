# tropical-conv-os-2023
[![DOI](https://zenodo.org/badge/654852706.svg)](https://zenodo.org/badge/latestdoi/654852706)

Code used in [Nugent and Bretherton (2023), _Geophysical Research Letters_](https://doi.org/10.1029/2023GL105083), on the observed distribution of convection that overshoots the cold point. 

The data used in this paper is not included here and must be [downloaded](https://github.com/jacnugent/tropical-conv-os-2023/#data-download) and [processed](https://github.com/jacnugent/tropical-conv-os-2023/#data-processing) following the steps listed below. The code used for the [analysis and figures](https://github.com/jacnugent/tropical-conv-os-2023/#analysis-and-figures) is primarily in Jupyter notebooks. Note that the scripts are not designed to be downloaded and run immediately; they must first be edited as specified below.

The South Pacific Convergence Zone (SPC) region crosses the International Date Line. For much of the data processing and analysis, this region is split into two halves to avoid any issues with the +/-180° longitude. In the code, “SPC1” refers to the eastern half (165°-180° E) and “SPC2” refers to the western half (180°-145° W), while “SPC” refers to the entire region (i.e., SPC1 and SPC2 concatenated by longitude).

## Data Download
### GPM_MERGIR
Downloads one season/year/region at a time. You must first install the CDS API key/client (see [here]( https://cds.climate.copernicus.eu/api-how-to) for instructions).
1. For each season/year/region, download the list of file links from [NASA GES DISC]( https://disc.gsfc.nasa.gov/datasets/GPM_MERGIR_1/summary). Rename these files with convention “subset_MMMYYYY_RRR.txt” (e.g., “subset_DJF2009_AMZ.txt”). 
2. Edit [get_mergir.sh](shell_scripts/get_mergir.sh) to change the file paths and season/year/
3. Run [get_mergir.sh](shell_scripts/get_mergir.sh) to download each file via wget. This loops through all regions in one year/season and concatenates the files into a single .nc4 file for each region.
4. Repeat for the next season/year.


### ERA5
Downloads one season/year at a time and then subsets into regions. You must first install the CDS API key/client (see [here]( https://cds.climate.copernicus.eu/api-how-to) for instructions). Note that this procedure can take over 24 hours to download three months of data because the files must be retrieved off of tapes.
1.	Edit [get_era5_ml.sh](shell_scripts/get_era5_ml.sh) to change the file paths and the season/year. 
2.	Run [get_era5_ml.sh](shell_scripts/get_era5_ml.sh) to download monthly files for the global tropics and split by variable. This runs the python script [get_era5_climo_ml.py](python_scripts/get_era5_climo_ml.py).
3.	Edit [process_era5_ml.sh](shell_scripts/process_era5_ml.sh) to change the file paths and specify the season/year.
4.	Run [process_era5_ml.sh](shell_scripts/process_era5_ml.sh) to subset the files into regions, compute the geopotential at each model level (runs [compute_geopotential_on_ml.py](python_scripts/compute_geopotential_on_ml.py), written by [ECMWF]( https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height#heading-Geopotentialonmodellevels)), convert files into netcdfs, and concatenate the monthly files into one file for that season. 50 GB memory is recommended for this step.
5.	Repeat for the next season/year.

### DARDAR
Downloads one season/year at a time and then subsets into regions. You must first register with the [AERIS/ICARE Data and Services Center](https://www.icare.univ-lille.fr/) to access the data archive.
1.	Edit [get_dardar_v3.sh](shell_scripts/get_dardar_v3.sh) with your username/password, change the file paths, and specify the season and year to download.
2.	Run [get_dardar_v3.sh](shell_scripts/get_dardar_v3.sh) to download the global files for the season/year via lftp.
3.	Edit [process_dardar_v3.sh](shell_scripts/process_dardar_v3.sh) to change the file paths and specify the year/season.
4.	Run [process_dardar_v3.sh](shell_scripts/get_dardar_v3.sh)  to subset the global file into one file for each region. This runs the python script [process_dardar_v3.py](python_scripts/process_dardar_v3.py). 25 GB memory is recommended for this step.
5.	Delete the global files and repeat for the next season/year.

## Data Processing
### Regridding and Binning
Regrids the GPM_MERGIR brightness temperature and ERA5 temperature/height data onto the DARDAR grid and bins the IWC by brightness temperature at cold point-relative levels. This processes all regions in one season/year at a time.
1.	Edit [regrid_bin_plot.sh]( shell_scripts/regrid_bin_plot.sh) to change the file paths and the season/year. 
2.	Run [regrid_bin_plot.sh](shell_scripts/regrid_bin_plot.sh). This runs the python script [regrid_data_cp.py](python_scripts/regrid_data_cp.py) to regrid and then [bin_obs_overshoot.py](python_scripts/bin_obs_overshoot.py) to do the binning. 50 GB memory is recommended for this step. This also saves “standard” versions of the plots for each year/season; note that these are NOT the same plots used in the paper. 
3.	Repeat for the next season/year.

### Joint Histograms
Calculates the joint brightness temperature-cold point histograms for each season/year. Saves a python dictionary with the histogram counts, bins, etc. to a pickle file.
1.	Edit [get_Tb-cpT_hist.sh](shell_scripts/get_Tb-cpT_hist.sh) to change the file paths and specify the season. 
2.	Run [get_Tb-cpT_hist.sh](shell_scripts/get_Tb-cpT_hist.sh). This runs the python script [cold_point_reindex.py](python_scripts/cold_point_reindex.py) to regrid the cold point temperature file _for SPC only_ (needed because of the date line issue) and [biv_hist.py](python_scripts/biv_hist.py) to calculate and save the histogram dictionary. 50 GB memory is recommended for this step.
3.	Repeat for the next season.


## Analysis and Figures
#### Figures 1, 2, S2, and S3
1. Get the counts for the $T_b-T_{cp}$ bins in [obs_climo_Tb-cpT_hists.ipynb](jupyter_notebooks/obs_climo_Tb-cpT_hists.ipynb).
2. Bin by $T_b-T_{cp}$ and make the figures in [obs_climo_paper_diffs_binned_plots.ipynb](jupyter_notebooks/obs_climo_paper_diffs_binned_plots.ipynb).

#### Figures 3 and S4
1. The joint histogram counts are saved as pickle files in [biv_hist.py](python_scripts/biv_hist.py).
2. Calculate the conditional probabilities of overshoots and make the figures in [obs_climo_cond_prob_joint_hists_plot.ipynb](jupyter_notebooks/obs_climo_cond_prob_joint_hists_plot.ipynb).

#### Figure 4
1. Calculate the frequencies of cold point overshoots and save as netcdf files in [obs_climo_calc_os_freqs.ipynb](jupyter_notebooks/obs_climo_calc_os_freqs.ipynb).
2. Make the figures and calculate the mean frequencies (including with alternate thresholds) in [obs_climo_overshooting_heatmaps.ipynb](jupyter_notebooks/obs_climo_overshooting_heatmaps.ipynb).

#### Figure S1
1. Make the time-mean cold point temperature and height files in [get_obs_mean_cpT.ipynb](jupyter_notebooks/get_obs_mean_cpT.ipynb).
2. Make the figure in [obs_climo_time_mean_cp_maps.ipynb](jupyter_notebooks/obs_climo_time_mean_cp_maps.ipynb).

#### Convective vs. other stratospheric cirrus fractions
Calculate the fractions in [obs_climo_cirrus_fracs.ipynb](jupyter_notebooks/obs_climo_cirrus_fracs.ipynb).

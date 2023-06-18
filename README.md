# tropical-conv-os-2023
Code used in Nugent and Bretherton 2023 (in preparation) on the observed distribution of convection that overshoots the cold point. 

The data used in this paper is not included here and must be downloaded (link) and processed (link) following the steps listed below. The code used for the analysis and figures (link) is primarily in Jupyter notebooks. Note that the scripts are not designed to be downloaded and run immediately; they must first be edited as specified below.

The South Pacific Convergence Zone (SPC) region crosses the International Date Line. For much of the data processing and analysis, this region is split into two halves to avoid any issues with the $\pm$180° longitude. In the code, “SPC1” refers to the eastern half (165°-180° E) and “SPC2” refers to the western half (180°-145° W), while “SPC” refers to the entire region (i.e., SPC1 and SPC2 concatenated by longitude).

## Data Download
### GPM_MERGIR
Downloads one season/year/region at a time.
(_Description coming soon._)

### ERA5
Downloads one season/year at a time and then subsets into regions. You must first install the CDS API key/client (see [here]( https://cds.climate.copernicus.eu/api-how-to) for instructions).
(_Description coming soon_)

### DARDAR
Downloads one season/year at a time and then subsets into regions. You must first register with the [AERIS/ICARE Data and Services Center](https://www.icare.univ-lille.fr/) to access the data archive.
1.	Edit [get_dardar_v3.sh](shell_scripts/get_dardar_v3.sh) with your username/password, change the file paths, and specify the season and year to download.
2.	Run [get_dardar_v3.sh](shell_scripts/get_dardar_v3.sh) to download the global files for the season/year via lftp.
3.	Edit [process_dardar_v3.sh](shell_scripts/process_dardar_v3.sh) to change the file paths and specify the year/season.
4.	Run [get_dardar_v3.sh](shell_scripts/get_dardar_v3.sh)  to subset the global file into one file for each region. This runs the python script [process_dardar_v3.py](python_scripts/process_dardar_v3.py). 25 GB memory is recommended for this step.
5.	Delete the global files and repeat for the next season/year.

## Data Processing
### Regridding and Binning
Regrids the GPM_MERGIR brightness temperature and ERA5 temperature/height data onto the DARDAR grid and bins the IWC by brightness temperature at cold point-relative levels. This processes all regions in one season/year at a time.
1.	Edit [regrid_bin_plot.sh]( shell_scripts/regrid_bin_plot.sh) to change the file paths and the season/year. 
2.	Run [regrid_bin_plot.sh](shell_scripts/regrid_bin_plot.sh). This runs the python script [regrid_data_cp.py](python_scripts/regrid_data_cp.py) to regrid and then [bin_obs_overshoot.py](python_scripts/bin_obs_overshoot.py) to do the binning. 50 GB memory is recommended for this step.
  a.	This also saves “standard” versions of the plots for each year/season; note that these are NOT the same plots used in the paper. 
  b.	(_Description of the netcdf files this produces coming soon_).
3.	Repeat for the next season/year.

### Joint Histograms
Calculates the joint brightness temperature-cold point histograms for each season/year. Saves a python dictionary with the histogram counts, bins, etc. to a pickle file.
1.	Edit [get_Tb-cpT_hist.sh](shell_scripts/get_Tb-cpT_hist.sh) to change the file paths and specify the season. 
2.	Run [get_Tb-cpT_hist.sh](shell_scripts/get_Tb-cpT_hist.sh). This runs the python script [cold_point_reindex.py](python_scripts/cold_point_reindex.py) to regrid the cold point temperature file _for SPC only_ (needed because of the date line issue) and [biv_hist.py](python_scripts/cold_point_reindex.py) to calculate and save the histogram dictionary.
3.	Repeat for the next season.


## Analysis and Figures
(_Description coming soon._)

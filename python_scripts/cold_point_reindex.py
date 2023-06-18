"""
cold_point_reindex.py

Script to get the cold point temperature and
reindex it onto the MERGIR brightness temperature
grid. Do for one year/region/season at a time.

-----------------------
usage: cold_point_reindex.py [-h] -y YEAR -f FILE_PATH -r REGION -m MONTHS [-o OUT_PATH]

optional arguments:
  -h, --help            show this help message and exit
  -y YEAR, --year YEAR  Year
  -f FILE_PATH, --file_path FILE_PATH
                        Input path for files
  -r REGION, --region REGION
                        Region abbreviation
  -m MONTHS, --months MONTHS
                        Months abbreviation (e.g., DJF)
  -o OUT_PATH, --out_path OUT_PATH
                        Output path to save files (default=file_path)
-----------------------

Assumes the following file name convention:
* temperature: ERA5_T_0.25deg_ml_12-20km_{m}{y}_{r}.nc
* brightness temp: MERGIR_Tb_4km_{m}{y}_{r}.nc4

Saves the output file as
    ERA5_cpT_reindexed_{m}{y}_{r}.nc
"""
import argparse
import dask
import xarray as xr


# horizontal resolution (degrees) of temperature grid
CPT_RES = 0.25


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
    parser.add_argument("-o", "--out_path", help="Output path to save files  (default=file_path)")

    return parser.parse_args()


def cold_point_temp(temp):
    """ Return data array of the cold point temperature
    """
    cpT_inds_chunked = temp.chunk("auto").argmin(dim="level")
    cpT_inds = cpT_inds_chunked.compute()
    cpT = temp.isel(level=cpT_inds)
    return cpT


def reindex_cpT(cpT, tb):
    """ Reindex cold point temperature to match brightness temp grid
    """
    cpT_ri_chunked = cpT.chunk({"time": 1}).reindex({"lat": tb.lat, "lon": tb.lon}, method="nearest")
    cpT_ri = cpT_ri_chunked.compute()
    return cpT_ri


def save_cpT(cpT_ri, out_path, months, year, region, cpT_res=CPT_RES):
    """ Save the reindexed cold point file as a netcdf
    """
    out_name = "ERA5_cpT_reindexed_{m}{y}_{r}.nc".format(m=months, y=year, r=region)
    cpT_ri.attrs["cold_point_grid"] = "{} deg".format(cpT_res)
    cpT_ri.to_netcdf(out_path + out_name)


def main():
    """ Calculate cold point, reindex grid, and save the netcdf
    """
    args = parse_args()
    print(args)

    year = args.year
    file_path = args.file_path
    region = args.region
    months = args.months
    out_path = args.out_path
    
    if file_path[-1] != "/":
        file_path = file_path + "/"
    if out_path is None:
        out_path = file_path
    else:
        if out_path[-1] != "/":
            out_path = out_path + "/"

    temp_file = "ERA5_T_0.25deg_ml_12-20km_{m}{y}_{r}.nc".format(m=months, r=region, y=year)
    tb_file = "MERGIR_Tb_4km_{m}{y}_{r}.nc4".format(m=months, y=year, r=region)
    temp = xr.open_dataset(file_path + temp_file)["t"]
    tb = xr.open_dataset(file_path + tb_file)["Tb"]
    
    cpT = cold_point_temp(temp).rename({"latitude": "lat", "longitude": "lon"})
    cpT_ri = reindex_cpT(cpT, tb)
    save_cpT(cpT_ri, out_path, months, year, region)    
    

if __name__ == "__main__":
    main()         

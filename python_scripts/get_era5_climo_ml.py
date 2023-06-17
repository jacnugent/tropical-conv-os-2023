"""
get_era5_climo_ml.py

This gets [temperature (12-20km)] and/or [temperature & spec hum (all levs), spec hum (all levs) and z & log of sfc pressure] for some region and some time frame.

------------------------------------------------------
usage: get_era5_climo_ml.py [-h] -s START -e END -o OUT_PATH -c COORDS
                            [COORDS ...] [-r REGION] [-v VARIABLE]

optional arguments:
  -h, --help            show this help message and exit
  -s START, --start START
                        start date
  -e END, --end END     end date
  -o OUT_PATH, --out_path OUT_PATH
                        path to output files
  -c COORDS [COORDS ...], --coords COORDS [COORDS ...]
                        coordinates in the format "N, W, S, E", e.g., 30 -180
                        -30 180
  -r REGION, --region REGION
                        region name
  -v VARIABLE, --variable VARIABLE
                        variable name (combined, temp, or combined_and_temp)
------------------------------------------------------
                       
Will raise an exception if the time frame spans more than one month.

Get info from here to make changes to the request: https://apps.ecmwf.int/data-catalogues/era5/?class=ea

See instructions here to set up the .cdsapirc file: 
https://cds.climate.copernicus.eu/api-how-to#install-the-cds-api-key
"""

import cdsapi
import sys
import argparse

from datetime import datetime


def parse_args():
    """ Parse command-line arguments
    """
    parser = argparse.ArgumentParser()

    # required
    parser.add_argument("-s", "--start", help="start date", required=True)
    parser.add_argument("-e", "--end", help="end date", required=True)
    parser.add_argument("-o", "--out_path", help="path to output files", required=True)
    parser.add_argument("-c", "--coords", nargs="+", help="coordinates in the format \"N, W, S, E\", e.g., 30 -180 -30 180", type=float, required=True)

    # optional
    parser.add_argument("-r", "--region", default="GT", help="region name")
    parser.add_argument("-v", "--variable", default="combined_and_temp", help="variable name (combined, temp, temp_and_qv, or combined_and_temp)")

    args = parser.parse_args()

    return args


def retrieve_temp(year, dates, region, coords, months, out_dir):
    """ Send the cdsapi request for temp, 12-20 km
    """
    c = cdsapi.Client()
    out_name = out_dir + 'ERA5_T_0.25deg_ml_12-20km_{m}{y}_{r}.nc'.format(y=year, m=months, r=region)

    # 49-78 is ~12-20 km
    c.retrieve('reanalysis-era5-complete', {
        'class': 'ea',
        'date': dates,
        'expver': '1',
        'levelist': '49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78',
        'levtype': 'ml',
        'param': '130',
        'stream': 'oper',
        'time': '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00',
        'type': 'an',
        'format': 'netcdf',
        'grid': '0.25/0.25',
        'area': coords,
   }, out_name )


def retrieve_tq(year, dates, region, coords, months, out_dir):
    """ Send the cdsapi request for temp & qv, all levels
    """
    c = cdsapi.Client()
    out_name = out_dir + 'tq_{m}{y}_{r}.grib'.format(y=year, m=months, r=region)

    c.retrieve('reanalysis-era5-complete', {
        'class': 'ea',
        'date': dates,
        'expver': '1',
        'levelist': '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137',
        'levtype': 'ml',
        'param': '130/133',
        'stream': 'oper',
        'time': '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00',
        'type': 'an',
        'grid': '0.25/0.25',
        'format': 'grib',
        'area': coords,
    }, out_name )


def retrieve_combined(year, dates, region, coords, months, out_dir):
    """ 
    Send the cdsapi request for temp, qv, and zlnsp, all levels
    **Will need to subset in cdo after (zlnsp vars on level 137 only, Tq in sep file)
    """
    c = cdsapi.Client()
    out_name = out_dir + 'tq_zlnsp_{m}{y}_{r}.grib'.format(y=year, m=months, r=region)

    c.retrieve('reanalysis-era5-complete', {
        'class': 'ea',
        'date': dates,
        'expver': '1',
        'levelist': '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137',
        'levtype': 'ml',
        'param': '130/133/129/152',
        'stream': 'oper',
        'time': '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00',
        'type': 'an',
        'grid': '0.25/0.25',
        'format': 'grib',
        'area': coords,
    }, out_name )


def retrieve_zlnsp(year, dates, region, coords, months, out_dir):
    """ Send the cdsapi request for log of surface pressure
    """
    c = cdsapi.Client()
    out_name = out_dir + 'zlnsp_{m}{y}_{r}.grib'.format(y=year, r=region, m=months)

    c.retrieve('reanalysis-era5-complete', {
        'class': 'ea',
        'date': dates,
        'expver': '1',
        'levelist': '1',
        'levtype': 'ml',
        'param': '129/152',
        'stream': 'oper',
        'time': '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00',
        'type': 'an',
        'grid': '0.25/0.25',
        'format': 'grib',
        'area': coords,
    }, out_name )


def main():
    """ Get year/dates/region from command line and run the request
    """
    args = parse_args()
    print(args)

    start = args.start
    end = args.end
    out_path = args.out_path
    coords = args.coords
    region = args.region
    variable = args.variable
    date_string = str(start) + "/to/" + str(end)

    # check that the requested time frame is one month or less
    start_mo = datetime.strptime(start, "%Y-%m-%d").strftime("%b")
    end_mo = datetime.strptime(end, "%Y-%m-%d").strftime("%b")
    start_y = datetime.strptime(start, "%Y-%m-%d").strftime("%Y")
    end_y = datetime.strptime(end, "%Y-%m-%d").strftime("%Y")
    if (start_mo != end_mo) or (start_y != end_y):
        raise Exception("Must get data for only one month in one year at a time!")
    else:
        months = start_mo
        year = start_y

    # submit the request(s)
    if variable == "combined":
        retrieve_combined(year, date_string, region, coords, months, out_path)
    elif variable == "combined_and_temp":
        retrieve_combined(year, date_string, region, coords, months, out_path)
        retrieve_temp(year, date_string, region, coords, months, out_path)
    elif variable == "temp": 
        retrieve_temp(year, date_string, region, coords, months, out_path)
    elif variable == "temp_and_qv": 
        retrieve_tq(year, date_string, region, coords, months, out_path)
    else:
        raise Exception("Variable \'{v}\' unknown! Must be \'temp\', \'combined\', \'temp_and_qv\', or \'combined_and_temp\'".format(v=variable))


if __name__ == "__main__":
    main()

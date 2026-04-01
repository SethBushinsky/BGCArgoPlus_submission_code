import os
import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats
import warnings

def month_grid(mon_dir, csv_output_dir, var, bins, targetz, targety, targetx, description, encoding, gridded_dir):
    if mon_dir.__contains__('DS'):
        return
    year = int(mon_dir.split('_')[0])
    month = int(mon_dir.split('_')[1])
    # year = 2020
    # month = 1

    # find all float files for that month
    month_path = csv_output_dir + f"{year}_{month:02d}" + '/'

    float_month_files = [f for f in os.listdir(month_path) if f.endswith('.csv')]
    print(len(float_month_files), 'float files found for', year, month)

    # initialize list, then loop through, read in each csv, and finally concatenate into one dataframe 
    float_month_all = []
    for file in float_month_files:
        mon_n = pd.read_csv(month_path + file)
        if not mon_n.empty:  # Only append non-empty DataFrames
            float_month_all.append(mon_n)
    
    # Only concatenate if we have non-empty DataFrames
    if float_month_all:
        float_month_all = pd.concat(float_month_all, ignore_index=True)
    else:
        print(f"No valid data for {year}-{month}, skipping...")
        return


    float_data = float_month_all.copy()
    spatial_data = [float_data['PRES'].values.astype(float), 
                    float_data.LATITUDE.values.astype(float), 
                    float_data.LONGITUDE.values.astype(float)]
    if len(spatial_data[0])==0:
        return
    to_bin = [float_data['TEMP'].values.astype(float), 
            float_data['PSAL'].values.astype(float), 
            float_data[var].values.astype(float)]
    
    # Suppress the "Mean of empty slice" warning
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        float_statistic_corrected, bin_edges, binnumber = stats.binned_statistic_dd(sample = spatial_data, values=to_bin, bins=bins, statistic=np.nanmean)

    # Create data_vars dict dynamically
    data_vars = {
        "Temperature": (["z", "lat", "lon"], float_statistic_corrected[0]),
        "Salinity": (["z", "lat", "lon"], float_statistic_corrected[1]),
        var: (["z", "lat", "lon"], float_statistic_corrected[2])
    }

    ds = xr.Dataset(
        data_vars=data_vars,
        coords=dict(
            z=targetz,
            lat=targety,
            lon=targetx
        ),
        attrs=dict(
            description=description
        )
    )
    ds = ds.expand_dims(time=[pd.Timestamp(year=year, month=month, day=15)])
    ds.to_netcdf(f"{gridded_dir}{var}_Gridded_{year}_{month}.nc", encoding=encoding)
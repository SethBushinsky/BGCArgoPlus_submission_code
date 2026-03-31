import xarray as xr
import pandas as pd
import numpy as np
import os

def float_month_separation(var, var_suffix, file, main_argo_dir, intermediate_gridding_dir):
    # open a processed float file, look for variables in a passed list, save into monthly files for each variable
    # file = '4903622_Sprof_BGCArgoPlus_full.nc'  # test file

    argo_n = xr.open_dataset(main_argo_dir + file)

    # check for data and validity
    var_list = argo_n.data_vars.keys()

    if var in var_list: # check if variable exists in dataset
        if np.isnan(argo_n[var].values).all():
            print(f'no data, skipping {var} for {file}')
            return
    else:
        print(f'{var} not found in {file}')
        return

    # get years and months present in data file
    time_vals = pd.to_datetime(argo_n['JULD'].values)
    years = np.unique(time_vals.year)

    # for each year, for each month with data in that year, save out a pandas dataframe as a csv with the data needed for gridding
    for year in years:
        months = np.unique(time_vals[time_vals.year == year].month)
        for month in months:
            # print(month)
            month_inds = np.where((time_vals.year == year) & (time_vals.month == month))[0]
            # month_inds
            
            if len(month_inds) == 0:
                continue
            argo_month = argo_n.isel(N_PROF=argo_n['N_PROF'][month_inds])
            month_df = pd.DataFrame()
            # for each profile in argo_month, create a dataframe, then concatenate all dataframes together
            for p in range(len(argo_month['N_PROF'])):
                if np.all(np.isnan(argo_month['PRES' + var_suffix].values[p,:])):
                    # print('Skipping profile', p, 'due to missing PRES' + var_suffix)
                    continue
                df_dict = {'Float Number': str(argo_month['WMO_ID'].values), 
                    'Datetime': pd.to_datetime(argo_month['JULD'].values[p]).strftime('%Y-%m-%d %H:%M:%S'),
                    'LATITUDE': argo_month['LATITUDE'].values[p],
                    'LONGITUDE': argo_month['LONGITUDE'].values[p],
                    'PRES': argo_month['PRES' + var_suffix].values[p,:],
                    'TEMP': argo_month['TEMP' + var_suffix].values[p,:],
                    'PSAL': argo_month['PSAL' + var_suffix].values[p,:],
                    'sigma': argo_month['sigma0'].values[p,:],
                    'gamma': argo_month['gamma'].values[p,:],
                    var: argo_month[var].values[p,:]
                }
                profile_df = pd.DataFrame(data=df_dict)
                month_df = pd.concat((month_df, profile_df), axis=0, ignore_index=True)

            # if month_df is empty, skip
            if month_df.empty:
                print(f'No valid data for {year}-{month} for {file}, skipping...')
                continue
            # drop rows with NaNs in the variable of interest
            month_df = month_df.dropna(subset=[var])
            # save to csv
            out_filename = f"{year}_{month:02d}_{var}_{str(argo_n['WMO_ID'].values)}.csv"
            # month_df.to_csv(out_filename, index=False)
            print(f'Saved {out_filename}')
            # if year_month_dir does not exist, create it
            month_path = intermediate_gridding_dir + var + '/' + f"{int(year)}_{month:02d}" + '/'
            if not os.path.exists(month_path):
                os.makedirs(month_path)

            month_df.to_csv(month_path + out_filename, index=False)
    return
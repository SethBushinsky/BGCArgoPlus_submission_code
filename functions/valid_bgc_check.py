import xarray as xr
import numpy as np
import pandas as pd
import datetime

def check_for_valid_bgc(output_dir, file, bgc_vars, local_outlier_dir):
    """
    Check if there is any valid BGC data in the filtered Argo files
    If not, saves out an outlier detection .csv file with 'XXX' as the researcher name 

    """
    argo_n = xr.open_dataset(output_dir + file)
    bgc_present=False
    researcher = 'XXX'

    for var in bgc_vars:
        if var not in argo_n.keys():
            continue
        else: # Check if there are any valid bgc data. If not, don't bother doing outlier detection. 
            if np.sum(~np.isnan(argo_n[var]))>0:
                bgc_present = True
    
    if bgc_present==False:
        # print('No valid bgc data in variables checked, skipping outlier removal')
        empty_df = pd.DataFrame(columns = ["Status"], data=['No valid bgc'])
        savename = local_outlier_dir + 'outliers_' + file[0:-3] + '_' + researcher
        current_time = datetime.datetime.now()
        current_time_str = str(current_time.year) + '_' + str(current_time.month) + '_' + str(current_time.day) + '_' + str(current_time.hour)
        empty_df.to_csv(savename + '_' + current_time_str + '.csv', mode='w', index=False, header=True)

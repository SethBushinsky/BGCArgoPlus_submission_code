#General Imports
import numpy as np
import xarray as xr
import pandas as pd
import warnings

def create_float_df_2(float_dir, argo_file, use_outlier_removed = 0, use_corrected = 0, use_bgc_argo_plus = 0):
    if use_outlier_removed == 0 and use_corrected == 0:
        PRES = "PRES_ADJUSTED"
        TEMP = "TEMP_ADJUSTED"
        PSAL = "PSAL_ADJUSTED"
        DOXY = "DOXY_ADJUSTED"
        NITRATE = "NITRATE_ADJUSTED"
        PH = "PH_IN_SITU_TOTAL_ADJUSTED"
        ADD_ON = "_ADJUSTED"
    elif use_outlier_removed == 0 and use_corrected == 1:
        PRES = "PRES_ADJUSTED"
        TEMP = "TEMP_ADJUSTED"
        PSAL = "PSAL_ADJUSTED"
        DOXY = "DOXY_ADJUSTED_C"
        NITRATE = "NITRATE_ADJUSTED"
        PH = "PH_IN_SITU_TOTAL_ADJUSTED"
        ADD_ON = "_ADJUSTED"
    elif use_outlier_removed == 1 and use_corrected == 0:
        PRES = "PRES_ADJUSTED_RO"
        TEMP = "TEMP_ADJUSTED_RO"
        PSAL = "PSAL_ADJUSTED_RO"
        DOXY = "DOXY_ADJUSTED_RO"
        NITRATE = "NITRATE_ADJUSTED_RO"
        PH = "PH_IN_SITU_TOTAL_ADJUSTED_RO"
        ADD_ON = "_ADJUSTED_RO"
    else: # outlier = 1, corrected = 1
        PRES = "PRES_ADJUSTED_RO"
        TEMP = "TEMP_ADJUSTED_RO"
        PSAL = "PSAL_ADJUSTED_RO"
        DOXY = "DOXY_ADJUSTED_ROC"
        NITRATE = "NITRATE_ADJUSTED_RO"
        PH = "PH_IN_SITU_TOTAL_ADJUSTED_RO"
        ADD_ON = "_ADJUSTED_RO"
    if use_bgc_argo_plus == 1:
        PRES = "PRES_ADJUSTED_BGCArgoPlus"
        TEMP = "TEMP_ADJUSTED_BGCArgoPlus"
        PSAL = "PSAL_ADJUSTED_BGCArgoPlus"
        DOXY = "DOXY_ADJUSTED_BGCArgoPlus"
        NITRATE = "NITRATE_ADJUSTED_BGCArgoPlus"
        PH = "PH_IN_SITU_TOTAL_ADJUSTED_BGCArgoPlus"
        ADD_ON = "_ADJUSTED_BGCArgoPlus"


    float_profile_data = pd.DataFrame(columns=['Float Number', 'Latitude', 'Longitude', 'Profile Number', 'Datetime', 'O2 Present', 'Pressure', 'Temperature', 'Salinity', 'Oxygen'])
    wmo = argo_file[0:7]
    print("Processing float " + wmo)
    with xr.open_dataset(float_dir+argo_file) as float_data:
        # flag_fp = '/Users/znachod/UHM_Ocean_BGC_Group Dropbox/Datasets/Data_Products/BGC_ARGO_GLOBAL/2025_01_24/processed/'
        # flag_file = str(wmo) + '_Sprof_flags_mode_only.nc'
        # flag_data = xr.open_dataset(flag_fp + flag_file)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            for profile_num, data in float_data.groupby('N_PROF'):
                if b'DOXY                                                            ' in data.PARAMETER.values[0]:
                    O2Present = "Y"
                else:
                    O2Present = "N"
                lon = data.LONGITUDE.values
                lat = data.LATITUDE.values
                t = pd.to_datetime(str(np.datetime64(data.JULD.values)))
                # flag_data_prof = flag_data.sel(N_PROF=profile_num)
                try:
                    timestring = t.strftime('%Y.%m.%d')
                    # print(data.JULD.values())
                    # date = pd.to_datetime(data.JULD, origin='julian', unit='D')
                    if O2Present == "N":
                        new_row = pd.DataFrame(data={'Float Number':wmo, 'Latitude': lat, 'Longitude': lon, 'Profile Number':profile_num, 'Datetime':t, 'O2 Present': O2Present, 'Pressure':[data[PRES].values], 'Temperature':[data[TEMP].values], 'Salinity':[data[PSAL].values], 'Oxygen':[np.nan] * len(data[PSAL].values)}, index=[0])
                    else:
                        new_row = pd.DataFrame(data={'Float Number':wmo, 'Latitude': lat, 'Longitude': lon, 'Profile Number':profile_num, 'Datetime':t, 'O2 Present': O2Present, 'Pressure':[data[PRES].values], 'Temperature':[data[TEMP].values], 'Salinity':[data[PSAL].values], 'Oxygen':[data[DOXY].values]}, index=[0])
                    float_profile_data = pd.concat((float_profile_data, new_row), axis=0, ignore_index=True)
                except ValueError:
                    pass
            # flag_data.close()
    return float_profile_data
from datetime import datetime
from tabnanny import verbose

import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import os.path
import gsw 
import traceback
import os

def append_flag(argo_n, flag_var, suffix):
    # centralizes the append+assign_attrs so long_name can't be dropped by a later suffix append
    long_name = argo_n[flag_var].long_name
    argo_n[flag_var] = argo_n[flag_var].values.item() + suffix
    argo_n[flag_var] = argo_n[flag_var].assign_attrs(long_name=long_name)
    return argo_n

def plot_processed_files(file, dir, processed_fig_dir):

    argo_n = xr.open_dataset(dir + file)

    argo_n['decimal_year'] = (['N_PROF'],np.empty(argo_n.PRES_ADJUSTED.shape[0])) #nprof 
    argo_n.decimal_year[:] = np.nan
    date_time = pd.to_datetime(argo_n.JULD.values)
    year = date_time.year
    decimal_year = year + (date_time.day_of_year - 1) / 365.25
    argo_n.decimal_year[:] = decimal_year

    data_proj = ccrs.PlateCarree(central_longitude=0)
    map_proj = ccrs.Robinson(central_longitude=np.nanmean(argo_n['LONGITUDE']))

    # repeated_array = np.tile(argo_n.JULD, (len(argo_n.N_LEVELS), 1))

    # # Create an xarray.DataArray with the repeated array
    # data_array = xr.DataArray(repeated_array.T, dims=('time', 'N_LEVELS'))

    list_var = list(argo_n.keys())
    plot_vars = [field for field in list_var if field.endswith('_ADJUSTED')]

    # plot_vars = []
    # list_var = list(argo_n.keys())
    # for var in plot_vars_all:
    #     if var in list_var:
    #         plot_vars.append(var)

    # print('variables in: ' + file)
    # print(plot_vars)
    # Then, plot all kind of fields to see if there are outliers
    f = plt.figure(figsize=(40,5*len(plot_vars)))
    gs = f.add_gridspec(1+len(plot_vars), 2)

    ax0 = f.add_subplot(gs[0], projection=map_proj)
    ax0.set_global()
    ax0.coastlines()
    map = ax0.scatter(argo_n.LONGITUDE.values, argo_n.LATITUDE.values, c=argo_n.N_PROF, transform=data_proj, cmap='cool')
    plt.colorbar(map, label='Profile')
    
    x_lims = [np.nan, np.nan]

    x_lims = [np.min(decimal_year), np.max(decimal_year)]
    if x_lims[0] == x_lims[1]:
        x_lims[1] = x_lims[1]+10/365
   
    plt.rcParams.update({'font.size': 14})  # Change 14 to the desired font size

    for idx, var in enumerate(plot_vars):
        plot_exist = 0
        ax = f.add_subplot(gs[idx*2+2])
        if idx==0:
            plt.title(file + '_' + var)
        else:
            plt.title(var)

        if np.isnan(argo_n[var]).all():
            continue

        if var=='DOXY_ADJUSTED':
            color_map = 'plasma_r'
        elif var=='NITRATE_ADJUSTED':
            color_map = 'magma_r'
        elif var=='TEMP_ADJUSTED':
            color_map = 'cool'

        elif var=='PH_IN_SITU_TOTAL_ADJUSTED':
            color_map = 'cividis'
        elif var=='PH_25C_TOTAL_ADJUSTED':
            color_map = 'cividis'
        elif var=='PSAL_ADJUSTED':
            color_map = 'Wistia'
        elif var=='DIC':
            color_map = 'copper_r'
        elif var=='CHLA_ADJUSTED':
            color_map = 'spring'
        elif var=='BBP700_ADJUSTED':
            color_map = 'autumn'
        elif var=='CDOM_ADJUSTED':
            color_map = 'bone'
        else:
            color_map = 'viridis_r'
        #set color limits according to min/max of MLD properties
        if 'MLD' not in list_var:
            if ~np.isnan(argo_n[var].where(argo_n.PRES_ADJUSTED<150)).all():
                c_limit = [np.nanmin(argo_n[var].where(argo_n.PRES_ADJUSTED<200)), np.nanmax(argo_n[var].where(argo_n.PRES_ADJUSTED<200))]
                y_limit = [200,0]
            else:
                c_limit = [0, 1]
                y_limit = [100,0] 
        elif np.logical_and(np.nansum(~np.isnan(argo_n.MLD))==0, ~np.isnan(argo_n[var].where(argo_n.PRES_ADJUSTED<150)).all()):
            c_limit = [np.nanmin(argo_n[var].where(argo_n.PRES_ADJUSTED<150)), np.nanmax(argo_n[var].where(argo_n.PRES_ADJUSTED<150))]
            y_limit = [200, 0]
        elif ~np.isnan(argo_n[var].where(argo_n.PRES_ADJUSTED<150)).all(): 
            c_limit = [np.nanmin(argo_n[var].where(argo_n.PRES_ADJUSTED<np.nanmax(argo_n.MLD))), np.nanmax(argo_n[var].where(argo_n.PRES_ADJUSTED<np.nanmax(argo_n.MLD)))]
            y_limit = [np.nanmax(argo_n.MLD)+50, 0]
        else:
            c_limit = [0, 1]
            y_limit = [100,0]

        for p in range(0, len(argo_n.N_PROF)):
            p_p = argo_n.PRES_ADJUSTED[p,np.logical_and(~np.isnan(argo_n[var][p,:]), ~np.isnan(argo_n.PRES_ADJUSTED[p,:]))].values

            t_p = argo_n.decimal_year[p:p+2].values
            if t_p.size==1:
                t_p = np.tile(t_p, (2,1))
                t_p[1] = t_p[1] + (argo_n.decimal_year[p] - argo_n.decimal_year[p-1]).values
            elif np.isnan(t_p).any(): # if any values in t_p are nans
                if np.isnan(t_p).all(): # if all are nans, continue
                    continue
                elif np.isnan(t_p[0]):
                    t_p[0] = t_p[1] - 10/365
                else:
                    t_p[1] = t_p[0] + 10/365
            xl,yl = np.meshgrid(t_p, p_p)


            c = argo_n[var][p,np.logical_and(~np.isnan(argo_n[var][p,:]), ~np.isnan(argo_n.PRES_ADJUSTED[p,:]))].values
            if c.size==0:
                continue
            c = np.tile(c, (2,1))
            c = c.T
            try:
                plt.pcolormesh(xl, yl, c[0:-1,0:-1], cmap=color_map, shading='flat')
            except:
                print(file + ' ' + var + ' failed to plot')
                
            
            plt.clim(c_limit)
            plot_exist = 1
        plt.ylim(y_limit)
        if plot_exist==1:
            plt.colorbar(label=var)
        if 'MLD' in list_var:
            plt.plot(argo_n.decimal_year, argo_n.MLD, 'm')

        plt.xlim(x_lims)

        ax = f.add_subplot(gs[idx*2+3])
        #set color limits according to the min / max water below 500m

        if ~np.isnan(argo_n['PRES_ADJUSTED'].where(argo_n.PRES_ADJUSTED>500)).all():
            c_limit = [np.nanmin(argo_n[var].where(argo_n.PRES_ADJUSTED>500)), np.nanmax(argo_n[var].where(argo_n.PRES_ADJUSTED>500))]

        # axs1
        for p in range(0, len(argo_n.N_PROF)):
            p_p = argo_n.PRES_ADJUSTED[p,np.logical_and(~np.isnan(argo_n[var][p,:]), ~np.isnan(argo_n.PRES_ADJUSTED[p,:]))].values

            t_p = argo_n.decimal_year[p:p+2].values
            if t_p.size==1:
                t_p = np.tile(t_p, (2,1))
                t_p[1] = t_p[1] + (argo_n.decimal_year[p] - argo_n.decimal_year[p-1]).values
            elif np.isnan(t_p).any(): # if any values in t_p are nans
                if np.isnan(t_p).all(): # if all are nans, continue
                    continue
                elif np.isnan(t_p[0]):
                    t_p[0] = t_p[1] - 10/365
                else:
                    t_p[1] = t_p[0] + 10/365
            xl,yl = np.meshgrid(t_p, p_p)


            c = argo_n[var][p,np.logical_and(~np.isnan(argo_n[var][p,:]), ~np.isnan(argo_n.PRES_ADJUSTED[p,:]))].values
            if c.size==0:
                continue
            c = np.tile(c, (2,1))
            c = c.T
            try:
                plt.pcolormesh(xl, yl, c[0:-1,0:-1], cmap=color_map, shading='flat')    
            except:
                print(file + ' ' + var + ' failed to plot')

            plt.clim(c_limit)
        if plot_exist==1:
            plt.colorbar(label=var)
        if 'MLD' in list_var:
            plt.plot(argo_n.decimal_year, argo_n.MLD, 'm')
        if ~np.isnan(argo_n['PRES_ADJUSTED']).all():
            plt.ylim([np.nanmax(argo_n.PRES_ADJUSTED), 0])
        plt.xlim(x_lims)
    plt.tight_layout()
    plot_filename = file[0:-3]
    plt.savefig(f'{processed_fig_dir}{plot_filename}_v2.png')
    plt.close()
    plt.clf()

    argo_n.close()        

    return

def apply_sensor_info(argo_n, file, sprof_path, argo_index, verbose=False):
    argo_n["WMO_ID"] = argo_n.PLATFORM_NUMBER.values.astype(int)[0]
    argo_n["WMO_ID"] = argo_n["WMO_ID"].astype(dtype='object')
    argo_n["WMO_ID"] = argo_n["WMO_ID"].assign_attrs(long_name='WMO ID of the float, from the "PLATFORM_NUMBER" field in the Sprof file')
    if verbose: print(f'Loading meta file for argo_n wmo: {argo_n["WMO_ID"].values.item()}')
    meta_n = xr.open_dataset(sprof_path + file[0:7] + '_meta.nc')
    if verbose: print(f'Loaded meta file for argo_n wmo: {argo_n["WMO_ID"].values.item()}')
    # parameters_n = meta_n.PARAMETER.values
    air_cal_list = ['in air', 'in-air']
    confusing_air_cal_list = ['no in air', 'no in-air']

    # parameters_n = meta_n.PARAMETER.values
    # print(parameters_n)

    sensors = meta_n.SENSOR.values # using sensors instead of parameters, because parameters includes DOXY Temperature and length differs from sensor_models
    # print(sensors)

    sensor_models = meta_n.SENSOR_MODEL.values
    # print(sensor_models)
    if verbose: print(sensor_models)
    for i in range(len(sensors)):
        # print(sensors[i])
        # print(sensor_models[i])
        argo_n[sensors[i].decode('utf-8').strip() + '_model'] = sensor_models[i].decode('utf-8').strip()
        argo_n[sensors[i].decode('utf-8').strip() + '_model'] = argo_n[sensors[i].decode('utf-8').strip() + '_model'].astype(dtype='object')
        argo_n[sensors[i].decode('utf-8').strip() + '_model'] = argo_n[sensors[i].decode('utf-8').strip() + '_model'].assign_attrs(long_name=f'Sensor model from the "SENSOR_MODEL" field in the meta.nc file')
    # load calibration information from Sprof file 
    n_prof = argo_n.sizes['N_PROF']
    # gets order of sensors to extract calibration comments
    if verbose: 
        if 'DOXY' in argo_n.keys():
            print('Oxygen is present')
        else:
            print('Oxygen is NOT present')
    if 'DOXY' in argo_n.keys(): # only look for oxygen calibration comments if oxygen sensor is present
        # some profiles might be missing sensor name (not sure why) so loop through looking
        o2_ind_all = np.full(n_prof, np.nan)
        # finds where o2 calibration comment is in each profile (at least for one float it changes from the first profile to the rest)
        for p in range(0, n_prof):
            cal_str = argo_n.STATION_PARAMETERS.values.astype(str)[p]
            # print(cal_str)
            for i, param in enumerate(cal_str):
                if 'DOXY' in param:
                    o2_ind_all[p] = i
                    break
        if verbose: 
            print(f'o2_ind_all: {o2_ind_all}')
            print(np.sum(~np.isnan(o2_ind_all)))

        o2_cal_full = []
        for idx, o2_ind in enumerate(o2_ind_all):
            if ~np.isnan(o2_ind):
                o2_cal_full.append(argo_n.SCIENTIFIC_CALIB_COMMENT.values[idx,-1,np.int32(o2_ind)])

        data_comment_eq = pd.DataFrame({'o2_cal_full': o2_cal_full})
        # Drop duplicate comments
        unique_comments = data_comment_eq.drop_duplicates()
        # Create an empty list to store DataFrames
        data_frames = []
        # print(o2_cal_unique)
        # might have multiple unique comments, so save out each one:
        for i in unique_comments.index:
            o2_cal_i = unique_comments.o2_cal_full[i].decode("utf-8")

            # print(o2_cal_i)
            # save with wmo only for now
            new_data_cal_info = pd.DataFrame({ 'o2_cal_comment': [o2_cal_i]})
            # Append the new_data_cal_info DataFrame to the list
            data_frames.append(new_data_cal_info)
            if verbose: print(new_data_cal_info)

        cal_type = 'not_air'
        # assume if there are any air calibration comments for a float then I can label the entire float as "air_cal"
        for d in range(len(data_frames)):
            for air_str in air_cal_list:
                air_test = (data_frames[d].o2_cal_comment.str.contains(air_str))
                for conf_str in confusing_air_cal_list:
                    confusing_test = (data_frames[d].o2_cal_comment.str.contains(conf_str))
                if air_test[0] and not confusing_test[0]:
                    cal_type = 'air_cal'
                    break
                
        argo_n['O2_cal_type'] = cal_type
        argo_n['O2_cal_type'] = argo_n['O2_cal_type'].astype(dtype='object')
        argo_n['O2_cal_type'] = argo_n['O2_cal_type'].assign_attrs(long_name='Type of oxygen calibration performed on the float, from the "SCIENTIFIC_CALIB_COMMENT" field in the Sprof file. "air_cal" = float was calibrated in air, "not_air" = float was not calibrated in air')

    if verbose: 
        print('Looking up ocean basin info')
        print(np.sum(argo_index['WMO_ID']==str(argo_n["WMO_ID"].values)))
    if np.sum(argo_index['WMO_ID']==str(argo_n["WMO_ID"].values))==0:
        if verbose: print(f'WMO_ID {argo_n["WMO_ID"].values.item()} not found in argo_index')
        argo_n['ocean'] = 'not found'
        argo_n['profiler_type'] = 'not found'
    else:
        # try to get ocean basin information from the synthetic argo profile file:
        profile_ocean = argo_index.loc[argo_index['WMO_ID']==str(argo_n["WMO_ID"].values), 'ocean']
        ocean_list = []
        # check if one value is nan and remove if so
        for oc in profile_ocean.unique():
            if type(oc)==str:
                ocean_list.append(oc)
                    
        argo_n['ocean'] = ', '.join(ocean_list)
        argo_n['ocean'] = argo_n['ocean'].astype(dtype='object')

        if verbose: print('Looking up profiler type')
        # profiler_type
        profiler_type = argo_index.loc[argo_index['WMO_ID']==str(argo_n["WMO_ID"].values), 'profiler_type']
        if verbose: print(profiler_type)
        argo_n['profiler_type'] = profiler_type.iloc[0]
        argo_n['profiler_type'] = argo_n['profiler_type'].astype(dtype='object')
           
    argo_n['ocean'] = argo_n['ocean'].assign_attrs(long_name='Ocean basin(s) where the float has been deployed, from the synthetic Argo profile index file')
    argo_n['profiler_type'] = argo_n['profiler_type'].assign_attrs(long_name='Type of profiler used by the float, from the synthetic Argo profile index file. Table of Argo profiler types can be found at: https://vocab.nerc.ac.uk/collection/R08/current/')
    # add O2 sat
    return argo_n 

def apply_flags(argo_n, flags_to_remove, verbose):
    plot_flags = False

    #Load sprof file 
    # argo_n = xr.open_dataset(sprof_path + file)
    
    
    pres_data = argo_n['PRES_ADJUSTED'].values
    nprof_n = argo_n.sizes['N_PROF']

    if verbose:
        print('setting flagged data to nan')
    for key in argo_n.keys():
        # if verbose: print(key)
        if key.endswith('_QC'):
            if key.startswith('PROFILE_') or key=='JULD_QC':  #or key=='POSITION_QC'
                continue

            # if key=='POSITION_QC':
            #     argo_n[key[:-3] + '_BGCArgoPlus_flag'] = ''
            #     qc_val = argo_n[key].values.astype('float')
            #     if verbose:
            #         print('In flag removal loop ' + key[:-3])
            #     for flag in flags_to_remove:
            #     # argo_n[key[:-3]].values[qc_val==flag]=np.nan
            #     argo_n[key[:-3] + '_BGCArgoPlus'].values[qc_val==flag]=np.nan
            #     # print(argo_n[key[:-3] + '_BGCArgoPlus_flag'].values)

            #     argo_n[key + '_n_' + str(flag) + '_removed'] = np.sum(qc_val==flag) # save out the number of data points removed for a given flag 
            # argo_n[key[:-3] + '_BGCArgoPlus_flag'] = argo_n[key[:-3] + '_BGCArgoPlus_flag'].values.item() + 'F_'


            # need to check position QC separately from the rest and apply to lat/lon
            if np.logical_or(key=='POSITION_QC', key.__contains__('ADJUSTED')==True): # only apply flag filtering to ADJUSTED variables or Position
                # if verbose: print('here')

                # make a copy of variable, w/ "BGCArgoPlus" appended - Keeps original "Var_ADJUSTED" unchanged
                if key!='POSITION_QC':
                    argo_n[key[:-3] + '_BGCArgoPlus'] = argo_n[key[:-3]].copy()
                    argo_n[key[:-3] + '_BGCArgoPlus'] = argo_n[key[:-3] + '_BGCArgoPlus'].assign_attrs(comment='This variable is originally copied from ' + key[:-3] + ' with additional QC and outlier detection applied. See Bushinsky et al. (2026) ESSD for more details')
                # create an empty string variable to contain information on what processing has been performed on different variables
                argo_n[key[:-3] + '_BGCArgoPlus_flag'] = ''
                # if verbose: print('here')
                # argo_n[key[:-3] + '_RO'] = argo_n[key[:-3]].copy()

                # print(key) if key.endswith('_QC') else None# 
                qc_val = argo_n[key].values.astype('float')
                if verbose:
                    print('In flag removal loop ' + key[:-3])
                # if key.startswith('PH'):
                #     break
                # sets data to nan (currently both raw and ADJUSTED, maybe can skip)
                flag_removed_str = ''
                for flag in flags_to_remove:
                    if verbose: print(f'removing flag {flag} from {key[:-3]}')
                    flag_removed_str = flag_removed_str + str(flag) + '_'
                    if verbose:
                        print('removing flag ' + str(flag) + ' from ' + key[:-3] + ' (' + str(np.sum(qc_val==flag)) + ' values set to nan)')
                    if key=='POSITION_QC':
                        argo_n['LATITUDE'].values[qc_val==flag]=np.nan
                        argo_n['LONGITUDE'].values[qc_val==flag]=np.nan
                    else:
                        argo_n[key[:-3] + '_BGCArgoPlus'].values[qc_val==flag]=np.nan
                    # print(argo_n[key[:-3] + '_BGCArgoPlus_flag'].values)

                    argo_n[key + '_n_' + str(flag) + '_removed'] = np.sum(qc_val==flag) # save out the number of data points removed for a given flag 
                    argo_n[key + '_n_' + str(flag) + '_removed'] = argo_n[key + '_n_' + str(flag) + '_removed'].assign_attrs(long_name=f'Number of data points removed from {key[:-3]} due to Argo QC flag {flag} in {key}')
                    argo_n[key[:-3] + '_BGCArgoPlus_flag'] = argo_n[key[:-3] + '_BGCArgoPlus_flag'].values.item() + 'F_'
                    argo_n[key[:-3] + '_BGCArgoPlus_flag'] = argo_n[key[:-3] + '_BGCArgoPlus_flag'].assign_attrs(long_name=f'A string listing each of the processing steps that have been applied to the variable. F = data with Argo QC flag(s) {flag_removed_str} set to NaN, S = surface data set to NaN, Inv = density inversions set to NaN, DMonly = non-delayed mode data removed')
            else:
                continue

            if plot_flags==1:
                qc_var_data = argo_n[key[:-3]].values
                plot_flag_filtering(argo_n, key, pres_data, qc_val, qc_var_data, nprof_n)
    
    #Finding and removing all non-delayed mode data
    # sometimes parameters are missing from profiles - 
    # need to loop through all profiles and check which parameters are present
    # assumes that mode applies to all levels of a profile
    parameter_array = argo_n.STATION_PARAMETERS.values.astype(str)
    all_parameters = []

    if verbose:
        print('checking profile data modes')
    for idx in range(len(parameter_array)): # loop through all profiles, idx is profile index 
        prof_parameters = parameter_array[idx] # get parameters present for each profile 
        # if verbose:
        #     print(idx)

        # loop through each paramter in the profile 
        for var in prof_parameters:
            # if verbose:
            #     print(var)
            var_str = var.strip()
            if var_str=='':
                continue
            all_parameters.append(var_str)
            # if verbose: print(var_str + '_profile_mode')
            if var_str + '_profile_mode' not in argo_n: # if profile removed count has not been initialized, do so here 
                # if verbose: print('initializing profile removal count')
                temp_var = {var_str + '_profile_mode':(['N_PROF'], np.zeros(len(argo_n['N_PROF'])))}
                argo_n = argo_n.assign(temp_var)
                argo_n[var_str + '_profile_mode'][:] = np.nan
                argo_n[var_str + '_profile_mode'] = argo_n[var_str + '_profile_mode'].assign_attrs(long_name='Numerical representation of the profile mode for each variable. 0 = Real Time, 1 = Real Time Adjusted, 2 = Delayed Mode. From the "PARAMETER_DATA_MODE" or "DATA_MODE" field in the Sprof file.',
                    missing_value=np.nan)
            # if np.logical_or(len(var_str)==0, var_str=='PRES'): # only proceed if the variable exists and is not pressure
            #     continue
            # if np.logical_or(var_str=='TEMP', var_str=='PSAL'): # don't set TEMP or PSAL to nan for non-delayed mode data
            #     continue
            # if var_str + '_profile_removed_not_D' not in argo_n: # if profile removed count has not been initialized, do so here 
            #     argo_n[var_str + '_profile_removed_not_D'] = 0
            # if verbose: 
            #     print(prof_parameters)
            #     print(var_str)

            if 'PARAMETER_DATA_MODE' in argo_n.keys():
                var_ind = [p_idx for p_idx, s in enumerate(prof_parameters) if s.strip()== var_str]

                data_mode_variable_name = 'PARAMETER_DATA_MODE'
                var_data_mode = argo_n[data_mode_variable_name][idx,var_ind].values
            else:
                # if verbose: print('Core Argo float, using "DATA MODE" instead of "PARAMETER_DATA_MODE"')
                data_mode_variable_name = 'DATA_MODE'
                var_data_mode = argo_n[data_mode_variable_name][idx].values

            # if verbose:
            #     print(var_ind)
            #     print(argo_n)
            #     print(argo_n[data_mode_variable_name])
            # get parameter data mode values for that profile / variable
            # if verbose: print(var_data_mode)
            # print(var_str)
            decoded_arr = np.array([elem.decode() if isinstance(elem, bytes) else np.nan for elem in var_data_mode.flatten()])
            # if verbose:
                # print(decoded_arr)
            if decoded_arr=='R':
                # if verbose: print('Real time')
                argo_n[var_str + '_profile_mode'][idx] = 0
            elif  decoded_arr=='A':
                # if verbose: print('Real Time Adjusted')
                argo_n[var_str + '_profile_mode'][idx] = 1
            elif  decoded_arr=='D':
                # if verbose: print('Delayed Mode')
                argo_n[var_str + '_profile_mode'][idx] = 2

            # # print(decoded_arr)
            # result = np.where(decoded_arr == 'D', False, True) # true whereever mode is not delayed
            # # print(result)
            # if result:
            #     argo_n[var_str +'_ADJUSTED_RO'][idx,:] = np.nan
            #     argo_n[var_str + '_profile_removed_not_D']+=1 # count the number of profiles removed b/c mode does not equal "D"
    if verbose:
        print('setting non-delayed mode profiles to nans')
    # set non delayed mode data to nan for now
    for var in np.unique(all_parameters):
        if var in {'NITRATE', 'DOXY', 'PH_IN_SITU_TOTAL'}:
            argo_n[var + '_ADJUSTED_BGCArgoPlus'][argo_n[var + '_profile_mode']!=2]=np.nan
            # argo_n[var + '_ADJUSTED_BGCArgoPlus_flag'] = argo_n[var + '_ADJUSTED_BGCArgoPlus_flag'].values.item() + 'DMonly_'
            argo_n = append_flag(argo_n, var + '_ADJUSTED_BGCArgoPlus_flag', 'DMonly_')
    return argo_n
def surface_removal(argo_n, surface_cutoff, outlier_df, verbose):
    
    pres_var = 'PRES' # using the least QCd pressure in case Pressure isn't marked as good but there is valid BGC that that needs to be removed... 
    wmo = str(argo_n['WMO_ID'].values)

    n_prof = len(argo_n.N_PROF)
    plot_vars_all=[ 'TEMP_ADJUSTED_BGCArgoPlus',
                'PSAL_ADJUSTED_BGCArgoPlus',
                'NITRATE_ADJUSTED_BGCArgoPlus',
                'DOXY_ADJUSTED_BGCArgoPlus',
                'PH_IN_SITU_TOTAL_ADJUSTED_BGCArgoPlus',]
    for p in range(n_prof):
        # find the presence of samples shallower than surface cutoff 
        surf_index = argo_n[pres_var][p,:]<surface_cutoff

        if any(surf_index): # if there are any true values
            # cycle through variables:
            var_removal_count = 0

            for var in plot_vars_all:
                if var in argo_n: # only do this for variables that are present
                    # first count number of non-nan instances of the variable in the surface
                    var_surf_index = np.logical_and(surf_index, ~np.isnan(argo_n[var][p,:]))

                    if any(var_surf_index): # if valid data needs to be removed
                        # print(var)
                        # add count to existing count for this variable:
                        if var + '_surface_samples_removed' not in argo_n: # if count has not been initialized, do so here 
                            argo_n[var + '_surface_samples_removed'] = np.sum(var_surf_index).values
                        else:
                            argo_n[var + '_surface_samples_removed'] = argo_n[var + '_surface_samples_removed'].values + np.sum(var_surf_index).values
                        argo_n[var + '_surface_samples_removed'] = argo_n[var + '_surface_samples_removed'].assign_attrs(long_name=f'Number of surface samples removed from {var} due to being shallower than {surface_cutoff} dbar')

                        juld_val = pd.to_datetime(argo_n['JULD'][p].values)
                        date_str = juld_val.strftime('%Y-%m-%d %H:%M:%S') if not pd.isna(juld_val) else 'NaT'

                        if verbose: 
                            print([date_str]*np.sum(var_surf_index).values.item())
                        outlier_df_n = pd.DataFrame({                 # creates a new "outlier_df" w wmo, variable name, profile, date, level, pressure, deletion reason 
                            "Float Number":[wmo]*np.sum(var_surf_index).values.item(),
                            "Variable":[var]*np.sum(var_surf_index).values.item(),
                            "N_PROF":[p]*np.sum(var_surf_index).values.item(),
                            "Date":[date_str]*np.sum(var_surf_index).values.item(),
                            "N_LEVELS":argo_n['N_LEVELS'][var_surf_index].values,
                            'Pressure (dbar)':argo_n['PRES'][p,var_surf_index].values,
                            'Deletion reason':['Surface Removal']*np.sum(var_surf_index).values.item()
                        })   
                        if len(outlier_df)==0:
                            outlier_df = outlier_df_n
                        else:
                            outlier_df = pd.concat([outlier_df, outlier_df_n])

                        # set surface obs to nan:
                        argo_n[var][p,var_surf_index] = np.nan
    for var in plot_vars_all:
        if var in argo_n: # only do this for variables that are present
            # argo_n[var[:-12] + '_BGCArgoPlus_flag'] = argo_n[var[:-12] + '_BGCArgoPlus_flag'].values.item() + 'S_'
            argo_n = append_flag(argo_n, var[:-12] + '_BGCArgoPlus_flag', 'S_')

    # print('Finishing Surface Removal')
    return argo_n, outlier_df

def sigma0(salinity,temperature,lon,lat,pressure):
    SA = gsw.SA_from_SP(salinity,
                        pressure,
                        lon,
                        lat)

    CT = gsw.CT_from_t(SA,
                       temperature,
                       pressure)

    sigma = gsw.sigma0(SA,CT)
    
    return sigma

def density_inversion(argo_n, processed_fig_dir, outlier_df, verbose):
    dens_thresh = -0.025 # acceptable density decrease

    # which temperature, salinity, and pressure are being checked? 
    temp_var = 'TEMP_ADJUSTED_BGCArgoPlus'
    sal_var = 'PSAL_ADJUSTED_BGCArgoPlus'
    pres_var = 'PRES_ADJUSTED_BGCArgoPlus'

    # argo_n = copy.deepcopy(argo_n)
    argo_n = argo_n.set_coords((pres_var, 'LONGITUDE', 'LATITUDE', 'JULD'))
    argo_n['Sigma_theta_gsw'] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_var].shape))
    argo_n.Sigma_theta_gsw[:] = np.nan
    argo_n['Sigma_theta_gsw'] = argo_n['Sigma_theta_gsw'].assign_attrs(long_name='Temporary density calculated for the purpose of determining density inversions. Potential density anomaly (sigma-theta) calculated using the Gibbs SeaWater (GSW) Oceanographic Toolbox for TEOS-10.')
    n_prof = len(argo_n.N_PROF)
    wmo = str(argo_n['WMO_ID'].values)

    params_to_del = [sal_var]

    for par in params_to_del:
        # argo_n[par + '_flag'] = argo_n[par + '_flag'].values.item() + 'Inv_'
        argo_n = append_flag(argo_n, par + '_flag', 'Inv_')
    # create copies of starting T and S to use for figures if a density inversion is found
    sal_orig = argo_n[sal_var].copy()
    temp_orig = argo_n[temp_var].copy()

    bad_data = 0
    removed_profile_record = []
    removed_pressure_record = []
    removed_var_record = []

    if verbose:
        print('Searching each profile for density inversions')
    for p in range(n_prof):
        # if verbose:
            # print('Starting profile ' + str(p))
        exit=0
        while exit==0:
            # count = 0
            
            pres_p = argo_n[pres_var][p,:].sortby(pres_var).values.copy()
            sal = argo_n[sal_var][p,:].sortby(pres_var).values
            temp = argo_n[temp_var][p,:].sortby(pres_var).values

            sigma_t_p = sigma0(sal,
                                        temp,
                                        argo_n.LONGITUDE[p].values,
                                        argo_n.LATITUDE[p].values,
                                        pres_p)
            
            sigma_t_p_nonan = sigma_t_p[~np.isnan(sigma_t_p)]
            pres_p_nonan = pres_p[~np.isnan(sigma_t_p)]
            sal_p_nonan = sal[~np.isnan(sigma_t_p)]

            # Check if there is a density profile
            if len(sigma_t_p_nonan) < 2:
                exit=1
            else:
                # Find most negative relative ([x+1]-[x]) and absolute (x - x[0]) difference
                min_diff_rel = min(np.diff(sigma_t_p_nonan))
                min_diff_tot = min(sigma_t_p_nonan-sigma_t_p_nonan[0])
                abs_min = min(min_diff_rel,min_diff_tot) 
                if abs_min < dens_thresh:  # if something trips the threshold, remove, otherwise go to the next profile 
                    #print(f"bad density, profile no. {p}")

                    bad_data += 1
                    if min_diff_rel < min_diff_tot:
                        outlier_index = np.argmin(np.diff(sigma_t_p_nonan))
                    else:
                        outlier_index = np.argmin(sigma_t_p_nonan-sigma_t_p_nonan[0])                
                        
                    window_size = 3 
                    df = pd.DataFrame({'Pressure':pres_p_nonan, 'Density':sigma_t_p_nonan})
                    df['Density_rollmean'] = df['Density'].rolling(window=window_size, center=True).mean()
                    if outlier_index+1==len(df):
                        exit=1 
                        continue
                    if pres_p_nonan[outlier_index].item()>600: # don't worry about density inversions deeper than 600db
                        exit=1
                        continue

                    diff1 = abs(df.loc[outlier_index, 'Density'] - df.loc[outlier_index, 'Density_rollmean'])
                    diff2 = abs(df.loc[outlier_index+1, 'Density'] - df.loc[outlier_index+1, 'Density_rollmean'])

                    if diff1 > diff2:   
                        outlier_pres = [pres_p_nonan[outlier_index]]
                        outlier_sal = [sal_p_nonan[outlier_index]]

                    elif diff1 < diff2:
                        outlier_pres = [pres_p_nonan[outlier_index+1]]
                        outlier_sal = [sal_p_nonan[outlier_index+1]]

                    else:
                        outlier_pres = [ pres_p_nonan[outlier_index],pres_p_nonan[outlier_index+1] ]
                        outlier_sal = [ sal_p_nonan[outlier_index],sal_p_nonan[outlier_index+1] ]
                    if verbose:
                        print(f'wmo {wmo}, profile {p}, outlier_pres {outlier_pres[0].item()}, salinity {outlier_sal[0].item()}')
                    for idx, i in enumerate(outlier_pres):
                        for par in params_to_del:
                            try:
                                mask = argo_n[pres_var].isel(N_PROF=p).data == i
                                # print(mask)
                                argo_n[par].data[p,mask] = np.nan
                                removed_profile_record.append(p)
                                removed_pressure_record.append(i)
                                removed_var_record.append(outlier_sal[idx])
                                outlier_df_n = pd.DataFrame({                 # creates a new "outlier_df" w wmo, variable name, profile, date, level, pressure, deletion reason 
                                                "Float Number":[wmo],
                                                "Variable":[par],
                                                "N_PROF":[p],
                                                "Date":[pd.to_datetime(argo_n['JULD'][p].values).strftime('%Y-%m-%d %H:%M:%S')],
                                                "N_LEVELS":[argo_n['N_LEVELS'][mask].values.item()],
                                                'Pressure (dbar)':[i],
                                                'Deletion reason':['Density Inversion']
                                            })   
                                if len(outlier_df)==0:
                                    outlier_df = outlier_df_n
                                else:
                                    outlier_df = pd.concat([outlier_df, outlier_df_n])  
                                # argo_n[par][p,:][argo_n[pres_var][p,:]==i] = np.nan
                            except:
                                pass # if parameter not present in argo_n
                else:
                    exit=1

    # save plots of bad profiles if they exist
    if bad_data>0:
        def autoscale_plot(x, y, x_min, x_max, y_margin=0.1):
            mask = (x >= x_min) & (x <= x_max)
            y_data = y[mask]
            if np.isnan(y_data).all():
                raise ValueError("No valid data within specified x-range.")
            y_min = np.nanmin(y_data)
            y_max = np.nanmax(y_data)
            y_range = y_max - y_min
            y_lims = [y_min - y_range * y_margin, y_max + y_range * y_margin]
            return y_lims

        # check if a Processed_Figures/[WMO] directory exists:
        
        float_plot_dir = os.path.join(processed_fig_dir, str(wmo))
        if not os.path.isdir(float_plot_dir):
            os.mkdir(float_plot_dir)
        
        removed_profile_record = np.array(removed_profile_record)
        removed_pressure_record = np.array(removed_pressure_record)
        removed_var_record = np.array(removed_var_record)

        unique_profiles = np.unique(removed_profile_record)
        for p in unique_profiles:
            if verbose:
                print('Plotting and saving profile ' + str(p))
            plot_filename = '999_' + str(argo_n['WMO_ID'].values) + '_Density_Inversions_Removed_Profile_' + str(p)

            selected_profile_index = removed_profile_record==p
            profile_pressures_removed = removed_pressure_record[selected_profile_index]
            profile_sal_removed = removed_var_record[selected_profile_index]

            # min_pres = outlier_pres[0] - 10
            # max_pres = outlier_pres[0] + 10

            shallow_pressure = np.min(profile_pressures_removed)-20
            deep_pressure = np.max(profile_pressures_removed)+20
            fig = plt.figure(figsize=(14,12))
            ax = fig.add_subplot(2,3,4) 
            ax.plot(temp_orig[p,:], argo_n[pres_var][p,:],'x-', label='Orig')
            ax.set_autoscalex_on

            axb = fig.add_subplot(2,3,1) # upper 500 m, not just zoomed in
            axb.plot(temp_orig[p,:], argo_n[pres_var][p,:],'x-', label='Orig')
            axb.set_autoscalex_on

            ax1 = fig.add_subplot(2,3,5)
            ax1.plot(sal_orig[p,:], argo_n[pres_var][p,:],'-x', label='Orig')

            ax1b = fig.add_subplot(2,3,2)
            ax1b.plot(sal_orig[p,:], argo_n[pres_var][p,:],'-x', label='Orig')


            for idx, i in enumerate(profile_pressures_removed):
                # ax1.plot(sal_orig[p,argo_n[pres_var][p,:]==i], i, 'ro')
                ax1.plot(profile_sal_removed[idx], i, 'ro')
                ax1b.plot(profile_sal_removed[idx], i, 'ro')




            ax2 = fig.add_subplot(2,3,6)
            ax2b = fig.add_subplot(2,3,3)

            orig_sigma_t_p = sigma0(sal_orig[p,:],
                                        temp_orig[p,:],
                                        argo_n.LONGITUDE[p].values,
                                        argo_n.LATITUDE[p].values,
                                        argo_n[pres_var][p,:])
            new_sigma_t_p = sigma0(argo_n[sal_var][p,:],
                                        argo_n[temp_var][p,:],
                                        argo_n.LONGITUDE[p].values,
                                        argo_n.LATITUDE[p].values,
                                        argo_n[pres_var][p,:])

            # 
            ax2.plot(orig_sigma_t_p, argo_n[pres_var][p,:],'-')
            ax2.plot(new_sigma_t_p, argo_n[pres_var][p,:],'m--')

            ax2b.plot(orig_sigma_t_p, argo_n[pres_var][p,:],'-')
            ax2b.plot(new_sigma_t_p, argo_n[pres_var][p,:],'m--', marker='^')


            ax.set_ylim(deep_pressure, shallow_pressure)
            ax1.set_ylim(deep_pressure, shallow_pressure)
            ax2.set_ylim(deep_pressure, shallow_pressure)

            zoom_out = [deep_pressure+250, np.max([shallow_pressure-250,0]).item()]
            axb.set_ylim(zoom_out[0], zoom_out[1])
            ax1b.set_ylim(zoom_out[0], zoom_out[1])
            ax2b.set_ylim(zoom_out[0], zoom_out[1])

            auto_x_lim = autoscale_plot(argo_n[pres_var][p,:].values,  temp_orig[p,:].values,   shallow_pressure, deep_pressure)
            if auto_x_lim[0]!=auto_x_lim[1]:
                ax.set_xlim(auto_x_lim)
            auto_x_lim = autoscale_plot(argo_n[pres_var][p,:].values,  temp_orig[p,:].values,  zoom_out[1], zoom_out[0])
            if auto_x_lim[0]!=auto_x_lim[1]:
                axb.set_xlim(auto_x_lim)

            auto_x_lim = autoscale_plot(argo_n[pres_var][p,:].values,  sal_orig[p,:].values,   shallow_pressure, deep_pressure)
            if auto_x_lim[0]!=auto_x_lim[1]:
                ax1.set_xlim(auto_x_lim)
            auto_x_lim = autoscale_plot(argo_n[pres_var][p,:].values,  sal_orig[p,:].values,  zoom_out[1], zoom_out[0])
            if auto_x_lim[0]!=auto_x_lim[1]:
                ax1b.set_xlim(auto_x_lim)

            auto_x_lim = autoscale_plot(argo_n[pres_var][p,:].values,  orig_sigma_t_p,   shallow_pressure, deep_pressure)
            if auto_x_lim[0]!=auto_x_lim[1]:
                ax2.set_xlim(auto_x_lim)
            auto_x_lim = autoscale_plot(argo_n[pres_var][p,:].values,  orig_sigma_t_p,  zoom_out[1], zoom_out[0])
            if auto_x_lim[0]!=auto_x_lim[1]:
                ax2b.set_xlim(auto_x_lim)

            ax.set_ylabel('Pressure')
            ax.set_xlabel('T')
            ax1.set_xlabel('S')
            ax2.set_xlabel('PDENS')

            axb.set_ylabel('Pressure')
            axb.set_xlabel('T')
            ax1b.set_xlabel('S')
            ax2b.set_xlabel('PDENS')

            axb.set_title(str(argo_n['WMO_ID'].values) + ' Profile: ' + str(p))
            fig.tight_layout()
            if verbose: print(f'saving figure to: {float_plot_dir}/{plot_filename}.png') 
            plt.savefig(f'{float_plot_dir}/{plot_filename}.png')
            plt.close()
            plt.clf()
    if verbose:
        print('Finishing Density Inversion')
    return argo_n, outlier_df

def plot_flag_filtering(argo_n, test_key, pres_data, qc_val, qc_var_data, nprof_n):
    time_2D = np.tile(argo_n['JULD'].values,(len(argo_n['N_LEVELS']),1)).T

    # plot all profiles and QC data
    var_data = argo_n[test_key[:-3]].values

    # Function for repetitive parts of plot
    def config_depth_section(cbar_label,title):
        plt.gca().invert_yaxis()
        plt.ylabel('Pressure (dbar)')
        plt.colorbar(label=cbar_label)
        plt.title(title)

    plt.figure(figsize=(16,6))
    plt.subplot(1,3,1)

    for p in range(nprof_n):    
        plt.scatter(time_2D[p,:],  pres_data[0,:],c=var_data[0,:],s=1)

    config_depth_section(argo_n[test_key[:-3]].units, test_key[:-3])

    plt.subplot(1,3,2)
    for p in range(nprof_n):    
        plt.scatter(time_2D[p,:], pres_data[0,:], c=qc_val[0,:],cmap='Set1',vmin=0.5,vmax=9.5)

    config_depth_section(argo_n[test_key].long_name, test_key)

    # remove data according to list of flags, plot remaining data:
    flags_to_remove = [3,4]

    plt.subplot(1,3,3)
    for p in range(nprof_n):    
        plt.scatter(time_2D[p,:],  pres_data[0,:],c=qc_var_data[0,:],s=1)

    config_depth_section(argo_n[test_key].long_name, test_key[:-3] + ' filtered')
    # plt.pcolor(time_2D,data['PRES_ADJUSTED'].values,data['NITRATE_QC'].values.astype(float),cmap='Set1',vmin=0.5,vmax=9.5)
    # config_depth_section('QC flag',"Parameter: 'NITRATE_QC'")
    plt.tight_layout()


def bottom_oxygen_check(argo_n, group_outlier_dir, outlier_df, verbose):
    deepest_meas = argo_n['PRES_ADJUSTED_BGCArgoPlus'].max().values.item()
    min_depth = deepest_meas-300 # change to deepest depth minus 300m 
    # argo_n = xr.open_dataset(output_dir + argo_file)
    nprof_n = len(argo_n['N_PROF'])
    deriv = {}
    o2_stored = {}

    if verbose:
        print('Getting NaN oxygen and pressure values, calculating derivative of oxygen data')
    for p in range(nprof_n):
        if len(argo_n['PRES_ADJUSTED_BGCArgoPlus'][p,argo_n['PRES_ADJUSTED_BGCArgoPlus'][p,:]>=min_depth])==0:
            continue
        # get data below min_depth db - this corresponds to the bottom 300m of the profile
        pres = argo_n['PRES_ADJUSTED_BGCArgoPlus'][p].where(argo_n['PRES_ADJUSTED_BGCArgoPlus'][p,:]>=min_depth,drop=True)
        o2 = argo_n['DOXY_ADJUSTED_BGCArgoPlus'][p].where(argo_n['PRES_ADJUSTED_BGCArgoPlus'][p,:]>=min_depth,drop=True)

        pres_nonan = pres.where(~np.isnan(pres) & ~np.isnan(o2), drop=True)
        o2_nonan = o2.where(~np.isnan(pres) & ~np.isnan(o2), drop=True)
    
        if len(o2_nonan)<15:
            continue
        # calculate the derivative of the bottom 300 db oxygen
        deriv[p] = o2_nonan.differentiate(coord='PRES_ADJUSTED_BGCArgoPlus')
        o2_stored[p] = o2_nonan
    # if verbose:
    #     print('Here')
    if len(deriv)==0:
        return argo_n, outlier_df
    
    if 'deriv_df' in locals():
        del(deriv_df)
    
    # assemble all derivatives here
    for p in deriv.keys():
        pres_out = deriv[p]['PRES_ADJUSTED_BGCArgoPlus'].values
        o2 = o2_stored[p].values

        o2_deriv_out = deriv[p].values
        # o2_deriv_out


        deriv_df_p = pd.DataFrame({            
                "p":len(pres_out)*[p],
                "pres":pres_out,
                "o2_conc":o2,
                "o2_deriv":o2_deriv_out,
            })    
        if 'deriv_df' in locals():
        
            deriv_df = pd.concat([deriv_df, deriv_df_p])
        else:
            deriv_df = deriv_df_p
    # if verbose:
    #     print('Here2')
    deriv_df = deriv_df.reset_index()

    wmo = str(argo_n['WMO_ID'].values)
    # wmo = argo_file[:7]

    float_RO_profile_dir = group_outlier_dir + '/' + wmo + '/'
    
    if verbose:
        print('Checking for extreme changes in dO2/dz at the bottom of each profile')
    for p in deriv_df['p'].unique(): # :#  range(47,50): #
        profile_deriv =  deriv_df['o2_deriv'].where(deriv_df['p']==p)
        if np.sum(~np.isnan(profile_deriv))==0: # move to next profile if no valid data
            continue
        profile_pres = deriv_df['pres'].where(deriv_df['p']==p)

        last_index = profile_deriv.last_valid_index()

        # calculate the mean and standard without considering that last few samples. 
        # Sometimes the large values at the end drive the mean and std so high that there is a significant bias

        deriv_mean = profile_deriv[:last_index-2].mean()
        deriv_std = profile_deriv[:last_index-2].std()

        # check whether the last three O2 derivatives are 3 times greater than the SD of the deep derivatives. If it is less, continue to the next profile 
        last_deriv_check = np.abs(profile_deriv[last_index] - deriv_mean)>deriv_std*3
        second_to_last_deriv_check = np.abs(profile_deriv[last_index-1] - deriv_mean)>deriv_std*3
        third_to_last_deriv_check = np.abs(profile_deriv[last_index-2] - deriv_mean)>deriv_std*3

        
        # check whether there is a large pressure gap between the last two valid points:
        last_pressure = profile_pres[last_index]
        second_to_last_pressure = profile_pres[last_index-1]
        third_to_last_pressure = profile_pres[last_index-2]

        last_pressure_gap = last_pressure - second_to_last_pressure
        second_to_last_pressure_gap = second_to_last_pressure - third_to_last_pressure
        mean_pressure_gap = profile_pres[:last_index-3].diff().mean()


        if third_to_last_deriv_check:
            outlier_index = [last_index-2, last_index-1, last_index]
            if verbose:
                print(f'Third to last deriv check: {third_to_last_deriv_check}, diff between third to last deriv and mean: {np.abs(profile_deriv[last_index-2] - deriv_mean)}, 3*std: {deriv_std*3}')
                        
        elif second_to_last_deriv_check:
            outlier_index = [last_index-1, last_index]
            if verbose:
                print(f'Second to last deriv check: {second_to_last_deriv_check}, diff between second to last deriv and mean: {np.abs(profile_deriv[last_index-1] - deriv_mean)}, 3*std: {deriv_std*3}')
                   
        elif last_deriv_check:
            outlier_index = last_index
            if verbose:
                print(f'Last deriv check: {last_deriv_check}, diff between last deriv and mean: {np.abs(profile_deriv[last_index] - deriv_mean)}, 3*std: {deriv_std*3}')
        elif last_pressure_gap>mean_pressure_gap*2:
            if second_to_last_pressure_gap>mean_pressure_gap*2:
                outlier_index = [last_index-1, last_index]
            else:
                outlier_index = last_index
        else:
            outlier_index = []

        if np.any(outlier_index):
            if not os.path.isdir(float_RO_profile_dir):
                os.mkdir(float_RO_profile_dir)

            plot_filename = wmo + '_bottom_oxygen_profile_' +  str(p)

            pres = argo_n['PRES_ADJUSTED_BGCArgoPlus'][p,argo_n['PRES_ADJUSTED_BGCArgoPlus'][p,:]>=min_depth]
            o2 = argo_n['DOXY_ADJUSTED_BGCArgoPlus'][p,argo_n['PRES_ADJUSTED_BGCArgoPlus'][p,:]>=min_depth]
            temp = argo_n['TEMP_ADJUSTED_BGCArgoPlus'][p,argo_n['PRES_ADJUSTED_BGCArgoPlus'][p,:]>=min_depth]
            sal = argo_n['PSAL_ADJUSTED_BGCArgoPlus'][p,argo_n['PRES_ADJUSTED_BGCArgoPlus'][p,:]>=min_depth]

            pres_nonan = pres.where(~np.isnan(pres) & ~np.isnan(o2), drop=True)
            o2_nonan = o2.where(~np.isnan(pres) & ~np.isnan(o2), drop=True)
            temp_nonan = temp.where(~np.isnan(pres) & ~np.isnan(o2), drop=True)
            sal_nonan = sal.where(~np.isnan(pres) & ~np.isnan(o2), drop=True)

            max_pres = max(pres_nonan) + 10
            fig = plt.figure(figsize=(16,6))

            ax = fig.add_subplot(2,4,1)
            ax.plot(temp_nonan, pres_nonan, 'k-', marker='.')
            ax.set_ylim([max_pres, min_depth])
            ax.set_title(f'a. WMO {wmo}, Profile {p}', loc='left')
            ax.set_ylabel('Pressure (dbar)')
            temp_units = argo_n["TEMP_ADJUSTED_BGCArgoPlus"].attrs.get('units', ' ')
            # replace degree_Celcius w/ degree symbol
            temp_units = temp_units.replace('degree_Celsius', '°C')
            ax.set_xlabel(f'Temperature ({temp_units})')
            ax.tick_params(axis='x', labelrotation=45)

            # second row zoomed in on bottom Xm
            zoom_depth_range = 50
            ax = fig.add_subplot(2,4,5)
            ax.plot(temp_nonan, pres_nonan, 'k-', marker='.', markersize=3)
            ax.set_ylim([max_pres, max_pres-zoom_depth_range])
            ax.set_title(f'e. Bottom {zoom_depth_range}m', loc='left')
            ax.set_ylabel('Pressure (dbar)')
            ax.set_xlabel(f'Temperature ({temp_units})')
            x_lim_min = min(temp_nonan[pres_nonan>max_pres-zoom_depth_range])
            x_lim_max = max(temp_nonan[pres_nonan>max_pres-zoom_depth_range])
            # expand by 5 % on either side
            x_lim_range = (x_lim_max - x_lim_min) * 0.05
            if x_lim_range==0:
                x_lim_range = 0.05 * x_lim_max
            ax.set_xlim(x_lim_min - x_lim_range, x_lim_max + x_lim_range)
            ax.tick_params(axis='x', labelrotation=45)

            ax = fig.add_subplot(2,4,2)
            ax.plot(sal_nonan, pres_nonan, 'k-', marker='.')
            ax.set_title('b.', loc='left')
            ax.set_ylim([max_pres, min_depth])
            ax.set_ylabel('Pressure (dbar)')
            ax.set_xlabel(f'Salinity ({argo_n["PSAL_ADJUSTED_BGCArgoPlus"].units})')
            ax.ticklabel_format(useOffset=False, style='plain')
            ax.tick_params(axis='x', labelrotation=45)

            ax = fig.add_subplot(2,4,6)
            ax.plot(sal_nonan, pres_nonan, 'k-', marker='.', markersize=3)
            ax.set_title('f.', loc='left')
            ax.set_ylim([max_pres, max_pres-zoom_depth_range])
            ax.set_ylabel('Pressure (dbar)')
            ax.set_xlabel(f'Salinity ({argo_n["PSAL_ADJUSTED_BGCArgoPlus"].units})')
            x_lim_min = min(sal_nonan[pres_nonan>max_pres-zoom_depth_range])
            x_lim_max = max(sal_nonan[pres_nonan>max_pres-zoom_depth_range])
            # expand by 5 % on either side
            x_lim_range = (x_lim_max - x_lim_min) * 0.05
            if x_lim_range==0:
                x_lim_range = 0.05 * x_lim_max
            ax.set_xlim(x_lim_min - x_lim_range, x_lim_max + x_lim_range)
            ax.ticklabel_format(useOffset=False, style='plain')
            ax.tick_params(axis='x', labelrotation=45)

            axo = fig.add_subplot(2,4,3)
            axo.plot(o2_nonan, pres_nonan, 'k-', marker='.')
            # axo.set_title('c.', loc='left')
            axo.set_ylim([max_pres, min_depth])
            axo.set_ylabel('Pressure (dbar)')
            oxy_units = argo_n["DOXY_ADJUSTED_BGCArgoPlus"].attrs.get('units', ' ')
            # replace micro w/ mu symbol
            oxy_units = oxy_units.replace('micro', 'μ')
            axo.set_xlabel(f'Oxygen ({oxy_units})')
            axo.tick_params(axis='x', labelrotation=45)

            axo2 = fig.add_subplot(2,4,7)
            axo2.plot(o2_nonan, pres_nonan, 'k-', marker='.', markersize=3)
            # axo.set_title('g.', loc='left')
            axo2.set_ylim([max_pres, max_pres-zoom_depth_range])
            axo2.set_ylabel('Pressure (dbar)')
            oxy_units = argo_n["DOXY_ADJUSTED_BGCArgoPlus"].attrs.get('units', ' ')
            # replace micro w/ mu symbol
            oxy_units = oxy_units.replace('micro', 'μ')
            axo2.set_xlabel(f'Oxygen ({oxy_units})')
            x_lim_min = min(o2_nonan[pres_nonan>max_pres-zoom_depth_range])
            x_lim_max = max(o2_nonan[pres_nonan>max_pres-zoom_depth_range])
            # expand by 5 % on either side
            x_lim_range = (x_lim_max - x_lim_min) * 0.05
            if x_lim_range==0:
                x_lim_range = 0.05 * x_lim_max
            axo2.set_xlim(x_lim_min - x_lim_range, x_lim_max + x_lim_range)
            axo2.tick_params(axis='x', labelrotation=45)

            ax = fig.add_subplot(2,4,4)
            ax.plot(profile_deriv, profile_pres, 'b-', marker='x')
            ax.set_ylim([max_pres, min_depth])
            ax.set_title('d.', loc='left')
            ax.vlines(deriv_mean, ymin=max_pres, ymax=min_depth, colors='k')
            ax.vlines(deriv_mean-deriv_std, ymin=max_pres, ymax=min_depth, colors='k', linestyles='--')
            ax.vlines(deriv_mean+deriv_std, ymin=max_pres, ymax=min_depth, colors='k', linestyles='--')
            ax.vlines(deriv_mean-deriv_std*2, ymin=max_pres, ymax=min_depth, colors='k', linestyles='-.')
            ax.vlines(deriv_mean+deriv_std*2, ymin=max_pres, ymax=min_depth, colors='k', linestyles='-.')
            ax.vlines(deriv_mean-deriv_std*3, ymin=max_pres, ymax=min_depth, colors='k', linestyles='dotted')
            ax.vlines(deriv_mean+deriv_std*3, ymin=max_pres, ymax=min_depth, colors='k', linestyles='dotted')
            ax.set_ylabel('Pressure (dbar)')
            ax.set_xlabel(f'dO2/dz ({oxy_units}/dbar)')
            ax.tick_params(axis='x', labelrotation=45)

            ax = fig.add_subplot(2,4,8)
            ax.plot(profile_deriv, profile_pres, 'b-', marker='x', markersize=3)
            ax.set_ylim([max_pres, max_pres-zoom_depth_range])
            ax.set_title('h.', loc='left')
            ax.vlines(deriv_mean, ymin=max_pres, ymax=min_depth, colors='k')
            ax.vlines(deriv_mean-deriv_std, ymin=max_pres, ymax=min_depth, colors='k', linestyles='--')
            ax.vlines(deriv_mean+deriv_std, ymin=max_pres, ymax=min_depth, colors='k', linestyles='--')
            ax.vlines(deriv_mean-deriv_std*2, ymin=max_pres, ymax=min_depth, colors='k', linestyles='-.')
            ax.vlines(deriv_mean+deriv_std*2, ymin=max_pres, ymax=min_depth, colors='k', linestyles='-.')
            ax.vlines(deriv_mean-deriv_std*3, ymin=max_pres, ymax=min_depth, colors='k', linestyles='dotted')
            ax.vlines(deriv_mean+deriv_std*3, ymin=max_pres, ymax=min_depth, colors='k', linestyles='dotted')
            ax.set_ylabel('Pressure (dbar)')
            ax.set_xlabel(f'dO2/dz ({oxy_units}/dbar)')
            # x_lim_min = min(profile_deriv[profile_pres>max_pres-zoom_depth_range])
            # x_lim_max = max(profile_deriv[profile_pres>max_pres-zoom_depth_range])
            # # expand by 5 % on either side
            # x_lim_range = (x_lim_max - x_lim_min) * 0.05
            # ax.set_xlim(x_lim_min - x_lim_range, x_lim_max + x_lim_range)
            ax.tick_params(axis='x', labelrotation=45)

            if np.any(outlier_index):
                pres_n = deriv_df['pres'][outlier_index]
                o2_n = deriv_df['o2_conc'][outlier_index]
                axo.plot(o2_n, pres_n, 'ro', markerfacecolor="None")
                axo.set_title('c. Red = Removed points')
                axo2.plot(o2_n, pres_n, 'ro', markerfacecolor="None")
                axo2.set_title('g. Red = Removed points')

                # remove points based on this function 
                if o2_n.size==1: # if only one point
                    o2_index = argo_n['DOXY_ADJUSTED_BGCArgoPlus'][p,:]==o2_n
                    pres_index = argo_n['PRES_ADJUSTED_BGCArgoPlus'][p,:]==pres_n
                    argo_n['DOXY_ADJUSTED_BGCArgoPlus'][p,o2_index & pres_index] = np.nan
                    if verbose:
                        print('Storing bottom oxygen removal points for profile ' + str(p))
                    outlier_df_n = pd.DataFrame({                 # creates a new "outlier_df" w wmo, variable name, profile, date, level, pressure, deletion reason 
                        "Float Number":[wmo],
                        "Variable":['DOXY_ADJUSTED_BGCArgoPlus'],
                        "N_PROF":[p],
                        "Date":[pd.to_datetime(argo_n['JULD'][p].values).strftime('%Y-%m-%d %H:%M:%S')],
                        "N_LEVELS":[argo_n['N_LEVELS'][o2_index & pres_index].values.item()],
                        'Pressure (dbar)':[pres_n],
                        'Deletion reason':['Bottom_Oxygen_Check']
                    })   
                    if len(outlier_df)==0:
                        outlier_df = outlier_df_n
                    else:
                        outlier_df = pd.concat([outlier_df, outlier_df_n])

                else: # if multiple points
                    if verbose:
                        print('Storing bottom oxygen removal points for profile ' + str(p))
                                           
                    for i in o2_n.index:
                        o2_index = argo_n['DOXY_ADJUSTED_BGCArgoPlus'][p,:]==o2_n[i]
                        pres_index = argo_n['PRES_ADJUSTED_BGCArgoPlus'][p,:]==pres_n[i]
                        argo_n['DOXY_ADJUSTED_BGCArgoPlus'][p,o2_index & pres_index] = np.nan
                        outlier_df_n = pd.DataFrame({                 # creates a new "outlier_df" w wmo, variable name, profile, date, level, pressure, deletion reason 
                                                "Float Number":[wmo],
                                                "Variable":['DOXY_ADJUSTED_BGCArgoPlus'],
                                                "N_PROF":[p],
                                                "Date":[pd.to_datetime(argo_n['JULD'][p].values).strftime('%Y-%m-%d %H:%M:%S')],
                                                "N_LEVELS":[argo_n['N_LEVELS'][o2_index & pres_index].values.item()],
                                                'Pressure (dbar)':[pres_n[i]],
                                                'Deletion reason':['Bottom_Oxygen_Check']
                                            })   
                        if len(outlier_df)==0:
                            outlier_df = outlier_df_n
                        else:
                            outlier_df = pd.concat([outlier_df, outlier_df_n])  

            plt.tight_layout()

            plt.savefig(f'{float_RO_profile_dir}{plot_filename}.png')
            plt.close(fig)
            
    # if len(outlier_df)>0:
    #     savename = group_outlier_dir + 'outliers_'+ wmo + '_Sprof_filtered_bottomoxygen_XXXX_X_X_X.csv'
    #     outlier_df.to_csv(savename, mode='w', index=False, header=True)
    return argo_n, outlier_df


def sensor_flag_wrapper(sprof_path, file, flags_to_remove, output_dir, argo_index, dataset_version, verbose=False):
    if verbose: print('Loading Sprof file: ' + file)

    try:
        # if verbose:
        #     print('Loading Sprof file: ' + file)
        try:
            argo_n = xr.open_dataset(sprof_path + file)
        except Exception as e:
            msg = f"{file} failed to open: {e}"
            print(msg, flush=True)
            return {
                "file": file,
                "status": "error",
                "message": str(e),
            }
            
        # load sensor information from meta file, add information about air calibrated or not, add information on sensor types
        if verbose:
            print('Adding sensor info')
        try:
            argo_n = apply_sensor_info(argo_n, file, sprof_path, argo_index, verbose)
        except Exception as e:
            print('Failed to add sensor info for ' + str(file) + ' check if meta.nc file exists: ' + str(e))
        
        # apply existing GDAC flags
        if verbose:
            print('Applying GDAC QC Flags')
        argo_n = apply_flags(argo_n, flags_to_remove, verbose)

        outlier_df = pd.DataFrame({
            "Float Number": pd.Series(dtype="str"),
            "Variable": pd.Series(dtype="str"),
            "N_PROF": pd.Series(dtype="int"),
            "Date": pd.Series(dtype="str"),  # or "datetime64[ns]" if preferred
            "N_LEVELS": pd.Series(dtype="int"),
            "Pressure (dbar)": pd.Series(dtype="float"),
            "Deletion reason": pd.Series(dtype="str")
        })
        # Remove surface data (everything above Xm - seems to be contaminated with in-air data
        if verbose:
            print('Checking for and removing any surface data')
        argo_n, outlier_df = surface_removal(argo_n, 2, outlier_df, verbose)

        # Find and remove density inversions
        if verbose:
            print('Checking for and removing any significant density inversions')
        argo_n, outlier_df  = density_inversion(argo_n, sprof_path + '../processed/Processed_Figures/', outlier_df , verbose)

        # check to see if bottom oxygen is off and remove 
        if 'DOXY_ADJUSTED_BGCArgoPlus' in argo_n.keys() and np.any(~np.isnan(argo_n.DOXY_ADJUSTED_BGCArgoPlus)):
            if verbose:
                print('Looking for jumps in bottom oxygen data')
            argo_n, outlier_df = bottom_oxygen_check(argo_n, output_dir + '../outlier_file_collection/', outlier_df, verbose)

        if len(outlier_df)>0:
            wmo = str(argo_n['WMO_ID'].values)

            savename = output_dir + '../outlier_file_collection/' + 'outliers_'+ wmo + '_Sprof_filtered_automaticoutliers_XXXX_X_X_X.csv'
            if verbose:
                print(savename)
            outlier_df.to_csv(savename, mode='w', index=False, header=True)
    except Exception as e:
        print(f'Error processing {file}: {e}')
        traceback.print_exc()
        return {
            "file": file,
            "status": "error",
            "message": str(e),
        }
    # save as a filtered file
    filtered_file_name = str(file[:-3])+'_filtered.nc'
    path_name = os.path.join(output_dir, filtered_file_name)
    # def find_bad_variable(ds, outdir="debug_nc"):
    #     os.makedirs(outdir, exist_ok=True)
    
    #     for var in ds.variables:
    #         test_path = os.path.join(outdir, f"{var}.nc")
    #         try:
    #             ds[[var]].to_netcdf(test_path)
    #             print(f"OK: {var}")
    #         except Exception as e:
    #             print(f"BAD: {var} -> {e}")
    #             traceback.print_exc()
    if os.path.isfile(path_name): # delete file if it already exists
        os.remove(path_name)
    history_vars = [
            "HISTORY_INSTITUTION",
            "HISTORY_STEP",
            "HISTORY_SOFTWARE",
            "HISTORY_SOFTWARE_RELEASE",
            "HISTORY_REFERENCE",
            "HISTORY_DATE",
            "HISTORY_ACTION",
            "HISTORY_PARAMETER",
            "HISTORY_QCTEST",
        ]
    
    argo_n.encoding.pop("unlimited_dims", None)
    
    if "N_HISTORY" not in argo_n.dims or argo_n.sizes.get("N_HISTORY", 0) == 0:
        drop_vars = [v for v in history_vars if v in argo_n.variables]
        argo_n = argo_n.drop_vars(drop_vars)

    # Add to the history attribute of the dataset to include dataset version and processing date
    current_history = argo_n.attrs.get("history", "")
    processing_date = datetime.now().strftime("%Y-%m-%d")
    new_history_entry = f"BGC-Argo+ {dataset_version} processing on {processing_date}"
    if current_history:
        argo_n.attrs["history"] = current_history + "; " + new_history_entry
    else:
        argo_n.attrs["history"] = new_history_entry

    if verbose:
        print(f'Saving out filtered file {filtered_file_name}')
    try:

        # Strip char_dim_name encoding to avoid STRING length mismatch warnings
        for var in argo_n.variables:
            argo_n[var].encoding.pop('char_dim_name', None)
        # find_bad_variable(argo_n)
        # for var in argo_n.variables:
        #     da = argo_n[var]
        #     print(var, da.dtype)
        #     if da.dtype == object:
        #         print("OBJECT DTYPE:", var, da.values.flat[:10])
        #     if da.dtype.kind in "iuf":
        #         for k in ["_FillValue", "missing_value", "fill_value"]:
        #             if k in da.attrs:
        #                 print("ATTR", var, k, repr(da.attrs[k]), type(da.attrs[k]))
        #             if k in da.encoding:
        #                 print("ENC", var, k, repr(da.encoding[k]), type(da.encoding[k]))

        argo_n.to_netcdf(path_name)
    except Exception as e:
        print(f'First attempt to save netcdf failed: {e}')
        traceback.print_exc()        # try to save again after removing encoding of "char_dim_name" which sometimes causes problems:
        print('Removing some encoding and trying again')

        try:
            if os.path.isfile(path_name): # delete file if it already exists
                os.remove(path_name)
            encode_var = "char_dim_name"
            for var in argo_n.variables:
                enc = argo_n[var].encoding
                if encode_var in enc:
                    del argo_n[var].encoding[encode_var]
            print('encoding removed, trying to save')

            argo_n.to_netcdf(path_name)
        except Exception as e2:
            print(f'{filtered_file_name} failed to save: {e2}')
            traceback.print_exc()
    finally:
        argo_n.close()
    return


def flag_wrapper_only(sprof_path, file, flags_to_remove, output_dir, verbose=False):
    try:
        if verbose:
            print('Loading Sprof file: ' + file)
        try:
            argo_n = xr.open_dataset(sprof_path + file)
        except:
            print(str(file) + ' failed to open')
            return
        # # load sensor information from meta file, add information about air calibrated or not, add information on sensor types
        # if verbose:
        #     print('Adding sensor info')
        # argo_n = apply_sensor_info(argo_n, file, sprof_path)
        
        # apply existing GDAC flags
        if verbose:
            print('Applying GDAC QC Flags')
        argo_n = apply_flags(argo_n, flags_to_remove, verbose)

        # outlier_df = pd.DataFrame({
        #     "Float Number": pd.Series(dtype="str"),
        #     "Variable": pd.Series(dtype="str"),
        #     "N_PROF": pd.Series(dtype="int"),
        #     "Date": pd.Series(dtype="str"),  # or "datetime64[ns]" if preferred
        #     "N_LEVELS": pd.Series(dtype="int"),
        #     "Pressure (dbar)": pd.Series(dtype="float"),
        #     "Deletion reason": pd.Series(dtype="str")
        # })
        # # Remove surface data (everything above Xm - seems to be contaminated with in-air data
        # if verbose:
        #     print('Checking for and removing any surface data')
        # argo_n, outlier_df = surface_removal(argo_n, 2, outlier_df)

        # # Find and remove density inversions
        # if verbose:
        #     print('Checking for and removing any significant density inversions')
        # argo_n, outlier_df  = density_inversion(argo_n, sprof_path + '../processed/Processed_Figures/', outlier_df , verbose)

        # # check to see if bottom oxygen is off and remove 
        # if verbose:
        #     print('Looking for jumps in bottom oxygen data')
        # argo_n, outlier_df = bottom_oxygen_check(argo_n, output_dir + '../outlier_file_collection/', outlier_df, verbose)

        # if len(outlier_df)>0:
        #     wmo = str(argo_n['WMO_ID'].values)

        #     savename = output_dir + '../outlier_file_collection/' + 'outliers_'+ wmo + '_Sprof_filtered_automaticoutliers_XXXX_X_X_X.csv'
        #     if verbose:
        #         print(savename)
        #     outlier_df.to_csv(savename, mode='w', index=False, header=True)
    except:
        print('Error in ' + str(file))
    # save as a filtered file
    filtered_file_name = str(file[:-3])+'_flags_mode_only.nc'

    if verbose:
        print('Saving out flag_removed file')
    try:

        if os.path.isfile(output_dir + filtered_file_name): # delete file if it already exists
            os.remove(output_dir + filtered_file_name)

        # Strip char_dim_name encoding to avoid STRING length mismatch warnings
        for var in argo_n.variables:
            argo_n[var].encoding.pop('char_dim_name', None)
        argo_n.to_netcdf(output_dir+filtered_file_name)
    except:
        # try to save again after removing encoding of "char_dim_name" which sometimes causes problems:
        print('First attempt to save netcdf failed, removing some encoding and trying again')
        try:
            if os.path.isfile(output_dir + filtered_file_name): # delete file if it already exists
                os.remove(output_dir + filtered_file_name)
            encode_var = "char_dim_name"
            for var in argo_n.variables:
                enc = argo_n[var].encoding
                if encode_var in enc:
                    del argo_n[var].encoding[encode_var]
            print('encoding removed, trying to save')
            argo_n.to_netcdf(output_dir+filtered_file_name)
        except:
            print(str(file[:-3]) + ' failed to save')
    argo_n.close()
    return
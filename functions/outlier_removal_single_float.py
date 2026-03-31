# attempting a python script for individual float plotting... 

import numpy as np
import pandas as pd
import xarray as xr
import os
import functions.find_wiggles_peak_version as wiggle
# import functions.find_density_inversions as func_inv
import random
import datetime
from io import StringIO
import gsw
import matplotlib as mpl
import time 

# import matplotlib.pyplot as plt
# import scipy.stats as stats
# import functions.float_download_sprof_meta as fl_download
#import functions.derived_parameter_utilities as fl_calcs
# import functions.float_z_interpolation_2 as fl_interp

#import functions.find_bad_oxy_bottom as func_badoxybottom
# from datetime import date
# from multiprocessing import Pool
# import cartopy.crs as ccrs
# import math
# import pickle
# import sys

#from plot_profiles_gdac import plot_gdac_profiles as pgp

#Plotly Packages
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Dash Packages
from dash import Dash, dcc, html, dash_table, ctx
import dash_bootstrap_components as dbc
import copy
from dash.dependencies import Input, Output, State
import json
import dash
# %matplotlib inline
suffix = 'BGCArgoPlus'
plot_vars_all=[ 'TEMP_ADJUSTED' + '_' + suffix,
                'PSAL_ADJUSTED' + '_' + suffix,
                'NITRATE_ADJUSTED' + '_' + suffix,
                'DOXY_ADJUSTED' + '_' + suffix,
                'PH_IN_SITU_TOTAL_ADJUSTED' + '_' + suffix,]
# stop_script = False  # Global flag
class GracefulExit(Exception):
    pass

# Function to load and serialize dataset
def load_dataset(output_dir, file, plot_vars_all, verbose):
    if verbose==True:
        print('initial dataset load')
    
    argo_n = xr.open_dataset(output_dir + file)

    # first, check to see if there are valid bgc data before progressing further 
    plot_vars = plot_vars_all.copy()
    # print(plot_vars)

    # t_s_present = True
    bgc_present=False
    for idx, var in enumerate(plot_vars_all):
        if var not in argo_n.keys():
            plot_vars.remove(var)
        else: # Check if there are any valid bgc data. If not, don't bother doing outlier detection. 
            if var in ['NITRATE_ADJUSTED' + '_' + suffix,
                    'DOXY_ADJUSTED' + '_' + suffix,
                    'PH_IN_SITU_TOTAL_ADJUSTED' + '_' + suffix,]:
                if np.sum(~np.isnan(argo_n[var]))>0:
                    bgc_present = True
            # also check that there is valid T/S data
            # else: # for now, stop checking if there is "_RO" temp or sal data, show anyway
            #     if var in ['TEMP_ADJUSTED_RO',
            #         'PSAL_ADJUSTED_RO']:
            #         if np.sum(~np.isnan(argo_n[var]))==0:
            #             t_s_present = False
    if bgc_present==False:
        print('No valid bgc data in variables checked, skipping outlier removal')
        raise GracefulExit  # Raises a custom exception
        return None, None, None, None, [None, None, None]
    # if t_s_present==False:
    #     print('Missing either valid T or valid S, skipping outlier removal')
    #     raise GracefulExit  # Raises a custom exception
    #     return None, None, None, None, [None, None, None]
    # check both nitrate and oxygen for wiggly profiles
    if verbose:
        print('Checking for wiggles')
        tic = time.time()

    vars_to_check_for_wiggles = ['DOXY_ADJUSTED' + '_' + suffix, 'NITRATE_ADJUSTED' + '_' + suffix, 'PH_IN_SITU_TOTAL_ADJUSTED' + '_' + suffix]
    for var in vars_to_check_for_wiggles:
        if var not in argo_n.keys():
            vars_to_check_for_wiggles.remove(var)
    for idx, var in enumerate(vars_to_check_for_wiggles):
        if var in argo_n.keys():
            # print(var)
            problem_profiles_var = wiggle.peak_detect_wiggles(output_dir, argo_n, var=var)
            # print(problem_profiles_var)
            if idx==0:
                problem_profiles = problem_profiles_var
            elif len(problem_profiles_var)>0:
                problem_profiles = np.concatenate((problem_profiles, problem_profiles_var))
            # print(problem_profiles)
    if verbose:
        toc = time.time()
        # print(str(np.round(toc-tic,2)) + 'sec Elapsed')

    problem_profiles = np.unique(problem_profiles)

    # # Identify density inversions
    # if verbose==True:
    #     print('Dataset loaded, starting density inversions')

    # levels_inv, prof_inv, points_dens_inversion = func_inv.density_inversion(argo_n, verbose)
    # if verbose==True:
    #     print('Density inversion function finished')

    # Identify oxygen increase in last profile point
    #levels_bo, prof_bo, points_bad_oxy = func_badoxybottom.badoxy_bottom(argo_n)

    argo_vars = argo_n[plot_vars+['PRES_ADJUSTED' + '_' + suffix,'CYCLE_NUMBER','JULD','LATITUDE','LONGITUDE', 'PRES_ADJUSTED', 'TEMP_ADJUSTED', 'PSAL_ADJUSTED','TEMP', 'PSAL', 'PRES']] # adding TEMP_ADJUSTED and PSAL_ADJUSTED so we can use them if data is missing
    print(plot_vars)

    df = argo_vars.to_dataframe()
    if verbose:
        print('argo_vars in dataframe')
    df = df.reset_index()
    
    # Set number of profiles to plot per figure as well as group profiles into ranges
    num_profiles_per_fig = 5
    finishing_prof = 0
    label_lists = []
    for i, data in df.groupby("N_PROF"):
        if i%num_profiles_per_fig == 0:
            starting_prof = i
        finishing_prof = starting_prof + 4
        label = str(starting_prof) + "-" + str(finishing_prof)
        label_list = [label] * len(data)
        label_lists = label_lists + label_list
    df["Profile Range"] = label_lists
    if verbose:
        print('Here')
    # # assign integer number tied to profile range i.e. 1:'0-4', 2:'5-9'
    marker_dict = {}
    for i, profiles in enumerate(df['Profile Range'].unique()):
        profiles_expanded = range(int(profiles.split('-')[0]), int(profiles.split('-')[1]))
        color = '#000000'
        for problem_profile in problem_profiles:
            if problem_profile in profiles_expanded:
                color = '#eb3434'
        marker_dict[i] = {
            'label':profiles,
            'style':{'color': color}
        }

    # try setting pressure range for each time series slider range to keep similar data together
    pressure_ranges_for_ts =[20, 40, 70, 100, 150, 200, 250, 350, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000]

    pressure_labels_list = []

    pressure_labels_list.append( '%i-%i'%(0, pressure_ranges_for_ts[0]) )
    for idx, lev in enumerate(pressure_ranges_for_ts[0:-1]):
        pressure_labels_list.append( '%i-%i'%(pressure_ranges_for_ts[idx]+1, pressure_ranges_for_ts[idx+1]) )
    pressure_labels_list

    pres_name = 'PRES_ADJUSTED' + '_' + suffix
    if np.sum(~np.isnan(df[pres_name]))==0: 
        pres_name = 'PRES_ADJUSTED'
        if np.sum(~np.isnan(df[pres_name]))==0: 
            pres_name = 'PRES'

    idx = np.digitize(df[pres_name], pressure_ranges_for_ts[0:-1])
    idx, len(pressure_labels_list)

    df['Pressure Range'] = [pressure_labels_list[i] for i in idx]

    # # assign integer number tied to level range i.e. 1:'0-4', 2:'5-9'
    press_ranges_sorted = df['Pressure Range'].unique()
    # Sorting function based on the lower bound of the range
    def sort_key(range_str):
        return int(range_str.split('-')[0])

    sorted_arr = sorted(press_ranges_sorted, key=sort_key)

    pressures_dict = {i: levels for i, levels in enumerate(sorted_arr)}
    if verbose:
        print('Dataset Loading Finished')
    return argo_n, df, marker_dict, pressures_dict#, [levels_inv, prof_inv, points_dens_inversion]

def create_map(argo_n, profile_lat_lon):
    n_prof = len(argo_n.N_PROF)

    # Create Map Figure using plotly express
    # Create scatter geo object
    map_fig = px.scatter_geo(
            lon = argo_n['LONGITUDE'],
            lat = argo_n['LATITUDE'],
            hover_name=argo_n['CYCLE_NUMBER'],
            color = argo_n['CYCLE_NUMBER'],
            )
    # map_fig.add_scattergeo(lon = profile_lat_lon[2:3], lat = profile_lat_lon[0:1]) # having trouble adding to this map

    # go.Scatter(x=[0,1,2,0], y=[0,2,0,0],)
               
    # Update figure to show grids and start zoom to where plotted points are
    mid_lon = argo_n['LONGITUDE'].values[int(len(argo_n['LONGITUDE'].values)/2)]
    mid_lat = argo_n['LATITUDE'].values[int(len(argo_n['LONGITUDE'].values)/2)]
    map_fig.update_geos(#fitbounds="locations",
                    lataxis_showgrid=True,
                    lonaxis_showgrid=True,
                    center=dict(lon=mid_lon, lat = mid_lat),
                    lataxis_range=[mid_lat-30, mid_lat+30], lonaxis_range=[mid_lon-40, mid_lon+40])
    
    # add information about variables and removed data
    sensor_removal_text=''
    var_names = ['DOXY', 'NITRATE', 'PH_IN_SITU_TOTAL'] #'TEMP', 'PSAL', 
    for var in var_names:
        if var+ '_ADJUSTED' + '_' + suffix in argo_n.keys():
            per_data_remaining = np.sum(~np.isnan(argo_n[var + '_ADJUSTED' + '_' + suffix])) / np.sum(~np.isnan(argo_n[var + '_ADJUSTED']))*100

            sensor_removal_text= sensor_removal_text + '<br>' + var + ' Adjusted QC flags removed: ' + str(argo_n[var +'_ADJUSTED_QC_n_4_removed'].values + argo_n[var +'_ADJUSTED_QC_n_3_removed'].values) \
                + '<br>' + var #+ ' Adjusted non Delayed-Mode profiles removed: ' + str(argo_n[var +'_profile_removed_not_D'].values) + ' of ' + str(len(argo_n['N_PROF'])) \
                # + "<br>Percent  <b>" + var + "</b> 'ADJUSTED' data remaining: <b>" + str(np.round(per_data_remaining.values,2)) + ' %</b>'

            model_var = [varname for varname in argo_n.data_vars if varname.endswith(var +'_model')]
            if model_var:
                sensor_removal_text = sensor_removal_text + '<br>     ' + model_var[0] + ': ' + str(argo_n[model_var[0]].values)
                
    title_text = (
        'Map of Float Profiles: ' + 
        'WMO: ' + str(argo_n.WMO_ID.values) + ', ' + str(n_prof) + ' Profiles, ' +
        str(argo_n['JULD'][0].dt.year.values) + '/' + str(argo_n['JULD'][0].dt.month.values) + ' to ' + str(argo_n['JULD'][-1].dt.year.values) + '/' + \
            str(argo_n['JULD'][-1].dt.month.values) + \
            ', Project: ' + argo_n['PROJECT_NAME'].values[0].decode('utf-8').strip() + \
                ', PI: ' + argo_n['PI_NAME'].values[0].decode('utf-8').strip() + \
                ', Float Type: ' + argo_n['PLATFORM_TYPE'].values[0].decode('utf-8').strip() + \
                                                          sensor_removal_text    
    )
    # Update Figure title text
    map_fig.update_layout(title_text=title_text)
    return map_fig

# Create outlier table function passing in outlier point dataframe
def create_table(df):
    # Using dash datatable, reading in a pandas dataframe
    t = dash_table.DataTable(df.to_dict('records'), [{"name": i, "id": i} for i in df.columns], id="Outlier_Table",
                             # Enable sorting
                         sort_action="native",
                         page_size=10,
                         style_table={"overflowX": "auto"},)
    return t

def pre_load_data(output_dir, file, researcher, port_num, verbose): # this is directly called from Float_Processing 
    # filtered_list = file
    # global stop_script
    if verbose==True:
        print('Starting dataset load')
    try:
        # print(output_dir)
        print(file)
        argo_n, df, marker_dict, pressures_dict \
            = load_dataset(output_dir, file, plot_vars_all, verbose)
        success = True
    except GracefulExit:
        print('Exiting gracefully')
        success = False
        return success
    
    
    # if stop_script:
    #     return
    
    def outlier_removal_plot(output_dir, file, argo_n, df, marker_dict, pressures_dict, researcher):
        # Select vairables to plot
        
        # function to create figures 
        def create_figure(df, profile_slider, skip_points=[], deleted_variables=[], y_scale=0):

            colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
            
            plot_vars = plot_vars_all.copy()
            for var in plot_vars_all:
                if not var in df.keys():
                    plot_vars.remove(var)
            #Set the width fraction of each subplot based on how many variables are present
            width = [1/len(plot_vars)] * len(plot_vars)
            
            dfs = df
            for dvar in np.unique(deleted_variables):
                set_nan = []
                for i in range(len(skip_points)):
                    if deleted_variables[i] == dvar:
                        set_nan.append(skip_points[i])
                # dfs[dvar][set_nan] = np.nan
                dfs.loc[set_nan, dvar] = np.nan
    
            # Filter float df to only include profile range that profile slider selected
            filtered_df = dfs[dfs["Profile Range"] == marker_dict[profile_slider]['label']]

            if y_scale==0:
                y_range=[6000,-5]
            elif y_scale==1:
                y_range=[10, -5]
            elif y_scale==2:
                y_range=[50,-5]
            elif y_scale==3:
                y_range=[2010, 1490]
            elif y_scale==4:
                y_range=[200,-5]
            elif y_scale==5:
                y_range=[500,-5]
            elif y_scale==-1:
                y_range=[np.nanmax(filtered_df['PRES'].values)+10, np.nanmax(filtered_df['PRES'].values)-210]
            elif y_scale==-2:
                y_range=[np.nanmax(filtered_df['PRES'].values)+10, np.nanmax(filtered_df['PRES'].values)-510]

            # Create sub plot figure object
            variable_fig = make_subplots(rows=1, cols=len(plot_vars),
                            column_width= width, row_heights=[1], shared_yaxes=True)
            variable_plot_position = {}
            # Loop through profiles available in filtered float df
            profile_dec_years = []
            profile_lats = []
            profile_lons = []

            profile_count = 0
            for prof, data_unused in filtered_df.groupby("N_PROF"):
                # Select data from one profile
                prof_df = filtered_df[filtered_df["N_PROF"] == prof]

                pres_name = 'PRES_ADJUSTED' + '_' + suffix
                pres = prof_df[pres_name]
                if np.sum(~np.isnan(pres))==0: # if PRES data is missing, might mean that only non-delayed mode is available, use "TEMP_ADJUSTED" instead
                    pres_name = 'PRES_ADJUSTED'
                    pres = prof_df[pres_name]
                    if np.sum(~np.isnan(pres))==0: # if PRES data is missing, might mean that only non-delayed mode is available, use "TEMP_ADJUSTED" instead
                        pres_name = 'PRES'
                        pres = prof_df[pres_name]

                temp_name = 'TEMP_ADJUSTED' + '_' + suffix
                temp_prof = prof_df[temp_name]
                if np.sum(~np.isnan(temp_prof))==0: # if temp data is missing, might mean that only non-delayed mode is available, use "TEMP_ADJUSTED" instead
                    temp_name = 'TEMP_ADJUSTED'
                    temp_prof = prof_df[temp_name]
                    if np.sum(~np.isnan(temp_prof))==0: # check if TEMP_ADJUSTED is also missing and use "TEMP" instead
                        temp_name = 'TEMP'
                        temp_prof = prof_df[temp_name]

                psal_name = 'PSAL_ADJUSTED' + '_' + suffix
                psal_prof = prof_df[psal_name]
                if np.sum(~np.isnan(psal_prof))==0: # if temp data is missing, might mean that only non-delayed mode is available, use "TEMP_ADJUSTED" instead
                    psal_name = 'PSAL_ADJUSTED'
                    psal_prof = prof_df[psal_name]
                    if np.sum(~np.isnan(psal_prof))==0: # check if TEMP_ADJUSTED is also missing and use "TEMP" instead
                        psal_name = 'PSAL'
                        psal_prof = prof_df[psal_name]
                prof_meta = prof_df['N_PROF']
                levels_meta = prof_df['N_LEVELS']
                cycle_meta = prof_df['CYCLE_NUMBER']
                juld_meta = [date_obj.strftime('%Y%m%d') for date_obj in prof_df['JULD']]

                # calculate decimal year of profiles so that you can pass the value to the time series plots
                decimal_year = prof_df['JULD'].dt.year.values + (prof_df['JULD'].dt.dayofyear.values)/365.25
                profile_dec_years.append(decimal_year[0])
       
                lat_meta = prof_df['LATITUDE']
                lon_meta = prof_df['LONGITUDE']
                psal_meta = prof_df['PSAL_ADJUSTED' + '_' + suffix]
                profile_lats.append(lat_meta.values[0])  # store to add to map
                profile_lons.append(lon_meta.values[0])  # store to add to map

                index_meta = prof_df.index

                for i, var in enumerate(plot_vars,1):
                    if var not in variable_plot_position:
                        variable_plot_position[var] = i

                    # plots better and seems to overall work without the mask and dropping nans. 
                    # Probably another way to fix things too, so just commenting out code above
                    symbol = 'circle-dot' 
                    line_width = 2
                    marker_size=4

                    if var=='TEMP_ADJUSTED' + '_' + suffix:
                        if temp_name !='TEMP_ADJUSTED' + '_' + suffix:
                            symbol = 'x-thin' 
                            line_width=2
                            marker_size=5
                        var = temp_name
                        marker_size=3
                    if var=='PSAL_ADJUSTED' + '_' + suffix:
                        if psal_name !='PSAL_ADJUSTED' + '_' + suffix:
                            symbol = 'x-thin' 
                            line_width=2
                            marker_size=5
                        var = psal_name
                        marker_size=3

                            
                    # if var=='PSAL_ADJUSTED_RO':
                    #     marker_size=3

                    #     if np.sum(~np.isnan(data))==0: # if temp data is missing, might mean that only non-delayed mode is available, use "TEMP_ADJUSTED" instead
                    #         var = 'PSAL_ADJUSTED'
                    #         data = prof_df[var]
                    #         if np.sum(~np.isnan(data))==0: # check if TEMP_ADJUSTED is also missing and use "TEMP" instead
                    #             var = 'PSAL'
                    #             data = prof_df[var]
                    #         symbol = 'x-thin'
                    #         line_width=2
                    #         marker_size=5
                    data = prof_df[var]

                    # Create meta data lists to include in hover data on figure
                    
                   
                    # doxy_meta = prof_df['DOXY_ADJUSTED_RO']

                    var_meta = len(prof_df.index) * [var]
                    n_meta = len(prof_df.index) * [i]
                    # Plot or "add" trace to figure with scatter plot function establishing custom data and hover template as well as name of trace for legend
                    variable_fig.add_trace(go.Scattergl(x=data, y=pres, name=f"Profile {prof}", line={'color':colors[profile_count], 'width':1}, 
                                                        mode='markers+lines', marker_line_width = line_width, marker_symbol=symbol, 
                                                        marker_line_color=colors[profile_count], 
                                                        marker_size = marker_size,
                                                        connectgaps=True, customdata=np.stack((
                        index_meta, prof_meta, levels_meta, cycle_meta, juld_meta, np.round(lat_meta,2), np.round(lon_meta,2), var_meta, n_meta, pres
                    ), axis=1),
                    hovertemplate='<br>'.join(['pressure: %{y}',
                                        'variable: %{x}',
                                            'index: %{customdata[0]}',
                                            'profile: %{customdata[1]}',
                                            'level: %{customdata[2]}',
                            #  #                                  'cycle: %{customdata[3]}',
                                            'juld: %{customdata[4]}',
                                            'lat: %{customdata[5]}',
                                            'lon: %{customdata[6]}'])
                                                ), row=1, col=i)
                    variable_fig.update_xaxes(title_text=var, showgrid=True, row=1, col=i)
                    variable_fig.update_yaxes(title_text=pres_name, showgrid=True, autorange="reversed", row=1, col=i)
                
                    if var=='DOXY_ADJUSTED' + '_' + suffix:
                        surf_psal = np.mean(prof_df[psal_name].where(prof_df[pres_name]<10))
                        # if np.isnan(surf_psal):
                        #     surf_psal = np.mean(prof_df['PSAL_ADJUSTED'].where(prof_df[pres_name]<10))
                        surf_temp = np.mean(prof_df[temp_name].where(prof_df[pres_name]<10))
                        # if np.isnan(surf_temp):
                        #     surf_temp = np.mean(prof_df[temp_name].where(prof_df[pres_name]<10))

                        o2_sat_conc = gsw.O2sol_SP_pt(surf_psal, surf_temp)
                        variable_fig.add_trace(go.Scattergl(y=[-2, -30], x=[o2_sat_conc, o2_sat_conc], name=f"O2 Sat {prof}", 
                                                            line={'color':colors[profile_count], 'width':3}, mode='lines'),
                                                              row=1, col=i)
                    def autoscale_plot(x, y, x_min, x_max, y_margin=0.1): # auto x-limit calc
                        mask = (x >= x_min) & (x <= x_max)
                        y_data = y[mask]
                        if np.isnan(y_data).all():
                            raise ValueError("No valid data within specified x-range.")
                        y_min = np.nanmin(y_data)
                        y_max = np.nanmax(y_data)
                        y_range = y_max - y_min
                        y_lims = [y_min - y_range * y_margin, y_max + y_range * y_margin]
                        return y_lims
                    
                    if y_range[0] ==6000:
                        variable_fig.update_yaxes(autorange="reversed")
                    
                    else:
                        x_range = autoscale_plot(filtered_df[pres_name], filtered_df[var], y_range[1], y_range[0]) 
                        variable_fig.update_yaxes(range=y_range, autorange=False)
                        variable_fig.update_xaxes(range=x_range, autorange=False, row=1, col=i)
                    
                profile_count = profile_count + 1  
                # Change marker size
                # variable_fig.update_traces(marker=dict(size=4))
                # print(marker_dict[profile_slider])
            # Create axis titles and enable grid also reverse y axis
            # print(juld_meta[0])
            variable_fig.update_layout(title_text = 'Profiles ' + str(marker_dict[profile_slider]['label']) + 
                                        '. If PSAL or TEMP have "x" symbols, non "RO" data is being used. Date of last profile shown: ' +
                                        juld_meta[0], height = 1000)
            profile_dec_year_range = [np.round(min(profile_dec_years),3), np.round(max(profile_dec_years),3)]
            profile_lat_lon = [np.round(min(profile_lats),3), np.round(max(profile_lats),3), np.round(min(profile_lons),3), np.round(max(profile_lons),3)]
            # print(profile_lat_lon)

            # print(profile_dec_year_range)

            return variable_fig, variable_plot_position, profile_dec_year_range, profile_lat_lon
        
        def create_figure_ts(df, level_slider, prof_dec_year_range, skip_points=[], deleted_variables=[]):
            # print(prof_dec_year_range[0])
            plot_vars = plot_vars_all.copy()
            for var in plot_vars_all:
                if not var in argo_n.keys():
                    plot_vars.remove(var)
            width = [1/len(plot_vars)] * len(plot_vars)
    
            # create a df for the profiles
            # Drop skip points list from dataframe
            dfs = df
            for dvar in np.unique(deleted_variables):
                set_nan = []
                for i in range(len(skip_points)):
                    if deleted_variables[i] == dvar:
                        set_nan.append(skip_points[i])
                # dfs[dvar][set_nan] = np.nan
                dfs.loc[set_nan, dvar] = np.nan
            # Filter float df to only include levels range that level slider selected
            # filtered_df = dfs[dfs["Level Range"] == levels_dict[level_slider]]
           
            # Filter float df to only include pressure range that level slider selected
            filtered_df = dfs[dfs["Pressure Range"] == pressures_dict[level_slider]]

            # print(filtered_df)
            # Create subplot figure object
            variable_fig = make_subplots(rows=len(plot_vars), cols=1,
                            column_width= [1], row_heights=width, shared_xaxes=True,
                            vertical_spacing=0.05)
    
            n_press_levs = 5
            cmap = mpl.colormaps["tab10"]

            colors_ts = [(mpl.colors.rgb2hex(cmap(n))) for n in range(n_press_levs)]

            # num_levels_per_fig = len(np.arange(int(pressure_range_n[0]), int(pressure_range_n[1]), press_step))
            # colors_ts = lambda n: ["#%06x" % random.randint(0, 0xFFFFFF) for _ in range(n)]
            # colors_ts = colors_ts(num_levels_per_fig)
            pres_name = 'PRES_ADJUSTED' + '_' + suffix
            if np.sum(~np.isnan(df[pres_name]))==0: 
                pres_name = 'PRES_ADJUSTED'
                if np.sum(~np.isnan(df[pres_name]))==0: 
                    pres_name = 'PRES'

            pressure_range_n = pressures_dict[level_slider].split('-')

            pressure_range_ts = np.linspace(int(pressure_range_n[0]), int(pressure_range_n[1]), n_press_levs)
            level_count = 0
            # Loop through pressure ranges available in filtered float df
            for idx_p, press_n in enumerate(pressure_range_ts[:-1]): 
                # print(lev)
                # print(len(filtered_df.N_LEVELS))
                # Select data from one level
                level_df = filtered_df[np.logical_and(filtered_df[pres_name]>=pressure_range_ts[idx_p]-1, filtered_df[pres_name]<pressure_range_ts[idx_p+1]+1)] # added -1 and +1 to account for the space between integer start/end values
                
                decimal_year = level_df['JULD'].dt.year.values + (level_df['JULD'].dt.dayofyear.values)/365.25
                time = decimal_year
                # time = [ date_obj.strftime('%Y%m%d') for date_obj in level_df['JULD'].where(mask).dropna()]
                # Create meta data lists to include in hover data on figure
                # pres_name = 'PRES_ADJUSTED_RO'
                pres = level_df[pres_name]
                # if np.sum(~np.isnan(pres))==0: # if PRES data is missing, might mean that only non-delayed mode is available, use "TEMP_ADJUSTED" instead
                #     pres_name = 'PRES_ADJUSTED'
                #     pres = level_df[pres_name]
                #     if np.sum(~np.isnan(pres))==0: # if PRES data is missing, might mean that only non-delayed mode is available, use "TEMP_ADJUSTED" instead
                #         pres_name = 'PRES'
                #         pres = level_df[pres_name]

                prof_meta = level_df['N_PROF']
                levels_meta = level_df['N_LEVELS']
                cycle_meta = level_df['CYCLE_NUMBER']
                juld_meta = [date_obj.strftime('%Y%m%d') if pd.notna(date_obj) else np.nan for date_obj in level_df['JULD']]
                lat_meta = level_df['LATITUDE']
                lon_meta = level_df['LONGITUDE']
                # temp_meta = level_df['TEMP_ADJUSTED_RO']
                # psal_meta = level_df['PSAL_ADJUSTED_RO']
                # doxy_meta = level_df['DOXY_ADJUSTED_RO']
                # pres_meta = level_df['PRES_ADJUSTED_RO']
                index_meta = level_df.index
                
                for i, var in enumerate(plot_vars,1):
                    
                    symbol = 'circle-dot'             
                    data = level_df[var]
                    if var=='TEMP_ADJUSTED' + '_' + suffix:
                        if np.sum(~np.isnan(data)) < np.sum(~np.isnan(pres)): # if temp data is missing, might mean that only non-delayed mode is available, use "TEMP_ADJUSTED" instead
                            var = 'TEMP_ADJUSTED'
                            data = level_df[var]
                            if np.sum(~np.isnan(data))==0: # check if TEMP_ADJUSTED is also missing and use "TEMP" instead
                                var = 'TEMP'
                                data = level_df[var]
                            symbol = 'x-thin' 

                    if var=='PSAL_ADJUSTED' + '_' + suffix:
                        if np.sum(~np.isnan(data)) < np.sum(~np.isnan(pres)): # if temp data is missing, might mean that only non-delayed mode is available, use "TEMP_ADJUSTED" instead
                            var = 'PSAL_ADJUSTED'
                            data = level_df[var]
                            if np.sum(~np.isnan(data))==0: # check if TEMP_ADJUSTED is also missing and use "TEMP" instead
                                var = 'PSAL'
                                data = level_df[var]
                            symbol = 'x-thin' 

                    var_meta = len(level_df.index) * [var]
                    n_meta = len(level_df.index) * [i]
                    # print(colors_ts[level_count])
                    # Plot or "add" trace to figure with scatter plot function establishing custom data and hover template as well as name of trace for legend
                    variable_fig.add_trace(go.Scattergl(x=time, y=data, name=f"Pressure {str(press_n)}", line={'color':colors_ts[level_count], 'width':1}, 
                                                        mode='markers', marker_symbol = symbol, marker_line_color=colors_ts[level_count], marker_line_width=2,  marker_size=5,
                                                        customdata=np.stack((
                        index_meta, prof_meta, levels_meta, cycle_meta, juld_meta, np.round(lat_meta,2), np.round(lon_meta,2), var_meta, n_meta, pres
                    ), axis=1),
                    hovertemplate='<br>'.join(['variable: %{y}',
                                            'juld: %{customdata[4]}',
        # #                                   'variable: %{customdata[7]}',
                                            'index: %{customdata[0]}',
                                            'profile: %{customdata[1]}',
                                            'level: %{customdata[2]}',
        # #                                   'cycle: %{customdata[3]}',
                                            
                                            'lat: %{customdata[5]}',
                                            'lon: %{customdata[6]}',
                                              'pressure: %{customdata[9]}'])
                                                ), row=i, col=1)
                    # Create axis titles and enable grid also reverse y axis
                    variable_fig.update_yaxes(title_text=var, showgrid=True, row=i, col=1, autorange=True)
                    if i==len(plot_vars):
                        variable_fig.update_xaxes(title_text='Date', showgrid=True, row=i, col=1, matches='x')
                    # variable_fig.update_layout(title_text = 'Levels ' + str(levels_dict[level_slider]) + " Pressures " + str(np.round(min(filtered_df['PRES_ADJUSTED_RO']),1)) + "-" + str(np.round(np.max(filtered_df['PRES_ADJUSTED_RO']),1)))
                variable_fig.update_layout(title_text = 'Levels ' + str(np.round(min(filtered_df['N_LEVELS']),1)) + "-" + str(np.round(np.max(filtered_df['N_LEVELS']),1))
                                                                                        + " Pressures " + str(np.round(min(filtered_df[pres_name]),1)) + "-" + str(np.round(np.max(filtered_df[pres_name]),1)) + 
                                                                                        '. If PSAL or TEMP have "x" symbols, non "RO" data is being used')

                level_count = level_count + 1  
                # Change marker size
                # variable_fig.update_traces(marker=dict(size=4))

            # add lines indicating the location of the profiles plotted above 
            for i, var in enumerate(plot_vars,1):
                variable_fig.add_vline(x=prof_dec_year_range[0],line_dash="dash")
                variable_fig.add_vline(x=prof_dec_year_range[1],line_dash="dash")

            return variable_fig
    
    
        # plotting elements
        dumdf = pd.DataFrame(columns = ["Float Number","Variable","N_PROF","Date","N_LEVELS",'Pressure (dbar)','Deletion reason'])
        # if verbose:
        #     print('Populating table with initial points from density inversion')
        # if len(points_dens_inversion)>0: # Add the detected density inversions to the table
        #     dumdf = pd.concat( [
        #         pd.DataFrame({
        #             "Float Number":[str(each.WMO_ID.values) for each in points_dens_inversion],
        #             "Variable":[var]*len(points_dens_inversion),
        #             "N_PROF":prof_inv,
        #             "Date":[each.JULD.values for each in points_dens_inversion],
        #             "N_LEVELS":levels_inv,
        #             'Pressure (dbar)':[float(each.PRES_ADJUSTED_RO.values) for each in points_dens_inversion],
        #             'Deletion reason':['Density inversion']*len(points_dens_inversion)
        #         })
        #         for var in ['PSAL_ADJUSTED_RO', 'TEMP_ADJUSTED_RO']], ignore_index=True )
        
        dtable = create_table(dumdf)
    
        profile_fig, variable_plot_position, prof_dec_year_range, profile_lat_lon = create_figure(df, 0)
        # print(profile_lat_lon)
        map_fig = create_map(argo_n, profile_lat_lon) # moving map figure here to add profile locations. 

        # print(prof_dec_year_range)
        ts_fig = create_figure_ts(df, 0, prof_dec_year_range)
    
        external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
        app = Dash(__name__, external_stylesheets=external_stylesheets)
        app.layout = html.Div([
                    dcc.Store(id='file_index', data=file),
                    # Header
                    html.H1([
                        "Float Outlier Removal"
                    ]),
                    # Place map figure in this row
                    dbc.Row(
                        [
                            dbc.Col(html.Div([dcc.Graph(id='graph_map', figure=map_fig)]))
                        ]
                    ),
                    # Row for download button
                    dbc.Row(
                        [
                            dbc.Col(html.Div(
                                [
                                    html.Button("Download CSV", id="btn_csv"),
                                    dcc.Download(id="download-dataframe-csv"),
                                ]
                            ))
                        ]
                    ),
                    # Place table in this row
                    dbc.Row(
                        [
                            dbc.Col([html.Div([dtable]), html.Div([dcc.Store(id='download_table'), dcc.Store(id='download_flag', data=0)]),], width = 12)
                        ]
                    ),
                    # Place outlier profile figure in this row, establish buttons and slider as well
                    dbc.Row(
                        [
                            dbc.Col(html.Div([
                                        html.Button('Remove Deepest measurements profiles shown:', id='bottom_remove'),
                                        html.Button('Remove deepest measurements for all profiles for variable:', id='bottom_remove_all'),
                                        dcc.RadioItems(id = 'bottom_remove_choice',
                                        options = [dict(label = 'Nitrate', value = 1),
                                                 dict(label = 'Oxygen', value = 2), 
                                                 dict(label = 'pH', value = 3)],
                                                 value=0),
                                        html.Button('Flag as outliers', id='delete'),
                                        html.Button('Clear Selection', id='clear'),
                                        dcc.RadioItems(id = 'y_scale',
                                        options = [dict(label = 'y-axis max: 10m', value = 1),
                                                 dict(label = 'y-axis max: 50m', value = 2), 
                                                 dict(label = 'y-axis max: 200m', value = 4), 
                                                 dict(label = 'y-axis max: 500m', value = 5), 
                                                 dict(label = 'y-axis max: Full depth', value = 0),
                                                 dict(label = 'y-axis max: Bottom 200m', value = -1),
                                                 dict(label = 'y-axis range: Bottom 500m', value = -2)],
                                                 value=0),
                                        
                                        dcc.Graph(id='graph_with_slider', figure=profile_fig),
                                        dcc.Slider(
                                            step=1, # for some reason changing to 1 allows arrow key usage instead of the len of marker_dict
                                            value=0,
                                            marks=marker_dict,
                                            id='profile_slider'
                                        ),
                                        dcc.Store(id='prev_profile_slider', data=0),
                                        #html.Div('selected:'),
                                        html.Div(id='selected_points', style={'display': 'none'}),
                                        #html.Div('deleted:'),
                                        html.Div(id='deleted_points', style={'display': 'none'})
                                    ])
                                    , width = 12)
                        ]
                    ),
                    # Place time series figure in this row: 
                    dbc.Row(
                        [
                            dbc.Col(html.Div([
                                        # html.Button('Flag as outliers', id='delete_ts'),
                                        html.Button('Clear Selection', id='clear_ts'),
                                        dcc.Graph(id='graph_with_slider_ts', figure=ts_fig, style={'height' :'800px'}),
                                        dcc.Slider(      
                                            step=1, # for some reason changing to 1 allows arrow key usage instead of the len of marker_dict
                                            value=0,
                                            marks=pressures_dict,
                                            id='level_slider'
                                        ),
                                        dcc.Store(id='prev_level_slider', data=0),
                                        dcc.Store(id='prof_dec_year_range', data=prof_dec_year_range),
                                        #html.Div('selected:'),
                                        html.Div(id='selected_points_ts', style={'display': 'none'}),
                                        #html.Div('deleted:'),
                                        html.Div(id='deleted_points_ts', style={'display': 'none'})
                                    ])
                                    , width = 12
                                    , style={"height": "800px"}
                                    )
                        ]
                    ),
        ]
        )
        # app.callback updates specific parts of the app
        # Input triggers callback which will perform the function below and return the updated output, state hands off data to the function
        #Handle deleting points
        
        
        # If Flag as outlier (delete) button is clicked then perform this callback # combined time series and profile into one button - simpler, no reason not to (I hope)
        # Inputs: 
        # - number of clicks for delete button (not used)
        # - selected points (state)
        # - deleted points (state)
        # - selected_points_ts (state)
        # - deleted_points_ts (state)
        # Loads "delete_points" json if it exists, otherwise creates an empty list
        # Adds selected points to deleted_points, then clears selection
        # Returns new list of deleted points 
        @app.callback([Output('deleted_points', 'children'), 
                       Output('deleted_points_ts', 'children')],
                    [Input('delete', 'n_clicks')],
                    [State('selected_points', 'children'),
                    State('deleted_points', 'children'), 
                    State('selected_points_ts', 'children'),
                    State('deleted_points_ts', 'children')])
        def delete_points(n_clicks, selected_points, delete_points, selected_points_ts, delete_points_ts):
            # If there are selected points, load them, otherwise create empty list
            if selected_points:
                selected_points = json.loads(selected_points)
            else:
                selected_points = []
            # print(selected_points)
            # If there are deleted points, load them, otherwise create empty list
            if delete_points:
                deleted_points = json.loads(delete_points)
            else:
                deleted_points = []
            # Create list of selected points indexes in the float dataframe
            new_deleted = selected_points
            # Expand deleted points list to include newly deleted points
            deleted_points.extend(new_deleted)
            selected_points = [] # empty for the next round of selection
            
            # same as above, but for time series
            if selected_points_ts:
                selected_points_ts = json.loads(selected_points_ts)
            else:
                selected_points_ts = []
            # If there are deleted points, load them, otherwise create empty list
            if delete_points_ts:
                deleted_points_ts = json.loads(delete_points_ts)
            else:
                deleted_points_ts = []
        
            # Create list of selected points indexes in the float dataframe
            new_deleted_ts = selected_points_ts
            # Expand deleted points list to include newly deleted points
            deleted_points_ts.extend(new_deleted_ts)
            selected_points_ts = []
            # print(deleted_points)
            # print(type(deleted_points))
            # Return deleted points in json object
            return json.dumps(deleted_points), json.dumps(deleted_points_ts)
        
        
        # Handle selecting points 
        # If the clear button or the graph performs any selection this perform this callback
        @app.callback(Output('selected_points', 'children', allow_duplicate=True),
                    [Input('graph_with_slider', 'clickData'),
                     Input('graph_with_slider', 'selectedData'),
                        Input('deleted_points', 'children'),
                        Input('clear', 'n_clicks')],
                    [State('selected_points', 'children')],
                    prevent_initial_call=True)
        def select_point(clickData, selectedData, deleted_points,  clear, selected_points):
            # Load in callback events into a list
            ctx = dash.callback_context
            #print(ctx.triggered)
            ids = [c['prop_id'] for c in ctx.triggered]
        
            # If selected points exists, load it, otherwise make empty results list
            if selected_points:
                results = json.loads(selected_points)
            else:
                results = []
        
            # If clickData callback event is in list then append the clicked points to results
            if 'graph_with_slider.clickData' in ids:
                if clickData:
                    for p in clickData['points']:
                        if p not in results:
                            results.append(p)
            # If selectData callback event is in list then append the clicked points to results
            if 'graph_with_slider.selectedData' in ids:
                if selectedData:
                    for p in selectedData['points']:
                        if p not in results:
                            results.append(p)
            # If there are deleted points present in ids or clear button is clicked then clear selected points
            if 'deleted_points.children' in ids or  'clear.n_clicks' in ids:
                results = []
                selected_points = []
            # Return selected points list
            #print(results)
            results = json.dumps(results)
            return results
        
        # Same for time series
        @app.callback(Output('selected_points_ts', 'children'),
                    [Input('graph_with_slider_ts', 'clickData'),
                     Input('graph_with_slider_ts', 'selectedData'),
                        Input('deleted_points_ts', 'children'),
                        Input('clear_ts', 'n_clicks')],
                    [State('selected_points_ts', 'children')])
        
        def select_point_ts(clickData, selectedData, deleted_points_ts, clear_ts, selected_points_ts):
            # Load in callback events into a list
            ctx = dash.callback_context
            #print(ctx.triggered)
            ids = [c['prop_id'] for c in ctx.triggered]
        
            # If selected points exists, load it, otherwise make empty results list
            if selected_points_ts:
                results = json.loads(selected_points_ts)
            else:
                results = []
        
            # If clickData callback event is in list then append the clicked points to results
            if 'graph_with_slider_ts.clickData' in ids:
                if clickData:
                    for p in clickData['points']:
                        if p not in results:
                            results.append(p)
            # If selectData callback event is in list then append the clicked points to results
            if 'graph_with_slider_ts.selectedData' in ids:
                if selectedData:
                    for p in selectedData['points']:
                        if p not in results:
                            results.append(p)
        
            # If there are deleted points present in ids or clear button is clicked then clear selected points
            if 'deleted_points_ts.children' in ids or  'clear_ts.n_clicks' in ids:
                results = []
                selected_points_ts = [] # empty for the next round of selection

            # Return selected points list
            results = json.dumps(results)
            return results
        
        # Handle updating figure based on selection and deletion of points
        # If profile slider is altered, perform this callback. If selected and deleted points are changed also perform this callback
        @app.callback(
            Output('graph_with_slider', 'figure'), Output('prev_profile_slider', 'data'), Output('prof_dec_year_range', 'data'),
            Input('file_index', 'data'),
            Input('profile_slider', 'value'),
                    [Input('selected_points', 'children'),
                   Input('y_scale', 'value')],
                    [State('deleted_points', 'children'),
                    State('prev_profile_slider', 'data'), ])
        def update_figure(file_idx, profile_slider, selected_points, y_scale, deleted_points_state, prev_profile_slider):

            # Load and serialize dataset
            # argo_n, df, marker_dict, level_dict = load_dataset(output_dir, file_idx, plot_vars_all)
            # print(deleted_points_state)
            # Load deleted points
            deleted_points = json.loads(deleted_points_state) if deleted_points_state else []
            # print(deleted_points)
            deleted_indexes = []
            deleted_variables = []
            for dp in deleted_points:
                deleted_indexes.append(int(float(dp['customdata'][0])))
                deleted_variables.append(dp['customdata'][7])
            # print(deleted_variables)
            # Create new figure based on profile slider state and deleted points
            f, junk, prof_dec_year_range, profile_lat_lon = create_figure(df, profile_slider, deleted_indexes, deleted_variables, y_scale)
            # Create a red point on every selected point
            selected_points = json.loads(selected_points) if selected_points else []
            if (len(selected_points)!=0) & (profile_slider == prev_profile_slider):
                highlight_plot_vars = {}
                for p in selected_points:
                    if p['customdata'][7] not in highlight_plot_vars:
                        highlight_plot_vars[p['customdata'][7]] = [p]
                    else:
                        highlight_plot_vars[p['customdata'][7]].append(p)
                for selected_variable in  highlight_plot_vars.keys():

                    v = variable_plot_position[selected_variable]
                    #print(highlight_plot_vars.keys())
                    # Add variables that are available
                    
                    # if selected_variable == 'TEMP_ADJUSTED_RO': v = 1
                    # elif selected_variable == 'PSAL_ADJUSTED_RO': v = 2
                    # elif selected_variable == 'DOXY_ADJUSTED_RO': v = 3
                    # elif selected_variable == 'NITRATE_ADJUSTED_RO': v = 4
                    # else:
                    #     v = 0
                    points = highlight_plot_vars[selected_variable]
                    f.add_trace(
                        go.Scattergl(
                            mode='markers',
                            x=[p['x'] for p in points],
                            y=[p['y'] for p in points],
                            marker=dict(
                                color='red',
                                size=5,
                                line=dict(
                                    color='red',
                                    width=2
                                )
                            ),
                            showlegend=False,   
                        ),
                    row=1, col=v
                    )
       
        
            # Keep same zoom characteristics with figure updates
            if profile_slider == prev_profile_slider : 
                f.update_layout({"uirevision": "foo"}, overwrite=True)
            # Include both event and selection clicks
            f.update_layout(clickmode='event+select')
            # Artificially make transition time, makes it seem 'smoother'
            f.update_layout(transition_duration=500)
            return f, profile_slider, prof_dec_year_range
        
        # Same for level slider
        @app.callback(
            Output('graph_with_slider_ts', 'figure'), Output('prev_level_slider', 'data'),
            Input('file_index', 'data'),
            Input('level_slider', 'value'),
                    [Input('selected_points_ts', 'children'),
                    Input('prof_dec_year_range', 'data')],
                    [State('deleted_points_ts', 'children'),
                     State('prev_level_slider', 'data')])
        def update_figure_ts(file_idx, level_slider, selected_points_ts, prof_dec_year_range, deleted_points_ts_state, prev_level_slider):
            # Load and serialize dataset
            # argo_n, df, marker_dict, level_dict = load_dataset(output_dir, file_idx, plot_vars_all)
            
            # Load deleted points
            deleted_points_ts = json.loads(deleted_points_ts_state) if deleted_points_ts_state else []
            deleted_indexes = []
            deleted_variables = []
            for dp in deleted_points_ts:
                deleted_indexes.append(int(float(dp['customdata'][0])))
                deleted_variables.append(dp['customdata'][7])
            f = create_figure_ts(df, level_slider, prof_dec_year_range, deleted_indexes, deleted_variables)
            # Create a red point on every selected point
            selected_points_ts = json.loads(selected_points_ts) if selected_points_ts else []
            if (len(selected_points_ts)!=0) & (level_slider == prev_level_slider):
                
                f.add_trace(
                    go.Scattergl(
                        mode='markers',
                        x=[p['x'] for p in selected_points_ts],
                        y=[p['y'] for p in selected_points_ts],
                        marker=dict(
                            color='red',
                            size=5,
                            line=dict(
                                color='red',
                                width=2
                            )
                        ),
                        showlegend=False
                    ),
                    col = 1,
                    row = int(selected_points_ts[-1]['customdata'][8])
                )
            
            # Keep same zoom characteristics with figure updates
            # only if did not move on the slider
            if level_slider == prev_level_slider : 
                f.update_layout({"uirevision": "foo"}, overwrite=True)
                
            # Include both event and selection clicks
            f.update_layout(clickmode='event+select')
            # Artificially make transition time, makes it seem 'smoother'
            f.update_layout(transition_duration=500)
            
            return f, level_slider
        
        # Update table based on what points have been removed
        # Inputs: 
        # - file_idx 
        # - selected_points 
        # - selected_points_ts 
        # - deleted_points 
        # - deleted_points_ts 
        # - csv_button_n_clicks 
        # - deleted_points 
        # - deleted_points_ts
        # - download_flag 
        # Loads deleted points from profile and time series if they exist, then combines them
        # creates an outlier_df using the deleted data index / vriables, plus data from "df"
        # Output:
        # - "Outlier_Table" (data1) which is a list-like dictionary of outlier dataframe
        # - "download_table" (json_outlier_df) which is a json of the outlier dataframe 
        @app.callback(Output('Outlier_Table','data',allow_duplicate=True),
                      Output('download_table', 'data'),
                    [Input('file_index', 'data'),
                    Input('selected_points', 'children'), Input('selected_points_ts', 'children'),
                    Input('deleted_points', 'children'), Input('deleted_points_ts', 'children'), Input('btn_csv','n_clicks')
                    ],
                    [State('deleted_points', 'children'), State('deleted_points_ts', 'children'), State('download_flag', 'data')
                    ],prevent_initial_call=True,)
        def update_table(file_idx, selected_points, selected_points_ts, deleted_points_input, deleted_points_input_ts, n_clicks, deleted_points_state, deleted_points_state_ts, download_flag):

            # If initial display, display density inversions - note that this currently will not display until more points are selected 
            # if n_clicks == None:
            #     return dumdf.to_dict('records'), dumdf.to_json(orient='table')
                
            # Load and serialize dataset
            # argo_n, df, marker_dict, levels_dict = load_dataset(output_dir, file_idx, plot_vars_all)

            # Load deleted points
            deleted_points = json.loads(deleted_points_state) if deleted_points_state else []
            deleted_points_ts = json.loads(deleted_points_state_ts) if deleted_points_state_ts else []
            combined_deleted_points = deleted_points + deleted_points_ts
            source_points = ['Profile'] * len(deleted_points) + ['Time series'] * len(deleted_points_ts)
            
            # If there are deleted points, make a new outlier df, populating the dataframe otherwise create an empty dataframe
            if combined_deleted_points:
                deleted_indexes = []
                deleted_variables = []
                # print(combined_deleted_points)
                for dp in combined_deleted_points: # step through deleted points, appending indexes and variable names
                    deleted_indexes.append(int(float(dp['customdata'][0])))
                    deleted_variables.append(dp['customdata'][7])
                outlier_data = df.loc[deleted_indexes]   # gets "outlier_data" from "df", using the indexes supplied by the "customdata"
                # print(outlier_data)
                pres_name = 'PRES_ADJUSTED' + '_' + suffix
                if np.sum(~np.isnan(df[pres_name]))==0: 
                    pres_name = 'PRES_ADJUSTED'
                    if np.sum(~np.isnan(df[pres_name]))==0: 
                        pres_name = 'PRES'

                outlier_df = pd.DataFrame({                 # creates a new "outlier_df" w wmo, variable name, profile, date, level, pressure, deletion reason 
                    "Float Number":len(outlier_data["N_PROF"])*[str(argo_n.WMO_ID.values)],
                    "Variable":deleted_variables,
                    "N_PROF":outlier_data["N_PROF"],
                    "Date":outlier_data["JULD"],
                    "N_LEVELS":outlier_data["N_LEVELS"],
                    'Pressure (dbar)':outlier_data[pres_name],
                    'Deletion reason':source_points
                })   
                # if download_flag == 0: # if download flag hasn't been clicked, concatenates the dummy dataframe and the new outlier dataframe 
                    # print(dumdf)
                    # print(outlier_df)
                    # outlier_df = pd.concat([dumdf, outlier_df]) # commenting this out as I no longer think it is necessary w/ no initally removed points
                data1 = outlier_df.to_dict('records')   # creates "data1" which is a list-like dictionary of the outlier dataframe
                json_outlier_df = outlier_df.to_json(orient='table') # also saves the outlier_df to a json
                return data1, json_outlier_df  # returns both dictionary and json of outlier_df 
            else:
                outlier_df = pd.DataFrame(columns = ["Float Number","Variable","N_PROF","Date","N_LEVELS",'Pressure (dbar)','Deletion reason'])
                if download_flag == 0:
                    outlier_df = dumdf
                data1 = outlier_df.to_dict('records')
                json_outlier_df = outlier_df.to_json(orient='table')
                return data1, json_outlier_df
        
        
        # Called when the download_button is clicked 
        # Inputs:
        # - number of clicks for csv button (not used)
        # - download_table - comes from "update_table" as "df_for_saving"
        # - file_index as "file_index"
        # Loads df_for_saving json
        # saves to a csv, either writing a new one or appending, depending on whether a csv with that name already exists
        # resets the outlier_df, data1, json_outlier_df, resets some of the deleted points, but not all. 
        # Outputs:
        # - download-dataframe-csv, value of "none" - still not sure what this is/does
        # - Outlier_Table - data1, should be empty
        # - download_table - json_outlier_df, should be empty
        # - deleted_points - deleted_points_ts - should be empty now 
        # - deleted_points_ts - deleted_points_ts - should be empty now 
        # - selected_points - results, also should be empty
        # - selected_points_ts - results, also should be empty

        # - download_flag, now set to 1

        @app.callback(
            [Output("download-dataframe-csv", "data"), 
            Output('Outlier_Table','data',allow_duplicate=True),
            Output('download_table', 'data',allow_duplicate=True),
            Output('deleted_points', 'children', allow_duplicate=True),
            Output('deleted_points_ts', 'children', allow_duplicate=True),
            Output('selected_points', 'children', allow_duplicate=True),
            Output('selected_points_ts', 'children', allow_duplicate=True),
            Output('download_flag', 'data')],
            Input("btn_csv", "n_clicks"),
            [State('download_table', 'data'),
            State('file_index', 'data')],
            prevent_initial_call=True,
        )
        def func(n_clicks, df_for_saving, file_index):
            # read_outlier_df = pd.read_json(df, orient='table')
            read_outlier_df = pd.read_json(StringIO(df_for_saving), orient='table')
            # print(read_outlier_df)
            current_time = datetime.datetime.now()
            current_time_str = str(current_time.year) + '_' + str(current_time.month) + '_' + str(current_time.day) + '_' + str(current_time.hour)
            
            savename = '../outlier_files/outliers_' + file_index[0:-3] + '_' + researcher + '_' + current_time_str + '.csv'
        
            if os.path.exists(savename):
                read_outlier_df.to_csv(savename, mode='a', index=False, header=False)
            else:
                read_outlier_df.to_csv(savename, mode='w', index=False, header=True)
        
            # reset outlier table
            outlier_df = pd.DataFrame(columns = ["Float Number","Variable","N_PROF","Date","N_LEVELS",'Pressure (dbar)','Deletion reason'])
            data1 = outlier_df.to_dict('records')
            json_outlier_df = outlier_df.to_json(orient='table')
        
            # reset list of deleted and selected points
            deleted_points_ts = []
            results = []
            results = json.dumps(results)

            download_flag = 1

            return None, data1, json_outlier_df, json.dumps(deleted_points_ts), json.dumps(deleted_points_ts), results, results, download_flag
        
        @app.callback([Output('deleted_points', 'children', allow_duplicate=True)],
                    [Input('bottom_remove', 'n_clicks')],
                    [State('profile_slider', 'value'),
                     State('bottom_remove_choice', 'value'),
                     State('selected_points', 'children'),
                    State('deleted_points', 'children')],
                    prevent_initial_call=True)
        def remove_bottom_points(n_clicks, profile_slider, bottom_remove_choice, selected_points, delete_points, ):
                       # If there are selected points, load them, otherwise create empty list
            if selected_points:
                selected_points = json.loads(selected_points)
            else:
                selected_points = []
            # If there are deleted points, load them, otherwise create empty list
            if delete_points:
                deleted_points = json.loads(delete_points)
            else:
                deleted_points = []

            dfs = df
            filtered_df = dfs[dfs["Profile Range"] == marker_dict[profile_slider]['label']]
            # print(filtered_df.keys())
            # print(filtered_df)
            match_key = ''
            match=0
            if bottom_remove_choice==1:
                match_var = 'NITRATE_ADJUSTED' + '_' + suffix
            elif bottom_remove_choice==2:
                match_var = 'DOXY_ADJUSTED' + '_' + suffix
            elif bottom_remove_choice==3:
                match_var = 'PH_IN_SITU_TOTAL_ADJUSTED' + '_' + suffix

            plot_vars = plot_vars_all.copy()
            for var in plot_vars_all:
                if not var in argo_n.keys():
                    plot_vars.remove(var)
                    
            for i_n, var_n in enumerate(plot_vars,1):
                # find "i"
                if match_var==var_n:
                    i = i_n
                    var = match_var
                    # print(i)
                    # print(var)
                    break
           
            # find each deepest point that is still valid
            for prof, data in filtered_df.groupby("N_PROF"):
                # Select data from one profile
                prof_df = filtered_df[filtered_df["N_PROF"] == prof]

                data = prof_df[var]
                # print(data.last_valid_index())
                bottom_index = data.last_valid_index()
                # print(bottom_index)
                if not bottom_index:
                    continue
                # print(prof_df)
                prof_meta = prof_df['N_PROF']
                prof_meta = prof_meta[bottom_index]
                # print(prof_meta)
                levels_meta = prof_df['N_LEVELS']
                levels_meta = levels_meta[bottom_index]
                
                cycle_meta = prof_df['CYCLE_NUMBER']
                cycle_meta = cycle_meta[bottom_index]

                index_meta = bottom_index
                # index_meta = index_meta[bottom_index]

                var_meta = var
                
                n_meta = i
                
                juld_meta = [date_obj.strftime('%Y%m%d') for date_obj in prof_df['JULD']]
                juld_meta = juld_meta[0]
                
                lat_meta = prof_df['LATITUDE']
                # print(lat_meta)
                lat_meta = lat_meta[bottom_index]
                
                lon_meta = prof_df['LONGITUDE']
                lon_meta = lon_meta[bottom_index]
                
                pres = prof_df['PRES_ADJUSTED' + '_' + suffix]
                if np.sum(~np.isnan(pres))==0: # if PRES data is missing, might mean that only non-delayed mode is available, use "TEMP_ADJUSTED" instead
                    pres = prof_df['PRES_ADJUSTED']
                    if np.sum(~np.isnan(pres))==0: # if PRES data is missing, might mean that only non-delayed mode is available, use "TEMP_ADJUSTED" instead
                        pres = prof_df['PRES']
                
                # print(data.last_valid_index())
                # bottom_index = data.last_valid_index()
                # print(pres[bottom_index])

                # print(data[bottom_index])
                
                # customdata=np.stack((
                #     index_meta, prof_meta, levels_meta, cycle_meta, 
                #     juld_meta, np.round(lat_meta,2), np.round(lon_meta,2), var_meta, n_meta, str(pres[bottom_index])))
                # print(customdata)
                selected_points = [{'curveNumber': '0', 'pointNumber': str(levels_meta), 'pointIndex': str(levels_meta), 'x': str(data[bottom_index]), 'y': str(pres[bottom_index]), 'customdata': [str(index_meta), str(prof_meta), str(levels_meta), str(cycle_meta), 
                    juld_meta, str(np.round(lat_meta,2)), str(np.round(lon_meta,2)), var_meta, str(n_meta), str(pres[bottom_index])]}] 
                
                new_deleted = selected_points
            # Expand deleted points list to include newly deleted points
                deleted_points.extend(new_deleted)
            # print(deleted_points)
            # print(type(deleted_points))
            
            selected_points = [] # empty for the next round of selection
            

            return json.dumps(deleted_points),

        @app.callback([Output('deleted_points', 'children', allow_duplicate=True)],
                    [Input('bottom_remove_all', 'n_clicks')],
                    [State('profile_slider', 'value'),
                     State('bottom_remove_choice', 'value'),
                     State('selected_points', 'children'),
                    State('deleted_points', 'children')],
                    prevent_initial_call=True)
        def remove_all_bottom_points(n_clicks, profile_slider, bottom_remove_choice, selected_points, delete_points, ):
                       # If there are selected points, load them, otherwise create empty list
            if selected_points:
                selected_points = json.loads(selected_points)
            else:
                selected_points = []
            # If there are deleted points, load them, otherwise create empty list
            if delete_points:
                deleted_points = json.loads(delete_points)
            else:
                deleted_points = []

            dfs = df
            filtered_df = dfs # instead of selecting profiles, loop through all profiles in the float and remove bottom points for the variable selected
            match_key = ''
            match=0
            if bottom_remove_choice==1:
                match_var = 'NITRATE_ADJUSTED' + '_' + suffix
            elif bottom_remove_choice==2:
                match_var = 'DOXY_ADJUSTED' + '_' + suffix
            elif bottom_remove_choice==3:
                match_var = 'PH_IN_SITU_TOTAL_ADJUSTED' + '_' + suffix

            plot_vars = plot_vars_all.copy()
            for var in plot_vars_all:
                if not var in argo_n.keys():
                    plot_vars.remove(var)
                    
            for i_n, var_n in enumerate(plot_vars,1):
                # find "i"
                if match_var==var_n:
                    i = i_n
                    var = match_var
                    # print(i)
                    # print(var)
                    break
           
            # find each deepest point that is still valid
            for prof, data in filtered_df.groupby("N_PROF"):
                # Select data from one profile
                prof_df = filtered_df[filtered_df["N_PROF"] == prof]

                data = prof_df[var]
                # print(data.last_valid_index())
                bottom_index = data.last_valid_index()
                # print(bottom_index)
                if not bottom_index:
                    continue
                # print(prof_df)
                prof_meta = prof_df['N_PROF']
                prof_meta = prof_meta[bottom_index]
                # print(prof_meta)
                levels_meta = prof_df['N_LEVELS']
                levels_meta = levels_meta[bottom_index]
                
                cycle_meta = prof_df['CYCLE_NUMBER']
                cycle_meta = cycle_meta[bottom_index]

                index_meta = bottom_index
                # index_meta = index_meta[bottom_index]

                var_meta = var
                
                n_meta = i
                
                juld_meta = [date_obj.strftime('%Y%m%d') for date_obj in prof_df['JULD']]
                juld_meta = juld_meta[0]
                
                lat_meta = prof_df['LATITUDE']
                # print(lat_meta)
                lat_meta = lat_meta[bottom_index]
                
                lon_meta = prof_df['LONGITUDE']
                lon_meta = lon_meta[bottom_index]
                
                pres = prof_df['PRES_ADJUSTED' + '_' + suffix]
                if np.sum(~np.isnan(pres))==0: # if PRES data is missing, might mean that only non-delayed mode is available, use "TEMP_ADJUSTED" instead
                    pres = prof_df['PRES_ADJUSTED']
                    if np.sum(~np.isnan(pres))==0: # if PRES data is missing, might mean that only non-delayed mode is available, use "TEMP_ADJUSTED" instead
                        pres = prof_df['PRES']
                
                # print(data.last_valid_index())
                # bottom_index = data.last_valid_index()
                # print(pres[bottom_index])

                # print(data[bottom_index])
                
                # customdata=np.stack((
                #     index_meta, prof_meta, levels_meta, cycle_meta, 
                #     juld_meta, np.round(lat_meta,2), np.round(lon_meta,2), var_meta, n_meta, str(pres[bottom_index])))
                # print(customdata)
                selected_points = [{'curveNumber': '0', 'pointNumber': str(levels_meta), 'pointIndex': str(levels_meta), 'x': str(data[bottom_index]), 'y': str(pres[bottom_index]), 'customdata': [str(index_meta), str(prof_meta), str(levels_meta), str(cycle_meta), 
                    juld_meta, str(np.round(lat_meta,2)), str(np.round(lon_meta,2)), var_meta, str(n_meta), str(pres[bottom_index])]}] 
                
                new_deleted = selected_points
            # Expand deleted points list to include newly deleted points
                deleted_points.extend(new_deleted)
            # print(deleted_points)
            # print(type(deleted_points))
            
            selected_points = [] # empty for the next round of selection
            

            return json.dumps(deleted_points),
            
        # if __name__ == "__main__":
        app.run(debug=True, port=port_num, jupyter_mode="external")
        # app.run_server(debug=True, port=port_num)

    outlier_removal_plot(output_dir, file, argo_n, df, marker_dict, pressures_dict, researcher)
    
    return success

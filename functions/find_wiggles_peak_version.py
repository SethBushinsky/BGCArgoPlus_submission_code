import numpy as np
import xarray as xr
from scipy.signal import find_peaks

def peak_detect_wiggles(output_dir, argo_n, var):
    """
    ID profiles with + and - wiggles by finding the 
    # of + and - peaks per profile. ID peaks using the 
    scipy find peaks function using the peak prominence method
    "https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html"
    If the number of peaks 
    in the variable with no corresponding peak in temperature
    exceeds 4, then the cycle number is 
    returned in wiggle_cycle_list"""
    wiggle_cycle_list = np.array([])
    # file_n = output_dir + file
    # argo_n = xr.open_dataset(file_n)

    if var not in argo_n.keys():
        return wiggle_cycle_list
    min_plot_prof = 0
    max_plot_prof = len(argo_n.N_PROF)
    
    for p in np.arange(min_plot_prof, max_plot_prof):
        var_p = argo_n[var][p,:].values
        press_p = argo_n.PRES_ADJUSTED[p,:].values
        temp_p = argo_n.TEMP_ADJUSTED[p,:].values
        var_nonan = var_p[np.logical_and(~np.isnan(var_p), ~np.isnan(press_p))]
        nlevels_nonan=argo_n.N_LEVELS[np.logical_and(~np.isnan(var_p), ~np.isnan(press_p))]
        if len(var_nonan)==0:
            continue
        press_nonan = press_p[np.logical_and(~np.isnan(var_p), ~np.isnan(press_p))]
        tempc_nonan = temp_p[np.logical_and(~np.isnan(var_p), ~np.isnan(press_p))]
        #find max and min temperature peaks
        Temp_max_peaks = find_peaks(tempc_nonan, prominence =1)
        T_max_peak_positions = nlevels_nonan[Temp_max_peaks[0]] 
        invert_temp_p = tempc_nonan*-1
        Temp_min_peaks = find_peaks(invert_temp_p, prominence=1)
        T_min_peak_positions = nlevels_nonan[Temp_min_peaks[0]]
        
        #Find variable peaks  
        #picked prominence thresholds that are around or greater than sensor accuracy, seem to get peaks 
        if var.__contains__('NITRATE'): thresh=1.5
        elif var.__contains__('DOXY'): thresh=10
        elif var.__contains__('PH_IN_SITU_TOTAL'): thresh=0.01
        peaks = find_peaks(var_nonan, prominence =thresh, width=(0,3))
        maxima_peak_values=var_nonan[peaks[0]]
        max_peak_pos = nlevels_nonan[peaks[0]] 
        invert_var = var_nonan*-1
        minima = find_peaks(invert_var, prominence=thresh, width=(0,3))
        min_peak_pos = nlevels_nonan[minima[0]]
        minima_peak_values = invert_var[minima[0]] 
        
        #position of all peaks
        all_var_peaks=np.append(max_peak_pos, min_peak_pos)
        all_temp_peaks=np.append(T_max_peak_positions, T_max_peak_positions)
        remove_temp_peaks = [x for x in all_var_peaks if x not in all_temp_peaks and (x - 1) not in all_temp_peaks and (x + 1) not in all_temp_peaks]
        
        if len(remove_temp_peaks)>=4:
            wiggle_cycle_list = np.append(wiggle_cycle_list, argo_n['CYCLE_NUMBER'][p].values)
    return wiggle_cycle_list